%        $Id: $
%      usage: trackpos
%         by: Josh Ryu
%       date: 05/15/2019
%    purpose: 

function myscreen = trackpos_hierarchical(varargin)
 

%getArgs(varargin,{'subjectID=s999','centerX=10','centerY=0','diameter=16'}); getArgs(varargin,{'subjectID=-1'});
% set up screen
if isempty(mglGetSID)
    myscreen.subjectID  = -1;
else
    myscreen.subjectID  = mglGetSID;
end

%myscreen.displayName = 'joshipad2'; %
%myscreen.displayName = 'windowed';
% myscreen.screenWidth = 860; myscreen.screenHeight = 600; 
% myscreen.hideCursor = 1;

myscreen                = initScreen(myscreen);
% mglSetGammaTable(0,1,1,0,1,1,0,1,1)
% t = mglGetGammaTable; figure; hold on; plot(t.redTable,'r'); plot(t.blueTable,'b'); plot(t.greenTable,'g');


%% parameters
global stimulus; stimulus = struct;

% Experimenter parameters
experimenter = struct();
experimenter.noeye           = 1; % 1 if no eyetracking (mouse for eye); 0 if there is eye tracking `
experimenter.showmouse       = 1; 
experimenter.grabframe       = 0; % grab frame
experimenter.fixateCenter    = 1;
experimenter.phasescrambleOn = 1; % consider different noise? (i.e. pink noise? https://www.mathworks.com/help/audio/ref/pinknoise.html)

experimenter.precompute_path = '~/proj/grustim/trackpos_hierarchical/testrun/'; %required if using precompute
% 1: load precomputed stimulus (background and stim)
% 2: load precomputed background and stimulus position only (during background trials)
experimenter.precompute      = 0; 
experimenter.precompute_gen  = 0; %1; % generate and save stimulus
experimenter.downsample_spatRes  = 2; % downsample spatially
if ~exist(experimenter.precompute_path), mkdir(experimenter.precompute_path);,end

%% Basic task design (todo: move into configureExperiment function?)
% S1: Stimulus (30s) 
% S2: Fixation (3s)
design_block = struct();
design_block.time_stim = 30;
design_block.time_fix = 3;
design_block.nTrials_train = 1; % learning period
design_block.nTrials_track = 1; % testing period; with noise

task{1}{1}.segmin = [design_block.time_stim, design_block.time_fix]; % fixation for shorter bc of the segment start takes time.
task{1}{1}.segmax = [design_block.time_stim, design_block.time_fix]; 
task{1}{1}.numTrials = design_block.nTrials_train;
task{1}{1}.getResponse = [1 0]; %segment to get response.
task{1}{1}.waitForBacktick = 0; %wait for backtick before starting each trial 

% task parameters for adaptation conditions
if experimenter.phasescrambleOn == 1
    task{1}{1}.parameter.phasescrambleOn    = 1;
    task{1}{1}.parameter.backLum            = 96; % todo: check if it works with colored noise
    task{1}{1}.parameter.noiseLum           = 32; % todo: this is not changing right now . fix this.
else 
    task{1}{1}.parameter.backLum = 32;  % background luminance; units: fraction of full luminance 
end
task{1}{1}.parameter.stimLum = 255 - task{1}{1}.parameter.backLum;  % stimulus luminance (out of 255)

% calculated parameters
task{1}{1}.randVars.calculated.randomSeed = 0;
task{1}{1}.randVars.calculated.trackResp = nan(ceil(task{1}{1}.segmax(1)*myscreen.framesPerSecond),2);
task{1}{1}.randVars.calculated.trackStim = nan; % defined later in the task configurations
task{1}{1}.randVars.calculated.trackEye  = nan(ceil(task{1}{1}.segmax(1)*myscreen.framesPerSecond),2);
task{1}{1}.randVars.calculated.trackTime = nan(1,ceil(task{1}{1}.segmax(1)*myscreen.framesPerSecond));
task{1}{1}.randVars.calculated.trackEyeTime = nan(1,ceil(task{1}{1}.segmax(1)*myscreen.framesPerSecond)); % for referencing edf file

%% Experimental Design
% predefined packages
selected_packages = {'balancedTree', 'independent5', 'singleblob'};
design = load_packages(myscreen, selected_packages);

% design.offset = 0; % if there are additional tasks before this, add the number of tasks here.
task = configureExperiment(task, myscreen, design, design_block);

%% initialize
% intiailize task
disp(' Initializing Task....')

for phaseN = 1:length(task{1})
    [task{1}{phaseN} myscreen] = initTask(task{1}{phaseN},myscreen,...
        @startSegmentCallback,@screenUpdateCallback,@responseCallback,@initTrialCallback);
end

% initialize stimulus
disp(' Initializing Stimulus....') 

myscreen = initStimulus('stimulus',myscreen); % what does this do???

% save parameters in the stimulus struct too
stimulus.experimenter       = experimenter;
stimulus.design             = design; % also save design
stimulus.design_block       = design_block;

if stimulus.experimenter.grabframe
    global frame
    frame = {};
end

%% Eye calibration
if ~stimulus.experimenter.noeye
    disp(' Calibrating Eye ....')
    myscreen = eyeCalibDisp(myscreen);
    
    % let the user know
    disp('(trackpos) Eye calibration finished...');
end

%% generate all the stimuli, if needed
if experimenter.precompute_gen
    if ~isempty(dir(fullfile(experimenter.precompute_path,'*.mat')))
        r = [];
        while isempty(r)
          question = 'There are mat files in precompute_path. Delete and regenerate?';
          r = input(sprintf('%s (y/n)? ',question),'s');
		     % make sure we got a valid answer
		     if (lower(r) == 'n')
		       r = 0;
		     elseif (lower(r) == 'y')
		       r = 1;
		     else
		       r =[];
             end
        end
        if r == 1
            precompute_gen_script(task, myscreen, stimulus)
        end
    else
        precompute_gen_script(task, myscreen, stimulus)
    end
end

%% run the task
disp(' Running Task....'); stimulus.t0 = mglGetSecs; % 

% let the experimentee know too...
mglClearScreen(task{1}{1}.parameter.backLum/255);
mglTextDraw('task (trackpos) starting... ', [0 0.5])
mglTextDraw('Track the brightest point of the center dot with the red mouse cursor',[0 -0.5]);
mglFlush

if ~stimulus.experimenter.showmouse, mglDisplayCursor(0);, end %hide cursor

% todo: put everything into phases instead of tasknumbers
phaseNum = 1;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    stimulus.condNum = floor(phaseNum/2)+1; % todo: make more flexible
    stimulus.blockNum = mod(phaseNum,2)+1;  % todo: make blocknum more flexible
    [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);     % update the task
    myscreen = tickScreen(myscreen,task);     % flip screen
end

%% End task
disp(' End Task....')

myscreen = endTask(myscreen,task);
mglClose 
endScreen(myscreen); mglDisplayCursor(1) %show cursor

% save the last frame
%if stimulus.experimenter.grabframe 
%    save(fullfile(stimulus.experimenter.precompute_path,'grabframe.mat'), 'frame')
%end

end

%% Experiment configurations (block design etc...)
function task = configureExperiment(task, myscreen, design, design_block)
    offset = 0; 
    ncond = length(design);
    % all blocks have an training period and a tracking period
    for condNum = (offset+1):(offset+ncond)
        phaseNum = 2*(condNum-1)+1;

        % training period
        task{1}{phaseNum}                                = task{1}{1};  
        task{1}{phaseNum}.numTrials                      = design_block.nTrials_train;
        task{1}{phaseNum}.parameter.phasescrambleOn      = 0;
        task{1}{phaseNum}.parameter.noiseLum             = 0;
        task{1}{phaseNum}.randVars.calculated.trackStim  = nan(ceil(task{1}{1}.segmax(1)*myscreen.framesPerSecond),2);
         
        % tracking period
        task{1}{phaseNum+1}                                = task{1}{phaseNum}; % inherit other parameters  
        task{1}{phaseNum+1}.numTrials                      = design_block.nTrials_track;
        task{1}{phaseNum+1}.parameter.phasescrambleOn      = 1;
        task{1}{phaseNum+1}.parameter.noiseLum             = 32;
    end
end

%% Initialize trials 
function [task, myscreen] = initTrialCallback(task, myscreen)
    global stimulus    
    disp(['(trackpos_hierarchical) running stimulus ', stimulus.design(stimulus.condNum).name])
    
    % todo: use these
    % task.thistrial.thisphase
    % task.thistrial.thisseg
    % task.blockTrialnum
    % task.trailnum
    % task.blocknum 
    
    % todo: this part is called similarly twice by the precompute and
    % the inittrial. try to resolve this difference as we change things
    % up. 

    % find current design parameters
    currDesign = stimulus.design(stimulus.condNum);
    stimulus.L                  = currDesign.L; % motion structure
    stimulus.start_pos          = currDesign.start_pos;
    stimulus.colors             = currDesign.colors;
    stimulus.sourceStd          = currDesign.sourceStd; % std/s ... we convert to std per frame later
    stimulus.stimSize           = currDesign.stimSize; 

    % other parameters (todo: take these out of thistrial)
    stimulus.phasescrambleOn    = task.thistrial.phasescrambleOn;
    stimulus.stimLum            = task.thistrial.stimLum;
    stimulus.backLum            = task.thistrial.backLum;
    stimulus.noiseLum           = task.thistrial.noiseLum;

    % fixation cross
    stimulus.fixColor = [1 1 1];
    
    % precomputed stuff here
    savedfile = [stimulus.experimenter.precompute_path, ...
        sprintf('cond%d_block%d_trial%d.mat',stimulus.condNum, stimulus.blockNum, task.trialnum)];
    
    if stimulus.experimenter.precompute == 1 % load stimulus and background.     
        if ~exist(savedfile), error('Cannot find precomputed image file');, end
        
        img = load(savedfile, 'img');
        stimulus.img = img.img;
        s = load(savedfile, 'stimulus');
        s = s.stimulus;
        
        % todo: do I need to load anything else from the stimulus?
        % I think the stimulus and s should be the same..
        stimulus.source_velx    = s.source_velx;
        stimulus.source_vely    = s.source_vely;
        stimulus.posx           = s.posx;
        stimulus.posy           = s.posy;
        stimulus.objects        = s.objects;
    elseif stimulus.experimenter.precompute == 2 % precompute background only        
        if ~exist(savedfile), error('Cannot find precomputed image file');, end

        img = load(savedfile, 'backgroundnoise_rgb');
        stimulus.backgroundnoise_rgb = img.backgroundnoise_rgb;
        s = load(savedfile, 'stimulus');
        s = s.stimulus;

        stimulus.source_velx    = s.source_velx;
        stimulus.source_vely    = s.source_vely;
        stimulus.posx           = s.posx;
        stimulus.posy           = s.posy;
        stimulus.objects        = s.objects;
    else
        stimulus = InitGroupStim(stimulus,myscreen,task);
    end
    
    task.thistrial.framecount = 0;
end

%% Start segment
function [task myscreen] = startSegmentCallback(task, myscreen)
% S1: Stimulus (30s)
% S2: Fixation (3s)
global stimulus 

if task.thistrial.thisseg == 1           
    % center mouse on first object
    x_img = stimulus.posx(1,1);  y_img = stimulus.posy(1,1); 
    x_screen = x_img*myscreen.screenWidth/myscreen.imageWidth + myscreen.screenWidth/2;
    y_screen = y_img*myscreen.screenHeight/myscreen.imageHeight + myscreen.screenHeight/2;
    mglSetMousePosition(ceil(x_screen),floor(y_screen), myscreen.screenNumber); % correct for screen resolution???
    if ~stimulus.experimenter.showmouse, mglDisplayCursor(0);, end %hide cursor
        
    %{ 
    see the stimulus
    for idx = 1:length(stimulus.backnoise)
        mglClearScreen(0);
        mglBltTexture(stimulus.backnoise{idx},[0 0 myscreen.imageWidth myscreen.imageHeight])
        mglBltTexture(stimulus.gaussian,stimulus.position);
        mglFlush
        stimulus        = updateTarget(stimulus,myscreen,task); % update position.
        pause(1/myscreen.framesPerSecond)
    end
    %}   
    % frame counter 
    task.thistrial.framecount = 0;
    
    % save images of the screen.
    if stimulus.experimenter.grabframe
        global frame
        % save frames every task segment
        savefile = [stimulus.experimenter.precompute_path, ...
            sprintf('grabframe_task%d_block%d_trial%d.mat',...
                    stimulus.condNum, stimulus.blockNum, task.trialnum)];
        save(savefile, 'frame','-v7.3')
        
        % delete and intialize new frames
        frame = {}; frame{stimulus.design_block.time_stim*myscreen.framesPerSecond} = [];
    end
    
    % time debugging
    % stimulus.timedebug = nan(15,ceil(task.segmax(1)*myscreen.framesPerSecond)); %nanmean(diff(stimulus.timedebug),2)
else %intertrial interval
    % *** fixation cross.

end    

end

%% screen update
function [task myscreen] = screenUpdateCallback(task, myscreen)
% S1: Stimulus (30s)
% S2: Fixation (3s)
% todo: make timedebugging simpler
% stimulus.timedebug(9,task.thistrial.framecount+1) = mglGetSecs(stimulus.t0); % takes ~0.0087425s

%% Update Screen
global stimulus 

% for debugging gamma:
%gam = 0.5;
%mglSetGammaTable((0:1/256:1).^gam); % set table to linear, so we can just clear Screen

mglClearScreen(stimulus.backLum/256);

if (task.thistrial.thisseg== 1)
    task.thistrial.framecount = task.thistrial.framecount + 1;
    
    % **** display stimulus and background
    if stimulus.experimenter.precompute == 1
        img = mglCreateTexture(stimulus.img(:,:,:,task.thistrial.framecount));
        mglBltTexture(img,[0 0 myscreen.imageWidth myscreen.imageHeight]) %strecth
    elseif stimulus.experimenter.precompute == 2
        img = mglCreateTexture(stimulus.backgroundnoise_rgb(:,:,:,task.thistrial.framecount));
        mglBltTexture(img,[0 0 myscreen.imageWidth myscreen.imageHeight]) %strecth
        for obj = 1:length(stimulus.objects)
            stimobj = stimulus.objects{obj};

            mglBltTexture(stimobj.gaussian,...
                [stimobj.position(1,task.thistrial.framecount), ...
                 stimobj.position(2,task.thistrial.framecount)]);

            % check stimuli, flush onto screen
            %{
            % figure;plot(stimobj.position(1,:))
            mglClearScreen(0)
            mglFlush
            %}
            % todo: blip background
        end
    else
        for obj = 1:length(stimulus.objects)
            stimobj = stimulus.objects{obj};

            mglBltTexture(stimobj.gaussian,...
                [stimobj.position(1,task.thistrial.framecount), ...
                 stimobj.position(2,task.thistrial.framecount)]);

            % check stimuli, flush onto screen
            %{
            % figure;plot(stimobj.position(1,:))
            mglClearScreen(0)
            mglFlush
            %}
            % todo: blip background
        end
    end

    % *** display mouse position
    mInfo = mglGetMouse(myscreen.screenNumber);
    mimg_x = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
    mimg_y = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;
    mglGluDisk(mimg_x, mimg_y, 0.1, [1 0 0])

    % *** record stimulus position and mouse position
    task.thistrial.trackStim(task.thistrial.framecount,:) = [stimulus.posx(1,task.thistrial.framecount), ...
                                                             stimulus.posy(1,task.thistrial.framecount)];
    task.thistrial.trackResp(task.thistrial.framecount,:) = [mimg_x, mimg_y];
    task.thistrial.trackTime(task.thistrial.framecount)   = mglGetSecs(stimulus.t0);
    
    % change fixation
    if stimulus.experimenter.fixateCenter == 1
        mglGluAnnulus(0,0,0.2,0.3,stimulus.fixColor,60,1);
        mglGluDisk(0,0,0.1,rand(1,3),60,1);
    end
    
elseif (task.thistrial.thisseg == 2) % fixation segment. 
    if stimulus.experimenter.fixateCenter == 1 % stop the flashing
        rng(task.thistrial.randomSeed,'twister');
        mglGluDisk(0,0,0.1,rand(1,3),60,1);
    end
    
    mglGluAnnulus(0,0,0.2,0.3,stimulus.fixColor,60,1);
end

% fixation cross for all tasks. 

% stimulus.timedebug(6,task.thistrial.framecount+1) = mglGetSecs(stimulus.t0); % takes ~2.48e-5 s
% mglFlush
% stimulus.timedebug(7,task.thistrial.framecount+1) = mglGetSecs(stimulus.t0); % takes ~0.007 s ..i.e. compensating for the other steps

%% eye tracking
% track for task
% *** track for fixation???
if (~stimulus.experimenter.noeye) && any(task.thistrial.thisseg==[1])
    % mouse version for testing with no eyetracker
    if stimulus.eyemousedebug
        mInfo = mglGetMouse(myscreen.screenNumber);
        degx = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
        degy = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;

        pos = [degx, degy];
    else  % check eye pos
        [pos,postime] = mglEyelinkGetCurrentEyePos; % is this in image coordinates?
    end
        
    task.thistrial.trackEye(task.thistrial.framecount,:)  = pos;
    task.thistrial.trackEyeTime(task.thistrial.framecount) = postime;
end

% stimulus.timedebug(8,task.thistrial.framecount+1) = mglGetSecs(stimulus.t0); % takes ~2.86661e-5 s

if stimulus.experimenter.grabframe && (task.thistrial.thisseg== 1)
    global frame; frame{task.thistrial.framecount} = mglFrameGrab;
end

end

%% Get response; do nothing. 
function [task myscreen] = responseCallback(task, myscreen)

global stimulus
 
end

%% Initialize stimulus
function stimulus = InitGroupStim(stimulus,myscreen,task)
    [nobject, nsource] = size(stimulus.L);
    
    % initialize background variables
    stimulus = InitStimulus(stimulus,myscreen);
    
    % initialize each object
    s = struct();
    stimulus.objects = {};
    for n = 1:nobject
        s.color = stimulus.colors(n,:); % todo: indicate color
        s.stimSize = stimulus.stimSize(n);
        s.stimLum = stimulus.stimLum;
        stimulus.objects{n} = myInitObjects(s,myscreen); 
    end
    
    % generate path
    stimulus = generateStimPath(stimulus,myscreen);
end

function stimulus = InitStimulus(stimulus,myscreen)
    % background noise
    if ~isfield(stimulus,'noiseLum'), stimulus.noiseLum = 122;,end; % unit: luminance
    
    % background luminance
    if ~isfield(stimulus,'backLum'), stimulus.backLum = 32;,end; % unit: luminance
end

function obj = myInitObjects(obj,myscreen)  
    % generate all timesteps here. **** 
    
    % stimulus size
    if ~isfield(obj,'stimSize'), obj.stimSize = 0.5;,end %unit: imageX, in deg. 
    obj.patchsize = min(6*obj.stimSize,min(myscreen.imageWidth,myscreen.imageHeight));
        
    % stimulus luminance
    if ~isfield(obj,'stimLum'), obj.stimLum = 122;,end %unit: luminance
            
    % stimulus color
    if ~isfield(obj,'color'), obj.color = [1 1 1]';, end
    if isrow(obj.color), obj.color = obj.color';, end
    
    % for loading the objects
    if isfield(obj,'gaussian'), mglDeleteTexture(obj.gaussian);, end 
    gaussian               =  mglMakeGaussian(obj.patchsize,obj.patchsize,...
                                              obj.stimSize,obj.stimSize)*(obj.stimLum);
    % gaussian_rgb = 255*ones(4,size(gaussian,2),size(gaussian,1),'uint8');
    gaussian_rgb           = 255*repmat([obj.color; 1], 1, size(gaussian,2),size(gaussian,1));
    gaussian_rgb(4,:,:)    = round(gaussian');
    gaussian_rgb           = uint8(gaussian_rgb);
    obj.gaussian           = mglCreateTexture(gaussian_rgb); % pre-generate texture here
end

function stimulus = generateStimPath(stimulus,myscreen)
    disp('(genstim) Generating stimulus trajectories... looping until a suitable one is found...')

    tic
    [nobject, nsource] = size(stimulus.L);
    nframes            = myscreen.framesPerSecond*stimulus.design_block.time_stim + 100; % 30s; %/downsample_timeRes;
    
    source_stds = stimulus.sourceStd/sqrt(myscreen.framesPerSecond); % change to deg/frame

    % todo: additional stimulus check here
    if isrow(source_stds), source_stds = source_stds';, end
    %if isrow(stimulus.start_pos), stimulus.start_pos = stimulus.start_pos';, end
    if isrow(stimulus.stimSize), stimulus.stimSize = stimulus.stimSize';,end
    
    source_velx      = nan(nsource,nframes); % source velocity
    source_vely      = nan(nsource,nframes); % source velocity
    
    posx             = nan(nobject,nframes); % object x 
    posy             = nan(nobject,nframes); % object x 
    
    maxiters = 10000;
    for iters = 1:maxiters
        source_velx = repmat(source_stds, [1, nframes]) .* normrnd(0, 1, nsource,nframes);
        source_vely = repmat(source_stds, [1, nframes]) .* normrnd(0, 1, nsource,nframes);
             
        velx = stimulus.L*source_velx;
        vely = stimulus.L*source_vely;
        velx(:,1) = 0; % to make sure the stimulus starts at initial position
        vely(:,1) = 0;
        
        posx = repmat(stimulus.start_pos(:,1), [1, nframes]) + cumsum(velx,2);
        posy = repmat(stimulus.start_pos(:,2), [1, nframes]) + cumsum(vely,2);
        
        % check overlap between stimuli (3 std boundary?) 
        % true if there are at least one overlap        
        if nobject > 1
            check_overlap = false(nchoosek(nobject,2),nframes);
            obj_dist = zeros(nchoosek(nobject,2),nframes);
            idx_combo = 1;
            for obj = 1:(nobject-1)
                for otherobj = (obj+1):nobject
                    xdiff = (posx(obj,:) - posx(otherobj,:)).^2;
                    ydiff = (posy(obj,:) - posy(otherobj,:)).^2;
                    obj_dist(idx_combo,:) = sqrt(xdiff + ydiff);
                    
                    thresh = stimulus.stimSize(obj) + stimulus.stimSize(otherobj);
                    check_overlap(idx_combo,:) = obj_dist(idx_combo,:) < 3*repmat(thresh, [1, nframes]);
                    
                    idx_combo = idx_combo +1;
                end
            end
        else
            check_overlap = false(nobject,nframes);
        end
        
        %{
        % check distances
        figure; 
        subplot(5,1,1);hold on; title('X position')
        subplot(5,1,2);hold on; title('Y position')
        subplot(5,1,3);hold on; title('Pairwise distances')
        subplot(5,1,4);hold on; title('Overlap')
        cmap = colormap(hsv(nobject));
        for obj = 1:nobject
            subplot(5,1,1); plot(posx(obj,:), 'color',cmap(obj,:));
            %plot(x_ub(obj,:),':', 'color',cmap(obj,:));
            %plot(x_lb(obj,:),':', 'color',cmap(obj,:));
            subplot(5,1,2); plot(posy(obj,:), 'color',cmap(obj,:))
            %plot(y_ub(obj,:),':', 'color',cmap(obj,:));
            %plot(y_lb(obj,:),':', 'color',cmap(obj,:));
        end
        for idx_combo = 1:nchoosek(nobject,2)
            subplot(5,1,3); plot(obj_dist(idx_combo,:))
            subplot(5,1,4); plot(check_overlap(idx_combo,:))
        end
        
        % check which combination is weird
        any(~check_overlap,2)
        %}
            
        
        % check if stimulus is out of bounds
        % true if there are at least one out of bound
        x_ub = posx + 3*repmat(stimulus.stimSize, [1, nframes]);
        x_lb = posx - 3*repmat(stimulus.stimSize, [1, nframes]);
        y_ub = posy + 3*repmat(stimulus.stimSize, [1, nframes]);
        y_lb = posy - 3*repmat(stimulus.stimSize, [1, nframes]);

        check1 = x_ub > myscreen.imageWidth/2;
        check2 = y_ub > myscreen.imageHeight/2;
        check3 = x_lb < -myscreen.imageWidth/2;
        check4 = y_lb < -myscreen.imageHeight/2;
        check_bounds = check1 | check2 | check3 | check4;
        
        %check = check_bounds | check_overlap;
        % check = false; % without any checks
        
        if mod(iters, 1000) == 0
            disp('(genstim) still checking... ')
        end
                    
        if ~any(check_bounds(:)) && ~any(check_overlap(:)) % if there arent any "true"
            disp('(genstim) suitable stimulus found!')
            break
        end
    end
    
    if any(check_bounds(:)) || any(check_overlap(:))
        error('The stimulus too big/fast and cannot be constrained within the screen. Please reset parameters.')
    end

    % save variables
    stimulus.source_velx    = source_velx;
    stimulus.source_vely    = source_vely;
    stimulus.posx           = posx;
    stimulus.posy           = posy;
    
    for obj = 1:nobject
        stimulus.objects{obj}.position = [posx(obj,:); posy(obj,:)];
    end
    
    tt = toc;
    disp(['(genstim) Stimulus trajectories generated! Took ', num2str(tt), 'secs'])
end

%% predefined stimuli
function design = load_packages(myscreen, selected_packages)
% selected_packages: cell of strings indicating which packages to use

design(length(selected_packages)) = struct('name',[],...
    'L',[], 'start_pos', [], 'colors',[], 'sourceStd', [], 'stimSize', []);

%% possible predefined stimuli
% one blob
singleblob = struct();
singleblob.name = 'singleblob';
singleblob.L = 1;
singleblob.start_pos = [0,0];
singleblob.colors = [1, 1, 1]; 
singleblob.sourceStd = 0.5; 
singleblob.stimSize = 0.5;

sum3std = 0.5;
% 9 motion sources (7 effective -- 3 for the center collapses to one
% source, but to balance each of the object effective variance..)
% two branched tree; 1 => body; 2,3=> layer2; 4-5 => layer 3
balancedTree = struct();
balancedTree.name = 'balancedTree';
balancedTree.L = [1, 0, 1, 0, 0 ,0 ,1, 0, 0; ...
                 1, 1, 0, 0, 1 ,0, 0 , 0, 0; ...
                 1, 1, 0, 0, 0 ,1, 0 , 0, 0; ... 
                 1, 0, 0, 1, 0 , 0 , 0 , 1, 0;...
                 1, 0, 0, 1, 0 , 0 , 0 , 0, 1]; 
start_pos = [0,0;-1,1;-1,-1;1,1;1,-1]; % scale by 1/8 of the screen width in degrees
balancedTree.start_pos  = start_pos .* repmat([myscreen.imageWidth/4, myscreen.imageHeight/4], [size(start_pos,1),1]);
balancedTree.colors       = [1, 1, 1;...
                             1, 0, 0;...
                             1, 0, 0;...
                             0, 0, 1;...
                             0, 0, 1]; 
std1 = sqrt(sum3std^2/3);
balancedTree.sourceStd    = ones(9,1)*std1;
balancedTree.stimSize     = [0.5; 0.5; 0.5; 0.5; 0.5];

% 5 independent objects
independent5 = struct();
independent5.name = 'independent5';
independent5.L = eye(5); 
start_pos = [0,0;-1,1;-1,-1;1,1;1,-1]; % scale by 1/8 of the screen width in degrees
independent5.start_pos  = start_pos .* repmat([myscreen.imageWidth/4, myscreen.imageHeight/4], [size(start_pos,1),1]);
independent5.colors       = [1, 1, 1;...
                             1, 0, 0;...
                             1, 0, 0;...
                             0, 0, 1;...
                             0, 0, 1]; 
independent5.sourceStd    = ones(5,1) * sum3std;
independent5.stimSize     = [0.5; 0.5; 0.5; 0.5; 0.5];

% some predefined structures todo: make these into "packages"
L1 = 1;
L2 = diag(ones(5,1)); ... % 5 independent motion
L3 = [1,  0, 0;... 
      1,  1, 0; ...
      1, -1, 0; ...
      1,  0, 1; ...
      1,  0, -1]; % motion sources; contrasting motion
L4 = [1,  0, 0, 0, 0;... 
      1,  1, 0, 0, 0; ...
      1,  0, 1, 0, 0; ...
      1,  1, 0, 1, 0; ...
      1,  0, 1, 0, 1]; ... % 5 motion sources. two branched tree; 1 => body; 2,3=> layer2; 4-5 => layer 3
L5 = [1, 0, 1, 0, 0 ,0 ,1, 0, 0; ...
      1, 1, 0, 0, 1 ,0, 0 , 0, 0; ...
      1, 1, 0, 0, 0 ,1, 0 , 0, 0; ... 
      1, 0, 0, 1, 0 , 0 , 0 , 1, 0;...
      1, 0, 0, 1, 0 , 0 , 0 , 0, 1]; % 7 motion sources. two branched tree; 1 => body; 2,3=> layer2; 4-5 => layer 3
  
%% add all the stimuli
for n = 1:length(selected_packages)
    eval(['design(', num2str(n), ') = ', selected_packages{n} ';'])
end

% check that the dimensions are correct
for n = 1:length(selected_packages)
    [nobject, nsources] = size(design(n).L);
    
    [a, c2] = size(design(n).start_pos);
    if a ~= nobject, error('The start_pos should have nobject elements');,end
    if c2 ~= 2, error('The start_pos should have 2 positions (x,y)');,end
    [a, c3] = size(design(n).colors);
    if a ~= nobject, error('The colors should have nobject elements');,end
    if c3 ~= 3, error('The colors should have 3 colors elements (rgb)');,end
    b = size(design(n).sourceStd);
    if b ~= nsources, error('The sourceStd should have nsources colors elements');,end
    a = size(design(n).stimSize);
    if a ~= nobject, error('The stimSize should have nobject elements');,end
end

end

%% Utility
function [stimx, stimy] = convertNearestPixel(myscreen,x_img,y_img)
    x_screen = x_img*myscreen.screenWidth/myscreen.imageWidth + myscreen.screenWidth/2;
    y_screen = y_img*myscreen.screenHeight/myscreen.imageHeight + myscreen.screenHeight/2;
    
    stimx = (ceil(x_screen)-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
    stimy = (floor(y_screen)-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight; % what is imagewidth???
end

function img = add_rgba_uint8(rgba1, rgba2)
% adding two uint8 rgba images (4 x x x y)
% https://stackoverflow.com/questions/28900598/how-to-combine-two-colors-with-varying-alpha-values

    if size(rgba1) ~= size(rgba2)
        error('size must be the same for adding two rgba images')
    end

    img = zeros(size(rgba1),'uint8');
    
    if length(size(rgba1)) == 3
        weight1 = (1-double(rgba1(4,:,:)/255)).*double(rgba2(4,:,:))/255;
        weight2 = double(rgba1(4,:,:))/255;
        img(4,:,:) = uint8(255*(weight1 + weight2)); % get alpha channel
        for idx = 1:3
            img(idx,:,:) = uint8((weight1 .* double(rgba2(idx,:,:)) + ...
                weight2 .* double(rgba1(idx,:,:))) ./ (weight1+weight2));
        end
    elseif length(size(rgba1)) == 4
        weight1 = (1-double(rgba1(4,:,:,:)/255)).*double(rgba2(4,:,:,:))/255;
        weight2 = double(rgba1(4,:,:,:))/255;
        img(4,:,:,:) = uint8(255*(weight1 + weight2)); % get alpha channel
        for idx = 1:3
            img(idx,:,:,:) = uint8((weight1 .* double(rgba2(idx,:,:,:)) + ...
                weight2 .* double(rgba1(idx,:,:,:))) ./ (weight1+weight2));
        end
    end
    
    
end

function precompute_gen_script(task, myscreen, stimulus)
% stimulus is not global, contained inside the script
% note that the precomputeStimulus saves and loads the stimulus struct!

    for condNum = 1:length(stimulus.design)
    for blockNum = 1:2
        % todo: make the blockNumber per condition more flexible?
        phaseNum = 2*(condNum-1) + blockNum;
    for trialNum = 1:task{1}{phaseNum}.numTrials
        stimulus.condNum = condNum;
        stimulus.blockNum = blockNum;
        stimulus.trialNum = trialNum; % keep track of trial number in stimulus for precomputing
    
        % find current design parameters
        currDesign                  = stimulus.design(stimulus.condNum);
        stimulus.L                  = currDesign.L; % motion structure
        stimulus.start_pos          = currDesign.start_pos;
        stimulus.colors             = currDesign.colors;
        stimulus.sourceStd          = currDesign.sourceStd; % std/s ... we convert to std per frame later
        stimulus.stimSize           = currDesign.stimSize; 

        % other parameters (todo: take these out of thistrial)
        stimulus.phasescrambleOn    = task{1}{phaseNum}.parameter.phasescrambleOn;
        stimulus.stimLum            = task{1}{phaseNum}.parameter.stimLum;
        stimulus.backLum            = task{1}{phaseNum}.parameter.backLum;
        stimulus.noiseLum           = task{1}{phaseNum}.parameter.noiseLum;

        % fixation cross
        stimulus.fixColor = [1 1 1];
        stimulus = InitGroupStim(stimulus,myscreen,task);
        img = precomputeStimulus(myscreen,stimulus);
    end
    end
    end
    
    done = 1;
end

function img = precomputeStimulus(myscreen,stimulus)
% note that the precomputeStimulus saves and loads the stimulus struct!
% generate a noisy stimulus
% Approximate generation time on csnl:
%     -(15s trial, 60 hz, ds= 2) stimulus generation: 380s
%     -(15s trial, 60 hz, ds= 2) background generation: 661s
%     -(15s trial, 60 hz, ds= 2) saving: 200s for stimulus
disp(['(genstim_precompute) Generating images for: ', ...
          'cond ', num2str(stimulus.condNum), ...
       ', block ', num2str(stimulus.blockNum), ...
       ', trial ', num2str(stimulus.trialNum), '...'])
   
clearvars img stim_rgb backgroundnoise_rgb

timer_allgen = tic;
savefile = [stimulus.experimenter.precompute_path, ...
    sprintf('cond%d_block%d_trial%d.mat',stimulus.condNum, stimulus.blockNum, stimulus.trialNum)];

downsample_spatRes  = stimulus.experimenter.downsample_spatRes;
nframes             = myscreen.framesPerSecond*stimulus.design_block.time_stim +100; % 30s; %/downsample_timeRes; 
xsize_deg           = round(myscreen.imageWidth/downsample_spatRes);
ysize_deg           = round(myscreen.imageHeight/downsample_spatRes);

% make a Gaussian for the stimulus and background intialization
backgaussian = mglMakeGaussian(xsize_deg,ysize_deg,...
        stimulus.objects{1}.stimSize/downsample_spatRes,...
        stimulus.objects{1}.stimSize/downsample_spatRes)*255;
backgaussianFFT = getHalfFourier(backgaussian);

x_npixels = size(backgaussian,2);
y_npixels = size(backgaussian,1);

%% generate objects
stim_rgb  = zeros(4,x_npixels,y_npixels,nframes,'uint8');

tic
for idx1 = 1:nframes 
    for obj = 1:length(stimulus.objects)
        stimobj = stimulus.objects{obj};
        obj_img = zeros(4,x_npixels,y_npixels,'uint8'); % object image for a single time
        
        % generate Gaussian 
        % todo: flip xsize and ysize??
        gaussian = mglMakeGaussian(xsize_deg,ysize_deg,...
                    stimobj.stimSize/downsample_spatRes, stimobj.stimSize/downsample_spatRes,...
                    stimobj.position(1,idx1)/downsample_spatRes, stimobj.position(2,idx1)/downsample_spatRes);
                
        %stim_rgb(1:3,:,:,idx1)  = stim_rgb(1:3,:,:,idx1) + ...
        %          uint8(255*repmat(ccc, [1, x_npixels,y_npixels]) ...
        %            .* permute(cat(3,gaussian',gaussian',gaussian')/max(gaussian(:)),[3,1,2])) ;
        obj_img(1:3,:,:)    = uint8(255*repmat(stimobj.color, [1, x_npixels,y_npixels]));
        obj_img(4,:,:)      = uint8(gaussian'/max(gaussian(:)) * stimobj.stimLum); % alpha channel
        
        % add to the ground image
        stim_rgb(:,:,:,idx1) = add_rgba_uint8(stim_rgb(:,:,:,idx1), obj_img);
        
        % check stimuli
        % figure;plot(stimobj.position(1,:))
        % imagesc the stimuli
        %{
        sum(sum(sum(stim_rgb(:,:,:,idx1) == obj_img)))
        
        figure;
        subplot(2,2,1);imagesc(squeeze(double(stim_rgb(1,:,:,idx1)) - double(obj_img(1,:,:))));colorbar;
        subplot(2,2,2);imagesc(squeeze(double(stim_rgb(2,:,:,idx1)) - double(obj_img(2,:,:))));colorbar;
        subplot(2,2,3);imagesc(squeeze(double(stim_rgb(3,:,:,idx1)) - double(obj_img(3,:,:))));colorbar;
        subplot(2,2,4);imagesc(squeeze(double(stim_rgb(4,:,:,idx1)) - double(obj_img(4,:,:))));colorbar;
        
        figure;
        subplot(4,2,1);imagesc(squeeze(double(stim_rgb(1,:,:,idx1))));colorbar;
        subplot(4,2,2);imagesc(squeeze(double(obj_img(1,:,:))));colorbar;
        subplot(4,2,3);imagesc(squeeze(double(stim_rgb(2,:,:,idx1))));colorbar;
        subplot(4,2,4);imagesc(squeeze(double(obj_img(2,:,:))));colorbar;
        subplot(4,2,5);imagesc(squeeze(double(stim_rgb(3,:,:,idx1))));colorbar;
        subplot(4,2,6);imagesc(squeeze(double(obj_img(3,:,:))));colorbar;
        subplot(4,2,7);imagesc(squeeze(double(stim_rgb(4,:,:,idx1))));colorbar;
        subplot(4,2,8);imagesc(squeeze(double(obj_img(4,:,:))));colorbar;
        %}               
        % flush onto screen
        %{
        mglClearScreen(0)
        txt = mglCreateTexture(obj_img);
        mglBltTexture(txt,[0 0 myscreen.imageWidth myscreen.imageHeight]) %strecth again.
        mglFlush
        
        mglClearScreen(0)
        txt = mglCreateTexture(stim_rgb(:,:,:,idx1));
        mglBltTexture(txt,[0 0 myscreen.imageWidth myscreen.imageHeight]) %strecth again.
        mglFlush
        %}
    end   
    
    %{
    mglClearScreen(0)
    txt = mglCreateTexture(stim_rgb(:,:,:,idx1));
    mglBltTexture(txt,[0 0 myscreen.imageWidth myscreen.imageHeight]) %strecth again.
    mglFlush
    %}
end
t_stimgen = toc;
disp(['(genstim_precompute) Took ', num2str(t_stimgen), ' sec to generate object images'])

%% noise
% background noise
tic
backgroundnoise_rgb  = zeros(4,x_npixels,y_npixels,nframes,'uint8');

if stimulus.noiseLum ~= 0 
    noise_alpha = stimulus.noiseLum;
    % for random color noise distribution:
    %{
    colorstd = 5; %stimulus.objects{1}.stimSize
    backgaussian_color = mglMakeGaussian(xsize_deg,ysize_deg,...
            colorstd/downsample_spatRes,...
            colorstd/downsample_spatRes)*255;
    backgaussian_color_FFT = getHalfFourier(backgaussian_color);
    %}

    for idx1 = 1:nframes %nframes
        % draw noise from the stimulus colors
        for obj = 1:length(stimulus.objects)
            stimobj = stimulus.objects{obj};
            back_img = zeros(4,x_npixels,y_npixels,'uint8'); % object image for a single time

            back                = backgaussianFFT; %0.02s
            back.phase          = rand(size(back.mag))*2*pi; % scramble phase % 0.02s
            backgroundnoise     = round(reconstructFromHalfFourier(back));   %0.04s
            back_img(4,:,:)     = uint8(backgroundnoise'/max(backgroundnoise(:))*noise_alpha);  % normalize contrast %0.025s        
            back_img(1:3,:,:)   = uint8(255*repmat(stimobj.color, [1, x_npixels,y_npixels]));

            % add to the ground image
            backgroundnoise_rgb(:,:,:,idx1) = add_rgba_uint8(backgroundnoise_rgb(:,:,:,idx1), ...
                                                             back_img);
            %{                                             
            mglClearScreen(0)
            txt = mglCreateTexture(back_img);
            mglBltTexture(txt,[0 0 myscreen.imageWidth myscreen.imageHeight]) %strecth again.
            mglFlush

            mglClearScreen(0)
            txt = mglCreateTexture(backgroundnoise_rgb(:,:,:,idx1));
            mglBltTexture(txt,[0 0 myscreen.imageWidth myscreen.imageHeight]) %strecth again.
            mglFlush
            %}
        end    

        % draw from random color distribution
        %{
        for dim = 1:4
            if dim == 4 % alpha channels => Gaussian noise
                back                        = backgaussianFFT; %0.02s
                back.phase                  = rand(size(back.mag))*2*pi; % scramble phase % 0.02s
                backgroundnoise             = round(reconstructFromHalfFourier(back));   %0.04s
                backgroundnoise_rgb(4,:,:,idx1)  = uint8(backgroundnoise'/max(backgroundnoise(:))*noise_alpha);  % normalize contrast %0.025s
            else % color channels => define a color distribution?.
                back                        = backgaussian_color_FFT; %0.02s
                back.phase                  = rand(size(back.mag))*2*pi; % scramble phase % 0.02s
                backgroundnoise             = round(reconstructFromHalfFourier(back));   %0.04s            
                backgroundnoise_rgb(dim,:,:,idx1)  = uint8(backgroundnoise'/max(backgroundnoise(:))*255);
            end
        end 
        %}

        % check the generated textures
        %{
        mglClearScreen(0)
        txt = mglCreateTexture(backgroundnoise_rgb(:,:,:,idx1));
        mglBltTexture(txt,[0 0 myscreen.imageWidth myscreen.imageHeight]) %strecth again.
        mglFlush

        figure;
        subplot(2,2,1);imagesc(squeeze(double(backgroundnoise_rgb(1,:,:,idx1))));colorbar;
        subplot(2,2,2);imagesc(squeeze(double(backgroundnoise_rgb(2,:,:,idx1))));colorbar;
        subplot(2,2,3);imagesc(squeeze(double(backgroundnoise_rgb(3,:,:,idx1))));colorbar;
        subplot(2,2,4);imagesc(squeeze(double(backgroundnoise_rgb(4,:,:,idx1))));colorbar;
        %}
    end
end
t_backgen = toc;
disp(['(genstim_precompute) Took ', num2str(t_backgen), ' sec to generate background images'])

img = add_rgba_uint8(stim_rgb, backgroundnoise_rgb);
save(savefile, 'stimulus', 'img', 'stim_rgb','backgroundnoise_rgb', '-v7.3')

ttt = toc(timer_allgen);
disp(['(genstim_precompute) Trial generation complete. Took ', num2str(ttt), ' sec'])
end
