%        $Id: $
%      usage: trackpos
%         by: Josh Ryu
%       date: 05/05/2021
%    purpose: 

% todo: make noiseLum more flexible

function myscreen = trackpos(varargin)
%getArgs(varargin,{'subjectID=s999','centerX=10','centerY=0','diameter=16'}); getArgs(varargin,{'subjectID=-1'});
% set up screen
myscreen = struct();
if isempty(mglGetSID)
    myscreen.subjectID  = -1;
else
    myscreen.subjectID  = mglGetSID;
end

%myscreen.displayName = 'debug'; myscreen.screenNumber = 1; 
%myscreen.screenWidth = 860; myscreen.screenHeight = 600; 
%myscreen.hideCursor = 1;
myscreen                = initScreen(myscreen);

%% parameters
% Experimenter parameters
%todo:  check these throughout the code!!
exp.debug               = 0; % debug code
exp.noeye               = 0; % 1 if no eyetracking (mouse for eye); 0 if there is eye tracking `
exp.showmouse           = 0; % show mouse during everything
exp.fixateCenter        = 1; %
exp.usejoystick 		= 1; % use joystick. Need "simulink 3D animation" package downloaded. 

exp.backprecompute      = '/Users/gru/proj/grustim/trackpos/trackpos.mat'; % precomputed background
exp.downsample_timeRes  = 1; % downsample temporal resolution of background noise the by this factor.
exp.phasescrambleOn     = 1; %
exp.whitenoiseOn        = 0; % 1: white noise; 2: 

exp.grabframe           = 0; % capture frames

global stimulus; stimulus = struct;
stimulus.exp = exp;

task = {}; 
% specify task design
sb      = stillblob();
phase1  = sb.configureExperiment(stimulus,task,myscreen);
task{1} = {phase1};

task{1} = add_calculated_params(task{1});
totaldur = approximate_total_task_dur(task);

%% initialize
% intiailize task
disp(' Initializing Task....')

for phaseN = 1:length(task{1})
    [task{1}{phaseN} myscreen] = initTask(task{1}{phaseN},myscreen,...
        @startSegmentCallback,@screenUpdateCallback,@responseCallback,@initTrialCallback);
end

% initialize stimulus
myscreen = initStimulus('stimulus',myscreen); % save the stimulus

%% load background
if stimulus.exp.phasescrambleOn == 1;
    disp('Loading phase scrambled background noise...')

    tic
    savefile = stimulus.exp.backprecompute;
    % savefile            = '/Users/joshua/data/trackpos_2afc/trackpos.mat'; % just use noise 1 and permute
    if ~exist(savefile,'file')
        error('need background file')
    end

    load(savefile,'backgroundnoise_rgb');

    if isfield(stimulus,'backnoise')
        for idx = 1:length(stimulus.backnoise)
            mglDeleteTexture(stimulus.backnoise{idx});
        end
    end

    % create all background textures and then load them later
    if stimulus.exp.debug 
        nnn = 200;
    else
        nnn = size(backgroundnoise_rgb,4);
    end
    for idx = 1:nnn %too big?? memory?
        stimulus.backnoise{idx} = mglCreateTexture(backgroundnoise_rgb(:,:,:,idx));
    end

    clearvars('backgroundnoise_rgb')
    toc
end

%% grabframe
if stimulus.exp.grabframe
    global frame
    frame = {};
end

%% Eye calibration and check joystick
if ~stimulus.exp.noeye && ~stimulus.exp.debug
    disp(' Calibrating Eye ....')
    myscreen = eyeCalibDisp(myscreen);
    
    % let the user know
    disp(sprintf('(trackpos) Starting Run...'));
end

if stimulus.exp.usejoystick
global joy; joy = vrjoystick(1); % use simulink 3d animation to load joystick object
if isempty(joy)
    stimulus.exp.usejoystick = 0;
    exp = stimulus.exp;
    disp(' FAILED TO FIND JOYSTICK! MAKE SURE SIMULINK 3D ANIMATION PACKAGE IS INSTALLED AND THE JOYSTICK IS PROPERLY CONNECTED');
    disp(' USING MOUSE FOR TRACKING');
else
    joy_params              = struct();
    joy_params.maxv         = 0.2;
    joy_params.deadzone     = 0.02;
    joy_params.sensitivity  = 2;
    joy_params.poly_order   = 1.2;
    stimulus.joy_params     = joy_params;
end
end

%% run the task
disp(' Running Task....'); stimulus.t0 = mglGetSecs; % 

% let the experimentee know too...
% mglClearScreen(task{1}{1}.parameter.backLum/255);
mglTextDraw('task starting... ', [0 0.5])
mglTextDraw('Track the target with the red pointer',[0 -0.5]);
mglTextDraw('When you are ready, press backtick to go to next trial',[0 -1.5]);
mglFlush

if ~exp.showmouse, mglDisplayCursor(0);, end %hide cursor

phaseNum = 1;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);     % update the task
    myscreen = tickScreen(myscreen,task);     % flip screen
end

%% End task
disp(' End Task....')

if isfield(stimulus,'backnoise')
    for idx = 1:length(stimulus.backnoise)
        mglDeleteTexture(stimulus.backnoise{idx});
    end
end

myscreen = endTask(myscreen,task);
mglClose 
endScreen(myscreen); mglDisplayCursor(1) %show cursor

if stimulus.exp.grabframe
    save('/Users/joshryu/Dropbox/GardnerLab/Stanford/Current/FYP/FYP talk/trackposTask.mat', 'frame')
end

end

%% Initialize trials 
function [task myscreen] = initTrialCallback(task, myscreen)
    global stimulus    
       
    stimulus.stepStd    = task.thistrial.stimStep/sqrt(myscreen.framesPerSecond); % convert to std per frame; i.e. variance of shz gaussian RVs.
    stimulus.stimLum    = task.thistrial.stimLum;
    
    stimulus.phasescrambleOn    = task.thistrial.phasescrambleOn;
    stimulus.backLum            = task.thistrial.backLum;
    stimulus.noiseLum           = task.thistrial.noiseLum;
    
    stimulus = trackposInitStimulus(stimulus,myscreen); 
    
    % save seed for generating random textures
    rng('shuffle','twister'); 
    s=rng; 
    task.thistrial.randomSeed = s.Seed;
    rng(task.thistrial.randomSeed,'twister');
    
    %% noise
    if stimulus.exp.phasescrambleOn == 1
        nframes = myscreen.framesPerSecond*task.segmax(1) + 20;%/downsample_timeRes; 
        task.thistrial.bgpermute(1:nframes) = randi(length(stimulus.backnoise),nframes,1);
    end
    
    if mod(task.trialnum,ceil(task.numTrials/20)) == 1
        disp(['(trackpos) '  num2str(task.trialnum/task.numTrials) ...
            '% finished: Trial ' num2str(task.trialnum) ' / ' num2str(task.numTrials)]);
    end
end

%% Start segment
function [task myscreen] = startSegmentCallback(task, myscreen)
% S1: Stimulus (30s)
% S2: Fixation (3s)
global stimulus 

if task.thistrial.thisseg == 1       
    %% set stimulus
    % after the noise bc we want to recreate the texture with the seed.
    stimulus = trackposInitStimulus(stimulus,myscreen); %centerX,Y, diameter called by getArgs.
    
    if ~stimulus.exp.usejoystick
        % set mouse position to the stimulus direction. 
        x_img = stimulus.position(1);  y_img = stimulus.position(2);
        x_screen = x_img*myscreen.screenWidth/myscreen.imageWidth + myscreen.screenWidth/2;
        y_screen = y_img*myscreen.screenHeight/myscreen.imageHeight + myscreen.screenHeight/2;
        mglSetMousePosition(ceil(x_screen),floor(y_screen), myscreen.screenNumber); % correct for screen resolution???
    end

    if ~stimulus.exp.showmouse, mglDisplayCursor(0);, end %hide cursor if needed
    
    task.thistrial.initStim = stimulus.position;% [stimx, stimy];
    
    %% 
    task.thistrial.framecount = 0;
    task.thistrial.movecursor = false; % need to set to true for pointer to move

    
    if stimulus.exp.grabframe
        global frame
        %save('/Users/jryu/Dropbox/Stanford/Current/FYP/FYP talk/figures/trackposTask_nonoise_160back.mat', 'frame','-v7.3')
        %save('/Users/jryu/Dropbox/Stanford/Current/FYP/FYP talk/figures/trackposTask_noise_160back.mat', 'frame','-v7.3')
        frame = {};
        frame{task.segmax(1)*myscreen.framesPerSecond} = [];
    end
    
    % time debugging
    % stimulus.timedebug = nan(15,ceil(task.segmax(1)*myscreen.framesPerSecond)); %nanmean(diff(stimulus.timedebug),2)
else %intertrial interval
    % *** fixation cross.
end    

end

%% screen update
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus

task  = update(obj, task, myscreen, stimulus);       
task.thistrial.movecursor
stimulus.pointer




    
    

%% Update Screen

if stimulus.exp.usejoystick, global joy;,end
% stimulus.timedebug(9,task.thistrial.framecount+1) = mglGetSecs(stimulus.t0); % takes ~0.0087425s
mglClearScreen(stimulus.backLum/255);

if (task.thistrial.thisseg== 1)
    task.thistrial.framecount = task.thistrial.framecount + 1;
    
    stimulus        = updateTarget(stimulus,myscreen,task); % update position.
    
    if task.thistrial.phasescrambleOn == 1 && mod(task.thistrial.framecount, stimulus.exp.downsample_timeRes) == 0
        idx = task.thistrial.bgpermute(task.thistrial.framecount);
        mglBltTexture(stimulus.backnoise{idx},...
            [0 0 myscreen.imageWidth myscreen.imageHeight])
    end
    
    % draw stimulus
    mglBltTexture(stimulus.gaussian,stimulus.position);
    
    % stimulus.timedebug(3,task.thistrial.framecount+1) = mglGetSecs(stimulus.t0); % takes ~0.000495s

    % **&display mouse position
    if ~stimulus.exp.usejoystick
        mInfo = mglGetMouse(myscreen.screenNumber);
        stimulus.pointer(1) = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
        stimulus.pointer(2) = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;
    else
        [vx, vy] = joy2vel(joy, stimulus.joy_params, myscreen);
        stimulus = update_pointer(stimulus, [vx, vy], myscreen);
        task.thistrial.trackJoy(task.thistrial.framecount,:)  = axis(joy);
    end
    
    mglGluDisk(stimulus.pointer(1), stimulus.pointer(2), 0.1, [1 0 0])
    % stimulus.timedebug(4,task.thistrial.framecount+1) = mglGetSecs(stimulus.t0); %takes ~9.2678e-5 s

    % ***record stimulus position and mouse position  
    task.thistrial.trackStim(task.thistrial.framecount,:) = stimulus.position;
    task.thistrial.trackResp(task.thistrial.framecount,:) = [mimg_x, mimg_y];
    task.thistrial.trackTime(task.thistrial.framecount)   = mglGetSecs(stimulus.t0);
    
    % stimulus.timedebug(5,task.thistrial.framecount+1) = mglGetSecs(stimulus.t0); %takes ~7.67331e-5 s

    if stimulus.exp.fixateCenter == 1
        mglGluAnnulus(0,0,0.2,0.3,stimulus.fixColor.response,60,1);
        mglGluDisk(0,0,0.1,rand(1,3),60,1);
    end
    
elseif (task.thistrial.thisseg == 2)
    if stimulus.exp.fixateCenter == 1 % stop the flashing
        rng(task.thistrial.randomSeed,'twister');
        mglGluDisk(0,0,0.1,rand(1,3),60,1);
    end
    
    mglGluAnnulus(0,0,0.2,0.3,stimulus.fixColor.response,60,1);
end

% fixation cross for all tasks. 


% stimulus.timedebug(6,task.thistrial.framecount+1) = mglGetSecs(stimulus.t0); % takes ~2.48e-5 s
% mglFlush
% stimulus.timedebug(7,task.thistrial.framecount+1) = mglGetSecs(stimulus.t0); % takes ~0.007 s ..i.e. compensating for the other steps

%% eye tracking
% track for task
% *** track for fixation???
if (~stimulus.exp.noeye) && any(task.thistrial.thisseg==[1])
    % mouse version for testing with no eyetracker
    if stimulus.exp.eyemousedebug
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

if stimulus.exp.grabframe && (task.thistrial.thisseg== 1)
    global frame; frame{task.thistrial.framecount} = mglFrameGrab;
end

end

%% Get response; do nothing. 
function [task myscreen] = responseCallback(task, myscreen)

global stimulus
 
end

function stimulus = update_pointer(stimulus, vel,myscreen)
    pos = stimulus.pointer;
    [horz_out, vert_out] = check_oob(pos + vel, myscreen, stimulus);
    stimulus.pointer(1) = stimulus.pointer(1) + (1-horz_out)*vel(1);
    stimulus.pointer(2) = stimulus.pointer(2) + (1-vert_out)*vel(1);
end

function stimulus = updateTarget(stimulus,myscreen,task)
    % convert steps to direction gradients
    xstep = normrnd(0,stimulus.stepStd);
    ystep = normrnd(0,stimulus.stepStd);
    
    %[xstep, ystep] = convertNearestPixel(myscreen,xstep,ystep);
    
    stimulus.position = stimulus.position + [xstep, ystep];
    
    % if thre circle goes out of the image 
    if stimulus.position(1)+3*stimulus.stimStd > myscreen.imageWidth/2 ...
            || stimulus.position(1)-3*stimulus.stimStd < -myscreen.imageWidth/2,
        stimulus.position(1) = stimulus.position(1) - 2*xstep;
    end
        
    if stimulus.position(2)+3*stimulus.stimStd > myscreen.imageHeight/2 ...
            || stimulus.position(2)-3*stimulus.stimStd < -myscreen.imageHeight/2,
        stimulus.position(2) = stimulus.position(2) - 2*ystep;
    end
    
    % if eye fixation is enforced
    if stimulus.exp.fixateCenter == 1
        if norm(stimulus.position) < 1 % center of the circle is never in the middle
            stimulus.position(1) = stimulus.position(1) - 2*xstep;
            stimulus.position(2) = stimulus.position(2) - 2*ystep;
        end
    end
    
    %disp(['Stimulus position (x,y): ' num2str(stimulus.position)])
end

function [stimx, stimy] = convertNearestPixel(myscreen,x_img,y_img)
    x_screen = x_img*myscreen.screenWidth/myscreen.imageWidth + myscreen.screenWidth/2;
    y_screen = y_img*myscreen.screenHeight/myscreen.imageHeight + myscreen.screenHeight/2;
    
    stimx = (ceil(x_screen)-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
    stimy = (floor(y_screen)-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight; % what is imagewidth???
end


function precomputBackground(myscreen,stimulus)
% generate a lot of noise
% takes about 12 minutes to generate on the stimulus test computer. 
savefile = '/Users/gru/data/trackpos/trackpos.mat';

%% noise
% background noise
downsample_timeRes  = 1;
downsample_spatRes  = 10;

%ntrials             = 1; %7*9; %7 trials * 9 conditions
nframes             = 10000; %myscreen.framesPerSecond*30; % 30s; %/downsample_timeRes; 

%for idx = 1:ntrials
tic
if ~isfield(stimulus,'gaussianFFT')
    xsize_deg          = round(myscreen.imageWidth/downsample_spatRes);
    ysize_deg          = round(myscreen.imageHeight/downsample_spatRes);

    backgaussian = mglMakeGaussian(xsize_deg,ysize_deg,...
        stimulus.stimStd/downsample_spatRes,stimulus.stimStd/downsample_spatRes)*255;
    stimulus.gaussianFFT = getHalfFourier(backgaussian);
end 

for idx2 = 1:nframes
    back                        = stimulus.gaussianFFT; %0.02s
    back.phase                  = rand(size(back.mag))*2*pi; % scramble phase % 0.02s
    backgroundnoise             = round(reconstructFromHalfFourier(back));   %0.04s
    if idx2 == 1 % to save time; only allocate memory once.
        backgroundnoise_rgb         = 255*ones(4,size(backgroundnoise,2),size(backgroundnoise,1),nframes,'uint8'); %0.1165 s
    end % otherwise just overwrite.
    backgroundnoise_rgb(4,:,:,idx2)  = backgroundnoise'/max(max(backgroundnoise))*stimulus.noiseLum;  % normalize contrast %0.025s
end
backgroundnoise_rgb = uint8(backgroundnoise_rgb);
%eval(['backgroundnoise_rgb =uint8(backgroundnoise_rgb);']); %0.02s 
toc
%end
save(savefile, 'backgroundnoise_rgb*','-v7.3')

end

%% utility
function totaldur = approximate_total_task_dur(task)
    % count trials
    numTrials        = 0; 
    totaltime        = 0;
    for phaseNum = 1:length(task{1})
        numTrials = numTrials + task{1}{phaseNum}.numTrials;
        totaltime = totaltime + task{1}{phaseNum}.numTrials * sum(task{1}{phaseNum}.segmax);
    end
    totaldur = totaltime/60/60;
    disp(['Approx task duration = ' num2str(totaldur) ' hours']);
end

function task = add_calculated_params(task)

    for phaseNum = 1:length(task)
        maxframes = ceil(task{phaseNum}.segmax(1)*myscreen.framesPerSecond) + 20;
        task{phaseNum}.randVars.calculated.randomSeed   = nan;
        task{phaseNum}.randVars.calculated.bgpermute    = nan(1,maxframes); % nframes x 1 for the background
        task{phaseNum}.randVars.calculated.perm         = nan(maxframes,1);
        task{phaseNum}.randVars.calculated.initStim     = [nan nan];
        task{phaseNum}.randVars.calculated.trackStim    = nan(maxframes,2);
        task{phaseNum}.randVars.calculated.trackResp    = nan(maxframes,2);
        task{phaseNum}.randVars.calculated.trackEye     = nan(maxframes,2);
        task{phaseNum}.randVars.calculated.trackJoy     = nan(maxframes,4);
        task{phaseNum}.randVars.calculated.trackTime    = nan(1,maxframes);
        task{phaseNum}.randVars.calculated.trackEyeTime = nan(1,maxframes); % for referencing edf file
    end
end