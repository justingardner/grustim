%
%        $Id: $
%      usage: trackpos_est
%         by: Josh Ryu
%       date: 05/06/2021
%    purpose: estimate position of blob
%
% staircase not implemented yet => takes staircase from afc runs.
% consider: adding another noise period after the stimulus

% S1: wait time while the main task calls this subtask.
% S2: random period of fixation (random ~0.5s)
% S3: stimulus period (stimdur s)
% S4: delay to mask
% S5: mask
% S6: repsonse period (inf)
% S7: feedback (1s)

% fixed position difference or samples from continuous distributions?

function [task, myscreen] = trackpos_sub_est(myscreen, params, exp)
%% task parameters
% stimulus and background
task{1}.random               = 1;

fieldnames = fields(params.task);
for fidx = 1:length(fieldnames) 
    fieldname = fieldnames{fidx};
    task{1}.parameter.(fieldname) = params.task.(fieldname);
end
task{1}.parameter.stimright  = [0,1]; % or polar angle?

% todo: change for no feedback and no mask
if mglIsFile(exp.noise_mask)
    task{1}.segmin           = [inf 0.2 0.5 inf inf inf inf 1];
    task{1}.segmax           = [inf 0.2 0.5 inf inf inf inf 1]; 
    task{1}.getResponse      = [0 0 0 0 0 0 1 0]; %segment to get response.
else
    % basically skip mask period
    task{1}.segmin           = [0.1 0.4 inf inf 1];
    task{1}.segmax           = [0.1 0.8 inf inf 1]; 
    task{1}.getResponse      = [0 0 0 1 0]; %segment to get response.
end

global stimulus
if strcmp(stimulus.exp.est.presSched, 'gaussian')
    task{1}.randVars.calculated.posDiff  = nan; 
elseif strcmp(exp.est.presSched, 'fixed')
    task{1}.parameter.posDiff            = params.task.posDiff; 
else
    task{1}.parameter.posDiff            = params.task.posDiff; 
end

if exp.block_design % with some overflow.
    task{1}.numBlocks        = trialpercond; % dont count stimright as condition %with some overflow
else
    task{1}.numTrials        = params.numTrials; % dont count stimright as condition %with some overflow
end

% calculated variables
% add response period?
trialdur = task{1}.segmax(2)+task{1}.segmax(3)+max(params.task.stimDur);
% if mglIsFile(exp.noise_mask)
%     trialdur = trialdur + max(params.mask_TOff2MOn) + max(params.maskDur);
% end
maxframes = ceil(trialdur)*myscreen.framesPerSecond+10; % with some additional overflow

task{1}.randVars.calculated.stimDur0     = nan; % stimulus duratino prescribed by the experimenter
task{1}.randVars.calculated.pos_estim    = nan(1,2); % nframes x 2 for position estimate

task{1}.randVars.calculated.bgpermute    = nan(1,maxframes); % nframes x 1 for the background
task{1}.randVars.calculated.stimON       = nan(1,maxframes); % nframes x 1 for the stimulus
   
task{1}.randVars.calculated.segTime     = nan(1,length(task{1}.segmin));

task{1}.randVars.calculated.trackTime    = nan(1,maxframes);
task{1}.randVars.calculated.trackEye     = nan(maxframes,2);
task{1}.randVars.calculated.trackEyeTime = nan(1,maxframes); % for referencing edf file

%% initializing task...
disp(' Initializing Task....')

for phaseN = 1:length(task)
    [task{phaseN}, myscreen] = initTask(task{phaseN},myscreen,...
        @startSegmentCallback,@screenUpdateCallback,@responseCallback,@initTrialCallback);
end

function [task, myscreen] = initTrialCallback(task, myscreen)

global stimulus

backLum     = task.thistrial.backLum;
noiseLum    = task.thistrial.noiseLum;
stimLum     = task.thistrial.stimLum;
stimDur     = task.thistrial.stimDur;
stimStd     = task.thistrial.stimStd;
stimColor   = task.thistrial.stimColor;

task.thistrial.seglen(4) = task.thistrial.stimDur;

if strcmp(stimulus.exp.est.presSched, 'staircase')
    idx         = findCondIdx(stimulus.staircaseTable,task.thistrial);
    stimulus.staircaseIdx = idx;
    [s, staircase_no_update] = doStaircase('gettestvalue',stimulus.staircaseTable.staircase{idx});
    task.thistrial.posDiff = s;
elseif strcmp(stimulus.exp.est.presSched, 'gaussian')
    task.thistrial.posDiff = normrnd(stimulus.params_est.prior.mean,stimulus.params_est.prior.std);
end

% noise mask
if mglIsFile(stimulus.exp.noise_mask)   
    task.thistrial.seglen(5) = task.thistrial.maskDur;
    task.thistrial.seglen(6) = task.thistrial.mask_TOff2MOn;
    maskLum = task.thistrial.maskLum;
    nframes = myscreen.framesPerSecond*task.thistrial.seglen(6) + 20; %/downsample_timeRes; 
    stimulus.noise_mask_trial = randi(size(stimulus.noise_mask.backgroundnoise_rgb,4),nframes,1); % sample with replacement
    
    % delete texture
    if isfield(stimulus,'noise_mask_texture') 
       for idx = 1:length(stimulus.noise_mask_texture)
           mglDeleteTexture(stimulus.noise_mask_texture{idx});
       end
    end

    % create texture
    for idx = 1:nframes
        midx    = stimulus.noise_mask_trial(idx);
        maskimg = stimulus.noise_mask.backgroundnoise_rgb(:,:,:,midx);
        maskimg(4,:,:) = maskLum * maskimg(4,:,:);
        % maskimg(4,:,:) = 0 * maskimg(4,:,:);
        maskimg = permute(maskimg,[2,3,1]);
        stimulus.noise_mask_texture{idx} = mglMetalCreateTexture(maskimg);
    end
end

%% Start segment
function [task, myscreen] = startSegmentCallback(task, myscreen)

global stimulus

% set flushMode based on noiseLum
if stimulus.exp.phasescrambleOn == 1 && stimulus.exp.backprecompute == 1 && task.thistrial.noiseLum >0
    myscreen.flushMode = 0;
else
    myscreen.flushMode = 1;
end

if task.thistrial.thisseg == 1
    % if we come back to beginning, means the previous trial finished. 
    % wait for the next task 
    myscreen.flushMode = 0; % update screen
    stimulus.currtask = 'done'; 
elseif task.thistrial.thisseg == 2
    if ~stimulus.exp.noeye
        myscreen.flushMode = 0;
    end
    % start the task.
    if ~stimulus.exp.showmouse, mglDisplayCursor(0);, end 
        
    stimulus.lum        = task.thistrial.stimLum;
    stimulus.std        = task.thistrial.stimStd;    
    stimulus.color      = task.thistrial.stimColor;    
    stimulus.backLum    = task.thistrial.backLum;
    stimulus.noiseLum   = task.thistrial.noiseLum;
    
    task.thistrial.framecount = 0;
    task.thistrial.stimDur0    = task.thistrial.seglen(4); % stimulus period
    
    stimulus.target = trackposInitStimulus(stimulus,myscreen); %centerX,Y, diameter called by getArgs.

    if stimulus.exp.phasescrambleOn == 1 && stimulus.exp.backprecompute == 1 && stimulus.noiseLum >0
        task.thistrial.framecount = 0;
        nframes = length(task.thistrial.bgpermute);
        task.thistrial.bgpermute(1:nframes) = randi(length(stimulus.backnoise),nframes,1);
    end
elseif task.thistrial.thisseg == 4
    if ~stimulus.exp.noeye
        myscreen.flushMode = 0;
    end
elseif task.thistrial.thisseg == 5
    task.thistrial.framecount = 0; % restart framecount
    task.thistrial.stimDur = task.thistrial.seglen(4); % stimulus period recorded by updateTask; make sure is same as Dur0
elseif task.thistrial.thisseg == 6
    task.thistrial.framecount = 0; % restart framecount
    myscreen.flushMode = 0; % has to be 0 for mask
elseif task.thistrial.thisseg == 7
    myscreen.flushMode = 0; % has to be 0 for response
    % response period
    % center the mouse 
    x_img = 0;  y_img = 0;
    x_screen = x_img*myscreen.screenWidth/myscreen.imageWidth + myscreen.screenWidth/2;
    y_screen = y_img*myscreen.screenHeight/myscreen.imageHeight + myscreen.screenHeight/2;
    mglSetMousePosition(x_screen,y_screen, myscreen.screenNumber);
end

% blt screen once before screenUpdates loops
if task.thistrial.thisseg > 1
    [task, myscreen] = screenUpdateCallback(task, myscreen);
    mglFlush;
    if task.thistrial.thisseg == 4
        stimulus.start = mglGetSecs;
    elseif task.thistrial.thisseg == 5
        stimulus.length = mglGetSecs - stimulus.start;
        disp(['Segment duration error: ', num2str(stimulus.length - task.thistrial.stimDur)])
    end
end

task.thistrial.segTime(task.thistrial.thisseg) = mglGetSecs;
    
%% screen update
function [task, myscreen] = screenUpdateCallback(task, myscreen)
global stimulus  % call stimulus

if task.thistrial.thisseg== 1
    %% waiting for task to start
    if strcmp(stimulus.currtask,'est')
        stimulus.currtask = 'running est';
        task = jumpSegment(task); 
    end
else
    %% do the task
    % set background luminance
    if task.thistrial.backLum > 1
        mglClearScreen(task.thistrial.backLum/255);
    else
        mglClearScreen(task.thistrial.backLum);
    end

    task.thistrial.framecount = task.thistrial.framecount + 1;
    task.thistrial.stimON(task.thistrial.framecount) = 0; %count stimulus

    % inject noise, track time
    if any(task.thistrial.thisseg == [2,3,4]) 
        if stimulus.exp.phasescrambleOn == 1 
            idx = task.thistrial.bgpermute(task.thistrial.framecount);
            mglBltTexture(stimulus.backnoise{idx},...
                [0 0 myscreen.imageWidth myscreen.imageHeight])
        end

        task.thistrial.trackTime(task.thistrial.framecount) = mglGetSecs(stimulus.t0);
    end

    % draw blob or mask 
    if task.thistrial.thisseg == 4 % stimulus
        stim_pos = (2*task.thistrial.stimright-1)*task.thistrial.posDiff;
        task.thistrial.stimON(task.thistrial.framecount) = 1;
        mglMetalBltTexture(stimulus.target.img,[stim_pos 0]);
    elseif task.thistrial.thisseg == 6
        mglMetalBltTexture(stimulus.noise_mask_texture{task.thistrial.framecount},[0 0]);
    end

    % add fixation
    if task.thistrial.thisseg == 2
        mglMetalArcs([0;0;0], [stimulus.fixColors.stim'; 0.5], [0.3;0.5],[0;2*pi], 1);
    end

    if any(task.thistrial.thisseg == [2,3,4])
        if stimulus.exp.colorfix
            % changing fixation colors
            % mglGluAnnulus(0,0,0.2,0.3,stimulus.fixColor,60,1);
            % mglGluDisk(0, 0, 0.1, 0.5+0.5*rand(1,3),60,1); 
            mglMetalDots([0;0;0], [0.5+0.5*rand(3,1);1], [0.1;0.1], 1, 1);
        else
            % mglGluDisk(0, 0, 0.1, stimulus.fixColors.stim,60,1); 
            mglMetalDots([0;0;0], [stimulus.fixColors.stim';1], [0.1;0.1], 1, 1);
        end
    elseif any(task.thistrial.thisseg == [5,6,7,8])
        % fixation indicating estimation task
        % mglGluDisk(0, 0, 0.1, stimulus.fixColors.est,60,1); 
        mglMetalDots([0;0;0], [stimulus.fixColors.est';1], [0.1;0.1], 1, 1);
    end


    % response and feedback
    if task.thistrial.thisseg == 7
        % show cursor
        mInfo = mglGetMouse(myscreen.screenNumber);
        degx = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
        degy = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;

        if ~stimulus.exp.estim_horiz, degx = 0; end
        if ~stimulus.exp.estim_verti, degy = 0; end
        % mglGluDisk(degx, degy, 0.1, [1 0 0]);
        mglMetalDots([degx;degy;0], [1;0;0;1], [0.1;0.1], 1, 1);
    elseif stimulus.exp.feedback && task.thistrial.thisseg == 8% feedback period
        stim_pos = (2*task.thistrial.stimright-1)*task.thistrial.posDiff;
        % mglGluDisk(stim_pos, 0, 0.1, [1 0 0]) ;    % draw center of blob
        mglMetalDots([stim_pos;0;0], [stimulus.fixColors.fb';1], [0.1;0.1], 1, 1);
    end

    % track eye
    if (~stimulus.exp.noeye) && any(task.thistrial.thisseg==[2,3,4]) 
        % mouse version for testing with no eyetracker
        if stimulus.exp.eyemousedebug
            mInfo = mglGetMouse(myscreen.screenNumber);
            degx = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
            degy = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;

            pos = [degx, degy];
        else  % check eye pos
            [pos,postime] = mglEyelinkGetCurrentEyePos; % is this in image coordinates?
        end

        task.thistrial.trackEye(task.thistrial.framecount,:)    = pos;
        task.thistrial.trackEyeTime(task.thistrial.framecount)  = postime;
    end
end


%% Get response 
function [task, myscreen] = responseCallback(task, myscreen)

global stimulus

% if the button 3 is pressed, record position
if task.thistrial.whichButton == 3
    % record position of the mouse
    mInfo = mglGetMouse(myscreen.screenNumber);
    degx = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
    degy = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;
    
    if ~stimulus.exp.estim_horiz, degx = 0; end
    if ~stimulus.exp.estim_verti, degy = 0; end
    task.thistrial.pos_estim = [degx, degy]; % save position estimates
    
    stim_pos = task.thistrial.posDiff;

    disp(['(subtask_est) Stimulus Hor. Pos.: ' num2str(stim_pos) '; ' ...
            'Response: ' num2str(degx)  ...
            '; Error: ' num2str(degx - stim_pos)]);
    
    % go to next segment
    task = jumpSegment(task);
end
