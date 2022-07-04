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
sb      = stillblob(myscreen);
phase1  = sb.configureExperiment(stimulus,task,myscreen);
task{1} = {phase1};
stimulus.task = {sb};

task{1} = add_calculated_params(task{1});
totaldur = approximate_total_task_dur(task);

%% initialize
% intiailize task
disp(' Initializing Task....')

for phaseN = 1:length(task{1})
    [task{1}{phaseN} myscreen] = initTask(task{1}{phaseN},myscreen,...
        @startSegmentCallback,@screenUpdateCallback,@responseCallback,@initTrialCallback);
end

myscreen = initStimulus('stimulus',myscreen); % save the stimulus into stimfile

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
    myscreen = eyeCalibDisp(myscreen); % calibrate eye every time.
    
    % let the user know
    disp(sprintf('(trackpos) Starting Run...'));
end

if stimulus.exp.usejoystick
stimulus.joy = vrjoystick(1); % use simulink 3d animation to load joystick object
if isempty(stimulus.joy)
    stimulus.exp.usejoystick = 0;
    exp = stimulus.exp;
    disp(' FAILED TO FIND JOYSTICK! MAKE SURE SIMULINK 3D ANIMATION PACKAGE IS INSTALLED AND THE JOYSTICK IS PROPERLY CONNECTED');
    disp(' USING MOUSE FOR TRACKING ...');
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
    
    % save seed for generating random textures
    rng('shuffle','twister'); 
    s=rng; 
    task.thistrial.randomSeed = s.Seed;
    rng(task.thistrial.randomSeed,'twister');
    
    % noise
    if stimulus.exp.phasescrambleOn == 1
        nframes = myscreen.framesPerSecond*task.segmax(1) + 20;%/downsample_timeRes; 
        task.thistrial.bgpermute(1:nframes) = randi(length(stimulus.backnoise),nframes,1);
    end
    
    if mod(task.trialnum,ceil(task.numTrials/20)) == 1
        disp(['(trackpos) '  num2str(task.trialnum/task.numTrials) ...
            '% finished: Trial ' num2str(task.trialnum) ' / ' num2str(task.numTrials)]);
    end
    
    %% task initTrial
    % todo: check phasenum variable
    stimulus.task{task.phaseNum}.initTrial(task, myscreen, stimulus);
end

%% Start segment
function [task myscreen] = startSegmentCallback(task, myscreen)
% S1: Stimulus (30s)
% S2: Fixation (3s)
global stimulus 

% todo: check phasenum variable
stimulus.task{task.phaseNum}.startSegment(task, myscreen, stimulus);

if stimulus.exp.grabframe
    global frame
    %save('/Users/jryu/Dropbox/Stanford/Current/FYP/FYP talk/figures/trackposTask_nonoise_160back.mat', 'frame','-v7.3')
    %save('/Users/jryu/Dropbox/Stanford/Current/FYP/FYP talk/figures/trackposTask_noise_160back.mat', 'frame','-v7.3')
    frame = {};
    frame{task.segmax(1)*myscreen.framesPerSecond} = [];
end


end

%% screen update
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus
% update framecount
task.thistrial.framecount = task.thistrial.framecount + 1;

% update task stimulus
task  = stimulus.task{task.phaseNum}.update(obj, task, myscreen, stimulus);

% move cursor        
if task.thistrial.movecursor
    % **&display mouse position
    if ~stimulus.exp.usejoystick
        mInfo = mglGetMouse(myscreen.screenNumber);
        stimulus.pointer(1) = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
        stimulus.pointer(2) = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;
    else
        [vx, vy] = joy2vel(stimulus.joy, stimulus.joy_params, myscreen);
        stimulus = update_pointer(stimulus, [vx, vy], myscreen);
        task.thistrial.trackJoy(task.thistrial.framecount,:)  = axis(stimulus.joy);
    end
    
    mglGluDisk(stimulus.pointer(1), stimulus.pointer(2), 0.1, [1 0 0])

    task.thistrial.trackStim(task.thistrial.framecount,:) = stimulus.position;
    task.thistrial.trackResp(task.thistrial.framecount,:) = [mimg_x, mimg_y];
    task.thistrial.trackTime(task.thistrial.framecount)   = mglGetSecs(stimulus.t0); 
end

% eye tracking
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
        
    task.thistrial.trackEye(task.thistrial.framecount,:)    = pos;
    task.thistrial.trackEyeTime(task.thistrial.framecount)  = postime;
end

% grabframe
if stimulus.exp.grabframe && (task.thistrial.thisseg== 1)
    global frame; frame{task.thistrial.framecount} = mglFrameGrab;
end

end

%% Get response; do nothing. 
function [task myscreen] = responseCallback(task, myscreen)

global stimulus
 
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