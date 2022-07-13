%        $Id: $
%      usage: trackpos
%         by: Josh Ryu
%       date: 05/05/2021
%    purpose: 

% todo: make noiseLum more flexible
% todo: currently only segment 1 is the tracking period (calculated
% variables stored only here including eye tracking, stimulus and response
% etc... make this more flexible? framecount is a bit awkward too.
% indicate tracking period or 
% todo: only load background once
% todo: move cursor into the task stimuli 
% todo: move things to startsegment. clear screen at inittrial

function myscreen = trackpos(varargin)
%getArgs(varargin,{'subjectID=s999','centerX=10','centerY=0','diameter=16'}); getArgs(varargin,{'subjectID=-1'});
% set up screen
myscreen = struct();
if isempty(mglGetSID)
    myscreen.subjectID  = -1;
else
    myscreen.subjectID  = mglGetSID;
end

%myscreen.screenWidth = 860; myscreen.screenHeight = 600; 
%myscreen.hideCursor = 1;
myscreen                = initScreen(myscreen);

%% experiment parameters
% Experimenter parameters
%todo:  check these throughout the code!!
exp.debug               = 0; % debug code
exp.noeye               = 1; % 1 if no eyetracking (mouse for eye); 0 if there is eye tracking `
exp.showMouse           = 0; % show mouse during everything

exp.fixateCenter        = 1; %
exp.dispPointer         = 1; % display pointer
exp.useJoystick 		= 0; % use joystick. Need "simulink 3D animation" package downloaded. 

exp.downsample_timeRes  = 1; % downsample temporal resolution of background noise the by this factor.
exp.phasescrambleOn     = 1; % load background, if specified by task

exp.grabframe           = 0; % capture frames; todo: change

global stimulus; stimulus = struct;
stimulus.exp = exp;

task = {}; 

%% specify task design
% sb      = stillblob(myscreen, 'backLum', 0, 'maxtrialtime',30,'steady_thresh_frame',50);
% phase1  = sb.configureExperiment(stimulus,task,myscreen);
stimsize = 1;
circ1   = circular(myscreen, 'noiseLum', 32, 'stimLum', 96, 'stimStd', stimsize, 'rStd', [1]);
circ2   = circular(myscreen, 'noiseLum', 32, 'stimLum', 96, 'stimStd', stimsize, 'r_logSpace', true, 'rStd', 0.3);
circ3   = circular(myscreen, 'noiseLum', 32, 'stimLum', 96, 'stimStd', stimsize, 'r_logSpace', false, 'rStd', 0,...
                   'thetaStep',0, 'thetaStd0', pi/10);
circ4   = circular(myscreen, 'noiseLum', 32, 'stimLum', 96, 'stimStd', stimsize, 'r_logSpace', false, 'rStd', 0,...
                   'thetaStep',0, 'thetaStd1', (pi/30)^2);

cp1     = brownian(myscreen, 'noiseLum', 0, 'stimLum', 92, 'stimStd', stimsize);
cp2     = brownian(myscreen, 'noiseLum', 32, 'stimLum', 96, 'stimStd', stimsize);
cp3     = brownian(myscreen, 'noiseLum', 32, 'stimLum', 96, 'stimStd', 0,  'pointLum', 96, 'pointStd', stimsize);
cp4     = brownian(myscreen, 'noiseLum', 32, 'stimLum', 64, 'stimStd', stimsize);
cp5     = brownian(myscreen, 'noiseLum', 32, 'stimLum', 64, 'stimStd', 0,  'pointLum', 64, 'pointStd', stimsize);
cp6     = brownian(myscreen, 'noiseLum', 32, 'stimLum', 0, 'stimStd', 0);
%stimulus.task = {circ1,circ2};
%stimulus.task = {cp2, cp3, cp4, cp5}; %, cp3, cp4, cp5, cp6};

% run the task at multiple radii.
stimulus.task = {};
thetaStd0_r1 = pi/3; % angular velocity at r=1; 
for r = [2,4,6,8,10]
    thetaStd0 = thetaStd0_r1/r; % keep linear velocity constant
    circ   = circular(myscreen, 'noiseLum', 32, 'stimLum', 80, 'stimStd', stimsize, ...
                      'r', [r], 'thetaStd0', thetaStd0, ...
                      'maxtrials',5);
    stimulus.task{end+1} = circ;                   
end

%% configure task
task{1} = cell(length(stimulus.task),1);
for ts = 1:length(stimulus.task)
    task{1}{ts} = stimulus.task{ts}.configureExperiment(task,myscreen,stimulus);
end
task{1} = add_calculated_params(task{1}, myscreen);
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
% todo: make stimulus-specific background
% this is loaded in the beginning and not during task initialization
for phaseN = 1 %:length(task{1})
    currtask = stimulus.task{phaseN};
    if ~isempty(currtask.bgfile) && exp.phasescrambleOn
        disp('Loading phase scrambled background noise...')

        savefile = currtask.bgfile;
        % savefile            = '/Users/joshua/data/trackpos_2afc/trackpos.mat'; % just use noise 1 and permute
        if ~exist(savefile,'file')
            disp(' THE SPECIFIED BACKGROUND FILE DOES NOT EXIST')
            someinput = input('press any button to continue');
            continue
        end

        load(savefile,'backgroundnoise_rgb');

        % create all background textures and then load them later
        if stimulus.exp.debug 
            nnn = 200;
        else
            nnn = size(backgroundnoise_rgb,4);
        end
        for idx = 1:nnn %too big?? memory?
            stimulus.backnoise{phaseN}{idx} = mglCreateTexture(backgroundnoise_rgb(:,:,:,idx));
        end

        clearvars('backgroundnoise_rgb')
        disp('Finished loading phase scrambled background noise...')
    end
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

if stimulus.exp.useJoystick
stimulus.joy = vrjoystick(1); % use simulink 3d animation to load joystick object
if isempty(stimulus.joy)
    stimulus.exp.useJoystick = 0;
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

if ~exp.showMouse, mglDisplayCursor(0);, end %hide cursor

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
    phaseNum = task.thistrial.thisphase; % phaseNum?
       
    % save seed for generating random textures
    rng('shuffle','twister'); 
    s=rng; 
    task.thistrial.randomSeed = s.Seed;
    rng(task.thistrial.randomSeed,'twister');
    
    % noise
    % permute background and save
    if stimulus.exp.phasescrambleOn == 1 && ~isempty(stimulus.task{phaseNum}.bgfile)
        nframes = myscreen.framesPerSecond*task.segmax(1) + 20;%/downsample_timeRes; 
        task.thistrial.bgpermute(1:nframes) = randi(length(stimulus.backnoise{1}),nframes,1);
    end
    
    if mod(task.trialnum,ceil(task.numTrials/20)) == 1
        disp(['(trackpos) '  num2str(task.trialnum/task.numTrials) ...
            '% finished: Trial ' num2str(task.trialnum) ' / ' num2str(task.numTrials)]);
    end
    
    % task initTrial
    stimulus.task{phaseNum}.initTrial(task, myscreen,stimulus);
    
    % initialize position to the stimulus position
    % for joystick and mouse tracking
    stimulus.pointer = stimulus.task{phaseNum}.positions{1}(1:2);
        
    % count frames
    task.thistrial.framecount = 0;
end

%% Start segment
function [task myscreen] = startSegmentCallback(task, myscreen)
% S1: Stimulus (30s)
% S2: Fixation (3s)
global stimulus 
phaseNum = task.thistrial.thisphase;

stimulus.task{phaseNum}.startSegment(task, myscreen, stimulus);

if stimulus.exp.grabframe
    global frame
    %save('/Users/jryu/Dropbox/Stanford/Current/FYP/FYP talk/figures/trackposTask_nonoise_160back.mat', 'frame','-v7.3')
    %save('/Users/jryu/Dropbox/Stanford/Current/FYP/FYP talk/figures/trackposTask_noise_160back.mat', 'frame','-v7.3')
    frame = {};
    frame{task.segmax(1)*myscreen.framesPerSecond} = [];
end

end

%% screen update
function [task, myscreen] = screenUpdateCallback(task, myscreen)

t0 = tic;
%if stimulus.exp.debug, t1 = toc(t0); end

global stimulus
phaseNum = task.thistrial.thisphase;
    
% move cursor        
if stimulus.task{phaseNum}.movecursor 
    % **&display mouse position
    if ~stimulus.exp.useJoystick
        mInfo = mglGetMouse(myscreen.screenNumber);
        stimulus.pointer(1) = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
        stimulus.pointer(2) = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;
    else
        [vx, vy] = joy2vel(stimulus.joy, stimulus.joy_params, myscreen);
        stimulus = update_pointer(stimulus, [vx, vy], myscreen);
    end
end

% update task
task  = stimulus.task{phaseNum}.update(task, myscreen, stimulus);

% display cursor
if stimulus.exp.dispPointer
    mglGluDisk(stimulus.pointer(1), stimulus.pointer(2), 0.2, [1 0 0])
end
if ~stimulus.exp.showMouse, mglDisplayCursor(0);, end %hide cursor

% if we are in the tracking period,  save tracking variables
if any(task.thistrial.thisseg==[1])
    % update framecount
    task.thistrial.framecount = task.thistrial.framecount + 1;

    if stimulus.exp.useJoystick
        task.thistrial.trackJoy(task.thistrial.framecount,:)  = axis(stimulus.joy);
    end

    task.thistrial.trackResp(task.thistrial.framecount,:) = stimulus.pointer;
    task.thistrial.trackStim(task.thistrial.framecount,:) = stimulus.task{phaseNum}.positions{1};
    task.thistrial.trackTime(task.thistrial.framecount)   = mglGetSecs(stimulus.t0); 
    
    % eye tracking
    if (~stimulus.exp.noeye)
        % mouse version for testing with no eyetracker
        [pos,postime] = mglEyelinkGetCurrentEyePos; % is this in image coordinates?
        task.thistrial.trackEye(task.thistrial.framecount,:)    = pos;
        task.thistrial.trackEyeTime(task.thistrial.framecount)  = postime;
    end
end

% grabframe
if stimulus.exp.grabframe && any(task.thistrial.thisseg==[1])
    global frame; frame{task.thistrial.framecount} = mglFrameGrab;
end

end

%% Get response; do nothing. 
function [task myscreen] = responseCallback(task, myscreen)

global stimulus
 
end

%% joystick update
function stimulus = update_pointer(stimulus, vel, myscreen)
    pos = stimulus.pointer;
    
    [horz_out, vert_out] = check_oob(pos + vel, myscreen, stimulus);
    stimulus.pointer(1) = stimulus.pointer(1) + (1-horz_out)*vel(1);
    stimulus.pointer(2) = stimulus.pointer(2) + (1-vert_out)*vel(2);
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


function task = add_calculated_params(task, myscreen)
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