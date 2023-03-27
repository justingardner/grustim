%        $Id: $
%      usage: 
%         by: Josh Ryu
%       date: 03/07/2023
%    purpose: 

% Tracking task for circular motion around an circle/ellipse
% Run in mglmetal 

% task.thistrial.trackStim and trackResp is in dva: r * polar_angle 

function myscreen = trackpos_circular(varargin)
%getArgs(varargin,{'subjectID=s999','centerX=10','centerY=0','diameter=16'}); getArgs(varargin,{'subjectID=-1'});

myscreen = setup_screen_jryu(); 
myscreen = initScreen(myscreen);
mglMetalSetViewColorPixelFormat(4);     % set to argb2101010 pixel format
rng(0, 'twister'); % set seed

%% experiment parameters
% Experimenter parameters

exp.debug               = 1; % debug code
exp.trackEye            = 1; % 0 if no eyetracking; 1 if there is eye tracking `
exp.showMouse           = 0; % show mouse during everything

exp.showRing            = 0; % show ring
exp.fixateCenter        = 1; % fixate center
exp.controlMethod       = 'wheel'; % available: wheel

exp.grabframe           = 0; % capture frames. specify save directory

global stimulus; stimulus = struct;
stimulus.exp = exp;

%% specify task design

% no noise run
cps                 = {};
stimStdList         = [1]; %[0.5, 1 ,2]; % size of gaussian blob
stim_noiseStdList   = [0]; % in dva per second
stimLums            = [0.2, 0.8]; %[0.1, 0.2, 0.5]; 
% backLum             = 0.7;

ecc_r_list          = [3, 7, 10]; % eccentricity
% ecc_a             = 1; % major axis
% ecc_b             = 1; % minor axis

stim_dyngroup       = [10]; % noise order, same size as stimStdList % 10: constant velocity
stim_vel            = [10]; 

pointStd            = 0.2; stimulus.pointerR = pointStd;
point_noiseStd      = 0;

ntrial_learn        = 3;  % learning phase at full luminance, not analyzed
ntrials             = 15; % trials per condition
nblocks             = 3;  % should divide ntrials, divide trial into blocks

maxtrialtime        = 25; % seconds

if exp.debug, ntrial_learn= 1; ntrials = 1; nblocks = 1; maxtrialtime=5; end

for ecc_r = ecc_r_list
for s = 1:length(stim_noiseStdList)
    stim_noiseStd = stim_noiseStdList(s);
    stimdyngroup = stim_dyngroup(s);

    for stimStd = stimStdList
        % learning phase
        cps{end+1} = circular(myscreen, 'numTrials', ntrial_learn, 'maxtrialtime', maxtrialtime, 'ecc_r', ecc_r, ...
            'stimLum', 1, 'stimStd', stimStd, 'stim_dyngroup', stimdyngroup, 'stim_noiseStd', stim_noiseStd, ...
            'stim_vel', stim_vel,...
            'pointLum',1, 'pointStd', pointStd, 'point_noiseStd', point_noiseStd);
    
        % tracking
        ntrials_phase = length(stimLums) * ceil(ntrials/nblocks);
        for b = 1:nblocks
            cps{end+1} = circular(myscreen, 'numTrials', ntrials_phase, 'maxtrialtime', maxtrialtime, 'ecc_r', ecc_r, ...
                'stimLum', stimLums, 'stimStd', stimStd, 'stim_dyngroup', stimdyngroup, 'stim_noiseStd', stim_noiseStd, ...
                'stim_vel', stim_vel,...
                'pointLum',1, 'pointStd', pointStd, 'point_noiseStd', point_noiseStd);
        end
    end
end
end  
stimulus.task = cps;

%% configure task
task{1} = cell(length(stimulus.task),1);
for ts = 1:length(stimulus.task)
    task{1}{ts} = stimulus.task{ts}.configureExperiment(task,myscreen,stimulus);
end
task{1}     = add_calculated_params(task{1}, myscreen);
totaldur    = approximate_total_task_dur(task);

%% initialize
% intiailize task
disp(' Initializing Task....')

for phaseN = 1:length(task{1})
    [task{1}{phaseN} myscreen] = initTask(task{1}{phaseN},myscreen,...
        @startSegmentCallback,@screenUpdateCallback,@responseCallback,@initTrialCallback);
end

myscreen = initStimulus('stimulus',myscreen); % save the stimulus into stimfile

%% grabframe
if stimulus.exp.grabframe
    global frame
    frame = {};
end

%% Eye calibration and check joystick, wheels
if strcmp(exp.controlMethod,'eye') && ~stimulus.exp.trackEye
    disp(' Need to track eye..  setting  exp.trackEye to true')
    stimulus.exp.trackEye = true;
    exp.trackEye = true;
end

if strcmp(exp.controlMethod,'eye')
    [positions_target, positions_eye] = calibrateEye(myscreen, stimulus, true);
    eyecalib = struct();
    eyecalib.target{1}      = positions_target;
    eyecalib.eye{1}         = positions_eye;
    stimulus.eyecalib       = eyecalib;
elseif stimulus.exp.trackEye
    disp(' Calibrating Eye ....')
    myscreen  = eyeCalibDisp(myscreen); % calibrate eye every time.
end

if strcmp(stimulus.exp.controlMethod, 'wheel')
    stimulus = calibrateWheel(myscreen, stimulus);
end

if strcmp(stimulus.exp.controlMethod, 'joystick')
    stimulus = calibrateJoy(myscreen, stimulus);
end

if ~strcmp(stimulus.exp.controlMethod, 'mouse') && ...
        ~strcmp(stimulus.exp.controlMethod, 'eye') && ...
        ~strcmp(stimulus.exp.controlMethod, 'joystick') && ...
        ~strcmp(stimulus.exp.controlMethod, 'wheel')
    
    disp('control method not found! using mouse for tracking')
    stimulus.exp.controlMethod = 'mouse';
end

exp = stimulus.exp;

%% run the task
disp(' Running Task....'); stimulus.t0 = mglGetSecs; % 

% let the experimentee know too...
% mglClearScreen(task{1}{1}.parameter.backLum/255);
%mglFlush;

mglBltTexture(mglText('task starting... '), [0 1.5]);
mglBltTexture(mglText('1. Track the target with the red pointer'),[0 0.5]);
if stimulus.exp.fixateCenter
    mglBltTexture(mglText('2. Fixate in center. You should be able to see the dot changing colors'),[0 -0.5]);    
end
mglBltTexture(mglText('When you are ready, press backtick to go to next trial'),[0 -1.5]);
mglFlush(); myscreen.flushMode = -1;

if ~exp.showMouse, mglDisplayCursor(0);, end %hide cursor

phaseNum = 1;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    [task{1}, myscreen, newphaseNum] = updateTask(task{1},myscreen,phaseNum); % update the task
    if newphaseNum ~= phaseNum
        mglClearScreen(0.5);
        
        if strcmp(exp.controlMethod,'eye')
            [positions_target, positions_eye]       = calibrateEye(myscreen, stimulus, true);
            stimulus.eyecalib.target{newphaseNum}   = positions_target;
            stimulus.eyecalib.eye{newphaseNum}       = positions_eye;
        end
        
        mglBltTexture(mglText('When you are ready, press backtick to go to next trial'),[0 -1.5]);
        mglFlush(); myscreen.flushMode = -1;
    end
    phaseNum = newphaseNum;
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
endScreen(myscreen); mglDisplayCursor(1) % show cursor

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
    
    % task status report
    if mod(task.trialnum,ceil(task.numTrials/20)) == 1
        disp(['(trackpos) '  num2str(task.trialnum/task.numTrials) ...
            '% finished: Trial ' num2str(task.trialnum) ' / ' num2str(task.numTrials)]);
    end  
    
    % these should be updated by the dynamics function 
    stimulus.target = struct();
    stimulus.pointer = struct();
    stimulus.otherObjs = {};
    
    % task initTrial
    [task, stimulus] = stimulus.task{phaseNum}.initTrial(task, myscreen, stimulus);
end


%% Start segment
function [task myscreen] = startSegmentCallback(task, myscreen)
    myscreen.flushMode = 0;
    global stimulus 
    phaseNum = task.thistrial.thisphase;
   
    stimulus = stimulus.task{phaseNum}.startSegment(task, myscreen, stimulus);
    
    if stimulus.task{phaseNum}.doTrack
        task.thistrial.framecount = 1;
    else
        task.thistrial.framecount = [];
    end

    if stimulus.exp.grabframe
        global frame
        frame = {};
        frame{task.segmax(1)*myscreen.framesPerSecond} = [];
    end
end


%% screen update
function [task, myscreen] = screenUpdateCallback(task, myscreen)
    global stimulus
    phaseNum = task.thistrial.thisphase;
    
    % set background luminance
    if task.thistrial.backLum > 1
        mglClearScreen(task.thistrial.backLum/255);
    else
        mglClearScreen(task.thistrial.backLum);
    end

    % draw blue ring for the trajectory path
    if stimulus.exp.showRing
        mglMetalRing_wlines(task.thistrial.ecc_r, stimulus.pointerR/2, [0,0,1,1], 600)
    end

    % eye tracking
    if stimulus.exp.trackEye && stimulus.task{phaseNum}.doTrack
        % mouse version for testing with no eyetracker
        [pos,postime] = mglEyelinkGetCurrentEyePos; % is this in image coordinates?
        task.thistrial.trackEye(task.thistrial.framecount,:)    = pos;
        task.thistrial.trackEyeTime(task.thistrial.framecount)  = postime;
    end

    % update task
    [task, stimulus]  = stimulus.task{phaseNum}.update(task, myscreen, stimulus);
    
    % if we are in the tracking period,  save tracking variables
    if stimulus.task{phaseNum}.doTrack
        % update framecount
        task.thistrial.framecount = task.thistrial.framecount + 1;

        % display target
        if ~isempty(stimulus.target)
            mglBltTexture(stimulus.target.img, stimulus.target.position);
        end

        % display other objects
        for ij = 1:length(stimulus.otherObjs)
            mglBltTexture(stimulus.otherObjs{ij}.img, stimulus.otherObjs{ij}.position);
        end
    end

        
    % display pointer
    if ~isempty(stimulus.pointer)
        if isfield(stimulus.pointer, 'img') && isempty(stimulus.pointer.img)
            mglBltTexture(stimulus.pointer.img, stimulus.pointer.position);
        else
            mglMetalDots([stimulus.pointer.position(1);stimulus.pointer.position(2);0], ...
             [1;0;0;1], [stimulus.pointer.std; stimulus.pointer.std], 1, 1);
        end
    end
    
    % display fixation
    if stimulus.exp.fixateCenter == 1 && stimulus.task{phaseNum}.displayFix % fixation below others.
        mglMetalArcs([0;0;0], [1;1;1; 1], [stimulus.pointerR+0.1;stimulus.pointerR+0.3],[0;2*pi], 1);
        mglMetalDots([0;0;0], [0.5+0.5*rand(3,1);1], [stimulus.pointerR;stimulus.pointerR], 1, 1);
    end
    
    if ~stimulus.exp.showMouse, mglDisplayCursor(0);, end %hide cursor

    % grabframe
    if stimulus.exp.grabframe && any(task.thistrial.thisseg==stimulus.task{phaseNum}.doTrack)
        global frame; frame{task.thistrial.framecount} = mglFrameGrab;
    end
end

%% Get response; do nothing. 
function [task myscreen] = responseCallback(task, myscreen)
    global stimulus
 
end

function task = add_calculated_params(task, myscreen)
    for phaseNum = 1:length(task)
        maxframes = ceil(task{phaseNum}.segmax(1)*myscreen.framesPerSecond) + 20;
        task{phaseNum}.randVars.calculated.randomSeed   = nan;
        task{phaseNum}.randVars.calculated.bgpermute    = nan(1,maxframes); % nframes x 1 for the background
        task{phaseNum}.randVars.calculated.initStim     = [nan nan];
        task{phaseNum}.randVars.calculated.trackStim    = nan(1,maxframes);
        task{phaseNum}.randVars.calculated.trackResp    = nan(1,maxframes);
        task{phaseNum}.randVars.calculated.trackEye     = nan(maxframes,2);
        task{phaseNum}.randVars.calculated.trackTime    = nan(1,maxframes);
        task{phaseNum}.randVars.calculated.trackEyeTime = nan(1,maxframes); % for referencing edf file
    end
end