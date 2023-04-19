%        $Id: $
%      usage: 
%         by: Josh Ryu
%       date: 03/07/2023
%    purpose: 

% Tracking task for circular motion around an circle/ellipse
% Run in mglmetal 

% task.thistrial.trackStim and trackResp is in dva: r * polar_angle 

function myscreen = trackpos_circular_ar(varargin)
%getArgs(varargin,{'subjectID=s999','centerX=10','centerY=0','diameter=16'}); getArgs(varargin,{'subjectID=-1'});

myscreen = setup_screen_jryu(); 
myscreen = initScreen(myscreen);
mglMetalSetViewColorPixelFormat(4);     % set to argb2101010 pixel format
rng(0, 'twister'); % set seed

%% experiment parameters
% Experimenter parameters

exp.debug               = 0; % debug code
exp.trackEye            = 1; % 0 if no eyetracking; 1 if there is eye tracking `
exp.showMouse           = 0; % show mouse during everything

exp.showRing            = 0; % show ring
exp.fixateCenter        = 1; % fixate center
exp.controlMethod       = 'wheel'; %'wheel'; % available: wheel

exp.grabframe           = 0; % capture frames. specify save directory

% for loading wheel
exp.lastStimFile        = ''; %'/Users/jryu/Dropbox/GardnerLab/data/trackpos_circular/test/230329_stim01.mat';

global stimulus; stimulus = struct;
stimulus.exp = exp;

%% specify task design

% cps = load_experiment(myscreen, exp, 1); % AR structure
cps = load_experiment(myscreen, exp, 2); % eccentricity
stimulus.task = cps;
stimulus.fixation_size = 0.4;

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
    eyecalib                = struct();
    eyecalib.target{1}      = positions_target;
    eyecalib.eye{1}         = positions_eye;
    stimulus.eyecalib       = eyecalib;
elseif stimulus.exp.trackEye
    disp(' Calibrating Eye ....')
    % http://sr-research.jp/support/manual/EyeLink%20Programmers%20Guide.pdf
    myscreen  = eyeCalibDisp(myscreen); % calibrate eye every time.
end 

if strcmp(stimulus.exp.controlMethod, 'wheel') || strcmp(stimulus.exp.controlMethod, 'mouse_circ')
    if mglIsFile(exp.lastStimFile)
        a = load(exp.lastStimFile);
        if isfield(a.stimulus,'wheel_params')
            stimulus.wheel_params = a.stimulus.wheel_params;
        end
    else
        disp('Wheel calibration information not found. Reinitializing... ')
    end
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

%% configure task
task{1} = cell(length(stimulus.task),1);
for ts = 1:length(stimulus.task)
    task{1}{ts} = stimulus.task{ts}.configureExperiment(task,myscreen,stimulus);
end
totaldur    = approximate_total_task_dur(task);

%% initialize
% intiailize task
disp(' Initializing Task....')

for phaseN = 1:length(task{1})
    [task{1}{phaseN} myscreen] = initTask(task{1}{phaseN},myscreen,...
        @startSegmentCallback,@screenUpdateCallback,@responseCallback,@initTrialCallback);
end

myscreen = initStimulus('stimulus',myscreen); % save the stimulus into stimfile

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
        
        mglBltTexture(mglText('When you are ready, press backtick to go to next trial'),[0 -1.5]);
        mglFlush(); myscreen.flushMode = -1;

        % todo make this work
%         if myscreen.userHitSpace
%             disp(' Calibrating Eye ....')
%             myscreen  = eyeCalibDisp(myscreen); % calibrate eye again
%         end

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
        mglMetalRing_wlines(task.thistrial.ecc_r, stimulus.pointer.std/2, [0,0,1,1], 600)
    end

    % update task
    [task, stimulus]  = stimulus.task{phaseNum}.update(task, myscreen, stimulus);
    
    % if we are in the tracking period,  save tracking variables
    if stimulus.task{phaseNum}.doTrack
        % update framecount
        task.thistrial.framecount = task.thistrial.framecount + 1;

        % display target
        if ~isempty(stimulus.target) 
            if isfield(stimulus.target, 'img') && ~isempty(stimulus.target.img)
                mglBltTexture(stimulus.target.img, stimulus.target.position);
            else
                mglMetalDots([stimulus.target.position(1);stimulus.target.position(2);0], ...
                    [stimulus.target.color;1], [stimulus.target.std; stimulus.target.std], 1, 1);
            end
        end

        % display other objects
        for ij = 1:length(stimulus.otherObjs)
            mglBltTexture(stimulus.otherObjs{ij}.img, stimulus.otherObjs{ij}.position);
        end
    end

        
    % display pointer
    if ~isempty(stimulus.pointer)
        if isfield(stimulus.pointer, 'img') && ~isempty(stimulus.pointer.img)
            mglBltTexture(stimulus.pointer.img, stimulus.pointer.position);
        else
            mglMetalDots([stimulus.pointer.position(1);stimulus.pointer.position(2);0], ...
             [stimulus.pointer.color;1], [stimulus.pointer.std; stimulus.pointer.std], 1, 1);
        end
    end
    
    % display fixation
    if stimulus.exp.fixateCenter == 1 && stimulus.task{phaseNum}.displayFix % fixation below others.
        mglMetalArcs([0;0;0], [1;1;1; 1], [stimulus.fixation_size+0.1;stimulus.fixation_size+0.3],[0;2*pi], 1);
        mglMetalDots([0;0;0], [0.5+0.5*rand(3,1);1], [stimulus.fixation_size;stimulus.fixation_size], 1, 1);
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

function cps = load_experiment(myscreen, exp, exp_num)
    cps = {};
    ntrial_learn        = 4;  % learning phase at full luminance, not analyzed
    ntrials             = 15; % trials per condition
    trials_per_block    = 5;  % should divide ntrials, divide trial into blocks
    maxtrialtime        = 20; % seconds

    if exp.debug, ntrial_learn= 0; ntrials = 1; trials_per_block = 1; maxtrialtime=10; end

    if exp_num == 1 % effect of pointer dynamics
        experiment          = {'mn'};
        experiment_paramset = [6, 5, 4, 3, 2, 1]; %., 9, 8, 7, 12,11,10];
        Nconds = length(experiment_paramset);

        for epset = experiment_paramset
            % learning phase -- not analyzed
            cps{end+1} = circular_ar(myscreen, 'numTrials', ntrial_learn, 'maxtrialtime', maxtrialtime, ...
                'experiment', experiment, 'experiment_paramset', epset);

            % tracking
            nblocks = ceil(ntrials/trials_per_block);
            for b = 1:nblocks
                cps{end+1} = circular_ar(myscreen, 'numTrials', trials_per_block, 'maxtrialtime', maxtrialtime, ...
                    'experiment', experiment, 'experiment_paramset', epset);
            end
        end
    elseif exp_num == 2 % effect of eccentricity
        experiment          = {'ecc'};
        
        experiment_paramset = [4,1,5,2,6,3]; %., 9, 8, 7, 12,11,10];
        Nconds = length(experiment_paramset);
        
        for epset = experiment_paramset
            % learning phase -- not analyzed
            cps{end+1} = circular_ar(myscreen, 'numTrials', ntrial_learn, 'maxtrialtime', maxtrialtime, ...
                'experiment', experiment, 'experiment_paramset', epset);

            % tracking
            nblocks = ceil(ntrials/trials_per_block);
            for b = 1:nblocks
                cps{end+1} = circular_ar(myscreen, 'numTrials', trials_per_block, 'maxtrialtime', maxtrialtime, ...
                    'experiment', experiment, 'experiment_paramset', epset);
            end
        end
    end
end