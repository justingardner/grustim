%        $Id: $
%      usage: 
%         by: Josh Ryu
%       date: 03/07/2023
%    purpose: 

% Tracking task for circular motion around an circle/ellipse
% Run in mglmetal 

% updates task and displays object on screen

% task.thistrial.trackStim and trackResp is in dva: r * polar_angle 

function myscreen = trackpos_corrdyn_sl(varargin)
%getArgs(varargin,{'subjectID=s999','centerX=10','centerY=0','diameter=16'}); getArgs(varargin,{'subjectID=-1'});

myscreen = setup_screen_jryu(); 
myscreen = initScreen(myscreen);
mglMetalSetViewColorPixelFormat(4);     % set to argb2101010 pixel format
rng(0, 'twister'); % set seed

%% experiment parameters
% Experimenter parameters

exp.debug               = 0; % debug code
exp.trackEye            = 0; % 0 if no eyetracking; 1 if there is eye tracking `
exp.trackEye_calibtrial = 1;
exp.showMouse           = 0; % show mouse during everything

exp.showRing            = 0; % show ring
exp.fixateCenter        = 1; % fixate center
exp.controlMethod       = 'wheel'; %'wheel'; % available: wheel

exp.grabframe           = 0; % capture frames. specify save directory

% for loading wheel and dynamics noise stimulus if needed
exp.lastStimFile        = ''; %'/Users/jryu/Dropbox/GardnerLab/data/trackpos_circular/test/230329_stim01.mat';

% define set of random colors
exp.randColorsFile      = '/Users/gru/proj/grustim/trackpos/util/labcolors.mat'; 
%'/Users/jryu/proj/grustim/trackpos/util/labcolors.mat';
%'/Users/gru/proj/grustim/trackpos/util/labcolors.mat'; 

global stimulus; stimulus = struct;
stimulus.exp = exp;

if mglIsFile(exp.randColorsFile)
    stimulus.randcolors = load(exp.randColorsFile);
end

stimulus.fixcolor = [1;0;0];

%% specify task design

cps = {};
maxtrialtime        = 20; % seconds

experiment          = 'tau';
versionnum          = 3;

nblocks_learn       = 1;
ntrial_learn        = 4;  % learning phase at full luminance, not analyzed
nblocks             = 12;  % number of same blocks for each condition
trials_per_block    = 5;  % number of trials per block

% shuffle_set         = true;

if exp.debug, nblocks_learn=0; ntrial_learn= 1; nblocks=1; trials_per_block = 1; maxtrialtime=20; end
% if exp.debug
%     shuffle_set = false;
% end

experiment_paramset = [3]; 

% if shuffle_set
%     experiment_paramset = experiment_paramset(randperm(length(experiment_paramset)));
% end
Nconds = length(experiment_paramset);

for epset = experiment_paramset
    % learning phase 
    for b =1:nblocks_learn
        cps{end+1} = circular_ar(myscreen, 'numTrials', ntrial_learn, 'maxtrialtime', maxtrialtime, ...
            'dyn_noise_phase', length(cps)+1,  'switch_tpnoise', false, ...
            'experiment', {experiment}, 'experiment_paramset', epset);
    end

    % tracking
    for b = 1:nblocks
        cps{end+1} = circular_ar(myscreen, 'numTrials', trials_per_block, 'maxtrialtime', maxtrialtime, ...
            'dyn_noise_phase', length(cps)+1,  'switch_tpnoise', false, ...
            'experiment', {experiment}, 'experiment_paramset', epset);
    end
end

dynnoisefile = fullfile(find_root_dir, 'proj/grustim/trackpos/noise', ...
    sprintf('circular_ar_%s_%s.mat', experiment, num2str(versionnum)));

stimulus = check_and_load_dynnoise(dynnoisefile, myscreen, stimulus, cps);

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

mglBltTexture(mglText('A pointer object will appear on the screen and will be cued with an large ring.'),[0 2.5]);
mglBltTexture(mglText('You will be controlling the pointer.'),[0 1.5]);
% mglBltTexture(mglText('With the cue, two arrows will appear. This is the tracking axis. Turn the wheel clockwise to move towards the red arrow.'),[0 1.5]);
mglBltTexture(mglText('The target will appear briefly after.'),[0 0.5]);
mglBltTexture(mglText('1. Use the pointer object to track the target.'),[0 -1]);
if stimulus.exp.fixateCenter
    mglBltTexture(mglText('2. Fixate in center. You should be able to see the dot changing colors'),[0 -2]);    
end
mglBltTexture(mglText('When you are ready, press space to go to next trial'),[0 -4]);
mglFlush(); myscreen.flushMode = -1;

keystate = mglGetKeys;
while ~keystate(50)
    keystate = mglGetKeys;
end

if ~exp.showMouse, mglDisplayCursor(0);, end %hide cursor

phaseNum = 1;trialNum=0;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    [task{1}, myscreen, newphaseNum] = updateTask(task{1},myscreen,phaseNum); % update the task
    
    if (newphaseNum > length(task{1}))
        break
    end
    
    if  task{1}{newphaseNum}.trialnum ~= trialNum
        fprintf("Phase %s/%s, trial %s/%s \n",num2str(newphaseNum), num2str(length(task{1})),...
            num2str(task{1}{newphaseNum}.trialnum), num2str(task{1}{newphaseNum}.numTrials));        
        trialNum = task{1}{newphaseNum}.trialnum;
    end
    
    if newphaseNum ~= phaseNum        
        mglClearScreen(0.5);
        mglBltTexture(mglText('When you are ready, press space to go to next trial'),[0 -1.5]);
        mglFlush;
        
        mglClearScreen(0.5);
        mglBltTexture(mglText('When you are ready, press space to go to next trial'),[0 -1.5]);
        mglFlush;
        
        % wait for user to press space
        keystate = mglGetKeys;
        while ~keystate(50)
            keystate = mglGetKeys;
        end
        
        % check eyetracker
        if stimulus.exp.trackEye && stimulus.exp.trackEye_calibtrial
            [pos,postime] = mglEyelinkGetCurrentEyePos; % is this in image coordinates?
            if any(isnan(pos))
                disp('eye not found... recalibrating eye');
                beep on; beep;
                % load gong.mat; sound(y)
                mglClearScreen(0.5);
                myscreen  = eyeCalibDisp(myscreen); % calibrate eye every time.
            end
        end
        
        % give control back to the subject
        mglClearScreen(0.5);
        mglBltTexture(mglText('When you are ready, press space to go to next trial'),[0 -1.5]);
        mglFlush;
        
        mglClearScreen(0.5);
        mglBltTexture(mglText('When you are ready, press space to go to next trial'),[0 -1.5]);
        mglFlush;
        
        % wait for user to press space
        keystate = mglGetKeys;
        while ~keystate(50)
            keystate = mglGetKeys;
        end
        
        fprintf("Space detected. Beginning trial\n");
        
        % since we used up the first segment's time to do this, extend the
        % segment length. add 1 second to first segment to transition into the task.
        thisseg = 1;
        thistrial = task{1}{newphaseNum}.thistrial;
        task{1}{newphaseNum}.thistrial.seglen(thisseg) = mglGetSecs-thistrial.segstart + 1;
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
    save('savedir', 'frame')
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
    
    % these should be updated by the task function 
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
    
    if ~isfield(stimulus, 'targetfirst') || ~stimulus.targetfirst
        %% if we are in the tracking period,  save tracking variables
        if stimulus.task{phaseNum}.doTrack
            % update framecount
            task.thistrial.framecount = task.thistrial.framecount + 1;

            % display target
            % todo: another function for displaying objects
            if ~isempty(stimulus.target) 
                if isfield(stimulus.target, 'img') && ~isempty(stimulus.target.img)
                    mglBltTexture(stimulus.target.img, stimulus.target.position);
                else
                    if strcmp(stimulus.target.color,'*') && isfield(stimulus,'randcolors')
                        targcolors = stimulus.randcolors.colors(randi(size(stimulus.randcolors.colors,1)),:)';
                    elseif strcmp(stimulus.target.color,'*')
                        targcolors = [0.5+0.5*rand(3,1)];
                    else
                        targcolors = stimulus.target.color;
                    end

                    mglMetalDots([stimulus.target.position(1);stimulus.target.position(2);0], ...
                        [targcolors;1], [stimulus.target.std; stimulus.target.std], 1, 1);
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
                if strcmp(stimulus.pointer.color,'*') && isfield(stimulus,'randcolors')
                    pointercolors = stimulus.randcolors.colors(randi(size(stimulus.randcolors.colors,1)),:)';
                elseif strcmp(stimulus.pointer.color,'*')
                    pointercolors = [0.5+0.5*rand(3,1)];
                else
                    pointercolors = stimulus.pointer.color;
                end

                mglMetalDots([stimulus.pointer.position(1);stimulus.pointer.position(2);0], ...
                    [pointercolors; 1], [stimulus.pointer.std; stimulus.pointer.std], 1, 1);
            end
        end

    else % display pointer first
        % display pointer
        if ~isempty(stimulus.pointer)
            if isfield(stimulus.pointer, 'img') && ~isempty(stimulus.pointer.img)
                mglBltTexture(stimulus.pointer.img, stimulus.pointer.position);
            else 
                if strcmp(stimulus.pointer.color,'*') && isfield(stimulus,'randcolors')
                    pointercolors = stimulus.randcolors.colors(randi(size(stimulus.randcolors.colors,1)),:)';
                elseif strcmp(stimulus.pointer.color,'*')
                    pointercolors = [0.5+0.5*rand(3,1)];
                else
                    pointercolors = stimulus.pointer.color;
                end

                mglMetalDots([stimulus.pointer.position(1);stimulus.pointer.position(2);0], ...
                    [pointercolors; 1], [stimulus.pointer.std; stimulus.pointer.std], 1, 1);
            end
        end
            
        if stimulus.task{phaseNum}.doTrack
            % update framecount
            task.thistrial.framecount = task.thistrial.framecount + 1;

            % display target
            % todo: another function for displaying objects
            if ~isempty(stimulus.target) 
                if isfield(stimulus.target, 'img') && ~isempty(stimulus.target.img)
                    mglBltTexture(stimulus.target.img, stimulus.target.position);
                else
                    if strcmp(stimulus.target.color,'*') && isfield(stimulus,'randcolors')
                        targcolors = stimulus.randcolors.colors(randi(size(stimulus.randcolors.colors,1)),:)';
                    elseif strcmp(stimulus.target.color,'*')
                        targcolors = [0.5+0.5*rand(3,1)];
                    else
                        targcolors = stimulus.target.color;
                    end

                    mglMetalDots([stimulus.target.position(1);stimulus.target.position(2);0], ...
                        [targcolors;1], [stimulus.target.std; stimulus.target.std], 1, 1);
                end
            end

            % display other objects
            for ij = 1:length(stimulus.otherObjs)
                mglBltTexture(stimulus.otherObjs{ij}.img, stimulus.otherObjs{ij}.position);
            end
        end

    end


        % display fixation
        if stimulus.exp.fixateCenter == 1 && stimulus.task{phaseNum}.displayFix % fixation below others.
            mglMetalArcs([0;0;0], [1; 0; 0; 1], [stimulus.fixation_size+0.1;stimulus.fixation_size+0.2],[0;2*pi], 1);
            if isfield(stimulus, 'fixcolor') && ~(isstring(stimulus.fixcolor) && stimulus.fixcolor == "*")
                fixcolors = stimulus.fixcolor;
            elseif isfield(stimulus,'randcolors')
                fixcolors = stimulus.randcolors.colors(randi(size(stimulus.randcolors.colors,1)),:)';
            else
                fixcolors = [0.5+0.5*rand(3,1)];    
            end

            mglMetalDots([0;0;0], [fixcolors; 1], [stimulus.fixation_size;stimulus.fixation_size], 1, 1);
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


function stimulus = check_and_load_dynnoise(dynnoisefile, myscreen, stimulus, cps)
% check noise file
allgood = 1;
if mglIsFile(dynnoisefile)
    dyn_noise = load(dynnoisefile).dyn_noise;
    if length(cps) ~= length(dyn_noise)
        allgood=0;
    else
        for n = 1:length(cps)
            % check that the conditions match
            thistask    = cps{n};

            if dyn_noise(n).trialN < thistask.numTrials
                allgood = 0;
            end

            T = ceil(thistask.maxtrialtime * myscreen.framesPerSecond);
            if dyn_noise(n).T < T
                allgood = 0;
            end

            if dyn_noise(n).framesPerSecond ~= myscreen.framesPerSecond
                allgood = 0;
            end

            paramset = thistask.parameter_set(thistask.experiment, thistask.experiment_paramset);

            for pname = {'stim_noiseStd', 'stim_noiseTau', 'stim_vel', 'point_noiseStd', 'point_noiseTau', 'point_vel'}
               if dyn_noise(n).(pname{1}) ~= paramset.(pname{1})
                   allgood = 0;
               end
            end 
        end
    end
else
    allgood = 0;
end

% generate noise (but not save)
if ~allgood
   % regenerate noise 
   disp('(trackpos_circular_ar) dynamic noise file does not match. regenerating dynamics noise')
   dyn_noise = ar_gen_noise('myscreen', myscreen, 'taskcfg',cps,'savenoise',false);
end

stimulus.dyn_noise = dyn_noise;
    
end








































