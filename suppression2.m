% Suppression 2
% V1 Aug 5 2015: Setup basic experiment parameters
% V2 July 8 2016: memory task: location cued centrally; if probe appears at
% the same location, make a saccade to it; if it appears elsewhere - keep
% fixating


function myscreen = suppression2

% check arguments
if ~any(nargin == [0])
    help taskTemplate
    return
end

% Initalize the screenx
myscreen = initScreen;
% myscreen = initScreen('fMRIproj32');
myscreen.background = myscreen.gray;

% MRI
task{1}.waitForBacktick = 0;

% Durations
task{1}.segmin = [1 0.5 10 0.5 3];
task{1}.segmax = [1 0.5 10 0.5 11];
task{1}.getResponse = [0 0 0 0 1];
task{1}.synchToVol = [0 0 0 0 0];

task{1}.random = 1;
task{1}.seglenPrecompute = 1;

% Parameteres
numrep1 = 6; % How many time to repeat trials per MRI session
t0=repmat([45,135,225,315],[1,numrep1]);
t1=randperm(length(t0)); t2=t0(t1);
task{1}.parameter.memory_cue_location = t2;

% Random settings
task{1}.randVars.uniform.t2_location = [-180];  % Relative to memory
task{1}.randVars.uniform.response_target = [0 1];
task{1}.randVars.uniform.peripheral_trial = [0 0];

% Exp duration
task{1}.numBlocks=1;
% task{1}.numTrials=20;


%% Stimulus parameters

% initialize the task
for phaseNum = 1:length(task)
    [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end

% init the stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);

% fix: you will change the funciton myInitStimulus
% to initialize the stimulus for your experiment.
stimulus = myInitStimulus(stimulus,myscreen);

% Fixation coordinates
stimulus.fixation_position = [0,0];
stimulus.fixation_size = [0.3, 0.3];
stimulus.fixation_color = [0.2, 0.2, 0.2];
stimulus.fixation_isi_color = [0.2, 0.2, 0.8];

% Memory target (line)
stimulus.memory_cue_color = [0.2, 0.2, 0.2];
stimulus.memory_cue_width = 5;
stimulus.memory_cue_length = 0.4;

% Memory target (peripheral)
stimulus.memory_cue_peripheral_size = [1, 1];
stimulus.memory_cue_peripheral_color = [0.8, 0.2, 0.2];

% Frames
stimulus.frames_positions = [45:90:359];
stimulus.frames_size = [1.5 1.5];
stimulus.frames_color = [0.45, 0.45, 0.45];

% Response targets
stimulus.target_radius = 8;
stimulus.t1_size = [1.5, 1.5];
stimulus.t1_color = [0.2, 0.2, 0.2];
stimulus.t2_size = [1.5, 1.5];
stimulus.t2_color = [0.2, 0.2, 0.2];

% mglClose;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phaseNum = 1;
while (phaseNum <= length(task)) && ~myscreen.userHitEsc
    % update the task
    [task myscreen phaseNum] = updateTask(task,myscreen,phaseNum);
    % flip screen
    myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;

% fix: do anything that needs to be done at the beginning
% of a segment (like for example setting the stimulus correctly
% according to the parameters etc).
if task.thistrial.thisseg==1
    if task.thistrial.response_target == 1
        disp('fixation trial');
    else
        disp('saccade trial');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus

% fix: display your stimulus here, for this code we just display
% a fixation cross that changes color depending on the segment
% we are on.
fixation_x = stimulus.fixation_position(1);
fixation_y = stimulus.fixation_position(2);

% Calculate coordinates of placeholders
posarc = stimulus.frames_positions;
radius = stimulus.target_radius;
[xc, yc] = pol2cart(posarc*pi/180, radius);
frames_x = xc; frames_y = yc;

% Calculate coordinates of the memory target cue (line)
posarc = task.thistrial.memory_cue_location;
radius = stimulus.memory_cue_length;
[xc, yc] = pol2cart(posarc*pi/180, radius);
memory_cue_x = xc;
memory_cue_y = yc;

% Calculate coordinates of the memory target (peripheral)
posarc = task.thistrial.memory_cue_location;
radius = stimulus.target_radius;
[xc, yc] = pol2cart(posarc*pi/180, radius);
memory_cue_peripheral_x = xc;
memory_cue_peripheral_y = yc;

% Calculate coordinates of T1
posarc = task.thistrial.memory_cue_location;
radius = stimulus.target_radius;
[xc, yc] = pol2cart(posarc*pi/180, radius);
target1_x = xc;
target1_y = yc;

% Calculate coordinates of T2
posarc = task.thistrial.memory_cue_location + task.thistrial.t2_location;
radius = stimulus.target_radius;
[xc, yc] = pol2cart(posarc*pi/180, radius);
target2_x = xc;
target2_y = yc;

mglClearScreen(0.5);
% mglBltTexture(stximulus.texture,[0 0,24,24]);

% Trial sequence
if (task.thistrial.thisseg == 1) % Fixation
    mglFillOval(fixation_x, fixation_y, stimulus.fixation_size,  stimulus.fixation_color);
    mglFillRect(frames_x, frames_y, stimulus.frames_size,  stimulus.frames_color);
elseif (task.thistrial.thisseg == 2) % Memory cue
    mglFillOval(fixation_x, fixation_y, stimulus.fixation_size,  stimulus.fixation_color);
    mglFillRect(frames_x, frames_y, stimulus.frames_size,  stimulus.frames_color);
    if task.thistrial.peripheral_trial==0
        mglLines2(fixation_x, fixation_y, memory_cue_x, memory_cue_y, stimulus.memory_cue_width, stimulus.memory_cue_color);
    elseif task.thistrial.peripheral_trial==1
        mglFillRect(memory_cue_peripheral_x, memory_cue_peripheral_y, stimulus.memory_cue_peripheral_size, stimulus.memory_cue_peripheral_color);
    end
elseif (task.thistrial.thisseg == 3) % Delay
    mglFillOval(fixation_x, fixation_y, stimulus.fixation_size,  stimulus.fixation_color);
    mglFillRect(frames_x, frames_y, stimulus.frames_size,  stimulus.frames_color);
elseif (task.thistrial.thisseg == 4) % Response targets
    mglFillOval(fixation_x, fixation_y, stimulus.fixation_size,  stimulus.fixation_color);
    mglFillRect(frames_x, frames_y, stimulus.frames_size,  stimulus.frames_color);
    if task.thistrial.response_target==0
        mglFillOval(target1_x, target1_y, stimulus.t1_size,  stimulus.t1_color);
    elseif task.thistrial.response_target==1
        mglFillOval(target2_x, target2_y, stimulus.t2_size,  stimulus.t2_color);
    end
elseif task.thistrial.thisseg == 5 % Inter-trial interval
    mglFillOval(fixation_x, fixation_y, stimulus.fixation_size,  stimulus.fixation_color);
    mglFillRect(frames_x, frames_y, stimulus.frames_size,  stimulus.frames_color);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)

global stimulus

% mglDeleteTexture(task.texture)

% fix: add the code you want to use to process subject response

% here, we just check whether this is the first time we got a response
% this trial and display what the subject's response was and the reaction time
if task.thistrial.gotResponse < 1
    disp(sprintf('Subject response: %i Reaction time: %0.2fs',task.thistrial.whichButton,task.thistrial.reactionTime));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen)

stimulus.texture = mglCreateTexture(round(rand(500,500)*255));



