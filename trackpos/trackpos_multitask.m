%
%        $Id: $
%      usage: trackpos_xtasks
%         by: Josh Ryu
%       date: 05/19/2021
%    purpose: estimate position of blob. Alternates between tasks.


% S1: random period of fixation (random ~0.5s)
% S2: stimulus period (stimdur s)

% task 1: estimation
% S3: repsonse period. move mouse (until response)
% S4: feedback (1s)

% task 2: 2afc
% S3: repsonse period (until response)
% S4: feedback (1s) position or corr/incorr?

% fixed position difference or samples from continuous distributions?

function myscreen = trackpos_multitask(varargin)

%% set up screen and experiment
% set input arguments
if isempty(mglGetSID)
    myscreen.subjectID  = -1;
else
    myscreen.subjectID  = mglGetSID;
    myscreen.saveData = 1;
end

myscreen.saveData       = 1; % save stimfile to data directory
myscreen = initScreen(myscreen);

% Experimenter parameters
exp                 = struct();
exp.noeye           = true;
exp.eyemousedebug   = false;
exp.showmouse       = false;
exp.debug           = true;
exp.phasescrambleOn = true;
exp.backprecompute  = true;
exp.feedback        = true; 
exp.estim_horiz     = true;  % do hoiztonal estimation
exp.estim_verti     = false; % do vertical estimation
exp.colorfix        = false;

%% task parameters
% stimulus and background
task{1}{1}.random = 1; 

params.backLum    = 90; %32;  % background luminance; units: fraction of full luminance 
params.noiseLum   = 32; % noise luminance, if there is one.

% main task parameters
tasks2run         = {'est', '2c'};
teststimLum       = [1, 1.5] * params.noiseLum; %SNR
teststimDur       = [2/60, 5/60, 10/60]; %[2/60 5/60 10/60 15/60]; %frames/hz
posDiff           = logspace(log(0.05)/log(10),log(0.5)/log(10),8); % in degs; minimum and maximum offset from fixation
trialpercond      = 12;
if exp.debug, trialpercond = 1; end

task{1}{1}.parameter.currtask   = tasks2run; % forst fixed values
params.posDiff   = posDiff; % forst fixed values
params.stimLum 	= teststimLum;
params.stimDur 	= teststimDur; % teststimDur is also saved under stimulus
params.numTrials = length(tasks2run) * trialpercond*length(teststimDur) * ...
    length(teststimLum)*2*length(posDiff);

task{1}{1}.segmin           = [inf]; % for running other tasks
task{1}{1}.segmax           = [inf]; % jumpsegment if the other task is finished
% task{1}{1}.synchToVol       = [1]; % wait for backtick before going onto next trial

task{1}{1}.waitForBacktick  = 1;

%tasks * ntrials x stimDur x stimLum x posDiff
task{1}{1}.numTrials        = params.numTrials; 
taskdur = (0.5 + max(teststimDur) + 1 + 1) * ...
    task{1}{1}.numTrials/60/60; % approximate duration in hours
disp(['Approx task duration = ' num2str(taskdur) ' hours']);

%% initialize stimulus
global stimulus;
stimulus = [];

stimulus.exp = exp;
stimulus.teststimDur        = teststimDur; % not saved in the task.

stimulus = trackposInitStimulus(stimulus,myscreen);
myscreen = initStimulus('stimulus',myscreen); % what does this do???

if stimulus.exp.phasescrambleOn == 1
    disp('Loading phase scrambled background noise...')

    tic
    if stimulus.exp.backprecompute == 1
        savefile = '/Users/gru/data/trackpos/trackpos.mat';
        %savefile = '/Users/jryu/data/trackpos/trackpos.mat'; 
        % savefile = '/Users/joshua/data/trackpos_2afc/trackpos.mat'; % just use noise 1 and permute
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
    end
    toc
end

%% Eye calibration
if ~stimulus.exp.noeye && ~ stimulus.exp.debug
    disp(' Calibrating Eye ....')
    myscreen = eyeCalibDisp(myscreen);
    
    % let the experimenter know
    disp(sprintf('(trackpos) Starting Run...'));
end


%% task blocks. 
% initializing task...
disp(' Initializing Task....')

for phaseN = 1:length(task{1})
    [task{1}{phaseN}, myscreen] = initTask(task{1}{phaseN},myscreen,...
        @startSegmentCallback,@screenUpdateCallback,@responseCallback,@initTrialCallback);
end

% estimation subtask
[task{2}, myscreen] = trackpos_sub_est(myscreen,params,exp); 
% 2AFC subtask
[task{3}, myscreen] = trackpos_sub_2c(myscreen,params,exp); 

%% run the task
stimulus.t0 = mglGetSecs; % 

% explain task.
mglDisplayCursor(0); %hide cursor
% mglClearScreen(task{1}{1}.parameter.backLum/255);
mglTextDraw('task (trackpos_multitask) starting... ', [0 1])
mglTextDraw('After the stimulus is presented, you will be asked to perform one of the two tasks, depending on the fixation color',[0 0]);
mglTextDraw('Estimation task (red fixation): move the mouse to the center of stimulus. Press 3 when done.',[0 -1],...
    'fontColor', [1 0 0]);
mglTextDraw('2AFC task (blue fixation): press 1 if the stimulus is to the right of fixation. 2 otherwise',[0 -2],...
    'fontColor', [0 0 1]);
mglTextDraw('press backtick to go to next trial',[0 -3]);
mglFlush

if ~exp.showmouse, mglDisplayCursor(0);, end %hide cursor

phaseNum = 1;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);     % update the task
    [task{2}, myscreen] = updateTask(task{2},myscreen,1);
    [task{3}, myscreen] = updateTask(task{3},myscreen,1);
    myscreen = tickScreen(myscreen,task);     % flip screen
end

myscreen = endTask(myscreen,task);
mglClose; endScreen(myscreen); mglDisplayCursor(1) %show cursor

end

%% Initialize trials;
function [task, myscreen] = initTrialCallback(task, myscreen)
    % nan out the parameters so that we don't analyze them (does this work?)
    task.thistrial.posDiff      = nan; % forst fixed values
    task.thistrial.stimLum      = nan;
    task.thistrial.stimDur      = nan; % teststimDur is also saved under stimulusnan
    
    % print trial number every 5%. 
    if mod(task.trialnum,ceil(task.numTrials/20)) == 1
        disp(['Trial ' num2str(task.trialnum) ' / ' num2str(task.numTrials)]);
    end
end

%% Start segment
function [task, myscreen] = startSegmentCallback(task, myscreen)
    global stimulus
    
    % select which task to run and save it in the stimulus
    stimulus.currtask = task.thistrial.currtask;
end

%% screen update
function [task, myscreen] = screenUpdateCallback(task, myscreen)
    global stimulus
    if strcmp(stimulus.currtask,'done')
        stimulus.currtask = 'initializing new task';
        task = jumpSegment(task);
    end
end