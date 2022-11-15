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

% /Users/JRyu/github/mgl/task/displays, 0001_dn0a22167c_220912.mat
myscreen.displayName    = 'vpixx';
myscreen.calibType      = 'Specify particular calibration';
myscreen.calibFilename  = '/Users/gru/proj/mgl/task/displays/0001_dn0a221834_221005.mat';%'0001_dn0a221834_221005.mat'; % '0001_dn0a22167c_220912.mat';
myscreen.saveData       = 1; % save stimfile to data directory
myscreen.datadir        = '/Users/gru/data/';

myscreen = initScreen(myscreen);

% set to argb2101010 pixel format
mglMetalSetViewColorPixelFormat(4);

% Experimenter parameters
exp                     = struct();
exp.debug               = false;
exp.noeye               = true;
exp.eyemousedebug       = false;
exp.showmouse           = false;
exp.phasescrambleOn     = false;
exp.backprecompute      = false;
exp.feedback            = false; 
exp.afc.feedback_center = false;  % feedback about the exact center

exp.estim_horiz     = true;  % do hoiztonal estimation
exp.estim_verti     = false; % do vertical estimation
exp.colorfix        = false; % colored fixation

exp.afc.presSched    = 'staircase'; %'/Users/JRyu/Dropbox/GardnerLab/data/trackpos_2afc_staircase/s374/220929_stim05_staircase.mat';
exp.est.presSched    = 'gaussian';

exp.block_design    = false; % in each block, present all combinations of parameters
exp.noise_mask      = '/Users/gru/proj/grustim/trackpos/noise/grating.mat'; 

%% task parameters
% stimulus and background
task{1}{1}.random = 1; 

params            = struct();
params.task       = struct();
params.task.backLum    = 0.7;%0.4; %32;  % background luminance; units: fraction of full luminance 
params.task.noiseLum   = 0; % noise luminance, if there is one.

% main task parameters
tasks2run                   = {'2afc'}; %{'est', '2afc'};
params.task.stimLum         = [0.1, 0.2, 0.4, 0.8]; %, 0.1, 0.2, 0.4]; %[0.1, 1]; %[0.05, 0.1, 0.2, 0.4]; % [0.1,0.2,0.5] % [16,32,48,96]
params.task.stimDur         = [2/60, 4/60, 6/60, 10/60, 15/60, 30/60]; %[2/60, 3/60, 4/60, 6/60, 10/60, 15/60]; %[2/60, 4/60, 6/60, 10/60, 15/60, 30/60]; 
params.task.stimStd         = [1]; % [1,1.5]
params.task.stimColor       = 'k';
params.trialpercond         = 40;

% mask parameters
if mglIsFile(exp.noise_mask)
    params.task.maskDur          = 3/60; %[0]; %4/60, 8/60]; % mask duration
    params.task.mask_TOff2MOn    = 0; % 0, 4/60, 8/60]; %, 2/60, 5/60]; % 3/60, 5/60]; % stimulus offset to mask onset (Neisser 1967)
    params.task.maskLum          = [0.6]; %0.7]; %[0.05, 0.7];
end

if exp.debug
    params.task.stimLum         = [0.05, 0.4]; % [0.1,0.2,0.5] % [16,32,48,96]
    params.task.stimDur         = [2/60]; %[2/60 5/60 10/60 15/60]; %frames/hz
    params.task.mask_TOff2MOn   = [1/60]; % stimulus offset to mask onset (Neisser 1967)
    params.task.maskLum         = [0.2];
    params.trialpercond         = 2; 
end

% afc parameters
params_afc = params;
params_afc.task.pointerOffset = [0]; % [-10,-5,-2,0,2,5,10];
if exp.afc.presSched == 'staircase'
    params_afc.presSched    = 'staircase';
    params_afc.staircase                    = struct();
    params_afc.staircase.initThreshold      = 0.3;
    params_afc.staircase.initThresholdSd    = 0.3;
end

if exp.debug
    params_afc.task.pointerOffset   = [0, 5];
end

% est parameters
params_est = params;
if exp.est.presSched    == 'gaussian'
    params_est.presSched    = 'gaussian';
    params_est.prior        = struct();
    params_est.prior.std    = 2;
    params_est.prior.mean   = 0;
elseif exp.est.presSched == 'staircase'
    params_est.staircase                    = struct();
    params_est.staircase.initThreshold      = 0.3;
    params_est.staircase.initThresholdSd    = 0.3;
end

% count conditions - 2AFC
params.numTrials            = 0; % try to match number of trials

if any(cellfun(@(x) strcmp(x,'2afc'),tasks2run))
    [afc_fields, afc_vals] = countconditions(params_afc.task);
    afc_comb = allcomb(afc_vals{:});
    disp(['Number of conditions (afc) = ' num2str(size(afc_comb,1))])
    
    params_afc.numTrials        = params_afc.trialpercond * size(afc_comb,1);
    params.numTrials            = params.numTrials + params_afc.numTrials; % main task
else
    params_afc.numTrials        = 0;
end

if any(cellfun(@(x) strcmp(x,'est'),tasks2run))
    [est_fields, est_vals] = countconditions(params_est.task);
    est_comb = allcomb(est_vals{:});
    disp(['Number of conditions (est) = ' num2str(size(est_comb,1))])

    params_est.numTrials        = params_afc.trialpercond * size(est_comb,1);
    params.numTrials            = params.numTrials + params_est.numTrials; % main task
else
    params_est.numTrials = 0;
end


task{1}{1}.parameter.currtask   = tasks2run; % for fixed values
task{1}{1}.segmin           = [inf]; % for running other tasks
task{1}{1}.segmax           = [inf]; % jumpsegment if the other task is finished
% task{1}{1}.synchToVol     = [1]; % wait for backtick before going onto next trial

task{1}{1}.waitForBacktick  = 1;

%tasks * ntrials x stimDur x stimLum x posDiff
task{1}{1}.numTrials        = params.numTrials; 

trialdur = 1.5 + max(params.task.stimDur) + max(params.task.mask_TOff2MOn) + params.task.maskDur + 1;
taskdur = trialdur * params.numTrials / 60/60; % approximate duration in hours
disp(['Approx task duration = ' num2str(taskdur) ' hours']);

%% initialize stimulus
global stimulus;
stimulus = [];

stimulus.exp            = exp;
stimulus.params         = params; % not saved in the task.
stimulus.params_afc     = params_afc;
stimulus.params_est     = params_est;
stimulus.target         = trackposInitStimulus(stimulus,myscreen);

stimulus.fixColors.stim     = [1 0 0]; % red
stimulus.fixColors.est      = [0 1 0]; % fixation color at response
stimulus.fixColors.afc      = [0 0 1]; % afc response period 
stimulus.fixColors.fb       = [1 1 1]; % position feedback

stimulus.t0 = mglGetSecs; % keeps track of trackTime

myscreen = initStimulus('stimulus',myscreen); % what does this do???

if stimulus.exp.phasescrambleOn == 1
    disp('Loading phase scrambled background noise...')

    tic
    if stimulus.exp.backprecompute == 1
        savefile = '/Users/gru/proj/grustim/trackpos/trackpos.mat';
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

% noise mask
if mglIsFile(stimulus.exp.noise_mask)
    stimulus.noise_mask = load(stimulus.exp.noise_mask);
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

% 2AFC subtask
[task{2}, myscreen] = trackpos_sub_2afc(myscreen,params_afc,exp); 
% estimation subtask
[task{3}, myscreen] = trackpos_sub_est(myscreen,params_est,exp); 

%% run the task
stimulus.t0 = mglGetSecs; % 

% explain task.
mglDisplayCursor(0); %hide cursor
mglClearScreen(params.task.backLum);
mglTextDraw('task (trackpos_multitask) starting... ', [0 3])
% mglTextDraw('After the stimulus is presented, you will be asked to perform one of the two tasks, depending on the fixation color',[0 1]);
if any(cellfun(@(x) strcmp(x,'est'),tasks2run))
    mglTextDraw('Estimation task (red fixation): move the mouse to the center of stimulus. Press 3 when done.',[0 1]);
end
if any(cellfun(@(x) strcmp(x,'2afc'),tasks2run))
    mglTextDraw('2AFC task (white fixation): press 1 if the stimulus is to the left of red reference. 2 otherwise',[0 -1]);
end
mglBltTexture(mglText('When you are ready, press backtick to go to the first trial'),[0 -4]);
mglFlush(); myscreen.flushMode = -1;

if ~exp.showmouse, mglDisplayCursor(0);, end %hide cursor

phaseNum{1} = 1; phaseNum{2}=1; phaseNum{3}=1;
while (phaseNum{1} <= length(task{1})) && ~myscreen.userHitEsc && ...
        (phaseNum{2} <= length(task{2}) || phaseNum{3} <= length(task{3}))
    [task{1}, myscreen, phaseNum{1}]   = updateTask(task{1},myscreen,phaseNum{1});     % update the main task
    
    if phaseNum{2} <= length(task{2}) % run 2afc
        [task{2}, myscreen, phaseNum{2}]  = updateTask(task{2},myscreen,phaseNum{2});
    end
    if phaseNum{3} <= length(task{3}) % run estimation
        [task{3}, myscreen, phaseNum{3}]  = updateTask(task{3},myscreen,phaseNum{3});
    end
    myscreen                        = tickScreen(myscreen,task);     % flip screen
end

myscreen = endTask(myscreen,task);
mglClose; endScreen(myscreen); mglDisplayCursor(1) %show cursor

if stimulus.exp.afc.presSched    == 'staircase'
    staircase = stimulus.staircaseTable;
    save([myscreen.stimfile(1:end-4),'_staircase.mat'], 'staircase');
end


end

%% Initialize trials;
function [task, myscreen] = initTrialCallback(task, myscreen)
    % nan out the parameters so that we don't analyze them (does this work?)
    task.thistrial.posDiff      = nan; % forst fixed values
    task.thistrial.stimLum      = nan;
    task.thistrial.stimDur      = nan;
    task.thistrial.stimStd      = nan; 
    
    % print trial number every 5%. 
    if mod(task.trialnum,ceil(task.numTrials/20)) == 1
        disp(['(trackpos_multitask) '  num2str(task.trialnum/task.numTrials) ...
            '% finished: Trial ' num2str(task.trialnum) ' / ' num2str(task.numTrials)]);
    end
end

%% Start segment
function [task, myscreen] = startSegmentCallback(task, myscreen)
    global stimulus
    
    % select which task to run and save it in the stimulus
    stimulus.currtask = task.thistrial.currtask;
    myscreen.flushMode = 0; % start updating screen for the subtasks to detect in screenupdate.
end

%% screen update
function [task, myscreen] = screenUpdateCallback(task, myscreen)
    global stimulus
    % code up a break => set flushMode to -1
    if strcmp(stimulus.currtask,'done')
        stimulus.currtask = 'initializing new task';
        task = jumpSegment(task);
    end
end

function [task, myscreen] = responseCallback(task, myscreen)
    global stimulus
    if task.thistrial.whichButton == 0
        % go to next segment
        task = jumpSegment(task);
    end
end