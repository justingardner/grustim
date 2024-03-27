%        $Id: $
%      usage: trackpos_2afc_masked_sc
%         by: Josh Ryu
%       date: 04/13/2021
%    purpose: 2 choice task for position of one blob against fixation.

% S1: wait time while the main task calls this subtask.
% S2: cue (1s)
% S3: random period of fixation (random ~0.5s)
% S4: stimulus period (stimdur s)
% S5: delay to mask
% S6: mask
% S7: repsonse period (inf)
% S8: feedback (1s)

function myscreen = trackpos_2afc_masked_sc(varargin)

%getArgs(varargin,{'subjectID=s999','centerX=10','centerY=0','diameter=16'}); getArgs(varargin,{'subjectID=-1'});

myscreen        = setup_screen_jryu(); 
myscreen        = initScreen(myscreen);
mglMetalSetViewColorPixelFormat(4);     % set to argb2101010 pixel format
rootdir         = find_root_dir();
rng(0, 'twister'); % set seed

%% Experimenter parameters
exp                     = struct();
exp.debug               = false;
exp.trackEye            = true;
exp.enforceFixThresh    = Inf;
exp.showmouse           = false;

exp.feedback            = false;  % correct or incorrect in 2afc
exp.feedback_center     = false;  % feedback about the exact center

exp.colorfix            = true; % colored fixation
exp.colorref            = true; % colored fixation
exp.displacement_type   = 'circular'; %'tangential'; % otherwise specify polar angle and displacement angle
exp.respDirArrow        = true;

exp.block_design        = false; % in each block, present all combinations of parameters
exp.noise_mask          = fullfile(rootdir, 'proj/grustim/trackpos/noise/grating.mat'); 
exp.staircase_init      = fullfile(rootdir,'data/trackpos_2afc_masked_sc/', mglGetSID, ...
                                   '230330_stim03_staircase.mat');


%% task parameters
% stimulus and background
task{1}{1}.random = 1; 

params                  = struct();

params.presSched        = 'staircase';
params.trialpercond     = 100;

params.task             = struct();
params.task.backLum     = 0.7;%0.4; %32;  % background luminance; units: fraction of full luminance 
params.task.noiseLum    = 0; % noise luminance, if there is one.

% main task parameters
params.task.stimLum         = 0.4; %[0.2, 0.8]; %, 0.1, 0.2, 0.4]; %[0.1, 1]; %[0.05, 0.1, 0.2, 0.4]; % [0.1,0.2,0.5] % [16,32,48,96]
params.task.stimDur         = [2/60, 4/60, 10/60, 20/60]; %, 30/60]; %[2/60, 3/60, 4/60, 6/60, 10/60, 15/60]; %[2/60, 4/60, 6/60, 10/60, 15/60, 30/60]; 
params.task.stimStd         = [1]; % [1,1.5]
params.task.stimColor       = 'k';

if strcmp(exp.displacement_type, 'tangential')
    params.task.angleSet    = 1:8; % polar angles
elseif strcmp(exp.displacement_type, 'circular')
    params.task.angleSet        = -1; %[1,2,3]; % polar angles
    params.task.displ_type      = {'circular'};
else
    params.task.polarAngle = 0;
    params.task.displAngle = pi/2;
end

% mask parameters
if mglIsFile(exp.noise_mask)
    params.task.maskDur          = 15/60; %[0]; %4/60, 8/60]; % mask duration
    params.task.mask_TOff2MOn    = 0; % 0, 4/60, 8/60]; %, 2/60, 5/60]; % 3/60, 5/60]; % stimulus offset to mask onset (Neisser 1967)
    params.task.maskLum          = [0.6]; %0.7]; %[0.05, 0.7];
end

params.task.pointerOffset       = [10, 15, 20]; %[3, 7, 10]; %3,7,10 % [-10,-5,-2,0,2,5,10];

% staircase parameters
params.staircase                    = struct();
% thresh = params.task.pointerOffset(1)*1.7/10 + 0.1;
params.staircase.initThreshold      = 0.3; %0.3;
params.staircase.initThresholdSd    = 0.3; %0.3;
params.staircase.threshstd_thresh   = 0.1; % 0.01;
params.staircase.staircase_init     = exp.staircase_init;

if exp.debug
    params = load_debug_params(params);
end

[afc_fields, afc_vals] = countconditions(params.task);
if isempty(afc_vals)
    afc_comb = {1};
else
    afc_comb = allcomb(afc_vals{:});
end
disp(['Number of conditions (afc) = ' num2str(size(afc_comb,1))])
params.numTrials        = params.trialpercond * size(afc_comb,1);

task{1}{1}.segmin           = [inf]; % for running other tasks
task{1}{1}.segmax           = [inf]; % jumpsegment if the other task is finished

task{1}{1}.waitForBacktick  = 1;

%tasks * ntrials x stimDur x stimLum x posDiff
task{1}{1}.numTrials        = params.numTrials; 

%% initialize stimulus
global stimulus;
stimulus = [];

stimulus.exp            = exp;
stimulus.params         = params; % not saved in the task.
stimulus.target         = trackposInitStimulus(stimulus,myscreen);

stimulus.fixColors.stim     = [1 0 0]; % red
stimulus.fixColors.est      = [0 1 0]; % fixation color at response
stimulus.fixColors.afc      = [0 0 1]; % afc response period 
stimulus.fixColors.fb       = [1 1 1]; % position feedback

stimulus.pointerR           = 0.4;

stimulus.reference  = struct();
if exp.colorref
    stimulus.reference.color    = '*';
else
    stimulus.reference.color    = 'r';
end

stimulus.t0 = mglGetSecs; % keeps track of trackTime

myscreen = initStimulus('stimulus', myscreen); % what does this do???

% phase Scrambled background
if isfield(stimulus.exp, 'phasescrambleOn') && (stimulus.exp.phasescrambleOn == 1)
    disp('Loading phase scrambled background noise...')

    tic
    if stimulus.exp.backprecompute == 1
        savefile = fullfile(rootdir, 'proj/grustim/trackpos/trackpos.mat');
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

% random colors
exp.randColorsFile      = '/Users/gru/proj/grustim/trackpos/util/labcolors.mat'; %'/Users/jryu/proj/grustim/trackpos/util/labcolors.mat'; % 
if mglIsFile(exp.randColorsFile)
    stimulus.randcolors = load(exp.randColorsFile);
else
    disp('random color data does not exist')
end


%% Eye calibration
if stimulus.exp.trackEye && ~ stimulus.exp.debug
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
[task{2}, myscreen] = trackpos_sub_2afc(myscreen,params,exp); 

%% run the task
stimulus.t0 = mglGetSecs; % 

% explain task.
mglDisplayCursor(0); %hide cursor
mglClearScreen(params.task.backLum);
mglTextDraw('task (trackpos_2afc_masked_sc) starting... ', [0 3])
mglTextDraw('Press 1 if the blob was in the direction of the blue arrow from the red reference dot. 2 otherwise',[0 1]);

mglBltTexture(mglText('When you are ready, press backtick to go to the first trial'),[0 -3]);
mglFlush(); myscreen.flushMode = -1;

if ~exp.showmouse, mglDisplayCursor(0);, end %hide cursor

phaseNum{1} = 1; phaseNum{2}=1; phaseNum{3}=1;
while (phaseNum{1} <= length(task{1})) && ~myscreen.userHitEsc && ...
        (phaseNum{2} <= length(task{2}))
    [task{1}, myscreen, phaseNum{1}]        = updateTask(task{1},myscreen,phaseNum{1});     % update the main task
    
    if phaseNum{2} <= length(task{2}) % run 2afc
        [task{2}, myscreen, phaseNum{2}]    = updateTask(task{2},myscreen,phaseNum{2});
    end
    
    myscreen                        = tickScreen(myscreen,task);     % flip screen
end

% save temporary staircase first
if isfield(task{2}{1}, 'private') && isfield(task{2}{1}.private,'staircaseTable')
    staircase = task{2}{1}.private.staircaseTable;
    save(fullfile(myscreen.datadir,'temp_staircase.mat'), 'staircase');
end

% end task
myscreen = endTask(myscreen,task);
mglClose; endScreen(myscreen); mglDisplayCursor(1) %show cursor

% save staircase to the correct stimfile
if isfield(task{2}{1}, 'private') && isfield(task{2}{1}.private,'staircaseTable')
    staircase = task{2}{1}.private.staircaseTable;
    for idx = 1:size(staircase,1)
        staircase.trialnum(idx) = staircase.staircase{idx,1}.trialNum;
        staircase.quest_thresh_std(idx) = QuestSd(staircase.staircase{idx,1}.s);
    end
    save([myscreen.stimfile(1:end-4),'_staircase.mat'], 'staircase');
end

end

%% Initialize trials;
function [task, myscreen] = initTrialCallback(task, myscreen)
    % nan out the parameters so that we don't analyze them (does this work?)
    task.thistrial.posDiff      = nan; % for fixed values
    task.thistrial.stimLum      = nan;
    task.thistrial.stimDur      = nan;
    task.thistrial.stimStd      = nan; 
    
    % print trial number every 5%. 
    if mod(task.trialnum,ceil(task.numTrials/20)) == 1
        disp(['(trackpos_multitask) '  num2str(task.trialnum/task.numTrials) ...
            '% finished: Trial ' num2str(task.trialnum) ' / ' num2str(task.numTrials)]);
    end

    global stimulus
    if stimulus.exp.trackEye
        [pos,postime] = mglEyelinkGetCurrentEyePos; 
        while norm(pos) > stimulus.exp.enforceFixThresh  %|| any(isnan(pos))
            mglMetalDots([0;0;0], [0.5+0.5*rand(3,1);1], [stimulus.pointerR; stimulus.pointerR], 1, 1);
            mglTextDraw('Please fixate on the middle of the screen !!', [0 2]);
            mglTextDraw('Please fixate on the middle of the screen !!', [0 -2]);
            mglFlush;
            [pos,postime] = mglEyelinkGetCurrentEyePos; 
        end
    end
end

%% Start segment
function [task, myscreen] = startSegmentCallback(task, myscreen)
    global stimulus
    
    % select which task to run and save it in the stimulus
    stimulus.currtask = '2afc';
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
    global stimulus;
    if task.thistrial.whichButton == 0
        % go to next segment
        task = jumpSegment(task);
    end
end



function params = load_debug_params(params)
    params.trialpercond         = 3; 
end