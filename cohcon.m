% cohCon
%
%      usage: myscreen=cohcon()
%         by: daniel birman
%       date: 11/10/14
%    purpose: contrast change detection with cued selective attention.
%
%        use: call flowAwe() to initialize. The first
%             run takes significantly longer due to loading stimuli.
%
%      flags: stimFileNum (-1/#) - Load a specific stimfile from a
%             subject's folder. Defaults to the last stimfile. Warning:
%             Only the first stimfile is saved with the image file data,
%             subsequent stimfiles only contain run data.
%             unattended (0/1) - If 1, runs a fixation task while the main
%             task just runs on idle in the background (no inputs)
%             plots (0/1) - Displays staircase plots (and estimated
%             psychophysic functions)
%             overrideTask (1/2) - Specifies the task to run: 1 =
%             coherence, 2 = contrast
%             projector (0/1) - Masks stimuli using the default projector
%             mask.
%             scan (0/1) - Scanner timing
%
%
%   TR .75 = 560 volumes (7:00 total)
%   TR 1.4 = 300 volumes (7:00 total)

function [myscreen] = cohcon(varargin)

global stimulus
clear fixStimulus
%% Initialize Variables

% add arguments later
stimFileNum = [];
plots = [];
overrideTask = 0;
scan = [];
nocatch = [];
stablecon = 0;
stablecoh = 0;
constant = [];
getArgs(varargin,{'stimFileNum=-1','nocatch=0',...
    'plots=0','overrideTask=0','scan=0','constant=1','stablecon=0','stablecoh=0'});
stimulus.scan = scan;
stimulus.plots = plots;
stimulus.nocatch = nocatch;
stimulus.stablecon = stablecon;
stimulus.stablecoh = stablecoh;
stimulus.constant = constant; % new param, keeps stimulus on screen at all times with 0% coherence

if stimulus.scan && ~stimulus.nocatch
    warning('You didn''t force nocatch during a scan, Auto-setting');
    stimulus.nocatch = 1;
end
if stimulus.stablecon && ~overrideTask==1
    warning('You didn''t specify task = motion, but you stabilized contrast. Auto-setting');
    overrideTask = 1;
end
if stimulus.stablecoh && ~overrideTask==2
    warning('You didn''t specify task = contrast, but you stabilized coherence. Auto-setting');
    overrideTask = 2;
end

if stimulus.scan && ~mglGetParam('ignoreInitialVols')==16 && ~mglGetParam('ignoreInitialVols')==4
    warning('ignoreInitialVols is set to %i.',mglGetParam('ignoreInitialVols'));
    if ~strcmp('y',input('Is this correct? [y/n]'))
        mglSetParam('ignoreInitialVols',input('Please input the correct value (mux8 = 16, mux2 = 4): '));
    end
end

stimulus.counter = 1; % This keeps track of what "run" we are on.

%% Setup Screen

if stimulus.scan
    myscreen = initScreen('fMRIprojFlex');
else
    myscreen = initScreen('VPixx');
end

myscreen.background = 0.5;

%% Open Old Stimfile
stimulus.initStair = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/cohcon/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/cohcon/%s/1*mat',mglGetSID));

    if length(files) >= 1
        if stimFileNum == -1
            if length(files) > 1
                warning('Multiple stimfiles found, loading last one. If you wanted a different functionality use stimFileNum=#');
            end
            fname = files(end).name;
        else
            fname = files(stimFileNum).name;
        end
        s = load(sprintf('~/data/cohcon/%s/%s',mglGetSID,fname));
        % copy staircases and run numbers
        stimulus.staircases = s.stimulus.staircases;
        stimulus.counter = s.stimulus.counter + 1;

        % load blocks too
        stimulus.runs = s.stimulus.runs;
        stimulus.runs.loaded = 1;

        clear s;
        stimulus.initStair = 0;
        disp(sprintf('(cohcon) Data file: %s loaded.',fname));
    end
end
disp(sprintf('(cohcon) This is run #%i',stimulus.counter));

if stimulus.plots==2
    dispInfoNum(stimulus);
    dispInfo(stimulus);
    return
end

%% Initialize Stimulus

myscreen = initStimulus('stimulus',myscreen);

if stimulus.scan
    stimulus.responseKeys = [1 2]; % corresponds to LEFT - RIGHT
else
    stimulus.responseKeys = [1 2]; % corresponds to LEFT - RIGHT
end

%% Block setup
if isfield(stimulus,'runs') && isfield(stimulus.runs,'loaded')
    % We already have our blocks
    stimulus.runs = rmfield(stimulus.runs,'loaded'); % remove the load field, otherwise it gets saved across runs
    if stimulus.counter > length(stimulus.runs.taskList)
        % double the length (maintains order)
        stimulus.runs.taskList  = [stimulus.runs.taskList stimulus.runs.taskBuilder]; %repmat(stimulus.runs.taskList,1,2);
    end
else
    % This is the first run, build up the blocks.
    stimulus.runs = struct;
    stimulus.runs.taskOpts = [1 2];
    if stimulus.scan
        stimulus.runs.taskBuild = {[1 1],[2 2]};
    else
        stimulus.runs.taskBuild = {[1 1 -1 1 -1] [2 2 -2 2 -2]};
    end
    stimulus.runs.taskOptsText = {'Motion','Contrast'};
    stimulus.runs.taskBuilder = [stimulus.runs.taskBuild{stimulus.runs.taskOpts(randperm(2))}];
    stimulus.runs.taskList = [stimulus.runs.taskBuilder];
end


%% Task Override
if overrideTask > 0
    stimulus.runs.curTask = overrideTask;
    if overrideTask>0
        stimulus.nocatch = 1;
        disp('(cohcon) Auto-setting nocatch!!');
    end
    % insert the override into the taskList
    pre = stimulus.runs.taskList(1:stimulus.counter-1);
    post = stimulus.runs.taskList(stimulus.counter:end);
    stimulus.runs.taskList = [pre overrideTask post];
else
    cT = stimulus.runs.taskList(stimulus.counter);
    if cT>0
        stimulus.nocatch = 1;
        disp('(cohcon) Auto-run setup is doing a no-catch run...');
    else
        disp('(cohcon) Auto-run setup is doing a main+catch run...');
    end
    stimulus.runs.curTask = abs(stimulus.runs.taskList(stimulus.counter));
end


%% Useful stimulus stuff

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% CONTROL BASE CONTRAST %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what contrast/coherence should run constantly in the background?
stimulus.baseCon = 0.25;
stimulus.baseCoh = 0;

% for initStair to know how many staircases to make
stimulus.stairInfo.catchP = 1;
stimulus.stairInfo.mainP = 1;

% for the real psychophysics experiment test we will run a discrimination
% task, i.e. "which side goes higher". Contrast will jump to 40% and
% coherence to 25%, but one side will go up more than the other side.
stimulus.stairInfo.pedestals.contrast = 0.4;
stimulus.stairInfo.pedestals.coherence = 0.3;
% for our distractors we will use these random increments:
stimulus.stairInfo.increments.coherence = [0 exp(-3:.33:-.6)];
stimulus.stairInfo.increments.contrast = [0 exp(-4:.33:-1.6)];
% initial thresholds to use in staircases
stimulus.stairInfo.initThresh.contrast = 0.25;
stimulus.stairInfo.initThresh.coherence = 0.85;

stimulus.stairInfo.pedOpts = {'coherence','contrast'};

if stimulus.scan
    stimulus.stairInfo.nocatchP = 4;
    % we are scanning, add more pedestals so we get the full range (these are now the same as for nocatch runs)
    stimulus.stairInfo.pedestals.contrast = [0.325 0.4 0.55 0.85];
    stimulus.stairInfo.pedestals.coherence = [0.15 0.3 0.45 0.6];
    
    if ~stimulus.nocatch
        disp('(cohcon) Auto-setting nocatch for scan run.');
        stimulus.nocatch = 1;
    end
    if stimulus.plots
        disp('(cohcon) Auto-setting no plots for scan run.');
        stimulus.plots =0 ;
    end
elseif stimulus.nocatch
    stimulus.stairInfo.nocatchP = 4;
    % we aren't scanning, so we can still run staircases, but we are doing
    % nocatch runs where there will be a LOT of runs. So let's add some
    % more pedestals in so we get a better estimate of the psychometric
    % function.
    stimulus.stairInfo.pedestals.contrast = [0.325 0.4 0.55 0.85];
    stimulus.stairInfo.pedestals.coherence = [0.15 0.3 0.45 0.6];
end

%% Colors
stimulus.colors.rmed = 127.5;

% We're going to add an equal number of reserved colors to the top and
% bottom, to try to keep the center of the gamma table stable.
stimulus.colors.reservedBottom = [0 0 0; .6 .6 .6]; % fixation cross colors
stimulus.colors.reservedTop = [.6 0 0; 0 .6 0]; % correct/incorrect colors
stimulus.colors.black = 0/255; stimulus.colors.white = 1/255;
stimulus.colors.red = 254/255; stimulus.colors.green = 255/255;
stimulus.colors.nReserved = 2; % this is /2 the true number, because it's duplicated
stimulus.colors.nUnreserved = 256-(2*stimulus.colors.nReserved);

stimulus.colors.mrmax = stimulus.colors.nReserved - 1 + stimulus.colors.nUnreserved;
stimulus.colors.mrmin = stimulus.colors.nReserved;

%% Everything else
stimulus.dots.xcenter = 0;
stimulus.dots.ycenter = 0;
stimulus.dots.dotsize = 3;
stimulus.dots.density = 21;
stimulus.dots.speed = 6;
stimulus.dots.centerOffset = 2;

stimulus.dotsR = stimulus.dots;
stimulus.dotsR.mult = 1;
stimulus.dotsL = stimulus.dots;
stimulus.dotsL.mult = -1;
stimulus = rmfield(stimulus,'dots');

stimulus.dotsR = initDotsRadial(stimulus.dotsR,myscreen);
stimulus.dotsL = initDotsRadial(stimulus.dotsL,myscreen);

if ~length(stimulus.baseCoh)==length(stimulus.baseCon)
    disp('(cohcon) Contrast and coherence have different pedestal #, this can cause errors');
end

%% Gamma Table Initialization

% get gamma table
if ~isfield(myscreen,'gammaTable')
  stimulus.linearizedGammaTable = mglGetGammaTable;
  disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
  disp(sprintf('(cuecon:initGratings) No gamma table found in myscreen. Contrast displays like this'));
  disp(sprintf('         should be run with a valid calibration made by moncalib for this monitor.'));
  disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
end
stimulus.linearizedGammaTable = myscreen.initScreenGammaTable;
    
%% Character textures
mglTextSet('Helvetica',32,stimulus.colors.white,0,0,0,0,0,0,0);
stimulus.text.mTexW = mglText('M');
stimulus.text.cTexW = mglText('C');

mglTextSet('Helvetica',32,stimulus.colors.black,0,0,0,0,0,0,0);
stimulus.text.mTexK = mglText('M');
stimulus.text.cTexK = mglText('C');
%% Setup Task

% This is the contrast change detection task
task{1}{1}.waitForBacktick = 1;

stimulus.seg.pre = 1;
stimulus.seg.stim = 2;
stimulus.seg.mask = 3;
stimulus.seg.ISI = 4;
stimulus.seg.resp = 5;
stimulus.seg.ITI = 6;
task{1}{1}.segmin = [0.2 0.5 0 .5 1 .2];
task{1}{1}.segmax = [0.2 0.5 0 1 1 .4];

if stimulus.scan
    task{1}{1}.segmin(stimulus.seg.ITI) = 2;
    task{1}{1}.segmax(stimulus.seg.ITI) = 11;
end

task{1}{1}.synchToVol = [0 0 0 0 0 0];
if stimulus.scan
    task{1}{1}.synchToVol(stimulus.seg.ITI) = 1;
end
task{1}{1}.getResponse = [0 0 0 0 0 0]; task{1}{1}.getResponse(stimulus.seg.resp)=1;
task{1}{1}.parameter.conSide = [1 2]; % 1 = left, 2 = right, the side will be the one with con/flow + delta (From staircase)
task{1}{1}.parameter.cohSide = [1 2];
task{1}{1}.parameter.dir = [-1 1];
task{1}{1}.parameter.conPedestal = 1:length(stimulus.stairInfo.pedestals.contrast); % target contrast
task{1}{1}.parameter.cohPedestal = 1:length(stimulus.stairInfo.pedestals.coherence); % target flow coherence
task{1}{1}.parameter.catch = [1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; % 15% chance of being a catch trial
task{1}{1}.random = 1;
task{1}{1}.numTrials = 65;

if stimulus.nocatch || stimulus.scan
    task{1}{1}.parameter.catch = -1;
end

%% Run variables

task{1}{1}.randVars.calculated.coherence = nan;
task{1}{1}.randVars.calculated.contrast = nan;
task{1}{1}.randVars.calculated.cohDelta = nan;
task{1}{1}.randVars.calculated.conDelta = nan;
%% Tracking

% these are variables that we want to track for later analysis.
task{1}{1}.randVars.calculated.task = nan; % Current task (calc per run)
task{1}{1}.randVars.calculated.correct = nan;
task{1}{1}.randVars.calculated.trialNum = nan;
task{1}{1}.randVars.calculated.lCoh = nan;
task{1}{1}.randVars.calculated.rCoh = nan;
task{1}{1}.randVars.calculated.lCon = nan;
task{1}{1}.randVars.calculated.rCon = nan;

stimulus.curTrial = 0;

if stimulus.scan
    % when scanning we add a 
    task{1}{1}.numTrials = inf;
    task{1}{2} = task{1}{1};
    task{1}{2}.waitForBacktick = 0;
    task{1}{1}.parameter.conPedestal = 1;
    task{1}{1}.parameter.cohPedestal = 1;
    task{1}{1}.numTrials = 1;
    task{1}{1}.getResponse = [0 0 0 0];
    task{1}{1}.segmin = [0 0 0 0 30];
    task{1}{1}.segmax = [0 0 0 0 30];
    task{1}{1}.parameter.catch = 0;
else
    % when scanning we add a 
    task{1}{2} = task{1}{1};
    task{1}{2}.waitForBacktick = 0;
    task{1}{1}.parameter.conPedestal = 1;
    task{1}{1}.parameter.cohPedestal = 1;
    task{1}{1}.numTrials = 1;
    task{1}{1}.getResponse = [0 0 0 0];
    task{1}{1}.segmin = [0 0 0 0 4];
    task{1}{1}.segmax = [0 0 0 0 4];
    task{1}{1}.parameter.catch = 0;
end
%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init staircase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if stimulus.initStair
    % We are starting our staircases from scratch
    disp(sprintf('(cohcon) Initializing staircases'));
    stimulus = initStaircase(stimulus);
else
    disp('(cohcon) Re-using staircase from previous run...');
    % Reset staircase if necessary
    checkStaircaseStop();
end

%% EYE CALIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if ~stimulus.scan
    myscreen = eyeCalibDisp(myscreen);
%end

%% Get Ready...
% clear screen    
setGammaTable_flowMax(1);
mglWaitSecs(1);
mglClearScreen(0.5);
if stimulus.scan        
    mglTextDraw('DO NOT MOVE',[0 1.5]);
    mglTextDraw(stimulus.runs.taskOptsText{stimulus.runs.curTask},[0 0]);
else

    mglTextDraw(stimulus.runs.taskOptsText{stimulus.runs.curTask},[0 0]);
end
mglFlush

tasktypes = {'Main+Catch','Nocatch'};
% let the user know
disp(sprintf('(cohcon) Starting run number: %i. Current task: %s, Type: %s',stimulus.counter,stimulus.runs.taskOptsText{stimulus.runs.curTask},tasktypes{stimulus.nocatch+1}));
% if stimulus.unattended
myscreen.flushMode = 1;
% end

%% Main Task Loop

phaseNum = 1;
% Again, only one phase.
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    % update the task
    [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);
    % flip screen
    myscreen = tickScreen(myscreen,task);
end

% task ended
mglClearScreen(0.5);
mglTextDraw('Run complete... please wait.',[0 0]);
mglFlush
myscreen.flushMode = 1;
mglWaitSecs(1);

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

if ~stimulus.scan
    dispInfoNum(stimulus);
    if stimulus.plots
        disp('(cohcon) Displaying plots');
        dispInfo(stimulus);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startTrialCallback(task,myscreen)
%%

global stimulus

if task.thistrial.thisphase==2 && stimulus.scan
    task.thistrial.seglen(end) = 1.05^(rand*30+20);
end

stimulus.curTrial = stimulus.curTrial + 1;

%  Set the current task

if task.thistrial.catch > 0
    switchTasks = [2 1];
    task.thistrial.task = switchTasks(stimulus.runs.curTask);
    % edit seglen
    task.thistrial.seglen(stimulus.seg.mask) = .75;
    task.thistrial.seglen(stimulus.seg.ISI) = .75;
    task.thistrial.seglen(stimulus.seg.resp) = 2.5;
    disp('(cohcon) Catch trial.');
else
    task.thistrial.task = stimulus.runs.curTask;
end

% Set the missing thistrial vars
task.thistrial.coherence = stimulus.stairInfo.pedestals.coherence(task.thistrial.cohPedestal);
task.thistrial.contrast = stimulus.stairInfo.pedestals.contrast(task.thistrial.conPedestal);
task.thistrial.trialNum = stimulus.curTrial;

% Get the pedestals
[cohTh, conTh, stimulus] = getDeltaPed(task,stimulus);


% Reduce if pedestals are too large
if (task.thistrial.coherence + cohTh) > 1
    cohTh = 1 - task.thistrial.coherence;
end
if (task.thistrial.contrast + conTh) > 1
    conTh = 1 - task.thistrial.contrast;
end

if stimulus.stablecon
    conTh = 0;
end
if stimulus.stablecoh
    cohTh = 0;
end
% Save info
task.thistrial.conDelta = conTh;
task.thistrial.cohDelta = cohTh;


if task.thistrial.conSide==1
    task.thistrial.lCon = task.thistrial.contrast+task.thistrial.conDelta;
    task.thistrial.rCon = task.thistrial.contrast;
else
    task.thistrial.rCon = task.thistrial.contrast+task.thistrial.conDelta;
    task.thistrial.lCon = task.thistrial.contrast;
end

if task.thistrial.cohSide==1
    task.thistrial.lCoh = task.thistrial.coherence+task.thistrial.cohDelta;
    task.thistrial.rCoh = task.thistrial.coherence;
else
    task.thistrial.rCoh = task.thistrial.coherence+task.thistrial.cohDelta;
    task.thistrial.lCoh = task.thistrial.coherence;
end

disp(sprintf('(cohcon) Trial %i starting. Coherence: L %.02f; R %.02f Contrast: L %.02f; R %.02f',task.thistrial.trialNum,...
    task.thistrial.lCoh,task.thistrial.rCoh,...
    task.thistrial.lCon,task.thistrial.rCon));

% set the gammaTable for this trial
setGammaTable_flowMax(task.thistrial.contrast + task.thistrial.conDelta);


% set directions
stimulus.dotsL.dir = task.thistrial.dir;
stimulus.dotsR.dir = task.thistrial.dir;

myscreen.flushMode = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)
%%

global stimulus

stimulus.live.mt = 0;
switch task.thistrial.thisseg
    case stimulus.seg.pre
        stimulus.live.dots = 0;
        stimulus.live.fixColor = stimulus.colors.white;
        stimulus.live.catchFix = 0;
    case stimulus.seg.ITI
        stimulus.live.dots = 0;
        stimulus.live.fixColor = stimulus.colors.black;
        stimulus.live.catchFix = 0;
    case stimulus.seg.stim
        stimulus.live.dots = 1;
        stimulus.live.fixColor = stimulus.colors.black;
        stimulus.live.catchFix = 0;
    case stimulus.seg.mask
        stimulus.live.dots = 0;
        stimulus.live.fixColor = stimulus.colors.black;
        stimulus.live.catchFix = 0;
    case stimulus.seg.ISI
        stimulus.live.dots = 0;
        stimulus.live.fixColor = stimulus.colors.black;
        stimulus.live.catchFix = 1;
    case stimulus.seg.resp
        stimulus.live.dots = 0;
        stimulus.live.fixColor = stimulus.colors.white;
        stimulus.live.catchFix = 1;
end

if stimulus.constant
    stimulus.live.dots=1;
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus


mglClearScreen(0.5);

if stimulus.live.dots==1, stimulus = upDots(task,stimulus,myscreen); end
upFix(task,stimulus);

function upFix(task,stimulus)
%%

if ~(task.thistrial.catch > 0) || stimulus.live.catchFix == 0
    mglFixationCross(1.5,1.5,stimulus.live.fixColor);
else
    if stimulus.runs.curTask==2
        if stimulus.live.fixColor==stimulus.colors.white
            mglBltTexture(stimulus.text.mTexW,[0 0]);
        else
            mglBltTexture(stimulus.text.mTexK,[0 0]);
        end
    else
        if stimulus.live.fixColor==stimulus.colors.white
            mglBltTexture(stimulus.text.cTexW,[0 0]);
        else
            mglBltTexture(stimulus.text.cTexK,[0 0]);
        end
    end
end

function stimulus = upDots(task,stimulus,myscreen)

if task.thistrial.thisseg==stimulus.seg.stim
    lcoh = task.thistrial.lCoh;
    rcoh = task.thistrial.rCoh;
    lcon = task.thistrial.lCon;
    rcon = task.thistrial.rCon;
else
    lcoh = stimulus.baseCoh;
    rcoh = stimulus.baseCoh;
    lcon = stimulus.baseCon;
    rcon = stimulus.baseCon;
end

%% Old update code start here
stimulus.dotsL = updateDotsRadial(stimulus.dotsL,lcoh,myscreen,true);
stimulus.dotsR = updateDotsRadial(stimulus.dotsR,rcoh,myscreen,true);

% Correct values for gamma table adjustments
lcon = lcon / stimulus.curMaxContrast;
rcon = rcon / stimulus.curMaxContrast;

% dotsR
% update +contrast

rcon = adjustConToTable(rcon,stimulus)/2;
lcon = adjustConToTable(lcon,stimulus)/2;

mglPoints2(stimulus.dotsR.xdisp(stimulus.dotsR.con==1),stimulus.dotsR.ydisp(stimulus.dotsR.con==1),...
    stimulus.dotsR.dotsize,[.5 .5 .5] - rcon);
% update - contrast
mglPoints2(stimulus.dotsR.xdisp(stimulus.dotsR.con==2),stimulus.dotsR.ydisp(stimulus.dotsR.con==2),...
    stimulus.dotsR.dotsize,[.5 .5 .5] + rcon);
% dotsL
% update +contrast
mglPoints2(stimulus.dotsL.xdisp(stimulus.dotsL.con==1),stimulus.dotsL.ydisp(stimulus.dotsL.con==1),...
    stimulus.dotsL.dotsize,[.5 .5 .5] - lcon);
% update - contrast
mglPoints2(stimulus.dotsL.xdisp(stimulus.dotsL.con==2),stimulus.dotsL.ydisp(stimulus.dotsL.con==2),...
    stimulus.dotsL.dotsize,[.5 .5 .5] + lcon);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Adjust contrast to the gamma table %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function conValue = adjustConToTable(conValue,stimulus)
conValue = conValue * stimulus.colors.nUnreserved / 256;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus

responseText = {'Incorrect','Correct'};
responsePos = {'Left','Right'};
fixColors = {stimulus.colors.red,stimulus.colors.green};

if any(task.thistrial.whichButton == stimulus.responseKeys)
    if task.thistrial.gotResponse == 0
        if task.thistrial.task==1
            cSide = task.thistrial.cohSide;
        else
            cSide = task.thistrial.conSide;
        end
        task.thistrial.correct = task.thistrial.whichButton == stimulus.responseKeys(cSide);
        % Store whether this was correct
        stimulus.live.fixColor = fixColors{task.thistrial.correct+1};
        disp(sprintf('(cohcon) Response %s: %s',responseText{task.thistrial.correct+1},responsePos{find(stimulus.responseKeys==task.thistrial.whichButton)}));
        if ~(task.thistrial.catch > 0)
            if ~stimulus.nocatch
                stimulus.staircases.main{task.thistrial.task,curPedValue(task,false)} = ...
                    doStaircase('update',stimulus.staircases.main{task.thistrial.task,curPedValue(task,false)},task.thistrial.correct);
            else
                stimulus.staircases.nocatch{task.thistrial.task,curPedValue(task,false)} = ...
                    doStaircase('update',stimulus.staircases.nocatch{task.thistrial.task,curPedValue(task,false)},task.thistrial.correct);
            end
        else
            stimulus.live.fixColor = stimulus.colors.black; % we never show information about catch trials
            stimulus.live.catchFix = 0;
            stimulus.staircases.catch{task.thistrial.task,curPedValue(task,true)} = ...
                doStaircase('update',stimulus.staircases.catch{task.thistrial.task,curPedValue(task,true)},task.thistrial.correct);
        end
    else
        disp(sprintf('(cohcon) Subject responded multiple times: %i',task.thistrial.gotResponse+1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              HELPER FUNCTIONS                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
%    getDeltaPed       %
%%%%%%%%%%%%%%%%%%%%%%%%

function [cohPed, conPed, stimulus] = getDeltaPed(task,stimulus)
%%
if stimulus.runs.curTask == 1
    % COHERENCE MAIN TASK    
    if ~stimulus.nocatch
        [cohPed, stimulus.staircases.main{1,curPedVal(task,1)}] = doStaircase('testValue',stimulus.staircases.main{1,curPedVal(task,1)});
    else
        [cohPed, stimulus.staircases.nocatch{1,curPedVal(task,1)}] = doStaircase('testValue',stimulus.staircases.nocatch{1,curPedVal(task,1)});
    end
    
    if task.thistrial.catch > 0
        [conPed, stimulus.staircases.catch{2,curPedVal(task,2)}] = doStaircase('testValue',stimulus.staircases.catch{2,curPedVal(task,2)});
    else
        conPed = stimulus.stairInfo.increments.contrast(randi(length(stimulus.stairInfo.increments.contrast)));
    end
else
    % CONTRAST MAIN TASK
    if ~stimulus.nocatch
        [conPed, stimulus.staircases.main{2,curPedVal(task,2)}] = doStaircase('testValue',stimulus.staircases.main{2,curPedVal(task,2)});
    else
        [conPed, stimulus.staircases.nocatch{2,curPedVal(task,2)}] = doStaircase('testValue',stimulus.staircases.nocatch{2,curPedVal(task,2)});
    end
    
    if task.thistrial.catch > 0
        [cohPed, stimulus.staircases.catch{1,curPedVal(task,1)}] = doStaircase('testValue',stimulus.staircases.catch{1,curPedVal(task,1)});
    else
        cohPed = stimulus.stairInfo.increments.coherence(randi(length(stimulus.stairInfo.increments.coherence)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%
%    curPedValue       %
%%%%%%%%%%%%%%%%%%%%%%%%

function ped = curPedValue(task,iscatch)
%%
switcher = [2 1];
if iscatch
    usetask = switcher(task.thistrial.task);
else
    usetask = task.thistrial.task;
end

ped = curPedVal(task,usetask);

function ped = curPedVal(task,taskNum)
%%
if taskNum==1
    ped = task.thistrial.cohPedestal;
else
    ped = task.thistrial.conPedestal;
end


%%%%%%%%%%%%%%%%%%%%%%%%
%    initStaircase     %
%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initStaircase(stimulus)
% we're going to be fucking organized this time and put all the staircases
% in one place.... duh.
stimulus.staircases = struct;
%%
stimulus.staircases.catch = cell(2,stimulus.stairInfo.catchP); % task first, pedestal second
stimulus.staircases.main = cell(2,stimulus.stairInfo.mainP);
stimulus.staircases.nocatch = cell(2,stimulus.stairInfo.nocatchP);

if size(stimulus.staircases.catch,2)>1, error('Staircase Initialization Failure'); end
% Catch && Main staircases
for task = 1:2
    stimulus.staircases.catch{task,1} = doStaircase('init','upDown',...
        'initialThreshold',stimulus.stairInfo.initThresh.(stimulus.stairInfo.pedOpts{task}),...
        'initialStepsize',stimulus.stairInfo.initThresh.(stimulus.stairInfo.pedOpts{task})/3,...
        'minThreshold=0.001','maxThreshold=2','stepRule','pest',...
        'nTrials=50','maxStepsize=0.2','minStepsize=0.001');
    stimulus.staircases.main{task,1} = doStaircase('init','upDown',...
        'initialThreshold',stimulus.stairInfo.initThresh.(stimulus.stairInfo.pedOpts{task}),...
        'initialStepsize',stimulus.stairInfo.initThresh.(stimulus.stairInfo.pedOpts{task})/3,...
        'minThreshold=0.001','maxThreshold=2','stepRule','pest',...
        'nTrials=50','maxStepsize=0.2','minStepsize=0.001');
end


% NoCatch staircases: Warning, these have different sizes in scan sessions

% motion first then contrast
for task = 1:2
    stimulus.staircases.nocatch{task,1} = doStaircase('init','upDown',...
        'initialThreshold',stimulus.stairInfo.initThresh.(stimulus.stairInfo.pedOpts{task}),...
        'initialStepsize',stimulus.stairInfo.initThresh.(stimulus.stairInfo.pedOpts{task})/3,...
        'minThreshold=0.001','maxThreshold=2','stepRule','pest',...
        'nTrials=50','maxStepsize=0.2','minStepsize=0.001');
    for p = 2:size(stimulus.staircases.nocatch,2)
        stimulus.staircases.nocatch{task,p} = stimulus.staircases.nocatch{task,1};
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%
%    dispInfo    %
%%%%%%%%%%%%%%%%%%%%%%%

function dispInfoNum(stimulus)
%%
trials = 0;
nmain = 0;
ncatch = 0;
ncontrol = 0;
for t = 1:2
    cmain = stimulus.staircases.main{t};
    for i = 1:length(cmain)
        trials = trials + cmain(i).trialNum;
        nmain = nmain + cmain(i).trialNum;
    end
    ccatch = stimulus.staircases.catch{t};
    for i = 1:length(ccatch)
        trials = trials + ccatch(i).trialNum;
        ncatch = ncatch + ccatch(i).trialNum;
    end
    
    for n = 1:4
        cno = stimulus.staircases.nocatch{t,n};
        for i = 1:length(cno)
            trials = trials + cno(i).trialNum;
            ncontrol = ncontrol + cno(i).trialNum;
        end
    end
end
disp(sprintf('(dispInfo) Subject %s has completed %i trials so far. %i main, %i catch, %i control.',mglGetSID,trials,nmain,ncatch,ncontrol));
function dispInfo(stimulus)
%%
try
    %% No-Catch Performance
    nocatch = zeros(2,4);
    nocatchs = zeros(2,4);
    
    for task = 1:size(nocatch,1)
        for ped = 1:size(nocatch,2)
            try
            out = doStaircase('threshold',stimulus.staircases.nocatch{task,ped},'type','weibull','dispFig=0');
            
%             out = doStaircase('threshold',stimulus.staircases.nocatch{task,ped});
            nocatch(task,ped) = out.threshold;
            nocatchs(task,ped) = out.thresholdSTE;
            catch % probably a missing staircase, whatever
            end
        end
    end
    if any(nocatch(:)>10)
        warning('removing stairs >10');
        nocatch(nocatch>10) = 0;
    end
    %% Main + Catch
    main = zeros(2,1);
    mains = zeros(2,1);
    catch_ = zeros(2,1);
    catch_s = zeros(2,1);
    
    for task = 1:2
        try
            out = doStaircase('threshold',stimulus.staircases.main{task},'type','weibull','dispFig=0');
%             out = doStaircase('threshold',stimulus.staircases.main{task});
            main(task) = out.threshold;
            mains(task) = out.thresholdSTE;
        catch % missing a staircase
        end
        try 
            out = doStaircase('threshold',stimulus.staircases.catch{task},'type','weibull','dispFig=0');
%             out = doStaircase('threshold',stimulus.staircases.catch{task});
            catch_(task) = out.threshold;
            catch_s(task) = out.thresholdSTE;
        catch
        end
    end
        %% The Plot!
    map = brewermap(6,'PuOr');
    figure, hold on
    h2 = plot([0.15 0.3 0.45 0.6],nocatch(1,:),'-','Color',map(1,:));
%     errbar([0.15 0.3 0.45 0.6],nocatch(1,:),nocatchs(1,:),'Color',map(1,:));
    h6 = plot([0.325 0.4 0.55 0.85],nocatch(2,:),'-','Color',map(6,:));
%     errbar([0.325 0.4 0.55 0.85],nocatch(2,:),nocatchs(2,:),'Color',map(6,:));
    
    
    h1 = plot([0.15 0.3 0.45 0.6],nocatch(1,:),'o','MarkerSize',15);
    set(h1(1),'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',map(1,:),'LineWidth',1.5);
    
    h3 = plot(0.305,main(1),'o','MarkerSize',15);
    set(h3(1),'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',map(2,:),'LineWidth',1.5);
%     errbar(0.305,main(1),mains(1),'Color',map(2,:));
    h7 = plot(0.305,catch_(1),'o','MarkerSize',15);
%     errbar(0.305,catch_(1),catch_s(1),'Color',map(3,:));
    set(h7(1),'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',map(3,:),'LineWidth',1.5);
    h5 = plot([0.325 0.4 0.55 0.85],nocatch(2,:),'o','MarkerSize',15);
    set(h5(1),'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',map(6,:),'LineWidth',1.5);
    h9 = plot(0.4,main(2),'o','MarkerSize',15);
%     errbar(0.4,main(2),mains(2),'Color',map(5,:));
    set(h9(1),'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',map(5,:),'LineWidth',1.5);
    
    h4 = plot(0.4,catch_(2),'o','MarkerSize',15);
%     errbar(0.4,catch_(2),catch_s(2),'Color',map(4,:));
    set(h4(1),'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',map(4,:),'LineWidth',1.5);
    %    ;
    
    
    legend([h1,h3,h7,h5,h9,h4],{'Coherence: Control','Coherence: Attended','Coherence: Unattended','Contrast: Control','Contrast: Attended','Contrast: Unattended'});
    
    title('Psychometric Functions for Cohcon');
    xlabel('Contrast/Coherence (%)');
    ylabel('Threshold (%)');
    a = axis();
    if a(4)>1
        axis([0 1 0 1]);
    else
        axis([0 1 0 a(4)]);
    end
    drawPublishAxis
%     set(h(1),'MarkerEdgeColor','r','MarkerFaceColor','none')
catch
    disp('(cohcon) Staircase figure was not generated successfully.');
end

%% checkStaircaseStop
function checkStaircaseStop()
global stimulus

for task = 1:2
    stimulus.staircases.catch{task,1} = resetStair(stimulus.staircases.catch{task,1});
end
% Check both staircases
for task = 1:2
    stimulus.staircases.main{task,1} = resetStair(stimulus.staircases.main{task,1});
end
% Check nocatch staircases

for task = 1:2
    for p = 1:size(stimulus.staircases.nocatch,2)
        stimulus.staircases.nocatch{task,p} = resetStair(stimulus.staircases.nocatch{task,p});
    end
end

function s = resetStair(s)

if isempty(s)
    return
end

if doStaircase('stop',s)
    % this is a bit of a pain... you can't pass an initialThreshold
    % argument do doStaircase('init',s, ...), it ignores everything and
    % resets using the calculated threshold. Because you can't override it
    [args, vals, ~] = getArgs(s(1).initArgs);
    threshPos = -1;
    stepPos = -1;
    for i = 1:length(args)
        switch args{i}
            case 'initialThreshold'
                threshPos = i;
            case 'initialStepsize'
                stepPos = i;
        end
    end
    out = doStaircase('threshold',s);
    in = input(sprintf('Resetting Staircase... Estimate is: %1.2f. Reset ([Y]/[C]ustom/[O]riginal): ',out.threshold),'s');
    switch in
        case 'Y'
            vals{threshPos} = out.threshold;
            vals{stepPos} = out.threshold / 3;
        case 'C'
            disp('Original values:');
            disp(sprintf('%s: %0.2f',args{threshPos},num2str(vals{threshPos})));
            val = str2double(input('New threshold value: ','s'));
            vals{threshPos} = val;
            vals{stepPos} = val / 3;
        case 'O'
    end
    if ~length(args) == 8
        disp('Args incorrect length...');
        keyboard
    end
    %             stimulus.staircase{task,ped}(end+1) = doStaircase('init',s,'initialThreshold',vals{threshPos},'initialStepsize',vals{stepPos});
    s(end+1) = doStaircase('init','upDown',args{1},vals{1},args{2},vals{2},args{3},vals{3},args{4},vals{4},args{5},vals{5},args{6},vals{6},args{7},vals{7},args{8},vals{8});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create dots for optic flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initDotsRadial(dots,~)

% maximum depth of points
dots.minX = 3.5;
dots.maxX = 12;
dots.minY = -7;
dots.maxY = 7;

dots.dir = 1;

area = (dots.maxX-dots.minX)*(dots.maxY-dots.minY);

dots.n = area * dots.density;

% make a some points
% dots.n = 500*dots.density;
% make sure it's an even number
dots.n = dots.n + mod(dots.n,2);

% set half to white and half to black
dots.con = repmat([1 2],1,dots.n/2);

dots.x = rand(1,dots.n)*(dots.maxX-dots.minX)+dots.minX;
dots.y = rand(1,dots.n)*abs(dots.maxY-dots.minY)+dots.minY;

dots.xdisp = dots.mult*dots.x;
dots.ydisp = dots.y;

% set incoherent dots to 0
dots.coherency = 1;
dots.incoherent = rand(1,dots.n) > dots.coherency;
dots.incoherentn = sum(dots.incoherent);
dots.coherent = ~dots.incoherent;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots for Radial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updateDotsRadial(dots,coherence,myscreen,repick)

% stuff to compute median speed
dots.oldx = dots.x;
dots.oldy = dots.y;

% get the coherent and incoherent dots
if repick
    dots.incoherent = rand(1,dots.n) > coherence;
    dots.incoherentn = sum(dots.incoherent);
    dots.coherent = ~dots.incoherent;
    dots.coherency = coherence;
elseif dots.coherency ~= coherence
    cohDiff = coherence - dots.coherency;
    numDots = round(abs(cohDiff) * dots.n); % actual number of dots to flip
    if numDots > dots.n, numDots = dots.n; end
    if cohDiff > 0
        % we need to add more coherent dots
        flipDots = [zeros(1,numDots) ones(1,sum(dots.incoherent)-numDots)];
        dots.incoherent(dots.incoherent) = flipDots(randperm(length(flipDots)));
    else
        % we need to add more incoherent dots
        flipDots = [ones(1,numDots) zeros(1,sum(dots.coherent)-numDots)];
        dots.incoherent(dots.coherent) = flipDots(randperm(length(flipDots)));
    end
    dots.incoherentn = sum(dots.incoherent);
    dots.coherent = ~dots.incoherent;
    dots.coherency = sum(dots.coherent)/dots.n;
end
dots.coherentn = dots.n-dots.incoherentn;

freq_factor = dots.speed/myscreen.framesPerSecond;

% move coherent dots
dots.x(dots.coherent) = dots.x(dots.coherent) + dots.dir*freq_factor;

% these are for flipping into the other quadrants
xmat = repmat([1 1 -1 -1],1,dots.incoherentn+4-mod(dots.incoherentn,4));
ymat = repmat([1 -1 1 -1],1,dots.incoherentn+4-mod(dots.incoherentn,4));
perms = randperm(dots.incoherentn);

% move incoherent dots
% get random vectors
dots.rX = rand(1,dots.incoherentn);
dots.rY = sqrt(1-(dots.rX.^2));
dots.rX = (dots.rX .* xmat(perms)) .* (freq_factor*(1+randn(1,dots.incoherentn)/3)); % rescale to match the velocity
dots.rY = (dots.rY .* ymat(perms)) .* (freq_factor*(1+randn(1,dots.incoherentn)/3));
% dots.rX = (dots.rX .* xmat) * freq_factor .* ((1.75*rand(1,dots.incoherentn)).^2); % rescale to match the velocity
% dots.rY = (dots.rY .* ymat) * freq_factor .* ((1.75*rand(1,dots.incoherentn)).^2);
dots.x(dots.incoherent) = dots.x(dots.incoherent) + dots.rX;
dots.y(dots.incoherent) = dots.y(dots.incoherent) + dots.rY;

offscreen = dots.x > dots.maxX;
dots.x(offscreen) = dots.x(offscreen) - abs(dots.maxX - dots.minX);
offscreen = dots.x < dots.minX;
dots.x(offscreen) = dots.x(offscreen) + abs(dots.maxX - dots.minX);

offscreen = dots.y > dots.maxY;
dots.y(offscreen) = dots.y(offscreen) - abs(dots.maxY - dots.minY);
offscreen = dots.y < dots.minY;
dots.y(offscreen) = dots.y(offscreen) + abs(dots.maxY - dots.minY);

dots.xdisp = dots.mult*dots.x;
dots.ydisp = dots.y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sets the gamma table so that we can have
% finest possible control over the stimulus contrast.
%
% stimulus.colors.reservedColors should be set to the reserved colors (for cue colors, etc).
% maxContrast is the maximum contrast you want to be able to display.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setGammaTable_flowMax(maxContrast)

global stimulus;

% set the bottom
gammaTable(1:size(stimulus.colors.reservedBottom,1),1:size(stimulus.colors.reservedBottom,2)) = stimulus.colors.reservedBottom;

% set the gamma table
if maxContrast == 1
    % create the rest of the gamma table
    cmax = 1;cmin = 0;
    luminanceVals = cmin:((cmax-cmin)/(stimulus.colors.nUnreserved-1)):cmax;

    % now get the linearized range
    redLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,luminanceVals,'linear');
    greenLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,luminanceVals,'linear');
    blueLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,luminanceVals,'linear');
elseif maxContrast > 0
    % create the rest of the gamma table
    cmax = 0.5+maxContrast/2;cmin = 0.5-maxContrast/2;
    luminanceVals = cmin:((cmax-cmin)/(stimulus.colors.nUnreserved-1)):cmax;

    % now get the linearized range
    redLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,luminanceVals,'linear');
    greenLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,luminanceVals,'linear');
    blueLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,luminanceVals,'linear');
else
    % if we are asked for 0 contrast then simply set all the values to gray
    redLinearized = repmat(.5,1,stimulus.colors.nUnreserved);
    greenLinearized = repmat(.5,1,stimulus.colors.nUnreserved);
    blueLinearized = repmat(.5,1,stimulus.colors.nUnreserved);
end

% add to the table!
gammaTable((stimulus.colors.mrmin:stimulus.colors.mrmax)+1,:)=[redLinearized;greenLinearized;blueLinearized]';

% set the top
gammaTable = [gammaTable; stimulus.colors.reservedTop];

if size(gammaTable,1)~=256
    disp('(setGammaTable) Failure: Incorrect number of colors in gamma table produced');
    keyboard
end

% set the gamma table
mglSetGammaTable(gammaTable);

% remember what the current maximum contrast is that we can display
stimulus.curMaxContrast = maxContrast;
