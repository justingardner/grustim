
% cohCon
%
%      usage: myscreen=coherentContrast()
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
%             mtloc (0/1) - Runs an mt localizer instead of the actual task
%             scan (0/1) - Scanner timing

function [myscreen] = coherentContrast(varargin)

global stimulus
clear fixStimulus
global fixStimulus
%% Initialize Variables

% add arguments later
stimFileNum = [];
unattended = [];
plots = [];
overrideTask = [];
projector = [];
mtloc = [];
scan = [];
getArgs(varargin,{'stimFileNum=-1','unattended=0', 'mtloc=0'...
    'plots=1','overrideTask=0','projector=0','scan=0'});
stimulus.projector = projector;
stimulus.mtloc = mtloc;
stimulus.unattended = unattended;
stimulus.scan = scan;
stimulus.plots = plots;

% if (stimulus.projector && ~stimulus.scan) || (stimulus.scan && ~ stimulus.projector)
%     warning('Running in scan mode or projector mode without the other... are you sure that''s what you wanted?');
%     keyboard
% end

if stimulus.mtloc && (~stimulus.unattended || ~stimulus.scan)
    warning('Running mtlocalizer attended or without scan... are you sure that''s what you wanted?');
    keyboard
end

if stimulus.mtloc
    stimulus.mt = struct;
end

stimulus.counter = 1; % This keeps track of what "run" we are on.
%% Setup Screen

myscreen = initScreen();

if stimulus.projector, stimulus.stencil = mglProjStencil(); end

%% Open Old Stimfile
stimulus.initStair = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/coherentContrast/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/coherentContrast/%s/1*mat',mglGetSID));

    if length(files) >= 1
        if stimFileNum == -1
            if length(files) > 1
                warning('Multiple stimfiles found, loading last one. If you wanted a different functionality use stimFileNum=#');
            end
            fname = files(end).name;
        else
            fname = files(stimFileNum).name;
        end
        s = load(sprintf('~/data/coherentContrast/%s/%s',mglGetSID,fname));
        stimulus.staircase = s.stimulus.staircase;
        stimulus.stairCatch = s.stimulus.stairCatch;
        stimulus.counter = s.stimulus.counter + 1;

        % load blocks too
        stimulus.runs = s.stimulus.runs;
        stimulus.runs.loaded = 1;

        clear s;
        stimulus.initStair = 0;
        disp(sprintf('(cohCon) Data file: %s loaded, this is run #%i',fname,stimulus.counter));
    end
end

%% Initialize Stimulus

myscreen = initStimulus('stimulus',myscreen);

stimulus.responseKeys = [9 10]; % corresponds to LEFT - RIGHT

% Sigmoid
x = 0:.015:1;

y = 1 ./ (1 + exp(-8 * (x - 0.5)));
y = y - min(y);
y = y ./ sum(y);
y = y ./ max(y);

y = [zeros(1,101-length(y)) y];

stimulus.sigmoid = y;
stimulus.sigmoidMu = mean(y);

%% Colors
stimulus.colors.rmed = 127.75;

% We're going to add an equal number of reserved colors to the top and
% bottom, to try to keep the center of the gamma table stable.
stimulus.colors.reservedBottom = [1 1 1; 0 0 0]; % fixation cross colors
stimulus.colors.reservedTop = [1 0 0; 0 1 0]; % correct/incorrect colors
stimulus.colors.black = 1/255; stimulus.colors.white = 0/255;
stimulus.colors.red = 254/255; stimulus.colors.green = 255/255;
stimulus.colors.nReserved = 2; % this is /2 the true number, because it's duplicated
stimulus.colors.nUnreserved = 256-(2*stimulus.colors.nReserved);

stimulus.colors.mrmax = stimulus.colors.nReserved - 1 + stimulus.colors.nUnreserved;
stimulus.colors.mrmin = stimulus.colors.nReserved;

%% Everything else
stimulus.dots.xcenter = 0;
stimulus.dots.ycenter = 0;
stimulus.dots.dotsize = 4;
stimulus.dots.density = 3;
stimulus.dots.speed = 3.25;
stimulus.dots.centerOffset = 2;

if stimulus.mtloc
    stimulus.dots.dotsize = 4;
    stimulus.dots.density = 3;
end

stimulus.dotsR = stimulus.dots;
stimulus.dotsR.mult = 1;
stimulus.dotsL = stimulus.dots;
stimulus.dotsL.mult = -1;
stimulus = rmfield(stimulus,'dots');

stimulus.pedestals.pedOpts = {'coherence','contrast'};

stimulus.pedestals.coherence = [.05 .125 .25 .45];
stimulus.pedestals.contrast = exp(-1.75:(1.25/3):-.5);

stimulus.pedestals.initThresh.coherence = .5;
stimulus.pedestals.initThresh.contrast = .2;

stimulus.pedestals.catch.coherence = [.05 .1 .175 .275 .4 .5];
stimulus.pedestals.catch.contrast = exp([-3.3 -3 -2.7 -2.4 -2.1 -1.8]);

if stimulus.mtloc
   stimulus.pedestals.coherence = [0 1];
   stimulus.pedestals.contrast = 1;
end

stimulus.dotsR = initDotsRadial(stimulus.dotsR,myscreen);
stimulus.dotsL = initDotsRadial(stimulus.dotsL,myscreen);

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

stimulus.seg.ITI = 1; % the ITI is either 20s (first time) or 1s
stimulus.seg.rampUP = 2;
stimulus.seg.stim = 3;
stimulus.seg.rampDOWN = 4;
stimulus.seg.ISI = 5;
stimulus.seg.resp = 6;
task{1}{1}.segmin = [.3 .55 .4 .55 .1 1];
task{1}{1}.segmax = [.7 .55 .4 .55 .4 1];

if stimulus.unattended
    task{1}{1}.segmin(stimulus.seg.ITI) = 1;
    task{1}{1}.segmax(stimulus.seg.ITI) = 2;
    % remove response
    task{1}{1}.segmin(stimulus.seg.ISI) = 0;
    task{1}{1}.segmax(stimulus.seg.ISI) = 0;
    task{1}{1}.segmin(stimulus.seg.resp) = 0;
    task{1}{1}.segmax(stimulus.seg.resp) = 0;
end
if stimulus.scan
    task{1}{1}.segmin(stimulus.seg.ITI) = 1;
    task{1}{1}.segmax(stimulus.seg.ITI) = 12;
end
if stimulus.mtloc
    % change stim length
    task{1}{1}.segmin(stimulus.seg.stim) = 12;
    task{1}{1}.segmax(stimulus.seg.stim) = 12;
    task{1}{1}.segmin(stimulus.seg.ITI) = 0;
    task{1}{1}.segmax(stimulus.seg.ITI) = 0;
end

task{1}{1}.synchToVol = [0 0 0 0 0 0];
task{1}{1}.getResponse = [0 0 0 0 0 1];
task{1}{1}.parameter.side = [1 2]; % 1 = left, 2 = right, the side will be the one with con/flow + delta (From staircase)
task{1}{1}.parameter.dir = [-1 1];
task{1}{1}.parameter.conPedestal = [1 2 3 4]; % target contrast
task{1}{1}.parameter.cohPedestal = [1 2 3 4]; % target flow coherence
task{1}{1}.parameter.catch = [1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; % 15% chance of being a catch trial
task{1}{1}.random = 1;
task{1}{1}.numTrials = 120;

if stimulus.scan
    task{1}{1}.numTrials = 40;
end
if stimulus.unattended
    task{1}{1}.getResponse = [0 0 0 0 0 0];
end
if stimulus.mtloc
    task{1}{1}.parameter.conPedestal = 1;
    task{1}{1}.parameter.cohPedestal = [1 2];
    task{1}{1}.parameter.catch = 0;
    task{1}{1}.numTrials = 22;
end

%% Run variables
task{1}{1}.randVars.calculated.task = nan; % Current task (calc per trial)
task{1}{1}.randVars.calculated.deltaPed = nan; % Current task (calc per trial)
task{1}{1}.randVars.calculated.coherence = nan;
task{1}{1}.randVars.calculated.contrast = nan;
task{1}{1}.randVars.calculated.trialNum = nan;
task{1}{1}.randVars.calculated.avgCohL = nan;
task{1}{1}.randVars.calculated.avgCohR = nan;
task{1}{1}.randVars.calculated.avgConL = nan;
task{1}{1}.randVars.calculated.avgConR = nan;

%% Tracking

% these are variables that we want to track for later analysis.
task{1}{1}.randVars.calculated.correct = nan;
task{1}{1}.randVars.calculated.trialNum = nan;

stimulus.curTrial = 0;

%% Block setup
if isfield(stimulus,'runs') && isfield(stimulus.runs,'loaded')
    % We already have our blocks
    stimulus.runs = rmfield(stimulus.runs,'loaded'); % remove the load field, otherwise it gets saved across runs
    if stimulus.counter > length(stimulus.runs.taskList)
        stimulus.runs.taskList  = repmat(stimulus.runs.taskList,1,2);
    end
else
    % This is the first run, build up the blocks.
    stimulus.runs = struct;
    stimulus.runs.taskOpts = [1 2];
    stimulus.runs.taskOptsText = {'Motion','Contrast'};
    stimulus.runs.taskList = stimulus.runs.taskOpts(randperm(2));
end
if overrideTask > 0
    stimulus.runs.curTask = overrideTask;
else
    stimulus.runs.curTask = stimulus.runs.taskList(stimulus.counter);
end

%% Unattended Mode
if stimulus.unattended
    fixStimulus.diskSize = 0;
    fixStimulus.stimColor = [.9 .9 0];
    fixStimulus.responseColor = stimulus.colors.white;
    fixStimulus.interColor = stimulus.colors.black;
    fixStimulus.correctColor = stimulus.colors.green;
    fixStimulus.incorrectColor = stimulus.colors.red;
    [task{2}, myscreen] = fixStairInitTask(myscreen);
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
    disp(sprintf('(cohCon) Initializing staircases'));
    stimulus = initStaircase(stimulus);
else
    disp('(cohCon) Re-using staircase from previous run...');
    % Reset staircase if necessary
    checkStaircaseStop();
end

%% EYE CALIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

%% Get Ready...
% clear screen    
mglWaitSecs(2);
setGammaTable_flowMax(1);
mglClearScreen(0.5);
if ~stimulus.unattended
    mglTextDraw(stimulus.runs.taskOptsText{stimulus.runs.curTask},[0 0]);
end
mglFlush

% let the user know
disp(sprintf('(cohCon) Starting run number: %i',stimulus.counter));
if ~stimulus.unattended
    myscreen.flushMode = 1;
end

%% Main Task Loop

phaseNum = 1;
% Again, only one phase.
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    % update the task
    [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);
    % update the fixation task
    if stimulus.unattended
        [task{2}, myscreen] = updateTask(task{2},myscreen,1);
    end
    % flip screen
    myscreen = tickScreen(myscreen,task);
end

stimulus.ended = mglGetSecs;

disp(sprintf('(cohCon) Run ending... Elapsed time: %3.0f s',stimulus.ended-stimulus.started));

% task ended
mglClearScreen(0.5);
mglTextDraw('Run complete... please wait.',[0 0]);
mglFlush
myscreen.flushMode = 1;

if stimulus.plots
    disp('(cohCon) Displaying plots');
    dispStaircase(stimulus);
    dispStaircaseCatch(stimulus);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

% if this was 'unattended' or 'mtloc' mode, we should copy the file we saved
% to a different folder.
dFolder = fullfile('~/data/coherentContrast/',mglGetSID);
files = dir(dFolder);
cFile = files(end);
if stimulus.mtloc 
    nFolder = fullfile('~/data/coherentContrast/',mglGetSID,'unattended');
    if ~isdir(nFolder), mkdir(nFolder); end
    s = movefile(fullfile(dFolder,cFile.name),fullfile(nFolder,cFile.name));
elseif stimulus.unattended
    nFolder = fullfile('~/data/coherentContrast/',mglGetSID,'unattended');
    if ~isdir(nFolder), mkdir(nFolder); end
    s = movefile(fullfile(dFolder,cFile.name),fullfile(nFolder,cFile.name));
else
    s = 1;
end

if ~s
    disp('File copy failed for some reason...');
    keyboard
end

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startTrialCallback(task,myscreen)

global stimulus

if stimulus.curTrial == 0
    stimulus.started = mglGetSecs;
    if stimulus.scan
      task.thistrial.seglen(stimulus.seg.ITI) = 7.5;
    end
end

if stimulus.mtloc
    % track the trial start time
    stimulus = mtRandDirs(stimulus);
end

stimulus.curTrial = stimulus.curTrial + 1;

%  Set the current task
if stimulus.unattended
    task.thistrial.task = 3;
else
    if task.thistrial.catch
        switchTasks = [2 1];
        task.thistrial.task = switchTasks(stimulus.runs.curTask);
        % edit seglen
        task.thistrial.seglen(stimulus.seg.ISI) = .5;
        task.thistrial.seglen(stimulus.seg.resp) = 2;
        disp('(cohCon) Catch trial.');
    else
        task.thistrial.task = stimulus.runs.curTask;
    end
end

% Set the missing thistrial vars
task.thistrial.coherence = stimulus.pedestals.coherence(task.thistrial.cohPedestal);
if stimulus.mtloc
    % fix the coherence so it switches
    task.thistrial.coherence = stimulus.pedestals.coherence(mod(stimulus.curTrial,2)+1);
end
task.thistrial.contrast = stimulus.pedestals.contrast(task.thistrial.conPedestal);
task.thistrial.trialNum = stimulus.curTrial;
if ~stimulus.unattended
    [task.thistrial.deltaPed, stimulus] = getDeltaPed(task,stimulus,task.thistrial.task);
else
    task.thistrial.deltaPed = 0;
end

% Assign the deltaPed to the correct locations
if task.thistrial.task==1
    % speed
    stimulus.live.cohDelta = task.thistrial.deltaPed;
    stimulus.live.conDelta = 0;
    if (task.thistrial.coherence + stimulus.live.cohDelta) > 0.95
        stimulus.live.cohDelta = .95 - task.thistrial.coherence;
    end
elseif task.thistrial.task==2
    % contrast
    stimulus.live.cohDelta = 0;
    stimulus.live.conDelta = task.thistrial.deltaPed;
    if (task.thistrial.contrast + stimulus.live.conDelta) > 1
        stimulus.live.conDelta = 1 - task.thistrial.contrast;
    end
elseif stimulus.unattended
    % unattended
    stimulus.live.cohDelta = 0;
    stimulus.live.conDelta = 0;
else
    warning('Never should get here... debug me');
    keyboard
end

ramps = task.thistrial.seglen(stimulus.seg.rampUP)*2;
main = task.thistrial.seglen(stimulus.seg.stim);
total = ramps+main;
if task.thistrial.side==1
    % left
    task.thistrial.avgCohL = task.thistrial.coherence + ramps/total*stimulus.sigmoidMu*stimulus.live.cohDelta + main/total*stimulus.live.cohDelta;
    task.thistrial.avgCohR = task.thistrial.coherence;
    task.thistrial.avgConL = task.thistrial.contrast + ramps/total*stimulus.sigmoidMu*stimulus.live.conDelta + main/total*stimulus.live.conDelta;
    task.thistrial.avgConR = task.thistrial.contrast;
    disp(sprintf('(cohCon) Trial %i starting. Coherence: L %.02f; R %.02f Contrast: L %.02f; R %.02f',task.thistrial.trialNum,...
        task.thistrial.coherence+stimulus.live.cohDelta,task.thistrial.coherence,...
        task.thistrial.contrast+stimulus.live.conDelta,task.thistrial.contrast));
else
    % right
    task.thistrial.avgCohR = task.thistrial.coherence + ramps/total*stimulus.sigmoidMu*stimulus.live.cohDelta + main/total*stimulus.live.cohDelta;
    task.thistrial.avgCohL = task.thistrial.coherence;
    task.thistrial.avgConR = task.thistrial.contrast + ramps/total*stimulus.sigmoidMu*stimulus.live.conDelta + main/total*stimulus.live.conDelta;
    task.thistrial.avgConL = task.thistrial.contrast;
    disp(sprintf('(cohCon) Trial %i starting. Coherence: L %.02f; R %.02f Contrast: L %.02f; R %.02f',task.thistrial.trialNum,...
        task.thistrial.coherence,task.thistrial.coherence+stimulus.live.cohDelta,...
        task.thistrial.contrast,task.thistrial.contrast+stimulus.live.conDelta));
end

% set the gammaTable for this trial
if ~stimulus.unattended
    setGammaTable_flowMax(task.thistrial.contrast + stimulus.live.conDelta);
else
    setGammaTable_flowMax(1);
end

% set directions
stimulus.dotsL.dir = task.thistrial.dir;
stimulus.dotsR.dir = task.thistrial.dir;

function ped = curPedValue(task)
if task.thistrial.task==1
    ped = task.thistrial.cohPedestal;
else
    ped = task.thistrial.conPedestal;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)

myscreen.flushMode = 0;

global stimulus

switch task.thistrial.thisseg
    case stimulus.seg.ITI
        stimulus.live.dots = 0;
        stimulus.live.fixColor = stimulus.colors.black;
        stimulus.live.catchFix = 0;
    case stimulus.seg.rampUP
        stimulus.live.dots = 1;
        stimulus.live.dotRampDir = 1;
        stimulus.live.fixColor = stimulus.colors.black;
        stimulus.live.catchFix = 0;
        stimulus.live.rampStart = mglGetSecs;
    case stimulus.seg.stim
        stimulus.live.dots = 1;
        stimulus.live.dotRampDir = 0;
        stimulus.live.fixColor = stimulus.colors.black;
        stimulus.live.catchFix = 0;
    case stimulus.seg.rampDOWN
        stimulus.live.dots = 1;
        stimulus.live.dotRampDir = -1;
        stimulus.live.fixColor = stimulus.colors.black;
        stimulus.live.catchFix = 0;
        stimulus.live.rampStart = mglGetSecs;
    case stimulus.seg.ISI
        stimulus.live.dots = 0;
        stimulus.live.fixColor = stimulus.colors.black;
        stimulus.live.catchFix = 1;
    case stimulus.seg.resp
        stimulus.live.dots = 0;
        stimulus.live.fixColor = stimulus.colors.white;
        stimulus.live.catchFix = 1;
end

if stimulus.unattended
    stimulus.live.dotRampDir = 0;
end

function value = calcPerc(stimulus,perc)
  value = stimulus.sigmoid(int8(round(perc,2)*100+1));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
global stimulus

if stimulus.mtloc && (mglGetSecs > stimulus.mt.nextFlip)
    stimulus = mtRandDirs(stimulus);
end

if stimulus.projector
    mglClearScreen(stimulus.colors.black);
    mglStencilSelect(stimulus.stencil);
    mglFillRect(0,0,[50 50],[.5 .5 .5]);
else
    mglClearScreen(0.5);
end

if stimulus.live.dots==1, stimulus = upDots(task,stimulus,myscreen); end
if ~stimulus.unattended, upFix(task,stimulus); end

if stimulus.projector, mglStencilSelect(0); end

%%
function upFix(task,stimulus)

if ~task.thistrial.catch || stimulus.live.catchFix == 0
    mglFixationCross(1.5,1.5,stimulus.live.fixColor);
else
    if task.thistrial.task==1
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

% update the dots
repick = logical(stimulus.mtloc);

tCoh = task.thistrial.coherence;
tCon = task.thistrial.contrast / stimulus.curMaxContrast;

switch stimulus.live.dotRampDir
    case 0
        perc = 1;
    case 1
        % we are ramping UP
        perc = (mglGetSecs-stimulus.live.rampStart) / task.thistrial.seglen(task.thistrial.thisseg);
    case -1
        % we are ramping DOWN        
        perc = 1-((mglGetSecs-stimulus.live.rampStart) / task.thistrial.seglen(task.thistrial.thisseg));
end

perc = calcPerc(stimulus,perc);

if task.thistrial.side==1
    lCohDel = perc*stimulus.live.cohDelta;
    rCohDel = 0;
    lConDel = perc*stimulus.live.conDelta;
    rConDel = 0;
else
    lCohDel = 0;
    rCohDel = perc*stimulus.live.cohDelta;
    lConDel = 0;
    rConDel = perc*stimulus.live.conDelta;
end

%% Old update code start here
stimulus.dotsL = updateDotsRadial(stimulus.dotsL,tCoh+lCohDel,myscreen,repick);
stimulus.dotsR = updateDotsRadial(stimulus.dotsR,tCoh+rCohDel,myscreen,repick);

% Correct values for gamma table adjustments
rConDel = rConDel / stimulus.curMaxContrast;
lConDel = lConDel / stimulus.curMaxContrast;

% dotsR
% update +contrast

mglPoints2(stimulus.dotsR.xdisp(stimulus.dotsR.con==1),stimulus.dotsR.ydisp(stimulus.dotsR.con==1),...
    stimulus.dotsR.dotsize,[.5 .5 .5] - adjustConToTable(tCon + rConDel,stimulus)/2);
% update - contrast
mglPoints2(stimulus.dotsR.xdisp(stimulus.dotsR.con==2),stimulus.dotsR.ydisp(stimulus.dotsR.con==2),...
    stimulus.dotsR.dotsize,[.5 .5 .5] + adjustConToTable(tCon + rConDel,stimulus)/2);
% dotsL
% update +contrast
mglPoints2(stimulus.dotsL.xdisp(stimulus.dotsL.con==1),stimulus.dotsL.ydisp(stimulus.dotsL.con==1),...
    stimulus.dotsL.dotsize,[.5 .5 .5] - adjustConToTable(tCon + lConDel,stimulus)/2);
% update - contrast
mglPoints2(stimulus.dotsL.xdisp(stimulus.dotsL.con==2),stimulus.dotsL.ydisp(stimulus.dotsL.con==2),...
    stimulus.dotsL.dotsize,[.5 .5 .5] + adjustConToTable(tCon + lConDel,stimulus)/2);

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
        task.thistrial.correct = task.thistrial.whichButton == stimulus.responseKeys(task.thistrial.side);
        % Store whether this was correct
        stimulus.live.fixColor = fixColors{task.thistrial.correct+1};
        disp(sprintf('(cohCon) Response %s: %s',responseText{task.thistrial.correct+1},responsePos{find(stimulus.responseKeys==task.thistrial.whichButton)}));
        if ~task.thistrial.catch
            stimulus.staircase{task.thistrial.task,curPedValue(task)} = ...
                doStaircase('update',stimulus.staircase{task.thistrial.task,curPedValue(task)},task.thistrial.correct);
        else
            stimulus.live.fixColor = stimulus.colors.black; % we never show information about catch trials
            stimulus.live.catchFix = 0;
            stimulus.stairCatch{task.thistrial.task} = ...
                doStaircase('update',stimulus.stairCatch{task.thistrial.task},task.thistrial.correct);
        end
    else
        disp(sprintf('(cohCon) Subject responded multiple times: %i',task.thistrial.gotResponse+1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              HELPER FUNCTIONS                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% getDeltaPed

function [deltaPed, stimulus] = getDeltaPed(task,stimulus,taskNum)
if task.thistrial.catch
    [deltaPed, stimulus.stairCatch{taskNum}] = doStaircase('testValue',stimulus.stairCatch{taskNum});
else
    [deltaPed, stimulus.staircase{taskNum,curPedValue(task)}] = doStaircase('testValue',stimulus.staircase{taskNum,curPedValue(task)});
end


%%%%%%%%%%%%%%%%%%%%%%%%
%    initStaircase     %
%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initStaircase(stimulus)

stimulus.stairCatch = cell(1,2);
stimulus.staircase = cell(2,length(stimulus.pedestals.contrast));

stimulus.stairCatch{1} = doStaircase('init','fixed',...
    'fixedVals',stimulus.pedestals.catch.coherence,'nTrials=25');
stimulus.stairCatch{2} = doStaircase('init','fixed',...
    'fixedVals',stimulus.pedestals.catch.contrast,'nTrials=25');
for task = 1:2
    stimulus.staircase{task,1} = doStaircase('init','upDown',...
        'initialThreshold',stimulus.pedestals.initThresh.(stimulus.pedestals.pedOpts{task}),...
        'initialStepsize',stimulus.pedestals.initThresh.(stimulus.pedestals.pedOpts{task})/3,...
        'minThreshold=0.001','maxThreshold=0.5','stepRule','pest', ...
        'nTrials=50','maxStepsize=.2','minStepsize=.001');
end

for i = 2:length(stimulus.pedestals.coherence)
    stimulus.staircase{1,i} = stimulus.staircase{1,1};
    stimulus.staircase{2,i} = stimulus.staircase{2,1};
end

function dispStaircaseCatch(stimulus)

try
    taskOpts = {'coherence','contrast'};
    drawing = {'--r' '-r'};
    
    figure
    hold on
    a1 = 1;
    a2 = 0;
    for task = 1:2
        pedSuccess = [];
        pedCount = [];
	pedPos = [];
	testV = [];
	resp = [];
	for i = 1:length(stimulus.stairCatch{task})
	  testV = [testV stimulus.stairCatch{task}(i).testValues];
	  resp = [resp stimulus.stairCatch{task}(i).response];
	end
        for i = 1:length(testV)
	  index = find(testV(i)==pedPos);
	  if isempty(index)
	    pedPos(end+1) = testV(i);
	    pedSuccess(end+1) = 0;
	    pedCount(end+1) = 0;
	    index = length(pedPos);
	  end
            pedSuccess(index) = pedSuccess(index) + resp(i);
            pedCount(index) = pedCount(index) + 1;
        end
        success = pedSuccess ./ pedCount;
	[pedPos is] = sort(pedPos);
	success = success(is);
        plot(pedPos,success,drawing{task});
	a1 = min(a1,min(stimulus.pedestals.catch.(taskOpts{task})));
	a2 = max(a2,max(stimulus.pedestals.catch.(taskOpts{task})));
        axis([a1 a2 -.05 1.05]);
    end
    legend(taskOpts);
catch
    disp('(cohCon) Catch figures did not display correctly...');
end

%%%%%%%%%%%%%%%%%%%%%%%
%    dispStaircase    %
%%%%%%%%%%%%%%%%%%%%%%%
function dispStaircase(stimulus)

try
    taskOpts = {'coherence','contrast'};
    
    plotting = zeros(2,4);
    
    drawing = {'-r' '-g' '-b' '-y'
                '--r' '--g' '--b' '--y'};
    for task = 1:2
        figure % this is the 'staircase' figure
        title(sprintf('%s, Staircase plot (R->G->B->Y high)',taskOpts{task}));
        hold on
        for ped = 1:4
            try
                testV = [];
                for i = 1:length(stimulus.staircase{task,ped})
                    testV = [testV stimulus.staircase{task,ped}(i).testValues];
                end
                plot(testV,drawing{task,ped});
            catch
            end
            try
                out = doStaircase('threshold',stimulus.staircase{task,ped},'type','weibull'); % noise, 1 cue, lowest
                plotting(task,ped) = out.threshold;
            catch
                plotting(task,ped) = -1;
            end
        end
    end
    hold off
    figure
    hold on
    title(sprintf('%s, R->G->B High',taskOpts{task}));
    plot(stimulus.pedestals.(taskOpts{task})(1:4),plotting(1,:),'-r');
    plot(stimulus.pedestals.(taskOpts{task})(1:4),plotting(2,:),'--r');
    legend(taskOpts);
    axis([stimulus.pedestals.(taskOpts{task})(1) stimulus.pedestals.(taskOpts{task})(4) 0 .5]);
    hold off

catch
    disp('(cohCon) Figures were not generated successfully.');
end

%% checkStaircaseStop
function checkStaircaseStop()
global stimulus
taskOpts = {'coherence','contrast'};
for task = 1:2
    s = stimulus.stairCatch{task};
    if doStaircase('stop',s)
        stimulus.stairCatch{task}(end+1) = doStaircase('init','fixed',...
            'fixedVals',stimulus.pedestals.catch.(taskOpts{task}),'nTrials=25');
    end
end
% Check both staircases
for task = 1:2
    for ped = 1:4
        s = stimulus.staircase{task,ped};
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
            stimulus.staircase{task,ped}(end+1) = doStaircase('init',s,'initialThreshold',vals{threshPos},'initialStepsize',vals{stepPos});
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create dots for optic flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initDotsRadial(dots,~)

global stimulus

% maximum depth of points
dots.minX = 3;
dots.maxX = 10;
dots.minY = -5;
dots.maxY = 5;

if stimulus.mtloc
    dots.minX = 2.5;
    dots.maxX = 11;
    dots.minY = -6;
    dots.maxY = 6;
end

dots.dir = 1;

% make a some points
dots.n = 500*dots.density;
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
freq_factor = dots.speed/myscreen.framesPerSecond;

% move coherent dots
dots.x(dots.coherent) = dots.x(dots.coherent) + dots.dir*freq_factor;

% these are for flipping into the other quadrants
xmat = [1 1 -1 -1];
ymat = [1 -1 1 -1];
xmat = repmat(xmat,1,dots.incoherentn+4-mod(dots.incoherentn,4));
ymat = repmat(ymat,1,dots.incoherentn+4-mod(dots.incoherentn,4));
xmat = xmat(randperm(dots.incoherentn));
ymat = ymat(randperm(dots.incoherentn));

% move incoherent dots
% get random vectors
dots.rX = rand(1,dots.incoherentn); % get a random number -1 -> 1 for the x offset
dots.rY = sqrt(1-dots.rX.^2); % solve for y
dots.rX = dots.rX .* xmat;
dots.rY = dots.rY .* ymat;
dots.rX = dots.rX * freq_factor; % rescale to match the velocity
dots.rY = dots.rY * freq_factor;
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

function stimulus = mtRandDirs(stimulus)

stimulus.mt.lastFlip = mglGetSecs;
stimulus.mt.nextFlip = stimulus.mt.lastFlip + 1.5;

dirOpts = [1 0 -1];
curDir = dirOpts(stimulus.dotsL.dir+2);
stimulus.dotsL.dir = curDir;
stimulus.dotsR.dir = curDir;

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