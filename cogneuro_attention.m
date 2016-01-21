
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
%   TR .5 = 296 volumes (10 * 14 * 2 + 16)

function [myscreen] = cogneuro_attention(varargin)

global stimulus
clear fixStimulus
%% Initialize Variables

disp('****************************************');
disp('** NEPR Into Cog Neuro Attention Task **');
disp('****************************************');
% add arguments later
stimFileNum = [];
plots = [];
overrideTask = [];
scan = [];
testing = [];
getArgs(varargin,{'stimFileNum=-1','testing=0' ...
    'plots=1','overrideTask=0','scan=0'});
stimulus.plots = plots;
stimulus.scan = scan;
stimulus.testing = testing;

if stimulus.scan && ~mglGetParam('ignoreInitialVols')==16 && ~mglGetParam('ignoreInitialVols')==4
    warning('ignoreInitialVols is set to %i.',mglGetParam('ignoreInitialVols'));
    if ~strcmp('y',input('Is this correct? [y/n]'))
        mglSetParam('ignoreInitialVols',input('Please input the correct value (mux8 = 16, mux2 = 4): '));
    end
end

stimulus.counter = 1; % This keeps track of what "run" we are on.

%% Useful stimulus stuff

stimulus.pedestals.coherence = .51;

stimulus.pedestals.initThresh.coherence = .4;

%% Setup Screen

if stimulus.scan
    myscreen = initScreen('fMRIprojFlex');
else
    myscreen = initScreen('VPixx');
end

%% Open Old Stimfile
stimulus.initStair = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/cogneuro_attention/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/cogneuro_attention/%s/1*mat',mglGetSID));
    
    if length(files) >= 1
        if stimFileNum == -1
            if length(files) > 1
                warning('Multiple stimfiles found, loading last one. If you wanted a different functionality use stimFileNum=#');
            end
            fname = files(end).name;
        else
            fname = files(stimFileNum).name;
        end
        s = load(sprintf('~/data/cogneuro_attention/%s/%s',mglGetSID,fname));
        stimulus.staircase = s.stimulus.staircase;
        stimulus.counter = s.stimulus.counter + 1;
        
        clear s;
        stimulus.initStair = 0;
        disp(sprintf('(cogneuro_att) Data file: %s loaded.',fname));
    end
end
disp(sprintf('(cogneuro_att) This is run #%i',stimulus.counter));

%% Initialize Stimulus

myscreen = initStimulus('stimulus',myscreen);

stimulus.responseKeys = [1 2]; % corresponds to LEFT - RIGHT


%% Everything else
stimulus.dots.xcenter = 0;
stimulus.dots.ycenter = 0;
stimulus.dots.dotsize = 4;
stimulus.dots.density = 1.4;
stimulus.dots.speed = 5;
stimulus.dots.centerOffset = 2;

stimulus.dotsR = stimulus.dots;
stimulus.dotsR.mult = 1;
stimulus.dotsL = stimulus.dots;
stimulus.dotsL.mult = -1;
stimulus = rmfield(stimulus,'dots');

stimulus.dotsR = initDotsRadial(stimulus.dotsR,myscreen);
stimulus.dotsL = initDotsRadial(stimulus.dotsL,myscreen);

%% Setup Task

% This is the contrast change detection task
task{1}{1}.waitForBacktick = 1;

stimulus.seg.cue = 1;
stimulus.seg.ITI = [4,7,10]; % the ITI is either 20s (first time) or 1s
stimulus.seg.stim = [2,5,8];
stimulus.seg.resp = [3, 6, 9];

l_b = 2;
l_s = .5;
l_r = 1;
l_i = 0;
task{1}{1}.seglen = [l_b l_s l_r l_i l_s l_r l_i l_s l_r l_i];

task{1}{1}.synchToVol = [0 0 0 0 0 0 0 0 0 0];
if stimulus.scan
    task{1}{1}.synchToVol(end) = 1;
end

task{1}{1}.getResponse = [0 0 1 0 0 1 0 0 1 0];
task{1}{1}.parameter.task = [1 2];
task{1}{1}.parameter.cohPedestal = 1; % target flow coherence
task{1}{1}.random = 0;
task{1}{1}.numTrials = 10;

%% Tracking

% these are variables that we want to track for later analysis.
task{1}{1}.randVars.calculated.correct = nan;
task{1}{1}.randVars.calculated.trialNum = nan;
task{1}{1}.randVars.calculated.lCoh = nan;
task{1}{1}.randVars.calculated.rCoh = nan;

stimulus.curTrial = 0;

%% Block setup
if isfield(stimulus,'runs') && isfield(stimulus.runs,'loaded')
    % We already have our blocks
    stimulus.runs = rmfield(stimulus.runs,'loaded'); % remove the load field, otherwise it gets saved across runs
else
    % This is the first run, build up the blocks.
    stimulus.runs = struct;
    stimulus.runs.taskOpts = [1 2];
    stimulus.runs.taskOptsText = {'Attend Left','Attend Right'};
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
    disp(sprintf('(cogneuro_att) Initializing staircases'));
    stimulus = initStaircase(stimulus);
else
    disp('(cogneuro_att) Re-using staircase from previous run...');
    % Reset staircase if necessary
    checkStaircaseStop();
end

%% Get Ready...
% clear screen
mglWaitSecs(2);
mglClearScreen(0.5);
mglTextDraw('DO NOT MOVE',[0 1.5]);
mglFlush

% let the user know
disp(sprintf('(cogneuro_att) Starting run number: %i',stimulus.counter));
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

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

if stimulus.plots
    disp('(cohCon) Displaying plots');
    dispStaircase(stimulus);
end

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startTrialCallback(task,myscreen)
%%

global stimulus

stimulus.curTrial = stimulus.curTrial + 1;
task.thistrial.trialNum = stimulus.curTrial;
disp(sprintf('(cogneuro_att) Block %i starting.',stimulus.curTrial));

stimulus.live.group = 0;
stimulus.live.upgroup = 0;

stimulus.live.dirsL = randi(2,1,3);
stimulus.live.dirsR = randi(2,1,3);
disp(sprintf('(cogneuro_att) Current task: %s',stimulus.runs.taskOptsText{task.thistrial.task}));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lCoh, rCoh, stimulus] = getCurrentThresholds(task,stimulus)
% Set the missing thistrial vars
task.thistrial.coherence = stimulus.pedestals.coherence(task.thistrial.cohPedestal);

[lCoh, stimulus] = getDeltaPed(task,1,stimulus);
[rCoh, stimulus] = getDeltaPed(task,2,stimulus);

if (task.thistrial.coherence + lCoh) > 0.95
    lCoh = 0.95 - task.thistrial.coherence;
end
if (task.thistrial.coherence + rCoh) > 0.95
    rCoh = 0.95 - task.thistrial.coherence;
end

lCoh = task.thistrial.coherence + lCoh;
rCoh = task.thistrial.coherence + rCoh;
%%%%%%%%%%%%%%%%%%%%%%%%
%    getDeltaPed       %
%%%%%%%%%%%%%%%%%%%%%%%%

function [cohPed, stimulus] = getDeltaPed(task,ct, stimulus)
[cohPed, stimulus.staircase{ct}] = doStaircase('testValue',stimulus.staircase{ct});

function [task, myscreen] = startSegmentCallback(task, myscreen)
%%
myscreen.flushMode = 0;
global stimulus

flipR = [-1 1];
flipL = [1 -1];
dirT = {'<-','->'};
stimulus.live.cue = 0;
stimulus.live.dots = 0;
stimulus.live.fixColor = 0;
if any(stimulus.seg.cue==task.thistrial.thisseg)
    stimulus.live.cue = 1;
end
if any(stimulus.seg.stim==task.thistrial.thisseg)
    if ~stimulus.live.upgroup
        stimulus.live.group = stimulus.live.group + 1;
        stimulus.live.upgroup = 1;
        stimulus.live.cdirL = stimulus.live.dirsL(stimulus.live.group);
        stimulus.live.cdirR = stimulus.live.dirsR(stimulus.live.group);
        stimulus.dotsL.dir = flipR(stimulus.live.cdirL);
        stimulus.dotsR.dir = flipL(stimulus.live.cdirR);
        [lc, rc, s] = getCurrentThresholds(task,stimulus);
        stimulus = s;
        stimulus.live.lCoh = lc;
        stimulus.live.rCoh = rc;
        disp(sprintf('(cogneuro_att) Group: %i',stimulus.live.group));
        disp(sprintf('(cogneuro_att) Left dir: %s, Right dir: %s',dirT{stimulus.live.cdirL},dirT{stimulus.live.cdirR}));
        disp(sprintf('(cogneuro_att) Left coherence: %0.2f, Right coherence: %0.2f',stimulus.live.lCoh,stimulus.live.rCoh));
    end
    stimulus.live.dots = 1;
end
if any(stimulus.seg.resp==task.thistrial.thisseg)
    stimulus.live.fixColor = 1;
    stimulus.live.upgroup = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus
mglClearScreen(0.5);
if stimulus.live.dots==1, stimulus = upDots(task,stimulus,myscreen); end
if stimulus.live.cue, upCue(task,stimulus); end
upFix(stimulus);

function upCue(task,stimulus)
%%
mglTextDraw(stimulus.runs.taskOptsText{task.thistrial.task},[0 1.5]);

function upFix(stimulus)
%%
mglFixationCross(1.5,1.5,stimulus.live.fixColor);


function stimulus = upDots(task,stimulus,myscreen)

% update the dots

stimulus.dotsL = updateDotsRadial(stimulus.dotsL,stimulus.live.lCoh,myscreen,false);
stimulus.dotsR = updateDotsRadial(stimulus.dotsR,stimulus.live.rCoh,myscreen,false);

% dotsR
% update +contrast

mglPoints2(stimulus.dotsR.xdisp(stimulus.dotsR.con==1),stimulus.dotsR.ydisp(stimulus.dotsR.con==1),...
    stimulus.dotsR.dotsize,[1 1 1]);
% update - contrast
mglPoints2(stimulus.dotsR.xdisp(stimulus.dotsR.con==2),stimulus.dotsR.ydisp(stimulus.dotsR.con==2),...
    stimulus.dotsR.dotsize,[0 0 0]);
% dotsL
% update +contrast
mglPoints2(stimulus.dotsL.xdisp(stimulus.dotsL.con==1),stimulus.dotsL.ydisp(stimulus.dotsL.con==1),...
    stimulus.dotsL.dotsize,[1 1 1]);
% update - contrast
mglPoints2(stimulus.dotsL.xdisp(stimulus.dotsL.con==2),stimulus.dotsL.ydisp(stimulus.dotsL.con==2),...
    stimulus.dotsL.dotsize,[0 0 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus

responseText = {'Incorrect','Correct'};
responsePos = {'Left','Right'};
fixColors = {[1 0 0],[0 1 0]};

if any(task.thistrial.whichButton == stimulus.responseKeys)
    if task.thistrial.gotResponse < 3
        if task.thistrial.task==1
            cSide = stimulus.live.cdirL;
        else
            cSide = stimulus.live.cdirR;
        end
        corr = task.thistrial.whichButton == stimulus.responseKeys(cSide);
        if isnan(task.thistrial.correct), task.thistrial.correct=0; end
        task.thistrial.correct = task.thistrial.correct + 10^(stimulus.live.group-1)*corr;
        % Store whether this was correct
        stimulus.staircase{task.thistrial.task} = doStaircase('update',stimulus.staircase{task.thistrial.task},corr);
        stimulus.live.fixColor = fixColors{corr+1};
        disp(sprintf('(cogneuro_att) Response %s: %s',responseText{corr+1},responsePos{find(stimulus.responseKeys==task.thistrial.whichButton)}));
    else
        disp(sprintf('(cogneuro_att) Subject responded multiple times: %i',task.thistrial.gotResponse+1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              HELPER FUNCTIONS                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%
%    initStaircase     %
%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initStaircase(stimulus)
%%
stimulus.staircase = cell(2,1);

% Main staircases
for task = 1:2
    stimulus.staircase{task} = doStaircase('init','upDown',...
        'initialThreshold',stimulus.pedestals.initThresh.coherence,...
        'initialStepsize',stimulus.pedestals.initThresh.coherence/3,...
        'minThreshold=0.001','maxThreshold=0.4','stepRule','pest', ...
        'nTrials=50','maxStepsize=.2','minStepsize=.001');
end
disp('(cogneuro_att) New staircase: nTrials=50, max=.4');
%%
%%%%%%%%%%%%%%%%%%%%%%%
%    dispStaircase    %
%%%%%%%%%%%%%%%%%%%%%%%
function dispStaircase(stimulus)
disp('(cogneuro_att) Todo: Implement');

%% checkStaircaseStop
function checkStaircaseStop()
global stimulus

% Check both staircases
for task = 1:2
    stimulus.staircase{task,1} = resetStair(stimulus.staircase{task,1});
end

function s = resetStair(s)

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
dots.maxX = 11;
dots.minY = -5;
dots.maxY = 5;

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
dots.coherency = 0.5;
dots.moveright = rand(1,dots.n) < dots.coherency;
dots.moverightn = sum(dots.moveright);
dots.moveleft = ~dots.moveright;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots for Radial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updateDotsRadial(dots,coherence,myscreen,repick)

% stuff to compute median speed
dots.oldx = dots.x;
dots.oldy = dots.y;

% get the coherent and incoherent dots
if repick
    dots.moveright = rand(1,dots.n) < dots.coherency;
    dots.moverightn = sum(dots.moveright);
    dots.moveleft = ~dots.moveright;
    dots.coherency = 0.5;
elseif dots.coherency ~= coherence
    cohDiff = coherence - dots.coherency;
    numDots = round(abs(cohDiff) * dots.n); % actual number of dots to flip
    if numDots > dots.n, numDots = dots.n; end
    if cohDiff > 0
        % we need to add more coherent dots
        flipDots = [zeros(1,numDots) ones(1,sum(dots.moveright)-numDots)];
        dots.moveright(dots.moveright) = flipDots(randperm(length(flipDots)));
    else
        % we need to add more incoherent dots
        flipDots = [ones(1,numDots) zeros(1,sum(dots.moveleft)-numDots)];
        dots.moveright(dots.moveleft) = flipDots(randperm(length(flipDots)));
    end
    dots.moverightn = sum(dots.moveright);
    dots.moveleft = ~dots.moveright;
    dots.coherency = sum(dots.moveleft)/dots.n;
end

freq_factor = dots.speed/myscreen.framesPerSecond;

% move coherent dots
dots.x(dots.moveright) = dots.x(dots.moveright) + dots.dir*freq_factor;
dots.x(dots.moveleft) = dots.x(dots.moveleft) - dots.dir*freq_factor;

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