
% cohCon
%
%      usage: myscreen=cohcon()
%         by: daniel birman
%       date: 11/10/14
%    purpose: Coherence test for cohcon experiment.
%
%        use: cohcon(flags)
%
%       desc: runs the cohcon experiment in coherence test mode. This has a
%       few options, either you can run it with the constant background (0%
%       coherence), or you can have it run with trials. There are also
%       adjustments that can be made for contrast (default 10%) and
%       coherence values (default 0/25/50/75/100). If task is enabled you
%       will do the regular coherence discrimination task from cohcon, if
%       task is disabled a fixation staircase task will run. Dots are
%       white/black on a grey background to maintain constant luminance
%       while contrast is adjustable.
%
%      flags:
%
%
%
%   TR .75 = 560 volumes (7:00 total)
%   TR 1.4 = 300 volumes (7:00 total)

function [myscreen] = cohcon_localizer(varargin)

global stimulus
global fixStimulus
%% Initialize Variables

% add arguments later
stimFileNum = [];
scan = [];
stablecon = [];
task = [];
constant = [];
test = 0;
timing = 0;
getArgs(varargin,{'stimFileNum=-1','scan=1', ...
    'stablecon=0','task=2','constant=1','test=0',...
    'timing=0'});
stimulus.scan = scan;
stimulus.stablecon = stablecon;
stimulus.task = task; clear task
stimulus.constant = constant;
stimulus.timing = timing;

stimulus.counter = 1; % This keeps track of what "run" we are on.

%% Setup Screen

if stimulus.scan
    myscreen = initScreen('fMRIprojFlex');
else
    myscreen = initScreen('VPixx');
end

%% Setup task
if stimulus.task==2
    fixStimulus.diskSize = 0;
    fixStimulus.stimColor = [0 .6 .6];
    fixStimulus.responseColor = [.6 .6 0];
    fixStimulus.interColor = [0 .6 .6];
    fixStimulus.correctColor = [0 0.6 0];
    fixStimulus.incorrectColor = [0.6 0 0];
    [task{2}, myscreen] = fixStairInitTask(myscreen);
end

%% Open Old Stimfile
stimulus.initStair = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/cohcon_localizer/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/cohcon_localizer/%s/1*mat',mglGetSID));
    
    if length(files) >= 1
        if stimFileNum == -1
            if length(files) > 1
                warning('Multiple stimfiles found, loading last one. If you wanted a different functionality use stimFileNum=#');
            end
            fname = files(end).name;
        else
            fname = files(stimFileNum).name;
        end
        s = load(sprintf('~/data/cohcon_localizer/%s/%s',mglGetSID,fname));
        stimulus.counter = s.stimulus.counter + 1;
        
        % load blocks too
        stimulus.runs = s.stimulus.runs;
        stimulus.runs.loaded = 1;
        
        clear s;
        stimulus.initStair = 0;
        disp(sprintf('(cohCon) Data file: %s loaded.',fname));
    end
end
disp(sprintf('(cohCon) This is run #%i',stimulus.counter));

%% Initialize Stimulus

myscreen = initStimulus('stimulus',myscreen);

if stimulus.scan
    stimulus.responseKeys = [1 2]; % corresponds to LEFT - RIGHT
else
    stimulus.responseKeys = [1 2]; % corresponds to LEFT - RIGHT
end

%% Colors
stimulus.colors.rmed = 127.5;

% We're going to add an equal number of reserved colors to the top and
% bottom, to try to keep the center of the gamma table stable.
% % stimulus.colors.reservedBottom = [0 0 0; .6 .6 .6]; % fixation cross colors
% % stimulus.colors.reservedTop = [.5 0 0; 0 .5 0]; % correct/incorrect colors
% % stimulus.colors.black = 0/255; stimulus.colors.white = 1/255;
% % stimulus.colors.red = 254/255; stimulus.colors.green = 255/255;
% % stimulus.colors.nReserved = 2; % this is /2 the true number, because it's duplicated
% % stimulus.colors.nUnreserved = 256-(2*stimulus.colors.nReserved);
% %
% % stimulus.colors.mrmax = stimulus.colors.nReserved - 1 + stimulus.colors.nUnreserved;
% % stimulus.colors.mrmin = stimulus.colors.nReserved;
stimulus.colors.black = [0 0 0];
stimulus.colors.white = [.6 .6 .6];
stimulus.colors.red = [.5 0 0];
stimulus.colors.green = [0 .5 0];

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

%% Gamma Table Initialization

% get gamma table
% % if ~isfield(myscreen,'gammaTable')
% %   stimulus.linearizedGammaTable = mglGetGammaTable;
% %   disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
% %   disp(sprintf('(cuecon:initGratings) No gamma table found in myscreen. Contrast displays like this'));
% %   disp(sprintf('         should be run with a valid calibration made by moncalib for this monitor.'));
% %   disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
% % end
% % stimulus.linearizedGammaTable = myscreen.initScreenGammaTable;

%% Setup Task

% This is the contrast change detection task
task{1}{1}.waitForBacktick = 1;

stimulus.seg.ITI = 5; % the ITI is either 20s (first time) or 1s
stimulus.seg.stim = 1;
stimulus.seg.mask = 2;
stimulus.seg.ISI = 3;
stimulus.seg.resp = 4;
task{1}{1}.segmin = [2.5 0 .2 1 .2];
task{1}{1}.segmax = [2.5 0 1 1 .4];

if stimulus.scan
    task{1}{1}.segmin(stimulus.seg.ITI) = 2;
    if test
        task{1}{1}.segmax(stimulus.seg.ITI) = 3;
    else
        task{1}{1}.segmax(stimulus.seg.ITI) = 11;
    end
    task{1}{1}.segmin(stimulus.seg.ISI) = .2;
    task{1}{1}.segmax(stimulus.seg.ISI) = 1;
end

if stimulus.task==2
    % remove the response and ISI segments
    task{1}{1}.segmin(stimulus.seg.ISI) = 0;
    task{1}{1}.segmax(stimulus.seg.ISI) = 0;
    task{1}{1}.segmin(stimulus.seg.resp) = 0;
    task{1}{1}.segmax(stimulus.seg.resp) = 0;
    % everything is now zero except for the ITI of 2-11
end

task{1}{1}.synchToVol = [0 0 0 0 0];
if stimulus.scan
    if test
        task{1}{1}.synchToVol(stimulus.seg.ITI) = 0;
    else
        task{1}{1}.synchToVol(stimulus.seg.ITI) = 1;
    end
end
task{1}{1}.getResponse = [0 0 0 1 0];
task{1}{1}.parameter.cohSide = [1 2];
task{1}{1}.parameter.dir = [-1 1];
task{1}{1}.parameter.contrast = [0.25 0.5 0.75 1]; % contrast starts at 25%
task{1}{1}.parameter.coherence = [0 0.25 0.5 0.75 1]; % coherence starts at 0%
task{1}{1}.random = 1;
task{1}{1}.numTrials = 100;

stimulus.baseCoh = 0;
stimulus.baseCon = 0.25;

if stimulus.scan
    task{1}{1}.numTrials = inf;
    task{1}{2} = task{1}{1};
    task{1}{2}.waitForBacktick = 0;
    task{1}{1}.parameter.contrast = stimulus.baseCon;
    task{1}{1}.parameter.coherence = stimulus.baseCoh;
    task{1}{1}.numTrials = 1;
    if test
        task{1}{1}.segmin = [0 0 0 0 3];
        task{1}{1}.segmax = [0 0 0 0 3];
    else
        task{1}{1}.segmin = [0 0 0 0 30];
        task{1}{1}.segmax = [0 0 0 0 30];
    end
    if stimulus.timing
        disp('(cohcon_localizer) Freezing contrast, coherence 25/100%, timing .25 .5 1 2 4');
        task{1}{2}.parameter.timing = [0.250 0.500 1.00 2.00 4.00];
        task{1}{2}.parameter.contrast = stimulus.baseCon;
        task{1}{2}.parameter.coherence = [0.25 1];
    end
end

%% Tracking

% these are variables that we want to track for later analysis.
task{1}{1}.randVars.calculated.trialNum = nan;

stimulus.curTrial = 0;

%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{1})
    if stimulus.task==1
        [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);
    else
        [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,[],@startTrialCallback,[],[]);
    end 
end

%% EYE CALIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~stimulus.scan
    myscreen = eyeCalibDisp(myscreen);
end

%% Get Ready...
% clear screen
mglWaitSecs(2);
% setGammaTable_flowMax(1);
mglClearScreen(0.5);
if stimulus.scan
    mglTextDraw('DO NOT MOVE',[0 1.5]);
end
if stimulus.task==1
    mglTextDraw('Motion',[0 0]);
else
    %mglTextDraw('Fixation',[0 0]);
end
mglFlush

stimulus.runs.taskOptsText = {'Motion','Fixation'};
% let the user know
disp(sprintf('(cohCon) Starting run number: %i. Current task: %s',stimulus.counter,stimulus.runs.taskOptsText{stimulus.task}));
% if stimulus.unattended
myscreen.flushMode = 1;
% end

%% Main Task Loop

phaseNum = 1;
% Again, only one phase.
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    % update the task
    [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);
    % update fixation
    [task{2}, myscreen] = updateTask(task{2},myscreen,1);
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

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startTrialCallback(task,myscreen)
%%

global stimulus

if task.thistrial.thisphase==2
    task.thistrial.seglen(end) = 1.05^(rand*30+20);

    if stimulus.timing==1
        task.thistrial.seglen(stimulus.seg.stim) = task.thistrial.timing;
    end
end

%%

stimulus.curTrial = stimulus.curTrial + 1;
task.thistrial.trialNum = stimulus.curTrial;

disp(sprintf('(cohCon) Trial %i starting. Coherence: %.02f Contrast: %.02f. Length %1.1f s ITI %1.1f s',task.thistrial.trialNum,...
    task.thistrial.coherence,task.thistrial.contrast,task.thistrial.seglen(stimulus.seg.stim), task.thistrial.seglen(end)));

% set directions
stimulus.dotsL.dir = task.thistrial.dir;
stimulus.dotsR.dir = task.thistrial.dir;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)
%%

myscreen.flushMode = 0;

global stimulus

stimulus.live.mt = 0;
switch task.thistrial.thisseg
    case stimulus.seg.ITI
        stimulus.live.dots = 0;
        stimulus.live.fixColor = stimulus.colors.black;
        stimulus.live.catchFix = 0;
        stimulus.live.mask = 0;
    case stimulus.seg.stim
        stimulus.live.dots = 1;
        stimulus.live.fixColor = stimulus.colors.black;
        stimulus.live.catchFix = 0;
        stimulus.live.mask = 0;
    case stimulus.seg.mask
        stimulus.live.dots = 0;
        stimulus.live.fixColor = stimulus.colors.black;
        stimulus.live.catchFix = 0;
        stimulus.live.mask = 1;
        stimulus.live.maskOn = mglGetSecs;
    case stimulus.seg.ISI
        stimulus.live.dots = 0;
        stimulus.live.fixColor = stimulus.colors.black;
        stimulus.live.catchFix = 1;
        stimulus.live.mask = 0;
    case stimulus.seg.resp
        stimulus.live.dots = 0;
        stimulus.live.fixColor = stimulus.colors.white;
        stimulus.live.catchFix = 1;
        stimulus.live.mask = 0;
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
if stimulus.task==1, upFix(task,stimulus); end


function upFix(task,stimulus)
%%
mglFixationCross(1.5,1.5,stimulus.live.fixColor);


function stimulus = upDots(task,stimulus,myscreen)

% update the dots

if task.thistrial.thisseg==stimulus.seg.stim
    coh = task.thistrial.coherence;
    con = task.thistrial.contrast;
else
    coh = stimulus.baseCoh;
    con = stimulus.baseCon;
end

%% Old update code start here
stimulus.dotsL = updateDotsRadial(stimulus.dotsL,coh,myscreen,true);
stimulus.dotsR = updateDotsRadial(stimulus.dotsR,coh,myscreen,true);

mglPoints2(stimulus.dotsR.xdisp(stimulus.dotsR.con==1),stimulus.dotsR.ydisp(stimulus.dotsR.con==1),...
    stimulus.dotsR.dotsize,[.5 .5 .5] - con/2);
% update - contrast
mglPoints2(stimulus.dotsR.xdisp(stimulus.dotsR.con==2),stimulus.dotsR.ydisp(stimulus.dotsR.con==2),...
    stimulus.dotsR.dotsize,[.5 .5 .5] + con/2);
% dotsL
% update +contrast
mglPoints2(stimulus.dotsL.xdisp(stimulus.dotsL.con==1),stimulus.dotsL.ydisp(stimulus.dotsL.con==1),...
    stimulus.dotsL.dotsize,[.5 .5 .5] - con/2);
% update - contrast
mglPoints2(stimulus.dotsL.xdisp(stimulus.dotsL.con==2),stimulus.dotsL.ydisp(stimulus.dotsL.con==2),...
    stimulus.dotsL.dotsize,[.5 .5 .5] + con/2);

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
        cSide = task.thistrial.cohSide;
        
        task.thistrial.correct = task.thistrial.whichButton == stimulus.responseKeys(cSide);
        % Store whether this was correct
        stimulus.live.fixColor = fixColors{task.thistrial.correct+1};
        disp(sprintf('(cohCon) Response %s: %s',responseText{task.thistrial.correct+1},responsePos{find(stimulus.responseKeys==task.thistrial.whichButton)}));
        stimulus.staircase = doStaircase('update',stimulus.staircase,task.thistrial.correct);
    else
        disp(sprintf('(cohCon) Subject responded multiple times: %i',task.thistrial.gotResponse+1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              HELPER FUNCTIONS                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
