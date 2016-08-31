function [ myscreen ] = berlin_experiment( varargin )
%CAT_AWE Testing category discrimination at low visibility conditions
%   This is dan's experiment for Berlin (fall 2016). Version 3! Now this is
%   a motion experiment. On each trial (0.3 s stim + 8 s ISI + 2s response
%   + 8 s ITI) there is a brief blast of motion/shape/category which is
%   then followed by a constant motion stimulus during the response period.
%   The task is *always* match/non-match, but can be varied to be either
%   match-to-sample (direction/shape) or match-to-category (up-left vs.
%   down-right). The underlying patch is moving dots. Shapes are defined by
%   luminance. Categories are defined by an arbitrary SW->NE dividing line.
%   There are six directions used at equal pi/3 spacing. By default the
%   task is MS on motion, obviously the category information is technically
%   present. You have to enable the shape stimulus directly.
%
%   The default task is invisible.
%   Staircasing mode runs a staircase over mask strengths to d'=1
%   performance.
%   Localizer mode runs the task with scanner timing (8s delays) vs.
%   psychophysics timing (0.5-1s delays).
%
%   OPTIONS
%
%   shape=1 enable shape stimulus (changes the frame from circular to
%   either circular or cross)
%   task=1,2,3 (default=1) DMS direction, DMS shape, DMC category
%   staircase=1 runs in staircase mode (variable mask strength)
%   localizer=1 runs with scanner timing and fixed mask strengths
%
%   noeye=0 (default=1) skips the eye movement detection script
%   plots=0,1,2 show plots, or show plots and then return
%
%   NOTES
%
%   The code has independent staircases for each mode. In invisible mode it
%   runs fixed value staircases for the different tasks, depending on
%   whether shapes are present. No staircases in localizer mode. In
%   staircase mode it runs a staircase specifically to determine necessary
%   mask luminance to block vision.


global stimulus

%% OVERRIDES (for testing)

% run: `

% set all to -1 when running:
stimulus.contrastOverride = 0/255;
stimulus.lowOverride = 5/255;
stimulus.timingOverride = -1;

% actual mask contrasts for scanning
stimulus.maskContrasts = [0 32];

%% Initialize Variables

% add arguments later
localizer = 0;
staircase = 0;
scan = 0;
plots = 0; task = 0;
noeye = 0; shape = 0;
getArgs(varargin,{'localizer=0','staircase=0','scan=0','plots=0','category=0','noeye=1','constant=1','shape=0','task=1'});
stimulus.shape = shape;
stimulus.task = task;
stimulus.localizer = localizer;
stimulus.scan = scan;
stimulus.staircasing = staircase;
stimulus.plots = plots;
stimulus.category = category; % use match to category rules instead of match to sample rules
stimulus.noeye = noeye;
stimulus.constant = constant;
clear localizer invisible scan category noeye task

if stimulus.task==2 && stimulus.shape==0
    disp('(berlin) Auto-setting shape for shape task run');
    stimulus.shape = 1;
end
if stimulus.category==1
    disp('(berlin) Not setup for categories yet');
    return
end
if stimulus.localizer==1 && stimulus.shape==0
    disp('(berlin) Auto-setting shape for localizer task');
    stimulus.shape = 1;
end
if stimulus.staircasing && stimulus.localizer
    disp('(berlin) Cannot run invisible and localizer simultaneously');
    return
end

stimulus.counter = 1; % This keeps track of what "run" we are on.

%% Setup Screen

if stimulus.scan
    myscreen = initScreen('fMRIprojFlex');
else
    myscreen = initScreen('VPixx');
end

myscreen.background = 0;

% stimulus.linearizedGammaTable = mglGetGammaTable;
% setGT(myscreen,stimulus);
% mglWaitSecs(1);

%% Open Old Stimfile
stimulus.initStair = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/berlin_experiment/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/berlin_experiment/%s/1*mat',mglGetSID));

    if length(files) >= 1
        fname = files(end).name;
        
        s = load(sprintf('~/data/berlin_experiment/%s/%s',mglGetSID,fname));
        % copy staircases and run numbers
        stimulus.staircase = s.stimulus.staircase;
        stimulus.istaircase = s.stimulus.istaircase;
        stimulus.run = s.stimulus.run;
        stimulus.counter = s.stimulus.counter + 1;
        stimulus.run = s.stimulus.run;

        clear s;
        stimulus.initStair = 0;
        disp(sprintf('(berlin) Data file: %s loaded.',fname));
        
    end
end
disp(sprintf('(berlin) This is run #%i',stimulus.counter));


if stimulus.plots==2
    dispInfo(stimulus);
    return
end

stimulus.run.points = 0; % reset points for this run

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% init staircase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if stimulus.initStair
    % We are starting our staircases from scratch
    disp(sprintf('(berlin) Initializing staircase'));
    stimulus = initStaircase(stimulus);
    
    stimulus = initRuns(stimulus);
else
    disp('(berlin) Checking staircase resets from previous run...');
    stimulus = resetStair(stimulus);
end

%% Initialize Stimulus

myscreen = initStimulus('stimulus',myscreen);

if stimulus.scan
    stimulus.responseKeys = [2 1]; % corresponds to NOMATCH, MATCH
else
    stimulus.responseKeys = [2 1]; % corresponds to  NOMATCH, MATCH
end

stimulus.colors.black = [0 0 0];
stimulus.colors.white = [0.2 0.2 0.2];
stimulus.colors.green = [0 0.2 0];
stimulus.colors.red = [0.2 0 0];
% stimulus.colors.black = 0;
% stimulus.colors.white = 1/255;
% stimulus.colors.red = 2/255;
% stimulus.colors.green = 3/255;

stimulus.ring.inner = 3.5;
stimulus.ring.outer = 7; %
stimulus.area = 3.14159265358979*((stimulus.ring.outer/2)^2-(stimulus.ring.inner/2)^2);

if stimulus.staircasing
    stimulus.counter = stimulus.counter+1;
    if stimulus.counter>length(stimulus.run.stimCon)
        stimulus.run.stimCon = repmat(stimulus.run.stimCon,1,2);
        stimulus.run.stimLengths = repmat(stimulus.run.stimLengths,1,2);
    end
    stimulus.lowCon = stimulus.run.stimCon(stimulus.counter)/255;
else
    stimulus.lowCon = 1/255; % minimum possible contrast
end
if stimulus.lowOverride>=0
    stimulus.lowCon = stimulus.lowOverride;
end

%% Generate stencils

% make the outer stencil (big, circular)
mglStencilCreateBegin(999);
mglFillOval(0,0,repmat(stimulus.ring.outer+3,1,2),1);
mglStencilCreateEnd;
mglClearScreen(0);
% make the circular stencil
mglStencilCreateBegin(99);
mglFillOval(0,0,repmat(stimulus.ring.outer,1,2),1);
mglStencilCreateEnd;
mglClearScreen(0);
% make the cross stencil, we try to approximate the area to be ~216 degrees
% squared, so it has the same total size as the circle. This means each of
% the five "parts" needs to be 6.6 degrees wide and high
mult = sqrt((stimulus.area+pi*(stimulus.ring.inner/2)^2)/5);
mglStencilCreateBegin(100);
mglFillRect([0 0 0 -mult mult],[0 -mult mult 0 0],[mult mult],[1 1 1]);
mglStencilCreateEnd;
mglClearScreen(0);


%% Setup Task
task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;

stimulus.curTrial = 0;

% calculate timing
stimulus.seg.mask1 = [];
stimulus.seg.mask2 = [];
stimulus.seg.stim = [];
task{1}{1}.seglen = [];

stimulus.stimRepeats = 8;
for i = 1:stimulus.stimRepeats
    stimulus.seg.mask1(end+1) = (i-1)*3+1;
    stimulus.seg.stim(end+1) = (i-1)*3+2;
    stimulus.seg.mask2(end+1) = (i-1)*3+3;
    task{1}{1}.seglen(end+1:end+3) = [0.05 0.1 0.05];
end
task{1}{1}.seglen(1) = 0.1;
task{1}{1}.seglen(end) = 0.1;

stimulus.seg.delay = (i-1)*3+4;
stimulus.seg.resp = (i-1)*3+5;
stimulus.seg.ITI = (i-1)*3+6;
task{1}{1}.seglen(end+1:end+3) = [0.5 2 1];

task{1}{1}.synchToVol = zeros(size(task{1}{1}.seglen));
task{1}{1}.getResponse = zeros(size(task{1}{1}.seglen)); task{1}{1}.getResponse(stimulus.seg.resp)=1;
task{1}{1}.numTrials = 55;
task{1}{1}.random = 1;
task{1}{1}.parameter.task = stimulus.task;

if stimulus.localizer
    % longer delay, necessary for scanning
    if stimulus.timingOverride<0
        task{1}{1}.seglen(stimulus.seg.delay) = 8;
        task{1}{1}.seglen(stimulus.seg.ITI) = 8;
    end 
    task{1}{1}.synchToVol(stimulus.seg.ITI) = 1;
    task{1}{1}.numTrials = Inf;
    task{1}{1}.parameter.contrast = stimulus.maskContrasts/255;
end
if stimulus.staircasing
    task{1}{1}.seglen(stimulus.seg.stim) = stimulus.run.stimLengths(stimulus.counter)/1000;
end

%% Tracking

% these are variables that we want to track for later analysis.
task{1}{1}.randVars.calculated.correct = nan;
task{1}{1}.randVars.calculated.dir1 = nan; % will be 0->2*pi
task{1}{1}.randVars.calculated.dir2 = nan; % will be 0->2*pi
task{1}{1}.randVars.calculated.shape1 = nan;
task{1}{1}.randVars.calculated.shape2 = nan;
task{1}{1}.randVars.calculated.match = nan;

if ~stimulus.localizer
    task{1}{1}.randVars.calculated.contrast = nan; % will be 0->100%
end

%% Add dead phase

task{1}{2} = task{1}{1};
task{1}{2}.waitForBacktick = 0;
task{1}{1}.numTrials = 1;
task{1}{1}.parameter.contrast = 0;
if stimulus.localizer
    task{1}{1}.seglen = [0 0 0 0 0 9.9];
else
    task{1}{1}.seglen = [0 0 0 0 0 4.9];
end

%% Dots
stimulus.dots.xcenter = 0;
stimulus.dots.ycenter = 0;
stimulus.dots.dotsize = 1;
stimulus.dots.density = 20;
stimulus.dots.speed = 3;
stimulus.idots = initDotsRadial(stimulus.dots);
dots = {};
for i = 1:3
    dots{i} = initDotsRadial(stimulus.dots);
end
stimulus.dots = dots;
stimulus.dot = 1;
%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);
end
%% EYE CALIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~stimulus.localizer
    myscreen = eyeCalibDisp(myscreen);
end

%% Get Ready...
% clear screen    
mglWaitSecs(1);
mglClearScreen(0);
mglFixationCross(1.5,1.5,stimulus.colors.white);
if stimulus.scan        
    mglTextDraw('DO NOT MOVE',[0 1.5]);
end
mglFlush

% let the user know
disp(sprintf('(berlin) Starting run number: %i.',stimulus.counter));
if stimulus.staircasing
    disp(sprintf('(berlin) Staircasing with stimulus contrast: %i/255 and timing: %3.0f',stimulus.lowCon*255,stimulus.run.stimLengths(stimulus.counter)));
end
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
mglClearScreen(0);
mglTextSet([],32,stimulus.colors.white);
mglTextDraw('Run complete...',[0 0]);
if ~stimulus.localizer && ~stimulus.staircasing
    mglTextDraw(sprintf('$%2.2f earned this run',stimulus.run.points/100),[0 2]);
end
mglFlush
myscreen.flushMode = 1;
mglWaitSecs(3);

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

if stimulus.plots
    disp('(berlin) Displaying plots');
    dispInfo(stimulus);
end

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startTrialCallback(task,myscreen)
%%

global stimulus

stimulus.curTrial = stimulus.curTrial + 1;

myscreen.flushMode = 0;

if stimulus.localizer
    % pass, should be already set
elseif stimulus.staircasing
    [task.thistrial.contrast, stimulus.staircase{end}] = doStaircase('testValue',stimulus.staircase{end});
    task.thistrial.contrast = 1-task.thistrial.contrast; % flip the coherence to be mask, instead of visibility
else
    [task.thistrial.contrast, stimulus.istaircase] = doStaircase('testValue',stimulus.istaircase);
end

% for testing
if stimulus.contrastOverride>=0
    task.thistrial.contrast = stimulus.contrastOverride;
end

% setup mask for this trial
stimulus.live.masktex = mglCreateTexture(task.thistrial.contrast*255*(rand(500,500)>0.5));

% set whether this trial matches or not
task.thistrial.match = randi(2)-1;

% okay if we are running task==1 we are doing direction judgments, but we
% can also have changing shapes
% if we are running task==2 we are doing shape judgments, but we can also
% have changing directions
if task.thistrial.task==1
    % DIRECTIONS task
    if stimulus.localizer
        opts = [0 pi/2];
        choice = randi(2);
        task.thistrial.dir1 = opts(choice);
        if task.thistrial.match
            task.thistrial.dir2 = task.thistrial.dir1;
        else
            opts = opts([2 1]);
            task.thistrial.dir2 = opts(choice);
        end
    else
        task.thistrial.dir1 = rand*2*pi;
        task.thistrial.dir2 = task.thistrial.dir1;
        if ~task.thistrial.match
            if randi(2)==1
                task.thistrial.dir2 = task.thistrial.dir1 + pi/2;
            else
                task.thistrial.dir2 = task.thistrial.dir1 - pi/2;
            end
        end
    end
    if stimulus.shape==1
        % randomize shapes
        task.thistrial.shape1 = randi(2); task.thistrial.shape2 = randi(2);
    else
        % same shapes (circles)
        task.thistrial.shape1 = 1; task.thistrial.shape2 = 1;
    end
else
    % SHAPE task
    opts = [1 2]; flip = [2 1]; choice = randi(2);
    if ~task.thistrial.match
        task.thistrial.shape1 = opts(choice);
        task.thistrial.shape2 = flip(choice);
    else
        task.thistrial.shape1 = opts(choice);
        task.thistrial.shape2 = opts(choice);
    end
    task.thistrial.dir1 = rand*2*pi;
    if randi(2)==1
        task.thistrial.dir2 = task.thistrial.dir1;
    else
        if randi(2)==1
            task.thistrial.dir2 = task.thistrial.dir1 + pi/2;
        else
            task.thistrial.dir2 = task.thistrial.dir1 - pi/2;
        end
    end
end
    
matches = {'No','Match'};
disp(sprintf('Trial %i Dir1: %i Dir2: %i Mask: %3.0f/255 Matching: %s',stimulus.curTrial,round(180/pi*task.thistrial.dir1),round(180/pi*task.thistrial.dir2),task.thistrial.contrast*255,matches{task.thistrial.match+1}));

stimulus.live.eyeCount = 0;
stimulus.dead = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)
%%

global stimulus

stimulus.live.dots = 0;
stimulus.live.resp = 0;
stimulus.live.mask = 0;
stimulus.live.fixColor = stimulus.colors.white;
stimulus.live.dotColor = stimulus.lowCon;
stimulus.live.coherence = 0;

if task.thistrial.thisphase==1
    return
end

if any([stimulus.seg.mask1 stimulus.seg.mask2]==task.thistrial.thisseg)
    stimulus.live.mask = 1;
elseif any(stimulus.seg.stim==task.thistrial.thisseg)
    stimulus.live.dots = 1;
    stimulus.dot = stimulus.dot+1; if stimulus.dot>3, stimulus.dot=1; end
    stimulus.dots{stimulus.dot}.dir = task.thistrial.dir1;
    stimulus.live.coherence = 1;
    stimulus.live.shape = task.thistrial.shape1;
elseif stimulus.seg.resp==task.thistrial.thisseg
    stimulus.live.dots = 1;
    stimulus.dot = stimulus.dot+1; if stimulus.dot>3, stimulus.dot=1; end
    stimulus.dots{stimulus.dot}.dir = task.thistrial.dir2;
    stimulus.live.coherence = 1;
    stimulus.live.shape = task.thistrial.shape2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus
mglClearScreen(0);
% check eye pos
if ~stimulus.noeye && task.thistrial.thisseg~=stimulus.seg.ITI && ~stimulus.localizer
    [pos,time] = mglEyelinkGetCurrentEyePos;
    if ~any(isnan(pos))
        dist = hypot(pos(1),pos(2));
        if dist > stimulus.ring.inner && stimulus.live.eyeCount > 30
            mglTextSet([],32,stimulus.colors.red);
            disp('Eye movement detected!!!!');
            mglTextDraw('Eye Movement Detected',[0 0]);
            mglFlush
            myscreen.flushMode = 1;
            stimulus.dead = 1;
            return
        elseif dist > stimulus.ring.inner-1
            stimulus.live.eyeCount = stimulus.live.eyeCount + 1;
        end
    end
end
% stimulus
% first draw the outer (always incoherent) dots
if stimulus.constant && ~stimulus.live.mask
	mglStencilSelect(999);
    stimulus = upDotsInc(stimulus,myscreen);
    mglStencilSelect(0);
end
% block the inside for the sometimes not incoherent dots
% now draw everything else
if stimulus.live.dots
    if stimulus.live.shape==2 % cross
        mglStencilSelect(100);
    else % circle
        mglStencilSelect(99);
    end
    mglFillOval(0,0,repmat(stimulus.ring.outer+3,1,2),0);
    stimulus = upDots(stimulus,myscreen);
    mglStencilSelect(0);
end
if stimulus.live.mask
    mglStencilSelect(999);
    mglBltTexture(stimulus.live.masktex,[0 0]);
    mglStencilSelect(0);
else
    mglFillOval(0,0,repmat(stimulus.ring.inner,1,2),0);
end
upFix(stimulus);

function upFix(stimulus)
%%
if ~stimulus.localizer && ~stimulus.staircasing && all(stimulus.live.fixColor==stimulus.colors.green)
    mglTextSet([],32,stimulus.live.fixColor);
    mglTextDraw('+7',[0 0]);
else
    mglFixationCross(1.5,1.5,stimulus.live.fixColor);
end

function stimulus = upDots(stimulus,myscreen)

stimulus.dots{stimulus.dot} = updateDotsRadial(stimulus.dots{stimulus.dot},stimulus.live.coherence,myscreen,true);

mglPoints2(stimulus.dots{stimulus.dot}.x,stimulus.dots{stimulus.dot}.y,...
    stimulus.dots{stimulus.dot}.dotsize,stimulus.live.dotColor);

function stimulus = upDotsInc(stimulus,myscreen)
stimulus.idots = updateDotsRadial(stimulus.idots,0,myscreen,true);
mglPoints2(stimulus.idots.x,stimulus.idots.y,...
    stimulus.idots.dotsize,stimulus.live.dotColor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [task, myscreen] = getResponseCallback(task, myscreen)


global stimulus
if stimulus.dead, return; end
responseText = {'Incorrect','Correct'};
fixColors = {stimulus.colors.red,stimulus.colors.green};

if stimulus.localizer
    
    if any(task.thistrial.whichButton == stimulus.responseKeys)
        if task.thistrial.gotResponse == 0
            task.thistrial.correct = task.thistrial.whichButton == stimulus.responseKeys(task.thistrial.match+1);
            disp(sprintf('Subject pressed %i: %s',task.thistrial.whichButton,responseText{task.thistrial.correct+1}));
            stimulus.live.fixColor = fixColors{task.thistrial.correct+1};
        else
            disp(sprintf('(berlin) Subject responded multiple times: %i',task.thistrial.gotResponse+1));
        end
        stimulus.live.resp = 0;
        stimulus.live.dots = 0;
    end
else

    if any(task.thistrial.whichButton == stimulus.responseKeys)
        if task.thistrial.gotResponse == 0
            task.thistrial.correct = task.thistrial.whichButton == stimulus.responseKeys(task.thistrial.match+1);
            disp(sprintf('Subject pressed %i: %s',task.thistrial.whichButton,responseText{task.thistrial.correct+1}));
            stimulus.live.fixColor = fixColors{task.thistrial.correct+1};
            % Store whether this was correct

            if stimulus.staircasing
                stimulus.staircase{end} = doStaircase('update',stimulus.staircase{end},task.thistrial.correct);
            else
                stimulus.istaircase = doStaircase('update',stimulus.istaircase,task.thistrial.correct);
                if task.thistrial.correct
                    stimulus.run.points = stimulus.run.points + 7;
                end
            end
        else
            disp(sprintf('(berlin) Subject responded multiple times: %i',task.thistrial.gotResponse+1));
        end
        stimulus.live.dots = 0;
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
stimulus.staircase = {};
stimulus.staircase{1} = doStaircase('init','upDown',...
        'initialThreshold',1,...
        'initialStepsize',0.1,...
        'minThreshold=0.001','maxThreshold=1','stepRule','pest',...
        'nTrials=50','maxStepsize=0.33','minStepsize=0.001');
stimulus.istaircase = doStaircase('init','fixed','fixedVals',32/255,'nTrials=50');

function stimulus = resetStair(stimulus)

if stimulus.staircasing
    disp('(berlin) Initializing new staircase');
    stimulus.staircase{end+1} = doStaircase('init','upDown',...
        'initialThreshold',1,...
        'initialStepsize',0.33,...
        'minThreshold=0.001','maxThreshold=1','stepRule','pest',...
        'nTrials=50','maxStepsize=0.33','minStepsize=0.001');
end
if doStaircase('stop',stimulus.istaircase)
    disp('(berlin) Resetting invisible staircase');
    stimulus.istaircase(end+1) = doStaircase('init','fixed','fixedVals',32/255,'nTrials=50');
end

function stimulus = initRuns(stimulus)
% initialize the run info for the staircasing mode
stimLengths = [100];
stimCon = [2];

sl = []; sc = [];

for i = 1:length(stimLengths)
    for j = 1:length(stimCon)
        sl(end+1) = stimLengths(i);
        sc(end+1) = stimCon(j);
    end
end

idxs = randperm(length(sl));
sl = sl(idxs); sc = sc(idxs);

stimulus.counter = 0;
stimulus.run.stimLengths = sl;
stimulus.run.stimCon = sc;

%%%%%%%%%%%%%%%%%%%%%%%
%    dispInfo    %
%%%%%%%%%%%%%%%%%%%%%%%
function dispInfo(stimulus)
%%

if ~stimulus.localizer && ~stimulus.staircasing
    disp(sprintf('Participant %s has earned $%2.2f',mglGetSID,stimulus.run.points/100));
end
% load the luminance table
% % % load(myscreen.calibFullFilename)
% % % luminance = interp1(calib.tableCorrected.outputValues,calib.tableCorrected.luminance,0:1/255:255);
if stimulus.staircasing
    %%
    notstaircase = stimulus.staircase;
    thresholds = zeros(size(stimulus.run.stimLengths));
    for i = 1:length(stimulus.staircase)
        out = doStaircase('threshold',notstaircase{i},'type','weibull','dispFig=0');
        thresholds(i) = out.threshold;
    end
    % reorganize into matrix
    stimCons = unique(stimulus.run.stimCon);
    stimCons = sort(stimCons);
    stimLengths = unique(stimulus.run.stimLengths);
    stimLengths = sort(stimLengths);
    datamat = nan(length(stimCons),length(stimLengths),5);
    for ci = 1:length(stimCons)
        for li = 1:length(stimLengths)
            idxs = logical((stimulus.run.stimLengths==stimLengths(li)) .* (stimulus.run.stimCon==stimCons(ci)));
            datamat(ci,li,1:sum(idxs)) = thresholds(idxs);
        end
    end
    if any(thresholds<0) || any(thresholds>1)
        % remove errant thresholds
        warning('should remove some thresholds...');
    end
    datamat(datamat>1) = NaN;
    datamat(datamat<=0) = NaN;
    %%
    datamu = nanmean(datamat,3);
    datamu(datamu==0) = NaN;
    datamu = round((1-datamu)*255);
    datasd = nanstd(datamat,[],3);
    datasd(datasd==0) = NaN;
    %% plot
    cmap = brewermap(length(stimCons)+1,'Purples');
    cmap = cmap(2:end,:);
    figure, hold on
    legs = {};
    for i = 1:length(stimCons)
        plot(stimLengths,datamu(i,:),'o','MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',[1 1 1],'MarkerSize',10);
        errbar(stimLengths,datamu(i,:),datasd(i,:),'-','Color',cmap(i,:));
        legs{end+1} = sprintf('Stimulus luminance: %i/255',stimCons(i));
    end
    a = axis;
    axis([50 100 0 a(4)]);
    legend(legs)
    xlabel('Stimulus length (ms)');
    ylabel('Mask contrast at just noticeable difference (% luminance)');
    drawPublishAxis
elseif stimulus.localizer
else
%     perf = zeros(size(stimulus.istaircase));
%     for i = 1:length(stimulus.istaircase)
%         perf(i) = mean(stimulus.istaircase(i).response);
%     end
%     figure
%     plot(1:length(perf),perf,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
%     set(gca,'XAxisTick',1:length(perf));
%     drawPublishAxis
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create dots for optic flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initDotsRadial(dots,~)

% maximum depth of points
dots.minX = -10;
dots.maxX = 10;
dots.minY = -10;
dots.maxY = 10;

dots.dir = 0;

area = (dots.maxX-dots.minX)*(dots.maxY-dots.minY);

dots.n = area * dots.density;

dots.x = rand(1,dots.n)*(dots.maxX-dots.minX)+dots.minX;
dots.y = rand(1,dots.n)*abs(dots.maxY-dots.minY)+dots.minY;

dots.x = dots.x;
dots.y = dots.y;

% set incoherent dots to 0
dots.coherency = 1;
dots.incoherent = rand(1,dots.n) > dots.coherency;
dots.incoherentn = sum(dots.incoherent);
dots.coherent = ~dots.incoherent;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots for Radial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updateDotsRadial(dots,coherence,myscreen,repick)

% get the coherent and incoherent dots
if coherence==1
    dots.coherent = true(size(dots.x));
    dots.incoherent = false(size(dots.x));
    dots.incoherentn = 0;
elseif repick
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

if dots.coherentn>0
    % move coherent dots
    % dots.dir is an angle, so compute x/y based on this
    dots.y(dots.coherent) = dots.y(dots.coherent) + (freq_factor * sin(dots.dir));
    dots.x(dots.coherent) = dots.x(dots.coherent) + (freq_factor * cos(dots.dir));
end

if dots.incoherentn>0
    % these are for flipping into the other quadrants
    rdir = rand(1,dots.incoherentn)*pi*2;
    dots.y(dots.incoherent) = dots.y(dots.incoherent) + (freq_factor * sin(rdir));
    dots.x(dots.incoherent) = dots.x(dots.incoherent) + (freq_factor * cos(rdir));
end

offscreen = dots.x > dots.maxX;
dots.x(offscreen) = dots.x(offscreen) - abs(dots.maxX - dots.minX);
offscreen = dots.x < dots.minX;
dots.x(offscreen) = dots.x(offscreen) + abs(dots.maxX - dots.minX);

offscreen = dots.y > dots.maxY;
dots.y(offscreen) = dots.y(offscreen) - abs(dots.maxY - dots.minY);
offscreen = dots.y < dots.minY;
dots.y(offscreen) = dots.y(offscreen) + abs(dots.maxY - dots.minY);


function setGT(myscreen,stimulus)

% multipliers relative
max = 95;
% load the calibration
load(myscreen.calibFullFilename);

% get the max value from the calibration file
localMax = calib.tableCorrected.luminance(end);

% get an approximate factor
factor = round(localMax/max);

disp(sprintf('Setting gamma table to use %02.2f%% of available space.',1/factor*100));
% if not 1, set the gammaTable accordingly

if factor > 1
    newTable = interp1(linspace(0,1,256),stimulus.linearizedGammaTable.redTable,linspace(0,1/factor,256));
    succ = mglSetGammaTable(repmat(newTable',1,3));

    if ~succ
        warning('Gamma table set failure');
        keyboard
    end
end
