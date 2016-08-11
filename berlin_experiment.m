function [ myscreen ] = berlin_experiment( varargin )
%CAT_AWE Testing category discrimination at low visibility conditions
%   This is dan's experiment for Berlin (fall 2016). Version 2! Now this is
%   a motion experiment. On each trial (15 s) a mask/stimulus pair occilate
%   on off every 250 ms for 6 s. The last 9 s are a response and ITI time.
%
%   OPTIONS
%
%   default (localizer=0, invisible=0) runs in staircase mode. The stimulus
%   is shown at 50% coherence and 5% contrast. The staircase alters the
%   contrast of the mask (0% coherence).
%
%   invisible=1 using staircases estimated by the default mode this
%   interpolates a mask contrast value to guarantee that the motion percept
%   is invisible.
%
%   localizer=1 runs a localizer mode. Uses either a mask contrast of 0%,
%   equal to dots (1/255), or invisible (15/255). The subject performs a
%   horizontal vs. vertical task and both motion go in the same direction.

global stimulus

%% Initialize Variables

% add arguments later
localizer = 0;
staircase = 0;
scan = 0;
plots = 0;
getArgs(varargin,{'localizer=0','staircase=0','scan=0','plots=0'});
stimulus.localizer = localizer;
stimulus.scan = scan;
stimulus.staircasing = staircase;
stimulus.plots = plots;
clear localizer invisible scan

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
        stimulus.counter = s.stimulus.counter + 1;

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

%% Initialize Stimulus

myscreen = initStimulus('stimulus',myscreen);

if stimulus.scan
    stimulus.responseKeys = [2 1]; % corresponds to NOMATCH, MATCH
else
    stimulus.responseKeys = [2 1]; % corresponds to  NOMATCH, MATCH
end

stimulus.colors.black = [0 0 0];
stimulus.colors.white = [0.25 0.25 0.25];
stimulus.colors.green = [0 0.25 0];
stimulus.colors.red = [0.25 0 0];

stimulus.ring.inner = 4;
stimulus.ring.outer = 8; %
stimulus.area = 3.14159265358979*(stimulus.ring.outer^2-stimulus.ring.inner^2);

stimulus.lowCon = 1/255; % minimum possible contrast
stimulus.contrastOverride = 0;

%% Generate stencils
mglStencilCreateBegin(99);
mglFillOval(0,0,repmat(stimulus.ring.outer,1,2),1);
mglStencilCreateEnd;

%% Setup Task

task{1}{1}.waitForBacktick = 1;

stimulus.curTrial = 0;

% calculate timing
stimulus.seg.mask1 = 1;
stimulus.seg.mask2 = 3;
stimulus.seg.mask3 = 5;
stimulus.seg.stim1 = 2;
stimulus.seg.stim2 = 4;
stimulus.seg.delay = -1; % not using yet (delay before resp)
stimulus.seg.resp = 6;
stimulus.seg.ITI = 7;
task{1}{1}.segmin = [0.100 0.100 0.100 0.100 0.100 2 1];
task{1}{1}.segmax = [0.100 0.100 0.100 0.100 0.100 2 3];
if stimulus.localizer
%     stimulus.seg.delay = 5;
%     stimulus.seg.resp = stimulus.seg.resp+1; stimulus.seg.ITI = stimulus.seg.ITI+1;
    stimulus.seg.mask3 = -1;
    stimulus.seg.delay = 4; 
    stimulus.seg.stim2 = 5;
    task{1}{1}.segmin = [0.100 0.100 0.100 3 0.200 2 1];
    task{1}{1}.segmax = [0.100 0.100 0.100 3 0.200 2 1];
end

task{1}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
if stimulus.scan
    task{1}{1}.synchToVol(stimulus.seg.ITI) = 1;
end
task{1}{1}.getResponse = zeros(size(task{1}{1}.segmin)); task{1}{1}.getResponse(stimulus.seg.resp)=1;

task{1}{1}.numTrials = 60;

if stimulus.localizer
    task{1}{1}.parameter.horiz = [0 1];
    task{1}{1}.parameter.contrast = 1-[0 1/255 15/255];
    task{1}{1}.random = 1;
    task{1}{1}.numTrials = Inf;
end

%% Tracking

% these are variables that we want to track for later analysis.
task{1}{1}.randVars.calculated.correct = nan;
task{1}{1}.randVars.calculated.dir1 = nan; % will be 0->2*pi
task{1}{1}.randVars.calculated.dir2 = nan; % will be 0->2*pi
if ~stimulus.localizer
    task{1}{1}.randVars.calculated.contrast = nan; % will be 0->100%
    task{1}{1}.randVars.calculated.match = nan;
end

%% Dots
stimulus.dots.xcenter = 0;
stimulus.dots.ycenter = 0;
stimulus.dots.dotsize = 2;
stimulus.dots.density = 15;
stimulus.dots.speed = 3;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init staircase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if stimulus.initStair
    % We are starting our staircases from scratch
    disp(sprintf('(berlin) Initializing staircase'));
    stimulus = initStaircase(stimulus);
else
    disp('(berlin) Checking staircase resets from previous run...');
    stimulus = resetStair(stimulus);
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
mglWaitSecs(1);
mglClearScreen(0);
if stimulus.scan        
    mglTextDraw('DO NOT MOVE',[0 1.5]);
end
mglFlush

% let the user know
disp(sprintf('(berlin) Starting run number: %i.',stimulus.counter));
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
mglTextDraw('Run complete...',[0 0]);
mglFlush
myscreen.flushMode = 1;
mglWaitSecs(1);

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

if stimulus.scan
    task.thistrial.seglen(end) = 1.05^(rand*30+15);
end

stimulus.curTrial = stimulus.curTrial + 1;

myscreen.flushMode = 0;

% set the current image
task.thistrial.dir1 = rand*2*pi;
task.thistrial.dir2 = task.thistrial.dir1;
if stimulus.localizer
    % pass
elseif stimulus.staircasing
    [task.thistrial.contrast, stimulus.staircase] = doStaircase('testValue',stimulus.staircase);
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

if ~stimulus.localizer
    task.thistrial.match = randi(2)-1;
    if ~task.thistrial.match
        if randi(2)==1
            task.thistrial.dir2 = task.thistrial.dir1 + pi/2;
        else
            task.thistrial.dir2 = task.thistrial.dir1 - pi/2;
        end
    end
    matches = {'No','Match'};
    disp(sprintf('Trial %i Dir1: %i Dir2: %i Contrast %3.0f/255 Match %s',stimulus.curTrial,round(180/pi*task.thistrial.dir1),round(180/pi*task.thistrial.dir2),task.thistrial.contrast*255,matches{task.thistrial.match+1}));
elseif task.thistrial.horiz==1
    % horizontal
    if randi(2)==1
        task.thistrial.dir1 = 0; task.thistrial.dir2 = task.thistrial.dir1;
    else
        task.thistrial.dir1 = pi; task.thistrial.dir2 = task.thistrial.dir1;
    end
    horizs = {'No','Yes'};
    disp(sprintf('Trial %i Dir1: %i Dir2: %i Contrast %3.0f/255 Horizontal %s',stimulus.curTrial,round(180/pi*task.thistrial.dir1),round(180/pi*task.thistrial.dir2),task.thistrial.contrast*255,horizs{task.thistrial.horiz+1}));
else
    % vertical
    if randi(2)==1
        task.thistrial.dir1 = pi/2; task.thistrial.dir2 = task.thistrial.dir1;
    else
        task.thistrial.dir1 = pi*1.5; task.thistrial.dir2 = task.thistrial.dir1;
    end
    horizs = {'No','Yes'};
    disp(sprintf('Trial %i Dir1: %i Dir2: %i Contrast %3.0f/255 Horizontal %s',stimulus.curTrial,round(180/pi*task.thistrial.dir1),round(180/pi*task.thistrial.dir2),task.thistrial.contrast*255,horizs{task.thistrial.horiz+1}));
end

if stimulus.localizer
    if randi(2)==1
        stimulus.live.respPos = [-5 5];
        stimulus.responseKeys = [2 1];
    else
        stimulus.live.respPos = [5 -5];
        stimulus.responseKeys = [1 2];
    end
end
    
stimulus.dots{1}.dir = task.thistrial.dir1;
stimulus.dots{2}.dir = task.thistrial.dir2;
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

if any([stimulus.seg.mask1 stimulus.seg.mask2 stimulus.seg.mask3]==task.thistrial.thisseg)
    stimulus.live.mask = 1;
elseif stimulus.seg.stim1==task.thistrial.thisseg
    stimulus.live.dots = 1;
    stimulus.dot = 1;
    stimulus.live.coherence = 1;
elseif stimulus.seg.stim2==task.thistrial.thisseg
    stimulus.live.dots = 1;
    stimulus.dot = 2;
    stimulus.live.coherence = 1;
elseif stimulus.seg.resp==task.thistrial.thisseg
    stimulus.live.resp = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus

mglClearScreen(0);
mglStencilSelect(99);
if stimulus.live.dots
    stimulus = upDots(stimulus,myscreen);
end
if stimulus.live.mask
    mglBltTexture(stimulus.live.masktex,[0 0]);
end
mglFillOval(0,0,repmat(stimulus.ring.inner,1,2),0);
upFix(stimulus);
mglStencilSelect(0);
if stimulus.live.resp && stimulus.localizer
    mglTextSet([],32,stimulus.live.fixColor);
    mglTextDraw('Horizontal',[stimulus.live.respPos(1) 0]);
    mglTextDraw('Vertical',[stimulus.live.respPos(2) 0]);
end

function upFix(stimulus)
%%
mglFixationCross(1.5,1.5,stimulus.live.fixColor);

function stimulus = upDots(stimulus,myscreen)

stimulus.dots{stimulus.dot} = updateDotsRadial(stimulus.dots{stimulus.dot},stimulus.live.coherence,myscreen,true);

mglPoints2(stimulus.dots{stimulus.dot}.xdisp,stimulus.dots{stimulus.dot}.ydisp,...
    stimulus.dots{stimulus.dot}.dotsize,stimulus.live.dotColor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus
responseText = {'Incorrect','Correct'};
fixColors = {stimulus.colors.red,stimulus.colors.green};

if stimulus.localizer
    
    if any(task.thistrial.whichButton == stimulus.responseKeys)
        if task.thistrial.gotResponse == 0
            task.thistrial.correct = task.thistrial.whichButton == stimulus.responseKeys(task.thistrial.horiz+1);
            disp(sprintf('Subject pressed %i: %s',task.thistrial.whichButton,responseText{task.thistrial.correct+1}));
            stimulus.live.fixColor = fixColors{task.thistrial.correct+1};
        else
            disp(sprintf('(berlin) Subject responded multiple times: %i',task.thistrial.gotResponse+1));
        end
        stimulus.live.resp = 0;
    end
else

    if any(task.thistrial.whichButton == stimulus.responseKeys)
        if task.thistrial.gotResponse == 0
            task.thistrial.correct = task.thistrial.whichButton == stimulus.responseKeys(task.thistrial.match+1);
            disp(sprintf('Subject pressed %i: %s',task.thistrial.whichButton,responseText{task.thistrial.correct+1}));
            stimulus.live.fixColor = fixColors{task.thistrial.correct+1};
            % Store whether this was correct

            if stimulus.staircasing
                stimulus.staircase = doStaircase('update',stimulus.staircase,task.thistrial.correct);
            else
                stimulus.istaircase = doStaircase('update',stimulus.istaircase,task.thistrial.correct);
            end
        else
            disp(sprintf('(berlin) Subject responded multiple times: %i',task.thistrial.gotResponse+1));
        end
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

stimulus.staircase = doStaircase('init','upDown',...
        'initialThreshold',1,...
        'initialStepsize',0.33,...
        'minThreshold=0.001','maxThreshold=1','stepRule','pest',...
        'nTrials=50','maxStepsize=0.33','minStepsize=0.001');
stimulus.istaircase = doStaircase('init','fixed','fixedVals',15/255,'nTrials=50');

function stimulus = resetStair(stimulus)

if doStaircase('stop',stimulus.staircase)
    disp('(berlin) Resetting staircase');
    out = doStaircase('threshold',stimulus.staircase);
    stimulus.staircase(end+1) = doStaircase('init','upDown',...
        'initialThreshold',out.threshold,...
        'initialStepsize',out.threshold/3,...
        'minThreshold=0.001','maxThreshold=1','stepRule','pest',...
        'nTrials=50','maxStepsize=0.33','minStepsize=0.001');
end
if doStaircase('stop',stimulus.istaircase)
    disp('(berlin) Resetting invisible staircase');
    stimulus.istaircase(end+1) = doStaircase('init','fixed','fixedVals',15/255,'nTrials=50');
end
%%%%%%%%%%%%%%%%%%%%%%%
%    dispInfo    %
%%%%%%%%%%%%%%%%%%%%%%%
function dispInfo(stimulus)
%%
if stimulus.staircasing
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

dots.xdisp = dots.x;
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
end

offscreen = dots.x > dots.maxX;
dots.x(offscreen) = dots.x(offscreen) - abs(dots.maxX - dots.minX);
offscreen = dots.x < dots.minX;
dots.x(offscreen) = dots.x(offscreen) + abs(dots.maxX - dots.minX);

offscreen = dots.y > dots.maxY;
dots.y(offscreen) = dots.y(offscreen) - abs(dots.maxY - dots.minY);
offscreen = dots.y < dots.minY;
dots.y(offscreen) = dots.y(offscreen) + abs(dots.maxY - dots.minY);

dots.xdisp = dots.x;
dots.ydisp = dots.y;