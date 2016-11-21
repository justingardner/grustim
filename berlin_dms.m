function [ myscreen ] = berlin_dms( varargin )
%CAT_AWE Testing match-to-sample for motion
%   This is dan's experiment for Berlin (fall 2016). Version 4!
%
%   Delayed match-to-sample w/ motion stimulus or gratings
%
%   No constant background (but optional for motion, e.g. for scanning)

global stimulus

%% OVERRIDES (for testing)

% run: `

% set all to -1 when running:
stimulus.contrastOverride = 128/255;
% stimulus.timingOverride = 1;

%% Initialize Variables

% add arguments later
scan = 0;
plots = 0; windowed=0;
noeye = 0; framegrab=0;
feature = 0; fixate = 0; trigger=0;
getArgs(varargin,{'scan=0','plots=0','noeye=1','fixate=0','constant=0','framegrab=0','complex=0','feature=1','trigger=0','windowed=0'});
stimulus.framegrab = framegrab;
stimulus.fixate = fixate;
stimulus.scan = scan;
stimulus.plots = plots;
stimulus.noeye = noeye;
stimulus.constant = constant;
stimulus.trigger = trigger;
stimulus.feature = feature; % 1 = motion, 2 = gratings
clear localizer invisible scan category noeye task

stimulus.counter = 1; % This keeps track of what "run" we are on.

%% Setup Screen

if windowed
    myscreen = initScreen('windowed');
else
    if stimulus.scan
        myscreen = initScreen('fMRIprojFlex');
    else
        myscreen = initScreen('Rolfs_VPixx');
    end
end

myscreen.background = 0;

if stimulus.framegrab
    deg2pix = myscreen.screenWidth/myscreen.imageWidth;
    total = round(20*deg2pix); if mod(total,2)==1, total=total+1; end
    stimulus.frame.size = total;
    stimulus.frame.frames = zeros(total,total,1000);
    stimulus.frame.count = 1;
    x = myscreen.screenWidth/2;
    y = myscreen.screenHeight/2;
    stimulus.frame.coords = [x-total/2 y-total/2 total total];
end

%% Eye calib
myscreen = eyeCalibDisp(myscreen);

%% Open Old Stimfile
stimulus.initStair = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/berlin_dms/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/berlin_dms/%s/1*mat',mglGetSID));

    if length(files) >= 1
        fname = files(end).name;
        
        s = load(sprintf('~/data/berlin_dms/%s/%s',mglGetSID,fname));
        % copy staircases and run numbers
        stimulus.staircase = s.stimulus.staircase;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% init staircase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if stimulus.initStair
    % We are starting our staircases from scratch
    disp(sprintf('(berlin) Initializing staircase'));
    stimulus = initStaircase(stimulus);
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
if stimulus.contrastOverride>0

    stimulus.colors.white = [1 1 1];
    stimulus.colors.green = [0 1 0];
    stimulus.colors.red = [1 0 0];
end
% stimulus.colors.black = 0;
% stimulus.colors.white = 1/255;
% stimulus.colors.red = 2/255;
% stimulus.colors.green = 3/255;

stimulus.ring.inner = 3.5;
stimulus.ring.outer = 5; %21
% stimulus.area = 3.14159265358979*((stimulus.ring.outer/2)^2-(stimulus.ring.inner/2)^2);

%% Jitter (to avoid making shit too easy)
stimulus.jitter = 0.3;

%% Motion Parameters
stimulus.orientations = [0 135];

%% Contrast
stimulus.contrast = stimulus.contrastOverride;

%% Grating parameters

stimulus.grating.radius = 5;
stimulus.grating.targetLoc = [1 4];
stimulus.grating.sf = 2;
stimulus.grating.tf = 0.5;
stimulus.grating.width = 11;
stimulus.grating.height = 11;
% stimulus.grating.phases = [0 pi];
stimulus.grating.phases = [0:2*pi/30:2*pi 2*pi:-2*pi/30:0]; % note this makes the actual tf = tf*2
stimulus.grating.phase = 0;
stimulus.grating.windowType = 'gabor'; % should be gabor or thresh
stimulus.grating.sdx = stimulus.grating.width/7;
stimulus.grating.sdy = stimulus.grating.width/7;
stimulus.x = 0;

stimulus.y = 1.5;
stimulus.grating.width = 8.5;
stimulus.grating.height = 8.5;
stimulus.grating.windowType = 'thresh'; % should be gabor or thresh
stimulus.grating.sdx = stimulus.grating.width/2;
stimulus.grating.sdy = stimulus.grating.height/2;

stimulus = initGratings(stimulus,myscreen);

%% Generate stencils

% stimulus.posAngles = [0 45 135 180]*pi/180;
stimulus.posAngles = 0;
stimulus.posDist = 0;
% stimulus.posAngles = [0 180]*pi/180;
% stimulus.posDist = 5;
stimulus.posopts = stimulus.posDist * [cos(stimulus.posAngles)' sin(stimulus.posAngles)'];

for i = 1:length(stimulus.posAngles)
    mglStencilCreateBegin(i);
    mglFillOval(stimulus.posopts(i,1),stimulus.posopts(i,2),repmat(stimulus.ring.outer+1,1,2),1);
    mglStencilCreateEnd;
end
%% Setup Task
task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;

stimulus.curTrial = 0;

task{1}{1}.segmin = [3 2 3 2 2];
task{1}{1}.segmax = [3 2 3 2 4];
if stimulus.trigger
    % set the delay intervals to infinite, eye position will be used to
    % trigger
    task{1}{1}.segmin = [3 inf 3 inf];
    task{1}{1}.segmax = [3 inf 3 inf];
end

stimulus.seg.stim1 = 1;
stimulus.seg.delay = 2;
% stimulus.seg.stim2 = 3;
stimulus.seg.resp = 3;
stimulus.seg.ITI = 4;

task{1}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
task{1}{1}.getResponse = zeros(size(task{1}{1}.segmin)); task{1}{1}.getResponse(stimulus.seg.resp)=1;
task{1}{1}.numTrials = Inf;
task{1}{1}.random = 1;
task{1}{1}.parameter.orientation = 1:length(stimulus.orientations); % which orientation to use (0 deg or 135 deg)
task{1}{1}.parameter.rotation = [-1 1];
task{1}{1}.parameter.contrast = stimulus.contrast;

if stimulus.scan
    task{1}{1}.synchToVol(stimulus.seg.ITI) = 1;
end

% if stimulus.timingOverride
%     task{1}{1}.seglen(stimulus.seg.delay) = stimulus.timingOverride;
%     task{1}{1}.seglen(stimulus.seg.ITI) = stimulus.timingOverride;
% end 

%% Tracking

% these are variables that we want to track for later analysis.
task{1}{1}.randVars.calculated.correct = nan;
task{1}{1}.randVars.calculated.dir1 = nan; % will be 0->2*pi
task{1}{1}.randVars.calculated.dir2 = nan; % will be 0->2*pi
% task{1}{1}.randVars.calculated.match = nan;
task{1}{1}.randVars.calculated.posopt = nan;
task{1}{1}.randVars.calculated.x = nan;
task{1}{1}.randVars.calculated.y = nan;

%% Add dead phase (not necessary without adaptation period)

% task{1}{2} = task{1}{1};
% task{1}{2}.waitForBacktick = 0;
% task{1}{1}.numTrials = 1;
% task{1}{1}.parameter.contrast = 0;
% if stimulus.scan
%     task{1}{1}.seglen = [0 0 0 0 0 9.9];
% else
%     task{1}{1}.seglen = [0 0 0 0 0 4.9];
% end

%% Dots
stimulus.dots = struct;
stimulus.dots.xcenter = 0;
stimulus.dots.ycenter = 0;
stimulus.dots.dotsize = 1;
stimulus.dots.density = 20;
stimulus.dots.speed = 3;

stimulus.dots.minX = -10;
stimulus.dots.maxX = 10;
stimulus.dots.minY = -10;
stimulus.dots.maxY = 10;

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
if ~stimulus.noeye
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

stimulus.live.phaseNum = randi(62);

stimulus.curTrial = stimulus.curTrial + 1;

myscreen.flushMode = 0;

[rotation, stimulus.staircase{stimulus.feature}] = doStaircase('testValue',stimulus.staircase{stimulus.feature});

stimulus.live.dotColor = task.thistrial.contrast;
% set whether this trial matches or not
task.thistrial.dir1 = stimulus.orientations(task.thistrial.orientation)*pi/180 + randn*stimulus.jitter;
task.thistrial.dir2 = task.thistrial.dir1 + rotation*task.thistrial.rotation;

% set the x/y coordinates of the gratings
% % % % % % pick = randi(size(stimulus.posopts,1));
% % % % % % task.thistrial.x = stimulus.posopts(pick,1);
% % % % % % task.thistrial.y = stimulus.posopts(pick,2);
% % % % % % task.thistrial.posopt = pick;

disp(sprintf('Trial %i Dir1: %i Dir2: %i',stimulus.curTrial,round(180/pi*task.thistrial.dir1),round(180/pi*task.thistrial.dir2)));

stimulus.live.eyeCount = 0;
stimulus.dead = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)
%%

global stimulus

stimulus.live.triggerWaiting = 0;
if stimulus.trigger && any(task.thistrial.thisseg==[stimulus.seg.delay stimulus.seg.ITI])
    stimulus.live.triggerWaiting = 1;
    stimulus.live.centered = 0;
    stimulus.live.triggerTime = 0;
    stimulus.live.lastTrigger = -1;
end

stimulus.live.dots = stimulus.constant;
stimulus.live.resp = 0;
stimulus.live.fixColor = stimulus.colors.white;
stimulus.live.coherence = 1;
stimulus.live.x = 0;
stimulus.live.y = 0;
% stimulus.live.x = task.thistrial.x;
% stimulus.live.y = task.thistrial.y;
stimulus.live.fix = 1;

% use only for constant mode
% if task.thistrial.thisphase==1
%     return
% end

if task.thistrial.thisseg==stimulus.seg.stim1
    stimulus.live.dots = 1;
    stimulus.dot = stimulus.dot+1; if stimulus.dot>3, stimulus.dot=1; end
    stimulus.dots{stimulus.dot}.dir = task.thistrial.dir1;
    stimulus.live.dir = task.thistrial.dir1; % for gratings
    stimulus.live.fix = 0;
elseif task.thistrial.thisseg==stimulus.seg.resp
    stimulus.live.dots = 1;
    stimulus.dot = stimulus.dot+1; if stimulus.dot>3, stimulus.dot=1; end
    stimulus.dots{stimulus.dot}.dir = task.thistrial.dir2;
    stimulus.live.dir = task.thistrial.dir2; % for gratings
    stimulus.live.fix = 0;
% elseif task.thistrial.thisseg==stimulus.seg.resp
%     stimulus.live.fix = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus
mglClearScreen(0);
% check eye pos
if ~stimulus.noeye && task.thistrial.thisseg~=stimulus.seg.ITI && ~stimulus.scan
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
% if stimulus.constant
% 	mglStencilSelect(999);
%     stimulus = upDotsInc(stimulus,myscreen);
%     mglStencilSelect(0);
% end
% block the inside for the sometimes not incoherent dots
% now draw everything else
if stimulus.live.dots
%     mglStencilSelect(task.thistrial.posopt);
    mglStencilSelect(1);
%     mglFillOval(0,0,repmat(stimulus.ring.outer+3,1,2),0);
    if stimulus.feature==1
        stimulus = upDots(stimulus,myscreen);
    else
        upGrating(stimulus,task);
    end
    mglStencilSelect(0);
end
if stimulus.fixate || stimulus.live.fix || stimulus.live.resp==1
%     mglFillOval(0,0,repmat(stimulus.ring.inner,1,2),0);
    upFix(stimulus);
end

if stimulus.framegrab==1
    if stimulus.frame.count < size(stimulus.frame.frames,3)
        stimulus.frame.frames(:,:,stimulus.frame.count) = mean(mglFrameGrab(stimulus.frame.coords),3);
        stimulus.frame.count = stimulus.frame.count+1;
    end
end

if stimulus.live.triggerWaiting
    now = mglGetSecs;
    % check eye position, if 
    [pos,time] = mglEyelinkGetCurrentEyePos;
    if ~any(isnan(pos))
        dist = hypot(pos(1),pos(2));
        wasCentered = stimulus.live.centered;
        stimulus.live.centered = dist<2;
        mglBltTexture(mglText(stimulus.live.centered),[5 5]);
        if wasCentered && stimulus.live.centered && stimulus.live.lastTrigger>0
            stimulus.live.triggerTime = stimulus.live.triggerTime + now-stimulus.live.lastTrigger;
        end
        stimulus.live.lastTrigger = now;
    end
    if stimulus.live.triggerTime > 1.5 % not in ms dummy, wait 1.5 seconds (reasonable slow time)
        disp('Eye position centered');
        task = jumpSegment(task);
    end
end

function upGrating(stimulus,task)

% phaseNum = floor(length(stimulus.grating.phases)*rem(mglGetSecs(task.thistrial.trialstart)*stimulus.grating.tf,1)+1);
mglBltTexture(stimulus.tex(stimulus.live.phaseNum),[stimulus.live.x stimulus.live.y 10],0,0,stimulus.live.dir*180/pi+90); % we add 90 so that it's aligned with the motion

function upFix(stimulus)
%%
% for this experiment use a circle to indicate where participants can
% fixate inside of (rather than a cross which might arbitrarily enforce
% poisitioning
% mglGluAnnulus(0,0,1.5,1.55,stimulus.live.fixColor,64);
mglFixationCross(1.5,1.5,stimulus.live.fixColor);

function stimulus = upDots(stimulus,myscreen)

stimulus.dots{stimulus.dot} = updateDotsRadial(stimulus.dots{stimulus.dot},stimulus.live.coherence,myscreen,true);

mglPoints2(stimulus.dots{stimulus.dot}.x+stimulus.live.x,stimulus.dots{stimulus.dot}.y+stimulus.live.y,...
    stimulus.dots{stimulus.dot}.dotsize,stimulus.live.dotColor);

% function stimulus = upDotsInc(stimulus,myscreen)
% stimulus.idots = updateDotsRadial(stimulus.idots,0,myscreen,true);
% mglPoints2(stimulus.idots.x,stimulus.idots.y,...
%     stimulus.idots.dotsize,stimulus.live.dotColor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus

if stimulus.dead, return; end
responseText = {'Incorrect','Correct'};
sideText = {'Left','Right'};
fixColors = {stimulus.colors.red,stimulus.colors.green};
    
responses = [1 0 2];
if any(task.thistrial.whichButton == stimulus.responseKeys)
    if task.thistrial.gotResponse == 0
        task.thistrial.correct = task.thistrial.whichButton == stimulus.responseKeys(responses(task.thistrial.rotation+2));
        disp(sprintf('Subject pressed %i: %s %s',task.thistrial.whichButton,sideText{task.thistrial.whichButton},responseText{task.thistrial.correct+1}));
        stimulus.live.fixColor = fixColors{task.thistrial.correct+1};
        stimulus.staircase{stimulus.feature} = doStaircase('update',stimulus.staircase{stimulus.feature},task.thistrial.correct);
        stimulus.live.resp = 1;
        stimulus.live.dots = 0;
        stimulus.live.fix=1;
    else
        disp(sprintf('(berlin) Subject responded multiple times: %i',task.thistrial.gotResponse+1));
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
        'initialThreshold',0.3,... % radians of rotation (
        'initialStepsize',0.1,...
        'minThreshold=0.001','maxThreshold=1','stepRule','pest',...
        'nTrials=65','maxStepsize=0.1','minStepsize=0.001');
stimulus.staircase{2} = stimulus.staircase{1};

function stimulus = resetStair(stimulus)

for i = 1:2
    if doStaircase('stop',stimulus.staircase{i})
        disp('(berlin) Initializing new staircase...');
        stimulus.staircase{i}(end+1) = doStaircase('init',stimulus.staircase{i}(end));
    end
end

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

function stimulus = initGratings(stimulus,myscreen)

stimulus.maxIndex = 255;
disppercent(-inf,'Creating grating textures');

nContrasts = length(stimulus.contrast);
nPhases = length(stimulus.grating.phases);

gaussianWin = mglMakeGaussian(stimulus.grating.width,stimulus.grating.height,stimulus.grating.sdx,stimulus.grating.sdy);
if strcmp(stimulus.grating.windowType,'gabor')
  % a gaussian window
  win = stimulus.maxIndex-stimulus.maxIndex*gaussianWin;
else
  % a simple window
  win = stimulus.maxIndex-stimulus.maxIndex*(gaussianWin>exp(-1/2));
end
mask = ones(size(win,1),size(win,2),4)*myscreen.grayIndex;
mask(:,:,4) = round(win);
stimulus.mask = mglCreateTexture(mask);

% make each one of he called for gratings
for iPhase = 1:nPhases
  for iContrast = 1:nContrasts
    disppercent(calcPercentDone(iPhase,nPhases,iContrast,nContrasts));
    % get the phase and contast
    thisPhase = (stimulus.grating.phase+stimulus.grating.phases(iPhase))*180/pi;
    thisContrast = stimulus.contrast;
    % make the grating
    thisGrating = round(stimulus.maxIndex*((thisContrast*mglMakeGrating(stimulus.grating.width,nan,stimulus.grating.sf,0,thisPhase))+1)/2);
    % create the texture
    stimulus.tex(iContrast,iPhase) = mglCreateTexture(thisGrating);
  end
end
disppercent(inf);
stimulus.randMaskSize = [size(mask,1) size(mask,2)];
stimulus.randMask = mglCreateTexture(floor(stimulus.maxIndex*rand(stimulus.randMaskSize)));


for iAngle = 1:length(stimulus.orientations)
  % get center of patch
  thisAngle = stimulus.orientations(iAngle);
  centerX = stimulus.x + stimulus.grating.radius*cos(pi*thisAngle/180);
  centerY = stimulus.y + stimulus.grating.radius*sin(pi*thisAngle/180);
  % now get top and bottom point of grating
  thisOrientation = thisAngle+90;
  radius = sqrt((stimulus.grating.width/2).^2 +(stimulus.grating.height/2).^2)-0.5;
  topX = centerX + radius*cos(pi*thisOrientation/180);
  topY = centerY + radius*sin(pi*thisOrientation/180);
  thisOrientation = thisOrientation+180;
  bottomX = centerX + radius*cos(pi*thisOrientation/180);
  bottomY = centerY + radius*sin(pi*thisOrientation/180);
  % place points
  stimulus.grating.refPoints.x{iAngle} = [topX bottomX];
  stimulus.grating.refPoints.y{iAngle} = [topY bottomY];
end

stimulus.waitForBacktickText = mglText('Hit backtick (`) key to start');
