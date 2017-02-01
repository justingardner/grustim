function [ myscreen ] = freedman_dms( varargin )
%FREEDMAN_DMS 
%
% Delayed match to sample with motion 
%
% Based on http://www.nature.com/neuro/journal/v19/n1/full/nn.4168.html#methods
% 
% Manual release of spacebar to indicate same dir / same category
%
% Dots move at 12/s, stimuli either 9 deg or 6 deg, at 7-10 deg
% eccentricity. DMC task was always 7 deg ecc along horizontal axis.
%
% DMS uses 8 equal spaced motions, and match/non-match run on a staircase
%
% Timing is taken from: http://www.jneurosci.org/content/33/32/13157.long
%
% Gaze fixation (within 2.5 degs of fixation) for 500 ms, sample for 650 ms
% 1000 ms delay, then 650 ms test stimulus
%
% No second response required (since no physiological recordings are made)
%
%
% WARNING
% This is behavioral code only!!

global stimulus

%% Initialize Variables

% add arguments later
scan = 0;
plots = 0;
noeye = 0;
getArgs(varargin,{'scan=0','plots=0','noeye=0'});
stimulus.scan = scan;
stimulus.plots = plots;
stimulus.noeye = noeye;
clear localizer invisible scan category noeye task

stimulus.counter = 1; % This keeps track of what "run" we are on.

if stimulus.scan
    warning('Not setup for scanning');
end
%% Setup Screen

if stimulus.scan
    myscreen = initScreen('fMRIprojFlex');
else
    myscreen = initScreen('VPixx');
end

myscreen.background = 0;

%% Open Old Stimfile
stimulus.initStair = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/freedman_dms/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/freedman_dms/%s/1*mat',mglGetSID));

    if length(files) >= 1
        fname = files(end).name;
        
        s = load(sprintf('~/data/freedman_dms/%s/%s',mglGetSID,fname));
        % copy staircases and run numbers
        stimulus.counter = s.stimulus.counter + 1;
        stimulus.staircases = s.stimulus.staircases;

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

disp(sprintf('(berlin) Initializing staircase'));
stimulus = initStaircase(stimulus);

%% Initialize Stimulus
[myscreen] = initStimulus('stimulus',myscreen);

if stimulus.scan
    stimulus.responseKeys = [2 1]; % corresponds to NOMATCH, MATCH
else
    stimulus.responseKeys = [2 1]; % corresponds to  NOMATCH, MATCH
end

stimulus.colors.black = [0 0 0];
stimulus.colors.white = [1 1 1];
stimulus.colors.green = [0 1 0];
stimulus.colors.red = [1 0 0];

stimulus.motion.minecc = 6;
stimulus.motion.maxecc = 10;
stimulus.motion.minrad = 3;
stimulus.motion.maxrad = 4.5;

stimulus.motion.ecc = 8;
stimulus.motion.pos = rand*2*pi; % randomize position
stimulus.motion.rad = 3.75;

stimulus.motion.x = stimulus.motion.ecc * cos(stimulus.motion.pos);
stimulus.motion.y = stimulus.motion.ecc * sin(stimulus.motion.pos);

stimulus.direction.opts = 0:pi/4:(2*pi-pi/4);

%% Dots

stimulus.dots.dotsize = 3;
stimulus.dots.density = 21;
stimulus.dots.speed = 12;
stimulus.dots.centerOffset = 1;

stimulus.dots.minX = -stimulus.motion.rad;
stimulus.dots.maxX = stimulus.motion.rad;
stimulus.dots.minY = stimulus.dots.minX;
stimulus.dots.maxY = stimulus.dots.maxX;

stimulus.dots = initDots(stimulus.dots,myscreen);

%% Generate stencils

mglStencilCreateBegin(1);
mglFillOval(stimulus.motion.x,stimulus.motion.y,[stimulus.motion.rad*2 stimulus.motion.rad*2],1);
mglStencilCreateEnd;

%% Setup Task
task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;

stimulus.curTrial = 0;

task{1}{1}.segmin = [0.4 .650 1.000 .650 0.5];
task{1}{1}.segmax = [0.4 .650 1.000 .650 1.5];

stimulus.seg.ITI1 = 1;
stimulus.seg.stim1 = 2;
stimulus.seg.delay = 3;
stimulus.seg.stim2 = 4;
stimulus.seg.ITI2 = 5;

task{1}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
task{1}{1}.getResponse = [0 0 0 1 1];
task{1}{1}.numTrials = 60;
task{1}{1}.random = 1;
task{1}{1}.parameter.dir1 = stimulus.direction.opts; % which orientation to use (0 deg or 135 deg)
task{1}{1}.parameter.match = [0 1];

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
task{1}{1}.randVars.calculated.dir2 = nan; % will be 0->2*pi
task{1}{1}.randVars.calculated.pos = nan; % will be 0->2*pi

%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);
end
%% EYE CALIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

%% Get Ready...
% clear screen    
mglWaitSecs(1);
mglClearScreen(0);
mglFixationCross(0.2,0.2,stimulus.colors.white);
% if stimulus.scan        
%     mglTextDraw('DO NOT MOVE',[0 1.5]);
% end
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

% save staircases
if ~isfield(stimulus,'staircases')
    stimulus.staircases = {};
end
stimulus.staircases{end+1} = stimulus.staircase;

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

[rotation, stimulus.staircase] = doStaircase('testValue',stimulus.staircase);

% set whether this trial matches or not
rot = [-1 1];
if task.thistrial.match
    task.thistrial.dir2 = task.thistrial.dir1;
else
    task.thistrial.dir2 = task.thistrial.dir1+rot(randi(2))*rotation;
end

disp(sprintf('Trial %i Dir1: %i Dir2: %i',stimulus.curTrial,round(180/pi*task.thistrial.dir1),round(180/pi*task.thistrial.dir2)));

stimulus.live.eyeCount = 0;
stimulus.dead = 0;

task.thistrial.pos = stimulus.motion.pos;

stimulus.live.fixColor = stimulus.colors.white;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)
%%

global stimulus

stimulus.live.triggerWaiting = 0;
if any(task.thistrial.thisseg==[stimulus.seg.ITI1])
    stimulus.live.triggerWaiting = 1;
    stimulus.live.centered = 0;
    stimulus.live.triggerTime = 0;
    stimulus.live.lastTrigger = -1;
end

stimulus.live.stim = 0;
stimulus.live.resp = 0;
stimulus.live.coherence = 1;
stimulus.live.x = stimulus.motion.x;
stimulus.live.y = stimulus.motion.y;
stimulus.live.fix = 1;

stimulus.live.dotColor = stimulus.colors.white;

% use only for constant mode
% if task.thistrial.thisphase==1
%     return
% end

if task.thistrial.thisseg==stimulus.seg.stim1
    stimulus.live.stim=1;
    stimulus.dots.dir = task.thistrial.dir1; % for gratings
elseif task.thistrial.thisseg==stimulus.seg.stim2
    stimulus.live.stim=1;
    stimulus.dots.dir = task.thistrial.dir2; % for gratings
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
if ~stimulus.noeye
    [pos,~] = mglEyelinkGetCurrentEyePos;
    dist = hypot(pos(1),pos(2));
end

if ~stimulus.noeye && task.thistrial.thisseg~=stimulus.seg.ITI2 && ~stimulus.scan
    if ~any(isnan(pos))
        if dist > 2 && stimulus.live.eyeCount > 30
            mglTextSet([],32,stimulus.colors.red);
            disp('Eye movement detected!!!!');
            mglTextDraw('Eye Movement Detected',[0 0]);
            mglFlush
            myscreen.flushMode = 1;
            stimulus.dead = 1;
            return
        elseif dist > 2
            stimulus.live.eyeCount = stimulus.live.eyeCount + 1;
        end
    end
end

if stimulus.live.stim
    mglStencilSelect(1);
    stimulus = upDots(stimulus,myscreen);

    mglStencilSelect(0);
end
if stimulus.live.fix || stimulus.live.resp==1
    upFix(stimulus);
end

% if stimulus.framegrab==1
%     if stimulus.frame.count < size(stimulus.frame.frames,3)
%         stimulus.frame.frames(:,:,stimulus.frame.count) = mean(mglFrameGrab(stimulus.frame.coords),3);
%         stimulus.frame.count = stimulus.frame.count+1;
%     end
% end

if stimulus.live.triggerWaiting
    now = mglGetSecs;
    % check eye position, if 
    if ~any(isnan(pos))
        wasCentered = stimulus.live.centered;
        stimulus.live.centered = dist<2;
        if wasCentered && stimulus.live.centered && stimulus.live.lastTrigger>0
            stimulus.live.triggerTime = stimulus.live.triggerTime + now-stimulus.live.lastTrigger;
        end
        stimulus.live.lastTrigger = now;
    end
    if stimulus.live.triggerTime > 0.5 % not in ms dummy, wait 1.5 seconds (reasonable slow time)
        disp('Eye position centered');
        task = jumpSegment(task);
    end
end

function upFix(stimulus)
%%
% for this experiment use a circle to indicate where participants can
% fixate inside of (rather than a cross which might arbitrarily enforce
% poisitioning
% mglGluAnnulus(0,0,1.5,1.55,stimulus.live.fixColor,64);
mglFixationCross(1,1,stimulus.live.fixColor);

function stimulus = upDots(stimulus,myscreen)

stimulus.dots = updateDots(stimulus.dots,stimulus.live.coherence,myscreen,true);

mglPoints2(stimulus.dots.x+stimulus.live.x,stimulus.dots.y+stimulus.live.y,...
    stimulus.dots.dotsize,stimulus.live.dotColor);

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
    
if any(task.thistrial.whichButton == stimulus.responseKeys)
    if task.thistrial.gotResponse == 0
        task.thistrial.correct = task.thistrial.whichButton == stimulus.responseKeys(task.thistrial.match+1);
        disp(sprintf('Subject pressed %i: %s %s',task.thistrial.whichButton,sideText{task.thistrial.whichButton},responseText{task.thistrial.correct+1}));
        stimulus.live.fixColor = fixColors{task.thistrial.correct+1};
        stimulus.staircase = doStaircase('update',stimulus.staircase,task.thistrial.correct);
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
stimulus.staircase = doStaircase('init','upDown',...
        'initialThreshold',0.25,... % radians of rotation (
        'initialStepsize',0.05,...
        'minThreshold=0.001','maxThreshold=1.57','stepRule','pest',...
        'nTrials=50','maxStepsize=0.1','minStepsize=0.001');

%%%%%%%%%%%%%%%%%%%%%%%
%    dispInfo    %
%%%%%%%%%%%%%%%%%%%%%%%
function dispInfo(stimulus)
%%
disp('not implemented');
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create dots for optic flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initDots(dots,~)

% maximum depth of points
dots.dir = 0;

area = (dots.maxX-dots.minX)*(dots.maxY-dots.minY);

dots.n = round(area * dots.density);

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
function dots = updateDots(dots,coherence,myscreen,repick)

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
mask = ones(size(win,1),size(win,2),4)*1;%myscreen.grayIndex;
win = win / max(win(:)) * 255;
mask(:,:,4) = round(win);
stimulus.mask = mglCreateTexture(mask);

% make each one of he called for gratings
for iPhase = 1:nPhases
  for iContrast = 1:nContrasts
    disppercent(calcPercentDone(iPhase,nPhases,iContrast,nContrasts));
    % get the phase and contast
    thisPhase = (stimulus.grating.phase+stimulus.grating.phases(iPhase))*180/pi;
    thisContrast = stimulus.contrastOverride;
    % make the grating
    thisGrating = round(stimulus.maxIndex*((thisContrast*mglMakeGrating(stimulus.grating.width,nan,stimulus.grating.sf,0,thisPhase))+1)/2);
    % create the texture
%     thisGrating = repmat(thisGrating,453,1);
    
    stimulus.tex(iContrast,iPhase) = mglCreateTexture(thisGrating);
  end
end
disppercent(inf);
% stimulus.randMaskSize = [size(mask,1) size(mask,2)];
% stimulus.randMask = mglCreateTexture(floor(stimulus.maxIndex*rand(stimulus.randMaskSize)));


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
