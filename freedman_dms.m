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
diff = 0;
getArgs(varargin,{'scan=0','plots=0','noeye=0','diff=0'});
stimulus.scan = scan;
stimulus.plots = plots;
stimulus.noeye = noeye;
stimulus.diff = diff/360;
clear localizer invisible scan category noeye task diff

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
        stimulus.nonmatchOpts = s.stimulus.nonmatchOpts;

        clear s;
        disp(sprintf('(freedman) Data file: %s loaded.',fname));
        
    end
end

if isfield(stimulus,'nonmatchOpts') && stimulus.counter>length(stimulus.nonmatchOpts)
    stimulus.nonmatchOpts = [stimulus.nonmatchOpts stimulus.nonmatchOpts];
end

if stimulus.plots==2
    dispInfo(stimulus);
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% init staircase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disp(sprintf('(freedman) Initializing staircase'));
% stimulus = initStaircase(stimulus);

if ~isfield(stimulus,'nonmatchOpts')
    disp('WARNING: New non-match options set');
    opts = [2,4,8,16,32];
    stimulus.nonmatchOpts = opts(randperm(length(opts)));
end

disp(sprintf('(freedman) This is run %i/%i',stimulus.counter,length(stimulus.nonmatchOpts)));

%% Initialize Stimulus
[myscreen] = initStimulus('stimulus',myscreen);

if stimulus.scan
    stimulus.responseKeys = [50];
else
    stimulus.responseKeys = [50];
end

stimulus.colors.black = [0 0 0];
stimulus.colors.white = [1 1 1];
stimulus.colors.green = [0 1 0];
stimulus.colors.red = [1 0 0];

% stimulus.motion.minecc = 6;
% stimulus.motion.maxecc = 10;
% stimulus.motion.minrad = 3;
% stimulus.motion.maxrad = 4.5;

stimulus.motion.ecc = 8;
stimulus.motion.pos = rand*2*pi; % randomize position
stimulus.motion.rad = 3.75; % 7.5 degree stimulus diameter

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

task{1}{1}.segmin = [inf .650 1.000 .650 1 inf];
task{1}{1}.segmax = [inf .650 1.000 .650 1 inf];

stimulus.seg.ITI1 = 1;
stimulus.seg.stim1 = 2;
stimulus.seg.delay = 3;
stimulus.seg.stim2 = 4;
stimulus.seg.resp = 5;
stimulus.seg.ITI2 = 6;

task{1}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
task{1}{1}.getResponse = [0 0 0 0 0];
task{1}{1}.numTrials = 45;
task{1}{1}.random = 1;
task{1}{1}.parameter.dir1 = stimulus.direction.opts; % which orientation to use (0 deg or 135 deg)
task{1}{1}.parameter.match = [0 1];

if stimulus.scan
    task{1}{1}.synchToVol(stimulus.seg.ITI) = 1;
end

%% Tracking

% these are variables that we want to track for later analysis.
task{1}{1}.randVars.calculated.correct = nan;
task{1}{1}.randVars.calculated.dir2 = nan; % will be 0->2*pi
task{1}{1}.randVars.calculated.pos = nan; % will be 0->2*pi
task{1}{1}.randVars.calculated.nomatchResp = nan; % will be 0->2*pi

%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,[],@startTrialCallback,[],[]);
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
disp(sprintf('(freedman) Starting run number: %i.',stimulus.counter));
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
    disp('(freedman) Displaying plots');
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


% set whether this trial matches or not
rot = [-1 1];
% [rotation, stimulus.staircase] = doStaircase('testValue',stimulus.staircase);
rotation = stimulus.nonmatchOpts(stimulus.counter)*pi/180;

if stimulus.diff
    rotation = stimulus.diff;
end

if task.thistrial.match
    task.thistrial.dir2 = task.thistrial.dir1;
else
    task.thistrial.dir2 = task.thistrial.dir1+rot(randi(2))*rotation;
end

disp(sprintf('(freedman) Trial %i Dir1: %i Dir2: %i',stimulus.curTrial,round(180/pi*task.thistrial.dir1),round(180/pi*task.thistrial.dir2)));

stimulus.live.eyeCount = 0;
stimulus.dead = 0;
stimulus.live.responded = 0;

task.thistrial.pos = stimulus.motion.pos;

stimulus.live.fixColor = stimulus.colors.white;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)
%%

global stimulus
stimulus.live.fixColor = stimulus.colors.white;

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

stimulus.live.eyeDead =0 ;
stimulus.live.barDead =0 ;

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
elseif task.thistrial.thisseg==stimulus.seg.ITI2
    stimulus.live.fix =0 ;
    if task.thistrial.gotResponse==0 && ~stimulus.dead
        % if we didn't respond yet, check and see if we are right
        task.thistrial.gotResponse=1;
        task.thistrial.nomatchResp = 1; 
        task.thistrial.correct = ~task.thistrial.match;
        doResponse(task);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus

mglClearScreen(0);

if stimulus.dead && mglGetSecs(task.thistrial.segStartSeconds)>1
    jumpSegment(task,inf); stimulus.dead=0;
end

if stimulus.dead
    if stimulus.live.barDead
        mglTextSet([],32,stimulus.colors.red);
        mglTextDraw('Space bar lifted',[0 0]);
    end
    if stimulus.dead && stimulus.live.eyeDead
        mglTextSet([],32,stimulus.colors.red);
        mglTextDraw('Eye Movement Detected',[0 0]);
    end
    return
end

% check eye pos
if ~stimulus.noeye
    [pos,~] = mglEyelinkGetCurrentEyePos;
    dist = hypot(pos(1),pos(2));
end

if ~stimulus.noeye && ~any(task.thistrial.thisseg==[stimulus.seg.ITI1 stimulus.seg.ITI2 stimulus.seg.resp]) && ~stimulus.scan
    if ~any(isnan(pos))
        if dist > 3 && stimulus.live.eyeCount > 30
            disp('Eye movement detected!!!!');
            stimulus.dead = 1;
            stimulus.live.eyeDead=1;
            return
        elseif dist > 3
            stimulus.live.eyeCount = stimulus.live.eyeCount + 1;
        end
    end
end

if ~any(task.thistrial.thisseg==[stimulus.seg.ITI1 stimulus.seg.ITI2 stimulus.seg.resp])
    if ~stimulus.live.spaceDown
        disp('Space bar lift detected!!!!');
        stimulus.dead =1;
        stimulus.live.barDead=1;
        return
    end
end

stimulus.live.spaceDown = logical(mglGetKeys(50));

if task.thistrial.thisseg==stimulus.seg.ITI1
    if stimulus.noeye && stimulus.live.spaceDown
        disp('(freedman) Starting trial--space detected');
        task = jumpSegment(task);
    end
end

if task.thistrial.thisseg==stimulus.seg.ITI2
    if ~stimulus.live.spaceDown && mglGetSecs(task.thistrial.segStartSeconds)>0.5
        disp('(freedman) Ending trial--space lifted');
        task = jumpSegment(task);
    end
end


if (task.thistrial.thisseg==stimulus.seg.resp) && ~stimulus.live.spaceDown && ~stimulus.live.responded
    disp('(freedman) Match response detected');
    stimulus.live.responded = 1;
    [task,myscreen] = customResponseCallback(task,myscreen);
end

if stimulus.live.stim
    mglStencilSelect(1);
    stimulus = upDots(stimulus,myscreen);

    mglStencilSelect(0);
end
if stimulus.live.fix || stimulus.live.resp==1
    upFix(stimulus);
end

if ~stimulus.noeye && stimulus.live.triggerWaiting
    now = mglGetSecs;
    % check eye position, if 
    if ~any(isnan(pos))
        wasCentered = stimulus.live.centered;
        stimulus.live.centered = dist<3;
        if wasCentered && stimulus.live.centered && stimulus.live.lastTrigger>0
            stimulus.live.triggerTime = stimulus.live.triggerTime + now-stimulus.live.lastTrigger;
        end
        stimulus.live.lastTrigger = now;
    end
    if stimulus.live.spaceDown && stimulus.live.triggerTime > 0.5 % not in ms dummy, wait 1.5 seconds (reasonable slow time)
        disp('Starting trial--eye centered and space pressed.');
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

%%%%%
function [task, myscreen] = customResponseCallback(task,myscreen)

if task.thistrial.gotResponse == 0
    task.thistrial.gotResponse = task.thistrial.gotResponse+1;
    
    task.thistrial.nomatchResp = 0; 
    task.thistrial.correct = task.thistrial.match; % if match, then we are correct
    
    doResponse(task);
else
    disp(sprintf('(freedman) Subject responded multiple times: %i',task.thistrial.gotResponse+1));
end

function doResponse(task)
global stimulus 

if stimulus.dead, return; end

responseText = {'Incorrect','Correct'};
fixColors = {stimulus.colors.red,stimulus.colors.green};
matchText = {'Match','Non-Match'};

disp(sprintf('Subject responded: %s, %s',matchText{task.thistrial.nomatchResp+1},responseText{task.thistrial.correct+1}));
stimulus.live.fixColor = fixColors{task.thistrial.correct+1};
% if ~task.thistrial.match
%     stimulus.staircase = doStaircase('update',stimulus.staircase,task.thistrial.nomatchResp);
% else
%     stimulus.staircase = doStaircase('update',stimulus.staircase,task.thistrial.nomatchResp,0);
% end
stimulus.live.resp = 1;
stimulus.live.dots = 0;
stimulus.live.fix = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              HELPER FUNCTIONS                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
%    initStaircase     %
%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initStaircase(stimulus)
%%
if ~isfield(stimulus,'staircase')
    stimulus.staircase = doStaircase('init','fixed','fixedVals',[3 6 9 12 24 36]*pi/180);
end

%%%%%%%%%%%%%%%%%%%%%%%
%    dispInfo    %
%%%%%%%%%%%%%%%%%%%%%%%
function dispInfo(stimulus)
%%

% get the files list
files = dir(fullfile(sprintf('~/data/freedman_dms/%s/*.mat',mglGetSID)));

% load the files and pull out the data (long form)
%    1          2               3           4       5        6
%  run #    local trial     real trial   match     dir1     dir2 
%    7           8             9
%   pos      nomatchResp     correct
count = 1; data = zeros(10000,9);

for fi = 1:length(files)
    load(fullfile(sprintf('~/data/freedman_dms/%s/%s',mglGetSID,files(fi).name)));
    
    e = getTaskParameters(myscreen,task);
    e = e{1}; % why?!
    
    data(count:count+(e.nTrials-1),:) = [repmat(stimulus.counter,e.nTrials,1) (1:e.nTrials)' (count:count+(e.nTrials-1))' ...
        e.parameter.match' e.parameter.dir1' e.randVars.dir2' ...
        e.randVars.pos' e.randVars.nomatchResp' e.randVars.correct'];
    
    count = count+e.nTrials;
end

data = data(1:(count-1),:);

data = data(~isnan(data(:,9)),:);

data(:,end+1) = abs(data(:,5)-data(:,6)); % diff
data(:,10) = round(data(:,10)*180/pi);

udegs = unique(data(:,10));
udegs = udegs(udegs~=0);
%%
runs = unique(data(:,1));
degs = zeros(1,length(udegs));
dprime = zeros(1,length(udegs));
crit = dprime;
for ui = 1:length(udegs) % 5 runs?
    dat__ = sel(data,10,udegs(ui)); % find data with this, get all unique runs
    cruns = unique(dat__(:,1));
    
    dat_ = [];
    for ci = 1:length(cruns)
        dat_ = [dat_;sel(data,1,cruns(ci))];
    end
    
    % norminv(hits) - norminv(false alarms)
    
    % present = non-match
    dat_nm = sel(dat_,4,0);
    dat_m = sel(dat_,4,1);
    
    nm_ = mean(dat_nm(:,8));
    m_ = mean(dat_m(:,8));
    
    if nm_==1, nm_=1-eps; end
    if m_==0, m_=eps; end
    
    dprime(ui) = norminv(nm_) - norminv(m_);
    crit(ui) = -0.5 * (norminv(nm_)+norminv(m_));
    degs(ui) = max(dat_(:,10));
end

% dprime(dprime==Inf) = max(dprime(dprime~=Inf));

%%

h = figure; 
subplot(211);
hold on

plot(degs,dprime,'o','MarkerFaceColor','black','MarkerEdgeColor','white');
hline(1,'--r');

a = axis;
axis([0 35 0 max(a(4),3)]);

set(gca,'YTick',[0 1 2 3]);
set(gca,'XTick',[2 4 8 16 32]);

xlabel('Angle difference (degs)');
ylabel('d''');

drawPublishAxis

subplot(212);
hold on
plot(degs,crit,'o','MarkerFaceColor','black','MarkerEdgeColor','white');
hline(0,'--r');

axis([0 35 -1.1 1.1]);
set(gca,'XTick',[2 4 8 16 32]);

title('Criterion');
xlabel('Angle difference (degs)');
% ylabel('Criterion');

set(gca,'YTick',[-1 0 1],'YTickLabel',{'Always non-match','Unbiased','Always match'});

drawPublishAxis

% savepdf(h,fullfile(sprintf('~/data/freedman_dms/%s/%s_dprime.pdf',mglGetSID,mglGetSID)));

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

