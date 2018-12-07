function [ myscreen ] = afcon( varargin )
% AFCOM (attention field continuous)
% *** set 'noeye=1' to turn of the eye tracker***
%
%   This is code for a psychophysics experiment using continuous response
%   with feature-based selection (attention?). The task is a free-viewing
%   task 

%%

global stimulus fixStimulus

stimulus = struct;
fixStimulus = struct;

stimulus.powerwheelFactor = 90;

%% Initialize Variables

% add arguments later
scan = 0;
plots = 0;
noeye = 0;
debug = 0;
powerwheel = 0;
run = 0; 
eyewindow=0; 
mouse=0; 
practice=0; 
practiceType=-1;
cue=0;

getArgs(varargin,{'scan=0','cue=1','plots=0','noeye=0','powerwheel=1','eyewindow=2.5','practice=0','practiceType=-1','debug=0','run=0','build=0','mouse=0'});
stimulus.scan = scan;
stimulus.plots = plots;
stimulus.noeye = noeye;
stimulus.cue = cue; % cue = 1 means direction, cue = 2 means color
stimulus.practice = practice;
stimulus.practiceType = practiceType;
stimulus.mousedebug = mouse;
stimulus.powerwheel = powerwheel;
stimulus.eyewindow = eyewindow;
stimulus.debug = debug;
stimulus.overrideRun = run;

clear localizer invisible scan noeye task test2 attend build eyewindow mouse practice powerwheel cue

if stimulus.scan
    warning('Disabling eyewindow');
    stimulus.eyewindow=0;
end

%% Open Old Stimfile
stimulus.counter = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/afcon/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/afcon/%s/1*mat',mglGetSID));

    if length(files) >= 1
        fname = files(end).name;

        s = load(sprintf('~/data/afcon/%s/%s',mglGetSID,fname));
        % copy staircases and run numbers
        stimulus.counter = s.stimulus.counter + 1;
        stimulus.colors = s.stimulus.colors;
        stimulus.colorwheel = s.stimulus.colorwheel;
        clear s;
        disp(sprintf('(afcon) Data file: %s loaded.',fname));
    else
        warning('(afcon) Unable to load previous data files. If this is *not* the first run there is something wrong.');
    end
end

%% Display run info
disp('*************************');
disp(sprintf('(afcon) This is run #%i',stimulus.counter));
disp('*************************');

%% Setup Screen
myscreen = initScreen('VPixx');
% set background to black
myscreen.background = 0;

%% load the calib
if isfield(myscreen,'calibFullFilename')
    calib = load(fullfile(myscreen.calibFullFilename));
    stimulus.calib = calib.calib;
else
    stimulus.calib = []; % need this so that mglLab2rgb doesn't fail
end

%% Plot and return
if stimulus.plots==2
    dispInfo;
    return
end

%% Initialize Stimulus
myscreen.stimulusNames{1} = 'stimulus';

if stimulus.powerwheel
    stimulus.responseKeys = 1;
else
    stimulus.responseKeys = [1 2]; % left right
end

%% load the calib
if isfield(myscreen,'calibFullFilename')
    calib = load(fullfile(myscreen.calibFullFilename));
    stimulus.calib = calib.calib;
else
    stimulus.calib = []; % need this so that mglLab2rgb doesn't fail
end

%% Colors
if ~isfield(stimulus,'colors')
    stimulus.colors.white = [0.2 0.2 0.2]; stimulus.colors.red = [0.8 0 0];
    stimulus.colors.green = [0 0.8 0]; stimulus.colors.black = [0 0 0];
end

%% Setup patches and stencils

% there will be 12 possible locations where we can show dots. We will use 6
% patches of dots 

stimulus.dotScale = 0.3;
stimulus.cueScale = 3;

%% Create the cue patch
rain.maxX = myscreen.imageWidth;
rain.maxY = myscreen.imageHeight;
rain.density = 0.2;
rain.dotScale = 3;
rain.maxAlive = myscreen.framesPerSecond/4;
rain.speed = 2;
rain.dir = pi*1.5;

D = 40;
ang1 = rand*2*pi;
ang2 = rand*2*pi;

stimulus.relRain = initDots(rain);
stimulus.relRain.color = mglLab2rgb([60 D*cos(ang1) D*sin(ang1)],stimulus.calib);
stimulus.irrelRain = initDots(rain);
stimulus.irrelRain.color = mglLab2rgb([60 D*cos(ang2) D*sin(ang2)],stimulus.calib);

%% Setup Probe Task
task{1}{1} = struct;
task{1}{1}.seglen = 30;

task{1}{1}.waitForBacktick = 1;

task{1}{1}.getResponse = 0;

if stimulus.scan==1
    task{1}{1}.numTrials = Inf;
else
    task{1}{1}.numTrials = 40;
end

task{1}{1}.random = 1;
task{1}{1}.parameter.trialType = 1;
task{1}{1}.parameter.difficulty = 8; % how fast reward drops off -- 4 is quite steep, 8 is quite shallow
task{1}{1}.parameter.raincoherence = 0.5; % how coherent the rain is
task{1}{1}.parameter.delay = 200; % ms delay before rain has an effect

%% Difficulty test code:
% plot the drop-off curve of performance
% x = 0:.01:20;
% y = 0.5*exp(-x/8);
% plot(x,y);
% axis([0 20 0 0.5]);

%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,@endTrialCallback);
end

%% Initial stimulus position (zero)
stimulus.live.pos = 0;

% track the powerwheel/mouse position
stimulus.live.cAngle = getAngle(myscreen);

% success rate
stimulus.live.rate = 1; % / second
stimulus.live.success = {};

% calculate how many frames the dots precede the noise
stimulus.noise.ticks = 1.0 * myscreen.framesPerSecond;

% motion noise
stimulus.live.noises = {'slow','fast'};
stimulus.noise.slow = initNoise(3,10,stimulus.noise.ticks,'slow');
stimulus.noise.fast = initNoise(1,1,stimulus.noise.ticks,'fast');

%% EYE CALIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~stimulus.noeye
    myscreen = eyeCalibDisp(myscreen);
end
% let the user know
disp(sprintf('(afcon) Starting run number: %i.',stimulus.counter));

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
mglClearScreen;
mglTextSet([],32,stimulus.colors.white);
% get count
mglTextDraw('Please wait',[0 0]);
mglFlush
myscreen.flushMode = 1;

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

if stimulus.plots
    disp('(afcon) Displaying plots');
    dispInfo(stimulus);
end

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

function dispInfo()
%%
files = dir(fullfile('~/data/afcon/',mglGetSID,'*.mat'));

for fi = 1:length(files)
    load(fullfile('~/data/afcon/',mglGetSID,files(fi).name));
    exp = getTaskParameters(myscreen,task);
    e{fi} = exp{1};
end

function [task, myscreen] = startTrialCallback(task,myscreen)
%%

function [task, myscreen] = endTrialCallback(task,myscreen)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = startSegmentCallback(task,myscreen)
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function angle = getAngle(myscreen)
global stimulus

if stimulus.powerwheel
    mInfo = mglGetMouse(myscreen.screenNumber);
    angle = -(mInfo.x-myscreen.screenWidth/2)/stimulus.powerwheelFactor;
    
else
    mInfo = mglGetMouse(myscreen.screenNumber);
    degx = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
    degy = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;
    angle = atan2(degy,degx);
end

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus

mglClearScreen();
angle = getAngle(myscreen);
if angle ~= stimulus.live.cAngle
    delta = angdist(angle, stimulus.live.cAngle); %angle - stimulus.live.cAngle;
    if (angle-stimulus.live.cAngle)>0
        stimulus.live.pos = stimulus.live.pos - delta;
    else
        stimulus.live.pos = stimulus.live.pos + delta;
    end
        
    stimulus.live.cAngle = angle;
end

% update the stimulus position with noise
for ni = 1:length(stimulus.live.noises)
    n = stimulus.noise.(stimulus.live.noises{ni});
    n = updateNoise(n,myscreen);
    % use the past stimulus position
    stimulus.live.pos = stimulus.live.pos + n.vals(end) / myscreen.framesPerSecond;
    stimulus.noise.(stimulus.live.noises{ni}) = n;
end

% update dots
stimulus.relRain.dir = stimulus.noise.slow.dir;
stimulus.relRain.speed = stimulus.noise.slow.vals(1);

if rand < 0.001
    stimulus.irrelRain.dir = 0;
elseif rand < 0.002
    stimulus.irrelRain.dir = pi;
end

stimulus.relRain = updateDots(stimulus.relRain,task.thistrial.raincoherence,1);
drawDots(stimulus.relRain.x,stimulus.relRain.y,0.1,stimulus.relRain.color,myscreen);
% stimulus.irrelRain = updateDots(stimulus.irrelRain,task.thistrial.raincoherence,1);
% drawDots(stimulus.irrelRain.x,stimulus.irrelRain.y,0.1,stimulus.irrelRain.color,myscreen);

if stimulus.live.pos<(-myscreen.imageWidth/2)
    stimulus.live.pos = -myscreen.imageWidth/2;
end
if stimulus.live.pos>(myscreen.imageWidth/2)
    stimulus.live.pos = myscreen.imageWidth/2;
end

% draw a vertical line at zero
mglLines2(0,-5,0,5,1,.3*[1 1 1]);
mglLines2(-.25,0,.25,0,1,.3*[1 1 1]);

% draw a circle at the stimulus pos (in degrees)
mglGluDisk(stimulus.live.pos,0,[1 0 0]);

% set rate 
stimulus.live.rate = 0.5*exp(-abs(stimulus.live.pos)/task.thistrial.difficulty)+0.025; %0.5 - (0.5*abs(stimulus.live.pos)/20);

if rand < (stimulus.live.rate/myscreen.framesPerSecond)
    stimulus.live.success{end+1} = newSuccess();
end

% deal with points
kill = zeros(size(stimulus.live.success));
for si = 1:length(stimulus.live.success)
    [stimulus.live.success{si},kill(si)] = updateSuccess(stimulus.live.success{si});
end
stimulus.live.success = stimulus.live.success(~kill);


%     task.thistrial.respAngle = mod(task.thistrial.respAngle,2*pi);
    
%     stimulus.data.mouseTrack(task.trialnum,stimulus.data.mouseTick) = task.thistrial.respAngle;
%     stimulus.data.mouseTick = stimulus.data.mouseTick + 1;

function drawDots(x,y,scale,c,msc)
cFlag = size(c,1)==1;
% draw the dots one at a time with mglGluDisk
for di = 1:length(x)
    if cFlag
        mglGluDisk(x(di)-msc.imageWidth/2,y(di)-msc.imageHeight/2,scale,c);
    else
        mglGluDisk(x(di)-msc.imageWidth/2,y(di)-msc.imageHeight/2,scale,c(di,:));
    end
end

function s = newSuccess()
s = struct;
s.tick = 0;

function [s,kill] = updateSuccess(s)
global stimulus

s.tick = s.tick+1;
kill = s.tick>25;
mglGluAnnulus(stimulus.live.pos,0,s.tick/10,(s.tick+1)/10,(25-s.tick)/25*[1 1 1]);

% noise functions
function n = initNoise(hz,power,ticks,type)
n = struct;
n.hz = hz;
n.power = power;
n.val = 0;
n.dir = 0;
n.vals = zeros(1,ticks);
n.type = type;

function [n] = updateNoise(n,myscreen)

if rand < (1/n.hz/myscreen.framesPerSecond)
    % reset noise
    nDir = randsample([1 -1],1);
    opts = [0 nan pi];
    n.dir = opts(nDir+2);
    n.val = rand*n.power*nDir;
    disp(sprintf('Updating %s noise, new power: %1.2f',n.type,n.val));
end

n.vals = shift(n.vals,1);
n.vals(1) = n.val; % current value at the start, past values later

function [task, myscreen] = getResponseCallback(task, myscreen)

if task.thistrial.dead, return; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helpers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rgb = ang2rgb(ang)

global stimulus

a = stimulus.colorwheel.distanceLab*cos(ang)+stimulus.colorwheel.acenter;
b = stimulus.colorwheel.distanceLab*sin(ang)+stimulus.colorwheel.bcenter;

rgb = mglLab2rgb([stimulus.backgroundLab(1) a b],stimulus.calib);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create dots for horizontal motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initDots(dots)

dots.dir = 0;

area = dots.maxX*dots.maxY;

dots.n = round(area * dots.density);

% make a some points
% dots.n = 500*dots.density;
% make sure it's an even number
dots.n = dots.n + mod(dots.n,2);

% set half to white and half to black
dots.con = repmat([1 2],1,dots.n/2);

dots.x = rand(1,dots.n)*dots.maxX;
dots.y = rand(1,dots.n)*dots.maxY;

% Why replace dots? Because if you don't then peripheral overlapping dot
% patches will rival!!

dots.alive = randi(dots.maxAlive,1,dots.n); % set to random up to 200 ms

dots.xdisp = dots.x;
dots.ydisp = dots.y;

% set incoherent dots to 0
dots.coherency = 1;
dots.incoherent = rand(1,dots.n) > dots.coherency;
dots.incoherentn = sum(dots.incoherent);
dots.coherent = ~dots.incoherent;

% track time
dots.time = mglGetSecs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots for horizontal motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updateDots(dots,coherence,repick)

elapsed = mglGetSecs-dots.time;
dots.time = mglGetSecs;

dots.alive = dots.alive+1;
rIdx = dots.alive>dots.maxAlive;
replace = sum(rIdx);
dots.x(rIdx) = rand(1,replace)*dots.maxX;
dots.y(rIdx) = rand(1,replace)*dots.maxY;
dots.alive(rIdx) = 0;

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

vectorLength = dots.speed*elapsed;

% move coherent dots
dots.x(dots.coherent) = dots.x(dots.coherent) + vectorLength * cos(dots.dir);
dots.y(dots.coherent) = dots.y(dots.coherent) + vectorLength * sin(dots.dir);

dots.x(dots.incoherent) = dots.x(dots.incoherent) + vectorLength * cos(rand(1,dots.incoherentn)*2*pi);
dots.y(dots.incoherent) = dots.y(dots.incoherent) + vectorLength * sin(rand(1,dots.incoherentn)*2*pi);

offscreen = dots.x > dots.maxX;
dots.x(offscreen) = mod(dots.x(offscreen),dots.maxX);
offscreen = dots.x < 0;
dots.x(offscreen) = mod(dots.x(offscreen),dots.maxX);

offscreen = dots.y > dots.maxY;
dots.y(offscreen) = mod(dots.y(offscreen),dots.maxY);
offscreen = dots.y < 0;
dots.y(offscreen) = mod(dots.y(offscreen),dots.maxY);
