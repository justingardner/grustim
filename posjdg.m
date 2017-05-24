function [ myscreen ] = posjdg( varargin )
%
% POSITION JUDGMENTS 
%    Position judgment task with three attentional conditions. 
%
% Usage: posjdg(varargin)
% Authors: Dan Birman, Akshay Jagadeesh
% Date: 05/10/2017
%
% Parameters:
%
%   To run the endogenous attention experiment
%      posjdg('att=1') --> default
%
%   To run the exogenous attention experiment
%      posjdg('att=2')
%     
%   To run without eyetracking
%      posjdg('noeye=1') --> default 'noeye=0'
%   
%   To run without powerwheel
%      posjdg('powerwheel=0') --> default 'powerwheel=1'
%
%   To run in scanner
%      posjdg('scan=1') --> default 'scan=0'
%
%   To generate plots
%      posjdg('plots=1') --> default 'plots=0'

global stimulus

stimulus = struct;

%% Initialize Variables

% add arguments later
scan = 0;
plots = 0;
noeye = 0;
debug = 0;
noimp = 0;
training = 0;
powerwheel = 0;
test2 = 0;
att = 1;
getArgs(varargin,{'att=1', 'scan=0','plots=0','noeye=0','debug=0','training=0','powerwheel=1','test2=0'});
stimulus.att = att;
stimulus.scan = scan;
stimulus.powerwheel = powerwheel;
stimulus.plots = plots;
stimulus.test2 = test2;
stimulus.noeye = noeye;
stimulus.debug = debug;
stimulus.noimp = noimp;
stimulus.training = training;
clear localizer invisible scan noeye task test2

if stimulus.scan
    warning('Not setup for scanning');
end

%% Stimulus parameters 

% new prior for each run
if stimulus.att == 1
    stimulus.condition = 'Endo';
    stimulus.prior = rand*2*pi;
    while (abs(mod(stimulus.prior,pi/2)<0.15))
      stimulus.prior = rand*2*pi;
    end
    disp(sprintf('Prior chosen to be %02.0f, dist was: %02.2f',stimulus.prior*180/pi,abs(mod(stimulus.prior,pi/2))));
elseif stimulus.att == 2
    stimulus.condition = 'Exo';
elseif stimulus.att == 3
    stimulus.condition = 'Sacc';
end
stimulus.sd = pi/8;
stimulus.ecc = 6;

%% Open Old Stimfile
stimulus.counter = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/posjdg_%s/%s',stimulus.condition,mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/posjdg_%s/%s/1*mat',stimulus.condition,mglGetSID));

    if length(files) >= 1
        fname = files(end).name;
        
        s = load(sprintf('~/data/posjdg_%s/%s/%s',stimulus.condition,mglGetSID,fname));
        % copy staircases and run numbers
        stimulus.counter = s.stimulus.counter + 1;
        stimulus.staircase = s.stimulus.staircase;
        clear s;
        disp(sprintf('(posjdg) Data file: %s loaded.',fname));
    end
end
disp(sprintf('(posjdg) This is run #%i',stimulus.counter));

%% Setup Screen
myscreen = initScreen('VPixx');

% set background to grey
myscreen.background = 0.5;

%% Staircase
if ~isfield(stimulus,'staircase')
    disp('(posjdg) WARNING: New staircase');
    initStair();
else
    resetStair();
end

%% Setup missing initial variables

if ~isfield(stimulus,'counter')
    stimulus.counter = 1; % This keeps track of what "run" we are on.
end

%% Plot and return
if stimulus.plots==2
    dispInfo(stimulus);
    return
end

%% Initialize Stimulus

myscreen = initStimulus('stimulus',myscreen);

localInitStimulus();
    
if stimulus.powerwheel
    stimulus.responseKeys = 1; % detection key (when powerwheel not in use)
else
    stimulus.responseKeys = [1 2 3]; % 
end

%% Colors
initGammaTable(myscreen);
stimulus.colors.rmed = 127.5;

% We're going to add an equal number of reserved colors to the top and
% bottom, to try to keep the center of the gamma table stable.
stimulus.colors.reservedBottom = [1 0 0; 0 0 0]; % fixation cross colors
stimulus.colors.reservedTop = [1 1 1; 0 1 0]; % correct/incorrect colors
stimulus.colors.black = 1/255; stimulus.colors.white = 254/255;
stimulus.colors.red = 0/255; stimulus.colors.green = 255/255;
stimulus.colors.nReserved = 2; % this is /2 the true number, because it's duplicated
stimulus.colors.nUnreserved = 256-(2*stimulus.colors.nReserved);

stimulus.colors.mrmax = stimulus.colors.nReserved - 1 + stimulus.colors.nUnreserved;
stimulus.colors.mrmin = stimulus.colors.nReserved;

%% Setup Task

%%%%%%%%%%%%% PHASE ONE %%%%%%%%%%%%%%%%%
%%%%% PRIOR + ESTIMATE OF THRESHOLD %%%%%

stimulus.curTrial(1) = 0;

task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;

switch stimulus.att
    case 1 %Endo
        % task waits for fixation on first segment
        task{1}{1}.segmin = [inf 0.000 .200 0.800];
        task{1}{1}.segmax = [inf 1.000 .200 0.800];

        stimulus.seg = {};
        stimulus.seg{1}.ITI1 = 1; % waits for user input (button press + held) and eye fixation (within 2 degrees)
        stimulus.seg{1}.delay = 2;
        stimulus.seg{1}.stim = 3;
        stimulus.seg{1}.resp = 4;

        % Task trial parameters
        task{1}{1}.parameter.target = stimulus.prior;
    case 2 %Exo
        task{1}{1}.segmin = [inf 0.000 0.100 0.200 0.800]; % Fixate, delay, cue, stim, response
        task{1}{1}.segmax = [inf 1.000 0.100 0.200 0.800];

        stimulus.seg = {};
        stimulus.seg{1}.ITI1 = 1;
        stimulus.seg{1}.delay = 2;
        stimulus.seg{1}.cue = 3;
        stimulus.seg{1}.stim = 4;
        stimulus.seg{1}.resp = 5;

    case 3 %Sacc

        task{1}{1}.segmin = [inf 0.000 0.100 0.200 1.000];
        task{1}{1}.segmax = [inf 1.000 0.100 2.000 1.000];

        stimulus.seg{1}.ITI1 = 1;
        stimulus.seg{1}.delay = 2;
        stimulus.seg{1}.stim = 3; % stimulus appears and stays on
        stimulus.seg{1}.delay2 = 4; % variable length delay
        stimulus.seg{1}.saccade = 5; % Fixation cross turns off, subj saccades to stim.
    otherwise
        disp('Invalid attention condition. Quitting...');
        return
end

if stimulus.noeye==1
    task{1}{1}.segmin(1) = 0.5;
    task{1}{1}.segmax(1) = 0.5;
end

task{1}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
task{1}{1}.getResponse = zeros(size(task{1}{1}.segmin));
if stimulus.powerwheel
    task{1}{1}.getResponse(stimulus.seg{1}.resp)=3;
else
    task{1}{1}.getResponse(stimulus.seg{1}.resp)=1;
end
if stimulus.counter<=1
    task{1}{1}.numTrials = 65;
else
    task{1}{1}.numTrials = 15;
end
task{1}{1}.random = 1;

if stimulus.scan
    task{1}{1}.synchToVol(stimulus.seg.ITI) = 1;
end

% Task trial parameters
task{1}{1}.parameter.ecc = stimulus.ecc;
task{1}{1}.parameter.priorSTD = stimulus.sd;

% Task variables to be calculated later
task{1}{1}.randVars.calculated.startRespAngle = nan;
task{1}{1}.randVars.calculated.respAngle = nan;
task{1}{1}.randVars.calculated.angle = nan; % angle at which displayed, depends on attention mode
task{1}{1}.randVars.calculated.rotation = nan; % rotation of the grating
task{1}{1}.randVars.calculated.contrast = nan; % contrast of the grating
task{1}{1}.randVars.calculated.detected = 0; % did they see the grating
task{1}{1}.randVars.calculated.dead = 0;
task{1}{1}.randVars.calculated.visible = 1;
if stimulus.att==2
    task{1}{1}.randVars.calculated.target = nan; 
end

%%%%%%%%%%%%% PHASE TWO %%%%%%%%%%%%%%%%%
%%%%% POSITION JUDGMENT + RESPONSE %%%%%%

stimulus.curTrial(2) = 0;

task{1}{2} = struct;
task{1}{2}.waitForBacktick = 1;

switch stimulus.att
    case 1 %Endo
        % task waits for fixation on first segment
        task{1}{2}.segmin = [inf .200 0.500 5.000];
        task{1}{2}.segmax = [inf .200 0.500 5.000];

        stimulus.seg{2}.ITI1 = 1; % waits for user input (button press + held) and eye fixation (within 2 degrees)
        stimulus.seg{2}.stim = 2;
        stimulus.seg{2}.delay = 3;
        stimulus.seg{2}.resp = 4;

        task{1}{2}.parameter.target = stimulus.prior; % prior center

    case 2
        task{1}{2}.segmin = [inf 0.100 0.200 0.500 5.000]; %fixate, cue, stim, delay, response
        task{1}{2}.segmax = [inf 0.100 0.200 0.500 5.000];

        stimulus.seg{2}.ITI1 = 1; % waits for user input (button press + held) and eye fixation (within 2 degrees)
        stimulus.seg{2}.cue = 2;
        stimulus.seg{2}.stim = 3;
        stimulus.seg{2}.delay = 4;
        stimulus.seg{2}.resp = 5;

    otherwise
        disp('Invalid attention condition. Quitting...');
end

if stimulus.noeye==1
    task{1}{2}.segmin(1) = 0.5;
    task{1}{2}.segmax(1) = 0.5;
end

task{1}{2}.synchToVol = zeros(size(task{1}{2}.segmin));
task{1}{2}.getResponse = zeros(size(task{1}{2}.segmin));
if stimulus.powerwheel
    task{1}{2}.getResponse(stimulus.seg{2}.resp)=3;
else
    task{1}{2}.getResponse(stimulus.seg{2}.resp)=3;
end

task{1}{2}.numTrials = 100;
task{1}{2}.random = 1;

if stimulus.scan
    task{1}{2}.synchToVol(stimulus.seg.ITI) = 1;
end

%% Tracking

% Task trial parameters
task{1}{2}.parameter.ecc = stimulus.ecc; % eccentricity of display
task{1}{2}.parameter.priorSTD = stimulus.sd; % radians
task{1}{2}.parameter.contrastOpt = 1:length(stimulus.live.contrastPercs); %60/70/80/90

% Task variables to be calculated later 
task{1}{2}.randVars.calculated.startRespAngle = nan;
task{1}{2}.randVars.calculated.respAngle = nan;
task{1}{2}.randVars.calculated.angle = nan; % angle at which displayed, depends on attention mode
task{1}{2}.randVars.calculated.rotation = nan; % rotation of the grating
task{1}{2}.randVars.calculated.contrast = nan; % contrast of the grating
task{1}{2}.randVars.calculated.detected = 0; % did they see the grating
task{1}{2}.randVars.calculated.dead = 0; % did the trial get canceled
task{1}{2}.randVars.calculated.visible = 1; % were they shown a grating

if stimulus.att==2
    task{1}{2}.randVars.calculated.target = nan;
end

%% Testing 2
if stimulus.test2
    task{1}{1} = task{1}{2};
    task{1}(2) = [];
end

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

% let the user know
disp(sprintf('(posjdg) Starting run number: %i.',stimulus.counter));

%% Main Task Loop

setGammaTable(1);
mglClearScreen(0.5); mglFixationCross(1,1,stimulus.colors.white);

mglFlush
mglClearScreen(0.5); mglFixationCross(1,1,stimulus.colors.white);

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
mglTextSet([],32,stimulus.colors.white);
% get count
mglTextDraw('Please wait',[0 0]);
mglFlush
myscreen.flushMode = 1;

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

% Move files from posjdg folder to either posjdg_Exo or posjdg_Endo
%
if ~isdir(sprintf('~/data/posjdg_%s/%s',stimulus.condition,mglGetSID))
    disp('No stim directory found, so creating directory now.');
    mkdir(sprintf('~/data/posjdg_%s/%s',stimulus.condition,mglGetSID))
end
movefile(sprintf('~/data/posjdg/%s/*',mglGetSID), sprintf('~/data/posjdg_%s/%s/', stimulus.condition, mglGetSID));

if stimulus.plots
    disp('(posjdg) Displaying plots');
    dispInfo(stimulus);
end

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startTrialCallback(task,myscreen)
%%
if ~isempty(task.lasttrial)
    disp(sprintf('Last trial, target: %02.0f, stim: %02.0f, resp: %02.0f',task.lasttrial.target*180/pi,task.lasttrial.angle*180/pi,task.lasttrial.respAngle*180/pi));
end
global stimulus

if (~isempty(task.lasttrial)) && (task.lasttrial.detected==0) && ~task.lasttrial.dead && task.lasttrial.visible==1
    stimulus.staircase = doStaircase('update',stimulus.staircase,task.lasttrial.detected);
    disp(sprintf('Subject did not see %01.2f%% contrast',task.lasttrial.contrast*100));
elseif (~isempty(task.lasttrial)) && task.lasttrial.visible==0
    disp('No stimulus displayed on this trial.');
end

stimulus.live.gotResponse = 0;
stimulus.curTrial(task.thistrial.thisphase) = stimulus.curTrial(task.thistrial.thisphase) + 1;

% compute missing variables
task.thistrial.rotation = rand*2*pi;
task.thistrial.startRespAngle = rand*2*pi;
switch stimulus.att
    case 1 %Endo
        task.thistrial.angle = randn*task.thistrial.priorSTD; %normally distributed around the prior
    case 2 %Exo
        task.thistrial.angle = rand*2*pi; % stim angle is random
        task.thistrial.target = rand*2*pi; % cue angle is random
        if task.thistrial.thisphase == 1 && ~stimulus.test2
            task.thistrial.visible = (rand > 0.5);
        end
    case 3 %Sacc
        task.thistrial.angle = rand*2*pi;
        task.thistrial.target = rand*2*pi;
    otherwise
        disp('Invalid attention condition. Quitting...');
end

% contrast from staircase
if task.thistrial.thisphase==1 && ~stimulus.test2
    [task.thistrial.contrast, stimulus.staircase] = doStaircase('testValue',stimulus.staircase);
else
    task.thistrial.contrast = stimulus.live.contrastOpts(task.thistrial.contrastOpt);
end

mglSetMousePosition(960,540,1);
myscreen.flushMode = 0;

disp(sprintf('(posjdg) Trial (%i): angle: %02.0f, rotation: %02.0f, contrast: %02.0f%%',...
    task.trialnum,(task.thistrial.angle+task.thistrial.target)*180/pi,task.thistrial.rotation*180/pi,task.thistrial.contrast*100));
    
stimulus.live.eyeCount = 0;

setGammaTable(task.thistrial.contrast);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)
%%

global stimulus

stimulus.live.triggerWaiting = 0;
if any(task.thistrial.thisseg==[stimulus.seg{task.thistrial.thisphase}.ITI1])
    stimulus.live.triggerWaiting = 1;
    stimulus.live.centered = 0;
    stimulus.live.triggerTime = 0;
    stimulus.live.lastTrigger = -1;
end

stimulus.live.eyeDead = 0;
stimulus.live.resp = 0;
stimulus.live.fixColor = stimulus.colors.white;
stimulus.live.fix = 1;
stimulus.live.stim = 0;
stimulus.live.cue = 0;

if task.thistrial.thisseg==stimulus.seg{task.thistrial.thisphase}.stim
    if task.thistrial.visible == 1
        stimulus.live.stim = 1;
    end
elseif task.thistrial.thisseg==stimulus.seg{task.thistrial.thisphase}.resp
    stimulus.live.resp = 1;
    stimulus.live.angle = task.thistrial.startRespAngle;
    stimulus.live.anyAngleAdj = false;
    mInfo = mglGetMouse(myscreen.screenNumber);
    stimulus.live.trackingAngle = -mInfo.x/90;
    convertRespXY(task);
elseif stimulus.att == 2 && task.thistrial.thisseg==stimulus.seg{task.thistrial.thisphase}.cue
    stimulus.live.cue = 1;
end

for i = 1:2
    mglClearScreen(0.5);
    if stimulus.live.stim
        x = task.thistrial.ecc * cos(task.thistrial.angle+task.thistrial.target);
        y = task.thistrial.ecc * sin(task.thistrial.angle+task.thistrial.target);
        mglBltTexture(stimulus.live.grating,[x y],0,0,task.thistrial.rotation*180/pi);
    elseif stimulus.live.cue && i == 1
        x = task.thistrial.ecc * cos(task.thistrial.target);
        y = task.thistrial.ecc * sin(task.thistrial.target);
        mglPolygon([x-.15, x-.15, x+.15, x+.15], [y-.15, y+.15, y+.15, y-.15], stimulus.colors.white);
        mglFlush();
        mglClearScreen(0.5);
    end
    
    % resp is updated in screenUpdate
    
    if stimulus.live.fix
        upFix(stimulus);
    end

    mglFlush
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus

% jump to next trial if you are dead and 1 second has elapsed since eye
% movement
if task.thistrial.dead && mglGetSecs(task.thistrial.segStartSeconds)>1
    task = jumpSegment(task,inf);
end

% skip screen updates if you are already dead
if task.thistrial.dead
    if task.thistrial.dead && stimulus.live.eyeDead
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

% Eye movement detection code
if ~stimulus.noeye && ~any(task.thistrial.thisseg==[stimulus.seg{task.thistrial.thisphase}.ITI1]) && ~stimulus.scan
    if ~any(isnan(pos))
        if dist > 1.5 && stimulus.live.eyeCount > 30
            disp('Eye movement detected!!!!');
            task.thistrial.dead = 1;
            stimulus.live.eyeDead=1;
            return
        elseif dist > 1.5
            stimulus.live.eyeCount = stimulus.live.eyeCount + 1;
        end
    end
end

if (task.thistrial.thisseg==stimulus.seg{task.thistrial.thisphase}.resp) && stimulus.powerwheel
    mInfo = mglGetMouse(myscreen.screenNumber);
    curPos = -mInfo.x/90;
    stimulus.live.angle = stimulus.live.angle + curPos-stimulus.live.trackingAngle;
    if abs(curPos-stimulus.live.trackingAngle)>0
        stimulus.live.anyAngleAdj = true;
    end
    stimulus.live.trackingAngle = curPos;
    convertRespXY(task);
elseif task.thistrial.thisseg==stimulus.seg{task.thistrial.thisphase}.resp
    keys = find(mglGetKeys);
    if any(keys==19)
        stimulus.live.angle = stimulus.live.angle+0.01;
        stimulus.live.anyAngleAdj = true;
    elseif any(keys==20)
        stimulus.live.angle = stimulus.live.angle-0.01;
        stimulus.live.anyAngleAdj = true;
    end
    convertRespXY(task);
end

if stimulus.live.resp && ((task.thistrial.thisphase==2) || stimulus.test2)
    mglClearScreen(0.5);
    if stimulus.live.fix, upFix(stimulus); end
    if stimulus.live.resp, mglFillOval(stimulus.live.respx,stimulus.live.respy,[0.5 0.5],stimulus.colors.white); end
end

% Trial trigger on eye fixation code  
if ~stimulus.noeye && stimulus.live.triggerWaiting
    now = mglGetSecs;
    % check eye position, if 
    if ~any(isnan(pos))
        wasCentered = stimulus.live.centered;
        stimulus.live.centered = dist<2.5;
        if wasCentered && stimulus.live.centered && stimulus.live.lastTrigger>0
            stimulus.live.triggerTime = stimulus.live.triggerTime + now-stimulus.live.lastTrigger;
        end
        stimulus.live.lastTrigger = now;
    end
    if stimulus.live.triggerTime > 0.5 % not in ms dummy, wait 1.5 seconds (reasonable slow time)
        disp('Starting trial--eye centered and space pressed.');
        task = jumpSegment(task);
    end
end

function convertRespXY(task)
global stimulus

stimulus.live.respx = task.thistrial.ecc*cos(stimulus.live.angle+task.thistrial.target);
stimulus.live.respy = task.thistrial.ecc*sin(stimulus.live.angle+task.thistrial.target);

function upFix(stimulus)
%%
% for this experiment use a circle to indicate where participants can
% fixate inside of (rather than a cross which might arbitrarily enforce
% poisitioning
% mglGluAnnulus(0,0,1.5,1.55,stimulus.live.fixColor,64);
mglFixationCross(1,1,stimulus.live.fixColor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus

if task.thistrial.dead, return; end
    
if isfield(task.thistrial,'whichButton') && (task.thistrial.whichButton==stimulus.responseKeys(1))
    % subject didn't see anything
    task = jumpSegment(task,inf);
    task.thistrial.detected = -1;
    disp('Subject reported not seeing anything');
    return
end

if stimulus.powerwheel
    validResponse = task.thistrial.mouseButton == 1;
else
    validResponse = any(task.thistrial.whichButton == stimulus.responseKeys);
end

if validResponse
    if stimulus.live.gotResponse==0
        if (task.thistrial.thisphase==1 && (~stimulus.test2))
            % they saw it
            task.thistrial.detected = 1;
            stimulus.staircase = doStaircase('update',stimulus.staircase,task.thistrial.detected);
            disp(sprintf('Subject saw %01.2f%% contrast',task.thistrial.contrast*100));
            stimulus.live.fix = 0;
        elseif ((task.thistrial.thisphase==2) || stimulus.test2)
            % they saw it
%             if stimulus.live.anyAngleAdj == true
            task.thistrial.detected = 1;
            stimulus.staircase = doStaircase('update',stimulus.staircase,task.thistrial.detected);
%             end
            % they are actually reporting locations
            task.thistrial.respAngle = stimulus.live.angle;
            disp(sprintf('Subject reported %02.0f real %02.0f at %02.0f%% contrast',task.thistrial.respAngle*180/pi,task.thistrial.angle*180/pi,task.thistrial.contrast*100));
            stimulus.live.fix = 0;
            stimulus.live.resp = 0;
            task = jumpSegment(task,inf);
        end
    else
        disp(sprintf('Subject responded multiple times: %i',stimulus.live.gotResponse));
    end
    stimulus.live.gotResponse=stimulus.live.gotResponse+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              HELPER FUNCTIONS                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function initStair()
global stimulus

stimulus.staircase = doStaircase('init','upDown',...
            'initialThreshold',0.25,...
            'initialStepsize',0.25/3,...
            'minThreshold=0.0001','maxThreshold=0.4','stepRule','pest',...
            'nTrials=40','maxStepsize=0.2','minStepsize=0.0001');
        
function resetStair()

global stimulus

if doStaircase('stop',stimulus.staircase)
    disp('(posjdg) Staircase is being reset');
    stimulus.staircase(end+1) = doStaircase('init',stimulus.staircase(end));
    if stimulus.staircase(end).s.threshold>0.3
        disp('(posjdg) Bad staircase threshold: setting to 0.3');
        stimulus.staircase(end).s.threshold=0.3;
    elseif stimulus.staircase(end).s.threshold<0
        disp('(posjdg) Bad staircase threshold: setting to 0.05');
        stimulus.staircase(end).s.threshold=0.05;
    end
end

function [trials] = totalTrials()
%%

% Counts trials + estimates the threshold based on the last 500 trials

% get the files list
files = dir(fullfile(sprintf('~/data/posjdg/%s/17*stim*.mat',mglGetSID)));

trials = 0;

for fi = 1:length(files)
    load(fullfile(sprintf('~/data/posjdg/%s/%s',mglGetSID,files(fi).name)));
    
    e = getTaskParameters(myscreen,task);
    e = e{1}; % why?!
    trials = trials + e.nTrials;
end

%%%%%%%%%%%%%%%%%%%%%%%
%    dispInfo    %
%%%%%%%%%%%%%%%%%%%%%%%
function dispInfo(rstimulus)
%%

doStaircase('threshold',rstimulus.staircase,'dispFig=1','type=weibull');
% ctask = task; cscreen = myscreen; % save this incase we need them

% compute % correct for valid and invalid trials, display learning over
% time (including history from other runs)
% exp = getTaskParameters(task,myscreen);

% get the files list
files = dir(fullfile(sprintf('~/data/posjdg_%s/%s/17*stim*.mat',rstimulus.condition,mglGetSID)));

% load the files and pull out the data (long form)
%  rrun # counter #    local trial     real trial   angle     respAngle    
%     1       2             3              4           5           6
%  target    startRespAngle     contrast     detected      ecc    priorsd
%     7            8                9           10          11      12
%    rotation
%       13
count = 1; data = zeros(10000,13);

for fi = 1:length(files)
    load(fullfile(sprintf('~/data/posjdg_%s/%s/%s',rstimulus.condition,mglGetSID,files(fi).name)));
    
    e = getTaskParameters(myscreen,task);
    if e{1}(1).nTrials>1
        e = e{1}(2); % why?!
    
        run = stimulus.counter;

        if stimulus.att==2
            data(count:count+(e.nTrials-1),:) = [repmat(fi,e.nTrials,1) repmat(run,e.nTrials,1) (1:e.nTrials)' (count:count+(e.nTrials-1))' ...
                e.randVars.angle' e.randVars.respAngle' e.randVars.target' ...
                e.randVars.startRespAngle' e.randVars.contrast' e.randVars.detected' ...
                e.parameter.ecc' e.parameter.priorSTD' e.randVars.rotation'];
        else
            data(count:count+(e.nTrials-1),:) = [repmat(fi,e.nTrials,1) repmat(run,e.nTrials,1) (1:e.nTrials)' (count:count+(e.nTrials-1))' ...
                e.randVars.angle' e.randVars.respAngle' e.parameter.target' ...
                e.randVars.startRespAngle' e.randVars.contrast' e.randVars.detected' ...
                e.parameter.ecc' e.parameter.priorSTD' e.randVars.rotation'];
        end

        count = count+e.nTrials;
    end
end

data = data(1:(count-1),:);

if any(data(:,6)>pi), data(data(:,6)>pi,6) = data(data(:,6)>pi,6)-2*pi; end
%% Compute angle-target and respAngle-target plot
h = figure; hold on

low = [0 0.075 0.081 0.09 inf];

%for i = 1:4
    %subplot(4,1,i); hold on
    % data(:,6) = data(:,6)-data(:,5);

    % remove no-response trials
    data_ = data(~isnan(data(:,6)),:);
    %data_ = data_(logical((data_(:,9)>low(i)).*(data_(:,9)<low(i+1))),:);
    % find the trials where stimulus is - relative to the prior
    flip = data_(:,5)<0; flip = flip*1;
    flip(flip==1) = -1; flip(flip==0) = 1;
    % flip all the stimulus-target to be in the positive space
    data_(:,5:6) = data_(:,5:6).* repmat(flip,1,2);
    Y = data_(:,6);
    X = [ones(size(Y)) data_(:,5)];
    b = X\Y;
    c = [X Y];

    bci = bootci(1000,@(x) x(:,1:2)\x(:,3),c);

    % [p,s] = polyfit(data_(:,5),data_(:,6),1);


    plot(data_(:,5),data_(:,6),'*');
    % plot constant line
    plot([-1 1],[-1 1],'--r');
    % plot fit
    x = -1:1;
    plot(x,b(1)+b(2)*x,'--k');
    % compute SD of residuals? 
    % todo
    xlabel('Stimulus - Target (deg)');
    ylabel('Resp - Target (deg)');
    title(sprintf('Bias %01.2f [%01.2f %01.2f], slope %01.2f [%01.2f %01.2f], Steeper = resp away, shallower = resp toward',b(1),bci(1,1),bci(2,1),b(2),bci(1,2),bci(2,2)));

    axis([0 1 -1 1]);
    axis square

    set(gca,'XTick',0:.5:1,'XTickLabel',round((0:.5:1)*180/pi,2),'YTick',-1:.5:1,'YTickLabel',round((-1:.5:1)*180/pi,2));

    drawPublishAxis;
    % h = figure;
    % hist(data(:,5)-data(:,6));
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localInitStimulus()

global stimulus

% stimulus.live.pos = 
% stimulus.live.pos = logspace(log10(stimulus.cur_.isize),log10(stimulus.cur_.osize),stimulus.cur_.N);

% mglClearScreen(0.5)
% for gi = 1:stimulus.cur_.N
    % for each grating distance
    % calculate the center position to estimate the radius
%     crad = stimulus.live.pos(gi);
    % get total degrees around circle
%     degs = 2*pi*crad;
%     sz = degs/stimulus.cur_.num*0.6;
sz = 1.5;
% use total degs / num to compute size
grating = 251/2*mglMakeGrating(sz,sz,2,0) + 255/2;
gauss = mglMakeGaussian(sz,sz,sz/6,sz/6);
alphamask = repmat(grating,1,1,4);
alphamask(:,:,4) = gauss*255;

% we'll adjust the gamma table to control contrast
stimulus.live.grating  = mglCreateTexture(alphamask); % high contrast

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sets the gamma table so that we can have
% finest possible control over the stimulus contrast.
%
% stimulus.colors.reservedColors should be set to the reserved colors (for cue colors, etc).
% maxContrast is the maximum contrast you want to be able to display.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setGammaTable(maxContrast)

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
    redLinearized = repmat(interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,.5,'linear'),1,stimulus.colors.nUnreserved);
    greenLinearized = repmat(interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,.5,'linear'),1,stimulus.colors.nUnreserved);
    blueLinearized = repmat(interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,.5,'linear'),1,stimulus.colors.nUnreserved);
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
succ = mglSetGammaTable(gammaTable);

if ~succ
    warning('Gamma table set failure');
    keyboard
end

% remember what the current maximum contrast is that we can display
stimulus.curMaxContrast = maxContrast;


function initGammaTable(myscreen)
global stimulus
%% Gamma Table Initialization

% get gamma table
if ~isfield(myscreen,'gammaTable')
  stimulus.linearizedGammaTable = mglGetGammaTable;
  disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
  disp(sprintf('(cuecon:initGratings) No gamma table found in myscreen. Contrast displays like this'));
  disp(sprintf('         should be run with a valid calibration made by moncalib for this monitor.'));
  disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
else
  % check to make sure this matches the calibration file
  
  % get each channel table that should have been set by mglGetGammaTable
  redTable = myscreen.initScreenGammaTable.redTable(:);
  greenTable = myscreen.initScreenGammaTable.greenTable(:);
  blueTable = myscreen.initScreenGammaTable.blueTable(:);
  % get what the calibration structure says it should have been set to
  gammaTable = myscreen.gammaTable(:);
  % table values are only good to 10 bits
  redTable = round(redTable*1024)/1024;
  greenTable = round(greenTable*1024)/1024;
  blueTable = round(blueTable*1024)/1024;
  gammaTable = round(gammaTable*1024)/1024;
  % compare, ignoring nans
  if ~isequaln(mglGetGammaTable,myscreen.initScreenGammaTable) || ~isequaln(redTable,gammaTable) || ~isequaln(greenTable,gammaTable) || ~isequaln(blueTable,gammaTable)
    disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
    disp(sprintf('(curecon:initGrating) Gamma table does not match calibration'));
    disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
    keyboard
  end
end
stimulus.linearizedGammaTable = myscreen.initScreenGammaTable;
