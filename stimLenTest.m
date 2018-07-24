function [ myscreen ] = stimLenTest( varargin )
%
% Saccade task to measure compression 
%
%  Usage: stimLenTest(varargin)
%  Authors: Akshay Jagadeesh
%  Date: 02/27/2018
%

global stimulus

stimulus = struct;

%% Initialize Variables

% add arguments later
scan = 0;
plots = 0;
noeye = 1;
debug = 0;
getData = 0;
fixSide = 0;
getArgs(varargin,{'getData=0','scan=0','plots=0','noeye=1','debug=0', 'fixSide=0'});
stimulus.scan = scan;
stimulus.plots = plots;
stimulus.noeye = noeye;
stimulus.debug = debug;
stimulus.getData = getData;
stimulus.fixSide = fixSide;
clear localizer invisible scan noeye task test2

if stimulus.plots
  dispInfo(stimulus);
  myscreen = 0;
  return
end

%% Stimulus parameters 
%% Open Old Stimfile
stimulus.counter = 1;
if ~isempty(mglGetSID) && isdir(sprintf('~/data/saccadeRF/%s',mglGetSID))
  % Directory exists, check for a stimfile
  files = dir(sprintf('~/data/saccadeRF/%s/1*mat',mglGetSID));

  if length(files) >= 1
    fname = files(end).name;
    
    s = load(sprintf('~/data/saccadeRF/%s/%s',mglGetSID,fname));
    stimulus.counter = s.stimulus.counter + 1;
    clear s;
    disp(sprintf('(saccadeRF) Data file: %s loaded.',fname));
  end
end
disp(sprintf('(saccadeRF) This is run #%i',stimulus.counter));

%% Setup Screen
if stimulus.scan
  myscreen = initScreen('fMRIprojFlex');
else
  myscreen = initScreen('VPixx2');
end

% set background to grey
myscreen.background = 0.5; 

%% Setup missing initial variables
if ~isfield(stimulus,'counter')
  stimulus.counter = 1; % This keeps track of what "run" we are on.
end

%% Initialize Stimulus
myscreen = initStimulus('stimulus',myscreen);

%localInitStimulus();
  
% Set response keys
stimulus.responseKeys = [11 12 13 14];

% set colors
stimulus.colors.white = [1 1 1];
stimulus.colors.black = [0 0 0];
stimulus.colors.red = [1 0 0];
stimulus.colors.green = [0 1 0];
stimulus.colors.blue = [0 0 1];
stimulus.live.fixColor = stimulus.colors.blue;
stimulus.live.cueColor = stimulus.colors.black;

stimulus.curTrial(1) = 0;

% Number of Trials
nTrials = 10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% TASK 1: Left Side %%%%%%%%%%%%%%%%%%%%%%
task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;

%% Define stimulus timing
stimTime = 1.000;
ITI_min = 3.000; ITI_max = 12.00;
task{1}{1}.segmin = [stimTime ITI_min];
task{1}{1}.segmax = [stimTime ITI_max];
stimulus.seg = {};
stimulus.seg{1}.stim = 1;
stimulus.seg{1}.ITI = 2;

% Task important variables
stimulus.eccentricity = 10;
stimulus.stimSize = 2;
stimulus.stimLength = 2; % Number of frames to present stimulus for.

task{1}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
task{1}{1}.getResponse = zeros(size(task{1}{1}.segmin));
task{1}{1}.numTrials = nTrials;
task{1}{1}.random = 1;
if stimulus.scan
  task{1}{1}.synchToVol(stimulus.seg{1}.ITI) = 1;
end

% Task trial parameters

% Task variables to be calculated late
task{1}{1}.randVars.calculated.detected = 0;
task{1}{1}.randVars.calculated.dead = 0;
task{1}{1}.randVars.calculated.visible = 1;

task{1}{1}.randVars.calculated.stimTime = NaN;


%% Initialize stimulus
initCheckerboard(stimulus.eccentricity, stimulus.stimSize, 4);
%keyboard
%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{1})
  [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% TASK 2: Left Side %%%%%%%%%%%%%%%%%%%%%%
task{2}{1} = struct;
task{2}{1}.waitForBacktick = 1;

%% Define stimulus timing
task{2}{1}.segmin = [stimTime ITI_min];
task{2}{1}.segmax = [stimTime ITI_max];
stimulus.seg = {};
stimulus.seg{1}.stim = 1;
stimulus.seg{1}.ITI = 2;

% Trial parameters
task{2}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
task{2}{1}.getResponse = zeros(size(task{1}{1}.segmin));
task{2}{1}.numTrials = nTrials;
task{2}{1}.random = 1;
if stimulus.scan
  task{2}{1}.synchToVol(stimulus.seg{1}.ITI) = 1;
end

% Task variables to be calculated late
task{2}{1}.randVars.calculated.detected = 0;
task{2}{1}.randVars.calculated.dead = 0;
task{2}{1}.randVars.calculated.visible = 1;

task{2}{1}.randVars.calculated.stimTime = NaN;

%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{2})
  [task{2}{phaseNum}, myscreen] = initTask(task{2}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% TASK 3: Center Side %%%%%%%%%%%%%%%%%%%%%%
task{3}{1} = struct;
task{3}{1}.waitForBacktick = 1;

%% Define stimulus timing
task{3}{1}.segmin = [stimTime ITI_min];
task{3}{1}.segmax = [stimTime ITI_max];
stimulus.seg = {};
stimulus.seg{1}.stim = 1;
stimulus.seg{1}.ITI = 2;

% Trial parameters
%task{3}{1}.parameter.stimLength = [1,2,3]; % Number of frames to present stimulus for
%task{2}{1}.parameter.stimSize = [1,2,4];

task{3}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
task{3}{1}.getResponse = zeros(size(task{1}{1}.segmin));
task{3}{1}.numTrials = nTrials;
task{3}{1}.random = 1;
if stimulus.scan
  task{3}{1}.synchToVol(stimulus.seg{1}.ITI) = 1;
end

% Task variables to be calculated late
task{3}{1}.randVars.calculated.detected = 0;
task{3}{1}.randVars.calculated.dead = 0;
task{3}{1}.randVars.calculated.visible = 1;

task{3}{1}.randVars.calculated.stimTime = NaN;

%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{3})
  [task{3}{phaseNum}, myscreen] = initTask(task{3}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% TASK 4: Center Side %%%%%%%%%%%%%%%%%%%%%%
task{4}{1} = struct;
task{4}{1}.waitForBacktick = 1;

%% Define stimulus timing
task{4}{1}.segmin = [3.00 1.00];
task{4}{1}.segmax = [12.00 1.00];

stimulus.seg{2}.fix = 1;
stimulus.seg{2}.sacc = 2;

task{4}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
task{4}{1}.getResponse = zeros(size(task{1}{1}.segmin));
task{4}{1}.numTrials = nTrials;
task{4}{1}.random = 1;

% Task variables to be calculated late
task{4}{1}.randVars.calculated.greenSide = NaN;

%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{4})
  [task{4}{phaseNum}, myscreen] = initTask(task{4}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);
end

%% EYE CALIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

% let the user know
disp(sprintf('(saccadeRF) Starting run number: %i.',stimulus.counter));

%% Main Task Loop

mglClearScreen(0.5); 
upFix(stimulus);

mglFlush
mglClearScreen(0.5); 
upFix(stimulus);

phaseNum = 1;
% Again, only one phase.
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
  mglClearScreen(0.5);
  % update the left task
  [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);
  % update the right task
  [task{2}, myscreen, phaseNum] = updateTask(task{2}, myscreen, phaseNum);
  % update the central task
  [task{3}, myscreen, phaseNum] = updateTask(task{3}, myscreen, phaseNum);
  % update the crosses
  [task{4}, myscreen, phaseNum] = updateTask(task{4}, myscreen, phaseNum);
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

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startTrialCallback(task,myscreen)

global stimulus

task.thistrial.dead = 0;
task.thistrial.detected = 0;
task.thistrial.visible = 1;
task.thistrial.response = 0;

stimulus.live.gotResponse = 0;
if task.taskID == 1
  stimulus.live.stim1 = 0;
  stimulus.live.task1Dead = 0;
elseif task.taskID == 2
  stimulus.live.stim2 = 0;
  stimulus.live.task2Dead = 0;
elseif task.taskID == 3
  stimulus.live.stim3 = 0;
  stimulus.live.task3Dead = 0;
end
if task.taskID == 1
    stimulus.live.t1 = nan;
elseif task.taskID == 2
    stimulus.live.t2 = nan;
elseif task.taskID == 3
    stimulus.live.t3 = nan;
end

if task.taskID == 4
  task.thistrial.greenSide = mod(task.trialnum,2);
end

stimulus.curTrial(task.thistrial.thisphase) = stimulus.curTrial(task.thistrial.thisphase) + 1;

% Disp trial parameters each trial
if task.taskID == 1
  disp(sprintf('Trial %d (L) - numFrames: %d', task.trialnum, stimulus.stimLength));
end

% Reset mouse to center of screen at start of every trial
mglSetMousePosition(960,540,1);
myscreen.flushMode = 0;
stimulus.live.eyeCount = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)

global stimulus

stimulus.live.text = nan;
stimulus.live.triggerWaiting = 0;
stimulus.live.eyeDead = 0;
stimulus.live.fix = 1;

if task.thistrial.thisseg == stimulus.seg{1}.ITI && task.taskID == 1
  disp(sprintf('Drawing %d frames took %g secs', stimulus.stimLength, stimulus.live.t1));
  task.thistrial.stimTime = stimulus.live.t1;
elseif task.thistrial.thisseg == stimulus.seg{1}.ITI && task.taskID == 2
  task.thistrial.stimTime = stimulus.live.t2;
elseif task.thistrial.thisseg == stimulus.seg{1}.ITI && task.taskID == 3
  task.thistrial.stimTime = stimulus.live.t3;
elseif task.taskID == 4
  stimulus.live.greenCtr = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus

% Draw two fixation crosses
mglFixationCross(1,1,stimulus.colors.white, [-5,0]);
mglFixationCross(1,1,stimulus.colors.white, [5,0]);

% Code for drawing checkerboard on each frame
if task.taskID == 1 && task.thistrial.thisseg==stimulus.seg{1}.stim && ~stimulus.live.task1Dead% LEFT side task
  % set clock
  if stimulus.live.stim1 == 0
      stimulus.live.t1 = mglGetSecs;
  end
  
  %draw stim
  if stimulus.live.stim1 < stimulus.stimLength
    %drawCheckerboard(-stimulus.eccentricity, 0, stimulus.stimSize);
    drawCheckerboard(stimulus.live.leftStim);
    stimulus.live.stim1 = stimulus.live.stim1 + 1;
  else
    stimulus.live.t1 = mglGetSecs(stimulus.live.t1);
    stimulus.live.task1Dead = 1;
  end
elseif task.taskID == 2 && task.thistrial.thisseg==stimulus.seg{1}.stim && ~stimulus.live.task2Dead % RIGHT side task
  % set clock
  if stimulus.live.stim2 == 0
    stimulus.live.t2 = mglGetSecs;
  end
  
  %draw stim
  if stimulus.live.stim2 < stimulus.stimLength
    drawCheckerboard(stimulus.live.rightStim);
    %drawCheckerboard(stimulus.eccentricity, 0, stimulus.stimSize);
    stimulus.live.stim2 = stimulus.live.stim2 + 1;
  else
    stimulus.live.t2 = mglGetSecs(stimulus.live.t2);
    stimulus.live.task2Dead = 1;
  end
elseif task.taskID == 3 && task.thistrial.thisseg==stimulus.seg{1}.stim && ~stimulus.live.task3Dead % Center task
  if stimulus.live.stim3 == 0
    stimulus.live.t3 = mglGetSecs;
  end

  if stimulus.live.stim3 < stimulus.stimLength
    drawCheckerboard(stimulus.live.centerStim);
    stimulus.live.stim3 = stimulus.live.stim3 +1;
  else
    stimulus.live.t3 = mglGetSecs(stimulus.live.t3);
    stimulus.live.task3Dead = 1;
  end
end


% Change the color of the crosses to green every now and then.
fixSide = [-5, 5];
if task.taskID == 4
  if task.thistrial.thisseg == stimulus.seg{2}.sacc && stimulus.live.greenCtr < 10
    mglFixationCross(1,1,stimulus.colors.green,[fixSide(task.thistrial.greenSide+1), 0]);
    stimulus.live.greenCtr = stimulus.live.greenCtr + 1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus

if task.thistrial.dead, return; end

validResponse = any(task.thistrial.whichButton == stimulus.responseKeys);

if validResponse
  if stimulus.live.gotResponse==0
    task.thistrial.detected = 1;
    task.thistrial.response = task.thistrial.whichButton - 10;
    stimulus.live.fix = 0;
  else
    disp(sprintf('Subject responded multiple times: %i',stimulus.live.gotResponse));
  end
  stimulus.live.gotResponse=stimulus.live.gotResponse+1;
  task = jumpSegment(task);
else
  disp(sprintf('Invalid response key. Subject pressed %d', task.thistrial.whichButton));
  task.thistrial.response = -1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              HELPER FUNCTIONS                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
% Draws a circle at center of the screen of color fixColor
function upFix(stimulus, fixColor)
if ieNotDefined('fixColor')
  fixColor = stimulus.live.fixColor;
end
%mglGluAnnulus(0,0,0,.1,fixColor);
mglFixationCross(1,1,fixColor);


%%% 
% Draws a circular cue at location x,y
function drawCue(x,y, stimulus)
mglGluAnnulus(x,y, 0.75, 0.8, stimulus.live.cueColor, 64);

%%%
% Turns image into a texture
function tex = genTexFromIm(im)
r = flipud(im);
r(:,:,4) = 255;
rP = permute(r, [3 2 1]);
tex = mglCreateTexture(rP);  

function [trials] = totalTrials()
%%
% Counts trials + estimates the threshold based on the last 500 trials
% get the files list
files = dir(fullfile(sprintf('~/data/saccadeRF/%s/18*stim*.mat',mglGetSID)));
trials = 0;

for fi = 1:length(files)
    load(fullfile(sprintf('~/data/saccadeRF/%s/%s',mglGetSID,files(fi).name)));
    e = getTaskParameters(myscreen,task);
    e = e{1}; % why?!
    trials = trials + e.nTrials;
end

%%%%%%%%%%%%%%%%%%%%%%%
%    dispInfo    %
%%%%%%%%%%%%%%%%%%%%%%%
function data = dispInfo(rstimulus)
%%

% get the files list
files = dir(fullfile(sprintf('~/data/saccadeRF/%s/18*stim*.mat',mglGetSID)));

count = 1; 
data = struct('nTrials', 0, 'subj_resp', [], 'corr_resp', [], 'corr_trials', [],...
              'image', [], 'layer', [], 'ecc', [], 'reaction_time', [], 'nValTrials', 0, 'rf_size', [], 'imSz', []);

for fi = 1:length(files)
  load(fullfile(sprintf('~/data/saccadeRF/%s/%s',mglGetSID,files(fi).name)));
  
  e = getTaskParameters(myscreen,task);
  if e{1}.nTrials>1
    
    subj_resp = e{1}.response-10;
    corr_resp = e{1}.randVars.targetPosition;
    data.run = stimulus.counter;
    data.subj_resp = [data.subj_resp subj_resp];
    data.corr_resp = [data.corr_resp corr_resp];
    data.corr_trials = [data.corr_trials subj_resp==corr_resp];
    data.reaction_time = [data.reaction_time e{1}.reactionTime];
    data.nTrials = data.nTrials + e{1}.nTrials;
    data.imSz = [data.imSz e{1}.randVars.imSz];
    % Calculate number of valid trials by excluding eye movements and pool5
    data.nValTrials = data.nValTrials + sum(~isnan(e{1}.response)) - sum(e{1}.parameter.layer == 5);
    
    data.image = [data.image e{1}.randVars.targIm];
    data.layer = [data.layer e{1}.parameter.layer];
    data.rf_size = [data.rf_size e{1}.parameter.rfSize];
    data.ecc = [data.ecc e{1}.parameter.eccentricity];
    
  end
  count = count + 1;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% drawCheckerboard()
% Draws checkerboard to screen using mglQuad
% Args:
%   - stim: struct returned by generateCheckerboard
function drawCheckerboard(stim)

mglQuad(stim.x, stim.y, stim.rgb);

%sqSz = sz/4;

%mglQuad([x-sqSz, x, x-sqSz, x; x, x+sqSz, x, x+sqSz; x, x+sqSz, x, x+sqSz; x-sqSz, x, x-sqSz, x],...
%        [y+sqSz, y+sqSz, y, y; y+sqSz, y+sqSz, y, y; y, y, y-sqSz, y-sqSz; y, y, y-sqSz, y-sqSz],...
%        [0, 1, 1, 0; 0, 1, 1, 0; 0, 1, 1, 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% initCheckerboard()
%  Function to initialize left and right checkerboards for this task
%     and save them into stimulus.live.leftStim and rightStim
%  Args:
%    - ecc: eccentricity of left and right stimulus centers
%    - sz: size in degrees along a side
%    - nSq: number of squares along a side
function initCheckerboard(ecc, sz, nSq)
% ecc = 6 (x=+/-6, y=0), sz = 2, nSq = 4
global stimulus;

stimulus.live.leftStim = generateCheckerboard(-ecc,0, sz, nSq);
stimulus.live.rightStim = generateCheckerboard(ecc,0, sz, nSq);
stimulus.live.centerStim = generateCheckerboard(0,0, sz, nSq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generateCheckerboard()
%   Function to pre-generate checkerboards
% Args:
%   - center position (x0,y0)
%   - height/width of checkerboard (sz)
%   - number of squares per side.
% Returns: 
%   - struct containing x,y, and rgb fields which can be used as arguments for mglQuad
function checkerboard = generateCheckerboard(x0,y0, sz, nSq)
sqSz = sz / nSq;
lim = nSq/2;

x = nan(4, nSq * nSq); y = nan(4, nSq*nSq); rgb = nan(3, nSq*nSq);
x(1,:) = repmat(x0-lim*sqSz : sqSz : x0 + (lim-1)*sqSz, 1, nSq);
x(2,:) = repmat(x0-(lim-1)*sqSz : sqSz : x0 + lim*sqSz, 1, nSq);
x(3,:) = x(2,:);
x(4,:) = x(1,:);

y(1,:) = vectify(repmat(y0+lim*sqSz: -sqSz: y0-(lim-1)*sqSz, nSq, 1))';
y(2,:) = y(1,:);
y(3,:) = vectify(repmat(y0+(lim-1)*sqSz: -sqSz: y0-lim*sqSz, nSq, 1))';
y(4,:) = y(3,:);

a = repmat([1 0], 1, lim);
b = flip(a);
for i = 1:3
  rgb(i,:) = repmat([a, b], 1, lim);
end

% for each side, return an x and y array
checkerboard = struct('x', x, 'y', y, 'rgb', rgb);
return;

function v = vectify(vec)
v = vec(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localInitStimulus()

global stimulus
sz = 1.5;

% use total degs / num to compute size
grating = 251/2*mglMakeGrating(sz,sz,2,0) + 255/2;
gauss = mglMakeGaussian(sz,sz,sz/6,sz/6);
alphamask = repmat(grating,1,1,4);
alphamask(:,:,4) = gauss*255;

alphamask = ones(20,20);
alphamask(1:10, 1:10) = 0;
alphamask(10:20, 10:20) = 0;

% we'll adjust the gamma table to control contrast
stimulus.live.grating  = mglCreateTexture(alphamask); % high contrast

