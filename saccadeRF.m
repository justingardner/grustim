function [ myscreen ] = saccadeRF( varargin )
%
% Saccade task to measure compression 
%
%  Usage: saccadeRF(varargin)
%  Authors: Akshay Jagadeesh
%  Date: 06/11/2018
%

global stimulus

stimulus = struct;

%% Initialize Variables

% add arguments later
scan = 0;
plots = 0;
noeye = 1;
getArgs(varargin,{'scan=0','plots=0','noeye=1'});
stimulus.scan = scan;
stimulus.plots = plots;
stimulus.noeye = noeye;
clear localizer invisible scan noeye task test2

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

% set colors
stimulus.colors.white = [1 1 1];
stimulus.colors.black = [0 0 0];
stimulus.colors.red = [1 0 0];
stimulus.colors.green = [0 1 0];
stimulus.colors.blue = [0 0 1];

stimulus.curTrial(1) = 0;

% Number of Trials
nTrials = 10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% TASK 1: Left Side %%%%%%%%%%%%%%%%%%%%%%
task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;

% Define stimulus timing
task{1}{1}.segmin = [2.00 1.00];
task{1}{1}.segmax = [5.00 2.00];
stimulus.seg = {};
stimulus.seg{1}.fix = 1;
stimulus.seg{1}.sacc = 2;

% Task important variables
stimulus.stimLength = 2;
stimulus.stimSize = 2;
stimulus.fixEcc = 5;
stimulus.saccLatency = 0.180; % Time to wait after the cue before flashing stimulus.
stimulus.cueLength = 4; % Number of frames to change color of fixation cross

% Trial parameters
task{1}{1}.parameter.stimPos = [-10, 0, 10]; % x position 
task{1}{1}.parameter.stimPresent = [0,1,1,1];

task{1}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
task{1}{1}.getResponse = zeros(size(task{1}{1}.segmin));
task{1}{1}.numTrials = nTrials;
task{1}{1}.random = 1;
if stimulus.scan
  task{1}{1}.synchToVol(stimulus.seg{1}.ITI) = 1;
end

% Task trial parameters

% Task variables to be calculated late
task{1}{1}.randVars.calculated.whichSide = NaN; % -1 is left, 1 is right.
task{1}{1}.randVars.calculated.detected = 0;
task{1}{1}.randVars.calculated.dead = 0;
task{1}{1}.randVars.calculated.visible = 1;

task{1}{1}.randVars.calculated.stimTime = NaN;

%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{1})
  [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);
end
%% EYE CALIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

% let the user know
disp(sprintf('(saccadeRF) Starting run number: %i.',stimulus.counter));

%% Main Task Loop

mglClearScreen(0.5); 
%upFix(stimulus);

mglFlush
mglClearScreen(0.5); 
%upFix(stimulus);

phaseNum = 1;
% Again, only one phase.
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
  mglClearScreen(0.5);
  % update the left task
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

% Alternate side of the stimulus on each trial
task.thistrial.whichSide = mod(task.trialnum,2); % 0 for left, 1 for right
task.thistrial.saccTarget = mod(task.trialnum+1,2);

stimulus.live.gotResponse = 0;
stimulus.curTrial(task.thistrial.thisphase) = stimulus.curTrial(task.thistrial.thisphase) + 1;

% Init checkerboard
initCheckerboard(task.thistrial.stimPos, stimulus.stimSize, 4);

% Disp trial parameters each trial
disp(sprintf('Trial %d - Side: %d, Stim Present: %d, Stim Position: %d', task.trialnum, task.thistrial.whichSide, task.thistrial.stimPresent, task.thistrial.stimPos));

% Reset mouse to center of screen at start of every trial
%mglSetMousePosition(960,540,1);
myscreen.flushMode = 0;
stimulus.live.eyeCount = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)

global stimulus

stimulus.live.nFrames = 0;
stimulus.live.tCue = NaN; % timer for each cue
stimulus.live.cueFrames = 0; % counter for displaying cue
stimulus.live.nFrames = 0; % counter for displaying stimulus
stimulus.live.saccStart = 0;
stimulus.live.saccLat = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus

fixLocs = [-stimulus.fixEcc, stimulus.fixEcc];
mglFixationCross(1,1,stimulus.colors.white,[fixLocs(1),0]);
mglFixationCross(1,1,stimulus.colors.white,[fixLocs(2),0]);

%% Get eye position
pos = mglEyelinkGetCurrentEyePos; % NaN if noeye
xEye = pos(1);
distFromFix = abs(fixLocs(task.thistrial.whichSide+1) - xEye);

%% Saccade Segment
if task.thistrial.thisseg == stimulus.seg{1}.sacc 

  % Draw cue at the start of segment.
  if stimulus.live.cueFrames < stimulus.cueLength
    if stimulus.live.cueFrames == 0 % Start timer on cue onset
      stimulus.live.tCue = mglGetSecs;
    end
    % Turn current side fixation cross green for cueLength frames.
    mglFixationCross(1,1,stimulus.colors.green,[fixLocs(task.thistrial.whichSide+1),0]);
    stimulus.live.cueFrames = stimulus.live.cueFrames +1;
  end

  %% Check if saccade has been initiated
  if ~isnan(distFromFix) && distFromFix > 0.5 && stimulus.live.saccStart == 0
    stimulus.live.saccStart = 1;
    stimulus.live.saccLat = mglGetSecs(stimulus.live.tCue);
    disp(sprintf('Saccade initiated. Latency = %g seconds', stimulus.live.saccLat));
  end
    
  % Draw checkerboard 180ms after cue onset.
  if task.thistrial.stimPresent == 1 && mglGetSecs(stimulus.live.tCue) > stimulus.saccLatency
    if stimulus.live.nFrames < stimulus.stimLength
      drawCheckerboard(stimulus.live.stim);
      stimulus.live.nFrames = stimulus.live.nFrames + 1;
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus

if task.thistrial.dead, return; end

validResponse = any(task.thistrial.whichButton == stimulus.responseKeys);

if validResponse
  if stimulus.live.gotResponse==0
    task.thistrial.detected = 1;
    task.thistrial.response = task.thistrial.whichButton - 10;
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
%   - stim: struct returned by generateCheckerboard
function drawCheckerboard(stim)
mglQuad(stim.x, stim.y, stim.rgb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% initCheckerboard()
%  Function to initialize left and right checkerboards for this task
%     and save them into stimulus.live.leftStim and rightStim
function initCheckerboard(xPos, sz, nSq)
global stimulus;
stimulus.live.stim = generateCheckerboard(xPos, 0, sz, nSq);

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
