function [ myscreen ] = texTest( varargin )
%
% EVENT-RELATED TEXTURES
%  Experiment to map neural responses to various textures
%
%  Usage: texTest(varargin)
%  Authors: Akshay Jagadeesh
%  Date: 11/05/2018
%

global stimulus

stimulus = struct;

%% Initialize Variables

% add arguments later
scan = 0;
run = 0;
getArgs(varargin,{'scan=1', 'sType=1'}, 'verbose=1');
stimulus.scan = scan;
stimulus.type = sType;
clear scan run type;

%% Stimulus parameters 
%% Open Old Stimfile
stimulus.counter = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/texTest/%s',mglGetSID))
  % Directory exists, check for a stimfile
  files = dir(sprintf('~/data/texTest/%s/1*mat',mglGetSID));

  if length(files) >= 1
    fname = files(end).name;
    
    s = load(sprintf('~/data/texTest/%s/%s',mglGetSID,fname));
    stimulus.counter = s.stimulus.counter + 1;
    clear s;
    disp(sprintf('(texTest) Data file: %s loaded.',fname));
  end
end
trialTypes = {'FAST_no-baseline', 'SLOW_no-baseline', 'SLOW_baseline'};
disp(sprintf('-------------------------'));
disp(sprintf('(texTest) Run Type %i: %s',stimulus.type, trialTypes{stimulus.type}));
disp(sprintf('-------------------------'))

%% Setup Screen
if stimulus.scan
  myscreen = initScreen('fMRIprojFlex');
else
  myscreen = initScreen('VPixx2');
end

% set background to grey
myscreen.background = 0.5;

%% Initialize Stimulus
myscreen = initStimulus('stimulus',myscreen);
localInitStimulus();
  
% Set response keys
stimulus.responseKeys = [1 2 3 4]; 

% set colors
stimulus.colors.white = [1 1 1];
stimulus.colors.black = [0 0 0];
stimulus.colors.red = [1 0 0];
stimulus.colors.green = [0 1 0];
stimulus.colors.blue = [0 0 1];
stimulus.live.fixColor = stimulus.colors.blue;
stimulus.live.cueColor = stimulus.colors.black;

%%%%%%%%%%%%% SETUP TEXTURE TASK %%%%%%%%%%%%%%%%%

stimulus.curTrial(1) = 0;

%% Define stimulus timing
stimRate = 4.5; % hertz
blockLen = 9; % seconds
stimulus.smpLen = 1.0 / stimRate; % seconds
nSegs = stimRate * blockLen;

% define the amount of time the stimulus should be on and off.
stimulus.tStimOn = 0.800;
stimulus.tStimOff = 0.200;

% Task important variables
stimulus.imNames = {'bark', 'rocks', 'spikes'}; %, 'im13', 'glass', 'bricks', 'branch'};
  
stimulus.layerNames = {'pool4'};
stimulus.poolSize = '1x1_';

stimulus.nTexFams = length(stimulus.imNames);
stimulus.imSize = 12;
stimulus.stimXPos = 7;
stimulus.nTexSmps = 1;
stimulus.nNoiseSmps = 1;
stimulus.nBkgdSmps = 25;
stimulus.nSmpsPerSeg = nSegs;

%% Select the condition for this run
% Choose which image and which pooling layer to display on this run on each side
stimulus.stimDir = '~/proj/TextureSynthesis/stimuli/texER/tex2';
stimulus.noiseDir = '~/proj/TextureSynthesis/stimuli/texER/noise2';
stimulus.bkgdDir = '~/proj/TextureSynthesis/stimuli/texER/bkgd';

%% Preload images
mask = imread('~/proj/TextureSynthesis/stimuli/Flattop8.tif');
stimulus.live.tex = struct();
stimulus.live.noise = struct();
disppercent(-inf, 'Preloading images');

% load texture and noise samples
for i = 1:stimulus.nTexFams
  imName = stimulus.imNames{i};
  for li = 1:length(stimulus.layerNames)
    layerI = stimulus.layerNames{li};
    for j = 1:stimulus.nTexSmps
      sd = rgb2gray(imread(sprintf('%s/%s%s_%s_smp%i.png', stimulus.stimDir, stimulus.poolSize, layerI, imName, j)));
      stimulus.live.tex.(sprintf('%s_%s_smp%i', layerI, imName, j)) = genTexFromIm(sd, mask);
    end
    for k = 1:stimulus.nNoiseSmps
      nd = rgb2gray(imread(sprintf('%s/noise_%s%s_%s_smp%i.png', stimulus.noiseDir, stimulus.poolSize, layerI, imName, k)));
      stimulus.live.noise.(sprintf('%s_%s_smp%i', layerI, imName, k)) = genTexFromIm(nd, mask);
    end
  end
  disppercent(i / stimulus.nTexFams);
end

% load background samples
for i = 1:stimulus.nBkgdSmps
  nd = rgb2gray(imread(sprintf('%s/bkgd_smp%i.png', stimulus.bkgdDir, i)));
  stimulus.live.bkgd.(sprintf('smp%i', i)) = genTexFromIm(nd, mask);
end

disppercent(inf);
clear sd nd

%%%%%%%%%%%%% TASK %%%%%%%%%%%%%%%%%
task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;
if stimulus.type==1 % Type 1: FAST event related
  task{1}{1}.segmin = [4.0 0.00]; % 3.75 sec trials
  task{1}{1}.segmax = [4.0 0.00];
  stimulus.blank = 1;
else % Type 2 and 3: Slow with long ISI
  task{1}{1}.segmin = [4.0 4.0]; % 7.5 sec trials
  task{1}{1}.segmax = [4.0 4.0]; % 7.5 sec trials
  if stimulus.type==2 % type 2: blank ISI
    stimulus.blank = 1;
  else % type 3: ISI with baseline/background images.
    stimulus.blank = 0;
  end
end
stimulus.seg = {};
stimulus.seg.stim = 1;
stimulus.seg.ITI = 2;

% Trial parameters
task{1}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
task{1}{1}.getResponse = zeros(size(task{1}{1}.segmin));
task{1}{1}.numTrials = 70;
task{1}{1}.random = 1;

if stimulus.scan
  task{1}{1}.synchToVol(end) = 1;
  % Shorten the last segment to account for synchtovol
  task{1}{1}.segmin(end) = max(0, task{1}{1}.segmin(end) - 0.200);
  task{1}{1}.segmax(end) = max(0, task{1}{1}.segmax(end) - 0.200);
end

% Specify task parameters
task{1}{1}.parameter.layer = 1:length(stimulus.layerNames);
task{1}{1}.parameter.texFam = 1:length(stimulus.imNames);
%task{1}{1}.parameter.sampleIdx = 1:(stimulus.nTexSmps + stimulus.nNoiseSmps); % 1-2 are textures, 3 is noise.
% Make textures twice as likely as noise.
task{1}{1}.parameter.sampleIdx = [repmat(1:stimulus.nTexSmps, 1,2), (stimulus.nTexSmps+1):(stimulus.nTexSmps+stimulus.nNoiseSmps)];

% Task variables to be calculated later
task{1}{1}.randVars.calculated.noiseOrTex = {NaN};
task{1}{1}.randVars.calculated.tSegStart = {NaN};

for phaseNum = 1:length(task{1})
  [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the third task to be the fixation staircase task
[task{2} myscreen] = fixStairInitTask(myscreen);

%% EYE CALIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%myscreen = eyeCalibDisp(myscreen);

%% Main Task Loop

mglClearScreen(0.5); 
upFix(stimulus);
mglFlush;
mglClearScreen(0.5); 
upFix(stimulus);

phaseNum = 1;
% Again, only one phase.
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
  mglClearScreen;
  % update the task
  [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,1);

  % update fixation task
  [task{2}, myscreen, phaseNum] = updateTask(task{2},myscreen,1);
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

stimulus.live.gotResponse = 0;
stimulus.curTrial(task.thistrial.thisphase) = stimulus.curTrial(task.thistrial.thisphase) + 1;

% Choose the 45 stimuli in this block by randomly sampling with replacement.
task.thistrial.tSegStart = [];

% At the start of each trial, choose which image to display.
if task.thistrial.sampleIdx > stimulus.nTexSmps % 1-4 are textures, 5-7 is noise
  task.thistrial.noiseOrTex = 'noise';
  sampleIdx = task.thistrial.sampleIdx - stimulus.nTexSmps;
else
  task.thistrial.noiseOrTex = 'tex';
  sampleIdx = task.thistrial.sampleIdx;
end
trialTexFam = stimulus.imNames{task.thistrial.texFam};
trialLayer = stimulus.layerNames{task.thistrial.layer};
stimulus.live.trialStim = stimulus.live.(task.thistrial.noiseOrTex).(sprintf('%s_%s_smp%i', trialLayer, trialTexFam, sampleIdx));

disp(sprintf('Trial %i - TexFam: %s, %s %s s%i', task.trialnum, trialTexFam, trialLayer, task.thistrial.noiseOrTex, sampleIdx));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = startSegmentCallback(task, myscreen)

global stimulus

% Save segment start time;
task.thistrial.tSegStart(task.thistrial.thisseg) = mglGetSecs;

%for i = 1:2
%  mglClearScreen;
%
%  if task.thistrial.thisseg == stimulus.seg.stim
%    mglBltTexture(stimulus.live.trialStim, [stimulus.stimXPos, 0, stimulus.imSize, stimulus.imSize]);
%    mglBltTexture(stimulus.live.trialStim, [-stimulus.stimXPos, 0, stimulus.imSize, stimulus.imSize]);
%  end
%
%  mglFlush;
%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus

% Select which stimulus to display as a function of time since seg start
timeSinceSegStart = mglGetSecs(task.thistrial.tSegStart(task.thistrial.thisseg));

% Flash stimuli on and off.
cycleLen = stimulus.tStimOn + stimulus.tStimOff;
stimOn = mod(timeSinceSegStart, cycleLen) < stimulus.tStimOn;

stimIdx = ceil(timeSinceSegStart / stimulus.smpLen);
stimIdx = ceil(timeSinceSegStart / cycleLen);

if task.thistrial.thisseg== stimulus.seg.stim
  if stimOn %mod(stimIdx,2) == 1
    mglBltTexture(stimulus.live.trialStim, [stimulus.stimXPos, 0, stimulus.imSize, stimulus.imSize]);
    mglBltTexture(stimulus.live.trialStim, [-stimulus.stimXPos, 0, stimulus.imSize, stimulus.imSize]);
  end
elseif task.thistrial.thisseg==stimulus.seg.ITI
  if stimOn && ~stimulus.blank %% mod(stimIdx,2) == 1
    thisBkgd = stimulus.live.bkgd.(sprintf('smp%i', min(25,stimIdx)));
    mglBltTexture(thisBkgd, [stimulus.stimXPos, 0, stimulus.imSize, stimulus.imSize]);
    mglBltTexture(thisBkgd, [-stimulus.stimXPos, 0, stimulus.imSize, stimulus.imSize]);
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

%%% Draws a cross at center of the screen of color fixColor
function upFix(stimulus, fixColor)
if ieNotDefined('fixColor')
  fixColor = stimulus.live.fixColor;
end
mglFixationCross(1,1,fixColor);

%%% 
% Draws a circular cue at location x,y
function drawCue(x,y, stimulus)
mglGluAnnulus(x,y, 0.75, 0.8, stimulus.live.cueColor, 64);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Turns image into a texture
function tex = genTexFromIm(im, mask)
r = flipud(im);

% Resize images to 256
if size(r,1) ~= 256;
  r = imresize(r, 256/size(r,1));
end

% Make sure they have 3 dimensions (even if grayscale)
if size(r,3) == 1
  r = cat(3, r, r, r);
end

% If a mask is passed in, apply as an alpha mask.
if ieNotDefined('mask')
  r(:,:,4) = 255;
else
  r(:,:,4) = mask(:,:,1);
end
% mgl requires the first dimension to be the RGBA dimension
rP = permute(r, [3 2 1]);
tex = mglCreateTexture(rP);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% totalTrials -- computes # of total trials
%
function [trials] = totalTrials()
%%
% Counts trials + estimates the threshold based on the last 500 trials
% get the files list
files = dir(fullfile(sprintf('~/data/texTest/%s/18*stim*.mat',mglGetSID)));
trials = 0;

for fi = 1:length(files)
    load(fullfile(sprintf('~/data/texTest/%s/%s',mglGetSID,files(fi).name)));
    e = getTaskParameters(myscreen,task);
    e = e{1}; % why?!
    trials = trials + e.nTrials;
end

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

% we'll adjust the gamma table to control contrast
stimulus.live.grating  = mglCreateTexture(alphamask); % high contrast        

