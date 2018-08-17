function [ myscreen ] = textureLoc( varargin )
%
% TEXTURE LOCALIZER
%  Experiment to map neural responses to various textures
%
%  Usage: textureLoc(varargin)
%  Authors: Akshay Jagadeesh
%  Date: 08/15/2018
%

global stimulus

stimulus = struct;

%% Initialize Variables

% add arguments later
scan = 0;
run = 0;
getArgs(varargin,{'scan=1', 'run=1'});
stimulus.scan = scan;
stimulus.run = run;
clear scan run;

%% Stimulus parameters 
%% Open Old Stimfile
stimulus.counter = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/textureLoc/%s',mglGetSID))
  % Directory exists, check for a stimfile
  files = dir(sprintf('~/data/textureLoc/%s/1*mat',mglGetSID));

  if length(files) >= 1
    fname = files(end).name;
    
    s = load(sprintf('~/data/textureLoc/%s/%s',mglGetSID,fname));
    stimulus.counter = s.stimulus.counter + 1;
    clear s;
    disp(sprintf('(textureLoc) Data file: %s loaded.',fname));
  end
end
disp(sprintf('(textureLoc) This is run #%i',stimulus.run));

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

%% Setup Task

%%%%%%%%%%%%% PHASE ONE %%%%%%%%%%%%%%%%%

stimulus.curTrial(1) = 0;

task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;

%% Define stimulus timing
stimRate = 5; % hertz
blockLen = 9; % seconds
stimLen = 1.0 / stimRate; % seconds
nSegs = stimRate * blockLen;

% task waits for fixation on first segment
task{1}{1}.segmin = repmat(stimLen, [1 nSegs]);
task{1}{1}.segmax = repmat(stimLen, [1 nSegs]);
stimulus.seg = {};

% Task important variables
stimulus.imNames = {'bark', 'branch', 'bricks', 'cracks', 'drops', 'floor', 'glass', 'rocks', 'spikes', 'wood'};
stimulus.layerNames = {'pool1', 'pool4'};
stimulus.rfNames = {'1x1'};
stimulus.stimDir = '~/proj/TextureSynthesis/out_bw';
stimulus.noiseDir = '~/proj/TextureSynthesis/spectral_noise';
stimulus.imSize = 12;
stimulus.nSegsPerBlock = nSegs;

% Choose which image and which pooling layer to display on this run
[a,b] = meshgrid(1:length(stimulus.imNames), 1:length(stimulus.layerNames));
stimulus.conditions = reshape(cat(2,a',b'), [], 2);
imageIndex = stimulus.conditions(stimulus.run, 1);
layerIndex = stimulus.conditions(stimulus.run, 2);
stimulus.runImage = stimulus.imNames{imageIndex};
stimulus.runLayer = stimulus.layerNames{layerIndex};

disp(sprintf('(textureLoc) Run #%i: Image = %s, Layer = %s', stimulus.run, stimulus.runImage, stimulus.runLayer));

% Trial parameters
task{1}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
task{1}{1}.getResponse = zeros(size(task{1}{1}.segmin));
%task{1}{1}.getResponse(stimulus.seg{1}.search)=1;
task{1}{1}.numTrials = 20;

if stimulus.scan
  task{1}{1}.synchToVol(nSegs) = 1;
  % Shorten the last segment to account for synchtovol
  task{1}{1}.segmin(end) = 0.100;
  task{1}{1}.segmax(end) = 0.100;
end

% Task trial parameters

% Task variables to be calculated later
task{1}{1}.randVars.calculated.noiseOrTex = NaN;
task{1}{1}.randVars.calculated.whichStimVersions = {'v1', 'v2'};

%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{1})
  [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the second task to be the fixation staircase task
%[task{2} myscreen] = fixStairInitTask(myscreen);


%% EYE CALIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%myscreen = eyeCalibDisp(myscreen);

%% Main Task Loop

mglClearScreen(0.5); 
upFix(stimulus);

mglFlush
mglClearScreen(0.5); 
upFix(stimulus);

phaseNum = 1;
% Again, only one phase.
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
  % update the task
  [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);

  % update fixation task
%  [task{2} myscreen] = updateTask(task{2},myscreen,1);
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

% directories
stimDirs = {'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7', 'v8', 'v9', 'v10'};

% Choose the 45 stimuli in this block by randomly sampling with replacement.
task.thistrial.whichStimVersions = randi(length(stimDirs), 1, stimulus.nSegsPerBlock);
stimulus.live.stimVersions = {stimDirs{task.thistrial.whichStimVersions}};

% Select noise or texture and display trial parameters
if mod(task.trialnum,2)==0 % even trials = spectral noise
  stimulus.live.trialType = 'spectral_noise';
  task.thistrial.noiseOrTex = 0;
  disp(sprintf('Block %i - Type: %s, Image: %s', task.trialnum, stimulus.live.trialType, stimulus.runImage));
  stimulus.live.stimPaths = cellfun(@(c)[stimulus.noiseDir '/' c '/noise_' stimulus.runImage '.jpg'], stimulus.live.stimVersions, 'uni', false);
else % odd trials = textures
  stimulus.live.trialType = 'texture';
  task.thistrial.noiseOrTex = 1;
  disp(sprintf('Block %i - Type: %s, Image: %s, Layer: %s', task.trialnum, stimulus.live.trialType, stimulus.runImage, stimulus.runLayer));
  stimulus.live.stimPaths = cellfun(@(c)[stimulus.stimDir '/' c '/1x1_' stimulus.runLayer '_' stimulus.runImage '.jpg'], stimulus.live.stimVersions, 'uni', false);
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

im = imread(stimulus.live.stimPaths{task.thistrial.thisseg});
tex = genTexFromIm(im);
clear im

for i = 1:2
  mglClearScreen(0.5);
  mglBltTexture(tex, [0 0 stimulus.imSize stimulus.imSize]);

  upFix(stimulus);
  mglFlush;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus

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
if size(r,3) == 1
  r = cat(3, r, r, r);
end
r(:,:,4) = 255;
rP = permute(r, [3 2 1]);
tex = mglCreateTexture(rP);  

function [trials] = totalTrials()
%%
% Counts trials + estimates the threshold based on the last 500 trials
% get the files list
files = dir(fullfile(sprintf('~/data/textureLoc/%s/18*stim*.mat',mglGetSID)));
trials = 0;

for fi = 1:length(files)
    load(fullfile(sprintf('~/data/textureLoc/%s/%s',mglGetSID,files(fi).name)));
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

