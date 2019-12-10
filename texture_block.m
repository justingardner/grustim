function [ myscreen ] = texture_block( varargin )
%
% TEXTURE - Block Design
%  Stimulus code for presenting textures and phase-scrambles in a blocked design.
%
%  Usage: texture_block(varargin)
%  Args:
%     - 'image=IMAGENAME' or 'image=1': Specifies which image class to use on this trial.
%          --> Arg can either be an image name (e.g. 'rocks') or a number to index into stimulus.imageNames
%     - 'textureType': Specifies which type of texture to use 
%              (default: 'PS' -- portilla-simoncelli). Other options: 'pool1', 'pool2', 'pool4'
%     - 
%  Example Usage:
%      myscreen = texture_block('image=1');
%      myscreen = texture_block('image', 'rocks');
%      myscreen = texture_block('image=rocks', 'textureType=pool2');
%  Authors: Akshay Jagadeesh
%  Date: 08/15/2018
%

global stimulus

stimulus = struct;

%% Initialize Variables
scan=0;
img=1;
textureType='ps';
getArgs(varargin,{'scan=0', 'img=1', 'textureType=ps'});
stimulus.scan = scan;

%% Select which stimuli to show on this trial
stimulus.imNames = {'beans', 'blossoms', 'bubbly', 'clouds', 'crystals', 'dahlias', 'fronds', 'fur', 'glass', 'leaves', 'leopard', 'noodles', 'paisley', 'plant', 'rocks', 'scales', 'spikes', 'tiles', 'waves', 'worms'};
stimulus.layerNames = {'ps', 'pool2', 'pool4'};

% Set which image class we will be showing on this trial.
if isnumeric(img) && img>=1 && img<=length(stimulus.imNames)
  stimulus.runImage = stimulus.imNames{img};
elseif ischar(img)
  stimulus.runImage = img;
else
  disp(sprintf('Input for image is of the wrong type. It must be either a number (between 1 and %i) or a name of an image.', length(stimulus.imNames)));
  return;
end

% Set which texture type (Portilla-Simoncelli, Pool1, Pool2, or Pool4) we will be showing on this trial.
if ~any(strcmp(stimulus.layerNames, textureType))
  disp(['Input for textureType is invalid. Must be one of the following: ', sprintf('%s ', stimulus.layerNames)]);
  return
else
  stimulus.runTextureType = textureType;
end

disp(sprintf('---TEXTURE_BLOCK: Starting run - Texture Type = %s, Image = %s---', stimulus.runTextureType, stimulus.runImage));

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

%%%%%%%%%%%%% SETUP TEXTURE PRESENTATION %%%%%%%%%%%%%%%%%

stimulus.curTrial(1) = 0;

% Define stimulus timing vars
stimRate = 5; % hertz (frequency of changing samples)
blockLength = 16; % seconds (Length of each block)
stimulus.smpLen = 1.0 / stimRate; % seconds (length of time to flash each stimulus for)

% Define stimulus presentations vars
stimulus.imSize = 16;
stimulus.nSmps = 15; % number of unique samples that we have
stimulus.nSamplesPerBlock = stimRate * blockLength; % Number of samples that get shown in each block = number of seconds per block * number of samples per second

% Specify which directory to get the stimuli from 
if strcmp(stimulus.runTextureType, 'ps')
  stimulus.stimDir = '~/proj/TextureSynthesis/stimuli/fzs/ps_tex';
  stimulus.phasescrambleDir = '~/proj/TextureSynthesis/stimuli/fzs/ps_noise';
else
  stimulus.stimDir = '~/proj/TextureSynthesis/stimuli/fzs/tex';
  stimulus.phasescrambleDir = '~/proj/TextureSynthesis/stimuli/fzs/noise';
end

%% Preload images
mask = imread('~/proj/TextureSynthesis/stimuli/Flattop8.tif');
stimulus.live.texture = struct();
stimulus.live.phasescramble = struct();
disppercent(-inf, 'Preloading images');
for i = 1:stimulus.nSmps
  sd1 = imread(sprintf('%s/%s_%s_smp%i.png', stimulus.stimDir, stimulus.runTextureType, stimulus.runImage, i));
  nd1 = imread(sprintf('%s/noise_%s_%s_smp%i.png', stimulus.phasescrambleDir, stimulus.runTextureType, stimulus.runImage, i));

  stimulus.live.texture.(sprintf('smp%i', i)) = genTexFromIm(sd1, mask);
  stimulus.live.phasescramble.(sprintf('smp%i', i)) = genTexFromIm(nd1, mask);
  disppercent(i / 10);
end
disppercent(inf);
clear sd1 nd1 

%%%%%%%%%%%%% MGL TASK SETUP %%%%%%%%%%%%%%%%%
task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;
task{1}{1}.segmin = [blockLength .01];
task{1}{1}.segmax = [blockLength .01];
stimulus.seg = {};
stimulus.seg.stim = 1;
stimulus.seg.iti = 2;

% Trial parameters
task{1}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
task{1}{1}.getResponse = zeros(size(task{1}{1}.segmin));
task{1}{1}.numTrials = 20;

if stimulus.scan
  task{1}{1}.synchToVol(end) = 1;
  % Shorten the last segment to account for synchtovol
  task{1}{1}.segmin(1) = blockLength-0.2;
  task{1}{1}.segmax(1) = blockLength-0.2;
end

% Specify task parameters
task{1}{1}.parameter.stimXPos = 0;

% Task variables to be calculated later
task{1}{1}.randVars.uniform.trialType = [0,1,2]; % 0: blank, 1: texture, 2: phase-scramble.
task{1}{1}.randVars.calculated.trialTypeStr = {NaN};
task{1}{1}.randVars.calculated.trialSamples = {NaN};
task{1}{1}.randVars.calculated.tSegStart = {NaN};

for phaseNum = 1:length(task{1})
  [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize a second task to be the fixation staircase task
[task{2} myscreen] = fixStairInitTask(myscreen);


%% Main Task Loop
mglClearScreen(0.5); 
upFix(stimulus);

mglFlush
mglClearScreen(0.5); 
upFix(stimulus);

phaseNum = 1;
% Again, only one phase.
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
  mglClearScreen;
  % update the task
  [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);

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

%%%%%%%%%% END OF MAIN TASK CODE: Callback / Helper Functions Follow %%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% startTrialCallback %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
function [task, myscreen] = startTrialCallback(task,myscreen)

global stimulus

stimulus.live.gotResponse = 0;
stimulus.curTrial(task.thistrial.thisphase) = stimulus.curTrial(task.thistrial.thisphase) + 1;

% Choose the 45 stimuli in this block by randomly sampling with replacement.
task.thistrial.trialSamples = randi(stimulus.nSmps, 1, stimulus.nSamplesPerBlock);
task.thistrial.tSegStart = [];
task.thistrial.phasescrambleOrTex = {};

% Determine which trial type (0: blank, 1: texture, or 2: phase-scramble).
if task.thistrial.trialType == 1
  stimulus.live.trialStim = stimulus.live.texture;
  task.thistrial.trialTypeStr = 'texture';
elseif task.thistrial.trialType==2
  stimulus.live.trialStim = stimulus.live.phasescramble;
  task.thistrial.trialTypeStr = 'phasescramble';
elseif task.thistrial.trialType==0
  task.thistrial.trialTypeStr = 'blank';
end

disp(sprintf('Trial %i: %s', task.trialnum, task.thistrial.trialTypeStr));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% startSegmentCallback %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = startSegmentCallback(task, myscreen)

global stimulus
  
% Save segment start time;
task.thistrial.tSegStart(task.thistrial.thisseg) = mglGetSecs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% screenUpdateCallback %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called on every screen refresh %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = screenUpdateCallback(task, myscreen)

global stimulus

%% Select which stimulus to display as a function of time since seg start
if task.thistrial.thisseg == stimulus.seg.stim && ~strcmp(task.thistrial.trialTypeStr, 'blank')
  % How much time has passed since the start of the segment?
  timeFromBlockStart = mglGetSecs(task.thistrial.tSegStart(task.thistrial.thisseg));

  % Which sample are we showing here?
  stimIdx = min(ceil(timefromBlockStart / stimulus.smpLen), length(task.thistrial.trialSamples));
  smp = sprintf('smp%i', task.thistrial.trialSamples(stimIdx));
  trialStim = stimulus.live.trialStim.(smp);

  % Display the texture to the screen.
  mglBltTexture(trialStim, [task.thistrial.stimXPos, 0, stimulus.imSize, stimulus.imSize]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% getResponseCallback %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus

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

%%% Loads image and turns into a mgl texture, which can then be Blt'd quickly
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

