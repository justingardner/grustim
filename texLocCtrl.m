function [ myscreen ] = texLocCtrl( varargin )
%
% TEXTURE LOCALIZER
%  Experiment to map neural responses to various textures
%
%  Usage: texLocCtrl(varargin)
%  Authors: Akshay Jagadeesh
%  Date: 08/15/2018
%

global stimulus

stimulus = struct;

%% Initialize Variables

% add arguments later
scan = 0;
run = 0;
getArgs(varargin,{'scan=0', 'run=1'});
stimulus.scan = scan;
stimulus.run = run;
clear scan run;

%% Stimulus parameters 
%% Open Old Stimfile
stimulus.counter = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/texLocCtrl/%s',mglGetSID))
  % Directory exists, check for a stimfile
  files = dir(sprintf('~/data/texLocCtrl/%s/1*mat',mglGetSID));

  if length(files) >= 1
    fname = files(end).name;
    
    s = load(sprintf('~/data/texLocCtrl/%s/%s',mglGetSID,fname));
    stimulus.counter = s.stimulus.counter + 1;
    clear s;
    disp(sprintf('(texLocCtrl) Data file: %s loaded.',fname));
  end
end
disp(sprintf('(texLocCtrl) This is run #%i',stimulus.run));

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
stimRate = 5; % hertz
blockLen = 9; % seconds
stimulus.smpLen = 1.0 / stimRate; % seconds

% Task important variables
a = {13, 71, 52, 48, 60, 18, 23, 30, 327, 336, 38, 393, 402, 56, 99};
stimulus.imNames = cellfun(@(x) sprintf('im%i', x), a, 'UniformOutput', 0);
stimulus.layerNames = {'ps', 'pool2', 'pool4'};
stimulus.imSize = 16;
stimulus.nSmps = 15;
stimulus.nSmpsPerSeg = stimRate * blockLen;

%% Select the condition for this run
% Choose which image and which pooling layer to display on this run on each side
[a,b] = meshgrid(1:length(stimulus.imNames), 1:length(stimulus.layerNames));
stimulus.conditions = reshape(cat(2,a',b'), [], 2);
imageIndex = stimulus.conditions(stimulus.run, 1);
layerIndex = stimulus.conditions(stimulus.run, 2);

stimulus.runImage = stimulus.imNames{imageIndex};
stimulus.runLayer = stimulus.layerNames{layerIndex};

disp(sprintf('(texLocCtrl) Run #%i: Layer = %s, LeftImg = %s', stimulus.run, stimulus.runLayer, stimulus.runImage));
stimulus.stimDir = {'~/proj/TextureSynthesis/stimuli/fzs/ps_tex', '~/proj/TextureSynthesis/stimuli/fzs/tex'};
stimulus.noiseDir = {'~/proj/TextureSynthesis/stimuli/fzs/ps_noise', '~/proj/TextureSynthesis/stimuli/fzs/noise'};

%% Preload images
mask = imread('~/proj/TextureSynthesis/stimuli/Flattop8.tif');
mask = mask(:,:,1);
largemask = imread('~/proj/TextureSynthesis/stimuli/inverseMask.tif');

stimulus.live.stim = struct();
stimulus.live.anti = struct();
disppercent(-inf, 'Preloading images');
ctr = 1;
for i = 1:stimulus.nSmps
  for j = 1:length(stimulus.layerNames)
    layerJ = stimulus.layerNames{j};
    for k = 1:length(stimulus.imNames)
      imgK = stimulus.imNames{k};
      if strcmp(layerJ, 'ps')
        stimDir = stimulus.stimDir{1};
        noiseDir = stimulus.noiseDir{1};
      else
        stimDir = stimulus.stimDir{2};
        noiseDir = stimulus.noiseDir{2};
      end

      sd1 = imread(sprintf('%s/%s_%s_smp%i.png', stimDir, layerJ, imgK, i));
      nd1 = imread(sprintf('%s/noise_%s_%s_smp%i.png', noiseDir, layerJ, imgK, i));

      stimulus.live.stim.(sprintf('smp%i', ctr)) = genTexFromIm(sd1, mask);
      stimulus.live.stim.(sprintf('smp%i', ctr+1)) = genTexFromIm(nd1, mask);

      sd2 = [sd1 fliplr(sd1); flipud(sd1) fliplr(flipud(sd1))];
      nd2 = [nd1 fliplr(nd1); flipud(nd1) fliplr(flipud(nd1))];
      stimulus.live.anti.(sprintf('smp%i', ctr)) = genTexFromIm(sd2, largemask);
      stimulus.live.anti.(sprintf('smp%i', ctr+1)) = genTexFromIm(nd2, largemask);
      ctr = ctr + 2;
    end
  end
  disppercent(i/stimulus.nSmps);
end
disppercent(inf);
clear sd1 nd1 

%%%%%%%%%%%%% TASK %%%%%%%%%%%%%%%%%
task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;
task{1}{1}.segmin = [9.00 9.00];
task{1}{1}.segmax = [9.00 9.00];
stimulus.seg = {};
stimulus.seg.stim = 1;
stimulus.seg.blank = 2;

% Trial parameters
task{1}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
task{1}{1}.getResponse = zeros(size(task{1}{1}.segmin));
task{1}{1}.numTrials = 20;

if stimulus.scan
  task{1}{1}.synchToVol(end) = 1;
  % Shorten the last segment to account for synchtovol
  task{1}{1}.segmin(end) = 8.800;
  task{1}{1}.segmax(end) = 8.800;
end

% Specify task parameters
task{1}{1}.parameter.stimXPos = 0;

% Task variables to be calculated later
task{1}{1}.randVars.calculated.stimOrBlank = {NaN};
task{1}{1}.randVars.calculated.whichStimVersions = {NaN};
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

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startTrialCallback(task,myscreen)

global stimulus

stimulus.live.gotResponse = 0;
stimulus.curTrial(task.thistrial.thisphase) = stimulus.curTrial(task.thistrial.thisphase) + 1;

% Choose the 45 stimuli in this block by randomly sampling with replacement.
task.thistrial.whichStimVersions = randi(length(fields(stimulus.live.stim)), 1, stimulus.nSmpsPerSeg);
task.thistrial.tSegStart = [];
task.thistrial.stimOrBlank = {};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)

global stimulus

if task.thistrial.thisseg == stimulus.seg.stim
  task.thistrial.stimOrBlank{task.thistrial.thisseg} = 'stim';
else
  task.thistrial.stimOrBlank{task.thistrial.thisseg} = 'blank';
end  
disp(sprintf('%i: %s block', task.trialnum, task.thistrial.stimOrBlank{task.thistrial.thisseg}));
% Save segment start time;
task.thistrial.tSegStart(task.thistrial.thisseg) = mglGetSecs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus

% Select which stimulus to display as a function of time since seg start
timeSinceSegStart = mglGetSecs(task.thistrial.tSegStart(task.thistrial.thisseg));

stimIdx = min(ceil(timeSinceSegStart / stimulus.smpLen), length(task.thistrial.whichStimVersions));
smp = sprintf('smp%i', task.thistrial.whichStimVersions(stimIdx));

if task.thistrial.thisseg == stimulus.seg.stim
  thisTex = stimulus.live.stim.(smp);
  mglBltTexture(thisTex, [0, 0, stimulus.imSize, stimulus.imSize]);
else
  thisTex = stimulus.live.anti.(smp);
  mglBltTexture(thisTex, [0, 0, stimulus.imSize*2, stimulus.imSize*2]);
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
if size(r,1) ~= size(mask,1)
  r = imresize(r, size(mask,1)/size(r,1));
end

% Make sure they have 3 dimensions (even if grayscale)
if size(r,3) == 1
  r = cat(3, r, r, r);
end

% If a mask is passed in, apply as an alpha mask.
if ieNotDefined('mask')
  r(:,:,4) = 255;
else
  r(:,:,4) = mask;
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
files = dir(fullfile(sprintf('~/data/texLocCtrl/%s/18*stim*.mat',mglGetSID)));
trials = 0;

for fi = 1:length(files)
    load(fullfile(sprintf('~/data/texLocCtrl/%s/%s',mglGetSID,files(fi).name)));
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
