function [ myscreen ] = texobj( varargin )
%
% EVENT-RELATED TEXTURES
%  Experiment to map neural responses to various textures
%
%  Usage: texobj(varargin)
%  Authors: Akshay Jagadeesh
%  Date: 11/05/2018
%

global stimulus

stimulus = struct;

%% Initialize Variables

% add arguments later
getArgs(varargin,{'scan=0', 'testing=0'}, 'verbose=1');
stimulus.scan = scan;
stimulus.debug = testing;
clear scan testing;

%% Stimulus parameters 
%% Open Old Stimfile
stimulus.counter = 1;
%% Setup Screen
if stimulus.scan || true
  myscreen = initScreen('fMRIprojFlex');
else
  myscreen = initScreen('VPixx2');
end

% set background to grey
myscreen.background = 0.5;

%% Initialize Stimulus
myscreen = initStimulus('stimulus',myscreen);
 
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


stimulus.curTrial(1) = 0;

% define the amount of time the stimulus should be on and off.
stimulus.tStimOn = 0.800;
stimulus.tStimOff = 0.200;

% Task important variables
stimulus.imageNames = {'apple', 'bear', 'cat', 'cruiseship', 'dalmatian', 'ferrari',...
                    'greatdane', 'helicopter', 'horse', 'house', 'iphone', 'laptop',...
                    'quarterback', 'samosa', 'shoes', 'stephcurry', 'tiger', 'truck',...
                    'face', 'jetplane', 'elephant', 'sand', 'lawn', 'dirt', 'moss', 'leaves',...
                    'stars', 'tiles', 'worms', 'bumpy', 'spiky', 'clouds', 'crowd', 'forest', 'frills'};
stimulus.imageNames = {'bear', 'cat', 'horse', 'truck', 'face', 'helicopter', 'house', 'lawn', 'moss', 'dirt'};
stimulus.layerNames = {'pool4'};
stimulus.poolSize = '1x1_';

stimulus.nTexFams = length(stimulus.imageNames);
stimulus.imSize = 16;
stimulus.stimXPos = 10;
stimulus.num_samples = 2;

%% Select the condition for this run
% Choose which image and which pooling layer to display on this run on each side
stimulus.stimDir = '~/proj/texture_stimuli/color/texobj';
stimulus.origDir = '~/proj/texture_stimuli/color/originals';

%% Preload images
mask = imread('~/proj/TextureSynthesis/stimuli/Flattop8.tif');
stimulus.images.synths = struct();
stimulus.images.origs = struct();
disppercent(-inf, 'Preloading images');

% load texture and noise samples
for i = 1:stimulus.nTexFams
  imName = stimulus.imageNames{i};
  for li = 1:length(stimulus.layerNames)
    layerI = stimulus.layerNames{li};
    for j = 1:stimulus.num_samples
       sd = imread(sprintf('%s/%s%s_%s_smp%i.png', stimulus.stimDir, stimulus.poolSize, layerI, imName, j));
       stimulus.images.synths.(sprintf('%s_%s_%i', imName, layerI, j)) = genTexFromIm(sd);
    end
 
    % Load noise samples.
    orig = imread(sprintf('%s/%s.png', stimulus.origDir, imName));
    stimulus.images.origs.(imName) = genTexFromIm(orig);
  end
  
  disppercent(i / stimulus.nTexFams);
end

disppercent(inf);
clear sd nd

%%%%%%%%%%%%% TASK %%%%%%%%%%%%%%%%%
task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;

% Stimulus Timing: FAST event related
task{1}{1}.segmin = [4.0]; % 4 sec trials
task{1}{1}.segmax = [4.0];

if stimulus.debug==1
  task{1}{1}.segmin = [.10]; % 4 sec trials
  task{1}{1}.segmax = [.10];
end

stimulus.blank = 1;

stimulus.seg = {};
stimulus.seg.stim = 1;
%stimulus.seg.ITI = 2;

% Trial parameters
task{1}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
task{1}{1}.getResponse = zeros(size(task{1}{1}.segmin));
task{1}{1}.numTrials = 90;
task{1}{1}.random = 1;
task{1}{1}.seglenPrecomputeSettings.framePeriod=1.0;
if stimulus.scan
  task{1}{1}.synchToVol(end) = 1;
  % Shorten the last segment to account for synchtovol
  task{1}{1}.segmin(end) = max(0, task{1}{1}.segmin(end) - 0.200);
  task{1}{1}.segmax(end) = max(0, task{1}{1}.segmax(end) - 0.200);
end

% Initialize task parameters
task{1}{1}.parameter.imageClass = 1:length(stimulus.imageNames);
task{1}{1}.parameter.stimType = [-1, 0, 0, 1, 1]; % -1 = blank, 0 = orig, 1:N = synths layers.
task{1}{1}.parameter.stimXPos = stimulus.stimXPos;

task{1}{1}.randVars.calculated.blank = NaN;
task{1}{1}.randVars.calculated.layer = NaN;
task{1}{1}.randVars.calculated.tSegStart = {NaN};
task{1}{1}.randVars.calculated.sample = NaN;


%%% 
% Create a second task on the left side.
task{2}{1} = task{1}{1};
task{2}{1}.parameter.stimXPos = -stimulus.stimXPos;


for i = 1:length(task{1})
  [task{1}{i}, myscreen] = initTask(task{1}{i},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);
end

for i = 1:length(task{2})
  [task{2}{i}, myscreen] = initTask(task{2}{i},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the third task to be the fixation staircase task
[task{3} myscreen] = fixStairInitTask(myscreen);

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
  [task{3}, myscreen, phaseNum] = updateTask(task{3}, myscreen, 1);
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

if task.thistrial.stimXPos < 0
  trial_str = sprintf('Trial %i (LEFT)', task.trialnum);
  side = 'left';
else
  trial_str = sprintf('Trial %i (RIGHT)', task.trialnum);
  side = 'right';
end
image = stimulus.imageNames{task.thistrial.imageClass};

%% Determine whether this is a blank trial
if task.thistrial.stimType == -1
  task.thistrial.blank = 1;
  disp(sprintf('%s: Blank', trial_str));
elseif task.thistrial.stimType==0
  task.thistrial.blank = 0;
  task.thistrial.layer = -1;
  task.thistrial.sample = -1;
  stimulus.live.(sprintf('%s_trialStim', side)) = stimulus.images.origs.(image);
  disp(sprintf('%s: Original %s image', trial_str, image));
else
  task.thistrial.blank = 0;
  task.thistrial.layer = task.thistrial.stimType;
  task.thistrial.sample = randi(stimulus.num_samples);
  layer = stimulus.layerNames{task.thistrial.layer};
  stimulus.live.(sprintf('%s_trialStim', side)) = stimulus.images.synths.(sprintf('%s_%s_%i', image, layer, task.thistrial.sample));
  disp(sprintf('%s: %s %s sample %i synth', trial_str, image, layer, task.thistrial.sample));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = startSegmentCallback(task, myscreen)

global stimulus

% Save segment start time;
task.thistrial.tSegStart(task.thistrial.thisseg) = mglGetSecs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus

if task.thistrial.stimXPos < 0
  side = 'left';
else
  side = 'right';
end

% Select which stimulus to display as a function of time since seg start
timeSinceSegStart = mglGetSecs(task.thistrial.tSegStart(task.thistrial.thisseg));

% Flash stimuli on and off - 800ms on, 200ms off.
cycleLen = stimulus.tStimOn + stimulus.tStimOff;
stimOn = mod(timeSinceSegStart, cycleLen) < stimulus.tStimOn;

%stimIdx = ceil(timeSinceSegStart / stimulus.smpLen);
stimIdx = ceil(timeSinceSegStart / cycleLen);

% Draw the stimuli at the correct flicker rate.
if task.thistrial.thisseg== stimulus.seg.stim
  if stimOn && ~task.thistrial.blank
    mglBltTexture(stimulus.live.(sprintf('%s_trialStim', side)), [task.thistrial.stimXPos, 0, stimulus.imSize, stimulus.imSize]);
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
files = dir(fullfile(sprintf('~/data/texobj/%s/18*stim*.mat',mglGetSID)));
trials = 0;

for fi = 1:length(files)
    load(fullfile(sprintf('~/data/texobj/%s/%s',mglGetSID,files(fi).name)));
    e = getTaskParameters(myscreen,task);
    e = e{1}; % why?!
    trials = trials + e.nTrials;
end


