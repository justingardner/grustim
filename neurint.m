function [ myscreen ] = neurint(varargin)
%
% Neural Interpolations
%  How well can observers discriminate along a specific dimension through
%  V4, IT or dCNN representational space?
%
%  Usage: neurint(varargin)
%  Authors: Akshay Jagadeesh
%  Date: 11/05/2018
%

global stimulus

stimulus = struct;

%% Initialize Variables

% add arguments later
scan = 1;
getArgs(varargin,{'scan=0', 'testing=0'}, 'verbose=1');
stimulus.scan = scan;
stimulus.debug = testing;
clear scan testing;

%% Stimulus parameters 
%% Open Old Stimfile
stimulus.counter = 1;


%% Setup Screen
if stimulus.scan
  myscreen = initScreen('fMRIprojFlex');
else
  myscreen = initScreen('VPixx2');
end

%%
if stimulus.debug == 0
  myscreen.saveData = 1;
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

% Task important variables
stimulus.intervals = [0, 25, 50, 75, 100];
stimulus.imSize = 12;
stimulus.stimXPos = 8;
stimulus.num_samples = 4;
stimulus.feature_spaces = {'pool3-V1', 'pool4-IT'};
stimulus.image_pairs = {'elephant1_face1', 'horse_cat', 'rocks_leaves'};

%% Select the condition for this run
% Choose which image and which pooling layer to display on this run on each side
stimulus.stimDir = '~/proj/NeuralImageSynthesis/neurint/bw_new';
stimulus.stimDir = '~/proj/neurint/outputs/bw';

%% Preload images
mask = imread('~/proj/NeuralImageSynthesis/Flattop8.tif');
stimulus.images= struct();
disppercent(-inf, 'Preloading images');

% load texture and noise samples
for fi = 1:length(stimulus.feature_spaces)
  for ii = 1:length(stimulus.image_pairs)
    image_pair = stimulus.image_pairs{ii};
    feature_space = strrep(stimulus.feature_spaces{fi}, '-', '_');
    for i = 1:length(stimulus.intervals)
      interval = stimulus.intervals(i);
      for j = 1:stimulus.num_samples
        image = imread(sprintf('%s/%s_%i_%s_s%i.png', stimulus.stimDir, image_pair, interval, stimulus.feature_spaces{fi}, j));
        stimulus.images.(sprintf('%s_%s_interval%i_smp%i', image_pair, feature_space, interval, j)) = genTexFromIm(image, mask);
      end
    end
  end
  stimulus.feature_spaces{fi} = feature_space;
  disppercent(fi / length(stimulus.feature_spaces));
end
disppercent(inf);
clear image

%%%%%%%%%%%%% TASK %%%%%%%%%%%%%%%%%
task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;
task{1}{1}.segmin = [0.2, 0.2, 0.2, 5.0, 0.2];
task{1}{1}.segmax = [0.2, 0.2, 0.2, 5.0, 0.2];

stimulus.seg = {};
stimulus.seg.fix = 1;
stimulus.seg.stim1 = 2;
stimulus.seg.isi = 3;
stimulus.seg.stim2 = 4;
stimulus.seg.feedback = 5;

% Trial parameters
task{1}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
task{1}{1}.getResponse = zeros(size(task{1}{1}.segmin));
task{1}{1}.getResponse(stimulus.seg.stim2) = 1;

task{1}{1}.numTrials = 250;
task{1}{1}.random = 1;

if stimulus.scan
  task{1}{1}.synchToVol(end) = 1;
  % Shorten the last segment to account for synchtovol
  task{1}{1}.segmin(end) = max(0, task{1}{1}.segmin(end) - 0.200);
  task{1}{1}.segmax(end) = max(0, task{1}{1}.segmax(end) - 0.200);
end

% Initialize task parameters

% Assign layer, texture family, and sample index in random blocks.
task{1}{1}.parameter.img1Sample = 1:stimulus.num_samples;
task{1}{1}.parameter.img2Sample = 1:stimulus.num_samples;
task{1}{1}.parameter.centerInterval = stimulus.intervals;
task{1}{1}.parameter.image_pair = 1:length(stimulus.image_pairs);
task{1}{1}.parameter.leftImg = [1 2]; % 1 = img1, 2=img2;
task{1}{1}.randVars.uniform.feature_space = 1:length(stimulus.feature_spaces);
task{1}{1}.randVars.calculated.correct = NaN;
task{1}{1}.randVars.calculated.response = NaN;
task{1}{1}.randVars.calculated.imgChoice = NaN;
task{1}{1}.randVars.calculated.rightImg = NaN;
task{1}{1}.randVars.calculated.centerSample = NaN;
task{1}{1}.randVars.calculated.leftSample = NaN;
task{1}{1}.randVars.calculated.rightSample = NaN;

for phaseNum = 1:length(task{1})
  [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the third task to be the fixation staircase task
% [task{2} myscreen] = fixStairInitTask(myscreen);

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
  % [task{2}, myscreen, phaseNum] = updateTask(task{2},myscreen,1);
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

feature_space = stimulus.feature_spaces{task.thistrial.feature_space};

if task.thistrial.centerInterval == 0
  task.thistrial.centerSample = randsample(setdiff(1:stimulus.num_samples, task.thistrial.img1Sample), 1);
elseif task.thistrial.centerInterval == 100
  task.thistrial.centerSample = randsample(setdiff(1:stimulus.num_samples, task.thistrial.img2Sample), 1);
else
  task.thistrial.centerSample = randsample(1:stimulus.num_samples, 1);
end

task.thistrial.rightImg = setdiff(1:2, task.thistrial.leftImg);
image_pair = stimulus.image_pairs{task.thistrial.image_pair};

if task.thistrial.leftImg == 1
  stimulus.live.leftStim = stimulus.images.(sprintf('%s_%s_interval0_smp%i', image_pair, feature_space, task.thistrial.img1Sample));
  stimulus.live.rightStim = stimulus.images.(sprintf('%s_%s_interval100_smp%i', image_pair, feature_space, task.thistrial.img2Sample));
  task.thistrial.leftSample = task.thistrial.img1Sample;
  task.thistrial.rightSample = task.thistrial.img2Sample;
  disp(sprintf('--- Trial %i - Left (0%%) Sample %i -- Center (%i%%), sample %i -- Right (100%%), sample %i ---', task.trialnum, task.thistrial.leftSample, task.thistrial.centerInterval, task.thistrial.centerSample, task.thistrial.rightSample));
else
  stimulus.live.leftStim = stimulus.images.(sprintf('%s_%s_interval100_smp%i', image_pair, feature_space, task.thistrial.img2Sample));
  stimulus.live.rightStim = stimulus.images.(sprintf('%s_%s_interval0_smp%i', image_pair, feature_space, task.thistrial.img1Sample));
  task.thistrial.leftSample = task.thistrial.img2Sample;
  task.thistrial.rightSample = task.thistrial.img1Sample;
  disp(sprintf('--- Trial %i - Left (100%%) Sample %i -- Center (%i%%), sample %i -- Right (0%%), sample %i ---', task.trialnum, task.thistrial.leftSample, task.thistrial.centerInterval, task.thistrial.centerSample, task.thistrial.rightSample));
end
stimulus.live.centerStim = stimulus.images.(sprintf('%s_%s_interval%i_smp%i', image_pair, feature_space, task.thistrial.centerInterval, task.thistrial.centerSample));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = startSegmentCallback(task, myscreen)

global stimulus

% Save segment start time;
task.thistrial.tSegStart(task.thistrial.thisseg) = mglGetSecs;
stimulus.live.screenGrab = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus
mglClearScreen(0.5);

if task.thistrial.thisseg == stimulus.seg.stim1
  mglBltTexture(stimulus.live.centerStim, [0,0, stimulus.imSize, stimulus.imSize]);
  if stimulus.live.screenGrab && stimulus.screenGrab
      frame=mglScreenGrab();
      
  end
elseif task.thistrial.thisseg == stimulus.seg.stim2
  mglBltTexture(stimulus.live.leftStim, [-stimulus.stimXPos, 0, stimulus.imSize, stimulus.imSize]);
  mglBltTexture(stimulus.live.rightStim, [stimulus.stimXPos, 0, stimulus.imSize, stimulus.imSize]);
 end

if task.thistrial.thisseg == stimulus.seg.feedback
  if ~isnan(task.thistrial.correct) && task.thistrial.correct==1
    upFix(stimulus, stimulus.colors.green);
  else
    upFix(stimulus, stimulus.colors.red);
  end
else
  upFix(stimulus);
end

%%%%
function upCue(task)

global stimulus

if task.thistrial.thisseg ~= stimulus.seg.resp && task.thistrial.cueFocal == 0
	drawArrow(0);
	drawArrow(pi);
elseif task.thistrial.cueSide == 1
	drawArrow(0);
else
	drawArrow(pi);
end

%%%
function drawArrow(theta)
x=[1,.75,2,.75];
y=[0,.5,0,-.5];
mglPolygon(cos(theta)*x - sin(theta)*y, sin(theta)*x + cos(theta)*y, [0,0,0])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus


validResponse = any(task.thistrial.whichButton == stimulus.responseKeys);

if validResponse
  if stimulus.live.gotResponse==0
    task.thistrial.detected = 1;
    task.thistrial.response = task.thistrial.whichButton;
    disp(sprintf('Response: %i', task.thistrial.response));
    if (task.thistrial.response == 1 && task.thistrial.leftImg == 1) || (task.thistrial.response == 2 && task.thistrial.leftImg==2)
      task.thistrial.correct = task.thistrial.centerInterval < 50;
      task.thistrial.imgChoice = 1;
    else
      task.thistrial.correct = task.thistrial.centerInterval >=50;
      task.thistrial.imgChoice = 2;
    end
    if task.thistrial.correct
    	disp('Correct!')
    else
    	disp('Incorrect')
    end
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

function imout = maxContrast(im)

imout = (im - min(im(:)))*(255/max(im(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Turns image into a texture
function tex = genTexFromIm(im, mask, contrast)

if ~ieNotDefined('contrast')
  im=maxContrast(im);
end

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
files = dir(fullfile(sprintf('~/data/neurint/%s/18*stim*.mat',mglGetSID)));
trials = 0;

for fi = 1:length(files)
    load(fullfile(sprintf('~/data/neurint/%s/%s',mglGetSID,files(fi).name)));
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

