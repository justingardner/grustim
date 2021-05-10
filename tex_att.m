function [ myscreen ] = tex_att( varargin )
%
% ATTENTION TEXTURES
%  Same-different task with focal vs distributed attentional cue.
%
%  Usage: tex_att(varargin)
%  Authors: Akshay Jagadeesh
%  Date: 04/01/2021
%

global stimulus

stimulus = struct;

%% Initialize Variables

% add arguments later
scan = 1;
getArgs(varargin,{'scan=0', 'testing=0', 'noeye=0'}, 'verbose=1');
stimulus.scan = scan;
stimulus.debug = testing;
stimulus.noeye = noeye;
clear scan testing

%% Stimulus parameters 
%% Open Old Stimfile
stimulus.counter = 1;

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

% Task important variables
%stimulus.imageNames = {'lawn', 'moss', 'dirt'};
stimulus.imageNames = {'blossoms', 'bricks', 'cherries', 'fronds', 'fur', 'grass', 'noodles', 'worms'};
stimulus.layerNames = {'pool4'};
stimulus.poolSize = '1x1_';

stimulus.nTexFams = length(stimulus.imageNames);
stimulus.imSize = 8;
stimulus.stimXPos = 6;
stimulus.stimYPos = 4;
stimulus.num_samples = 1;

%% Select the condition for this run
% Choose which image and which pooling layer to display on this run on each side
stimulus.stimDir = '~/proj/texture_stimuli/color/textures';
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
      if isfile( sprintf('%s/%s%s_%s_smp%i.png', stimulus.stimDir, stimulus.poolSize, layerI, imName, j))
        sd = imread(sprintf('%s/%s%s_%s_smp%i.png', stimulus.stimDir, stimulus.poolSize, layerI, imName, j));
      elseif isfile(sprintf('%s/%s%s_%s_smp%i.jpg', stimulus.stimDir, stimulus.poolSize, layerI, imName, j))
        sd = imread(sprintf('%s/%s%s_%s_smp%i.jpg', stimulus.stimDir, stimulus.poolSize, layerI, imName, j));
      else
        disp(sprintf('%s%s_%s_smp%i not found', stimulus.poolsize, layerI, imName, j));
        keyboard
      end
      stimulus.images.synths.(sprintf('%s_%s_%i', imName, layerI, j)) = genTexFromIm(sd, mask);
    end
 
    % Load noise samples.
    orig = imread(sprintf('%s/%s.jpg', stimulus.origDir, imName));
    stimulus.images.origs.(imName) = genTexFromIm(orig, mask);
  end
  
  disppercent(i / stimulus.nTexFams);
end

disppercent(inf);
clear sd

%%%%%%%%%%%%% TASK %%%%%%%%%%%%%%%%%
task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;

task{1}{1}.segmin = [inf, 1.0, 1.0, 2.2, 0.2, 0.1];
task{1}{1}.segmax = [inf, 1.0, 1.0, 2.2, 0.2, 0.1];

stimulus.seg = {};
stimulus.seg.fix = 1;
stimulus.seg.cue = 2;
stimulus.seg.stim1 = 3;
stimulus.seg.resp = 4;
stimulus.seg.feedback = 5;
stimulus.seg.ITI = 6;
if stimulus.noeye==1
  task{1}{1}.segmin(1) = 0.2;
  task{1}{1}.segmin(1) = 0.2;
end

% Trial parameters
task{1}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
task{1}{1}.getResponse = zeros(size(task{1}{1}.segmin));
task{1}{1}.getResponse(stimulus.seg.resp) = 1;

task{1}{1}.numTrials = 60;
task{1}{1}.random = 1;

if stimulus.scan
  task{1}{1}.synchToVol(end) = 1;
  % Shorten the last segment to account for synchtovol
  task{1}{1}.segmin(end) = max(0, task{1}{1}.segmin(end) - 0.200);
  task{1}{1}.segmax(end) = max(0, task{1}{1}.segmax(end) - 0.200);
end

% Initialize task parameters

% Assign layer, texture family, and sample index in random blocks.
task{1}{1}.parameter.leftImgClass = 1:length(stimulus.imageNames);
task{1}{1}.parameter.rightImgClass = 1:length(stimulus.imageNames);
task{1}{1}.parameter.cueFocal = [0,1]; %0 = distributed, 1 = focal
task{1}{1}.parameter.cueSide = [1 2];
task{1}{1}.randVars.uniform.leftSame = [0,1];
task{1}{1}.randVars.uniform.rightSame = [0,1];
task{1}{1}.randVars.calculated.correct = NaN;
task{1}{1}.randVars.calculated.response = NaN;
task{1}{1}.randVars.calculated.dead = 0;

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
myscreen = eyeCalibDisp(myscreen);

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

task.thistrial.dead = 0;
stimulus.live.gotResponse = 0;
stimulus.curTrial(task.thistrial.thisphase) = stimulus.curTrial(task.thistrial.thisphase) + 1;


% At the start of each trial, choose which image to display.
% 1. Randomly select a texture family
leftImgClass = stimulus.imageNames{task.thistrial.leftImgClass};
rightImgClass = stimulus.imageNames{task.thistrial.rightImgClass};

sampleIdx=1;

% Load the stimulus for this trial.
if task.thistrial.leftSame == 0
	stimulus.live.leftStims = {stimulus.images.synths.(sprintf('%s_pool4_1', leftImgClass)), stimulus.images.origs.(leftImgClass)};
else
	if rand() > 0.5
		stimulus.live.leftStims = {stimulus.images.origs.(leftImgClass), stimulus.images.origs.(leftImgClass)};
	else
		stimulus.live.leftStims = {stimulus.images.synths.(sprintf('%s_pool4_1', leftImgClass)), stimulus.images.synths.(sprintf('%s_pool4_1', leftImgClass))};
	end
end
if task.thistrial.rightSame == 0
	stimulus.live.rightStims = {stimulus.images.synths.(sprintf('%s_pool4_1', rightImgClass)), stimulus.images.origs.(rightImgClass)};
else
	if rand() > 0.5
		stimulus.live.rightStims = {stimulus.images.origs.(rightImgClass), stimulus.images.origs.(rightImgClass)};
	else
		stimulus.live.rightStims = {stimulus.images.synths.(sprintf('%s_pool4_1', rightImgClass)), stimulus.images.synths.(sprintf('%s_pool4_1', rightImgClass))};
	end
end

% if task.thistrial.rightSame == 1
% 	stimulus.live.rightStim1 = stimulus.images.synths.(sprintf('%s_pool4_1', leftImgClass));
% 	stimulus.live.rightStim2 = stimulus.images.synths.(sprintf())
% end
% stimulus.live.trialStim = 
stimulus.live.eyeCount = 0;
cueSide = {'Right', 'Left'};
sameDiff = {'Different', 'Same'};
disp(sprintf('--- Trial %i - Cue Side: %s; Left: %s %s, Right: %s %s ---', task.trialnum, cueSide{task.thistrial.cueSide},...
                                   leftImgClass, sameDiff{task.thistrial.leftSame+1}, rightImgClass, sameDiff{1+task.thistrial.rightSame}));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = startSegmentCallback(task, myscreen)

global stimulus

% Save segment start time;
task.thistrial.tSegStart(task.thistrial.thisseg) = mglGetSecs;
stimulus.live.eyeDead = 0;
stimulus.live.triggerWaiting = 0;
if any(task.thistrial.thisseg==[stimulus.seg.fix])
  stimulus.live.triggerWaiting = 1;
  stimulus.live.centered = 0;
  stimulus.live.triggerTime = 0;
  stimulus.live.lastTrigger = -1;
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
if ~stimulus.noeye && ~any(task.thistrial.thisseg==[stimulus.seg.fix]) 
  if ~any(isnan(pos))
    if dist > 1.5 && stimulus.live.eyeCount > 20
      disp('Eye movement detected!!!!');
      task.thistrial.dead = 1;
      stimulus.live.eyeDead=1;
      return
    elseif dist > 1.5
      stimulus.live.eyeCount = stimulus.live.eyeCount + 1;
    end
  end
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
    disp('Starting trial--eye centered.');
    task = jumpSegment(task);
  end
end




mglClearScreen(0.5);
if ~any(task.thistrial.thisseg == [stimulus.seg.feedback,stimulus.seg.ITI])
    upCue(task);
end

% Draw the stimuli at the correct flicker rate.
if task.thistrial.thisseg == stimulus.seg.stim1
  	mglBltTexture(stimulus.live.leftStims{1}, [-stimulus.stimXPos, -stimulus.stimYPos, stimulus.imSize, stimulus.imSize]);
    mglBltTexture(stimulus.live.rightStims{1}, [stimulus.stimXPos, -stimulus.stimYPos, stimulus.imSize, stimulus.imSize]);
  	mglBltTexture(stimulus.live.leftStims{2}, [-stimulus.stimXPos, stimulus.stimYPos, stimulus.imSize, stimulus.imSize]);
    mglBltTexture(stimulus.live.rightStims{2}, [stimulus.stimXPos, stimulus.stimYPos, stimulus.imSize, stimulus.imSize]);
end

if task.thistrial.thisseg == stimulus.seg.feedback
	if task.thistrial.correct == 1
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
    % Button 1 = Same; Button 2 = different
    if task.thistrial.cueSide == 1
    	task.thistrial.correct = mod(task.thistrial.response,2) == task.thistrial.rightSame;
    else
    	task.thistrial.correct = mod(task.thistrial.response,2) == (task.thistrial.leftSame);
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
files = dir(fullfile(sprintf('~/data/tex_att/%s/18*stim*.mat',mglGetSID)));
trials = 0;

for fi = 1:length(files)
    load(fullfile(sprintf('~/data/tex_att/%s/%s',mglGetSID,files(fi).name)));
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

