function [ myscreen ] = texatt_2isd( varargin )
% TEXATT_2ISD: How does attention affect visual cortical representations of natural and synthesized images?
%
%  2-interval same-different task with focal vs distributed attentional cue.
%
%  Usage: texatt_2isd(varargin)
%  Authors: Akshay Jagadeesh
%  Date: 02/03/2022
%

global stimulus

stimulus = struct;

%% Initialize Variables

% add arguments later
scan = 1;
getArgs(varargin,{'scan=0','noeye=0', 'nTrials=90', 'testing=0'}, 'verbose=1');
stimulus.scan = scan;
stimulus.noeye = noeye;
stimulus.nTrials = nTrials;

%% Stimulus parameters 
%% Open Old Stimfile
stimulus.counter = 1;

%% Setup Screen
if(stimulus.scan)
  dispName = 'fMRIprojFlex';
elseif(testing==1)
  dispName = 'test';
else
  dispName = 'VPixx2';
end
myscreen = initScreen(dispName);

% set background to grey
myscreen.background = 0.5;

%% Initialize Stimulus
myscreen = initStimulus('stimulus',myscreen);
  
% Set response keys
stimulus.responseKeys = [1 2 0]; 
%stimulus.responseKeys = [2 3];
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
stimulus.imageNames = {'elephant', 'tiger', 'car', 'apple', 'bear', 'cat', 'jetplane', 'helicopter', 'face', 'grass', 'sand', 'rocks', 'dirt', 'lawn', 'bark', 'worms', 'bricks'};
stimulus.layerNames = {'pool4'};
stimulus.poolSize = '1x1_';

stimulus.nTexFams = length(stimulus.imageNames);
stimulus.imSize = 8;
stimulus.stimXPos = 8;
stimulus.stimYPos = 6;
stimulus.num_samples = 2;

%% Select the condition for this run
% Choose which image and which pooling layer to display on this run on each side
stimulus.stimDir = '~/proj/texture_stimuli/color/textures';
stimulus.origDir = '~/proj/texture_stimuli/color/originals';

%% Preload images
mask = imread('~/proj/TextureSynthesis/stimuli/Flattop8.tif');
stimulus.images = struct();
disppercent(-inf, 'Preloading images');

% load texture and noise samples
for i = 1:stimulus.nTexFams
  imName = stimulus.imageNames{i};

  for j = 1:stimulus.num_samples
    if isfile(    sprintf('%s/1x1_pool4_%s_smp%i.png', stimulus.stimDir, imName, j))
      sd = imread(sprintf('%s/1x1_pool4_%s_smp%i.png', stimulus.stimDir, imName, j));
    elseif isfile(sprintf('%s/1x1_pool4_%s_smp%i.jpg', stimulus.stimDir, imName, j))
      sd = imread(sprintf('%s/1x1_pool4_%s_smp%i.jpg', stimulus.stimDir, imName, j));
    else
      disp(sprintf('1x1_pool4_%s_smp%i not found', imName, j));
      keyboard
    end
    stimulus.images.(sprintf('%s_s%i', imName, j)) = genTexFromIm(sd, mask);
  end

  % Load noise samples.
  orig = imread(sprintf('%s/%s.jpg', stimulus.origDir, imName));
  stimulus.images.(sprintf('%s_s0', imName)) = genTexFromIm(orig, mask);
  
  disppercent(i / stimulus.nTexFams);
end
disppercent(inf);
clear sd

%%%%%%%%%%%%% TASK %%%%%%%%%%%%%%%%%
task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;

                   % fix  cue stim1 isi stim2 resp iti
task{1}{1}.segmin = [inf, 0.4, 0.6, 0.2, 0.6, 1.6, 0.4];
task{1}{1}.segmax = [inf, 0.4, 0.6, 0.2, 0.6, 1.6, 0.4];

stimulus.seg = {};
stimulus.seg.fix = 1;
stimulus.seg.cue = 2;
stimulus.seg.stim1 = 3;
stimulus.seg.isi = 4;
stimulus.seg.stim2 = 5;
stimulus.seg.resp = 6;
stimulus.seg.iti = 7;

% Trial parameters
task{1}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
task{1}{1}.getResponse = zeros(size(task{1}{1}.segmin));
task{1}{1}.getResponse(stimulus.seg.resp) = 1;

task{1}{1}.numTrials = stimulus.nTrials;
task{1}{1}.random = 1;

if stimulus.noeye==1
  task{1}{1}.segmin(1) = 0.1;
  task{1}{1}.segmax(1) = 0.1;
end

if stimulus.scan
  task{1}{1}.synchToVol(end) = 1;
  % Shorten the last segment to account for synchtovol
  task{1}{1}.segmin(end) = max(0, task{1}{1}.segmin(end) - 0.100);
  task{1}{1}.segmax(end) = max(0, task{1}{1}.segmax(end) - 0.100);
end

% Initialize task parameters

% Assign layer, texture family, and sample index in random blocks.
task{1}{1}.parameter.leftTopImgClass1 = 1:length(stimulus.imageNames);
task{1}{1}.parameter.rightTopImgClass1 = 1:length(stimulus.imageNames);
task{1}{1}.parameter.leftBottomImgClass1 = 1:length(stimulus.imageNames);
task{1}{1}.parameter.rightBottomImgClass1 = 1:length(stimulus.imageNames);

task{1}{1}.randVars.calculated.leftTopImgClass2 = NaN;
task{1}{1}.randVars.calculated.rightTopImgClass2 = NaN;
task{1}{1}.randVars.calculated.leftBottomImgClass2 = NaN;
task{1}{1}.randVars.calculated.rightBottomImgClass2 = NaN;

task{1}{1}.parameter.cueFocal = [0,1,1]; % 0 = distributed, 1 = focal
task{1}{1}.parameter.cueSide = [1 2 3 4]; % 1 = left top, 2 = right top, 3 = left bottom, 4 = right bottom
task{1}{1}.randVars.calculated.cueType = NaN; % 1 = left, 2 = right, 0 = distributed
task{1}{1}.randVars.uniform.leftTopSame = [0,1]; % 0 = different, 1 = same
task{1}{1}.randVars.uniform.leftBottomSame = [0,1]; % 0 = different, 1 = same
task{1}{1}.randVars.uniform.rightTopSame = [0,1];
task{1}{1}.randVars.uniform.rightBottomSame = [0,1];

task{1}{1}.randVars.uniform.leftTopSample1 = [0]; % 0 = original, 1 = synth sample 1, 2 = synth sample2
task{1}{1}.randVars.uniform.leftBottomSample1 = [0]; % 0 = original, 1 = synth sample 1, 2 = synth sample2
task{1}{1}.randVars.uniform.rightTopSample1 = [0]; 
task{1}{1}.randVars.uniform.rightBottomSample1 = [0]; 
task{1}{1}.randVars.uniform.leftTopSample2 = [0];
task{1}{1}.randVars.uniform.leftBottomSample2 = [0];
task{1}{1}.randVars.uniform.rightTopSample2 = [0];
task{1}{1}.randVars.uniform.rightBottomSample2 = [0];

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
% Select the stimulus for this trial
% Left side
locations = {'leftTop', 'rightTop', 'leftBottom', 'rightBottom'};
for li = 1:length(locations)
	loc = locations{li};
	imgClass1 = stimulus.imageNames{task.thistrial.(sprintf('%sImgClass1', loc))};
	if task.thistrial.(sprintf('%sSame', loc)) == 0
	  other_imgclasses = setdiff(1:length(stimulus.imageNames), task.thistrial.(sprintf('%sImgClass1', loc)));
	  task.thistrial.(sprintf('%sImgClass2', loc)) = other_imgclasses(randsample(length(other_imgclasses),1));
	else
    task.thistrial.(sprintf('%sImgClass2', loc)) = task.thistrial.(sprintf('%sImgClass1', loc));
	end
  imgClass2 = stimulus.imageNames{task.thistrial.(sprintf('%sImgClass2', loc))};

	stimulus.live.(sprintf('%sStim1', loc)) = stimulus.images.(sprintf('%s_s%i', imgClass1, task.thistrial.(sprintf('%sSample1', loc))));
	stimulus.live.(sprintf('%sStim2', loc)) = stimulus.images.(sprintf('%s_s%i', imgClass2, task.thistrial.(sprintf('%sSample2', loc))));
end

stimulus.live.eyeCount = 0;
cueType = {'Distributed', 'Focal'};
cueSide = locations;
sameDiff = {'Different', 'Same'};
disp(sprintf('--- Trial %i - %s Cue, %s Response; Correct Response: %s', task.trialnum, cueType{1+task.thistrial.cueFocal}, cueSide{task.thistrial.cueSide}, sameDiff{1+task.thistrial.(sprintf('%sSame', locations{task.thistrial.cueSide}))}));

upFix(stimulus);

if task.thistrial.cueFocal==1
  task.thistrial.cueType = task.thistrial.cueSide;
else
  task.thistrial.cueType = 0;
end

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

%%% CODE TO DRAW IMAGES
% Draw pre-cue and response cues 
if ~any(task.thistrial.thisseg == [stimulus.seg.fix, stimulus.seg.iti])
  upCue(task);
end
locations = {'leftTop', 'rightTop', 'leftBottom', 'rightBottom'};
x = stimulus.stimXPos; 
y = stimulus.stimYPos;
coordinates = [-x y; x y; -x -y; x -y];

% Draw the stimuli.
if task.thistrial.thisseg == stimulus.seg.stim1
	for li = 1:length(locations)
		loc = locations{li};
	  	mglBltTexture(stimulus.live.(sprintf('%sStim1', loc)), [coordinates(li,1), coordinates(li,2), stimulus.imSize, stimulus.imSize]);
    end
elseif task.thistrial.thisseg == stimulus.seg.stim2
	for li = 1:length(locations)
		loc = locations{li};
	  	mglBltTexture(stimulus.live.(sprintf('%sStim2', loc)), [coordinates(li,1), coordinates(li,2), stimulus.imSize, stimulus.imSize]);
    end
end

if task.thistrial.thisseg == stimulus.seg.iti && stimulus.live.gotResponse == 1
  % If subject has already responded, change cross color to green/red
	if task.thistrial.correct == 1
		upFix(stimulus, stimulus.colors.green);
  else
		upFix(stimulus, stimulus.colors.red);
	end
else
  % Draw neutral color cross.
	upFix(stimulus);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Draws attention cue %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function upCue(task)

global stimulus
x = stimulus.stimXPos; 
y = stimulus.stimYPos;
angles = [pi-atan(y/x); atan(y/x); pi+atan(y/x); atan(-y/x)];

if ~(task.thistrial.thisseg == stimulus.seg.resp) && task.thistrial.cueFocal == 0 % Draw distributed cue
	drawArrow(angles(1));
	drawArrow(angles(2));
	drawArrow(angles(3));
	drawArrow(angles(4));
elseif task.thistrial.cueSide == 1 % Cue left top
	drawArrow(angles(1));
elseif task.thistrial.cueSide == 2 % Cue right top 
	drawArrow(angles(2));
elseif task.thistrial.cueSide == 3 % Cue left bottom
	drawArrow(angles(3));
elseif task.thistrial.cueSide == 4 % Cue right bottom
	drawArrow(angles(4));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Draws arrow %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawArrow(theta)
x=[1,.75,2,.75]*.6;
y=[0,.5,0,-.5]*.6;
mglPolygon(cos(theta)*x - sin(theta)*y, sin(theta)*x + cos(theta)*y, [0,0,0])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus


validResponse = any(str2num(myscreen.responseKeys{task.thistrial.whichButton}) == stimulus.responseKeys);
locations = {'leftTop', 'rightTop', 'leftBottom', 'rightBottom'};
if validResponse
  if stimulus.live.gotResponse==0
    task.thistrial.detected = 1;
    task.thistrial.response = str2num(myscreen.responseKeys{task.thistrial.whichButton});
    disp(sprintf('Response: %i', task.thistrial.response));
    % Button 1 = Same; Button 0 = different
    task.thistrial.correct = mod(task.thistrial.response,2) == task.thistrial.(sprintf('%sSame', locations{task.thistrial.cueSide}));
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
function images = localInitStimulus()

global stimulus

mask = imread('~/proj/TextureSynthesis/stimuli/Flattop8.tif');
images = struct();
disppercent(-inf, 'Preloading images');

% load texture and noise samples
for i = 1:stimulus.nTexFams
  imName = stimulus.imageNames{i};

  for j = 1:stimulus.num_samples
    if isfile(    sprintf('%s/1x1_pool4_%s_smp%i.png', stimulus.stimDir, imName, j))
      sd = imread(sprintf('%s/1x1_pool4_%s_smp%i.png', stimulus.stimDir, imName, j));
    elseif isfile(sprintf('%s/1x1_pool4_%s_smp%i.jpg', stimulus.stimDir, imName, j))
      sd = imread(sprintf('%s/1x1_pool4_%s_smp%i.jpg', stimulus.stimDir, imName, j));
    else
      disp(sprintf('1x1_pool4_%s_smp%i not found', imName, j));
      keyboard
    end
    images.(sprintf('%s_s%i', imName, j)) = genTexFromIm(sd, mask);
  end

  % Load noise samples.
  orig = imread(sprintf('%s/%s.jpg', stimulus.origDir, imName));
  images.(sprintf('%s_s0', imName)) = genTexFromIm(orig, mask);
  
  disppercent(i / stimulus.nTexFams);
end
disppercent(inf);
clear sd


