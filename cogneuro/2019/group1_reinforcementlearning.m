function [ myscreen ] = rewardLearning( varargin )
%
% Reward Learning task
%  Reinforcement learning task
%
%  Usage: rewardLearning(varargin)
%  Authors: Akshay Jagadeesh
%  Created: 04/11/2019
%

global stimulus

stimulus = struct;

%% Initialize Variables

% add arguments later
plots = 0;
noeye = 0;
getArgs(varargin,{'plots=0','noeye=1'});
stimulus.plots = plots;
stimulus.noeye = noeye;
clear noeye plots

%% Stimulus parameters 
%% Open Old Stimfile
stimulus.counter = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/rewardLearning/%s',mglGetSID))
  % Directory exists, check for a stimfile
  files = dir(sprintf('~/data/rewardLearning/%s/1*mat',mglGetSID));
 
  if length(files) >= 1
    fname = files(end).name;
     
    s = load(sprintf('~/data/rewardLearning/%s/%s',mglGetSID,fname));
    stimulus.counter = s.stimulus.counter + 1;
    clear s;
    disp(sprintf('(rewardLearning) Data file: %s loaded.',fname));
  end
end
disp(sprintf('(rewardLearning) This is run #%i',stimulus.counter));

%% Setup Screen
myscreen = initScreen('VPixx2');

% set background to grey
myscreen.background = 0.5;

%% Initialize Stimulus
myscreen = initStimulus('stimulus',myscreen);
localInitStimulus();
  
% Set response keys
stimulus.responseKeys = [1 2];
%% [aj]: change to add response keys for whichever 3 keys we want.

% set colors
stimulus.colors.white = [1 1 1];
stimulus.colors.black = [0 0 0];
stimulus.colors.red = [1 0 0];
stimulus.colors.green = [0 1 0];
stimulus.colors.blue = [0 0 1];
stimulus.live.fixColor = stimulus.colors.blue;
stimulus.live.cueColor = stimulus.colors.black;

%%%%%%%%%%%%% SET TASK VARIABLES  %%%%%%%%%%%%%%%%%
stimulus.curTrial(1) = 0;
task{1} = struct;
task{1}.waitForBacktick = 1;

% Define stimulus timing
task{1}.segmin = [inf, 2.0, 2.0, 1.0, 2.0, 2.0];
task{1}.segmax = [inf, 2.0, 2.0, 1.0, 2.0, 4.0];
stimlus.seg = {};
stimulus.seg.fix = 1;
stimulus.seg.stim = 2; % stimuli are presented
stimulus.seg.response = 3; % fixation cross goes off and they have to respond (stimuli still present) 
stimulus.seg.wait = 4; % ISI
stimulus.seg.feedback = 5; % stimuli turn off, did you get a reward or not
stimulus.seg.ITI=6;
% Set fixation length to constant if not eye tracking.
if stimulus.noeye==1
  task{1}.segmin(1) = 0.1;
  task{1}.segmax(1) = 0.1;
end

%% Task important variables
% Set the images which should be used as unconditioned stimuli
stimulus.all_imNames = {'fronds', 'monarchs', 'paisley', 'rocks', 'zebras', 'bananas'};
counter = mod(stimulus.counter, 3)+1;
stimulus.imNames = {stimulus.all_imNames{2*counter - 1}, stimulus.all_imNames{2*counter}};
stimulus.stimDir = '~/proj/rewardLearning/stimuli';

% Set the range of values. Stimuli will be randomly assigned probabilities
%     at evenly spaced intervals between rewardRange(0) and rewardRange(1).
stimulus.rewardRange = [.2, .8];
stimulus.rewardProbs = stimulus.rewardRange(1):((stimulus.rewardRange(2)-stimulus.rewardRange(1))/(length(stimulus.imNames)-1)):stimulus.rewardRange(2);

% Set display parameters (image size, image eccentricity).
stimulus.imSize = 6;
stimulus.eccentricity = 10;
stimulus.live.mask = imread('~/proj/TextureSynthesis/stimuli/Flattop8.tif');

% Trial parameters
%task{1}.parameter.leftIm = 1:length(stimulus.imNames);

task{1}.synchToVol = zeros(size(task{1}.segmin));
task{1}.getResponse = zeros(size(task{1}.segmin));
task{1}.getResponse(stimulus.seg.response)=1;

% Make numTrials some multiple of number of TrialTypes 
task{1}.numTrials = 180; 
task{1}.random = 1;

% Task variables to be calculated later
task{1}.randVars.calculated.leftIm = 1:length(stimulus.imNames);
task{1}.randVars.calculated.rightIm=NaN;
task{1}.randVars.calculated.leftProb=NaN;
task{1}.randVars.calculated.rightProb=NaN;
task{1}.randVars.calculated.chosenSideRewardProb=NaN; % reward probability on the chosen side.

task{1}.randVars.calculated.chosenSide=NaN; % Which side (0: left or 1: right) did the subject choose.
task{1}.randVars.calculated.rewarded=NaN; % did they get a reward on this trial or not

task{1}.randVars.calculated.dead = 0;
task{1}.randVars.calculated.visible = 1;

%% Preload images
stims = struct();
disppercent(-inf, 'Preloading images');
for i = 1:length(stimulus.imNames)
  imName = stimulus.imNames{i};
  stims.(imName) = genTexFromIm(imread(sprintf('%s/%s.jpg', stimulus.stimDir, imName)), stimulus.live.mask);
  disp(sprintf('Assigning image %s to probability %g', imName, stimulus.rewardProbs(i)));
  disppercent(i / length(stimulus.imNames));
end
stims.win = genTexFromIm(imread(sprintf('%s/win.png', stimulus.stimDir)), stimulus.live.mask);
stims.lose = genTexFromIm(imread(sprintf('%s/lose.png', stimulus.stimDir)), stimulus.live.mask);
disppercent(inf);
stimulus.live.stims = stims;

%% Save
clear stims

%%%%%%%%%%%%%%%% MGL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialize task and begin update loop
[task{1}, myscreen] = initTask(task{1},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);

% Run the eye calibration
myscreen = eyeCalibDisp(myscreen);

% let the user know
disp(sprintf('(rewardLearning) Starting run number: %i.',stimulus.counter));

%% Main Task Loop
mglClearScreen(0.5); 
upFix(stimulus);

mglFlush
mglClearScreen(0.5); 
upFix(stimulus);

phaseNum = 1;
while (phaseNum <= length(task)) && ~myscreen.userHitEsc
  % update the task
  [task, myscreen, phaseNum] = updateTask(task,myscreen,phaseNum);
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
task.thistrial.response = NaN;

stimulus.live.gotResponse = 0;
stimulus.curTrial(task.thistrial.thisphase) = stimulus.curTrial(task.thistrial.thisphase) + 1;

% directories
texDir = stimulus.stimDir;

%% Load both images for this trial
trial_ims = randsample(1:length(stimulus.imNames), 2, false);
task.thistrial.leftIm = trial_ims(1);
task.thistrial.rightIm = trial_ims(2);
task.thistrial.leftProb = stimulus.rewardProbs(task.thistrial.leftIm);
task.thistrial.rightProb = stimulus.rewardProbs(task.thistrial.rightIm);

leftImName = stimulus.imNames{task.thistrial.leftIm};
rightImName = stimulus.imNames{task.thistrial.rightIm};
stimulus.live.leftStim = stimulus.live.stims.(leftImName);
stimulus.live.rightStim = stimulus.live.stims.(rightImName);
stimulus.live.trialProbs = [task.thistrial.leftProb, task.thistrial.rightProb];

%dispp trial parameters each trial
disp(sprintf('(rewardLearning) Trial %d - Left Side: %s (P=%g), Right Side: %s (P=%g)', task.trialnum, leftImName, task.thistrial.leftProb, rightImName, task.thistrial.rightProb));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)

global stimulus

stimulus.live.triggerWaiting = 0;
if any(task.thistrial.thisseg==[stimulus.seg.fix])
  stimulus.live.triggerWaiting = 1;
  stimulus.live.centered = 0;
  stimulus.live.triggerTime = 0;
  stimulus.live.lastTrigger = -1;
end

stimulus.live.eyeDead = 0;

% Select image parameters: size, eccentricity, and location
imSz = stimulus.imSize; % Size in Degrees of Visual Angle to display image
ecc = stimulus.eccentricity; % Eccentricity to display image at
stimulus.locations = [0, ecc; ecc*cosd(30), -ecc*sind(30); -ecc*cosd(30), -ecc*sind(30)];

if task.thistrial.thisseg == stimulus.seg.wait
  task.thistrial.chosenSide = task.thistrial.response; % Save which side they responded for
  if isnan(task.thistrial.chosenSide)
    task.thistrial.chosenSideRewardProb = -1;
  else
    task.thistrial.chosenSideRewardProb = stimulus.live.trialProbs(task.thistrial.chosenSide);
  end
  % Then flip a weighted coin to determine if they get a reward or not on this trial.
  task.thistrial.rewarded = (rand() < task.thistrial.chosenSideRewardProb); 
  
  if task.thistrial.rewarded
    disp(sprintf('(rewardLearning) Nice! Subject won $1 :)'));
  else
    disp(sprintf('(rewardLearning) Aww, no reward :('));
  end
end

for i = 1:2
  mglClearScreen(0.5);
  if any(task.thistrial.thisseg == [stimulus.seg.stim, stimulus.seg.response])
    mglBltTexture(stimulus.live.leftStim, [-ecc, 0, imSz, imSz]);
    mglBltTexture(stimulus.live.rightStim,[ecc, 0, imSz, imSz]);
    
    if task.thistrial.thisseg == stimulus.seg.stim
      upFix(stimulus, stimulus.colors.white); % turn off fixation to indicate response.
    end
  elseif any(task.thistrial.thisseg == [stimulus.seg.wait, stimulus.seg.feedback])
    % Draw the stimulus on the chosen side.
    if task.thistrial.chosenSide == 1
      mglBltTexture(stimulus.live.leftStim, [-ecc, 0, imSz, imSz]);
    elseif task.thistrial.chosenSide == 2
      mglBltTexture(stimulus.live.rightStim, [ecc, 0, imSz, imSz]);
    end
    
    if task.thistrial.thisseg == stimulus.seg.feedback
      % Show the reward
      if task.thistrial.rewarded
        mglBltTexture(mglText('Nice! You win $1 :)'), [0,3]);
        mglBltTexture(stimulus.live.stims.win, [0,0,4,4]);
        %upFix(stimulus, stimulus.colors.green);
      else
        mglBltTexture(mglText('Sorry! No reward :('), [0,3]);
        mglBltTexture(stimulus.live.stims.lose, [0,0,4,4]);
        %upFix(stimulus, stimulus.colors.red);
      end
    end
  else
    upFix(stimulus, stimulus.colors.black);
  end
  mglFlush
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
    if dist > 1.5 && stimulus.live.eyeCount > 30
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus

if task.thistrial.dead, return; end

% only allow a response if it is one of the pre-specified responses (in this case, y g or h)
validResponse = any(task.thistrial.whichButton == stimulus.responseKeys);
if validResponse
  if stimulus.live.gotResponse==0
    task.thistrial.detected = 1;
    task.thistrial.response = find(stimulus.responseKeys == task.thistrial.whichButton);
    stimulus.live.fix = 0;
  else
    disp(sprintf('Subject responded multiple times: %i',stimulus.live.gotResponse));
  end
  %disp('jumping segment');
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

