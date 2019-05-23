function [ myscreen ] = group1_reinforcementlearning( varargin )
%
% Reward Learning task
%  Reinforcement learning task
%
%  Usage: group1_reinforcementlearning(varargin)
%  Authors: Akshay Jagadeesh
%  Created: 04/11/2019
%

global stimulus

stimulus = struct;

%% Initialize Variables

% add arguments later
getArgs(varargin,{'plots=0', 'scan=0'});
stimulus.plots = plots;
stimulus.scan = scan;
clear plots

%% Stimulus parameters 
%% Open Old Stimfile
stimulus.counter = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/group1_reinforcementlearning/%s',mglGetSID))
  % Directory exists, check for a stimfile
  files = dir(sprintf('~/data/group1_reinforcementlearning/%s/1*mat',mglGetSID));
 
  if length(files) >= 1
    fname = files(end).name;
     
    s = load(sprintf('~/data/group1_reinforcementlearning/%s/%s',mglGetSID,fname));
    stimulus.counter = s.stimulus.counter + 1;
    clear s;
    disp(sprintf('(group1_reinforcementlearning) Data file: %s loaded.',fname));
  end
end
disp(sprintf('(group1_reinforcementlearning) This is run #%i',stimulus.counter));

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
stimulus.responseKeys = [1 2];

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
task{1}.segmin = [0.1, 2.0, 2.0, 2.0];
task{1}.segmax = [0.1, 2.0, 2.0, 5.0];
stimlus.seg = {};
stimulus.seg.fix = 1;
stimulus.seg.stim = 2; % stimuli are presented
%stimulus.seg.response = 3; % fixation cross goes off and they have to respond (stimuli still present) 
stimulus.seg.feedback = 3; % stimuli turn off, did you get a reward or not
stimulus.seg.ITI=4;

%% Task important variables
% Set the images which should be used as unconditioned stimuli
stimulus.all_imNames = {'fronds', 'monarchs', 'paisley', 'rocks', 'zebras', 'bananas'};
counter = mod(stimulus.counter, 3)+1;
stimulus.imNames = {stimulus.all_imNames{2*counter - 1}, stimulus.all_imNames{2*counter}};
stimulus.stimDir = '~/proj/rewardlearning/stimuli';

% Set the range of values. Stimuli will be randomly assigned probabilities
%     at evenly spaced intervals between rewardRange(0) and rewardRange(1).
stimulus.rewardRange = [.2, .8];
stimulus.rewardProbs = stimulus.rewardRange(1):((stimulus.rewardRange(2)-stimulus.rewardRange(1))/(length(stimulus.imNames)-1)):stimulus.rewardRange(2);

stimulus.initialRewardProbs = [.4 .6];

% Set display parameters (image size, image eccentricity).
stimulus.imSize = 6;
stimulus.eccentricity = 10;
stimulus.live.mask = imread('~/proj/TextureSynthesis/stimuli/Flattop8.tif');

% Trial parameters
task{1}.synchToVol = zeros(size(task{1}.segmin));
if stimulus.scan
  task{1}.synchToVol(end)=1;
end
task{1}.getResponse = zeros(size(task{1}.segmin));
task{1}.getResponse(stimulus.seg.stim)=1;

% Make numTrials some multiple of number of TrialTypes 
task{1}.numTrials = 180; 
task{1}.random = 1;

% Task variables to be calculated later
task{1}.randVars.calculated.leftIm = NaN; %1:length(stimulus.imNames);
task{1}.randVars.calculated.rightIm=NaN;
task{1}.randVars.calculated.imgRewardProbs = {NaN};
task{1}.randVars.calculated.leftProb=NaN;
task{1}.randVars.calculated.rightProb=NaN;
task{1}.randVars.calculated.chosenSideRewardProb=NaN; % reward probability on the chosen side.

task{1}.randVars.calculated.chosenSide=NaN; % Which side (0: left or 1: right) did the subject choose.
task{1}.randVars.calculated.rewarded=NaN; % did they get a reward on this trial or not

task{1}.randVars.calculated.visible = 1;

%% Preload images
stims = struct();
disppercent(-inf, 'Preloading images');
for i = 1:length(stimulus.imNames)
  imName = stimulus.imNames{i};
  stims.(imName) = genTexFromIm(imread(sprintf('%s/%s.jpg', stimulus.stimDir, imName)), stimulus.live.mask);
  disp(sprintf('Initializing image %s to probability %g', imName, stimulus.rewardProbs(i)));
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
disp(sprintf('(group1_reinforcementlearning) Starting run number: %i.',stimulus.counter));

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

task.thistrial.detected = 0;
task.thistrial.visible = 1;
task.thistrial.response = NaN;

stimulus.live.gotResponse = 0;
stimulus.curTrial(task.thistrial.thisphase) = stimulus.curTrial(task.thistrial.thisphase) + 1;

% directories
texDir = stimulus.stimDir;

%% Randomly select which image will be on the left and which on the right.
trial_ims = randsample(1:length(stimulus.imNames), 2, false);
task.thistrial.leftIm = trial_ims(1);
task.thistrial.rightIm = trial_ims(2);

% Set the reward probabilities for this trial
if task.trialnum == 1
  % Initialize reward probabilities on trial 1.
  task.thistrial.imgRewardProbs = stimulus.initialRewardProbs;
else
  % Change one of the images reward probability by a random amount.
  task.thistrial.rewardDelta = .03*randn();
  task.thistrial.imgRewardProbs(1) = task.lasttrial.imgRewardProbs(1) + task.thistrial.rewardDelta;
  % If the change would take you out of the range, then bounce back the other direction.
  if task.thistrial.imgRewardProbs(1) < stimulus.rewardRange(1) || task.thistrial.imgRewardProbs(1) > stimulus.rewardRange(2)
    task.thistrial.imgRewardProbs(1) = task.lasttrial.imgRewardProbs(1) - task.thistrial.rewardDelta;
    task.thistrial.rewardDelta = -1*task.thistrial.rewardDelta;
  end

  % Make the other images reward probability a mirror (1-n) of the first.
  task.thistrial.imgRewardProbs(2) = 1-task.thistrial.imgRewardProbs(1);
end
task.thistrial.leftProb = task.thistrial.imgRewardProbs(task.thistrial.leftIm);
task.thistrial.rightProb = task.thistrial.imgRewardProbs(task.thistrial.rightIm);

leftImName = stimulus.imNames{task.thistrial.leftIm};
rightImName = stimulus.imNames{task.thistrial.rightIm};
stimulus.live.leftStim = stimulus.live.stims.(leftImName);
stimulus.live.rightStim = stimulus.live.stims.(rightImName);
stimulus.live.trialProbs = [task.thistrial.leftProb, task.thistrial.rightProb];

%dispp trial parameters each trial
disp(sprintf('(group1_reinforcementlearning) Trial %d - Left Side: %s (P=%g), Right Side: %s (P=%g)', task.trialnum, leftImName, task.thistrial.leftProb, rightImName, task.thistrial.rightProb));

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

% After the response segment (at start of feedback segment), decide whether they get a reward or not.
if task.thistrial.thisseg == stimulus.seg.feedback
  task.thistrial.chosenSide = task.thistrial.response; % Save which side they responded for
  if isnan(task.thistrial.chosenSide) || task.thistrial.chosenSide==-1
    task.thistrial.chosenSideRewardProb = -1;
  else
    task.thistrial.chosenSideRewardProb = stimulus.live.trialProbs(task.thistrial.chosenSide);
  end
  % Then flip a weighted coin to determine if they get a reward or not on this trial.
  task.thistrial.rewarded = (rand() < task.thistrial.chosenSideRewardProb); 
  
  if task.thistrial.rewarded
    disp(sprintf('(group1_reinforcementlearning) Nice! Subject won $1 :)'));
  else
    disp(sprintf('(group1_reinforcementlearning) Aww, no reward :('));
  end
end

for i = 1:2
  mglClearScreen(0.5);
  if any(task.thistrial.thisseg == [stimulus.seg.stim])
    mglBltTexture(stimulus.live.leftStim, [-ecc, 0, imSz, imSz]);
    mglBltTexture(stimulus.live.rightStim,[ecc, 0, imSz, imSz]);
    
    if task.thistrial.thisseg == stimulus.seg.stim
      upFix(stimulus, stimulus.colors.white); % turn off fixation to indicate response.
    end
  elseif any(task.thistrial.thisseg == [stimulus.seg.feedback])
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus


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
  task.thistrial.response = NaN;
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

