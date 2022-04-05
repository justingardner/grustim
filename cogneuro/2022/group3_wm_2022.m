function [ myscreen ] = group3_wm_2022( varargin )
%
% Attention task
%  Contrast increment detection task.
%
%  Usage: group3_wm_2022(varargin)
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
task{1}.segmin = [0.1, 4.0, 10.0, 2.0, 0.2, 2.0];
task{1}.segmax = [0.1, 4.0, 10.0, 2.0, 0.2, 2.0];
stimulus.seg = {};
stimulus.seg.fix = 1;
stimulus.seg.stim = 2; % stimuli are presented
stimulus.seg.delay = 3;
stimulus.seg.response = 4;
stimulus.seg.feedback = 5; % stimuli turn off, did you get a reward or not
stimulus.seg.ITI=6;

%% Task important variables
stimulus.nImages = 20;
stimulus.arraySize = 8;

% Set display parameters (image size, image eccentricity).
stimulus.imSize = 6;
stimulus.eccentricity = 10;

% Trial parameters
task{1}.synchToVol = zeros(size(task{1}.segmin));
if stimulus.scan
  task{1}.synchToVol(end)=1;
end
task{1}.getResponse = zeros(size(task{1}.segmin));
task{1}.getResponse(stimulus.seg.response)=1;

% Make numTrials some multiple of number of TrialTypes 
task{1}.numTrials = 180; 
task{1}.random = 1;

% Task parameters
task{1}.parameter.numCued = [2,4,6];
task{1}.parameter.rotation = [0, 50]; % How much is the probe stim rotated?
task{1}.parameter.match  = [-1, 0, 1]; % 1: yes cued, 0: noncued but in stimulus array, -1: not presented in stim array

% Task variables to be calculated later
task{1}.randVars.calculated.detected = NaN;
task{1}.randVars.calculated.response=NaN;
task{1}.randVars.calculated.visible = 1;
task{1}.randVars.calculated.correct = NaN;
task{1}.randVars.calculated.stim_inds = {NaN};
task{1}.randVars.calculated.cued_stims = {NaN};
task{1}.randVars.calculated.probe_stim = NaN;

% Preload images
stim_dir = '~/proj/grustim/cogneuro/2022/task_images';
stimulus.nImages = 20;
for i = 1:stimulus.nImages
  stimulus.stims.(sprintf('image%i', i)) = genTexFromIm(imread(sprintf('%s/task_symbol%02i.png', stim_dir, i)));
end

%%%%%%%%%%%%%%%% MGL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialize task and begin update loop
[task{1}, myscreen] = initTask(task{1},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);

% Run the eye calibration
myscreen = eyeCalibDisp(myscreen);

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
task.thistrial.correct = 0;

stimulus.live.gotResponse = 0;
stimulus.curTrial(task.thistrial.thisphase) = stimulus.curTrial(task.thistrial.thisphase) + 1;

task.thistrial.stim_inds = randsample(stimulus.nImages, stimulus.arraySize);
task.thistrial.cued_stims = randsample(task.thistrial.stim_inds, task.thistrial.numCued);
if task.thistrial.match == 1 
  task.thistrial.probe_stim = randsample(task.thistrial.cued_stims, 1);
elseif task.thistrial.match == 0
  task.thistrial.probe_stim = randsample(setdiff(task.thistrial.stim_inds, task.thistrial.cued_stims), 1);
else
  task.thistrial.probe_stim = randsample(setdiff(1:20, task.thistrial.stim_inds), 1);
end

%disp trial parameters each trial
disp(sprintf('(group3_wm_2022) Trial %d - Cued %i stim; Rotation: %i', task.trialnum, task.thistrial.numCued, task.thistrial.rotation));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)

global stimulus

% Select image parameters: size, eccentricity, and location
ecc = stimulus.eccentricity;
imSz = stimulus.imSize;
cueDirs = [pi, 0];
colors = [stimulus.colors.red; stimulus.colors.green];
cue_colors = [stimulus.colors.black; stimulus.colors.blue];

% After the response segment (at start of feedback segment), decide whether they get a reward or not.
for i = 1:2
  mglClearScreen(0.5);
  upFix(stimulus, stimulus.colors.white);
  if task.thistrial.thisseg == stimulus.seg.stim
    for j = 1:stimulus.arraySize
      theta = (j-1)*(2*pi/stimulus.arraySize);
      mglQuad([ecc*cos(theta)-imSz/2; ecc*cos(theta)-imSz/2; ecc*cos(theta)+imSz/2; ecc*cos(theta)+imSz/2],...
              [ecc*sin(theta)-imSz/2; ecc*sin(theta)+imSz/2; ecc*sin(theta)+imSz/2; ecc*sin(theta)-imSz/2],...
              cue_colors(any(task.thistrial.cued_stims == task.thistrial.stim_inds(j))+1,:)');
      mglBltTexture(stimulus.stims.(sprintf('image%i', task.thistrial.stim_inds(j))), [ecc*cos(theta), ecc*sin(theta), imSz*.8, imSz*.8]);
    end 
  elseif task.thistrial.thisseg == stimulus.seg.response
    mglBltTexture(stimulus.stims.(sprintf('image%i', task.thistrial.probe_stim)), [0,0,imSz,imSz], 0, 0, task.thistrial.rotation);
  elseif task.thistrial.thisseg == stimulus.seg.feedback
    upFix(stimulus, colors(task.thistrial.correct+1,:));
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
    if task.thistrial.match == 1
      task.thistrial.correct = (task.thistrial.response == 1);
    else
      task.thistrial.correct = (task.thistrial.response == 2);
    end
    if task.thistrial.correct
      disp('Correct response');
    else
      disp('Incorrect response');
    end
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
% Draws an arrow cuing towards direction specified by theta
function drawCue(theta)
x=[1,.75,2,.75];
y=[0,.5,0,-.5];
mglPolygon(cos(theta)*x - sin(theta)*y, sin(theta)*x + cos(theta)*y, [0,0,0])


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

