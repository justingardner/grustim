function [ myscreen ] = group3_attention( varargin )
%
% Attention task
%  Contrast increment detection task.
%
%  Usage: group3_attention(varargin)
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

if ~isempty(mglGetSID) && isdir(sprintf('~/data/group3_attention/%s',mglGetSID))
  % Directory exists, check for a stimfile
  files = dir(sprintf('~/data/group3_attention/%s/1*mat',mglGetSID));
 
  if length(files) >= 1
    fname = files(end).name;
     
    s = load(sprintf('~/data/group3_attention/%s/%s',mglGetSID,fname));
    stimulus.counter = s.stimulus.counter + 1;
    clear s;
    disp(sprintf('(group3_attention) Data file: %s loaded.',fname));
  end
end
disp(sprintf('(group3_attention) This is run #%i',stimulus.counter));

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
task{1}.segmin = [0.1, 1.0, 2.0, 1.0, 1.0, 3.9];
task{1}.segmax = [0.1, 1.0, 2.0, 1.0, 1.0, 3.9];
stimulus.seg = {};
stimulus.seg.fix = 1;
stimulus.seg.cue = 2;
stimulus.seg.stim = 3; % stimuli are presented
stimulus.seg.increment = 4;
stimulus.seg.feedback = 5; % stimuli turn off, did you get a reward or not
stimulus.seg.ITI=6;

%% Task important variables
% Set contrast parameters
stimulus.contrastDelta = .10;
stimulus.baseContrast = .30;

% Set display parameters (image size, image eccentricity).
stimulus.imSize = 6;
stimulus.eccentricity = 10;

% Trial parameters
task{1}.synchToVol = zeros(size(task{1}.segmin));
if stimulus.scan
  task{1}.synchToVol(end)=1;
end
task{1}.getResponse = zeros(size(task{1}.segmin));
task{1}.getResponse(stimulus.seg.increment)=1;

% Make numTrials some multiple of number of TrialTypes 
task{1}.numTrials = 180; 
task{1}.random = 1;

% Task parameters
task{1}.parameter.leftStimPresent = [0 1];
task{1}.parameter.rightStimPresent = [0 1];
task{1}.parameter.cueSide = [1 2];

% Task variables to be calculated later
task{1}.randVars.calculated.leftContrastDelta = NaN;
task{1}.randVars.calculated.rightContrastDelta = NaN;
task{1}.randVars.calculated.detected = NaN;
task{1}.randVars.calculated.response=NaN;
task{1}.randVars.calculated.visible = 1;


%% Save
clear stims

%%%%%%%%%%%%%%%% MGL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialize task and begin update loop
[task{1}, myscreen] = initTask(task{1},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);

% Run the eye calibration
myscreen = eyeCalibDisp(myscreen);

% let the user know
disp(sprintf('(group3_attention) Starting run number: %i.',stimulus.counter));

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

% Make the gratings for this trial.
grating = mglMakeGrating(10,10,0.5,0,0);
gauss = mglMakeGaussian(10,10,2,2);

makeGrating = @(contrast) round(255*(contrast*grating.*gauss + 1) / 2);

% Select contrast on this trial.
contrast = stimulus.baseContrast;
task.thistrial.leftContrastDelta = task.thistrial.leftStimPresent*stimulus.contrastDelta;
task.thistrial.rightContrastDelta = task.thistrial.rightStimPresent*stimulus.contrastDelta;

% Left stimulus
stimulus.live.tex11 = mglCreateTexture(makeGrating(contrast));
newContrast = contrast + task.thistrial.leftContrastDelta;
stimulus.live.tex12 = mglCreateTexture(makeGrating(newContrast));

% Right stimulus
stimulus.live.tex21 = mglCreateTexture(makeGrating(contrast));
newContrast = contrast + task.thistrial.rightContrastDelta;
stimulus.live.tex22 = mglCreateTexture(makeGrating(newContrast));

%disp trial parameters each trial
sides = {'left', 'right'};
disp(sprintf('(group3_attention) Trial %d - Cue %s; Left: %d; Right %d', task.trialnum, sides{task.thistrial.cueSide}, task.thistrial.leftStimPresent, task.thistrial.rightStimPresent));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)

global stimulus

% Select image parameters: size, eccentricity, and location
ecc = stimulus.eccentricity;
imSz = stimulus.imSize;

% After the response segment (at start of feedback segment), decide whether they get a reward or not.
if task.thistrial.thisseg == stimulus.seg.feedback
  task.thistrial.chosenSide = task.thistrial.response; % Save which side they responded for
  if isnan(task.thistrial.response)
    task.thistrial.response = 1;
  end
  if task.thistrial.cueSide == 1
    task.thistrial.correct = (task.thistrial.response-1 == task.thistrial.leftStimPresent);
  else
    task.thistrial.correct = (task.thistrial.response-1 == task.thistrial.rightStimPresent);
  end
  if task.thistrial.correct
    disp(sprintf('(group3_attention) Correct response'));
  else
    disp(sprintf('(group3_attention) Incorrect response'));
  end
end

for i = 1:2
  mglClearScreen(0.5);
  upFix(stimulus, stimulus.colors.white);
  if task.thistrial.thisseg == stimulus.seg.stim
    mglBltTexture(stimulus.live.tex11, [-ecc, 0]);
    mglBltTexture(stimulus.live.tex21,[ecc, 0]);
  elseif task.thistrial.thisseg == stimulus.seg.increment
    mglBltTexture(stimulus.live.tex12, [-ecc, 0]);
    mglBltTexture(stimulus.live.tex22, [ecc,0]);
  elseif task.thistrial.thisseg == stimulus.seg.cue
    if task.thistrial.cueSide == 1
      mglPolygon([-1,-.75, -2, -.75], [0, .5, 0, -.5], [0,0,0]);
    else
      mglPolygon([1,.75,2,.75], [0,.5,0,-.5], [0,0,0]);
    end
  elseif task.thistrial.thisseg == stimulus.seg.feedback
    if task.thistrial.correct
      upFix(stimulus, stimulus.colors.green);
    else
      upFix(stimulus, stimulus.colors.red);
    end
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

