function [ myscreen ] = texAtt( varargin )
%
% TEXTURE SEARCH 
%  Visual search task using textures 
%
%  Usage: texAtt(varargin)
%  Authors: Akshay Jagadeesh
%  Date: 02/27/2018
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

if ~isempty(mglGetSID) && isdir(sprintf('~/data/texAtt/%s',mglGetSID))
  % Directory exists, check for a stimfile
  files = dir(sprintf('~/data/texAtt/%s/1*mat',mglGetSID));

  if length(files) >= 1
    fname = files(end).name;
    
    s = load(sprintf('~/data/texAtt/%s/%s',mglGetSID,fname));
    stimulus.counter = s.stimulus.counter + 1;
    clear s;
    disp(sprintf('(texAtt) Data file: %s loaded.',fname));
  end
end
disp(sprintf('(texAtt) This is run #%i',stimulus.counter));

%% Setup Screen
myscreen = initScreen('VPixx2');

% set background to grey
myscreen.background = 0.5;

%% Setup missing initial variables
if ~isfield(stimulus,'counter')
  stimulus.counter = 1; % This keeps track of what "run" we are on.
end


%% Initialize Stimulus
myscreen = initStimulus('stimulus',myscreen);
localInitStimulus();
  
% Set response keys
stimulus.responseKeys = [11 12 13 14];
stimulus.responseKeys = [1 2 3];

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

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% [akshay] (1) Edit the code below to specify the length 
%         and name of each segment
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Define stimulus timing
% inf, Cue, Stimulus, Pause, Response, Feedback
task{1}.segmin = [inf, 1.0, 0.6, 0.4, 1.2, .2];
task{1}.segmax = [inf, 1.0, 0.6, 0.4, 1.2, .2];
stimulus.seg = {};
stimulus.seg.fix = 1;
stimulus.seg.cue = 2;
stimulus.seg.stim = 3;
stimulus.seg.blank = 4;
stimulus.seg.response = 5;
stimulus.seg.feedback = 6;

if stimulus.noeye==1
  task{1}.segmin(1) = 0.1;
  task{1}.segmax(1) = 0.1;
end

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% [akshay] (2) Edit the code below to store some important variables in the 
%               stimulus struct. Remember that stimulus is just a place to store
%               useful variables -- it doesn't really interact with the mgl code.
%         - Important variables such as a list of the images you want to load, and the directory they're saved in
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% Task important variables
stimulus.imNames = {'im13', 'im18', 'im23', 'im30', 'im38', 'im48', 'im52', 'im56', 'im60', 'im71', 'im99', 'im327', 'im336', 'im393', 'im402'};
stimulus.layerNames = {'pool1', 'pool2', 'pool4'};
stimulus.stimDir = '~/proj/TextureSynthesis/stimuli/textures/tex_bw';
% Tess 12-13-18
stimulus.origDir = '~/proj/TextureSynthesis/stimuli/textures/orig_bw';
stimulus.imSize = 6;
stimulus.eccentricity = 8;
stimulus.poolSizes = {'1x1', '2x2', '3x3', '4x4'};
stimulus.cueEcc = 4;
stimulus.live.mask = imread('~/proj/TextureSynthesis/stimuli/Flattop8.tif');


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% [akshay] (3) Edit the code below to specify task trial parameters which will be selected in each
%         trial. For example, the image family or the layer can be randomly assigned.
%        -  MGL will attempt to assign trial parameters such that there are equal number of each trial type.
%        -  Remember that in later functions, these variables can be access from task.thistrial.layer or whatever.   
%        -  The getResponse is a list of length nSegs and is used to specify which segments MGL should expect to be 
%           listening for a response. 
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Trial parameters
task{1}.parameter.imgFam = 1:length(stimulus.imNames);
task{1}.parameter.layer = 1:length(stimulus.layerNames);
task{1}.parameter.poolSize = 1:length(stimulus.poolSizes);

task{1}.synchToVol = zeros(size(task{1}.segmin));
task{1}.getResponse = zeros(size(task{1}.segmin));
task{1}.getResponse(stimulus.seg.response)=1;

% Make numTrials some multiple of number of TrialTypes (in this case, # of attentional conditions x # of layers).
task{1}.numTrials = 144; %where does this number come from?
task{1}.random = 1;


%!!!!!!!!!!!!!!!!!!!!!!!!!/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% [akshay] (4) Edit the code below to define and initialize variables that you may want to 
%          calculate later on (e.g. in startTrialCallback or startSegmentCallback).
%       - you may want to create a variable here to keep track of which location is the cued location.
%       - or to keep track of what the response of the subject is.
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Task variables to be calculated later
task{1}.randVars.calculated.targetSide = NaN;
task{1}.randVars.calculated.isCueFocal = NaN;
task{1}.randVars.calculated.targetPosition = NaN;
task{1}.randVars.calculated.detected = 0; % did they see the grating
task{1}.randVars.calculated.dead = 0;
task{1}.randVars.calculated.visible = 1;
task{1}.randVars.calculated.correct = NaN;
task{1}.randVars.calculated.isTargetSame = NaN;

%%%%%%%%%%%%%%%% MGL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialize task and begin update loop
[task{1}, myscreen] = initTask(task{1},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);

% Run the eye calibration
myscreen = eyeCalibDisp(myscreen);

% let the user know
disp(sprintf('(texAtt) Starting run number: %i.',stimulus.counter));

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

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% [akshay] (5) Modify the code below to specify the directories from which we will be loading images.
%       - This is where we will preload all 8 images (4 images x 2 intervals) at the start of each trial.
%       - Use imread to load the images and then pass that into genTexFromIm, and store the output into a field of
%         the stimulus.live struct.
%
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% directories
texDir = stimulus.stimDir;
origDir = stimulus.origDir;

%% Load all 8 images for this trial
imName = stimulus.imNames{task.thistrial.imgFam};
layer = stimulus.layerNames{task.thistrial.layer};
poolSize = stimulus.poolSizes{task.thistrial.poolSize};

% Randomly select 4 images to display on this trial.
%smpls = randperm(15,8);
% smpls = randperm(10,8);
smpls = randperm(10,4);
% isSame = randi([0,1], 1, 4);


stimulus.live.left_oddball = genTexFromIm(imread(sprintf('%s/%s.png', origDir, imName)), stimulus.live.mask);
stimulus.live.right_oddball = genTexFromIm(imread(sprintf('%s/%s.png', origDir, imName)), stimulus.live.mask);

stimulus.live.dl1 = genTexFromIm(imread(sprintf('%s/%s_%s_%s_smp%i.png', texDir, poolSize, layer, imName, smpls(1))), stimulus.live.mask);
stimulus.live.dl2 = genTexFromIm(imread(sprintf('%s/%s_%s_%s_smp%i.png', texDir, poolSize, layer, imName, smpls(2))), stimulus.live.mask);
stimulus.live.dr1 = genTexFromIm(imread(sprintf('%s/%s_%s_%s_smp%i.png', texDir, poolSize, layer, imName, smpls(3))), stimulus.live.mask);
stimulus.live.dr2 = genTexFromIm(imread(sprintf('%s/%s_%s_%s_smp%i.png', texDir, poolSize, layer, imName, smpls(4))), stimulus.live.mask);



%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% [akshay] (6) Edit and add code below to set the values of the randVars.calculated variables which you had defined above.
%       - For example, you may want to use this space to choose which side
%       (right or left) will be the response-cued direction,
%       as well as whether this trial is a distributed cue or a focal cue trial.
%
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Select target position and target size.
% task.thistrial.isTargetSame = isSame(1);
task.thistrial.targetSide = randi(2); % 1 is left 2 is right
task.thistrial.isCueFocal = randi([0, 1]);
task.thistrial.targetLoc = randi(3);
task.thistrial.otherOddLoc = randi(3);

% Disp trial parameters each trial
disp(sprintf('Trial %d - Image %s, Layer %s, Pooling Region Size: %s, targetSide: %i, isCueFocal: %i', task.trialnum, imName, layer,poolSize, task.thistrial.targetSide, task.thistrial.isCueFocal));
if task.trialnum > 1
    disp(sprintf('--Target on Last Trial: %g, Response on last trial: %g, LastTrialCorrect?: %g', task.lasttrial.targetLoc, task.lasttrial.response, task.lasttrial.correct));
end
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

% If in feedback segment, change the color of the cross to green or red.
if task.thistrial.thisseg == stimulus.seg.feedback
  if task.thistrial.response == task.thistrial.targetLoc
    task.thistrial.correct = 1;
  else
    task.thistrial.correct = 0;
  end
end

% Select image parameters: size, eccentricity, and location
imSz = stimulus.imSize; % Size in Degrees of Visual Angle to display image
ecc = stimulus.eccentricity; % Eccentricity to display image at
cueX = stimulus.cueEcc/sqrt(2);
locations = [-ecc ecc; -ecc 0; -ecc -ecc; ecc ecc; ecc 0; ecc -ecc];
% locations = [ecc ecc; -ecc ecc; -ecc -ecc; ecc -ecc];

% if cue is left
if task.thistrial.targetSide == 1
    % left
    leftOddPos = task.thistrial.targetLoc; %targetLoc is left
    rightOddPos = task.thistrial.otherOddLoc+3;
elseif task.thistrial.targetSide == 2
    % right
    rightOddPos = task.thistrial.targetLoc + 3; %shift targetLoc by 3
    leftOddPos = task.thistrial.otherOddLoc;   
end

distLocations = setdiff(1:6, [rightOddPos, leftOddPos]);



%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% [akshay] (6) Edit the code below to specify what to draw onto the screen at the start of each segment
%          - Recall that MGL has a front buffer and a back buffer, so the code below loops through
%          both buffers, clears whatever was on the screen from the last segment, and draws something else.
%          - There may be some things that you want to draw on every segment (eg fixation cross), and other 
%          things that you may want to draw only on some segments (e.g. attentional focal/distributed cue),
%          and other things you may only want to draw on one segment (e.g. interval 1/2 stimuli).
%
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if task.thistrial.isCueFocal == 1 || task.thistrial.thisseg == stimulus.seg.feedback || task.thistrial.thisseg == stimulus.seg.response
  stimulus.live.draw4 = 0;
  if task.thistrial.thisseg == stimulus.seg.response
      trialColor = stimulus.colors.blue;
  elseif task.thistrial.thisseg == stimulus.seg.feedback
      if task.thistrial.correct
        trialColor = stimulus.colors.green;
      else
        trialColor = stimulus.colors.red;
      end
  else
      trialColor = stimulus.colors.black;
  end
else
  stimulus.live.draw4 = 1;
end

% cueXLocs = [cueX, cueX, cueX-1; -cueX, -cueX, -cueX+1; -cueX, -cueX, -cueX+1; cueX-1 cueX, cueX];
% cueYLocs = [cueX, cueX-1, cueX; cueX, cueX-1, cueX; -cueX+1, -cueX, -cueX; -cueX, -cueX, -cueX+1];

% cueXLocs = [cueX, cueX-.25, cueX-.75; -cueX, -cueX+.25, -cueX+.75; -cueX+.25, -cueX, -cueX+.75; cueX-.75 cueX, cueX-.25];
% cueYLocs = [cueX, cueX-.75, cueX-.25; cueX, cueX-.75, cueX-.25; -cueX+.75, -cueX, -cueX+.25; -cueX+.25, -cueX, -cueX+.75];

cueXLocs = [-cueX, -cueX, -cueX-1; cueX, cueX, cueX+1];
cueYLocs = [-.25, .25, 0; -.25, .25, 0];

stimulus.live.cueXLocs = cueXLocs;
stimulus.live.cueYLocs = cueYLocs;
for i = 1:2
  mglClearScreen(0.5);
  if task.thistrial.thisseg == stimulus.seg.stim
    mglBltTexture(stimulus.live.left_oddball,  [locations(leftOddPos, :) imSz imSz]);
    mglBltTexture(stimulus.live.right_oddball, [locations(rightOddPos, :) imSz imSz]);
    mglBltTexture(stimulus.live.dl1, [locations(distLocations(1), :) imSz imSz]);
    mglBltTexture(stimulus.live.dl2, [locations(distLocations(2), :) imSz imSz]);
    mglBltTexture(stimulus.live.dr1, [locations(distLocations(3), :) imSz imSz]);
    mglBltTexture(stimulus.live.dr2, [locations(distLocations(4), :) imSz imSz]);
    
    upFix(stimulus, stimulus.colors.white);
    
  else
    upFix(stimulus, stimulus.colors.black);
  end

  if stimulus.live.draw4
    %draw 2 triangles
    mglPolygon(cueXLocs(1,:), cueYLocs(1,:), stimulus.colors.black);
    mglPolygon(cueXLocs(2,:), cueYLocs(2,:), stimulus.colors.black);
  else 
    %draw 1 triangle
    mglPolygon(cueXLocs(task.thistrial.targetSide,:), cueYLocs(task.thistrial.targetSide,:), trialColor);
  end
  mglFlush
end
stimulus.live.firstTime = 0;

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

validResponse = any(task.thistrial.whichButton == stimulus.responseKeys);
if validResponse
  if stimulus.live.gotResponse==0
    task.thistrial.detected = 1;
    task.thistrial.response = task.thistrial.whichButton;
    stimulus.live.fix = 0;
  else
    disp(sprintf('Subject responded multiple times: %i',stimulus.live.gotResponse));
  end
  disp('jumping segment');
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

