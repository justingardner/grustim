% 
%
%        $Id: veins.m,v 1.9 2008/09/13 16:17:11 justin Exp $
%      usage: taskTemplateStaticStaircase
%         by: justin gardner
%       date: 01/27/07
%  copyright: (c) 2007 Justin Gardner (GPL see mgl/COPYING)
%    purpose: oriented grating stimulus
%
function myscreen = veins

% check arguments
if ~any(nargin == [0])
  help veins
  return
end

% check to see whether screen is still open
global stimulus;
if ~isfield(stimulus,'texturesCreated')
  disp(sprintf('(veins) Recreating stimulus'));
  stimulus = [];
  stimulus.texturesCreated = 0;
end
global MGL;
if ~isfield(MGL,'displayNumber') || (MGL.displayNumber == -1)
  % no screen open, force stimulus to recreate textures
  stimulus.texturesCreated = 0;
end

% initalize the screen
myscreen.background = 'gray';
myscreen.autoCloseScreen = 0;
myscreen = initScreen(myscreen);

% set the first task to be the fixation staircase task
[task{1} myscreen] = fixStairInitTask(myscreen);

% update stimulus every 250 ms
frameLen = 0.250;

% stimulus is on for 30 seconds
stimLen = 20;

% first phase, we show randomized orientation
%orientations = 15:60:135;
orientations = 11.25:22.5:180;
task{2}{1}.waitForBacktick =1;
task{2}{1}.seglen = frameLen * ones(1,(stimLen/10)/frameLen);
task{2}{1}.parameter.orientation = [orientations;orientations];
task{2}{1}.random = 1;
task{2}{1}.numTrials = 10;

stimLen = stimLen-frameLen;
% second phase actually show orientations
task{2}{2}.seglen = [frameLen * ones(1,stimLen/frameLen) frameLen/2];
task{2}{2}.synchToVol = zeros(1,(stimLen/frameLen)+1);
task{2}{2}.synchToVol(end) = 1;
task{2}{2}.parameter.orientation = [orientations;orientations];
task{2}{2}.random = 1;

% initialize the task
for phaseNum = 1:length(task)
  [task{2}{phaseNum} myscreen] = initTask(task{2}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end

% register stimulus in myscreen
myscreen = initStimulus('stimulus',myscreen);
% do our initialization which creates a grating
stimulus = myInitStimulus(stimulus,myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task)) && ~myscreen.userHitEsc
  % update the task
  [task{2} myscreen phaseNum] = updateTask(task{2},myscreen,phaseNum);
  % update the fixation task
  [task{1} myscreen] = updateTask(task{1},myscreen,1);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;

% randomize the current phase of the stimulus
stimulus.phaseNum = ceil(rand(1)*stimulus.numPhases);
if (task.thistrial.thisseg == 2)
  keyboard
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus;

% clear the screen
mglClearScreen;

% draw the left grating
mglStencilSelect(2);
mglBltTexture(stimulus.tex(stimulus.phaseNum),[-stimulus.centerX 0],0,0,task.thistrial.orientation(1));

% draw the right grating
mglStencilSelect(1);
mglBltTexture(stimulus.tex(stimulus.phaseNum),[stimulus.centerX 0],0,0,task.thistrial.orientation(2));

% return to unstenciled drawing
mglStencilSelect(0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets subject  response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task, myscreen)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen,task)

% spatial frequency
stimulus.sf = 1.5;

% which phases we will have
stimulus.numPhases = 16;
stimulus.phases = 0:(360-0)/stimulus.numPhases:360;

% use circular aperature or not
stimulus.circularAperture = 0;

% use single side
stimulus.singleAperture = 0;

% stimulus contrast
stimulus.contrast = 0.9;

% size of stimulus
stimulus.width = 0.5*floor(myscreen.imageHeight/0.5);
stimulus.height = stimulus.width;

% choose whether we want a gabor/simple cutout
stimulus.gabor = 0;

% chose a sin or square
stimulus.square = 1;

% initial phase number
stimulus.phaseNum = 1;

% check against old parameters to see if we have to
% recreate the stimulus
fieldnames = {'gabor','square','width','height','phases','sf'};
if isfield(stimulus,'old') && isfield(stimulus,'texturesCreated') && stimulus.texturesCreated
  for i = 1:length(fieldnames)
    if ~isequal(stimulus.(fieldnames{i}),stimulus.old.(fieldnames{i}))
      % set init to 0, meaning that we have to init the gratings again
      stimulus.init = 0;
    end
  end
end

if ~stimulus.texturesCreated
  disppercent(-inf,'Making gratings');
  for i = 1:length(stimulus.phases)
    % make a grating
    if stimulus.square
      grating(:,:,1) = 255*(stimulus.contrast*sign(makeGrating(stimulus.width,stimulus.height,stimulus.sf,0,stimulus.phases(i)))+1)/2;
    else
      grating(:,:,1) = 255*(stimulus.contrast*makeGrating(stimulus.width,stimulus.height,stimulus.sf,0,stimulus.phases(i))+1)/2;
    end
    % set g and b channels
    grating(:,:,2) = grating(:,:,1);
    grating(:,:,3) = grating(:,:,1);
    % either make a gabor
    if stimulus.gabor
      grating(:,:,4) = 255*makeGaussian(stimulus.width,stimulus.height,stimulus.width/6,stimulus.height/6);
      % or a simple cutout
    else
      grating(:,:,4) = 255*(makeGaussian(stimulus.width,stimulus.height,stimulus.width/2,stimulus.height/2)>exp(-1/2));
    end
    % make it into a texture
    stimulus.tex(i) = mglCreateTexture(grating);
    disppercent(i/length(stimulus.phases));
  end
  disppercent(inf);

  % create stencil for left and right hand side

  % stencil locations
  stimulus.stencilAngle = 10;
  stencilRadius = max(stimulus.width,stimulus.height)+1;

  % screen places
  screenRight = myscreen.imageWidth/2;
  screenTop = myscreen.imageHeight/2;
  screenBottom = -myscreen.imageWidth/2;
  screenLeft = -myscreen.imageHeight/2;
  
  % make right cutout
  x = [0 cos(d2r(90-stimulus.stencilAngle))*stencilRadius screenRight ...
       screenRight cos(d2r(270+stimulus.stencilAngle))*stencilRadius 0];
  y = [0 sin(d2r(90-stimulus.stencilAngle))*stencilRadius screenTop ...
       screenBottom sin(d2r(270+stimulus.stencilAngle))*stencilRadius 0];

  % or size of circular cutot
  if stimulus.circularAperture
    fixDiskSize = 1;
    distFromEdge = 0.5;
    circleSize = (myscreen.imageWidth/2) - fixDiskSize - distFromEdge;
    stimulus.centerX = fixDiskSize+circleSize/2;
    circleSize = [circleSize circleSize];
  else
    stimulus.centerX = 0;
  end
  
  % right stencil
  mglStencilCreateBegin(1);
  if stimulus.circularAperture
    mglFillOval(stimulus.centerX,0,circleSize,[1 1 1]);
    mglGluDisk(stimulus.centerX,0,circleSize(1)/2,[1 1 1],128);
  else
    % for the single grating
    if stimulus.singleAperture
      mglGluDisk(0,0,myscreen.imageWidth,[1 1 1],128);
    else
      mglPolygon(x,y,[1 1 1]);
    end
  end
  mglStencilCreateEnd;

  % make left stencil
  mglClearScreen;
  mglStencilCreateBegin(2);
  if stimulus.circularAperture
    mglGluDisk(-stimulus.centerX,0,circleSize(1)/2,[1 1 1],128);
  else
    if ~stimulus.singleAperture
      mglPolygon(-x,y,[1 1 1]);
    end
  end
  mglStencilCreateEnd;
  
end

% stimulus has textures now
stimulus.texturesCreated = 1;

% remember old settings so that we can test whether we will need to
for i = 1:length(fieldnames)
  stimulus.old.(fieldnames{i}) = stimulus.(fieldnames{i});
end
