% introCogNeuroAttention
%
%      usage: introCogNeuroAttention
%         by: justin gardner
%       date: 02/06/17
%  copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%    purpose: program for demonstration of attention tasks for psych 50
% 
%             First run with no arguments:
%             introCogNeuroAttention
%
%             This will show a single motion patch
%             on the left which moves inwards or outwards. Task is to hit 1 for inwards
%             and 2 for outwards. It will step down slowly in coherence - when you get
%             to a coherence where people start to get confused, call that threshold (it's
%             not so important to get exactly right)
%
%             Now you can run the full task with that coherence threshold. Say the coherence
%             threshold was 0.3
%             introCogNeuroAttention(.3)
%
%             This is the main task. Two patches come up. The fixation arm will point red
%             towards one or both of the patches. Tell them to attend to the location indicated
%             by the cue - if both arms are red, then they should distribute their attention
%             across both patches. After the patches disappear, one arm of the fixation cross
%             will turn green and they will have to report about that patch (same task - in or outward
%             motion).
%
%
%             To discuss about suppressing stimuli, you can add distractors:
%             introCogNeuroAttention(0.3,1);
% 
%             This is exactly the same task except there are 4 other distractor patches. They
%             should simply be ignored - the hope is that people will find this difficult - and
%             can talk about selection / biased competition models
%       
%
function myscreen = introCogNeuro_attention(thresholdCoherence,useDistractors)

if nargin == 0
  taskType = 'coherenceThreshold';
  thresholdCoherence = [];
elseif (nargin >= 1) && isscalar(thresholdCoherence)
  % allow values between 1 and 100 (or 0 - 1)
  if (thresholdCoherence > 1) thresholdCoherence = thresholdCoherence / 100;end
  if thresholdCoherence > 1,thresholdCoherence = 1;end
  taskType = 'attention';
else
  help introCogNeuro_attention
  return
end  
if nargin < 2,useDistractors = 0;end

% initalize the screen
myscreen.autoCloseScreen = 1;
myscreen.displayname = 'projector';
myscreen.background = 'black';
myscreen = initScreen(myscreen);

% task parameters
task{1}.seglen = [1 0.5 1 5];
task{1}.getResponse = [0 0 0 1];
task{1}.parameter.cue = [-1 1 0];
task{1}.parameter.coherence = thresholdCoherence;
task{1}.random = 1;
task{1}.fudgeLastVolume = 1;

% init the stimulus
clear global stimulus;
global stimulus;
myscreen = initStimulus('stimulus',myscreen);

% setup dot patches
global gStencilNumber;gStencilNumber = 1;
dots(1).type = 'Opticflow';
dots(1).mask = 1;
dots(1).xcenter = -12;
dots(1).ycenter = 0;
dots(1).coherence = 1;
dots(1).rmax = 9.5;
dots(2) = dots(1);
dots(2).xcenter = 12;

% size of fixation
stimulus.fixSize = 4;
stimulus.fixLineWidth = 3;

% initialize some other parameters
stimulus.responseText = [];
stimulus.taskType = taskType;

% if we are running coherence threshold then set things up differently
if strcmp(stimulus.taskType,'coherenceThreshold')
  task{1}.parameter.coherence = [1 1 1 0.8 0.5 0.5 0.4 0.4 0.3 0.3 0.2 0.2 0.1];
  task{1}.random = 0;
  task{1}.parameter.cue = -1;
  % only display left patch
  dots = dots(1);
end

% makes some distractors
if useDistractors
  eccentricity = 12;
  xcenter = cos(pi/3)*eccentricity;
  ycenter = sin(pi/3)*eccentricity;
  dots(1).rmax = 5.5;
  dots(2).rmax = 5.5;
  dots(3) = dots(2);
  dots(3).xcenter = -xcenter;
  dots(3).ycenter = -ycenter;

  dots(4) = dots(2);
  dots(4).xcenter = -xcenter;
  dots(4).ycenter = ycenter;

  dots(5) = dots(2);
  dots(5).xcenter = xcenter;
  dots(5).ycenter = -ycenter;

  dots(6) = dots(2);
  dots(6).xcenter = xcenter;
  dots(6).ycenter = ycenter;
end

% initialize dot patches
for iDots = 1:length(dots)
  stimulus.dots(iDots) = initDots(dots(iDots),myscreen);
end

% initialize task
[task{1} myscreen] = initTask(task{1},myscreen,@startSegmentCallback,@updateScreenCallback,@responseCallback);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
  % update the dots
  [task myscreen] = updateTask(task,myscreen,1);
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

% first segment
if (task.thistrial.thisseg == 1)
  % set the direction of the patches randomly
  nDots = length(stimulus.dots);
  stimulus.dotDirs = -1 + 2*(rand(1,nDots)>0.5);
  % actually set the speed
  for iDots = 1:nDots
    stimulus.dots(iDots) = setDotsSpeedOpticflow(stimulus.dots(iDots),abs(stimulus.dots(iDots).speed)*stimulus.dotDirs(iDots),myscreen);;
  end
  % decide which is the target to respond to
  if task.thistrial.cue == 0
    stimulus.responseSide = -1 + 2*(rand(1,1)>0.5);
  else
    stimulus.responseSide = task.thistrial.cue;
  end
  % init response text
  if ~isempty(stimulus.responseText)
    mglDeleteTexture(stimulus.responseText);
    stimulus.responseText = [];
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task, myscreen)

global stimulus
mglClearScreen;

% draw fixation
mglFixationCross(stimulus.fixSize,stimulus.fixLineWidth,1);

if any(task.thistrial.thisseg == [2 3])
  % draw cue
  switch task.thistrial.cue
    case {-1}
     mglLines2(0,0,-stimulus.fixSize/2,0,stimulus.fixLineWidth,[1 0 0]);
   case {1}
     mglLines2(0,0,stimulus.fixSize/2,0,stimulus.fixLineWidth,[1 0 0]);
   case {0}
     mglLines2(-stimulus.fixSize/2,0,stimulus.fixSize/2,0,stimulus.fixLineWidth,[1 0 0]);
  end
end

if task.thistrial.thisseg == 3
  % draw the dots in the stimulus interval
  for iDots = 1:length(stimulus.dots)

    % update dot positions
    stimulus.dots(iDots) = feval(sprintf('updateDots%s',stimulus.dots(iDots).type),stimulus.dots(iDots),task.thistrial.coherence,myscreen);
    
    % draw the dots (with stencil)
    if stimulus.dots(iDots).mask,mglStencilSelect(stimulus.dots(iDots).stencilNumber);end
    mglPoints2(stimulus.dots(iDots).x+stimulus.dots(iDots).xcenter,stimulus.dots(iDots).y+stimulus.dots(iDots).ycenter,stimulus.dots(iDots).dotsize,[1 1 1]);
    if stimulus.dots(iDots).mask,mglStencilSelect(0);end
  end
end

% draw response cue
if task.thistrial.thisseg == 4
  % draw responseCue
  switch stimulus.responseSide
    case {-1}
     mglLines2(0,0,-stimulus.fixSize/2,0,stimulus.fixLineWidth,[0 1 0]);
   case {1}
     mglLines2(0,0,stimulus.fixSize/2,0,stimulus.fixLineWidth,[0 1 0]);
  end
end

if task.thistrial.thisseg == 4
  if ~isempty(stimulus.responseText)
    mglBltTexture(stimulus.responseText,[0 stimulus.fixSize/2+0.5],0,1);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called when the subject responds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task, myscreen)

global stimulus;

if any(task.thistrial.whichButton == [1 2]);
  responseMap = [-1 1];
  % patch direction
  if stimulus.dotDirs(1+(stimulus.responseSide+1)/2) == -1
    targetDir = 1;
    targetDirStr = 'outward';
  else
    targetDir = -1;
    targetDirStr = 'inward';
  end
  % see if we are correct or not
  if responseMap(task.thistrial.whichButton) == targetDir
    correctIncorrect = 1;
  else
    correctIncorrect = 0;
  end
  if stimulus.responseSide == 1
    targetLocStr = 'Right';
  else
    targetLocStr = 'Left';
  end
  % display what happened
  if correctIncorrect
    mglTextSet('Helvetica',48,[0 1 0],0,0,0);
    % see which task we are doing
    if strcmp(stimulus.taskType,'coherenceThreshold')
      stimulus.responseText = mglText(sprintf('Correct: %s at %i%% coherence',targetDirStr,ceil(task.thistrial.coherence*100)));
    else
      stimulus.responseText = mglText(sprintf('Correct: %s patch was %s',targetLocStr,targetDirStr));
    end
  else
    mglTextSet('Helvetica',48,[1 0 0],0,0,0);
    % see which task we are doing
    if strcmp(stimulus.taskType,'coherenceThreshold')
      stimulus.responseText = mglText(sprintf('Wrong: %s at %i%% coherence',targetDirStr,ceil(task.thistrial.coherence*100)));
    else
      stimulus.responseText = mglText(sprintf('Wrong: %s patch was %s',targetLocStr,targetDirStr));
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initDots(dots,myscreen)

% convert the passed in parameters to real units
if ~isfield(dots,'rmax'), dots.rmax = 7.5;,end
if ~isfield(dots,'xcenter'), dots.xcenter = 0;,end
if ~isfield(dots,'ycenter'), dots.ycenter = 0;,end
if ~isfield(dots,'dotsize'), dots.dotsize = 4;,end
if ~isfield(dots,'density'), dots.density = 5;,end
if ~isfield(dots,'coherence'), dots.coherence = 1;,end
if ~isfield(dots,'speed'), dots.speed = 4;,end
if ~isfield(dots,'dir'), dots.dir = 0;,end
if ~isfield(dots,'mask'), dots.mask = 1;,end

% update the dots
dots = feval(sprintf('initDots%s',dots.type),dots,myscreen);

% set color
dots.color = ones(dots.n,1);

% create stencil
if dots.mask
  % give the dots a stencil number
  global gStencilNumber;
  dots.stencilNumber = gStencilNumber;
  gStencilNumber = gStencilNumber + 1;
  % and create the stencil mask
  mglClearScreen;
  mglStencilCreateBegin(dots.stencilNumber);
  % and draw that oval
  mglGluDisk(dots.xcenter,dots.ycenter,[dots.rmax dots.rmax],[1 1 1],60);
  mglStencilCreateEnd;
  mglClearScreen;
end

dots = feval(sprintf('updateDots%s',dots.type),dots,dots.coherence,myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the dots speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = setDotsSpeedLinear(dots,speed,myscreen)

% get the step size
dots.speed = speed;
dots.stepsize = speed/myscreen.framesPerSecond;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the dots direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = setDotsDirLinear(dots,direction,myscreen)

% get the step size
dots.dir = direction;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots for linear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updateDotsLinear(dots,coherence,myscreen)

% get the dots step
dots.xstep = cos(pi*dots.dir/180)*dots.stepsize;
dots.ystep = sin(pi*dots.dir/180)*dots.stepsize;

% pick a random set of dots
dots.coherent = rand(1,dots.n) < coherence;

% now move those dots in the right direction
dots.x(dots.coherent) = dots.x(dots.coherent)+dots.xstep;
dots.y(dots.coherent) = dots.y(dots.coherent)+dots.ystep;

% randomwalk rule
thisdir = rand(1,sum(~dots.coherent))*2*pi;
dots.x(~dots.coherent) = dots.x(~dots.coherent)+cos(thisdir)*dots.stepsize;
dots.y(~dots.coherent) = dots.y(~dots.coherent)+sin(thisdir)*dots.stepsize;

% movshon noise
%if dots.movshonNoise
  dots.x(~dots.coherent) = rand(1,sum(~dots.coherent))*dots.width;
  dots.y(~dots.coherent) = rand(1,sum(~dots.coherent))*dots.height;
%end

% make sure we haven't gone off the patch
% do the dots separately for left and right hand side
dots.x(dots.x < dots.xmin) = dots.x(dots.x < dots.xmin)+dots.width;
dots.x(dots.x > dots.xmax) = dots.x(dots.x > dots.xmax)-dots.width;
dots.y(dots.y < dots.ymin) = dots.y(dots.y < dots.ymin)+dots.height;
dots.y(dots.y > dots.ymax) = dots.y(dots.y > dots.ymax)-dots.height;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots for opticflow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updateDotsOpticflow(dots,coherence,myscreen)

% get the coherent and incoherent dots
%if (dots.coherency ~= coherence)
  dots.incoherent = rand(1,dots.n) > coherence;
  dots.incoherentn = sum(dots.incoherent);
  dots.coherent = ~dots.incoherent;
  dots.coherency = coherence;
  % generate a random transformation matrix for each incoherent point
  dots.randT = rand(3,dots.incoherentn)-0.5;
  % and normalize the transformation to have the same length
  % (i.e. speed) as the real transformation matrix
  dots.randT = sqrt(sum(dots.T.^2))*dots.randT./([1 1 1]'*sqrt(sum(dots.randT.^2)));
%end

% update relative position of dots in 3-space to observer
dots.X(dots.coherent) = dots.X(dots.coherent)-dots.T(1);
dots.Y(dots.coherent) = dots.Y(dots.coherent)-dots.T(2);
dots.Z(dots.coherent) = dots.Z(dots.coherent)-dots.T(3);

% now move the incoherent points according to the random trasnformation
dots.X(dots.incoherent) = dots.X(dots.incoherent)-dots.randT(1,:);
dots.Y(dots.incoherent) = dots.Y(dots.incoherent)-dots.randT(2,:);
dots.Z(dots.incoherent) = dots.Z(dots.incoherent)-dots.randT(3,:);

% get all points that have fallen off the screen
offscreen = dots.Z<dots.minZ;

% and put them at the furthest distance
dots.Z(offscreen) = dots.maxZ;

% get all points that have fallen out of view
offscreen = dots.Z>dots.maxZ;
% and move them to the front plane
dots.Z(offscreen) = dots.minZ;

% put points fallen off the X edge back
offscreen = dots.X < -dots.maxX;
dots.X(offscreen) = dots.X(offscreen)+2*dots.maxX;
offscreen = dots.X > dots.maxX;
dots.X(offscreen) = dots.X(offscreen)-2*dots.maxX;

% put points fallen off the Y edge back
offscreen = dots.Y < -dots.maxY;
dots.Y(offscreen) = dots.Y(offscreen)+2*dots.maxY;
offscreen = dots.Y > dots.maxY;
dots.Y(offscreen) = dots.Y(offscreen)-2*dots.maxY;

% project on to screen
dots.xproj = dots.f*dots.X./dots.Z;
dots.yproj = dots.f*dots.Y./dots.Z;

% stuff to compute median speed
dots.oldx = dots.x;
dots.oldy = dots.y;

% get actual screen coordinates
dots.x = dots.xproj*myscreen.imageWidth;
dots.y = dots.yproj*myscreen.imageHeight;

%medianSpeed = median(sqrt((dots.oldx-dots.x).^2+(dots.oldy-dots.y).^2)*myscreen.framesPerSecond);
%minSpeed = min(sqrt((dots.oldx-dots.x).^2+(dots.oldy-dots.y).^2)*myscreen.framesPerSecond);
%disp(sprintf('min: %f median: %f',minSpeed,medianSpeed));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create dots for linear2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initDotsLinear(dots,myscreen)

% actually a square patch of dots that get stenciled
% so calculate width and height
dots.width = dots.rmax*2;
dots.height = dots.rmax*2;

% get the number of dots
dots.n = round(dots.width*dots.height*dots.density);

% get max and min points for dots
dots.xmin = -dots.width/2;
dots.xmax = dots.width/2;
dots.ymin = -dots.height/2;
dots.ymax = dots.height/2;

% get initial position
dots.x = rand(1,dots.n)*dots.width;
dots.y = rand(1,dots.n)*dots.height;

% get the step size
dots.stepsize = dots.speed/myscreen.framesPerSecond;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the dots speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = setDotsSpeedOpticflow(dots,speed,myscreen)

% get the step size
dots.speed = speed;
dots.T = [0 0 dots.speed/myscreen.framesPerSecond];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the dots direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = setDotsDirOpticflow(dots,direction,myscreen)

% get the step size
dots.T = [0 0 direction*dots.speed/myscreen.framesPerSecond];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create dots for optic flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initDotsOpticflow(dots,myscreen)

% focal length to projection plane
% projection plane is defined to be 
% 1 unit wide and high, so with 
% this focal length, we are looking at
% a view of the world with a 90 deg fov
dots.f = .5;

% translation and rotation matrices
dots.T = [0 0 dots.speed/myscreen.framesPerSecond];
dots.R = [0 0 0];

% maximum depth of points
dots.maxZ = 10;dots.minZ = dots.f;
dots.maxX = 10;
dots.maxY = 10;

% make a brick of points
dots.n = round(myscreen.imageWidth*myscreen.imageHeight*dots.density);

% initial position of dots
dots.X = 2*dots.maxX*rand(1,dots.n)-dots.maxX;
dots.Y = 2*dots.maxY*rand(1,dots.n)-dots.maxY;
dots.Z = (dots.maxZ-dots.minZ)*rand(1,dots.n)+dots.minZ;

% get projection on to plane
dots.xproj = dots.f*dots.X./dots.Z;
dots.yproj = dots.f*dots.Y./dots.Z;

% put into screen coordinates
dots.x = dots.xproj*myscreen.imageWidth;
dots.y = dots.yproj*myscreen.imageHeight;

% set incoherent dots to 0
dots.coherency = 1;
dots.incoherent = rand(1,dots.n) > dots.coherency;
dots.incoherentn = sum(dots.incoherent);
dots.coherent = ~dots.incoherent;

dots.randT = zeros(3,dots.incoherentn);