% 
%
%        $Id: veins.m,v 1.8 2007/09/07 19:46:08 justin Exp $
%      usage: dotsdir
%         by: justin gardner
%       date: 01/27/07
%  copyright: (c) 2007 Justin Gardner (GPL see mgl/COPYING)
%    purpose: 2 direction patches
% Some options: 'taskType=dots' Can be either dots (the default task) or the block localizer ('taskType=blockLoc')
%               Change the position of the fixation cross: 'fixX=0','fixY=0'
%               Or the position of the stimulus: 'centerX=[],'centerY=[]'
%               Or the diameter of the stimulus: 'diameter=[]'
%               Or the difficulty of the fixation task: 'easyFixTask=1'
%
function myscreen = dotsdir(varargin)

taskType=[];
fixX = [];
fixY = [];
centerX = [];
centerY = [];
diameter = [];
easyFixTask = [];
getArgs(varargin,{'taskType=dots','fixX=0','fixY=0','centerX=[]','centerY=[]','diameter=[]','easyFixTask=1'});

% check to see whether screen is still open
global stimulus;
if ~isfield(stimulus,'texturesCreated')
  disp(sprintf('(dotsdir) Recreating stimulus'));
  stimulus = [];
  stimulus.texturesCreated = 0;
end

global MGL;
if ~isfield(MGL,'displayNumber') || (MGL.displayNumber == -1)
  % no screen open, force stimulus to recreate textures
  stimulus.texturesCreated = 0;
end

% initalize the screen
if isfield(stimulus,'type')
  if stimulus.type == 1
    myscreen.background = 'gray';
  else
    myscreen.background = 'black';
  end
else
  myscreen.background = 'black';
end

% init screen
myscreen = initScreen(myscreen);

% set the first task to be the fixation staircase task
global fixStimulus;
fixStimulus.pos = [fixX fixY];
if ~easyFixTask
  % default values
  fixStimulus.diskSize = 0.5;
  fixStimulus.fixWidth = 1;
  fixStimulus.fixLineWidth = 3;
  fixStimulus.stimTime = 0.4;
  fixStimulus.responseTime = 1;
else
  % make cross bigger and task slower
  fixStimulus.diskSize = 0.5;
  fixStimulus.fixWidth = 1+1*easyFixTask;
  fixStimulus.fixLineWidth = 3+2*easyFixTask;
  fixStimulus.stimTime = 0.4+0.4*easyFixTask;
  fixStimulus.responseTime = 1+1*easyFixTask;
end
[task{1} myscreen] = fixStairInitTask(myscreen);

stimLen = 12;

% first phase, we show randomized orientation
%Aorientations = 15:60:135;
orientations = 22.5:22.5:180;
orientations = 15:60:135;
orientations = 0:90:270;
task{2}{1}.waitForBacktick = 1;
task{2}{1}.seglen = [0 11.9];
task{2}{1}.synchToVol = [0 0];
task{2}{1}.parameter.orientation = [orientations;orientations];
task{2}{1}.random = 1;
task{2}{1}.numTrials = 1;

% set the task type
stimulus.taskType = taskType;
stimulus.onoff = 1;

if strcmp(stimulus.taskType,'blockLoc')
  % half the segments will have the stimulsu on, half will have it off
  task{2}{2}.seglen =     repmat(1.5,1,16);
  task{2}{2}.synchToVol = zeros(1,16);
  task{2}{2}.synchToVol(8) = 1;
  task{2}{2}.synchToVol(16) = 1;
  task{2}{2}.parameter.orientation = [orientations;orientations];
  task{2}{2}.random = 1;
  task{2}{2}.waitForBacktick = 1;
else
  % second phase actually show orientations
  task{2}{2}.seglen = [12 5.9];
  task{2}{2}.synchToVol = [0 1];
  task{2}{2}.parameter.orientation = [orientations;orientations];
  task{2}{2}.random = 1;
  task{2}{2}.waitForBacktick = 1;
end
  

% initialize the task
for phaseNum = 1:length(task)
  [task{2}{phaseNum} myscreen] = initTask(task{2}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback,@initTrialCallback);
end

% register stimulus in myscreen
myscreen = initStimulus('stimulus',myscreen);
% do our initialization which creates a grating
stimulus = myInitStimulus(stimulus,myscreen,task,centerX,centerY,diameter);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   initTrialCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = initTrialCallback(task, myscreen)

% get a random orientation sequence for block design
global stimulus;
if strcmp(stimulus.taskType,'blockLoc')
  orientations = task.parameter.orientation(1,:);
  rperm1 = randperm(size(orientations,2));
  rperm2 = randperm(size(orientations,2));
  stimulus.thisOrientationOrder = orientations([rperm1;rperm2]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;

if stimulus.type == 1
  % randomize the current phase of the stimulus
  stimulus.phaseNum = ceil(rand(1)*stimulus.numPhases);
end

% set the orientations in the tast.thistrial variable for blockLoc
if strcmp(stimulus.taskType,'blockLoc')
  % set the orientation
  onum = mod(task.thistrial.thisseg-1,8)+1;
  task.thistrial.orientation = stimulus.thisOrientationOrder(:,onum);
  % set on or off
  stimulus.onoff = (task.thistrial.thisseg > (length(task.seglen)/2));
end

if (task.thistrial.thisphase == 2) && (task.thistrial.thisseg == 1)
  disp(sprintf('Orientations: %s',num2str(task.thistrial.orientation','%i ')));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus;

% clear the screen
mglClearScreen;
if task.thistrial.thisseg == 2,return,end

% return if the stimulus is off
if stimulus.onoff == 0, return,end

% get coherence for this phase
coherence = (task.thistrial.thisphase == 2);

% draw the left grating
%mglStencilSelect(0);
if stimulus.type == 1
  mglBltTexture(stimulus.tex(stimulus.phaseNum),[-stimulus.centerX 0],0,0,task.thistrial.orientation(1));
else
  stimulus.dots{1}.dir = task.thistrial.orientation(1);
  stimulus.dots{1} = updateDots(stimulus.dots{1},coherence,myscreen);
end

% draw the right grating
%mglStencilSelect(0);
if stimulus.type == 1
  mglBltTexture(stimulus.tex(stimulus.phaseNum),[stimulus.centerX 0],0,0,task.thistrial.orientation(2));
else
  stimulus.dots{2}.dir = task.thistrial.orientation(2);
  stimulus.dots{2} = updateDots(stimulus.dots{2},coherence,myscreen);
end

% return to unstenciled drawing
%mglStencilSelect(0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets subject  response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task, myscreen)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen,task,centerX,centerY,diameter)

% 1 is for gratings. 2 is for dots
stimulus.type = 2;

% use circular aperature or not
stimulus.circularAperture = 1;

% use single side
stimulus.singleAperture = 0;

% size of stimulus
stimulus.width = 0.5*floor(myscreen.imageHeight/0.5);
stimulus.height = stimulus.width;

if stimulus.type == 1
  % spatial frequency
  stimulus.sf = 1.5;

  % which phases we will have
  stimulus.numPhases = 16;
  stimulus.phases = 0:(360-0)/stimulus.numPhases:360;

  % stimulus contrast
  stimulus.contrast = 0.9;
  
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
    % stimulus has textures now
    stimulus.texturesCreated = 1;
    % remember old settings so that we can test whether we will need to
    for i = 1:length(fieldnames)
      stimulus.old.(fieldnames{i}) = stimulus.(fieldnames{i});
    end
  end
end

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
  if ~isempty(diameter)
    circleSize = diameter;
  else
    circleSize = (myscreen.imageWidth/2) - fixDiskSize - distFromEdge;
  end
  if ~isempty(centerX)
    stimulus.centerX = centerX;
  else
    stimulus.centerX = fixDiskSize+circleSize/2;
  end
  if ~isempty(centerY)
    stimulus.centerY = centerY;
  else
    stimulus.centerY = 0;
  end
  circleSize = [circleSize circleSize];
else
  stimulus.centerX = 0;
end

% right stencil
mglClearScreen;
mglStencilCreateBegin(1);
if stimulus.circularAperture
  mglGluDisk(stimulus.centerX,stimulus.centerY,circleSize(1)/2,[1 1 1],128);
else
  if ~stimulus.singleAperture
    mglPolygon(x,y,[1 1 1]);
  end
end
mglStencilCreateEnd;

% make left stencil
mglClearScreen;
mglStencilCreateBegin(2);
if stimulus.circularAperture
  mglGluDisk(-stimulus.centerX,stimulus.centerY,circleSize(1)/2,[1 1 1],128);
else
  if ~stimulus.singleAperture
    mglPolygon(-x,y,[1 1 1]);
  end
end
mglStencilCreateEnd;

% right stencil
mglClearScreen;
mglStencilCreateBegin(3);
if stimulus.circularAperture
  mglGluDisk(stimulus.centerX,stimulus.centerY,circleSize(1)/2,[1 1 1],128);
else
  if ~stimulus.singleAperture
    mglPolygon(x,y,[1 1 1]);
  end
end
mglStencilCreateEnd;

% create dots stimulus here, easier to do once we have stencil info
if stimulus.type == 2
  % create dots stimulus.
  stimulus.dots = {};
  dots.rmax = circleSize(1)/2;
  dots.xcenter = -stimulus.centerX;
  dots.ycenter = stimulus.centerY;
  stimulus.dots{1} = initDots(myscreen,dots);
  dots.xcenter = stimulus.centerX;
  dots.ycenter = stimulus.centerY;
  stimulus.dots{2} = initDots(myscreen,dots);
  mglClearScreen(0);mglFlush;
  mglClearScreen(0);mglFlush;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initDots(myscreen,dots)

% init dots structure
if ieNotDefined('dots'),dots = [];end

% convert the passed in parameters to real units
if ~isfield(dots,'type'), dots.type = 'Linear';,end
if ~isfield(dots,'rmax'), dots.rmax = max(myscreen.imageWidth,myscreen.imageHeight)/2;,end
if ~isfield(dots,'xcenter'), dots.xcenter = 0;,end
if ~isfield(dots,'ycenter'), dots.ycenter = 0;,end
if ~isfield(dots,'dotsize'), dots.dotsize = 4;,end
if ~isfield(dots,'density'), dots.density = 5;,end
if ~isfield(dots,'coherence'), dots.coherence = 1;,end
if ~isfield(dots,'speed'), dots.speed = 6;,end
if ~isfield(dots,'dir'), dots.dir = 0;,end

% init the dots
dots = feval(sprintf('initDots%s',dots.type),dots,myscreen);

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
% step dots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updateDots(dots,coherence,myscreen)

% update the dots
dots = feval(sprintf('updateDots%s',dots.type),dots,coherence,myscreen);

% make circular aperture
x = dots.x;
y = dots.y;
good = find(sqrt(x.^2 + y.^2) <= dots.rmax);
x = x(good);
y = y(good);
% draw the dots
%mglPoints2(dots.x+dots.xcenter,dots.y+dots.ycenter,dots.dotsize,[1 1 1]);
mglPoints2(x+dots.xcenter,y+dots.ycenter,dots.dotsize,[1 1 1]);

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
%dots.x(~dots.coherent) = rand(1,sum(~dots.coherent))*dots.width;
%dots.y(~dots.coherent) = rand(1,sum(~dots.coherent))*dots.height;

% make sure we haven't gone off the patch
% do the dots separately for left and right hand side
dots.x(dots.x < dots.xmin) = dots.x(dots.x < dots.xmin)+dots.width;
dots.x(dots.x > dots.xmax) = dots.x(dots.x > dots.xmax)-dots.width;
dots.y(dots.y < dots.ymin) = dots.y(dots.y < dots.ymin)+dots.height;
dots.y(dots.y > dots.ymax) = dots.y(dots.y > dots.ymax)-dots.height;
