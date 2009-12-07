% dots.m
%
%      usage: [stimulus,task,myscreen]=dots(stimulus)
%         by: justin gardner
%         ed: ds
%       date: 04/15/06
%    purpose: block design moving dots localizer
%
%
%
function myscreen = dotsvsblank(tasktype)

% check arguments
if ~any(nargin == [0 1])
  help dotsnew
  return
end

if ~exist('tasktype') == 1
  tasktype = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other screen parameters
myscreen.autoCloseScreen = 1;
myscreen.allowpause = 0;
myscreen.eatkeys = 0;
myscreen.displayname = 'projector';
myscreen.background = 'black';

myscreen = initScreen(myscreen);

global stimulus
myscreen = initStimulus('stimulus',myscreen);
% optic flow stimuli:
stimulus.dots.type = 'randomwalk_linear2';
stimulus.dots.rmin = 1.5;
stimulus.dots.rmax = 8.5; 
stimulus.dots.dotsize = 2;
stimulus.dots.speed = 5;
stimulus.dots.coherence = [0.6 0.6];
stimulus.dots.n = 600;
stimulus.dots.dir = 0;
stimulus.dots.stencilPos = 8;
stimulus.dots.stencilWidth = myscreen.imageHeight;
stimulus.dots.stencilHeight = myscreen.imageHeight;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up baseline task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
triallen = 1;

% set the first task to be the fixation staircase task
[task{1} myscreen] = fixStairInitTask(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
task{2}{1}.parameter.coherence = [1];
task{2}{1}.waitForBacktick = 1;
task{2}{1}.seglen = 6;
task{2}{1}.private.onsegs = 0;
task{2}{1}.numBlocks = 1;

task{2}{2}.parameter.coherence = [1];
task{2}{2}.random = 0;
% set up segments of trials
task{2}{2}.seglen = ones(1,24);
%task{2}{1}.timeInSecs = 1;
task{2}{2}.private.onsegs = 12;
% set up which traces to write out
task{2}{2}.writeTrace{1}.tracenum = 1;
task{2}{2}.writeTrace{1}.tracevar{1} = 'coherence';

% make event related
if tasktype == 1
  task{2}{2} = rmfield(task{2}{2},'seglen');
  task{2}{2}.segmin = [1 1 5];
  task{2}{2}.segmax = [1 1 11];
  task{2}{2}.synchToVol = [0 0 1];
  task{2}{2}.private.onsegs = 2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialze tasks and stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stimulus = myInitStimulus(stimulus,task,myscreen);

% initialze tasks
for tasknum = 1:length(task{2})
  task{2}{tasknum} = initTask(task{2}{tasknum},myscreen,@startSegmentCallback,@trialStimulusCallback);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the tasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set which phase is active
tnum = 1;

while (tnum <= length(task{2})) && ~myscreen.userHitEsc
  % updatethe task
  [task{2} myscreen tnum] = updateTask(task{2},myscreen,tnum);
  % display fixation cross
  [task{1} myscreen] = updateTask(task{1},myscreen,1);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     STIMULUS Dots start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initDots(stimulus,myscreen,task)

% convert the passed in parameters to real units
if ~isfield(stimulus.dots,'rmin'), stimulus.dots.rmin = 0;,end
if ~isfield(stimulus.dots,'rmax'), stimulus.dots.rmax = 5;,end
if ~isfield(stimulus.dots,'dotsize'), stimulus.dots.dotsize = 3,end

% get the number of dots
if ~isfield(stimulus.dots,'density') 
  if ~isfield(stimulus.dots,'n')
    stimulus.dots.n = round(pi*(stimulus.dots.rmax^2)*stimulus.dots.density);
  else
    disp(sprintf('stimulus.dots.n already set to %i, ignoring stimulus.dots.density',stimulus.dots.n));
  end
end

% set the number of dots if we didn't get anything
if ~isfield(stimulus.dots,'n')
  stimulus.dots.n = 100;
  disp(sprintf('Set stimulus.dots.n=%i',stimulus.dots.n));
end

% get initial position
stimulus.dots.r = stimulus.dots.rmin+rand(1,stimulus.dots.n)*(stimulus.dots.rmax-stimulus.dots.rmin);
stimulus.dots.theta = rand(1,stimulus.dots.n)*2*pi;

% get the step size
stimulus.dots.stepsize = stimulus.dots.speed/myscreen.framesPerSecond;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to update dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = updateDots(stimulus,myscreen)

% now move those dots in the right direction
if (stimulus.dots.motion == 2) || (stimulus.dots.motion == 4)
  stimulus.dots.r = stimulus.dots.r+stimulus.dots.stepsize;
end

% make sure we haven't gone off the patch
stimulus.dots.r(stimulus.dots.r > stimulus.dots.rmax) = stimulus.dots.rmin;
stimulus.dots.r(stimulus.dots.r < stimulus.dots.rmin) = stimulus.dots.rmax;

% get cartesian position
stimulus.dots.x = stimulus.dots.r.*cos(stimulus.dots.theta);
stimulus.dots.y = stimulus.dots.r.*sin(stimulus.dots.theta);

% plot the points
if (stimulus.dots.motion > 0)
  mglPoints2(stimulus.dots.x,stimulus.dots.y,stimulus.dots.dotsize);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     STIMULUS initdots_linear2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init the dots stimuli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initdots_linear2(stimulus,myscreen,task)

% convert the passed in parameters to real units
if ~isfield(stimulus.dots,'rmin'), stimulus.dots.rmin = 0;,end
if ~isfield(stimulus.dots,'rmax'), stimulus.dots.rmax = 5;,end
if ~isfield(stimulus.dots,'dotsize'), stimulus.dots.dotsize = 3;,end
if ~isfield(stimulus.dots,'density'), stimulus.dots.density = 3;,end
if ~isfield(stimulus.dots,'width'), stimulus.dots.width = myscreen.imageWidth;,end
if ~isfield(stimulus.dots,'height'), stimulus.dots.height = myscreen.imageHeight;,end
if ~isfield(stimulus.dots,'coherence'), stimulus.dots.coherence = [1 1];,end

% get the number of dots
stimulus.dots.n = round(stimulus.dots.width*stimulus.dots.height*stimulus.dots.density);

% get max and min points for dots
stimulus.dots.xmin = -stimulus.dots.width/2;
stimulus.dots.xmax = stimulus.dots.width/2;
stimulus.dots.xhalf = 0;
stimulus.dots.ymin = -stimulus.dots.height/2;
stimulus.dots.ymax = stimulus.dots.height/2;
stimulus.dots.yhalf = 0;

% set direction of dots
if (length(stimulus.dots.dir) < 2)
  stimulus.dots.dir = [stimulus.dots.dir stimulus.dots.dir];
else
  stimulus.dots.dir = stimulus.dots.dir;
end
  
% get initial position
stimulus.dots.x = stimulus.dots.xmin+rand(1,stimulus.dots.n)*stimulus.dots.width;
stimulus.dots.y = stimulus.dots.ymin+rand(1,stimulus.dots.n)*stimulus.dots.height;

% get the step size
stimulus.dots.stepsize = stimulus.dots.speed/myscreen.framesPerSecond;
stimulus.dots.xstep = cos(stimulus.dots.dir)*stimulus.dots.stepsize;
stimulus.dots.ystep = sin(stimulus.dots.dir)*stimulus.dots.stepsize;

% create stencil
mglStencilCreateBegin(1);
% get position of first cutout
xpos = 0;
ypos = 0;
% and size of the oval
stencilSize(1) = stimulus.dots.stencilWidth;
stencilSize(2) = stimulus.dots.stencilHeight;
% and draw that oval
mglFillOval(xpos,ypos,stencilSize);
mglStencilCreateEnd;
mglClearScreen;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots randomwalk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = update_randomwalk_linear2(stimulus,myscreen)

% pick a random set of dots
stimulus.dots.coherent = rand(1,stimulus.dots.n) < stimulus.dots.coherence;

% now move those dots in the right direction
stimulus.dots.x(stimulus.dots.coherent) = stimulus.dots.x(stimulus.dots.coherent)+stimulus.dots.xstep;
stimulus.dots.y(stimulus.dots.coherent) = stimulus.dots.y(stimulus.dots.coherent)+stimulus.dots.ystep;

% randomwalk rule
%thisdir = rand(1,sum(~stimulus.dots.coherent))*2*pi;
%stimulus.dots.x(~stimulus.dots.coherent) = stimulus.dots.x(~stimulus.dots.coherent)+cos(thisdir)*stimulus.dots.stepsize;
%stimulus.dots.y(~stimulus.dots.coherent) = stimulus.dots.y(~stimulus.dots.coherent)+sin(thisdir)*stimulus.dots.stepsize;

% movshon noise
stimulus.dots.x(~stimulus.dots.coherent) = stimulus.dots.xmin+rand(1,sum(~stimulus.dots.coherent))*stimulus.dots.width;
stimulus.dots.y(~stimulus.dots.coherent) = rand(1,sum(~stimulus.dots.coherent))*stimulus.dots.height;

% make sure we haven't gone off the patch

stimulus.dots.x((stimulus.dots.x < stimulus.dots.xmin)) = stimulus.dots.x((stimulus.dots.x < stimulus.dots.xmin))+stimulus.dots.width;
stimulus.dots.x((stimulus.dots.x > stimulus.dots.xmax)) = stimulus.dots.x((stimulus.dots.x > stimulus.dots.xmax))-stimulus.dots.width;
stimulus.dots.y(stimulus.dots.y > stimulus.dots.ymax) = stimulus.dots.y(stimulus.dots.y > stimulus.dots.ymax)-stimulus.dots.height;
stimulus.dots.y(stimulus.dots.y < stimulus.dots.ymin) = stimulus.dots.y(stimulus.dots.y < stimulus.dots.ymin)+stimulus.dots.height;

% draw the dots
mglStencilSelect(1);
mglPoints2(stimulus.dots.x,stimulus.dots.y,stimulus.dots.dotsize,[1 1 1]);
mglStencilSelect(0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to set stimulus parameters at
% the beginning of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task,myscreen)

global stimulus;
% set the stimulus parameters
% set direction of dots
stimulus.dots.dir = rand(1,1)*2*pi;

% get the step size
stimulus.dots.xstep = cos(stimulus.dots.dir)*stimulus.dots.stepsize;
stimulus.dots.ystep = sin(stimulus.dots.dir)*stimulus.dots.stepsize;

% set the coherence
stimulus.dots.coherence = task.thistrial.coherence;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to display stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = trialStimulusCallback(task,myscreen)

global stimulus;
mglClearScreen;

if task.thistrial.thisseg <= task.private.onsegs
  % now update the dots, by calling update function
  stimulus = feval(stimulus.updateFunction,stimulus,myscreen);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stimulus task myscreen] = myInitStimulus(stimulus,task,myscreen)

% init dots of different flavors
if isfield(stimulus,'dots') 
  if isfield(stimulus.dots,'type') && strcmp(stimulus.dots.type, 'randomwalk_linear2')
    disp('using _random walk_ linear 2 dots...')
    stimulus = initdots_linear2(stimulus, myscreen, task); 
    stimulus.updateFunction = @update_randomwalk_linear2;
  else
    stimulus = initDots(stimulus,myscreen,task); % normal dots
    stimulus.updateFunction = @updateDots;
  end
end

