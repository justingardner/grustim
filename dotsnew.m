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
function myscreen = dotsnew

% check arguments
if ~any(nargin == [0])
  help dotsnew
  return
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
stimulus.dots.stencilWidth = 12;
stimulus.dots.stencilHeight = 12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up baseline task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
triallen = 1;

% set the first task to be the fixation staircase task
[task{1} myscreen] = fixStairInitTask(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
task{2}{1}.numBlocks = 1;
task{2}{1}.parameter.dir = [0;0];
task{2}{1}.parameter.coherence = [0;0];
task{2}{1}.random = 0;
task{2}{1}.seglen = 10;
task{2}{1}.waitForBacktick = 1;
task{2}{1}.timeInVols = 1;

task{2}{2}.parameter.dir = [15 135 255;15 135 255];
task{2}{2}.parameter.coherence = [0.5 1;0.5 1];
%task{2}{2}.parameter.coherence = [1;1];
task{2}{2}.random = 1;
% set up segments of trials
task{2}{2}.segmin = [2 3];
task{2}{2}.segmax = [2 7];
task{2}{2}.waitForBacktick = 0;
task{2}{2}.timeInVols = 1;
% set up which traces to write out
task{2}{2}.writeTrace{1}.tracenum = [1 2 3 4];
task{2}{2}.writeTrace{1}.tracevar{1} = 'dir';
task{2}{2}.writeTrace{1}.tracevar{2} = 'dir';
task{2}{2}.writeTrace{1}.tracevar{3} = 'coherence';
task{2}{2}.writeTrace{1}.tracevar{4} = 'coherence';
task{2}{2}.writeTrace{1}.tracerow = [1 2 1 2];
task{2}{2}.writeTrace{1}.usenum = [1 1 1 1];

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

% if we are using dots that appear on one half or the other screen
stimulus.dots.lefthalf = stimulus.dots.x < stimulus.dots.xhalf;
stimulus.dots.righthalf = stimulus.dots.x >= stimulus.dots.xhalf;

% get number of dots for each
stimulus.dots.leftN = sum(stimulus.dots.lefthalf);
stimulus.dots.rightN = sum(stimulus.dots.righthalf);

% create stencil
mglStencilCreateBegin(1);
% get position of first cutout
xpos = stimulus.dots.stencilPos;
ypos = 0;
% and size of the oval
stencilSize(1) = stimulus.dots.stencilWidth;
stencilSize(2) = stimulus.dots.stencilHeight;
% and draw that oval
mglFillOval(xpos,ypos,stencilSize);
% now shift over the oval to the otherside and draw it there too
xpos = -stimulus.dots.stencilPos;
mglFillOval(xpos,ypos,stencilSize);
mglStencilCreateEnd;
mglClearScreen;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots randomwalk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = update_randomwalk_linear2(stimulus,myscreen)

% pick a random set of dots
stimulus.dots.coherent(stimulus.dots.lefthalf) = rand(1,stimulus.dots.leftN) < stimulus.dots.coherence(1);
stimulus.dots.coherent(stimulus.dots.righthalf) = rand(1,stimulus.dots.rightN) < stimulus.dots.coherence(2);

% now move those dots in the right direction
stimulus.dots.x(stimulus.dots.coherent&stimulus.dots.lefthalf) = stimulus.dots.x(stimulus.dots.coherent&stimulus.dots.lefthalf)+stimulus.dots.xstep(1);
stimulus.dots.y(stimulus.dots.coherent&stimulus.dots.lefthalf) = stimulus.dots.y(stimulus.dots.coherent&stimulus.dots.lefthalf)+stimulus.dots.ystep(1);
stimulus.dots.x(stimulus.dots.coherent&stimulus.dots.righthalf) = stimulus.dots.x(stimulus.dots.coherent&stimulus.dots.righthalf)+stimulus.dots.xstep(2);
stimulus.dots.y(stimulus.dots.coherent&stimulus.dots.righthalf) = stimulus.dots.y(stimulus.dots.coherent&stimulus.dots.righthalf)+stimulus.dots.ystep(2);

% randomwalk rule
%thisdir = rand(1,sum(~stimulus.dots.coherent))*2*pi;
%stimulus.dots.x(~stimulus.dots.coherent) = stimulus.dots.x(~stimulus.dots.coherent)+cos(thisdir)*stimulus.dots.stepsize;
%stimulus.dots.y(~stimulus.dots.coherent) = stimulus.dots.y(~stimulus.dots.coherent)+sin(thisdir)*stimulus.dots.stepsize;

% movshon noise
stimulus.dots.x(~stimulus.dots.coherent & stimulus.dots.lefthalf) = stimulus.dots.xmin+rand(1,sum(~stimulus.dots.coherent & stimulus.dots.lefthalf))*stimulus.dots.width/2;
stimulus.dots.x(~stimulus.dots.coherent & stimulus.dots.righthalf) = stimulus.dots.xmin+rand(1,sum(~stimulus.dots.coherent & stimulus.dots.righthalf))*stimulus.dots.width/2+stimulus.dots.width/2;
stimulus.dots.y(~stimulus.dots.coherent) = rand(1,sum(~stimulus.dots.coherent))*stimulus.dots.height;

% make sure we haven't gone off the patch
% do the dots separately for left and right hand side
stimulus.dots.x(stimulus.dots.lefthalf & (stimulus.dots.x < stimulus.dots.xmin)) = stimulus.dots.x(stimulus.dots.lefthalf & (stimulus.dots.x < stimulus.dots.xmin))+stimulus.dots.width;
stimulus.dots.x(stimulus.dots.lefthalf & (stimulus.dots.x > stimulus.dots.xhalf)) = stimulus.dots.x(stimulus.dots.lefthalf & (stimulus.dots.x > stimulus.dots.xhalf))-stimulus.dots.width/2;
stimulus.dots.x(stimulus.dots.righthalf & (stimulus.dots.x < stimulus.dots.xhalf)) = stimulus.dots.x(stimulus.dots.righthalf & (stimulus.dots.x < stimulus.dots.xhalf))+stimulus.dots.width/2;
stimulus.dots.x(stimulus.dots.righthalf & (stimulus.dots.x > stimulus.dots.xmax)) = stimulus.dots.x(stimulus.dots.righthalf & (stimulus.dots.x > stimulus.dots.xmax))-stimulus.dots.width;
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
stimulus.dots.dir = d2r([task.thistrial.dir]);

% get the step size
stimulus.dots.xstep = cos(stimulus.dots.dir)*stimulus.dots.stepsize;
stimulus.dots.ystep = sin(stimulus.dots.dir)*stimulus.dots.stepsize;

stimulus.dots.coherence = task.thistrial.coherence;

if task.thistrial.thisseg == 2
  stimulus.dots.coherence = [0 0];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to display stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = trialStimulusCallback(task,myscreen)

global stimulus;
mglClearScreen;

% now update the dots, by calling update function
stimulus = feval(stimulus.updateFunction,stimulus,myscreen);


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

