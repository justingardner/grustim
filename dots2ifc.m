% dots2ifc - two interval forced choice version of dot task
%
%      usage: [  ] = dots2ifc( stimtype )
%         by: denis schluppeck
%       date: 2007-06-15
%        $Id:
%     inputs: stimtype
%    outputs: 
%
%    purpose: a two interval forced choice version of a dot task, as
%    described in the grant submitted to the MRC by DS (May
%    2007). comparision between stimulus in interval 1 and two,
%    separated by long and variable delays... trials separated by long
%    and variable intertrial intervals.
%
%    based on dotsvblank by justin g 

function [  ]=dots2ifc( stimtype )

% check arguments
if ~any(nargin == [0 1])
  help dots2ifc
  return
end

if ~exist('stimtype', 'var')
  % random walk dots
  stimtype = 0;
  % stimtype = 1: % movshon noise
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other screen parameters
myscreen.autoCloseScreen = 0;
myscreen.saveData = 1;
myscreen.allowpause = 0;
myscreen.eatkeys = 1;
myscreen.displayname = 'projector';
myscreen.background = 'black';
myscreen.TR = 1.5; %TR = 1500ms

myscreen = initScreen(myscreen);
myscreen.keyboard.backtick = mglCharToKeycode({'5'}); %TR = 1500ms
myscreen.keyboard.nums = mglCharToKeycode({'1' '2' '3' '4'    '6' '7' '8' '9' '0'});

global stimulus
myscreen = initStimulus('stimulus',myscreen);
% linear dots
stimulus.dots.type = 'randomwalk_linear2';
stimulus.dots.movshonNoise = stimtype;
stimulus.dots.rmin = 1.5;
stimulus.dots.rmax = 6.5; 
stimulus.dots.dotsize = 4;
stimulus.dots.speed = 8;
stimulus.dots.density = 8;
stimulus.dots.dir = 0;
stimulus.dots.stencilPos = 8;
stimulus.dots.stencilWidth = myscreen.imageHeight;
stimulus.dots.stencilHeight = myscreen.imageHeight;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up baseline task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
triallen = 1;

global fixStimulus
fixStimulus.responseTime = 2;
fixStimulus.stimTime = 0.6;
fixStimulus.interTime = 0.5; 

% set the first task to be the fixation staircase task
[task{1} myscreen] = fixStairInitTask(myscreen); % change this///

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
task{2}{1}.parameter.coherence = [1];
task{2}{1}.waitForBacktick = 1;
task{2}{1}.seglen = 4*myscreen.TR;
task{2}{1}.private.onsegs = 0;
task{2}{1}.numBlocks = 1;

task{2}{2}.parameter.coherence = [1/2 2/3 1];
task{2}{2}.random = 1;
% block design
% set up segments of trials
% task{2}{2}.seglen = ones(1,24);
% task{2}{1}.timeInSecs = 1;
% task{2}{2}.private.onsegs = 12;
% set up which traces to write out
% task{2}{2}.writeTrace{1}.tracenum = 1;
% task{2}{2}.writeTrace{1}.tracevar{1} = 'coherence';

% make event related
task{2}{2}.segmin = [1 1 2 1 3]*myscreen.TR;
task{2}{2}.segmax = [1 1 4 1 6]*myscreen.TR;
task{2}{2}.segquant = [0 0 1 0 1]*myscreen.TR;
task{2}{2}.private.onsegs = [2 4]; % are dots on the screen or not?

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
  stimulus.dots.n = 100;
  disp(sprintf('Set stimulus.dots.n=%i',stimulus.dots.n));
else
   stimulus.dots.n = round(pi*(stimulus.dots.rmax^2)*stimulus.dots.density);  
   disp(sprintf('stimulus.dots.n %i, from density: %i per deg^2',stimulus.dots.n,stimulus.dots.density));
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

% area of annulus is pi*(r_o^2-r_i^2)

stimulus.dots.n = round(    [stimulus.dots.rmax^2-stimulus.dots.rmin^2].*stimulus.dots.density  );
disp(sprintf('stimulus.dots.n %i, from density: %i per deg^2',stimulus.dots.n,stimulus.dots.density));

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

if stimulus.dots.movshonNoise
  % movshon noise
  stimulus.dots.x(~stimulus.dots.coherent) = stimulus.dots.xmin+rand(1,sum(~stimulus.dots.coherent))*stimulus.dots.width;
  stimulus.dots.y(~stimulus.dots.coherent) = rand(1,sum(~stimulus.dots.coherent))*stimulus.dots.height;
else
  % randomwalk rule
  thisdir = rand(1,sum(~stimulus.dots.coherent))*2*pi;
  stimulus.dots.x(~stimulus.dots.coherent) = stimulus.dots.x(~stimulus.dots.coherent)+cos(thisdir)*stimulus.dots.stepsize;
  stimulus.dots.y(~stimulus.dots.coherent) = stimulus.dots.y(~stimulus.dots.coherent)+sin(thisdir)*stimulus.dots.stepsize;
end

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

% set the coherence (0 for segments that are NOT ON)
if any(task.thistrial.thisseg == task.private.onsegs)
  stimulus.dots.coherence = task.thistrial.coherence;
else
  stimulus.dots.coherence = 0;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to display stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = trialStimulusCallback(task,myscreen)

global stimulus;
mglClearScreen;

if 1 %any(task.thistrial.thisseg == task.private.onsegs) %% changed to vector [ which segment is ON??]
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

