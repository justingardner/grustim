% saccadotopy - reimplementation of schluppeck et al (2006) stimuli
%
%      usage: [  ] = saccadotopy(direction)
%         by: denis schluppeck
%       date: 2007-03-19
%        $Id: saccadotopy.m,v 1.3 2007/03/19 18:01:31 ds Exp $:
%     inputs: direction
%    outputs: 
%
%    purpose: 
%
%       task:
%  
%        e.g:  
%             saccadotopy('cw')
%
function [  ]=saccadotopy( direction )

% check arguments
if ~any(nargin == [1])
  help saccadotopy
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
stimulus.fixSize = 0.1;
stimulus.fixColor = [1 1 1];
stimulus.targetSize = 0.1;
stimulus.targetColor = [1 1 0];
stimulus.dots.density = 0.5;
stimulus.dots.iradius = 7; % inner annulus radius
stimulus.dots.oradius = 9; % inner annulus radius
stimulus.dots.showdots = 0; % initialize to 0 

stimulus = initDots(stimulus,myscreen);

% ----------------------------------------------------
% add horizontal and vertical offset by changing the MODELVIEW matrix
% this is to make stimuli display correctly at Nott'm scanner setup
stimulus.xoffset = 0;
stimulus.yoffset = 0;
mglTransform('GL_MODELVIEW','glTranslate',stimulus.xoffset,stimulus.yoffset,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up baseline task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
triallen = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prescan
task{1}{1}.numBlocks = 1;
task{1}{1}.seglen = 0;
task{1}{1}.getResponse = 1;
task{1}{1}.waitForBacktick = 0;

% saccadotopy
% task{1}{2}.numBlocks = 1;
if strcmp(direction, 'cw')
  task{1}{2}.parameter.theta = 0:30:330;
elseif strcmp(direction, 'ccw')
  task{1}{2}.parameter.theta = 360-[0:30:330];
end

task{1}{2}.private.thetaJitter = 10; % cw/ccw angle offset (i.e. +- 5)
task{1}{2}.parameter.r = 8; %deg of visual angle
task{1}{2}.private.rJitter = 1; % radius offset (+- 0.5deg)
task{1}{2}.random = 0;
task{1}{2}.seglen = [.25 3 .25 1];
task{1}{2}.waitForBacktick = 1;
task{1}{2}.numBlocks = 1;

% task{1}{2}.timeInVols = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialze tasks and stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialze task
% pre-scan
task{1}{1} = initTask(task{1}{1},myscreen,@pre_startSegmentCallback,@pre_trialStimulusCallback);
% saccadotopy - phase encoded version
task{1}{2} = initTask(task{1}{2},myscreen,@startSegmentCallback,@trialStimulusCallback,@trialResponseCallback);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the tasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set which phase is active
tnum = 1;

while (tnum <= length(task{1})) && ~myscreen.userHitEsc
  % updatethe task
  [task{1} myscreen tnum] = updateTask(task{1},myscreen,tnum);
  
  % update the dots (if the flag is set in that segment)
  stimulus = updateDots(stimulus, myscreen);

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
% pre-scan callback functions... before experiment 
% starts...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = pre_startSegmentCallback(task,myscreen)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to display stimulus during
% pre-scan period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = pre_trialStimulusCallback(task,myscreen)

global stimulus;
mglClearScreen;

% draw a red marker for now
mglGluDisk(0, 0, stimulus.fixSize ,  stimulus.fixColor, 72, 2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to set stimulus parameters at
% the beginning of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task,myscreen)

global stimulus;
% set the stimulus parameters
% mydisp(sprintf('task.thistrial.thisseg: %i\n', task.thistrial.thisseg));
if task.thistrial.thisseg == 1
  task.thistrial.private.thetaJitter = task.private.thetaJitter.*(rand(1)-0.5);
  task.thistrial.private.rJitter = task.private.rJitter.*(rand(1)-0.5);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to display stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = trialStimulusCallback(task,myscreen)

global stimulus;
mglClearScreen;

% values for targets +- jitter
theta = task.thistrial.theta + task.thistrial.private.thetaJitter;
r = task.thistrial.r + task.thistrial.private.rJitter;
x = sind(theta).*r;
y = cosd(theta).*r;

% do some display logic here:
if task.thistrial.thisseg == 1 
  % fix + flash up target
  mglGluDisk(0, 0, stimulus.fixSize, stimulus.fixColor, 24, 2);
  mglGluDisk(x, y, stimulus.targetSize, stimulus.targetColor, 24, 2);
  stimulus.dots.showdots = 0;

elseif task.thistrial.thisseg == 2
  % fix + distractor dots? - update flashing dots...
  mglGluDisk(0, 0, stimulus.fixSize, stimulus.fixColor, 24, 2);
  stimulus.dots.showdots = 1;
    
elseif task.thistrial.thisseg == 3 
  stimulus.dots.showdots = 0;
  
elseif task.thistrial.thisseg == 4 
  % return to fixation
  mglGluDisk(0, 0, stimulus.fixSize, stimulus.fixColor, 24, 2);
  stimulus.dots.showdots = 0;
end

% now update the dots, by calling update function
% stimulus = feval(stimulus.updateFunction,stimulus,myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to display stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = trialResponseCallback(task,myscreen)

global stimulus;


% distractor stuff

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initDots(stimulus,myscreen)

% convert the passed in parameters to real units
if ~isfield(stimulus,'dots') || ~isfield(stimulus.dots,'rmax'), stimulus.dots.rmax = 15;,end
if ~isfield(stimulus.dots,'xcenter'), stimulus.dots.xcenter = 0;,end
if ~isfield(stimulus.dots,'ycenter'), stimulus.dots.ycenter = 0;,end
if ~isfield(stimulus.dots,'dotsize'), stimulus.dots.dotsize = 4;,end
if ~isfield(stimulus.dots,'density'), stimulus.dots.density = 0.5;,end
if ~isfield(stimulus.dots,'color'), stimulus.dots.color = [1 1 1]; end 
if ~isfield(stimulus.dots,'colorClock'), stimulus.dots.colorClock = 0; end % clock not set
if ~isfield(stimulus.dots,'colorClockLen'), stimulus.dots.colorClockLen = 10; end % used to be 20.

% actually a square patch of dots that get stenciled
% so calculate width and height
stimulus.dots.width = stimulus.dots.rmax*2;
stimulus.dots.height = stimulus.dots.rmax*2;

% get the number of dots
stimulus.dots.n = round(stimulus.dots.width*stimulus.dots.height*stimulus.dots.density);
% their lifetimes
stimulus.dots.lifetime = floor( rand(1, stimulus.dots.n).* stimulus.dots.colorClockLen);

% get max and min points for dots
stimulus.dots.xmin = -stimulus.dots.width/2;
stimulus.dots.xmax = stimulus.dots.width/2;
stimulus.dots.ymin = -stimulus.dots.height/2;
stimulus.dots.ymax = stimulus.dots.height/2;

% get initial position
stimulus.dots.x = rand(1,stimulus.dots.n)*stimulus.dots.width;
stimulus.dots.y = rand(1,stimulus.dots.n)*stimulus.dots.height;

% create stencil
mglClearScreen;

mglStencilCreateBegin(1);
mglClearScreen;
mglGluAnnulus(stimulus.dots.xcenter, stimulus.dots.ycenter, stimulus.dots.iradius, stimulus.dots.oradius, [0 0 0], 60, 2); 
mglStencilCreateEnd;
mglClearScreen;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to update dot positions and draw them to screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = updateDots(stimulus,myscreen)

if stimulus.dots.showdots == 1
  % if this is a segment when dots should be shown (i.e. during the memory interval)
  
  % stimulus dots are refreshed after stimulus.dots.liftime
  stimulus.dots.lifetime = stimulus.dots.lifetime-1; 
  stimulus.dots.replace = (stimulus.dots.lifetime < 1);
  stimulus.dots.nreplace = sum(stimulus.dots.replace);

  % reset liftime for dots where lifetime < 1
  stimulus.dots.lifetime(stimulus.dots.replace) = ones(1, stimulus.dots.nreplace)*stimulus.dots.colorClockLen;
 
  % refresh position of those dots
  stimulus.dots.x(stimulus.dots.replace) = rand(stimulus.dots.nreplace,1)*stimulus.dots.width;
  stimulus.dots.y(stimulus.dots.replace) = rand(stimulus.dots.nreplace,1)*stimulus.dots.height;
     
  % this forces things to be wrapped around
  % make sure we haven't gone off the patch
  % do the dots separately for left and right hand side
  stimulus.dots.x(stimulus.dots.x < stimulus.dots.xmin) = stimulus.dots.x(stimulus.dots.x < stimulus.dots.xmin)+stimulus.dots.width; 
  stimulus.dots.x(stimulus.dots.x > stimulus.dots.xmax) = stimulus.dots.x(stimulus.dots.x > stimulus.dots.xmax)-stimulus.dots.width; 
  stimulus.dots.y(stimulus.dots.y < stimulus.dots.ymin) = stimulus.dots.y(stimulus.dots.y < stimulus.dots.ymin)+stimulus.dots.height; 
  stimulus.dots.y(stimulus.dots.y > stimulus.dots.ymax) = stimulus.dots.y(stimulus.dots.y > stimulus.dots.ymax)-stimulus.dots.height; 
  
  % draw the dots
  mglStencilSelect(1); 
  mglPoints2(stimulus.dots.x,stimulus.dots.y,stimulus.dots.dotsize, stimulus.dots.color); 
  mglStencilSelect(0); 
else 
  % an interval in which we need to do nothing
end
