% morrone.m
%
%      usage: morrone()
%         by: justin gardner
%       date: 12/15/06
%    purpose: 
%
function myscreen = morroneloc(type)


% check arguments
if ~any(nargin == [0 1])
  help morrone
  return
end

if ~exist('type','var'), type = 0; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other screen parameters
myscreen.autoCloseScreen = 1;
myscreen.displayname = 'projector';
myscreen.background = 'gray';
myscreen.allowpause = 1;

myscreen = initScreen(myscreen);

global MGL;
clear global stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);
yoffset = 0;
% stimuli are different in size for localizer
stimulus.dots.windowSize = [myscreen.imageWidth myscreen.imageHeight];
stimulus.dots.lifetime = 18;
stimulus.dots.speed = 15;
stimulus.dots.n = 150;
stimulus.dots.size = 9/60;
stimulus.dots.v = 7;
stimulus.dots.phi = pi/4;
% make radius out to screen corner
stimulus.dots.radius = sqrt(((myscreen.imageWidth/2)^2) + ((myscreen.imageHeight/2)^2));
stimulus.dots.coherence = 0;
stimulus.dots.radial = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up segments of trials

task{2}{1}.parameter.phi = [0 nan -pi/2 nan pi/2 nan -3*pi/2 nan];
task{2}{1}.random = 0;
task{2}{1}.seglen = [1 1 1 1 1 1];
task{2}{1}.timeInVols = 1;
task{2}{1}.waitForBacktick = 1;
task{2}{1}.private.fixation = 0;
global fixStimulus;
fixStimulus.pos = [0 yoffset];

% set the first task to be the fixation staircase task
[task{1} myscreen] = fixStairInitTask(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialze tasks and stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initDots(task{2},myscreen);

% initialze tasks
task{2}{1} = initTask(task{2}{1},myscreen,@startSegmentCallback,@trialStimulusCallback);

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
  [myscreen task] = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init dots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initDots(task,myscreen)

global stimulus;

% calculate x and y positions
stimulus.dots.r = rand(1,stimulus.dots.n)*stimulus.dots.radius;
stimulus.dots.theta = rand(1,stimulus.dots.n)*2*pi;

% init dots to have an age of 0 or 1
stimulus.dots.age = round(rand(1,stimulus.dots.n) * (stimulus.dots.lifetime-1));

% set color randomly to white or black
stimulus.dots.color = rand(1,stimulus.dots.n) > 0.5;

% calculate x and y
stimulus.dots.x = stimulus.dots.r.*cos(stimulus.dots.theta);
stimulus.dots.y = stimulus.dots.r.*sin(stimulus.dots.theta);

% calculate size in pixels
global MGL;
stimulus.dots.pixsize = stimulus.dots.size*MGL.xDeviceToPixels;
stimulus.dots.pixsize = round(stimulus.dots.pixsize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to display dots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function displayDots(task,myscreen)

global stimulus;

% draw white points
mglGluDisk(stimulus.dots.x(stimulus.dots.color==1),stimulus.dots.y(stimulus.dots.color==1),stimulus.dots.size,1);
% draw black points
mglGluDisk(stimulus.dots.x(stimulus.dots.color==0),stimulus.dots.y(stimulus.dots.color==0),stimulus.dots.size,0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to update dot position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateDots(task,myscreen)

global stimulus;

% add one to the age
stimulus.dots.age = mod(stimulus.dots.age+1,stimulus.dots.lifetime);

% if age is 0 then the dot needs to be repositioned
dotsToReposition = find(stimulus.dots.age == 0);
% otherwise it has to be moved in one direction
dotsToMove = find(stimulus.dots.age ~= 0);

if stimulus.dots.radial
  % reposition dots
  stimulus.dots.r(dotsToReposition) = rand(1,length(dotsToReposition))*stimulus.dots.radius;
  stimulus.dots.theta(dotsToReposition) = rand(1,length(dotsToReposition))*2*pi;

  % move dots according to equation 1 of methods
  stimulus.dots.r(dotsToMove) = stimulus.dots.r(dotsToMove) + stimulus.dots.v*cos(stimulus.dots.phiall(dotsToMove))/myscreen.framesPerSecond;
  stimulus.dots.theta(dotsToMove) = stimulus.dots.theta(dotsToMove) + ((stimulus.dots.v*sin(stimulus.dots.phiall(dotsToMove)))/myscreen.framesPerSecond)./stimulus.dots.r(dotsToMove);

  % make sure we have not gone over radius
  stimulus.dots.r(stimulus.dots.r > stimulus.dots.radius) = stimulus.dots.r(stimulus.dots.r > stimulus.dots.radius)-stimulus.dots.radius;
  stimulus.dots.r(stimulus.dots.r < 0) = stimulus.dots.r(stimulus.dots.r < 0)+stimulus.dots.radius;
  
  % calculate x and y
  stimulus.dots.x = stimulus.dots.r.*cos(stimulus.dots.theta);
  stimulus.dots.y = stimulus.dots.r.*sin(stimulus.dots.theta);
else
  % reposition dots
  stimulus.dots.x(dotsToReposition) = rand(1,length(dotsToReposition))*stimulus.dots.windowSize(1)-stimulus.dots.windowSize(1)/2;
  stimulus.dots.y(dotsToReposition) = rand(1,length(dotsToReposition))*stimulus.dots.windowSize(2)-stimulus.dots.windowSize(2)/2;

  % update position
  stimulus.dots.x(dotsToMove) = stimulus.dots.x(dotsToMove)+(stimulus.dots.v.*cos(stimulus.dots.phiall(dotsToMove)))/myscreen.framesPerSecond;
  stimulus.dots.y(dotsToMove) = stimulus.dots.y(dotsToMove)+(stimulus.dots.v.*sin(stimulus.dots.phiall(dotsToMove)))/myscreen.framesPerSecond;

  % make sure we haven't gone out of window
  stimulus.dots.x(stimulus.dots.x > stimulus.dots.windowSize(1)/2) = stimulus.dots.x(stimulus.dots.x > stimulus.dots.windowSize(1)/2)-stimulus.dots.windowSize(1);
  stimulus.dots.y(stimulus.dots.y > stimulus.dots.windowSize(2)/2) = stimulus.dots.y(stimulus.dots.y > stimulus.dots.windowSize(2)/2)-stimulus.dots.windowSize(2);
  stimulus.dots.x(stimulus.dots.x < -stimulus.dots.windowSize(1)/2) = stimulus.dots.x(stimulus.dots.x < -stimulus.dots.windowSize(1)/2)+stimulus.dots.windowSize(1);
  stimulus.dots.y(stimulus.dots.y < -stimulus.dots.windowSize(2)/2) = stimulus.dots.y(stimulus.dots.y < -stimulus.dots.windowSize(2)/2)+stimulus.dots.windowSize(2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to set stimulus parameters at
% the beginning of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task,myscreen)

global stimulus;
if task.thistrial.thisseg == 1
  if ~isnan(task.thistrial.phi)
    if task.thistrial.phi < 0
      stimulus.dots.phiall = ones(1,stimulus.dots.n)*task.thistrial.phi;
      stimulus.dots.radial = 0;
    else
      stimulus.dots.phiall = ones(1,stimulus.dots.n)*task.thistrial.phi;
      stimulus.dots.radial = 1;
    end
  else
    stimulus.dots.phiall = rand(1,stimulus.dots.n)*2*pi;
  end
else
  % only reverse motion for radial
  if stimulus.dots.radial
    stimulus.dots.phiall = mod(stimulus.dots.phiall+pi,2*pi);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to display stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = trialStimulusCallback(task,myscreen)

global stimulus;

mglClearScreen;

% draw and update dots
displayDots(task,myscreen);
updateDots(task,myscreen);


