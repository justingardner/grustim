% morrone.m
%
%      usage: morrone()
%         by: justin gardner
%       date: 12/15/06
%    purpose: 
%
function myscreen = morrone(type)


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
yoffset = -3.5;
% stimuli are different in size for localizer
if type == 0
  stimulus.dots.windowSize = [4.5 4.5];
  stimulus.dots.lifetime = 1;
  stimulus.dots.speed = 15;
  stimulus.dots.n = 50;
  stimulus.dots.windowPosition = [0 yoffset];
  stimulus.dots.direction = 90;
  stimulus.dots.size = 9/60;
  % and for main experiment
else
  stimulus.dots.windowSize = [1 8];
  stimulus.dots.lifetime = 2;
  stimulus.dots.speed = 15;
  stimulus.dots.n = 48;
  stimulus.dots.windowPosition = [0 yoffset];
  stimulus.dots.direction = 90;
  stimulus.dots.size = 9/60;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up segments of trials

% can't actually fit the same dimensions on our display
% so we make smaller by a scale factor
scalefactor =  ((myscreen.imageWidth/2)-1)/15;

%scalefactor = 0.8;
disp(sprintf('Scalefactor=%0.2f 15=%0.2f 10=%0.2f 7.5=%0.2f 5=%0.2f',scalefactor,scalefactor*15,scalefactor*10,scalefactor*7.5,scalefactor*5));
if type == 0
  task{2}{1}.parameter.position = [-7.5 7.5]*scalefactor;
  task{2}{1}.random = 0;
  task{2}{1}.segmin = [9];
  task{2}{1}.segmax = [9];
  task{2}{1}.segquant = [0];
  task{2}{1}.synchToVol = [0 0];
  task{2}{1}.waitForBacktick = 1;
  task{2}{1}.private.fixation = 0;
  task{2}{1}.writeTrace{1}.tracenum = 1;
  task{2}{1}.writeTrace{1}.tracevar{1} = 'position';
  global fixStimulus;
  fixStimulus.pos = [0 yoffset];
else
  task{2}{1}.parameter.position = [-15 -5 5 15]*scalefactor;
  task{2}{1}.random = 1;
  task{2}{1}.segmin = [3 3];
  task{2}{1}.segmax = [3 9];
  task{2}{1}.segquant = [0 1.5];
  task{2}{1}.synchToVol = [0 0];
  task{2}{1}.waitForBacktick = 1;
  task{2}{1}.writeTrace{1}.tracenum = 1;
  task{2}{1}.writeTrace{1}.tracevar{1} = 'position';
  task{2}{1}.writeTrace{1}.usenum = 1;

  % set the fixation position
  global fixStimulus;
  fixationPositions = [-10 0 10]*scalefactor;
  fixStimulus.pos = [fixationPositions(type) yoffset];

  % write the fixation position
  myscreen = writeTrace(type,myscreen.stimtrace+2,myscreen);

end

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
stimulus.dots.x = rand(1,stimulus.dots.n)*stimulus.dots.windowSize(1);
stimulus.dots.y = rand(1,stimulus.dots.n)*stimulus.dots.windowSize(2);

% init dots to have an age of 0 or 1
stimulus.dots.age = round(rand(1,stimulus.dots.n) * (stimulus.dots.lifetime-1));

% based on the speed and direction calculate the x and y step of the dots
stimulus.dots.xStep = cos(stimulus.dots.direction*pi/180)*stimulus.dots.speed/myscreen.framesPerSecond;
stimulus.dots.yStep = sin(stimulus.dots.direction*pi/180)*stimulus.dots.speed/myscreen.framesPerSecond;

stimulus.dots.color = rand(1,stimulus.dots.n) > 0.5;

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
mglGluDisk(task.thistrial.position+stimulus.dots.x(stimulus.dots.color==1)+stimulus.dots.windowPosition(1)-stimulus.dots.windowSize(1)/2,stimulus.dots.y(stimulus.dots.color==1)+stimulus.dots.windowPosition(2)-stimulus.dots.windowSize(2)/2,stimulus.dots.size,1);
% draw black points
mglGluDisk(task.thistrial.position+stimulus.dots.x(stimulus.dots.color==0)+stimulus.dots.windowPosition(1)-stimulus.dots.windowSize(1)/2,stimulus.dots.y(stimulus.dots.color==0)+stimulus.dots.windowPosition(2)-stimulus.dots.windowSize(2)/2,stimulus.dots.size,0);

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

% reposition dots
stimulus.dots.x(dotsToReposition) = rand(1,length(dotsToReposition))*stimulus.dots.windowSize(1);
stimulus.dots.y(dotsToReposition) = rand(1,length(dotsToReposition))*stimulus.dots.windowSize(2);

% move dots
stimulus.dots.x(dotsToMove) = stimulus.dots.x(dotsToMove)+stimulus.dots.xStep;
stimulus.dots.y(dotsToMove) = stimulus.dots.y(dotsToMove)+stimulus.dots.yStep;

% make sure we haven't gone out of window
stimulus.dots.x(stimulus.dots.x > stimulus.dots.windowSize(1)) = stimulus.dots.x(stimulus.dots.x > stimulus.dots.windowSize(1))-stimulus.dots.windowSize(1);
stimulus.dots.y(stimulus.dots.y > stimulus.dots.windowSize(2)) = stimulus.dots.y(stimulus.dots.y > stimulus.dots.windowSize(2))-stimulus.dots.windowSize(2);
stimulus.dots.x(stimulus.dots.x < 0) = stimulus.dots.x(stimulus.dots.x < 0)+stimulus.dots.windowSize(1);
stimulus.dots.y(stimulus.dots.y < 0) = stimulus.dots.y(stimulus.dots.y < 0)+stimulus.dots.windowSize(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to set stimulus parameters at
% the beginning of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task,myscreen)

global stimulus;
if task.thistrial.thisseg == 1
  % randomly set direction up or down
  stimulus.dots.direction = (rand > 0.5)*180+90;
  % save it in the trace
  myscreen = writeTrace(stimulus.dots.direction,myscreen.stimtrace+1,myscreen);
  % and change the dot xstep and ystep
  stimulus.dots.xStep = cos(stimulus.dots.direction*pi/180)*stimulus.dots.speed/myscreen.framesPerSecond;
  stimulus.dots.yStep = sin(stimulus.dots.direction*pi/180)*stimulus.dots.speed/myscreen.framesPerSecond;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to display stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = trialStimulusCallback(task,myscreen)

global stimulus;

mglClearScreen;

% draw and update dots
if task.thistrial.thisseg == 1
  displayDots(task,myscreen);
  updateDots(task,myscreen);
end


