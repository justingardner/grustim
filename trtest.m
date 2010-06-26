% trtest
%
%        $Id: taskTemplate.m 217 2007-04-04 16:54:42Z justin $
%      usage: trtest
%         by: justin gardner
%       date: 06/26/2010
%  copyright: (c) 2010 Justin Gardner (GPL see mgl/COPYING)
%    purpose: show tr length
%
function myscreen = trtest

% check arguments
if ~any(nargin == [0])
  help taskTemplateSaccade
  return
end

% initialize the screen
myscreen.background = 0;
myscreen = initScreen(myscreen);

% setup task
task{1}.waitForBacktick = 1;
task{1}.seglen = 1;


% initialize the task
for phaseNum = 1:length(task)
  [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@stimStartSegmentCallback,@stimDrawStimulusCallback);
end

mglTextDraw('Hit ESC to quit',[5 4.5],1);
mglTextDraw('Waiting for scanner to start',[5 6],1);
mglFlush;
myscreen.flushMode = -1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task)) && ~myscreen.userHitEsc
  % update the task
  [task myscreen phaseNum] = updateTask(task,myscreen,phaseNum);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = stimStartSegmentCallback(task, myscreen)

mglClearScreen;

% get volume times
volTimes = myscreen.events.time(find(myscreen.events.tracenum == 1));

% display time scanner has been running for and volume number
mglTextDraw('Hit ESC to end',[5 4.5],1);
mglTextDraw(sprintf('Scan time: %0.0fs',mglGetSecs(volTimes(1))),[5 3],1);
mglTextDraw(sprintf('Volume number: %i',myscreen.volnum),[5 1.5],1);

% if we have more than one volume than display average framePeriod
if length(volTimes) > 1
  framePeriod = diff(volTimes);
  mglTextDraw(sprintf('Mean frame period: %0.4f',mean(framePeriod)),[5 0],1);
  mglTextDraw(sprintf('Min frame period: %0.4f',min(framePeriod)),[5 -1.5],1);
  mglTextDraw(sprintf('Max frame period: %0.4f',max(framePeriod)),[5 -3],1);
  mglTextDraw(sprintf('Last frame period: %0.4f',framePeriod(end)),[5 -4.5],1);
end

% don't redraw screen
myscreen.flushMode = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = stimDrawStimulusCallback(task, myscreen)

