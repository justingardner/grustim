%    dianaTestTriggers.m
%
%       
%       usage: dianaTestTriggers
%          by: josh wilson
%        date: October 2022
%       
%       Test stimulus for DIANA imaging.
%      
function myscreen = dianaTestTriggers(varargin)

% check arguments
getArgs(varargin);

% initilaize the screen
myscreen = initScreen('debug');

% clear the screen to gray
mglClearScreen(0.5);mglFlush

% wait for backtick to start experiment
task{1}.waitForBacktick = 1;

% init the stimulus
global stimulus;

% initStimulus
myscreen = initStimulus('stimulus',myscreen)
myInitStimulus(myscreen);

% set task parameters
task{1}.segmin = [0.1 0.5];
task{1}.segmax = [0.1 0.5];
task{1}.getResponse = [0 0];

% set number of trials to infinite
task{1}.numTrials = inf;

% synch to vol (note this will sync *after* the first segment is done in time, so that
% the stimulus presented in the 2nd segment will start just after the volume acquisition
% trigger comes.
task{1}.synchToVol = [1 0];

% initialize the task
for phaseNum = 1:length(task)
  [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback);
end

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
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;

% print out what we are doing and clear screen
if task.thistrial.thisseg == 1
  stimulus.endTime = mglGetSecs;
  if ~isnan(stimulus.startTime)
    disp(sprintf('(dianaTestTriggers) Stim segment lasted %0.1f ms: Ready to start trial %i',1000*(stimulus.endTime-stimulus.startTime),task.trialnum));
  else    
    disp(sprintf('(dianaTestTriggers) Ready to start trial %i',task.trialnum));
  end
  mglClearScreen(0.5);
elseif task.thistrial.thisseg == 2
  % keep starting tick
  stimulus.startTick = myscreen.tick;
  stimulus.startTime = mglGetSecs;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%resp%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen, checkboardBltTex)

global stimulus

% for the secound segment show a stimulus
if task.thistrial.thisseg == 2
  if 0
    % flashing box from white to black
    if isodd(myscreen.tick)
      mglQuad(stimulus.quadX,stimulus.quadY,[1;1;1]);
    else
      mglQuad(stimulus.quadX,stimulus.quadY,[0;0;0]);
    end
  else
    % flashing box from white to black
    if isodd(floor((myscreen.tick-stimulus.startTick)/stimulus.framesPerHalfCycle))
      mglBltTexture(stimulus.plaid1);
    else
      mglBltTexture(stimulus.plaid2);
    end
  end
else
  % otherwise just clear the screen to black
  mglClearScreen(0.5);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = getResponseCallback(task,myscreen)

%%%%%%%%%%%%%%%%%%%%%%%%
%    myInitStimulus    %
%%%%%%%%%%%%%%%%%%%%%%%%
function myInitStimulus(myscreen)

global stimulus

% set where the quad will be
stimulus.xPosMin = -5;
stimulus.xPosMax = 5;
stimulus.yPosMin = -5;
stimulus.yPosMax = 5;
% compute location of quad
stimulus.quadX = [stimulus.xPosMin stimulus.xPosMin stimulus.xPosMax stimulus.xPosMax]';
stimulus.quadY = [stimulus.yPosMin stimulus.yPosMax stimulus.yPosMax stimulus.yPosMin]';

% ok, now make a checkerboard stimulus
width = 32;
height = 32;
sf = 0.5;
orientation = 90;
phase = 0;
contrast = 0.5;
stimulus.tf = 2;

% gratings
g1 = contrast*mglMakeGrating(width,height,sf,orientation,phase);
g2 = contrast*mglMakeGrating(width,height,sf,orientation+90,phase);

% plaid
plaid = g1 + g2;

% make into checkerboard
plaid(plaid>0) = 1;
plaid(plaid<0) = -1;

stimulus.plaid1 = mglCreateTexture(plaid);
stimulus.plaid2 = mglCreateTexture(-plaid);

% compute frames per each cycle
stimulus.framesPerHalfCycle = round(((1/stimulus.tf) * myscreen.framesPerSecond)/2);

% set initial startTime
stimulus.startTime = nan;