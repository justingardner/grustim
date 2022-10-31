%    dianaTestTriggers.m
%
%       
%       usage: dianaTestTriggers
%          by: jw + jlg
%        date: October 2022
%       
%       Test stimulus for DIANA imaging.
%      
function myscreen = dianaTestTriggers(varargin)

% check arguments
getArgs(varargin);

% init the stimulus
global stimulus;

% set stimulus vaules, these are in ms, but will be
% translated into number of video frames and the
% actual values that you can expect will be displayed.
stimulus.onsetDelay = 50;
% this is in Hz. The stimulus is phase-reversing checkerboard
stimulus.tf = 2;
% number of cycles to display for
% size in degrees
stimulus.width = 78;
stimulus.height = 64;
% spatial frequency in cycles/deg
stimulus.sf = 0.5;
stimulus.orientation = 90;
stimulus.phase = 0;
stimulus.contrast = 0.5;
% number of cycles to run stimulus for
stimulus.numCycles = 1;
% use fix task (or just draw a fixation cross if not)
stimulus.useFixTask = 1;
stimulus.fixWidth = 1;

% initilaize the screen
%myscreen = initScreen('debug');
myscreen = initScreen('diana');

% clear the screen to gray
mglClearScreen(0.5);mglFlush;

% wait for backtick to start experiment (probably not necessary as we
% wait at the beginning of each segment anyway).
task{1}.waitForBacktick = 0;

% initStimulus
myscreen = initStimulus('stimulus',myscreen);
myInitStimulus(myscreen);

% set task parameters
task{1}.segmin = [0.1 stimulus.stimLen];
task{1}.segmax = [0.1 stimulus.stimLen];
task{1}.getResponse = [0 0];

% set number of trials to infinite (to stop stimulus hit ESC)
task{1}.numTrials = inf;

% synch to vol (note this will sync *after* the first segment is done in time, so that
% the stimulus presented in the 2nd segment will start just after the volume acquisition
% trigger comes.
task{1}.synchToVol = [1 0];

% setup fixation task
if stimulus.useFixTask
  global fixStimulus
  fixStimulus.fixWidth = stimulus.fixWidth
  [fixTask myscreen] = fixStairInitTask(myscreen);
end

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
  if stimulus.useFixTask
    % draw fixation task
    [fixTask myscreen] = updateTask(fixTask,myscreen,1);
  end
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
  % calculate how long it took to run last segment
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
  % flashing box from white to black

  %if isodd(floor((myscreen.tick-stimulus.startTick)/stimulus.framesPerHalfCycle))
  if floor((myscreen.tick-stimulus.startTick)/stimulus.framesPerHalfCycle) == 0 % edit here - set 0 instead of odd; otherwise flashes back
      
    mglBltTexture(stimulus.plaid1,[0 0]);
  else
    mglBltTexture(stimulus.plaid2,[0 0]);
  end
  % draw fixation cross, if not using task
  if ~stimulus.useFixTask
    mglFixationCross(stimulus.fixWidth,2,[0 1 1]);
  end
else
  % otherwise just clear the screen to black
  mglClearScreen(0.5);
  % draw fixation cross, if not using task
  if ~stimulus.useFixTask
    mglFixationCross(stimulus.fixWidth,2,[0 1 1]);
  end
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

% calculate number of frames for onestDelay
stimulus.onsetDelayFrames = round(stimulus.onsetDelay*myscreen.framesPerSecond/1000);
stimulus.onsetDelayActual = 1000*stimulus.onsetDelayFrames*(1/myscreen.framesPerSecond);
disp(sprintf('(dianaTestTriggers) onsetDelay for %i frames: %0.1f ms (desired=%0.1f ms)',stimulus.onsetDelayFrames,stimulus.onsetDelayActual,stimulus.onsetDelay));

% gratings
stimulus.g1 = stimulus.contrast*mglMakeGrating(stimulus.width,stimulus.height,stimulus.sf,stimulus.orientation,stimulus.phase);
stimulus.g2 = stimulus.contrast*mglMakeGrating(stimulus.width,stimulus.height,stimulus.sf,stimulus.orientation+90,stimulus.phase);

% plaid
stimulus.plaid = stimulus.g1 + stimulus.g2;

% make into checkerboard
stimulus.plaid(stimulus.plaid>0) = 1;
stimulus.plaid(stimulus.plaid<0) = -1;

stimulus.plaid1 = mglCreateTexture(stimulus.plaid);
stimulus.plaid2 = mglCreateTexture(-stimulus.plaid);

% compute frames per each cycle
stimulus.framesPerHalfCycle = round(((1/stimulus.tf) * myscreen.framesPerSecond)/2);
stimulus.cycleLenActual = 1000*(stimulus.framesPerHalfCycle*2)*(1/myscreen.framesPerSecond);
stimulus.tfActual = 1/(stimulus.cycleLenActual/1000);
disp(sprintf('(dianaTestTriggers) tf=%0.1f (desired: %0.1f) cycleLen: %i frames (%0.1f ms, desired: %0.1f ms)',stimulus.tfActual,stimulus.tf,stimulus.framesPerHalfCycle*2,stimulus.cycleLenActual,1000/stimulus.tf));

% compute length of time stimulus will be on for
stimulus.stimLen = stimulus.cycleLenActual*stimulus.numCycles/1000;
disp(sprintf('(dianaTestTriggers) Running %i cycles for %i frames (%0.1f ms)',stimulus.numCycles,stimulus.stimLen*myscreen.framesPerSecond,stimulus.stimLen));

% Display total length of stimulus
disp(sprintf('(dianaTestTriggers) Total stimulus time (including onset delay): %0.1f',stimulus.onsetDelayActual+stimulus.cycleLenActual*stimulus.numCycles));

% set initial startTime
stimulus.startTime = nan;