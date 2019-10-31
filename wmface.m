% wmface
%
%        $Id$
%      usage: wmface
%         by: justin gardner
%       date: 10/15/2019
%    purpose: working memory task based on Harrison & Tong with camera capture for Druckmann lab
%
function e = wmface(scan)

% check arguments
if ~any(nargin == [0 1])
  help cameraTest
  return
end

% default to non-scan mode
if nargin < 1,scan = 0;end

debugMode = 0;

% stimulus parameters
global stimulus;
stimulus.scan = scan;
if stimulus.scan
    waitForBacktick = 1;
else
    waitForBacktick = 0;
end
if scan
    stimulus.contrast = 0.175;
else
    stimulus.contrast = 0.1;
end
stimulus.innerWidth = 1.5;
stimulus.outerWidth = 10;
stimulus.sf = 1;
% set to false for sharp edged stimulus
stimulus.gabor = 0;
% fixation width
stimulus.fixWidth = 1.5;

% orientations to display
stimulus.orientations = [25 115];
stimulus.orientationJitter = 3;
stimulus.constantVals = [3 6];

% delay interval in seconds
stimulus.delayInterval = 11;


if debugMode
  stimulus.delayInterval = 5;
%   stimulus.constantVals = [20 30];
%   stimulus.constantVals = 25;
end
% initalize the screen
myscreen.background = 128/255;
myscreen = initScreen(myscreen);

% intialize the camera
dispHeader('(cameraTest) Starting Camera Thread');
mglCameraThread('init');
stimulus.cameraImages = {};
dispHeader;
% place to save data
% FIX, this should endup being put in myscreen.datadir;
stimulus.cameraDataDir = '~/Desktop/camera';
if ~isdir(stimulus.cameraDataDir),mkdir(stimulus.cameraDataDir);end
if isempty(myscreen.SID)
  stimulus.cameraFileStem = 'TEST';
else
  stimulus.cameraFileStem = sprintf('%s',myscreen.SID);
end
stimulus.cameraFileStem = sprintf('[%s_%s_%s]',stimulus.cameraFileStem,datestr(now,'yyyymmdd'),datestr(now,'hhmmss'));

% by waiting for the backtick key to be pressed before starting the experiment
% (for systems that use NI digital I/O, this will wait for the digital
% signal that the scanner has started collecting data)
task{1}.waitForBacktick = 1;
task{1}.segmin = [0.2 0.4 0.2 0.4 0.8 stimulus.delayInterval 0.5 2.5 2 1];
task{1}.segmax = [0.2 0.4 0.2 0.4 0.8 stimulus.delayInterval 0.5 2.5 2 1];
task{1}.getResponse = [0 0 0 0 0 0 0 1 0 0];
task{1}.synchToVol = [0 0 0 0 0 0 0 waitForBacktick];
% parameters
task{1}.randVars.uniform.cue = [1 2];
task{1}.randVars.uniform.orientationOrder = [1 2];
task{1}.randVars.uniform.clockwiseCounterclockwise = [-1 1];
task{1}.randVars.calculated.orientationJitter = [nan nan];
task{1}.randVars.calculated.orientationJitter = [nan nan];
task{1}.randVars.calculated.orientation = [nan nan];
task{1}.randVars.calculated.orientationThreshold = nan;
task{1}.random = 1;
task{1}.numTrials = 3;

% initialize the task
for phaseNum = 1:length(task)
  [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end

% init the stimulus
myscreen = initStimulus('stimulus',myscreen);

% to initialize the stimulus for your experiment.
stimulus = myInitStimulus(stimulus,myscreen);

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
  [task myscreen phaseNum] = updateTask(task,myscreen,phaseNum);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

% quit camera
dispHeader('(cameraTest) Ending Camera Thread');
mglCameraThread('quit');
dispHeader;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;
% set fixation white
if task.thistrial.thisseg == 1
  % start saving camera
  mglCameraThread('capture','timeToCapture',sum(task.segmax(1:(end-2))));
  
  stimulus.fixColor = [1 1 1];

  % set the orientation jitter of both
  task.thistrial.orientationJitter(1) = round(stimulus.orientationJitter*(2*(rand - 0.5)));
  task.thistrial.orientationJitter(2) = round(stimulus.orientationJitter*(2*(rand - 0.5)));

  % get the threshold value to test
  [task.thistrial.orientationThreshold stimulus.s] = doStaircase('testValue',stimulus.s);

  % set first and second orientation
  task.thistrial.orientation(1) = stimulus.orientations(task.thistrial.orientationOrder);
  task.thistrial.orientation(2) = stimulus.orientations(setdiff([1 2],task.thistrial.orientationOrder));

  % keep saved orientations to display
  stimulus.displayOrientation(1) = task.thistrial.orientation(1) + task.thistrial.orientationJitter(1);
  stimulus.displayOrientation(2) = task.thistrial.orientation(2) + task.thistrial.orientationJitter(2);
  stimulus.matchOrientation = task.thistrial.orientation(task.thistrial.cue)+task.thistrial.orientationJitter(task.thistrial.cue)+task.thistrial.orientationThreshold*task.thistrial.clockwiseCounterclockwise;

  % display what we have selected
  disp(sprintf('Trial %i: Orientation 1: %0.1f + %0.1f: %f',task.trialnum,task.thistrial.orientation(1),task.thistrial.orientationJitter(1),stimulus.displayOrientation(1)));
  disp(sprintf('Trial %i: Orientation 2: %0.1f + %0.1f: %f',task.trialnum,task.thistrial.orientation(2),task.thistrial.orientationJitter(2),stimulus.displayOrientation(2)));
  disp(sprintf('Trial %i: Match to %i (%0.1f): %0.1f',task.trialnum,task.thistrial.cue,task.thistrial.orientationThreshold*task.thistrial.clockwiseCounterclockwise,stimulus.matchOrientation));
elseif task.thistrial.thisseg == 9
  stimulus.fixColor = [1 1 1];
elseif task.thistrial.thisseg == 10
  % save camera images
%  stimulus.cameraImages{end+1} = mglCameraThread('get');
  % create filename for images
  saveCameraTime = mglGetSecs;
  filename = fullfile(stimulus.cameraDataDir,sprintf('%s-[%04i]',stimulus.cameraFileStem,task.trialnum));
  stimulus.cameraImages{end+1} = mglCameraThread('save','videoFilename',filename);
  disp(sprintf('(wmface) Save of camera file took %0.2fs',mglGetSecs(saveCameraTime)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus

mglClearScreen(0.5);
mglClearScreen;
if task.thistrial.thisseg == 1
  % put up grating
  mglBltTexture(stimulus.tex,[0 0],0,0,stimulus.displayOrientation(1));
elseif task.thistrial.thisseg == 3
  % put up grating
  mglBltTexture(stimulus.tex,[0 0],0,0,stimulus.displayOrientation(2));
elseif task.thistrial.thisseg == 7
  % put up match either clockwise or counterclockwise 
  mglBltTexture(stimulus.tex,[0 0],0,0,stimulus.matchOrientation);
elseif (task.thistrial.thisseg == 5)
  % put up text
  mglBltTexture(stimulus.text(task.thistrial.cue),[0 1.5],0,0,0);
end  

% put up fixation cross
%mglGluDisk(0,0,stimulus.innerWidth,0.5,48)
mglFixationCross(stimulus.fixWidth,1,stimulus.fixColor);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)

global stimulus

butts = [2 0 1];
corrButt = butts(task.thistrial.clockwiseCounterclockwise+2);
% here, we just check whether this is the first time we got a response
% this trial and display what the subject's response was and the reaction time
if task.thistrial.gotResponse < 1
  if task.thistrial.whichButton == corrButt
    disp(sprintf('Trial %i: Correct (reaction time: %0.2fs)',task.trialnum,task.thistrial.reactionTime));
    stimulus.fixColor = [0 1 0];
    stimulus.s = doStaircase('update',stimulus.s,1);
  else
    disp(sprintf('Trial %i: Incorrect (reaction time: %0.2fs)',task.trialnum,task.thistrial.reactionTime));
    stimulus.fixColor = [1 0 0];
    stimulus.s = doStaircase('update',stimulus.s,0);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen)

% compute the grating
grating = mglMakeGrating(stimulus.outerWidth,stimulus.outerWidth,stimulus.sf,0,0);

[gaussian xMesh yMesh] = mglMakeGaussian(stimulus.outerWidth,stimulus.outerWidth,stimulus.outerWidth/8,stimulus.outerWidth/8); 
% calculate radius of each point
r = sqrt(xMesh.^2 + yMesh.^2);
% now calculate gaussian centered on width of stimulus
gaussianRing = exp(-((r-(stimulus.outerWidth/4)).^2)/(stimulus.outerWidth/10)^2);
gaussianRing(gaussianRing < 0.01) = 0;
% hard edge
if ~stimulus.gabor
    gaussianRing(gaussianRing>0) = 1;
    gaussianRing(r<stimulus.innerWidth/2) = 0;
end
gabor = round(255*(stimulus.contrast*grating.*gaussianRing+1)/2);
stimulus.tex = mglCreateTexture(gabor);

% make text
mglTextSet('Helvetica',32,[1 1 1],0,0,0);
stimulus.text(1) = mglText('1');
stimulus.text(2) = mglText('2');

% start staircase
stimulus.s = doStaircase('init','fixed','fixedVals',stimulus.constantVals);