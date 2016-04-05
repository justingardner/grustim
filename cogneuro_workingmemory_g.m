% cogneuro_workingmemory.m.m
%
%      usage: cogneuro_workingmemory(scan) 
%             scan is a tf value which when set to 1 waits for
%             scanner sync pulses and sets the timing appropriately
%         by: justin gardner
%       date: 03/03/2016
%    purpose: Code for replicating Harrison & Tong, 2009
%       
function myscreen = cogneuro_workingmemory_g(scan)

% check arguments
if ~any(nargin == [0 1])
  help cogneuro_workingmemory.m
  return
end
if nargin < 1, scan = 0; end

% set for debugging - makes the delay period shorter among other things.
debugMode = 0;

% stimulus parameters get stored in a global variable
global stimulus;

% parameters that are different for when scanning or not
stimulus.scan = scan;
if stimulus.scan
  % wait for scanner acq pulse at beginning of experiment
  waitForBacktick = 1;
  % contrast of the stimulus
  stimulus.contrast = 0.175;
else
  % no wait for scanner acq pulse at beginning of experiment
  waitForBacktick = 0;
  % contrast of the stimulus
  stimulus.contrast = 0.1;
end

% inner and outer width of stimulus
stimulus.innerWidth = 1.5;
stimulus.outerWidth = 10;
% spatial frequency of grating stimulus
stimulus.sf = 1;
% set to false for sharp edged stimulus
stimulus.gabor = 0;
% fixation cross width
stimulus.fixWidth = 1.5;

% orientations to display (these are the two different orientations)
stimulus.orientations = [25 115];
% jitter of 3 degrees around those standard orientations
stimulus.orientationJitter = 3;
% constant vals are the differences in orientations
% that the oberver is asked to discriminate (i.e if set to 3 and 6 then
% the psychophysics is done as method of constant stimuli where
% the probe stimulus is 3 or 6 degrees different from the memorized
% grating.
stimulus.constantVals = [3 6];

% delay interval in seconds
stimulus.delayInterval = 11;

% only sets these in debug mode for testing
if debugMode
  % sets delay interval shorted
  stimulus.delayInterval = 5;
  % make the task easier to see the difference in orientation
%   stimulus.constantVals = [20 30];
end

% initalize the screen, set the background to gray
myscreen.background = 128/255;
myscreen = initScreen(myscreen);

% by waiting for the backtick key to be pressed before starting the experiment
% (for systems that use NI digital I/O, this will wait for the digital
% signal that the scanner has started collecting data)
task{1}.waitForBacktick = 1;

% length in seconds of different segments of the trial
task{1}.segmin = [0.2 0.4 0.8 stimulus.delayInterval 0.5 2.5];
task{1}.segmax = [0.2 0.4 0.8 stimulus.delayInterval 0.5 2.5];
% synchToVol sets to wait for a scanner sync pulse after the end of the 
% segment - in this case the last segment will last 2.5 seconds after
% which the program will wait for the scanner sync pulse. This will insure
% that the beginning of each trial will be synchronized to the acquisition
% of an imaging volume.
task{1}.synchToVol = [0 0 0 0 0 waitForBacktick];
% get responses sets whether to aquire a keyboard response from the subject 
% we only need to set it to 1 for the response interval at the end of the trial.
task{1}.getResponse = [0 0 0 0 0 1];
% parameters
% cue can be 1 or 2 (meaning which stimulus, the first or second the subject should remember)
% this is picked from a uniform distribution, so just varies randomly from trial to trial
task{1}.randVars.uniform.cue = [1 2];
% oreintationOrder specifies which orientation will be shown 1st and which 2nd
task{1}.randVars.uniform.orientationOrder = [1 2];
% sets whether the probe stimulus at the end will be clockwise or counterclockwise of the remembered stimulus
task{1}.randVars.uniform.clockwiseCounterclockwise = [-1 1];
% these values are calculated during the experiment
% orientationJitter is the small amount of jitter that is 
% added to each orientation (its an array of 2, one fore each orientation)
task{1}.randVars.calculated.orientationJitter = [nan nan];
% oreintation is the actual oreintation shown
task{1}.randVars.calculated.orientation = [nan nan];

task{1}.randVars.calculated.memOrientation = nan;
% threshold is the threshold orientation shown
task{1}.randVars.calculated.orientationThreshold = nan;
% random sets to randomize parameters
task{1}.random = 1;

% initialize the task
for phaseNum = 1:length(task)
  [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end

% init the stimulus
myscreen = initStimulus('stimulus',myscreen);

% creates the stimulus images
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;
% on the first segment
if task.thistrial.thisseg == 1
  % set fixation white
  stimulus.fixColor = [1 1 1];

  % set the orientation jitter of both stimuli
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

  task.thistrial.memOrientation = task.thistrial.orientation(task.thistrial.cue);
  % display what we have selected
  disp(sprintf('Trial %i: Orientation 1: %0.1f + %0.1f: %f',task.trialnum,task.thistrial.orientation(1),task.thistrial.orientationJitter(1),stimulus.displayOrientation(1)));
  disp(sprintf('Trial %i: Orientation 2: %0.1f + %0.1f: %f',task.trialnum,task.thistrial.orientation(2),task.thistrial.orientationJitter(2),stimulus.displayOrientation(2)));
  disp(sprintf('Trial %i: Match to %i (%0.1f): %0.1f',task.trialnum,task.thistrial.cue,task.thistrial.orientationThreshold*task.thistrial.clockwiseCounterclockwise,stimulus.matchOrientation));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus

% clear screen to gray
mglClearScreen;

if task.thistrial.thisseg == 1
  % put up gratings to the left and right of fixation
  mglBltTexture(stimulus.tex,[-1.5*stimulus.outerWidth 0],0,0,stimulus.displayOrientation(1));
  mglBltTexture(stimulus.tex,[1.5*stimulus.outerWidth 0],0,0,stimulus.displayOrientation(2));
elseif (task.thistrial.thisseg == 3) 
  % draw a red line towards the cued side - one is left, two is right
  mglLines2((task.thistrial.cue-1.5)*3,0,(task.thistrial.cue-1.5)*5,0,2,[1,0,0]);
elseif task.thistrial.thisseg == 5
  % put up match either clockwise or counterclockwise 
  mglBltTexture(stimulus.tex,[0 0],0,0,stimulus.matchOrientation);
end  

% put up fixation cross, gluDisk can be used to make
% sure that there is a blank area around the fixation cross
%mglGluDisk(0,0,stimulus.innerWidth,0.5,48)
mglFixationCross(stimulus.fixWidth,1,stimulus.fixColor);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)

global stimulus

% figure out which button is the correct button for the trial
butts = [2 0 1];
corrButt = butts(task.thistrial.clockwiseCounterclockwise+2);

% here, we just check whether this is the first time we got a response
if task.thistrial.gotResponse < 1
  % if they got it correct
  if task.thistrial.whichButton == corrButt
    % display that with the reaction time
    disp(sprintf('Trial %i: Correct (reaction time: %0.2fs)',task.trialnum,task.thistrial.reactionTime));
    % change fixation cross to green
    stimulus.fixColor = [0 1 0];
    % and update the staircase
    stimulus.s = doStaircase('update',stimulus.s,1);
  else
    % incorrect
    disp(sprintf('Trial %i: Incorrect (reaction time: %0.2fs)',task.trialnum,task.thistrial.reactionTime));
    % change fixation color to red
    stimulus.fixColor = [1 0 0];
    % update the staircase
    stimulus.s = doStaircase('update',stimulus.s,0);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen)

% compute the grating with the appropriate size and spatial frequency
grating = mglMakeGrating(stimulus.outerWidth,stimulus.outerWidth,stimulus.sf,0,0);

% make a gaussian 
[gaussian xMesh yMesh] = mglMakeGaussian(stimulus.outerWidth,stimulus.outerWidth,stimulus.outerWidth/8,stimulus.outerWidth/8); 

% calculate radius of each point in the grating
r = sqrt(xMesh.^2 + yMesh.^2);

% Now we make a circular gaussian ring based on the radius of each point
gaussianRing = exp(-((r-(stimulus.outerWidth/4)).^2)/(stimulus.outerWidth/10)^2);
% set small values to 0 to try to make sure that the stimulus goes back to the background color gracefully
gaussianRing(gaussianRing < 0.01) = 0;

% hard edge
if ~stimulus.gabor
  % just make the gaussian ring a hard cutoff (either 1 or 0)
  gaussianRing(gaussianRing>0) = 1;
  gaussianRing(r<stimulus.innerWidth/2) = 0;
end

% now compute the "gabor" by multiplying the grating by the ring stimulus.
% we multiply that by the desired contrast. Then make those numbers go from 0 to 255
% so that they can be display as a texture.
gabor = round(255*(stimulus.contrast*grating.*gaussianRing+1)/2);

% create the texture
stimulus.tex = mglCreateTexture(gabor);

% make text for 1 and 2 that are displayed to subjects
mglTextSet('Helvetica',32,[1 1 1],0,0,0);
stimulus.text(1) = mglText('1');
stimulus.text(2) = mglText('2');

% start staircase
stimulus.s = doStaircase('init','fixed','fixedVals',stimulus.constantVals);