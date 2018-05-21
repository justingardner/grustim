% cogneuro_chichi.m
%
%      usage: cogneuro_workingmemory(scan) 
%             scan is a tf value which when set to 1 waits for
%             scanner sync pulses and sets the timing appropriately
%         by: justin gardner
%       date: 03/03/2016
%    purpose: Code for replicating Harrison & Tong, 2009
%       
function myscreen = cogneuro_chichi(scan,numJitter)

% check arguments
if ~any(nargin == [0 1 2])
  help cogneuro_workingmemory.m
  return
end
if nargin < 1, scan = 0; end
if nargin < 1, numJitter = 4; end

% set for debugging - makes the delay period shorter among other things.
debugMode = 1;

% stimulus parameters get stored in a global variable
global stimulus;

% parameters that are different for when scanning or not
stimulus.scan = scan;
if stimulus.scan
  % wait for scanner acq pulse at beginning of experiment
  waitForBacktick = 1;
  % contrast of the stimulus
  stimulus.contrast = 0.9;
else
  % no wait for scanner acq pulse at beginning of experiment
  waitForBacktick = 0;
  % contrast of the stimulus
  stimulus.contrast = 0.9;
end

% inner and outer width of stimulus
stimulus.outerWidth = 10;
% spatial frequency of grating stimulus
stimulus.sf = 1;
% set to false for sharp edged stimulus
stimulus.gabor = 0;
% stimulus eccentricity
stimulus.eccentricity = 7.5;
% fixation cross width
stimulus.fixWidth = 1.5;

% orientation jitter values
stimulus.orientationJitterValues = [-5 -15 -30 5 15 30];
stimulus.orientationJitterValues = 90;%[-5 -15 -30 5 15 30];

% delay interval in seconds
stimulus.delayInterval = 10;
stimulus.itimin = 3;
stimulus.itimax = 12;

% only sets these in debug mode for testing
if debugMode
  % sets delay interval shorted
  stimulus.delayInterval = 1;
  % make the task easier to see the difference in orientation
  waitForBacktick = 0;
  % set stimulus iti
  stimulus.itimin = 2;
  stimulus.itimax = 2;
end

% initalize the screen, set the background to gray
myscreen.background = 128/255;
myscreen = initScreen(myscreen);

% by waiting for the backtick key to be pressed before starting the experiment
% (for systems that use NI digital I/O, this will wait for the digital
% signal that the scanner has started collecting data)
task{1}.waitForBacktick = 1;

% length in seconds of different segments of the trial
task{1}.segmin = [1 stimulus.delayInterval 1 2 stimulus.itimin];
task{1}.segmax = [1 stimulus.delayInterval 1 2 stimulus.itimax];
% synchToVol sets to wait for a scanner sync pulse after the end of the 
% segment - in this case the last segment will last 2.5 seconds after
% which the program will wait for the scanner sync pulse. This will insure
% that the beginning of each trial will be synchronized to the acquisition
% of an imaging volume.
task{1}.synchToVol = [0 0 0 0 waitForBacktick];
% get responses sets whether to aquire a keyboard response from the subject 
% we only need to set it to 1 for the response interval at the end of the trial.
task{1}.getResponse = [0 0 1 1 0];
% parameters
task{1}.parameter.numJitter = numJitter;
task{1}.randVars.uniform.match = [0 1];
% orientation is the actual orientation shown
task{1}.randVars.calculated.orientation = [nan nan nan nan];
% orientation is the actual orientation shown
task{1}.randVars.calculated.orientationJitter = [nan nan nan nan];
% base orientation from which everything was randomized
task{1}.randVars.calculated.baseOrientation = [nan nan nan nan];
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

  % get a random base orientation
  task.thistrial.baseOrientation = round(rand(1,4)*180);

  % set all four orientations
  task.thistrial.orientation = task.thistrial.baseOrientation;

  % second phase may (or may not) have jitter
  task.thistrial.orientationJitter = task.thistrial.baseOrientation;
  
  % if a different trial (sameDifferent) then add orientation jitter
  if ~task.thistrial.match

    % figure out which stimuli will be jittered
    stimulusNums = randperm(4);
    jitterNums = stimulusNums(1:task.thistrial.numJitter);

    % figure out jitter values
    for i = 1:length(jitterNums)
      jitterVals(i) = stimulus.orientationJitterValues(randi(length(stimulus.orientationJitterValues)));
    end

    % set orientation jittered values
    task.thistrial.orientationJitter(jitterNums) = mod(task.thistrial.orientation(jitterNums) + jitterVals,180);
  end

  % display orientations
  disp(sprintf('Trial %i: [%i %i %i %i] -> [%i %i %i %i] numJitter: %i match: %i',task.trialnum,task.thistrial.orientation(1),task.thistrial.orientation(2),task.thistrial.orientation(3),task.thistrial.orientation(4),task.thistrial.orientationJitter(1),task.thistrial.orientationJitter(2),task.thistrial.orientationJitter(3),task.thistrial.orientationJitter(4),task.thistrial.numJitter,task.thistrial.match));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus

% clear screen to gray
mglClearScreen;

if task.thistrial.thisseg == 1
  % put up four grating in first semgent
  mglBltTexture(stimulus.tex,stimulus.loc(1,:),0,0,task.thistrial.orientation(1));
  mglBltTexture(stimulus.tex,stimulus.loc(2,:),0,0,task.thistrial.orientation(2));
  mglBltTexture(stimulus.tex,stimulus.loc(3,:),0,0,task.thistrial.orientation(3));
  mglBltTexture(stimulus.tex,stimulus.loc(4,:),0,0,task.thistrial.orientation(4));
elseif task.thistrial.thisseg == 3 
  % put up match stimulus
  mglBltTexture(stimulus.tex,stimulus.loc(1,:),0,0,task.thistrial.orientationJitter(1));
  mglBltTexture(stimulus.tex,stimulus.loc(2,:),0,0,task.thistrial.orientationJitter(2));
  mglBltTexture(stimulus.tex,stimulus.loc(3,:),0,0,task.thistrial.orientationJitter(3));
  mglBltTexture(stimulus.tex,stimulus.loc(4,:),0,0,task.thistrial.orientationJitter(4));
elseif task.thistrial.thisseg == 5 
  stimulus.fixColor = [1 1 1];
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
if task.thistrial.match
  corrButt = 1;
else
  corrButt = 2;
end

% here, we just check whether this is the first time we got a response
if task.thistrial.gotResponse < 1
  % if they got it correct
  if task.thistrial.whichButton == corrButt
    % display that with the reaction time
    disp(sprintf('Trial %i: Correct (reaction time: %0.2fs)',task.trialnum,task.thistrial.reactionTime));
    % change fixation cross to green
    stimulus.fixColor = [0 1 0];
  else
    % incorrect
    disp(sprintf('Trial %i: Incorrect (reaction time: %0.2fs)',task.trialnum,task.thistrial.reactionTime));
    % change fixation color to red
    stimulus.fixColor = [1 0 0];
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


% now compute the gabor by multiplying the grating by the gaussian
% we multiply that by the desired contrast. Then make those numbers go from 0 to 255
% so that they can be display as a texture.
gabor = round(255*(stimulus.contrast*grating.*gaussian+1)/2);

% create the texture
stimulus.tex = mglCreateTexture(gabor);

% make text for 1 and 2 that are displayed to subjects
mglTextSet('Helvetica',32,[1 1 1],0,0,0);
stimulus.text(1) = mglText('1');
stimulus.text(2) = mglText('2');

% set stimulus locations
stimulus.loc(1,:) = round([cos(pi/4) sin(pi/4)]*stimulus.eccentricity);
stimulus.loc(2,:) = round([cos(3*pi/4) sin(3*pi/4)]*stimulus.eccentricity);
stimulus.loc(3,:) = round([cos(5*pi/4) sin(5*pi/4)]*stimulus.eccentricity);
stimulus.loc(4,:) = round([cos(7*pi/4) sin(7*pi/4)]*stimulus.eccentricity);
