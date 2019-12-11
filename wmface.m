% wmface
%
%        $Id$
%      usage: wmface
%         by: justin gardner
%       date: 10/15/2019
%    purpose: working memory task based on Harrison & Tong with camera capture for Druckmann lab
%
%             Has both an auditory and visual version. Default is to run the auditory experiment,
%             but this can be swapped by setting the taskType
% 
%             wmface('taskType=visual');
%
%             For the auditory task, two tones are presented and then the word one or two.
%             Subject is instructed to remember the first or second tone depending on the
%             spoken cue. After a delay period, the subject will hear a third tone and
%             report whether the third tone is higher (button 1) or lower (button 2) then the
%             remembered tone.
%
%             For the visual task, subjects are shown two gratings of different orientations in
%             quick succession followed by the number 1 or 2 which indicates which grating orientation
%             the subject should remember. After a delay period, the subject is shown a third grating
%             and asked to report if it was rotated clockwise (button 1) or counterclockwise (button 2)
%             from the remembered stimulus.
%
%             Other arguments:
%
%             scan=0: Set to 1 to add a synchToVol at end of trial
%             delayInterval=5: Set to the desired delay interval in seconds
%             useCamera=1: Set to 0 to not use camera (useful for debugging w/out camera connected)
%             
%
function e = wmface(varargin)

% check arguments
getArgs(varargin,{'scan=0','taskType=auditory','delayInterval=5','useCamera=1','loFreq=440','hiFreq=560','dispStaircaseFig=1','numTrials=60','soundPath=~/proj/grustim/sounds'});
% default to non-scan mode
if nargin < 1,scan = 0;end

% global variable containing stimulus information
global stimulus;

% special settings for scanner
stimulus.scan = scan;
if stimulus.scan,waitForBacktick = 1;else, waitForBacktick = 0;end

% check taskType and set settings appropriately
switch lower(taskType)
 case 'auditory'
  % set auditory type
  stimulus.taskType = lower(taskType);

  % sampling frequency
  stimulus.audio.samplesPerSecond = 22000;
  % amplitude of stimulus
  stimulus.audio.amplitude = 0.1;
  % set the frequencies
  stimulus.audio.stimFreq = [loFreq hiFreq];
  % range of comparison frequencies in percentage
  % e.g. if set to 10, will compute frequencies that are -10 to 10 % of the stimulus frequency
  stimulus.audio.comparisonPercentRange = 0.20;
  % number of comparison freqencies in the range specified above - make this
  % an odd number so that we always include the sample frequency
  stimulus.audio.nComparisonFreq = 301;
  % set the initial threshold guess in percent
  stimulus.audio.initialThreshold = 0.05;
  % and to whether to display staircase figure
  stimulus.audio.dispStaircaseFig = dispStaircaseFig;
  % set background color
  myscreen.background = 0.5;
  % length in seconds to display stimulus for
  stimulus.stimLen = 0.5;
  stimulus.audio.len = 0.5;
  % set path to audio stimuli (i.e. "one" and "two")
  stimulus.audio.soundPath = soundPath;
  % set to use the audio cue (set to 0 for visual cue)
  stimulus.audio.audioCue = 1;

 case 'visual'
  % set visual type
  stimulus.taskType = lower(taskType);

  % set the contrast higher for in the scanner
  if scan,stimulus.contrast = 0.175;else, stimulus.contrast = 0.1;end
  
  % size of grating stimuli
  stimulus.innerWidth = 1.5;
  stimulus.outerWidth = 10;

  % spatial frequency of grating
  stimulus.sf = 1;

  % set to false for sharp edged stimulus
  stimulus.gabor = 0;

  % orientations to display
  stimulus.orientations = [25 115];
  stimulus.orientationJitter = 3;
  stimulus.constantVals = [3 6];
  
  % length in seconds to display stimulus for
  stimulus.stimLen = 0.2;
  
  % set background color
  myscreen.background = 128/255;
 otherwise
  disp(sprintf('(wmface) Unknown taskType: %s',taskType));
  return
end

% fixation width
stimulus.fixWidth = 1.5;


% delay interval in seconds
stimulus.delayInterval = delayInterval;

% initalize the screen
myscreen = initScreen(myscreen);

% intialize the camera
stimulus.useCamera = useCamera;
if stimulus.useCamera
  dispHeader('(cameraTest) Starting Camera Thread');
  tf = mglCameraThread('init');
  dispHeader;
  if ~tf
    endScreen(myscreen);
    dispHeader('(wmface) Could not initialize camera');
    return
  end
  stimulus.cameraImages = {};

  % place to save data
  stimulus.cameraDataDir = myscreen.datadir;
  if ~isdir(stimulus.cameraDataDir),mkdir(stimulus.cameraDataDir);end
  if isempty(myscreen.SID)
    stimulus.cameraFileStem = 'TEST';
  else
    stimulus.cameraFileStem = sprintf('%s',myscreen.SID);
  end
  stimulus.cameraFileStem = sprintf('%s_%s_%s',stimulus.cameraFileStem,datestr(now,'yyyymmdd'),datestr(now,'HHMMSS'));
else
  dispHeader;
  dispHeader('(wmface) Not using camera: useCamera set to 0');
  dispHeader;
end  

% by waiting for the backtick key to be pressed before starting the experiment
% (for systems that use NI digital I/O, this will wait for the digital
% signal that the scanner has started collecting data)
task{1}.waitForBacktick = 1;
task{1}.segmin = [stimulus.stimLen 0.4 stimulus.stimLen 0.4 0.8 stimulus.delayInterval stimulus.stimLen 5 1 1 inf];
task{1}.segmax = [stimulus.stimLen 0.4 stimulus.stimLen 0.4 0.8 stimulus.delayInterval stimulus.stimLen 5 1 1 inf];
task{1}.getResponse = [0 0 0 0 0 0 0 1 0 0];
task{1}.synchToVol = [0 0 0 0 0 0 0 waitForBacktick];
% parameters
task{1}.randVars.uniform.cue = [1 2];
task{1}.randVars.uniform.stimulusOrder = [1 2];

if strcmp(stimulus.taskType,'visual')
  % values for visual experiment
  task{1}.randVars.uniform.clockwiseCounterclockwise = [-1 1];
  task{1}.randVars.calculated.orientationJitter = [nan nan];
  task{1}.randVars.calculated.orientationJitter = [nan nan];
  task{1}.randVars.calculated.orientation = [nan nan];
  task{1}.randVars.calculated.deltaOrientation = nan;
else
  % values for audio experiment
  % whether to have the delta be a higher or lower frequency
  task{1}.randVars.uniform.highOrLow = [-1 1];
  % whether to put the delta on the first interval (stimulus to be remembered)
  % or the second interval (at end of working memory period)
  task{1}.randVars.uniform.deltaInterval = [1 2];
  % what the actual frequency in each interval was
  task{1}.randVars.calculated.freq = [nan nan nan];
  % and the frequency number from our precomputed table
  task{1}.randVars.calculated.freqNum = [nan nan nan];
  % what the current delta is that we are testing
  task{1}.randVars.calculated.deltaFreq = nan;
  task{1}.randVars.calculated.deltaPercent = nan;
end

% set correct answer
task{1}.randVars.calculated.correctAnswer = nan;
% set whether the subject got it correct or not
task{1}.randVars.calculated.correct = nan;

task{1}.random = 1;
task{1}.numTrials = numTrials;

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

% put up fix cross
mglClearScreen;mglFixationCross(stimulus.fixWidth,1,[1 1 1]);mglFlush;
mglClearScreen;mglFixationCross(stimulus.fixWidth,1,[1 1 1]);mglFlush;

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
if stimulus.useCamera
  dispHeader('(cameraTest) Ending Camera Thread');
  mglCameraThread('quit');
  dispHeader;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;
% set fixation white
if task.thistrial.thisseg == 1
  % start saving camera
  if stimulus.useCamera
    mglCameraThread('capture','timeToCapture',sum(task.segmax(1:(end-3))));
  end
  
  stimulus.fixColor = [1 1 1];

  % visual or auditory task
  if strcmp(stimulus.taskType,'visual')
    
    % set the orientation jitter of both
    task.thistrial.orientationJitter(1) = round(stimulus.orientationJitter*(2*(rand - 0.5)));
    task.thistrial.orientationJitter(2) = round(stimulus.orientationJitter*(2*(rand - 0.5)));

    % get the threshold value to test
    [task.thistrial.deltaOrientation stimulus.s] = doStaircase('testValue',stimulus.s);

    % set first and second orientation
    task.thistrial.orientation(1) = stimulus.orientations(task.thistrial.stimulusOrder);
    task.thistrial.orientation(2) = stimulus.orientations(setdiff([1 2],task.thistrial.stimulusOrder));

    % keep saved orientations to display
    stimulus.displayOrientation(1) = task.thistrial.orientation(1) + task.thistrial.orientationJitter(1);
    stimulus.displayOrientation(2) = task.thistrial.orientation(2) + task.thistrial.orientationJitter(2);
    stimulus.matchOrientation = task.thistrial.orientation(task.thistrial.cue)+task.thistrial.orientationJitter(task.thistrial.cue)+task.thistrial.deltaOrientation*task.thistrial.clockwiseCounterclockwise;

    % get correct answer
    butts = [2 0 1];
    task.thistrial.correctAnswer = butts(task.thistrial.clockwiseCounterclockwise+2);
    
    % display what we have selected
    disp(sprintf('Trial %i: Orientation 1: %0.1f + %0.1f: %f',task.trialnum,task.thistrial.orientation(1),task.thistrial.orientationJitter(1),stimulus.displayOrientation(1)));
    disp(sprintf('Trial %i: Orientation 2: %0.1f + %0.1f: %f',task.trialnum,task.thistrial.orientation(2),task.thistrial.orientationJitter(2),stimulus.displayOrientation(2)));
    disp(sprintf('Trial %i: Match to %i (%0.1f): %0.1f',task.trialnum,task.thistrial.cue,task.thistrial.deltaOrientation*task.thistrial.clockwiseCounterclockwise,stimulus.matchOrientation));
  else
    % set what the base frequencies are
    task.thistrial.freq(1) = stimulus.audio.stimFreq(task.thistrial.stimulusOrder);
    task.thistrial.freq(2) = stimulus.audio.stimFreq(setdiff([1 2],task.thistrial.stimulusOrder));
    
    % get the deltas we are testing
    [task.thistrial.deltaPercent stimulus.s] = doStaircase('testValue',stimulus.s);
    task.thistrial.deltaFreq = task.thistrial.freq(task.thistrial.cue)*task.thistrial.deltaPercent;
    
    % set the comparison stimulus at end to be which one is cued
    task.thistrial.freq(3) = task.thistrial.freq(task.thistrial.cue);
    
    % and then add the delta (with correct polarity) to the right stimulus
    if task.thistrial.deltaInterval == 1
      deltaStimNum = task.thistrial.cue;
    else
      deltaStimNum = 3;
    end
    
    % set the delta on the correct frequency
    task.thistrial.freq(deltaStimNum) = task.thistrial.freq(deltaStimNum)+task.thistrial.highOrLow*task.thistrial.deltaFreq;
      
    % compute the correct answer (which is 1 if the stimulus at the last interval is lower than 
    % the remembered one or 2 if it is higher)
    % the 1.5 and 2's are all about converting between [1 2] and [-1 1] so that the fina
    % answer should be 1 or 2
    task.thistrial.correctAnswer = 1.5+((task.thistrial.highOrLow * (2*(task.thistrial.deltaInterval-1.5)))/2);

    % now convert the frequencies into stimulus nums so that we can play them
    for iFreq = 1:3
      % look up the frequency
      [freqDiff task.thistrial.freqNum(iFreq)] = min(abs(task.thistrial.freq(iFreq)-stimulus.audio.comparisonFreq(:)));
      % store the actual frequency that was presented
      task.thistrial.freq(iFreq) = stimulus.audio.comparisonFreq(task.thistrial.freqNum(iFreq));
    end
    
    % now display what was done
    if task.thistrial.cue==1
      disp(sprintf('Trial %i: *%0.1f* %0.1f ------ %0.1f (Correct Answer: %i)',task.trialnum,task.thistrial.freq(1),task.thistrial.freq(2),task.thistrial.freq(3),task.thistrial.correctAnswer));
    else
      disp(sprintf('Trial %i: %0.1f *%0.1f* ------ %0.1f (Correct Answer: %i)',task.trialnum,task.thistrial.freq(1),task.thistrial.freq(2),task.thistrial.freq(3),task.thistrial.correctAnswer));
    end
  end
elseif task.thistrial.thisseg == 9
  stimulus.fixColor = [1 1 1];
elseif task.thistrial.thisseg == 10
  if stimulus.useCamera
    % create filename for images
    stimulus.saveCameraTime = mglGetSecs;
    filename = fullfile(stimulus.cameraDataDir,sprintf('%s-%04i',stimulus.cameraFileStem,task.trialnum));
    % tell thread to save images
    stimulus.cameraImages{end+1} = mglCameraThread('save','videoFilename',filename);
    % now save the task variables (just in case we crash)
    save(fullfile(stimulus.cameraDataDir,stimulus.cameraFileStem),'myscreen','task','stimulus');
    disp(sprintf('(wmface) Setting to save and task variable save took: %0.2fs',mglGetSecs(stimulus.saveCameraTime)));
  end
elseif task.thistrial.thisseg == 11
  if stimulus.useCamera
    % wait till camera is done
    mglCameraThread('blockTillDone');
    % display how long it took
    disp(sprintf('(wmface) Save of camera files took %0.2fs',mglGetSecs(stimulus.saveCameraTime)));
    % then jump to next segment (effectively end of trial)
  end
  % end the segment
  task = jumpSegment(task);
end

% play audio stimulus
if ~strcmp(stimulus.taskType,'visual')
 if task.thistrial.thisseg == 1
   mglPlaySound(task.thistrial.freqNum(1));
 elseif task.thistrial.thisseg == 3
   mglPlaySound(task.thistrial.freqNum(2));
 elseif task.thistrial.thisseg == 7
   mglPlaySound(task.thistrial.freqNum(3));
 end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus

mglClearScreen(0.5);

if strcmp(stimulus.taskType,'visual')
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
else
  if (task.thistrial.thisseg == 5)
    if stimulus.audio.audioCue
      % put up text
      mglPlaySound(stimulus.audio.cue(task.thistrial.cue));
    else
      % put up text
      mglBltTexture(stimulus.text(task.thistrial.cue),[0 1.5],0,0,0);
    end
  end  
end

% put up fixation cross
%mglGluDisk(0,0,stimulus.innerWidth,0.5,48)
mglFixationCross(stimulus.fixWidth,1,stimulus.fixColor);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)

global stimulus

% here, we just check whether this is the first time we got a response
% this trial and display what the subject's response was and the reaction time
if task.thistrial.gotResponse < 1
  if task.thistrial.whichButton == task.thistrial.correctAnswer
    disp(sprintf('Trial %i: Correct (reaction time: %0.2fs)',task.trialnum,task.thistrial.reactionTime));
    stimulus.fixColor = [0 1 0];
    stimulus.s = doStaircase('update',stimulus.s,1);
    task.thistrial.correct = true;
  else
    disp(sprintf('Trial %i: Incorrect (reaction time: %0.2fs)',task.trialnum,task.thistrial.reactionTime));
    stimulus.fixColor = [1 0 0];
    stimulus.s = doStaircase('update',stimulus.s,0);
    task.thistrial.correct = false;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen)

% setup visual stimulus
if strcmp(stimulus.taskType,'visual')

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
  stimulus.s = doStaircase('init','fixed','fixedVals',stimulus.constantVals);% setup auditory stimulus
else
  % make sure we have an odd number of frequencies to make it so that we have
  % the sample frequncy in the compariosn set
  stimulus.audio.nComparisonFreq = floor(stimulus.audio.nComparisonFreq/2)*2+1;
  
  % set the sample stimulus number
  stimulus.audio.sampleNum = ceil(stimulus.audio.nComparisonFreq/2);
  
  % uninstall all sounds
  mglInstallSound;
  
  % make a time sequence that goes over 2 pi every second
  t = 0:(stimulus.audio.len*2*pi)/((stimulus.audio.len*stimulus.audio.samplesPerSecond)-1):(stimulus.audio.len*2*pi);

  % cycle over stimulus frequencies
  for iFreq = 1:length(stimulus.audio.stimFreq)
    waveform = stimulus.audio.amplitude * sin(stimulus.audio.stimFreq(iFreq)*t);
    stimulus.audio.stim(iFreq) = mglInstallSound(waveform,stimulus.audio.samplesPerSecond);
  end
  
  % now make sounds that are in range around these frequencies
  stimulus.audio.comparisonFreqPercentage = -stimulus.audio.comparisonPercentRange:(2*stimulus.audio.comparisonPercentRange)/(stimulus.audio.nComparisonFreq-1):stimulus.audio.comparisonPercentRange;

  % cycle over comparison frequncies
  for iComparisonFreq = 1:stimulus.audio.nComparisonFreq
    % make the comparison frequency for each stimulus frequency
    for iStimFreq = 1:length(stimulus.audio.stimFreq)
      % compute the comparisonfrequency
      stimulus.audio.comparisonFreq(iStimFreq,iComparisonFreq) = stimulus.audio.stimFreq(iStimFreq) + ((stimulus.audio.comparisonFreqPercentage(iComparisonFreq))*stimulus.audio.stimFreq(iStimFreq));
      % compute the waveform
      waveform = stimulus.audio.amplitude * sin(stimulus.audio.comparisonFreq(iStimFreq,iComparisonFreq)*t);
      % create the sounds
      stimulus.audio.comparison(iStimFreq,iComparisonFreq) = mglInstallSound(waveform,stimulus.audio.samplesPerSecond);
    end
  end

  % make text
  mglTextSet('Helvetica',32,[1 1 1],0,0,0);
  stimulus.text(1) = mglText('1');
  stimulus.text(2) = mglText('2');

  if stimulus.audio.audioCue
    % load sounds
    stimulus.audio.oneName = fullfile(stimulus.audio.soundPath,'one.m4a');
    stimulus.audio.twoName = fullfile(stimulus.audio.soundPath,'two.m4a');

    % load one
    if ~isfile(stimulus.audio.oneName)
      disp(sprintf('(wmface) Could not load sound: %s',stimulus.audio.oneName));
      stimulus.audio.audioCue = 0;
    else
      [y Fs] = audioread(stimulus.audio.oneName);
      stimulus.audio.cue(1) = mglInstallSound(y',Fs);
    end

    % load two
    if ~isfile(stimulus.audio.twoName)
      disp(sprintf('(wmface) Could not load sound: %s',stimulus.audio.twoName));
      stimulus.audio.audioCue = 0;
    else
      [y Fs] = audioread(stimulus.audio.twoName);
      stimulus.audio.cue(2) = mglInstallSound(y',Fs);
    end
  end
  
% staircase 
  % min stepsize is determined by the difference in comparison freq percentages we have pre-computed
  minStepsize = median(diff(stimulus.audio.comparisonFreqPercentage));
  % max stepsize is set to a fraction of the total range
  maxStepsize = stimulus.audio.comparisonPercentRange/5;
  % initialize the staircase
  stimulus.s = doStaircase('init','upDown','initialThreshold',stimulus.audio.initialThreshold,'stepRule=pest','minThreshold=0','maxThreshold',stimulus.audio.comparisonPercentRange,'minStepsize',minStepsize,'maxStepsize',maxStepsize,'initialStepsize',maxStepsize/3,'dispFig',stimulus.audio.dispStaircaseFig);
end
