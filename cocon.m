% cocon.m
%
%        $Id$
%      usage: cocon
%         by: justin gardner
%       date: 02/25/2020
%  copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%    purpose: 2AFC motion direction paradigm for testing COherence CONtext
%
%             getConfidence=1 (sets whether to run confidence measurement)
%             taskType=newsome (which task, can be newsome or dirdiff)
%             runType=staircase (first run as staircase to get threshold. Then run as low or high
%               for low or high cohernce context);
%
function myscreen = cocon(varargin)

% check arguments
getArgs(varargin,{'stimulusType=dots','getConfidence=1','scannerMode=0','taskType=newsome','runType=staircase'});

% initalize the screen
myscreen.background = 0.5;
myscreen = initScreen(myscreen.background);

% init the stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);

% init info about what kind of stimulus and task we are doing
stimulus.stimulusType = lower(stimulusType);
stimulus.getConfidence = getConfidence;
stimulus.taskType = lower(taskType);
stimulus.runType = lower(runType);

% set which variable to staircase
stimulus.staircase = 'coherence';
% set which variable over which each staircase is made
% for example, you might want to have one staricase for
% each direction difference
stimulus.staircaseOver = 'delta';

% stimulus parameters
stimulus.width = 10;
stimulus.speed = 2.5;
stimulus.eccentricity = [7.5];
stimulus.coherence = [0.25 0.5];
stimulus.delta = 15;

% select the range of values to use for low and high coherence. These are multiplicative of threshold
stimulus.lowVals = [0.5 1 1.5];
stimulus.highVals = [1.5 2 2.5];

% init the stimulus
if strcmp(stimulus.stimulusType,'dots')
  % init the dots
  stimulus = initDotsStimulus(stimulus, myscreen);
  stimulus.isGrating = false;
else
  disp(sprintf('(cocon) Unknown stimulus type: %s',stimulus.stimulusType));
  endScreen(myscreen);
  return
end

if strcmp(stimulus.runType,'staircase')
  % init the staircases
  initialThreshold = 25;
  initialStepsize = 5;
  % displays the staircase as trial data comes in a figure
  dispStaircaseFig = 1;
  % number of trials per staircase
  nTrialsPerStaircase = 40;
  % how many staircases to run until program ends
  stimulus.nStaircases = 2;
  % whether to try to start from the last computed threshold
  useLastThreshold = false;
  % call with the above parameters
  stimulus = initStaircases(stimulus, myscreen, initialThreshold, initialStepsize, nTrialsPerStaircase, dispStaircaseFig, useLastThreshold);
elseif any(strcmp(stimulus.runType,{'low','high'}))
  % this is a regular run, set number of trials
  stimulus.nTrials = 150;
  % get threshold
  stimulus.threshold = getThreshold(stimulus, myscreen);
  if isempty(stimulus.threshold)
    endScreen(myscreen);
    return
  end
  % now init constant stimulus
  stimulus = initConstantStimuli(stimulus);
end

% wait for backtick to start
task{1}.waitForBacktick = scannerMode;

% task timing and response
if stimulus.getConfidence
  % init the confidence display parameters (see function definition below for definition of parameters)
  stimulus = initConfidence(stimulus,0,0,3,8,2,[1 1 1],[0.3 0.3 0.3]);
  stimulus.feedback.segnum = 5;

  % the fifth segment is the one in which the confidence judgement comes up
  task{1}.segmin = [0.5 1 inf inf 0.5];
  task{1}.segmax = [0.5 1 inf inf 0.5];
  task{1}.getResponse = [0 1 1 0 0];
else
  % set feedback segment
  stimulus.feedback.segnum = [2 3];
  % not confidence judgement
  task{1}.segmin = [0.5 1 2];
  task{1}.segmax = [0.5 1 2];
  task{1}.getResponse = [0 1 1];
end

if ~strcmp(stimulus.runType,'staircase')
    % no feedback
    stimulus.feedback.segnum = nan;
end

% taks paraeters
task{1}.parameter.eccentricity = stimulus.eccentricity;
task{1}.randVars.calculated.confidnece = nan;
task{1}.randVars.calculated.correctIncorrect = nan;
task{1}.randVars.uniform.whichAnswer = [1 2];
task{1}.random = 1;
task{1}.randVars.calculated.direction = nan;

% now set the variable that we will staircase on
task{1}.randVars.calculated.(stimulus.staircase) = nan;
% and the variable that we will staircase over
task{1}.parameter.(stimulus.staircaseOver) = stimulus.(stimulus.staircaseOver);

% initialize the task
for phaseNum = 1:length(task)
  [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while ~all(stimulus.staircaseCompleted == stimulus.nStaircases) && ~myscreen.userHitEsc
  % update the task
  [task myscreen phaseNum] = updateTask(task,myscreen,phaseNum);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

if ~myscreen.userHitEsc
  % display thresholds
  for iStaircase = 1:length(stimulus.eccentricity);
    % compute threshold
    t = doStaircase('threshold',stimulus.s(iStaircase,:));
    % display it
    disp(sprintf('Eccentricity: %0.2f Threshold: %0.2f',stimulus.eccentricity(iStaircase),t.threshold));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;

if task.thistrial.thisseg == 1

  % get which staircase, note that we have kxn staircases - k coherences and n staircases
  stimulus.whichStaircase = find(stimulus.(stimulus.staircaseOver)==task.thistrial.(stimulus.staircaseOver));
  stimulus.staircaseNum = stimulus.staircaseCompleted(stimulus.whichStaircase)+1;

  % get the staircase value to test
  [staircaseVal stimulus.s(stimulus.whichStaircase,stimulus.staircaseNum)] = doStaircase('testValue',stimulus.s(stimulus.whichStaircase,stimulus.staircaseNum));

  % now set the values
  task.thistrial.(stimulus.staircase) = staircaseVal;

  % for dirdiff task, we set up two patches
  if strcmp(stimulus.taskType,'dirdiff')
    % make sure delta does not go below zero
    task.thistrial.delta = max(0,task.thistrial.delta);

    % set  direction
    task.thistrial.direction = round(rand*360);
    if task.thistrial.whichAnswer == 1
      leftDirection = mod(task.thistrial.direction+stimulus.delta,360);
      rightDirection = task.thistrial.direction;
    else
      leftDirection = task.thistrial.direction;
      rightDirection = mod(task.thistrial.direction+stimulus.delta,360);
    end

    % display what is going on
    disp(sprintf('%i: Coherence: %0.1f Side: %i delta: %0.2f (Direction %0.1f vs %0.1f Eccentricity: %0.1f)',task.trialnum,task.thistrial.coherence,task.thistrial.whichAnswer,stimulus.delta,leftDirection,rightDirection,task.thistrial.eccentricity));
  
    % set the coherence
    stimulus.dotsLeft = stimulus.dotsLeft.setCoherence(stimulus.dotsLeft,task.thistrial.coherence/100);
    stimulus.dotsRight = stimulus.dotsRight.setCoherence(stimulus.dotsRight,task.thistrial.coherence/100);

    % set the position
    stimulus.dotsLeft = stimulus.dotsLeft.setCenter(stimulus.dotsLeft,-task.thistrial.eccentricity,0);
    stimulus.dotsRight = stimulus.dotsRight.setCenter(stimulus.dotsRight,task.thistrial.eccentricity,0);

    % set the direction
    stimulus.dotsLeft = stimulus.dotsLeft.setDir(stimulus.dotsLeft,leftDirection);
    stimulus.dotsRight = stimulus.dotsRight.setDir(stimulus.dotsRight,rightDirection);
  elseif strcmp(stimulus.taskType,'newsome')
    % do a left/right discrimination task
    stimulus.dots = stimulus.dots.setCoherence(stimulus.dots,task.thistrial.coherence/100);
    stimulus.dots = stimulus.dots.setCenter(stimulus.dots,0,task.thistrial.eccentricity);
    % set the direction to left if answer is 1 and right for answer is 2
    stimulus.dots = stimulus.dots.setDir(stimulus.dots,180*(1-(task.thistrial.whichAnswer-1)));
    disp(sprintf('%i: Coherence: %0.1f Direction: %i',task.trialnum,task.thistrial.coherence,180*(1-(task.thistrial.whichAnswer-1))));
  end

  % set the fixation color
  stimulus.fixColor = stimulus.normalFixColor;

% confidence segment
elseif stimulus.getConfidence && (task.thistrial.thisseg == stimulus.confidence.segnum)
  % set starting confidnece
  task.thistrial.confidence = 0.5;
  scrollEvents = mglListener('getAllScrollEvents');
  mglListener('getAllMouseEvents');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus

% clear screen to gray
mglClearScreen(stimulus.backgroundColor);

if task.thistrial.thisseg == 2

  if strcmp(stimulus.taskType,'dirdiff')
    % update the dots
    stimulus.dotsLeft = stimulus.dotsLeft.update(stimulus.dotsLeft);
    stimulus.dotsRight = stimulus.dotsRight.update(stimulus.dotsRight);

    % draw the dots
    stimulus.dotsLeft = stimulus.dotsLeft.draw(stimulus.dotsLeft);
    stimulus.dotsRight = stimulus.dotsRight.draw(stimulus.dotsRight);
  else
    % update and draw single patch
    stimulus.dots = stimulus.dots.update(stimulus.dots);
    stimulus.dots = stimulus.dots.draw(stimulus.dots);
  end
    
elseif stimulus.getConfidence && (task.thistrial.thisseg == stimulus.confidence.segnum)
  % set the confidence
  [task.thistrial.confidence confidenceDone] = setConfidence(task.thistrial.confidence, stimulus);
  if confidenceDone
    task = jumpSegment(task);
    disp(sprintf('(cocon) Confidence: %0.2f',task.thistrial.confidence));
  end
end

% feedback segment
if any(task.thistrial.thisseg == stimulus.feedback.segnum)
  % if correct
  if isequal(task.thistrial.correctIncorrect,1)
    stimulus.fixColor = stimulus.correctFixColor;
  % if incorrect
  elseif isequal(task.thistrial.correctIncorrect,0)
    stimulus.fixColor = stimulus.incorrectFixColor;
  end
end

% draw fixation cross
mglFixationCross(1,1,stimulus.fixColor);

%%%%%%%%%%%%%%%%%%%%%%%%
%    initConfidence    %
%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initConfidence(stimulus,centerX,centerY,width,height,lineSize,lineColor,fillColor)

% set the segment in which the confidence judgement happens
stimulus.confidence.segnum = 4;
  
% set the dimensions of the confidence display
stimulus.confidence.width = width;
stimulus.confidence.height = height;
stimulus.confidence.centerX = centerX;
stimulus.confidence.centerY = centerY;

% set line size and color
stimulus.confidence.outlineSize = lineSize;
stimulus.confidence.outlineColor = lineColor;
stimulus.confidence.fillColor = fillColor;

% now compute the x and y of the outline
xL = stimulus.confidence.centerX-stimulus.confidence.width/2;
xR = stimulus.confidence.centerX+stimulus.confidence.width/2;
yB = stimulus.confidence.centerY-stimulus.confidence.height/2;
yT = stimulus.confidence.centerY+stimulus.confidence.height/2;

% dimensions of rectangle
stimulus.confidence.X0 = [xL xR xR xL];
stimulus.confidence.X1 = [xR xR xL xL];
stimulus.confidence.Y0 = [yB yB yT yT];
stimulus.confidence.Y1 = [yB yT yT yB];

% dimensions of fill
stimulus.confidence.fillX = [xL xR xR xL];
stimulus.confidence.fillY = [yB yB yT yT];

% values that need to be changed to reflect confidence level
stimulus.confidence.fillTop = [0 0 1 1];

% turn off scrolling
mglListener('eatscroll');
% read all pending scroll events
scrollEvents = mglListener('getAllScrollEvents');

%%%%%%%%%%%%%%%%%%%%%%%%
%    drawConfidence    %
%%%%%%%%%%%%%%%%%%%%%%%%
function drawConfidence(confidenceLevel, stimulus)

% draw confidence level as a filled bar.

% draw filled inside, compute top coordinate
fillY = stimulus.confidence.fillY;
fillY(find(stimulus.confidence.fillTop)) = stimulus.confidence.centerY+(-0.5+confidenceLevel)*stimulus.confidence.height;
% now draw as a filled polygon
mglPolygon(stimulus.confidence.fillX,fillY,stimulus.confidence.fillColor);

% draw outline
mglLines2(stimulus.confidence.X0,stimulus.confidence.Y0,stimulus.confidence.X1,stimulus.confidence.Y1,stimulus.confidence.outlineSize,stimulus.confidence.outlineColor);

%%%%%%%%%%%%%%%%%%%%%%%
%    setConfidence    %
%%%%%%%%%%%%%%%%%%%%%%%
function [confidence confidenceDone] = setConfidence(confidence, stimulus)

% get scroll
scrollEvents = mglListener('getAllScrollEvents');
if ~isempty(scrollEvents)
  % get the sum of the vertical and horizontal scrolls
  verticalScroll = -sum(scrollEvents.scrollVertical);
  horizontalScroll = sum(scrollEvents.scrollHorizontal);
  % set confidence
  confidence = confidence+verticalScroll/100;
else
  horizontalScroll = 0;
end

% check bounds
if confidence > 1,confidence = 1;end
if confidence < 0,confidence = 0;end
  
% draw the confidence
drawConfidence(confidence,stimulus);

% if mouse button down (or horizontal scroll is non-zero) then we are done setting confidence
mouse = mglGetMouse;
keyEvent = mglGetKeyEvent;

if ~isequal(mouse.buttons,0) || ~isempty(keyEvent)
  confidenceDone = 1;
else
  confidenceDone = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)

global stimulus

% check the response
if task.thistrial.gotResponse < 1
  % see if it is correct
  if isequal(task.thistrial.whichButton,task.thistrial.whichAnswer)
    % report answer
    disp(sprintf(' !! Correct !!. Reaction time: %0.2f',task.thistrial.reactionTime));
    % save it in task
    task.thistrial.correctIncorrect = 1;
    % change fixation color
    stimulus.fixColor = stimulus.responseFixColor;
    % and update staircase
    stimulus.s(stimulus.whichStaircase,stimulus.staircaseNum) = doStaircase('update',stimulus.s(stimulus.whichStaircase,stimulus.staircaseNum),1,task.thistrial.(stimulus.staircase));
  else
    % report answer
    disp(sprintf(' ++ Incorrect ++. Reaction time: %0.2f',task.thistrial.reactionTime));
    % save it in task
    task.thistrial.correctIncorrect = 0;
    % change fixation color
    stimulus.fixColor = stimulus.responseFixColor;
    % and update staircase
    stimulus.s(stimulus.whichStaircase,stimulus.staircaseNum) = doStaircase('update',stimulus.s(stimulus.whichStaircase,stimulus.staircaseNum),0,task.thistrial.(stimulus.staircase));
  end    
  % see if we are done
  if doStaircase('stop',stimulus.s(stimulus.whichStaircase,stimulus.staircaseNum))
    % update variable that says we are done
    stimulus.staircaseCompleted(stimulus.whichStaircase) = stimulus.staircaseCompleted(stimulus.whichStaircase) + 1;
    % and initialze the next staircase
    if stimulus.staircaseNum < stimulus.nStaircases
      stimulus.s(stimulus.whichStaircase,stimulus.staircaseNum+1) = doStaircase('init',stimulus.s(stimulus.whichStaircase,stimulus.staircaseNum));
    end
  end
  if stimulus.getConfidence
    % jump to confidence rating part of text  
    task = jumpSegment(task);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dots stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initDotsStimulus(stimulus,myscreen)

if strcmp(stimulus.taskType,'dirdiff')
  % init the dot patchs
  stimulus.dotsLeft = dotsInit('framesPerSecond',myscreen.framesPerSecond,'dir=0','width',stimulus.width,'speed',stimulus.speed);
  stimulus.dotsRight = dotsInit('framesPerSecond',myscreen.framesPerSecond,'dir=180','width',stimulus.width,'speed',stimulus.speed);
elseif strcmp(stimulus.taskType,'newsome')
  stimulus.dots = dotsInit('framesPerSecond',myscreen.framesPerSecond,'dir=180','width',stimulus.width,'speed',stimulus.speed);
else
  disp(sprintf('(cocon:initDotsStimulus) Unknown task type: %s',stimulus.taskType));
  keyboard
end

% set background color
stimulus.backgroundColor = 0.5;

% and fixation color
stimulus.normalFixColor = [1 1 1];
stimulus.correctFixColor = [0 1 0];
stimulus.incorrectFixColor = [1 0 0];
stimulus.responseFixColor = [1 1 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the staircases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initStaircases(stimulus,myscreen,initialThreshold,initialStepsize,nTrials,dispStaircaseFig,useLastThreshold)

if ~useLastThreshold
  % if we are not using last threshold, 
  % then don't bother trying to load last stimfile
  stimfile = [];
else
  % get the last stimfile
  stimfile = getLastStimfile(myscreen);
end

% check that the last stimfile had the same values
if ~isempty(stimfile)
  if ~isfield(stimfile,'stimulus') || ~isfield(stimfile.stimulus,stimulus.staircaseOver) || ~isequal(stimfile.stimulus.(stimulus.staircaseOver),stimulus.(stimulus.staircaseOver))
    % dump this stimfile, since it does not match current eccentricity
    disp(sprintf('(cocon) Found stimfile, but does not have %s match. Ignoring',stimulus.staircaseOver));
    stimfile = [];
  end
end

% get the variable to staircase over
staircaseOver = stimulus.(stimulus.staircaseOver);
nStaircase = length(staircaseOver);

% if no stimfile found
if isempty(stimfile)
  % then initialize
  for iStaircase = 1:nStaircase
    % print message of what we are doing
    disp(sprintf('(cocon) Initializing staircase for %s: %0.2f',stimulus.staircaseOver,staircaseOver(iStaircase)));
    % init staircase
    stimulus.s(iStaircase,1) = doStaircase('init','upDown','nup=1','ndown=2','initialStepsize',initialStepsize,'nTrials',nTrials,'initialThreshold',initialThreshold,'subplotCols',nStaircase,'subplotNum',iStaircase,'dispFig',dispStaircaseFig,'subplotName',sprintf('%s: %0.2f',stimulus.staircaseOver,staircaseOver(iStaircase)),'minThreshold',0,'stepRule=pest','maxStepsize',0.5,'minStepsize',0.005);
  end
else
  disp(sprintf('(cocon) Found stimfile'))
  % init using threshold from last stimfile
  for iStaircase = 1:nStaircase
    % print message of what we are doing
    threshold = doStaircase('threshold',stimfile.stimulus.s(iStaircase,:));
    disp(sprintf('(cocon) Initializing staircase for %s: %0.2f from last stimfile with threshold: %0.2f',stimulus.staircaseOver,staircaseOver(iStaircase),threshold.threshold));
    % init staircase
    stimulus.s(iStaircase,1) = doStaircase('init',stimfile.stimulus.s(iStaircase,end));
  end
end

% set that the staircases are not yet done
stimulus.staircaseCompleted = zeros(1,nStaircase);

% and clear threshold field
stimulus.threshold = [];
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to get threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function threshold = getThreshold(stimulus,myscreen)

% get the last stimfile
stimfile = getLastStimfile(myscreen);
if isempty(stimfile)
  disp(sprintf('(cocon:getThreshold) No stimfile found, you need to run a staircase first!'));
  threshold = [];
  return
end

disp(sprintf('(cocon) Found stimfile'));

% now see if we have the proper staircases run for these conditions
% check that the last stimfile had the same values
if ~isfield(stimfile,'stimulus') || ~isfield(stimfile.stimulus,stimulus.staircaseOver) || ~isequal(stimfile.stimulus.(stimulus.staircaseOver),stimulus.(stimulus.staircaseOver))
  % dump this stimfile, since it does not match current variable settings
  disp(sprintf('(cocon:getThreshold) Found stimfile, but does not have %s match. Need to run staircase again.',stimulus.staircaseOver));
  threshold = [];
  return
end

% see if threshold already exists
if isfield(stimulus,'threshold') && ~isempty(stimulus.threshold)
  disp(sprintf('(cocon:getThreshold) Found existing threshold: %s',num2str(stimulus.threshold.threshold)));
  threshold = stimulus.threshold;
  return
end
    
% get number of staircases
staircaseOver = stimulus.(stimulus.staircaseOver);
nStaircase = length(staircaseOver);
nRepeats = size(stimfile.stimulus.s,2);

% now get the threshold values
for iStaircase = 1:nStaircase
  staircases = {};
  % collect all staircases with more than 80% of trials run
  for iRepeat = 1:nRepeats
    % get the staircase
    thisStaircase = stimfile.stimulus.s(iStaircase,iRepeat);
    % get how many trials were run
    if (thisStaircase.trialNum/thisStaircase.stopCriterion) > 0.8
      disp(sprintf('(cocon:getThreshold) Using staircase for %s: %0.2f (%i/%i) has (%i/%i) trials',stimulus.staircaseOver,staircaseOver(iStaircase),iRepeat,nRepeats,thisStaircase.trialNum,thisStaircase.stopCriterion));
      staircases{end+1} = thisStaircase;
    else
      disp(sprintf('(cocon:getThreshold) IGNORING Staircase for %s: %0.2f (%i/%i) has (%i/%i) trials - TOO FEW TRIALS',stimulus.staircaseOver,staircaseOver(iStaircase),iRepeat,nRepeats,thisStaircase.trialNum,thisStaircase.stopCriterion));
    end
  end
  % no staircases found, return empty
  if isempty(staircases)
    disp(sprintf('(cocon:getThreshold) No staircases found'));
    threshold = [];
    return
  end
  % get threshold
  threshold(iStaircase) = doStaircase('threshold',cell2mat(staircases));
  % print message of what we are doing
  disp(sprintf('(cocon:getThreshold) Staircase for %s: %0.2f from last stimfile had threshold: %0.2f',stimulus.staircaseOver,staircaseOver(iStaircase),threshold.threshold));
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init constant stimuli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initConstantStimuli(stimulus,myscreen)

% get number of thresholds
nThreshold = length(stimulus.threshold);

% set doStaircase to run constant stimuli
for iThreshold = 1:nThreshold
  threshold = stimulus.threshold(iThreshold).threshold;
  % get the fixed values
  if strcmp(stimulus.runType,'low')
    fixedVals = threshold * stimulus.lowVals;
  else    
    fixedVals = threshold * stimulus.highVals;
  end
  % display what we are doing
  disp(sprintf('(cocon:initConstantStimuli) Running %s with values: %s (threshold=%f)',stimulus.runType,num2str(fixedVals),threshold));
  % set the staircase to run constant stimuli
  stimulus.s(iThreshold,1) = doStaircase('init','fixed','fixedVals',fixedVals,'nTrials',stimulus.nTrials,'dispFig',true);
end

% set this so the code knows wow many staircases have been run (in this case staircase is a set of constant stimuli)
stimulus.staircaseCompleted = zeros(1,nThreshold);
stimulus.nStaircases = 1;
