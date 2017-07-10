% liaoTemplate.m
%
%        $Id$
%      usage: liaoTemplate
%         by: justin gardner
%       date: 07/97/2017
%  copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%    purpose: Template for contrast discrimination using dots at different eccentricities
%
function myscreen = liaoTemplate

% check arguments
if ~any(nargin == [0])
  help liao
  return
end

% initalize the screen
myscreen = initScreen;

% init the stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);

% init the dots
stimulus.width = 5;
stimulus.eccentricity = [5 10];
stimulus = initDotsStimulus(stimulus, myscreen);

% init the staircases
stimulus.pedestalContrast = 0.5;
initialThreshold = 0.2;
initialStepsize = 0.05;
% displays the staircase as trial data comes in a figure
dispStaircaseFig = 1;
% number of trials per staircase
nTrialsPerStaircase = 10;
% how many staircases to run until program ends
stimulus.nStaircases = 2;
% call with the above parameters
stimulus = initStaircases(stimulus, myscreen, initialThreshold, initialStepsize, nTrialsPerStaircase, dispStaircaseFig);

% fix: set waitForBacktick if you want to synch with the scanner
% by waiting for the backtick key to be pressed before starting the experiment
% (for systems that use NI digital I/O, this will wait for the digital
% signal that the scanner has started collecting data)
task{1}.waitForBacktick = 0;

task{1}.segmin = [0.5 0.5 2];
task{1}.segmax = [0.5 0.5 2];
task{1}.getResponse = [0 1 1];
task{1}.parameter.eccentricity = stimulus.eccentricity;
task{1}.randVars.uniform.whichSide = [1 2];
task{1}.random = 1;

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
  % set the eccentricity of the dots
  stimulus.dotsLeft = stimulus.dotsLeft.setCenter(stimulus.dotsLeft,-task.thistrial.eccentricity,0);
  stimulus.dotsRight = stimulus.dotsLeft.setCenter(stimulus.dotsRight,task.thistrial.eccentricity,0);

  % get which staircase
  stimulus.thisStaircase = find(stimulus.eccentricity==task.thistrial.eccentricity);
  stimulus.thisStaircase(2) = stimulus.staircaseCompleted(stimulus.thisStaircase)+1;

  % get the delta contrast to test
  [delta stimulus.s(stimulus.thisStaircase(1),stimulus.thisStaircase(2))] = doStaircase('testValue',stimulus.s(stimulus.thisStaircase(1),stimulus.thisStaircase(2)));

  % make sure delta does not go below zero
  stimulus.delta = max(0,delta);
  
  % get contrast for left and right
  if task.thistrial.whichSide == 1
    leftContrast = stimulus.pedestalContrast+stimulus.delta;
    rightContrast = stimulus.pedestalContrast;
  else
    leftContrast = stimulus.pedestalContrast;
    rightContrast = stimulus.pedestalContrast+stimulus.delta;
  end
  
  % set the contrast
  stimulus.dotsLeft = stimulus.dotsLeft.setContrast(stimulus.dotsLeft,leftContrast);
  stimulus.dotsRight = stimulus.dotsRight.setContrast(stimulus.dotsRight,rightContrast);

  % display what is going on
  disp(sprintf('%i: Eccentricity: %0.1f Side: %i delta: %0.2f',task.trialnum,task.thistrial.eccentricity,task.thistrial.whichSide,stimulus.delta));
  
  % set fix color to white
  stimulus.fixColor = [1 1 1];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus

mglClearScreen(0.5);

if task.thistrial.thisseg == 2
  % update the dots
  stimulus.dotsLeft = stimulus.dotsLeft.update(stimulus.dotsLeft);
  stimulus.dotsRight = stimulus.dotsRight.update(stimulus.dotsRight);

  % draw the dots
  stimulus.dotsLeft = stimulus.dotsLeft.draw(stimulus.dotsLeft);
  stimulus.dotsRight = stimulus.dotsRight.draw(stimulus.dotsRight);
end

% draw fixation cross
mglFixationCross(1,1,stimulus.fixColor);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)

global stimulus

% check the response
if task.thistrial.gotResponse < 1
  % see if it is correct
  if isequal(task.thistrial.whichButton,task.thistrial.whichSide)
    % report answer
    disp(sprintf(' !! Correct !!. Reaction time: %0.2f',task.thistrial.reactionTime));
    % change fixation color
    stimulus.fixColor = [0 1 0];
    % and update staircase
    stimulus.s(stimulus.thisStaircase(1),stimulus.thisStaircase(2)) = doStaircase('update',stimulus.s(stimulus.thisStaircase(1),stimulus.thisStaircase(2)),1,stimulus.delta);
  else
    % report answer
    disp(sprintf(' ++ Incorrect ++. Reaction time: %0.2f',task.thistrial.reactionTime));
    % change fixation color
    stimulus.fixColor = [1 0 0];
    % and update staircase
    stimulus.s(stimulus.thisStaircase(1),stimulus.thisStaircase(2)) = doStaircase('update',stimulus.s(stimulus.thisStaircase(1),stimulus.thisStaircase(2)),0,stimulus.delta);
  end    
  % see if we are done
  if doStaircase('stop',stimulus.s(stimulus.thisStaircase(1),stimulus.thisStaircase(2)))
    % update variable that says we are done
    stimulus.staircaseCompleted(stimulus.thisStaircase(1)) = stimulus.staircaseCompleted(stimulus.thisStaircase(1)) + 1;
    % and initialze the next staircase
    if stimulus.thisStaircase(2) < stimulus.nStaircases
      stimulus.s(stimulus.thisStaircase(1),stimulus.thisStaircase(2)+1) = doStaircase('init',stimulus.s(stimulus.thisStaircase(1),stimulus.thisStaircase(2)));
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initDotsStimulus(stimulus,myscreen)

% init the dot patchs
stimulus.dotsLeft = dotsInit('framesPerSecond',myscreen.framesPerSecond,'dir=0');
stimulus.dotsRight = dotsInit('framesPerSecond',myscreen.framesPerSecond,'dir=180');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the staircases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initStaircases(stimulus, myscreen,initialThreshold,initialStepsize,nTrials,dispStaircaseFig)

% get the last stimfile
stimfile = getLastStimfile(myscreen);

% check that the last stimfile had the same eccentricites
if ~isempty(stimfile)
  if ~isfield(stimfile,'stimulus') || ~isfield(stimfile.stimulus,'eccentricity') || ~isequal(stimfile.stimulus.eccentricity,stimulus.eccentricity)
    % dump this stimfile, since it does not match current eccentricity
    disp(sprintf('(liaoTemplate) Found stimfile, but does not have eccentricty match. Ignoring'));
    stimfile = [];
  end
end

% if no stimfile found
if isempty(stimfile)
  % then initialize
  for iStaircase = 1:length(stimulus.eccentricity)
    % print message of what we are doing
    disp(sprintf('(liaoTemplate) Initializing staircase for eccentricity: %0.2f',stimulus.eccentricity(iStaircase)));
    % init staircase
    stimulus.s(iStaircase,1) = doStaircase('init','upDown','nup=1','ndown=2','initialStepsize',initialStepsize,'nTrials',nTrials,'initialThreshold',initialThreshold,'subplotCols',length(stimulus.eccentricity),'subplotNum',iStaircase,'dispFig',dispStaircaseFig,'subplotName',sprintf('Eccentricity: %0.2f',stimulus.eccentricity(iStaircase)),'minThreshold',0,'stepRule=pest','maxStepsize',0.5,'minStepsize',0.005);
  end
else
  disp(sprintf('(liaoTemplate) Found stimfile'))
  % init using threshold from last stimfile
  for iStaircase = 1:length(stimulus.eccentricity)
    % print message of what we are doing
    threshold = doStaircase('threshold',stimfile.stimulus.s(iStaircase,:));
    disp(sprintf('(liaoTemplate) Initializing staircase for eccentricity: %0.2f from last stimfile with threshold: %0.2f',stimulus.eccentricity(iStaircase),threshold.threshold));
    % init staircase
    stimulus.s(iStaircase,1) = doStaircase('init',stimfile.stimulus.s(iStaircase,end));
  end
end

% set that the staircases are not yet done
stimulus.staircaseCompleted = zeros(1,length(stimulus.eccentricity));
