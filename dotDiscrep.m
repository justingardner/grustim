% twoPatchDiscrep.m
%
%        $Id$
%      usage: twoPatchDiscrep
%         by: josh wilson (dots by austin kuo)
%       date: 2/26/21
%  copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%    purpose: cue combination with 2 parameterizable dot stimuli
%
% NOTE: parameters you (might) want to change for this task are coherence, direction, speed.
% There are 2 patches so that you can have multiple motion components; they are individually parameterizable 

function myscreen = dotDiscrep(varargin)

% check arguments
getArgs(varargin,'stimulusType=dots');

% initilaize the screen
myscreen = initScreen('hplp');
mglClearScreen
task{1}.waitForBacktick = 1;

global stimulus

% task parameters
task{1}.segmin = [3 1 6 1];
task{1}.segmax = [3 1 6 1];
task{1}.getResponse = [0 0 1 0];
task{1}.numTrials = 100;
task{1}.randVars.calculated.confidence = nan;
% stimulus values are constant throughout experiment
stimulus.leftEcc = 0;
stimulus.rightEcc = 0;
stimulus.width = 10;
stimulus.speed = 5;
stimulus.contrast = .8;
stimulus.density = 2.5
% parameters are called *every* trial
task{1}.parameter.coherenceL = [1];
task{1}.parameter.coherenceR = [1];
task{1}.parameter.dirOffset = [60]; % 'right' dots are shifted clockwise this # of degrees from "dirI"
task{1}.parameter.speedOffset = [0];

% initialize responses
task{1}.response.leftDir = [];
task{1}.response.rightDir = [];
task{1}.response.dirOffset = [];
task{1}.response.leftSpeed = [];
task{1}.response.rightSpeed = [];
task{1}.response.speedOffset = [];
task{1}.response.leftCoherence = [];
task{1}.response.rightCoherence = [];
tasl{1}.response.estimate = [];


% initialize the task
for phaseNum = 1:length(task)
    [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback);
end

% init the stimulus
myscreen = initStimulus('stimulus',myscreen);
stimulus = initDots(stimulus,myscreen);


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
global stimulus
if task.thistrial.thisseg == 1
        
    % set the coherence
    stimulus.dotsLeft = stimulus.dotsLeft.setCoherence(stimulus.dotsLeft,task.thistrial.coherenceL);
    stimulus.dotsRight = stimulus.dotsRight.setCoherence(stimulus.dotsRight,task.thistrial.coherenceR);
  
    % set the motion direction of left and right patches 
    dirI = randi(360);
    stimulus.dirL(task.trialnum) = dirI;
    stimulus.dirR(task.trialnum) = dirI+task.thistrial.dirOffset;
        if stimulus.dirR(task.trialnum) > 360; stimulus.dirR(task.trialnum) = stimulus.dirR(task.trialnum)-360;end;
    stimulus.dotsLeft = stimulus.dotsLeft.setDir(stimulus.dotsLeft,stimulus.dirL(task.trialnum));
    stimulus.dotsRight = stimulus.dotsRight.setDir(stimulus.dotsRight,stimulus.dirR(task.trialnum));
    
     % set the speed of left and right patches
    stimulus.speedL(task.trialnum) = stimulus.speed;
    stimulus.speedR(task.trialnum) = stimulus.speed+task.thistrial.speedOffset;
    stimulus.dotsLeft = stimulus.dotsLeft.setSpeed(stimulus.dotsLeft,stimulus.speedL(task.trialnum));
    stimulus.dotsRight = stimulus.dotsRight.setSpeed(stimulus.dotsRight,stimulus.speedR(task.trialnum));

    % set the fixation color
    stimulus.fixColor = stimulus.stimulatingFixColor;
    
elseif (task.thistrial.thisseg == 2 | task.thistrial.thisseg == 4)

    stimulus.dotsLeft = stimulus.dotsLeft.setCoherence(stimulus.dotsLeft,0);
    stimulus.dotsRight = stimulus.dotsRight.setCoherence(stimulus.dotsRight,0);
    
    % set the fixation color to indicate response period
    stimulus.fixColor = stimulus.responseFixColor;
    

elseif task.thistrial.thisseg == 3
  % set starting confidence
  task.thistrial.confidence = 0;
  task.thistrial.gotResponse = 0;
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

if task.thistrial.thisseg == 1
    
    % update the dots
    stimulus.dotsLeft = stimulus.dotsLeft.update(stimulus.dotsLeft);
    stimulus.dotsRight = stimulus.dotsRight.update(stimulus.dotsRight);
    
    % draw the dots
    stimulus.dotsLeft = stimulus.dotsLeft.draw(stimulus.dotsLeft);
    stimulus.dotsRight = stimulus.dotsRight.draw(stimulus.dotsRight);
    
elseif task.thistrial.thisseg == 2 | task.thistrial.thisseg == 4
    
    % update the dots
    stimulus.dotsLeft = stimulus.dotsLeft.update(stimulus.dotsLeft);
    stimulus.dotsRight = stimulus.dotsRight.update(stimulus.dotsRight);
    
    % draw the dots
    stimulus.dotsLeft = stimulus.dotsLeft.draw(stimulus.dotsLeft);
    stimulus.dotsRight = stimulus.dotsRight.draw(stimulus.dotsRight);
    
end

if task.thistrial.thisseg == 3
  % set the confidence
  [task.thistrial.confidence confidenceDone] = setConfidence(task.thistrial.confidence, stimulus);
  if confidenceDone
    task.thistrial.gotResponse = 1
    % save trial parameters in a convenient place for analysis (personal preference) IF we get response
    task.response.leftDir(task.trialnum) = stimulus.dirL(task.trialnum);
    task.response.rightDir(task.trialnum) = stimulus.dirR(task.trialnum);
    task.response.dirOffset(task.trialnum) = task.thistrial.dirOffset
    task.response.leftSpeed(task.trialnum) = stimulus.speedL(task.trialnum);
    task.response.rightSpeed(task.trialnum) = stimulus.speedR(task.trialnum);
    task.response.speedOffset(task.trialnum) = task.thistrial.speedOffset
    task.response.leftCoherence(task.trialnum) = task.thistrial.coherenceL;
    task.response.rightCoherence(task.trialnum) = task.thistrial.coherenceR;
    task.response.estimate(task.trialnum) = task.thistrial.confidence;
    task = jumpSegment(task);
    disp(sprintf('(estimation) Estimate: %0.2f',task.thistrial.confidence))
  end
end

% draw fixation cross
mglFixationCross(1,1,stimulus.fixColor);


%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = getResponseCallback(task,myscreen)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dots stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initDots(stimulus,myscreen,task)

width = strcat('width=',num2str(stimulus.width));
contrast = strcat('contrast=',num2str(stimulus.contrast));
speed = strcat('speed=',num2str(stimulus.speed));
leftEcc = strcat('xCenter=',num2str(stimulus.leftEcc));
rightEcc = strcat('xCenter=',num2str(stimulus.rightEcc));
dotDense = strcat('density=',num2str(stimulus.density));

% init the dot patches
stimulus.dotsLeft = dotsInit2Patch('framesPerSecond',myscreen.framesPerSecond,width,contrast,speed,leftEcc,dotDense);
stimulus.dotsRight = dotsInit2Patch('framesPerSecond',myscreen.framesPerSecond,width,contrast,speed,rightEcc,dotDense);

% set background color
stimulus.backgroundColor = 0.5;

% and fixation color
stimulus.stimulatingFixColor = [0 0 0];
stimulus.responseFixColor = [1 1 1];
stimulus.correctFixColor = [0 1 0];
stimulus.incorrectFixColor = [1 0 0];
stimulus.neutralFixColor = [0 0 1];


%%%%%%%%%%%%%%%%%%%%%%%%
%       estimate       %
%%%%%%%%%%%%%%%%%%%%%%%%
function [confidence confidenceDone] = setConfidence(confidence, stimulus)

% get scroll
scrollEvents = mglListener('getAllScrollEvents');
if ~isempty(scrollEvents)
  % get the sum of the vertical and horizontal scrolls
  verticalScroll = sum(scrollEvents.scrollVertical);
  horizontalScroll = -sum(scrollEvents.scrollHorizontal);
  % set confidence
  confidence = confidence+verticalScroll;
else
  horizontalScroll = 0;
end

% check bounds
if confidence > 360,confidence = 0;end
if confidence < 0,confidence = 360;end

% draw the confidence
drawConfidence(confidence);

% if mouse button down (or horizontal scroll is non-zero) then we are done setting confidence
mouse = mglGetMouse;
keyEvent = mglGetKeyEvent;

if ~isequal(mouse.buttons,0) || ~isempty(keyEvent)
  confidenceDone = 1;
else
  confidenceDone = 0;
end

function drawConfidence(confidenceLevel)
guessValue = sprintf('%.0f',confidenceLevel*100);
%mglTextDraw(guessValue,[0 0]);
mglLines2(0,0,7*cos(deg2rad(confidenceLevel)),7*sin(deg2rad(confidenceLevel)),3,[0 0 0])



