% twoPatchDiscrep.m
%
%        $Id$
%      usage: twoPatchDiscrep
%         by: josh wilson (dots by austin kuo)
%       date: 2/26/21
%  copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%    purpose: cue combination with 2 parameterizable dot stimuli
%
% NOTE: parameters you (might) want to change for this task are coherence, direction, density, speed.
% There are 2 patches so that you can have multiple motion components; they are individually parameterizable 

function myscreen = twoPatchDiscrep(varargin)

% check arguments
getArgs(varargin,'stimulusType=dots');

% initilaize the screen
myscreen = initScreen('hplp');
mglClearScreen
task{1}.waitForBacktick = 1;

global stimulus

% task parameters
task{1}.segmin = [3 2 10 2];
task{1}.segmax = [3 2 10 2];
task{1}.getResponse = [0 0 1 0];
task{1}.numTrials = 100;
task{1}.randVars.calculated.confidence = nan;
stimulus = initConfidence(stimulus,0,0,3,8,2,[1 1 1],[0.3 0.3 0.3]);
stimulus.feedback.segnum = 3;
% stimulus values are constant throughout experiment
stimulus.leftEcc = 0;
stimulus.rightEcc = 0;
stimulus.width = 10;
stimulus.speed = 5;
stimulus.contrast = .8;
stimulus.density = 2
% parameters are called *every* trial
task{1}.parameter.coherenceL = [.6 .8 1];
task{1}.parameter.coherenceR = [.6 .8 1];
task{1}.parameter.dirOffset = [100]; % 'right' dots are shifted clockwise this # of degrees from "dirI"
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


% initialize the task
for phaseNum = 1:length(task)
    [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
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
    
    % save trial parameters in a convenient place for analysis (personal preference)
    task.response.leftDir(task.trialnum) = stimulus.dirL(task.trialnum);
    task.response.rightDir(task.trialnum) = stimulus.dirR(task.trialnum);
    task.response.dirOffset(task.trialnum) = task.thistrial.dirOffset
    task.response.leftSpeed(task.trialnum) = stimulus.speedL(task.trialnum);
    task.response.rightSpeed(task.trialnum) = stimulus.speedR(task.trialnum);
    task.response.speedOffset(task.trialnum) = task.thistrial.speedOffset
    task.response.leftCoherence(task.trialnum) = task.thistrial.coherenceL;
    task.response.rightCoherence(task.trialnum) = task.thistrial.coherenceR;
    
    
    % set the fixation color
    stimulus.fixColor = stimulus.stimulatingFixColor;
    
elseif task.thistrial.thisseg == 2 | task.thistrial.thisseg == 4

    stimulus.dotsLeft = stimulus.dotsLeft.setCoherence(stimulus.dotsLeft,0);
    stimulus.dotsRight = stimulus.dotsRight.setCoherence(stimulus.dotsRight,0);
    
    % set the fixation color to indicate response period
    stimulus.fixColor = stimulus.responseFixColor;
    

elseif task.thistrial.thisseg == 3
  % set starting confidence
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
    task = jumpSegment(task);
    disp(sprintf('(estimation) Estimate: %0.2f',task.thistrial.confidence))
  end
end

% draw fixation cross
mglFixationCross(1,1,stimulus.fixColor);


%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dots stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initDots(stimulus,myscreen,task)

width = strcat('width=',num2str(stimulus.width));
contrast = strcat('contrast=',num2str(stimulus.contrast));
speed = strcat('speed=',num2str(stimulus.speed));
leftEcc = strcat('xCenter=',num2str(stimulus.leftEcc));
rightEcc = strcat('xCenter=',num2str(stimulus.rightEcc));
dotDensity = strcat('density=',num2str(stimulus.density));

% init the dot patches
stimulus.dotsLeft = dotsInit2Patch('framesPerSecond',myscreen.framesPerSecond,width,contrast,speed,leftEcc,dotDensity);
stimulus.dotsRight = dotsInit2Patch('framesPerSecond',myscreen.framesPerSecond,width,contrast,speed,rightEcc,dotDensity);

% set background color
stimulus.backgroundColor = 0.5;

% and fixation color
stimulus.stimulatingFixColor = [0 0 0];
stimulus.responseFixColor = [1 1 1];
stimulus.correctFixColor = [0 1 0];
stimulus.incorrectFixColor = [1 0 0];
stimulus.neutralFixColor = [0 0 1];









%%%%%%%%%%%%%%%%%%%%%%%%
%    initConfidence    %
%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initConfidence(stimulus,centerX,centerY,width,height,lineSize,lineColor,fillColor)

% set the segment in which the confidence judgement happens
stimulus.confidence.segnum = 8;
  
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
%fillY = stimulus.confidence.fillY;
%fillY(find(stimulus.confidence.fillTop)) = stimulus.confidence.centerY+(-0.5+confidenceLevel)*stimulus.confidence.height;
% now draw as a filled polygon
%mglPolygon(stimulus.confidence.fillX,fillY,stimulus.confidence.fillColor);

% draw outline
%mglLines2(stimulus.confidence.X0,stimulus.confidence.Y0,stimulus.confidence.X1,stimulus.confidence.Y1,stimulus.confidence.outlineSize,stimulus.confidence.outlineColor);

guessValue = sprintf('%.0f',confidenceLevel);
%mglTextDraw(guessValue,[0 0]);
mglLines2(0,0,7*cos(deg2rad(confidenceLevel)),7*sin(deg2rad(confidenceLevel)),3,[0 0 0])


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
  confidence = confidence+verticalScroll;
else
  horizontalScroll = 0;
end

% check bounds
if confidence > 360,confidence = 0;end
if confidence < 0,confidence = 360;end
  
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

function drawCorrect(task, stimulus)
cg = sprintf('%.0f',(task+20)*2.5);
mglTextDraw(cg,[0 0]);



