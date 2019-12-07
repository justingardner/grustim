% twoPatchMotionLocalizer.m
%
%        $Id$
%      usage: twoPatchMotionLocalizer
%         by: austin kuo
%       date: 11/15/2019
%  copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%    purpose: localizer for twoPatchMotionDirfMRI

function myscreen = twoPatchMotionLocalizer(varargin)

% check arguments
getArgs(varargin,'stimulusType=dots');

% initalize the screen
myscreen = initScreen();
mglClearScreen(0.5);
% fix: set waitForBacktick if you want to synch with the scanner
% by waiting for the backtick key to be pressed before starting the experiment
% (for systems that use NI digital I/O, this will wait for the digital
% signal that the scanner has started collecting data)
task{1}.waitForBacktick = 1;

% task parameters
task{1}.segmin = [6 6];
task{1}.segmax = [6 6];

% initialize the task
for phaseNum = 1:length(task)
    [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback);
end

% global stimulus parameters
global stimulus;

% if scanning, do synch to volume
stimulus.scan = 1;
if stimulus.scan
    task{1}.synchToVol = [0 1];
    task{1}.segmin = task{1}.segmin - [0 0.5];
    task{1}.segmax = task{1}.segmax - [0 0.5];
end

% set coherence
stimulus.coherence = 1; % 0-1

% set left and right eccentricities of L and R patches
stimulus.leftEcc = -6;
stimulus.rightEcc = 6;

stimulus.contrast = 0.8;
stimulus.speed = 5; % in deg/s
stimulus.width = 10; % diameter

% init the stimulus
myscreen = initStimulus('stimulus',myscreen);
stimulus = initDots(stimulus,myscreen);

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
clear global stimulus

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;

if task.thistrial.thisseg == 1 % set left side
        
    % set the coherence (1st value per cell of stimulus.conditions)
    stimulus.dotsLeft = stimulus.dotsLeft.setCoherence(stimulus.dotsLeft,stimulus.coherence);
    
    % set the motion direction of left and right patches (2nd value per cell of stimulus.conditions = dir offset)
    dirI = randi(360);
    stimulus.dirL(task.trialnum) = dirI;
    stimulus.dotsLeft = stimulus.dotsLeft.setDir(stimulus.dotsLeft,stimulus.dirL(task.trialnum));
    
    % set the fixation color
    stimulus.fixColor = stimulus.normalFixColor;
    
elseif task.thistrial.thisseg == 2 % set right side
    
    % set the coherence (1st value per cell of stimulus.conditions)
    stimulus.dotsRight = stimulus.dotsRight.setCoherence(stimulus.dotsRight,stimulus.coherence);
    
    % set the motion direction of left and right patches (2nd value per cell of stimulus.conditions = dir offset)
    dirI = randi(360);
    stimulus.dirR(task.trialnum) = dirI;
    stimulus.dotsRight = stimulus.dotsRight.setDir(stimulus.dotsRight,stimulus.dirR(task.trialnum));
        
    % set the fixation color
    stimulus.fixColor = stimulus.normalFixColor;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus

% clear screen to gray
mglClearScreen(stimulus.backgroundColor);

if task.thistrial.thisseg == 1 % left side
    
    % update the dots
    stimulus.dotsLeft = stimulus.dotsLeft.update(stimulus.dotsLeft);
    
    % draw the dots
    stimulus.dotsLeft = stimulus.dotsLeft.draw(stimulus.dotsLeft);
    
elseif task.thistrial.thisseg == 2 % right side
    
    % update the dots
    stimulus.dotsRight = stimulus.dotsRight.update(stimulus.dotsRight);
    
    % draw the dots
    stimulus.dotsRight = stimulus.dotsRight.draw(stimulus.dotsRight);
    
end

% draw fixation cross
mglFixationCross(1,1,stimulus.fixColor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dots stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initDots(stimulus,myscreen)

width = strcat('width=',num2str(stimulus.width));
contrast = strcat('contrast=',num2str(stimulus.contrast));
speed = strcat('speed=',num2str(stimulus.speed));
leftEcc = strcat('xCenter=',num2str(stimulus.leftEcc));
rightEcc = strcat('xCenter=',num2str(stimulus.rightEcc));

% init the dot patches
stimulus.dotsLeft = dotsInit2Patch('framesPerSecond',myscreen.framesPerSecond,width,contrast,speed,leftEcc);
stimulus.dotsRight = dotsInit2Patch('framesPerSecond',myscreen.framesPerSecond,width,contrast,speed,rightEcc);

% set background color
stimulus.backgroundColor = 0.5;

% and fixation color
stimulus.normalFixColor = [1 1 1];
stimulus.correctFixColor = [0 1 0];
stimulus.incorrectFixColor = [1 0 0];
stimulus.neutralFixColor = [0 0 0];
