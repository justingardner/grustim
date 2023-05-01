% twoPatchMotionDirfMRI.m
%
%        $Id$
%      usage: twoPatchMotionDirfMRI
%         by: austin kuo
%       date: 11/14/2019
%  copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%    purpose: obtain fMRI responses to patches of varying coherence and motion direction differences
%
% NOTE: parameters you (might) want to change for this task are - 
% numTrials, coherence, dirs, scan (contrast, speed, width, leftEcc, rightEcc)

function myscreen = fourPatchWM(varargin)

% check arguments
getArgs(varargin,'stimulusType=dots');

% set default parameters
if ieNotDefined('atScanner'),atScanner = 0;end
if ieNotDefined('saveParam'),saveParam = 0;end
if ieNotDefined('screenParam')
    myscreen.displayName = 'fMRIproj_akuo2';
else
    myscreen.displayName = screenParam;
end

% initalize the screen
myscreen.background = 'gray';
myscreen.autoCloseScreen = 0;
myscreen.allowpause = 1;
myscreen.saveData = saveParam;
myscreen = initScreen(myscreen);

% fix: set waitForBacktick if you want to synch with the scanner
% by waiting for the backtick key to be pressed before starting the experiment
% (for systems that use NI digital I/O, this will wait for the digital
% signal that the scanner has started collecting data)
task{1}{1}.waitForBacktick = 1;

% task
task{1}{1}.getResponse = [0 0 0 0 0 0 1 0];
task{1}{1}.collectEyeData = true;
task{1}{1}.seglen = [1, 0.5, 0.5, 1, 0.5, 1, 1.5, 2];

% task parameters
task{1}{1}.parameter.cue = [1 2 3 4 5 6 7 8];
task{1}{1}.parameter.speed = [1 10];
task{1}{1}.random = 1;

% number of direction differences * number of coherences * number of sets desired (7*2*x)
task{1}{1}.numTrials = 32; % 1 set = 1 trial per (number of direction differences * number of coherences)

% initialize the task
for phaseNum = 1:length(task)
    [task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end

% global stimulus parameters
global stimulus;

% if scanning, do synch to volume
task{1}{1}.synchToVol = zeros(size(task{1}{1}.seglen));
if atScanner
    task{1}{1}.fudgeLastVolume = 1;
    task{1}{1}.seglen(end) = task{1}{1}.seglen(end)-0.5;
    task{1}{1}.synchToVol(end) = 1;
end

% set coherence
stimulus.coherence = 0; % 0-1

% set directions
stimulus.contrast = 1;
stimulus.speed = 5; % in deg/s
stimulus.width = 5; % diameter

stimulus.response = zeros(1,task{1}{1}.numTrials);
stimulus.correctResponse = zeros(1,task{1}{1}.numTrials);

% init the stimulus
myscreen = initStimulus('stimulus',myscreen);
stimulus = initDots(stimulus,myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

%mglSimulateRun(1,1000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task)) && ~myscreen.userHitEsc
    % update the task
    [task{1} myscreen phaseNum] = updateTask(task{1},myscreen,phaseNum);
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

if task.thistrial.thisseg == 1
    
    % set cue color to neutral
    stimulus.cuecolor{1} = [0 0 0];
    stimulus.cuecolor{2} = [0 0 0];
    stimulus.cuecolor{3} = [0 0 0];
    stimulus.cuecolor{4} = [0 0 0];

elseif task.thistrial.thisseg == 2
    
    % set cue color if this is a cued condition
    if task.thistrial.cue < 5
        stimulus.cuecolor{task.thistrial.cue} = [1 1 1];
    end

elseif task.thistrial.thisseg == 3
    
    % reset cue color to neutral
    stimulus.cuecolor{1} = [0 0 0];
    stimulus.cuecolor{2} = [0 0 0];
    stimulus.cuecolor{3} = [0 0 0];
    stimulus.cuecolor{4} = [0 0 0];

elseif task.thistrial.thisseg == 4

    % reset speed of dots
    stimulus.dotsUpLeft = stimulus.dotsUpLeft.setSpeed(stimulus.dotsUpLeft,5);
    stimulus.dotsUpRight = stimulus.dotsUpRight.setSpeed(stimulus.dotsUpRight,5);
    stimulus.dotsDownLeft = stimulus.dotsDownLeft.setSpeed(stimulus.dotsDownLeft,5);
    stimulus.dotsDownRight = stimulus.dotsDownRight.setSpeed(stimulus.dotsDownRight,5);

elseif task.thistrial.thisseg == 6

    % change correct patch to move at a different speed
    if task.thistrial.cue > 5
        cuenum = task.thistrial.cue - 4;
    else
        cuenum = task.thistrial.cue;
    end

    if cuenum == 1
        stimulus.dotsUpLeft = stimulus.dotsUpLeft.setSpeed(stimulus.dotsUpLeft,task.thistrial.speed);
    elseif cuenum == 2
        stimulus.dotsUpRight = stimulus.dotsUpRight.setSpeed(stimulus.dotsUpRight,task.thistrial.speed);
    elseif cuenum == 3
        stimulus.dotsDownLeft = stimulus.dotsDownLeft.setSpeed(stimulus.dotsDownLeft,task.thistrial.speed);
    elseif cuenum == 4
        stimulus.dotsDownRight = stimulus.dotsDownRight.setSpeed(stimulus.dotsDownRight,task.thistrial.speed);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus

% clear screen to gray
mglClearScreen(stimulus.backgroundColor);

if task.thistrial.thisseg == 1

    % draw the fixation cross in neutral colors
    mglFixationCrossDiag(1,1);

elseif task.thistrial.thisseg == 2
    
    % draw the cue
    mglFixationCrossDiag(1,1,stimulus.cuecolor);
    
elseif task.thistrial.thisseg == 3

    % draw the fixation cross in neutral colors
    mglFixationCrossDiag(1,1);
    
elseif task.thistrial.thisseg == 4

    % draw the fixation cross in neutral colors
    mglFixationCrossDiag(1,1);
    
    % update the dots
    stimulus.dotsUpLeft = stimulus.dotsUpLeft.update(stimulus.dotsUpLeft);
    stimulus.dotsUpRight = stimulus.dotsUpRight.update(stimulus.dotsUpRight);
    stimulus.dotsDownLeft = stimulus.dotsDownLeft.update(stimulus.dotsDownLeft);
    stimulus.dotsDownRight = stimulus.dotsDownRight.update(stimulus.dotsDownRight);
    
    % draw the dots
    stimulus.dotsUpLeft = stimulus.dotsUpLeft.draw(stimulus.dotsUpLeft);
    stimulus.dotsUpRight = stimulus.dotsUpRight.draw(stimulus.dotsUpRight);
    stimulus.dotsDownLeft = stimulus.dotsDownLeft.draw(stimulus.dotsDownLeft);
    stimulus.dotsDownRight = stimulus.dotsDownRight.draw(stimulus.dotsDownRight);

elseif task.thistrial.thisseg == 5

    % draw the fixation cross in neutral colors
    mglFixationCrossDiag(1,1);

elseif task.thistrial.thisseg == 6

    % draw the fixation cross in neutral colors
    mglFixationCrossDiag(1,1);

    % update the dots
    stimulus.dotsUpLeft = stimulus.dotsUpLeft.update(stimulus.dotsUpLeft);
    stimulus.dotsUpRight = stimulus.dotsUpRight.update(stimulus.dotsUpRight);
    stimulus.dotsDownLeft = stimulus.dotsDownLeft.update(stimulus.dotsDownLeft);
    stimulus.dotsDownRight = stimulus.dotsDownRight.update(stimulus.dotsDownRight);
    
    % draw the dots
    stimulus.dotsUpLeft = stimulus.dotsUpLeft.draw(stimulus.dotsUpLeft);
    stimulus.dotsUpRight = stimulus.dotsUpRight.draw(stimulus.dotsUpRight);
    stimulus.dotsDownLeft = stimulus.dotsDownLeft.draw(stimulus.dotsDownLeft);
    stimulus.dotsDownRight = stimulus.dotsDownRight.draw(stimulus.dotsDownRight);

elseif task.thistrial.thisseg == 7

    % draw the fixation cross in neutral colors
    mglFixationCrossDiag(1,1);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)

global stimulus

fprintf('Response received : %g', task.thistrial.whichButton);

% check the response
% if task.thistrial.gotResponse < 1
%     % see if it is correct
%     if isequal(task.thistrial.whichButton,task.thistrial.whichSide)
%         % report answer
%         fprintf(' !! Correct !!. Reaction time: %0.2f',task.thistrial.reactionTime);
%         % change fixation color
%         stimulus.fixColor = stimulus.correctFixColor;
%         stimulus.correctResponse(task.trialnum) = 1;
%         if task.thistrial.whichButton == 1
%             stimulus.response(task.trialnum) = 0; % subject answered left
%         elseif task.thistrial.whichButton == 2
%             stimulus.response(task.trialnum) = 1; % subject answered right
%         end
%         % task = jumpSegment(task);
%     elseif ~task.thistrial.whichSide
%         % report answer
%         fprintf(' -- Same direction --. Reaction time: %0.2f',task.thistrial.reactionTime);
%         stimulus.fixColor = stimulus.neutralFixColor;
%         stimulus.correctResponse(task.trialnum) = nan;
%         if task.thistrial.whichButton == 1
%             stimulus.response(task.trialnum) = 0; % subject answered left
%         elseif task.thistrial.whichButton == 2
%             stimulus.response(task.trialnum) = 1; % subject answered right
%         end
%         % task = jumpSegment(task);
%     else
%         % report answer
%         fprintf(' ++ Incorrect ++. Reaction time: %0.2f',task.thistrial.reactionTime);
%         % change fixation color
%         stimulus.fixColor = stimulus.incorrectFixColor;
%         stimulus.correctResponse(task.trialnum) = 0;
%         if task.thistrial.whichButton == 1
%             stimulus.response(task.trialnum) = 0; % subject answered left
%         elseif task.thistrial.whichButton == 2
%             stimulus.response(task.trialnum) = 1; % subject answered right
%         end
%         % task = jumpSegment(task);
%     end
%     
% % no response
% else 
%     stimulus.correctResponse(task.trialnum) = nan;
%     stimulus.response(task.trialnum) = nan;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dots stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initDots(stimulus,myscreen)

aperwidth = strcat('aperwidth=',num2str(stimulus.width));
contrast = strcat('contrast=',num2str(stimulus.contrast));
speed = strcat('speed=',num2str(stimulus.speed));
coherence = strcat('coherence=',num2str(stimulus.coherence));

% init the dot patches
stimulus.dotsUpLeft = dotsInitNew('framesPerSecond',myscreen.framesPerSecond,aperwidth,contrast,speed,coherence,'xCenter=-5/sqrt(2)','yCenter=5/sqrt(2)');
stimulus.dotsUpRight = dotsInitNew('framesPerSecond',myscreen.framesPerSecond,aperwidth,contrast,speed,coherence,'xCenter=5/sqrt(2)','yCenter=5/sqrt(2)');
stimulus.dotsDownLeft = dotsInitNew('framesPerSecond',myscreen.framesPerSecond,aperwidth,contrast,speed,coherence,'xCenter=-5/sqrt(2)','yCenter=-5/sqrt(2)');
stimulus.dotsDownRight = dotsInitNew('framesPerSecond',myscreen.framesPerSecond,aperwidth,contrast,speed,coherence,'xCenter=5/sqrt(2)','yCenter=-5/sqrt(2)');

% set background color
stimulus.backgroundColor = 0.5;

% and fixation color
stimulus.cueColor = [1 1 1];
stimulus.noncueColor = [0 0 0];
