% twoPatchMotionDir.m
%
%        $Id$
%      usage: twoPatchMotionDir
%         by: austin kuo
%       date: 10/7/2019
%  copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%    purpose: motion direction discrimination between 2 L/R patches
%
function myscreen = twoPatchMotionDir(varargin)

% check arguments
getArgs(varargin,'stimulusType=dots');

% initalize the screen
myscreen = initScreen;

% fix: set waitForBacktick if you want to synch with the scanner
% by waiting for the backtick key to be pressed before starting the experiment
% (for systems that use NI digital I/O, this will wait for the digital
% signal that the scanner has started collecting data)
task{1}.waitForBacktick = 1;

% task parameters
task{1}.segmin = [0.5 1 inf];
task{1}.segmax = [0.5 1 inf];
task{1}.getResponse = [0 0 1];
task{1}.randVars.uniform.whichSide = [1 2];
task{1}.parameter.eccentricity = [6 6];

% number of direction differences * number of sets desired * number of coherences (7*6*5)
task{1}.numTrials = 210; % 1 set = 1 trial per (number of direction differences * number of coherences)


% initialize the task
for phaseNum = 1:length(task)
    [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end

% global stimulus parameters
global stimulus;
stimulus.contrast = 0.5;
stimulus.dots.speed = 5; % in deg/s
stimulus.dots.width = 10; % diameter
stimulus.dots.dirs = (-15 : 5 : 15); % difference in L/R motion directions (in deg)
stimulus.nConditions = length(stimulus.dots.dirs);
stimulus.dots.dirL = zeros(task{1}.numTrials/stimulus.nConditions,stimulus.nConditions);
stimulus.dots.dirR = zeros(task{1}.numTrials/stimulus.nConditions,stimulus.nConditions);
stimulus.dots.correctResponse = zeros(task{1}.numTrials/stimulus.nConditions,stimulus.nConditions);
stimulus.dots.correctResponse = zeros(task{1}.numTrials/stimulus.nConditions,stimulus.nConditions);
stimulus.setNumber = 0;
stimulus.trialNumber = 0;

% set coherence
stimulus.dots.coherence = [.05 .30 .55 .70 .95]; % 0-1
coherenceUnordered = repmat(stimulus.dots.coherence,1,task{1}.numTrials/length(stimulus.dots.coherence));
stimulus.dots.coherenceOrder = coherenceUnordered(randperm(length(coherenceUnordered)));

% init the stimulus
myscreen = initStimulus('stimulus',myscreen);

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

% do quick wrap up calcs
% left stim - right stim direction
stimulus.dots.LmR = stimulus.dots.dirL-stimulus.dots.dirR;

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;

if task.thistrial.thisseg == 1
    
    % shuffle the order of direction difference every stimulus.nConditions trials
    stimulus.nTrialMod = mod(task.trialnum,stimulus.nConditions);
    if stimulus.nTrialMod == 1
        stimulus.dots.dirs = stimulus.dots.dirs(randperm(length(stimulus.dots.dirs)));
        stimulus.dots.dirOffset = stimulus.dots.dirs(1);
        stimulus.setNumber = stimulus.setNumber+1;
    elseif stimulus.nTrialMod == 0
        stimulus.dots.dirOffset = stimulus.dots.dirs(stimulus.nConditions);
        stimulus.nTrialMod = stimulus.nConditions;
    else
        stimulus.dots.dirOffset = stimulus.dots.dirs(stimulus.nTrialMod);
    end
    
    % get contrast for left and right
    if task.thistrial.whichSide == 1
        leftContrast = stimulus.contrast;
        rightContrast = stimulus.contrast;
    else
        leftContrast = stimulus.contrast;
        rightContrast = stimulus.contrast;
    end
    
    % dots
    stimulus.trialNumber = stimulus.trialNumber+1;
    stimulus = initDots(stimulus,myscreen);
    
    % determine if left or right patch is more CW
    if stimulus.dots.dirL(stimulus.setNumber,stimulus.nTrialMod) < stimulus.dots.dirR(stimulus.setNumber,stimulus.nTrialMod)
        % left side is more CW
        task.thistrial.whichSide = 1;
    elseif stimulus.dots.dirL(stimulus.setNumber,stimulus.nTrialMod) > stimulus.dots.dirR(stimulus.setNumber,stimulus.nTrialMod)
        % right side is more CW
        task.thistrial.whichSide = 2;
    else
        % both are moving in the same direction
        task.thistrial.whichSide = 0;
    end
    
    % set the contrast
    stimulus.dotsLeft = stimulus.dotsLeft.setContrast(stimulus.dotsLeft,leftContrast);
    stimulus.dotsRight = stimulus.dotsRight.setContrast(stimulus.dotsRight,rightContrast);
    
    % set the position
    stimulus.dotsLeft = stimulus.dotsLeft.setCenter(stimulus.dotsLeft,-task.thistrial.eccentricity,0);
    stimulus.dotsRight = stimulus.dotsLeft.setCenter(stimulus.dotsRight,task.thistrial.eccentricity,0);
    
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
        stimulus.fixColor = stimulus.correctFixColor;
        stimulus.dots.correctResponse(stimulus.setNumber,stimulus.nTrialMod) = 1;
        if task.thistrial.whichButton == 1
            stimulus.dots.response(stimulus.setNumber,stimulus.nTrialMod) = 0; % subject answered left
        elseif task.thistrial.whichButton == 2
            stimulus.dots.response(stimulus.setNumber,stimulus.nTrialMod) = 1; % subject answered right
        end
        task = jumpSegment(task);
    elseif ~task.thistrial.whichSide
        % report answer
        disp(sprintf(' -- Same direction --. Reaction time: %0.2f',task.thistrial.reactionTime));
        stimulus.fixColor = stimulus.neutralFixColor;
        if task.thistrial.whichButton == 1
            stimulus.dots.response(stimulus.setNumber,stimulus.nTrialMod) = 0; % subject answered left
        elseif task.thistrial.whichButton == 2
            stimulus.dots.response(stimulus.setNumber,stimulus.nTrialMod) = 1; % subject answered right
        end
        task = jumpSegment(task);
    else
        % report answer
        disp(sprintf(' ++ Incorrect ++. Reaction time: %0.2f',task.thistrial.reactionTime));
        % change fixation color
        stimulus.fixColor = stimulus.incorrectFixColor;
        stimulus.dots.correctResponse(stimulus.setNumber,stimulus.nTrialMod) = 0;
        if task.thistrial.whichButton == 1
            stimulus.dots.response(stimulus.setNumber,stimulus.nTrialMod) = 0; % subject answered left
        elseif task.thistrial.whichButton == 2
            stimulus.dots.response(stimulus.setNumber,stimulus.nTrialMod) = 1; % subject answered right
        end
        task = jumpSegment(task);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dots stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initDots(stimulus,myscreen)

coherence = strcat('coherence=',num2str(stimulus.dots.coherenceOrder(stimulus.trialNumber)));
width = strcat('width=',num2str(stimulus.dots.width));
dirI = randi(360);
stimulus.dots.dirL(stimulus.setNumber,stimulus.nTrialMod) = dirI;
stimulus.dots.dirR(stimulus.setNumber,stimulus.nTrialMod) = dirI+stimulus.dots.dirOffset;
dirL = strcat('dir=',num2str(dirI));
dirR = strcat('dir=',num2str(dirI+stimulus.dots.dirOffset));
speed = strcat('speed=',num2str(stimulus.dots.speed));

% init the dot patches
stimulus.dotsLeft = dotsInit('framesPerSecond',myscreen.framesPerSecond,dirL,coherence,width,speed);
stimulus.dotsRight = dotsInit('framesPerSecond',myscreen.framesPerSecond,dirR,coherence,width,speed);

% set background color
stimulus.backgroundColor = 0.5;

% and fixation color
stimulus.normalFixColor = [1 1 1];
stimulus.correctFixColor = [0 1 0];
stimulus.incorrectFixColor = [1 0 0];
stimulus.neutralFixColor = [0 0 0];


