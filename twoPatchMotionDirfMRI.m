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

function myscreen = twoPatchMotionDirfMRI(varargin)

% check arguments
getArgs(varargin,'stimulusType=dots');

% initalize the screen
myscreen = initScreen('fMRIprojFlex');
mglClearScreen(0.5);
% fix: set waitForBacktick if you want to synch with the scanner
% by waiting for the backtick key to be pressed before starting the experiment
% (for systems that use NI digital I/O, this will wait for the digital
% signal that the scanner has started collecting data)
task{1}.waitForBacktick = 1;

% task parameters
task{1}.segmin = [1 7];
task{1}.segmax = [1 7];
task{1}.getResponse = [0 1];
task{1}.synchToVol = [0 0];

% number of direction differences * number of coherences * number of sets desired (7*2*x)
task{1}.numTrials = 112; % 1 set = 1 trial per (number of direction differences * number of coherences)

% initialize the task
for phaseNum = 1:length(task)
    [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
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
stimulus.coherence = [0.1 0.9]; % 0-1

% set left and right eccentricities of L and R patches
stimulus.leftEcc = -6;
stimulus.rightEcc = 6;

% set directions
stimulus.dirs = (-12 : 4 : 12); % difference in L/R motion directions (in deg)

% create randomized loading vector so each trial loads a unique condition
orderedCondsPerSet = cell(length(stimulus.coherence),length(stimulus.dirs));
for i = 1:length(stimulus.coherence)
    for j = 1:length(stimulus.dirs)
        orderedCondsPerSet{i,j} = [stimulus.coherence(i) stimulus.dirs(j)];
    end
end
orderedCondsPerSet = reshape(orderedCondsPerSet',[],1); % turn matrix into vector; coh1(dir1 dir2 ... dirN) ... cohN(dir1 ... dirN)
stimulus.nSets = task{1}.numTrials/(length(stimulus.coherence)*length(stimulus.dirs));
stimulus.orderedCondsTotal = repmat(orderedCondsPerSet,[stimulus.nSets 1]);
stimulus.conditions = stimulus.orderedCondsTotal(randperm(length(stimulus.orderedCondsTotal))); % completely randomized conditions across nSets

stimulus.contrast = 0.8;
stimulus.speed = 5; % in deg/s
stimulus.width = 10; % diameter

stimulus.dirL = zeros(1,task{1}.numTrials);
stimulus.dirR = zeros(1,task{1}.numTrials);
stimulus.response = zeros(1,task{1}.numTrials);
stimulus.correctResponse = zeros(1,task{1}.numTrials);

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
    [task myscreen phaseNum] = updateTask(task,myscreen,phaseNum);
    % flip screen
    myscreen = tickScreen(myscreen,task);
end

% do quick wrap up calcs
% left stim - right stim direction
stimulus.LmR = stimulus.dirL-stimulus.dirR;

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;

if task.thistrial.thisseg == 1
        
    % set the coherence (1st value per cell of stimulus.conditions)
    stimulus.dotsLeft = stimulus.dotsLeft.setCoherence(stimulus.dotsLeft,stimulus.conditions{task.trialnum}(1));
    stimulus.dotsRight = stimulus.dotsRight.setCoherence(stimulus.dotsRight,stimulus.conditions{task.trialnum}(1));
    
    % set the motion direction of left and right patches (2nd value per cell of stimulus.conditions = dir offset)
    dirI = randi(360);
    stimulus.dirL(task.trialnum) = dirI;
    stimulus.dirR(task.trialnum) = dirI+stimulus.conditions{task.trialnum}(2);
    stimulus.dotsLeft = stimulus.dotsLeft.setDir(stimulus.dotsLeft,stimulus.dirL(task.trialnum));
    stimulus.dotsRight = stimulus.dotsRight.setDir(stimulus.dotsRight,stimulus.dirR(task.trialnum));
    
    % determine if left or right patch is more CW
    if stimulus.dirL(task.trialnum) < stimulus.dirR(task.trialnum)
        % left side is more CW
        task.thistrial.whichSide = 1;
    elseif stimulus.dirL(task.trialnum) > stimulus.dirR(task.trialnum)
        % right side is more CW
        task.thistrial.whichSide = 2;
    else
        % both are moving in the same direction
        task.thistrial.whichSide = 0;
    end
    
    % set the fixation color
    stimulus.fixColor = stimulus.stimulatingFixColor;
    
elseif task.thistrial.thisseg == 2

    stimulus.dotsLeft = stimulus.dotsLeft.setCoherence(stimulus.dotsLeft,0);
    stimulus.dotsRight = stimulus.dotsRight.setCoherence(stimulus.dotsRight,0);
    
    % set the fixation color to indicate response period
    stimulus.fixColor = stimulus.responseFixColor;
    
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
    
elseif task.thistrial.thisseg == 2
    
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

disp(sprintf('Response received : %g', task.thistrial.whichButton));

% check the response
if task.thistrial.gotResponse < 1
    % see if it is correct
    if isequal(task.thistrial.whichButton,task.thistrial.whichSide)
        % report answer
        disp(sprintf(' !! Correct !!. Reaction time: %0.2f',task.thistrial.reactionTime));
        % change fixation color
        stimulus.fixColor = stimulus.correctFixColor;
        stimulus.correctResponse(task.trialnum) = 1;
        if task.thistrial.whichButton == 1
            stimulus.response(task.trialnum) = 0; % subject answered left
        elseif task.thistrial.whichButton == 2
            stimulus.response(task.trialnum) = 1; % subject answered right
        end
        % task = jumpSegment(task);
    elseif ~task.thistrial.whichSide
        % report answer
        disp(sprintf(' -- Same direction --. Reaction time: %0.2f',task.thistrial.reactionTime));
        stimulus.fixColor = stimulus.neutralFixColor;
        stimulus.correctResponse(task.trialnum) = nan;
        if task.thistrial.whichButton == 1
            stimulus.response(task.trialnum) = 0; % subject answered left
        elseif task.thistrial.whichButton == 2
            stimulus.response(task.trialnum) = 1; % subject answered right
        end
        % task = jumpSegment(task);
    else
        % report answer
        disp(sprintf(' ++ Incorrect ++. Reaction time: %0.2f',task.thistrial.reactionTime));
        % change fixation color
        stimulus.fixColor = stimulus.incorrectFixColor;
        stimulus.correctResponse(task.trialnum) = 0;
        if task.thistrial.whichButton == 1
            stimulus.response(task.trialnum) = 0; % subject answered left
        elseif task.thistrial.whichButton == 2
            stimulus.response(task.trialnum) = 1; % subject answered right
        end
        % task = jumpSegment(task);
    end
    
% no response
else 
    stimulus.correctResponse(task.trialnum) = nan;
    stimulus.response(task.trialnum) = nan;
end

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
stimulus.stimulatingFixColor = [0 0 0];
stimulus.responseFixColor = [1 1 1];
stimulus.correctFixColor = [0 1 0];
stimulus.incorrectFixColor = [1 0 0];
stimulus.neutralFixColor = [0 0 1];
