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

function myscreen = fourPatchAtt(varargin)

% check arguments
getArgs(varargin, [], 'verbose=0');

% set default parameters
if ieNotDefined('atScanner'),atScanner = 0;end
if ieNotDefined('saveParam'),saveParam = 0;end
if ieNotDefined('screenParam')
    myscreen.displayName = 'fMRIproj_akuo2';
else
    myscreen.displayName = screenParam;
end
if ieNotDefined('initStair'),initStair = 0;end

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
task{1}{1}.collectEyeData = false;
task{1}{1}.segmin = [0.5, 1, 0.5, 1, 0.5, 1, 2, 1.5];
task{1}{1}.segmax = [0.5, 1, 0.5, 1, 0.5, 1, 2, 6.5];

% task parameters
task{1}{1}.parameter.cue = [1 2 3 4 5 6];
task{1}{1}.parameter.change = [0 1];
task{1}{1}.randVars.uniform.speed1 = [2 4 6 8 10];
task{1}{1}.randVars.uniform.speed2 = [2 4 6 8 10];
task{1}{1}.randVars.uniform.speed3 = [2 4 6 8 10];
task{1}{1}.randVars.uniform.speed4 = [2 4 6 8 10];
task{1}{1}.randVars.calculated.trialType = nan;
task{1}{1}.randVars.calculated.cuenum = nan;
task{1}{1}.randVars.calculated.kWeb1 = nan;
task{1}{1}.randVars.calculated.kWeb2 = nan;
task{1}{1}.random = 1;

task{1}{1}.numTrials = 24; % 12 conditions

% initialize the task
for phaseNum = 1:length(task)
    [task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end

% global stimulus parameters
global stimulus;

% if scanning, do synch to volume
task{1}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
if atScanner
    task{1}{1}.fudgeLastVolume = 1;
    task{1}{1}.synchToVol(end) = 1;
end

% set coherence
stimulus.coherence = 0; % 0-1

% set directions
stimulus.contrast = 1;
stimulus.initSpeed = 5; % in deg/s
stimulus.width = 5; % diameter
stimulus.initStair = initStair;

stimulus.response = zeros(1,task{1}{1}.numTrials);
stimulus.correctResponse = zeros(1,task{1}{1}.numTrials);

% init the stimulus
myscreen = initStimulus('stimulus',myscreen);
stimulus = initDots(stimulus,myscreen);

% init the staircase
if stimulus.initStair
    fprintf('\n(fourPatchAtt) Initializing staircases from scratch...\n\n');
    stimulus = initStaircase(stimulus);
else
    fprintf('\n(fourPatchAtt) Re-using staircase from previous run...\n\n');
end

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
    if task.thistrial.cue < 4
        stimulus.cuecolor{task.thistrial.cue} = [1 1 1]; 
        task.thistrial.trialType = 1;
    else
        task.thistrial.trialType = 2;
    end

elseif task.thistrial.thisseg == 3
    
    % reset cue color to neutral
    stimulus.cuecolor{1} = [0 0 0];
    stimulus.cuecolor{2} = [0 0 0];
    stimulus.cuecolor{3} = [0 0 0];
    stimulus.cuecolor{4} = [0 0 0];

elseif task.thistrial.thisseg == 4

    % reset speed of dots
    stimulus.dotsUpLeft = stimulus.dotsUpLeft.setSpeed(stimulus.dotsUpLeft,task.thistrial.speed1);
    stimulus.dotsUpRight = stimulus.dotsUpRight.setSpeed(stimulus.dotsUpRight,task.thistrial.speed2);
    stimulus.dotsDownLeft = stimulus.dotsDownLeft.setSpeed(stimulus.dotsDownLeft,task.thistrial.speed3);
    stimulus.dotsDownRight = stimulus.dotsDownRight.setSpeed(stimulus.dotsDownRight,task.thistrial.speed4);

elseif task.thistrial.thisseg == 6

    % change correct patch to move at a different speed
    if task.thistrial.trialType == 2
        task.thistrial.cuenum = task.thistrial.cue - randi(4,1);
    else
        task.thistrial.cuenum = task.thistrial.cue;
    end

    % Get the new delta for this trial from the staircase
    if task.thistrial.trialType == 1
        [kWeb, stimulus.staircase1] = doStaircase('testValue',stimulus.staircase1);
        task.thistrial.kWeb1 = kWeb;
        fprintf('Current kWeb %d value: %0.2f \n\n', task.thistrial.trialType, kWeb)
    elseif task.thistrial.trialType == 2
        [kWeb, stimulus.staircase2] = doStaircase('testValue',stimulus.staircase2);
        task.thistrial.kWeb2 = kWeb;
        fprintf('Current kWeb %d value: %0.2f \n\n', task.thistrial.trialType, kWeb)
    end
    
    % adjust the speed of the dots by a value determined by Weber's law
    if task.thistrial.cuenum == 1 && task.thistrial.change == 1
        stimulus.dotsUpLeft = stimulus.dotsUpLeft.setSpeed(stimulus.dotsUpLeft,task.thistrial.speed1 + task.thistrial.speed1*kWeb);
    elseif task.thistrial.cuenum == 2 && task.thistrial.change == 1
        stimulus.dotsUpRight = stimulus.dotsUpRight.setSpeed(stimulus.dotsUpRight,task.thistrial.speed2 + task.thistrial.speed2*kWeb);
    elseif task.thistrial.cuenum == 3 && task.thistrial.change == 1
        stimulus.dotsDownLeft = stimulus.dotsDownLeft.setSpeed(stimulus.dotsDownLeft,task.thistrial.speed3 + task.thistrial.speed3*kWeb);
    elseif task.thistrial.cuenum == 4 && task.thistrial.change == 1
        stimulus.dotsDownRight = stimulus.dotsDownRight.setSpeed(stimulus.dotsDownRight,task.thistrial.speed4 + task.thistrial.speed4*kWeb);
    end
    
elseif task.thistrial.thisseg == 7

    % set cuecolor to indicate response period
    if task.thistrial.cuenum == 1
        stimulus.cuecolor{1} = [1 1 0];
    elseif task.thistrial.cuenum == 2
        stimulus.cuecolor{2} = [1 1 0];
    elseif task.thistrial.cuenum == 3
        stimulus.cuecolor{3} = [1 1 0];
    elseif task.thistrial.cuenum == 4
        stimulus.cuecolor{4} = [1 1 0];
    end

elseif task.thistrial.thisseg == 8

    % reset cue color to neutral
    stimulus.cuecolor{1} = [0 0 0];
    stimulus.cuecolor{2} = [0 0 0];
    stimulus.cuecolor{3} = [0 0 0];
    stimulus.cuecolor{4} = [0 0 0];

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

    % draw the fixation cross in respond colors
    mglFixationCrossDiag(1,1,stimulus.cuecolor);

elseif task.thistrial.thisseg == 8

    % draw the fixation cross in neutral colors
    mglFixationCrossDiag(1,1,stimulus.cuecolor);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)

global stimulus

% check the response
if task.thistrial.gotResponse < 1

    fprintf('Response received : %g\n', task.thistrial.whichButton);

    if ~any(task.thistrial.whichButton == [1 2]) % [1 2] should (hopefully) correspond to left/right buttons
        error('Check your button inputs')
    end
    
    % determine what the correct button should be
    if task.thistrial.change
        correctbutton = 1;
    elseif ~task.thistrial.change
        correctbutton = 2;
    end

    % see if it is correct
    if isequal(task.thistrial.whichButton,correctbutton) % correct
        corr = 1;
        % report answer
        fprintf('Trial %d. !! Correct !!. \n',task.trialnum);
        task.thistrial.correct = corr;

        % update staircases
        if task.thistrial.trialType == 1
            stimulus.staircase1 = doStaircase('update',stimulus.staircase1,corr);
        elseif task.thistrial.trialType == 2
            stimulus.staircase2 = doStaircase('update',stimulus.staircase2,corr);
        end
        
    else % incorrect
        corr = 0;
        % report answer
        fprintf('Trial %d. ++ Incorrect ++. \n',task.trialnum);
        task.thistrial.correct = corr;
        
        % update staircases
        if task.thistrial.trialType == 1
            stimulus.staircase1 = doStaircase('update',stimulus.staircase1,corr);
        elseif task.thistrial.trialType == 2
            stimulus.staircase2 = doStaircase('update',stimulus.staircase2,corr);
        end

    end

end

if task.thistrial.trialType == 1
    [kWeb, stimulus.staircase1] = doStaircase('testValue',stimulus.staircase1);
    fprintf('Previous kWeb %d value: %0.2f \n', task.thistrial.trialType, task.thistrial.kWeb1)
    fprintf('Next kWeb %d value: %0.2f \n\n', task.thistrial.trialType, kWeb)
else
    [kWeb, stimulus.staircase2] = doStaircase('testValue',stimulus.staircase2);
    fprintf('Previous kWeb %d value: %0.2f \n', task.thistrial.trialType, task.thistrial.kWeb2)
    fprintf('Next kWeb %d value: %0.2f \n\n', task.thistrial.trialType, kWeb)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dots stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initDots(stimulus,myscreen)

aperwidth = strcat('aperwidth=',num2str(stimulus.width));
contrast = strcat('contrast=',num2str(stimulus.contrast));
speed = strcat('speed=',num2str(stimulus.initSpeed));
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the staircase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initStaircase(stimulus)

stimulus.staircase1 = doStaircase('init','upDown','nup=1','ndown=1','initialThreshold=5','initialStepsize=0.05','nTrials=24','stepRule=levitt','maxStepsize=1','minStepsize=.05');
stimulus.staircase2 = doStaircase('init','upDown','nup=1','ndown=1','initialThreshold=5','initialStepsize=0.05','nTrials=24','stepRule=levitt','maxStepsize=1','minStepsize=.05');



