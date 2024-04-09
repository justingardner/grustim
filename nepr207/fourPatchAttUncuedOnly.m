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

function myscreen = fourPatchAttUncuedOnly(varargin)

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
if ieNotDefined('initStair'),initStair = 1;end

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
task{1}{1}.segmin = [0.5, 1, 0.5, 1, 0.5, 1, 2, 1.5]; % average time = 10.5s
task{1}{1}.segmax = [0.5, 1, 0.5, 1, 0.5, 1, 2, 6.5];

% task parameters
task{1}{1}.parameter.cue = 5;
task{1}{1}.randVars.uniform.faster1 = [1 2];
task{1}{1}.randVars.uniform.faster2 = [1 2];
task{1}{1}.randVars.uniform.faster3 = [1 2];
task{1}{1}.randVars.uniform.faster4 = [1 2];
task{1}{1}.randVars.uniform.speed1 = [2 4 6 8 10];
task{1}{1}.randVars.uniform.speed2 = [2 4 6 8 10];
task{1}{1}.randVars.uniform.speed3 = [2 4 6 8 10];
task{1}{1}.randVars.uniform.speed4 = [2 4 6 8 10];
task{1}{1}.randVars.calculated.trialType = nan;
task{1}{1}.randVars.calculated.cuenum = nan;
task{1}{1}.randVars.calculated.kWeb = nan;
task{1}{1}.random = 1;

task{1}{1}.numTrials = inf; % 12 conditions

% initialize the task
for phaseNum = 1:length(task)
    [task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end

% global stimulus parameters
global stimulusUncued;

% if scanning, do synch to volume
task{1}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
if atScanner
    task{1}{1}.fudgeLastVolume = 1;
    task{1}{1}.synchToVol(end) = 1;
end

% set coherence
stimulusUncued.coherence = 0; % 0-1

% set directions
stimulusUncued.contrast = 1;
stimulusUncued.initSpeed = 5; % in deg/s
stimulusUncued.width = 5; % diameter
stimulusUncued.initStair = initStair;

% init the stimulus
myscreen = initStimulus('stimulusUncued',myscreen);
stimulusUncued = initDots(stimulusUncued,myscreen);

% init the staircase
if stimulusUncued.initStair
    fprintf('\n(fourPatchAttCuedOnly) Initializing staircases from scratch...\n\n');
    stimulusUncued = initStaircase(stimulusUncued);
else
    fprintf('\n(fourPatchAttCuedOnly) Re-using staircase from previous run...\n\n');
end

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

global stimulusUncued;

if task.thistrial.thisseg == 1
    
    % set cue color to neutral
    stimulusUncued.cuecolor{1} = [0 0 0];
    stimulusUncued.cuecolor{2} = [0 0 0];
    stimulusUncued.cuecolor{3} = [0 0 0];
    stimulusUncued.cuecolor{4} = [0 0 0];

elseif task.thistrial.thisseg == 2
    
    stimulusUncued.cuecolor{1} = [1 1 1];
    stimulusUncued.cuecolor{2} = [1 1 1];
    stimulusUncued.cuecolor{3} = [1 1 1];
    stimulusUncued.cuecolor{4} = [1 1 1];

elseif task.thistrial.thisseg == 3
    
    % reset cue color to neutral
    stimulusUncued.cuecolor{1} = [0 0 0];
    stimulusUncued.cuecolor{2} = [0 0 0];
    stimulusUncued.cuecolor{3} = [0 0 0];
    stimulusUncued.cuecolor{4} = [0 0 0];

elseif task.thistrial.thisseg == 4

    % change correct patch to move at a different speed
    task.thistrial.cuenum = randi(4,1);

    % Get the new delta for this trial from the staircase
    [kWeb, stimulusUncued.staircase] = doStaircase('testValue',stimulusUncued.staircase);
    task.thistrial.kWeb = kWeb;
    fprintf('Current kWeb distributed value: %0.2f \n\n', kWeb)

    % set dot base speed
    stimulusUncued.dotsUpLeft = stimulusUncued.dotsUpLeft.setSpeed(stimulusUncued.dotsUpLeft,task.thistrial.speed1);
    stimulusUncued.dotsUpRight = stimulusUncued.dotsUpRight.setSpeed(stimulusUncued.dotsUpRight,task.thistrial.speed2);
    stimulusUncued.dotsDownLeft = stimulusUncued.dotsDownLeft.setSpeed(stimulusUncued.dotsDownLeft,task.thistrial.speed3);
    stimulusUncued.dotsDownRight = stimulusUncued.dotsDownRight.setSpeed(stimulusUncued.dotsDownRight,task.thistrial.speed4);

    % adjust the speed of the dots by a value determined by Weber's law
    if task.thistrial.faster1 == 1
        stimulusUncued.dotsUpLeft = stimulusUncued.dotsUpLeft.setSpeed(stimulusUncued.dotsUpLeft,task.thistrial.speed1 + task.thistrial.speed1*kWeb);
    end
    if task.thistrial.faster2 == 1
        stimulusUncued.dotsUpRight = stimulusUncued.dotsUpRight.setSpeed(stimulusUncued.dotsUpRight,task.thistrial.speed2 + task.thistrial.speed2*kWeb);
    end
    if task.thistrial.faster3 == 1
        stimulusUncued.dotsDownLeft = stimulusUncued.dotsDownLeft.setSpeed(stimulusUncued.dotsDownLeft,task.thistrial.speed3 + task.thistrial.speed3*kWeb);
    end
    if task.thistrial.faster4 == 1
        stimulusUncued.dotsDownRight = stimulusUncued.dotsDownRight.setSpeed(stimulusUncued.dotsDownRight,task.thistrial.speed4 + task.thistrial.speed4*kWeb);
    end

elseif task.thistrial.thisseg == 6

    % Get the new delta for this trial from the staircase
    [kWeb, stimulusUncued.staircase] = doStaircase('testValue',stimulusUncued.staircase);
    task.thistrial.kWeb = kWeb;

    % set dot base speed
    stimulusUncued.dotsUpLeft = stimulusUncued.dotsUpLeft.setSpeed(stimulusUncued.dotsUpLeft,task.thistrial.speed1);
    stimulusUncued.dotsUpRight = stimulusUncued.dotsUpRight.setSpeed(stimulusUncued.dotsUpRight,task.thistrial.speed2);
    stimulusUncued.dotsDownLeft = stimulusUncued.dotsDownLeft.setSpeed(stimulusUncued.dotsDownLeft,task.thistrial.speed3);
    stimulusUncued.dotsDownRight = stimulusUncued.dotsDownRight.setSpeed(stimulusUncued.dotsDownRight,task.thistrial.speed4);
    
    % adjust the speed of the dots by a value determined by Weber's law
    if task.thistrial.faster1 == 2
        stimulusUncued.dotsUpLeft = stimulusUncued.dotsUpLeft.setSpeed(stimulusUncued.dotsUpLeft,task.thistrial.speed1 + task.thistrial.speed1*kWeb);
    end
    if task.thistrial.faster2 == 2
        stimulusUncued.dotsUpRight = stimulusUncued.dotsUpRight.setSpeed(stimulusUncued.dotsUpRight,task.thistrial.speed2 + task.thistrial.speed2*kWeb);
    end
    if task.thistrial.faster3 == 2
        stimulusUncued.dotsDownLeft = stimulusUncued.dotsDownLeft.setSpeed(stimulusUncued.dotsDownLeft,task.thistrial.speed3 + task.thistrial.speed3*kWeb);
    end
    if task.thistrial.faster4 == 2
        stimulusUncued.dotsDownRight = stimulusUncued.dotsDownRight.setSpeed(stimulusUncued.dotsDownRight,task.thistrial.speed4 + task.thistrial.speed4*kWeb);
    end
    
elseif task.thistrial.thisseg == 7

    % set cuecolor to indicate response period
    if task.thistrial.cuenum == 1
        stimulusUncued.cuecolor{1} = [1 1 0];
    elseif task.thistrial.cuenum == 2
        stimulusUncued.cuecolor{2} = [1 1 0];
    elseif task.thistrial.cuenum == 3
        stimulusUncued.cuecolor{3} = [1 1 0];
    elseif task.thistrial.cuenum == 4
        stimulusUncued.cuecolor{4} = [1 1 0];
    end

elseif task.thistrial.thisseg == 8

    % reset cue color to neutral
    stimulusUncued.cuecolor{1} = [0 0 0];
    stimulusUncued.cuecolor{2} = [0 0 0];
    stimulusUncued.cuecolor{3} = [0 0 0];
    stimulusUncued.cuecolor{4} = [0 0 0];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulusUncued

% clear screen to gray
mglClearScreen(stimulusUncued.backgroundColor);

if task.thistrial.thisseg == 1

    % draw the fixation cross in neutral colors
    mglFixationCrossDiag(1,1);

elseif task.thistrial.thisseg == 2
    
    % draw the cue
    mglFixationCrossDiag(1,1,stimulusUncued.cuecolor);
    
elseif task.thistrial.thisseg == 3

    % draw the fixation cross in neutral colors
    mglFixationCrossDiag(1,1);
    
elseif task.thistrial.thisseg == 4

    % draw the fixation cross in neutral colors
    mglFixationCrossDiag(1,1);
    
    % update the dots
    stimulusUncued.dotsUpLeft = stimulusUncued.dotsUpLeft.update(stimulusUncued.dotsUpLeft);
    stimulusUncued.dotsUpRight = stimulusUncued.dotsUpRight.update(stimulusUncued.dotsUpRight);
    stimulusUncued.dotsDownLeft = stimulusUncued.dotsDownLeft.update(stimulusUncued.dotsDownLeft);
    stimulusUncued.dotsDownRight = stimulusUncued.dotsDownRight.update(stimulusUncued.dotsDownRight);
    
    % draw the dots
    stimulusUncued.dotsUpLeft = stimulusUncued.dotsUpLeft.draw(stimulusUncued.dotsUpLeft);
    stimulusUncued.dotsUpRight = stimulusUncued.dotsUpRight.draw(stimulusUncued.dotsUpRight);
    stimulusUncued.dotsDownLeft = stimulusUncued.dotsDownLeft.draw(stimulusUncued.dotsDownLeft);
    stimulusUncued.dotsDownRight = stimulusUncued.dotsDownRight.draw(stimulusUncued.dotsDownRight);

elseif task.thistrial.thisseg == 5

    % draw the fixation cross in neutral colors
    mglFixationCrossDiag(1,1);

elseif task.thistrial.thisseg == 6

    % draw the fixation cross in neutral colors
    mglFixationCrossDiag(1,1);

    % update the dots
    stimulusUncued.dotsUpLeft = stimulusUncued.dotsUpLeft.update(stimulusUncued.dotsUpLeft);
    stimulusUncued.dotsUpRight = stimulusUncued.dotsUpRight.update(stimulusUncued.dotsUpRight);
    stimulusUncued.dotsDownLeft = stimulusUncued.dotsDownLeft.update(stimulusUncued.dotsDownLeft);
    stimulusUncued.dotsDownRight = stimulusUncued.dotsDownRight.update(stimulusUncued.dotsDownRight);
    
    % draw the dots
    stimulusUncued.dotsUpLeft = stimulusUncued.dotsUpLeft.draw(stimulusUncued.dotsUpLeft);
    stimulusUncued.dotsUpRight = stimulusUncued.dotsUpRight.draw(stimulusUncued.dotsUpRight);
    stimulusUncued.dotsDownLeft = stimulusUncued.dotsDownLeft.draw(stimulusUncued.dotsDownLeft);
    stimulusUncued.dotsDownRight = stimulusUncued.dotsDownRight.draw(stimulusUncued.dotsDownRight);

elseif task.thistrial.thisseg == 7

    % draw the fixation cross in respond colors
    mglFixationCrossDiag(1,1,stimulusUncued.cuecolor);

elseif task.thistrial.thisseg == 8

    % draw the fixation cross in neutral colors
    mglFixationCrossDiag(1,1,stimulusUncued.cuecolor);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)

global stimulusUncued

% check the response
if task.thistrial.gotResponse < 1

    fprintf('Response received : %g\n', task.thistrial.whichButton);

    if ~any(task.thistrial.whichButton == [1 2]) % [1 2] should (hopefully) correspond to left/right buttons
        error('Check your button inputs')
    end
    
    % determine what the correct button should be
    if task.thistrial.cuenum == 1
        if task.thistrial.faster1 == 1
            correctbutton = 1;
        elseif task.thistrial.faster1 == 2
            correctbutton = 2;
        end
    elseif task.thistrial.cuenum == 2
        if task.thistrial.faster2 == 1
            correctbutton = 1;
        elseif task.thistrial.faster2 == 2
            correctbutton = 2;
        end
    elseif task.thistrial.cuenum == 3
        if task.thistrial.faster3 == 1
            correctbutton = 1;
        elseif task.thistrial.faster3 == 2
            correctbutton = 2;
        end
    elseif task.thistrial.cuenum == 4
        if task.thistrial.faster4 == 1
            correctbutton = 1;
        elseif task.thistrial.faster4 == 2
            correctbutton = 2;
        end
    end

    % see if it is correct
    if isequal(task.thistrial.whichButton,correctbutton) % correct
        corr = 1;
        % report answer
        fprintf('Trial %d. !! Correct !!. \n',task.trialnum);
        task.thistrial.correct = corr;

        % update staircases
        stimulusUncued.staircase = doStaircase('update',stimulusUncued.staircase,corr);

        stimulusUncued.cuecolor{1} = [0 1 0];
        stimulusUncued.cuecolor{2} = [0 1 0];
        stimulusUncued.cuecolor{3} = [0 1 0];
        stimulusUncued.cuecolor{4} = [0 1 0];
        
    else % incorrect
        corr = 0;
        % report answer
        fprintf('Trial %d. ++ Incorrect ++. \n',task.trialnum);
        task.thistrial.correct = corr;
        
        % update staircases
        stimulusUncued.staircase = doStaircase('update',stimulusUncued.staircase,corr);

        stimulusUncued.cuecolor{1} = [1 0 0];
        stimulusUncued.cuecolor{2} = [1 0 0];
        stimulusUncued.cuecolor{3} = [1 0 0];
        stimulusUncued.cuecolor{4} = [1 0 0];

    end

    [kWeb, stimulusUncued.staircase] = doStaircase('testValue',stimulusUncued.staircase);
    fprintf('Previous kWeb distributed value: %0.2f \n', task.thistrial.kWeb)
    fprintf('Next kWeb distributed value: %0.2f \n\n',kWeb)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dots stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulusUncued = initDots(stimulusUncued,myscreen)

aperwidth = strcat('aperwidth=',num2str(stimulusUncued.width));
contrast = strcat('contrast=',num2str(stimulusUncued.contrast));
speed = strcat('speed=',num2str(stimulusUncued.initSpeed));
coherence = strcat('coherence=',num2str(stimulusUncued.coherence));

% init the dot patches
[~,name] = system('hostname');
if any(strfind(name,'oban'))
    stimulusUncued.dotsUpLeft = dotsInitNew('framesPerSecond',myscreen.framesPerSecond,aperwidth,contrast,speed,coherence,'xCenter=-5/sqrt(2)','yCenter=5/sqrt(2)');
    stimulusUncued.dotsUpRight = dotsInitNew('framesPerSecond',myscreen.framesPerSecond,aperwidth,contrast,speed,coherence,'xCenter=5/sqrt(2)','yCenter=5/sqrt(2)');
    stimulusUncued.dotsDownLeft = dotsInitNew('framesPerSecond',myscreen.framesPerSecond,aperwidth,contrast,speed,coherence,'xCenter=-5/sqrt(2)','yCenter=-5/sqrt(2)');
    stimulusUncued.dotsDownRight = dotsInitNew('framesPerSecond',myscreen.framesPerSecond,aperwidth,contrast,speed,coherence,'xCenter=5/sqrt(2)','yCenter=-5/sqrt(2)');
else
    stimulusUncued.dotsUpLeft = dotsInitNew('framesPerSecond',myscreen.framesPerSecond,aperwidth,contrast,speed,coherence,'xCenter=-5/sqrt(2)','yCenter=5/sqrt(2)','dotSize=0.1');
    stimulusUncued.dotsUpRight = dotsInitNew('framesPerSecond',myscreen.framesPerSecond,aperwidth,contrast,speed,coherence,'xCenter=5/sqrt(2)','yCenter=5/sqrt(2)','dotSize=0.1');
    stimulusUncued.dotsDownLeft = dotsInitNew('framesPerSecond',myscreen.framesPerSecond,aperwidth,contrast,speed,coherence,'xCenter=-5/sqrt(2)','yCenter=-5/sqrt(2)','dotSize=0.1');
    stimulusUncued.dotsDownRight = dotsInitNew('framesPerSecond',myscreen.framesPerSecond,aperwidth,contrast,speed,coherence,'xCenter=5/sqrt(2)','yCenter=-5/sqrt(2)','dotSize=0.1');
end

% set background color
stimulusUncued.backgroundColor = 0.5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the staircase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulusUncued = initStaircase(stimulusUncued)

stimulusUncued.staircase = doStaircase('init','upDown','nup=1','ndown=2','initialThreshold=0.8','minThreshold=0','initialStepsize=0.05','nTrials=24','stepRule=levitt','maxStepsize=1','minStepsize=.001');



