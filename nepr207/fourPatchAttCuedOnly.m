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

function myscreen = fourPatchAttCuedOnly(varargin)

% check arguments
getArgs(varargin, [], 'verbose=0');

% set default parameters
if ieNotDefined('atScanner'),atScanner = 1;end
if ieNotDefined('saveParam'),saveParam = 1;end
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
task{1}{1}.parameter.cue = [1 2 3 4];
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
global stimulusCued;

% if scanning, do synch to volume
task{1}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
if atScanner
    task{1}{1}.fudgeLastVolume = 1;
    task{1}{1}.synchToVol(end) = 1;
end

% set coherence
stimulusCued.coherence = 0; % 0-1

% set directions
stimulusCued.contrast = 1;
stimulusCued.initSpeed = 5; % in deg/s
stimulusCued.width = 5; % diameter
stimulusCued.initStair = initStair;

% init the stimulus
myscreen = initStimulus('stimulusCued',myscreen);
stimulusCued = initDots(stimulusCued,myscreen);

% init the staircase
if stimulusCued.initStair
    fprintf('\n(fourPatchAttCuedOnly) Initializing staircases from scratch...\n\n');
    stimulusCued = initStaircase(stimulusCued);
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

global stimulusCued;

if task.thistrial.thisseg == 1
    
    % set cue color to neutral
    stimulusCued.cuecolor{1} = [0 0 0];
    stimulusCued.cuecolor{2} = [0 0 0];
    stimulusCued.cuecolor{3} = [0 0 0];
    stimulusCued.cuecolor{4} = [0 0 0];

elseif task.thistrial.thisseg == 2
    
    % set cue color if this is a cued condition
    stimulusCued.cuecolor{task.thistrial.cue} = [1 1 1];
    task.thistrial.trialType = task.thistrial.cue;

elseif task.thistrial.thisseg == 3
    
    % reset cue color to neutral
    stimulusCued.cuecolor{1} = [0 0 0];
    stimulusCued.cuecolor{2} = [0 0 0];
    stimulusCued.cuecolor{3} = [0 0 0];
    stimulusCued.cuecolor{4} = [0 0 0];

elseif task.thistrial.thisseg == 4

    % change correct patch to move at a different speed
    task.thistrial.cuenum = task.thistrial.cue;

    % Get the new delta for this trial from the staircase
    [kWeb, stimulusCued.staircase] = doStaircase('testValue',stimulusCued.staircase);
    task.thistrial.kWeb = kWeb;
    fprintf('Current kWeb focal value: %0.2f \n\n', kWeb)

    % set dot base speed
    stimulusCued.dotsUpLeft = stimulusCued.dotsUpLeft.setSpeed(stimulusCued.dotsUpLeft,task.thistrial.speed1);
    stimulusCued.dotsUpRight = stimulusCued.dotsUpRight.setSpeed(stimulusCued.dotsUpRight,task.thistrial.speed2);
    stimulusCued.dotsDownLeft = stimulusCued.dotsDownLeft.setSpeed(stimulusCued.dotsDownLeft,task.thistrial.speed3);
    stimulusCued.dotsDownRight = stimulusCued.dotsDownRight.setSpeed(stimulusCued.dotsDownRight,task.thistrial.speed4);

    % adjust the speed of the dots by a value determined by Weber's law
    if task.thistrial.faster1 == 1
        stimulusCued.dotsUpLeft = stimulusCued.dotsUpLeft.setSpeed(stimulusCued.dotsUpLeft,task.thistrial.speed1 + task.thistrial.speed1*kWeb);
    end
    if task.thistrial.faster2 == 1
        stimulusCued.dotsUpRight = stimulusCued.dotsUpRight.setSpeed(stimulusCued.dotsUpRight,task.thistrial.speed2 + task.thistrial.speed2*kWeb);
    end
    if task.thistrial.faster3 == 1
        stimulusCued.dotsDownLeft = stimulusCued.dotsDownLeft.setSpeed(stimulusCued.dotsDownLeft,task.thistrial.speed3 + task.thistrial.speed3*kWeb);
    end
    if task.thistrial.faster4 == 1
        stimulusCued.dotsDownRight = stimulusCued.dotsDownRight.setSpeed(stimulusCued.dotsDownRight,task.thistrial.speed4 + task.thistrial.speed4*kWeb);
    end

elseif task.thistrial.thisseg == 6

    % Get the new delta for this trial from the staircase
    [kWeb, stimulusCued.staircase] = doStaircase('testValue',stimulusCued.staircase);
    task.thistrial.kWeb = kWeb;

    % set dot base speed
    stimulusCued.dotsUpLeft = stimulusCued.dotsUpLeft.setSpeed(stimulusCued.dotsUpLeft,task.thistrial.speed1);
    stimulusCued.dotsUpRight = stimulusCued.dotsUpRight.setSpeed(stimulusCued.dotsUpRight,task.thistrial.speed2);
    stimulusCued.dotsDownLeft = stimulusCued.dotsDownLeft.setSpeed(stimulusCued.dotsDownLeft,task.thistrial.speed3);
    stimulusCued.dotsDownRight = stimulusCued.dotsDownRight.setSpeed(stimulusCued.dotsDownRight,task.thistrial.speed4);
    
    % adjust the speed of the dots by a value determined by Weber's law
    if task.thistrial.faster1 == 2
        stimulusCued.dotsUpLeft = stimulusCued.dotsUpLeft.setSpeed(stimulusCued.dotsUpLeft,task.thistrial.speed1 + task.thistrial.speed1*kWeb);
    end
    if task.thistrial.faster2 == 2
        stimulusCued.dotsUpRight = stimulusCued.dotsUpRight.setSpeed(stimulusCued.dotsUpRight,task.thistrial.speed2 + task.thistrial.speed2*kWeb);
    end
    if task.thistrial.faster3 == 2
        stimulusCued.dotsDownLeft = stimulusCued.dotsDownLeft.setSpeed(stimulusCued.dotsDownLeft,task.thistrial.speed3 + task.thistrial.speed3*kWeb);
    end
    if task.thistrial.faster4 == 2
        stimulusCued.dotsDownRight = stimulusCued.dotsDownRight.setSpeed(stimulusCued.dotsDownRight,task.thistrial.speed4 + task.thistrial.speed4*kWeb);
    end
    
elseif task.thistrial.thisseg == 7

    % set cuecolor to indicate response period
    if task.thistrial.cuenum == 1
        stimulusCued.cuecolor{1} = [1 1 0];
    elseif task.thistrial.cuenum == 2
        stimulusCued.cuecolor{2} = [1 1 0];
    elseif task.thistrial.cuenum == 3
        stimulusCued.cuecolor{3} = [1 1 0];
    elseif task.thistrial.cuenum == 4
        stimulusCued.cuecolor{4} = [1 1 0];
    end

elseif task.thistrial.thisseg == 8

    % reset cue color to neutral
    stimulusCued.cuecolor{1} = [0 0 0];
    stimulusCued.cuecolor{2} = [0 0 0];
    stimulusCued.cuecolor{3} = [0 0 0];
    stimulusCued.cuecolor{4} = [0 0 0];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulusCued

% clear screen to gray
mglClearScreen(stimulusCued.backgroundColor);

if task.thistrial.thisseg == 1

    % draw the fixation cross in neutral colors
    mglFixationCrossDiag(1,1);

elseif task.thistrial.thisseg == 2
    
    % draw the cue
    mglFixationCrossDiag(1,1,stimulusCued.cuecolor);
    
elseif task.thistrial.thisseg == 3

    % draw the fixation cross in neutral colors
    mglFixationCrossDiag(1,1);
    
elseif task.thistrial.thisseg == 4

    % draw the fixation cross in neutral colors
    mglFixationCrossDiag(1,1);
    
    % update the dots
    stimulusCued.dotsUpLeft = stimulusCued.dotsUpLeft.update(stimulusCued.dotsUpLeft);
    stimulusCued.dotsUpRight = stimulusCued.dotsUpRight.update(stimulusCued.dotsUpRight);
    stimulusCued.dotsDownLeft = stimulusCued.dotsDownLeft.update(stimulusCued.dotsDownLeft);
    stimulusCued.dotsDownRight = stimulusCued.dotsDownRight.update(stimulusCued.dotsDownRight);
    
    % draw the dots
    stimulusCued.dotsUpLeft = stimulusCued.dotsUpLeft.draw(stimulusCued.dotsUpLeft);
    stimulusCued.dotsUpRight = stimulusCued.dotsUpRight.draw(stimulusCued.dotsUpRight);
    stimulusCued.dotsDownLeft = stimulusCued.dotsDownLeft.draw(stimulusCued.dotsDownLeft);
    stimulusCued.dotsDownRight = stimulusCued.dotsDownRight.draw(stimulusCued.dotsDownRight);

elseif task.thistrial.thisseg == 5

    % draw the fixation cross in neutral colors
    mglFixationCrossDiag(1,1);

elseif task.thistrial.thisseg == 6

    % draw the fixation cross in neutral colors
    mglFixationCrossDiag(1,1);

    % update the dots
    stimulusCued.dotsUpLeft = stimulusCued.dotsUpLeft.update(stimulusCued.dotsUpLeft);
    stimulusCued.dotsUpRight = stimulusCued.dotsUpRight.update(stimulusCued.dotsUpRight);
    stimulusCued.dotsDownLeft = stimulusCued.dotsDownLeft.update(stimulusCued.dotsDownLeft);
    stimulusCued.dotsDownRight = stimulusCued.dotsDownRight.update(stimulusCued.dotsDownRight);
    
    % draw the dots
    stimulusCued.dotsUpLeft = stimulusCued.dotsUpLeft.draw(stimulusCued.dotsUpLeft);
    stimulusCued.dotsUpRight = stimulusCued.dotsUpRight.draw(stimulusCued.dotsUpRight);
    stimulusCued.dotsDownLeft = stimulusCued.dotsDownLeft.draw(stimulusCued.dotsDownLeft);
    stimulusCued.dotsDownRight = stimulusCued.dotsDownRight.draw(stimulusCued.dotsDownRight);

elseif task.thistrial.thisseg == 7

    % draw the fixation cross in respond colors
    mglFixationCrossDiag(1,1,stimulusCued.cuecolor);

elseif task.thistrial.thisseg == 8

    % draw the fixation cross in neutral colors
    mglFixationCrossDiag(1,1,stimulusCued.cuecolor);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)

global stimulusCued

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
        stimulusCued.staircase = doStaircase('update',stimulusCued.staircase,corr);

        stimulusCued.cuecolor{1} = [0 1 0];
        stimulusCued.cuecolor{2} = [0 1 0];
        stimulusCued.cuecolor{3} = [0 1 0];
        stimulusCued.cuecolor{4} = [0 1 0];
        
    else % incorrect
        corr = 0;
        % report answer
        fprintf('Trial %d. ++ Incorrect ++. \n',task.trialnum);
        task.thistrial.correct = corr;
        
        % update staircases
        stimulusCued.staircase = doStaircase('update',stimulusCued.staircase,corr);

        stimulusCued.cuecolor{1} = [1 0 0];
        stimulusCued.cuecolor{2} = [1 0 0];
        stimulusCued.cuecolor{3} = [1 0 0];
        stimulusCued.cuecolor{4} = [1 0 0];

    end

    [kWeb, stimulusCued.staircase] = doStaircase('testValue',stimulusCued.staircase);
    fprintf('Previous kWeb focal value: %0.2f \n', task.thistrial.kWeb)
    fprintf('Next kWeb focal value: %0.2f \n\n', kWeb)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dots stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulusCued = initDots(stimulusCued,myscreen)

aperwidth = strcat('aperwidth=',num2str(stimulusCued.width));
contrast = strcat('contrast=',num2str(stimulusCued.contrast));
speed = strcat('speed=',num2str(stimulusCued.initSpeed));
coherence = strcat('coherence=',num2str(stimulusCued.coherence));

% init the dot patches
[~,name] = system('hostname');
if any(strfind(name,'oban'))
    stimulusCued.dotsUpLeft = dotsInitNew('framesPerSecond',myscreen.framesPerSecond,aperwidth,contrast,speed,coherence,'xCenter=-5/sqrt(2)','yCenter=5/sqrt(2)');
    stimulusCued.dotsUpRight = dotsInitNew('framesPerSecond',myscreen.framesPerSecond,aperwidth,contrast,speed,coherence,'xCenter=5/sqrt(2)','yCenter=5/sqrt(2)');
    stimulusCued.dotsDownLeft = dotsInitNew('framesPerSecond',myscreen.framesPerSecond,aperwidth,contrast,speed,coherence,'xCenter=-5/sqrt(2)','yCenter=-5/sqrt(2)');
    stimulusCued.dotsDownRight = dotsInitNew('framesPerSecond',myscreen.framesPerSecond,aperwidth,contrast,speed,coherence,'xCenter=5/sqrt(2)','yCenter=-5/sqrt(2)');
else
    stimulusCued.dotsUpLeft = dotsInitNew('framesPerSecond',myscreen.framesPerSecond,aperwidth,contrast,speed,coherence,'xCenter=-5/sqrt(2)','yCenter=5/sqrt(2)','dotSize=0.1');
    stimulusCued.dotsUpRight = dotsInitNew('framesPerSecond',myscreen.framesPerSecond,aperwidth,contrast,speed,coherence,'xCenter=5/sqrt(2)','yCenter=5/sqrt(2)','dotSize=0.1');
    stimulusCued.dotsDownLeft = dotsInitNew('framesPerSecond',myscreen.framesPerSecond,aperwidth,contrast,speed,coherence,'xCenter=-5/sqrt(2)','yCenter=-5/sqrt(2)','dotSize=0.1');
    stimulusCued.dotsDownRight = dotsInitNew('framesPerSecond',myscreen.framesPerSecond,aperwidth,contrast,speed,coherence,'xCenter=5/sqrt(2)','yCenter=-5/sqrt(2)','dotSize=0.1');
end

% set background color
stimulusCued.backgroundColor = 0.5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the staircase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulusCued = initStaircase(stimulusCued)

stimulusCued.staircase = doStaircase('init','upDown','nup=1','ndown=2','initialThreshold=0.4','minThreshold=0','initialStepsize=0.05','nTrials=24','stepRule=levitt','maxStepsize=1','minStepsize=.001');


