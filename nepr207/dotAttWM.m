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

function myscreen = dotAttWM(varargin)

% check arguments
getArgs(varargin);

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
task{1}{1}.getResponse = [0 0 0 1 0];
task{1}{1}.collectEyeData = false;
task{1}{1}.segmin = [0.5, 1, nan, 2.5, 1]; % average trial time = 18s
task{1}{1}.segmax = [0.5, 1, nan, 2.5, 7];
task{1}{1}.segdur{3} = [8 12];

% task parameters
task{1}{1}.parameter.cue = [-1 -0.1 0.1 1]; % -1 0 1 calculated variable, put the distributed together
task{1}{1}.randVars.uniform.probePosition = [-1 1];
task{1}{1}.randVars.calculated.duration = nan; % this should be 0 or 1 for short/long
task{1}{1}.randVars.calculated.trialType = nan;
task{1}{1}.randVars.calculated.xOffset1 = nan;
task{1}{1}.randVars.calculated.xOffset2 = nan;
task{1}{1}.randVars.calculated.cueJitterLeft = nan;
task{1}{1}.randVars.calculated.cueJitterRight = nan;
task{1}{1}.random = 1;

task{1}{1}.numTrials = inf; % 16 conditions

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

% set dot parameters
stimulus.cuedotEcc = 5;
stimulus.cuedotX = [stimulus.cuedotEcc; -stimulus.cuedotEcc];
stimulus.cuedotY = [0; 0];
stimulus.dotSize = 0.5;
stimulus.dotColor = [1 1 1];
stimulus.probeDotY = 0;

% set cue parameters
width = 0.5;
stimulus.cueLeftX0 = 0;
stimulus.cueLeftX1 = -width;
stimulus.cueLeftY0 = 0;
stimulus.cueLeftY1 = 0;

stimulus.cueRightX0 = 0;
stimulus.cueRightX1 = width;
stimulus.cueRightY0 = 0;
stimulus.cueRightY1 = 0;

stimulus.cueVertX0 = 0;
stimulus.cueVertX1 = 0;
stimulus.cueVertY0 = width;
stimulus.cueVertY1 = -width;

% init the stimulus
myscreen = initStimulus('stimulus',myscreen);
stimulus = initDots(stimulus,myscreen);

% init the staircase
stimulus.initStair = initStair;
if stimulus.initStair
    fprintf('\n(dotAttWM) Initializing staircases from scratch...\n\n');
    stimulus = initStaircase(stimulus);
else
    fprintf('\n(dotAttWM) Re-using staircase from previous run...\n\n');
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

    % set trialType
    if task.thistrial.seglen(3) == 8
        if task.thistrial.cue == -1
            task.thistrial.trialType = 1;

        elseif task.thistrial.cue == -0.1 || task.thistrial.cue == 0.1
            task.thistrial.trialType = 2;

        elseif task.thistrial.cue == 1
            task.thistrial.trialType = 3;

        end

        task.thistrial.duration = 8;

    elseif task.thistrial.seglen(3) == 12
        if task.thistrial.cue == -1
            task.thistrial.trialType = 4;

        elseif task.thistrial.cue == -0.1 || task.thistrial.cue == 0.1
            task.thistrial.trialType = 5;

        elseif task.thistrial.cue == 1
            task.thistrial.trialType = 6;
            
        end

        task.thistrial.duration = 12;

    end

    % jitter
    task.thistrial.cueJitterLeft = rand - 0.5;
    task.thistrial.cueJitterRight = rand - 0.5;

elseif task.thistrial.thisseg == 4

    stimulus.feedbackColors{1} = [1 1 1];
    stimulus.feedbackColors{2} = [1 1 1];
    stimulus.feedbackColors{3} = [1 1 1];

    % Get the new delta for this trial from the staircase
    if task.thistrial.trialType < 4
        [xOffset, stimulus.staircase1] = doStaircase('testValue',stimulus.staircase1);
        task.thistrial.xOffset1 = xOffset;
        fprintf('Current xOffset value (%d second WM): %0.2f \n\n', task.thistrial.seglen(3), xOffset)
    elseif task.thistrial.trialType > 3
        [xOffset, stimulus.staircase2] = doStaircase('testValue',stimulus.staircase2);
        task.thistrial.xOffset2 = xOffset;
        fprintf('Current xOffset value (%d second WM): %0.2f \n\n', task.thistrial.seglen(3), xOffset)
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

    % draw cue
    if task.thistrial.cue == -1
        % Vert
        mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, [0 0 0] );
        % left
        mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [1 1 1]);
        % right
        mglLines2( stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, [0 0 0]);
    elseif round(task.thistrial.cue) == 0
        % Vert
        mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, [0 0 0] );
        % left
        mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [1 1 1]);
        % right
        mglLines2( stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, [1 1 1]);
    elseif task.thistrial.cue == 1
        % Vert
        mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, [0 0 0] );
        % left
        mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [0 0 0]);
        % right
        mglLines2( stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, [1 1 1]);
    end

elseif task.thistrial.thisseg == 2

    % draw cue
    if task.thistrial.cue == -1
        % Vert
        mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, [0 0 0] );
        % left
        mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [1 1 1]);
        % right
        mglLines2( stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, [0 0 0]);
    elseif round(task.thistrial.cue) == 0
        % Vert
        mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, [0 0 0] );
        % left
        mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [1 1 1]);
        % right
        mglLines2( stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, [1 1 1]);
    elseif task.thistrial.cue == 1
        % Vert
        mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, [0 0 0] );
        % left
        mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [0 0 0]);
        % right
        mglLines2( stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, [1 1 1]);
    end

    % draw the left/right dots
    mglBltTexture(stimulus.gauss,[stimulus.cuedotX(1)+task.thistrial.cueJitterLeft,stimulus.cuedotY(1)]);
    mglBltTexture(stimulus.gauss,[stimulus.cuedotX(2)+task.thistrial.cueJitterRight,stimulus.cuedotY(2)]);

elseif task.thistrial.thisseg == 3

    % memory period

    % draw fix
    % Vert
    mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, [0 0 0] );
    % left
    mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [0 0 0]);
    % right
    mglLines2( stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, [0 0 0]);

elseif task.thistrial.thisseg == 4

    % draw fix
    % Vert
    mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, stimulus.feedbackColors{1});
    % left
    mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, stimulus.feedbackColors{2});
    % right
    mglLines2( stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, stimulus.feedbackColors{3});

    % draw probe
    cursign = sign(task.thistrial.cue);
    if task.thistrial.trialType < 4
        if cursign < 0
            mglBltTexture(stimulus.gauss,...
                [stimulus.cuedotEcc*cursign + (task.thistrial.xOffset1*task.thistrial.probePosition+task.thistrial.cueJitterLeft),... % xpos
                stimulus.probeDotY]);                                                                   % ypos
        else
            mglBltTexture(stimulus.gauss,...
                [stimulus.cuedotEcc*cursign + (task.thistrial.xOffset1*task.thistrial.probePosition+task.thistrial.cueJitterRight),... % xpos
                stimulus.probeDotY]);
        end
    else
        if cursign < 0
        mglBltTexture(stimulus.gauss,...
            [stimulus.cuedotEcc*cursign + (task.thistrial.xOffset2*task.thistrial.probePosition+task.thistrial.cueJitterLeft),... % xpos
             stimulus.probeDotY]);                                                                   % ypos
        else
            mglBltTexture(stimulus.gauss,...
            [stimulus.cuedotEcc*cursign + (task.thistrial.xOffset2*task.thistrial.probePosition+task.thistrial.cueJitterRight),... % xpos
             stimulus.probeDotY]);  
        end
    end

elseif task.thistrial.thisseg == 5

    % draw fix
    % Vert
    mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, [0 0 0] );
    % left
    mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [0 0 0]);
    % right
    mglLines2( stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, [0 0 0]);

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
    if task.thistrial.probePosition == -1
        correctbutton = 1;
    elseif task.thistrial.probePosition == 1
        correctbutton = 2;
    end

    % see if it is correct
    if isequal(task.thistrial.whichButton,correctbutton) % correct
        corr = 1;
        % report answer
        fprintf('Trial %d. !! Correct !!. \n',task.trialnum);
        task.thistrial.correct = corr;

        % update staircases
        if task.thistrial.trialType < 4
            stimulus.staircase1 = doStaircase('update',stimulus.staircase1,corr);
        elseif task.thistrial.trialType > 3
            stimulus.staircase2 = doStaircase('update',stimulus.staircase2,corr);
        end
        
        stimulus.feedbackColors{1} = [0 1 0];
        stimulus.feedbackColors{2} = [0 1 0];
        stimulus.feedbackColors{3} = [0 1 0];

    else % incorrect
        corr = 0;
        % report answer
        fprintf('Trial %d. ++ Incorrect ++. \n',task.trialnum);
        task.thistrial.correct = corr;
        
        % update staircases
        if task.thistrial.trialType < 4
            stimulus.staircase1 = doStaircase('update',stimulus.staircase1,corr);
        elseif task.thistrial.trialType > 3
            stimulus.staircase2 = doStaircase('update',stimulus.staircase2,corr);
        end

        stimulus.feedbackColors{1} = [1 0 0];
        stimulus.feedbackColors{2} = [1 0 0];
        stimulus.feedbackColors{3} = [1 0 0];

    end

    if task.thistrial.trialType < 4
        [xOffset, stimulus.staircase1] = doStaircase('testValue',stimulus.staircase1);
        fprintf('Previous xOffset value (%d second WM): %0.2f \n', task.thistrial.seglen(3), task.thistrial.xOffset1)
        fprintf('Next xOffset value (%d second WM): %0.2f \n\n', task.thistrial.seglen(3), xOffset)
    else
        [xOffset, stimulus.staircase2] = doStaircase('testValue',stimulus.staircase2);
        fprintf('Previous xOffset value (%d second WM): %0.2f \n', task.thistrial.seglen(3), task.thistrial.xOffset2)
        fprintf('Next xOffset value (%d second WM): %0.2f \n\n', task.thistrial.seglen(3), xOffset)
    end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dots stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initDots(stimulus,myscreen)

stimulus.pixRes = min(myscreen.screenHeight/myscreen.imageHeight, myscreen.screenWidth/myscreen.imageWidth);

gauss = mglMakeGaussian(2,2,0.3,0.3,0,0,stimulus.pixRes,stimulus.pixRes);
gauss = 255*(gauss+1)/2;
stimulus.gauss = mglCreateTexture(gauss);

% set background color
stimulus.backgroundColor = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the staircase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initStaircase(stimulus)

stimulus.staircase1 = doStaircase('init','upDown','nup=2','ndown=1','initialThreshold=1','initialStepsize=0.05','minThreshold=0','nTrials=24','stepRule=levitt','maxStepsize=1','minStepsize=.01');
stimulus.staircase2 = doStaircase('init','upDown','nup=2','ndown=1','initialThreshold=1','initialStepsize=0.05','minThreshold=0','nTrials=24','stepRule=levitt','maxStepsize=1','minStepsize=.01');


