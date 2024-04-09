% compareGratingsWM.m
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

function myscreen = compareGratingsWM(varargin)

% check arguments
getArgs(varargin);

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
task{1}{1}.getResponse = [0 0 0 0 0 0 0 1 0];
task{1}{1}.segmin = [0.2, 0.4, 0.2, 0.4, 0.8, 11, 0.5, 2.5, 1]; % average trial time = 20s
task{1}{1}.segmax = [0.2, 0.4, 0.2, 0.4, 0.8, 11, 0.5, 2.5, 7];

% task parameters
task{1}{1}.parameter.sameStim = [1 2]; % 1 stim1, 2 stim2
task{1}{1}.parameter.offsetDir = [-1 1]; % 1 CW, -1 CCW
task{1}{1}.randVars.calculated.oriStim1 = nan;
task{1}{1}.randVars.calculated.oriStim2 = nan;
task{1}{1}.randVars.calculated.oriOffset1 = nan;
task{1}{1}.randVars.calculated.oriOffset2 = nan;
task{1}{1}.randVars.calculated.probeOri = nan;

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

% set grating parameters
stimulus.SF = 0.5;
stimulus.height = 20;
stimulus.width = 20;
stimulus.aperOuterHeight = 16; % minor axis diameter
stimulus.aperOuterWidth = 16; % major axis diameter
stimulus.outerHeightRatio = stimulus.aperOuterHeight/stimulus.height;
stimulus.outerWidthRatio = stimulus.aperOuterWidth/stimulus.width;
stimulus.aperInnerHeight = 2; % minor axis diameter
stimulus.aperInnerWidth = 2; % major axis diameter
stimulus.innerHeightRatio = stimulus.aperInnerHeight/stimulus.height;
stimulus.innerWidthRatio = stimulus.aperInnerWidth/stimulus.width;

% set fix parameters
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
stimulus = initGrating(stimulus,myscreen);

% init the staircase
stimulus.initStair = initStair;
if stimulus.initStair
    fprintf('\n(compareGratingsWM) Initializing staircases from scratch...\n\n');
    stimulus = initStaircase(stimulus);
else
    fprintf('\n(compareGratingsWM) Re-using staircase from previous run...\n\n');
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

global stimulus;

if task.thistrial.thisseg == 1

    % calculate the necessary offset for stimulus 1
    % Get the new delta for this trial from the staircase
    if task.thistrial.sameStim == 1
        [oriOffset1, stimulus.staircase1] = doStaircase('testValue',stimulus.staircase1);
        task.thistrial.oriOffset1 = oriOffset1;
        fprintf('Current oriOffset1 value: %0.2f \n\n', oriOffset1)
    elseif task.thistrial.sameStim == 2
        [oriOffset2, stimulus.staircase2] = doStaircase('testValue',stimulus.staircase2);
        task.thistrial.oriOffset2 = oriOffset2;
        fprintf('Current oriOffset2 value: %0.2f \n\n', oriOffset2)
    end

    % generate a random grating orientation for the first grating + make it orthogonal for the second grating
    task.thistrial.oriStim1 = randi(360);
    task.thistrial.oriStim2 = task.thistrial.oriStim1 + 90;

    % calculate the probe orientation
    if task.thistrial.sameStim == 1
        task.thistrial.probeOri = task.thistrial.oriStim1 + task.thistrial.offsetDir*task.thistrial.oriOffset1;
    elseif task.thistrial.sameStim == 2
        task.thistrial.probeOri = task.thistrial.oriStim2 + task.thistrial.offsetDir*task.thistrial.oriOffset2;
    end

elseif task.thistrial.thisseg == 7

    stimulus.feedbackColors{1} = [1 1 1];
    stimulus.feedbackColors{2} = [1 1 1];
    stimulus.feedbackColors{3} = [1 1 1];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus

% clear screen to gray
mglClearScreen(stimulus.backgroundColor);

if task.thistrial.thisseg == 1 % first stim presentation

    % draw the sample grating
    mglBltTexture(stimulus.grating, [0 0 stimulus.height stimulus.height], 0, 0, task.thistrial.oriStim1);
    mglBltTexture(stimulus.aperture, [0 0 stimulus.height*2 stimulus.height*2], 0, 0, 0);

    % draw fix
    % Vert
    mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, [0 0 0] );
    % left
    mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [0 0 0]);
    % right
    mglLines2( stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, [0 0 0]);

elseif task.thistrial.thisseg == 2 % blank

    % draw fix
    % Vert
    mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, [0 0 0] );
    % left
    mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [0 0 0]);
    % right
    mglLines2( stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, [0 0 0]);

elseif task.thistrial.thisseg == 3 % second stim presentation

    % draw the sample grating
    mglBltTexture(stimulus.grating, [0 0 stimulus.height stimulus.height], 0, 0, task.thistrial.oriStim2);
    mglBltTexture(stimulus.aperture, [0 0 stimulus.height*2 stimulus.height*2], 0, 0, 0);

    % draw fix
    % Vert
    mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, [0 0 0] );
    % left
    mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [0 0 0]);
    % right
    mglLines2( stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, [0 0 0]);

elseif task.thistrial.thisseg == 4 % blank

    % draw fix
    % Vert
    mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, [0 0 0] );
    % left
    mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [0 0 0]);
    % right
    mglLines2( stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, [0 0 0]);
    
elseif task.thistrial.thisseg == 5 % cue period

    % draw cue
    if task.thistrial.sameStim == 1
        mglBltTexture(stimulus.cueText1,[0 0],'left','top');
    elseif task.thistrial.sameStim == 2
        mglBltTexture(stimulus.cueText2,[0 0],'left','top');
    end

elseif task.thistrial.thisseg == 6 % working memory period

    % draw fix
    % Vert
    mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, [1 0.5 1] );
    % left
    mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [1 0.5 1]);
    % right
    mglLines2( stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, [1 0.5 1]);

elseif task.thistrial.thisseg == 7

    % draw the probe grating
    mglBltTexture(stimulus.grating, [0 0 stimulus.height stimulus.height], 0, 0, task.thistrial.probeOri);
    mglBltTexture(stimulus.aperture, [0 0 stimulus.height*2 stimulus.height*2], 0, 0, 0);

    % draw fix
    % Vert
    mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, [0 0 0] );
    % left
    mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [0 0 0]);
    % right
    mglLines2( stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, [0 0 0]);

elseif task.thistrial.thisseg == 8 % response period

    % draw fix
    % Vert
    mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, stimulus.feedbackColors{1});
    % left
    mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, stimulus.feedbackColors{2});
    % right
    mglLines2( stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, stimulus.feedbackColors{3});

elseif task.thistrial.thisseg == 9 % inter-trial interval

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
    if task.thistrial.offsetDir == 1
        correctbutton = 1;
    elseif task.thistrial.offsetDir == -1
        correctbutton = 2;
    end

    % see if it is correct
    if isequal(task.thistrial.whichButton,correctbutton) % correct
        corr = 1;
        % report answer
        fprintf('Trial %d. !! Correct !!. \n',task.trialnum);
        task.thistrial.correct = corr;

        % update staircases
        if task.thistrial.sameStim == 1
            stimulus.staircase1 = doStaircase('update',stimulus.staircase1,corr);
        elseif task.thistrial.sameStim == 2
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
        if task.thistrial.sameStim == 1
            stimulus.staircase1 = doStaircase('update',stimulus.staircase1,corr);
        elseif task.thistrial.sameStim == 2
            stimulus.staircase2 = doStaircase('update',stimulus.staircase2,corr);
        end

        stimulus.feedbackColors{1} = [1 0 0];
        stimulus.feedbackColors{2} = [1 0 0];
        stimulus.feedbackColors{3} = [1 0 0];

    end
    
    if task.thistrial.sameStim == 1
        [oriOffset1, stimulus.staircase1] = doStaircase('testValue',stimulus.staircase1);
        fprintf('Previous oriOffset1 value: %0.2f \n', task.thistrial.oriOffset1)
        fprintf('Next oriOffset1 value: %0.2f \n\n', oriOffset1)
    elseif task.thistrial.sameStim == 2
        [oriOffset2, stimulus.staircase2] = doStaircase('testValue',stimulus.staircase2);
        fprintf('Previous oriOffset2 value: %0.2f \n', task.thistrial.oriOffset2)
        fprintf('Next oriOffset2 value: %0.2f \n\n', oriOffset2)
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dots stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initGrating(stimulus,myscreen)

stimulus.pixRes = min(myscreen.screenHeight/myscreen.imageHeight, myscreen.screenWidth/myscreen.imageWidth);

% grating phase
phase = randi(360);

% create the grating texture (just need phase and sf, mglBltTexture will handle the rotations)
% make a grating  but now scale it
grating = mglMakeGrating(stimulus.width, stimulus.height, stimulus.SF, 0, phase, stimulus.pixRes, stimulus.pixRes);

% scale to range of display
grating = 255*(grating+1)/2;

% create a texture
stimulus.grating = mglCreateTexture(grating, [], 1);

% create the aperture texture - make a elliptical (circular) aperture
% grating = grating .*  mkDisc(size(grating), (length(grating)/2)-2, (size(grating)+1)/2, 1);
stimSize = size(grating,1);
apertureOuter = circStim([stimSize/2*stimulus.outerHeightRatio stimSize/2*stimulus.outerWidthRatio],[stimSize stimSize],[ceil(stimSize/2) ceil(stimSize/2)]);
apertureInner = ~circStim([stimSize/2*stimulus.innerHeightRatio stimSize/2*stimulus.innerWidthRatio],[stimSize stimSize],[ceil(stimSize/2) ceil(stimSize/2)]);
apertureAlpha = ~and(apertureOuter,apertureInner);
apertureAlpha = padarray(apertureAlpha,[(stimSize-1)/2 (stimSize-1)/2],1,'both'); % pad the array so that we block out the grating edges

% make rgb matrix (n x m x 3)
apertureRGB = 0.5 * ones(size(apertureAlpha,1),size(apertureAlpha,2),3);

% tack on the alpha values
aperture = 255 * cat(3,apertureRGB,apertureAlpha);
% keyboard
% create a texture
stimulus.aperture = mglCreateTexture(aperture, [], 1);

% build the cue texture
mglTextSet('Helvetica',56,[0 0 0],0,0,0);
stimulus.cueText1 = mglText('1');
stimulus.cueText2 = mglText('2');

% set background color
stimulus.backgroundColor = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the staircase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initStaircase(stimulus)

stimulus.staircase1 = doStaircase('init','upDown','nup=1','ndown=2','initialThreshold=15','initialStepsize=2','minThreshold=0','nTrials=24','stepRule=levitt','maxStepsize=1','minStepsize=.01');
stimulus.staircase2 = doStaircase('init','upDown','nup=1','ndown=2','initialThreshold=15','initialStepsize=2','minThreshold=0','nTrials=24','stepRule=levitt','maxStepsize=1','minStepsize=.01');

