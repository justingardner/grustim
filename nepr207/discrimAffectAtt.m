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

function myscreen = discrimAffectAtt(varargin)

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
task{1}{1}.getResponse = [0 0 0 1 0];
task{1}{1}.segmin = [0.5, 3, 1, 1.75, 1]; % average trial time = 9.5s
task{1}{1}.segmax = [0.5, 3, 1, 1.75, 7];

% task parameters
task{1}{1}.parameter.affectCue = [-1 0]; % -1 angry, 0 neutral
task{1}{1}.parameter.offsetDir = [-1 1]; % -1 CW, 1 CCW
task{1}{1}.randVars.uniform.uncuedDir = [-1 1]; % -1 CW, 1 CCW
task{1}{1}.randVars.uniform.angrySide = [-1 1]; % -1 left, 1 right
task{1}{1}.randVars.uniform.angryImage = 1:217;
task{1}{1}.randVars.uniform.neutralImage = 1:216;
task{1}{1}.randVars.calculated.cueDir = nan; % -1 left, 1 right
task{1}{1}.randVars.calculated.oriOffsetAngry = nan;
task{1}{1}.randVars.calculated.oriOffsetNeutral = nan;
task{1}{1}.randVars.calculated.cueAngle = nan;
task{1}{1}.randVars.calculated.disAngle = nan;
task{1}{1}.randVars.calculated.correctButton = nan;
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
stimulus.gratingEcc = 5;
stimulus.gaussStd = 0.5;

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
stimulus = initFaces(stimulus);

% init the staircase
stimulus.initStair = initStair;
if stimulus.initStair
    fprintf('\n(discrimAffectAtt) Initializing staircases from scratch...\n\n');
    stimulus = initStaircase(stimulus);
else
    fprintf('\n(discrimAffectAtt) Re-using staircase from previous run...\n\n');
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
    
    % choose both affect images and make their textures
    stimulus.currentAngryFace = mglMetalCreateTexture(cat(3,rescale(stimulus.angryFaces(:,:,:,task.thistrial.angryImage)),ones(150)));
    stimulus.currentNeutralFace = mglMetalCreateTexture(cat(3,rescale(stimulus.neutralFaces(:,:,:,task.thistrial.neutralImage)),ones(150)));

    % figure out cue left or cue right
    if (task.thistrial.affectCue == -1 && task.thistrial.angrySide == -1) || (task.thistrial.affectCue == 0 && task.thistrial.angrySide == 1)
        task.thistrial.cueDir = -1;
    elseif (task.thistrial.affectCue == -1 && task.thistrial.angrySide == 1) || (task.thistrial.affectCue == 0 && task.thistrial.angrySide == -1)
        task.thistrial.cueDir = 1;
    end

    % Get the new delta for this trial from the staircase - need to instantiate both angry and neutral, or else the uncued one has no texture
    [oriOffset, stimulus.staircase1] = doStaircase('testValue',stimulus.staircase1);
    task.thistrial.oriOffsetAngry = oriOffset;
    fprintf('Current oriOffset value (Angry face): %0.2f \n\n', oriOffset)

    [oriOffset, stimulus.staircase2] = doStaircase('testValue',stimulus.staircase2);
    task.thistrial.oriOffsetNeutral = oriOffset;
    fprintf('Current oriOffset value (Neutral face): %0.2f \n\n', oriOffset)

    if task.thistrial.affectCue == -1
        % create the grating texture (just need phase and sf, mglBltTexture will handle the rotations)
        % make a grating  but now scale it
        task.thistrial.cueAngle = 90 + task.thistrial.offsetDir*task.thistrial.oriOffsetAngry;
        task.thistrial.disAngle = 90 + task.thistrial.uncuedDir*task.thistrial.oriOffsetNeutral;
        gratingAngry = mglMakeGrating(5, 5, stimulus.SF, task.thistrial.cueAngle, randi(360));
        gratingNeutral = mglMakeGrating(5, 5, stimulus.SF, task.thistrial.disAngle, randi(360));

        if task.thistrial.cueAngle > 90
            task.thistrial.correctButton = 6;
        else
            task.thistrial.correctButton = 1;
        end

    elseif task.thistrial.affectCue == 0
        % create the grating texture (just need phase and sf, mglBltTexture will handle the rotations)
        % make a grating  but now scale it
        task.thistrial.cueAngle = 90 + task.thistrial.offsetDir*task.thistrial.oriOffsetNeutral;
        task.thistrial.disAngle = 90 + task.thistrial.uncuedDir*task.thistrial.oriOffsetAngry;
        gratingAngry = mglMakeGrating(5, 5, stimulus.SF, 90 + task.thistrial.uncuedDir*task.thistrial.oriOffsetAngry, randi(360));
        gratingNeutral = mglMakeGrating(5, 5, stimulus.SF, 90 + task.thistrial.offsetDir*task.thistrial.oriOffsetNeutral, randi(360));

        if task.thistrial.cueAngle > 90
            task.thistrial.correctButton = 6;
        else
            task.thistrial.correctButton = 1;
        end

    end
    
    angryGabor = cat(3,repmat((gratingAngry.*stimulus.gaussian+1)/2, [1 1 3]),ones(size(stimulus.gaussian)));
    % scale to range of display
    stimulus.angryTex = mglMetalCreateTexture(angryGabor);

    neutralGabor = cat(3,repmat((gratingNeutral.*stimulus.gaussian+1)/2, [1 1 3]),ones(size(stimulus.gaussian)));
    % scale to range of display
    stimulus.neutralTex = mglMetalCreateTexture(neutralGabor);

elseif task.thistrial.thisseg == 4

    stimulus.feedbackColors{1} = [0 0 0];
    if task.thistrial.cueDir == 1
        stimulus.feedbackColors{2} = [0 0 0];
        stimulus.feedbackColors{3} = [1 1 1];
    elseif task.thistrial.cueDir == -1
        stimulus.feedbackColors{2} = [1 1 1];
        stimulus.feedbackColors{3} = [0 0 0];
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus

% clear screen to gray
mglClearScreen(stimulus.backgroundColor);

if task.thistrial.thisseg == 1 % blank

    % draw fix
    % Vert
    mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, [0 0 0] );
    % left
    mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [0 0 0]);
    % right
    mglLines2( stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, [0 0 0]);

elseif task.thistrial.thisseg == 2 % faces
    
    if task.thistrial.affectCue == -1
        % draw the faces
        % mglMetalBltTexture([stimulus.currentAngryFace stimulus.currentNeutralFace], [task.thistrial.cueDir*stimulus.gratingEcc 0; -task.thistrial.cueDir*stimulus.gratingEcc 0]);
        mglMetalBltTexture(stimulus.currentAngryFace, [task.thistrial.cueDir*stimulus.gratingEcc 0]);
        mglMetalBltTexture(stimulus.currentNeutralFace, [-task.thistrial.cueDir*stimulus.gratingEcc 0]);
    elseif task.thistrial.affectCue == 0
        % draw the faces
        % mglMetalBltTexture([stimulus.currentAngryFace stimulus.currentNeutralFace], [task.thistrial.cueDir*stimulus.gratingEcc 0; -task.thistrial.cueDir*stimulus.gratingEcc 0]);
        mglMetalBltTexture(stimulus.currentAngryFace, [-task.thistrial.cueDir*stimulus.gratingEcc 0]);
        mglMetalBltTexture(stimulus.currentNeutralFace, [task.thistrial.cueDir*stimulus.gratingEcc 0]);
    end

    % draw fix
    % Vert
    mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, [0 0 0] );
    % left
    mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [0 0 0]);
    % right
    mglLines2( stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, [0 0 0]);

elseif task.thistrial.thisseg == 3 % grating period
    
    if task.thistrial.affectCue == -1
        % draw the cued grating on the angry side
        mglMetalBltTexture(stimulus.angryTex, [task.thistrial.cueDir*stimulus.gratingEcc 0]);
        mglMetalBltTexture(stimulus.neutralTex, [-task.thistrial.cueDir*stimulus.gratingEcc 0]);
    elseif task.thistrial.affectCue == 0
        % draw the cued grating on the neutral side 
        mglMetalBltTexture(stimulus.angryTex, [-task.thistrial.cueDir*stimulus.gratingEcc 0]);
        mglMetalBltTexture(stimulus.neutralTex, [task.thistrial.cueDir*stimulus.gratingEcc 0]);
    end

    % draw fix
    % Vert
    mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, [0 0 0] );
    % left
    mglLines2(stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [0 0 0]);
    % right
    mglLines2(stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, [0 0 0]);

elseif task.thistrial.thisseg == 4 % response period

    % draw fix
    % Vert
    mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, stimulus.feedbackColors{1});
    % left
    mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, stimulus.feedbackColors{2});
    % right
    mglLines2( stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, stimulus.feedbackColors{3});

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

    if ~any(task.thistrial.whichButton == [1 6]) % [1 6] should (hopefully) correspond to left/right buttons
        error('Check your button inputs')
    end
    
%     % determine what the correct button should be
%     if task.thistrial.offsetDir == -1
%         correctbutton = 6;
%     elseif task.thistrial.offsetDir == 1
%         correctbutton = 1;
%     end

    % see if it is correct
    if isequal(task.thistrial.whichButton,task.thistrial.correctButton) % correct
        corr = 1;
        % report answer
        fprintf('Trial %d. !! Correct !!. \n',task.trialnum);
        task.thistrial.correct = corr;

        % update staircases
        if task.thistrial.affectCue == -1
            stimulus.staircase1 = doStaircase('update',stimulus.staircase1,corr);
        elseif task.thistrial.affectCue == 0
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
        if task.thistrial.affectCue == -1
            stimulus.staircase1 = doStaircase('update',stimulus.staircase1,corr);
        elseif task.thistrial.affectCue == 0
            stimulus.staircase2 = doStaircase('update',stimulus.staircase2,corr);
        end

        stimulus.feedbackColors{1} = [1 0 0];
        stimulus.feedbackColors{2} = [1 0 0];
        stimulus.feedbackColors{3} = [1 0 0];

    end

    % print out the updated values for the staircase
    [oriOffset, stimulus.staircase1] = doStaircase('testValue',stimulus.staircase1);
    fprintf('Previous oriOffset value (Angry faces): %0.2f \n', task.thistrial.oriOffsetAngry)
    fprintf('Next oriOffset value (Angry faces): %0.2f \n\n', oriOffset)

    [oriOffset, stimulus.staircase2] = doStaircase('testValue',stimulus.staircase2);
    fprintf('Previous oriOffset value (Neutral faces): %0.2f \n', task.thistrial.oriOffsetNeutral)
    fprintf('Next oriOffset value (Neutral faces): %0.2f \n\n', oriOffset)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the grating stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initGrating(stimulus,myscreen)

% make gaussian
gaussian = mglMakeGaussian(5, 5, stimulus.gaussStd, stimulus.gaussStd);

% create a texture
stimulus.gaussian = gaussian;

% set background color
stimulus.backgroundColor = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the faces stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initFaces(stimulus)

% load angry and neutral faces
load('~/Documents/MATLAB/nepr207/neutralFaces.mat')
load('~/Documents/MATLAB/nepr207/angryFaces.mat')
stimulus.angryFaces = angryImgStack;
stimulus.neutralFaces = neutralImgStack;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the staircase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initStaircase(stimulus)

stimulus.staircase1 = doStaircase('init','upDown','nup=1','ndown=2','initialThreshold=20','initialStepsize=2','minThreshold=0','nTrials=24','stepRule=levitt','maxStepsize=1','minStepsize=.01');
stimulus.staircase2 = doStaircase('init','upDown','nup=1','ndown=2','initialThreshold=20','initialStepsize=2','minThreshold=0','nTrials=24','stepRule=levitt','maxStepsize=1','minStepsize=.01');


