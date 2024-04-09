% flashSeqWM.m
%
%        $Id$
%      usage: flashSeqWM(varargin)
%         by: austin kuo
%       date: 4/8/2024
%  copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%    purpose: nepr207 stimulus (2024)
%
% varargin - specify 'atScanner', 'saveParam', defaults = 0

function myscreen = flashSeqWM(varargin)

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

% global stimulus parameters
global stimulus

% task
stimulus.maxSequence = 6;
task{1}{1}.getResponse = [zeros(1,2*stimulus.maxSequence+1) 0 1 0];
task{1}{1}.collectEyeData = false;
% task structure: [cue, [flash1,blank1,...,flash9,blank9] wait, response, iti] - total: 22 segments
task{1}{1}.segmin = [1.1, repmat([0.4,0.2],1,stimulus.maxSequence), 8, 5, 1]; % average trial time = 23s
task{1}{1}.segmax = [1.1, repmat([0.4,0.2],1,stimulus.maxSequence), 8, 5, 6];

% task parameters
task{1}{1}.parameter.sequenceLength = linspace(0,stimulus.maxSequence,4);
task{1}{1}.parameter.cue = [-1 1]; % -1 left, 1 right
task{1}{1}.randVars.calculated.trialType = nan;
task{1}{1}.randVars.calculated.nResponses = 0;
task{1}{1}.randVars.calculated.noMoreResponses = 0;
task{1}{1}.random = 1;

task{1}{1}.fullSequence = {};
task{1}{1}.correctSequence = {};
task{1}{1}.distractorSequence = {};
task{1}{1}.sequenceResponse = {};

task{1}{1}.numTrials = inf;

% initialize the task
for phaseNum = 1:length(task)
    [task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end

% if scanning, do synch to volume
task{1}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
if atScanner
    task{1}{1}.fudgeLastVolume = 1;
    task{1}{1}.synchToVol(end) = 1;
end

% set dot parameters
stimulus.symbolEcc = 5;

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
stimulus = initShapes(stimulus,myscreen);

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
    switch task.thistrial.sequenceLength
        case task.parameter.sequenceLength(1)
            task.thistrial.trialType = 1;
            task.fullSequence{task.trialnum} = zeros(1,stimulus.maxSequence);
            task.distractorSequence{task.trialnum} = randi(3,[1,stimulus.maxSequence])-2;
            task.correctSequence{task.trialnum} = nan(1,stimulus.maxSequence);
        case task.parameter.sequenceLength(2)
            task.thistrial.trialType = 2;
            memSegs = sort(randsample(stimulus.maxSequence,task.thistrial.sequenceLength));
            memLR = zeros(1,stimulus.maxSequence);
            memLR(memSegs) = 1;
            memLR = memLR.*randsample([-1 1],stimulus.maxSequence,true);
            task.fullSequence{task.trialnum} = memLR;
            task.distractorSequence{task.trialnum} = memLR(randperm(length(memLR)));
            task.correctSequence{task.trialnum} = nonzeros(task.fullSequence{task.trialnum})';
        case task.parameter.sequenceLength(3)
            task.thistrial.trialType = 3;
            memSegs = sort(randsample(stimulus.maxSequence,task.thistrial.sequenceLength));
            memLR = zeros(1,stimulus.maxSequence);
            memLR(memSegs) = 1;
            memLR = memLR.*randsample([-1 1],stimulus.maxSequence,true);
            task.fullSequence{task.trialnum} = memLR;
            task.distractorSequence{task.trialnum} = memLR(randperm(length(memLR)));
            task.correctSequence{task.trialnum} = nonzeros(task.fullSequence{task.trialnum})';
        case task.parameter.sequenceLength(4)
            task.thistrial.trialType = 4;
            memSegs = sort(randsample(stimulus.maxSequence,task.thistrial.sequenceLength));
            memLR = zeros(1,stimulus.maxSequence);
            memLR(memSegs) = 1;
            memLR = memLR.*randsample([-1 1],stimulus.maxSequence,true);
            task.fullSequence{task.trialnum} = memLR;
            task.distractorSequence{task.trialnum} = memLR(randperm(length(memLR)));
            task.correctSequence{task.trialnum} = nonzeros(task.fullSequence{task.trialnum})';
    end

elseif task.thistrial.thisseg == length(task.segmin) % last seg
    
    if all(isnan(task.correctSequence{task.trialnum})) && (length(task.correctSequence) > length(task.sequenceResponse))
        corr = 1;
        % report answer
        fprintf('Trial %d. No response !! Correct !!. \n',task.trialnum);
        task.thistrial.correct = corr;

        stimulus.feedbackColors{1} = [0 1 0];
        stimulus.feedbackColors{2} = [0 1 0];
        stimulus.feedbackColors{3} = [0 1 0];
    else
        % reset the feedback colors for next trial
        stimulus.feedbackColors{1} = [1 1 1];
        stimulus.feedbackColors{2} = [1 1 1];
        stimulus.feedbackColors{3} = [1 1 1];
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus

% clear screen to gray
mglClearScreen(stimulus.backgroundColor);

if task.thistrial.thisseg == 1 % cue period

    % draw cue
    if task.thistrial.cue == -1
        % Vert
        mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, [0 0 0] );
        % left
        mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [1 1 1]);
        % right
        mglLines2( stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, [0 0 0]);
    elseif task.thistrial.cue == 1
        % Vert
        mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, [0 0 0] );
        % left
        mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [0 0 0]);
        % right
        mglLines2( stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, [1 1 1]);
    end

elseif ismember(task.thistrial.thisseg,2:2:stimulus.maxSequence*2) % flash symbol segments

    % draw cue
    if task.thistrial.cue == -1
        % Vert
        mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, [0 0 0] );
        % left
        mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [1 1 1]);
        % right
        mglLines2( stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, [0 0 0]);
    elseif task.thistrial.cue == 1
        % Vert
        mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, [0 0 0] );
        % left
        mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [0 0 0]);
        % right
        mglLines2( stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, [1 1 1]);
    end

    % draw the memorization shape
    shape = task.fullSequence{task.trialnum}(task.thistrial.thisseg/2);
    switch shape
        case -1 % respond left, draw a black square
            mglPolygon(stimulus.squareX * task.thistrial.cue, stimulus.squareY, stimulus.squareColor);
        case 1 % respond right, draw a white circle
            mglMetalArcs([stimulus.symbolEcc * task.thistrial.cue; 0; 0], stimulus.circleColor, stimulus.circleSize, [0; 2*pi], 0);
        case 0 % filler, draw an outlined triangle
            mglPolygon(stimulus.triangleX * task.thistrial.cue, stimulus.triangleY, stimulus.triangleColor);
    end

    % draw the distractor shape
    distractor = task.distractorSequence{task.trialnum}(task.thistrial.thisseg/2);
    switch distractor
        case -1 % respond left, draw a black square
            mglPolygon(stimulus.squareX * -task.thistrial.cue, stimulus.squareY, stimulus.squareColor);
        case 1 % respond right, draw a white circle
            mglMetalArcs([stimulus.symbolEcc * -task.thistrial.cue; 0; 0], stimulus.circleColor, stimulus.circleSize, [0; 2*pi], 0);
        case 0 % filler, draw an outlined triangle
            mglPolygon(stimulus.triangleX * -task.thistrial.cue, stimulus.triangleY, stimulus.triangleColor);
    end

elseif ismember(task.thistrial.thisseg,3:2:stimulus.maxSequence*2+1) % blank segments

    % draw cue
    if task.thistrial.cue == -1
        % Vert
        mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, [0 0 0] );
        % left
        mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [1 1 1]);
        % right
        mglLines2( stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, [0 0 0]);
    elseif task.thistrial.cue == 1
        % Vert
        mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, [0 0 0] );
        % left
        mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [0 0 0]);
        % right
        mglLines2( stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, [1 1 1]);
    end

elseif task.thistrial.thisseg == stimulus.maxSequence*2+2 % memory period

    % draw fix
    % Vert
    mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, [0 0 0] );
    % left
    mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [0 0 0]);
    % right
    mglLines2( stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, [0 0 0]);

elseif task.thistrial.thisseg == stimulus.maxSequence*2+3 % response period

    % draw fix
    % Vert
    mglLines2(stimulus.cueVertX0,stimulus.cueVertY0,stimulus.cueVertX1,stimulus.cueVertY1, 2, stimulus.feedbackColors{1});
    % left
    mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, stimulus.feedbackColors{2});
    % right
    mglLines2( stimulus.cueRightX0,stimulus.cueRightY0,stimulus.cueRightX1,stimulus.cueRightY1, 2, stimulus.feedbackColors{3});

elseif task.thistrial.thisseg == stimulus.maxSequence*2+4 % ITI

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

if ~task.thistrial.noMoreResponses

    global stimulus

    fprintf('Response received : %g\n', task.thistrial.whichButton);
    task.thistrial.nResponses = task.thistrial.nResponses + 1;

    if ~any(task.thistrial.whichButton == [1 6]) % [1 6] should (hopefully) correspond to left/right buttons
        error('Check your button inputs')
    end

    % determine what the current correct button should be
    expectedCorrectResponse = task.correctSequence{task.trialnum}(task.thistrial.nResponses);
    if expectedCorrectResponse == -1
        correctbutton = 1;
    elseif expectedCorrectResponse == 1
        correctbutton = 6;
    elseif isnan(expectedCorrectResponse)
        task.sequenceResponse{task.trialnum}(task.thistrial.nResponses) = nan;
        corr = 0;
        % report answer
        fprintf('Trial %d. ++ Incorrect ++. \n',task.trialnum);
        task.thistrial.correct = corr;
        task.thistrial.nResponses = -length(task.correctSequence{task.trialnum});
        correctbutton = nan;
    end

    % see if it is correct
    if isequal(task.thistrial.whichButton,correctbutton) % correct
        task.sequenceResponse{task.trialnum}(task.thistrial.nResponses) = expectedCorrectResponse;
        corr = 1;
        % report answer
        fprintf('Trial %d. Response %d of %d !! Correct !!. \n',task.trialnum,task.thistrial.nResponses,length(task.correctSequence{task.trialnum}));
        task.thistrial.correct = corr;

    else % incorrect
        task.sequenceResponse{task.trialnum}(task.thistrial.nResponses) = -expectedCorrectResponse;
        corr = 0;
        % report answer
        fprintf('Trial %d. Response %d of %d ++ Incorrect ++. \n',task.trialnum,task.thistrial.nResponses,length(task.correctSequence{task.trialnum}));
        task.thistrial.correct = corr;

    end

end

% if inputted all responses, stop taking additional responses and let participant know
if abs(task.thistrial.nResponses) == length(task.correctSequence{task.trialnum})
    
    % see if all responses are correct
    if isequal(task.sequenceResponse{task.trialnum},task.correctSequence{task.trialnum})
        fprintf('Trial %d. All correct!! \n',task.trialnum);

        stimulus.feedbackColors{1} = [0 1 0];
        stimulus.feedbackColors{2} = [0 1 0];
        stimulus.feedbackColors{3} = [0 1 0];

    else
        fprintf('Trial %d. At least one error. \n',task.trialnum);

        stimulus.feedbackColors{1} = [1 0 0];
        stimulus.feedbackColors{2} = [1 0 0];
        stimulus.feedbackColors{3} = [1 0 0];
    end

    task.thistrial.noMoreResponses = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the shapes stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initShapes(stimulus,myscreen)

stimulus.pixRes = min(myscreen.screenHeight/myscreen.imageHeight, myscreen.screenWidth/myscreen.imageWidth);

stimulus.symbolEcc = 5;

% square parameters
stimulus.squareColor = [1 1 1];
stimulus.squareX = [3 7 7 3];
stimulus.squareY = [2 2 -2 -2];

% circle parameters
stimulus.circleColor = [1 1 1 1]';
stimulus.circleSize = [0; 2];

% triangle parameters
stimulus.triangleColor = [1 1 1];
stimulus.triangleX = [5 7 3];
stimulus.triangleY = [2 -2 -2];

% set initial colors for feedback trial
stimulus.feedbackColors{1} = [1 1 1];
stimulus.feedbackColors{2} = [1 1 1];
stimulus.feedbackColors{3} = [1 1 1];

% set background color
stimulus.backgroundColor = 0.5;


