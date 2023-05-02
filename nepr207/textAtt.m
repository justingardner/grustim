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

function myscreen = textAtt(varargin)

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

% task
task{1}{1}.getResponse = [zeros(1,11) 1 0];
task{1}{1}.collectEyeData = false;
task{1}{1}.segmin = [2, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 2, 1]; % average trial time: 13s
task{1}{1}.segmax = [2, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 2, 7];

% task parameters
task{1}{1}.parameter.cue = [-1 0 1];
task{1}{1}.parameter.real = [1 0];
task{1}{1}.parameter.wordPresent = [1 0];
task{1}{1}.randVars.calculated.trialType = nan;
task{1}{1}.randVars.calculated.distributedSide = nan;
task{1}{1}.random = 1;

task{1}{1}.numTrials = 24;

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

% set words
stimulus.realList = {'OTHER','ABOUT','WHICH','MAYBE','LUNCH','SERVE','SHARP','STAND','STONE','EAGER','EARTH','PIZZA'};
stimulus.fakeList = {'TOHER','OBTUA','HCIWH','YMBAE','UNLHC','EERSV','RPHSA','TSNDA','ONTSE','REAEG','HRATE','ZIZPA'};
stimulus.allList = [stimulus.realList stimulus.fakeList];
stimulus.matchWords = cell(task{1}{1}.numTrials,1);
stimulus.presentedWords = cell(task{1}{1}.numTrials,2);
stimulus.matchRealOrder = stimulus.realList(randperm(length(stimulus.realList)));
stimulus.matchFakeOrder = stimulus.fakeList(randperm(length(stimulus.fakeList)));
stimulus.realmatchCounter = 1;
stimulus.fakematchCounter = 1;

stimulus.wordEcc = 5;
stimulus.response = zeros(1,task{1}{1}.numTrials);
stimulus.correctResponse = zeros(1,task{1}{1}.numTrials);

% init the stimulus
myscreen = initStimulus('stimulus',myscreen);
stimulus = initWords(stimulus,myscreen);

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
    
    % pull a word from a predetermined random order
    if task.thistrial.real == 1
        stimulus.currentTrial.matchWord = stimulus.matchRealOrder{stimulus.realmatchCounter};
        stimulus.realmatchCounter = stimulus.realmatchCounter + 1;
    else
        stimulus.currentTrial.matchWord = stimulus.matchFakeOrder{stimulus.fakematchCounter};
        stimulus.fakematchCounter = stimulus.fakematchCounter + 1;
    end

    % set when the match will appear
    if task.thistrial.wordPresent == 1
        stimulus.presentSegment = randi(10,1);
    else
        stimulus.presentSegment = 0;
    end

    % set trialType
    if task.thistrial.cue == 0
        if task.thistrial.real == 1
            if task.thistrial.wordPresent == 1
                task.thistrial.trialType = 1;
            elseif task.thistrial.wordPresent == 0
                task.thistrial.trialType = 2;
            end

        elseif task.thistrial.real == 0
            if task.thistrial.wordPresent == 1
                task.thistrial.trialType = 3;
            elseif task.thistrial.wordPresent == 0
                task.thistrial.trialType = 4;
            end
        end

    elseif task.thistrial.cue == -1
        if task.thistrial.real == 1
            if task.thistrial.wordPresent == 1
                task.thistrial.trialType = 5;
            elseif task.thistrial.wordPresent == 0
                task.thistrial.trialType = 6;
            end

        elseif task.thistrial.real == 0
            if task.thistrial.wordPresent == 1
                task.thistrial.trialType = 7;
            elseif task.thistrial.wordPresent == 0
                task.thistrial.trialType = 8;
            end
        end

    elseif task.thistrial.cue == 1
        if task.thistrial.real == 1
            if task.thistrial.wordPresent == 1
                task.thistrial.trialType = 9;
            elseif task.thistrial.wordPresent == 0
                task.thistrial.trialType = 10;
            end

        elseif task.thistrial.real == 0
            if task.thistrial.wordPresent == 1
                task.thistrial.trialType = 11;
            elseif task.thistrial.wordPresent == 0
                task.thistrial.trialType = 12;
            end
        end

    end
    
    % set cue color
    if task.thistrial.cue == -1
        stimulus.cuecolor{1} = [1 1 1];
        fprintf('\nTrial %d: Cue left', task.trialnum)
    elseif task.thistrial.cue == 1
        stimulus.cuecolor{2} = [1 1 1];
        fprintf('\nTrial %d: Cue right', task.trialnum)
    else
        stimulus.cuecolor{1} = [0 0 0];
        stimulus.cuecolor{2} = [0 0 0];
        stimulus.cuecolor{3} = [0 0 0];
        stimulus.cuecolor{4} = [0 0 0];
        fprintf('\nTrial %d: Distributed cue', task.trialnum)
    end

elseif task.thistrial.thisseg > 1 && task.thistrial.thisseg < 12
    
    % set the left and right words
    currentTextLeft = stimulus.allList{randi(length(stimulus.allList),1)};
    currentTextRight = stimulus.allList{randi(length(stimulus.allList),1)};
    while strcmp(currentTextLeft,stimulus.currentTrial.matchWord) || strcmp(currentTextRight,stimulus.currentTrial.matchWord)
        currentTextLeft = stimulus.allList{randi(length(stimulus.allList),1)};
        currentTextRight = stimulus.allList{randi(length(stimulus.allList),1)};
    end
    stimulus.currentTextLeft = currentTextLeft;
    stimulus.currentTextRight = currentTextRight;

    if task.thistrial.thisseg == stimulus.presentSegment
        task.thistrial.distributedSide = round(-1 + 2*rand);

        fprintf('\nMatch Segment %d \n', stimulus.presentSegment)
    end

elseif task.thistrial.thisseg == 12
    stimulus.cuecolor{1} = [1 1 0];
    stimulus.cuecolor{2} = [1 1 0];
    stimulus.cuecolor{3} = [1 1 0];
    stimulus.cuecolor{4} = [1 1 0];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus

% clear screen to gray
mglClearScreen(stimulus.backgroundColor);

if task.thistrial.thisseg == 1

    % draw the text with flankers
    mglTextSet('Helvetica', 32, [0 0 0], 0, 0, 0, 0, 0, 0, 0);
    matchWord = mglText(stimulus.currentTrial.matchWord);
    mglBltTexture(matchWord,[0 3],'left','top');
    mglTextSet('Helvetica', 32, [0 0 0], 0, 0, 0, 0, 0, 0, 0);
    leftText = mglText('XXXXX');
    mglBltTexture(leftText,[-stimulus.wordEcc 3],'left','top');
    mglTextSet('Helvetica', 32, [0 0 0], 0, 0, 0, 0, 0, 0, 0);
    rightText = mglText('XXXXX');
    mglBltTexture(rightText,[stimulus.wordEcc 3],'left','top');

    mglFixationCrossDiag(0.5,2,stimulus.cuecolor);

elseif task.thistrial.thisseg > 1 && task.thistrial.thisseg < 12
    
    % draw fixation
    mglFixationCrossDiag(0.5,2)

    % draw the cue text
    mglTextSet('Helvetica', 32, [0 0 0], 0, 0, 0, 0, 0, 0, 0);
    matchWord = mglText(stimulus.currentTrial.matchWord);
    mglBltTexture(matchWord,[0 3],'left','top');

    % draw the left/right words
    if task.thistrial.thisseg == stimulus.presentSegment
        if task.thistrial.cue == 0
            mglTextSet('Helvetica', 32, [0 0 0], 0, 0, 0, 0, 0, 0, 0);
            matchWord = mglText(stimulus.currentTrial.matchWord);
            mglBltTexture(matchWord,[task.thistrial.distributedSide*stimulus.wordEcc 3],'left','top');
            mglTextSet('Helvetica', 32, [0 0 0], 0, 0, 0, 0, 0, 0, 0);
            matchWord = mglText(stimulus.currentTextLeft);
            mglBltTexture(matchWord,[task.thistrial.distributedSide*(-1)*stimulus.wordEcc 3],'left','top');

        else
            mglTextSet('Helvetica', 32, [0 0 0], 0, 0, 0, 0, 0, 0, 0);
            matchWord = mglText(stimulus.currentTrial.matchWord);
            mglBltTexture(matchWord,[task.thistrial.cue*stimulus.wordEcc 3],'left','top');
            mglTextSet('Helvetica', 32, [0 0 0], 0, 0, 0, 0, 0, 0, 0);
            currentText = mglText(stimulus.currentTextLeft);
            mglBltTexture(currentText,[task.thistrial.cue*(-1)*stimulus.wordEcc 3],'left','top');

        end
    
    else
        mglTextSet('Helvetica', 32, [0 0 0], 0, 0, 0, 0, 0, 0, 0);
        leftText = mglText(stimulus.currentTextLeft);
        mglBltTexture(leftText,[-stimulus.wordEcc 3],'left','top');
        mglTextSet('Helvetica', 32, [0 0 0], 0, 0, 0, 0, 0, 0, 0);
        rightText = mglText(stimulus.currentTextRight);
        mglBltTexture(rightText,[stimulus.wordEcc 3],'left','top');
    end
    
elseif task.thistrial.thisseg == 12
    
    % fixation cross in yellow
    mglFixationCrossDiag(0.5,2,stimulus.cuecolor);

    % print the cue
    mglTextSet('Helvetica', 32, [0 0 0], 0, 0, 0, 0, 0, 0, 0);
    matchWord = mglText(stimulus.currentTrial.matchWord);
    mglBltTexture(matchWord,[0 3],'left','top');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)

global stimulus

if task.thistrial.gotResponse < 1

    fprintf('Response received : %g\n', task.thistrial.whichButton);

    if ~any(task.thistrial.whichButton == [1 2]) % [1 2] should (hopefully) correspond to left/right buttons
        error('Check your button inputs')
    end
    
    % determine what the correct button should be
    if task.thistrial.wordPresent
        correctbutton = 1;
    elseif ~task.thistrial.wordPresent
        correctbutton = 2;
    end

    % see if it is correct
    if isequal(task.thistrial.whichButton,correctbutton) % correct
        corr = 1;
        % report answer
        fprintf('Trial %d. !! Correct !!. \n',task.trialnum);
        task.thistrial.correct = corr;
        
    else % incorrect
        corr = 0;
        % report answer
        fprintf('Trial %d. ++ Incorrect ++. \n',task.trialnum);
        task.thistrial.correct = corr;
        
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the words stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initWords(stimulus,myscreen)

stimulus.xWord = mglText('XXXXX');

% set background color
stimulus.backgroundColor = 0.5;

