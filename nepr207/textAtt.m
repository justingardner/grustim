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
task{1}{1}.seglen = [2, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 2, 1];

% task parameters
task{1}{1}.parameter.cue = [-1 1];
task{1}{1}.parameter.real = [1 0];
task{1}{1}.parameter.wordPresent = [1 0];
task{1}{1}.random = 1;

task{1}{1}.numTrials = 24;

% initialize the task
for phaseNum = 1:length(task)
    [task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end

% global stimulus parameters
global stimulus;

% if scanning, do synch to volume
task{1}{1}.synchToVol = zeros(size(task{1}{1}.seglen));
if atScanner
    task{1}{1}.fudgeLastVolume = 1;
    task{1}{1}.seglen(end) = task{1}{1}.seglen(end)-0.5;
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

stimulus.wordEcc = 10;
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
    mglBltTexture(matchWord,[0 0],'left','top');
    if task.thistrial.cue == -1
        mglTextSet('Helvetica', 32, [1 1 1], 0, 0, 0, 0, 0, 0, 0);
        leftText = mglText('XXXXX');
        mglBltTexture(leftText,[-stimulus.wordEcc 0],'left','top');
        mglTextSet('Helvetica', 32, [0 0 0], 0, 0, 0, 0, 0, 0, 0);
        rightText = mglText('XXXXX');
        mglBltTexture(rightText,[stimulus.wordEcc 0],'left','top');
    elseif task.thistrial.cue == 1
        mglTextSet('Helvetica', 32, [0 0 0], 0, 0, 0, 0, 0, 0, 0);
        leftText = mglText('XXXXX');
        mglBltTexture(leftText,[-stimulus.wordEcc 0],'left','top');
        mglTextSet('Helvetica', 32, [1 1 1], 0, 0, 0, 0, 0, 0, 0);
        rightText = mglText('XXXXX');
        mglBltTexture(rightText,[stimulus.wordEcc 0],'left','top');
    end

elseif task.thistrial.thisseg > 1 && task.thistrial.thisseg < 12
    
    % draw the cue text
    mglTextSet('Helvetica', 32, [0 0 0], 0, 0, 0, 0, 0, 0, 0);
    matchWord = mglText(stimulus.currentTrial.matchWord);
    mglBltTexture(matchWord,[0 0],'left','top');

    % draw the left/right words
    if task.thistrial.thisseg == stimulus.presentSegment
        mglTextSet('Helvetica', 32, [0 0 0], 0, 0, 0, 0, 0, 0, 0);
        matchWord = mglText(stimulus.currentTrial.matchWord);
        mglBltTexture(matchWord,[task.thistrial.cue*stimulus.wordEcc 0],'left','top');
        mglTextSet('Helvetica', 32, [0 0 0], 0, 0, 0, 0, 0, 0, 0);
        currentText = mglText(stimulus.currentTextLeft);
        mglBltTexture(currentText,[task.thistrial.cue*(-1)*stimulus.wordEcc 0],'left','top');
        fprintf('Match Segment %d \n\n', stimulus.presentSegment)
    else
        mglTextSet('Helvetica', 32, [0 0 0], 0, 0, 0, 0, 0, 0, 0);
        leftText = mglText(stimulus.currentTextLeft);
        mglBltTexture(leftText,[-stimulus.wordEcc 0],'left','top');
        mglTextSet('Helvetica', 32, [0 0 0], 0, 0, 0, 0, 0, 0, 0);
        rightText = mglText(stimulus.currentTextRight);
        mglBltTexture(rightText,[stimulus.wordEcc 0],'left','top');
    end
    
elseif task.thistrial.thisseg == 12

    % print the cue in yellow to indicate response period
    mglTextSet('Helvetica', 32, [1 1 0], 0, 0, 0, 0, 0, 0, 0);
    matchWord = mglText(stimulus.currentTrial.matchWord);
    mglBltTexture(matchWord,[0 0],'left','top');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)

global stimulus

fprintf('Response received : %g', task.thistrial.whichButton);

% check the response
% if task.thistrial.gotResponse < 1
%     % see if it is correct
%     if isequal(task.thistrial.whichButton,task.thistrial.whichSide)
%         % report answer
%         fprintf(' !! Correct !!. Reaction time: %0.2f',task.thistrial.reactionTime);
%         % change fixation color
%         stimulus.fixColor = stimulus.correctFixColor;
%         stimulus.correctResponse(task.trialnum) = 1;
%         if task.thistrial.whichButton == 1
%             stimulus.response(task.trialnum) = 0; % subject answered left
%         elseif task.thistrial.whichButton == 2
%             stimulus.response(task.trialnum) = 1; % subject answered right
%         end
%         % task = jumpSegment(task);
%     elseif ~task.thistrial.whichSide
%         % report answer
%         fprintf(' -- Same direction --. Reaction time: %0.2f',task.thistrial.reactionTime);
%         stimulus.fixColor = stimulus.neutralFixColor;
%         stimulus.correctResponse(task.trialnum) = nan;
%         if task.thistrial.whichButton == 1
%             stimulus.response(task.trialnum) = 0; % subject answered left
%         elseif task.thistrial.whichButton == 2
%             stimulus.response(task.trialnum) = 1; % subject answered right
%         end
%         % task = jumpSegment(task);
%     else
%         % report answer
%         fprintf(' ++ Incorrect ++. Reaction time: %0.2f',task.thistrial.reactionTime);
%         % change fixation color
%         stimulus.fixColor = stimulus.incorrectFixColor;
%         stimulus.correctResponse(task.trialnum) = 0;
%         if task.thistrial.whichButton == 1
%             stimulus.response(task.trialnum) = 0; % subject answered left
%         elseif task.thistrial.whichButton == 2
%             stimulus.response(task.trialnum) = 1; % subject answered right
%         end
%         % task = jumpSegment(task);
%     end
%     
% % no response
% else 
%     stimulus.correctResponse(task.trialnum) = nan;
%     stimulus.response(task.trialnum) = nan;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the words stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initWords(stimulus,myscreen)

stimulus.xWord = mglText('XXXXX');

% set background color
stimulus.backgroundColor = 0.5;

