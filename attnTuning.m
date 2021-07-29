% taskTemplate.m
%
%        $Id$
%      usage: taskTemplate
%         by: justin gardner
%       date: 09/07/06
%  copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%    purpose: example program to show how to use the task structure
%
function myscreen = attnTuning(varargin)
clear global stimulus

% get arguments
staircase = 0;
getArgs(varargin,{'staircase=0'},'verbose=1');

% initalize the screen
myscreen.background = 0.5;
myscreen = initScreen(myscreen);

% init the stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);
stimulus.staircase = staircase;

% fix: set waitForBacktick if you want to synch with the scanner
% by waiting for the backtick key to be pressed before starting the experiment
% (for systems that use NI digital I/O, this will wait for the digital
% signal that the scanner has started collecting data)
task{1}.waitForBacktick = 1;
% fix: the task defined here has two segments, one that
% is 3 seconds long followed by another that is 
% 6-9 seconds (randomized in steps of 1.5 seconds)
% change this to what you want for your trial
task{1}.seglen = [1.75 2 1.5 2 1.5 2 1.5];
% task{1}.segquant = [0 1.5];
% fix: set which segment(s) you want to collect subject responses for
task{1}.getResponse = [0 0 1 0 1 0 1];
% fix: enter the parameter of your choice
% task{1}.parameter.myParameter = [0 30 90];
task{1}.random = 1;
if ~stimulus.staircase
    task{1}.parameter.offset = [0, 3, 5, 7, 15, 30, 90]; % reporter grating offsets
end

attentionCond = [1,2];
flickerSide = [1,2]; % left or right
offset = [0, 3, 5, 7, 15, 30, 90]; % reporter grating offsets
targetOffset = [-5, -2.5, -1, -0.5, 0.5, 1, 2.5, 5];
nTargetOffsets = length(targetOffset);
designMat(:,1) = [repmat(1, [nTargetOffsets*2,1]); repmat(2, [nTargetOffsets*2,1])];
designMat(:,2) = [repmat(1,[nTargetOffsets,1]); repmat(2,[nTargetOffsets,1]); repmat(1,[nTargetOffsets,1]); repmat(2,[nTargetOffsets,1])];
if stimulus.staircase
    designMat(:,3) = [repmat(offset',[4,1])];
else
    designMat(:,3) = [repmat(targetOffset',[4,1])];
end
designMat = Shuffle(designMat);%, 2);
nTrials = length(designMat);

task{1}.randVars.attentionCond = designMat(:,1);
task{1}.randVars.flickerSide = designMat(:,2);
if stimulus.staircase
    task{1}.randVars.offset = designMat(:,3);
else
    task{1}.randVars.targetOffset = designMat(:,3);
end

% task{1}.numBlocks = 1;
task{1}.numTrials = nTrials;

task{1}.randVars.calculated.cue1 = nan; %left or right 1 or 2
task{1}.randVars.calculated.cue2 = nan;
task{1}.randVars.calculated.cue3 = nan;

task{1}.randVars.calculated.flickerResp1 = nan; % same or different 1 or 2
task{1}.randVars.calculated.flickerResp2 = nan;
task{1}.randVars.calculated.flickerResp3 = nan;

task{1}.randVars.calculated.resp1 = nan;
task{1}.randVars.calculated.rt1 = nan;
task{1}.randVars.calculated.correct1 = nan;
task{1}.randVars.calculated.resp2 = nan;
task{1}.randVars.calculated.rt2 = nan;
task{1}.randVars.calculated.correct2 = nan;
task{1}.randVars.calculated.resp3 = nan;
task{1}.randVars.calculated.rt3 = nan;
task{1}.randVars.calculated.correct2 = nan;

task{1}.randVars.calculated.targetOrientation1 = nan; % 1 or 2 counter-clockwise or clockwise
task{1}.randVars.calculated.targetOrientation2 = nan;
task{1}.randVars.calculated.targetOrientation3 = nan;

task{1}.randVars.calculated.nontargetOrientation1 = nan; % 1 or 2 counter-clockwise or clockwise
task{1}.randVars.calculated.nontargetOrientation2 = nan;
task{1}.randVars.calculated.nontargetOrientation3 = nan;

task{1}.randVars.calculated.targetOffset1 = nan;
task{1}.randVars.calculated.targetOffset2 = nan;
task{1}.randVars.calculated.targetOffset3 = nan;

task{1}.randVars.calculated.nontargetOffset1 = nan;
task{1}.randVars.calculated.nontargetOffset2 = nan;
task{1}.randVars.calculated.nontargetOffset3 = nan;


% initialize the task
for phaseNum = 1:length(task)
  [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end


% Setup parameters for stimulus
display.dist = myscreen.displayDistance;
display.width = myscreen.displaySize(1);
display.resolution = myscreen.screenWidth;
% pix = angle2pix(display,1);
% stimulus.width = 400 / pix; % gabor
% stimulus.height = 6`2200 / pix;
stimulus.sf = 1;
stimulus.orientation1 = 90;
stimulus.orientation2 = 60;
stimulus.orientation3 = 90;
stimulus.orientation4 = 0;

stimulus.leftPos = -7.5;
stimulus.rightPos = 7.5;

stimulus.annulusInner = 2.5;
stimulus.annulusOuter = 5;
stimulus.gap = 0.25;
stimulus.targetSize = stimulus.annulusInner-stimulus.gap; %2.5-0.25 = 2.25

stimulus.framesPerSecond = 120;%myscreen.framesPerSecond;
stimulus.flickerDur = 1.75; 
stimulus.flickerFrames = stimulus.flickerDur * stimulus.framesPerSecond;

stimulus.framesPerCycle1 = 35; %210/6
stimulus.framesPerCycle2 = 30; %210/7
stimulus.f1 = stimulus.framesPerSecond / stimulus.framesPerCycle1;
stimulus.f2 = stimulus.framesPerSecond / stimulus.framesPerCycle2;

stimulus.stimArray1 = makeStimArray(stimulus.flickerFrames, stimulus.framesPerCycle1);
stimulus.stimArray2 = makeStimArray(stimulus.flickerFrames, stimulus.framesPerCycle2);

stimulus.cueSize = 0.75; % length
stimulus.cueLineWidth = 2;
stimulus.preCueColor = [1 1 1];
stimulus.postCueColor = [0 1 1];
stimulus.leftCueStart = -0.75;
stimulus.leftCueEnd = stimulus.leftCueStart - stimulus.cueSize;
stimulus.rightCueStart = 0.75;
stimulus.rightCueEnd = stimulus.rightCueStart + stimulus.cueSize;

stimulus.fixColor = [1 1 1];
stimulus.fixWhite = [1 1 1];
stimulus.fixCorrect = [0 1 0];
stimulus.fixIncorrect = [1 0 0];


% fix: you will change the funciton myInitStimulus
% to initialize the stimulus for your experiment.
% stimulus = myInitStimulus(stimulus,myscreen);
stimulus = initGrating(stimulus,myscreen);

mglStencilCreateBegin(1);
mglFillOval(stimulus.leftPos,0, [stimulus.annulusOuter stimulus.annulusOuter]);
mglFillOval(stimulus.rightPos, 0, [stimulus.annulusOuter stimulus.annulusOuter]);
% mglFillOval(0,0, [1 1]);
mglStencilCreateEnd;

mglStencilCreateBegin(2);
mglFillOval(stimulus.leftPos,0, [stimulus.targetSize stimulus.targetSize]);
mglFillOval(stimulus.rightPos, 0, [stimulus.targetSize stimulus.targetSize]);
mglStencilCreateEnd;
mglClearScreen;

if stimulus.staircase
stimulus.initialThreshold = 10;
stimulus.initialStepsize = 2.5;
stimulus.minThreshold = 0;
stimulus.maxThreshold = 30;
stimulus.minStepsize = 0.75;
stimulus.maxStepsize = 5;

stimulus = initStair(stimulus);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% myscreen = eyeCalibDisp(myscreen);

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

stimulus.grating1 = []; stimulus.grating2 = []; stimulus.grating3 = []; stimulus.grating4 = [];
stimulus.gratingSuperimposed1 = []; stimulus.gratingSuperimposed2 = []; 
stimulus.texSuperimposed1 = []; stimulus.texSuperimposed2 = []; stimulus.texTarget = []; stimulus.texTargetOffset=[]; stimulus.texNontargetOffset=[];

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;
if task.thistrial.thisseg == 1

    stimulus.fixColor = stimulus.fixWhite;
    stimulus.counter1 = 1;
    stimulus.counter2 = 1;
    
    task.thistrial.targetOrientation1 = randi(2);
    task.thistrial.targetOrientation2 = randi(2);
    task.thistrial.targetOrientation3 = randi(2);
    
    task.thistrial.nontargetOrientation1 = randi(2);
    task.thistrial.nontargetOrientation2 = randi(2);
    task.thistrial.nontargetOrientation3 = randi(2);
    
    task.thistrial.cue1 = randi(2);
    task.thistrial.cue2 = randi(2);
    task.thistrial.cue3 = randi(2);
    
    if task.thistrial.flickerSide == task.thistrial.cue1
        task.thistrial.flickerResp1 = 1;
    else
        task.thistrial.flickerResp1 = 2;
    end
    if task.thistrial.flickerSide == task.thistrial.cue2
        task.thistrial.flickerResp2 = 1;
    else
        task.thistrial.flickerResp2 = 2;
    end
    if task.thistrial.flickerSide == task.thistrial.cue3
        task.thistrial.flickerResp3 = 1;
    else
        task.thistrial.flickerResp3 = 2;
    end
    
    stimulus = initFlicker(stimulus,task,myscreen);
    
    if task.thistrial.flickerSide == 1 % left
        stimulus.flickerPos = stimulus.leftPos;
        stimulus.fixPos = stimulus.rightPos;
    else
        stimulus.flickerPos = stimulus.rightPos;
        stimulus.fixPos = stimulus.leftPos;
    end
    
    if stimulus.staircase
    [testValue, stimulus.stair{task.thistrial.attentionCond}{task.thistrial.flickerResp1}] = doStaircase('testValue', stimulus.stair{task.thistrial.attentionCond}{task.thistrial.flickerResp1});
    
    if task.thistrial.targetOrientation1 == 1
        sign_target = 1; % counter-clockwise
    else
        sign_target = -1; % clockwise
    end
    if task.thistrial.nontargetOrientation1 == 1
        sign_nontarget = 1; % counter-clockwise
    else
        sign_nontarget = -1; % clockwise
    end
    
    task.thistrial.targetOffset1 = testValue * sign_target;
    task.thistrial.nontargetOffset1 = testValue * sign_nontarget;
    
    else
        if task.thistrial.nontargetOrientation1 == 1
            sign_nontarget = 1; % counter-clockwise
        else
            sign_nontarget = -1; % clockwise
        end
        task.thistrial.targetOffset1 = task.thistrial.targetOffset;
        task.thistrial.nontargetOffset1 = task.thistrial.targetOffset * sign_nontarget;
        
    end
    stimulus.targetOffset = task.thistrial.targetOffset1;
    stimulus.nontargetOffset = task.thistrial.nontargetOffset1;
    
    stimulus = initTarget(stimulus,task,myscreen);
    
    disp(sprintf('Trial %i: target %0.2f reporter %0.2f ',task.trialnum, task.thistrial.targetOffset, task.thistrial.offset));
    
elseif task.thistrial.thisseg == 2
    stimulus.fixColor = stimulus.fixWhite;
    stimulus.tic = tic;
    if task.thistrial.cue1 == 1 % cue left
        stimulus.targetPos = stimulus.leftPos;
        stimulus.nontargetPos = stimulus.rightPos;
    else
        stimulus.targetPos = stimulus.rightPos;
        stimulus.nontargetPos = stimulus.leftPos;
    end
    
    if task.thistrial.cue1 == 1; %left
        stimulus.thisPostCueStart = stimulus.leftCueStart;
        stimulus.thisPostCueEnd = stimulus.leftCueEnd;
    else
        stimulus.thisPostCueStart = stimulus.rightCueStart;
        stimulus.thisPostCueEnd = stimulus.rightCueEnd;
    end
        
elseif task.thistrial.thisseg == 4
    stimulus.fixColor = stimulus.fixWhite;
    stimulus.tic = tic;
    
    if task.thistrial.cue2 == 1 % cue left
        stimulus.targetPos = stimulus.leftPos;
        stimulus.nontargetPos = stimulus.rightPos;
    else
        stimulus.targetPos = stimulus.rightPos;
        stimulus.nontargetPos = stimulus.leftPos;
    end
    
    if stimulus.staircase
    [testValue, stimulus.stair{task.thistrial.attentionCond}{task.thistrial.flickerResp2}] = doStaircase('testValue', stimulus.stair{task.thistrial.attentionCond}{task.thistrial.flickerResp2});
    
    if task.thistrial.targetOrientation2 == 1
        sign_target = 1; % counter-clockwise
    else
        sign_target = -1; % clockwise
    end
    if task.thistrial.nontargetOrientation2 == 1
        sign_nontarget = 1; % counter-clockwise
    else
        sign_nontarget = -1; % clockwise
    end
    
    task.thistrial.targetOffset2 = testValue * sign_target;
    task.thistrial.nontargetOffset2 = testValue * sign_nontarget;
    
    else
        if task.thistrial.nontargetOrientation2 == 1
            sign_nontarget = 1; % counter-clockwise
        else
            sign_nontarget = -1; % clockwise
        end
        task.thistrial.targetOffset2 = task.thistrial.targetOffset;
        task.thistrial.nontargetOffset2 = task.thistrial.targetOffset * sign_nontarget;
    end
    
    stimulus.targetOffset = task.thistrial.targetOffset2;
    stimulus.nontargetOffset = task.thistrial.nontargetOffset2;
    
    stimulus = initTarget(stimulus,task,myscreen);
    
    if task.thistrial.cue1 == 2; %left
        stimulus.thisPostCueStart = stimulus.leftCueStart;
        stimulus.thisPostCueEnd = stimulus.leftCueEnd;
    else
        stimulus.thisPostCueStart = stimulus.rightCueStart;
        stimulus.thisPostCueEnd = stimulus.rightCueEnd;
    end
    
elseif task.thistrial.thisseg == 6
    stimulus.fixColor = stimulus.fixWhite;
    stimulus.tic = tic;
    
    if task.thistrial.cue3 == 1 % cue left
        stimulus.targetPos = stimulus.leftPos;
        stimulus.nontargetPos = stimulus.rightPos;
    else % cue right
        stimulus.targetPos = stimulus.rightPos;
        stimulus.nontargetPos = stimulus.leftPos;
    end
    
    if stimulus.staircase
    [testValue, stimulus.stair{task.thistrial.attentionCond}{task.thistrial.flickerResp3}] = doStaircase('testValue', stimulus.stair{task.thistrial.attentionCond}{task.thistrial.flickerResp3});
    
    if task.thistrial.targetOrientation3 == 1
        sign_target = 1; % counter-clockwise
    else
        sign_target = -1; % clockwise
    end
    if task.thistrial.nontargetOrientation3 == 1
        sign_nontarget = 1; % counter-clockwise
    else
        sign_nontarget = -1; % clockwise
    end
    
    task.thistrial.targetOffset3 = testValue * sign_target;
    task.thistrial.nontargetOffset3 = testValue * sign_nontarget;
    else
        if task.thistrial.nontargetOrientation3 == 1
            sign_nontarget = 1; % counter-clockwise
        else
            sign_nontarget = -1; % clockwise
        end
        task.thistrial.targetOffset3 = task.thistrial.targetOffset;
        task.thistrial.nontargetOffset3 = task.thistrial.targetOffset * sign_nontarget;
    end
    stimulus.targetOffset = task.thistrial.targetOffset3;
    stimulus.nontargetOffset = task.thistrial.nontargetOffset3;
    
    stimulus = initTarget(stimulus,task,myscreen);
    
    if task.thistrial.cue3 == 1; %left
        stimulus.thisPostCueStart = stimulus.leftCueStart;
        stimulus.thisPostCueEnd = stimulus.leftCueEnd;
    else
        stimulus.thisPostCueStart = stimulus.rightCueStart;
        stimulus.thisPostCueEnd = stimulus.rightCueEnd;
    end
end
% fix: do anything that needs to be done at the beginning
% of a segment (like for example setting the stimulus correctly
% according to the parameters etc).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus

% fix: display your stimulus here, for this code we just display 
% a fixation cross that changes color depending on the segment
% we are on.
if task.thistrial.thisseg == 1

% mglClearScreen;
mglStencilSelect(1);
if stimulus.counter1 <= 210
    mglBltTexture(stimulus.texSuperimposed1(stimulus.counter1), [stimulus.flickerPos 0]);
end
mglBltTexture(stimulus.texSuperimposed2, [stimulus.fixPos 0]);
stimulus.counter1 = stimulus.counter1 + 1;

mglFillOval(stimulus.leftPos, 0, [stimulus.annulusInner stimulus.annulusInner], [0.5 0.5 0.5]);
mglFillOval(stimulus.rightPos, 0, [stimulus.annulusInner stimulus.annulusInner], [0.5 0.5 0.5]);

mglStencilSelect(2);
mglBltTexture(stimulus.texTarget, [stimulus.leftPos 0]);
mglBltTexture(stimulus.texTarget, [stimulus.rightPos 0]);

mglStencilSelect(0);

    mglFixationCross(1,1,stimulus.fixColor);
%   mglFixationCross(1,1,[0 1 1]);

% mglStencilSelect(0);
else
    mglClearScreen;
    mglStencilSelect(1);
%     if stimulus.counter1 <= 210
        mglBltTexture(stimulus.texSuperimposed1(stimulus.counter2), [stimulus.flickerPos 0]);
%     end
    mglBltTexture(stimulus.texSuperimposed2, [stimulus.fixPos 0]);
    if stimulus.counter2 == 210
        stimulus.counter2 = 0;
    end
    stimulus.counter2 = stimulus.counter2 + 1;
    
    mglFillOval(stimulus.leftPos, 0, [stimulus.annulusInner stimulus.annulusInner], [0.5 0.5 0.5]);
    mglFillOval(stimulus.rightPos, 0, [stimulus.annulusInner stimulus.annulusInner], [0.5 0.5 0.5]);

    mglStencilSelect(2);
    mglBltTexture(stimulus.texTarget, [stimulus.leftPos 0]);
    mglBltTexture(stimulus.texTarget, [stimulus.rightPos 0]);

    mglStencilSelect(0);
    
    mglFixationCross(1,1,stimulus.fixColor);
    
    if any(task.thistrial.thisseg == [2, 4, 6]) % precue + target
        if toc(stimulus.tic) >= 1 && toc(stimulus.tic) <= 1.2
            mglStencilSelect(2);
            mglBltTexture(stimulus.texTargetOffset, [stimulus.targetPos 0]);
            mglBltTexture(stimulus.texNontargetOffset, [stimulus.nontargetPos 0]);

            mglStencilSelect(0);
        end
        
        if task.thistrial.attentionCond == 1 % focal cue
            %mglLines(x0, y0, x1, y1,size,color)
            mglLines2(stimulus.thisPostCueStart, 0, stimulus.thisPostCueEnd, 0, stimulus.cueLineWidth, stimulus.preCueColor);
        else % divided cue
            mglLines2(stimulus.leftCueStart, 0, stimulus.leftCueEnd, 0, stimulus.cueLineWidth, stimulus.preCueColor);
            mglLines2(stimulus.rightCueStart, 0, stimulus.rightCueEnd, 0, stimulus.cueLineWidth, stimulus.preCueColor);
        end
                

    elseif any(task.thistrial.thisseg == [3, 5, 7]) % postcue + response
        
        % postcue
        mglLines2(stimulus.thisPostCueStart, 0, stimulus.thisPostCueEnd, 0, stimulus.cueLineWidth, stimulus.postCueColor);
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)

global stimulus

% fix: add the code you want to use to process subject response

% here, we just check whether this is the first time we got a response
% this trial and display what the subject's response was and the reaction time
% if task.thistrial.gotResponse < 1
    if task.thistrial.thisseg == 3
        if (task.thistrial.targetOffset1 > 0 && task.thistrial.whichButton == 1) || (task.thistrial.targetOffset1 < 0 && task.thistrial.whichButton == 2)
            task.thistrial.correct1 = 1;
            stimulus.fixColor = stimulus.fixCorrect;
        else
            task.thistrial.correct1 = 0;
            stimulus.fixColor = stimulus.fixIncorrect;
        end
        if stimulus.staircase
        stimulus.stair{task.thistrial.attentionCond}{task.thistrial.flickerResp1} = doStaircase('update', stimulus.stair{task.thistrial.attentionCond}{task.thistrial.flickerResp1},...
            task.thistrial.correct1);
        end
        
        task.thistrial.resp1 = task.thistrial.whichButton;
        task.thistrial.rt1 = task.thistrial.reactionTime;
        disp(sprintf('Subject response: %i Reaction time: %0.2fs',task.thistrial.whichButton,task.thistrial.reactionTime));
        
        task.thistrial.gotResponse = 0;
        task.thistrial.whichButton = nan;
        task.thistrial.reactionTime = nan;
        
    elseif task.thistrial.thisseg == 5
        if (task.thistrial.targetOffset2 > 0 && task.thistrial.whichButton == 1) || (task.thistrial.targetOffset2 < 0 && task.thistrial.whichButton == 2)
            task.thistrial.correct2 = 1;
            stimulus.fixColor = stimulus.fixCorrect;
        else
            task.thistrial.correct2 = 0;
            stimulus.fixColor = stimulus.fixIncorrect;
        end
        if stimulus.staircase
        stimulus.stair{task.thistrial.attentionCond}{task.thistrial.flickerResp2} = doStaircase('update', stimulus.stair{task.thistrial.attentionCond}{task.thistrial.flickerResp2},...
            task.thistrial.correct2);
        end
        
        task.thistrial.resp2 = task.thistrial.whichButton;
        task.thistrial.rt2 = task.thistrial.reactionTime;
        disp(sprintf('Subject response: %i Reaction time: %0.2fs',task.thistrial.whichButton,task.thistrial.reactionTime));
        
        task.thistrial.gotResponse = 0;
        task.thistrial.whichButton = nan;
        task.thistrial.reactionTime = nan;
        
    elseif task.thistrial.thisseg == 7
        if (task.thistrial.targetOffset3 > 0 && task.thistrial.whichButton == 1) || (task.thistrial.targetOffset3 < 0 && task.thistrial.whichButton == 2)
            task.thistrial.correct3 = 1;
            stimulus.fixColor = stimulus.fixCorrect;
        else
            task.thistrial.correct3 = 0;
            stimulus.fixColor = stimulus.fixIncorrect;
        end
        if stimulus.staircase
        stimulus.stair{task.thistrial.attentionCond}{task.thistrial.flickerResp3} = doStaircase('update', stimulus.stair{task.thistrial.attentionCond}{task.thistrial.flickerResp3},...
            task.thistrial.correct3);
        end
        
        task.thistrial.resp3 = task.thistrial.whichButton;
        task.thistrial.rt3 = task.thistrial.reactionTime;
        disp(sprintf('Subject response: %i Reaction time: %0.2fs',task.thistrial.whichButton,task.thistrial.reactionTime));
        
        task.thistrial.gotResponse = 0;
        task.thistrial.whichButton = nan;
        task.thistrial.reactionTime = nan;
    end
  
    
%     disp(sprintf('Subject response: %i Reaction time: %0.2fs',task.thistrial.whichButton,task.thistrial.reactionTime));
%     task.thistrial.resp = task.thistrial.whichButton;
%     task.thistrial.rt = task.thistrial.reactionTime;
  
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulusArray = makeStimArray(totalFrames,framesPerCycle)
stimulusArray = zeros(1,totalFrames);
rep = totalFrames/framesPerCycle;

if isodd(framesPerCycle)
    firsthalf = ceil(framesPerCycle/2);
    secondhalf = floor(framesPerCycle/2);
else
    firsthalf = framesPerCycle/2;
    secondhalf = firsthalf;
end

cycle = [repmat(1, [1,firsthalf]), repmat(2, [1,secondhalf])];
stimulusArray = repmat(cycle, [1, rep]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initGrating(stimulus,myscreen)
stimulus.grating1 = []; stimulus.grating2 = [];
stimulus.grating1{1} = mglMakeGrating(stimulus.annulusOuter, stimulus.annulusOuter, stimulus.sf, stimulus.orientation1,0);
stimulus.grating1{2} = mglMakeGrating(stimulus.annulusOuter, stimulus.annulusOuter, stimulus.sf, stimulus.orientation1,180);

stimulus.grating1{1} = 255*(stimulus.grating1{1}+1)/2;
stimulus.grating1{2} = 255*(stimulus.grating1{2}+1)/2;

stimulus.grating3 = mglMakeGrating(stimulus.annulusOuter, stimulus.annulusOuter, stimulus.sf, stimulus.orientation3,0);
stimulus.grating4 = mglMakeGrating(stimulus.annulusOuter, stimulus.annulusOuter, stimulus.sf, stimulus.orientation4,0);
stimulus.grating3 = 255*(stimulus.grating3+1)/2;
stimulus.grating4 = 255*(stimulus.grating4+1)/2;

stimulus.gratingSuperimposed2 = (stimulus.grating3 + stimulus.grating4) /2;

% stimulus.texSuperimposed(1) = mglCreateTexture(stimulus.gratingSuperimposed1);
stimulus.texSuperimposed2 = mglCreateTexture(stimulus.gratingSuperimposed2);

% fix: add stuff to initalize your stimulus
function stimulus = initFlicker(stimulus,task,myscreen)

% compute the grating
% reporter annuli
% mglMakeGrating(width,height,sf,angle,phase,<xDeg2pix>,<yDeg2pix>)

stimulus.grating2{1} = mglMakeGrating(stimulus.annulusOuter, stimulus.annulusOuter, stimulus.sf, stimulus.orientation1-task.thistrial.offset,0);
stimulus.grating2{2} = mglMakeGrating(stimulus.annulusOuter, stimulus.annulusOuter, stimulus.sf, stimulus.orientation1-task.thistrial.offset,180);

stimulus.grating2{1} = 255*(stimulus.grating2{1}+1)/2;
stimulus.grating2{2} = 255*(stimulus.grating2{2}+1)/2;

stimulus.gratingSuperimposed1 = []; 
for i = 1:stimulus.flickerFrames
    stimulus.gratingSuperimposed1{i} = (stimulus.grating1{stimulus.stimArray1(i)} + stimulus.grating2{stimulus.stimArray2(i)}) / 2;
    stimulus.texSuperimposed1(i) = mglCreateTexture(stimulus.gratingSuperimposed1{i});
end

% stimulus.gratingSuperimposed1 = (stimulus.grating1 + stimulus.grating2) /2;

function stimulus = initTarget(stimulus,task,myscreen)

stimulus.gratingTarget = mglMakeGrating(stimulus.targetSize, stimulus.targetSize, stimulus.sf, 90, 0);
stimulus.gratingTarget = 255*(stimulus.gratingTarget+1)/2;

stimulus.texTarget = mglCreateTexture(stimulus.gratingTarget);

targetOffset = mglMakeGrating(stimulus.targetSize, stimulus.targetSize, stimulus.sf, 90 + stimulus.targetOffset, 0);
nontargetOffset = mglMakeGrating(stimulus.targetSize, stimulus.targetSize, stimulus.sf, 90 + stimulus.nontargetOffset, 0);
targetOffset = 255*(targetOffset+1)/2;
nontargetOffset = 255*(nontargetOffset+1)/2;

stimulus.texTargetOffset = mglCreateTexture(targetOffset);
stimulus.texNontargetOffset = mglCreateTexture(nontargetOffset);


function stimulus = initStair(stimulus)
for attCond = 1:2
    for flickerCond = 1:2
        stimulus.stair{attCond}{flickerCond} = doStaircase('init', 'upDown', 'nup=1', 'ndown=2',...
            'initialThreshold', stimulus.initialThreshold, 'initialStepsize', stimulus.initialStepsize, ...
            'minStepsize', stimulus.minStepsize, 'maxStepsize', stimulus.maxStepsize, 'minThreshold', stimulus.minThreshold, 'maxThreshold', stimulus.maxThreshold,...
            'stepRule=pest','dispFig=1');
    end
end
% if stimulus.auditoryTrain || stimulus.visualTrain
% 	 stimulus.stair = doStaircase('init','upDown', 'nup=1','ndown=2',...
%         'initialThreshold', stimulus.initialThreshold, 'initialStepsize',stimulus.initialStepsize, ...
%     'minStepsize',stimulus.minStepsize,'maxStepsize',stimulus.maxStepsize,'minThreshold',stimulus.minThreshold,'maxThreshold', stimulus.maxThreshold,...
%      'stepRule=pest','dispFig=1');
% else

% 	condNames = {'vision','auditory','noOffset','posOffset','negOffset'};
%   relNames = {'high','low'};
% for rel = 1:2
% for cond = 1:5
% 	stimulus.stair{rel}{cond} = doStaircase('init','quest', 'initialThreshold', stimulus.initOffset, 'initialThresholdSd', stimulus.initOffsetSd, ...
% 		'pThreshold', 0.75,'dispFig=1','subplotRows=5','subplotCols=2','subplotNum',2*(cond-1)+rel,'subplotName',[relNames{rel},' ', condNames{cond}]);
% end
% end
   
% end