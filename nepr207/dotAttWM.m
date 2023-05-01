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
task{1}{1}.getResponse = [0 0 0 1 0];
task{1}{1}.collectEyeData = false;
task{1}{1}.seglen = [0.5, 1, 8, 2.5, 3];

% task parameters
task{1}{1}.parameter.cue = [-1 -0.1 0.1 1];
task{1}{1}.parameter.probePosition = [-1 1];
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

% set dot parameters
stimulus.cuedotEcc = 5;
stimulus.cuedotX = [stimulus.cuedotEcc; -stimulus.cuedotEcc];
stimulus.cuedotY = [0; 0];
stimulus.dotSize = 0.5;
stimulus.dotColor = [1 1 1];
stimulus.probeDotOffset = 2;
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
        mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [0 0 0]);
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
        mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [0 0 0]);
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

    % draw the left/right dots
    mglGluDisk(stimulus.cuedotX,stimulus.cuedotY,stimulus.dotSize,stimulus.dotColor);
    
elseif task.thistrial.thisseg == 3

    % memory period

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
        mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [0 0 0]);
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

elseif task.thistrial.thisseg == 4

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
        mglLines2( stimulus.cueLeftX0,stimulus.cueLeftY0,stimulus.cueLeftX1,stimulus.cueLeftY1, 2, [0 0 0]);
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

    % draw probe
    cursign = sign(task.thistrial.cue);
    mglGluDisk(stimulus.cuedotEcc*cursign + (stimulus.probeDotOffset*task.thistrial.probePosition),stimulus.probeDotY,stimulus.dotSize,stimulus.dotColor);

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
% function to init the dots stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initDots(stimulus,myscreen)

% set background color
stimulus.backgroundColor = 0.5;

