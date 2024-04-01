% fixStairInitTask.m
%
%        $Id$
%      usage: [fixTask myscreen] = fixStairInitTask(myscreen)
%         by: justin gardner
%       date: 09/07/06
%  copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%    purpose: Implements a fixation task. In this task, the fixation cross
%             starts out cyan and then darkens twice. The fixation cross then
%             turns yellow to indicate the response interval, and the subject
%             is required to press 1 or 2 to indicate in which interval the cross
%             appeared darker. The cross will then turn green or red to indicate
%             correct or incorrect responses. The dimness of the target is
%             controlled by a 2 down 1 up staircase to keep task difficulty
%             the same.
%
%             See testExperiment.m for how this is used in a task. If you want
%             to change parameters, before you call fixStairInitTask, set
%             appropriate fields of the global variable fixStimulus. e.g.:
%
%             global fixStimulus
%             fixStimulus.interTime = 1;
%
%             See the code, for a list of all parameters that can be changed.
%
function [task myscreen] = fixStairInitTaskLR_RT(myscreen)

% check arguments
if ~any(nargin == [1])
    help fixDispStairInitTask
    return
end

% create the stimulus for the experiment, use defaults if they are
% not already set
global fixStimulus;
myscreen = initStimulus('fixStimulus',myscreen);
if ~isfield(fixStimulus,'threshold') fixStimulus.threshold = 1; end
if ~isfield(fixStimulus,'pedestal') fixStimulus.pedestal = 0.4; end
if ~isfield(fixStimulus,'stairUp') fixStimulus.stairUp = 1; end
if ~isfield(fixStimulus,'stairDown') fixStimulus.stairDown = 2; end
if ~isfield(fixStimulus,'stairStepSize') fixStimulus.stairStepSize = 0.05; end
if ~isfield(fixStimulus,'stairUseLevitt') fixStimulus.stairUseLevitt = 0; end
if ~isfield(fixStimulus,'stimColor') fixStimulus.stimColor = [0 1 1]; end
if ~isfield(fixStimulus,'responseColor') fixStimulus.responseColor = [1 1 0]; end
if ~isfield(fixStimulus,'interColor') fixStimulus.interColor = [0 1 1]; end
if ~isfield(fixStimulus,'correctColor') fixStimulus.correctColor = [0 0.8 0]; end
if ~isfield(fixStimulus,'incorrectColor') fixStimulus.incorrectColor = [0.8 0 0]; end
if ~isfield(fixStimulus,'responseTime') fixStimulus.responseTime = 1; end
if ~isfield(fixStimulus,'stimTime') fixStimulus.stimTime = 0.4; end
if ~isfield(fixStimulus,'interTime') fixStimulus.interTime = 0.8; end
if ~isfield(fixStimulus,'diskSize') fixStimulus.diskSize = 1; end
if ~isfield(fixStimulus,'pos') fixStimulus.pos = [0 0]; end
if ~isfield(fixStimulus,'fixWidth') fixStimulus.fixWidth = 1; end
if ~isfield(fixStimulus,'fixLineWidth') fixStimulus.fixLineWidth = 0.3; end
if ~isfield(fixStimulus,'trainingMode') fixStimulus.trainingMode = 0;end
if ~isfield(fixStimulus,'verbose') fixStimulus.verbose = 1;end

% for trainingMode set text
if fixStimulus.trainingMode
    mglTextSet('Helvetica',64,[1 1 1],0,0,0,0,0,0,0);
end

% create a fixation task
task{1}.seglen = [fixStimulus.interTime, fixStimulus.stimTime+fixStimulus.responseTime, fixStimulus.interTime];
task{1}.getResponse = [0 1 0]; %[0 0 0 0 0 1];
[task{1} myscreen] = addTraces(task{1}, myscreen, 'segment', 'phase', 'response');

fixStimulus.leftArmStartX = fixStimulus.pos(1) - fixStimulus.fixWidth/2;
fixStimulus.leftArmEndX = fixStimulus.pos(1);
fixStimulus.leftArmStartY = fixStimulus.pos(2);
fixStimulus.leftArmEndY = fixStimulus.pos(2);

fixStimulus.rightArmStartX = fixStimulus.pos(1);
fixStimulus.rightArmEndX = fixStimulus.pos(1) + fixStimulus.fixWidth/2;
fixStimulus.rightArmStartY = fixStimulus.pos(2);
fixStimulus.rightArmEndY = fixStimulus.pos(2);

% init the staircase
fixStimulus.staircase = upDownStaircase(fixStimulus.stairUp,fixStimulus.stairDown,fixStimulus.threshold,fixStimulus.stairStepSize,fixStimulus.stairUseLevitt);
fixStimulus.staircase.minThreshold = 0;
fixStimulus.staircase.maxThreshold = 1;

% init the task
[task{1} myscreen] = initTask(task{1},myscreen,@fixStartSegmentCallback,@fixDrawStimulusCallback,@fixTrialResponseCallback,@fixTrialStartCallback);

[task{1} myscreen] = addTraces(task{1}, myscreen, 'fixStair');

% keep the correct and incorrect counts
task{1}.correct = 0;
task{1}.incorrect = 0;

% keep track of which side signal was presented on
task{1}.signalLR = [];
% keep track of subject response
task{1}.correctResponse = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each trial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = fixTrialStartCallback(task, myscreen)

global fixStimulus;
% choose stimulus interval
task.thistrial.sigLocation = 1+(rand > 0.5); % 1 for left, 2 for right
task.signalLR = [task.signalLR task.thistrial.sigLocation]; % record signal location
if fixStimulus.verbose
    disp(sprintf('sigint = %i threshold = %0.2f',task.thistrial.sigLocation,fixStimulus.threshold));
end

% update staircase from last trial's response (or lack thereof)
if task.trialnum > 1
    if task.correctResponse(end) > 0
        response = 1;
    else
        response = 0;
    end
    fixStimulus.staircase = upDownStaircase(fixStimulus.staircase,response);
    fixStimulus.threshold = fixStimulus.staircase.threshold;
end

% update the correct response array such that the current trial is coded as no response (-1)
task.correctResponse = [task.correctResponse -1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = fixStartSegmentCallback(task, myscreen)

global fixStimulus;

if isfield(fixStimulus,'displayText')
    % delete the texture
    if ~isempty(fixStimulus.displayText)
        mglDeleteTexture(fixStimulus.displayText);
    end
end
fixStimulus.displayText = [];

% choose what color the arm of the fixation cross will be
if task.thistrial.thisseg == 1
    fixStimulus.thisColor = fixStimulus.interColor; % reset from previous trial
elseif task.thistrial.thisseg == 2 % change
    fixStimulus.thisColor = fixStimulus.interColor;
    fixStimulus.thisArmColor = fixStimulus.interColor;
    fixStimulus.thisArmColor(1) = fixStimulus.threshold;
else
    fixStimulus.thisColor = fixStimulus.responseColor;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called every frame udpate to draw the fixation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = fixDrawStimulusCallback(task, myscreen)

global fixStimulus;

if fixStimulus.trainingMode,mglClearScreen;end

if ~isempty(fixStimulus.displayText)
    mglBltTexture(fixStimulus.displayText,fixStimulus.displayTextLoc);
end

% mglGluDisk(fixStimulus.pos(1),fixStimulus.pos(2),fixStimulus.diskSize*[1 1],myscreen.background,60);
mglPoints2(fixStimulus.pos(1),fixStimulus.pos(2),fixStimulus.diskSize*[1 1],myscreen.background, 1);

drawFixationCross(fixStimulus.fixWidth,fixStimulus.fixLineWidth,fixStimulus.thisColor',fixStimulus.pos);

if task.thistrial.thisseg == 2 && task.thistrial.gotResponse == 0
    if task.thistrial.sigLocation == 1
        % mglLines2(fixStimulus.leftArmStartX, fixStimulus.leftArmStartY, fixStimulus.leftArmEndX, fixStimulus.leftArmEndY, ...
        %     fixStimulus.fixLineWidth, fixStimulus.thisArmColor);
        mglMetalLines(fixStimulus.leftArmStartX, fixStimulus.leftArmStartY, fixStimulus.leftArmEndX, fixStimulus.leftArmEndY, ...
            fixStimulus.fixLineWidth, fixStimulus.thisArmColor');
    elseif task.thistrial.sigLocation == 2
        % mglLines2(fixStimulus.rightArmStartX, fixStimulus.rightArmStartY, fixStimulus.rightArmEndX, fixStimulus.rightArmEndY, ...
        %     fixStimulus.fixLineWidth, fixStimulus.thisArmColor);
        mglMetalLines(fixStimulus.rightArmStartX, fixStimulus.rightArmStartY, fixStimulus.rightArmEndX, fixStimulus.rightArmEndY, ...
            fixStimulus.fixLineWidth, fixStimulus.thisArmColor');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called when subject responds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = fixTrialResponseCallback(task, myscreen)

global fixStimulus;

if task.thistrial.gotResponse < 1
    % get correct or incorrect
    response = find(task.thistrial.buttonState) == task.thistrial.sigLocation;
    response = response(1);

    if response
        % for training mode
        if fixStimulus.trainingMode
            mglDeleteTexture(fixStimulus.displayText);
            fixStimulus.displayText = mglText('Correct!');
            fixStimulus.displayTextLoc = [fixStimulus.pos(1) fixStimulus.pos(2)+fixStimulus.diskSize*2];
        end
        % set to correct fixation color
        fixStimulus.thisColor = fixStimulus.correctColor;
        task.correct = task.correct+1;
        task.correctResponse(end) = 1;
    else
        % set to incorrect fixation color
        fixStimulus.thisColor = fixStimulus.incorrectColor;
        task.incorrect = task.incorrect+1;
        task.correctResponse(end) = 0;
    end
end

function drawFixationCross(width, linewidth, color, origin)
% horizontal line
x0=origin(1)-width/2;
x1=origin(1)+width/2;
y0=origin(2);
y1=y0;
mglMetalLines(x0,y0,x1,y1, linewidth, color);
% vertical line
x0=origin(1);
x1=x0;
y0=origin(2)-width/2;
y1=origin(2)+width/2;
mglMetalLines(x0,y0,x1,y1, linewidth, color);
