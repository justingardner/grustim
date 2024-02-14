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
function [task myscreen] = fixStairInitTaskLR(myscreen)

% check arguments
if ~any(nargin == [1])
  help fixDispStairInitTask
  return
end

% create the stimulus for the experiment, use defaults if they are
% not already set
global fixStimulus;
myscreen = initStimulus('fixStimulus',myscreen);
if ~isfield(fixStimulus,'threshold') fixStimulus.threshold = 0.5; end
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
task{1}.seglen = [fixStimulus.interTime, fixStimulus.stimTime, fixStimulus.interTime, fixStimulus.responseTime];
%[fixStimulus.interTime fixStimulus.stimTime fixStimulus.interTime fixStimulus.stimTime fixStimulus.interTime fixStimulus.responseTime];
task{1}.getResponse = [0 0 0 1]; %[0 0 0 0 0 1];
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = fixTrialStartCallback(task, myscreen)

global fixStimulus;
% choose stimulus interval
task.thistrial.sigLocation = 1+(rand > 0.5); % 1 for left, 2 for right
if fixStimulus.verbose
  disp(sprintf('sigint = %i threshold = %0.2f',task.thistrial.sigLocation,fixStimulus.threshold));
end

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
if task.thistrial.thisseg == 2 % change
    fixStimulus.thisColor = fixStimulus.interColor;
    fixStimulus.thisArmColor = fixStimulus.interColor;
    fixStimulus.thisArmColor(1) = fixStimulus.threshold;
    
    % write out what the strength is
  % myscreen = writeTrace(fixStimulus.thisArmColor(1),task.fixStairTrace,myscreen);
  % training mode text
  if fixStimulus.trainingMode
      fixStimulus.displayText = mglText('Change');
      fixStimulus.displayTextLoc = [fixStimulus.pos(1) fixStimulus.pos(2)+fixStimulus.diskSize*2];
  end
% if this is the response interval
elseif task.thistrial.thisseg == 4 
    fixStimulus.thisColor = fixStimulus.responseColor;
    if fixStimulus.trainingMode
        fixStimulus.displayText = mglText('Which side changed (1 or 2)?');
        fixStimulus.displayTextLoc = [fixStimulus.pos(1) fixStimulus.pos(2)+fixStimulus.diskSize*2];
    end
% if this is the inter stimulus interval
else
    % training mode, clear sceen here
  if fixStimulus.trainingMode,mglClearScreen;end
  fixStimulus.thisColor = fixStimulus.interColor;
  % myscreen = writeTrace(0,task.fixStairTrace,myscreen);
end

% % choose what color the fixation cross will be
% whichInterval = find(task.thistrial.thisseg == [2 4]);
% % if this is the signal interval
% if ~isempty(whichInterval)
%   if task.thistrial.sigInterval == whichInterval
%     fixStimulus.thisStrength = (1-(fixStimulus.pedestal+fixStimulus.threshold));
%   else
%     fixStimulus.thisStrength = (1-fixStimulus.pedestal);
%   end
%   fixStimulus.thisColor = fixStimulus.stimColor*fixStimulus.thisStrength;
%   % write out what the strength is
%   myscreen = writeTrace(fixStimulus.thisStrength,task.fixStairTrace,myscreen);
%   % training mode text
%   if fixStimulus.trainingMode
%     if whichInterval == 1
%       fixStimulus.displayText = mglText('Interval 1');
%       fixStimulus.displayTextLoc = [fixStimulus.pos(1) fixStimulus.pos(2)+fixStimulus.diskSize*2];
%     else
%       fixStimulus.displayText = mglText('Interval 2');
%       fixStimulus.displayTextLoc = [fixStimulus.pos(1) fixStimulus.pos(2)+fixStimulus.diskSize*2];
%     end
%   end
% % if this is the response interval
% elseif task.thistrial.thisseg == 6
%   fixStimulus.thisColor = fixStimulus.responseColor;
%   if fixStimulus.trainingMode
%     fixStimulus.displayText = mglText('Which interval was darker (1 or 2)?');
%     fixStimulus.displayTextLoc = [fixStimulus.pos(1) fixStimulus.pos(2)+fixStimulus.diskSize*2];
%   end
% % if this is the inter stimulus interval
% else
%   % training mode, clear sceen here
%   if fixStimulus.trainingMode,mglClearScreen;end
%   fixStimulus.thisColor = fixStimulus.interColor;
%   myscreen = writeTrace(0,task.fixStairTrace,myscreen);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called every frame udpate to draw the fixation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = fixDrawStimulusCallback(task, myscreen)

global fixStimulus;

if fixStimulus.trainingMode,mglClearScreen;end

if ~isempty(fixStimulus.displayText)
  mglBltTexture(fixStimulus.displayText,fixStimulus.displayTextLoc);
end
mglGluDisk(fixStimulus.pos(1),fixStimulus.pos(2),fixStimulus.diskSize*[1 1],myscreen.background,60);
% mglMetalDots([fixStimulus.pos(1),fixStimulus.pos(2),0],[myscreen.background, 1],fixStimulus.diskSize*[1 1]);

drawFixationCross(fixStimulus.fixWidth,fixStimulus.fixLineWidth,fixStimulus.thisColor',fixStimulus.pos);

if task.thistrial.thisseg == 2
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
  % set trace to 2 to indicate correct response
  % myscreen = writeTrace(2,task.fixStairTrace,myscreen);
  % and update correct count
  task.correct = task.correct+1;
else
  if fixStimulus.trainingMode
    mglDeleteTexture(fixStimulus.displayText);
    fixStimulus.displayText = mglText('Incorrect.');
    fixStimulus.displayTextLoc = [fixStimulus.pos(1) fixStimulus.pos(2)+fixStimulus.diskSize*2];
  end
  % set to incorrect fixation color
  fixStimulus.thisColor = fixStimulus.incorrectColor;
  % set trace to -2 to indicate incorrect response
  % myscreen = writeTrace(-2,task.fixStairTrace,myscreen);
  % and update incorrect count
  task.incorrect = task.incorrect+1;
end

% update staircase
fixStimulus.staircase = upDownStaircase(fixStimulus.staircase,response);
fixStimulus.threshold = fixStimulus.staircase.threshold;

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
