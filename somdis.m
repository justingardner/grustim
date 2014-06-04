% somdis.m
%
%        $Id:$ 
%      usage: somdis()
%         by: justin gardner
%       date: 06/03/14
%    purpose: Basic somato discrim task
%
function myscreen = somdis

% check arguments
if ~any(nargin == [0])
  help taskTemplate
  return
end

% initalize the screen
myscreen = initScreen;

% get the last stimfile
stimfile = getLastStimfile(myscreen);

% task parameters
task{1}.waitForBacktick = 1;
task{1}.segmin = [1 2 1];
task{1}.segmax = [1 2 9];
task{1}.synchToVol = [0 0 1];
task{1}.getResponse = [1 1 0];
task{1}.parameter.pedestal = [-1 0.5 0.25 0.125 0];
task{1}.randVars.uniform.side = [-1 1];
task{1}.randVars.calculated.threshold = nan;
task{1}.randVars.calculated.correct = nan;
task{1}.randVars.calculated.responseSide = nan;
task{1}.random = 1;

testing = false;
if testing
  task{1}.segmin = [1 3 1];
  task{1}.segmax = [1 3 1];
  task{1}.parameter.pedestal = [-1 0.5];
  task{1}.synchToVol = [0 0 0];
end

% initialize the task
for phaseNum = 1:length(task)
  [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end

% init the stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);
stimulus = myInitStimulus(stimulus,myscreen,task,stimfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%myscreen = eyeCalibDisp(myscreen);

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

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;

if task.thistrial.thisseg == 1
  % time, so that we can preciesly set the two stimulation intervals
  timeNow = mglGetSecs;
  % pedestal stimulus strength
  ped = task.thistrial.pedestal;
  % delta 
  staircaseNum = find(stimulus.pedestal==ped);
  [delta stimulus.s(staircaseNum)] = doStaircase('testValue',stimulus.s(staircaseNum));
  if (delta < 0) delta = 0;end
  % if this is ped == -1, then it means to give max difference randomly betwen sides
  if ped == -1
    ped = 0;
    delta = 1;
  end
  % add the delta to the correct side
  if task.thistrial.side == -1
    stimLeft = ped+delta;
    stimRight = ped;
  else
    stimLeft = ped;
    stimRight = ped+delta;
  end
  % make sure we do not go over the maximum amplitude
  if stimLeft > 1,stimLeft = 1;end
  if stimRight > 1,stimRight = 1;end
  % multiply to get the actual voltage amplitude
  stimLeft = stimLeft * stimulus.maxStim;
  stimRight = stimRight * stimulus.maxStim;
  % set up the two stimulation period
  mglDigIO('ao',timeNow+0.1,[0 1],250,[stimRight stimLeft],0.4,[],2);
  mglDigIO('ao',timeNow+0.6,[0 1],250,[stimRight stimLeft],0.4,[],2);
  % display what we are doing
  disp(sprintf('Trial %i: %0.3f vs %0.3f (pedestal: %0.2f delta: %0.2f side: %i)',task.trialnum,stimLeft,stimRight,ped,delta,task.thistrial.side));
  % store delta
  task.thistrial.delta = delta;
end

if task.thistrial.thisseg == 3
  % no response
  if isnan(task.thistrial.correct)
    task.thistrial.correct = false;
    % update appropriate staircase
    staircaseNum = find(stimulus.pedestal==task.thistrial.pedestal);
    stimulus.s(staircaseNum) = doStaircase('update',stimulus.s(staircaseNum),task.thistrial.correct);
    %  threshold = doStaircase('threshold',stimulus.s(staircaseNum));
    disp(sprintf('  No response. Delta: %0.2f',task.thistrial.delta));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus

% display fixation cross
mglClearScreen;
if any(task.thistrial.thisseg == [1 2])
%  mglFixationCross(1,1,[0 1 1]);
  if task.thistrial.side == -1
    leftColor = [0 1 0];
    rightColor = [1 0 0];
  else
    leftColor = [1 0 0];
    rightColor = [0 1 0];
  end
  % set colors only after the subject has responded
  if ~isnan(task.thistrial.correct)
    % display the chosen side in green or red contingent
    % on whether subject was correct or not
    if task.thistrial.responseSide == -1
      mglFillOval(-4,0,[0.5 0.5],leftColor);
    else
      mglFillOval(-4,0,[0.5 0.5],[1 1 1]);
    end
    if task.thistrial.responseSide == 1
      mglFillOval(4,0,[0.5 0.5],rightColor);
    else
      mglFillOval(4,0,[0.5 0.5],[1 1 1]);
    end
  else
    mglFillOval(-4,0,[0.5 0.5],[1 1 1]);
    mglFillOval(4,0,[0.5 0.5],[1 1 1]);
  end

end  

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)

global stimulus

if task.thistrial.gotResponse < 1
  % caluclate response
  task.thistrial.responseSide = (task.thistrial.whichButton-1.5)*2;
  if task.thistrial.responseSide == task.thistrial.side
    task.thistrial.correct = true;
  else
    task.thistrial.correct = false;
  end
  % update appropriate staircase
  staircaseNum = find(stimulus.pedestal==task.thistrial.pedestal);
  stimulus.s(staircaseNum) = doStaircase('update',stimulus.s(staircaseNum),task.thistrial.correct);
%  threshold = doStaircase('threshold',stimulus.s(staircaseNum));
  disp(sprintf('  Response: %i Correct: %i Delta: %0.2f',task.thistrial.responseSide,task.thistrial.correct,task.thistrial.delta));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen,task,stimfile)

% set the maximum amplitude
stimulus.maxStim = 2;

% determine pedestals
params = getTaskParameters(myscreen,task);
pedestal = params.originalTaskParameter.pedestal;
numPedestal = length(pedestal);
stimulus.pedestal = pedestal;
stimulus.numPedestal = numPedestal;

% first time, initialize the staircases,
if isempty(stimfile)
  for iStaircase = 1:numPedestal
    stimulus.s(iStaircase) = doStaircase('init','upDown','nup=1','ndown=2','initialThreshold=0.3','initialStepsize=0.1','nTrials=40','stepRule=levitt','minThreshold=0');
  end
% subsequent times, continue staircases from where we left off
else
  oldStimulus = stimfile.stimulus;
  % make sure all is the same
  if ~isfield(oldStimulus,'pedestal') || ~isfield(oldStimulus,'numPedestal') || ~isfield(oldStimulus,'s') || ~isequal(oldStimulus.pedestal,pedestal) || ~isequal(oldStimulus.numPedestal,numPedestal) || ~(length(oldStimulus.s) == numPedestal)
    disp(sprintf('(somdis) !!! Previous stimfile has different conditions than this run, so restarting staircases !!!'));
    stimulus = myInitStimulus(stimulus,myscreen,task,[]);
    return
  end
  % if all is the same, then just reinit the staircases with the previous ones
  for iStaircase = 1:numPedestal
    stimulus.s(iStaircase) = doStaircase('init',oldStimulus.s(iStaircase));
  end
end
  
  

