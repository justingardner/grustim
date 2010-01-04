% od
%
%      usage: myscreen=od(stimulus)
%         by: justin gardner
%       date: 04/15/06
%    purpose: orientation discrimination task
%
%
%
function myscreen = cpoolrule(varargin)

% check arguments
if ~any(nargin == [0 1])
  help orientationDiscrimination
  return
end

isLoc = [];
getArgs(varargin,{'isLoc=0'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other screen parameters
myscreen.autoCloseScreen = 1;
myscreen.allowpause = 0;
myscreen.eatkeys = 0;
myscreen.displayname = 'projector';
myscreen.background = 'gray';

myscreen = initScreen(myscreen);

global MGL;
clear global stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);

% set the contrasts of distractor and target
pedestalContrasts = [0.1 0.2 0.4];
deltaContrasts = [0 0.05 0.1 0.15 0.2];

pedestalContrasts = [0.1];
deltaContrasts = [0 0.05 0.1 0.2]/8;

% create an array of all contrasts that will need to be displayed
allContrasts = repmat(deltaContrasts,length(pedestalContrasts),1);
allContrasts = allContrasts(:)';
allContrasts = allContrasts + repmat(pedestalContrasts,1,length(deltaContrasts));
allContrasts = unique(allContrasts);

% parameters
stimulus.int1 = 2;
stimulus.int2 = 4;
stimulus.grating.radius = 4.5;
stimulus.grating.orientations = -45:90:359;
stimulus.grating.n = 4;
stimulus.grating.contrasts = allContrasts;
stimulus.grating.sf = 2;
stimulus.grating.tf = 2;
stimulus.grating.width = 6;
stimulus.grating.height = 6;
stimulus.grating.phases = [0 pi];
stimulus.grating.phases = [0:2*pi/10:2*pi 2*pi:-2*pi/10:0]; % note this makes the actual tf = tf*2
stimulus.grating.phase = 0;
stimulus.grating.windowType = 'gabor'; % should be gabor or thresh
stimulus.grating.sdx = stimulus.grating.width/7;
stimulus.grating.sdy = stimulus.grating.width/7;

stimulus.grating.width = 6;
stimulus.grating.height = 6;
stimulus.grating.windowType = 'thresh'; % should be gabor or thresh
stimulus.grating.sdx = 2.5;
stimulus.grating.sdy = 2.5;

stimulus.fixColor = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up fixation task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global fixStimulus;
easyFixTask = 1;
if ~easyFixTask
  % default values
  fixStimulus.diskSize = 0.5;
  fixStimulus.fixWidth = 1;
  fixStimulus.fixLineWidth = 3;
  fixStimulus.stimTime = 0.4;
  fixStimulus.responseTime = 1;
else
  % make cross bigger and task slower
  fixStimulus.diskSize = 0.5;
  fixStimulus.fixWidth = 1+1*easyFixTask;
  fixStimulus.fixLineWidth = 3+2*easyFixTask;
  fixStimulus.stimTime = 0.4+0.4*easyFixTask;
  fixStimulus.responseTime = 1+1*easyFixTask;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up segments of trials
if ~isLoc
  stimulus.isLocalizer = 0;
  task{1}{1}.random = 1;
  task{1}{1}.segmin = [0];
  task{1}{1}.segmax = [0];
  task{1}{1}.synchToVol = [0];
  task{1}{1}.getResponse = [0];
  task{1}{1}.numTrials = 1;
  task{1}{1}.waitForBacktick = 0;
  task{1}{2}.parameter.pedestalContrast = repmat(pedestalContrasts,stimulus.grating.n,1);
  task{1}{2}.parameter.deltaContrast = repmat(deltaContrasts,stimulus.grating.n,1);
  task{1}{2}.parameter.interval = repmat([stimulus.int1 stimulus.int2],stimulus.grating.n,1);
  task{1}{2}.random = 1;
  task{1}{2}.segmin = [1 0.6 0.3 0.6 1.5];
  task{1}{2}.segmax = [1 0.6 0.3 0.6 1.5];
  task{1}{2}.synchToVol = [0 0 0 0 0];
  task{1}{2}.getResponse = [0 0 0 0 1];
  task{1}{2}.waitForBacktick = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if this is localizer then change a few things
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
  stimulus.isLocalizer = 1;
  stimulus.grating.contrasts = [0 1];
  task{1}{1}.parameter.distractorContrast = distractorContrast;
  task{1}{1}.parameter.targetContrast = targetContrast;
  task{1}{1}.seglen = [12 12];
  task{1}{1}.synchToVol = [1 1];
  task{1}{1}.waitForBacktick=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialze tasks and stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stimulus = initGratings(stimulus,myscreen,task)

% initialze tasks
for phaseNum = 1:length(task{1})
  [task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@updateScreenCallback,@trialResponseCallback);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the tasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set which phase is active
tnum = 1;

phaseNum = 1;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
  % update the task
  [task{1} myscreen phaseNum] = updateTask(task{1},myscreen,phaseNum);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     STIMULUS Gratings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initGratings(stimulus,myscreen,task)

disppercent(-inf,'Creating grating textures');

nContrasts = length(stimulus.grating.contrasts);
nOrientations = length(stimulus.grating.orientations);
nPhases = length(stimulus.grating.phases);

% make each one of he called for gratings
for iPhase = 1:nPhases
  for iContrast = 1:nContrasts
    for iOrientation = 1:nOrientations
      pDone = calcPercentDone(iPhase,nPhases,iContrast,nContrasts,iOrientation,nOrientations);
      % display pDone to screen
      mglClearScreen;
      mglTextDraw(sprintf('%i%%',round(100*pDone)),[0 0]);
      myscreen = tickScreen(myscreen);
      disppercent(pDone);
      if myscreen.userHitEsc,mglClose;keyboard,end
      % get the orientation we want to create
      thisOrientation = stimulus.grating.orientations(iOrientation);
      % get the phase
      thisPhase = (stimulus.grating.phase+stimulus.grating.phases(iPhase))*180/pi;
      % make the grating
      thisGrating = 255*(mglMakeGrating(stimulus.grating.width,stimulus.grating.height,...
					stimulus.grating.sf,thisOrientation,...
					thisPhase)+1)/2;
  
      thisGaussian = mglMakeGaussian(stimulus.grating.width,stimulus.grating.height,...
				     stimulus.grating.sdx,stimulus.grating.sdy);

      % create an rgb/a matrix
      thisStimulus(:,:,1) = thisGrating;
      thisStimulus(:,:,2) = thisGrating;
      thisStimulus(:,:,3) = thisGrating;

      % make the gaussian window
      if strcmp(stimulus.grating.windowType,'gabor')
	% create an rgb/a matrix
	mask = stimulus.grating.contrasts(iContrast)*255*thisGaussian;
      else
	mask = stimulus.grating.contrasts(iContrast)*255*(thisGaussian>exp(-1/2));
      end

      thisStimulus(:,:,4) = mask;
      thatStimulus(:,:,4) = mask;
      
      % create the texture
      stimulus.tex(iContrast,iOrientation,iPhase) = mglCreateTexture(thisStimulus);
    end
  end
end
disppercent(inf);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to set stimulus parameters at
% the beginning of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task,myscreen)

global stimulus;

if any(task.thistrial.thisseg == [stimulus.int1 stimulus.int2])
  % set up the stimuli
  for i = 1:stimulus.grating.n
    % orientation of stimulus
    orientation = pi*stimulus.grating.orientations(i)/180;

    % get x,y position of grating
    stimulus.x(i) = stimulus.grating.radius*cos(orientation);
    stimulus.y(i) = stimulus.grating.radius*sin(orientation);

    % set the contrast num
    if (task.thistrial.thisseg == task.thistrial.interval(i))
      thisContrast = task.thistrial.pedestalContrast(i)+task.thistrial.deltaContrast(i);
    else
      thisContrast = task.thistrial.pedestalContrast(i);
    end
    stimulus.contrastNum(i) = find(thisContrast == stimulus.grating.contrasts);
    
  end
  % compute overall contrast in the interval
  if task.thistrial.thisseg == stimulus.int1
    stimulus.contrast1 = sum(stimulus.grating.contrasts(stimulus.contrastNum));
  else
    stimulus.contrast2 = sum(stimulus.grating.contrasts(stimulus.contrastNum));
  end
end

% set the fixation color
if any(task.thistrial.thisseg == [stimulus.int1 stimulus.int2])
  stimulus.fixColor = [1 1 0];
else
  stimulus.fixColor = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to display stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task,myscreen)

global stimulus;
mglClearScreen;

if any(task.thistrial.thisseg == [stimulus.int1 stimulus.int2]) 
  for i = 1:stimulus.grating.n
    % get which phase we are on
    phaseNum = floor(length(stimulus.grating.phases)*rem(mglGetSecs(task.thistrial.trialstart)*stimulus.grating.tf,1)+1);
    % reverse the phase every other grating
    if iseven(i)
      phaseNum = mod(phaseNum+length(stimulus.grating.phases)/2,length(stimulus.grating.phases))+1;
    end
    % blt texture
    mglBltTexture(stimulus.tex(stimulus.contrastNum(i),i,phaseNum),[stimulus.x(i) stimulus.y(i)]);
  end
end

mglFixationCross(0.5,1,stimulus.fixColor);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to set stimulus parameters at
% the beginning of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = trialResponseCallback(task,myscreen)

global stimulus;

if (stimulus.contrast1 > stimulus.contrast2) && (task.thistrial.whichButton == 1)
  disp(sprintf('correct'));
  stimulus.fixColor = [0 1 0];
elseif (stimulus.contrast1 < stimulus.contrast2) && (task.thistrial.whichButton == 2)
  disp(sprintf('correct'));
  stimulus.fixColor = [0 1 0];
elseif (stimulus.contrast1 > stimulus.contrast2) && (task.thistrial.whichButton == 2)
  disp(sprintf('incorrect'));
  stimulus.fixColor = [1 0 0];
elseif (stimulus.contrast1 < stimulus.contrast2) && (task.thistrial.whichButton == 1)
  disp(sprintf('incorrect'));
  stimulus.fixColor = [1 0 0];
end
  
