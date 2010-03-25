% spatcon
%
%      usage: myscreen=spatcon(stimulus)
%         by: justin gardner
%       date: 04/15/06
%    purpose: spatial contrast response function experiment
%
%
%
function myscreen = spatcon1d_broken_fixme(varargin)

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

% compute contrast, we want them evenly spaced on a log scale
% and we specify the middle, max and number of contrast
midContrast = 50;
maxContrast = 85;
nContrasts = 5;

% now calculate the contrasts we need
logContrastDifference = log(maxContrast)-log(midContrast);
targetContrast = exp(log(midContrast)+(-logContrastDifference:2*logContrastDifference/(nContrasts-1):logContrastDifference));
targetContrast = targetContrast/100;
% and display them
disp(sprintf('(spatcon) targetContrasts: %s',mynum2str(targetContrast)))

% set the contrasts of distractor and target
%distractorContrast = [0.1 1];
%targetContrast = [0.1 0.5 1];
%targetContrast = [0.25 0.375 0.5 0.75];
distractorContrast = [0.5];

% parameters
stimulus.grating.radius = 6.5;
stimulus.grating.targetLoc = [1 4];
stimulus.grating.orientations = 0:60:359;
stimulus.grating.contrasts = union(distractorContrast,targetContrast);
stimulus.grating.sf = 2;
stimulus.grating.tf = 2;
stimulus.grating.width = 9;
stimulus.grating.height = 9;
stimulus.grating.phases = [0 pi];
stimulus.grating.phases = [0:2*pi/30:2*pi 2*pi:-2*pi/30:0]; % note this makes the actual tf = tf*2
stimulus.grating.phase = 0;
stimulus.grating.windowType = 'gabor'; % should be gabor or thresh
stimulus.grating.sdx = stimulus.grating.width/7;
stimulus.grating.sdy = stimulus.grating.width/7;
stimulus.x = 0;

stimulus.y = 0.5;
stimulus.grating.width = 5.5;
stimulus.grating.height = 5.5;
stimulus.grating.windowType = 'thresh'; % should be gabor or thresh
stimulus.grating.sdx = stimulus.grating.width/2;
stimulus.grating.sdy = stimulus.grating.height/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up fixation task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global fixStimulus;
easyFixTask = 1;
fixStimulus.pos = [stimulus.x stimulus.y];
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
[task{1} myscreen] = fixStairInitTask(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up segments of trials
if ~isLoc
  stimulus.isLocalizer = 0;
  task{2}{1}.parameter.distractorContrast = 0;
  task{2}{1}.parameter.targetContrast = 0;
  task{2}{1}.parameter.targetLoc = 1;
  task{2}{1}.random = 1;
  task{2}{1}.segmin = [10];
  task{2}{1}.segmax = [10];
  task{2}{1}.synchToVol = [1];
  task{2}{1}.getResponse = [0];
  task{2}{1}.numTrials = 1;
  task{2}{1}.waitForBacktick = 1;

  task{2}{2}.parameter.distractorContrast = distractorContrast;
  task{2}{2}.parameter.targetContrast = targetContrast;
  task{2}{2}.parameter.targetLoc = stimulus.grating.targetLoc;
  task{2}{2}.random = 1;
  task{2}{2}.segmin = [3 4];
  task{2}{2}.segmax = [3 8];
  task{2}{2}.synchToVol = [0 1];
  task{2}{2}.getResponse = [0 0];
  task{2}{2}.waitForBacktick = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if this is localizer then change a few things
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
  stimulus.isLocalizer = 1;
  stimulus.grating.contrasts = [0 1];
  task{2}{1}.parameter.distractorContrast = distractorContrast;
  task{2}{1}.parameter.targetContrast = targetContrast;
  task{2}{1}.seglen = [12 12];
  task{2}{1}.synchToVol = [1 1];
  task{2}{1}.waitForBacktick=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialze tasks and stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stimulus = initGratings(stimulus,myscreen,task);

% initialze tasks
for phaseNum = 1:length(task{2})
  [task{2}{phaseNum} myscreen] = initTask(task{2}{phaseNum},myscreen,@startSegmentCallback,@updateScreenCallback,@trialResponseCallback);
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
while (phaseNum <= length(task{2})) && ~myscreen.userHitEsc
  % update the dots
  [task{2} myscreen phaseNum] = updateTask(task{2},myscreen,phaseNum);
  % update the fixation task
  [task{1} myscreen] = updateTask(task{1},myscreen,1);
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

maxIndex = 255;
disppercent(-inf,'Creating grating textures');

nContrasts = length(stimulus.grating.contrasts);
nOrientations = length(stimulus.grating.orientations);
nPhases = length(stimulus.grating.phases);

gaussianWin = mglMakeGaussian(stimulus.grating.width,stimulus.grating.height,stimulus.grating.sdx,stimulus.grating.sdy);
if strcmp(stimulus.grating.windowType,'gabor')
  % a gaussian window
  win = maxIndex-maxIndex*gaussianWin;
else
  % a simple window
  win = maxIndex-maxIndex*(gaussianWin>exp(-1/2));
end
mask = ones(size(win,1),size(win,2),4)*myscreen.grayIndex;
mask(:,:,4) = round(win);
stimulus.mask = mglCreateTexture(mask);

% make each one of he called for gratings
for iPhase = 1:nPhases
  for iContrast = 1:nContrasts
    disppercent(calcPercentDone(iPhase,nPhases,iContrast,nContrasts));
    % get the phase and contast
    thisPhase = (stimulus.grating.phase+stimulus.grating.phases(iPhase))*180/pi;
    thisContrast = stimulus.grating.contrasts(iContrast);
    % make the grating
    thisGrating = round(maxIndex*((thisContrast*mglMakeGrating(stimulus.grating.width,nan,stimulus.grating.sf,0,thisPhase))+1)/2);
    % create the texture
    stimulus.tex(iContrast,iPhase) = mglCreateTexture(thisGrating);
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
if task.thistrial.thisseg == 1
  disp(sprintf('Target: %0.2f Distractor: %0.2f',task.thistrial.targetContrast,task.thistrial.distractorContrast));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to display stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task,myscreen)

global stimulus;
mglClearScreen;

if (task.thistrial.thisseg == 1) || (stimulus.isLocalizer)
  % create a ring of stimuli
  for iAngle = 1:length(stimulus.grating.orientations)

    % get x,y position of grating
    thisAngle = stimulus.grating.orientations(iAngle);
    x = stimulus.x + stimulus.grating.radius*cos(pi*thisAngle/180);
    y = stimulus.y + stimulus.grating.radius*sin(pi*thisAngle/180);
    angleNum = iAngle;

    % get which phase we are on
    phaseNum = floor(length(stimulus.grating.phases)*rem(mglGetSecs(task.thistrial.trialstart)*stimulus.grating.tf,1)+1);
    % reverse the phase every other grating
    if iseven(iAngle)
      phaseNum = mod(phaseNum+length(stimulus.grating.phases)/2,length(stimulus.grating.phases))+1;
    end
    

    % get what contrast to display
    % Main experiment
    if ~stimulus.isLocalizer
      if iAngle == task.thistrial.targetLoc
	contrastNum = find(stimulus.grating.contrasts==task.thistrial.targetContrast);
	angleNum = find(stimulus.grating.orientations==0);
      else
	contrastNum = find(stimulus.grating.contrasts==task.thistrial.distractorContrast);
      end    

    % Localizer
    else
      if any(iAngle == stimulus.grating.targetLoc)
	% flip phase so that left and right gratings are opposite each other
	if (iAngle == stimulus.grating.targetLoc(1))
	  phaseNum = mod(phaseNum+length(stimulus.grating.phases)/2,length(stimulus.grating.phases))+1;
	end
	if task.thistrial.thisseg == 1
	  contrastNum = 2;
	else
	  contrastNum = 1;
	end
      else
	if task.thistrial.thisseg == 1
	  contrastNum = 1;
	else
	  contrastNum = 2;
	end
      end
    end

    % display the texture
    mglBltTexture(stimulus.tex(contrastNum,phaseNum),[x y stimulus.grating.height],0,0,thisAngle);
    mglBltTexture(stimulus.mask,[x y stimulus.grating.height],0,0,thisAngle);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to set stimulus parameters at
% the beginning of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = trialResponseCallback(task,myscreen)
v
global stimulus;

