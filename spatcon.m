% spatcon
%
%      usage: myscreen=spatcon(stimulus)
%         by: justin gardner
%       date: 04/15/06
%    purpose: spatial contrast response function experiment
%
%
%
function myscreen = spatcon(varargin)

taskType = [];
initStair = [];
threshold = [];
stepsize = [];
useLevittRule = [];
stimFile = [];
numBlocks = [];
targetContrast = [];
distractorContrast = [];
subjectID = [];
easyFixTask = [];
getArgs(varargin,{'taskType=1','initStair=1','threshold=6','stepsize=2','useLevittRule=1','stimFile=[]','numBlocks=100','targetContrast=[1 1/9]','distractorContrast=[1/3]','subjectID=[]','easyFixTask=0'});

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

% set data directory to be ~/data/subjectID
myscreen.datadir = fullfile(myscreen.datadir,'spatcon');
if ~isdir(myscreen.datadir)
  mkdir(myscreen.datadir);
end
if ~isempty(subjectID)
  myscreen.datadir = fullfile(myscreen.datadir,subjectID);
  if ~isdir(myscreen.datadir)
    mkdir(myscreen.datadir);
  end
end
disp(sprintf('(spatcon) Saving data into %s',myscreen.datadir));

global stimulus;
if initStair
  clear global stimulus
  global stimulus;
end

if ~isempty(stimFile)
  stimFile = setext(stimFile,'mat');
  if isfile(stimFile)
    s = load(stimFile);
  elseif isfile(fullfile(myscreen.datadir,stimFile))
    s = load(fullfile(myscreen.datadir,stimFile));
  else
    disp(sprintf('(spatcon) Could not find file %s',stimFile));
    endScreen(myscreen);
    return
  end
  stimulus = s.stimulus;
  numBlocks = numBlocks - s.task{2}{2}.blocknum + 1;
  clear s;
  initStair = 0;
end
myscreen = initStimulus('stimulus',myscreen);

% set the time to wait before starting task
initWaitTime = 0.1;

% set target contrasts
if taskType == 1
  % compute contrast, we want them evenly spaced on a log scale
  % and we specify the middle, max and number of contrast
  %midContrast = 50;
  %maxContrast = 100;
  %nContrasts = 3;
  % now calculate the contrasts we need
  %logContrastDifference = log(maxContrast)-log(midContrast);
  %stimulus.targetContrast = exp(log(midContrast)+(-logContrastDifference:2*logContrastDifference/(nContrasts-1):logContrastDifference));
  %stimulus.targetContrast = stimulus.targetContrast/100;
%  targetContrast = 0:0.2:1;
%  distractorContrast = [0.2 0.8];
  targetContrast = [0.1 0.5 1];
  distractorContrast = 0.5;
  
elseif any(taskType==[2 3])
%  stimulus.targetContrast = [1 0.25];
elseif taskType == 0
  targetContrast = 1;
  distractorContrast = 0;
end

% set the contrasts
stimulus.targetContrast = targetContrast;
stimulus.distractorContrast = distractorContrast;

% display contaststhem
disp(sprintf('(spatcon) targetContrasts: %s distractorContrast: %s',mynum2str(stimulus.targetContrast),mynum2str(stimulus.distractorContrast)));


% parameters
stimulus.grating.radius = 10;
stimulus.grating.targetLoc = [1 4];
stimulus.grating.angles = 0:60:359;
stimulus.grating.orientations = stimulus.grating.angles;
if any(taskType == [2 3])
  stimulus.grating.orientations(1) = 45;
  stimulus.grating.orientations(4) = -45;
end
stimulus.grating.contrasts = union(stimulus.distractorContrast,stimulus.targetContrast);
stimulus.grating.sf = 2;
stimulus.grating.tf = 2;
stimulus.grating.width = 11;
stimulus.grating.height = 11;
stimulus.grating.phases = [0 pi];
stimulus.grating.phases = [0:2*pi/30:2*pi 2*pi:-2*pi/30:0]; % note this makes the actual tf = tf*2
stimulus.grating.phase = 0;
stimulus.grating.windowType = 'gabor'; % should be gabor or thresh
stimulus.grating.sdx = stimulus.grating.width/7;
stimulus.grating.sdy = stimulus.grating.width/7;
stimulus.x = 0;

stimulus.y = 1.5;
stimulus.grating.width = 8.5;
stimulus.grating.height = 8.5;
stimulus.grating.windowType = 'thresh'; % should be gabor or thresh
stimulus.grating.sdx = stimulus.grating.width/2;
stimulus.grating.sdy = stimulus.grating.height/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up fixation task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global fixStimulus;
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
if taskType == 1
  disp(sprintf('(spatcon) Running fixation task main experiment.'));
  stimulus.isLoc = 0;
  stimulus.taskType = taskType; % central fixation
  task{2}{1}.parameter.distractorContrast = 0;
  task{2}{1}.parameter.targetContrast = 0;
  task{2}{1}.parameter.targetLoc = 1;
  task{2}{1}.random = 1;
  task{2}{1}.seglen = [initWaitTime];
  task{2}{1}.synchToVol = [1];
  task{2}{1}.getResponse = [0];
  task{2}{1}.numTrials = 1;
  task{2}{1}.waitForBacktick = 0;

  task{2}{2}.parameter.distractorContrast = stimulus.distractorContrast;
  task{2}{2}.parameter.targetContrast = stimulus.targetContrast;
  task{2}{2}.parameter.targetLoc = stimulus.grating.targetLoc;
  task{2}{2}.random = 1;
  task{2}{2}.segmin = [3 4];
  task{2}{2}.segmax = [3 8];
  task{2}{2}.synchToVol = [0 1];
  task{2}{2}.getResponse = [0 0];
  task{2}{2}.waitForBacktick = 0;

  stimulus.stimSegment = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up attentive task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up segments of trials
elseif any(taskType == [2 3])
  disp(sprintf('(spatcon) Running attentive task, main experiment.'));
  stimulus.isLoc = 0;
  stimulus.taskType = taskType; % attentive task
  task{2}{1}.parameter.distractorContrast = 0;
  task{2}{1}.parameter.targetContrast = 0;
  task{2}{1}.parameter.targetLoc = 1;
  task{2}{1}.random = 1;
  task{2}{1}.getResponse = [0];
  task{2}{1}.numTrials = 1;
  task{2}{1}.waitForBacktick = 1;

  task{2}{2}.parameter.distractorContrast = stimulus.distractorContrast;
  task{2}{2}.parameter.targetContrast = stimulus.targetContrast;
  task{2}{2}.parameter.targetLoc = stimulus.grating.targetLoc;
  task{2}{2}.random = 1;

  % length of time a stimulus will be on for
  stimLen = 0.5;
  
  % scanner and psychophysics room
  if taskType == 2
    % psychophyscis room
    task{2}{1}.seglen = 0;
    task{2}{1}.synchToVol = [1];
    task{2}{2}.waitForBacktick = 0;
    task{2}{1}.waitForBacktick = 0;

    task{2}{2}.seglen = [1 stimLen 1 0.5];
%    task{2}{2}.seglen = 3*task{2}{2}.seglen;
    task{2}{2}.getResponse = [0 0 1];
    task{2}{2}.numBlocks = numBlocks;
    disp(sprintf('(spatcon) Number of blocks: %i',numBlocks));
  else
    % scanner
    task{2}{1}.seglen = initWaitTime;
    task{2}{1}.synchToVol = 1;
    task{2}{2}.segmin = [1 stimLen 1 3];
    task{2}{2}.segmax = [1 stimLen 1 7];
    task{2}{2}.synchToVol = [0 0 0 1];
    task{2}{2}.getResponse = [0 0 1 0];
    task{2}{2}.waitForBacktick = 0;
  end

  stimulus.stimSegment = 2;
  stimulus.fixColor = [1 1 1];
  stimulus.fixWidth = 1;
  stimulus.fixLineWidth = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if this is localizer then change a few things
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif taskType == 0
  disp(sprintf('(spatcon) Running Localizer'));
  stimulus.isLoc = 1;
  stimulus.taskType = taskType;
  stimulus.grating.contrasts = [0 1];
  task{2}{1}.parameter.distractorContrast = stimulus.distractorContrast;
  task{2}{1}.parameter.targetContrast = stimulus.targetContrast;
  task{2}{1}.seglen = [12 12];
  task{2}{1}.synchToVol = [1 1];
  task{2}{1}.waitForBacktick=1;

  stimulus.stimSegment = [1 2];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialze tasks and stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stimulus = initGratings(stimulus,myscreen,task);

if initStair
  disp(sprintf('(spatcon) Initializing staircase with threshold: %f stepsize: %f useLevittRule: %i',threshold,stepsize,useLevittRule));
  if any(taskType == [2 3])
    stimulus = initStaircase(threshold,stimulus,stepsize,useLevittRule);
  end
else
  disp(sprintf('(spatcon) Continuing staircase from last run'));
  targetDistractorStr = {'distractor','target'};
  for iLoc = stimulus.grating.targetLoc
    for iContrast = 1:length(stimulus.targetContrast);
      for iTargetDistractor = 1:2
	threshold = stimulus.staircase{iLoc}{iContrast}{iTargetDistractor}.threshold;
	if isfield(stimulus.staircase{iLoc}{iContrast}{iTargetDistractor},'strength')
	  n = length(stimulus.staircase{iLoc}{iContrast}{iTargetDistractor}.strength);
	else
	  n = 0;
	end
	disp(sprintf('%s (Loc: %i Contrast: %i): %f (n=%i)',targetDistractorStr{iTargetDistractor},iLoc,iContrast,threshold,n));
      end
    end
  end
end

% initialze tasks
for phaseNum = 1:length(task{2})
  [task{2}{phaseNum} myscreen] = initTask(task{2}{phaseNum},myscreen,@startSegmentCallback,@updateScreenCallback,@trialResponseCallback);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen,'Hit <space> to calibrate eye tracker, <return> to skip.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the tasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set which phase is active
tnum = 1;

phaseNum = 1;
while (phaseNum <= length(task{2})) && ~myscreen.userHitEsc
  % update the dotspXU
  [task{2} myscreen phaseNum] = updateTask(task{2},myscreen,phaseNum);
  if any(stimulus.taskType == [0 1])
    % update the fixation task
    [task{1} myscreen] = updateTask(task{1},myscreen,1);
  end
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

if any(stimulus.taskType == [2 3])
  if task{2}{2}.trialnum > 10
    spatconPsycho(stimulus);
  end
end


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

stimulus.maxIndex = 255;
disppercent(-inf,'Creating grating textures');

nContrasts = length(stimulus.grating.contrasts);
nPhases = length(stimulus.grating.phases);

gaussianWin = mglMakeGaussian(stimulus.grating.width,stimulus.grating.height,stimulus.grating.sdx,stimulus.grating.sdy);
if strcmp(stimulus.grating.windowType,'gabor')
  % a gaussian window
  win = stimulus.maxIndex-stimulus.maxIndex*gaussianWin;
else
  % a simple window
  win = stimulus.maxIndex-stimulus.maxIndex*(gaussianWin>exp(-1/2));
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
    thisGrating = round(stimulus.maxIndex*((thisContrast*mglMakeGrating(stimulus.grating.width,nan,stimulus.grating.sf,0,thisPhase))+1)/2);
    % create the texture
    stimulus.tex(iContrast,iPhase) = mglCreateTexture(thisGrating);
  end
end
disppercent(inf);
stimulus.randMaskSize = [size(mask,1) size(mask,2)];
stimulus.randMask = mglCreateTexture(floor(stimulus.maxIndex*rand(stimulus.randMaskSize)));


for iAngle = 1:length(stimulus.grating.angles)
  % get center of patch
  thisAngle = stimulus.grating.angles(iAngle);
  centerX = stimulus.x + stimulus.grating.radius*cos(pi*thisAngle/180);
  centerY = stimulus.y + stimulus.grating.radius*sin(pi*thisAngle/180);
  % now get top and bottom point of grating
  thisOrientation = stimulus.grating.orientations(iAngle)+90;
  radius = sqrt((stimulus.grating.width/2).^2 +(stimulus.grating.height/2).^2)-0.5;
  topX = centerX + radius*cos(pi*thisOrientation/180);
  topY = centerY + radius*sin(pi*thisOrientation/180);
  thisOrientation = thisOrientation+180;
  bottomX = centerX + radius*cos(pi*thisOrientation/180);
  bottomY = centerY + radius*sin(pi*thisOrientation/180);
  % place points
  stimulus.grating.refPoints.x{iAngle} = [topX bottomX];
  stimulus.grating.refPoints.y{iAngle} = [topY bottomY];
end

stimulus.waitForBacktickText = mglText('Hit backtick (`) key to start');

%%%%%%%%%%%%%%%%%%%%%%%%
%    startStaircase    %
%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initStaircase(threshold,stimulus,stepsize,useLevittRule)

for iLoc = stimulus.grating.targetLoc
  for iContrast = 1:length(stimulus.targetContrast);
    for iTargetDistractor = 1:2
      stimulus.staircase{iLoc}{iContrast}{iTargetDistractor} = upDownStaircase(1,2,threshold,stepsize,useLevittRule);
      stimulus.staircase{iLoc}{iContrast}{iTargetDistractor}.minThreshold = 0;
      stimulus.staircase{iLoc}{iContrast}{iTargetDistractor}.maxThreshold = 90;
    end
  end
end

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
  if stimulus.taskType > 0
    disp(sprintf('Target: %0.2f Distractor: %0.2f TargetLoc: %i',task.thistrial.targetContrast,task.thistrial.distractorContrast,task.thistrial.targetLoc));
  end
end

% set the orientation of the targets
if (stimulus.taskType == 2) && (task.thistrial.thisphase == 2)
  if (task.thistrial.thisseg == 1)
    % choose which direction is target
    stimulus.task.cuedLoc(task.trialnum) = stimulus.grating.targetLoc(ceil(length(stimulus.grating.targetLoc)*rand));
    % set a random orientation of the target at the current threshold estimate
    stimulus.task.direction{task.trialnum} = nan(1,length(stimulus.grating.orientations));
    stimulus.task.orientation{task.trialnum} = nan(1,length(stimulus.grating.orientations));
    % set these to decide which staircase to use. staircases are sorted by targetContrast
    % and whether the cuedLoc is the same as the target
    stimulus.thisContrastNum = find(task.thistrial.targetContrast == stimulus.targetContrast);
    stimulus.thisTargetDistractor = (stimulus.task.cuedLoc(task.trialnum) == task.thistrial.targetLoc)+1;
    for iTarget = stimulus.grating.targetLoc
      stimulus.task.direction{task.trialnum}(iTarget) = round((rand>0.5)*2-1);
      stimulus.task.orientation{task.trialnum}(iTarget) = stimulus.staircase{iTarget}{stimulus.thisContrastNum}{stimulus.thisTargetDistractor}.threshold*stimulus.task.direction{task.trialnum}(iTarget);
    end
    stimulus.fixColor = [1 1 1];
    % set up how to draw fixation lines
    o = d2r(stimulus.grating.angles(stimulus.task.cuedLoc(task.trialnum)));
    cueAngle = d2r(30);
    cueWidth = stimulus.fixWidth;
    stimulus.cueX = [cos(o)*cueWidth/2 -cos(o+cueAngle)*cueWidth/2 cos(o)*cueWidth/2 -cos(o-cueAngle)*cueWidth/2];
    stimulus.cueY = [sin(o)*cueWidth/2 -sin(o+cueAngle)*cueWidth/2 sin(o)*cueWidth/2 -sin(o-cueAngle)*cueWidth/2];
    disp(sprintf('Cued targer: %i',stimulus.task.cuedLoc(task.trialnum)));
  elseif (task.thistrial.thisseg == 2)
      % init a new randMask
      mglDeleteTexture(stimulus.randMask);
      stimulus.randMask = mglCreateTexture(floor(stimulus.maxIndex*rand(stimulus.randMaskSize)));
  elseif (task.thistrial.thisseg == 3)
    stimulus.fixColor = [1 1 0];
  elseif (task.thistrial.thisseg == 4)
    stimulus.fixColor = [1 1 1];
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to display stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task,myscreen)

global stimulus;
mglClearScreen;

% put up text to tell subject to hit backtick to start
if (stimulus.taskType == 2) && (task.thistrial.thisphase == 1)
  mglBltTexture(stimulus.waitForBacktickText,[0 0]);
  return;
end

if any(task.thistrial.thisseg == stimulus.stimSegment)
  % create a ring of stimuli
  for iAngle = 1:length(stimulus.grating.orientations)

    % get x,y position of grating
    thisAngle = stimulus.grating.angles(iAngle);
    x = stimulus.x + stimulus.grating.radius*cos(pi*thisAngle/180);
    y = stimulus.y + stimulus.grating.radius*sin(pi*thisAngle/180);
    angleNum = iAngle;
    thisOrientation = stimulus.grating.orientations(iAngle);

    % get which phase we are on
    phaseNum = floor(length(stimulus.grating.phases)*rem(mglGetSecs(task.thistrial.trialstart)*stimulus.grating.tf,1)+1);
    % reverse the phase every other grating
    if iseven(iAngle)
      phaseNum = mod(phaseNum+length(stimulus.grating.phases)/2,length(stimulus.grating.phases))+1;
    end
    
    % get what contrast to display
    % Main experiment
    if ~stimulus.isLoc
      if iAngle == task.thistrial.targetLoc
	contrastNum = find(stimulus.grating.contrasts==task.thistrial.targetContrast);
	angleNum = find(stimulus.grating.orientations==0);
      else
	contrastNum = find(stimulus.grating.contrasts==task.thistrial.distractorContrast);
      end    
      % if this is the attentive task then we put threshold level changes in orientation
      % on the target gratings
      if any(stimulus.taskType == [2 3])
	if any(iAngle ==stimulus.grating.targetLoc)
	  thisOrientation = thisOrientation+stimulus.task.orientation{task.trialnum}(iAngle);
	end
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
    mglBltTexture(stimulus.tex(contrastNum,phaseNum),[x y stimulus.grating.height],0,0,thisOrientation);
    mglBltTexture(stimulus.mask,[x y stimulus.grating.height],0,0,thisOrientation);
  end
end


% draw task reference points 
if any(stimulus.taskType == [2 3]) && any(task.thistrial.thisseg == [stimulus.stimSegment])
  for i = stimulus.grating.targetLoc
    % draw reference points
%    mglPoints2(stimulus.grating.refPoints.x{i},stimulus.grating.refPoints.y{i},4,[1 1 1]);
  end
end

% put up a mask durin response interval
if any(stimulus.taskType == [2 3]) && any(task.thistrial.thisseg ~= stimulus.stimSegment)
  for i = 1:length(stimulus.grating.angles)
    thisAngle = stimulus.grating.angles(i);
    x = stimulus.x + stimulus.grating.radius*cos(pi*thisAngle/180);
    y = stimulus.y + stimulus.grating.radius*sin(pi*thisAngle/180);
    mglBltTexture(stimulus.randMask,[x y stimulus.grating.height],0,0,0);
    mglBltTexture(stimulus.mask,[x y stimulus.grating.height],0,0,0);
  end
end
% put up fixtaion cross
if stimulus.taskType == 2
  if (task.thistrial.thisphase == 2) && any(task.thistrial.thisseg == [1 2]) 
    % put up cue
    mglLines2(stimulus.cueX(1),stimulus.cueY(1),stimulus.cueX(2),stimulus.cueY(2),stimulus.fixLineWidth,stimulus.fixColor,1);
    mglLines2(stimulus.cueX(3),stimulus.cueY(3),stimulus.cueX(4),stimulus.cueY(4),stimulus.fixLineWidth,stimulus.fixColor,1);
    
  else
    % put up fixation cross
    mglFixationCross(stimulus.fixWidth,stimulus.fixLineWidth,stimulus.fixColor);
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to set stimulus parameters at
% the beginning of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = trialResponseCallback(task,myscreen)

global stimulus;
if ~task.thistrial.gotResponse
  % get correct or incorrect
  whichButton = -1*((task.thistrial.whichButton-1)*2-1);
  thisDirection = stimulus.task.direction{task.trialnum}(stimulus.task.cuedLoc(task.trialnum));

  dispstr = sprintf('Unknown button: %i',whichButton);

  % get correct or incorrect
  if (thisDirection == whichButton)
    % correct
    stimulus.fixColor = [0 1 0];
    % update staircase
    stimulus.staircase{stimulus.task.cuedLoc(task.trialnum)}{stimulus.thisContrastNum}{stimulus.thisTargetDistractor} = upDownStaircase(stimulus.staircase{stimulus.task.cuedLoc(task.trialnum)}{stimulus.thisContrastNum}{stimulus.thisTargetDistractor},1);
    dispstr = 'correct';
  else
    % incorrect
    stimulus.fixColor = [1 0 0];
    % update staircase
    stimulus.staircase{stimulus.task.cuedLoc(task.trialnum)}{stimulus.thisContrastNum}{stimulus.thisTargetDistractor} = upDownStaircase(stimulus.staircase{stimulus.task.cuedLoc(task.trialnum)}{stimulus.thisContrastNum}{stimulus.thisTargetDistractor},0);
    dispstr = 'incorrect';
  end
else
  % observer gave multiple responses
  dispstr = sprintf('multiple responses: %i',task.thistrial.gotResponse);
end

% display threshold
disp(sprintf('(%s) Threshold: %f',dispstr,stimulus.staircase{stimulus.task.cuedLoc(task.trialnum)}{stimulus.thisContrastNum}{stimulus.thisTargetDistractor}.threshold));




