% cuecon
%
%      usage: myscreen=cuecon(stimulus)
%         by: justin gardner
%       date: 04/15/06
%    purpose: contrast discrimination cuing different numbers of stimuli
%
%
%  if training, set to 1, 2, or 3.
%  Training set will always have large threshold (0.5) and the pedestals defined below, UNLESS SPECIFIED.
%  Must specify cueCondOneOnly=1 OR cueCondOneFour=1.
%               training=1  long seglens, pedestalContrasts = [0.5]
%               training=2  long seglens, pedestalContrasts = [0.75 0.50 0.25]
%               training=3  standard seglens, pedestalContrasts = [0.75 0.50 0.25]
%  defaults:
%           threshold=0.2
%           stepsize=0.1
%           useLevittRule=1
%           stimFile=[]
%           numBlocks=12
%           pedestalContrasts=[0.0625 0.125 0.25]
%           subjectID=default
%           training=0
%           cueCondOneFour=0
%           cueCondOneOnly=4
%               (default cueConditions = {'one','two_leftRightHemi','two_upperLowerHemi','two_kittyCorners','four'})

            
function myscreen = cuecon(varargin)

taskType = [];
initStair = [];
threshold = [];
stepsize = [];
useLevittRule = [];
stimFile = [];
numBlocks = [];
pedestalContrasts = [];
subjectID = [];
cueConditions = [];
training = [];
cueCondOneFour = [];
cueCondOneOnly = [];
cueCondFourOnly = [];
getArgs(varargin,{'taskType=1','initStair=1','threshold=0.2','stepsize=0.1','useLevittRule=1','stimFile=[]','numBlocks=12','pedestalContrasts=[0.0625 0.125 0.25]','subjectID=default','training=0','cueCondOneFour=0','cueCondOneOnly=0','cueCondFourOnly=0'});

%if this is a training set...
if training > 0
    if threshold == 0.2 %and threshold is not set (is default)
        threshold = 0.5;
    end
    if isequal(pedestalContrasts,[0.0625 0.125 0.25]) %and pedestals are not set (are default)
        if training == 1 %training 1
            pedestalContrasts = [0.50];
        else %training 2 and 3
            pedestalContrasts = [0.75 0.50 0.25];
        end
    end
    if isequal(cueCondOneOnly,0) && isequal(cueCondOneFour,0) %we will not run training if cueconditions are not specified
        disp('cueCondition not set for this training set. please choose between cueCondOneOnly=1 or cueCondOneFour=1');
        return
    end
end

global stimulus;
if initStair
  clear global stimulus
  global stimulus;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other screen parameters
myscreen.background = 'gray';
myscreen = initScreen(myscreen);

% set data directory to be ~/data/subjectID
myscreen.datadir = fullfile(myscreen.datadir,'cuecon');
if ~isdir(myscreen.datadir)
  mkdir(myscreen.datadir);
end
if ~isempty(subjectID)
  myscreen.datadir = fullfile(myscreen.datadir,subjectID);
  if ~isdir(myscreen.datadir)
    mkdir(myscreen.datadir);
  end
end
disp(sprintf('(cuecon) Saving data into %s',myscreen.datadir));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open old stimfile to continue on, if asked for
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(stimFile)
  if isfile(setext(stimFile,'mat'))
    s = load(stimFile);
  elseif isfile(fullfile(myscreen.datadir,setext(stimFile,'mat')))
    s = load(fullfile(myscreen.datadir,setext(stimFile,'mat')));
  else
    disp(sprintf('(cuecon) Could not find file %s',stimFile));
    return
  end
  stimulus = s.stimulus;
  numBlocks = numBlocks - s.task{1}{2}.blocknum + 1;
  clear s;
  initStair = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = initStimulus('stimulus',myscreen);

% set the contrasts of distractor and target
stimulus.pedestalContrasts = pedestalContrasts;

% parameters
stimulus.int1 = 2;
stimulus.int2 = 4;

% grating parameters
stimulus.grating.n = 4;
stimulus.grating.orientationOfFirstGrating = 45;
stimulus.grating.radius = 6;
stimulus.colors.reservedColors = [0 0 0; 1 1 1; 0 1 0; 1 0 0;0.2 0.3 0.7];
stimulus.grating.sf = 2;
stimulus.grating.tf = 1;
stimulus.grating.width = 6;
stimulus.grating.height = 6;
stimulus.grating.nPhases = 36;

%stimulus.grating.windowType = 'gabor'; % should be gabor or thresh
%stimulus.grating.sdx = stimulus.grating.width/7;
%stimulus.grating.sdy = stimulus.grating.width/7;

stimulus.grating.windowType = 'thresh'; % should be gabor or thresh
stimulus.grating.sdx = stimulus.grating.width/2;
stimulus.grating.sdy = stimulus.grating.height/2;

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

% descriptive names for cue conditions. Add a name here for a new cue condition
% Yuko add here!!!
%stimulus.cueConditions = {'one','four'};
%stimulus.cueConditions = {'one','two_LeftRightHemi','four'};
%stimulus.cueConditions = {'one','two_leftRightHemi','two_upperLowerHemi','two_kittyCorners','four'};
if cueCondOneOnly == 1
    stimulus.cueConditions = {'one'};
elseif cueCondFourOnly == 1
    stimulus.cueConditions = {'four'};
elseif cueCondOneFour == 1
    stimulus.cueConditions = {'one', 'four'};
else
    stimulus.cueConditions = {'one','two_leftRightHemi','two_upperLowerHemi','two_kittyCorners','four'};
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up segments of trials
if taskType == 1
  stimulus.isLocalizer = 0;
  task{1}{1}.random = 1;
  task{1}{1}.segmin = [0];
  task{1}{1}.segmax = [0];
  task{1}{1}.synchToVol = [0];
  task{1}{1}.getResponse = [0];
  task{1}{1}.numTrials = 1;
  task{1}{1}.waitForBacktick = 0;
  task{1}{2}.parameter.pedestalContrast = repmat(stimulus.pedestalContrasts,stimulus.grating.n,1);
  task{1}{2}.parameter.interval = [stimulus.int1 stimulus.int2];
  task{1}{2}.parameter.targetLoc = 1:stimulus.grating.n;
  task{1}{2}.parameter.cueCondition = 1:length(stimulus.cueConditions);
  task{1}{2}.random = 1;
  if training > 0
      if training == 1 || training == 2 % training 1 and 2 is slow
        task{1}{2}.segmin = [1 1 0.3 1 2.5 1];
        task{1}{2}.segmax = [1 1 0.3 1 2.5 1];
      else% training 3 is regular speed
        task{1}{2}.segmin = [1 0.6 0.3 0.6 1.5 1];
        task{1}{2}.segmax = [1 0.6 0.3 0.6 1.5 1];
      end
  else % not training set
    task{1}{2}.segmin = [1 0.6 0.3 0.6 1.5 1];
    task{1}{2}.segmax = [1 0.6 0.3 0.6 1.5 1];
  end
  task{1}{2}.synchToVol = [0 0 0 0 0];
  task{1}{2}.getResponse = [0 0 0 0 1];
  task{1}{2}.waitForBacktick = 0;
  task{1}{2}.numBlocks = numBlocks;

  disp(sprintf('(cuecon) Number of blocks: %i',numBlocks));
  
  
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stimulus = initGratings(stimulus,myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init staircase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if initStair
  disp(sprintf('(cuecon) Initializing staircase with threshold: %f stepsize: %f useLevittRule: %i',threshold,stepsize,useLevittRule));
  stimulus = initStaircase(threshold,stimulus,stepsize,useLevittRule);
else
  disp(sprintf('(cuecon) Continuing staircase from last run'));
  dispStaircase(stimulus);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialze tasks 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

dispStaircase(stimulus,1);

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

% get number of things
stimulus.nPedestals = length(stimulus.pedestalContrasts);
stimulus.nCueConditions = length(stimulus.cueConditions);
stimulus.thisCue = 1;
% set maximum color index (for 24 bit color we have 8 bits per channel, so 255)
maxIndex = 255;

% get gamma table
if ~isfield(myscreen,'gammaTable')
  stimulus.linearizedGammaTable = mglGetGammaTable;
  disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
  disp(sprintf('(cuecon:initGratings) No gamma table found in myscreen. Contrast displays like this'));
  disp(sprintf('         should be run with a valid calibration made by moncalib for this monitor.'));
  disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
end
stimulus.linearizedGammaTable = myscreen.initScreenGammaTable;
disppercent(-inf,'Creating grating textures');

% calculate all the orientations we need
stimulus.grating.orientations = stimulus.grating.orientationOfFirstGrating:360/stimulus.grating.n:stimulus.grating.orientationOfFirstGrating+359;

% calculate all the phases we are going to compute
stimulus.grating.centerPhase = 0;
stimulus.grating.phases = [0:2*pi/(stimulus.grating.nPhases):2*pi];
stimulus.grating.phases = stimulus.grating.phases(1:end-1);
nPhases = length(stimulus.grating.phases);

% this gives the phase order to go back and forth
stimulus.grating.phaseIndex = [1:nPhases 1 nPhases:-1:2];
stimulus.grating.phaseIndexLen = length(stimulus.grating.phaseIndex);

% calculate some colors information
%  number of reserved colors
stimulus.colors.nReservedColors = size(stimulus.colors.reservedColors,1);

% number of colors possible for gratings, make sure that we 
% have an odd number
stimulus.colors.nGratingColors = maxIndex+1-stimulus.colors.nReservedColors;
if iseven(stimulus.colors.nGratingColors)
  stimulus.colors.nGratingColors = stimulus.colors.nGratingColors-1;
end
% min, mid and max index of gratings colors (index values are 0 based)
stimulus.colors.minGratingIndex = maxIndex+1-stimulus.colors.nGratingColors;
stimulus.colors.midGratingIndex = stimulus.colors.minGratingIndex+floor(stimulus.colors.nGratingColors/2);
stimulus.colors.maxGratingIndex = maxIndex;
% number of contrasts we can display (not including 0 contrast)
stimulus.colors.nDisplayContrasts = floor(stimulus.colors.nGratingColors/2);

% get the color value for gray (i.e. the number between 0 and 1 that corresponds to the midGratingIndex)
stimulus.colors.grayColor = stimulus.colors.midGratingIndex/maxIndex;

% set the reserved colors - this gives a convenient value between 0 and 1 to use the reserved colors with
for i = 1:stimulus.colors.nReservedColors
  stimulus.colors.reservedColor(i) = (i-1)/maxIndex;
end

stimulus.fixColor = stimulus.colors.reservedColor(1);
% make the window through with the gratings will be displayed
gaussianWin = mglMakeGaussian(stimulus.grating.width,stimulus.grating.height,stimulus.grating.sdx,stimulus.grating.sdy);
if strcmp(stimulus.grating.windowType,'gabor')
  % a gaussian window
  win = maxIndex-maxIndex*gaussianWin;
else
  % a simple window
  win = maxIndex-maxIndex*(gaussianWin>=exp(-1/2));
end
mask = ones(size(win,1),size(win,2),4)*stimulus.colors.midGratingIndex;
mask(:,:,4) = win;
stimulus.mask = mglCreateTexture(mask);

% make all the 1D gratings. We compute all phases and all possible contrast values given the
% range of indexes available to us. The 1st texture is gray the nth texture is full
% contrast for the current gamma setting
for iPhase = 1:nPhases
  for iContrast = 0:stimulus.colors.nDisplayContrasts
    pDone = calcPercentDone(iPhase,nPhases,iContrast,stimulus.colors.nDisplayContrasts);
    disppercent(pDone);
    if myscreen.userHitEsc,mglClose;keyboard,end
    % get the phase
    thisPhase = (stimulus.grating.centerPhase+stimulus.grating.phases(iPhase))*180/pi;
    % make the grating
    thisGrating = round(iContrast*mglMakeGrating(stimulus.grating.width,nan,stimulus.grating.sf,0,thisPhase)+stimulus.colors.midGratingIndex);
    if min(thisGrating(:)) < stimulus.colors.minGratingIndex
      keyboard
    end
    % create the texture
    stimulus.tex(iContrast+1,iPhase) = mglCreateTexture(thisGrating);
  end
end
disppercent(inf);

% set up cuelines
for i = 1:stimulus.grating.n
  stimulus.grating.cueLines(i,1) = cos(d2r(stimulus.grating.orientations(i)))*0.5;
  stimulus.grating.cueLines(i,2) = sin(d2r(stimulus.grating.orientations(i)))*0.5;
end

% set up reference lines
stimulus.grating.refLines.x1 = [];
stimulus.grating.refLines.y1 = [];
stimulus.grating.refLines.x2 = [];
stimulus.grating.refLines.y2 = [];

for iAngle = 1:length(stimulus.grating.orientations)
  % get center of patch
  thisAngle = stimulus.grating.orientations(iAngle);
  centerX = stimulus.grating.radius*cos(pi*thisAngle/180);
  centerY = stimulus.grating.radius*sin(pi*thisAngle/180);
  % get radius
  radius = sqrt(((stimulus.grating.width/2)^2)+((stimulus.grating.height/2)^2))-0.5;
  % get left/right top/bottom;
  left = centerX-stimulus.grating.width/2;
  right = centerX+stimulus.grating.width/2;
  top = centerY-stimulus.grating.height/2;
  bottom = centerY+stimulus.grating.height/2;
  % square reference points 
  %stimulus.grating.refLines.x1 = [stimulus.grating.refLines.x1 left left right right];
  %stimulus.grating.refLines.y1 = [stimulus.grating.refLines.y1 top bottom bottom top];
  %stimulus.grating.refLines.x2 = [stimulus.grating.refLines.x2 left right right left];
  %stimulus.grating.refLines.y2 = [stimulus.grating.refLines.y2 bottom bottom top top];
  
  % circular reference lines
  d = 0:1:360;
  for dIndex = 1:length(d)-1
    stimulus.grating.refLines.x1 = [stimulus.grating.refLines.x1 centerX+radius*cos(d2r(d(dIndex)))];
    stimulus.grating.refLines.y1 = [stimulus.grating.refLines.y1 centerY+radius*sin(d2r(d(dIndex)))];
    stimulus.grating.refLines.x2 = [stimulus.grating.refLines.x2 centerX+radius*cos(d2r(d(dIndex+1)))];
    stimulus.grating.refLines.y2 = [stimulus.grating.refLines.y2 centerY+radius*sin(d2r(d(dIndex+1)))];
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

if (task.thistrial.thisphase == 2) && (task.thistrial.thisseg == 1)
  % get current threshold from staircase
  stimulus.pedestalNum = find(task.thistrial.pedestalContrast(task.thistrial.targetLoc)==stimulus.pedestalContrasts);
  stimulus.cueNum = task.thistrial.cueCondition;
  stimulus.deltaContrast(task.trialnum) = stimulus.staircase{stimulus.pedestalNum}{stimulus.cueNum}.threshold;
  % check to see if we go over 100% contrast
  if (task.thistrial.pedestalContrast(task.thistrial.targetLoc)+stimulus.deltaContrast(end)) > 1
    stimulus.deltaContrast(end) = 1-task.thistrial.pedestalContrast(task.thistrial.targetLoc);
    disp(sprintf('(cuecon) !!!! Delta contast (%f) + Pedestal Contrast (%f) is greater than 100% contrast',stimulus.deltaContrast(end),task.thistrial.pedestalContrast(task.thistrial.targetLoc)));
  end
  % set the maximum contrast we can display
  setGammaTableForMaxContrast(max([task.thistrial.pedestalContrast; task.thistrial.pedestalContrast(task.thistrial.targetLoc)+stimulus.deltaContrast(end)]));
  % choose which cue to show
  % Yuko hey!!!
  switch stimulus.cueConditions{task.thistrial.cueCondition} 
    case {'one'}
      stimulus.thisCue = task.thistrial.targetLoc;
    case {'four'}
      stimulus.thisCue = [1 2 3 4];
    %upR=1, upL=2, loL=3, loR=4
    case {'two_leftRightHemi'} %left vs. right (made a mistake and named this cueCondition two_sameHemi for 100415_stim01 and stim02.
      if task.thistrial.targetLoc == 1
        stimulus.thisCue = [1 4];
      elseif task.thistrial.targetLoc == 2
        stimulus.thisCue = [2 3];
      elseif task.thistrial.targetLoc == 3
        stimulus.thisCue = [2 3];
      elseif task.thistrial.targetLoc == 4
        stimulus.thisCue = [1 4];
      end
    case {'two_upperLowerHemi'} %upper vs. lower
      if task.thistrial.targetLoc == 1
        stimulus.thisCue = [1 2];
      elseif task.thistrial.targetLoc == 2
        stimulus.thisCue = [1 2];
      elseif task.thistrial.targetLoc == 3
        stimulus.thisCue = [3 4];
      elseif task.thistrial.targetLoc == 4
        stimulus.thisCue = [3 4];
      end      
    case {'two_kittyCorners'} %opp diagonal corners
      if task.thistrial.targetLoc == 1
        stimulus.thisCue = [1 3];
      elseif task.thistrial.targetLoc == 2
        stimulus.thisCue = [2 4];
      elseif task.thistrial.targetLoc == 3
        stimulus.thisCue = [1 3];
      elseif task.thistrial.targetLoc == 4
        stimulus.thisCue = [2 4];
      end
  end
end

if any(task.thistrial.thisseg == [stimulus.int1 stimulus.int2])
  % get the contrast for the stimuli in this interval
  stimulusContrast = task.thistrial.pedestalContrast;
  deltaContrast = stimulus.deltaContrast(end).*(task.thistrial.interval==task.thistrial.thisseg);
  stimulusContrast(task.thistrial.targetLoc) = stimulusContrast(task.thistrial.targetLoc)+deltaContrast;

  % set up the stimuli
  for i = 1:stimulus.grating.n
    % orientation of stimulus
    stimulus.o(i) = stimulus.grating.orientations(i);

    % get x,y position of grating
    stimulus.x(i) = stimulus.grating.radius*cos(d2r(stimulus.o(i)));
    stimulus.y(i) = stimulus.grating.radius*sin(d2r(stimulus.o(i)));
    
    % get the contrasts
    stimulus.contrastIndex(i) = getContrastIndex(stimulusContrast(i));
  end
end

% set the fixation color
if any(task.thistrial.thisseg == [stimulus.int1 stimulus.int2])
  stimulus.fixColor = stimulus.colors.reservedColor(1);
else
  stimulus.fixColor = stimulus.colors.reservedColor(2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to display stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task,myscreen)

global stimulus;
mglClearScreen(stimulus.colors.grayColor);

if any(task.thistrial.thisseg == [stimulus.int1 stimulus.int2]) 
  for i = 1:stimulus.grating.n
    % get which phase we are on, note that the phaseIndex goes through approximately 2 cycles
    % of the stimulus so we adjust the tf accordingly
    cycleTime = rem(mglGetSecs(task.thistrial.trialstart)*stimulus.grating.tf*(stimulus.grating.nPhases+1)/(stimulus.grating.phaseIndexLen+1),1);
    phaseIndex = floor(stimulus.grating.phaseIndexLen*cycleTime)+1;

    % reverse the phase every other grating
    if iseven(i)
      phaseIndex = mod(phaseIndex+stimulus.grating.phaseIndexLen/2,stimulus.grating.phaseIndexLen)+1;
    end
    phaseIndex = stimulus.grating.phaseIndex(phaseIndex);
    
    % blt texture
    mglBltTexture(stimulus.tex(stimulus.contrastIndex(i),phaseIndex),[stimulus.x(i) stimulus.y(i) stimulus.grating.height],0,0,stimulus.o(i));
    mglBltTexture(stimulus.mask,[stimulus.x(i) stimulus.y(i)],0,0,stimulus.o(i));
  end
end

% draw reference points
mglLines2(stimulus.grating.refLines.x1,stimulus.grating.refLines.y1,stimulus.grating.refLines.x2,stimulus.grating.refLines.y2,1,stimulus.colors.reservedColor(2));

if task.thistrial.thisseg == 5
  % draw target
  mglLines2(0,0,stimulus.grating.cueLines(task.thistrial.targetLoc,1),stimulus.grating.cueLines(task.thistrial.targetLoc,2),1,stimulus.colors.reservedColor(3));
elseif (task.thistrial.thisphase == 2) && (task.thistrial.thisseg < 5)
  % draw cues
  mglLines2(zeros(length(stimulus.thisCue),1),zeros(length(stimulus.thisCue),1),stimulus.grating.cueLines(stimulus.thisCue,1),stimulus.grating.cueLines(stimulus.thisCue,2),1,stimulus.colors.reservedColor(2));
end
mglFixationCross(0.5,1,stimulus.fixColor);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to set stimulus parameters at
% the beginning of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = trialResponseCallback(task,myscreen)

global stimulus;

whichInterval = find(task.thistrial.interval == task.parameter.interval);
if (task.thistrial.whichButton == whichInterval)
  correctIncorrect = 'correct';
  stimulus.fixColor = stimulus.colors.reservedColor(3);
  keyboard
  stimulus.staircase{stimulus.pedestalNum}{stimulus.cueNum} = upDownStaircase(stimulus.staircase{stimulus.pedestalNum}{stimulus.cueNum},1);
else
  correctIncorrect = 'incorrect';
  stimulus.fixColor = stimulus.colors.reservedColor(4);
  stimulus.staircase{stimulus.pedestalNum}{stimulus.cueNum} = upDownStaircase(stimulus.staircase{stimulus.pedestalNum}{stimulus.cueNum},0);
end
disp(sprintf('Cue: %s pedestal: %f deltaC: %f (%s)',stimulus.cueConditions{task.thistrial.cueCondition},task.thistrial.pedestalContrast(task.thistrial.targetLoc),stimulus.deltaContrast(end),correctIncorrect));
  
%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getContrastIndex    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function contrastIndex = getContrastIndex(desiredContrast,verbose)

if nargin < 2,verbose = 0;end
if desiredContrast < 0, desiredContrast = 0;end

global stimulus;

% now find closest matching contrast we can display with this gamma table
contrastIndex = min(round(stimulus.colors.nDisplayContrasts*desiredContrast/stimulus.currentMaxContrast),stimulus.colors.nDisplayContrasts);

% display the desired and actual contrast values if verbose is set
if verbose
  actualContrast = stimulus.currentMaxContrast*(contrastIndex/stimulus.colors.nDisplayContrasts);
  disp(sprintf('(getContrastIndex) Desired contrast: %0.4f Actual contrast: %0.4f Difference: %0.4f',desiredContrast,actualContrast,desiredContrast-actualContrast));
end

% out of range check
if round(stimulus.colors.nDisplayContrasts*desiredContrast/stimulus.currentMaxContrast)>stimulus.colors.nDisplayContrasts
 disp(sprintf('(getContrastIndex) Desired contrast (%0.9f) out of range max contrast : %0.9f',desiredContrast,stimulus.currentMaxContrast));
 keyboard
end

% 1 based indexes (0th index is gray, nDisplayContrasts+1 is full contrast)
contrastIndex = contrastIndex+1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sets the gamma table so that we can have
% finest possible control over the stimulus contrast.
%
% stimulus.colors.reservedColors should be set to the reserved colors (for cue colors, etc).
% maxContrast is the maximum contrast you want to be able to display.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setGammaTableForMaxContrast(maxContrast)

global stimulus;

% set the reserved colors
gammaTable(1:size(stimulus.colors.reservedColors,1),1:size(stimulus.colors.reservedColors,2))=stimulus.colors.reservedColors;

% set the gamma table
if maxContrast > 0
  % create the rest of the gamma table
  cmax = 0.5+maxContrast/2;cmin = 0.5-maxContrast/2;
  luminanceVals = cmin:((cmax-cmin)/(stimulus.colors.nGratingColors-1)):cmax;

  % now get the linearized range
  redLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,luminanceVals,'linear');
  greenLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,luminanceVals,'linear');
  blueLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,luminanceVals,'linear');
  
  % add these values to the table


  gammaTable((stimulus.colors.minGratingIndex:stimulus.colors.maxGratingIndex)+1,:)=[redLinearized;greenLinearized;blueLinearized]';
else
  % if we are asked for 0 contrast then simply set all the values to gray
  gammaTable((stimulus.colors.minGratingIndex:stimulus.colors.maxGratingIndex)+1,1)=interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,0.5,'linear');
  gammaTable((stimulus.colors.minGratingIndex:stimulus.colors.maxGratingIndex)+1,2)=interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,0.5,'linear');
  gammaTable((stimulus.colors.minGratingIndex:stimulus.colors.maxGratingIndex)+1,3)=interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,0.5,'linear');
end

% set the gamma table
mglSetGammaTable(gammaTable);

% remember what the current maximum contrast is that we can display
stimulus.currentMaxContrast = maxContrast;

%%%%%%%%%%%%%%%%%%%%%%%%
%    startStaircase    %
%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initStaircase(threshold,stimulus,stepsize,useLevittRule)

for iPedestal = 1:stimulus.nPedestals
  for iCueCondition = 1:stimulus.nCueConditions
    stimulus.staircase{iPedestal}{iCueCondition} = upDownStaircase(1,2,threshold,stepsize,useLevittRule);
    stimulus.staircase{iPedestal}{iCueCondition}.minThreshold = 0;
    stimulus.staircase{iPedestal}{iCueCondition}.maxThreshold = 1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%
%    dispStaircase    %
%%%%%%%%%%%%%%%%%%%%%%%
function dispStaircase(stimulus,makeplot)

if nargin < 2,makeplot = 0;end
%if makeplot,smartfig('cuecon','reuse');end

for iPedestal = 1:stimulus.nPedestals
  for iCueCondition = 1:stimulus.nCueConditions
    s = stimulus.staircase{iPedestal}{iCueCondition};
    if isfield(s,'strength')
      n = length(s.strength);
    else
      n = 0;
    end
    disp(sprintf('(Contrast: %f CueCondition: %s): %f (n=%i)',stimulus.pedestalContrasts(iPedestal),stimulus.cueConditions{iCueCondition},s.threshold,n));
  end
end
