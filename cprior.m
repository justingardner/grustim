% cprior
%
%      usage: myscreen=cprior
%         by: justin gardner
%       date: 01/11/2011
%    purpose: contrast prior experiment. Code based con cuecon
%
function myscreen = cprior(varargin)

taskType = [];
initStair = [];
stimFile = [];
numBlocks = [];
pedestalContrasts = [];
subjectID = [];
cueConditions = [];
training = [];
getArgs(varargin,{'taskType=1','initStair=1','stimFile=[]','numBlocks=12','pedestalContrasts=[0.0625 0.125 0.25]','subjectID=[]','training=0'});

%i Training run
if training > 0
  disp(sprintf('(cprior) Training not implemented yet!'));
  keyboard
end

global stimulus;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other screen parameters
myscreen.background = 'gray';
myscreen.subjectID = subjectID;
myscreen = initScreen(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open old stimfile to continue on, if asked for
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(stimFile)
  if isfile(setext(stimFile,'mat'))
    s = load(stimFile);
  elseif isfile(fullfile(myscreen.datadir,setext(stimFile,'mat')))
    s = load(fullfile(myscreen.datadir,setext(stimFile,'mat')));
  else
    disp(sprintf('(cprior) Could not find file %s',stimFile));
    return
  end
  disp(sprintf('(cprior) Loading from stimfile not implemented yet!'));
  keyboard
  %  stimulus.staircase = s.stimulus.staircase;
%  numBlocks = numBlocks - s.task{1}{2}.blocknum + 1;
%  clear s;
%  initStair = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = initStimulus('stimulus',myscreen);

% set the contrasts of distractor and target
stimulus.pedestalContrasts = pedestalContrasts;

% grating parameters
stimulus.grating.n = 2;
stimulus.grating.orientationOfFirstGrating = 0;
stimulus.grating.radius = 10;
stimulus.colors.reservedColors = [0 0 0; 1 1 1; 0 1 0; 1 0 0;0.2 0.3 0.7];
stimulus.grating.sf = 2;
stimulus.grating.tf = 1;
stimulus.grating.width = 10;
stimulus.grating.height = 10;
stimulus.grating.nPhases = 36;
stimulus.grating.windowType = 'thresh'; % should be gabor or thresh
stimulus.grating.sdx = stimulus.grating.width/2;
stimulus.grating.sdy = stimulus.grating.height/2;

% other stimulus parameters
stimulus.interval = 2;
stimulus.priorProb = [0.9 0.6];
stimulus.cueColor = [3 4];
stimulus.nProb = length(stimulus.priorProb);
stimulus.nPedestals = length(stimulus.pedestalContrasts);

% init the stimulus with these settings
stimulus = initGratings(stimulus,myscreen);

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

  task{1}{2}.parameter.pedestalContrast = stimulus.pedestalContrasts;
  task{1}{2}.parameter.targetLoc = 1:stimulus.grating.n;
  task{1}{2}.parameter.priorProb = stimulus.priorProb;
  task{1}{2}.randVars.calculated.validity = nan;
  task{1}{2}.randVars.calculated.deltaContrast = nan;
  task{1}{2}.randVars.calculated.distractorLoc = nan;
  task{1}{2}.random = 1;
  task{1}{2}.segmin = [0.5 3 1.5];
  task{1}{2}.segmax = [0.5 3 1.5];
  task{1}{2}.synchToVol = [0 0 0];
  task{1}{2}.getResponse = [0 0 1];
  task{1}{2}.waitForBacktick = 0;
  task{1}{2}.numBlocks = numBlocks;

  disp(sprintf('(cprior) Number of blocks: %i',numBlocks));
  
  
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
% init staircase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if initStair
  disp(sprintf('(cprior) Initializing staircase'));
  stimulus.staircase = [];
  % one staircase for each pedestal and each probability condition and
  % for whether the cue was valid or invalid.
  for iPedestal = 1:stimulus.nPedestals
    for iProb = 1:stimulus.nProb
      for iValidity = 1:2
	stimulus.staircase{iPedestal}{iProb}{iValidity} = doStaircase('init','fixed','fixedVals=[0 0.01 0.02 0.04 0.08 0.16]');
      end
    end
  end
else
  disp(sprintf('(cprior) Continuing staircase from last run'));
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
phaseNum = 1;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
  % update the task
  [task{1} myscreen phaseNum] = updateTask(task{1},myscreen,phaseNum);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);
mglSetGammaTable(mglGetParam('initialGammaTable'));

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
function stimulus = initGratings(stimulus,myscreen)

% set maximum color index (for 24 bit color we have 8 bits per channel, so 255)
maxIndex = 255;

% get gamma table
if ~isfield(myscreen,'gammaTable')
  stimulus.linearizedGammaTable = mglGetGammaTable;
  disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
  disp(sprintf('(cprior:initGratings) No gamma table found in myscreen. Contrast displays like this'));
  disp(sprintf('         should be run with a valid calibration made by moncalib for this monitor.'));
  disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
end
stimulus.linearizedGammaTable = myscreen.initScreenGammaTable;

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
% remember the GLContext, this is usefult to test whether we need to rebuild textures or not
stimulus.grating.GLContext = mglGetParam('GLContext');
if isfield(stimulus,'oldGrating') && isequal(stimulus.grating.GLContext,stimulus.oldGrating.GLContext)
  sameContextAsLastRun = true;
else
  sameContextAsLastRun = false;
end
  
% decide whether to make textures or not
createGratings = false;
deleteGratings = false;
if ~isfield(stimulus,'oldGrating') || ~isfield(stimulus,'tex') || isempty(stimulus.tex)
  % if oldGrating has not been set in stimulus then need to create
  createGratings = true;
else
  % otherwise check that none of the following have changed, otherwise we need to remake
  fieldsThatMustBeTheSame = {'GLContext','orientationOfFirstGrating','sf','tf','width','height','nPhases'};
  for i = 1:length(fieldsThatMustBeTheSame)
    if ~isequal(stimulus.grating.(fieldsThatMustBeTheSame{i}),stimulus.oldGrating.(fieldsThatMustBeTheSame{i}))
      createGratings = true;
      deleteGratings = true;
      break;
    end
  end
end

if deleteGratings
  % delete old gratings if they exist
  if sameContextAsLastRun && isfield(stimulus,'tex')
    disppercent(-inf,'(cprior) Deleting old grating textures');
    for i = 1:length(stimulus.tex(:))
      disppercent(i/length(stimulus.tex(:)));
      mglDeleteTexture(stimulus.tex(i));
    end
  end  
  disppercent(inf);
end

if createGratings
  disppercent(-inf,'(cprior) Creating grating textures');
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
else
  disp(sprintf('(cprior) Using gratings from last run'));
end

% remove old mask if it is there
if sameContextAsLastRun && isfield(stimulus,'mask') && ~isempty(stimulus.mask)
  mglDeleteTexture(stimulus.mask);
end

% make the window through with the gratings will be displayed
gaussianWin = mglMakeGaussian(stimulus.grating.width+0.5,stimulus.grating.height+0.5,stimulus.grating.sdx,stimulus.grating.sdy);
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

% save current settings
stimulus.oldGrating = stimulus.grating;

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
  % set whether cue is valid or not on this try
  task.thistrial.validity = (task.thistrial.priorProb > rand);
  % get the distractor location (a random location out of the n), for the case of 2 targets
  % then is always the other location
  task.thistrial.distractorLoc = setxor(task.thistrial.targetLoc,randperm(1:stimulus.grating.n));
  task.thistrial.distractorLoc = task.thistrial.distractorLoc(first(randperm(length(task.thistrial.distractorLoc))));
  % get deltaContrast to test from relevant staircase
  stimulus.pedestalNum = find(task.thistrial.pedestalContrast==stimulus.pedestalContrasts);
  stimulus.probNum = find(task.thistrial.priorProb==task.parameter.priorProb);
  stimulus.validityNum = task.thistrial.validity+1;
  % now get the delta contrast for this trial
  [task.thistrial.deltaContrast stimulus.staircase{stimulus.pedestalNum}{stimulus.probNum}{stimulus.validityNum}] = doStaircase('testValue',stimulus.staircase{stimulus.pedestalNum}{stimulus.probNum}{stimulus.validityNum});
  % check to see if we go over 100% contrast
  if (task.thistrial.pedestalContrast+task.thistrial.deltaContrast) > 1
    task.thistrial.deltaContrast = 1-task.thistrial.pedestalContrast;
    disp(sprintf('(cprior) !!!! Delta contast (%f) + Pedestal Contrast (%f) is greater than 100% contrast',task.thistrial.deltaContrast,task.thistrial.pedestalContrast));
  end
  % set the maximum contrast we can display
  setGammaTableForMaxContrast(task.thistrial.pedestalContrast+task.thistrial.deltaContrast);
  % display some stuff
  disp(sprintf('(cprior) pedestal: %i:%0.3f priorProb: %i:%0.3f validity: %i targetLoc: %i distractorLoc: %i deltaContrast:%f',stimulus.pedestalNum,100*task.thistrial.pedestalContrast,stimulus.probNum,task.thistrial.priorProb,task.thistrial.validity,task.thistrial.targetLoc,task.thistrial.distractorLoc,task.thistrial.deltaContrast));
end

if any(task.thistrial.thisseg == stimulus.interval)
  % get the contrast for the stimuli in this interval
  stimulusContrast(1:stimulus.grating.n) = task.thistrial.pedestalContrast;
  stimulusContrast(task.thistrial.targetLoc) = stimulusContrast(task.thistrial.targetLoc)+task.thistrial.deltaContrast;

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
if any(task.thistrial.thisseg == [stimulus.interval])
  stimulus.fixColor = stimulus.colors.reservedColor(1);
else
  stimulus.fixColor = stimulus.colors.reservedColor(2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to display stimulus (on every refresh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task,myscreen)

global stimulus;
mglClearScreen(stimulus.colors.grayColor);

% for stimulus interval
if any(task.thistrial.thisseg == stimulus.interval) 
  for i = 1:stimulus.grating.n
    % get which phase we are on, note that the phaseIndex goes through approximately 2 cycles
    % of the stimulus so we adjust the tf accordingly
    cycleTime = rem(mglGetSecs(task.thistrial.trialstart)*stimulus.grating.tf*(stimulus.grating.nPhases+1)/(stimulus.grating.phaseIndexLen+1),1);
    phaseIndex = floor(stimulus.grating.phaseIndexLen*cycleTime)+1;

    % reverse the phase every other grating
%    if iseven(i),phaseIndex = mod(phaseIndex+stimulus.grating.phaseIndexLen/2,stimulus.grating.phaseIndexLen)+1;,end
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
elseif (task.thistrial.thisphase == 2) && (task.thistrial.thisseg < 5)
  if task.thistrial.validity
    % draw the valid cue
    mglLines2(0,0,stimulus.grating.cueLines(task.thistrial.targetLoc,1),stimulus.grating.cueLines(task.thistrial.targetLoc,2),1,stimulus.colors.reservedColor(3));
  else
    % draw the invalid cue
    mglLines2(0,0,stimulus.grating.cueLines(task.thistrial.distractorLoc,1),stimulus.grating.cueLines(task.thistrial.distractorLoc,2),1,stimulus.colors.reservedColor(3));
  end
end
mglFixationCross(1,1,stimulus.fixColor);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function response callback 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = trialResponseCallback(task,myscreen)

global stimulus;
if task.thistrial.gotResponse == 0
  whichInterval = find(task.thistrial.interval == task.parameter.interval);
  if (task.thistrial.whichButton == whichInterval)
    correctIncorrect = 'correct';
    stimulus.fixColor = stimulus.colors.reservedColor(3);
    stimulus.staircase{stimulus.pedestalNum}{stimulus.staircaseNum} = upDownStaircase(stimulus.staircase{stimulus.pedestalNum}{stimulus.staircaseNum},1);
  else
    correctIncorrect = 'incorrect';
    stimulus.fixColor = stimulus.colors.reservedColor(4);
    stimulus.staircase{stimulus.pedestalNum}{stimulus.staircaseNum} = upDownStaircase(stimulus.staircase{stimulus.pedestalNum}{stimulus.staircaseNum},0);
  end
    disp(sprintf('Cue: %s pedestal: %f deltaC: %f (%s)',stimulus.cueConditions{task.thistrial.cueCondition},task.thistrial.pedestalContrast(task.thistrial.targetLoc),stimulus.deltaContrast(task.trialnum),correctIncorrect));
else
  disp(sprintf('Subject responded multiple times: %i',task.thistrial.gotResponse+1));
end
  
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


%%%%%%%%%%%%%%%%%%%%%%%
%    dispStaircase    %
%%%%%%%%%%%%%%%%%%%%%%%
function dispStaircase(stimulus,makeplot)

if nargin < 2,makeplot = 0;end

%if makeplot,smartfig('cuecon','reuse');end

for iPedestal = 1:stimulus.nPedestals
  for iProb = 1:stimulus.nProb
    for iValid = 1:2
      s = stimulus.staircase{iPedestal}{iProb}{iValid};
%      disp(sprintf('(Contrast: %f CueCondition: %i): %f
%      (n=%i)',stimulus.pedestalContrasts(iPedestal),stimulus.priorProb(iStaircase),s.threshold,n));
    end
  end
end
