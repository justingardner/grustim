% liaoTemplate.m
%
%        $Id$
%      usage: liaoTemplate
%         by: justin gardner
%       date: 07/97/2017
%  copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%    purpose: Template for contrast discrimination using dots at different eccentricities
%
function myscreen = liaoTemplate(varargin)

% check arguments
getArgs(varargin,'stimulusType=dots');

% initalize the screen
myscreen = initScreen;

% init the stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);
stimulus.stimulusType = lower(stimulusType);

% stimulus parameters
stimulus.width = 5;
stimulus.eccentricity = [5 10];

% init the stimulus
if strcmp(stimulus.stimulusType,'dots')
  % init the dots
  stimulus = initDotsStimulus(stimulus, myscreen);
  stimulus.isGrating = false;
elseif strcmp(stimulus.stimulusType,'grating')
  % init the grating
  stimulus = initGratingStimulus(stimulus, myscreen);
  stimulus.isGrating = true;
else
  disp(sprintf('(liaoTemplate) Unknown stimulus type: %s',stimulus.stimulusType));
  endScreen;
  return
end

% init the staircases
stimulus.pedestalContrast = 0.5;
initialThreshold = 0.2;
initialStepsize = 0.05;
% displays the staircase as trial data comes in a figure
dispStaircaseFig = 1;
% number of trials per staircase
nTrialsPerStaircase = 10;
% how many staircases to run until program ends
stimulus.nStaircases = 2;
% whether to try to start from the last computed threshold
useLastThreshold = false;
% call with the above parameters
stimulus = initStaircases(stimulus, myscreen, initialThreshold, initialStepsize, nTrialsPerStaircase, dispStaircaseFig, useLastThreshold);

% fix: set waitForBacktick if you want to synch with the scanner
% by waiting for the backtick key to be pressed before starting the experiment
% (for systems that use NI digital I/O, this will wait for the digital
% signal that the scanner has started collecting data)
task{1}.waitForBacktick = 0;

% task parameters
task{1}.segmin = [0.5 0.5 2];
task{1}.segmax = [0.5 0.5 2];
task{1}.getResponse = [0 1 1];
task{1}.parameter.eccentricity = stimulus.eccentricity;
task{1}.randVars.uniform.whichSide = [1 2];
task{1}.random = 1;

% initialize the task
for phaseNum = 1:length(task)
  [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while ~all(stimulus.staircaseCompleted == stimulus.nStaircases) && ~myscreen.userHitEsc
  % update the task
  [task myscreen phaseNum] = updateTask(task,myscreen,phaseNum);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

if ~myscreen.userHitEsc
  % display thresholds
  for iStaircase = 1:length(stimulus.eccentricity);
    % compute threshold
    t = doStaircase('threshold',stimulus.s(iStaircase,:));
    % display it
    disp(sprintf('Eccentricity: %0.2f Threshold: %0.2f',stimulus.eccentricity(iStaircase),t.threshold));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;

if task.thistrial.thisseg == 1
  % set the eccentricity of the dots
  % get which staircase
  stimulus.thisStaircase = find(stimulus.eccentricity==task.thistrial.eccentricity);
  stimulus.thisStaircase(2) = stimulus.staircaseCompleted(stimulus.thisStaircase)+1;

  % get the delta contrast to test
  [delta stimulus.s(stimulus.thisStaircase(1),stimulus.thisStaircase(2))] = doStaircase('testValue',stimulus.s(stimulus.thisStaircase(1),stimulus.thisStaircase(2)));

  % make sure delta does not go below zero
  stimulus.delta = max(0,delta);
  
  % get contrast for left and right
  if task.thistrial.whichSide == 1
    leftContrast = stimulus.pedestalContrast+stimulus.delta;
    rightContrast = stimulus.pedestalContrast;
  else
    leftContrast = stimulus.pedestalContrast;
    rightContrast = stimulus.pedestalContrast+stimulus.delta;
  end
  
  % display what is going on
  disp(sprintf('%i: Eccentricity: %0.1f Side: %i delta: %0.2f',task.trialnum,task.thistrial.eccentricity,task.thistrial.whichSide,stimulus.delta));
  
  % grating
  if stimulus.isGrating
    % set the gamma table
    setGammaTableForMaxContrast(stimulus.pedestalContrast+stimulus.delta);
    % set the grating indexes
    stimulus.leftContrastIndex = getContrastIndex(leftContrast);
    stimulus.rightContrastIndex = getContrastIndex(rightContrast);
  % dots
  else
    % set the contrast
    stimulus.dotsLeft = stimulus.dotsLeft.setContrast(stimulus.dotsLeft,leftContrast);
    stimulus.dotsRight = stimulus.dotsRight.setContrast(stimulus.dotsRight,rightContrast);

    % set the position
    stimulus.dotsLeft = stimulus.dotsLeft.setCenter(stimulus.dotsLeft,-task.thistrial.eccentricity,0);
    stimulus.dotsRight = stimulus.dotsLeft.setCenter(stimulus.dotsRight,task.thistrial.eccentricity,0);
  end

  % set the fixation color
  stimulus.fixColor = stimulus.normalFixColor;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus

% clear screen to gray
mglClearScreen(stimulus.backgroundColor);

if task.thistrial.thisseg == 2
  % grating stimulus
  if stimulus.isGrating
    % blt gratings
    mglBltTexture(stimulus.tex(stimulus.leftContrastIndex),[-task.thistrial.eccentricity 0 stimulus.grating.height]);
    mglBltTexture(stimulus.tex(stimulus.rightContrastIndex),[task.thistrial.eccentricity 0 stimulus.grating.height]);

    % and mask with the gaussian
    mglBltTexture(stimulus.mask,[-task.thistrial.eccentricity 0]);
    mglBltTexture(stimulus.mask,[task.thistrial.eccentricity 0]);

  % dots stimulus
  else
    % update the dots
    stimulus.dotsLeft = stimulus.dotsLeft.update(stimulus.dotsLeft);
    stimulus.dotsRight = stimulus.dotsRight.update(stimulus.dotsRight);

    % draw the dots
    stimulus.dotsLeft = stimulus.dotsLeft.draw(stimulus.dotsLeft);
    stimulus.dotsRight = stimulus.dotsRight.draw(stimulus.dotsRight);
  end
end

% draw fixation cross
mglFixationCross(1,1,stimulus.fixColor);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)

global stimulus

% check the response
if task.thistrial.gotResponse < 1
  % see if it is correct
  if isequal(task.thistrial.whichButton,task.thistrial.whichSide)
    % report answer
    disp(sprintf(' !! Correct !!. Reaction time: %0.2f',task.thistrial.reactionTime));
    % change fixation color
    stimulus.fixColor = stimulus.correctFixColor;
    % and update staircase
    stimulus.s(stimulus.thisStaircase(1),stimulus.thisStaircase(2)) = doStaircase('update',stimulus.s(stimulus.thisStaircase(1),stimulus.thisStaircase(2)),1,stimulus.delta);
  else
    % report answer
    disp(sprintf(' ++ Incorrect ++. Reaction time: %0.2f',task.thistrial.reactionTime));
    % change fixation color
    stimulus.fixColor = stimulus.incorrectFixColor;
    % and update staircase
    stimulus.s(stimulus.thisStaircase(1),stimulus.thisStaircase(2)) = doStaircase('update',stimulus.s(stimulus.thisStaircase(1),stimulus.thisStaircase(2)),0,stimulus.delta);
  end    
  % see if we are done
  if doStaircase('stop',stimulus.s(stimulus.thisStaircase(1),stimulus.thisStaircase(2)))
    % update variable that says we are done
    stimulus.staircaseCompleted(stimulus.thisStaircase(1)) = stimulus.staircaseCompleted(stimulus.thisStaircase(1)) + 1;
    % and initialze the next staircase
    if stimulus.thisStaircase(2) < stimulus.nStaircases
      stimulus.s(stimulus.thisStaircase(1),stimulus.thisStaircase(2)+1) = doStaircase('init',stimulus.s(stimulus.thisStaircase(1),stimulus.thisStaircase(2)));
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dots stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initDotsStimulus(stimulus,myscreen)

% init the dot patchs
stimulus.dotsLeft = dotsInit('framesPerSecond',myscreen.framesPerSecond,'dir=0','width',stimulus.width);
stimulus.dotsRight = dotsInit('framesPerSecond',myscreen.framesPerSecond,'dir=180','width',stimulus.width);

% set background color
stimulus.backgroundColor = 0.5;

% and fixation color
stimulus.normalFixColor = [1 1 1];
stimulus.correctFixColor = [0 1 0];
stimulus.incorrectFixColor = [1 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the staircases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initStaircases(stimulus, myscreen,initialThreshold,initialStepsize,nTrials,dispStaircaseFig,useLastThreshold)

if ~useLastThreshold
  % if we are not using last threshold, 
  % then don't bother trying to load last stimfile
  stimfile = [];
else
  % get the last stimfile
  stimfile = getLastStimfile(myscreen);
end

% check that the last stimfile had the same eccentricites
if ~isempty(stimfile)
  if ~isfield(stimfile,'stimulus') || ~isfield(stimfile.stimulus,'eccentricity') || ~isequal(stimfile.stimulus.eccentricity,stimulus.eccentricity)
    % dump this stimfile, since it does not match current eccentricity
    disp(sprintf('(liaoTemplate) Found stimfile, but does not have eccentricty match. Ignoring'));
    stimfile = [];
  end
end

% if no stimfile found
if isempty(stimfile)
  % then initialize
  for iStaircase = 1:length(stimulus.eccentricity)
    % print message of what we are doing
    disp(sprintf('(liaoTemplate) Initializing staircase for eccentricity: %0.2f',stimulus.eccentricity(iStaircase)));
    % init staircase
    stimulus.s(iStaircase,1) = doStaircase('init','upDown','nup=1','ndown=2','initialStepsize',initialStepsize,'nTrials',nTrials,'initialThreshold',initialThreshold,'subplotCols',length(stimulus.eccentricity),'subplotNum',iStaircase,'dispFig',dispStaircaseFig,'subplotName',sprintf('Eccentricity: %0.2f',stimulus.eccentricity(iStaircase)),'minThreshold',0,'stepRule=pest','maxStepsize',0.5,'minStepsize',0.005);
  end
else
  disp(sprintf('(liaoTemplate) Found stimfile'))
  % init using threshold from last stimfile
  for iStaircase = 1:length(stimulus.eccentricity)
    % print message of what we are doing
    threshold = doStaircase('threshold',stimfile.stimulus.s(iStaircase,:));
    disp(sprintf('(liaoTemplate) Initializing staircase for eccentricity: %0.2f from last stimfile with threshold: %0.2f',stimulus.eccentricity(iStaircase),threshold.threshold));
    % init staircase
    stimulus.s(iStaircase,1) = doStaircase('init',stimfile.stimulus.s(iStaircase,end));
  end
end

% set that the staircases are not yet done
stimulus.staircaseCompleted = zeros(1,length(stimulus.eccentricity));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the grating stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initGratingStimulus(stimulus,myscreen)

% these are the reserved colors, if you need them later
% you can display them by setting your color to the appropriate
% index in stimulus.colors.reservedColor e.g. to get the
% second color, in this case white, you would do
% mglClearScreen(stimulus.colors.reservedColor(2));
stimulus.colors.reservedColors = [1 1 1; 0 1 0; 1 0 0];

% grating parameters
stimulus.grating.n = 4;
stimulus.grating.sf = 2;
stimulus.grating.tf = 2;
stimulus.grating.width = stimulus.width;
stimulus.grating.height = stimulus.width;
stimulus.grating.windowType = 'gabor'; % should be gabor or thresh
stimulus.grating.sdx = stimulus.grating.width/7;
stimulus.grating.sdy = stimulus.grating.width/7;

% init the gratins
stimulus = initGratings(stimulus,myscreen);

% set background color
stimulus.backgroundColor = stimulus.colors.grayColor;

% and fixation color
stimulus.normalFixColor = stimulus.colors.reservedColor(1);
stimulus.correctFixColor = stimulus.colors.reservedColor(2);
stimulus.incorrectFixColor = stimulus.colors.reservedColor(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     initGratings - from taskTemplateContrast10bit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initGratings(stimulus,myscreen)

% set maximum color index (for 24 bit color we have 8 bits per channel, so 255)
maxIndex = 255;

% get gamma table
if ~isfield(myscreen,'gammaTable')
  stimulus.linearizedGammaTable = mglGetGammaTable;
  disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
  disp(sprintf('(taskTemplateContrast10bit:initGratings) No gamma table found in myscreen. Contrast'));
  disp(sprintf('         displays like this should be run with a valid calibration made by moncalib'));
  disp(sprintf('         for this monitor.'));
  disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
end
stimulus.linearizedGammaTable = myscreen.initScreenGammaTable;

disppercent(-inf,'Creating grating textures');

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

% make the window through with the gratings will be displayed
gaussianWin = mglMakeGaussian(stimulus.grating.width,stimulus.grating.height,stimulus.grating.sdx,stimulus.grating.sdy);
if strcmp(stimulus.grating.windowType,'gabor')
  % a gaussian window
  win = maxIndex-maxIndex*gaussianWin;
else
  % a simple window
  win = maxIndex-maxIndex*(gaussianWin>exp(-1/2));
end
mask = ones(size(win,1),size(win,2),4)*stimulus.colors.midGratingIndex;
mask(:,:,4) = win;
stimulus.mask = mglCreateTexture(mask);

% make all the 1D gratings. We compute all possible contrast values given the
% range of indexes available to us. The 1st texture is gray the nth texture is full
% contrast for the current gamma setting
for iContrast = 0:stimulus.colors.nDisplayContrasts
  disppercent(iContrast/stimulus.colors.nDisplayContrasts);
  if myscreen.userHitEsc,mglClose;keyboard,end
  % make the grating
  thisGrating = round(iContrast*mglMakeGrating(stimulus.grating.width,nan,stimulus.grating.sf,0,0)+stimulus.colors.midGratingIndex);
  % create the texture
  stimulus.tex(iContrast+1) = mglCreateTexture(thisGrating);
end
disppercent(inf);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sets the gamma table so that we can have
% finest possible control over the stimulus contrast.
%
% stimulus.reservedColors should be set to the reserved colors (for cue colors, etc).
% maxContrast is the maximum contrast you want to be able to display.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setGammaTableForMaxContrast(maxContrast)

global stimulus;
% if you just want to show gray, that's ok, but to make the
% code work properly we act as if you want to display a range of contrasts
if maxContrast <= 0,maxContrast = 0.01;end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getContrastIndex    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function contrastIndex = getContrastIndex(desiredContrast,verbose)

if nargin < 2,verbose = 0;end

global stimulus;
if desiredContrast < 0, desiredContrast = 0;end

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

