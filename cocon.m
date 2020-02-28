% cocon.m
%
%        $Id$
%      usage: cocon
%         by: justin gardner
%       date: 02/25/2020
%  copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%    purpose: 2AFC motion direction paradigm for testing COherence CONtext
%
function myscreen = cocon(varargin)

% check arguments
getArgs(varargin,'stimulusType=dots');

% initalize the screen
myscreen = initScreen;

% init the stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);
stimulus.stimulusType = lower(stimulusType);

% stimulus parameters
stimulus.width = 10;
stimulus.eccentricity = [7.5];
stimulus.coherence = [0.25 0.5];

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
  disp(sprintf('(cocon) Unknown stimulus type: %s',stimulus.stimulusType));
  endScreen;
  return
end

% init the confidence display parameters (see function definition below for definition of parameters)
stimulus = initConfidence(stimulus,0,0,3,8,2,[1 1 1],[0.3 0.3 0.3]);
stimulus.feedback.segnum = 5;

% init the staircases
initialThreshold = 30;
initialStepsize = 5;
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
task{1}.segmin = [0.5 1 inf inf 0.5];
task{1}.segmax = [0.5 1 inf inf 0.5];
task{1}.getResponse = [0 1 1];
task{1}.parameter.eccentricity = stimulus.eccentricity;
task{1}.parameter.coherence = stimulus.coherence;
task{1}.randVars.uniform.whichSide = [1 2];
task{1}.randVars.calculated.confidnece = nan;
task{1}.randVars.calculated.correctIncorrect = nan;
task{1}.randVars.calculated.direction = nan;
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

  % get which staircase, note that we have kxn staircases - k coherences and n staircases
  stimulus.thisStaircase = find(stimulus.coherence==task.thistrial.coherence);
  stimulus.thisStaircase(2) = stimulus.staircaseCompleted(stimulus.thisStaircase)+1;

  % get the delta direction to test
  [delta stimulus.s(stimulus.thisStaircase(1),stimulus.thisStaircase(2))] = doStaircase('testValue',stimulus.s(stimulus.thisStaircase(1),stimulus.thisStaircase(2)));

  % make sure delta does not go below zero
  stimulus.delta = max(0,delta);
  
  % set  direction
  task.thistrial.direction = round(rand*360);
  if task.thistrial.whichSide == 1
    leftDirection = mod(task.thistrial.direction+stimulus.delta,360);
    rightDirection = task.thistrial.direction;
  else
    leftDirection = task.thistrial.direction;
    rightDirection = mod(task.thistrial.direction+stimulus.delta,360);
  end

  % display what is going on
  disp(sprintf('%i: Coherence: %0.1f Side: %i delta: %0.2f (Direction %0.1f vs %0.1f Eccentricity: %0.1f)',task.trialnum,task.thistrial.coherence,task.thistrial.whichSide,stimulus.delta,leftDirection,rightDirection,task.thistrial.eccentricity));
  
  % grating
  if stimulus.isGrating
    disp(sprintf('NOT YET IMPLEMENTED'));
    mglClose;
    keyboard
    % set the gamma table
    setGammaTableForMaxContrast(stimulus.pedestalContrast+stimulus.delta);
    % set the grating indexes
    stimulus.leftContrastIndex = getContrastIndex(leftContrast);
    stimulus.rightContrastIndex = getContrastIndex(rightContrast);
  % dots
  else
    % set the coherence
    stimulus.dotsLeft = stimulus.dotsLeft.setCoherence(stimulus.dotsLeft,task.thistrial.coherence);
    stimulus.dotsRight = stimulus.dotsRight.setCoherence(stimulus.dotsRight,task.thistrial.coherence);

    % set the position
    stimulus.dotsLeft = stimulus.dotsLeft.setCenter(stimulus.dotsLeft,-task.thistrial.eccentricity,0);
    stimulus.dotsRight = stimulus.dotsRight.setCenter(stimulus.dotsRight,task.thistrial.eccentricity,0);

    % set the direction
    stimulus.dotsLeft = stimulus.dotsLeft.setDir(stimulus.dotsLeft,leftDirection);
    stimulus.dotsRight = stimulus.dotsRight.setDir(stimulus.dotsRight,rightDirection);
  end

  % set the fixation color
  stimulus.fixColor = stimulus.normalFixColor;
elseif task.thistrial.thisseg == stimulus.confidence.segnum
  % set starting confidnece
  task.thistrial.confidence = 0.5;
  scrollEvents = mglListener('getAllScrollEvents');
  mglListener('getAllMouseEvents');
elseif task.thistrial.thisseg == stimulus.feedback.segnum
  if isequal(task.thistrial.correctIncorrect,1)
    stimulus.fixColor = stimulus.correctFixColor;
  elseif isequal(task.thistrial.correctIncorrect,0)
    stimulus.fixColor = stimulus.incorrectFixColor;
  end
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
elseif task.thistrial.thisseg == stimulus.confidence.segnum
  % set the confidence
  [task.thistrial.confidence confidenceDone] = setConfidence(task.thistrial.confidence, stimulus);
  if confidenceDone
    task = jumpSegment(task);
    disp(sprintf('(cocon) Confidence: %0.2f',task.thistrial.confidence));
  end
end

% draw fixation cross
mglFixationCross(1,1,stimulus.fixColor);

%%%%%%%%%%%%%%%%%%%%%%%%
%    initConfidence    %
%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initConfidence(stimulus,centerX,centerY,width,height,lineSize,lineColor,fillColor)

% set the segment in which the confidence judgement happens
stimulus.confidence.segnum = 4;
  
% set the dimensions of the confidence display
stimulus.confidence.width = width;
stimulus.confidence.height = height;
stimulus.confidence.centerX = centerX;
stimulus.confidence.centerY = centerY;

% set line size and color
stimulus.confidence.outlineSize = lineSize;
stimulus.confidence.outlineColor = lineColor;
stimulus.confidence.fillColor = fillColor;

% now compute the x and y of the outline
xL = stimulus.confidence.centerX-stimulus.confidence.width/2;
xR = stimulus.confidence.centerX+stimulus.confidence.width/2;
yB = stimulus.confidence.centerY-stimulus.confidence.height/2;
yT = stimulus.confidence.centerY+stimulus.confidence.height/2;

% dimensions of rectangle
stimulus.confidence.X0 = [xL xR xR xL];
stimulus.confidence.X1 = [xR xR xL xL];
stimulus.confidence.Y0 = [yB yB yT yT];
stimulus.confidence.Y1 = [yB yT yT yB];

% dimensions of fill
stimulus.confidence.fillX = [xL xR xR xL];
stimulus.confidence.fillY = [yB yB yT yT];

% values that need to be changed to reflect confidence level
stimulus.confidence.fillTop = [0 0 1 1];

% turn off scrolling
mglListener('eatscroll');
% read all pending scroll events
scrollEvents = mglListener('getAllScrollEvents');

%%%%%%%%%%%%%%%%%%%%%%%%
%    drawConfidence    %
%%%%%%%%%%%%%%%%%%%%%%%%
function drawConfidence(confidenceLevel, stimulus)

% draw confidence level as a filled bar.

% draw filled inside, compute top coordinate
fillY = stimulus.confidence.fillY;
fillY(find(stimulus.confidence.fillTop)) = stimulus.confidence.centerY+(-0.5+confidenceLevel)*stimulus.confidence.height;
% now draw as a filled polygon
mglPolygon(stimulus.confidence.fillX,fillY,stimulus.confidence.fillColor);

% draw outline
mglLines2(stimulus.confidence.X0,stimulus.confidence.Y0,stimulus.confidence.X1,stimulus.confidence.Y1,stimulus.confidence.outlineSize,stimulus.confidence.outlineColor);

%%%%%%%%%%%%%%%%%%%%%%%
%    setConfidence    %
%%%%%%%%%%%%%%%%%%%%%%%
function [confidence confidenceDone] = setConfidence(confidence, stimulus)

% get scroll
scrollEvents = mglListener('getAllScrollEvents');
if ~isempty(scrollEvents)
  % get the sum of the vertical and horizontal scrolls
  verticalScroll = -sum(scrollEvents.scrollVertical);
  horizontalScroll = sum(scrollEvents.scrollHorizontal);
  % set confidence
  confidence = confidence+verticalScroll/100;
else
  horizontalScroll = 0;
end

% check bounds
if confidence > 1,confidence = 1;end
if confidence < 0,confidence = 0;end
  
% draw the confidence
drawConfidence(confidence,stimulus);

% if mouse button down (or horizontal scroll is non-zero) then we are done setting confidence
mouse = mglGetMouse;

if ~isequal(mouse.buttons,0) % || ~isequal(horizontalScroll,0)
  confidenceDone = 1;
else
  confidenceDone = 0;
end

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
    % save it in task
    task.thistrial.correctIncorrect = 1;
    % change fixation color
    stimulus.fixColor = stimulus.responseFixColor;
    % and update staircase
    stimulus.s(stimulus.thisStaircase(1),stimulus.thisStaircase(2)) = doStaircase('update',stimulus.s(stimulus.thisStaircase(1),stimulus.thisStaircase(2)),1,stimulus.delta);
  else
    % report answer
    disp(sprintf(' ++ Incorrect ++. Reaction time: %0.2f',task.thistrial.reactionTime));
    % save it in task
    task.thistrial.correctIncorrect = 0;
    % change fixation color
    stimulus.fixColor = stimulus.responseFixColor;
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
  % jump to confidence rating part of text  
  task = jumpSegment(task);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dots stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initDotsStimulus(stimulus,myscreen)

% init the dot patchs
stimulus.dotsLeft = dotsInit('framesPerSecond',myscreen.framesPerSecond,'dir=0','width',stimulus.width,'speed=2.5');
stimulus.dotsRight = dotsInit('framesPerSecond',myscreen.framesPerSecond,'dir=180','width',stimulus.width,'speed=2.5');

% set background color
stimulus.backgroundColor = 0.5;

% and fixation color
stimulus.normalFixColor = [1 1 1];
stimulus.correctFixColor = [0 1 0];
stimulus.incorrectFixColor = [1 0 0];
stimulus.responseFixColor = [1 1 1];
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

% check that the last stimfile had the same coherences
if ~isempty(stimfile)
  if ~isfield(stimfile,'stimulus') || ~isfield(stimfile.stimulus,'coherence') || ~isequal(stimfile.stimulus.coherence,stimulus.coherence)
    % dump this stimfile, since it does not match current eccentricity
    disp(sprintf('(cocon) Found stimfile, but does not have coherence match. Ignoring'));
    stimfile = [];
  end
end

% if no stimfile found
if isempty(stimfile)
  % then initialize
  for iStaircase = 1:length(stimulus.coherence)
    % print message of what we are doing
    disp(sprintf('(cocon) Initializing staircase for coherence: %0.2f',stimulus.coherence(iStaircase)));
    % init staircase
    stimulus.s(iStaircase,1) = doStaircase('init','upDown','nup=1','ndown=2','initialStepsize',initialStepsize,'nTrials',nTrials,'initialThreshold',initialThreshold,'subplotCols',length(stimulus.coherence),'subplotNum',iStaircase,'dispFig',dispStaircaseFig,'subplotName',sprintf('Coherence: %0.2f',stimulus.coherence(iStaircase)),'minThreshold',0,'stepRule=pest','maxStepsize',0.5,'minStepsize',0.005);
  end
else
  disp(sprintf('(cocon) Found stimfile'))
  % init using threshold from last stimfile
  for iStaircase = 1:length(stimulus.coherence)
    % print message of what we are doing
    threshold = doStaircase('threshold',stimfile.stimulus.s(iStaircase,:));
    disp(sprintf('(cocon) Initializing staircase for coherence: %0.2f from last stimfile with threshold: %0.2f',stimulus.coherence(iStaircase),threshold.threshold));
    % init staircase
    stimulus.s(iStaircase,1) = doStaircase('init',stimfile.stimulus.s(iStaircase,end));
  end
end

% set that the staircases are not yet done
stimulus.staircaseCompleted = zeros(1,length(stimulus.coherence));

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

