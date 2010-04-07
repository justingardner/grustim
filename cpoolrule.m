% cpoolrule
%
%      usage: myscreen=cpoolrule(stimulus)
%         by: justin gardner
%       date: 04/15/06
%    purpose: orientation discrimination task
%
%
%
function myscreen = cpoolrule(varargin)

% check arguments
if ~any(nargin == [0 1])
  help cpoolrule
  return
end

isLoc = [];
getArgs(varargin,{'isLoc=0'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other screen parameters
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
deltaContrasts = [-0.2 -0.1 -0.05 0 0.05 0.1 0.2]/4;

%pedestalContrasts = [0.5];
%deltaContrasts = [-0.5 0.5];

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

stimulus.grating.windowType = 'gabor'; % should be gabor or thresh
stimulus.grating.sdx = stimulus.grating.width/7;
stimulus.grating.sdy = stimulus.grating.width/7;

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
  task{1}{2}.numBlocks = 10;

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
stimulus = initGratings(stimulus,myscreen,task);

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

% set maximum color index (for 24 bit color we have 8 bits per channel, so 255)
maxIndex = 255;

% get gamma table
if ~isfield(myscreen,'gammaTable')
  stimulus.linearizedGammaTable = mglGetGammaTable;
  disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
  disp(sprintf('(cpoolrule:initGratings) No gamma table found in myscreen. Contrast displays like this'));
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
  % set the maximum contrast we can display
  setGammaTableForMaxContrast(max([task.thistrial.pedestalContrast; task.thistrial.pedestalContrast+task.thistrial.deltaContrast]));
end

if any(task.thistrial.thisseg == [stimulus.int1 stimulus.int2])
  % get the contrast for the stimuli in this interval
  stimulusContrast = task.thistrial.pedestalContrast+task.thistrial.deltaContrast.*(task.thistrial.interval==task.thistrial.thisseg);

 
  % get the overall contrast
  if task.thistrial.thisseg == stimulus.int1
    stimulus.contrast1= sum(stimulusContrast);
  else
    stimulus.contrast2= sum(stimulusContrast);
  end

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
  
    mglLines2(0,0,stimulus.grating.cueLines(i,1),stimulus.grating.cueLines(i,2),1,stimulus.colors.reservedColor(2));
  end
  
  
end

% draw reference points
mglLines2(stimulus.grating.refLines.x1,stimulus.grating.refLines.y1,stimulus.grating.refLines.x2,stimulus.grating.refLines.y2,1,stimulus.colors.reservedColor(2));


mglFixationCross(0.5,1,stimulus.fixColor);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to set stimulus parameters at
% the beginning of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = trialResponseCallback(task,myscreen)

global stimulus;

fprintf(1,sprintf('difference: %f ',stimulus.contrast1-stimulus.contrast2));
if (stimulus.contrast1 > stimulus.contrast2) && (task.thistrial.whichButton == 1)
  disp(sprintf('correct'));
  stimulus.fixColor = stimulus.colors.reservedColor(3);
elseif (stimulus.contrast1 < stimulus.contrast2) && (task.thistrial.whichButton == 2)
  disp(sprintf('correct'));
  stimulus.fixColor = stimulus.colors.reservedColor(3);
elseif (stimulus.contrast1 > stimulus.contrast2) && (task.thistrial.whichButton == 2)
  disp(sprintf('incorrect'));
  stimulus.fixColor = stimulus.colors.reservedColor(4);
elseif (stimulus.contrast1 < stimulus.contrast2) && (task.thistrial.whichButton == 1)
  disp(sprintf('incorrect'));
  stimulus.fixColor = stimulus.colors.reservedColor(4);
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getContrastIndex    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function contrastIndex = getContrastIndex(desiredContrast,verbose)

if nargin < 2,verbose = 0;end

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

