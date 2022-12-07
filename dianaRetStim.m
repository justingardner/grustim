%    dianaTestTriggers.m
%
%       
%       usage: dianaTestTriggers
%          by: jw + jlg
%        date: October 2022
%       
%       Test stimulus for DIANA imaging.
%      
function myscreen = dianaTestTriggers(varargin)

% check arguments
getArgs(varargin);

% initilaize the screen
%myscreen = initScreen('debug');
myscreen = initScreen('diana');

% init the stimulus
global stimulus;

% onset after cardiac pulse to stimuls presentation in seconds.
stimulus.retinotopyOnsetDelay = 0.250;
% each hemifield will be stimulated for this long, so overall
% the stimulus will be on for 2x this value.
stimulus.retinotopyStimLen = 0.250;

% set stimulus vaules, these are in ms, but will be
% translated into number of video frames and the
% actual values that you can expect will be displayed.
stimulus.onsetDelay = 50;
% this is in Hz. The stimulus is phase-reversing checkerboard
stimulus.tf = 2;
% number of cycles to display for
% size in degrees
stimulus.width = 78;
stimulus.height = 64;
% spatial frequency in cycles/deg
stimulus.sf = 0.5;
stimulus.orientation = 90;
stimulus.phase = 0;
stimulus.contrast = 0.5;
% number of cycles to run stimulus for
stimulus.numCycles = 1;
% use fix task (or just draw a fixation cross if not)
stimulus.useFixTask = 1;
stimulus.fixWidth = 1;

% CODE FROM RETINOTOPY (NEEDS TO BE CLEANED UP)
stimulus.imageWidth = myscreen.imageWidth;
stimulus.imageHeight = myscreen.imageHeight;
if ieNotDefined('fixedRandom') stimulus.fixedRandom = 0;else stimulus.fixedRandom = fixedRandom;end
if ieNotDefined('blanks'),stimulus.blanks = false;end
if ieNotDefined('direction'),stimulus.direction = 1;end
if ieNotDefined('dutyCycle'),stimulus.dutyCycle = 0.5;end
if ieNotDefined('stepsPerCycle'),stimulus.stepsPerCycle = 24;end
if ieNotDefined('stimulusPeriod'),stimulus.stimulusPeriod = 24;end
if ieNotDefined('numCycles'),stimulus.numCycles = 10;end
if ieNotDefined('doEyeCalib'),stimulus.doEyeCalib = -1;end
if ieNotDefined('initialHalfCycle') stimulus.initialHalfCycle = 1;end
if ieNotDefined('easyFixTask'),stimulus.easyFixTask = 1;end
if ieNotDefined('dispText'),stimulus.dispText = '';end
if ieNotDefined('barWidth'),stimulus.barWidth = 3;end
if ieNotDefined('elementAngle'),stimulus.elementAngle = 'parallel';end
if ieNotDefined('barSweepExtent'),stimulus.barSweepExtent = 'min';end
if ieNotDefined('elementSize'),stimulus. elementSize = [];end
if ieNotDefined('barStepsMatchElementSize'),stimulus.barStepsMatchElementSize=true;end
if ieNotDefined('synchToVolEachCycle'),stimulus.synchToVolEachCycle=false;end
% settings that are used to adjust the position on the screen
% the stimuli are shown in - for cases when the subject can
% not see the whole screen
if ieNotDefined('xOffset'),stimulus.xOffset = 0;end
if ieNotDefined('yOffset'),stimulus.yOffset = 0;end
% settings for dots
if ieNotDefined('dotsDensity'),stimulus.dotsDensity = 25;end
if ieNotDefined('dotsSpeed'),stimulus.dotsSpeed = 3;end
if ieNotDefined('dotsSpace'),stimulus.dotsSpace = 0.5;end % space between three patches
if ieNotDefined('dotsSize'),stimulus.dotsSize = 4;end
% parameters for barsTask
if ieNotDefined('startCoherence'),stmulus.startCoherence = 0.5;end
if ieNotDefined('stepCoherence'),stimulus.stepCoherence = 0.05;end
if ieNotDefined('minCoherence'),stimulus.minCoherence = 0.025;end
if ieNotDefined('maxCoherence'),stimulus.maxCoherence = 0.975;end
if ieNotDefined('minStepCoherence'),stimulus.minStepCoherence = 0.005;end
if ieNotDefined('fixFeedback'),stimulus.fixFeedback = true;end % give task feedback on fixation cross
% set the parameters of the stimulus
% whether to display wedges or rings or bars
% 1 for wedges 2 for rings 3 for bars
stimulus.stimulusType = 1;
% min/max radius is the size of the stimulus
stimulus.minRadius = 0.5;
stimulus.maxRadius = min(stimulus.imageWidth/2,stimulus.imageHeight/2);
stimulus.barContrast = 1;
%stimulus.maxRadius = stimulus.imageWidth;
% direction of stimulus
%stimulus.direction = direction;
% the duty cycle of the stimulus.
% angle size is the size in degrees of the
% elements of the wedge that slide against each other
stimulus.elementAngleSize = 5;
% radius is the radial length of these elements
stimulus.elementRadiusSize = 2;
% radial speed of elements moving each other
stimulus.elementRadialVelocity = 7.5;
% element size parameters for bars (non-radial pattern)
stimulus.elementWidth = 3;
stimulus.elementHeight = 3;
stimulus.elementVelocity = 6;
stimulus = initRetinotopyStimulus(stimulus,myscreen);
stimulus.cycleTime = mglGetSecs;
% CODE FROM RETINOTOPY END (NEEDS TO BE CLEANED UP)

% clear the screen to gray
mglClearScreen(0.5);mglFlush;

% wait for backtick to start experiment (probably not necessary as we
% wait at the beginning of each segment anyway).
task{1}.waitForBacktick = 0;

% initStimulus
myscreen = initStimulus('stimulus',myscreen);
myInitStimulus(myscreen);

% set task parameters
task{1}.segmin = [0.1 stimulus.retinotopyOnsetDelay stimulus.retinotopyStimLen stimulus.retinotopyStimLen 0.1];
task{1}.segmax = [0.1 stimulus.retinotopyOnsetDelay stimulus.retinotopyStimLen stimulus.retinotopyStimLen 0.1];
task{1}.getResponse = [0 0 0 0 0];

% set number of trials to infinite (to stop stimulus hit ESC)
task{1}.numTrials = inf;

% synch to vol (note this will sync *after* the first segment is done in time, so that
% the stimulus presented in the 2nd segment will start just after the volume acquisition
% trigger comes.
task{1}.synchToVol = [1 0 0 0 0];

% setup fixation task
if stimulus.useFixTask
  global fixStimulus
  fixStimulus.fixWidth = stimulus.fixWidth
  [fixTask myscreen] = fixStairInitTask(myscreen);
end

% initialize the task
for phaseNum = 1:length(task)
  [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task)) && ~myscreen.userHitEsc
  % update the task
  [task myscreen phaseNum] = updateTask(task,myscreen,phaseNum);
  if stimulus.useFixTask
    % draw fixation task
    [fixTask myscreen] = updateTask(fixTask,myscreen,1);
  end
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

% print out what we are doing and clear screen
if task.thistrial.thisseg == 1
  % calculate how long it took to run last two segment
  %stimulus.endTime = mglGetSecs;
  if ~isnan(stimulus.startTime)
    disp(sprintf('(dianaTestTriggers) Stim segment lasted %0.1f ms: Ready to start trial %i',1000*(stimulus.endTime-stimulus.startTime),task.trialnum));
  else    
    disp(sprintf('(dianaTestTriggers) Ready to start trial %i',task.trialnum));
  end
  mglClearScreen(0.5);
elseif task.thistrial.thisseg == 3
  % keep starting tick
  stimulus.startTick = myscreen.tick;
  stimulus.startTime = mglGetSecs;
  stimulus.currentMask = 7;
elseif task.thistrial.thisseg == 4
  stimulus.currentMask = 19;
elseif task.thistrial.thisseg == 5
  stimulus.endTime = mglGetSecs;
  mglClearScreen(0.5);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%resp%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen, checkboardBltTex)

global stimulus

% for the secound segment show a stimulus
if any(task.thistrial.thisseg == [3 4])
  % flashing box from white to black

  %if isodd(floor((myscreen.tick-stimulus.startTick)/stimulus.framesPerHalfCycle))
  %if floor((myscreen.tick-stimulus.startTick)/stimulus.framesPerHalfCycle) == 0 % edit here - set 0 instead of odd; otherwise flashes back
      
  %  mglBltTexture(stimulus.plaid1,[0 0]);
  %else
  %  mglBltTexture(stimulus.plaid2,[0 0]);
  %end
  % draw retinotopy stimulus  
  stimulus = updateRetinotopyStimulus(stimulus,myscreen);

 
 % draw fixation cross, if not using task
  if ~stimulus.useFixTask
    mglFixationCross(stimulus.fixWidth,2,[0 1 1]);
  end
else
  % otherwise just clear the screen to black
  mglClearScreen(0.5);
  % draw fixation cross, if not using task
  if ~stimulus.useFixTask
    mglFixationCross(stimulus.fixWidth,2,[0 1 1]);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = getResponseCallback(task,myscreen)

%%%%%%%%%%%%%%%%%%%%%%%%
%    myInitStimulus    %
%%%%%%%%%%%%%%%%%%%%%%%%
function myInitStimulus(myscreen)

global stimulus

% calculate number of frames for onestDelay
stimulus.onsetDelayFrames = round(stimulus.onsetDelay*myscreen.framesPerSecond/1000);
stimulus.onsetDelayActual = 1000*stimulus.onsetDelayFrames*(1/myscreen.framesPerSecond);
disp(sprintf('(dianaTestTriggers) onsetDelay for %i frames: %0.1f ms (desired=%0.1f ms)',stimulus.onsetDelayFrames,stimulus.onsetDelayActual,stimulus.onsetDelay));

% gratings
stimulus.g1 = stimulus.contrast*mglMakeGrating(stimulus.width,stimulus.height,stimulus.sf,stimulus.orientation,stimulus.phase);
stimulus.g2 = stimulus.contrast*mglMakeGrating(stimulus.width,stimulus.height,stimulus.sf,stimulus.orientation+90,stimulus.phase);

% plaid
stimulus.plaid = stimulus.g1 + stimulus.g2;

% make into checkerboard
stimulus.plaid(stimulus.plaid>0) = 1;
stimulus.plaid(stimulus.plaid<0) = -1;

stimulus.plaid1 = mglCreateTexture(stimulus.plaid);
stimulus.plaid2 = mglCreateTexture(-stimulus.plaid);

% compute frames per each cycle
stimulus.framesPerHalfCycle = round(((1/stimulus.tf) * myscreen.framesPerSecond)/2);
stimulus.cycleLenActual = 1000*(stimulus.framesPerHalfCycle*2)*(1/myscreen.framesPerSecond);
stimulus.tfActual = 1/(stimulus.cycleLenActual/1000);
disp(sprintf('(dianaTestTriggers) tf=%0.1f (desired: %0.1f) cycleLen: %i frames (%0.1f ms, desired: %0.1f ms)',stimulus.tfActual,stimulus.tf,stimulus.framesPerHalfCycle*2,stimulus.cycleLenActual,1000/stimulus.tf));

% compute length of time stimulus will be on for
stimulus.stimLen = stimulus.cycleLenActual*stimulus.numCycles/1000;
disp(sprintf('(dianaTestTriggers) Running %i cycles for %i frames (%0.1f ms)',stimulus.numCycles,stimulus.stimLen*myscreen.framesPerSecond,stimulus.stimLen));

% Display total length of stimulus
disp(sprintf('(dianaTestTriggers) Total stimulus time (including onset delay): %0.1f',stimulus.onsetDelayActual+stimulus.cycleLenActual*stimulus.numCycles));

% set initial startTime
stimulus.startTime = nan;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the retinotopy stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initRetinotopyStimulus(stimulus,myscreen)

% round to nearest quarter of a degree, this reduces
% some edge effects
stimulus.maxRadius = floor(stimulus.maxRadius/.25)*.25;
disp(sprintf('Stimulus radius = [%0.2f %0.2f] degrees',stimulus.minRadius,stimulus.maxRadius));
% calculate some parameters
% size of wedges
stimulus.wedgeAngle = 360*stimulus.dutyCycle;
% how much to step the wedge angle by
stimulus.wedgeStepSize = 360/stimulus.stepsPerCycle;
% Based on the duty cycle calculate the ringSize.
% Note that this is not just maxRadius-minRadius * dutyCycle
% because we start the rings off the inner edge and end
% them off the outer edge (that is the screen will be blank for
% a time, rathter than showing a ring). We also have to
% compensate in our dutyCycle for this, since the effective
% duty cycle is reduced by the time that we are offscreen
% that is for two periods
dutyCycle = stimulus.dutyCycle*stimulus.stepsPerCycle/(stimulus.stepsPerCycle-2);
stimulus.ringSize = (dutyCycle*(stimulus.maxRadius-stimulus.minRadius))/(1-dutyCycle);
% get the radii for the rings
minRadius = stimulus.minRadius-stimulus.ringSize;
maxRadius = stimulus.maxRadius+stimulus.ringSize;
% now we have the inner and outer ring radius that will be used
% add a little fudge factor so that we don't get any rings
% with only a small bit showing
epsilon = 0.1;
stimulus.ringRadiusMin = max(0,minRadius:(epsilon+stimulus.maxRadius-minRadius)/(stimulus.stepsPerCycle-1):stimulus.maxRadius+epsilon);
stimulus.ringRadiusMax = stimulus.minRadius:(maxRadius-stimulus.minRadius)/(stimulus.stepsPerCycle-1):maxRadius;

% set some parameters for bars
stimulus.barHeight = stimulus.imageWidth*1.5;
stimulus.barMaskWidth = stimulus.imageWidth*1.5;

% dots stimulus
if any(stimulus.stimulusType == [4 5])
  % parameters of dots
  stimulus.dots.dir = pi/2;
  stimulus.dots.initialCoherence = 1;

  % make dots for each of the three segments
  stimulus.dots.length = (stimulus.imageHeight-stimulus.dots.space*4)/3;
  stimulus.dots.upper = initDots(myscreen,0,stimulus.dots.length+stimulus.dots.space,stimulus.barWidth,stimulus.dots.length,stimulus.dots.dir,stimulus.dots.speed,stimulus.dots.size,stimulus.dots.density,stimulus.dots.initialCoherence);
  stimulus.dots.middle = initDots(myscreen,0,0,stimulus.barWidth,stimulus.dots.length,stimulus.dots.dir,stimulus.dots.speed,stimulus.dots.size,stimulus.dots.density,stimulus.dots.initialCoherence);
  stimulus.dots.lower = initDots(myscreen,0,-stimulus.dots.length-stimulus.dots.space,stimulus.barWidth,stimulus.dots.length,stimulus.dots.dir,stimulus.dots.speed,stimulus.dots.size,stimulus.dots.density,stimulus.dots.initialCoherence);

  % see if there was a previous staircase
  s = getLastStimfile(myscreen);
  if ~isempty(s)
    if isfield(s,'stimulus') && isfield(s.stimulus,'stair')
      % check if we are doing barsTaskEasy
      if ~stimulus.barsTaskEasy
	% and that the old stimfile was not doing barsTaskEasy
	if (~isfield(s.stimulus,'barsTaskEasy') || ~s.stimulus.barsTaskEasy)
	  % set starting thershold to staircase value
	  stimulus.barsTask.startCoherence = s.stimulus.stair.threshold;
	  % display what we are doing
	  disp(sprintf('(mglRetinotopy) Setting barTask starting thershold based on last stimfile to: %f',stimulus.barsTask.startCoherence));
	end
      end
    end
  end
  % initialize  staircase
  stimulus.stair = upDownStaircase(1,2,stimulus.barsTask.startCoherence,[stimulus.barsTask.stepCoherence stimulus.barsTask.minStepCoherence stimulus.barsTask.stepCoherence],'pest');
  stimulus.stair.minThreshold = stimulus.barsTask.minCoherence;
  stimulus.stair.maxThreshold = stimulus.barsTask.maxCoherence;
end

  % we only need to recompute the mglQuad points of the elements if something has
  % changed in the stimulus. This is for the radial element pattern
  if ~isfield(stimulus,'last') || ~isfield(stimulus,'x') || ...
	(stimulus.elementAngleSize ~= stimulus.last.elementAngleSize) || ...
	(stimulus.elementRadiusSize ~= stimulus.last.elementRadiusSize) || ...
	(stimulus.elementRadialVelocity ~= stimulus.last.elementRadialVelocity) || ...
	(stimulus.maxRadius ~= stimulus.last.maxRadius) || ...
	(stimulus.minRadius ~= stimulus.last.minRadius)
    % all the angles that the elements will be made up of
    allAngles = (0:stimulus.elementAngleSize:(360-stimulus.elementAngleSize));
    % all the phases. The phase refers to the radial position of the
    % black and white pattern (the pattern that is seen as moving
    % in the stimulus). There are two sets here since the wedges slide
    % against each other. That is every other sector will go in a 
    % different direction. 
    allPhases1 = 0:(stimulus.elementRadialVelocity/myscreen.framesPerSecond):(stimulus.elementRadiusSize*2);
    allPhases2 = fliplr(allPhases1);
    disppercent(-inf,'(mglRetinotopy) Calculating coordinates of elements in stimulus pattern');
    for phaseNum = 1:length(allPhases1)
      stimulus.x{phaseNum} = [];stimulus.y{phaseNum} = [];stimulus.c{phaseNum} = [];
      for angleNum = 1:length(allAngles)
	% get the angle
	angle = allAngles(angleNum);
	% choose which phase we are going to be
	if isodd(angleNum)
	  thisMinRadius = stimulus.minRadius-allPhases1(phaseNum);
	else
	  thisMinRadius = stimulus.minRadius-allPhases2(phaseNum);
	end
	% all the radiuses
	allRadius = thisMinRadius:stimulus.elementRadiusSize:stimulus.maxRadius;
	% now create all the quads for this wedge
	for radiusNum = 1:length(allRadius)
	  radius = allRadius(radiusNum);
	  if (radius+stimulus.elementRadiusSize) >= stimulus.minRadius
	    radius1 = max(radius,stimulus.minRadius);
	    radius2 = min(radius+stimulus.elementRadiusSize,stimulus.maxRadius);
	    % calculate in polar angle coordinates the corners of this quad
	    r = [radius1 radius1 radius2 radius2];
	    a = [angle angle+stimulus.elementAngleSize angle+stimulus.elementAngleSize angle];
	    % convert into rectilinear coordinates and save in array
	    stimulus.x{phaseNum}(:,end+1) = r.*cos(d2r(a));
	    stimulus.y{phaseNum}(:,end+1) = r.*sin(d2r(a));
	    % also calculate what color we ant
	    stimulus.c{phaseNum}(:,end+1) = [1 1 1]*(isodd(radiusNum+isodd(angleNum)));
	  end
	end
      end
      disppercent(phaseNum/length(allPhases1));
    end
    disppercent(inf);
    stimulus.n = length(allPhases1);
    stimulus.phaseNum = 1;
  else
    disp(sprintf('(mglRetinotopy) Using precomputed stimulus pattern'));
  end

  % we only need to recompute the mglQuad points of the elements if something has
  % changed in the stimulus. This is for the rectilinear elements (i.e. for the bars)
  if ~isfield(stimulus,'last') || ~isfield(stimulus,'x') || ...
	(stimulus.elementWidth ~= stimulus.last.elementWidth) || ...
	(stimulus.elementHeight ~= stimulus.last.elementHeight) || ...
	(stimulus.elementVelocity ~= stimulus.last.elementVelocity) || ...
    (isfield(stimulus,'last') && (stimulus.barContrast ~= stimulus.last.barContrast))
    maxDim = ceil(max(stimulus.imageWidth,stimulus.imageHeight)/stimulus.elementWidth)*stimulus.elementWidth;
    minRect = -maxDim/2-2*stimulus.elementWidth;
    maxRect = maxDim/2;
    % all the angles that the elements will be made up of
    allY = minRect:stimulus.elementHeight:maxRect;
    % all the phases. The phase refers to the radial position of the
    % black and white pattern (the pattern that is seen as moving
    % in the stimulus). There are two sets here since the wedges slide
    % against each other. That is every other sector will go in a 
    % different direction. 
    allPhases1 = 0:(stimulus.elementVelocity/myscreen.framesPerSecond):(stimulus.elementWidth*2);
    allPhases2 = fliplr(allPhases1);
    disppercent(-inf,'(mglRetinotopy) Calculating coordinates of elements in stimulus pattern for rectilinear patterns (i.e. ones for bars)');
    for phaseNum = 1:length(allPhases1)
      stimulus.xRect{phaseNum} = [];stimulus.yRect{phaseNum} = [];stimulus.cRect{phaseNum} = [];
      for yNum = 1:length(allY)
	% get the y
	y = allY(yNum)+stimulus.elementHeight/2;
	% choose which phase we are going to be
	if isodd(yNum)
	  thisMinX = minRect-allPhases1(phaseNum);
	else
	  thisMinX = minRect-allPhases2(phaseNum);
	end
	% all the X
	allX = thisMinX:stimulus.elementWidth:maxRect;
	% now create all the quads for this wedge
	for xNum = 1:length(allX)
	  x = allX(xNum);
	  % calculate element
	  stimulus.xRect{phaseNum}(:,end+1) = [x x x+stimulus.elementWidth x+stimulus.elementWidth];
	  stimulus.yRect{phaseNum}(:,end+1) = [y y+stimulus.elementHeight y+stimulus.elementHeight y];
	  % also calculate what color we want
	  stimulus.cRect{phaseNum}(:,end+1) = ([1 1 1]*(isodd(xNum+isodd(yNum)))-0.5) * stimulus.barContrast + 0.5;
	end
      end
      disppercent(phaseNum/length(allPhases1));
    end
    disppercent(inf);
    stimulus.nRect = length(allPhases1);
    stimulus.phaseNumRect = 1;
    % later below when we are calculating the optimal placement of bars
    % over the elements we may wish to offset the location of the elements
    % we start off with no offset.
    stimulus.xRectOffset = 0;
  else
    disp(sprintf('(mglRetinotopy) Using precomputed stimulus pattern'));
  end

  % remember these parameters, so that we can know whether we
  % need to recompute
  stimulus.last.elementRadiusSize = stimulus.elementRadiusSize;
  stimulus.last.elementAngleSize = stimulus.elementAngleSize;
  stimulus.last.elementRadialVelocity = stimulus.elementRadialVelocity;
  stimulus.last.maxRadius = stimulus.maxRadius;
  stimulus.last.minRadius = stimulus.minRadius;
  stimulus.last.elementWidth = stimulus.elementWidth;
  stimulus.last.elementHeight = stimulus.elementHeight;
  stimulus.last.elementVelocity = stimulus.elementVelocity;
  stimulus.last.barContrast = stimulus.barContrast;

% new we calculate the masks that cover the stimulus so that we can
% have either rings or wedges, we start by making a set of wedge masks
angles = (0:stimulus.wedgeStepSize:(360-stimulus.wedgeStepSize))+90+stimulus.wedgeAngle/2;
% create masks for wedges
for angleNum = 1:length(angles)
  angle = angles(angleNum);
  % Make the mask from two overlapping semicircles, offset by wedgeAngle.
  % This wayt the mask polygons will be convex and triangle-strip-able.
  maskRadius = stimulus.maxRadius + 1;
  topRadians = d2r(angle:angle + 180);
  topX = maskRadius * cos(topRadians);
  topY = maskRadius * sin(topRadians);

  bottomAngle = angle - stimulus.wedgeAngle;
  bottomRadians = d2r(bottomAngle-180:bottomAngle);
  bottomX = maskRadius * cos(bottomRadians);
  bottomY = maskRadius * sin(bottomRadians);
  stimulus.maskWedgeX{angleNum} = {topX, bottomX};
  stimulus.maskWedgeY{angleNum} = {topY, bottomY};
end
stimulus.wedgeN = length(angles);

% now we will make the masks for the rings. We will
% make an inner and outer set of ring masks so that we
% can cover up everything but one ring of the stimulus
for radiusNum = 1:length(stimulus.ringRadiusMin)
  % create the inner mask
  stimulus.maskInnerX{radiusNum} = [];
  stimulus.maskInnerY{radiusNum} = [];
  % compute in radial coordinates
  r = [];a = [];
  for angle = 0:stimulus.elementAngleSize:360
    r(end+1) = stimulus.ringRadiusMin(radiusNum);
    a(end+1) = angle;
  end
  %r(end+1) = 0;a(end+1) = 0;
  % now convert to rectilinear
  stimulus.maskInnerX{radiusNum}(:,end+1) = r.*cos(d2r(a));
  stimulus.maskInnerY{radiusNum}(:,end+1) = r.*sin(d2r(a));
  % create the outer mask, this will be 
  % a set of quads that make a torus
  stimulus.maskOuterX{radiusNum} = [];
  stimulus.maskOuterY{radiusNum} = [];
  allAngles = 0:stimulus.elementAngleSize:360;
  for angleNum = 1:length(allAngles)
    angle = allAngles(angleNum);
    r = stimulus.ringRadiusMax(radiusNum);
    a = angle;
    r(end+1) = stimulus.maxRadius+1;
    a(end+1) = angle;
    r(end+1) = stimulus.maxRadius+1;
    a(end+1) = angle+stimulus.elementAngleSize;
    r(end+1) = stimulus.ringRadiusMax(radiusNum);
    a(end+1) = angle+stimulus.elementAngleSize;
    % convert to rectilinear
    stimulus.maskOuterX{radiusNum}(:,angleNum) = r.*cos(d2r(a));
    stimulus.maskOuterY{radiusNum}(:,angleNum) = r.*sin(d2r(a));
    stimulus.maskOuterC{radiusNum}(:,angleNum) = [0.5 0.5 0.5];
  end
end
stimulus.ringN = length(stimulus.ringRadiusMin);

% now make masks for bars
if any(stimulus.stimulusType == [3 4 5])
  % get x and y of bar center (note that this is before we apply the
  % rotation to coordinates, so we are making the coordinates as if
  % we are going to make horizontally sweeping bars - we will later
  % rotate the coordinates around the center to get the other sweep angels

  % start by figuring out over what the extent of visual angle
  % that we want to sweep the bars over
  if isempty(stimulus.barSweepExtent),stimulus.barSweepExtent = 'max';end
  if isstr(stimulus.barSweepExtent)
    if strcmp(stimulus.barSweepExtent,'min')
      sweepExtent = min(stimulus.imageWidth,stimulus.imageHeight);
    elseif strcmp(stimulus.barSweepExtent,'max')
      barDirVec = abs(stimulus.maskBarRotMatrix*[1 0]');
      sweepExtent = max(barDirVec(1)*stimulus.imageWidth,barDirVec(2)*stimulus.imageHeight);
    else
      disp(sprintf('(mglRetinotopy) Unknown Bar sweep extent; %s',stimulus.barSweepExtent));
      keyboard
    end
  else
    sweepExtent = stimulus.barSweepExtent;
  end
  % calculate the stepSize that we will have with this sweepExtent
  stepSize = sweepExtent/(stimulus.stepsPerCycle-1);
  % now round that stepSize to the size of half an element
  % so that the bars all show at the same part of the underlying
  % checkerboard
  if stimulus.barStepsMatchElementSize
    % set the granularity - that is, stepSize will need
    % to be an integer multiple of stepSizeGranularity.
    % If this is just elementWidth then the stepSize
    % will never have only part of element in it.
    stepSizeGranularity = stimulus.elementWidth/4;
    % keep original step size for comparison, to
    % tell user what we are doing
    originalStepSize = stepSize;
    stepSize = round(stepSize/stepSizeGranularity)*stepSizeGranularity;
    % make sure stepSize is greater than 0
    stepSize = max(stepSize,stepSizeGranularity);
    % if the stepSize is not evenly divisible into the elementWidth
    % then we add an offset to the start position, so that the stimulus
    % is always centered on the element size
    barAndElementNeedAlignment = stepSize/stimulus.elementWidth;
    barAndElementNeedAlignment = (barAndElementNeedAlignment - floor(barAndElementNeedAlignment)) ~= 0;
    % tell user what is going on if the stepSize is greater than 10% of desired
    if (stepSize >= 1.05*originalStepSize) ||  (stepSize <= 0.95*originalStepSize) 
      disp(sprintf('(mglRetinotopy) Bar step size has been set to %0.2f (from ideal %0.2f) which changes the coverage',stepSize,originalStepSize));
      disp(sprintf('                of the bars from the desired barSweepExtent of %0.2f to %0.2f ',-sweepExtent/2,sweepExtent/2));
      disp(sprintf('                (see barCenter setting below). This is done to match the underlying'));
      disp(sprintf('                pattern better. If you want to have the bars exactly cover the '));
      disp(sprintf('                barSweepExtent, set barStepsMatchElementSize to false or set the'));
      disp(sprintf('                elementSize (%0.2fx%0.2f) such that the bar step size is an integer multiple',stimulus.elementWidth,stimulus.elementHeight));
      disp(sprintf('                of that elementSize'));
    end
  else
    barAndElementNeedAlignment = false;
  end
  % now create the steps
  stimulus.barCenter = [];
  if isodd(stimulus.stepsPerCycle)
    % for odd number of steps cover the center of the screen
    stimulus.barCenter(:,1) = -stepSize*(stimulus.stepsPerCycle-1)/2:stepSize:stepSize*(stimulus.stepsPerCycle-1)/2;
  else
    % for even number of steps, will be symmetric around the center
    stimulus.barCenter(:,1) = -stepSize*stimulus.stepsPerCycle/2+stepSize/2:stepSize:stepSize*stimulus.stepsPerCycle/2-stepSize/2;
  end
  stimulus.barCenter(:,2) = 0;
  % check to see if we need to align bars or elements to get a match (we want the element pattern to be
  % centered on bars)
  if barAndElementNeedAlignment
    adjustmentSize = min(abs(stimulus.barCenter(:,1)));
    switch (lower(stimulus.barElementAlignmentCorrection))
     case {'movebars'}
      stimulus.barCenter(:,1) = stimulus.barCenter(:,1)-adjustmentSize;
      % remove any offset on the elements
      for i = 1:length(stimulus.xRect)
	stimulus.xRect{i} = stimulus.xRect{i}-stimulus.xRectOffset;
      end
      stimulus.xRectOffset = 0;
      disp(sprintf('(mglRetinotopy) Moving bar centers by %f to match to underlying elements. Check stimulus.barElementAlignmentCorrection if you want something different.',adjustmentSize));
     case {'moveelements'}
      for i = 1:length(stimulus.xRect)
	stimulus.xRect{i} = stimulus.xRect{i}-stimulus.xRectOffset+adjustmentSize;
      end
      stimulus.xRectOffset = adjustmentSize;
      disp(sprintf('(mglRetinotopy) Moving element centers by %f to match bar locations. Check stimulus.barElementAlignmentCorrection if you want something different.',adjustmentSize));
     case {'none'}
      % remove any offset on the elements
      for i = 1:length(stimulus.xRect)
	stimulus.xRect{i} = stimulus.xRect{i}-stimulus.xRectOffset;
      end
      stimulus.xRectOffset = 0;
      disp(sprintf('(mglRetinotopy) Detected misalignment of bar and elements center, but not doing any adjustment, if you want to adjust the bar location set stimulus.barElementAlignmentCorrection to moveBars. If you want to adjust the elements to match then set to moveElements'));
    end
  end

  % display the bar centers
  disp(sprintf('(mglRetinotopy) barCenter: %s',num2str(stimulus.barCenter(:,1)','%0.02f ')));

  % now make the left hand bar mask - remember we are making the inverse
  % of the bar so we need to mask out everything to the left and to the
  % right of the bar.
  stimulus.maskBarLeft(1,:) = [-stimulus.barWidth/2 -stimulus.barHeight/2]';
  stimulus.maskBarLeft(2,:) = [-stimulus.barWidth/2-stimulus.barMaskWidth -stimulus.barHeight/2]';
  stimulus.maskBarLeft(3,:) = [-stimulus.barWidth/2-stimulus.barMaskWidth stimulus.barHeight/2]';
  stimulus.maskBarLeft(4,:) = [-stimulus.barWidth/2 stimulus.barHeight/2]';
  % now make the right hand bar mask
  stimulus.maskBarRight(1,:) = [stimulus.barWidth/2 -stimulus.barHeight/2]';
  stimulus.maskBarRight(2,:) = [stimulus.barWidth/2+stimulus.barMaskWidth -stimulus.barHeight/2]';
  stimulus.maskBarRight(3,:) = [stimulus.barWidth/2+stimulus.barMaskWidth stimulus.barHeight/2]';
  stimulus.maskBarRight(4,:) = [stimulus.barWidth/2 stimulus.barHeight/2]';
end

% set the current mask that will be displayed
if stimulus.direction == 1
  stimulus.currentMask = 0;
else
  stimulus.currentMask = 1;
end

% if we are supposed to start halfway through
if stimulus.initialHalfCycle
  stimulus.currentMask = stimulus.currentMask+round(stimulus.stepsPerCycle/2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to draw retinotopy stimulus to screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = updateRetinotopyStimulus(stimulus,myscreen)

if any(stimulus.stimulusType == [3 4 5])

  % if doing standard bars, then draw sliding elements
  if stimulus.stimulusType == 3
    % update the phase of the sliding wedges
    stimulus.phaseNumRect = 1+mod(stimulus.phaseNumRect,stimulus.nRect);

    % draw the whole stimulus pattern, rotate to the element angle
    x = stimulus.xRect{stimulus.phaseNumRect};
    y = stimulus.yRect{stimulus.phaseNumRect};
    coords(1:2,:) = stimulus.elementRotMatrix*[x(1,:);y(1,:)];
    coords(3:4,:) = stimulus.elementRotMatrix*[x(2,:);y(2,:)];
    coords(5:6,:) = stimulus.elementRotMatrix*[x(3,:);y(3,:)];
    coords(7:8,:) = stimulus.elementRotMatrix*[x(4,:);y(4,:)];
    mglQuad(coords(1:2:8,:)+stimulus.xOffset,coords(2:2:8,:)+stimulus.yOffset,stimulus.cRect{stimulus.phaseNumRect},1);

    % compute the center of the bar
    barCenter = repmat(stimulus.barCenter(stimulus.currentMask,:),size(stimulus.maskBarLeft,1),1);
    % compute the left and right masks (covering up everything except the bar)
    % by shifting by the barCenter and rotating the coordinates for the angle we want
    maskBarLeft = stimulus.maskBarRotMatrix*(barCenter+stimulus.maskBarLeft)';
    maskBarRight = stimulus.maskBarRotMatrix*(barCenter+stimulus.maskBarRight)';

    % draw the bar masks
    mglPolygon(maskBarLeft(1,:)+stimulus.xOffset,maskBarLeft(2,:)+stimulus.yOffset,0.5);
    mglPolygon(maskBarRight(1,:)+stimulus.xOffset,maskBarRight(2,:)+stimulus.yOffset,0.5);
  elseif stimulus.stimulusType == 4
    % doing dots task, so draw dots, first clear screen
    mglClearScreen;
    % draw the dots
    xOffset = stimulus.xOffset+stimulus.barCenter(stimulus.currentMask,1);
    yOffset = stimulus.yOffset+stimulus.barCenter(stimulus.currentMask,2);
    drawDots(stimulus.dots.middle,stimulus.maskBarRotMatrix,xOffset,yOffset);
    drawDots(stimulus.dots.upper,stimulus.maskBarRotMatrix,xOffset,yOffset);
    drawDots(stimulus.dots.lower,stimulus.maskBarRotMatrix,xOffset,yOffset);
    % update the dots
    stimulus.dots.upper = updateDots(stimulus.dots.upper);
    stimulus.dots.middle = updateDots(stimulus.dots.middle);
    stimulus.dots.lower = updateDots(stimulus.dots.lower);
    % draw the fixation cross
    global fixStimulus;
    mglGluDisk(fixStimulus.pos(1),fixStimulus.pos(2),fixStimulus.diskSize*[1 1],myscreen.background,60);
    mglFixationCross(fixStimulus.fixWidth,fixStimulus.fixLineWidth,stimulus.fixColor,fixStimulus.pos);
  elseif stimulus.stimulusType == 5
       % doing dots task, so draw dots, first clear screen
    mglClearScreen;
    % draw the dots
    xOffset = stimulus.xOffset+stimulus.barCenter(stimulus.currentMask,1);
    yOffset = stimulus.yOffset+stimulus.barCenter(stimulus.currentMask,2);
    drawDots(stimulus.dots.middle,stimulus.maskBarRotMatrix,xOffset,yOffset);
    drawDots(stimulus.dots.upper,stimulus.maskBarRotMatrix,xOffset,yOffset);
    drawDots(stimulus.dots.lower,stimulus.maskBarRotMatrix,xOffset,yOffset);
    % update the dots
    stimulus.dots.upper = updateDots(stimulus.dots.upper);
    stimulus.dots.middle = updateDots(stimulus.dots.middle);
    stimulus.dots.lower = updateDots(stimulus.dots.lower);
      
  end
  
else
  % update the phase of the sliding wedges
  stimulus.phaseNum = 1+mod(stimulus.phaseNum,stimulus.n);
  % draw the whole stimulus pattern
  mglQuad(stimulus.x{stimulus.phaseNum}+stimulus.xOffset,stimulus.y{stimulus.phaseNum}+stimulus.yOffset,stimulus.c{stimulus.phaseNum},1);
  
  % mask out to get a wedge
  if stimulus.stimulusType == 1
    % For the mask is made of two overlapping semicircles.
    topX = stimulus.maskWedgeX{stimulus.currentMask}{1};
    topY = stimulus.maskWedgeY{stimulus.currentMask}{1};
    mglPolygon(topX+stimulus.xOffset,topY+stimulus.yOffset,0.5);

    bottomX = stimulus.maskWedgeX{stimulus.currentMask}{2};
    bottomY = stimulus.maskWedgeY{stimulus.currentMask}{2};
    mglPolygon(bottomX+stimulus.xOffset,bottomY+stimulus.yOffset,0.5);
    % or mask out to get a ring
  else
    mglPolygon(stimulus.maskInnerX{stimulus.currentMask}+stimulus.xOffset,stimulus.maskInnerY{stimulus.currentMask}+stimulus.yOffset,0.5);
    mglQuad(stimulus.maskOuterX{stimulus.currentMask}+stimulus.xOffset,stimulus.maskOuterY{stimulus.currentMask}+stimulus.yOffset,stimulus.maskOuterC{stimulus.currentMask});
  end
end
