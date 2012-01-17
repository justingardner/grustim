% sigdetect.m
%
%        $Id: sigdetect.m,v 1.6 2009/02/02 22:54:37 justin Exp $
%      usage: sigdetect
%         by: justin gardner & yuko
%       date: 11/25/09
%    purpose: sdt stimulus program using optic flow dots / faces
%    sigdetect(varargin)
%       'numTrials=10' for 10 trials.
%       'doEyeCalib' Default is -1 (at begining of each block)
%       'presentButton=1' for present to be 1 (default)
%       'absentButton=2' for absent to be 2 (default)
%       'feedback=0' gives correct/incorrect feedback for as many trials as it is set for.
%                    To always give feedback set to inf.
%       'presentProb=0.5' The probability for which a target is presetn
%       'subjectID=999' Saves in directory with that subjectID
%
%        sigdetect('staircase=1','staircaseType=quest');
%
%        sigdetect('staircase=1','staircaseType=fixed','fixedValues=[0.3 0.4 0.5 0.6 0.7]);
%
%        % To run signal detection with dots and coherence of 0.3
%        sigdetect('staircase=0','stimulusType=dots','strength=0.3')
% 
%        Same thing, but does sdt for faces
%        sigdetect('staircase=0','stimulusType=faces,'strength=0.3')
%
%        Note that for displaying images, a clipRange has been set arbitrarily
%        for the images in grustim/images/facesWithTransparentBackground.
%        It scales the image by 0.7 and then clips to 0 and 255. This seems
%        to work find, but it might be better to calculate the scale factor
%        on the fly for the set of images when they are loaded using loadNormalizedImages
%  nTrials                number of trials to run before quiting
%  doEyeCalib             run eye calibration
%  presentButtion      button to press when stimulus is present
%  absentButton        button to press when stimulus is absernt
%  correctSound           correct sound
%  incorrectSound         incorrect sound
%  stimSound              sound for when the stimulus period occurs
%  soundDir               where the sounds are to be found
%  strength               strength is the signal present strength for sigdetect trials
%  pedestals              pedestals for which to run staircase on
%  feedback               give feedback 
%  presentProb            change probability of stimulus present for sigdetect trials
%  staircase              run a staircase instead of a detection experiment
%  continueStaircase      continues staircases from one run to the other
%  stimulusType           dots or faces
%  staircaseType          quest or fixed
%  dispParams             display params, sometimes useful to turn off so that subejct's don't view what is going on
%  blankIntertrial        set to true to have no stimulus in between trials (FIX - does this work for images?)
%  nStaircaseTrials       number of trials in each staircase
%  scanner                scanner run. Sets ITIs, startDelay and waitForBacktick etc for parameters for magnet runs
%  startDelay             time to wait before starting experiment, useful for fMRI experiments
%  waitForBacktick        wait for backtick before starting, useful for scanner runs
function [myscreen stimulus] = sigdetect(varargin)

% parse arguments, note that stimulus variable gets
% set because it is a global
myscreen = parseArgs(varargin);
global stimulus;

% initalize the screen
myscreen.background = 'black';
myscreen = initScreen(myscreen);

% display some settings
if stimulus.runParams.dispParams,doDispParams(myscreen,stimulus);end

% first phase is just for having a start delay
task{1}.waitForBacktick = stimulus.runParams.waitForBacktick;
task{1}.seglen = stimulus.runParams.startDelay;
task{1}.numTrials = 1;
task{1}.parameter.strength = 0;

% set up the actual task
if ~stimulus.runParams.scanner
  % setup trials for either a staircase (with two intervals) or a sigdetect task
  if stimulus.params.staircase
    task{2}.seglen = [0.5 1 0.4 1 1];
    task{2}.getResponse = [0 0 0 0 1];
    task{2}.parameter.pedestal = stimulus.params.pedestals;
  else
    task{2}.seglen = [0.5 1 1.5];
    task{2}.getResponse = [0 1 1];
  end
else
  % non-psychophysics
  task{2}.segmin = [1 1 3];
  task{2}.segmax = [5 1 3]; % randomize fixation between 4 and 8 seconds
  task{2}.synchToVol = [1 0 0]; % sync at end of fixation
  task{2}.getResponse = [0 0 1];
  keyboard % FIX these
end
task{2}.random = 1;
task{2}.numTrials = stimulus.runParams.nTrials;
task{2}.waitForBacktick = 0;
task{2}.randVars.calculated.response = nan;
if stimulus.params.staircase
  task{2}.randVars.calculated.whichInterval = nan;
else
  task{2}.randVars.calculated.responseType = nan;
end
task{2}.randVars.calculated.strength = nan;

% initialize our task
[task{1} myscreen] = initTask(task{1},myscreen,@startSegmentCallback,@updateScreenCallback,@responseCallback);
[task{2} myscreen] = initTask(task{2},myscreen,@startSegmentCallback,@updateScreenCallback,@responseCallback);

% set up stimulus
if strcmp(stimulus.params.stimulusType ,'dots')
  % type can be Opticflow or Linear
  stimulus.dots.type = 'Opticflow'; 
  % set to 1 to make a circular mask around stimulus
  stimulus.dots.mask = 1;
  stimulus.dots.xcenter = 0;
  stimulus.dots.ycenter = 0;
  stimulus.dots.rmax = myscreen.imageHeight;
  stimulus = initDots(stimulus,myscreen);
elseif strcmp(stimulus.params.stimulusType ,'grating')
  % grating parameters
  stimulus.grating.radius = 9;
  stimulus.grating.n = 4;
  stimulus.grating.sf = 2;
  stimulus.grating.tf = 2;
  stimulus.grating.width = 12;
  stimulus.grating.height = 12;
  stimulus.grating.windowType = 'gabor'; % should be gabor or thresh
  stimulus.grating.sdx = stimulus.grating.width/7;
  stimulus.grating.sdy = stimulus.grating.width/7;
  % these are the reserved colors, if you need them later
  % you can display them by setting your color to the appropriate
  % index in stimulus.colors.reservedColor e.g. to get the
  % second color, in this case white, you would do
  % mglClearScreen(stimulus.colors.reservedColor(2));
  stimulus.colors.reservedColors = [0 0 0; 1 1 1; 1 1 0; 0 1 1;1 0 0;0 1 0;1 0 1];
  % init the grating stimulus.
  stimulus = initGratings(stimulus,myscreen);
  myscreen.background = stimulus.colors.grayColor;
  mglClearScreen(stimulus.colors.grayColor);
else
  stimulus = initFaces(stimulus,myscreen,stimulus.runParams.imageDir,320,240);
end

% sounds
mglInstallSound(stimulus.sounds.soundDir);

% other stimulus settings
stimulus.display = 1;

% size of image to display
stimulus.imageWidth = 24;
stimulus.imageHeight = 32;

% see if we can load staircase from last run
stimulus.s = [];
% see if we can continue from last run
if stimulus.runParams.continueStaircase
  % load last stimfile
  lastStimfile = getLastStimfile(myscreen);
  % see if we have a match
  if ~isempty(lastStimfile)
    % match, then grab that staircase
    if isfield(lastStimfile.stimulus,'params') && isequal(lastStimfile.stimulus.params,stimulus.params)
      disp(sprintf('(sigdetect) Continuing staircases from last run'));
      stimulus.s = lastStimfile.stimulus.s;
    else
      % no match, warn and go on
      dispHeader('',80);
      disp(sprintf('(sigdetect) !!! Last stimfile params does not match, starting new staircase !!!'));
      dispHeader('',80);
    end
  end
end

% init staircase
if isempty(stimulus.s)
  if stimulus.params.staircase
    nPedestals = length(stimulus.params.pedestals);
    % initialize staircase if we haven't been able to continue
    for i = 1:nPedestals
      if strcmp(stimulus.params.staircaseType,'quest')
	stimulus.s{i} = doStaircase('init','quest','initialThreshold=0.8','tGuessSd=4','nTrials',stimulus.runParams.nStaircaseTrials,'dispFig=1','subplotRows',nPedestals,'subplotNum',i);
      elseif strcmp(stimulus.params.staircaseType,'pest')
	stimulus.s{i} = doStaircase('init','upDown','stepRule=pest','initialThreshold=0.8','initialStepsize=0.2','nTrials',stimulus.runParams.nStaircaseTrials,'dispFig=1','minThreshold=0','maxThreshold=1','subplotRows',nPedestals,'subplotNum',i);
      elseif strcmp(stimulus.params.staircaseType,'fixed')
	stimulus.s{i} = doStaircase('init','fixed','fixedVals',stimulus.params.fixedValues,'nTrials',stimulus.runParams.nStaircaseTrials,'dispFig=1','subplotRows',nPedestals,'subplotNum',i);
      else
	disp(sprintf('(sigdetect) Unknown staircase type: %s',stimulus.params.staircaseType));
	keyboard
      end
    end
  else
    % sdt trials, use doStaircase for running sdt
    stimulus.s = doStaircase('init','sdt','strength',stimulus.params.strength,'dispFig=1','p',stimulus.params.presentProb);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task)) && ~myscreen.userHitEsc
  % update the dots
  [task myscreen phaseNum] = updateTask(task,myscreen,phaseNum);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

% display the data plots
if stimulus.params.staircase
  for i = 1:length(stimulus.params.pedestals)
    t(i) = doStaircase('threshold',stimulus.s{i});
  end
  smartfig('sigdetect','reuse');clf;
  plot(stimulus.params.pedestals,[t.threshold],'ko-','MarkerFaceColor','k','MarkerSize',9);
  yaxis(0,1);
  xlabel('pedestal');
  ylabel('threshold');
else
  t = doStaircase('threshold',stimulus.s,'dispFig=1');
end

% save log
if ~isempty(myscreen.stimfile)
  fid = fopen(fullfile(myscreen.datadir,'log.txt'),'a');
  if stimulus.params.staircase
    fprintf(fid,sprintf('%s: nTrials:%i staircaseType: %s threshold: %s',myscreen.stimfile,task{2}.trialnum,stimulus.params.staircaseType,mynum2str([t.threshold])));
  else
    fprintf(fid,sprintf('%s: nTrials:%i d-prime: %s',myscreen.stimfile,task{2}.trialnum,mynum2str(t.dprime,'sigfigs=3','compact=1')));
  end
  fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;
stimulus.fixColor = [1 1 1];
% when to show the stimulus (i.e. set strength high)
if stimulus.params.staircase
  % find which staircase we are working on.
  if isfield(task.thistrial,'pedestal')
    stimulus.snum = find(task.thistrial.pedestal == task.parameter.pedestal);
  else
    stimulus.strength = nan;
    stimulus.sigStrength = inf;
    return;
  end
  % staircase to find a strength threshold
  % set which interval the signal occurs in.
  if task.thistrial.thisseg == 1
    task.thistrial.whichInterval = (rand > 0.5)+1;
    % now set the strength for each trial
    [task.thistrial.strength stimulus.s{stimulus.snum}] = doStaircase('testValue',stimulus.s{stimulus.snum});
    % for faces create the signal and noise images
    if strcmp(stimulus.params.stimulusType,'faces')
      if ~isempty(stimulus.sigTex) mglDeleteTexture(stimulus.sigTex);end
      stimulus.sigTex = createScrambledFace(task.thistrial.strength+task.thistrial.pedestal);
      stimulus.sigStrength = task.thistrial.pedestal+task.thistrial.strength;
      if ~isempty(stimulus.noiseTex) mglDeleteTexture(stimulus.noiseTex);end
      stimulus.noiseTex = createScrambledFace(task.thistrial.pedestal);
    end
  end
  % set which strength should be displayed in which segment
  if stimulus.params.blankIntertrial,stimulus.display = true;end
  if (task.thistrial.thisseg == 2) && (task.thistrial.whichInterval == 1)
    stimulus.strength = task.thistrial.strength+task.thistrial.pedestal;
  elseif (task.thistrial.thisseg == 4) && (task.thistrial.whichInterval == 2)
    stimulus.strength = task.thistrial.strength+task.thistrial.pedestal;
  elseif any(task.thistrial.thisseg == [2 4])
    stimulus.strength = task.thistrial.pedestal;
  else
    if stimulus.params.blankIntertrial,stimulus.display = false;end
    stimulus.strength = [];
    stimulus.fixColor = [1 1 1];
  end
  % set fixation color
  if any(task.thistrial.thisseg == [2 4])
    stimulus.fixColor = [1 1 0];
    mglPlaySound(stimulus.sounds.stimSound);
  elseif task.thistrial.thisseg == 5
    stimulus.fixColor = [0 1 1];
  end
else
  % randomize whether stimulus will be presented or not
  if task.thistrial.thisseg == 1
    [task.thistrial.strength stimulus.s] = doStaircase('testValue',stimulus.s);
    if strcmp(stimulus.params.stimulusType,'faces')
      % remove old textures
      if ~isempty(stimulus.sigTex) mglDeleteTexture(stimulus.sigTex);end
      if ~isempty(stimulus.noiseTex) mglDeleteTexture(stimulus.noiseTex);end
      % create new ones
      stimulus.sigTex = createScrambledFace(task.thistrial.strength);
      stimulus.noiseTex = createScrambledFace(0);
    end
  end
  
  % show the stimulus in segment 2 when it is present
  if task.thistrial.thisseg == 2
    stimulus.strength = task.thistrial.strength;
    if stimulus.params.blankIntertrial,stimulus.display = true;end
  else
    stimulus.strength = nan;
    if stimulus.params.blankIntertrial,stimulus.display = false;end
  end
  % set fixation color
  if any(task.thistrial.thisseg == [2])
    stimulus.fixColor = [1 1 0];
    mglPlaySound(stimulus.sounds.stimSound);
  elseif task.thistrial.thisseg == 3
    stimulus.fixColor = [0 1 1];
  end
end

% set speed
if strcmp(stimulus.params.stimulusType,'dots')
  stimulus.dots = feval(sprintf('setDotsSpeed%s',stimulus.dots.type),stimulus.dots,stimulus.speed,myscreen);
  stimulus.dots = feval(sprintf('setDotsDir%s',stimulus.dots.type),stimulus.dots,2*mod(task.thistrial.thisseg,2)-1,myscreen);
end

% special stuff for gratings that deals with the gamma table
if strcmp(stimulus.params.stimulusType,'grating')
  setGammaTableForMaxContrast(stimulus.strength);
  % find the color in the reserved color array and get its index for display
  [matchColor reservedColorIndex] = intersect(stimulus.colors.reservedColors,stimulus.fixColor,'rows');
  if isempty(matchColor),disp(sprintf('(sigdetect) Missing fixation color in reserved colors'));end
  stimulus.fixColor = stimulus.colors.reservedColor(reservedColorIndex);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task, myscreen)

global stimulus
mglClearScreen;

% display the dots
if strcmp(stimulus.params.stimulusType,'dots')
  if isempty(stimulus.strength) || isnan(stimulus.strength), stimulus.strength = 0;end
  % update the dots
  stimulus.dots = feval(sprintf('updateDots%s',stimulus.dots.type),stimulus.dots,stimulus.strength,myscreen);

  % draw the dots
  if stimulus.display
    if stimulus.dots.mask,mglStencilSelect(1);end
    mglPoints2(stimulus.dots.x,stimulus.dots.y,stimulus.dots.dotsize,[1 1 1]);
    if stimulus.dots.mask,mglStencilSelect(0);end
  end
% display the grating stimulus
elseif strcmp(stimulus.params.stimulusType,'grating')
  if isempty(stimulus.strength) || isnan(stimulus.strength), stimulus.strength = 0;end
  if stimulus.display && stimulus.strength > 0
    % get the contrast index
    contrastIndex = getContrastIndex(stimulus.strength);
    
    % blt texture
    mglBltTexture(stimulus.tex(contrastIndex),[0 0 stimulus.grating.height]);

    % and mask with the gaussian
    mglBltTexture(stimulus.mask,[0 0]);

  end
% display the face images
else
  if stimulus.display
    if stimulus.strength == stimulus.sigStrength
      mglBltTexture(stimulus.sigTex,[0 0 stimulus.imageWidth stimulus.imageHeight]);
    else
      mglBltTexture(stimulus.noiseTex,[0 0 stimulus.imageWidth stimulus.imageHeight]);
    end
  end
end

% draw fixation cross
fixCrossSize = 2;
mglGluDisk(0,0,repmat(fixCrossSize*.6,1,2),myscreen.background,60);
mglFixationCross(fixCrossSize,4,stimulus.fixColor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    createScreambledFace    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function imageTex = createScrambledFace(scrambleFactor)

global stimulus;

% select which image to show
stimulus.faceNum(end+1) = ceil(stimulus.face.n*rand);

% get its half fourier representation
faceIm  = stimulus.face.halfFourier{stimulus.faceNum(end)};

% now create a hybrid image which contains scambleFactor percent of the 
% phase information from the place image
scrambledPhases = randperm(faceIm.n);
scrambledPhases = scrambledPhases(1:floor(faceIm.n*(1-scrambleFactor)));
% scramble the phases
faceIm.phase(scrambledPhases) = rand(1,length(scrambledPhases))*2*pi;

% get the average
faceIm.mag = stimulus.averageMag;

% create current image
imageTex = mglCreateTexture(flipud(clipRange(reconstructFromHalfFourier(faceIm))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called when subject responds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)

global stimulus

% only record the first response
if task.thistrial.gotResponse == 0
  % for staircase just get correct, incorrect and update staircase
  if stimulus.params.staircase
    task.thistrial.response = (task.thistrial.whichInterval == task.thistrial.whichButton);
    stimulus.s{stimulus.snum} = doStaircase('update',stimulus.s{stimulus.snum},task.thistrial.response);
    % if at the end of the staircase, reinitialize
    if doStaircase('stop',stimulus.s{stimulus.snum})
      stimulus.s{stimulus.snum}(end+1) = doStaircase('init',stimulus.s{stimulus.snum}(end));
    end
    if task.thistrial.response
      disp(sprintf('Trial %i: pedestal %s strength %f interval %i correct',task.trialnum,mynum2str(task.thistrial.pedestal),task.thistrial.strength,task.thistrial.whichInterval));
      stimulus.fixColor = [0 1 0];
    else
      disp(sprintf('Trial %i: pedestal %s strength %f interval %i incorrect',task.trialnum,mynum2str(task.thistrial.pedestal),task.thistrial.strength,task.thistrial.whichInterval));
      stimulus.fixColor = [1 0 0];
    end
    % deal with gamma table for gratings
    if strcmp(stimulus.params.stimulusType,'grating')
      if length(stimulus.fixColor) > 1
	% find the color in the reserved color array and get its index for display
	[matchColor reservedColorIndex] = intersect(stimulus.colors.reservedColors,stimulus.fixColor,'rows');
	if isempty(matchColor),disp(sprintf('(sigdetect) Missing fixation color in reserved colors'));end
	stimulus.fixColor = stimulus.colors.reservedColor(reservedColorIndex);
      end
    end
    return
  end
  %classify the respones into 4 SDT categories
  if task.thistrial.whichButton == stimulus.runParams.presentButton
    % HIT trial
    if task.thistrial.strength > 0
      task.thistrial.response = 1;
      % feedback
      if stimulus.runParams.feedback
	mglPlaySound(stimulus.sounds.correctSound); 
	stimulus.fixColor = [0 1 0];
      end
      responseType = 'Hit';
      % FALSE ALARM, incorrect
    else
      task.thistrial.response = 0;
      % feedback
      if stimulus.runParams.feedback
	mglPlaySound(stimulus.sounds.incorrectSound); 
	stimulus.fixColor = [1 0 0];
      end
      responseType = 'False Alarm';
    end
  %response absent
  elseif task.thistrial.whichButton == stimulus.runParams.absentButton
    % CORRECT REJECT, correct
    if task.thistrial.strength == 0
      task.thistrial.response = 1;
      % feedback
      if stimulus.runParams.feedback
	mglPlaySound(stimulus.sounds.correctSound); 
	stimulus.fixColor = [0 1 0];
      end
      responseType = 'Correct Reject';
      %MISS, incorrect
    else
      task.thistrial.response = 0;
      % feedback
      if stimulus.runParams.feedback
	mglPlaySound(stimulus.sounds.incorrectSound); 
	stimulus.fixColor = [1 0 0];
      end
      responseType = 'Miss';
    end
  end
  if any(task.thistrial.whichButton == [stimulus.runParams.presentButton stimulus.runParams.absentButton])
    %update staircase
    stimulus.s = doStaircase('update',stimulus.s,task.thistrial.response);;
    % keep response type
    task.thistrial.responseType = responseType(1);
    disp(sprintf('(sigdetect) %s strength=%f',responseType,task.thistrial.strength));
    if ~stimulus.runParams.feedback,stimulus.fixColor = [1 0 1];end
  end
  % turn off feedback after stimulus.feedback number of trials
  stimulus.runParams.feedback = max(stimulus.runParams.feedback-1,0);
end

% deal with gamma table for gratings
if strcmp(stimulus.params.stimulusType,'grating')
  if length(stimulus.fixColor) > 1
    % find the color in the reserved color array and get its index for display
    [matchColor reservedColorIndex] = intersect(stimulus.colors.reservedColors,stimulus.fixColor,'rows');
    if isempty(matchColor),disp(sprintf('(sigdetect) Missing fixation color in reserved colors'));end
    stimulus.fixColor = stimulus.colors.reservedColor(reservedColorIndex);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initDots(stimulus,myscreen)

stimulus.speed = 7;
% convert the passed in parameters to real units
if ~isfield(stimulus,'dots') || ~isfield(stimulus.dots,'rmax'), stimulus.dots.rmax = min(myscreen.imageWidth,myscreen.imageHeight);,end
if ~isfield(stimulus.dots,'xcenter'), stimulus.dots.xcenter = 0;,end
if ~isfield(stimulus.dots,'ycenter'), stimulus.dots.ycenter = 0;,end
if ~isfield(stimulus.dots,'dotsize'), stimulus.dots.dotsize = 4;,end
if ~isfield(stimulus.dots,'density'), stimulus.dots.density = 5;,end
if ~isfield(stimulus.dots,'coherence'), stimulus.dots.coherence = 1;,end
if ~isfield(stimulus.dots,'speed'), stimulus.dots.speed = stimulus.speed;,end
if ~isfield(stimulus.dots,'dir'), stimulus.dots.dir = 0;,end
if ~isfield(stimulus.dots,'mask'), stimulus.dots.mask = 0;,end

% update the dots
stimulus.dots = feval(sprintf('initDots%s',stimulus.dots.type),stimulus.dots,myscreen);

% set color
stimulus.dots.color = ones(stimulus.dots.n,1);
%stimulus.dots.color(rand(1,stimulus.dots.n)>0.5) = 1;

% create stencil
if stimulus.dots.mask
  mglClearScreen;
  mglStencilCreateBegin(1);
  % and draw that oval
  mglGluDisk(stimulus.dots.xcenter,stimulus.dots.ycenter,[stimulus.dots.rmax stimulus.dots.rmax]/2,[1 1 1],60);
  mglStencilCreateEnd;
  mglClearScreen;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the dots speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = setDotsSpeedLinear(dots,speed,myscreen)

% get the step size
dots.speed = speed;
dots.stepsize = speed/myscreen.framesPerSecond;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the dots direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = setDotsDirLinear(dots,direction,myscreen)

% get the step size
dots.dir = direction;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots for linear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updateDotsLinear(dots,coherence,myscreen)

% get the dots step
dots.xstep = cos(pi*dots.dir/180)*dots.stepsize;
dots.ystep = sin(pi*dots.dir/180)*dots.stepsize;

% pick a random set of dots
dots.coherent = rand(1,dots.n) < coherence;

% now move those dots in the right direction
dots.x(dots.coherent) = dots.x(dots.coherent)+dots.xstep;
dots.y(dots.coherent) = dots.y(dots.coherent)+dots.ystep;

% randomwalk rule
thisdir = rand(1,sum(~dots.coherent))*2*pi;
dots.x(~dots.coherent) = dots.x(~dots.coherent)+cos(thisdir)*dots.stepsize;
dots.y(~dots.coherent) = dots.y(~dots.coherent)+sin(thisdir)*dots.stepsize;

% movshon noise
%dots.x(~dots.coherent) = rand(1,sum(~dots.coherent))*dots.width;
%dots.y(~dots.coherent) = rand(1,sum(~dots.coherent))*dots.height;

% make sure we haven't gone off the patch
% do the dots separately for left and right hand side
dots.x(dots.x < dots.xmin) = dots.x(dots.x < dots.xmin)+dots.width;
dots.x(dots.x > dots.xmax) = dots.x(dots.x > dots.xmax)-dots.width;
dots.y(dots.y < dots.ymin) = dots.y(dots.y < dots.ymin)+dots.height;
dots.y(dots.y > dots.ymax) = dots.y(dots.y > dots.ymax)-dots.height;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots for opticflow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updateDotsOpticflow(dots,coherence,myscreen)

% get the coherent and incoherent dots
%if (dots.coherency ~= coherence)
  dots.incoherent = rand(1,dots.n) > coherence;
  dots.incoherentn = sum(dots.incoherent);
  dots.coherent = ~dots.incoherent;
  dots.coherency = coherence;
  % generate a random transformation matrix for each incoherent point
  dots.randT = rand(3,dots.incoherentn)-0.5;
  % and normalize the transformation to have the same length
  % (i.e. speed) as the real transformation matrix
  dots.randT = sqrt(sum(dots.T.^2))*dots.randT./([1 1 1]'*sqrt(sum(dots.randT.^2)));
%end

% update relative position of dots in 3-space to observer
dots.X(dots.coherent) = dots.X(dots.coherent)-dots.T(1);
dots.Y(dots.coherent) = dots.Y(dots.coherent)-dots.T(2);
dots.Z(dots.coherent) = dots.Z(dots.coherent)-dots.T(3);

% now move the incoherent points according to the random trasnformation
dots.X(dots.incoherent) = dots.X(dots.incoherent)-dots.randT(1,:);
dots.Y(dots.incoherent) = dots.Y(dots.incoherent)-dots.randT(2,:);
dots.Z(dots.incoherent) = dots.Z(dots.incoherent)-dots.randT(3,:);

% get all points that have fallen off the screen
offscreen = dots.Z<dots.minZ;

% and put them at the furthest distance
dots.Z(offscreen) = dots.maxZ;

% get all points that have fallen out of view
offscreen = dots.Z>dots.maxZ;
% and move them to the front plane
dots.Z(offscreen) = dots.minZ;

% put points fallen off the X edge back
offscreen = dots.X < -dots.maxX;
dots.X(offscreen) = dots.X(offscreen)+2*dots.maxX;
offscreen = dots.X > dots.maxX;
dots.X(offscreen) = dots.X(offscreen)-2*dots.maxX;

% put points fallen off the Y edge back
offscreen = dots.Y < -dots.maxY;
dots.Y(offscreen) = dots.Y(offscreen)+2*dots.maxY;
offscreen = dots.Y > dots.maxY;
dots.Y(offscreen) = dots.Y(offscreen)-2*dots.maxY;

% project on to screen
dots.xproj = dots.f*dots.X./dots.Z;
dots.yproj = dots.f*dots.Y./dots.Z;

% stuff to compute median speed
%dots.oldx = dots.x;
%dots.oldy = dots.y;

% get actual screen coordinates
dots.x = dots.xproj*myscreen.imageWidth;
dots.y = dots.yproj*myscreen.imageHeight;

%medianSpeed = median(sqrt((dots.oldx-dots.x).^2+(dots.oldy-dots.y).^2)*myscreen.framesPerSecond);
%minSpeed = min(sqrt((dots.oldx-dots.x).^2+(dots.oldy-dots.y).^2)*myscreen.framesPerSecond);
%disp(sprintf('min: %f median: %f',minSpeed,medianSpeed));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create dots for linear2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initDotsLinear(dots,myscreen)

% actually a square patch of dots that get stenciled
% so calculate width and height
dots.width = dots.rmax*2;
dots.height = dots.rmax*2;

% get the number of dots
dots.n = round(dots.width*dots.height*dots.density);

% get max and min points for dots
dots.xmin = -dots.width/2;
dots.xmax = dots.width/2;
dots.ymin = -dots.height/2;
dots.ymax = dots.height/2;

% get initial position
dots.x = rand(1,dots.n)*dots.width;
dots.y = rand(1,dots.n)*dots.height;

% get the step size
dots.stepsize = dots.speed/myscreen.framesPerSecond;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the dots speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = setDotsSpeedOpticflow(dots,speed,myscreen)

% get the step size
dots.speed = speed;
dots.T = [0 0 dots.speed/myscreen.framesPerSecond];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the dots direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = setDotsDirOpticflow(dots,direction,myscreen)

% get the step size
dots.T = [0 0 direction*dots.speed/myscreen.framesPerSecond];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create dots for optic flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initDotsOpticflow(dots,myscreen)

% focal length to projection plane
% projection plane is defined to be 
% 1 unit wide and high, so with 
% this focal length, we are looking at
% a view of the world with a 90 deg fov
dots.f = .5;

% translation and rotation matrices
dots.T = [0 0 dots.speed/myscreen.framesPerSecond];
dots.R = [0 0 0];

% maximum depth of points
dots.maxZ = 10;dots.minZ = dots.f;
dots.maxX = 10;
dots.maxY = 10;

% make a brick of points
dots.n = round(myscreen.imageWidth*myscreen.imageHeight*dots.density);

% initial position of dots
dots.X = 2*dots.maxX*rand(1,dots.n)-dots.maxX;
dots.Y = 2*dots.maxY*rand(1,dots.n)-dots.maxY;
dots.Z = (dots.maxZ-dots.minZ)*rand(1,dots.n)+dots.minZ;

% get projection on to plane
dots.xproj = dots.f*dots.X./dots.Z;
dots.yproj = dots.f*dots.Y./dots.Z;

% put into screen coordinates
dots.x = dots.xproj*myscreen.imageWidth;
dots.y = dots.yproj*myscreen.imageHeight;

% set incoherent dots to 0
dots.coherency = 1;
dots.incoherent = rand(1,dots.n) > dots.coherency;
dots.incoherentn = sum(dots.incoherent);
dots.coherent = ~dots.incoherent;

dots.randT = zeros(3,dots.incoherentn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initFaces(stimulus,myscreen,faceDir,width,height)

if isfield(stimulus,'imagesLoaded') && (stimulus.imagesLoaded)
  disp(sprintf('(facesdt) Stimulus already initialized'));
  return
end

% face dir
stimulus.face = loadNormalizedImages(faceDir,'width',width,'height',height);
if isempty(stimulus.face)
  disp(sprintf('(facesdt) Could not load face images; %s',faceDir));
  keyboard
end

% get average mag
stimulus.averageMag = stimulus.face.averageMag;

stimulus.imagesLoaded = 1;
stimulus.faceNum = [];
stimulus.sigTex = [];
stimulus.noiseTex = [];

%%%%%%%%%%%%%%%%%%%
%    clipRange    %
%%%%%%%%%%%%%%%%%%%
function im = clipRange(im)

% new version (July 2011) to keep luminance and rms contrast
% constant across images. Scaling factor 0.7 chosen for a specific set
% of stimuli (the set of images that was in grustim/images/ObjLocImages/
% in June 2011)

im = im*0.7;
im(find(im>255)) = 255;
im(find(im<0)) = 0;

%%%%%%%%%%%%%%%%%%%%%%
%%   doDispParams   %%
%%%%%%%%%%%%%%%%%%%%%%
function doDispParams(myscreen,stimulus)

dispHeader(sprintf('subjectID: %s dataDir: %s',myscreen.subjectID,myscreen.datadir),80);
if ~stimulus.params.staircase
  disp(sprintf('(sigdetect) Running signal detection task: %s',stimulus.params.stimulusType));
  disp(sprintf('(sigdetect) Testing strength: %s',mynum2str(stimulus.params.strength)));
  disp(sprintf('(sigdetect) Probability: %s',mynum2str(stimulus.params.presentProb)));
else
  disp(sprintf('(sigdetect) Running 2AFC discrimination task: %s',stimulus.params.stimulusType));
  disp(sprintf('(sigdetect) Pedestals: %s',mynum2str(stimulus.params.pedestals)));
  disp(sprintf('(sigdetect) Staircase type: %s',stimulus.params.staircaseType));
end
disp(sprintf('(sigdetect) startDelay: %s waitForBacktick: %i doEyeCalib: %i',mynum2str(stimulus.runParams.startDelay),stimulus.runParams.waitForBacktick,stimulus.runParams.doEyeCalib));
disp(sprintf('(sigdetect) soundDir: %s correctSound: %s incorrectSound: %s stimSound: %s',stimulus.sounds.soundDir,stimulus.sounds.correctSound,stimulus.sounds.incorrectSound,stimulus.sounds.stimSound));
disp(sprintf('(sigdetect) blankIntertrial: %i',stimulus.params.blankIntertrial));
dispHeader('',80);

%%%%%%%%%%%%%%%%%%%%
%%   paraseArgs   %%
%%%%%%%%%%%%%%%%%%%%
function [myscreen stimulus] = parseArgs(args)

% evaluate the arguments
nTrials = [];        % number of trials to run before quiting
doEyeCalib = [];     % run eye calibration
presentButton=[];    % button to press when stimulus is present
absentButton = [];   % button to press when stimulus is absernt
correctSound = [];   % correct sound
incorrectSound = []; % incorrect sound
stimSound = [];      % sound for when the stimulus period occurs
soundDir = [];       % where the sounds are to be found
strength = [];       % strength is the signal present strength for sigdetect trials
pedestals = [];      % pedestals for which to run staircase on
feedback = [];       % give feedback 
presentProb = [];    % change probability of stimulus present for sigdetect trials
staircase = [];      % run a staircase instead of a detection experiment
continueStaircase=[];% continues staircases from one run to the other
stimulusType = [];   % dots or faces
staircaseType=[];    % quest or fixed
dispParams=[];       % display params, sometimes useful to turn off so that subejct's don't view what is going on
blankIntertrial=[];  % set to true to have no stimulus in between trials (FIX - does this work for images?)
nStaircaseTrials=[]; % number of trials in each staircase
scanner = [];        % scanner run. Sets ITIs, startDelay and waitForBacktick etc for parameters for magnet runs
startDelay = [];     % time to wait before starting experiment, useful for fMRI experiments
waitForBacktick=[];  % wait for backtick before starting, useful for scanner runs
getArgs(args,{'startDelay=[]','nTrials=500','doEyeCalib=1','presentButton=1','absentButton=2','correctSound=Pop','incorrectSound=Basso','stimSound=stimsound','soundDir=~/proj/grustim/sounds','strength=1','feedback=0','presentProb=0.5','staircase=0','stimulusType=faces','imageDir=~/proj/grustim/images/facesWithTransparentBackground','subjectID=999','staircaseType=quest','fixedValues=[0.1 0.2 0.3 0.4 0.5 0.6 0.7]','dispParams=1','blankIntertrial=[]','nStaircaseTrials=50','scanner=0','waitForBacktick=[]','pedestals=0','continueStaircase=1'});

% init the stimulus with current settings
clear global stimulus;
global stimulus;
stimulus.params.stimulusType = stimulusType;
stimulus.params.staircase = staircase;
stimulus.params.blankIntertrial = blankIntertrial;
if stimulus.params.staircase
  stimulus.params.staircaseType = staircaseType;
  stimulus.params.fixedValues = fixedValues;
  stimulus.params.pedestals = pedestals;
  myscreen.subjectFolder = 'staircase';
  if isempty(stimulus.params.blankIntertrial),stimulus.params.blankIntertrial = true;end
else
  stimulus.params.strength = strength;
  stimulus.params.presentProb = presentProb;
  myscreen.subjectFolder = 'sdt';
  if isempty(stimulus.params.blankIntertrial),stimulus.params.blankIntertrial = false;end
end

% change the settings based on whether we are in the scanner or not
if scanner
  if isempty(waitForBacktick),waitForBacktick = true;end
  if isempty(startDelay),startDelay = 10;end
else
  if isempty(waitForBacktick),waitForBacktick = false;end
  if isempty(startDelay),startDelay = 0.1;end
end

% set run params
stimulus.runParams.startDelay = startDelay;
stimulus.runParams.nTrials = nTrials;
stimulus.runParams.doEyeCalib = doEyeCalib;
stimulus.runParams.presentButton = presentButton;
stimulus.runParams.absentButton = absentButton;
stimulus.runParams.feedback = feedback;
stimulus.runParams.imageDir = imageDir;
stimulus.runParams.dispParams = dispParams;
stimulus.runParams.nStaircaseTrials = nStaircaseTrials;
stimulus.runParams.scanner = scanner;
stimulus.runParams.waitForBacktick = waitForBacktick;
stimulus.runParams.continueStaircase = continueStaircase;

% set sound parameters
stimulus.sounds.soundDir = soundDir;
stimulus.sounds.correctSound = correctSound;
stimulus.sounds.incorrectSound = incorrectSound;
stimulus.sounds.stimSound = stimSound;

% create subject directory
if isscalar(subjectID),subjectID = sprintf('s%03i',subjectID);end
myscreen.subjectID = subjectID;

% register stimulus name
myscreen = initStimulus('stimulus',myscreen);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     STIMULUS Gratings
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
gaussianWin = mglMakeGaussian(stimulus.grating.width+1,stimulus.grating.height+1,stimulus.grating.sdx,stimulus.grating.sdy);
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

