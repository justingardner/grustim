% alaisburr.m
%
%      usage: myscreen = alaisburr()
%         by: minyoung lee
%       date: 
%    purpose: replication of Alais & Burr, 2004
%
%            Set for visual / auditory or bimodal run:
%            alaisburr('visual=1');
%            alaisburr('auditory=1');
%            alaisburr('bimodal=1');
%
%            set to run a test of the gamma settings - this is important
%            because new versions of the operating system do not seem
%            to set the gamma correctly. It will put up a screen with
%            reserved colors at the top and on the bottom the gradation
%            of colors used for the gaussian
%
%            alaisburr('doGammaTest=1');
%
%            To run with a noise background with SNR of 2.5 that updates 4.5 every sec
%
%            alaisBurr('SNR=2.5','backgroundFreq=4.5');
%
function myscreen = alaisburr(varargin)
 
clear global stimulus
mglEatKeys('12`');
global stimulus
 
% get arguments
bimodal = 0;
getArgs(varargin,{'width=6','visual=0','auditory=0','bimodal=0','dispPlots=1','auditoryTrain=0','visualTrain=0','tenbit=1','doGammaTest=0','stimulusContrast=0.5','SNR=3','doTestSNR=0','backgroundFreq=4.5','doTestStimSize=0'},'verbose=1');

% close screen if open - to make sure that gamma gets sets correctly
mglClose;

% set task
if visual || visualTrain
    stimulus.task = 1;
elseif auditory || auditoryTrain
    stimulus.task = 2;
else
    stimulus.task = 3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% set arguments in stimulus variable
%%%%%%%%%%%%%%%%%%%%%%%%%
stimulus.disp = dispPlots;
stimulus.visual=visual;
stimulus.auditory=auditory;
stimulus.bimodal=bimodal;
stimulus.auditoryTrain = auditoryTrain;
stimulus.visualTrain = visualTrain;
stimulus.tenbit = tenbit;
stimulus.SNR = SNR;
%%%%%%%%%%%%%%%%%%%%%%%%%
stimulus.width = width;
stimulus.stimDur = 0.015; % 15ms
stimulus.gaussainDur = 0.015; % 15ms
stimulus.clickDur = 0.0015; % 1.5ms
stimulus.samplesPerSecond = 44100;
stimulus.ISI = .500; % 500ms
stimulus.trialnum = 0;
if stimulus.tenbit
    if stimulus.width >= 64
        stimulus.contrast = .002;%.01; 
    else
       stimulus.contrast = .005;
    end
else
  stimulus.contrast = 0.5;
end

% set the stimulus contrast
stimulus.contrast = stimulusContrast;

stimulus.interval = [2 4];
% fixation cross
stimulus.fixWidth = 1;
stimulus.fixColor = [1 1 1];
stimulus.colors.reservedColors = [1 1 1; 0.3 0.3 0.3; 0 1 0;1 0 0; 0 1 1];

if stimulus.auditoryTrain || stimulus.visualTrain
    stimulus.initialThreshold = 10;
    stimulus.initialStepsize = 2.5;
    stimulus.minThreshold = 0;
    stimulus.maxThreshold = 15;
    stimulus.minStepsize = 0.75;
    stimulus.maxStepsize = 5;
end

screenParams = mglGetScreenParams;
stimulus.displayDistance = screenParams{1}.displayDistance*.01;

% initalize the screen
myscreen.background = 0;  %black
myscreen = initScreen;

%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%% 
task{1}{1}.waitForBacktick = 1;
% trial: Fixation + Stim1 (.015s/.0015s) + ISI (.5s) + Stim2 (.015s/.0015s) + Resp + ITI
task{1}{1}.segmin = [1 stimulus.stimDur stimulus.ISI stimulus.stimDur 1.5 1];
task{1}{1}.segmax = [1 stimulus.stimDur stimulus.ISI stimulus.stimDur 1.5 1];
task{1}{1}.getResponse = [0 0 0 0 1 0];
if stimulus.bimodal
  task{1}{1}.numBlocks = 2;
elseif stimulus.visual || stimulus.auditory
  task{1}{1}.numBlocks = 16;
else
  task{1}{1}.randVars.uniform.sign = [1,-1];
end
% parameters & randomization
task{1}{1}.parameter.centerWhich = [1 2]; % centered in which interval
task{1}{1}.random = 1;
task{1}{1}.parameter.posDiff = [-18 -15 -12 -10 -8 -6 -4 -2.5 -1.25 0 1.25 2.5 4 6 8 10 12 15 18]; 
if stimulus.task == 3
  task{1}{1}.parameter.displacement = [-5 0 5];
end

task{1}{1}.randVars.calculated.resp = nan;
task{1}{1}.randVars.calculated.correct = nan;
task{1}{1}.randVars.calculated.diff = nan;
task{1}{1}.randVars.calculated.rt = nan;
task{1}{1}.randVars.calculated.centerInt = nan;
task{1}{1}.randVars.calculated.displ = nan;
task{1}{1}.randVars.calculated.left = nan;
 
% initialize the task
for phaseNum = 1:length(task{1})
  [task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end
 
% init the stimulus
myscreen = initStimulus('stimulus',myscreen);

% to initialize the stimulus for your experiment.
stimulus = initGaussian(stimulus,myscreen);

% init auditory stimulus
stimulus = initClick(stimulus,myscreen);
if stimulus.auditoryTrain || stimulus.visualTrain
  stimulus = initStair(stimulus);
end

% check gamma if called for
if doGammaTest && tenbit
  tf = testGammaTable(stimulus,myscreen);
  if ~tf,mglClose,return,end;
end

if doTestStimSize
  % clear screen
  mglClearScreen(stimulus.colors.black);
  % draw gaussian
  mglBltTexture(stimulus.tex,[0 0]);
  % draw fixation cross
  mglFixationCross(1,1,stimulus.colors.white,[ 0 0]);
  % draw circle around gaussian width
  x0 = stimulus.width * cos(0);
  y0 = stimulus.width * sin(0);
  for iDeg = 1:360
    x1 = stimulus.width * cos(d2r(iDeg));
    y1 = stimulus.width * sin(d2r(iDeg));
    mglLines2(x0, y0, x1, y1, 2, stimulus.colors.red);
    x0 = x1;
    y0 = y1;
  end
  mglFlush;
  disp(sprintf('(alaisBurr) Screen width: %0.1f Screen height: %0.1f Gaussian width: %0.1f',myscreen.imageWidth, myscreen.imageHeight, width));
  if ~askuser('Is this ok'), return,end
end

% if test then display several levels of SNR
if doTestSNR
  % set to test these SNR levels
  stimulus.SNR = [1];
  % setup background
  stimulus = initBackgroundNoise(stimulus, myscreen);
  stimulus = setBackgroundNoise(stimulus,myscreen,task,stimulus.SNR,1);
  % now cycle through each and display
  for iSNR = 1:length(stimulus.SNR);
    % create texture
    stimulus = setStimulusOnBackground(stimulus,0,0,1,1,stimulus.SNR(iSNR));
    % and blt texture
    mglBltTexture(stimulus.stimTexture(1),[0 0]);mglFlush;
    if ~askuser(sprintf('(alaisburr) Testing SNR of %0.2f',stimulus.SNR(iSNR)))
      break;
    end
  end
  mglClose;return
end

% init the background noise if we have to
if ~isinf(stimulus.SNR)
  % init background noise
  stimulus = initBackgroundNoise(stimulus, myscreen);
  % and set them
  stimulus = setBackgroundNoise(stimulus,myscreen,task,stimulus.SNR,backgroundFreq);
end

% put up display string
mglWaitSecs(1);
mglClearScreen(stimulus.colors.black);
mglTextSet([],32,stimulus.colors.white);
mglTextDraw('Press ` key to start when you are ready',[0 0]);
mglFlush;
mglClearScreen(stimulus.colors.black);
mglTextDraw('Press ` key to start when you are ready',[0 0]);
mglFlush;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
stimulus.endflag = 0;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc && ~stimulus.endflag
  % update the task
  [task{1} myscreen phaseNum] = updateTask(task{1},myscreen,phaseNum);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end
 
% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

% display psychometric functions
if stimulus.disp
  % first check if there are any trials with responses
  e = getTaskParameters(myscreen,task);
  if any(~isnan(e{end}.response))
    dispPsychometric(task{1}{1},stimulus);
  else
    disp(sprintf('(alaisburr) No subject responses to plot psychometric function with.'));
  end
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)
global stimulus

% set fixation cross color
stimulus.fixColor = stimulus.colors.white;

if task.thistrial.thisseg == 1
  % set random jittering between -0.5 and 0.5 deg
  task.thistrial.jitter = rand - 0.5; 
  % horizontal position of first, second stim
  if stimulus.auditoryTrain || stimulus.visualTrain
    % get test value
    [testValue, stimulus.stair] = doStaircase('testValue', stimulus.stair);
    task.thistrial.diff = testValue * task.thistrial.sign;

    % Which stimulus should be displayed in the center (i.e. reference position)
    if task.thistrial.centerWhich == 1
      task.thistrial.xposV = [task.thistrial.jitter, task.thistrial.diff + task.thistrial.jitter];
      task.thistrial.xposA = task.thistrial.xposV;
      task.thistrial.centerInt = 1;
    else
      task.thistrial.xposV = [task.thistrial.jitter + task.thistrial.diff, task.thistrial.jitter];
      task.thistrial.xposA = task.thistrial.xposV;
      task.thistrial.centerInt = 2;
    end
  else
    if stimulus.task == 3
      if task.thistrial.centerWhich == 1
	task.thistrial.xposV = [task.thistrial.jitter, task.thistrial.posDiff + task.thistrial.jitter + task.thistrial.displacement];
	task.thistrial.xposA = [task.thistrial.jitter, task.thistrial.posDiff + task.thistrial.jitter - task.thistrial.displacement];
	task.thistrial.centerInt = 1;
      else
	task.thistrial.xposV = [task.thistrial.jitter + task.thistrial.posDiff + task.thistrial.displacement, task.thistrial.jitter];
	task.thistrial.xposA = [task.thistrial.jitter + task.thistrial.posDiff - task.thistrial.displacement, task.thistrial.jitter];
	task.thistrial.centerInt = 2;
      end
      task.thistrial.displ = task.thistrial.displacement;
    else
      if task.thistrial.centerWhich == 1
	task.thistrial.xposV = [task.thistrial.jitter, task.thistrial.posDiff + task.thistrial.jitter];
	task.thistrial.xposA = task.thistrial.xposV;
	task.thistrial.centerInt = 1;
      else
	task.thistrial.xposV = [task.thistrial.jitter + task.thistrial.posDiff, task.thistrial.jitter];
	task.thistrial.xposA = task.thistrial.xposV;
	task.thistrial.centerInt = 2;
      end
    end
    task.thistrial.diff = task.thistrial.posDiff;
  end
    
  % auditory or bimodal condition
  if stimulus.task ~= 1 
    % create the necessary ITD for the sound position
    for int = 1:2
      stimulus.sound(int) = createITD(stimulus,task.thistrial.xposA(int));
    end
  end
  
  % setup background noise
  if ~isinf(stimulus.SNR)
    % get the order of the background noise images to display
    stimulus.background.frameOrder = randperm(stimulus.background.n);
    stimulus.background.frameNum = 0;
    % set the frame timer
    stimulus.background.frameStart = -inf;
    % figure out which frame of noise that the stimulus will be imbeded on
    stimulus.background.stim1Frame = ceil(sum(task.segmax(1:2))/stimulus.background.frameTime)-1;
    stimulus.background.stim2Frame = ceil(sum(task.segmax(1:4))/stimulus.background.frameTime)-1;
    % and create the stimuli on the background that we guess to be the
    % one that should be being presented
    stimulus = setStimulusOnBackground(stimulus,task.thistrial.xposV(1),0,1,stimulus.background.frameOrder(stimulus.background.stim1Frame),stimulus.SNR);
    stimulus = setStimulusOnBackground(stimulus,task.thistrial.xposV(2),0,2,stimulus.background.frameOrder(stimulus.background.stim2Frame),stimulus.SNR);
  end

elseif any(task.thistrial.thisseg == [2 4])
  % turn fixation color gray
  stimulus.fixColor = stimulus.colors.grey;
end
if task.thistrial.thisseg == 6
  if exist('task.thistrial.reactionTime', 'var')
    task.thistrial.rt = task.thistrial.reactionTime;
  end
  stimulus.trialnum = stimulus.trialnum+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)
global stimulus

% clear screen and put up fixation cross
mglClearScreen(stimulus.colors.black);

% visual or bimodal condition
if stimulus.task ~= 2 
  % display visual stimulus
  % on a blank background if snr is inf
  if isinf(stimulus.SNR)
    if task.thistrial.thisseg == stimulus.interval(1)
      mglBltTexture(stimulus.tex, [task.thistrial.xposV(1), 1]);
    elseif task.thistrial.thisseg == stimulus.interval(2)
      mglBltTexture(stimulus.tex, [task.thistrial.xposV(2), 1]);
    end
  else
    % otherwise display on noise background
    if task.thistrial.thisseg == stimulus.interval(1)
      mglBltTexture(stimulus.stimTexture(1),[0 0]);
      if stimulus.background.frameNum ~= stimulus.background.stim1Frame
	disp(sprintf('!!! (alaisburr) Stimulus is being displayed on a different noisy background then what is currently being presented. You should adjust the backgroundFreq until this no longer happens'));
      end
    elseif task.thistrial.thisseg == stimulus.interval(2)
      mglBltTexture(stimulus.stimTexture(2),[0 0]);
      if stimulus.background.frameNum ~= stimulus.background.stim2Frame
	disp(sprintf('!!! (alaisburr) Stimulus is being displayed on a different noisy background then what is currently being presented. You should adjust the backgroundFreq until this no longer happens'));
      end
    else
      % display background
      % see if we need to update frame number
      if mglGetSecs(stimulus.background.frameStart) > stimulus.background.frameTime
	% update the count
	stimulus.background.frameNum = mod(stimulus.background.frameNum,stimulus.background.n)+1;
	% reset timer
	stimulus.background.frameStart = mglGetSecs;
      end
      % draw the background texture
      mglBltTexture(stimulus.backTexture(stimulus.background.frameOrder(stimulus.background.frameNum)),[0 0]);
    end
  end
end

% auditory or bimodal condition
if stimulus.task ~= 1 
  if task.thistrial.thisseg == stimulus.interval(1)
    mglPlaySound(stimulus.sound(1));
  elseif task.thistrial.thisseg == stimulus.interval(2)
    mglPlaySound(stimulus.sound(2));
  end
end

mglFixationCross(stimulus.fixWidth,1.5,stimulus.fixColor);

% %draw fixation cross
% if task.thistrial.thisseg == 5 || task.thistrial.thisseg == 6
%     mglFixationCross(stimulus.fixWidth,1.5,stimulus.fixColor + stimulus.colors.reservedColor(2));
% else
%     mglFixationCross(stimulus.fixWidth,1.5,stimulus.fixColor);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)
global stimulus
% here, we just check whether this is the first time we got a response
if ~task.thistrial.gotResponse
    % which one seemed more to the LEFT
    % centerWhich (1/2) first or second one more eccentric
    if (task.thistrial.centerWhich == 1 && task.thistrial.whichButton == 2) || ...
    (task.thistrial.centerWhich == 2 && task.thistrial.whichButton == 1)
    task.thistrial.left = 'L';
    else
    task.thistrial.left = 'R';
    end
    if (task.thistrial.centerWhich == 1 && ((task.thistrial.diff > 0 && task.thistrial.whichButton == 1) || ...
      (task.thistrial.diff < 0 && task.thistrial.whichButton == 2))) || ...
    (task.thistrial.centerWhich == 2 && ((task.thistrial.diff > 0 && task.thistrial.whichButton == 2) || ...
      (task.thistrial.diff < 0 && task.thistrial.whichButton == 1)))
        % correct
        task.thistrial.correct = 1;
        if stimulus.auditoryTrain || stimulus.visualTrain
        % feeback
        stimulus.fixColor = stimulus.colors.green;%[0 1 0];
        end
        if ~stimulus.bimodal
            disp(sprintf('(alaisburr) Trial %i: %0.2f %c correct centerInt %i resp %i', ...
             task.trialnum, task.thistrial.diff, task.thistrial.left, task.thistrial.centerWhich, task.thistrial.whichButton))
        else
            disp(sprintf('(alaisburr) Trial %i: %0.2f %i %c correct centerInt %i resp %i', ...
             task.trialnum, task.thistrial.diff, task.thistrial.displacement, task.thistrial.left, task.thistrial.centerWhich, task.thistrial.whichButton))
        end
    else
        % incorrect
        task.thistrial.correct = 0;
        if stimulus.auditoryTrain || stimulus.visualTrain
        stimulus.fixColor = stimulus.colors.red;%[1 0 0];
        end
        if ~stimulus.bimodal
            disp(sprintf('(alaisburr) Trial %i: %0.2f %c incorrect centerInt %i resp %i', ...
             task.trialnum, task.thistrial.diff, task.thistrial.left, task.thistrial.centerWhich, task.thistrial.whichButton))
        else
            disp(sprintf('(alaisburr) Trial %i: %0.2f %i %c incorrect centerInt %i resp %i', ...
             task.trialnum, task.thistrial.diff, task.thistrial.displacement, task.thistrial.left, task.thistrial.centerWhich, task.thistrial.whichButton))
        end
    end
        
    task.thistrial.resp = task.thistrial.whichButton;
    task.thistrial.rt = task.thistrial.reactionTime;

    if ~(stimulus.auditoryTrain || stimulus.visualTrain)
      stimulus.fixColor = stimulus.colors.cyan;
    else
      stimulus.stair = doStaircase('update', stimulus.stair, task.thistrial.correct);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initStair(stimulus)

   stimulus.stair = doStaircase('init','upDown', 'nup=1','ndown=2',...
        'initialThreshold', stimulus.initialThreshold, 'initialStepsize',stimulus.initialStepsize, ...
    'minStepsize',stimulus.minStepsize,'maxStepsize',stimulus.maxStepsize,'minThreshold',stimulus.minThreshold,'maxThreshold', stimulus.maxThreshold,...
     'stepRule=pest','dispFig=1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initGaussian(stimulus,myscreen)

global stimulus;
if stimulus.tenbit
  % set maximum color index (for 24 bit color we have 8 bits per channel, so 255)
  maxIndex = 255;

  % get gamma table
  if ~isfield(myscreen,'gammaTable')
    stimulus.linearizedGammaTable = mglGetGammaTable;
    disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
    disp(sprintf('(alaisburr:initGratings) No gamma table found in myscreen. Contrast'));
    disp(sprintf('         displays like this should be run with a valid calibration made by moncalib'));
    disp(sprintf('         for this monitor.'));
    disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
  end
  stimulus.linearizedGammaTable = myscreen.initScreenGammaTable;

  % disppercent(-inf,'Creating gaussian textures');

  % calculate some colors information
  %  number of reserved colors
  stimulus.colors.nReservedColors = size(stimulus.colors.reservedColors,1);
  % number of colors possible for gratings, make sure that we 
  % have an odd number
  stimulus.colors.nGaussianColors = maxIndex+1-stimulus.colors.nReservedColors;
  % if iseven(stimulus.colors.nGaussianColors)
  %   stimulus.colors.nGaussianColors = stimulus.colors.nGaussianColors-1;
  % end

  % min,mid,max index of gaussian colors
  stimulus.colors.minGaussianIndex = maxIndex+1 - stimulus.colors.nGaussianColors;
  stimulus.colors.midGaussianIndex = stimulus.colors.minGaussianIndex + floor(stimulus.colors.nGaussianColors/2);
  stimulus.colors.maxGaussianIndex = maxIndex;
  % number of contrasts we can display (not including 0 contrast)
  stimulus.colors.nDisplayContrasts = floor(stimulus.colors.nGaussianColors-1);

  % set the reserved colors - this gives a convenient value between 0 and 1 to use the reserved colors with
  for i = 1:stimulus.colors.nReservedColors
    stimulus.colors.reservedColor(i) = (i-1)/maxIndex;
  end

  setGammaTableForMaxContrast(stimulus.contrast);
  contrastIndex = getContrastIndex(stimulus.contrast,1);

  % make all the 1D gaussians. We compute all possible contrast values given the
  % range of indexes available to us. The 1st texture is black the nth texture is full
  % contrast for the current gamma setting
  stimulus.gaussian = mglMakeGaussian(7*stimulus.width, 7*stimulus.width, stimulus.width,stimulus.width);
  stimulus.colors.gaussRange = contrastIndex-1;
  % make the grating
  thisGaussian = round(stimulus.colors.gaussRange*stimulus.gaussian + stimulus.colors.minGaussianIndex);
  % create the texture
  stimulus.tex = mglCreateTexture(thisGaussian);

  % get the color value for black (i.e. the number between 0 and 1 that corresponds to the minGaussianIndex)
  stimulus.colors.black = stimulus.colors.minGaussianIndex/maxIndex;
  [1 1 1; 0.3 0.3 0.3; 0 1 0;1 0 0; 0 1 1];
  % get the color values (i.e. reserved color)
  stimulus.colors.white = stimulus.colors.reservedColor(1);
  stimulus.colors.red = stimulus.colors.reservedColor(4);
  stimulus.colors.green = stimulus.colors.reservedColor(3);
  stimulus.colors.grey = stimulus.colors.reservedColor(2);
  stimulus.colors.cyan = stimulus.colors.reservedColor(5);
  
else

  gauss = mglMakeGaussian(7*stimulus.width, 7*stimulus.width, stimulus.width, stimulus.width);
  gaussian = zeros(size(gauss,1), size(gauss,2), 4);
  for i = 1:3
      gaussian(:,:,i) = 255*ones(size(gauss,1), size(gauss,2));
  end
  gaussian(:,:,4) = 255*gauss*stimulus.contrast;
  stimulus.tex = mglCreateTexture(gaussian);
  % get the color value for black (i.e. the number between 0 and 1 that corresponds to the minGaussianIndex)
  stimulus.colors.black = 0;
  % get the color values (i.e. reserved color)
  stimulus.colors.white = 1;
  stimulus.colors.grey = 0.3;
  stimulus.colors.green = [0 1 0];
  stimulus.colors.red = [1 0 0];
  stimulus.colors.cyan = [0 1 1];
end

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
%   cmax = 0.5+maxContrast/2;cmin = 0.5-maxContrast/2;
  cmin = 0;
  cmax = maxContrast;
  luminanceVals = cmin:((cmax-cmin)/(stimulus.colors.nGaussianColors-1)):cmax;

  % replace NaN in gamma tables with zero
  stimulus.linearizedGammaTable.redTable(isnan(stimulus.linearizedGammaTable.redTable)) = 0;
  stimulus.linearizedGammaTable.greenTable(isnan(stimulus.linearizedGammaTable.greenTable)) = 0;
  stimulus.linearizedGammaTable.blueTable(isnan(stimulus.linearizedGammaTable.blueTable)) = 0;

  % now get the linearized range
  redLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,luminanceVals,'linear');
  greenLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,luminanceVals,'linear');
  blueLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,luminanceVals,'linear');
  
  % add these values to the table
  gammaTable((stimulus.colors.minGaussianIndex:stimulus.colors.maxGaussianIndex)+1,:)=[redLinearized;greenLinearized;blueLinearized]';
else
  % if we are asked for 0 contrast then simply set all the values to BLACK
  gammaTable((stimulus.colors.minGaussianIndex:stimulus.colors.maxGaussianIndex)+1,1)=interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,0,'linear');
  gammaTable((stimulus.colors.minGaussianIndex:stimulus.colors.maxGaussianIndex)+1,2)=interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,0,'linear');
  gammaTable((stimulus.colors.minGaussianIndex:stimulus.colors.maxGaussianIndex)+1,3)=interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,0,'linear');
end

% set the gamma table
mglSetGammaTable(gammaTable);

% keep the gamma table
stimulus.gammaTable = gammaTable;

% remember what the current maximum contrast is that we can display
stimulus.currentMaxContrast = maxContrast;

%%%%%%%%%%%%%%%%%%%%%%%%
%    testGammaTable    %
%%%%%%%%%%%%%%%%%%%%%%%%
function tf = testGammaTable(stimulus, myscreen)

% clear screen to black
mglClearScreen(stimulus.colors.nReservedColors/255);

% setup text
mglTextSet('Helvetica',32,0.8,0,0,0);

% setup rect dimensions
rectHeight = myscreen.imageHeight/2;
rectY = rectHeight/2;
rectWidth = myscreen.imageWidth/stimulus.colors.nReservedColors;

% make nReservedColors rectangles with the reserved colors
for iColor = 1:stimulus.colors.nReservedColors
  % get color index
  colorIndex = stimulus.colors.reservedColor(iColor);
  % make color square
  rectX = -(myscreen.imageWidth/2) + rectWidth * (iColor-1)+rectWidth/2;
  mglFillRect(rectX,rectY,[rectWidth rectHeight],colorIndex);
  % draw text
  mglTextDraw(sprintf('Color: %i',iColor),[rectX,rectY]);
end

disp(sprintf('(alaisburr:testGammaTable) Top row should be reserved colors'));

% setup rect dimensions
rectHeight = myscreen.imageHeight/2;
rectY = -rectHeight/2;
rectWidth = myscreen.imageWidth/(256-stimulus.colors.nReservedColors);

% make nReservedColors rectangles with the reserved colors
for iColor = stimulus.colors.nReservedColors:255
  % get color index
  colorIndex = iColor/255;
  % make color square
  rectX = -(myscreen.imageWidth/2) + rectWidth * (iColor-1)+rectWidth/2;
  mglFillRect(rectX,rectY,[rectWidth rectHeight],colorIndex);
end

disp(sprintf('(alaisburr:testGammaTable) Bottom row should be stimulus colors'));
mglFlush;

tf = askuser('Continue');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sound
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initClick(stimulus,myscreen)
% sampling frequency (samples per sec)
duration = stimulus.clickDur;
fs = stimulus.samplesPerSecond-1;
t = 0:1/fs:duration;
% frequency of signal in hz
% hz = 440;
% amplitude = 0.5;
% stimulus.wav = amplitude * sin(2*pi*hz*t);

% wav = 0.5 * randn(1,length(t));
% stimulus.wav = wav;

wav = 0.5 * randn(1,length(t));

fc = 2000; % cutoff frequency
% 5th order Butterworth filter
% [b,a] = butter(5, fc/(stimulus.tone.samplesPerSecond/2));
a = [1  -4.07876493416512 6.72527084144657  -5.59474636818042 2.34559680959441  -0.396133028715511];
b = [3.82287493728255e-05 0.000191143746864127  0.000382287493728255  0.000382287493728255  0.000191143746864127  3.82287493728255e-05];
wav = randn(1,length(t));
wavFiltered = filter(b,a,wav);
stimulus.wav = wavFiltered;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to update the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = createITD(stimulus,theta)
% theta: angle to sound source (-90 to left and +90 to right)
% %%% interaural time difference
% % radius of head in meter (assuming spherical)
% r = 0.0875;
% % speed of sound at room temperature (m/s)
% c = 346;
% fs = stimulus.samplesPerSecond - 1;
% % distance from monitor
% d = stimulus.displayDistance;
% % for low frequency sound
% % td = (2*sind(theta)) * r / c;
% % left - right
% td = (sqrt((d*tand(theta)-r).^2 + (d+r)^2) - sqrt((d*tand(theta)+r).^2 + (d+r)^2))./c;
% td_a = 0:1/fs:abs(td);
% clear waveform s
% if td > 0
%     waveform(1,:) = [stimulus.wav, zeros(1,length(td_a))];
%     waveform(2,:) = [zeros(1,length(td_a)), stimulus.wav];
% elseif td < 0
%     waveform(2,:) = [stimulus.wav, zeros(1,length(td_a))];
%     waveform(1,:) = [zeros(1,length(td_a)), stimulus.wav];
% else
%     waveform(1,:) = stimulus.wav;
%     waveform(2,:) = stimulus.wav;
% end
len = length(stimulus.wav);
% load impulse response
elev = 0;
if mod(theta, 5) == 0
    if theta >= 0
        temp = readhrtf(elev,theta,'L');
        IR_L = temp(1,:);
        IR_R = temp(2,:);
    else
        temp = readhrtf(elev,-theta,'L');
        IR_L = temp(2,:);
        IR_R = temp(1,:);
    end
    % fft
    yl = fft(IR_L, len);
    yr = fft(IR_R, len);
else
    if theta >= 0
        temp{1} = readhrtf(elev, floor(theta/5) * 5, 'L');
        temp{2} = readhrtf(elev, ceil(theta/5) * 5, 'L');
        for i = 1:2
            IR_L{i} = temp{i}(1,:);
            IR_R{i} = temp{i}(2,:);
        end
    else
        theta = -theta;
        temp{1} = readhrtf(elev, floor(theta/5) * 5, 'L');
        temp{2} = readhrtf(elev, ceil(theta/5) * 5, 'L');
        for i = 1:2
            IR_L{i} = temp{i}(2,:);
            IR_R{i} = temp{i}(1,:);
        end
    end
    for i = 1:2
        YL{i} = fft(IR_L{i},len);
        YR{i} = fft(IR_R{i},len);
    end
    
    yl = ((ceil(theta/5)*5 - theta)/5 * YL{1}) + ((theta - floor(theta/5)*5)/5 * YL{2});
    yr = ((ceil(theta/5)*5 - theta)/5 * YR{1}) + ((theta - floor(theta/5)*5)/5 * YR{2});
end
yw = fft(stimulus.wav,len);
if size(yw,1) ~= size(yl,1);
    yl = yl';
    yr = yr';
end

leftFunc = yw .* yl;
rightFunc = yw .* yr;

leftwav = ifft(leftFunc, len);
rightwav = ifft(rightFunc,len);
clear waveform
if size(leftwav,1) > size(leftwav,2)
    waveform(1,:) = leftwav';
    waveform(2,:) = rightwav';
else
    waveform(1,:) = leftwav;
    waveform(2,:) = rightwav;
end
s = mglInstallSound(waveform, stimulus.samplesPerSecond);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display psychometric functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispPsychometric(task,stimulus)
posDiff = task.randVars.diff;
allresp = task.randVars.resp;
ind = ~isnan(posDiff)&~isnan(allresp);
posDiff =  posDiff(ind);
allresp =  allresp(ind);
testVal = unique(posDiff);

if ~stimulus.bimodal

n = zeros(1,length(testVal)); k = zeros(1,length(testVal));
isLeft = cell(1); resp= cell(1); whichInt = cell(1);

for i = 1:length(testVal)
    resp{i} = allresp(posDiff == testVal(i));
    whichInt{i} = task.randVars.centerInt(posDiff == testVal(i));
    for j = 1:length(resp{i})
      switch whichInt{i}(j)
        case 1
          if resp{i}(j) == 2
            % If probe seen Left
            isLeft{i}(j) = 1;
          elseif resp{i}(j) == 1
            % probe seen Right
            isLeft{i}(j) = 0;
          end
        case 2
          if resp{i}(j) == 1
            % probe seen Left
            isLeft{i}(j) = 1;
          elseif resp{i}(j) == 2
            % probe seen Right
            isLeft{i}(j) = 0;
          end
        end
     end
    n(i) = sum(resp{i} == 1 | resp{i}==2);
    % k(i) = sum(resp{i} == 2);
    k(i) = sum(isLeft{i} == 1);
end
percent = k./n;

mlrSmartfig('alaisburr','reuse');clf;
h = plot(testVal, percent, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
ylabel('Proportion of trials probe seen "left"');
xlabel('Displacement of probe (degs)');
axis([-20 20 0 1]); box off;

else

displ = task.randVars.displ;
displ = displ(ind);
testVal = unique(posDiff);
displacement = unique(displ);
whichint = task.randVars.centerInt;
whichint = whichint(ind);
isLeft = cell(1);
mlrSmartfig('alaisburr','reuse');clf;
    brewer = brewermap(5, 'Set1');
    for d = 1:length(displacement)
      posdiffConflict{d} = posDiff(displ==displacement(d));
      whichIntConflict{d} = whichint(displ==displacement(d));
      respConflict{d} = allresp(displ==displacement(d));

      for j = 1:length(respConflict{d})
        switch whichIntConflict{d}(j)
               case 1
                if respConflict{d}(j) == 2
                  isLeft{d}(j) = 1;
                else
                  isLeft{d}(j) = 0;
                end
               case 2
                if respConflict{d}(j) == 1
                  isLeft{d}(j) = 1;
                else
                  isLeft{d}(j) = 0;
                end
         end
      end
    nTotal{d} = arrayfun(@(x) sum(posdiffConflict{d}==x), testVal);
    nLeft{d} = arrayfun(@(x) sum(posdiffConflict{d}==x & isLeft{d}==1), testVal); 
    pLeft = nLeft{d}./nTotal{d};
          hold on;
      myerrorbar(testVal, pLeft, 'MarkerSize',6, 'Color', brewer(d,:), 'Symbol',getsymbol(d));

        if d == length(displacement)
        title(sprintf('bimodal: Width=%d \n', stimulus.width));
        ylabel('Proportion probe "left"');
        xlabel('Displacement of probe (degs)');
        axis([-20 20 0 1]); box off;
        
        l = cellfun(@num2str, num2cell(displacement),'UniformOutput',0);
        mylegend(l, {{getsymbol(1), brewer(1,:)},{getsymbol(2), brewer(2,:)},{getsymbol(3), brewer(3,:)},{getsymbol(4), brewer(4,:)},{getsymbol(5), brewer(5,:)}});

        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    initBackgroundNoise    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initBackgroundNoise(stimulus,myscreen)

% get the gaussian full size of the screen
[stimulus.background.gaussian stimulus.background.x stimulus.background.y] = mglMakeGaussian(myscreen.imageWidth, myscreen.imageHeight, stimulus.width,stimulus.width);

% now make sure that dimensions are odd numbers of pixels
oddWidth = 2*floor(myscreen.screenWidth/2)+1;
oddHeight = 2*floor(myscreen.screenHeight/2)+1;

% resize everything to odd
stimulus.background.gaussian = stimulus.background.gaussian(1:oddHeight,1:oddWidth);
stimulus.background.x = stimulus.background.x(1:oddHeight,1:oddWidth);
stimulus.background.y = stimulus.background.y(1:oddHeight,1:oddWidth);

% get the fourier transform
stimulus.background.gaussianTransform = getHalfFourier(stimulus.background.gaussian);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    setBackgroundNoise    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = setBackgroundNoise(stimulus, myscreen, task, SNR, backgroundFreq)

% figure out how many frames we should have
% for the length of the trial so that we can
% randomly show at least one per frame
trialTime = sum(task{end}{end}.segmax);

% make a few times more than necessary so that there will be more randomization
numBackgrounds = round(trialTime * backgroundFreq * 2);

% set the noise maximum (i.e. maximum luminance of noise)
% set it so that we can achieve the maximu SNR that is asked for.
stimulus.background.noiseMax = 1 / (max(SNR) + 1);
% now set all the various SNR levels
stimulus.background.sigMax = SNR * stimulus.background.noiseMax;

disppercent(-inf,'(alaisburr:setBackgroundNoise) Precomputing background noise images');

% delete old textures
if isfield(stimulus,'backTexture') && ~isempty(stimulus.backTexture)
  for iBackground = 1:numBackgrounds
    mglDeleteTexture(stimulus.backTexture(iBackground));
  end
end

for iBackground = 1:numBackgrounds
  % randomize phase and reconstruct
  stimulus.background.gaussianTransform.phase = (rand(1,stimulus.background.gaussianTransform.n)*2*pi - pi);
  im = reconstructFromHalfFourier(stimulus.background.gaussianTransform);

  % scale from 0 to noise max
  maxIm = max(im(:));
  minIm = min(im(:));
  stimulus.background.im(iBackground,:,:) = stimulus.background.noiseMax * (im - minIm) / (maxIm-minIm);

  % make into texture
  stimulus.backTexture(iBackground) = mglCreateTexture(round(stimulus.colors.gaussRange*squeeze(stimulus.background.im(iBackground,:,:)) + stimulus.colors.minGaussianIndex));
  
  % update disppercent
  disppercent(iBackground/numBackgrounds);
end
disppercent(inf);

% set how many textures we have
stimulus.background.n = numBackgrounds;

% compute how long to show each frame for
stimulus.background.frameTime = 1/backgroundFreq;
stimulus.backgroundFreq = backgroundFreq;


%%%%%%%%%%%%%%%%%%%%%%%%%
%    setStimuliOnNoise  %
%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = setStimulusOnBackground(stimulus, xPos, yPos, stimulusNum, backgroundNum, SNR)

% get signal max for this trial
sigMax = stimulus.background.sigMax(find(stimulus.SNR == SNR));

% make the gaussian at the xPos, yPos and scale to signal max
im = squeeze(stimulus.background.im(backgroundNum,:,:)) + sigMax * exp(-(((stimulus.background.x-xPos).^2)/(2*(stimulus.width^2))+((stimulus.background.y-yPos).^2)/(2*(stimulus.width^2))));

% delete any existing texture
if isfield(stimulus,'stimTexture') && (length(stimulus.stimTexture) >= stimulusNum)
  mglDeleteTexture(stimulus.stimTexture(stimulusNum));
end

% scale and set texture
stimulus.stimTexture(stimulusNum) = mglCreateTexture(round(stimulus.colors.gaussRange*im + stimulus.colors.minGaussianIndex));


