function myscreen = spatialBisection(varargin)

clear global stimulus
mglEatKeys('12`');
global stimulus

% get arguments
high = 0; low = 0; med = 0; tenbit = 1; practice = 0; auditoryTrain = 0; visualTrain=0;
getArgs(varargin,{'high=0','low=0','med=0','tenbit=1','practice=0','auditoryTrain=0','visualTrain=0'},'verbose=1');

if high
    stimulus.gaussian.diameter = 4;
    stimulus.gaussian.contrast = 0.005;
elseif low
    stimulus.gaussian.diameter = 240;
    stimulus.gaussian.contrast = 0.002;
elseif med
    stimulus.gaussian.diameter = 32;
    stimulus.gaussian.contrast = 0.002;
else 
    if auditoryTrain
        tenbit =0;
        stimulus.gaussian.diameter = nan;
    else
        return
    end
end
stimulus.high = high;
stimulus.low = low;
stimulus.med = med;
stimulus.tenbit = tenbit;
stimulus.practice = practice;
stimulus.auditoryTrain = auditoryTrain;
stimulus.visualTrain = visualTrain;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stimulus.gaussian.diameter = 14;
stimulus.gaussian.sd = stimulus.gaussian.diameter/7;
stimulus.gaussian.duration = .015;% .025;%1/60; % one(or two) frame 
    if ~stimulus.tenbit
        stimulus.gaussian.contrast = .1;
    end
stimulus.colors.reservedColors = [1 1 1; 0.3 0.3 0.3; 0 1 0;1 0 0; 0 1 1];

stimulus.tone.samplesPerSecond = 44100;
stimulus.tone.duration = .0015;

stimulus.pos1 = -17.5;
stimulus.pos3 = 17.5;
stimulus.midPoint = (stimulus.pos1 + stimulus.pos3)/2;

% stimulus.pos1 = -7.5;
% stimulus.pos3 = stimulus.pos1+30;
% stimulus.midPoint = (stimulus.pos1 + stimulus.pos3)/2; %7.5 deg from center
stimulus.delta=2.5;

% fixation cross
stimulus.fixWidth = 1;
stimulus.fixColor = [1 1 1];

stimulus.initOffset = 2;
stimulus.initOffsetSd = 2.5;
if stimulus.auditoryTrain || stimulus.visualTrain
    stimulus.initialThreshold = 10;
    stimulus.initialStepsize = 2.5;
    stimulus.minThreshold = 0;
    stimulus.maxThreshold = 15;
    stimulus.minStepsize = 0.75;
    stimulus.maxStepsize = 5;
end

% initalize the screen
myscreen = initScreen;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%% 
task{1}{1}.waitForBacktick = 1;
stimDur = stimulus.gaussian.duration;
task{1}{1}.segmin = [1 0.5 stimDur 0.5-stimDur stimDur 0.5-stimDur stimDur 1.5 1];
task{1}{1}.segmax = task{1}{1}.segmin;
task{1}{1}.getResponse = [0 0 0 0 0 0 0 1 0];

% parameters & randomization
if stimulus.practice
    task{1}{1}.parameter.condition = {'vision','auditory','noOffset'};
    task{1}{1}.numTrials = 10*length(task{1}{1}.parameter.condition);
    task{1}{1}.randVars.uniform.closeTo = [1 3];
elseif stimulus.auditoryTrain
%     task{1}{1}.numTrials = 30;
    task{1}{1}.parameter.condition = {'auditory'};
%     task{1}{1}.parameter.offset = [7.5 10];
    task{1}{1}.parameter.closeTo = [1 3];
    % figure;
    % xlabel('trial'); ylabel('correct'); 
    % xaxis([0 30]); yaxis([-1 2]);
elseif stimulus.visualTrain
	task{1}{1}.parameter.condition = {'vision'};
    task{1}{1}.parameter.closeTo = [1 3];
else
    task{1}{1}.parameter.condition = {'vision','auditory','noOffset','posOffset','negOffset'};
	task{1}{1}.numTrials = 25*length(task{1}{1}.parameter.condition);
    task{1}{1}.randVars.uniform.closeTo = [1 3];
end

task{1}{1}.random = 1;

task{1}{1}.randVars.calculated.resp = nan;
task{1}{1}.randVars.calculated.correct = nan;
task{1}{1}.randVars.calculated.probeOffset = nan; % from midpoint
task{1}{1}.randVars.calculated.noise = nan;
task{1}{1}.randVars.calculated.rt = nan;
task{1}{1}.randVars.calculated.condNum = nan;
task{1}{1}.randVars.calculated.pos1 = nan;
% task{1}{1}.randVars.calculated.hz = nan;
 
% initialize the task
for phaseNum = 1:length(task{1})
  [task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end
 
% init the stimulus
myscreen = initStimulus('stimulus',myscreen);

% init the staircase
stimulus = initStair(stimulus);
% to initialize the stimulus for your experiment.
stimulus = initGaussian(stimulus,myscreen);
stimulus = initClick(stimulus,task);

mglWaitSecs(1);
mglClearScreen(stimulus.colors.black);
mglTextSet([],32,stimulus.colors.white);
mglTextDraw('Press ` key to start when you are ready',[0 0]);
mglFlush;
mglWaitSecs(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% if stimulus.disp
dispPsychometric(task{1}{1},stimulus);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)
global stimulus
if task.thistrial.thisseg == 1
	stimulus.fixColor = stimulus.colors.white;
	task.thistrial.condNum = find(strcmp(char(task.thistrial.condition),{'vision','auditory','noOffset','posOffset','negOffset'}));
    if stimulus.auditoryTrain || stimulus.visualTrain
    	[testValue, stimulus.stair] = doStaircase('testValue', stimulus.stair);
%         testValue = task.thistrial.offset;
        task.thistrial.noise = 0;

    else
        % get Test Value
        [testValue, stimulus.stair{task.thistrial.condNum}] = doStaircase('testValue', stimulus.stair{task.thistrial.condNum});
        if testValue > 10.5
            testValue = 10.5;
        end
        task.thistrial.noise = 4 * randn(1); % random number from a gaussian distribution with a std of 4 deg
        while (stimulus.midPoint + (testValue + task.thistrial.noise) <= stimulus.pos1+stimulus.delta ) || (stimulus.midPoint + (testValue + task.thistrial.noise) >= stimulus.pos3-stimulus.delta)
            task.thistrial.noise = 4 * randn(1);
        end
    end

	if task.thistrial.closeTo == 1
		sign = -1;
	else
		sign = 1;
    end
		
	task.thistrial.probeOffset = sign * (testValue + task.thistrial.noise);
 	switch char(task.thistrial.condition)
		case 'vision'
			task.thistrial.xposV = [stimulus.pos1 stimulus.midPoint+task.thistrial.probeOffset stimulus.pos3];
			task.thistrial.xposA = [nan nan nan];
		case 'auditory'
			task.thistrial.xposV = [nan nan nan];
			task.thistrial.xposA = [stimulus.pos1 stimulus.midPoint+task.thistrial.probeOffset stimulus.pos3];
		case 'noOffset'
			task.thistrial.xposV = [stimulus.pos1 stimulus.midPoint+task.thistrial.probeOffset stimulus.pos3];
			task.thistrial.xposA = [stimulus.pos1 stimulus.midPoint+task.thistrial.probeOffset stimulus.pos3];
		case 'posOffset' % 1st/3rd tone shifted LEFTward & disc shifted RIGHTward
			task.thistrial.xposV = [stimulus.pos1+stimulus.delta stimulus.midPoint+task.thistrial.probeOffset stimulus.pos3+stimulus.delta];
			task.thistrial.xposA = [stimulus.pos1-stimulus.delta stimulus.midPoint+task.thistrial.probeOffset stimulus.pos3-stimulus.delta];
		case 'negOffset' % 1st/3rd tone shifted LEFTward & disc shifted RIGHTward
			task.thistrial.xposV = [stimulus.pos1-stimulus.delta stimulus.midPoint+task.thistrial.probeOffset stimulus.pos3-stimulus.delta];
			task.thistrial.xposA = [stimulus.pos1+stimulus.delta stimulus.midPoint+task.thistrial.probeOffset stimulus.pos3+stimulus.delta];
	end

	if task.thistrial.condNum ~= 1 % if not vision condition
		% stimulus = initClick(stimulus,task);
        for int = 1:3
            stimulus.sound(int) = createITD(stimulus,task.thistrial.xposA(int));
        end
        % task.thistrial.hz = stimulus.thisHz;
    end

% elseif task.thistrial.thisseg == 9
% 	stimulus.fixColor = stimulus.colors.darkgrey;

elseif task.thistrial.thisseg == 8
	stimulus.fixColor = stimulus.colors.darkgrey;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)
global stimulus
mglClearScreen(stimulus.colors.black);
mglFixationCross(stimulus.fixWidth,1.5,stimulus.fixColor);

if any(task.thistrial.thisseg == [3 5 7])
  switch char(task.thistrial.condition)
	case 'vision'
		mglBltTexture(stimulus.tex, [task.thistrial.xposV(floor(task.thistrial.thisseg/2)), 0]);
	case 'auditory'
		mglPlaySound(stimulus.sound(floor(task.thistrial.thisseg/2)));
	otherwise
		mglPlaySound(stimulus.sound(floor(task.thistrial.thisseg/2)));
		mglBltTexture(stimulus.tex, [task.thistrial.xposV(floor(task.thistrial.thisseg/2)), 0]);

  end		
% elseif any(task.thistrial.thisseg == [1 8 9])
% 	mglFixationCross(stimulus.fixWidth,1.5,stimulus.fixColor);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)
global stimulus
% here, we just check whether this is the first time we got a response
if ~task.thistrial.gotResponse

	if (task.thistrial.probeOffset < 0 && task.thistrial.whichButton == 1) || ... % closer to the first
		(task.thistrial.probeOffset > 0 && task.thistrial.whichButton == 2)  % closer to the third 
		% correct
		task.thistrial.correct = 1;
		% if any(task.thistrial.condNum == [1 2 3])	
		% % feedback
		if stimulus.practice || stimulus.auditoryTrain || stimulus.visualTrain
			stimulus.fixColor = stimulus.colors.green;
		end
		disp(sprintf('(spatialBisection) %i:%s offset %0.3f resp %i correct', ...
            task.trialnum, char(task.thistrial.condition), task.thistrial.probeOffset, task.thistrial.whichButton))

	else
		% incorrect
		task.thistrial.correct = 0;
		% if any(task.thistrial.condNum == [1 2 3])
		if stimulus.practice || stimulus.auditoryTrain || stimulus.visualTrain
			stimulus.fixColor = stimulus.colors.red;
		end
		disp(sprintf('(spatialBisection) %i:%s offset %0.3f resp %i incorrect', ...
            task.trialnum, char(task.thistrial.condition), task.thistrial.probeOffset, task.thistrial.whichButton))
	end
	if ~(stimulus.practice || stimulus.auditoryTrain || stimulus.visualTrain)
		stimulus.fixColor = stimulus.colors.cyan;
    end
    if ~(stimulus.auditoryTrain||stimulus.visualTrain)
        
        if strcmp(char(task.thistrial.condition),'vision') || strcmp(char(task.thistrial.condition),'auditory')
            stimulus.stair{task.thistrial.condNum} = doStaircase('update', stimulus.stair{task.thistrial.condNum}, task.thistrial.correct, ...
            abs(task.thistrial.probeOffset));        
        else
          randomnumbers = [1 1 1 1 1 1 1 1 1 0];
          thisrandom = randomnumbers(randi(length(randomnumbers),1));
           stimulus.stair{task.thistrial.condNum} = doStaircase('update', stimulus.stair{task.thistrial.condNum}, thisrandom, ...
            abs(task.thistrial.probeOffset));
        end
    else
        stimulus.stair = doStaircase('update', stimulus.stair, task.thistrial.correct);
    end
	task.thistrial.resp = task.thistrial.whichButton;
    task.thistrial.rt = task.thistrial.reactionTime;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initStair(stimulus)
if stimulus.auditoryTrain || stimulus.visualTrain
	 stimulus.stair = doStaircase('init','upDown', 'nup=1','ndown=2',...
        'initialThreshold', stimulus.initialThreshold, 'initialStepsize',stimulus.initialStepsize, ...
    'minStepsize',stimulus.minStepsize,'maxStepsize',stimulus.maxStepsize,'minThreshold',stimulus.minThreshold,'maxThreshold', stimulus.maxThreshold,...
     'stepRule=pest','dispFig=1');
else

	condNames = {'vision','auditory','noOffset','posOffset','negOffset'};
for cond = 1:5
	stimulus.stair{cond} = doStaircase('init','quest', 'initialThreshold', stimulus.initOffset, 'initialThresholdSd', stimulus.initOffsetSd, ...
		'pThreshold', 0.75,'dispFig=1','subplotRows=5','subplotCols=1','subplotNum',cond,'subplotName',condNames{cond});
end
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

setGammaTableForMaxContrast(stimulus.gaussian.contrast);
contrastIndex = getContrastIndex(stimulus.gaussian.contrast,1);

% make all the 1D gaussians. We compute all possible contrast values given the
% range of indexes available to us. The 1st texture is black the nth texture is full
% contrast for the current gamma setting
gaussian = mglMakeGaussian(stimulus.gaussian.diameter, stimulus.gaussian.diameter, stimulus.gaussian.sd,stimulus.gaussian.sd);
iContrast = contrastIndex-1;
% for iContrast = 0:stimulus.colors.nDisplayContrasts
  % disppercent(iContrast/stimulus.colors.nDisplayContrasts);
  % if myscreen.userHitEsc,mglClose;keyboard,end
  % make the grating
  thisGaussian = round(iContrast*gaussian + stimulus.colors.minGaussianIndex);
  % create the texture
  % stimulus.tex(iContrast+1) = mglCreateTexture(thisGaussian);
  stimulus.tex = mglCreateTexture(thisGaussian);
% end
% disppercent(inf);

% get the color value for black (i.e. the number between 0 and 1 that corresponds to the minGaussianIndex)
stimulus.colors.black = stimulus.colors.minGaussianIndex/maxIndex;
% get the color values (i.e. reserved color)
stimulus.colors.white = stimulus.colors.reservedColor(1);
stimulus.colors.darkgrey = stimulus.colors.reservedColor(2);
stimulus.colors.green = stimulus.colors.reservedColor(3);
stimulus.colors.red = stimulus.colors.reservedColor(4);
stimulus.colors.cyan = stimulus.colors.reservedColor(5);

% % compute the guassian
% gauss = mglMakeGaussian(stimulus.width,stimulus.width, stimulus.width/8,stimulus.width/8);

% gaussian = zeros(size(gauss,1), size(gauss,2), 4);
% for i = 1:3
%     gaussian(:,:,i) = 255*ones(size(gauss,1), size(gauss,2));
% end
%     gaussian(:,:,4) = 255*gauss*stimulus.contrast;
    
% %create texture
% stimulus.tex = mglCreateTexture(gaussian);

else
	gauss = mglMakeGaussian(stimulus.gaussian.diameter, stimulus.gaussian.diameter, stimulus.gaussian.sd,stimulus.gaussian.sd);
	gaussian = zeros(size(gauss,1), size(gauss,2), 4);
	for i = 1:3
    	gaussian(:,:,i) = 255*ones(size(gauss,1), size(gauss,2));
	end
    gaussian(:,:,4) = 255*gauss*stimulus.gaussian.contrast;
    stimulus.tex = mglCreateTexture(gaussian);
    % get the color value for black (i.e. the number between 0 and 1 that corresponds to the minGaussianIndex)
stimulus.colors.black = 0;
% get the color values (i.e. reserved color)
stimulus.colors.white = 1;
stimulus.colors.darkgrey = 0.3;
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

% remember what the current maximum contrast is that we can display
stimulus.currentMaxContrast = maxContrast;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sound
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initClick(stimulus,myscreen)
% sampling frequency (samples per sec)
fs = stimulus.tone.samplesPerSecond-1;
t = 0:1/fs:stimulus.tone.duration;

wav = 0.5 * randn(1,length(t));

fc = 2000; % cutoff frequency
% 5th order Butterworth filter
% [b,a] = butter(5, fc/(stimulus.tone.samplesPerSecond/2));
a = [1	-4.07876493416512	6.72527084144657	-5.59474636818042	2.34559680959441	-0.396133028715511];
b = [3.82287493728255e-05	0.000191143746864127	0.000382287493728255	0.000382287493728255	0.000191143746864127	3.82287493728255e-05];
wav = randn(1,length(t));
wavFiltered = filter(b,a,wav);
stimulus.wav = wavFiltered;

% stimulus.wav = wav;

% stimulus.thisHz = stimulus.tone.hz + (stimulus.tone.hz/10)*randn;

% wav = sin(2*pi*stimulus.thisHz*t);
% if isodd(length(t))
% 	mu = t((1+length(t))/2);
% else
% 	mu = t(round((1+length(t))/2));
% end
% sigma = stimulus.tone.duration;
% envelope = normpdf(t,mu,sigma);
% %normalize
% envelope = envelope / max(envelope);
% waveform = wav .* envelope;
% stimulus.wav = waveform;

% frequency of signal in hz
% stimulus.wav = sin(2*pi*stimulus.tone.hz*t);
stimulus.sound(1) = mglInstallSound(stimulus.wav, stimulus.tone.samplesPerSecond);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to update the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = createITD(stimulus,theta)
% theta: angle to sound source (-90 to left and +90 to right)

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
s = mglInstallSound(waveform, stimulus.tone.samplesPerSecond);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display psychometric functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispPsychometric(task,stimulus)
uni = figure; 
if ~(stimulus.auditoryTrain || stimulus.visualTrain)
bi = figure;
condNames = {'vision','auditory','noOffset','posOffset','negOffset'};
nCon = 5;
else
	nCon = 2;
	
end

% condition = {};
% for b = 1:length(task.block)
%   condition = [condition, task.block(b).parameter.condition];
% end
condition = task.randVars.condNum;

for cond = 1:nCon
  thisTrials = find(condition == cond);%find(strcmp(condition, condNames{cond}));
  thisOffset = task.randVars.probeOffset(thisTrials);
  thisResp = task.randVars.resp(thisTrials);
  isThird = (thisResp == 2);

	binCenter = -14.5:1:14.5; space = max(diff(binCenter));
  for b = 1:length(binCenter)
    switch b
      case 1
        binnedTrial{b} = find(thisOffset<= min(binCenter)+space/2); % -200 ~ -150 ms
      case length(binCenter)
        binnedTrial{b} = find(thisOffset > max(binCenter)-space/2);
      otherwise
        binnedTrial{b} = find((thisOffset > min(binCenter) + space*(b-2)) & (thisOffset <= min(binCenter) + space*(b-1)));
    end
  end
  % binCenter = -12.5:1:12.5;%-0.175:.05:0.175;
  nTotalBin = cellfun(@(x) length(x), binnedTrial);
  binnedTrial(nTotalBin==0) = [];
  binCenter(nTotalBin==0) = [];
  nTotalBin(nTotalBin==0) = [];

  nThirdBin = cellfun(@(x) sum(isThird(x)), binnedTrial);
  pThirdBin = nThirdBin./nTotalBin;

  % xs = -0.2:.001:0.2;
  % fit{cond} = fitCG(binCenter, nTotalBin, nLateBin);

  % mu = fit{cond}.fitparams(1);
  % sigma = fit{cond}.fitparams(2);
  % lambda = fit{cond}.fitparams(3);
  % cgfunc = fit{cond}.cgfunc;
  % fitVal = arrayfun(@(x) cgfunc(mu,sigma,lambda,x), xs);
  % pfitValBin = arrayfun(@(x) cgfunc(mu,sigma,lambda,x), binCenter);
  % error = sqrt(pfitValBin.*(1-pfitValBin)./nTotalBin);

  switch cond
    case {1,2} 
      figure(uni);
      hold on;
      myerrorbar(binCenter, pThirdBin,'Symbol',getsymbol(cond),'Color',getcolor(cond), 'MarkerSize',6);
      % myerrorbar(binCenter+(rand-0.5)*0.002, pLateBin,'yError',error, 'Symbol',getsymbol(cond),'Color',getcolor(cond), 'MarkerSize',6);
      % plot(xs, fitVal, 'Color',getcolor(cond),'LineWidth',0.8);
    otherwise
    	figure(bi);
     hold on;
      myerrorbar(binCenter, pThirdBin,'Symbol',getsymbol(cond),'Color',getcolor(cond), 'MarkerSize',6);
      % myerrorbar(binCenter+(rand-0.5)*0.002, pLateBin,'yError',error, 'Symbol',getsymbol(cond),'Color',getcolor(cond), 'MarkerSize',6);
      % plot(xs, fitVal, 'Color',getcolor(cond),'LineWidth',0.8);
  end

end

if ~(stimulus.auditoryTrain || stimulus.visualTrain)
figure(bi); mylegend({'No offset','Pos offset','Neg offset'}, {{getcolor(3)},{ getcolor(4)},{ getcolor(5)}});
end
figure(uni); mylegend({'vision','auditory'},{{getcolor(1)},{getcolor(2)}});

% 
% function dispPerformance(task)
% hold on;
% if task.thistrial.correct
% plot(task.trialnum, task.thistrial.correct, 'go', 'markerFaceColor',[0 1 0]);
% else
%     plot(task.trialnum, task.thistrial.correct, 'ro', 'markerFaceColor',[1 0 0]);
% end
% drawnow;
% if any(task.trialnum == [10 20 30])
%     disp(sprintf('(spatialBisection) Percent Correct of Last 10 Trials: %0.2f%', sum(task.randVars.correct(task.trialnum-9:task.trialnum))/task.trialnum*100));
% end