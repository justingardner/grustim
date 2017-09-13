% alaisburr.m
%
%      usage: myscreen = alaisburr()
%         by: minyoung lee
%       date: 
%    purpose: replication of Alais & Burr, 2004

function myscreen = alaisburr(varargin)
 
clear global stimulus
mglEatKeys('12`');
global stimulus
 
% get arguments
width = 32; visual = 0; auditory = 0; bimodal = 0; disp = 0;
getArgs(varargin,{'width=32','visual=0','auditory=0','bimodal=0','disp=0'},'verbose=1');

if sum([visual,auditory,bimodal]) > 1
    warning('(alaisburr) More than one task type detected.');
    return
elseif sum([visual,auditory,bimodal]) == 0
    warning('(alaisburr) Task type unspecified. Running visual task...')
    return
end
if visual
    stimulus.task = 1;
elseif auditory
    stimulus.task = 2;
else
    stimulus.task = 3;
end
stimulus.disp = disp;
%%%%%%%%%%%%%%%%%%%%%%%%%
stimulus.width = width;
stimulus.stimDur = .015; % 15ms
stimulus.gaussainDur = .015; % 15ms
stimulus.clickDur = 0.0015; % 1.5ms
stimulus.samplesPerSecond = 44100;
stimulus.ISI = .500; % 500ms
stimulus.contrast = .05; % 10% contrast
stimulus.interval = [2 4];
% fixation cross
stimulus.fixWidth = 1;
stimulus.fixColor = [1 1 1];
stimulus.colors.reservedColors = [0 0 0; 1 1 1;0.5 0.5 0.5; 0.5 0 0; 1 0 0; 0 0.5 0; 0 1 0];

screenParams = mglGetScreenParams;
stimulus.displayDistance = screenParams{1}.displayDistance*.01;

% initalize the screen
% myscreen.background = 0;  %black
myscreen = initScreen;

%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%% 
task{1}{1}.waitForBacktick = 1;
% trial: Fixation + Stim1 (.015s/.0015s) + ISI (.5s) + Stim2 (.015s/.0015s) + Resp + ITI
task{1}{1}.segmin = [1 stimulus.stimDur stimulus.ISI stimulus.stimDur 1.5 1];
task{1}{1}.segmax = [1 stimulus.stimDur stimulus.ISI stimulus.stimDur 1.5 1];
task{1}{1}.getResponse = [0 0 0 0 1 0];
if stimulus.task ~= 3
  task{1}{1}.numBlocks = 5;
else
  task{1}{1}.numBlocks = 1;
end
% parameters & randomization
task{1}{1}.parameter.centerWhich = [1 2]; % centered in which interval
task{1}{1}.random = 1;
task{1}{1}.parameter.posDiff = [-20 -15 -10 -5 -2.5 -1.25 0 1.25 2.5 5 10 15 20]; 
task{1}{1}.parameter.displacement = [-5 -2.5 0 2.5 5];

task{1}{1}.randVars.calculated.resp = nan;
task{1}{1}.randVars.calculated.correct = nan;
task{1}{1}.randVars.calculated.diff = nan;
task{1}{1}.randVars.calculated.rt = nan;
task{1}{1}.randVars.calculated.centerint = nan;
 
% initialize the task
for phaseNum = 1:length(task{1})
  [task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end
 
% init the stimulus
myscreen = initStimulus('stimulus',myscreen);

% to initialize the stimulus for your experiment.
stimulus = initGaussian(stimulus,myscreen);
stimulus = initClick(stimulus,myscreen);
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

if stimulus.disp
dispPsychometric(task{1}{1});
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)
global stimulus
if task.thistrial.thisseg == 1
    stimulus.fixColor = stimulus.colors.white;%[1 1 1];
    task.thistrial.jitter = rand - 0.5; %random jittering between -0.5 and 0.5 deg
    % horizontal position of first, second stim
    if stimulus.task == 3
        if task.thistrial.centerWhich == 1
            task.thistrial.xposV = [task.thistrial.jitter, task.thistrial.posDiff + task.thistrial.jitter + task.thistrial.displacement];
            task.thistrial.xposA = [task.thistrial.jitter, task.thistrial.posDiff + task.thistrial.jitter - task.thistrial.displacement];
            task.thistrial.centerint = 1;
        else
            task.thistrial.xposV = [task.thistrial.jitter - task.thistrial.posDiff + task.thistrial.displacement, task.thistrial.jitter];
            task.thistrial.xposA = [task.thistrial.jitter - task.thistrial.posDiff - task.thistrial.displacement, task.thistrial.jitter];
            task.thistrial.centerint = 2;
        end
    else
      if ~task.thistrial.posDiff
        task.thistrial.xposV = [task.thistrial.jitter, task.thistrial.jitter];
        task.thistrial.xposA = [task.thistrial.jitter, task.thistrial.jitter];
      else
        if task.thistrial.centerWhich == 1
            task.thistrial.xposV = [task.thistrial.jitter, task.thistrial.posDiff + task.thistrial.jitter];
            task.thistrial.xposA = task.thistrial.xposV;
            task.thistrial.centerint = 1;
        else
            task.thistrial.xposV = [task.thistrial.jitter - task.thistrial.posDiff, task.thistrial.jitter];
            task.thistrial.xposA = task.thistrial.xposV;
            task.thistrial.centerint = 2;
        end
      end
    end

    task.thistrial.diff = task.thistrial.posDiff;
    
    if stimulus.task ~= 1 %auditory or bimodal condition
        for int = 1:2
            stimulus.sound(int) = createITD(stimulus,task.thistrial.xposA(int));
        end
    end

end
if task.thistrial.thisseg == 6
    if exist('task.thistrial.reactionTime', 'var')
        task.thistrial.rt = task.thistrial.reactionTime;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)
global stimulus
mglClearScreen(stimulus.colors.black);

if stimulus.task ~= 2 %visual or bimodal condition
    if task.thistrial.thisseg == stimulus.interval(1)
        mglBltTexture(stimulus.tex, [task.thistrial.xposV(1), 0]);
    elseif task.thistrial.thisseg == stimulus.interval(2)
        mglBltTexture(stimulus.tex, [task.thistrial.xposV(2), 0]);
    end
end
if stimulus.task ~= 1 %auditory or bimodal condition
    if task.thistrial.thisseg == stimulus.interval(1)
        mglPlaySound(stimulus.sound(1));
    elseif task.thistrial.thisseg == stimulus.interval(2)
        mglPlaySound(stimulus.sound(2));
    end
end

%draw fixation cross
if task.thistrial.thisseg == 5 || task.thistrial.thisseg == 6
    mglFixationCross(stimulus.fixWidth,1.5,stimulus.fixColor + stimulus.colors.reservedColor(2));
else
    mglFixationCross(stimulus.fixWidth,1.5,stimulus.fixColor);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)
global stimulus
% here, we just check whether this is the first time we got a response
if ~task.thistrial.gotResponse
    % which one seemed more to the LEFT
    % centerWhich (1/2) first or second one more eccentric
    if (task.thistrial.posDiff > 0 && task.thistrial.whichButton == 1) || ...
            (task.thistrial.posDiff < 0 && task.thistrial.whichButton == 2)
        % correct
        task.thistrial.correct = 1;
        % feeback
        stimulus.fixColor = stimulus.colors.green;%[0 1 0];
        disp(sprintf('(alaisburr) Trial %i: %0.4f resp %i correct', ...
            task.trialnum, task.thistrial.posDiff, task.thistrial.whichButton))
    else
        % incorrect
        task.thistrial.correct = 0;
        stimulus.fixColor = stimulus.colors.red;%[1 0 0];
        disp(sprintf('(alaisburr) Trial %i: %0.4f resp %i incorrect', ...
            task.trialnum, task.thistrial.posDiff, task.thistrial.whichButton))
    end
        
    task.thistrial.resp = task.thistrial.whichButton;
    task.thistrial.rt = task.thistrial.reactionTime;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initGaussian(stimulus,myscreen)
global stimulus;
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
gaussian = mglMakeGaussian(stimulus.width, stimulus.width, stimulus.width/8,stimulus.width/8);
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
stimulus.colors.white = stimulus.colors.reservedColor(2);
stimulus.colors.red = stimulus.colors.reservedColor(4);
stimulus.colors.green = stimulus.colors.reservedColor(6);
% % compute the guassian
% gauss = mglMakeGaussian(stimulus.width,stimulus.width, stimulus.width/8,stimulus.width/8);

% gaussian = zeros(size(gauss,1), size(gauss,2), 4);
% for i = 1:3
%     gaussian(:,:,i) = 255*ones(size(gauss,1), size(gauss,2));
% end
%     gaussian(:,:,4) = 255*gauss*stimulus.contrast;
    
% %create texture
% stimulus.tex = mglCreateTexture(gaussian);

 %stim centers
% [stimulus.x, stimulus.y] = pol2cart(0*pi/180,stimulus.eccentricity);
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
duration = stimulus.clickDur;
fs = stimulus.samplesPerSecond-1;
t = 0:1/fs:duration;
% frequency of signal in hz
% hz = 440;
% amplitude = 0.5;
% stimulus.wav = amplitude * sin(2*pi*hz*t);

wav = 0.5 * randn(1,length(t));
stimulus.wav = wav;

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
function dispPsychometric(task)
posDiff = task.parameter.posDiff;
n = zeros(1,length(posDiff)); k = zeros(1,length(posDiff));
% percent interval 2 (probe "Left")
for i = 1:length(posDiff)
    resp{i} = task.randVars.resp(task.randVars.diff == posDiff(i));
    n(i) = sum(resp{i} == 1 | resp{i}==2);
    k(i) = sum(resp{i} == 2);
end
percent = k./n;

figure;
h = plot(posDiff, percent, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
ylabel('Proportion of trials probe seen "left"');
xlabel('Displacement of probe (degs)');
axis([-20 20 0 1]); box off;

