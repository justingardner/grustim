% NOTES:
% (1) Uses Contrast (Instead of filter) and adjusts pixel range to [-1, 1]
% (2) Is not functionated
% (3) This does not run with a visually correct stimulus on any monitor but
% can be used for reference for the following:
    % 1. Code to use contrast instead of filter
    % 2. Code to adjust pixel range to [-1,1]
    % 3. Original global variables that Josh used for creating 1/f noise (here
    % they are commented out)

function myscreen = zz(varargin)
% check arguments
getArgs(varargin);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initilaize the screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = initScreen(); mglClearScreen; task{1}.waitForBacktick = 1; 
contrast = 0.25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code specific to generating noise
%global Colors; global stimulus; global randomLocations;`
%Colors.reservedColors = [1 1 1; 0.3 0.3 0.3; 0 1 0;1 0 0; 0 1 1];
%Colors.nReservedColors = size(Colors.reservedColors,1);
%maxIndex = 255;
%Colors.nGaussianColors = maxIndex+1-Colors.nReservedColors;
%Colors.minGaussianIndex = maxIndex+1 - Colors.nGaussianColors;
%Colors.maxGaussianIndex = maxIndex;
%stimulus.linearizedGammaTable = mglGetGammaTable;
%Colors.nDisplayContrasts = floor(Colors.nGaussianColors-1);
%setGammaTableForMaxContrast(contrast);
%contrastIndex = getContrastIndex(contrast,1);
%Colors.gaussRange = contrastIndex-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mglClearScreen;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set task parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
task{1}.segmin = [2.5 5 0.5 3 inf]; task{1}.segmax = [2.5 5 0.5 3 inf]; task{1}.getResponse = [1 0 0 0 1]; 
task{1}.numTrials = 6;
task{1}.random=1; % each trial pulls random values from the parameters below 
task{1}.parameter.targetContrast = [0.5];
task{1}.parameter.whichSegment = [1 2];
% intialize response arrays %
task{1}.response.correct = []
% initialize locations arrays 
locations = [0 0; 0 2; 0 4; 0 6; 0 8; 
                  2 2; 4 4; 6 6; 8 8;
                  2 0; 4 0; 6 0; 8 0;
                  2 -2; 4 -4; 6 -6; 8 -8; 
                  0 -2; 0 -4; 0 -6; 0 -8;
                  -2 -2; -4 -4; -6 -6; -8 -8;
                  -2 0; -4 0; -6 0; -8 0;
                  -2 2; -4 4; -6 6; -8 8;];
randomLocations = locations(randperm(size(locations, 1)), :);

% initialize the task
for phaseNum = 1:length(task)
    [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop (STANDARD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task)) && ~myscreen.userHitEsc
    % update the task
    [task myscreen phaseNum] = updateTask(task,myscreen,phaseNum);
    % flip screen
    myscreen = tickScreen(myscreen,task);
end
% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)
global randomLocations;
% Draw the stimulus in segments that contain the stimulus
if task.thistrial.thisseg == 2 | task.thistrial.thisseg == 4
    
    % (1) Making The Target 
    % (1.1) The Stencil
    mglStencilCreateBegin(1);
    mglFillOval(0,0,[17 17]);
    mglStencilCreateEnd;
    mglStencilSelect(1);
    % (1.2) The Gausian
    % Determining the location of the target based on the session (trial number)
    if task.trialnum > 0 & task.trialnum <= 2
        locationVector = randomLocations(1,:);
        x = locationVector(1)
        y = locationVector(2)
    end
    if task.trialnum > 2 & task.trialnum <= 4
        locationVector = randomLocations(2,:);
        x = locationVector(1)
        y = locationVector(2)
    end
    if task.trialnum > 4 & task.trialnum <= 6
        locationVector = randomLocations(3,:);
        x = locationVector(1);
        y = locationVector(2);
    end
    pixX = 38.8567214157064*x;
    pixY = 31.9291779098311*y;
    gaussian = mglMakeGaussian(60,60,0.1,0.1); [h w] = size(gaussian); gaussian = gaussian((h/2-400+pixY):(h/2+400+pixY),(w/2-400+pixX):(w/2+400+pixX)); 
    % (1.3) The Grating
    grating = mglMakeGrating(60,60,4,45,0); [h w] = size(grating); grating = grating((h/2-400+pixY):(h/2+400+pixY),(w/2-400+pixX):(w/2+400+pixX));
    % Setting the contrast of the graing (i.e. the target contrast)
    targetContrast = task.thistrial.targetContrast;
    grating = grating * targetContrast;
    
    % (2) Making The Background Noise
    % (2.1) Generating 1/f noise
    [noiseImage] = makestim(myscreen);
    % (2.2) Subtract the mean to center around 0 and multiply by 2 to get [-1, 1] range
    noiseImageMean = mean(noiseImage(:));
    noiseImage = noiseImage - noiseImageMean;
    noiseImage = 2 * noiseImage;
    % (2.3) Setting the RMS contrast
    sumOfSquares = sum(sum(noiseImage.^2));
    n = numel(noiseImage);     
    backgroundRmsContrast = 0.25;  
    rmsAdjust = sqrt(sumOfSquares/(n*(backgroundRmsContrast)^2)); 
    noiseImage = noiseImage / rmsAdjust;
    
    % (3) Making the Final image (Target embedded in background noise)
    % (3.1) Multiplying grating with background noise so that it blends into The final Image 
    grating = grating.*noiseImage;   
    % (3.2) Assembling the grating windowed by the gaussian (i.e. a gabor) and the background noise windowed by the reverse of the gaussian
    stimImage = grating.*gaussian + noiseImage.*(1-gaussian);
    % (3.3) Actually creating the image through mgl 
    tex = mglCreateTexture(stimImage*255);
    mglBltTexture(tex,[0 0]);
    
elseif task.thistrial.thisseg == 1 | task.thistrial.thisseg == 3
    mglFixationCross;
elseif task.thistrial.thisseg == 5
    mglClearScreen();
end
myscreen.flushMode = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%resp%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame (STANDARD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% responseCallback  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = getResponseCallback(task,myscreen)
if task.thistrial.whichButton == 1 & task.thistrial.thisseg == 5
 task = jumpSegment(task);
end
mglClearScreen();
if task.thistrial.whichButton == 2 & task.thistrial.thisseg == 5
 task = jumpSegment(task);
end

if task.thistrial.whichButton && task.thistrial.thisseg ==1
    task = jumpSegment(task)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating 1/f noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [noiseImage] = makestim(myscreen);
global Colors
width = 2;
SNR = 1;
maxSNR = 0;
xPos = 10; yPos = 5;

% get the gaussian full size of the screen
[Background.gaussian Background.x Background.y] = mglMakeGaussian(myscreen.imageWidth, myscreen.imageHeight, width, width);

% now make sure that dimensions are odd numbers of pixels
oddWidth = 2*floor(myscreen.screenWidth/2)+1;
oddHeight = 2*floor(myscreen.screenHeight/2)+1;

% resize everything to odd
%Background.gaussian = Background.gaussian(1:oddHeight,1:oddWidth);
Background.gaussian = imread('pic03.png'); Background.gaussian = imresize(Background.gaussian,[oddHeight oddWidth]);
Background.x = Background.x(1:oddHeight,1:oddWidth);
Background.y = Background.y(1:oddHeight,1:oddWidth);
  
% get the fourier transform
Background.gaussianTransform = getHalfFourier(Background.gaussian);
  
% pull out magnitude and dc for averaging
mag = Background.gaussianTransform.mag;
dc = Background.gaussianTransform.dc;

% make average transform
Background.averageGaussianTransform = Background.gaussianTransform;
Background.averageGaussianTransform.dc = dc;
Background.averageGaussianTransform.mag = mag;

% max noise and signal
Background.noiseMax = 1 / (maxSNR + 1);
Background.sigMax = SNR * Background.noiseMax;

% randomize phase and reconstruct
Background.averageGaussianTransform.phase = (rand(1,Background.averageGaussianTransform.n)*2*pi - pi);
im = reconstructFromHalfFourier(Background.averageGaussianTransform);

% scale from 0 to noise max
maxIm = max(im(:));
minIm = min(im(:));
Background.im(:,:) = Background.noiseMax * (im - minIm) / (maxIm-minIm);

% make into texture // currently unused
%BackTexture = mglCreateTexture(round(Colors.gaussRange*squeeze(Background.im(:,:)) + Colors.minGaussianIndex));

% Images`
fullImage = squeeze(Background.im(:,:)) + Background.sigMax * exp(-((((Background.x-xPos).^2) + (Background.y+yPos).^2))/(2*(width^2)));
noiseImage = Background.im(100:900,100:900);
stimulusImage = Background.gaussian;

% scale and set on texture // unnecessary right now
stimulus.stimTexture = mglCreateTexture(round(Colors.gaussRange*im + Colors.minGaussianIndex));

% getContrastIndex 
function contrastIndex = getContrastIndex(desiredContrast,verbose)

if nargin < 2,verbose = 0;end

global stimulus; global Colors;
if desiredContrast < 0, desiredContrast = 0;end

% now find closest matching contrast we can display with this gamma table
contrastIndex = min(round(Colors.nDisplayContrasts*desiredContrast/stimulus.currentMaxContrast),Colors.nDisplayContrasts);

% display the desired and actual contrast values if verbose is set
if verbose
  actualContrast = stimulus.currentMaxContrast*(contrastIndex/Colors.nDisplayContrasts);
  disp(sprintf('(getContrastIndex) Desired contrast: %0.4f Actual contrast: %0.4f Difference: %0.4f',desiredContrast,actualContrast,desiredContrast-actualContrast));
end

% out of range check
if round(Colors.nDisplayContrasts*desiredContrast/stimulus.currentMaxContrast)>Colors.nDisplayContrasts
 disp(sprintf('(getContrastIndex) Desired contrast (%0.9f) out of range max contrast : %0.9f',desiredContrast,stimulus.currentMaxContrast));
 keyboard
end

% 1 based indexes (0th index is gray, nDisplayContrasts+1 is full contrast)
contrastIndex = contrastIndex+1;

% setGammaTableForMaxContrast 
function setGammaTableForMaxContrast(maxContrast)

global Colors; global stimulus;
% if you just want to show gray, that's ok, but to make the
% code work properly we act as if you want to display a range of contrasts
if maxContrast <= 0,maxContrast = 0.01;end

% set the reserved colors
gammaTable(1:size(Colors.reservedColors,1),1:size(Colors.reservedColors,2))=Colors.reservedColors;

% set the gamma table
if maxContrast > 0
  % create the rest of the gamma table
  cmin = 0;
  cmax = maxContrast;
  luminanceVals = cmin:((cmax-cmin)/(Colors.nGaussianColors-1)):cmax;

  % replace NaN in gamma tables with zero
  stimulus.linearizedGammaTable.redTable(isnan(stimulus.linearizedGammaTable.redTable)) = 0;
  stimulus.linearizedGammaTable.greenTable(isnan(stimulus.linearizedGammaTable.greenTable)) = 0;
  stimulus.linearizedGammaTable.blueTable(isnan(stimulus.linearizedGammaTable.blueTable)) = 0;

  % now get the linearized range
  redLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,luminanceVals,'linear');
  greenLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,luminanceVals,'linear');
  blueLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,luminanceVals,'linear');
  
  % add these values to the table
  gammaTable((Colors.minGaussianIndex:Colors.maxGaussianIndex)+1,:)=[redLinearized;greenLinearized;blueLinearized]';
else
  % if we are asked for 0 contrast then simply set all the values to BLACK
  gammaTable((Colors.minGaussianIndex:Colors.maxGaussianIndex)+1,1)=interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,0,'linear');
  gammaTable((Colors.minGaussianIndex:Colors.maxGaussianIndex)+1,2)=interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,0,'linear');
  gammaTable((Colors.minGaussianIndex:Colors.maxGaussianIndex)+1,3)=interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,0,'linear');
end

% set the gamma table
mglSetGammaTable(gammaTable);

% keep the gamma table
stimulus.gammaTable = gammaTable;

% remember what the current maximum contrast is that we can display
stimulus.currentMaxContrast = maxContrast;