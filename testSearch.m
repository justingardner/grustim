function myscreen = testSearch(varargin)
% check arguments
getArgs(varargin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initilaize the screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = initScreen(); mglClearScreen; task{1}.waitForBacktick = 1; 
contrast = .2
global Colors; global stimulus;
Colors.reservedColors = [1 1 1; 0.3 0.3 0.3; 0 1 0;1 0 0; 0 1 1];
Colors.nReservedColors = size(Colors.reservedColors,1);
maxIndex = 255;
Colors.nGaussianColors = maxIndex+1-Colors.nReservedColors;
Colors.minGaussianIndex = maxIndex+1 - Colors.nGaussianColors;
Colors.maxGaussianIndex = maxIndex;
stimulus.linearizedGammaTable = mglGetGammaTable;
Colors.nDisplayContrasts = floor(Colors.nGaussianColors-1);
setGammaTableForMaxContrast(contrast);
contrastIndex = getContrastIndex(contrast,1);
Colors.gaussRange = contrastIndex-1;
mglClearScreen;

% Code specific to generating noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set task parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
task{1}.segmin = [2.5 1 0.5 1 inf]; task{1}.segmax = [2.5 1 0.5 1 inf]; task{1}.getResponse = [1 0 0 0 1]; 
task{1}.numTrials = 500; % number of trials you want it to run 
task{1}.random=1; % each trial pulls random values from the parameters below 
task{1}.parameter.x = [-3 3]
task{1}.parameter.y = [-3 3]
% intialize response arrays %1
task{1}.response.correct = []

% initialize the task
for phaseNum = 1:length(task)
 [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
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
% Set the variables for drawing a picture in trials that have a picture

% Draw the picture in trials that have a picture
if task.thistrial.thisseg == 2 | task.thistrial.thisseg == 4
 mglStencilCreateBegin(1);
 mglFillOval(0,0,[15 15]);
 mglStencilCreateEnd;
 mglStencilSelect(1);
 [noiseImage] = makestim(myscreen);
 
 
 %gabor
 x = -task.thistrial.x
 y = -task.thistrial.y
 pixX = 38.8567214157064*x;
 pixY = 31.9291779098311*y;
 gaussian = mglMakeGaussian(60,60,.1,.1); [h w] = size(gaussian); gaussian = gaussian((h/2-400+pixY):(h/2+400+pixY),(w/2-400+pixX):(w/2+400+pixX));
 grating = mglMakeGrating(60,60,2,45,0); [h w] = size(grating); grating = grating((h/2-400+pixY):(h/2+400+pixY),(w/2-400+pixX):(w/2+400+pixX));
 grating = grating.*noiseImage;
 
 stimImage = grating.*gaussian   +    noiseImage.*(1-gaussian);
 
 
 %%%% edits this %%%%%%
%create a bitmap of -1:1 then set it from 0:255
 
 
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
Background.gaussian = imread('pic01.png'); Background.gaussian = imresize(Background.gaussian,[oddHeight oddWidth]);
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

% Images
fullImage = squeeze(Background.im(:,:)) + Background.sigMax * exp(-((((Background.x-xPos).^2) + (Background.y+yPos).^2))/(2*(width^2)));
noiseImage = Background.im(100:900,100:900);
stimulusImage = Background.gaussian;

% scale and set on texture // unnecessary right now
stimulus.stimTexture = mglCreateTexture(round(Colors.gaussRange*im + Colors.minGaussianIndex));






%%%%%%%%%%%%%%%%%%%%%%
%% getContrastIndex %%
%%%%%%%%%%%%%%%%%%%%%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setGammaTableForMaxContrast %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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