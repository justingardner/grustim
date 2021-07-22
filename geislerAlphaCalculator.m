function alpha = calculateAlpha()
myscreen = initScreen;

[gaussian grating] = makeGrating(myscreen);
gabor = grating.*gaussian;
normalizer = dot(gabor(:), gabor(:))

templateResponses = [];

for i=1:10000
    
    stimBackground = makeStimBackground(myscreen);
    
    templateResponse = dot(stimBackground(:), gabor(:));
   
    templateResponses = [templateResponses templateResponse];
    
end

normalizedTemplateResponses = templateResponses / normalizer;
alpha = std(normalizedTemplateResponses)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Making background noise 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimBackground = makeStimBackground(myscreen)
% (1) Generating 1/f noise
noiseImage = makestim(myscreen);
% (2) Subtracting the mean to center around 
noiseImageMean = mean(noiseImage(:));
noiseImage = noiseImage - noiseImageMean;
% (3) Adjusting the range to [-1,1] (Written by Justin)
noiseImage = 2*noiseImage / (max(noiseImage(:))-min(noiseImage(:)));
% (4) Setting the RMS contrast
sumOfSquares = sum(noiseImage(:).^2);
n = numel(noiseImage);     
backgroundRmsContrast = 0.25;  
rmsAdjust = sqrt(sumOfSquares/(n*(backgroundRmsContrast)^2)); 
stimBackground = noiseImage / rmsAdjust;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Making grating and gaussian for the target
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gaussian grating] = makeGrating(task,myscreen)

% Determining the location of the target 
x = 0;
y = 0;
pixX = 38.8567214157064*x;
pixY = 31.9291779098311*y;
% NOTE: the parameters are set s.t. the FWHM of the Gaussian is equal to 1/cpd of the grating
gaussian = mglMakeGaussian(60,60,0.16,0.16); [h w] = size(gaussian); gaussian = gaussian((h/2-400+pixY):(h/2+400+pixY),(w/2-400+pixX):(w/2+400+pixX)); 
grating = mglMakeGrating(60,60,2.65413,45,0); [h w] = size(grating); grating = grating((h/2-400+pixY):(h/2+400+pixY),(w/2-400+pixX):(w/2+400+pixX));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating 1/f noise 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function noiseImage = makestim(myscreen)
global Colors
width = 2;
SNR = 1;
maxSNR = .3;
xPos = 10; yPos = 5;

% get the gaussian full size of the screen
[Background.gaussian Background.x Background.y] = mglMakeGaussian(myscreen.imageWidth, myscreen.imageHeight, width, width);

% now make sure that dimensions are odd numbers of pixels
oddWidth = 2*floor(myscreen.screenWidth/2)+1;
oddHeight = 2*floor(myscreen.screenHeight/2)+1;

% resize everything to odd
% Background.gaussian = Background.gaussian(1:oddHeight,1:oddWidth);
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
% BackTexture = mglCreateTexture(round(Colors.gaussRange*squeeze(Background.im(:,:)) + Colors.minGaussianIndex));

% Images
fullImage = squeeze(Background.im(:,:)) + Background.sigMax * exp(-((((Background.x-xPos).^2) + (Background.y+yPos).^2))/(2*(width^2)));
noiseImage = Background.im(100:900,100:900);
stimulusImage = Background.gaussian;

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

% setGammaTableForMaxContrast %
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




