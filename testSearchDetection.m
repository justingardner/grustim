% NOTES:
% (1) Uses filter instead of Contrast
% (2) Is functionated
% (3) Saves the locations into a task variable, making it easy to access in the SearchAnalysis script
% PROBLEM: how to randomize locations. Right now, if you put the locations
% as a task parameter, mgl does not pick columns at a time. It picks a
% random x and a random y
% NOTE: last updated on Dubonnet July 6


function myscreen = testSearch(varargin)
% check arguments
getArgs(varargin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initilaize the screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen.saveData = 1; myscreen = initScreen(myscreen); mglClearScreen; task{1}.waitForBacktick = 1; doEye =1;
contrast = 1;
global Colors; global stimulus;
mglClearScreen;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set task and stimulus parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
task{1}.seglen = [2.5 inf];
task{1}.getResponse = [0 1]; 
task{1}.numTrials = 18; % 3 locations x 3 target contrasts x 20 trials for each location, target contrast
task{1}.random=1; % each trial pulls random values from the parameters below 
% filters
filters = [0.2, 0.5, 0.9];
global randFilters; 
randFilters = filters(randperm(length(filters)));
task{1}.filters = randFilters;
% locations
task{1}.parameter.locations = [0 0 2 ;
                               5 8 2 ;];
% intialize response arrays %
task{1}.response.correct = [];
task{1}.response.filter = [];
% intialize arrays to hold (x,y) position 
task{1}.locations.x = [];
task{1}.locations.y = [];

% initialize the task
for phaseNum = 1:length(task)
    [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if doEye == 1; myscreen = eyeCalibDisp(myscreen); end;


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
% Start Segment Callback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)
global randFilters;
% (1) Making The Background Noise
stimBackground = makeStimBackground(myscreen);
    
if task.thistrial.thisseg == 2
        
    % (2) Making The Stencil
    mglStencilCreateBegin(1);
    mglFillOval(0,0,[17 17]);
    mglStencilCreateEnd;
    mglStencilSelect(1);
        
    % (3) Making the Target
    % (3.1) Dtermine the x,y position of the target 
    x = -task.thistrial.locations(1)
    y = -task.thistrial.locations(2)
    pixX = 38.8567214157064*x;
    pixY = 31.9291779098311*y;
    % Save the x,y position of the trial in a task variable
    task.locations.x = [task.locations.x -x];
    task.locations.y = [task.locations.y -y];
    % (3.2) Make the gaussian and the grating
    gaussian = mglMakeGaussian(60,60,0.1,0.1); [h w] = size(gaussian); gaussian = gaussian((h/2-400+pixY):(h/2+400+pixY),(w/2-400+pixX):(w/2+400+pixX)); 
    grating = mglMakeGrating(60,60,2,45,0); [h w] = size(grating); grating = grating((h/2-400+pixY):(h/2+400+pixY),(w/2-400+pixX):(w/2+400+pixX));
    % Setting the target contrast (i.e. the contrast of the grating)
    targetContrast = 1;
    grating = grating * targetContrast;
    
    % (4) Making the Final image (target embedded in background noise)
    % (4.1) Multiplying grating with background noise so that it blends into The final Image 
    grating = grating.*stimBackground;
    
    % (4.2) Assembling the grating windowed by the gaussian (i.e. a gabor) and the background noise windowed by the opposite of the gaussian
    % Pick a contrast based on the block (i.e trial number)
    if task.trialnum > 0 & task.trialnum <= 6
        filter = randFilters(1)
        stimImage = (grating.*gaussian)*filter + stimBackground;
        % (4.3) Actually creating the image through mgl 
        tex = mglCreateTexture(stimImage*255);
        mglBltTexture(tex,[0 0]);
    end
    if task.trialnum > 6 & task.trialnum <= 12
        filter = randFilters(2)
        stimImage = (grating.*gaussian)*filter + stimBackground;
        % (4.3) Actually creating the image through mgl 
        tex = mglCreateTexture(stimImage*255);
        mglBltTexture(tex,[0 0]);
    end
    if task.trialnum > 12 & task.trialnum <= 18
        filter = randFilters(3)
        stimImage = (grating.*gaussian)*filter + stimBackground;
        % (4.3) Actually creating the image through mgl 
        tex = mglCreateTexture(stimImage*255);
        mglBltTexture(tex,[0 0]);
    end
end 

if task.thistrial.thisseg == 1 
    mglFixationCross;
end

myscreen.flushMode = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame (STANDARD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% responseCallback  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = getResponseCallback(task,myscreen)

if task.thistrial.thisseg == 2 & task.thistrial.whichButton == 2
    task = jumpSegment(task);
end

mglClearScreen();


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Making background noise 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimBackground = makeStimBackground(myscreen)
% (1) Generating 1/f noise
noiseImage = makestim(myscreen);
% (2) Setting the RMS contrast
sumOfSquares = sum(sum(noiseImage.^2));
n = numel(noiseImage);     
backgroundRmsContrast = 0.39;  
rmsAdjust = sqrt(sumOfSquares/(n*(backgroundRmsContrast)^2)); 
stimBackground = noiseImage / rmsAdjust;
