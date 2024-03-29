% WRITTEN BY:
% (1) Yehia Elkersh 
% (2) Josh Wilson (generating 1/f noise)

% DESCRIPTION: 
% This script is set up to run the detection task in the Najemnik & Geisler 2005 Nature paper. 

% THE TASK: This is a basic 2AFC task, where ach "block" refers to a set of trials where the target appears in one particular location, 
% which should be indicated to the subject at the beginning of each block. Within each block (i.e. at each location) the target contrasts is varied randomly. 

% FAULTS:
% As of July 15, 2021, this scripts has a few faults
% (1) It only runs on two locations, whereas N&G use 25, but that should be easy to fix.  
% (2) More importantly, it uses one set of target contrasts for all the locations, whereas N&G use use a different set of contrasts for each location. 
% (3) The way the script "blocks" the trials is by using the trial number (e.g. for the first 200 trials, it puts the target in location x1, 
% for the next 200 it puts it in location x2, etc.). This is not ideal because it does not guarantee that the target contrasts gets randomized properly over in each block. 
% For instance, targetContrat might be 0.5 more often in location x1 than in location x2. 
% (4) The script does not inform the subject where the location is going to be for the proceeding block of trials.

% SOLUTION: 
% To fix the stated faults, one should split up the task into phases, such that phase{1} is text (or illustrations) telling the subject that the target
% will be at location x1, phase{2} is running the block of trials at location x1, phase(3) is informing the subject that the target will be at location x2, etc. 
% This will also make it much relatively easy to have different contrasts for each location, and will ensure that the contrasts are properly randomized.

% NOTES: 
% (1) This scripts saves the locations into a task variable (task{1}.locations) in order to make it easy to access in geislerDetectionAnalysisMultipleLocs

    
function myscreen = testSearch(varargin)
% check arguments
getArgs(varargin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initilaize the screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen.saveData = 1; myscreen = initScreen(myscreen); mglClearScreen(0.5); task{1}.waitForBacktick = 1; 
contrast = 1;
global Colors; global stimulus;
mglClearScreen(0.5);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set task and stimulus parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
task{1}.seglen = [2.5 0.25 0.5 0.25 inf 0.5];
task{1}.getResponse = [0 0 0 0 1 0]; 
task{1}.numTrials = 544; 
task{1}.random=1; % each trial pulls random values from the parameters below 
task{1}.parameter.contrast = [0 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.225 0.25 0.275 0.3 0.35 0.4 0.5 0.7];
% Determines which segent to embed the target in
% For instance, if whichSegmemt = 1, then embed the target in the first segment
task{1}.parameter.whichSegment = [1 2];
% intialize response arrays 
task{1}.response.correct = [];
task{1}.response.contrast = [];

% initialize locations array and save it in a task variable, making it easy to access in the analysis script
locations = [0 0;
             0 4.5;];
global randomLocations;
randomLocations = locations(randperm(size(locations, 1)), :); 
task{1}.locations = randomLocations;

% Used to change the color of the cross on the 5th segment (i.e  to give feedback)
global correct;


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
% Start Segment Callback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)
global correct

% (1) Making The Background Noise
stimBackground = makeStimBackground(myscreen);
    
if task.thistrial.thisseg == 2  
    % if the temporally first segment (seg 2) has the target, then embed the target in the background noise
    if task.thistrial.whichSegment == 1
        
        % (2) Making The Stencil
        mglStencilCreateBegin(1);
        mglFillOval(0,0,[15 15]);
        mglStencilCreateEnd;
        mglStencilSelect(1);
        
        % (3) Making the Target
        [gaussian grating] = makeGrating(task,myscreen);
    
        % (4) Making the Final image (target embedded in background noise)
        % (4.1) Making the Gabor and adjusting its contrast
        gabor = grating.*gaussian; 
        gabor = gabor * task.thistrial.contrast;
        % (4.2) Adding the Gabor to the background and clipping values outside [-1, 1] range
        stimImage = gabor + stimBackground;
        stimImage(stimImage > 1) = 1;
        stimImage(stimImage < -1) = -1;
        % (4.3) Actually creating the image through mgl 
        tex = mglCreateTexture(((stimImage+1)/2)*255);
        mglBltTexture(tex,[0 0]);
        
        mglFixationCross(0.4, 1, [0 0 0]);
    end
    % if it does not have the target, then just present background noise
    if task.thistrial.whichSegment == 2
        
        % (1) Making The Stencil
        mglStencilCreateBegin(1);
        mglFillOval(0,0,[15 15]);
        mglStencilCreateEnd;
        mglStencilSelect(1);
        
        % (2) Clipping values outside [-1, 1] range 
        stimBackground(stimBackground > 1) = 1;
        stimBackground(stimBackground < -1) = -1;
        
        % (3) Actually creating the image through mgl
        tex = mglCreateTexture(((stimBackground+1)/2)*255);
        mglBltTexture(tex,[0 0]);
        
        mglFixationCross(0.4, 1, [0 0 0]);
    end
end
    
if task.thistrial.thisseg == 4
    % if the temporally second segment (seg 4) has the target, then embed the target in the background noise
    if task.thistrial.whichSegment == 2
        
        % (2) Making The Stencil
        mglStencilCreateBegin(1);
        mglFillOval(0,0,[15 15]);
        mglStencilCreateEnd;
        mglStencilSelect(1);
        
        % (3) Making the Target
        [gaussian grating] = makeGrating(task,myscreen);
    
        % (4) Making the Final image (target embedded in background noise)
        % (4.1) Making the Gabor and adjusting its contrast
        gabor = grating.*gaussian; 
        gabor = gabor * task.thistrial.contrast;
        % (4.2) Adding the Gabor to the background and clipping values outside [-1, 1] range
        stimImage = gabor + stimBackground;
        stimImage(stimImage > 1) = 1;
        stimImage(stimImage < -1) = -1;
        % (4.3) Actually creating the image through mgl 
        tex = mglCreateTexture(((stimImage+1)/2)*255);
        mglBltTexture(tex,[0 0]);
        
        mglFixationCross(0.4, 1, [0 0 0]);
    end 
    % if it does not have the target, then just present background noise
    if task.thistrial.whichSegment == 1
        
       % (1) Making The Stencil
        mglStencilCreateBegin(1);
        mglFillOval(0,0,[15 15]);
        mglStencilCreateEnd;
        mglStencilSelect(1);
        
        % (2) Clipping values outside [-1, 1] range 
        stimBackground(stimBackground > 1) = 1;
        stimBackground(stimBackground < -1) = -1;
        
        % (3) Actually creating the image through mgl
        tex = mglCreateTexture(((stimBackground+1)/2)*255);
        mglBltTexture(tex,[0 0]);
        
        mglFixationCross(0.4, 1, [0 0 0]);
    end
end

if task.thistrial.thisseg == 1 | task.thistrial.thisseg == 3 | task.thistrial.thisseg == 5 
    mglClearScreen(0.5);
    mglFixationCross(0.4, 1, [0 0 0]);
end

if task.thistrial.thisseg == 6
    if correct == 1
        % if the subject was correct, display a green cross 
        mglFixationCross(0.4, 2, [0 255 0]);
    end
    if correct == 0
        % if the subject was incorrect, display a red cross 
        mglFixationCross(0.4, 2, [255 0 0]);
    end
    
    
end
myscreen.flushMode = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% responseCallback  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = getResponseCallback(task,myscreen)
global correct;

if task.thistrial.whichButton == 1 & task.thistrial.thisseg == 5
    % if the button is 1 and the 1st segment had the target, mark as correct
    if task.thistrial.whichSegment == 1 
        task.response.correct = [task.response.correct 1];
        correct = 1;
    end
    % if the button is 1 and the 2nd segment had the target, mark as incorrect
    if task.thistrial.whichSegment == 2
        task.response.correct = [task.response.correct 0];
        correct = 0;
    end
    % Adding the contrast to the reponse struct
    task.response.contrast = [task.response.contrast task.thistrial.contrast];
    task = jumpSegment(task);
end
mglClearScreen(0.5);

if task.thistrial.whichButton == 2 & task.thistrial.thisseg == 5
    % if the button is 2 and the 1st segment had the target, mark as incorrect
    if task.thistrial.whichSegment == 1 
        task.response.correct = [task.response.correct 0];
        correct = 0;
    end
    % if the button is 2 and the 2nd segment had the target, mark as correct
    if task.thistrial.whichSegment == 2
        task.response.correct = [task.response.correct 1];
        correct = 1;
    end
    % Adding the contrast to the reponse struct
    task.response.contrast = [task.response.contrast task.thistrial.contrast];
    
    task = jumpSegment(task);
end

if task.thistrial.whichButton && task.thistrial.thisseg == 1
    task = jumpSegment(task);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Making background noise 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimBackground = makeStimBackground(myscreen)
% (1) Generating 1/f noise
noiseImage = makestim(myscreen);
% (2) Subtracting the mean to center around 0
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
% Making the grating and gaussian for the target
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gaussian grating] = makeGrating(task,myscreen)
global randomLocations;

% Determining the location of the target based on the session (trial number). For now, we are only doing two location with 544 trials for each location
if task.trialnum > 0 & task.trialnum <= 544
    locationVector = randomLocations(1,:);
    x = -locationVector(1);
    y = -locationVector(2);
end
if task.trialnum > 544 & task.trialnum <= 1088
    locationVector = randomLocations(2,:);
    x = -locationVector(1);
    y = -locationVector(2);
end
pixX = 38.8567214157064*x;
pixY = 31.9291779098311*y;
% NOTE: the parameters of mglMakeGaussian are set s.t. the FWHM of the Gaussian is equal to 1/cpd of the grating
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




