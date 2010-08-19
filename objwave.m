% objwave.m
%
%        $Id:$ 
%      usage: objwave(<'faceDir=faceDir'>,<'placeDir=placeDir'>)
%         by: justin gardner
%       date: 06/27/10
%    purpose: object localizer scan
%
% options:
% objloc('imageDir=somedir') 
%   the directory which stores images
% objloc('categories',{'faces','houses') % categories
%   of images, should be directories of image under imageDir for everything
% objloc('keepAspectRatio=0')
%    defaults to zero - keeps the aspect ratio of the original
% objloc('widthPix=180','heightPix=180','widthDeg','heightDeg')
%    image sizes
% objloc('repeatFreq=0.1')
%    frequency with which to repeat images
% objloc('waitForBacktick=1')
%    whether to wait for backtick (for running in scanner)
function retval = objwave(varargin)

% get arguments
categories = [];
imageDir = [];
dispLoadFig = [];
categoryWeight = [];
keepAspectRatio = [];
repeatFreq = [];
waitForBacktick = [];
widthPix = [];
heightPix = [];
widthDeg = [];
heightDeg = [];
getArgs(varargin,{'categories',{'human_face','building'},'imageDir=~/proj/grustim/images/ObjLocImages','dispLoadFig=0','keepAspectRatio=0','repeatFreq=0.1','waitForBacktick=1','widthPix=180','heightPix=180','widthDeg=18','heightDeg=18'});

% initalize the screen
myscreen.background = 'gray';
myscreen = initScreen(myscreen);
% init the stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);

% set up task
task{1}.waitForBacktick = waitForBacktick;
task{1}.seglen = repmat([0.75 0.25],1,26);
if waitForBacktick,task{1}.seglen(end) = 0.1;end
task{1}.getResponse = ones(1,length(task{1}.seglen));
task{1}.getResponse(1:2) = 0;
task{1}.synchToVol = zeros(1,length(task{1}.seglen));
task{1}.synchToVol(end) = waitForBacktick;
task{1}.numBlocks = 100;
task{1}.random = 1;

% initialize the task
for phaseNum = 1:length(task)
  [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@stimStartSegmentCallback,@stimDrawStimulusCallback,@responseCallback);
end

% load images
stepsPerCycle = length(task{1}.seglen)/2;
range = .75;midPoint = 0.4;
oneCycle = cos(0:2*pi/(stepsPerCycle):2*pi);
oneCycle = (range/2)*oneCycle(1:end-1)+midPoint;
scrambleFactors = oneCycle;
stimulus.widthPix = widthPix;
stimulus.heightPix = heightPix;
stimulus.widthDeg = widthDeg;
stimulus.heightDeg = heightDeg;
stimulus.repeatFreq = repeatFreq;
stimulus = myInitStimulus(stimulus,myscreen,scrambleFactors,categories,imageDir,dispLoadFig,keepAspectRatio);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

mglFixationCross(1,2,[0 0 0]);
myscreen.flushMode = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task)) && ~myscreen.userHitEsc
  % update the task
  [task myscreen phaseNum] = updateTask(task,myscreen,phaseNum);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% delete texture
for i  = 1:length(stimulus.tex(:))
  mglDeleteTexture(stimulus.tex(i));
end
stimulus = rmfield(stimulus,'tex');

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = stimStartSegmentCallback(task, myscreen)

global stimulus;

if isodd(task.thistrial.thisseg)
  % clear screen
  mglClearScreen;
  if task.thistrial.thisseg == 1
    disp(sprintf('%i: %0.2f volnum: %i',task.trialnum,mglGetSecs(stimulus.trialStart),myscreen.volnum));
    stimulus.trialStart = mglGetSecs;
  end
  % display a random image
  % see if we want to do a repeat
  if (task.thistrial.thisseg > 1) && (rand < stimulus.repeatFreq);
    stimulus.isRepeat = 1;
  else
    iPhase = (task.thistrial.thisseg+1)/2;
    stimulus.thisTex = stimulus.tex(task.trialnum,iPhase);
    stimulus.isRepeat = 0;
  end
  
  % draw the textures
  mglBltTexture(stimulus.thisTex,[0 0 stimulus.widthDeg stimulus.heightDeg]);
  % fixation cross
  mglFixationCross(1,2,[0 0 0]);
  % tell mgl/task not to flush after displaying image
  myscreen.flushMode = 1;
else
  % just display fixation cross
  myscreen.flushMode = 1;
  mglClearScreen;
  mglFixationCross(1,2,[0 0 0]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = stimDrawStimulusCallback(task, myscreen)

global stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called every response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task, myscreen)

global stimulus;

fixColor = [1 0 1];
mglClearScreen;
% redisplay image
if isodd(task.thistrial.thisseg)
  mglBltTexture(stimulus.thisTex,[0 0 stimulus.widthDeg stimulus.heightDeg]);;
end
% check correct
if stimulus.isRepeat
  disp('(objloc) Correct');
  fixColor = [0 1 0];
else
  disp('(objloc) Incorrect')
  fixColor = [1 0 0];
end

% change fixation color
mglFixationCross(1,1,fixColor);
myscreen.flushMode = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen,scrambleFactors,categories,imageDir,dispFig,keepAspectRatio)

if ~isfield(stimulus,'objwave'),stimulus.imagesLoaded = 0;end
stimulus.objwave = 1;
stimulus.thisTex = [];

% make sure widht and height are odd
if iseven(stimulus.widthPix), stimulus.widthPix = stimulus.widthPix-1;end
if iseven(stimulus.heightPix), stimulus.heightPix = stimulus.heightPix-1;end

% check whether images are loaded
averageN = 0;
if ~isfield(stimulus,'imagesLoaded') || (~stimulus.imagesLoaded) || ~isequal(stimulus.categories,categories) || ~isequal(stimulus.imageDir,imageDir) || dispFig
  stimulus.nCategories = length(categories);

  % keep the averageMag and averageDc so that we can normalize images
  stimulus.averageMag = 0;
  stimulus.averageDC = 0;

  for i = 1:stimulus.nCategories
    if ~any(strcmp(categories{i},{'scramble','blank','gray'}))
      % load images
      stimulus.raw{i} = loadNormalizedImages(fullfile(imageDir,categories{i}),'width',stimulus.widthPix,'height',stimulus.heightPix,'dispFig',dispFig,'keepAspectRatio',keepAspectRatio);

      % make sure we opened ok
      if isempty(stimulus.raw{i})
	disp(sprintf('(objloc) Could not load images; %s',categories{i}));
	keyboard
      end
      % get average mag
      stimulus.averageMag = stimulus.averageMag + stimulus.raw{i}.averageMag;
      stimulus.averageDC = stimulus.averageDC +stimulus.raw{i}.averageDC;
      averageN = averageN + 1;
    else
      stimulus.raw{i}.n = 0;
    end
  end
  % compute average
  stimulus.averageMag = stimulus.averageMag/averageN;
  stimulus.averageDC = stimulus.averageDC/averageN;
  % and set that we have loaded
  stimulus.imageDir = imageDir;
  stimulus.imagesLoaded = 1;

else
  disp(sprintf('(objloc) Stimulus already initialized'));
end

% count how many images
nImages = [];
for i = 1:stimulus.nCategories
  if stimulus.raw{i}.n > 0
    % keep number of images so that we can make scrambles with equal number of images
    nImages(end+1) = stimulus.raw{i}.n;
  end
end

stimulus.scrambleFactors = scrambleFactors;
stimulus.nCycles = 12;

disppercent(-inf,'(objloc) Crateating mixture images');
for iCycle = 1:stimulus.nCycles
  for iScrambleFactor = 1:length(stimulus.scrambleFactors)
    % otherwise get a random image number
    randImageNum1 = ceil(rand*stimulus.raw{1}.n);
    randImageNum2 = ceil(rand*stimulus.raw{2}.n);
    image1 = stimulus.raw{1}.halfFourier{randImageNum1};
    image2 = stimulus.raw{2}.halfFourier{randImageNum2};
    % create the mixture
    stimulus.tex(iCycle,iScrambleFactor) = makeMixImage(stimulus.scrambleFactors(iScrambleFactor),image1,image2,stimulus.averageDC,stimulus.averageMag);
    disppercent(calcPercentDone(iCycle,stimulus.nCycles,iScrambleFactor,length(stimulus.scrambleFactors)));
  end
end
disppercent(inf);

stimulus.trialStart = mglGetSecs;
stimulus.categories = categories;
%%%%%%%%%%%%%%%%%%%%%%%%
%    getHalfFourier    %
%%%%%%%%%%%%%%%%%%%%%%%%
function d = getHalfFourier(im)

% make sure there are an odd number of pixels
if iseven(size(im,1)),im = im(1:end-1,:);end
if iseven(size(im,2)),im = im(:,1:end-1);end

% take fourier transform of image
imf = fft2(im);

% get input dimensions
d.originalDims = size(im);

% get one half of fourier image
imfh = fftshift(imf);
imfh = imfh(1:d.originalDims(1),1:ceil(d.originalDims(2)/2));

% extract dc form half fourier image
d.dc = imfh(ceil(d.originalDims(1)/2),end);
halfFourier = imfh(1:(prod(size(imfh))-ceil(d.originalDims(1)/2)));

d.mag = abs(halfFourier);
d.phase = angle(halfFourier);
d.n = length(d.phase);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    reconstructFromHalfFourier    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function im = reconstructFromHalfFourier(d)

d.halfFourier = d.mag.*(cos(d.phase)+i*sin(d.phase));

% first make the last column of the half fourier space which includes
% the dc and should have the frequency components replicated corectly
halfFourier = [d.halfFourier d.dc];
halfFourier(end+1:end+floor(d.originalDims(1)/2)) = conj(d.halfFourier(end:-1:end-floor(d.originalDims(1)/2)+1));
halfFourier = reshape(halfFourier,d.originalDims(1),ceil(d.originalDims(2)/2));

% replicate the frequency components to make the negative frequencies which
% are the complex conjugate of the positive frequncies
halfFourier2 = fliplr(flipud(conj(halfFourier(:,1:floor(d.originalDims(2)/2)))));
im = ifft2(ifftshift([halfFourier halfFourier2]));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    contrastNormalize    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function im = contrastNormalize(im)

% image max/min
immax = max(im(:));
immin = min(im(:));

% normalize to range of 0:1
im = 255*(im-immin)/(immax-immin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    loadNormalizedImages    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = loadNormalizedImages(dirname,varargin)

if nargin == 0
  help loadNormalizedImages;
  return
end

% get variable arguments
width=[];
height=[];
dispFig=[];
keepAspectRatio=[];
getArgs(varargin,{'height=320','width=240','dispFig=0','keepAspectRatio=0'});

% check directory
d = [];
if ~isdir(dirname)
  disp(sprintf('(loadNormalizedImages) Could not find directory %s',dirname));
  return
end

% size that image will be resampled to
d.width = width;
d.height = height;

% get a listing of directory
d.dirName = dirname;
d.dir = dir(dirname);
d.n = 0;

% load each image
if dispFig,smartfig('loadNormalizedImages','reuse');end
disppercent(-inf,sprintf('(loadNormalizedImages) Loading images for %s',dirname));
d.im = zeros(height,width,length(d.dir));
d.averageMag = 0;d.averageDC = 0;
for i = 1:length(d.dir)
  % get filename
  thisFilename = fullfile(d.dirName,d.dir(i).name);
  % and load if it exists
  if isfile(thisFilename) && ~isempty(imformats(getext(thisFilename)))
    d.n = d.n + 1;
    % read the image
    im = imread(thisFilename);
    % display if called for
    if dispFig,clf;subplot(1,2,1);imagesc(im);axis equal; axis off;end
    % normalize to grayscale and same width height
    im = flipud(imageNormalize(im,d.width,d.height,keepAspectRatio));
    if dispFig,subplot(1,2,2);imagesc(im);colormap(gray);axis equal; axis off;drawnow;end
    % save
    d.im(1:height,1:width,d.n) = im;
    d.filenames{d.n} = thisFilename;
    % get its half fourier image
    d.halfFourier{d.n} = getHalfFourier(d.im(:,:,d.n));
    d.averageMag = d.averageMag + d.halfFourier{d.n}.mag;
    d.averageDC = d.averageDC + d.halfFourier{d.n}.dc;
  end
  disppercent(i/length(d.dir));
end
disppercent(inf);
d.im = d.im(:,:,1:d.n);

% now get average magnitude
d.averageMag = d.averageMag/d.n;
d.averageDC = d.averageDC/d.n;

%%%%%%%%%%%%%%%%%%%%%%%%
%    imageNormalize    %
%%%%%%%%%%%%%%%%%%%%%%%%
function im = imageNormalize(im,width,height,keepAspectRatio)

% get image dimensions
imdim = size(im);

% first convert to grayscale
if length(imdim > 2)
  im = mean(im,3);
end

% get the image coordinates of the image
[x y] = meshgrid(0:1/(imdim(2)-1):1,0:1/(imdim(1)-1):1);


if keepAspectRatio
  % get aspect ratios
  aspectRatioIn = imdim(2)/imdim(1);
  aspectRatioOut = (height/width);

  % now resample to the same dimensions, making sure to keep the same
  % aspect ratio (this will cause some cropping of the image in the
  % appropriate dimension if your aspect ratios do not match).

  % set the image coordinates of the output
  if (aspectRatioOut > aspectRatioIn)
    minX = (1-aspectRatioIn/aspectRatioOut)/2;
    maxX = 1-minX;
    [xi yi] = meshgrid(minX:(maxX-minX)/(width-1):maxX,0:1/(height-1):1);
  else
    minY = (1-aspectRatioOut/aspectRatioIn)/2;
    maxY = 1-minY;
    [xi yi] = meshgrid(0:1/(width-1):1,minY:(maxY-minY)/(height-1):maxY);
  end
else
  [xi yi] = meshgrid(0:1/(width-1):1,0:1/(height-1):1);
end
  
% interpolate the image
im = interp2(x,y,im,xi,yi,'cubic');


%%%%%%%%%%%%%%%%%%%%%%
%%   makeMixImage   %%
%%%%%%%%%%%%%%%%%%%%%%
function tex = makeMixImage(scrambleFactor,image1,image2,dc,mag)

% create the mixture image
mixImage.originalDims = image1.originalDims;
mixImage.dc = dc;
mixImage.mag = mag;
phase1 = image1.phase;
phase2 = image2.phase;

% scramble the phases together
scrambleNum = floor(scrambleFactor*length(phase1));
phaseValues = randperm(length(phase1));
phaseValues1 = phaseValues(1:scrambleNum);
phaseValues2 = phaseValues(scrambleNum+1:end);
mixImage.phase = phase1;
mixImage.phase(phaseValues2) = phase2(phaseValues2);

% create the texture
tex = mglCreateTexture(contrastNormalize(reconstructFromHalfFourier(mixImage)));

