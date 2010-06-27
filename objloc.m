% objloc.m
%
%        $Id:$ 
%      usage: objloc(<'faceDir=faceDir'>,<'placeDir=placeDir'>)
%         by: justin gardner
%       date: 06/25/10
%    purpose: object localizer scan
%
% options:
% objloc('imageDir=somedir') 
%   the directory which stores images
% objloc('categories',{'faces','houses','scramble','gray'}) % categories
%   of images, should be directories of image under imageDir for everything
%   except scramble and gray which are created by the program
% objloc('categoryWeight=[4 2 1 1]')
%    This is an array that specifies how you want to weight each category in
%    terms of how frequently it will be shown. In this case, category 1 will
%    be shown 4 times for every time category 2 is shown 2 and category 3
%    and 4 shown 1 time.
% objloc('keepAspectRatio=0')
%    defaults to zero - keeps the aspect ratio of the original
% objloc('widthPix=180','heightPix=180','widthDeg','heightDeg')
%    image sizes
% objloc('repeatFreq=0.1')
%    frequency with which to repeat images
% objloc('waitForBacktick=1')
%    whether to wait for backtick (for running in scanner)
function retval = objloc(varargin)

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
getArgs(varargin,{'categories',{'faces','houses','scramble','gray'},'imageDir=~/proj/faceplace/FaceHouseStim','dispLoadFig=0','categoryWeight=[]','keepAspectRatio=0','repeatFreq=0.1','waitForBacktick=0','widthPix=180','heightPix=180','widthDeg=18','heightDeg=18'});

% initalize the screen
myscreen.background = 'gray';
myscreen = initScreen(myscreen);
% init the stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);

% setup categoryWeight
if isempty(categoryWeight)
  categoryWeight = ones(1,length(categories));
elseif length(categoryWeight) ~= length(categories)
  disp(sprintf('(objloc) Length of categories: %i does not mach categoryWeight length: %i',length(categories),length(categoryWeight)));
  return
end
categoryNums = [];
for i = 1:length(categoryWeight)
  categoryNums = [categoryNums repmat(i,1,categoryWeight(i))];
end

% load images
stimulus.widthPix = widthPix;
stimulus.heightPix = heightPix;
stimulus.widthDeg = widthDeg;
stimulus.heightDeg = heightDeg;
stimulus.repeatFreq = repeatFreq;
stimulus = myInitStimulus(stimulus,myscreen,categories,imageDir,dispLoadFig,keepAspectRatio);

% set up task
task{1}.waitForBacktick = waitForBacktick;
task{1}.seglen = repmat([0.75 0.25],1,12);
if waitForBacktick,task{1}.seglen(end) = 0.1;end
task{1}.getResponse = ones(1,length(task{1}.seglen));
task{1}.getResponse(1:2) = 0;
task{1}.synchToVol = zeros(1,length(task{1}.seglen));
task{1}.synchToVol(end) = waitForBacktick;
% fix: enter the parameter of your choice
task{1}.parameter.categoryNum = categoryNums;
task{1}.numBlocks = 100;
task{1}.random = 1;

% initialize the task
for phaseNum = 1:length(task)
  [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@stimStartSegmentCallback,@stimDrawStimulusCallback,@responseCallback);
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
  % display what block we are on at the first segment only
  if task.thistrial.thisseg == 1
    disp(sprintf('%i: %s',task.trialnum,stimulus.categories{task.thistrial.categoryNum}));
  end
  % display a random image
  if stimulus.raw{task.thistrial.categoryNum}.n > 0
    % see if we want to do a repeat
    if (task.thistrial.thisseg > 1) && (rand < stimulus.repeatFreq);
      randImageNum = stimulus.randImageNum{task.trialnum}((task.thistrial.thisseg-1)/2);
    else
      % otherwise get a random image number
      randImageNum = ceil(rand*stimulus.raw{task.thistrial.categoryNum}.n);
    end
    stimulus.randImageNum{task.trialnum}((task.thistrial.thisseg+1)/2) = randImageNum;
    mglBltTexture(stimulus.tex{task.thistrial.categoryNum}(randImageNum),[0 0 stimulus.widthDeg stimulus.heightDeg]);;
  end
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
if isodd(task.thistrial.thisseg) && (stimulus.raw{task.thistrial.categoryNum}.n ~= 0)
  randImageNum = stimulus.randImageNum{task.trialnum}(floor((task.thistrial.thisseg+1)/2));
  mglBltTexture(stimulus.tex{task.thistrial.categoryNum}(randImageNum),[0 0 stimulus.widthDeg stimulus.heightDeg]);;
  % check correct
  if randImageNum == stimulus.randImageNum{task.trialnum}(floor((task.thistrial.thisseg-1)/2));
    disp('(objloc) Correct');
    fixColor = [0 1 0];
  else
    disp('(objloc) Incorrect')
    fixColor = [1 0 0];
  end
end
% change fixation color
mglFixationCross(1,1,fixColor);
myscreen.flushMode = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen,categories,imageDir,dispFig,keepAspectRatio)

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

% create textures of all images
disppercent(-inf,'(objloc) Converting images to textures');
for iCategory = 1:stimulus.nCategories
  if ~strcmp(categories{iCategory},'scramble')
    for iImage = 1:stimulus.raw{iCategory}.n
      thisImage = stimulus.raw{iCategory}.halfFourier{iImage};
      thisImage.mag = stimulus.averageMag;
      thisImage.dc = stimulus.averageDC;
      thisImage = reconstructFromHalfFourier(thisImage);
      stimulus.tex{iCategory}(iImage) = mglCreateTexture(contrastNormalize(flipud(thisImage)));
      disppercent(calcPercentDone(iCategory,stimulus.nCategories,iImage,stimulus.raw{iCategory}.n));
    end
  end
end
disppercent(inf);

% create scramble
doScramble = find(strcmp('scramble',categories));
if ~isempty(doScramble)
  disppercent(-inf,'(objloc) Creating scrambles');
  % how many screamble to make
  nScrambles = median(nImages);
  for i = 1:length(doScramble)
    stimulus.raw{doScramble(i)}.n = nScrambles;
    for j = 1:stimulus.raw{doScramble(i)}.n
      thisImage = [];
      thisImage.dc = stimulus.averageDC;
      thisImage.mag = stimulus.averageMag;
      thisImage.phase = rand(size(thisImage.mag))*2*pi;
      thisImage.originalDims = [stimulus.heightPix stimulus.widthPix];
      thisImage = reconstructFromHalfFourier(thisImage);
      stimulus.tex{doScramble(i)}(j) = mglCreateTexture(contrastNormalize(thisImage));
      disppercent(calcPercentDone(i,length(doScramble),j,stimulus.raw{doScramble(i)}.n));
    end
  end
end
disppercent(inf);

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
    im = imageNormalize(im,d.width,d.height,keepAspectRatio);
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




