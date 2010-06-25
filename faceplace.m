% facePlace.m
%
%        $Id:$ 
%      usage: facePlace(<'faceDir=faceDir'>,<'placeDir=placeDir'>)
%         by: justin gardner
%       date: 06/14/10
%    purpose: get psychometric function for faces vs places
%
function retval = faceplace(varargin)

% get arguments
faceDir = [];
placeDir = [];
imageDir = [];
stimfile = [];
getArgs(varargin,{'faceDir=faces','placeDir=houses','imageDir=~/proj/faceplace/FaceHouseStim','stimfile=[]'});

% just display data if called with stimfile name
if ~isempty(stimfile),dispStimfile(stimfile);return;end

% initalize the screen
myscreen.background = 'gray';
myscreen = initScreen(myscreen)
% init the stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);

% load images
stimulus = myInitStimulus(stimulus,myscreen,faceDir,placeDir,imageDir);

% set up task
task{1}.waitForBacktick = 0;
task{1}.segmin = [0.5 2];
task{1}.segmax = [0.5 2];
task{1}.getResponse = [0 1];
% fix: enter the parameter of your choice
task{1}.parameter.scrambleFactor = 0.4:0.05:0.8;%0:0.1:1;%[0 0.4 0.5 0.6 1];
task{1}.numBlocks = 40;
%task{1}.parameter.scrambleFactor = [0 1];
%task{1}.numBlocks = 10;
%task{1}.parameter.scrambleFactor = [0.7 0.8];
task{1}.parameter.scrambleFactor = [0.6];
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

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

dispPsychometricFunction(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%
%    dispStimfile    %
%%%%%%%%%%%%%%%%%%%%%%
function dispStimfile(filename)

filename = setext(filename,'mat');
if ~isfile(filename)
  disp(sprintf('(faceplace) Could not find file %s',filename));
  return
end

stimfile = load(filename);
if ~isfield(stimfile,'myscreen') || ~isfield(stimfile,'task')
  disp(sprintf('(faceplace) File %s does not contain myscreen/task',filename));
  return
end

% display the data
dispPsychometricFunction(stimfile.myscreen,stimfile.task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    dispPsychometricFunction    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispPsychometricFunction(myscreen,task)

if ~isfield(task{1},'trialnum') || (task{1}.trialnum <= length(task{1}.parameter.scrambleFactor))
  disp(sprintf('(faceplace:dispPsychometricFunction) Not enough data to display'));
  return
end

e = getTaskparameters(myscreen,task);
p = task{1}.parameter.scrambleFactor;
for i = 1:length(p)
  thisResponse = denan(e.response(e.parameter.scrambleFactor == p(i)));
  n(i) = length(thisResponse);
  nChoice2(i) = sum(thisResponse==2);
  percentChoices2(i) = nChoice2(i)/n(i);
  percentChoices2Error(i) = sqrt((percentChoices2(i) * (1-percentChoices2(i)))/n(i));
end
weibull = fitweibull(p,nChoice2,n,0,0,0,0);
weibull.bias = interp1(weibull.y,weibull.x,0.5);

smartfig('faceplace','reuse');
clf;
myerrorbar(p,percentChoices2,'yError',percentChoices2Error,'Symbol=o');
hold on
plot(weibull.x,weibull.y,'k-');
hline(0.5);vline(weibull.bias);
xlabel('face <- scramble factor -> place');
ylabel('percent choices place');
title(sprintf('%s: %i blocks\n(bias: %f beta: %f lambda: %f)',myscreen.starttime,task{1}.blocknum,weibull.bias,weibull.fitparams(2),weibull.fitparams(3)));
keyboard

% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = stimStartSegmentCallback(task, myscreen)

global stimulus;

if task.thistrial.thisseg == 1

  % select which image to show
  stimulus.faceNum(task.trialnum) = ceil(stimulus.face.n*rand);
  stimulus.placeNum(task.trialnum) = ceil(stimulus.place.n*rand);

  % get its half fourier representation
  faceIm  = stimulus.face.halfFourier{stimulus.faceNum(task.trialnum)};
  placeIm  = stimulus.place.halfFourier{stimulus.placeNum(task.trialnum)};

  % now create a hybrid image which contains scambleFactor percent of the 
  % phase information from the place image
  scrambledPhases = randperm(faceIm.n);
  scrambledPhases = scrambledPhases(1:floor(faceIm.n*task.thistrial.scrambleFactor));

  faceIm.phase(scrambledPhases) = placeIm.phase(scrambledPhases);
  %faceIm.phase(scrambledPhases) = rand(1,length(scrambledPhases))*2*pi;

  % get the average
  faceIm.mag = stimulus.averageMag;

  % delete old image
  if ~isempty(stimulus.imageTex)
    mglDeleteTexture(stimulus.imageTex);
  end

  % get the current image
  stimulus.imageTex = mglCreateTexture(flipud(contrastNormalize(reconstructFromHalfFourier(faceIm))));
  mglBltTexture(stimulus.imageTex,[0 0 16 24]);
  mglFixationCross(1,1,[0 1 1]);
  myscreen.flushMode = 1;
else
  myscreen.flushMode = 1;
  mglClearScreen;
  mglFixationCross(1,1,[0 1 1]);
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

mglClearScreen;
mglFixationCross(1,1,[1 0 1]);
myscreen.flushMode = 1;
disp(sprintf('scrambleFactor: %f choice: %i',task.thistrial.scrambleFactor,task.thistrial.whichButton));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen,faceDir,placeDir,imageDir)

% make fully qualified paths
faceDir = fullfile(imageDir,faceDir);
placeDir = fullfile(imageDir,placeDir);

% fix: add stuff to initalize your stimulus
if isfield(stimulus,'imagesLoaded') && (stimulus.imagesLoaded)
  disp(sprintf('(faceplace) Stimulus already initialized'));
  return
end

% face dir
stimulus.face = loadNormalizedImages(faceDir);
if isempty(stimulus.face)
  disp(sprintf('(facePlace) Could not load face images; %s',faceDir));
  keyboard
end

% place dir
stimulus.place = loadNormalizedImages(placeDir);
if isempty(stimulus.place)
  disp(sprintf('(facePlace) Could not load face images; %s',faceDir));
  keyboard
end

% get average mag
stimulus.averageMag = (stimulus.place.averageMag + stimulus.face.averageMag)/2;

stimulus.imagesLoaded = 1;
stimulus.faceNum = [];
stimulus.placeNum = [];
stimulus.imageTex = [];

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
getArgs(varargin,{'height=240','width=320','dispFig=0'});

% check directory
d = [];
if ~isdir(dirname)
  disp(sprintf('(loadNormalizedImages) Could not find directory %s',dirname));
  return
end

% size that image will be resampled to
d.width = 320;
d.height = 240;

% get a listing of directory
d.dirName = dirname;
d.dir = dir(dirname);
d.n = 0;

% load each image
if dispFig,smartfig('loadNormalizedImages','reuse');end
disppercent(-inf,sprintf('(loadNormalizedImages) Loading images for %s',dirname));
d.im = zeros(width,height,length(d.dir));
d.averageMag = 0;
for i = 1:length(d.dir)
  % get filename
  thisFilename = fullfile(d.dirName,d.dir(i).name);
  % and load if it exists
  if isfile(thisFilename) && ~isempty(imformats(getext(thisFilename)))
    d.n = d.n + 1;
    % read the image
    im = imread(thisFilename);
    % normalize to grayscale and same width height
    im = imageNormalize(im,d.width,d.height);
    if dispFig,clf;imagesc(im);drawnow;colormap(gray);axis equal; axis off;end
    % save
    d.im(1:width,1:height,d.n) = im;
    d.filenames{d.n} = thisFilename;
    % get its half fourier image
    d.halfFourier{d.n} = getHalfFourier(d.im(:,:,d.n));
    d.averageMag = d.averageMag + d.halfFourier{d.n}.mag;
  end
  disppercent(i/length(d.dir));
end
disppercent(inf);
d.im = d.im(:,:,1:d.n);

% now get average magnitude
d.averageMag = d.averageMag/d.n;

%%%%%%%%%%%%%%%%%%%%%%%%
%    imageNormalize    %
%%%%%%%%%%%%%%%%%%%%%%%%
function im = imageNormalize(im,width,height)

% get image dimensions
imdim = size(im);

% first convert to grayscale
if length(imdim > 2)
  im = mean(im,3);
end

% now resample to the same dimensions
[x y] = meshgrid(0:1/(imdim(2)-1):1,0:1/(imdim(1)-1):1);
[xi yi] = meshgrid(0:1/(height-1):1,0:1/(width-1):1);
im = interp2(x,y,im,xi,yi,'cubic');

