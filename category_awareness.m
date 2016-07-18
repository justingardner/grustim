function [ myscreen ] = category_awareness( varargin )
%CAT_AWE Testing category discrimination at low visibility conditions
%   This is dan's experiment for Berlin (fall 2016). Three object
%   categories are shown across trials and participants simply identify
%   the current category. Each stimulus is 500 ms (50 ms stim + 450 ms
%   mask). 
%
%   There is also a localizer mode ('localizer=1') which just runs a
%   category localizer with the fixation task.
%
%   There is a staircase mode ('staircase=1') which runs three independent
%   staircases for each category identifying what % phase randomization
%   gets the subj to 3/6, 4/6, and 5/6 % correct (2/6 being chance for
%   3AFC).
%
%   If no staircases were calculated the function will default to some
%   arbitrary values.
%
%   Stimuli are shown around fixation in a 10x10 degree view. The remaining
%   screen background is phase randomized noise. 

global stimulus images

%% Initialize Variables

% add arguments later
localizer = 0;
staircase = 0;
scan = 0;
plots = 0;
getArgs(varargin,{'localizer=0','staircase=0','scan=0','plots=0'});
stimulus.localizer = localizer;
stimulus.dostaircase = staircase;
stimulus.scan = scan;
stimulus.plots = plots;
clear localizer staircase scan

stimulus.counter = 1; % This keeps track of what "run" we are on.

%% Setup Screen

if stimulus.scan
    myscreen = initScreen('fMRIprojFlex');
else
    myscreen = initScreen('VPixx');
end

myscreen.background = 0.5;

%% Open Old Stimfile
stimulus.initStair = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/category_awareness/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/category_awareness/%s/1*mat',mglGetSID));

    if length(files) >= 1
        fname = files(end).name;
        
        s = load(sprintf('~/data/category_awareness/%s/%s',mglGetSID,fname));
        % copy staircases and run numbers
        stimulus.staircase = s.stimulus.staircase;
        stimulus.counter = s.stimulus.counter + 1;

        clear s;
        stimulus.initStair = 0;
        disp(sprintf('(cohcon) Data file: %s loaded.',fname));
    end
end
disp(sprintf('(cohcon) This is run #%i',stimulus.counter));

%% Initialize Stimulus

myscreen = initStimulus('stimulus',myscreen);

if stimulus.scan
    stimulus.responseKeys = [2 1]; % corresponds to NOMATCH, MATCH
else
    stimulus.responseKeys = [2 1]; % corresponds to  NOMATCH, MATCH
end

stimulus.colors.black = [0 0 0];
stimulus.colors.white = [1 1 1];
stimulus.colors.green = [0 1 0];
stimulus.colors.red = [1 0 0];

%% Setup Task

task{1}{1}.waitForBacktick = 1;

stimulus.curTrial = 0;

stimulus.seg.warn = 1;
stimulus.seg.stim = 2;
stimulus.seg.mask = 3;
stimulus.seg.ISI = 4;
stimulus.seg.resp = 5;
stimulus.seg.ITI = 6;
task{1}{1}.segmin = [0.15 0.050 0.450 0.4 1.2 .5];
task{1}{1}.segmax = [0.15 0.050 0.450 0.4 1.2 .5];

task{1}{1}.synchToVol = [0 0 0 0 0];
if stimulus.scan
    task{1}{1}.synchToVol(stimulus.seg.ITI) = 1;
end
task{1}{1}.getResponse = [0 0 0 0 0]; task{1}{1}.getResponse(stimulus.seg.resp)=1;
task{1}{1}.parameter.category = [1 2 3]; % which category is shown on this trial
task{1}{1}.parameter.match = [0 1]; % whether this is a "match" trial or not
task{1}{1}.random = 1;
task{1}{1}.numTrials = inf;

%% Tracking

% these are variables that we want to track for later analysis.
task{1}{1}.randVars.calculated.correct = nan;
task{1}{1}.randVars.calculated.phase = nan; % will be 1->5 for 0% to 100% phase

%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init staircase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stimulus.nostairperf = [3/6 4/6 5/6]; % these are the performances we want to use when no staircase is running
if stimulus.dostaircase
    % use a range of phases
    images.phases = repmat(0:.1:.9,3,1);
    if stimulus.initStair
        % We are starting our staircases from scratch
        disp(sprintf('(cohcon) Initializing staircases'));
        stimulus = initStaircase(stimulus);
    else
        disp('(cohcon) Re-using staircase from previous run...');
    end
else
    if stimulus.initStair
        disp('(cat_awe) No staircases present, please run staircase=1 mode first');
        return
    else
        disp('(cat_awe) Loading staircases and predicting phases from weibull fits.');
        out(1) = doStaircase('threshold',stimulus.staircase{1},'type','weibull','dispFig=1','gamma=1/2');
        out(2) = doStaircase('threshold',stimulus.staircase{2},'type','weibull','dispFig=0','gamma=1/2');
        out(3) = doStaircase('threshold',stimulus.staircase{3},'type','weibull','dispFig=0','gamma=1/2');
        testPhases = zeros(3,length(stimulus.nostairperf));
        for cat = 1:3
            for pi = 1:length(stimulus.nostairperf)
                x = out(cat).fit.x;
                y = out(cat).fit.y;
                pos = find(y>stimulus.nostairperf(pi),1);
                if ~isempty(pos)
                    testPhases(cat,pi) = interp1(1:size(stimulus.phases,2),stimulus.phases,x(pos));
                else
                    testPhases(cat,pi) = 1;
                end
            end
        end
        stop = 1;
        stimulus.phases = testPhases;
    end
end

%% Load images
images = myInitStimulus(images);
stimulus.categories = images.categories;
stimulus.phases = images.phases;

%% EYE CALIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if ~stimulus.scan
%     myscreen = eyeCalibDisp(myscreen);
% end

%% Get Ready...
% clear screen    
mglWaitSecs(1);
mglClearScreen(0.5);
if stimulus.scan        
    mglTextDraw('DO NOT MOVE',[0 1.5]);
end
mglFlush

% let the user know
disp(sprintf('(cohcon) Starting run number: %i.',stimulus.counter));
% if stimulus.unattended
myscreen.flushMode = 1;
% end

%% Main Task Loop

phaseNum = 1;
% Again, only one phase.
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    % update the task
    [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);
    % flip screen
    myscreen = tickScreen(myscreen,task);
end

% task ended
mglClearScreen(0.5);
mglTextDraw('Run complete...',[0 0]);
mglFlush
myscreen.flushMode = 1;
mglWaitSecs(1);

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

if stimulus.plots
    disp('(cohcon) Displaying plots');
    dispInfo(stimulus);
end

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startTrialCallback(task,myscreen)
%%

global stimulus images

if stimulus.scan
    task.thistrial.seglen(end) = 1.05^(rand*30+15);
end

stimulus.curTrial = stimulus.curTrial + 1;

myscreen.flushMode = 0;

% set the current image

[task.thistrial.phase, stimulus.staircase{task.thistrial.category}] = doStaircase('testValue',stimulus.staircase{task.thistrial.category});

stimulus.live.imgNum = randi(50);

disp(sprintf('Trial %i Category: %s Phase %2.2f',stimulus.curTrial,stimulus.categories{task.thistrial.category},images.phases(task.thistrial.phase)));


oneof = 1:3;
oneof = oneof(oneof~=task.thistrial.category);
stimulus.live.matchCompare = oneof(randi(2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)
%%

global stimulus

stimulus.live.mask=0;
stimulus.live.obj=0;
stimulus.live.match=0;

switch task.thistrial.thisseg
    case stimulus.seg.ITI
        stimulus.live.fixColor = stimulus.colors.black;
    case stimulus.seg.warn
        stimulus.live.fixColor = stimulus.colors.white;
    case stimulus.seg.stim
        stimulus.live.fixColor = stimulus.colors.black;
        stimulus.live.obj = 1;
    case stimulus.seg.mask
        stimulus.live.fixColor = stimulus.colors.black;
        stimulus.live.mask = 1;
    case stimulus.seg.ISI
        stimulus.live.fixColor = stimulus.colors.black;
    case stimulus.seg.resp
        stimulus.live.match = 1;
        stimulus.live.fixColor = stimulus.colors.white;
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus images

mglClearScreen(0.5);

if ~stimulus.live.match ||  task.thistrial.gotResponse>0
    upFix(stimulus);
else
    upMatch(task,stimulus);
end
if stimulus.live.mask
    mglBltTexture(images.tex{task.thistrial.category}(stimulus.live.imgNum,1),[0 0]);
elseif stimulus.live.obj
    mglBltTexture(images.tex{task.thistrial.category}(stimulus.live.imgNum,task.thistrial.phase),[0 0]);
end

function upMatch(task,stimulus)
mglTextSet('Helvetica',32,[1 1 1 1],0,0,0,0,0,0,0);
if task.thistrial.match
    mglTextDraw(stimulus.categories{task.thistrial.category},[0 0]);
else
    mglTextDraw(stimulus.categories{stimulus.live.matchCompare},[0 0]);
end
mglFlush;


function upFix(stimulus)
%%
mglFixationCross(1.5,1.5,stimulus.live.fixColor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus

responseText = {'Incorrect','Correct'};
fixColors = {stimulus.colors.red,stimulus.colors.green};

if any(task.thistrial.whichButton == stimulus.responseKeys)
    if task.thistrial.gotResponse == 0
        task.thistrial.correct = task.thistrial.whichButton == stimulus.responseKeys(task.thistrial.match+1);
        disp(sprintf('Subject pressed %i: %s',task.thistrial.whichButton,responseText{task.thistrial.correct+1}));
        % Store whether this was correct
        stimulus.live.fixColor = fixColors{task.thistrial.correct+1};
        
        stimulus.staircase{task.thistrial.category} = doStaircase('update',stimulus.staircase{task.thistrial.category},task.thistrial.correct);
    else
        disp(sprintf('(cohcon) Subject responded multiple times: %i',task.thistrial.gotResponse+1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              HELPER FUNCTIONS                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
%    initStaircase     %
%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initStaircase(stimulus)
global images
%%
% we're going to be fucking organized this time and put all the staircases
% in one place.... duh.
stimulus.staircase = cell(1,3);

stimulus.staircase{1} = doStaircase('init','fixed','fixedVals',1:length(images.phases));
stimulus.staircase{2} = stimulus.staircase{1};
stimulus.staircase{3} = stimulus.staircase{1};

%%%%%%%%%%%%%%%%%%%%%%%
%    dispInfo    %
%%%%%%%%%%%%%%%%%%%%%%%
function dispInfo(stimulus)
%%
out = doStaircase('threshold',stimulus.staircase{1},'type','weibull','dispFig=1','gamma=1/2');
disp(sprintf('Threshold for category %s = %2.2f%%',stimulus.categories{1},interp1(1:size(stimulus.phases,2),100*stimulus.phases,out.threshold)));
%%
out = doStaircase('threshold',stimulus.staircase{2},'type','weibull','dispFig=1','gamma=1/2');
disp(sprintf('Threshold for category %s = %0.2f%%',stimulus.categories{2},interp1(1:size(stimulus.phases,2),100*stimulus.phases,out.threshold)));

%%
out = doStaircase('threshold',stimulus.staircase{3},'type','weibull','dispFig=1','gamma=1/2');
disp(sprintf('Threshold for category %s = %0.2f%%',stimulus.categories{3},interp1(1:size(stimulus.phases,2),100*stimulus.phases,out.threshold)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus)

imageDir = '~/proj/grustim/images/category_awareness/';
categories = {'building','car','flower'};
stimulus.categories = categories;
keepAspectRatio=1;
stimulus.widthPix = 512;
stimulus.heightPix = 512;
% make sure widht and height are odd
if iseven(stimulus.widthPix), stimulus.widthPix = stimulus.widthPix-1;end
if iseven(stimulus.heightPix), stimulus.heightPix = stimulus.heightPix-1;end

% check whether images are loaded
averageN = 0;
if ~isfield(stimulus,'imagesLoaded') || (~stimulus.imagesLoaded) 
  stimulus.nCategories = length(categories);

  % keep the averageMag and averageDc so that we can normalize images
  stimulus.averageMag = 0;
  stimulus.averageDC = 0;

  for i = 1:stimulus.nCategories
    if ~any(strcmp(categories{i},{'scramble','blank','gray'}))
      % load images
      stimulus.raw{i} = loadNormalizedImages(fullfile(imageDir,categories{i}),'width',stimulus.widthPix,'height',stimulus.heightPix,'keepAspectRatio',keepAspectRatio);

      % make sure we opened ok
      if isempty(stimulus.raw{i})
	disp(sprintf('(cat_awe) Could not load images; %s',categories{i}));
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
  disp(sprintf('(cat_awe) Stimulus already initialized'));
end


% count how many images
nImages = [];
for i = 1:stimulus.nCategories
  if stimulus.raw{i}.n > 0
    % keep number of images so that we can make scrambles with equal number of images
    nImages(end+1) = stimulus.raw{i}.n;
  end
end
%%
% create textures of all images
disppercent(-inf,'(cat_awe) Converting images to textures');
for iCategory = 1:stimulus.nCategories
    for iImage = 1:stimulus.raw{iCategory}.n
      thisImage = stimulus.raw{iCategory}.halfFourier{iImage};
      thisImage.mag = stimulus.averageMag;
      thisImage.dc = stimulus.averageDC;
      original = reconstructFromHalfFourier(thisImage);
      original = flipud(original);
      thisImage.phase = rand(size(thisImage.mag))*2*pi;
      thisImage.originalDims = [stimulus.heightPix stimulus.widthPix];
      scramble = reconstructFromHalfFourier(thisImage);
      for z = 1:size(stimulus.phases,2)
          percSignal = stimulus.phases(iCategory,z);
          stimulus.tex{iCategory}(iImage,z) = mglCreateTexture(contrastNormalize(original .* percSignal + scramble .* (1-percSignal)));
      end
      disppercent(calcPercentDone(iCategory,stimulus.nCategories,iImage,stimulus.raw{iCategory}.n));
    end
end
disppercent(inf);

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

d.halfFourier = d.mag.*(cos(d.phase)+1i*sin(d.phase));

% first make the last column of the half fourier space which includes
% the dc and should have the frequency components replicated corectly
halfFourier = [d.halfFourier d.dc];
halfFourier(end+1:end+floor(d.originalDims(1)/2)) = conj(d.halfFourier(end:-1:end-floor(d.originalDims(1)/2)+1));
halfFourier = reshape(halfFourier,d.originalDims(1),ceil(d.originalDims(2)/2));

% replicate the frequency components to make the negative frequencies which
% are the complex conjugate of the positive frequncies
halfFourier2 = rot90(conj(halfFourier(:,1:floor(d.originalDims(2)/2))),2);
im = ifft2(ifftshift([halfFourier halfFourier2]));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    contrastNormalize    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function im = contrastNormalize(im)

% new version (July 2011) to keep luminance and rms contrast
% constant across images. Scaling factor 0.7 chosen for a specific set
% of stimuli (the set of images that was in grustim/images/ObjLocImages/
% in June 2011)

 im = im*0.7;
 im(find(im>255)) = 255;
 im(find(im<0)) = 0;

% original version of contrastNormalize
%
% image max/min
% immax = max(im(:));
% immin = min(im(:));
%
% normalize to range of 0:1
% im = 255*(im-immin)/(immax-immin);

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
getArgs(varargin,{'height=512','width=512','dispFig=0','keepAspectRatio=0'});

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
  
% uncomment the following line if you want to rescale the aspect ratio
% [xi yi] = meshgrid(0:1/(width-1):1,0:1/(height-1):1);

% interpolate the image
im = interp2(x,y,im,xi,yi,'cubic');

