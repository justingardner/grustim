function [ myscreen ] = fLoc( varargin )
%
% Functional Localizer to define Category-Selective Cortical Regions
%  
%  Presents characters, bodies, faces, places, and objects in a block design.
%
%   Citation
%       Stigliani, A., Weiner, K. S., & Grill-Spector, K. (2015). 
%       Temporal processing capacity in high-level visual cortex is 
%       domain specific. Journal of Neuroscience, 35(36), 12412-12424.
%
%   For more info: http://vpnl.stanford.edu/fLoc/
%
%   Download stimuli from above URL and store at ~/proj/fLoc_stimuli
% 
%  Usage: fLoc(varargin)
%  Authors: Akshay Jagadeesh
%  Date: 09/05/2020
%

global stimulus

stimulus = struct;

%% Initialize Variables

% add arguments later
scan = 0;
getArgs(varargin,{'scan=0', 'oban=1'}, 'verbose=1');
stimulus.scan = scan;
stimulus.oban = oban;
clear scan;

%% Setup Screen
if stimulus.scan || stimulus.oban;
  myscreen = initScreen('fMRIprojFlex');
else
  myscreen = initScreen('VPixx2');
end

%% Initialize Stimulus
% set background to grey
myscreen.background = 0.5;
myscreen = initStimulus('stimulus',myscreen);
  
% Set response keys
stimulus.responseKeys = [1 2 3 4]; 

% set colors
stimulus.colors.white = [1 1 1];
stimulus.colors.black = [0 0 0];
stimulus.colors.red = [1 0 0];
stimulus.colors.green = [0 1 0];
stimulus.colors.blue = [0 0 1];
stimulus.live.fixColor = stimulus.colors.blue;
stimulus.live.cueColor = stimulus.colors.black;

%%%%%%%%%%%%% SETUP TASK %%%%%%%%%%%%%%%%%

stimulus.curTrial(1) = 0;

stimulus.imSize = 10;

%% Define stimulus timing
% Task important variables
stimulus.categories = {'characters', 'bodies', 'faces', 'places', 'objects'};
stimulus.subcategories = struct();
stimulus.subcategories.characters = {'word', 'number'};
stimulus.subcategories.bodies = {'body', 'limb'};
stimulus.subcategories.faces = {'adult', 'child'};
stimulus.subcategories.places = {'corridor', 'house'};
stimulus.subcategories.objects = {'car', 'instrument'};

stimulus.num_samples = 20;

%% Define stimulus directory
stimulus.stimDir = '~/proj/fLoc_stimuli';

%% Preload images
mask = imread('~/proj/TextureSynthesis/stimuli/Flattop8.tif');
stimulus.images = struct();

% load texture and noise samples
disppercent(-inf, 'Preloading images');
for i = 1:length(stimulus.categories)
  cat = stimulus.categories{i};
  for j = 1:2
    subcat = stimulus.subcategories.(cat){j};
    for k = 1:stimulus.num_samples
      img = imread(sprintf('%s/%s/%s-%i.jpg', stimulus.stimDir, subcat, subcat, k));
      stimulus.images.(sprintf('%s_%i', subcat, k)) = genTexFromIm(img);
    end
  end
  disppercent(i / length(stimulus.categories));
end
disppercent(inf);
clear img

%%%%%%%%%%%%% TASK %%%%%%%%%%%%%%%%%
task{1} = struct;
task{1}.waitForBacktick = 1;
task{1}.segmin = zeros(1,8) + 0.5;
task{1}.segmax = zeros(1,8) + 0.5;
stimulus.seg = {};
for i = 1:8
  stimulus.seg.(sprintf('stim%i', i)) = i;
end

% Trial parameters
task{1}.synchToVol = zeros(size(task{1}.segmin));
task{1}.getResponse = zeros(size(task{1}.segmin));
task{1}.numTrials = 75;
task{1}.random = 1;

if stimulus.scan
  task{1}.synchToVol(end) = 1;
  % Shorten the last segment to account for synchtovol
  task{1}.segmin(end) = max(0, task{1}.segmin(end) - 0.100);
  task{1}.segmax(end) = max(0, task{1}.segmax(end) - 0.100);
end

% Select trials
task{1}.parameter.category = 1:(length(stimulus.categories)+1); % 1: characters, 2: bodies, 3: faces, 4: places, 5: objects, 6: blank
task{1}.parameter.subcategory = 1:2;
task{1}.randVars.calculated.category_str = NaN;
task{1}.randVars.calculated.subcategory_str = NaN;
task{1}.randVars.calculated.blank = NaN;

%% Init task
[task{1}, myscreen] = initTask(task{1}, myscreen, @startSegmentCallback, @screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main Task Loop

mglClearScreen(0.5); 
upFix(stimulus);
mglFlush;
mglClearScreen(0.5); 
upFix(stimulus);

phaseNum = 1;
% Again, only one phase.
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
  % update the task
  [task, myscreen, phaseNum] = updateTask(task, myscreen, phaseNum);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% task ended
mglClearScreen(0.5);
mglTextSet([],32,stimulus.colors.white);
% get count
mglTextDraw('Please wait',[0 0]);
mglFlush
myscreen.flushMode = 1;

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = startTrialCallback(task,myscreen)

global stimulus

if task.thistrial.category <= length(stimulus.categories)
  task.thistrial.category_str = stimulus.categories{task.thistrial.category};
  task.thistrial.subcategory_str = stimulus.subcategories.(task.thistrial.category_str){task.thistrial.subcategory};
  task.thistrial.blank = 0;
else
  task.thistrial.blank = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = startSegmentCallback(task, myscreen)

global stimulus

if task.thistrial.category <= length(stimulus.categories)
  % If non-blank trial
  seg_sample = randi(stimulus.num_samples);
  seg_image = stimulus.images.(sprintf('%s_%i', task.thistrial.subcategory_str, seg_sample));
  for i = 1:2
    mglClearScreen(0.5);
    mglBltTexture(seg_image, [0, 0, stimulus.imSize, stimulus.imSize]);
    upFix(stimulus, stimulus.colors.black);
    mglFlush;
  end
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus

% Select which stimulus to display as a function of time since seg start

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus

validResponse = any(task.thistrial.whichButton == stimulus.responseKeys);

if validResponse
  if stimulus.live.gotResponse==0
    task.thistrial.response = task.thistrial.whichButton - 10;
  else
    disp(sprintf('Subject responded multiple times: %i',stimulus.live.gotResponse));
  end
  %task = jumpSegment(task);
else
  disp(sprintf('Invalid response key. Subject pressed %d', task.thistrial.whichButton));
  task.thistrial.response = -1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              HELPER FUNCTIONS                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Draws a cross at center of the screen of color fixColor
function upFix(stimulus, fixColor)
if ieNotDefined('fixColor')
  fixColor = stimulus.live.fixColor;
end
mglFixationCross(1,1,fixColor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Turns image into a texture
function tex = genTexFromIm(im, mask)
r = flipud(im);

% Resize images to 256
if size(r,1) ~= 256;
  r = imresize(r, 256/size(r,1));
end

% Make sure they have 3 dimensions (even if grayscale)
if size(r,3) == 1
  r = cat(3, r, r, r);
end

% If a mask is passed in, apply as an alpha mask.
if ieNotDefined('mask')
  r(:,:,4) = 255;
else
  r(:,:,4) = mask(:,:,1);
end
% mgl requires the first dimension to be the RGBA dimension
rP = permute(r, [3 2 1]);
tex = mglCreateTexture(rP);  

