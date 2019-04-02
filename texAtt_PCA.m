function [ myscreen ] = texAtt_PCA( varargin )
%
% TEXTURE SEARCH 
%  Visual search task using textures 
%
%  Usage: texAtt_PCA(varargin)
%  Authors: Akshay Jagadeesh
%  Date: 02/27/2018
%

global stimulus

stimulus = struct;

%% Initialize Variables

% add arguments later
plots = 0;
noeye = 0;
getArgs(varargin,{'plots=0','noeye=1', 'analyze=0'});
stimulus.plots = plots;
stimulus.noeye = noeye;
clear noeye plots

if analyze
    analyzeData;
    return;
end

%% Stimulus parameters 
%% Open Old Stimfile
stimulus.counter = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/texAtt_PCA/%s',mglGetSID))
  % Directory exists, check for a stimfile
  files = dir(sprintf('~/data/texAtt_PCA/%s/1*mat',mglGetSID));

  if length(files) >= 1
    fname = files(end).name;
    
    s = load(sprintf('~/data/texAtt_PCA/%s/%s',mglGetSID,fname));
    stimulus.counter = s.stimulus.counter + 1;
    clear s;
    disp(sprintf('(texAtt_PCA) Data file: %s loaded.',fname));
  end
end
disp(sprintf('(texAtt_PCA) This is run #%i',stimulus.counter));

%% Setup Screen
myscreen = initScreen('VPixx2');

% set background to grey
myscreen.background = 0.5;

%% Setup missing initial variables
if ~isfield(stimulus,'counter')
  stimulus.counter = 1; % This keeps track of what "run" we are on.
end


%% Initialize Stimulus
myscreen = initStimulus('stimulus',myscreen);
localInitStimulus();
  
% Set response keys
% stimulus.responseKeys = [11 12 13 14];
stimulus.responseKeys = [11 12 13];


% set colors
stimulus.colors.white = [1 1 1];
stimulus.colors.black = [0 0 0];
stimulus.colors.red = [1 0 0];
stimulus.colors.green = [0 1 0];
stimulus.colors.blue = [0 0 1];
stimulus.live.fixColor = stimulus.colors.blue;
stimulus.live.cueColor = stimulus.colors.black;

%%%%%%%%%%%%% SET TASK VARIABLES  %%%%%%%%%%%%%%%%%
stimulus.curTrial(1) = 0;
task{1} = struct;
task{1}.waitForBacktick = 1;

% Define stimulus timing
% inf, Cue, Stimulus, Pause, Response, Feedback
task{1}.segmin = [inf, 1.0, 0.6, 1.2, .2];
task{1}.segmax = [inf, 1.0, 0.6, 1.2, .2];
stimulus.seg = {};
stimulus.seg.fix = 1;
stimulus.seg.cue = 2;
stimulus.seg.stim = 3;
stimulus.seg.response = 4;
stimulus.seg.feedback = 5;

if stimulus.noeye==1
  task{1}.segmin(1) = 0.1;
  task{1}.segmax(1) = 0.1;
end


% Task important variables
%stimulus.imNames = {'balls', 'beansalad', 'biryani','bubbles', 'cherries', 'clouds', 'crowd', 'dahlias',...
%    'fireworks', 'forest', 'fronds', 'leaves', 'noodles', 'paneer','rocks',  'tulips', 'worms', 'zebras'};
% stimulus.imNames = {'balls', 'beansalad', 'biryani', 'bubbles', 'cherries', 'clouds', ...
%         'crowd', 'dahlias', 'fireworks', 'leaves', 'noodles', 'rocks', 'tulips', 'worms', ...
%         'zebras', 'bananas', 'bark', 'bison', 'blossoms', 'blotch', 'braids', 'bricks',...
%         'bubbly', 'bumpy', 'crystals', 'dalmatians', 'ducks', 'face', 'frills', 'fur', ...
%         'galaxy', 'gourds', 'grass', 'honeycomb', 'lace', 'marbled', 'marbles', 'monarchs',...
%         'paisley', 'pears', 'phlox', 'rorschach', 'spiky', 'splotchy', 'stars', 'succulent', 'tiles'};
stimulus.imNames = {'bark', 'bricks', 'drops', 'face', 'fireworks', 'fronds', 'glass', 'paisley', 'rocks', 'tulips'};
% stimulus.layerNames = {'pool1', 'pool2', 'pool4'};
stimulus.layerNames = {'pool2'};
stimulus.stimDir = '~/proj/TextureSynthesis/stimuli/textures/tex_pca';
stimulus.imSize = 6;
stimulus.eccentricity = 12;
% Stimulus and oddball poolSize is same for this task -- var leftover from
% previous task
%stimulus.obPoolSize = '4x4';

stimulus.PCs = {'PC0', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5'};
stimulus.PCDists = {'0'};
stimulus.ob_PCDists = {'-10', '-5', '-2', '2', '5', '10'};
% stimulus.poolSizes = {'1x1','1.25x1.25', '1.5x1.5','1.75x1.75', '2x2', '3x3', '4x4'};
stimulus.poolSizes = {'1x1'};
stimulus.cueEcc = 4;
stimulus.live.mask = imread('~/proj/TextureSynthesis/stimuli/Flattop8.tif');


% Trial parameters
task{1}.parameter.layer = 1:length(stimulus.layerNames);
task{1}.parameter.poolSize = 1:length(stimulus.poolSizes);
task{1}.parameter.PCs = stimulus.PCs;
task{1}.parameter.PCDists = stimulus.PCDists;
task{1}.parameter.ob_PCDists = stimulus.ob_PCDists;
task{1}.parameter.isCueFocal = [0 1]; % Is the cue distributed or focal?

task{1}.synchToVol = zeros(size(task{1}.segmin));
task{1}.getResponse = zeros(size(task{1}.segmin));
task{1}.getResponse(stimulus.seg.response)=1;

% Make numTrials some multiple of number of TrialTypes (in this case, # of attentional conditions x # of layers).
task{1}.numTrials = 144; 
task{1}.random = 1;

% Task variables to be calculated later
task{1}.randVars.calculated.cueSide = NaN; % which side is the target on
task{1}.randVars.calculated.cueside_imFam = NaN; % which image are we showing on the target side
task{1}.randVars.calculated.cueside_oddpos = NaN; % which of the 3 positions is the oddball on.
task{1}.randVars.calculated.cueside_poolSize = NaN;
task{1}.randVars.calculated.cueside_smp = NaN;
task{1}.randVars.calculated.cueside_PC = NaN;

task{1}.randVars.calculated.otherside_smp = NaN;
task{1}.randVars.calculated.otherside_layer = NaN;
task{1}.randVars.calculated.otherside_poolSize = NaN;
task{1}.randVars.calculated.otherside_imFam = NaN; % which image are we showing on the non-target side
task{1}.randVars.calculated.otherside_oddpos = NaN; 
task{1}.randVars.calculated.otherside_PC = NaN;

task{1}.randVars.calculated.ob_PCdist = NaN;
task{1}.randVars.calculated.d_PCdist = NaN;

task{1}.randVars.calculated.correct = NaN;

task{1}.randVars.calculated.detected = 0; % did they see the grating
task{1}.randVars.calculated.dead = 0;
task{1}.randVars.calculated.visible = 1;


%%%%%%%%%%%%%%%% MGL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialize task and begin update loop
[task{1}, myscreen] = initTask(task{1},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);

% Run the eye calibration
myscreen = eyeCalibDisp(myscreen);

% let the user know
disp(sprintf('(texAtt_PCA) Starting run number: %i.',stimulus.counter));

%% Main Task Loop
mglClearScreen(0.5); 
upFix(stimulus);

mglFlush
mglClearScreen(0.5); 
upFix(stimulus);

phaseNum = 1;
while (phaseNum <= length(task)) && ~myscreen.userHitEsc
  % update the task
  [task, myscreen, phaseNum] = updateTask(task,myscreen,phaseNum);
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

task.thistrial.dead = 0;
task.thistrial.detected = 0;
task.thistrial.visible = 1;
task.thistrial.response = NaN;

stimulus.live.gotResponse = 0;
stimulus.curTrial(task.thistrial.thisphase) = stimulus.curTrial(task.thistrial.thisphase) + 1;

% directories
texDir = stimulus.stimDir;

%% Load all 6 images for this trial

% Making imFam & PC distance between oddball and distractors same for both sides
% all of these are the same between cueside and other side, only difference
% is cueing

%Set Cue side task parameters
task.thistrial.cueside_imFam = randi(length(stimulus.imNames));
task.thistrial.cueside_layer = randi(length(stimulus.layerNames));
task.thistrial.cueside_poolSize = randi(length(stimulus.poolSizes));
task.thistrial.cueside_PC = randi(length(stimulus.PCs));

% Get image, layer, and poolSize for cue side
cueside_PC = stimulus.PCs{task.thistrial.cueside_PC};
cueside_layer = stimulus.layerNames{task.thistrial.cueside_layer};
cueside_poolSize = stimulus.poolSizes{task.thistrial.cueside_poolSize};
cueside_imFam = stimulus.imNames{task.thistrial.cueside_imFam};

% PC Distance is the same for both cueside and otherside.
task.thistrial.ob_PCDist = randi(length(stimulus.ob_PCDists));
task.thistrial.d_PCDist = randi(length(stimulus.PCDists));
ob_PCDist = stimulus.ob_PCDists{task.thistrial.ob_PCDist};
d_PCDist = stimulus.PCDists{task.thistrial.d_PCDist};

% Get image, layer, and poolSize for other side 
sidesAreDifferent = 0;
if sidesAreDifferent
    task.thistrial.otherside_imFam = randi(length(stimulus.imNames));
    task.thistrial.otherside_layer = randi(length(stimulus.layerNames));
    task.thistrial.otherside_poolSize = randi(length(stimulus.poolSizes));
    task.thistrial.otherside_PC = randi(length(stimulus.PCs));
else
    task.thistrial.otherside_imFam = task.thistrial.cueside_imFam;
    task.thistrial.otherside_layer = task.thistrial.cueside_layer;
    task.thistrial.otherside_poolSize = task.thistrial.cueside_poolSize;
    task.thistrial.otherside_PC = task.thistrial.cueside_PC;
end

otherside_imFam = stimulus.imNames{task.thistrial.otherside_imFam};
otherside_layer = stimulus.layerNames{task.thistrial.otherside_layer};
otherside_poolSize = stimulus.poolSizes{task.thistrial.otherside_poolSize};
otherside_PC = stimulus.PCs{task.thistrial.otherside_PC};

% Randomly select 3 images to display on each side for this trial.
% cueside_smps = randperm(3,3); otherside_smps = randperm(3,3);
% task.thistrial.cueside_smp = cueside_smps(1);  % Store which sample was the oddball on each trial.
% task.thistrial.otherside_smp = otherside_smps(1);
oddball_smps = randperm(4, 2); dist_smps = randperm(4, 4);
task.thistrial.cueside_smp = oddball_smps(1); % Store which sample was the oddball on each trial.
task.thistrial.otherside_smp = oddball_smps(2);

% Load the two oddball images.
stimulus.live.cueside_odd = genTexFromIm(imread(sprintf('%s/%s_%s_%s_%s_%s_smp%i.png', texDir, cueside_poolSize, cueside_layer, cueside_imFam, cueside_PC, ob_PCDist, oddball_smps(1))), stimulus.live.mask);
stimulus.live.otherside_odd = genTexFromIm(imread(sprintf('%s/%s_%s_%s_%s_%s_smp%i.png', texDir, otherside_poolSize, otherside_layer, otherside_imFam, otherside_PC, ob_PCDist, oddball_smps(2))), stimulus.live.mask);

% Load the 4 distractor images.
stimulus.live.cueside_dist1 = genTexFromIm(imread(sprintf('%s/%s_%s_%s_%s_%s_smp%i.png', texDir, cueside_poolSize, cueside_layer, cueside_imFam, cueside_PC, d_PCDist, dist_smps(1))), stimulus.live.mask);
stimulus.live.cueside_dist2 = genTexFromIm(imread(sprintf('%s/%s_%s_%s_%s_%s_smp%i.png', texDir, cueside_poolSize, cueside_layer, cueside_imFam, cueside_PC, d_PCDist, dist_smps(2))), stimulus.live.mask);
stimulus.live.otherside_dist1 = genTexFromIm(imread(sprintf('%s/%s_%s_%s_%s_%s_smp%i.png', texDir, otherside_poolSize, otherside_layer, otherside_imFam, otherside_PC, d_PCDist, dist_smps(3))), stimulus.live.mask);
stimulus.live.otherside_dist2 = genTexFromIm(imread(sprintf('%s/%s_%s_%s_%s_%s_smp%i.png', texDir, otherside_poolSize, otherside_layer, otherside_imFam, otherside_PC, d_PCDist, dist_smps(4))), stimulus.live.mask);

%% Select cueside position and target size.
task.thistrial.cueSide = randi(2); % 1 is left 2 is right
task.thistrial.cueside_oddpos = randi(3);
task.thistrial.otherside_oddpos = randi(3);

% Disp trial parameters each trial
disp(sprintf('Trial %d - %s, %s, PC: %s, PC_Dist: %s, PoolSize: %s, cueSide: %i, isCueFocal: %i', task.trialnum, cueside_imFam, cueside_layer,cueside_PC, ob_PCDist, cueside_poolSize, task.thistrial.cueSide, task.thistrial.isCueFocal));
if task.trialnum > 1
    disp(sprintf('--Target on Last Trial: %g, Response on last trial: %g, LastTrialCorrect?: %g', task.lasttrial.cueside_oddpos, task.lasttrial.response, task.lasttrial.correct));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)

global stimulus

stimulus.live.triggerWaiting = 0;

if any(task.thistrial.thisseg==[stimulus.seg.fix])
  stimulus.live.triggerWaiting = 1;
  stimulus.live.centered = 0;
  stimulus.live.triggerTime = 0;
  stimulus.live.lastTrigger = -1;
end

stimulus.live.eyeCount = 0;
stimulus.live.eyeDead = 0;

% Set the fixation cross color.
if task.thistrial.thisseg == stimulus.seg.feedback
  if task.thistrial.response == task.thistrial.cueside_oddpos
    task.thistrial.correct = 1;
    fixColor = stimulus.colors.green;
  else
    task.thistrial.correct = 0;
    fixColor = stimulus.colors.red;
  end
elseif task.thistrial.thisseg == stimulus.seg.stim
  fixColor = stimulus.colors.white;
else
  fixColor = stimulus.colors.black;
end

% Select image parameters: size, eccentricity, and location
imSz = stimulus.imSize; % Size in Degrees of Visual Angle to display image
ecc = stimulus.eccentricity; % Eccentricity to display image at

% Define possible stimulus locations.
locations = zeros(2, 3, 2);
locations(1,:,:) = [-ecc/sqrt(2) ecc/sqrt(2); -ecc 0; -ecc/sqrt(2) -ecc/sqrt(2)];
locations(2,:,:) = [ecc/sqrt(2) ecc/sqrt(2); ecc 0; ecc/sqrt(2) -ecc/sqrt(2)];

% Define possible cue positions
cueX = stimulus.cueEcc/sqrt(2);
cueXLocs = [-cueX, -cueX, -cueX-1; cueX, cueX, cueX+1];
cueYLocs = [-.25, .25, 0; -.25, .25, 0];
cueColor = stimulus.colors.blue;
cueSide = task.thistrial.cueSide;
otherSide = setdiff([1 2], cueSide);

% Select which location each stimulus will go in.
cueside_oddpos = task.thistrial.cueside_oddpos;
otherside_oddpos = task.thistrial.otherside_oddpos;
cueside_distpos = setdiff(1:3, cueside_oddpos);
otherside_distpos = setdiff(1:3, otherside_oddpos);

for i = 1:2
  mglClearScreen(0.5);
  if task.thistrial.thisseg == stimulus.seg.stim
    mglBltTexture(stimulus.live.cueside_odd, [squeeze(locations(cueSide, cueside_oddpos, :))' imSz imSz]);
    mglBltTexture(stimulus.live.otherside_odd, [squeeze(locations(otherSide, otherside_oddpos, :))' imSz imSz]);
    
    mglBltTexture(stimulus.live.cueside_dist1, [squeeze(locations(cueSide, cueside_distpos(1), :))' imSz imSz]);
    mglBltTexture(stimulus.live.cueside_dist2, [squeeze(locations(cueSide, cueside_distpos(2), :))' imSz imSz]);
    
    mglBltTexture(stimulus.live.otherside_dist1, [squeeze(locations(otherSide, otherside_distpos(1), :))', imSz, imSz]);
    mglBltTexture(stimulus.live.otherside_dist2, [squeeze(locations(otherSide, otherside_distpos(2), :))', imSz, imSz]);
  
  end
  upFix(stimulus, fixColor);  
  
  if task.thistrial.isCueFocal == 1 || any(task.thistrial.thisseg == [stimulus.seg.feedback stimulus.seg.response])
    %draw 1 triangle
    mglPolygon(cueXLocs(task.thistrial.cueSide,:), cueYLocs(task.thistrial.cueSide,:), stimulus.colors.black);
  else 
    %draw 2 triangles
    mglPolygon(cueXLocs(1,:), cueYLocs(1,:), stimulus.colors.black);
    mglPolygon(cueXLocs(2,:), cueYLocs(2,:), stimulus.colors.black);

  end
  mglFlush
end
stimulus.live.firstTime = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus

% jump to next trial if you are dead and 1 second has elapsed since eye
% movement
if task.thistrial.dead && mglGetSecs(task.thistrial.segStartSeconds)>1
  task = jumpSegment(task,inf);
end

% skip screen updates if you are already dead
if task.thistrial.dead
  if stimulus.live.eyeDead
    mglTextSet([],32,stimulus.colors.red);
    mglTextDraw('Eye Movement Detected',[0 0]);
  end
  return
end

% check eye pos
if ~stimulus.noeye
  [pos,~] = mglEyelinkGetCurrentEyePos;
  dist = hypot(pos(1),pos(2));
end

% Eye movement detection code
if ~stimulus.noeye && ~any(task.thistrial.thisseg==[stimulus.seg.fix]) 
  if ~any(isnan(pos))
    if dist > 1.5 && stimulus.live.eyeCount > 30
      disp('Eye movement detected!!!!');
      task.thistrial.dead = 1;
      stimulus.live.eyeDead=1;
      return
    elseif dist > 1.5
      stimulus.live.eyeCount = stimulus.live.eyeCount + 1;
    end
  end
end

% Trial trigger on eye fixation code  
if ~stimulus.noeye && stimulus.live.triggerWaiting
  now = mglGetSecs;
  % check eye position, if 
  if ~any(isnan(pos))
    wasCentered = stimulus.live.centered;
    stimulus.live.centered = dist<2.5;
    if wasCentered && stimulus.live.centered && stimulus.live.lastTrigger>0
      stimulus.live.triggerTime = stimulus.live.triggerTime + now-stimulus.live.lastTrigger;
    end
    stimulus.live.lastTrigger = now;
  end
  if stimulus.live.triggerTime > 0.5 % not in ms dummy, wait 1.5 seconds (reasonable slow time)
    disp('Starting trial--eye centered.');
    task = jumpSegment(task);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus

if task.thistrial.dead, return; end

validResponse = any(task.thistrial.whichButton == stimulus.responseKeys);
if validResponse
  if stimulus.live.gotResponse==0
    task.thistrial.detected = 1;
    task.thistrial.response = task.thistrial.whichButton-10;
    stimulus.live.fix = 0;
  else
    disp(sprintf('Subject responded multiple times: %i',stimulus.live.gotResponse));
  end
  stimulus.live.gotResponse=stimulus.live.gotResponse+1;
  task = jumpSegment(task);
else
  disp(sprintf('Invalid response key. Subject pressed %d', task.thistrial.whichButton));
  task.thistrial.response = -1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              HELPER FUNCTIONS                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
% Draws a circle at center of the screen of color fixColor
function upFix(stimulus, fixColor)
if ieNotDefined('fixColor')
  fixColor = stimulus.live.fixColor;
end
%mglGluAnnulus(0,0,0,.1,fixColor);
mglFixationCross(1,1,fixColor);


%%% 
% Draws a circular cue at location x,y
function drawCue(x,y, stimulus)
mglGluAnnulus(x,y, 0.75, 0.8, stimulus.live.cueColor, 64);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Turns image into a texture
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localInitStimulus()

global stimulus
sz = 1.5;

% use total degs / num to compute size
grating = 251/2*mglMakeGrating(sz,sz,2,0) + 255/2;
gauss = mglMakeGaussian(sz,sz,sz/6,sz/6);
alphamask = repmat(grating,1,1,4);
alphamask(:,:,4) = gauss*255;

% we'll adjust the gamma table to control contrast
stimulus.live.grating  = mglCreateTexture(alphamask); % high contrast        

%%%%%%%%%%%%%%%%%%%%%%%
%    analyzeData %
%%%%%%%%%%%%%%%%%%%%%%%
function data = analyzeData()
%%

% get the files list
files = dir(fullfile(sprintf('~/data/texAtt_PCA/%s/19*stim*.mat',mglGetSID)));

count = 1; 
data = struct('response', [], 'reaction_time', [], 'nTrials', 0, 'nValTrials', 0, 'accByRuns', []);

for fi = 1:length(files)
  load(fullfile(sprintf('~/data/texAtt_PCA/%s/%s',mglGetSID,files(fi).name)));
  
  e = getTaskParameters(myscreen,task);
  if e.numTrials>1
    
    f = fields(e.parameter);
    for i = 1:length(f)
        if ~isfield(data, f{i})
            data.(f{i}) = [];
        end
        data.(f{i}) = [data.(f{i}) e.parameter.(f{i})];
    end
    f = fields(e.randVars);
    for i = 1:length(f)
        if ~isfield(data, f{i})
            data.(f{i}) = [];
        end
        data.(f{i}) = [data.(f{i}) e.randVars.(f{i})];
    end
    
    data.response = [data.response e.response-10];

    data.reaction_time = [data.reaction_time e.reactionTime];
    data.nTrials = data.nTrials + e.nTrials;

    % Calculate number of valid trials by excluding eye movements
    data.nValTrials = data.nValTrials + sum(~isnan(e.response));
    
    data.accByRuns = [data.accByRuns nanmean(e.randVars.correct)];
    
  end
  count = count + 1;
end

%%
% all poolSizes are 1x1, also recorded incorrectly for stimfiles 1-7
%all_pools = unique(data.cueside_poolSize);
all_cueTypes = unique(data.isCueFocal);

all_PCDists = unique(data.ob_PCDists);
accs = nan(length(all_PCDists), length(all_cueTypes));
SEs = nan(length(all_PCDists), length(all_cueTypes));
for i = 1:length(all_PCDists)
    for j = 1:length(all_cueTypes)
        ct = data.correct(data.ob_PCDist==all_PCDists(i) & data.isCueFocal==all_cueTypes(j));
        accs(i,j) = nanmean(ct);
        SEs(i,j) = 1.96*nanstd(ct) / length(ct);
    end
end

%%
%poolsizes = [1, 1.25, 1.5, 1.75, 2, 3, 4];
PCDists = [-10, -5, -2, 2, 5, 10];
%x = 6./poolsizes;
x = PCDists;
figure;
h1 = myerrorbar(x, accs(:,1), 'yError', SEs(:,1), 'Color', 'g', 'Symbol=o-');
h2 = myerrorbar(x, accs(:,2), 'yError', SEs(:,2), 'Color', 'b', 'Symbol=o-'); hold on;
xlim([0 length(all_PCDists)+1]);

hline(1/3, ':');
legend([h1,h2], {'Distributed', 'Focal'});
set(gca, 'XTick', sort(x));
set(gca, 'XTickLabel', round(sort(x),2));
xlabel('Oddball Pooling Size (degrees)');
ylabel('Accuracy (% correct)');
title(sprintf('Oddity Task: nTrials = %i', data.nTrials));

%%
all_PCs = unique(data.PCs);
accs = nan(length(all_PCDists), length(all_cueTypes), length(all_PCs));
SEs = nan(length(all_PCDists), length(all_cueTypes), length(all_PCs));
for i = 1:length(all_PCDists)
  for j = 1:length(all_cueTypes)
    for k = 1:length(all_PCs)
      ct = data.correct(data.ob_PCDist==all_pools(i) & data.isCueFocal==all_cueTypes(j) & data.cueside_imFam==all_PCs(k));
      accs(i,j,k) = nanmean(ct);
      SEs(i,j,k) = 1.96*nanstd(ct) / length(ct);
    end
  end
end

%%
figure;
for i = 1:5
  for j = 1:4
    ind = 4*(i-1)+j;
    
    if ind < length(all_PCs)
        subplot(5,4, ind);
        h1 = myerrorbar(x, accs(:,1,ind), 'yError', SEs(:,1,ind), 'Color', 'g', 'Symbol=o');
        h2 = myerrorbar(x, accs(:,2,ind), 'yError', SEs(:,2,ind), 'Color', 'b', 'Symbol=o'); hold on;
        xlim([0 length(all_pools)+1]);
        ylim([-.2 1.2]);

        hline(1/3, ':');
        if ind == length(all_PCs)-1
          legend([h1,h2], {'Distributed', 'Focal'});
        end
        set(gca, 'XTick', sort(x));
        set(gca, 'XTickLabel', round(sort(x),2));
        xlabel('Oddball Pooling Size (degrees)');
        ylabel('Accuracy (% correct)');
        title(sprintf('%s: nTrials = %i', stimulus.imNames{all_PCs(ind)}, data.nTrials));
    end
  end
end

%%
keyboard

