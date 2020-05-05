function [ myscreen ] = texOdd( varargin )
%
% TEXTURE ODDBALL TASK
%  Oddball task using textures 
%
%  Usage: texOdd(varargin)
%  Authors: Akshay Jagadeesh
%  Date: 02/27/2018
%

global stimulus

stimulus = struct;

%% Initialize Variables

% add arguments later
plots = 0;
noeye = 0;
getArgs(varargin,{'periph=0','plots=0','noeye=1', 'analyze=0', 'training=0'}, 'verbose=1');
stimulus.plots = plots;
stimulus.noeye = noeye;
stimulus.training = training;
stimulus.analyze = analyze;
stimulus.periph = periph;
clear noeye plots training analyze

if stimulus.analyze
  analyzeData();
  return
end

%% Stimulus parameters 
%% Open Old Stimfile
stimulus.counter = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/texOdd/%s',mglGetSID))
  % Directory exists, check for a stimfile
  files = dir(sprintf('~/data/texOdd/%s/1*mat',mglGetSID));

  if length(files) >= 1
    fname = files(end).name;
    
    s = load(sprintf('~/data/texOdd/%s/%s',mglGetSID,fname));
    stimulus.counter = s.stimulus.counter + 1;
    clear s;
    disp(sprintf('(texOdd) Data file: %s loaded.',fname));
  end
end
disp(sprintf('(texOdd) This is run #%i',stimulus.counter));

%% Setup Screen
myscreen = initScreen('VPixx2');

% set background to grey
myscreen.background = 0.5;

%% Initialize Stimulus
myscreen = initStimulus('stimulus',myscreen);
localInitStimulus();
  
% Set response keys
stimulus.responseKeys = [11 14 12];
%% [aj]: change to add response keys for whichever 3 keys we want.

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
task{1}.segmin = [inf, 0.5, 2.0, 0.2];
task{1}.segmax = [inf, 0.5, 2.0, 0.2];
stimulus.seg = {};
stimulus.seg.fix = 1;
stimulus.seg.stim = 2;
stimulus.seg.response = 3;
stimulus.seg.feedback = 4;
% Set fixation length to constant if not eye tracking.
if stimulus.noeye==1
  task{1}.segmin(1) = 0.2;
  task{1}.segmax(1) = 0.2;
end

% Set when to synchtovol and getResponse
task{1}.synchToVol = zeros(size(task{1}.segmin));
task{1}.getResponse = zeros(size(task{1}.segmin));
task{1}.getResponse(stimulus.seg.response)=1;

%%% Stimulus Variables
stimulus.imNames = {'face', 'jetplane', 'elephant', 'sand', 'lawn', 'tulips', 'dirt', 'fireworks', 'bananas'};
%stimulus.imNames = {'beans', 'blossoms', 'bubbly', 'clouds', 'crystals', 'dahlias',...
%          'fronds', 'fur', 'glass', 'leaves', 'leopard', 'noodles', 'paisley', 'plant',...
%          'rocks', 'scales', 'spikes', 'tiles', 'waves', 'worms'};
stimulus.layerNames = {'pool1', 'pool2', 'pool4'};
stimulus.stimDir = '~/proj/texture_stimuli/color/textures';
stimulus.origImDir = '~/proj/texture_stimuli/color/originals';
stimulus.imSize = 4;
if stimulus.periph
  stimulus.eccentricity = 10;
else
  stimulus.eccentricity = 3;
end
stimulus.poolSizes = {'1x1', '2x2', '3x3', '4x4'}; % Add back 2x2, 3x3, and 4x4 later.
stimulus.cueEcc = 4;
stimulus.live.mask = imread('~/proj/TextureSynthesis/stimuli/Flattop8.tif');

% Trial parameters
task{1}.parameter.layer = 1:length(stimulus.layerNames);
task{1}.parameter.poolSize = 1:length(stimulus.poolSizes);

% Make numTrials some multiple of number of TrialTypes 
task{1}.numTrials = 240;
task{1}.random = 1;

%%% Task Variables
% Keep track of texturefamily, oddball layer/poolsize, standard layer/poolsize, and target position.
task{1}.randVars.uniform.imgFam = 1:length(stimulus.imNames);
task{1}.randVars.uniform.oddball_layer = [0]; % oddball is always original
task{1}.randVars.uniform.distractor_layer = 1:length(stimulus.layerNames);
task{1}.randVars.uniform.distractor_poolsize = 1:length(stimulus.poolSizes);
task{1}.randVars.uniform.targetPosition = 1:3; % Track target position

if stimulus.training
  task{1}.randVars.uniform.oddball_layer = 0;
  task{1}.randVars.uniform.distractor_layer = 1:2;
  task{1}.numTrials = 50;
end
% Task variables to keep track of the status of each trial
task{1}.randVars.calculated.detected = 0; % did they see the grating
task{1}.randVars.calculated.dead = 0;
task{1}.randVars.calculated.correct = NaN;

%% Preload images
stimulus.nSamples = 2; % preload 3 samples of each kind.
%if ~exist(presavedStimLoc) % on the first time, load each image, convert to mgl texture, and save it to a struct.
stims = struct();
disppercent(-inf, 'Preloading images');
for i = 1:length(stimulus.imNames)
  imName = stimulus.imNames{i};
  orig = imread(sprintf('%s/%s.jpg', stimulus.origImDir, imName));
  stims.(imName) = genTexFromIm(orig, stimulus.live.mask);
  for j = 1:length(stimulus.layerNames)
    for k = 1:length(stimulus.poolSizes)
      for l = 1:stimulus.nSamples
        ps = stimulus.poolSizes{k}; ln = stimulus.layerNames{j};
        if strcmp(ln, 'PS'), ps = '1x1'; end
        filename = sprintf('%s/%s_%s_%s_smp%i.%s', stimulus.stimDir, ps, ln, imName, l, 'png');
        if ~isfile(filename)
            filename = sprintf('%s/%s_%s_%s_smp%i.%s', stimulus.stimDir, ps, ln, imName, l, 'jpg');
        end
        smp = imread(filename);
        stims.(sprintf('%s_%s_%s_smp%i', imName, ps, ln, l)) = genTexFromIm(smp, stimulus.live.mask);
      end
    end
    disppercent( ( (i-1) + (j/length(stimulus.layerNames))) / length(stimulus.imNames));
  end
  disppercent(i / length(stimulus.imNames));
end
disppercent(inf);
%  save(presavedStimLoc, '-struct', 'stims');
%else % after first time, just load the saved struct containing the mgl textures.
%  disp('Presaved stimulus structure found. Loading now...');
%  stims = load(presavedStimLoc);
%end
stimulus.live.stims = stims;

%% Save
clear smp smpName orig stims

%%%%%%%%%%%%%%%% MGL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialize task and begin update loop
[task{1}, myscreen] = initTask(task{1},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);

% Run the eye calibration
myscreen = eyeCalibDisp(myscreen);

% let the user know
disp(sprintf('(texOdd) Starting run number: %i.',stimulus.counter));

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
task.thistrial.response = NaN;

stimulus.live.gotResponse = 0;
stimulus.curTrial(task.thistrial.thisphase) = stimulus.curTrial(task.thistrial.thisphase) + 1;

% directories
texDir = stimulus.stimDir;

%% Load all 3 images for this trial
% Load all 3 images for this trial
imName = stimulus.imNames{task.thistrial.imgFam};

std_layer = stimulus.layerNames{task.thistrial.distractor_layer};
std_poolsize = stimulus.poolSizes{task.thistrial.distractor_poolsize};
if strcmp(std_layer, 'PS')
    std_poolsize = '1x1';
    task.thistrial.distractor_poolsize = 1;
end
stimulus.live.d1 = stimulus.live.stims.(sprintf('%s_%s_%s_smp1', imName, std_poolsize, std_layer));
stimulus.live.d2 = stimulus.live.stims.(sprintf('%s_%s_%s_smp2', imName, std_poolsize, std_layer));

stimulus.live.targetImg = stimulus.live.stims.(imName);
oddball_text = 'original';

stimulus.live.eyeCount = 0;

% Disp trial parameters each trial
fprintf('Trial %d - %s: Oddball = %s, Standard = %s %s, TargetLoc: %i - ', task.trialnum, imName, oddball_text, std_poolsize, std_layer, task.thistrial.targetPosition);

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

stimulus.live.eyeDead = 0;

% Select image parameters: size, eccentricity, and location
imSz = stimulus.imSize; % Size in Degrees of Visual Angle to display image
ecc = stimulus.eccentricity; % Eccentricity to display image at
stimulus.locations = [0, ecc; ecc*cosd(30), -ecc*sind(30); -ecc*cosd(30), -ecc*sind(30)];

targetLoc = stimulus.locations(task.thistrial.targetPosition,:);
distLocs = stimulus.locations(setdiff(1:3, task.thistrial.targetPosition), :);

for i = 1:2
  mglClearScreen(0.5);
  if task.thistrial.thisseg == stimulus.seg.stim
    mglBltTexture(stimulus.live.targetImg, [targetLoc(1), targetLoc(2), imSz, imSz]);
    mglBltTexture(stimulus.live.d1, [distLocs(1,:) imSz imSz]);
    mglBltTexture(stimulus.live.d2, [distLocs(2,:) imSz, imSz]);
    upFix(stimulus, stimulus.colors.black);
  elseif task.thistrial.thisseg == stimulus.seg.feedback
    if task.thistrial.response == task.thistrial.targetPosition
      task.thistrial.correct = 1;
      upFix(stimulus, stimulus.colors.green);
      if i==1, fprintf('Correct! \n'); end
    else
      task.thistrial.correct = 0;
      upFix(stimulus, stimulus.colors.red);
      if i==1, fprintf('Incorrect :( \n'); end
    end
  else
    upFix(stimulus, stimulus.colors.black);
  end
  mglFlush
end
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
  if task.thistrial.dead && stimulus.live.eyeDead
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
    if dist > 1.5 && stimulus.live.eyeCount > 20
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

% only allow a response if it is one of the pre-specified responses (in this case, y g or h)
validResponse = any(task.thistrial.whichButton == stimulus.responseKeys);
if validResponse
  if stimulus.live.gotResponse==0
    task.thistrial.detected = 1;

    % Set thistrial response to index of response key (i.e. 1,2, or 3) rather than the button id.
    task.thistrial.response = find(stimulus.responseKeys == task.thistrial.whichButton);
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
files = dir(fullfile(sprintf('~/data/texOdd/%s/19*stim*.mat',mglGetSID)));

count = 1; 
data = struct('response', [], 'reaction_time', [], 'nTrials', 0, 'nValTrials', 0, 'accByRuns', []);
for fi = 1:length(files)
  load(fullfile(sprintf('~/data/texOdd/%s/%s',mglGetSID,files(fi).name)));
  
  e = getTaskParameters(myscreen,task);
  if e.nTrials>1
    
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
disp(sprintf('SUBJECT %s: Found %i runs with a total of %i trials', mglGetSID, length(data.accByRuns), data.nTrials));

%% Get average accuracies across each condition.
standard_pools = unique(data.distractor_poolsize);
distractor_layers = unique(data.distractor_layer);
odd_pools = unique(data.oddball_poolsize);
odd_layers = unique(data.oddball_layer);

accs = nan(length(standard_pools)*length(distractor_layers), length(odd_pools)*length(odd_layers));
SEs = nan( length(standard_pools)*length(distractor_layers), length(odd_pools)*length(odd_layers));
Ns = nan( length(standard_pools)*length(distractor_layers), length(odd_pools)*length(odd_layers));
x = 0; 
std_labels = {}; odd_labels = {};
for i = 1:length(standard_pools)
  for j = 1:length(distractor_layers)
   x = x+1; y = 0;
   std_labels{x} = sprintf('%s %s', stimulus.poolSizes{standard_pools(i)}, stimulus.layerNames{distractor_layers(j)});
   for k = 1:length(odd_pools)
      for l = 1:length(odd_layers)
        y = y+1;
        ct = data.correct(data.distractor_poolsize==standard_pools(i) & data.distractor_layer==distractor_layers(j) & data.oddball_poolsize==odd_pools(k) & data.oddball_layer==odd_layers(l));
        accs(x,y) = nanmean(ct);
        SEs(x,y) = 1.96*nanstd(ct) / sqrt(length(ct));
        Ns(x,y) = length(ct);
        if odd_layers(l)==0
            odd_labels{y} = 'Original';
        elseif odd_layers(l)==1
            odd_labels{y} = 'P-S';
        else
            odd_labels{y} = sprintf('%s %s', stimulus.poolSizes{odd_pools(k)}, stimulus.layerNames{odd_layers(l)});
        end
      end
    end
  end
end

%% Plot heatmap of accuracies in each condition (distractor conditions x target conditions)
figure;
%subplot(1,2,1);
imagesc(accs);
for i=1:size(accs,1)
    for j =1:size(accs,2)
        text(j-0.5,i, sprintf('%.2f', accs(i,j)));
    end
end
xlabel('Oddball Texture');
ylabel('Standard (non-oddball) Texture');
set(gca, 'XTick', 1:length(odd_pools)*length(odd_layers));
set(gca, 'YTick', 1:length(standard_pools)*length(distractor_layers));
set(gca, 'XTickLabel', odd_labels);
set(gca, 'XTickLabelRotation', 45);
set(gca, 'YTickLabel', std_labels);
title('Proportion correct')
colormap('Hot');
colorbar;

% Plot heatmap of Ns
%subplot(1,2,2);
figure;
imagesc(Ns);
for i=1:size(Ns,1)
    for j =1:size(Ns,2)
        text(j,i, sprintf('%i', Ns(i,j)));
    end
end
title('Number of trials');
xlabel('Oddball Texture');
ylabel('Standard (non-oddball) Texture');
set(gca, 'XTick', 1:length(odd_pools)*length(odd_layers));
set(gca, 'YTick', 1:length(standard_pools)*length(distractor_layers));
set(gca, 'XTickLabel', odd_labels);
set(gca, 'XTickLabelRotation', 45);
set(gca, 'YTickLabel', std_labels);
colormap('Hot');
colorbar;

%% When oddball is ORIGINAL`
accs = nan(length(standard_pools), length(distractor_layers));
SEs = nan( length(standard_pools), length(distractor_layers));
Ns = nan( length(standard_pools), length(distractor_layers));
x = 0;
for i = 1:length(standard_pools)
  for j = 1:length(distractor_layers)
    if strcmp(stimulus.layerNames{distractor_layers(j)}, 'PS')
        ct = data.correct(data.distractor_layer==distractor_layers(j) & data.oddball_layer==0);
        if i~=1, ct = NaN; end
    else
        ct = data.correct(data.distractor_poolsize==standard_pools(i) & data.distractor_layer==distractor_layers(j) & data.oddball_layer==0);
    end
    accs(i,j) = nanmean(ct);
    SEs(i,j) = 1.96*nanstd(ct) / sqrt(length(ct));
    Ns(i,j) = length(ct);
  end
end

figure;
colors = brewermap(4, 'Set1');
handles = [];
for i = 1:size(accs,2)
  h = myerrorbar(1:size(accs,1), accs(i,:), 'yError', SEs(i,:), 'Symbol', 'o', 'Color', colors(i,:)); hold on;
  handles = [handles h];
end
xlim([0, size(accs,1)+1]);
hline(1/3, ':k');
legend(handles, stimulus.poolSizes);
title('Comparison to original');
set(gca, 'XTick', 1:size(accs,1));
set(gca, 'XTickLabel', stimulus.layerNames);
ylabel('Accuracy (proportion correct)')
xlabel('Model layer to which distractors were feature-matched')

%%
keyboard
figure; bar(accs');
set(gca, 'XTick', 1:length(odd_pools)*length(odd_layers));
set(gca, 'YTick', 1:length(standard_pools)*length(distractor_layers));
set(gca, 'XTickLabel', ['Original' stimulus.layerNames]);
legend(stimulus.layerNames);

%%
figure;
subplot(2,1,1);
myerrorbar(1:length(all_layers), mean(accs, 1), 'yError', mean(SEs,1), 'Symbol', 'o');
hold on; 
xlim([0 length(all_layers)+1]); ylim([0 1]);
hline(1/3, ':k');
xlabel('Layer'); ylabel('Accuracy');
set(gca, 'XTick', 1:length(all_layers));
set(gca, 'XTickLabel', stimulus.layerNames);

subplot(2,1,2);
myerrorbar(6./(1:length(all_pools)), mean(accs,2), 'yError', mean(SEs,2), 'Symbol', 'o');
hold on; 
xlim([0, 7]); ylim([0 1]);
hline(1/3, ':k');
xlabel('Pooling Region Size (degs)'); ylabel('Accuracy');
set(gca, 'XTick', sort(6./(1:length(all_pools)), 'ascend'));
%set(gca, 'XTickLabel', stimulus.poolSizes);


%%
figure;
handles = []; cols = brewermap(length(all_pools)+1, 'Blues');
for i = 1:length(all_pools)
  h = myerrorbar(1:length(all_layers), accs(i,:), 'yError', SEs(i,:), 'Symbol', 'o', 'Color', cols(i+1,:)); 
  hold on;
  handles = [handles h];
end
xlim([0 length(all_layers)+1]); ylim([0 1]);
hline(1/3, ':k');
xlabel('Layer'); ylabel('Accuracy');
legend(handles, stimulus.poolSizes);
set(gca, 'XTick', 1:length(all_layers));
set(gca, 'XTickLabel', stimulus.layerNames);
box off;
