function [ myscreen ] = texCrowd( varargin )
%
% TEXTURE ODDBALL TASK
%  Oddball task using textures 
%
%  Usage: texCrowd(varargin)
%  Authors: Akshay Jagadeesh
%  Date: 02/27/2018
%

global stimulus

stimulus = struct;

%% Initialize Variables

% add arguments later
plots = 0;
noeye = 0;
getArgs(varargin,{'periph=1','plots=0','noeye=1', 'analyze=0', 'training=0', 'crowders_present=1'}, 'verbose=1');
stimulus.plots = plots;
stimulus.noeye = noeye;
stimulus.training = training;
stimulus.analyze = analyze;
stimulus.periph = periph;
stimulus.crowders_present = crowders_present;
clear noeye plots training analyze

if stimulus.analyze
  analyzeData();
  return
end

%% Stimulus parameters 
%% Open Old Stimfile
stimulus.counter = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/texCrowd/%s',mglGetSID))
  % Directory exists, check for a stimfile
  files = dir(sprintf('~/data/texCrowd/%s/1*mat',mglGetSID));

  if length(files) >= 1
    fname = files(end).name;
    
    s = load(sprintf('~/data/texCrowd/%s/%s',mglGetSID,fname));
    stimulus.counter = s.stimulus.counter + 1;
    clear s;
    disp(sprintf('(texCrowd) Data file: %s loaded.',fname));
  end
end
disp(sprintf('(texCrowd) This is run #%i',stimulus.counter));

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
%stimulus.imNames = {'face', 'tulips', 'fur', 'bricks', 'blotch', 'bark'};
stimulus.imNames = {'face', 'jetplane', 'elephant', 'car', ...
                    'crowd', 'rocks', 'tulips', 'bananas',...
                    'sand', 'dirt', 'lawn', 'succulent'}; 
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
stimulus.poolSizes = {'1x1', '2x2', '3x3', '4x4'}; 
stimulus.live.mask = imread('~/proj/TextureSynthesis/stimuli/Flattop8.tif');

% Trial parameters
task{1}.parameter.crowders_present = [stimulus.crowders_present];
task{1}.parameter.distractor_layer = 1:length(stimulus.layerNames);
task{1}.parameter.distractor_poolsize = 1:length(stimulus.poolSizes);

% Make numTrials some multiple of number of TrialTypes 
task{1}.numTrials = 240;
task{1}.random = 1;

%%% Task Variables
% Keep track of texturefamily, oddball layer/poolsize, standard layer/poolsize, and target position.
task{1}.randVars.uniform.imgFam = 1:length(stimulus.imNames);
task{1}.randVars.uniform.targetPosition = 1:3; % Track target position

% Crowder parameters
if stimulus.crowders_present == 1
    task{1}.randVars.uniform.crowder_imgfam = 1:length(stimulus.imNames);
    task{1}.randVars.uniform.crowder_layer = 1:length(stimulus.layerNames);
else % if stimulus.crowders_present==2, assign imgfam and layer to same as target at startTrial.
    task{1}.randVars.calculated.crowder_imgfam = NaN;
    task{1}.randVars.calculated.crowder_layer = NaN;
end
task{1}.randVars.uniform.crowder_poolsize = 1:length(stimulus.poolSizes);

if stimulus.training
  task{1}.randVars.uniform.distractor_layer = 1:2;
  task{1}.numTrials = 50;
end

% Task variables to keep track of the status of each trial
task{1}.randVars.calculated.detected = 0; % did they see the grating
task{1}.randVars.calculated.dead = 0;
task{1}.randVars.calculated.correct = NaN;

%% Preload images
stimulus.nSamples = 3; % preload 3 samples of each kind.
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

stimulus.live.stims = stims;

%% Save
clear smp smpName orig stims

%%%%%%%%%%%%%%%% MGL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialize task and begin update loop
[task{1}, myscreen] = initTask(task{1},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);

% Run the eye calibration
myscreen = eyeCalibDisp(myscreen);

% let the user know
disp(sprintf('(texCrowd) Starting run number: %i.',stimulus.counter));

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

%% Load both crowder images for this trial
if task.thistrial.crowders_present == 2
    task.thistrial.crowder_imgfam = task.thistrial.imgFam;
    task.thistrial.crowder_layer = task.thistrial.distractor_layer;
end
crowder_imgfam = stimulus.imNames{task.thistrial.crowder_imgfam};
crowder_layer = stimulus.layerNames{task.thistrial.crowder_layer};
crowder_poolsize = stimulus.poolSizes{task.thistrial.crowder_poolsize};
stimulus.live.c1 = stimulus.live.stims.(sprintf('%s_%s_%s_smp1', crowder_imgfam, crowder_poolsize, crowder_layer));
stimulus.live.c2 = stimulus.live.stims.(sprintf('%s_%s_%s_smp2', crowder_imgfam, crowder_poolsize, crowder_layer));

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

    % Draw crowders
    if task.thistrial.crowders_present
      mglBltTexture(stimulus.live.c1, [targetLoc(1)-4, targetLoc(2), imSz, imSz]);
      mglBltTexture(stimulus.live.c2, [targetLoc(1)+4, targetLoc(2), imSz, imSz]);
      mglBltTexture(stimulus.live.c1, [distLocs(1,1)-4, distLocs(1,2), imSz, imSz]);
      mglBltTexture(stimulus.live.c2, [distLocs(1,1)+4, distLocs(1,2), imSz, imSz]);
      mglBltTexture(stimulus.live.c1, [distLocs(2,1)-4, distLocs(2,2), imSz, imSz]);
      mglBltTexture(stimulus.live.c2, [distLocs(2,1)+4, distLocs(2,2), imSz, imSz]);
    end
    
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
%    fixResponse %
%  fix the response mapping which I messed up like an idiot
%%%%%%%%%%%%%%%%%%%%%%%
function resp = fixResponse(response)

resp(response==11) = 1;
resp(response==14) = 2;
resp(response==12) = 3;


%%%%%%%%%%%%%%%%%%%%%%%
%    analyzeData %
%%%%%%%%%%%%%%%%%%%%%%%
function data = analyzeData()
%%

% get the files list
files = dir(fullfile(sprintf('~/data/texCrowd/%s/20*stim*.mat',mglGetSID)));

count = 1; 
data = struct('response', [], 'reaction_time', [], 'periphery',[], 'eccentricity', [], ...
              'nTrials', 0, 'nValTrials', 0, 'accByRuns', []);
for fi = 1:length(files)
  load(fullfile(sprintf('~/data/texCrowd/%s/%s',mglGetSID,files(fi).name)));
  
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
    
    data.response = [data.response fixResponse(e.response)];
    data.reaction_time = [data.reaction_time e.reactionTime];
    data.nTrials = data.nTrials + e.nTrials;
    
    data.periphery = [data.periphery ones(1,e.nTrials)*stimulus.periph];
    data.eccentricity = [data.eccentricity ones(1,e.nTrials)*stimulus.eccentricity];

    % Calculate number of valid trials by excluding eye movement trials and no-response trials.
    data.nValTrials = data.nValTrials + sum(~isnan(e.response));
    
    data.accByRuns = [data.accByRuns nanmean(e.randVars.correct)];
    
  end
  count = count + 1;
end
data.imSize = stimulus.imSize;

disp(sprintf('SUBJECT %s: Found %i runs with a total of %i trials', mglGetSID, length(data.accByRuns), data.nTrials));

keyboard

%% First, plot results across both eccentricities.
distractor_pools = unique(data.distractor_poolsize);
distractor_layers = unique(data.distractor_layer);
accs = nan(length(distractor_pools), length(distractor_layers));
SEs = nan( length(distractor_pools), length(distractor_layers));
Ns = nan( length(distractor_pools), length(distractor_layers));
x = 0;
for i = 1:length(distractor_pools)
  for j = 1:length(distractor_layers)
    ct = data.correct(data.distractor_poolsize==distractor_pools(i) & data.distractor_layer==distractor_layers(j));
    accs(i,j) = nanmean(ct);
    SEs(i,j) = 1.96*nanstd(ct) / sqrt(length(ct));
    Ns(i,j) = length(ct);
  end
end

figure;
colors = brewermap(4, 'Set1');
handles = [];
for i = 1:size(accs,1)
  h = myerrorbar(1:size(accs,2), accs(i,:), 'yError', SEs(i,:), 'Symbol', 'o', 'Color', colors(i,:)); hold on;
  handles = [handles h];
end
xlim([0, length(distractor_layers)+1]);
hline(1/3, ':k');
legend(handles, stimulus.poolSizes);
title('Overall results');
set(gca, 'XTick', 1:length(distractor_layers));
set(gca, 'XTickLabel', stimulus.layerNames);
ylabel('Accuracy (proportion correct)')
xlabel('Model layer to which distractors were feature-matched')

%% Next, split out results by eccentricity.
eccentricities = unique(data.periphery);
accs = nan(length(distractor_pools), length(distractor_layers), length(eccentricities));
SEs  = nan(length(distractor_pools), length(distractor_layers), length(eccentricities));
Ns   = nan(length(distractor_pools), length(distractor_layers), length(eccentricities));
for i = 1:length(distractor_pools)
  for j = 1:length(distractor_layers)
    for k = 1:length(eccentricities)
      ct = data.correct(data.distractor_poolsize==distractor_pools(i) & data.distractor_layer==distractor_layers(j) & data.periphery==eccentricities(k));
      accs(i,j,k) = nanmean(ct);
      SEs(i,j,k) = 1.96*nanstd(ct) / sqrt(length(ct));
      Ns(i,j,k) = length(ct);
    end
  end
end

figure;
colors = brewermap(4, 'Set2');
subplot(1,3,1);
handles = [];
for i = 1:size(accs,1)
  h = myerrorbar(1:length(distractor_layers), accs(i,:,1), 'yError', SEs(i,:,1), 'Symbol', 'o', 'Color', colors(i,:)); hold on;
  handles = [handles h];
end
xlim([0, length(distractor_layers)+1]); ylim([0.25, 1]);
hline(1/3, ':k');
legend(handles, stimulus.poolSizes);
title('Foveal (ecc=3deg)');
set(gca, 'XTick', 1:length(distractor_layers));
set(gca, 'XTickLabel', stimulus.layerNames);
ylabel('Accuracy (proportion correct)')
xlabel('Model layer to which distractors were feature-matched')

subplot(1,3,2);
handles = [];
for i = 1:size(accs,1)
  h = myerrorbar(1:length(distractor_layers), accs(i,:,2), 'yError', SEs(i,:,2), 'Symbol', 's', 'Color', colors(i,:)); hold on;
  handles = [handles h];
end
xlim([0, length(distractor_layers)+1]); ylim([0.25, 1]);
hline(1/3, ':k');
legend(handles, stimulus.poolSizes);
title('Peripheral (ecc=10deg)');
set(gca, 'XTick', 1:length(distractor_layers));
set(gca, 'XTickLabel', stimulus.layerNames);
ylabel('Accuracy (proportion correct)')
xlabel('Model layer to which distractors were feature-matched')

subplot(1,3,3);
delta_acc = accs(:,:,1) - accs(:,:,2);
bar(delta_acc'); colormap(colors);
title('Difference in performance between foveal and peripheral');
ylabel('Delta Performance (foveal - peripheral)');
legend(stimulus.poolSizes);
set(gca, 'XTick', 1:length(distractor_layers));
set(gca, 'XTickLabel', stimulus.layerNames);
box off;

%%

%% Finally, split out results by individual images.
imageClasses = unique(data.imgFam);
accs = nan(length(distractor_pools), length(distractor_layers), length(eccentricities), length(imageClasses));
SEs  = nan(length(distractor_pools), length(distractor_layers), length(eccentricities), length(imageClasses));
Ns   = nan(length(distractor_pools), length(distractor_layers), length(eccentricities), length(imageClasses));
for i = 1:length(distractor_pools)
  for j = 1:length(distractor_layers)
    for k = 1:length(eccentricities)
      for l = 1:length(imageClasses)
        ct = data.correct(data.distractor_poolsize==distractor_pools(i) & data.distractor_layer==distractor_layers(j) ...
                          & data.periphery==eccentricities(k) & data.imgFam==imageClasses(l));
        accs(i,j,k,l) = nanmean(ct);
        SEs(i,j,k,l) = 1.96*nanstd(ct) / sqrt(length(ct));
        Ns(i,j,k,l) = length(ct);
      end
    end
  end
end

figure;
for l = 1:length(imageClasses)
  subplot(1,length(imageClasses),l);
  delta_acc = accs(:,:,1,l) - accs(:,:,2,l);
  bar(delta_acc'); colormap(colors);
  title(sprintf('Image Class: %s', stimulus.imNames{l}));
  ylabel('Delta Performance (foveal - peripheral)');
  legend(stimulus.poolSizes, 'Location', 'northwest');
  set(gca, 'XTick', 1:length(distractor_layers));
  set(gca, 'XTickLabel', stimulus.layerNames);
  box off;
end

%% Plot mean accuracies across conditions, split by individual images.
periph_acc = squeeze(accs(:,:,2,:));
periph_Ns = squeeze(Ns(:,:,2,:));
periph_SEs = squeeze(SEs(:,:,2,:));

figure;
for l = 1:length(imageClasses)
  subplot(1,length(imageClasses),l);
  
  handles = [];
  for i = 1:size(accs,1)
    h = myerrorbar(1:length(distractor_layers), periph_acc(i,:,l), 'yError', periph_SEs(i,:,2), 'Symbol', 'o', 'Color', colors(i,:)); hold on;
    handles = [handles h];
  end
  xlim([0, length(distractor_layers)+1]); ylim([0.25, 1]);
  hline(1/3, ':k');
  legend(handles, stimulus.poolSizes);
  title(sprintf('%s', stimulus.imNames{l}));
  set(gca, 'XTick', 1:length(distractor_layers));
  set(gca, 'XTickLabel', stimulus.layerNames);
  ylabel('Accuracy (proportion correct)')
  xlabel('Model layer to which distractors were feature-matched')
end

%% Save results
data.accuracy = accs;
data.numTrials = Ns;
data.stdErr = SEs;
data.layerNames = stimulus.layerNames;
data.imNames = stimulus.imNames;
data.poolSizes = stimulus.poolSizes;
data.origImDir = stimulus.origImDir;
data.stimDir = stimulus.stimDir;

save(sprintf('~/proj/texOdd/data/%s_behavior.mat', mglGetSID), '-struct', 'data');
keyboard
