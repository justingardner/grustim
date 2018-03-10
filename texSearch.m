function [ myscreen ] = texSearch( varargin )
%
% TEXTURE SEARCH 
%  Visual search task using textures 
%
%  Usage: texSearch(varargin)
%  Authors: Akshay Jagadeesh
%  Date: 02/27/2018
%

global stimulus

stimulus = struct;

%% Initialize Variables

% add arguments later
scan = 0;
plots = 0;
noeye = 0;
debug = 0;
getData = 0;
getArgs(varargin,{'getData=0', 'scan=0','plots=0','noeye=0','debug=0'});
stimulus.scan = scan;
stimulus.plots = plots;
stimulus.noeye = noeye;
stimulus.debug = debug;
stimulus.getData = getData;
clear localizer invisible scan noeye task test2


if stimulus.plots
    dispInfo(stimulus);
    myscreen = 0;
    return
end

if stimulus.scan
  warning('Not setup for scanning');
end

%% Stimulus parameters 
%% Open Old Stimfile
stimulus.counter = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/texSearch/%s',mglGetSID))
  % Directory exists, check for a stimfile
  files = dir(sprintf('~/data/texSearch/%s/1*mat',mglGetSID));

  if length(files) >= 1
    fname = files(end).name;
    
    s = load(sprintf('~/data/texSearch/%s/%s',mglGetSID,fname));
    stimulus.counter = s.stimulus.counter + 1;
    clear s;
    disp(sprintf('(texSearch) Data file: %s loaded.',fname));
  end
end
disp(sprintf('(texSearch) This is run #%i',stimulus.counter));

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
stimulus.responseKeys = [1 2 3 4]; % 
stimulus.responseKeys = [11 12 13 14];

% set colors
stimulus.colors.white = [1 1 1];
stimulus.colors.black = [0 0 0];
stimulus.colors.red = [1 0 0];
stimulus.colors.green = [0 1 0];
stimulus.colors.blue = [0 0 1];
stimulus.live.fixColor = stimulus.colors.blue;
stimulus.live.cueColor = stimulus.colors.black;

%% Setup Task

%%%%%%%%%%%%% PHASE ONE %%%%%%%%%%%%%%%%%

stimulus.curTrial(1) = 0;

task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;

%% Define stimulus timing
tTarg = 0.500;
isi = 0.500;
searchTime = 2.00;
stimDirectory = '~/proj/TextureSynthesis/stimuli';

% task waits for fixation on first segment
task{1}{1}.segmin = [inf tTarg isi searchTime .200];
task{1}{1}.segmax = [inf tTarg isi searchTime .200];
stimulus.seg = {};
stimulus.seg{1}.fix = 1;
stimulus.seg{1}.target = 2;
stimulus.seg{1}.isi = 3;
stimulus.seg{1}.search = 4;
stimulus.seg{1}.feedback = 5;

if stimulus.noeye==1
  task{1}{1}.segmin(1) = 0.5;
  task{1}{1}.segmax(1) = 0.5;
end

% Task important variables
task{1}{1}.imNames = {'rocks', 'tulips', 'leaves', 'fronds', 'cherries', 'clouds', 'bubbles', 'balls', 'forest', 'worms'};
task{1}{1}.layerNames = {'pool1', 'pool2', 'pool3', 'pool4'};
task{1}{1}.stimDir = stimDirectory;

% Trial parameters
task{1}{1}.parameter.targIm = 1:10;
task{1}{1}.parameter.layer = [1 2 3 4];
task{1}{1}.parameter.eccentricity = [5 8 11];

task{1}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
task{1}{1}.getResponse = zeros(size(task{1}{1}.segmin));
task{1}{1}.getResponse(stimulus.seg{1}.search)=1;
task{1}{1}.numTrials = 120;
task{1}{1}.random = 1;

if stimulus.scan
  task{1}{1}.synchToVol(stimulus.seg.ITI) = 1;
end

% Task trial parameters

% Task variables to be calculated later
task{1}{1}.randVars.calculated.targetPosition = NaN;
task{1}{1}.randVars.calculated.detected = 0; % did they see the grating
task{1}{1}.randVars.calculated.dead = 0;
task{1}{1}.randVars.calculated.visible = 1;

%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{1})
  [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);
end

%% EYE CALIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

% let the user know
disp(sprintf('(texSearch) Starting run number: %i.',stimulus.counter));

%% Main Task Loop

mglClearScreen(0.5); 
upFix(stimulus);

mglFlush
mglClearScreen(0.5); 
upFix(stimulus);

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

%%
%if ~isempty(task.lasttrial)
%  disp(sprintf('Last trial - image: %d, scaling: 0.%d, correct response: %d, subject''s response: %d', task.lasttrial.image, task.lasttrial.scaling, task.lasttrial.correctResponse, task.lasttrial.response));
%end
global stimulus

task.thistrial.dead = 0;
task.thistrial.detected = 0;
task.thistrial.visible = 1;
task.thistrial.response = 0;

stimulus.live.gotResponse = 0;
stimulus.curTrial(task.thistrial.thisphase) = stimulus.curTrial(task.thistrial.thisphase) + 1;

% directories
targetDir = '~/proj/TextureSynthesis/orig_ims';
distDir = '~/proj/TextureSynthesis/stimuli';

%% Load all 4 images for this trial
imName = task.imNames{task.thistrial.targIm};
layer = task.layerNames{task.thistrial.layer};

stimulus.live.target_image = genTexFromIm(imread(sprintf('%s/%s.jpg', targetDir, imName)));
stimulus.live.d1 = genTexFromIm(imread(sprintf('%s/v1/%s_%s_step_10000.jpg', distDir, layer, imName)));
stimulus.live.d2 = genTexFromIm(imread(sprintf('%s/v2/%s_%s_step_10000.jpg', distDir, layer, imName)));
stimulus.live.d3 = genTexFromIm(imread(sprintf('%s/v3/%s_%s_step_10000.jpg', distDir, layer, imName)));

% Select target position
task.thistrial.targetPosition = randi(4, 1);

% set response text
stimulus.live.responseText = mglText('1 or 2?');

% Disp trial parameters each trial
disp(sprintf('Trial %d - Image: %s, Layer: %s, Ecc: %d', task.trialnum, imName, layer, task.thistrial.eccentricity));

% Reset mouse to center of screen at start of every trial
mglSetMousePosition(960,540,1);
myscreen.flushMode = 0;
stimulus.live.eyeCount = 0;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)
%%

global stimulus

stimulus.live.triggerWaiting = 0;
if any(task.thistrial.thisseg==[stimulus.seg{task.thistrial.thisphase}.fix])
  stimulus.live.triggerWaiting = 1;
  stimulus.live.centered = 0;
  stimulus.live.triggerTime = 0;
  stimulus.live.lastTrigger = -1;
end

stimulus.live.eyeDead = 0;
stimulus.live.fix = 1;
stimulus.live.target = 0;
stimulus.live.search = 0;
stimulus.live.feedback = 0;

if task.thistrial.thisseg == stimulus.seg{task.thistrial.thisphase}.target
  stimulus.live.target = 1;
elseif task.thistrial.thisseg == stimulus.seg{task.thistrial.thisphase}.search
  stimulus.live.search = 1;
elseif task.thistrial.thisseg == stimulus.seg{task.thistrial.thisphase}.feedback
  stimulus.live.feedback = 1;
  stimulus.live.fix = 0;
end

% Select image parameters: size, eccentricity, and location
imSz = 5; % Size in Degrees of Visual Angle to display image
ecc = task.thistrial.eccentricity / sqrt(2); % 
locations = [-ecc ecc; ecc ecc; ecc -ecc; -ecc -ecc];
distLocations = setdiff(1:4, task.thistrial.targetPosition);

for i = 1:2
  mglClearScreen(0.5);
  if stimulus.live.target
    mglBltTexture(stimulus.live.target_image, [0 0 imSz imSz]);
  elseif stimulus.live.search
    mglBltTexture(stimulus.live.d1, [locations(distLocations(1), :), imSz, imSz]);
    mglBltTexture(stimulus.live.d2, [locations(distLocations(2), :), imSz, imSz]);
    mglBltTexture(stimulus.live.d3, [locations(distLocations(3), :), imSz, imSz]);
    mglBltTexture(stimulus.live.target_image, [locations(task.thistrial.targetPosition,:), imSz, imSz]);
  elseif stimulus.live.feedback
    if task.thistrial.response == task.thistrial.targetPosition
      if i == 1
        disp(sprintf('Correct! You responded: %i, Correct answer: %i', task.thistrial.response, task.thistrial.targetPosition));
      end
      upFix(stimulus, stimulus.colors.green);
    else
      if i == 1
        disp(sprintf('Incorrect! You responded: %i, correct was: %i', task.thistrial.response, task.thistrial.targetPosition));
      end
      upFix(stimulus, stimulus.colors.red);
    end
  end

  if stimulus.live.fix
    upFix(stimulus, stimulus.colors.blue);
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
if ~stimulus.noeye && ~any(task.thistrial.thisseg==[stimulus.seg{task.thistrial.thisphase}.fix]) && ~stimulus.scan
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
    disp('Starting trial--eye centered and space pressed.');
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
    task.thistrial.response = task.thistrial.whichButton - 10;
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

%%%
% Turns image into a texture
function tex = genTexFromIm(im)
r = flipud(im);
r(:,:,4) = 255;
rP = permute(r, [3 2 1]);
tex = mglCreateTexture(rP);  

function [trials] = totalTrials()
%%
% Counts trials + estimates the threshold based on the last 500 trials
% get the files list
files = dir(fullfile(sprintf('~/data/texSearch/%s/18*stim*.mat',mglGetSID)));
trials = 0;

for fi = 1:length(files)
    load(fullfile(sprintf('~/data/texSearch/%s/%s',mglGetSID,files(fi).name)));
    e = getTaskParameters(myscreen,task);
    e = e{1}; % why?!
    trials = trials + e.nTrials;
end

%%%%%%%%%%%%%%%%%%%%%%%
%       getTrialData       %
%%%%%%%%%%%%%%%%%%%%%%%
function data = getTrialData(rstimulus)
%%
files = dir(fullfile(sprintf('~/data/texSearch/%s/18*stim*.mat',mglGetSID)));
data = struct('subjResp', [], 'corrResp', [], 'scaling', [], 'image', [], 'synthPairs', [], 'whichRun', [], 'nTrials', []);
data.nRuns = length(files);
scaleIdx = [NaN NaN 1 2 3 4 5 NaN NaN 6];

for fi = 1:length(files)
  load(fullfile(sprintf('~/data/texSearch/%s/%s',mglGetSID,files(fi).name)));
  e = getTaskParameters(myscreen,task);
  if e{1}.nTrials>1
    % Set data for all trials
    data.subjResp = [data.subjResp e{1}.response];
    data.corrResp = [data.corrResp e{1}.parameter.correctResponse];
    data.synthPairs = [data.synthPairs e{1}.parameter.synthPair];
    data.whichRun = [data.whichRun fi*ones(1,e{1}.nTrials)];
    data.nTrials = [data.nTrials e{1}.nTrials];
    % Remap images
    im = e{1}.parameter.image;
    im(im==0) = 1; im(im==5) = 4;
    data.image = [data.image im];
    % Remap scalings
    data.scaling = [data.scaling scaleIdx(e{1}.parameter.scaling)];
  end
end

%%%%%%%%%%%%%%%%%%%%%%%
%    dispInfo    %
%%%%%%%%%%%%%%%%%%%%%%%
function dispInfo(rstimulus)
%%

% get the files list
files = dir(fullfile(sprintf('~/data/texSearch/%s/18*stim*.mat',mglGetSID)));

count = 1; 
data = struct('nTrials', 0, 'subj_resp', [], 'corr_resp', [], 'corr_trials', [],...
              'image', [], 'layer', [], 'ecc', [], 'reaction_time', [], 'nValTrials', 0);

for fi = 1:length(files)
  load(fullfile(sprintf('~/data/texSearch/%s/%s',mglGetSID,files(fi).name)));
  
  e = getTaskParameters(myscreen,task);
  if e{1}.nTrials>1
    
    subj_resp = e{1}.response-10;
    corr_resp = e{1}.randVars.targetPosition;
    data.run = stimulus.counter;
    data.subj_resp = [data.subj_resp subj_resp];
    data.corr_resp = [data.corr_resp corr_resp];
    data.corr_trials = [data.corr_trials subj_resp==corr_resp];
    data.reaction_time = [data.reaction_time e{1}.reactionTime];
    data.nTrials = data.nTrials + e{1}.nTrials;
    % Calculate number of valid trials by excluding eye movements and pool5
    data.nValTrials = data.nValTrials + sum(~isnan(e{1}.response)) - sum(e{1}.parameter.layer == 5);
    
    data.image = [data.image e{1}.parameter.targIm];
    data.layer = [data.layer e{1}.parameter.layer];
    data.ecc = [data.ecc e{1}.parameter.eccentricity];
    
  end
  count = count + 1;
end

%%
figure;
subplot(2,1,1);
all_eccs = unique(data.ecc);
ct = data.corr_trials;
colors = brewermap(length(all_eccs), 'Dark2');
y = [];
for i = 1:length(all_eccs)
  ei = all_eccs(i);
  y(i,:) = [nansum(ct(data.ecc==ei & data.layer==1)), nansum(ct(data.ecc==ei & data.layer==2)), nansum(ct(data.ecc==ei & data.layer==3)), nansum(ct(data.ecc==ei & data.layer==4))]/(data.nTrials/12);
  plot(1:4, y(i,:), '.', 'MarkerSize', 15, 'Color', colors(i,:)); hold on;
end
plot(1:4, nanmean(y,1), '.k', 'MarkerSize', 20);

se = @(x) 1.96*nanstd(x) / sqrt(length(x));
eb = [se(ct(data.layer==1)), se(ct(data.layer==2)), se(ct(data.layer==3)), se(ct(data.layer==4))];
errorbar(1:4, mean(y,1), eb, '.k');
legend('5 degrees', '8 degrees', '11 degrees', 'Mean');
title(sprintf('Accuracy as a function of distractor layer. nTrials=%i', data.nValTrials), 'FontSize', 18);
xlim([0 5]);ylim([0 1.2]);
xlabel('CNN Layer from which distractors were generated', 'FontSize', 16);
ylabel('Identification Accuracy', 'FontSize', 16);
hline(0.25, ':');
set(gca, 'FontSize', 14);
set(gca, 'XTick', 1:4);
set(gca, 'XTickLabel', {'Pool1', 'Pool2', 'Pool3', 'Pool4'});


subplot(2,1,2);
all_eccs = unique(data.ecc);
rt = data.reaction_time;
y2 = [];
for i = 1:length(all_eccs)
  ei = all_eccs(i);
  y2(i,:) = [nansum(rt(data.ecc==ei & data.layer==1)), nansum(rt(data.ecc==ei & data.layer==2)), nansum(rt(data.ecc==ei & data.layer==3)), nansum(rt(data.ecc==ei & data.layer==4))]/(data.nTrials/4);
  plot(1:4, y2(i,:), '.', 'MarkerSize', 15, 'Color', colors(i,:)); hold on;
end
plot(1:4, nanmean(y2,1), '.k', 'MarkerSize', 20);

eb = [se(rt(data.layer==1)), se(rt(data.layer==2)), se(rt(data.layer==3)), se(rt(data.layer==4))];
errorbar(1:4, mean(y2,1), eb, '.k');
legend('5 degrees', '8 degrees', '11 degrees', 'Mean');
title('Reaction Time as a function of distractor layer', 'FontSize', 18);
xlim([0 5]);
ylim([.1 .3]);
xlabel('CNN Layer from which distractors were generated', 'FontSize', 16);
ylabel('Reaction Time', 'FontSize', 16);
set(gca, 'FontSize', 14);
set(gca, 'XTick', 1:4);
set(gca, 'XTickLabel', {'Pool1', 'Pool2', 'Pool3', 'Pool4'});

%%

keyboard

mglClose;


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

