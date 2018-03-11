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
plotEye = 0;
stairImSz = 0;
varyImSz = 0;
analyzeStair = 0;
getArgs(varargin,{'varyImSz=0', 'analyzeStair=0', 'stairImSz=1','getData=0', 'plotEye=0', 'scan=0','plots=0','noeye=0','debug=0'});
stimulus.scan = scan;
stimulus.plots = plots;
stimulus.noeye = noeye;
stimulus.debug = debug;
stimulus.getData = getData;
stimulus.plotEye= plotEye;
stimulus.stairImSz= stairImSz;
stimulus.varyImSz = varyImSz;
stimulus.analyzeStair = analyzeStair;
clear localizer invisible scan noeye task test2

if stimulus.analyzeStair
    analyzeStaircase();
    return
end

if stimulus.plots
    dispInfo(stimulus);
    myscreen = 0;
    return
end

if stimulus.plotEye
  plotEyetraces(stimulus);
  return
end

if stimulus.stairImSz
  disp('Staircasing image size');
  stimulus.strcs = {};
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

%
stimulus.allSizes = [3, 4, 5, 6, 7, 8, 9, 10];

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
task{1}{1}.randVars.calculated.imSz = NaN;
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

% Staircase image size.
if stimulus.varyImSz
  task.thistrial.imSz = stimulus.allSizes(randi(length(stimulus.allSizes), 1));
elseif stimulus.stairImSz
  task.thistrial.imSz = stairImSize(task);
else
  task.thistrial.imSz = 5;
end

% set response text
stimulus.live.responseText = mglText('1 or 2?');

% Disp trial parameters each trial
disp(sprintf('Trial %d - Image: %s, Layer: %s, Ecc: %d, Size, %g', task.trialnum, imName, layer, task.thistrial.eccentricity, task.thistrial.imSz));

% Reset mouse to center of screen at start of every trial
mglSetMousePosition(960,540,1);
myscreen.flushMode = 0;
stimulus.live.eyeCount = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Run Staircase over imsize for each condition %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newImSz = stairImSize(task)

global stimulus

lastN = 4;
threshAcc = 0.75;
increment = 0.25;

eccs = [5, 8, 11];

% first figure out which trial index it is
trLayer = task.thistrial.layer;
trEcc = find(eccs==task.thistrial.eccentricity);
idxMtx = reshape(1:12, [3,4]);
trIdx = idxMtx(trEcc, trLayer);

% then update that element of stimulus.strcs
cnm = sprintf('cond%d_acc', trIdx);
fnm = sprintf('cond%d', trIdx);
if ~isfield(stimulus.strcs, fnm)
  stimulus.strcs.(fnm) = [];
  stimulus.strcs.(cnm) = [];
end


% Update the last trial's condition
if ~isempty(task.lasttrial)
  lastIdx = idxMtx(find(eccs==task.lasttrial.eccentricity), task.lasttrial.layer);
  lastAcc = (task.lasttrial.response == task.lasttrial.targetPosition);
  lnm = sprintf('cond%d_acc', lastIdx);
  % Update the correct condition's accuracy depending on if they got the last trial right or wrong.
  stimulus.strcs.(lnm) = [stimulus.strcs.(lnm) lastAcc];
end

% for each of 12 conditions, store a running log of accuracies, and imSizes
if length(stimulus.strcs.(fnm)) < 1
  % Load previous runs
  f = dir(fullfile(sprintf('~/data/texSearch/%s/18*.mat', mglGetSID)));
  lsz = 5;
  if length(f) > 0
    l = load(fullfile(sprintf('~/data/texSearch/%s/%s', mglGetSID, f(end).name)));
    if isfield(l.stimulus, 'strcs')
      lsz = l.stimulus.strcs.(fnm)(end);
      disp(sprintf('Previous runs threshold found.. Starting image size at %g', lsz));
    end
  end
  newImSz = lsz;
else
  % Get the image size on the last trial of this condition
  prevSz = stimulus.strcs.(fnm)(end);

  % Check a running average accuracy over the last 4 trials
  % If the running average is less than threshold (0.75), increment imsize, otherwise decrement.
  st = max(length(stimulus.strcs.(cnm)) - lastN, 1);
  if mean(stimulus.strcs.(cnm)(st:end)) < threshAcc
    newImSz = prevSz + increment;
  else
    newImSz = prevSz - increment;
  end
end

stimulus.strcs.(fnm) = [stimulus.strcs.(fnm) newImSz];
  
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
imSz = task.thistrial.imSz; % Size in Degrees of Visual Angle to display image
ecc = task.thistrial.eccentricity / sqrt(2); % 
locations = [-ecc ecc; ecc ecc; ecc -ecc; -ecc -ecc];
distLocations = setdiff(1:4, task.thistrial.targetPosition);

for i = 1:2
  mglClearScreen(0.5);
  if stimulus.live.target
    mglBltTexture(stimulus.live.target_image, [0 0 5 5]);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% analyze staircase %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function analyzeStaircase()

files = dir(fullfile(sprintf('~/data/texSearch/%s/18*.mat', mglGetSID)));

strDt = struct('nTrials', 0, 'subj_resp', [], 'corr_resp', [], 'corr_trials', [],...
              'image', [], 'layer', [], 'ecc', [], 'reaction_time', [], 'nValTrials', 0);

for fi = 1:length(files)
  l = load(fullfile(sprintf('~/data/texSearch/%s/%s',mglGetSID,files(fi).name)));
   
  e = getTaskParameters(l.myscreen, l.task);
  strcs = l.stimulus.strcs;
  
  for ci = 1:12
      cnm = sprintf('cond%d', ci);
      anm = sprintf('cond%d_acc', ci);

      if ~isfield(strDt, cnm)
          strDt.(cnm) = strcs.(cnm);
          strDt.(anm) = strcs.(anm);
      else
          strDt.(cnm) = [strDt.(cnm) strcs.(cnm)];
          strDt.(anm) = [strDt.(anm) strcs.(anm)];
      end
  end
  
  if e{1}.nTrials>1
    
    subj_resp = e{1}.response-10;
    corr_resp = e{1}.randVars.targetPosition;
    strDt.run = l.stimulus.counter;
    strDt.subj_resp = [strDt.subj_resp subj_resp];
    strDt.corr_resp = [strDt.corr_resp corr_resp];
    strDt.corr_trials = [strDt.corr_trials subj_resp==corr_resp];
    strDt.reaction_time = [strDt.reaction_time e{1}.reactionTime];
    strDt.nTrials = strDt.nTrials + e{1}.nTrials;
    % Calculate number of valid trials by excluding eye movements and pool5
    strDt.nValTrials = strDt.nValTrials + sum(~isnan(e{1}.response)) - sum(e{1}.parameter.layer == 5);
    
    strDt.image = [strDt.image e{1}.parameter.targIm];
    strDt.layer = [strDt.layer e{1}.parameter.layer];
    strDt.ecc = [strDt.ecc e{1}.parameter.eccentricity];
    
  end
  
end

%% Plot threshold staircase changes over time
conds = reshape(1:12, [3, 4]);
figure;
colors = brewermap(4, 'Dark2');

thresh = nan(3,4);

for ci = 1:12
    [eccI, layI] = find(conds==ci);
    
    condi = sprintf('cond%d', ci);
    plot(1:length(strDt.(condi)), strDt.(condi), '*-', 'Color', colors(layI,:)); hold on;
    thresh(eccI, layI) = mean(strDt.(condi)(end-4:end));
end
title('Threshold changes over the course of all trials', 'FontSize', 18);
xlabel('Trial count');
ylabel('Threshold');

%% Plot thresholds as a function of eccentricity.
figure;
all_layers = 1:4;
colors = brewermap(length(all_layers), 'Dark2');
y = [];
ct = strDt.corr_trials;
for i = 1:length(all_layers)
  plot([5 8 11], thresh(:,i)', '.', 'MarkerSize', 15, 'Color', colors(i,:)); hold on;
  y(i,:) = [nansum(ct(strDt.ecc==5 & strDt.layer==i)), nansum(ct(strDt.ecc==8 & strDt.layer==i)), nansum(ct(strDt.ecc==11 & strDt.layer==i))]/(strDt.nValTrials/12);
end
for i = 1:length(all_layers)
   X = [ones(3,1) [5 8 11]'];
   b = X \ thresh(:,i);
   yPred = X*b;
   plot([5 8 11], yPred, '-', 'Color', colors(i,:));
end
    
legend('Pool1', 'Pool2', 'Pool3', 'Pool4');
title(sprintf('Size threshold as a function of eccentricity. nTrials=%i', strDt.nValTrials), 'FontSize', 18);
xlim([3 13]);ylim([0 15]);
xlabel('Eccentricity', 'FontSize', 16);
ylabel('Size Threshold', 'FontSize', 16);
set(gca, 'FontSize', 14);

%%
keyboard


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% plot Eye Traces %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotEyetraces(stimulus)

cdir = pwd;

cd(sprintf('~/data/texSearch/%s', mglGetSID));

% get files list
files = dir(fullfile(sprintf('~/data/texSearch/%s/18*.mat', mglGetSID)));

idata = struct('xPos', {}, 'yPos', {}, 'pupil', {}, 'time', {}, 'isVal', [], 'respLoc', [], 'imLoc', []);
idata(1).xPos{1} = 0;

for fi = 1:length(files)
  fnm = files(fi).name(1:end-4);

  trace = getTaskEyeTraces(fnm);

  %idata.xPos = [idata.xPos; trace.eye.xPos];
  %idata.yPos = [idata.yPos; trace.eye.yPos];
  idata.xPos{fi} = trace.eye.xPos;
  idata.yPos{fi} = trace.eye.yPos;
  idata.pupil{fi} = trace.eye.pupil;
  idata.time{fi} = trace.eye.time;
  
  %idata.pupil = [idata.pupil; trace.eye.pupil];
  %idata.time = [idata.time; trace.eye.time];

  idata.isVal = [idata.isVal; ~isnan(trace.response)];

  idata.respLoc = [idata.respLoc; trace.response-10];
  idata.imLoc = [idata.imLoc; trace.randVars.targetPosition];

end

% Now, plot the eye traces by condition maybe?

keyboard

%%
midX = nan(120, length(files));
midY = nan(120, length(files));

for ri = 1:length(files)
  respLocs = idata.respLoc(ri,:);
  imLocs = idata.imLoc(ri,:);

  midX(:,ri) = median(idata.xPos{ri}, 2);
  midY(:,ri) = median(idata.yPos{ri}, 2);
end


%%%%%%%%%%%%%%%%%%%%%%%
%    dispInfo    %
%%%%%%%%%%%%%%%%%%%%%%%
function data = dispInfo(rstimulus)
%%

% get the files list
files = dir(fullfile(sprintf('~/data/texSearch/%s/18*stim*.mat',mglGetSID)));

count = 1; 
data = struct('nTrials', 0, 'subj_resp', [], 'corr_resp', [], 'corr_trials', [],...
              'image', [], 'layer', [], 'ecc', [], 'reaction_time', [], 'nValTrials', 0, 'imSz', []);

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
    data.imSz = [data.imSz e{1}.randVars.imSz];
    % Calculate number of valid trials by excluding eye movements and pool5
    data.nValTrials = data.nValTrials + sum(~isnan(e{1}.response)) - sum(e{1}.parameter.layer == 5);
    
    data.image = [data.image e{1}.parameter.targIm];
    data.layer = [data.layer e{1}.parameter.layer];
    data.ecc = [data.ecc e{1}.parameter.eccentricity];
    
  end
  count = count + 1;
end

%% Plot Accuracy and Reaction Time as a function of distractor layer.
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

%% Accuracy Vs Image Size (for different eccentricities)
all_layers = 1:4;
allSz = unique(data.imSz);
colors = brewermap(3, 'Dark2');
figure;
y4 = [];
ctd = @(x) nanmean(ct(data.imSz==x));
ctd2 = @(x,y) nanmean(ct(data.imSz==x & data.ecc==all_eccs(y)));
plot(allSz, [ctd(3) ctd(4) ctd(5) ctd(6) ctd(7) ctd(8) ctd(9) ctd(10)], '.k', 'MarkerSize', 15); hold on;
for i = 1:3
  y4(i,:) = [ctd2(3,i) ctd2(4,i) ctd2(5,i) ctd2(6,i) ctd2(7,i) ctd2(8,i) ctd2(9,i) ctd2(10,i)];
  plot(allSz, y4(i,:), '.', 'MarkerSize', 15, 'Color', colors(i,:)); hold on;
end
for i = 1:3
   X = [ones(length(allSz),1) allSz'];
   b = X \ y4(i,:)';
   yPred = X*b;
   plot(allSz, yPred, '-', 'Color', colors(i,:));
end
  
xlim([2, 11]);
ylim([.5, 1]);
legend({'Mean', '5 degrees', '8 degrees', '11 degrees'});
title(sprintf('Performance at different eccentricities as a function of image size. N=%d', data.nValTrials), 'FontSize', 16);
xlabel('Image Size (degrees)', 'FontSize', 14)
ylabel('Accuracy (% correct)', 'FontSize', 14);
set(gca, 'XTick', allSz);


%% Accuracy Vs Image Size (for different distractor layers)
all_layers = 1:4;
allSz = unique(data.imSz);
figure;
y4 = [];
ctd = @(x) nanmean(ct(data.imSz==x));
ctd2 = @(x,y) nanmean(ct(data.imSz==x & data.layer==y));
plot(allSz, [ctd(3) ctd(4) ctd(5) ctd(6) ctd(7) ctd(8) ctd(9) ctd(10)], '.k', 'MarkerSize', 15); hold on;
for i = 1:4
  y4(i,:) = [ctd2(3,i) ctd2(4,i) ctd2(5,i) ctd2(6,i) ctd2(7,i) ctd2(8,i) ctd2(9,i) ctd2(10,i)];
  plot(allSz, y4(i,:), '.', 'MarkerSize', 15, 'Color', colors(i,:)); hold on;
end
for i = 1:4
   X = [ones(length(allSz),1) allSz'];
   b = X \ y4(i,:)';
   yPred = X*b;
   plot(allSz, yPred, '-', 'Color', colors(i,:));
  end
  
xlim([2, 11]);
ylim([.5, 1.2]);
legend({'Mean', 'Pool1', 'Pool2', 'Pool3', 'Pool4'});
title(sprintf('Performance at different distractor layers as a function of image size. N=%d', data.nValTrials), 'FontSize', 16);
xlabel('Image Size (degrees)', 'FontSize', 14)
ylabel('Accuracy (% correct)', 'FontSize', 14);
set(gca, 'XTick', allSz);
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

