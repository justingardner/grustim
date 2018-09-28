function [ myscreen ] = texAtt( varargin )
%
% TEXTURE SEARCH 
%  Visual search task using textures 
%
%  Usage: texAtt(varargin)
%  Authors: Akshay Jagadeesh
%  Date: 02/27/2018
%

global stimulus

stimulus = struct;

%% Initialize Variables

% add arguments later
plots = 0;
noeye = 0;
getArgs(varargin,{'plots=0','noeye=0'});
stimulus.plots = plots;
stimulus.noeye = noeye;
clear noeye plots

%% Stimulus parameters 
%% Open Old Stimfile
stimulus.counter = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/texAtt/%s',mglGetSID))
  % Directory exists, check for a stimfile
  files = dir(sprintf('~/data/texAtt/%s/1*mat',mglGetSID));

  if length(files) >= 1
    fname = files(end).name;
    
    s = load(sprintf('~/data/texAtt/%s/%s',mglGetSID,fname));
    stimulus.counter = s.stimulus.counter + 1;
    clear s;
    disp(sprintf('(texAtt) Data file: %s loaded.',fname));
  end
end
disp(sprintf('(texAtt) This is run #%i',stimulus.counter));

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
stimulus.responseKeys = [11 12 13 14];

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

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% [akshay] (1) Edit the code below to specify the length 
%         and name of each segment
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Define stimulus timing
task{1}.segmin = [inf, 0.2, 0.5, 1.00, .200];
task{1}.segmax = [inf, 0.2, 0.5, 1.00, .200];
stimulus.seg = {};
stimulus.seg.fix = 1;
stimulus.seg.target = 2;
stimulus.seg.isi = 3;
stimulus.seg.search = 4;
stimulus.seg.feedback = 5;

if stimulus.noeye==1
  task{1}.segmin(1) = 0.5;
  task{1}.segmax(1) = 0.5;
end

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% [akshay] (2) Edit the code below to store some important variables in the 
%               stimulus struct. Remember that stimulus is just a place to store
%               useful variables -- it doesn't really interact with the mgl code.
%         - Important variables such as a list of the images you want to load, and the directory they're saved in
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% Task important variables
stimulus.imNames = {'balls', 'beansalad', 'biryani', 'bubbles', 'cherries', 'clouds', 'crowd', 'dahlias', 'fireworks', 'leaves', 'noodles', 'rocks', 'tulips', 'worms', 'zebras', 'bananas', 'bark', 'bison', 'blossoms', 'blotch', 'braids', 'bricks', 'bubbly', 'bumpy', 'crystals', 'dalmatians', 'ducks', 'face', 'frills', 'fur', 'galaxy', 'gourds', 'grass', 'honeycomb', 'lace', 'marbled', 'marbles', 'monarchs', 'paisley', 'pears', 'phlox', 'rorschach', 'spiky', 'splotchy', 'stars', 'succulent', 'tiles'};
stimulus.layerNames = {'pool1', 'pool2', 'pool4'};
stimulus.stimDir = '~/proj/TextureSynthesis/stimuli/rf_stim';
stimulus.imSize = 6;
stimulus.eccentricity = 10;


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% [akshay] (3) Edit the code below to specify task trial parameters which will be selected in each
%         trial. For example, the image family or the layer can be randomly assigned.
%        -  MGL will attempt to assign trial parameters such that there are equal number of each trial type.
%        -  Remember that in later functions, these variables can be access from task.thistrial.layer or whatever.   
%        -  The getResponse is a list of length nSegs and is used to specify which segments MGL should expect to be 
%           listening for a response. 
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Trial parameters
task{1}.parameter.imgFam = 1:length(stimulus.imNames);
task{1}.parameter.layer = 1:length(stimulus.layerNames);

task{1}.synchToVol = zeros(size(task{1}.segmin));
task{1}.getResponse = zeros(size(task{1}.segmin));
task{1}.getResponse(stimulus.seg.search)=1;

task{1}.numTrials = 144;
task{1}.random = 1;


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% [akshay] (4) Edit the code below to define and initialize variables that you may want to 
%          calculate later on (e.g. in startTrialCallback or startSegmentCallback).
%       - you may want to create a variable here to keep track of which location is the cued location.
%       - or to keep track of what the response of the subject is.
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Task variables to be calculated later
task{1}.randVars.calculated.targetPosition = NaN;
task{1}.randVars.calculated.detected = 0; % did they see the grating
task{1}.randVars.calculated.dead = 0;
task{1}.randVars.calculated.visible = 1;

%%%%%%%%%%%%%%%% MGL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialize task and begin update loop
[task{1}, myscreen] = initTask(task{1},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);

% Run the eye calibration
myscreen = eyeCalibDisp(myscreen);

% let the user know
disp(sprintf('(texAtt) Starting run number: %i.',stimulus.counter));

%% Main Task Loop
mglClearScreen(0.5); 
upFix(stimulus);

mglFlush
mglClearScreen(0.5); 
upFix(stimulus);

phaseNum = 1;
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

global stimulus

task.thistrial.dead = 0;
task.thistrial.detected = 0;
task.thistrial.visible = 1;
task.thistrial.response = 0;

stimulus.live.gotResponse = 0;
stimulus.curTrial(task.thistrial.thisphase) = stimulus.curTrial(task.thistrial.thisphase) + 1;

keyboard

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% [akshay] (5) Modify the code below to specify the directories from which we will be loading images.
%       - This is where we will preload all 8 images (4 images x 2 intervals) at the start of each trial.
%       - Use imread to load the images and then pass that into genTexFromIm, and store the output into a field of
%         the stimulus.live struct.
%
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% directories
targetDir = '~/proj/TextureSynthesis/stimuli/orig_ims';
distDir = stimulus.stimDir;

% Select image
task.thistrial.targIm = randi(length(stimulus.imNames),1);

%% Load all 4 images for this trial
imName = stimulus.imNames{task.thistrial.targIm};
layer = stimulus.layerNames{task.thistrial.layer};
rfSz = stimulus.rfNames{task.thistrial.rfSize};

stimulus.live.target_image = genTexFromIm(imread(sprintf('%s/%s.jpg', targetDir, imName)));
stimulus.live.d1 = genTexFromIm(imread(sprintf('%s/s1/%s_%s_%s_step_10000.jpg', distDir, rfSz, layer, imName)));
stimulus.live.d2 = genTexFromIm(imread(sprintf('%s/s2/%s_%s_%s_step_10000.jpg', distDir, rfSz, layer, imName)));
stimulus.live.d3 = genTexFromIm(imread(sprintf('%s/s3/%s_%s_%s_step_10000.jpg', distDir, rfSz, layer, imName)));

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% [akshay] (6) Edit and add code below to set the values of the randVars.calculated variables which you had defined above.
%       - For example, you may want to use this space to choose which of the 4 locations will be the response-cued direction,
%       as well as whether this trial is a distributed cue or a focal cue trial.
%
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Select target position and target size.
task.thistrial.targetPosition = randi(4, 1);






% Disp trial parameters each trial
disp(sprintf('Trial %d - Image %s, Layer %s', task.trialnum, imName, layer));

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
locations = [-ecc ecc; ecc ecc; ecc -ecc; -ecc -ecc];
distLocations = setdiff(1:4, task.thistrial.targetPosition);

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% [akshay] (6) Edit the code below to specify what to draw onto the screen at the start of each segment
%          - Recall that MGL has a front buffer and a back buffer, so the code below loops through
%
%
%
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
for i = 1:2
  mglClearScreen(0.5);
  if task.thistrial.thisseg == stimulus.seg.target
    mglBltTexture(stimulus.live.target_image, [0 0 imSz imSz]);
  elseif task.thistrial.thisseg == stimulus.seg.search
    mglBltTexture(stimulus.live.d1, [locations(distLocations(1), :), imSz, imSz]);
    mglBltTexture(stimulus.live.d2, [locations(distLocations(2), :), imSz, imSz]);
    mglBltTexture(stimulus.live.d3, [locations(distLocations(3), :), imSz, imSz]);
    mglBltTexture(stimulus.live.target_image, [locations(task.thistrial.targetPosition,:), imSz, imSz]);
  end

  % If in feedback segment, change the color of the cross to green or red.
  if task.thistrial.thisseg == stimulus.seg.feedback
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
  else % If not in feedback segment, just draw a blue fixation cross.
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
files = dir(fullfile(sprintf('~/data/texAtt/%s/18*stim*.mat',mglGetSID)));
trials = 0;

for fi = 1:length(files)
    load(fullfile(sprintf('~/data/texAtt/%s/%s',mglGetSID,files(fi).name)));
    e = getTaskParameters(myscreen,task);
    e = e{1}; % why?!
    trials = trials + e.nTrials;
end

%%%%%%%%%%%%%%%%%%%%%%%
%       getTrialData       %
%%%%%%%%%%%%%%%%%%%%%%%
function data = getTrialData(rstimulus)
%%
files = dir(fullfile(sprintf('~/data/texAtt/%s/18*stim*.mat',mglGetSID)));
data = struct('subjResp', [], 'corrResp', [], 'scaling', [], 'image', [], 'synthPairs', [], 'whichRun', [], 'nTrials', []);
data.nRuns = length(files);
scaleIdx = [NaN NaN 1 2 3 4 5 NaN NaN 6];

for fi = 1:length(files)
  load(fullfile(sprintf('~/data/texAtt/%s/%s',mglGetSID,files(fi).name)));
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
%%%%%%% plot Eye Traces %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotEyetraces(stimulus)

cdir = pwd;

cd(sprintf('~/data/texAtt/%s', mglGetSID));

% get files list
files = dir(fullfile(sprintf('~/data/texAtt/%s/18*.mat', mglGetSID)));

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
files = dir(fullfile(sprintf('~/data/texAtt/%s/18*stim*.mat',mglGetSID)));

count = 1; 
data = struct('nTrials', 0, 'subj_resp', [], 'corr_resp', [], 'corr_trials', [],...
              'image', [], 'layer', [], 'ecc', [], 'reaction_time', [], 'nValTrials', 0, 'rf_size', [], 'imSz', [], 'accByRuns', []);

for fi = 1:length(files)
  load(fullfile(sprintf('~/data/texAtt/%s/%s',mglGetSID,files(fi).name)));
  
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
    
    data.image = [data.image e{1}.randVars.targIm];
    data.layer = [data.layer e{1}.parameter.layer];
    data.rf_size = [data.rf_size e{1}.parameter.rfSize];
    data.ecc = [data.ecc e{1}.parameter.eccentricity];
    
    data.accByRuns = [data.accByRuns nanmean(subj_resp==corr_resp)];
        
  end
  count = count + 1;
end

% fitline function
fitLine = @(x,y) polyval(polyfit(x,y,1),x);

%% Plot performance by runs
figure; plot(data.accByRuns, '.', 'MarkerSize', 25); hold on;
plot(1:length(data.accByRuns), fitLine(1:length(data.accByRuns), data.accByRuns), '-');
title(sprintf('%s: Performance across runs', mglGetSID), 'FontSize', 18); 
ylabel('Accuracy', 'FontSize', 16); xlabel('Run Number', 'FontSize', 16);

%% Plot accuracy and reaction time as a function of distractor layer & RF Size
figure;
set(gcf, 'Position', [436, 460, 1038, 624]);

all_eccs = unique(data.ecc);
all_RFs = unique(data.rf_size);
all_layers = unique(data.layer);

nLayers = length(unique(data.layer));
x1 = 1:nLayers;

ct = data.corr_trials;

colors = brewermap(length(all_RFs), 'Dark2');
for k = 1:length(all_eccs)
  subplot(2,3,k);
  y = [];
  for i = 1:length(all_RFs)
    ei = all_RFs(i);
    for j = 1:nLayers
        y(i,j) = nanmean(ct(data.rf_size==ei & data.layer == j & data.ecc==all_eccs(k)));
    end
    plot(x1, y(i,:), '.', x1, fitLine(x1, y(i,:)), '-', 'Color', colors(i,:)); hold on;
  end

  h = findobj(gca, 'Type', 'line');
  legend(h(length(h)-1:-2:1), stimulus.rfNames, 'Location', 'southwest');
  nTr = sum(~isnan(ct(data.ecc==all_eccs(k))));
  if k ~= 1
    title(sprintf('Ecc=%i deg. nTrials=%i', all_eccs(k), nTr), 'FontSize', 18);
  else
    title(sprintf('%s: Ecc = %i deg.', mglGetSID, all_eccs(k)), 'FontSize', 18);
  end
  xlim([0 nLayers+1]); xlabel('CNN Layer', 'FontSize', 16);
  ylim([0 1]); ylabel('Identification Accuracy', 'FontSize', 16);
  hline(0.25, ':');
  set(gca, 'FontSize', 14);
  set(gca, 'XTick', 1:nLayers); set(gca, 'XTickLabel', stimulus.layerNames);
end

% calculate RF size in degrees
all_RF_deg = data.imSz(1) ./ cellfun(@(x) str2num(x(1)), stimulus.rfNames);

colors = brewermap(length(all_layers), 'Dark2');
for k = 1:length(all_eccs)
  subplot(2,3,k+length(all_eccs));
  y = [];
  for i = 1:length(all_layers)
      li = all_layers(i);
      for j = 1:length(all_RFs)
          ei = all_RFs(j);
          y(i,j) = nanmean(ct(data.rf_size==ei & data.layer == li & data.ecc==all_eccs(k)));
      end
      plot(all_RF_deg, y(i,:), '.', 'MarkerSize', 15, 'Color', colors(i,:)); hold on;
      plot(all_RF_deg, fitLine(all_RF_deg, y(i,:)), '-', 'Color', colors(i,:)); hold on;
  end

  h = findobj(gca, 'Type', 'line');
  legend(h(length(h)-1:-2:1), stimulus.layerNames, 'Location', 'southeast');

  xlim([0 max(all_RF_deg)+1]); ylim([0 1]);
  xlabel('Gram RF Size (dva)', 'FontSize', 16);
  ylabel('Accuracy', 'FontSize', 16);
  if k~=1
    title(sprintf('RF Size vs Accuracy - Ecc=%i deg', all_eccs(k)), 'FontSize', 18);
  else
    title(sprintf('Ecc = %i deg', all_eccs(k)), 'FontSize', 18);
  end
  hline(0.25, ':');
  set(gca, 'FontSize', 14);
end

%saveas(gcf, sprintf('~/proj/TextureSynthesis/Figures/%s_results.png', mglGetSID));

%% Eccentricity Plots
figure; set(gcf, 'Position', [974, 344, 389, 738]);
subplot(2,1,1);

%  First plot all RF Sizes
y = []; colors = brewermap(length(all_RFs), 'Dark2');
for i = 1:length(all_RFs)
  for j = 1:length(all_eccs)
    y(i,j) = nanmean(ct(data.rf_size==all_RFs(i) & data.ecc==all_eccs(j)));
  end
  plot(all_eccs, y(i,:), '.', 'MarkerSize', 15, 'Color', colors(i,:)); hold on;
  plot(all_eccs, fitLine(all_eccs, y(i,:)), '-', 'Color', colors(i,:));
end
h = findobj(gca, 'Type', 'line');
legend(h(length(h)-1:-2:1), stimulus.rfNames, 'Location', 'southeast');
xlim([min(all_eccs)-1, max(all_eccs)+1]); ylim([0 1]);
xlabel('Eccentricity (degs)', 'FontSize', 16);
ylabel('Accuracy', 'FontSize', 16);
title(sprintf('%s: Pooling Region Size', mglGetSID), 'FontSize', 18);
hline(0.25, ':');
set(gca, 'FontSize', 14);

subplot(2,1,2);

%  Then, plot Sizes
y = []; colors = brewermap(length(all_layers), 'Dark2');
for i = 1:length(all_layers)
  for j = 1:length(all_eccs)
    y(i,j) = nanmean(ct(data.layer==all_layers(i) & data.ecc==all_eccs(j)));
  end
  plot(all_eccs, y(i,:), '.', 'MarkerSize', 15, 'Color', colors(i,:)); hold on;
  plot(all_eccs, fitLine(all_eccs, y(i,:)), '-', 'Color', colors(i,:));
end
h = findobj(gca, 'Type', 'line');
legend(h(length(h)-1:-2:1), stimulus.layerNames, 'Location', 'southeast');
xlim([min(all_eccs)-1, max(all_eccs)+1]); ylim([0 1]);
xlabel('Eccentricity (degs)', 'FontSize', 16);
ylabel('Accuracy', 'FontSize', 16);
title(sprintf('Feature Complexity'), 'FontSize', 18);
hline(0.25, ':');
set(gca, 'FontSize', 14);



%%
keyboard
%% Plot individual images -- Accuracy vs Gram RF Size
figure;
set(gcf, 'Position', [680, 268, 1005, 830]);
all_ims = unique(data.image);
imnames = stimulus.imNames;

% Calculate image slopes
imslopes = nan(length(all_ims), length(all_layers));

all_RFs = unique(data.rf_size);
all_layers = unique(data.layer);
ct = data.corr_trials;
colors = brewermap(length(all_layers), 'Dark2');

% Labels for x axis and legend
xLabels = stimulus.rfNames;
lnLabels = stimulus.layerNames;

% calculate RF size in degrees
all_RF_deg = data.imSz(1) ./ cellfun(@(x) str2num(x(1)), stimulus.rfNames);

for imi = 1:length(all_ims)
  im = all_ims(imi);
  subplot(7,7,imi);

  y = [];
  for i = 1:length(all_layers)
    li = all_layers(i);
    for j = 1:length(all_RFs)
        ei = all_RFs(j);
        y(i,j) = nanmean(ct(data.rf_size==ei & data.layer == li & data.image == im));
    end
    plot(all_RF_deg, y(i,:), '.', 'MarkerSize', 15, 'Color', colors(i,:)); hold on;
    plot(all_RF_deg, fitLine(all_RF_deg, y(i,:)), '-', 'Color', colors(i,:));

    p = polyfit(all_RF_deg, y(i,:), 1);
    imslopes(imi,i) = p(1);
  end

  nTrials = sum(~isnan(data.subj_resp(data.image==im)));
  if imi == length(all_ims)
    h = findobj(gca, 'Type', 'line');
    legend(h(length(h)-1:-2:1), lnLabels, 'Location', 'bestoutside');
  end
  title(sprintf('%s: nTrials=%i', imnames{imi}, nTrials), 'FontSize', 16);
  xlim([0 max(all_RF_deg)+1]);ylim([0 1]);
  xlabel('Gram RF Size (degrees)', 'FontSize', 12);
  ylabel('Accuracy', 'FontSize', 12);
  hline(0.25, ':');
  set(gca, 'FontSize', 12);

end

%slopes.imnames = stimulus.imnames; 
%slopes.imslopes = imslopes;
%save('~/proj/TextureSynthesis/tmp/imslopes.mat', '-struct', 'slopes');

%% Plot individual images -- Accuracy vs Layer
figure; 
set(gcf, 'Position', [680, 268, 1005, 830]);

% Calculate image slopes
imslopes_layer = nan(length(all_ims), length(all_layers));

for imi = 1:length(all_ims)
  im = all_ims(imi);
  subplot(4,4,imi);

  all_RFs = unique(data.rf_size);
  nLayers = length(unique(data.layer));
  ct = data.corr_trials;
  colors = brewermap(length(all_RFs), 'Dark2');
  x1 = 1:nLayers;
  y = [];
  for i = 1:length(all_RFs)
    ei = all_RFs(i);
    for j = 1:nLayers
        y(i,j) = nanmean(ct(data.rf_size==ei & data.layer == j & data.image ==im));
    end
    plot(x1, y(i,:), '.', x1, fitLine(x1, y(i,:)), '-', 'Color', colors(i,:)); hold on;
    p = polyfit(x1, y(i,:), 1);
    imslopes_layer(imi,i) = p(1);
  end
  %plot(1:4, nanmean(y,1), '.k', 'MarkerSize', 20);

  if imi == length(all_ims) 
    h = findobj(gca, 'Type', 'line');
    legend(h(length(h)-1:-2:1), stimulus.rfNames, 'Location', 'bestoutside');
  end
  
  title(sprintf('%s: nTrials=%i', stimulus.imNames{imi}, data.nValTrials), 'FontSize', 16);
  xlim([0 nLayers+1]);ylim([0 1]);
  xlabel('CNN Layer Matched', 'FontSize', 12);
  ylabel('Identification Accuracy', 'FontSize', 12);
  hline(0.25, ':');
  set(gca, 'FontSize', 12);
  set(gca, 'XTick', 1:nLayers);
  set(gca, 'XTickLabel', stimulus.layerNames);
end

slopes = struct();
slopes.layernames = stimulus.layerNames;
slopes.imnames = stimulus.imNames; 
slopes.imslopes = imslopes_layer;
%save('~/proj/TextureSynthesis/tmp/imslopes_layer.mat', '-struct', 'slopes');

%%
keyboard
return


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

