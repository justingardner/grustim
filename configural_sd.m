function [ myscreen ] = configural( varargin )
% 
%  2-AFC task using configural shape stimuli
%
%  Usage: configural(varargin)
%  Authors: Akshay Jagadeesh
%  Date: 10/19/2021
%

global stimulus

stimulus = struct;

%% Initialize Variables
getArgs(varargin,{'noeye=1', 'analyze=0', 'fixate=0', 'smallStim=0', 'fastStim=0'}, 'verbose=1');
stimulus.noeye = noeye;
stimulus.analyze = analyze;
stimulus.fixate = fixate;
stimulus.smallStim = smallStim;
stimulus.fastStim = fastStim;
clear noeye analyze fixate smallStim fastStim

if stimulus.analyze
  analyzeData();
  return
end

%% Stimulus parameters 
if stimulus.fastStim == 1
  stimulus.stimLen = 0.200;
else
  stimulus.stimLen = 2.0;
end

if stimulus.smallStim == 1
  stimulus.stimSize = 2;
else
  stimulus.stimSize = 6;
end
stimulus.eccentricity = 4;
stimulus.base_seeds = {'008'};
stimulus.surfaces = [0,1,2];
stimulus.structures = [0,1];
stimulus.partlocations = [0,1];
stimulus.partstructures = [0,1];
stimulus.views = [0,1,2];
stimulus.conditions = {'surface', 'structure', 'partstructure', 'partlocation'};

%% Setup Screen
myscreen = initScreen('VPixx2');

% set background to grey
myscreen.background = 0.5;

%% Initialize Stimulus
myscreen = initStimulus('stimulus',myscreen);
  
% Set response keys
stimulus.responseKeys = [11 14 12];

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
task{1}.segmin = [inf, stimulus.stimLen, 0.2, stimulus.stimLen, 0.2];
task{1}.segmax = [inf, stimulus.stimLen, 0.2, stimulus.stimLen, 0.2];
stimulus.seg = {};
stimulus.seg.fix = 1;
stimulus.seg.stim1 = 2;
stimulus.seg.ISI = 3;
stimulus.seg.stim2 = 4;
stimulus.seg.feedback = 5;
% Set fixation length to constant if not eye tracking.
if stimulus.noeye==1
  task{1}.segmin(1) = 0.2;
  task{1}.segmax(1) = 0.2;
end

% Set when to synchtovol and getResponse
task{1}.synchToVol = zeros(size(task{1}.segmin));
task{1}.getResponse = zeros(size(task{1}.segmin));
task{1}.getResponse(stimulus.seg.stim2)=1;

%%% Stimulus Variables
stimulus.live.mask = imread('~/proj/TextureSynthesis/stimuli/Flattop8.tif');

% Make numTrials some multiple of number of TrialTypes 
task{1}.numTrials = 240;
task{1}.random = 1;

%%% Task Variables
% Trial parameters, which will be preassigned to be block randomized.
task{1}.parameter.base_seed = 1:length(stimulus.base_seeds);
task{1}.parameter.sample_surface = stimulus.surfaces;
task{1}.parameter.sample_structure = stimulus.structures;
task{1}.parameter.sample_partlocation = stimulus.partlocations;
task{1}.parameter.sample_partstructure = stimulus.partstructures;
task{1}.parameter.same = [0,1]; % 1: same, 0: different

task{1}.randVars.uniform.sample_view = stimulus.views;

% Trial variables, which will be calculated at the start of each trial.
task{1}.randVars.calculated.sample2_view = NaN;
task{1}.randVars.uniform.nonmatch_dimension = [1,2,3,4]; % which dimension does nonmatch vary in - 1:surface, 2:structure, 3:partlocation, 4:partstructure

task{1}.randVars.calculated.sample2_position = NaN;
task{1}.randVars.calculated.sample2_surface = NaN;
task{1}.randVars.calculated.sample2_structure = NaN;
task{1}.randVars.calculated.sample2_partlocation = NaN;
task{1}.randVars.calculated.sample2_partstructure = NaN;

% Task variables to keep track of the status of each trial
task{1}.randVars.calculated.detected = 0; % did they see the grating
task{1}.randVars.calculated.dead = 0;
task{1}.randVars.calculated.correct = NaN;

%% Preload images
%if ~exist(presavedStimLoc) % on the first time, load each image, convert to mgl texture, and save it to a struct.
stims = struct();
disppercent(-inf, sprintf('Preloading %i images', length(stimulus.base_seeds)*length(stimulus.surfaces)*length(stimulus.structures)*length(stimulus.partlocations)*length(stimulus.partstructures)*length(stimulus.views)));
stim_dir = '~/proj/oddity/pilot1_hybrid';
stimulus.stimDir = stim_dir;
for bsi = 1:length(stimulus.base_seeds)
  bs = stimulus.base_seeds{bsi};
  for sui = stimulus.surfaces
    for sti = stimulus.structures
      for pli = stimulus.partlocations
        for psi = stimulus.partstructures
          for vi = stimulus.views
            filename = sprintf('%s/%s_%i%i%i%i_view%i.png', stim_dir, bs, sui, sti, pli, psi, vi);
            [image, map, alpha] = imread(filename);
            stims.(sprintf('img%s_%i%i%i%i_v%i', bs, sui, sti, pli, psi, vi)) = genTexFromIm(cat(3, image, alpha), stimulus.live.mask);
          end
        end
      end
    end
  end
  disppercent(bsi / length(stimulus.base_seeds));
end
disppercent(inf);
stimulus.live.stims = stims;

clear stims

%%%%%%%%%%% MGL %%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialize task and begin update loop
[task{1}, myscreen] = initTask(task{1},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);

% Run the eye calibration
myscreen = eyeCalibDisp(myscreen);

% let the user know
disp(sprintf('(configural) Starting run'));

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

task.thistrial.sample2_view = randsample(repmat(setdiff(stimulus.views, task.thistrial.sample_view), 1,2), 1); % choose a different view for the match image

% Initialize the nonmatch properties to all be the same, then change one of them
task.thistrial.sample2_surface = task.thistrial.sample_surface;
task.thistrial.sample2_structure = task.thistrial.sample_structure;
task.thistrial.sample2_partlocation = task.thistrial.sample_partlocation;
task.thistrial.sample2_partstructure = task.thistrial.sample_partstructure;

if task.thistrial.same == 0
  switch task.thistrial.nonmatch_dimension
    case 1
      task.thistrial.sample2_surface = randsample(repmat(setdiff(stimulus.surfaces, task.thistrial.sample_surface), 1,2), 1);
    case 2
      task.thistrial.sample2_structure = randsample(repmat(setdiff(stimulus.structures, task.thistrial.sample_structure), 1, 2), 1);
    case 3
      task.thistrial.sample2_partlocation = randsample(repmat(setdiff(stimulus.partlocations, task.thistrial.sample_partlocation), 1,2), 1);
    case 4
      task.thistrial.sample2_partstructure = randsample(repmat(setdiff(stimulus.partstructures, task.thistrial.sample_partstructure), 1,2), 1);
  end
end

%% Get all 3 images for this trial
base_seed = stimulus.base_seeds{task.thistrial.base_seed};
stimulus.live.sample_image = stimulus.live.stims.(sprintf('img%s_%i%i%i%i_v%i',base_seed, task.thistrial.sample_surface, task.thistrial.sample_structure, task.thistrial.sample_partlocation, task.thistrial.sample_partstructure, task.thistrial.sample_view));
stimulus.live.sample2_image = stimulus.live.stims.(sprintf('img%s_%i%i%i%i_v%i', base_seed, task.thistrial.sample2_surface, task.thistrial.sample2_structure, task.thistrial.sample2_partlocation, task.thistrial.sample2_partstructure, task.thistrial.sample2_view));

stimulus.live.eyeCount = 0;

% Disp trial parameters each trial
sd = {'Different', 'Same'};
fprintf('Trial %d - %s: BaseSeed = %s. Sample = %i%i%i%i, Non-Match = %i%i%i%i', task.trialnum, sd{task.thistrial.same+1}, base_seed, stimulus.conditions{task.thistrial.nonmatch_dimension}, task.thistrial.sample_surface, task.thistrial.sample_structure,...
                                                                                 task.thistrial.sample_partlocation, task.thistrial.sample_partstructure, task.thistrial.sample2_surface, task.thistrial.sample2_structure, task.thistrial.sample2_partlocation, task.thistrial.sample2_partstructure);

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
imSz = stimulus.stimSize; % Size in Degrees of Visual Angle to display image
ecc = stimulus.eccentricity; % Eccentricity to display image at

stimulus.choice_locations = [-ecc, 0; ecc,0];

for i = 1:2
  mglClearScreen(0.5);
  if task.thistrial.thisseg == stimulus.seg.stim1
    mglBltTexture(stimulus.live.sample_image, [0,0, imSz, imSz]);
  elseif task.thistrial.thisseg == stimulus.seg.stim2
    mglBltTexture(stimulus.live.sample2_image, [0,0, imSz, imSz]);
    upFix(stimulus, stimulus.colors.black);
  elseif task.thistrial.thisseg == stimulus.seg.feedback
    if task.thistrial.response == (2-task.thistrial.same) % response 1 (left) = same,  2 (right) = different
      task.thistrial.correct = 1;
      if i == 1, fprintf('Correct! \n'); end
      upFix(stimulus, stimulus.colors.green);
    else
      task.thistrial.correct = 0;
      if i == 1, fprintf('Incorrect :( \n'); end
      upFix(stimulus, stimulus.colors.red);
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
if ~stimulus.noeye && ~any(task.thistrial.thisseg==[stimulus.seg.fix, stimulus.seg.stim2]) 
  if ~any(isnan(pos))
    if dist > 1.5 && stimulus.live.eyeCount > 20 && stimulus.fixate == 1
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

% Draws a circle at center of the screen of color fixColor
function upFix(stimulus, fixColor)
if ieNotDefined('fixColor')
  fixColor = stimulus.live.fixColor;
end
%mglGluAnnulus(0,0,0,.1,fixColor);
mglFixationCross(1,1,fixColor);


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
elseif size(r,3)==4)
  % do nothing.
else
  r(:,:,4) = mask(:,:,1);
end
% mgl requires the first dimension to be the RGBA dimension
rP = permute(r, [3 2 1]);
tex = mglCreateTexture(rP); 

%%%%%%%%%%%%%%%%%%%%%%%
%    analyzeData %
%%%%%%%%%%%%%%%%%%%%%%%
function data = analyzeData()
%%

% get the files list
files = dir(fullfile(sprintf('~/data/configural/%s/20*stim*.mat',mglGetSID)));

count = 1; 
data = struct('response', [], 'reaction_time', [], 'periphery',[], 'eccentricity', [], ...
              'nTrials', 0, 'nValTrials', 0, 'accByRuns', []);
for fi = 1:length(files)
  load(fullfile(sprintf('~/data/configural/%s/%s',mglGetSID,files(fi).name)));
  
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
    
    data.response = [data.response e.response];
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
data.imSize = stimulus.stimSize;

disp(sprintf('SUBJECT %s: Found %i runs with a total of %i trials', mglGetSID, length(data.accByRuns), data.nTrials));

