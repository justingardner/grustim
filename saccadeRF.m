function [ myscreen ] = saccadeRF( varargin )
%
% Saccade task to measure compression 
%
%  Usage: saccadeRF(varargin)
%  Authors: Akshay Jagadeesh
%  Date: 02/27/2018
%

global stimulus

stimulus = struct;

%% Initialize Variables

% add arguments later
scan = 0;
plots = 0;
noeye = 1;
debug = 0;
getData = 0;
getArgs(varargin,{'getData=0','scan=0','plots=0','noeye=1','debug=0'});
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

%% Stimulus parameters 
%% Open Old Stimfile
stimulus.counter = 1;
if ~isempty(mglGetSID) && isdir(sprintf('~/data/saccadeRF/%s',mglGetSID))
  % Directory exists, check for a stimfile
  files = dir(sprintf('~/data/saccadeRF/%s/1*mat',mglGetSID));

  if length(files) >= 1
    fname = files(end).name;
    
    s = load(sprintf('~/data/saccadeRF/%s/%s',mglGetSID,fname));
    stimulus.counter = s.stimulus.counter + 1;
    clear s;
    disp(sprintf('(saccadeRF) Data file: %s loaded.',fname));
  end
end
disp(sprintf('(saccadeRF) This is run #%i',stimulus.counter));

%% Setup Screen
myscreen = initScreen('VPixx2');

% set background to grey
myscreen.background = 0; % change back to 0.5;

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

stimulus.curTrial(1) = 0;

% Number of Trials
nTrials = 10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% TASK 1: Left Side %%%%%%%%%%%%%%%%%%%%%%
task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;

%% Define stimulus timing
stimTime = 1.000;
ITI_min = 3.000; ITI_max = 12.00;
task{1}{1}.segmin = [stimTime ITI_min];
task{1}{1}.segmax = [stimTime ITI_max];
stimulus.seg = {};
stimulus.seg{1}.stim = 1;
stimulus.seg{1}.ITI = 2;

% Task important variables
stimulus.eccentricity = 6;
stimulus.stimSize = 2;

% Trial parameters
task{1}{1}.parameter.stimLength = [1,2,3]; % Number of frames to present stimulus for
task{1}{1}.parameter.stimSize = [1,2,4];

task{1}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
task{1}{1}.getResponse = zeros(size(task{1}{1}.segmin));
task{1}{1}.numTrials = nTrials;
task{1}{1}.random = 1;
if stimulus.scan
  task{1}{1}.synchToVol(stimulus.seg{1}.ITI) = 1;
end

% Task trial parameters

% Task variables to be calculated late
task{1}{1}.randVars.calculated.detected = 0;
task{1}{1}.randVars.calculated.dead = 0;
task{1}{1}.randVars.calculated.visible = 1;

task{1}{1}.randVars.calculated.stimTime = NaN;

%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{1})
  [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% TASK 2: Left Side %%%%%%%%%%%%%%%%%%%%%%
task{2}{1} = struct;
task{2}{1}.waitForBacktick = 1;

%% Define stimulus timing
task{2}{1}.segmin = [stimTime ITI_min];
task{2}{1}.segmax = [stimTime ITI_max];
stimulus.seg = {};
stimulus.seg{1}.stim = 1;
stimulus.seg{1}.ITI = 2;

% Trial parameters
task{2}{1}.parameter.stimLength = [1,2,3]; % Number of frames to present stimulus for
task{2}{1}.parameter.stimSize = [1,2,4];

task{2}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
task{2}{1}.getResponse = zeros(size(task{1}{1}.segmin));
task{2}{1}.numTrials = nTrials;
task{2}{1}.random = 1;
if stimulus.scan
  task{2}{1}.synchToVol(stimulus.seg{1}.ITI) = 1;
end

% Task variables to be calculated late
task{2}{1}.randVars.calculated.detected = 0;
task{2}{1}.randVars.calculated.dead = 0;
task{2}{1}.randVars.calculated.visible = 1;

task{2}{1}.randVars.calculated.stimTime = NaN;

%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{2})
  [task{2}{phaseNum}, myscreen] = initTask(task{2}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);
end

%% EYE CALIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

% let the user know
disp(sprintf('(saccadeRF) Starting run number: %i.',stimulus.counter));

%% Main Task Loop

%% change back to 0.5
mglClearScreen(0); 
upFix(stimulus);

mglFlush
mglClearScreen(0); 
upFix(stimulus);

phaseNum = 1;
% Again, only one phase.
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
  mglClearScreen(0);
  % update the left task
  [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);
  % update the right task
  [task{2}, myscreen, phaseNum] = updateTask(task{2}, myscreen, phaseNum);
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
if task.taskID == 1
  stimulus.live.stim1 = 0;
  stimulus.live.task1Dead = 0;
elseif task.taskID == 2
  stimulus.live.stim2 = 0;
  stimulus.live.task2Dead = 0;
end
if task.taskID == 1
    stimulus.live.t1 = nan;
elseif task.taskID == 2
    stimulus.live.t2 = nan;
end

stimulus.curTrial(task.thistrial.thisphase) = stimulus.curTrial(task.thistrial.thisphase) + 1;

% Disp trial parameters each trial
if task.taskID == 1
    disp(sprintf('Trial %d (L) - numFrames: %d', task.trialnum, task.thistrial.stimLength));
end

% Reset mouse to center of screen at start of every trial
mglSetMousePosition(960,540,1);
myscreen.flushMode = 0;
stimulus.live.eyeCount = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)

global stimulus

stimulus.live.text = nan;
stimulus.live.triggerWaiting = 0;
stimulus.live.eyeDead = 0;
stimulus.live.fix = 1;

if task.thistrial.thisseg == stimulus.seg{1}.ITI && task.taskID == 1
  disp(sprintf('Drawing %d frames took %g secs', task.thistrial.stimLength, stimulus.live.t1));
  stimulus.live.text = mglText(sprintf('numFrames: %d', task.thistrial.stimLength));
end

%for i = 1:2
  %if ~(stimulus.live.text && task.taskID == 2)
  %  mglClearScreen(0); %% change back
  %end
  %if stimulus.live.text && task.taskID == 1
  %  mglTextDraw(sprintf('Previous numFrames: %d', task.thistrial.stimLength), [-stimulus.eccentricity, 10]);
  %end
  
  %if stimulus.live.fix
  %  upFix(stimulus, stimulus.colors.blue);
  %end
  %mglFlush
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus

% % jump to next trial if you are dead and 1 second has elapsed since eye
% % movement
% if task.thistrial.dead && mglGetSecs(task.thistrial.segStartSeconds)>1
%   task = jumpSegment(task,inf);
% end
% 
% % skip screen updates if you are already dead
% if task.thistrial.dead
%   if task.thistrial.dead && stimulus.live.eyeDead
%     mglTextSet([],32,stimulus.colors.red);
%     mglTextDraw('Eye Movement Detected',[0 0]);
%   end
%   return
% end
% 
% % Eye movement detection Code
% if ~stimulus.noeye
%   [pos,~] = mglEyelinkGetCurrentEyePos;
%   dist = hypot(pos(1),pos(2));
%   if ~any(isnan(pos))
%     if dist > 1.5 && stimulus.live.eyeCount > 30
%       disp('Eye movement detected!!!!');
%       task.thistrial.dead = 1;
%       stimulus.live.eyeDead=1;
%       return
%     elseif dist > 1.5
%       stimulus.live.eyeCount = stimulus.live.eyeCount + 1;
%     end
%   end
% end

% Code for drawing checkerboard on each frame
if task.taskID == 1 && task.thistrial.thisseg==stimulus.seg{1}.stim && ~stimulus.live.task1Dead% LEFT side task
  % set clock
  if stimulus.live.stim1 == 0
      stimulus.live.t1 = mglGetSecs;
  end
  
  %draw stim
  if stimulus.live.stim1 < task.thistrial.stimLength
    drawCheckerboard(-stimulus.eccentricity, 0, stimulus.stimSize);
    stimulus.live.stim1 = stimulus.live.stim1 + 1;
  else
    stimulus.live.t1 = mglGetSecs(stimulus.live.t1);
    stimulus.live.task1Dead = 1;
  end
elseif task.taskID == 2 && task.thistrial.thisseg==stimulus.seg{1}.stim && ~stimulus.live.task2Dead % RIGHT side task
  % set clock
  if stimulus.live.stim2 == 0
    stimulus.live.t2 = mglGetSecs;
  end
  
  %draw stim
  if stimulus.live.stim2 < task.thistrial.stimLength
    drawCheckerboard(stimulus.eccentricity, 0, stimulus.stimSize);
    stimulus.live.stim2 = stimulus.live.stim2 + 1;
  else
    stimulus.live.t2 = mglGetSecs(stimulus.live.t2);
    stimulus.live.task2Dead = 1;
  end
end

% Code for drawing text during ITI
%if task.taskID == 1 && task.thistrial.thisseg == stimulus.seg{1}.ITI
  %mglBltTexture(stimulus.live.text, [-stimulus.eccentricity, 10]);
%end

% % Trial trigger on eye fixation code  
% if ~stimulus.noeye && stimulus.live.triggerWaiting
%   now = mglGetSecs;
%   % check eye position, if 
%   if ~any(isnan(pos))
%     wasCentered = stimulus.live.centered;
%     stimulus.live.centered = dist<2.5;
%     if wasCentered && stimulus.live.centered && stimulus.live.lastTrigger>0
%       stimulus.live.triggerTime = stimulus.live.triggerTime + now-stimulus.live.lastTrigger;
%     end
%     stimulus.live.lastTrigger = now;
%   end
%   if stimulus.live.triggerTime > 0.5 % not in ms dummy, wait 1.5 seconds (reasonable slow time)
%     disp('Starting trial--eye centered.');
%     task = jumpSegment(task);
%   end
% end

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
files = dir(fullfile(sprintf('~/data/saccadeRF/%s/18*stim*.mat',mglGetSID)));
trials = 0;

for fi = 1:length(files)
    load(fullfile(sprintf('~/data/saccadeRF/%s/%s',mglGetSID,files(fi).name)));
    e = getTaskParameters(myscreen,task);
    e = e{1}; % why?!
    trials = trials + e.nTrials;
end

%%%%%%%%%%%%%%%%%%%%%%%
%    dispInfo    %
%%%%%%%%%%%%%%%%%%%%%%%
function data = dispInfo(rstimulus)
%%

% get the files list
files = dir(fullfile(sprintf('~/data/saccadeRF/%s/18*stim*.mat',mglGetSID)));

count = 1; 
data = struct('nTrials', 0, 'subj_resp', [], 'corr_resp', [], 'corr_trials', [],...
              'image', [], 'layer', [], 'ecc', [], 'reaction_time', [], 'nValTrials', 0, 'rf_size', [], 'imSz', []);

for fi = 1:length(files)
  load(fullfile(sprintf('~/data/saccadeRF/%s/%s',mglGetSID,files(fi).name)));
  
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
    
  end
  count = count + 1;
end

return

function drawCheckerboard(x,y,sz)

sqSz = sz/2;

mglQuad([x-sqSz, x, x-sqSz, x; x, x+sqSz, x, x+sqSz; x, x+sqSz, x, x+sqSz; x-sqSz, x, x-sqSz, x],...
        [y+sqSz, y+sqSz, y, y; y+sqSz, y+sqSz, y, y; y, y, y-sqSz, y-sqSz; y, y, y-sqSz, y-sqSz],...
        [0, 1, 1, 0; 0, 1, 1, 0; 0, 1, 1, 0]);


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

alphamask = ones(20,20);
alphamask(1:10, 1:10) = 0;
alphamask(10:20, 10:20) = 0;

% we'll adjust the gamma table to control contrast
stimulus.live.grating  = mglCreateTexture(alphamask); % high contrast

