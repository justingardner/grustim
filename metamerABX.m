function [ myscreen ] = metamerABX( varargin )
%
% METAMER DISCRIMINATION 
%  Metamer discrimination task 
%
%  Usage: metamerABX(varargin)
%  Authors: Akshay Jagadeesh
%  Date: 10/31/2017
%

global stimulus

stimulus = struct;

%% Initialize Variables

% add arguments later
scan = 0;
plots = 0;
noeye = 0;
debug = 0;
noimp = 0;
training = 0;
test2 = 0;
att = 0;
mod = 0;
long = 1;
getArgs(varargin,{'long=1', 'att=0', 'mod=0', 'scan=0','plots=0','noeye=0','debug=0','training=0','test2=0'});
stimulus.mod = mod;
stimulus.att = att;
stimulus.scan = scan;
stimulus.plots = plots;
stimulus.test2 = test2;
stimulus.noeye = noeye;
stimulus.debug = debug;
stimulus.noimp = noimp;
stimulus.training = training;
stimulus.long = long;
clear localizer invisible scan noeye task test2


if stimulus.plots
    dispInfo(stimulus);
end

if stimulus.scan
  warning('Not setup for scanning');
end

%% Stimulus parameters 
%% Open Old Stimfile
stimulus.counter = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/metamerABX/%s',mglGetSID))
  % Directory exists, check for a stimfile
  files = dir(sprintf('~/data/metamerABX/%s/1*mat',mglGetSID));

  if length(files) >= 1
    fname = files(end).name;
    
    s = load(sprintf('~/data/metamerABX/%s/%s',mglGetSID,fname));
    stimulus.counter = s.stimulus.counter + 1;
    clear s;
    disp(sprintf('(metamerABX) Data file: %s loaded.',fname));
  end
end
disp(sprintf('(metamerABX) This is run #%i',stimulus.counter));

%% Setup Screen
myscreen = initScreen('VPixx2');

% set background to grey
myscreen.background = 0.5;

%% Setup missing initial variables

if ~isfield(stimulus,'counter')
  stimulus.counter = 1; % This keeps track of what "run" we are on.
end

%% Plot and return
if stimulus.plots==2
  dispInfo(stimulus);
  return
end

%% Initialize Stimulus

myscreen = initStimulus('stimulus',myscreen);

localInitStimulus();
  
% Set response keys
stimulus.responseKeys = [1 2]; % 

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

if stimulus.long
  st = 0.400;
else
  st = 0.200;
end

if ~stimulus.att
  % task waits for fixation on first segment
  task{1}{1}.segmin = [inf st .500 st 1.00 st 2.00 .200];
  task{1}{1}.segmax = [inf st .500 st 1.00 st 2.00 .200];

  stimulus.seg = {};
  stimulus.seg{1}.fix = 1;
  stimulus.seg{1}.stim1 = 2;
  stimulus.seg{1}.ISI1 = 3;
  stimulus.seg{1}.stim2 = 4;
  stimulus.seg{1}.ISI2 = 5;
  stimulus.seg{1}.stim3 = 6;
  stimulus.seg{1}.resp = 7;
  stimulus.seg{1}.feedback = 8;
else
  task{1}{1}.segmin = [inf .200 .200 .500 .200 .200 1.00 .200 .200 2.00];
  task{1}{1}.segmax = [inf .200 .200 .500 .200 .200 1.00 .200 .200 2.00];
  stimulus.seg = {};
  stimulus.seg{1}.fix = 1;
  stimulus.seg{1}.cue1 = 2;
  stimulus.seg{1}.stim1 = 3;
  stimulus.seg{1}.ISI1 = 4;
  stimulus.seg{1}.cue2 = 5;
  stimulus.seg{1}.stim2 = 6;
  stimulus.seg{1}.ISI2 = 7;
  stimulus.seg{1}.cue3 = 8;
  stimulus.seg{1}.stim3 = 9;
  stimulus.seg{1}.resp = 10;
  stimulus.seg{1}.cue = [2 5 8];
end

if stimulus.noeye==1
  task{1}{1}.segmin(1) = 0.5;
  task{1}{1}.segmax(1) = 0.5;
end

% Task trial parameters
task{1}{1}.parameter.image = [0 2 3 5];
task{1}{1}.parameter.scaling = [3 4 5 6 7 10];
task{1}{1}.parameter.synthPair = [1 2 3]; % 1: (1,1); 2:(1,2); 3:(2,3)
task{1}{1}.parameter.stimOrder = [1 2];
task{1}{1}.parameter.correctResponse = [1 2];

task{1}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
task{1}{1}.getResponse = zeros(size(task{1}{1}.segmin));
task{1}{1}.getResponse(stimulus.seg{1}.resp)=1;
task{1}{1}.numTrials = 288;
task{1}{1}.random = 1;

if stimulus.scan
  task{1}{1}.synchToVol(stimulus.seg.ITI) = 1;
end

% Task trial parameters

% Task variables to be calculated later
if stimulus.att
  task{1}{1}.randVars.calculated.cuePeriphery = 0; % is the cue central or in the periphery.
  task{1}{1}.randVars.calculated.cueAngle = nan;
end
if stimulus.mod
  task{1}{1}.randVars.calculated.mod1 = nan;
  task{1}{1}.randVars.calculated.mod2 = nan;
end
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
disp(sprintf('(metamerABX) Starting run number: %i.',stimulus.counter));

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

if stimulus.plots
  disp('(metamerABX) Displaying plots');
  dispInfo(stimulus);
end

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startTrialCallback(task,myscreen)
imDir = '~/proj/v2model/output';
%%
if ~isempty(task.lasttrial)
  disp(sprintf('Last trial - image: %d, scaling: 0.%d, correct response: %d, subject''s response: %d', task.lasttrial.image, task.lasttrial.scaling, task.lasttrial.correctResponse, task.lasttrial.response));
end
global stimulus

task.thistrial.dead = 0;
task.thistrial.detected = 0;
task.thistrial.visible = 1;
task.thistrial.response = 0;

stimulus.live.gotResponse = 0;
stimulus.curTrial(task.thistrial.thisphase) = stimulus.curTrial(task.thistrial.thisphase) + 1;

%% Load two images for this trial
switch task.thistrial.synthPair
  case 1
    synth1 = 1; synth2 = 2;
  case 2
    synth1 = 1; synth2 = 3;
  case 3
    synth1 = 2; synth2 = 3;
end
% Load images according to synthpair and stimOrder
if stimulus.mod
  task.thistrial.mod1 = randi(2);
  task.thistrial.mod2 = randi(2);
  imDir = '~/proj/v2model/modified';
  if task.thistrial.stimOrder == 1
    im1 = imread(sprintf('%s/im%d_s%d_synth%d_m1.png', imDir, task.thistrial.image, task.thistrial.scaling, synth1));
    im2 = imread(sprintf('%s/im%d_s%d_synth%d_m1.png', imDir, task.thistrial.image, task.thistrial.scaling, synth2));
  else
    im1 = imread(sprintf('%s/im%d_s%d_synth%d_m1.png', imDir, task.thistrial.image, task.thistrial.scaling, synth2));
    im2 = imread(sprintf('%s/im%d_s%d_synth%d_m1.png', imDir, task.thistrial.image, task.thistrial.scaling, synth1));
  end
else
  imDir = '~/proj/v2model/output';
  if task.thistrial.stimOrder == 1
    im1 = imread(sprintf('%s/im%d_s%d_synth%d.png', imDir, task.thistrial.image, task.thistrial.scaling, synth1));
    im2 = imread(sprintf('%s/im%d_s%d_synth%d.png', imDir, task.thistrial.image, task.thistrial.scaling, synth2));
  else
    im1 = imread(sprintf('%s/im%d_s%d_synth%d.png', imDir, task.thistrial.image, task.thistrial.scaling, synth2));
    im2 = imread(sprintf('%s/im%d_s%d_synth%d.png', imDir, task.thistrial.image, task.thistrial.scaling, synth1));
  end
end
% Use images to generate textures
tex1 = generateTextureFromImage(im1);
tex2 = generateTextureFromImage(im2);
clear('im1', 'im2');
stimulus.live.tex1 = tex1;
stimulus.live.tex2 = tex2;

% Set tex3 according to task.thistrial.correctResponse
if task.thistrial.correctResponse == 1
  stimulus.live.tex3 = tex1;
else
  stimulus.live.tex3 = tex2;
end

% set response text
stimulus.live.responseText = mglText('1 or 2?');


%% For attention task, determine cue location
if stimulus.att
  if rand > 0.5
    task.thistrial.cuePeriphery = 1;
  end
  task.thistrial.cueAngle = rand*2*pi;
end

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
stimulus.live.resp = 0;
stimulus.live.fix = 1;
stimulus.live.stim1 = 0;
stimulus.live.stim2 = 0;
stimulus.live.stim3 = 0;
stimulus.live.cue = 0;
stimulus.live.feedback = 0;

if task.thistrial.thisseg==stimulus.seg{task.thistrial.thisphase}.stim1
  stimulus.live.stim1 = 1;
elseif task.thistrial.thisseg==stimulus.seg{task.thistrial.thisphase}.stim2
  stimulus.live.stim2 = 1;
elseif task.thistrial.thisseg==stimulus.seg{task.thistrial.thisphase}.stim3
  stimulus.live.stim3 = 1;
elseif task.thistrial.thisseg==stimulus.seg{task.thistrial.thisphase}.resp
  stimulus.live.resp = 1;
elseif task.thistrial.thisseg==stimulus.seg{task.thistrial.thisphase}.feedback
  stimulus.live.fix = 0;
  stimulus.live.feedback = 1;
elseif stimulus.att && any(task.thistrial.thisseg==stimulus.seg{task.thistrial.thisphase}.cue)
  stimulus.live.cue = 1;
end

for i = 1:2
  mglClearScreen(0.5);
  if stimulus.live.stim1 
    mglBltTexture(stimulus.live.tex1,[0 0 26 26]);
  elseif stimulus.live.stim2
    mglBltTexture(stimulus.live.tex2,[0 0 26 26]);
  elseif stimulus.live.stim3
    mglBltTexture(stimulus.live.tex3, [0 0 26 26]);
  elseif stimulus.live.cue
    if ~task.thistrial.cuePeriphery, ecc = 0;
    else ecc = 3.5; end
    x = ecc * cos(task.thistrial.cueAngle);
    y = ecc * sin(task.thistrial.cueAngle);
    drawCue(x,y, stimulus);
  elseif stimulus.live.resp
    mglDeleteTexture(stimulus.live.tex1);
    mglDeleteTexture(stimulus.live.tex2);
  end
  
  if stimulus.live.fix
    upFix(stimulus);
  elseif stimulus.live.feedback
    if task.thistrial.response == task.thistrial.correctResponse
      upFix(stimulus, stimulus.colors.green);
    else
      upFix(stimulus, stimulus.colors.red);
    end
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
    task.thistrial.response = task.thistrial.whichButton;
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
mglGluAnnulus(0,0,0,.2,fixColor);


%%% 
% Draws a circular cue at location x,y
function drawCue(x,y, stimulus)
mglGluAnnulus(x,y, 0.75, 0.8, stimulus.live.cueColor, 64);

%%%
% Turns image into a texture
function tex = generateTextureFromImage(im)
r = flipud(repmat(im, 1, 1, 3));
r(:,:,4) = 255;
rP = permute(r, [3 2 1]);
tex = mglCreateTexture(rP);  

function [trials] = totalTrials()
%%

% Counts trials + estimates the threshold based on the last 500 trials

% get the files list
files = dir(fullfile(sprintf('~/data/metamerABX/%s/17*stim*.mat',mglGetSID)));

trials = 0;

for fi = 1:length(files)
    load(fullfile(sprintf('~/data/metamerABX/%s/%s',mglGetSID,files(fi).name)));
    
    e = getTaskParameters(myscreen,task);
    e = e{1}; % why?!
    trials = trials + e.nTrials;
end

%%%%%%%%%%%%%%%%%%%%%%%
%    dispInfo    %
%%%%%%%%%%%%%%%%%%%%%%%
function dispInfo(rstimulus)
%%

% ctask = task; cscreen = myscreen; % save this incase we need them

% compute % correct for valid and invalid trials, display learning over
% time (including history from other runs)
% exp = getTaskParameters(task,myscreen);

% get the files list
files = dir(fullfile(sprintf('~/data/metamerABX/%s/17*stim*.mat',mglGetSID)));

% load the files and pull out the data (long form)
%  rrun # counter #    local trial     real trial   angle     respAngle    
%     1       2             3              4           5           6
%  target    startRespAngle     contrast     detected      ecc    priorsd
%     7            8                9           10          11      12
%    rotation
%       13
count = 1; data = zeros(10000, 5,6);

ims = [0, 0, 2, 3, 5];

for fi = 1:length(files)
    load(fullfile(sprintf('~/data/metamerABX/%s/%s',mglGetSID,files(fi).name)));
    
    e = getTaskParameters(myscreen,task);
    if e{1}.nTrials>1
    
        run = stimulus.counter;
        corrects = e{1}.response == e{1}.parameter.correctResponse;
        scaling = e{1}.parameter.scaling;
        image = e{1}.parameter.image;
        data(count,1, :) = [sum(corrects(scaling==3)), sum(corrects(scaling==4)), sum(corrects(scaling==5)), sum(corrects(scaling==6)), sum(corrects(scaling==7)), sum(corrects(scaling==10))]/48;
        for imi = 2:5
            data(count,imi, :) = [sum(corrects(scaling==3 & image == ims(imi))), sum(corrects(scaling==4 & image == ims(imi))), sum(corrects(scaling==5 & image == ims(imi))), sum(corrects(scaling==6 & image == ims(imi))), sum(corrects(scaling==7 & image == ims(imi))), sum(corrects(scaling==10 & image == ims(imi)))] / 12;
        end
    end
    count = count + 1;
end

data = data(1:(count-1),:,:);

mglClose;

%% fit parameters
respAcc = squeeze(mean(data(:,1,:),1))';
x0 = [1, 0.5];
options = optimset('Display','iter','PlotFcns',@optimplotfval);
[x,fval,exitflag,output] = fminsearch(@(x) objectivefcn(x,respAcc), x0, options);

scales = [0.3, 0.4, 0.5, 0.6, 0.7, 1.0];
for si = 1:length(scales)
  pc(si) = probCorrect(scales(si), x(1), x(2));
end

figure;
x1 = [0.3 0.4 0.5 0.6 0.7 1.0];
y = squeeze(mean(data(:,1,:),1));
plot(x1, squeeze(mean(data(:,1,:),1)), '*-k', 'LineWidth', 4, 'MarkerSize', 5); hold on;
plot(x1, squeeze(mean(data(:,2,:),1)), '*-y', 'LineWidth', 0.5, 'MarkerSize',1); hold on;
plot(x1, squeeze(mean(data(:,3,:),1)), '*-g', 'LineWidth', 0.5, 'MarkerSize',1); hold on;
plot(x1, squeeze(mean(data(:,4,:),1)), '*-b', 'LineWidth', 0.5, 'MarkerSize',1); hold on;
plot(x1, squeeze(mean(data(:,5,:),1)), '*-m', 'LineWidth', 0.5, 'MarkerSize',1); hold on;

plot([0.3 0.4 0.5 0.6 0.7 1.0], pc, '*-r', 'LineWidth', 4, 'MarkerSize', 5); hold on;
errorbar(x1,y, 1.96*squeeze(std(data(:,1,:),1))/sqrt(length(y)), 'k');

legend('All', 'Fountain', 'Railroad', 'City', 'Waterfall', 'BestFit');
title(sprintf('Scaling vs Accuracy - Gain: %g; Scaling: %g', x(1), x(2)), 'FontSize', 18);
xlabel('Scaling Constant', 'FontSize', 14);
ylabel('Accuracy', 'FontSize', 14);

keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions for analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pc = probCorrect(s, a0, s0)
pc = normcdf(d(s,a0,s0) / sqrt(2)) * normcdf(d(s,a0,s0) / 2) + normcdf(-d(s,a0,s0) / sqrt(2))*normcdf(-d(s,a0,s0) / 2);

function dist = d(s, a0, s0)
if s > s0
  dist = a0 * (1 - (s0^2)/(s^2));
else
  dist = 0;
end

function f = objectivefcn(x, data)

%data = [0.50 0.604166666666667 0.651041666666667 0.625 0.59375 0.723958333333333];
%data = [.48 .51 .51 .6 .7 .8];

a0 = x(1); s0 = x(2);
scales = [0.3 0.4 0.5 0.6 0.7 1.0];
for i = 1:6
  probs(i) = probCorrect(scales(i), a0, s0);
end
f = sum((probs - data).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localInitStimulus()

global stimulus

% stimulus.live.pos = 
% stimulus.live.pos = logspace(log10(stimulus.cur_.isize),log10(stimulus.cur_.osize),stimulus.cur_.N);

% mglClearScreen(0.5)
% for gi = 1:stimulus.cur_.N
    % for each grating distance
    % calculate the center position to estimate the radius
%     crad = stimulus.live.pos(gi);
    % get total degrees around circle
%     degs = 2*pi*crad;
%     sz = degs/stimulus.cur_.num*0.6;
sz = 1.5;
% use total degs / num to compute size
grating = 251/2*mglMakeGrating(sz,sz,2,0) + 255/2;
gauss = mglMakeGaussian(sz,sz,sz/6,sz/6);
alphamask = repmat(grating,1,1,4);
alphamask(:,:,4) = gauss*255;

% we'll adjust the gamma table to control contrast
stimulus.live.grating  = mglCreateTexture(alphamask); % high contrast        

