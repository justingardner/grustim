function [ myscreen ] = afcom( varargin )
% AFCOM (attention field color mapping)
% *** set 'noeye=1' to turn of the eye tracker***
%
%   This is code for a psychophysics experiment interested in finding a
%   measure that can be compared across spatial and feature-based cueing.
%   We do this by having a shared report (either color angle or motion
%   direction) that is cued by either space or the other feature (motion
%   direction, or color angle, respectively). 
%
%   The stimuli are two pairs of overlapped dot patches on the left and
%   right of fixation. The pairs share feature-cues which allows us to 
%   either cue one side of the screen or one feature, and then post-cue to
%   resolve down to one patch. Observers always report about a single patch
%   of dots, reporting either the motion direction or the color.
%
%   Control conditions: there are two types of control trials which are
%   presented 10% of the time each. The first control is a double cue, i.e.
%   we directly cue the dot patch that will need to be reported. The second
%   control is no-cue, i.e. observers have to remember all four dot
%   patches. 
%
%   We expect to see that feature cueing is not as powerful as spatial
%   cueing and that the sum of feature/spatial cueing is more than the
%   individuals. We also expect to see weird effects when you are asked to
%   remember two similar features, vs. remembering two dissimilar features.
%
%   Scanning version: in the scanner we expect to see that any cueing
%   condition causes tuning shifts, because of selection. E.g. a V4 or MT
%   response to color/motion is a post-selection response, hence there are
%   apparent tuning shifts.

% TODO:
 % - add the color cue report direction variation
 % - add control conditions
%%

global stimulus fixStimulus

stimulus = struct;
fixStimulus = struct;

%% Initialize Variables

% add arguments later
scan = 0;
plots = 0;
noeye = 0;
debug = 0;
replay = 0;
powerwheel = 0;
run = 0; 
eyewindow=0; 
mouse=0; 
practice=0; 
cue=0;

getArgs(varargin,{'scan=0','cue=1','plots=0','noeye=0','powerwheel=1','eyewindow=1.5','practice=0','debug=0','replay=0','run=0','build=0','mouse=0'});
stimulus.scan = scan;
stimulus.plots = plots;
stimulus.noeye = noeye;
stimulus.cue = cue; % cue = 1 means direction, cue = 2 means color
stimulus.practice = practice;
stimulus.mousedebug = mouse;
stimulus.powerwheel = powerwheel;
stimulus.eyewindow = eyewindow;
stimulus.debug = debug;
stimulus.replay = replay;
stimulus.overrideRun = run;

clear localizer invisible scan noeye task test2 attend build eyewindow mouse practice powerwheel cue

if stimulus.scan
    warning('Disabling eyewindow');
    stimulus.eyewindow=0;
end

%% Open Old Stimfile
if ~stimulus.replay
    stimulus.counter = 1;
    
    if ~isempty(mglGetSID) && isdir(sprintf('~/data/afcom/%s',mglGetSID))
        % Directory exists, check for a stimefile
        files = dir(sprintf('~/data/afcom/%s/1*mat',mglGetSID));
        
        if length(files) >= 1
            fname = files(end).name;
            
            s = load(sprintf('~/data/afcom/%s/%s',mglGetSID,fname));
            % copy staircases and run numbers
            stimulus.counter = s.stimulus.counter + 1;
            stimulus.colors = s.stimulus.colors;
            stimulus.colorwheel = s.stimulus.colorwheel;
            clear s;
            disp(sprintf('(afcom) Data file: %s loaded.',fname));
        else
            warning('(afcom) Unable to load previous data files. If this is *not* the first run there is something wrong.');
        end
    end
end

%% Display run info
stimulus.counter = -1;
if ~stimulus.replay
    disp('*************************');
    disp(sprintf('(afcom) This is scan #%i',stimulus.counter));
    disp('*************************');
end

%% Setup Screen
if stimulus.replay
    myscreen = initScreen('replayScreen');
else
    myscreen = initScreen('VPixx');
end
% set background to black
myscreen.background = 0;


%% Plot and return
if stimulus.plots==2
    dispInfo;
    return
end

%% Initialize Stimulus

if ~stimulus.replay
    myscreen.stimulusNames{1} = 'stimulus';
    
    if stimulus.powerwheel
        stimulus.responseKeys = 1;
    else
        stimulus.responseKeys = [1 2]; % left right
    end
else
end

%% load the calib
if isfield(myscreen,'calibFullFilename')
    calib = load(fullfile(myscreen.calibFullFilename));
    stimulus.calib = calib.calib;
else
    stimulus.calib = []; % need this so that mglLab2rgb doesn't fail
end

%% Colors
if ~isfield(stimulus,'colors')
    stimulus.colors.white = [0.8 0.8 0.8]; stimulus.colors.red = [0.8 0 0];
    stimulus.colors.green = [0 1 0]; stimulus.colors.black = [0 0 0];
end

% available range of color/direction increments
stimulus.theta_ = pi/64; % increment size
stimulus.thetas = 0:stimulus.theta_:(2*pi);
    
if 1 %~isfield(stimulus,'colorwheel')
    % get the lab space rgb values
    stimulus.backgroundLab = rgb2lab([0.5 0.5 0.5]);

    % setup color wheel
    stimulus.colorwheel.acenter = 0;%stimulus.backgroundLab(2);
    stimulus.colorwheel.bcenter = 0;%stimulus.backgroundLab(3);

    % compute the ranges around 0 and pi for the colorwheel

    D = 60;
    stimulus.colorwheel.distanceLab = D;

    stimulus.colorwheel.rgb = zeros(length(stimulus.thetas),3);
    for ti = 1:length(stimulus.thetas)
        theta = stimulus.thetas(ti);
        rgb = ang2rgb(theta);
    %     rgb = lab2rgb([stimulus.backgroundLab(1) a b]);
        stimulus.colorwheel.rgb(ti,:) = rgb;
    end

    % if any values are outside RGB space just cut them off
    stimulus.colorwheel.rgb(stimulus.colorwheel.rgb<0) = 0;
    stimulus.colorwheel.rgb(stimulus.colorwheel.rgb>1) = 1;
end

stimulus.colors.mean = [1 1 1]*mean(stimulus.colorwheel.rgb(:));

%% Sizes
stimulus.fixWidth = 0.5;
stimulus.targetWidth = 10;
stimulus.patchEcc = 8;

%% Setup patches and stencils

% there will be 12 possible locations where we can show dots. We will use 6
% patches of dots 

stimulus.dotScale = 0.3;
stimulus.cueScale = 0.1;

stimulus.dotDirs = [0.5 1.5 0.5 1.5]*pi; % when cue=1 we use these to set the dot directions
stimulus.dotColors = [0.5 1.5 0.5 1.5]*pi; % when cue=2 we use these to set the color
% horizontal directions:
% stimulus.dotDirs = [0 1 0 1]*pi;
 
dots = struct;

dots.density = 0.2;
dots.speed = 3.5;
dots.maxAlive = myscreen.framesPerSecond/4;
dots.maxX = stimulus.targetWidth;
dots.maxY = stimulus.targetWidth;

stimulus.dotThetas = [0 0 pi pi];

for di = 1:4
    stimulus.patches{di} = struct;
    
    % patch dots
    stimulus.patches{di}.dots = initDots(dots);

    % color
    if stimulus.cue==1
        stimulus.patches{di}.color = stimulus.colors.mean;
        stimulus.patches{di}.dots.dir = stimulus.dotDirs(di);
    else
        stimulus.patches{di}.color = ang2rgb(stimulus.dotColors(di));
        stimulus.patches{di}.dots.dir = 0;
    end
    
    % location
    stimulus.patches{di}.theta = stimulus.dotThetas(di);
    stimulus.patches{di}.ecc = stimulus.patchEcc;
    stimulus.patches{di}.xcenter = stimulus.patches{di}.ecc * cos(stimulus.patches{di}.theta);
    stimulus.patches{di}.ycenter = stimulus.patches{di}.ecc * sin(stimulus.patches{di}.theta);
end

% stencils
mglClearScreen(0);
mglStencilCreateBegin(1);
for di = [1 3]
    mglFillOval(stimulus.patches{di}.xcenter,stimulus.patches{di}.ycenter,[stimulus.targetWidth, stimulus.targetWidth],[1 1 1]);
end
mglFillOval(0,0,[stimulus.targetWidth stimulus.targetWidth]/4,[1 1 1]);
mglFlush;
mglStencilCreateEnd;

%% Extra stuff
stimulus.live.trackingAngle = 0;
stimulus.eyeFrames = myscreen.framesPerSecond * 0.300; % eye movements occur when for 300 ms someone moves out of the fixation region

%% Create the cue patch

dots.maxX = stimulus.targetWidth/4;
dots.maxY = stimulus.targetWidth/4;
dots.density = 2;
dots.dotScale = 3;
dots.maxAlive = 1000;
stimulus.cueDots = initDots(dots);

%% Setup Probe Task

task{1}{1} = struct;
% task waits for fixation on first segment
stimulus.seg.iti = 1;
stimulus.seg.cue = 2;
stimulus.seg.isi = 3;
stimulus.seg.stim = 4;
stimulus.seg.delay = 5;
stimulus.seg.resp = 6;
stimulus.seg.feedback = 7;

task{1}{1}.segmin = [0.5 0.5 0.5 inf 1 4 0.75];
task{1}{1}.segmax = [2.5 0.5 0.5 inf 1 4 0.75];

if stimulus.practice
    task{1}{1}.segmin = [1 2 1 4 1 inf 2];
    task{1}{1}.segmax = [2 2 1 4 1 inf 2];
end

task{1}{1}.waitForBacktick = 1;

task{1}{1}.getResponse = zeros(1,length(task{1}{1}.segmin));
task{1}{1}.getResponse(stimulus.seg.resp) = 1;

if stimulus.scan==1
    task{1}{1}.numTrials = Inf;
else
    task{1}{1}.numTrials = 40;
end

task{1}{1}.random = 1;

task{1}{1}.parameter.trialType = [1 1 1 2 2 2 0 0 3 3 4]; % 1 = spatial, 2 = feature, 0 = no cue, 3 = exact cue (1+2), 4 = target only
task{1}{1}.parameter.target = [1 2 3 4]; % which patch is the target
task{1}{1}.parameter.duration = [0.25];% 1.0]; % bump to 0.25/0.50/1.00 for full task? 
warning('Reset to 0.25/1.0'); 

task{1}{1}.parameter.cue = stimulus.cue; % which cue condition, 1=direction cues, 2=color cues

if ~stimulus.replay && stimulus.scan
    task{1}{1}.synchToVol = zeros(1,length(task{1}{1}.segmin));
    task{1}{1}.synchToVol(end) = 1;
end

% feature target
task{1}{1}.randVars.calculated.dead = nan;
task{1}{1}.randVars.calculated.targetAngle = nan; % angle of the target
task{1}{1}.randVars.calculated.distractorAngle = nan; % angle of the other thing you had to attend
task{1}{1}.randVars.calculated.distractor = nan;
task{1}{1}.randVars.calculated.featdist = nan; % number of the matched feature
task{1}{1}.randVars.calculated.sidedist = nan; % number of the matched side
task{1}{1}.randVars.calculated.distdist = nan; % number of the one you ignored (not matched side or matched feature)
task{1}{1}.randVars.calculated.angle1 = nan;
task{1}{1}.randVars.calculated.angle2 = nan;
task{1}{1}.randVars.calculated.angle3 = nan;
task{1}{1}.randVars.calculated.angle4 = nan;
task{1}{1}.randVars.calculated.respAngle = nan;
task{1}{1}.randVars.calculated.respDistance = nan;
task{1}{1}.randVars.calculated.distDistance = nan;
task{1}{1}.randVars.calculated.cwOffset = nan; % colorwheel offset rotation

%% Mouse movement storage data

% average reaction time is ~300, but the matrix will get filled with zeros
% (bad) if we don't pre-fill it with nan
stimulus.data.mouseTrack = nan(min(task{1}{1}.numTrials,50),500);
stimulus.data.mouseTick = 1;

%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,@endTrialCallback,[]);
end

%% EYE CALIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~stimulus.replay && ~stimulus.noeye
    myscreen = eyeCalibDisp(myscreen);
    
    % let the user know
    disp(sprintf('(afcom) Starting run number: %i.',stimulus.counter));
end

%% Main Task Loop

% setGammaTable(1);
mglClearScreen;
mglFlush;
mglClearScreen;
mglFlush;

phaseNum = 1;
% Again, only one phase.
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    % update the task
    [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);
    % flip screen
    myscreen = tickScreen(myscreen,task);
end

% task ended
mglClearScreen;
mglTextSet([],32,stimulus.colors.white);
% get count
mglTextDraw('Please wait',[0 0]);
mglFlush
myscreen.flushMode = 1;

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

if ~stimulus.replay && stimulus.plots
    disp('(afcom) Displaying plots');
    dispInfo(stimulus);
end

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

function dispInfo()
%%
files = dir(fullfile('~/data/afcom/',mglGetSID,'*.mat'));

maxTrackLength = 0;
for fi = 1:length(files)
    load(fullfile('~/data/afcom/',mglGetSID,files(fi).name));
    exp = getTaskParameters(myscreen,task);
    e{fi} = exp{1};
    mt{fi} = stimulus.data.mouseTrack(1:e{fi}.nTrials,:);
    mt{fi}(mt{fi}==0) = nan;
    maxTrackLength = max(maxTrackLength,size(mt{fi},2));
end

%% concatenate all trials
pvars = {'target','trialType','cue'};
rvars = {'dead','targetAngle','distractorAngle','angle1','angle2','angle3',...
    'angle4','respAngle','respDistance','distDistance'};
runs = [];

for pii = 1:length(pvars)
    eval(sprintf('%s = [];',pvars{pii}));
end
for ri = 1:length(rvars)
    eval(sprintf('%s = [];',rvars{ri}));
end

runcount = [0 0];
for run = 1:length(e)
    runs = [runs ones(1,e{run}.nTrials)];
    runcount(e{run}.parameter.cue(1)) = runcount(e{run}.parameter.cue(1)) + 1;
    for pii = 1:length(pvars)
        eval(sprintf('%s = [%s e{run}.parameter.%s];',pvars{pii},pvars{pii},pvars{pii}));
    end
    for ri = 1:length(rvars)
        eval(sprintf('%s = [%s e{run}.randVars.%s];',rvars{ri},rvars{ri},rvars{ri}));
    end
end

%% concatenate mouse tracks
amt = nan(length(target),maxTrackLength);
start = 1;
for run = 1:length(e)
    stop = (start+e{run}.nTrials-1);
    amt(start:stop,1:size(mt{run},2)) = mt{run};
    start = stop + 1;
end

%% go backward through mouseTracks and fix jumps
% assume that you end near zero, so if you jump -pi you need to -pi the
% earlier section, etc
amt = fliplr(amt);
for ai = 1:size(amt,1)
    track = amt(ai,:);
    dtrack = diff(track);
    posidx = find(dtrack>5);
    negidx = find(dtrack<-5);
    for pii = 1:length(posidx)
        idx = posidx(pii)+1;
        track(idx:end) = track(idx:end)-2*pi;
    end
    for nii = 1:length(negidx)
        idx = negidx(nii)+1;
        track(idx:end) = track(idx:end)+2*pi;
    end
    dtrack = diff(track);
    amt(ai,:) = track;
end
amt = fliplr(amt);

%% create one giant matrix, but just of a few variables that matter
data = [cue' runs' trialType' respDistance' distDistance'];
keepIdxs = ~any(isnan(data(:,4)),2);
data = data(keepIdxs,:);
amt = amt(keepIdxs,:);

%% print out information
disp(sprintf('Runs so far: %i cue direction (cue=1), %i cue color (cue=2)',runcount(1),runcount(2)));

%% plot mousetracks
% step 1: rotate mousetracks so that they are relative to the target
amt_ = amt - repmat(targetAngle(keepIdxs)',1,size(amt,2));
% test plot the average mousetrack
figure; hold on
plot(amt_','-k');
hline(0,'--r');
xlabel('Time from response window start');
ylabel('Rotation (rad)');
drawPublishAxis;

%% plot

% build one figure for each task
titles = {'Cue direction','Cue color'};
bins = pi/32:pi/16:pi;
blabels = {};
for bi = 0:(length(bins)-1)
    blabels{bi+1} = sprintf('%i/16',bi);
end

for cue = 1:2
    cdata = data(data(:,1)==cue,:);
    
    disp(sprintf('Trials of: %s so far %i',titles{cue},size(cdata,1)));
    
    all = cdata(cdata(:,3)==0,:);
    disp(sprintf('Type all: %i',size(all,1)));
    spatial = cdata(cdata(:,3)==1,:);
    disp(sprintf('Type spatial: %i',size(spatial,1)));
    feature = cdata(cdata(:,3)==2,:);
    disp(sprintf('Type feature: %i',size(feature,1)));
    target = cdata(cdata(:,3)==3,:);
    disp(sprintf('Type target: %i',size(target,1)));
    base = cdata(cdata(:,3)==4,:);
    disp(sprintf('Type baseline: %i',size(base,1)));

    figure;

    a_h = hist(all(:,4),bins);
    s_h = hist(spatial(:,4),bins);
    f_h = hist(feature(:,4),bins);
    t_h = hist(target(:,4),bins);
    b_h = hist(base(:,4),bins);
    % normalize
    a_h = a_h/sum(a_h);
    s_h = s_h/sum(s_h);
    f_h = f_h/sum(f_h);
    t_h = t_h/sum(t_h);
    b_h = b_h/sum(b_h);
    
    bar(bins,[a_h' s_h' f_h' t_h' b_h']);
    legend({'All','Spatial','Feature','Target','Baseline'});
    ylabel('Proportion (%)');
    xlabel('Response distance from target (target=0');
    set(gca,'XTick',bins,'XTickLabel',blabels);
    title(titles{cue});
    drawPublishAxis;
end

function [task, myscreen] = startTrialCallback(task,myscreen)
global stimulus

% swap seglen in
task.thistrial.seglen(stimulus.seg.stim) = task.thistrial.duration;

if stimulus.powerwheel
    mglSetMousePosition(myscreen.screenWidth/2,myscreen.screenHeight/2,1);
else
    mglSetMousePosition(myscreen.screenWidth/2,myscreen.screenHeight/2,2);
end

if stimulus.cue==1
    stimulus.cueDots.dir = stimulus.patches{task.thistrial.target}.dots.dir;
else
    stimulus.cueDots.dir = 0; % doesn't matter, dots are incoherent
end

if (task.thistrial.trialType==0) || (task.thistrial.trialType==4)
    task.thistrial.distractor = nan; % no cue, so everything is a potential distractor
else
    if task.thistrial.trialType==2
        distractors = [3 4 1 2];
    else % this could be trial type 1 or 3, but in both cases this is a spatial distractor
        distractors = [2 1 4 3];
    end
    task.thistrial.distractor = distractors(task.thistrial.target);
end

% set the angles of the patches
for di = 1:length(stimulus.patches)
    ctheta = randsample(stimulus.thetas,1);
    
    if di==task.thistrial.target
        task.thistrial.targetAngle = ctheta;
    end
    if di==task.thistrial.distractor
        task.thistrial.distractorAngle = ctheta;
    end
    
    if stimulus.cue==1
        % direction cue, so set the colors to be different
        stimulus.patches{di}.color = ang2rgb(ctheta);
    else
        % color cue, so set the directions to be different
        stimulus.patches{di}.dots.dir = ctheta;
    end
    
    task.thistrial.(sprintf('angle%i',di)) = ctheta;
end

% colorwheel random rotation
task.thistrial.cwOffset = rand*2*pi;
task.thistrial.respAngle = -task.thistrial.cwOffset;

if stimulus.cue==1
    trialTypes = {'nocue','spatial','direction','target','baseline'};
else
    trialTypes = {'nocue','spatial','color','target','baseline'};
end
disp(sprintf('(afcom) Starting trial %i. Attending %s',task.trialnum,trialTypes{task.thistrial.trialType+1}));

% eye tracking 
task.thistrial.dead = 0;
stimulus.live.eyeCount=0;

% mouse tracking
stimulus.data.mouseTick = 1;

function [task, myscreen] = endTrialCallback(task,myscreen)

if task.thistrial.dead, return; end

respType = {'timeout','click','multiclick','multiclick','multiclick'};
if isnan(task.thistrial.respDistance)
    task.thistrial.respDistance = angdist(task.thistrial.respAngle,task.thistrial.targetAngle);
    task.thistrial.distDistance = angdist(task.thistrial.respAngle,task.thistrial.distractorAngle);
    disp(sprintf('Recorded: %s. angle of %1.2f true %1.2f: %1.2f distance',respType{task.thistrial.gotResponse+1},task.thistrial.respAngle,task.thistrial.targetAngle,task.thistrial.respDistance));
end

function d = angdist(t1,t2)
d = acos(cos(t1)*cos(t2)+sin(t1)*sin(t2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = startSegmentCallback(task,myscreen)
%%
% global stimulus

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Drawing functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function drawCue(task)

global stimulus

switch task.thistrial.trialType
    case 0
        cues = 0;
    case 1
        cues = 1;
    case 2
        cues = 2;
    case 3
        cues = [1 2];
    case 4
        cues = [1 2];
end

if any(cues==0)
    % draw lines to both sides
    dotDirs = unique(stimulus.dotThetas);
    for di = 1:length(dotDirs)
        x = 1.5*stimulus.fixWidth * cos(dotDirs(di));
        y = 1.5*stimulus.fixWidth * sin(dotDirs(di));
        mglLines2(x,y,2*x,2*y,2,stimulus.colors.white);    
    end
end

if any(cues==1)
    % spatial - draw lines to attended locations
        
    % draw the line from fixWidth to 2*fixWidth
    x = 1.5*stimulus.fixWidth * cos(stimulus.patches{task.thistrial.target}.theta);
    y = 1.5*stimulus.fixWidth * sin(stimulus.patches{task.thistrial.target}.theta);
    mglLines2(x,y,2*x,2*y,2,stimulus.colors.white);
end
if any(cues==2)
    % feature - draw the motion direction or the color

    if stimulus.cue==1
        coherence = 1;
        color = stimulus.colors.white;
    elseif stimulus.cue==2
        coherence = 0;
        color = ang2rgb(stimulus.dotColors(task.thistrial.target));
    end
    % cue dots version
    stimulus.cueDots = updateDots(stimulus.cueDots,coherence,false);
    
    mglStencilSelect(1);
    afPoints(stimulus.cueDots.x-stimulus.cueDots.maxX/2,stimulus.cueDots.y-stimulus.cueDots.maxY/2,stimulus.cueScale,color);
    mglStencilSelect(0);
end

function drawStim(task,stimSeg)

global stimulus

mglStencilSelect(1);
% update and collapse x/y coordinates for drawing
n = stimulus.patches{1}.dots.n;

x = nan(1,n*length(stimulus.patches));
y = x;
r = stimulus.colors.mean(1)*ones(1,n*length(stimulus.patches));
g = r;
b = r;

for di = 1:length(stimulus.patches)
    if task.thistrial.trialType~=4 || ~stimSeg || (stimSeg && (di==task.thistrial.target))
        if stimulus.cue==1 || (stimulus.cue==2 && stimSeg)
            % if this is the actual stim seg and using motion, update
            % coherently
            stimulus.patches{di}.dots = updateDots(stimulus.patches{di}.dots,1,false);
        else
            % otherwise use incoherent motion
            stimulus.patches{di}.dots = updateDots(stimulus.patches{di}.dots,0,false);
        end

        offX = stimulus.patches{di}.xcenter - stimulus.patches{di}.dots.maxX/2;
        offY = stimulus.patches{di}.ycenter - stimulus.patches{di}.dots.maxY/2;

        x(((di-1)*n+1):(di*n)) = offX + stimulus.patches{di}.dots.x;
        y(((di-1)*n+1):(di*n)) = offY + stimulus.patches{di}.dots.y;
        if stimulus.cue==2 || (stimulus.cue==1 && stimSeg) % if this is the actual stimulus segment and we are using color, show the colors
            r(((di-1)*n+1):(di*n)) = stimulus.patches{di}.color(1);
            g(((di-1)*n+1):(di*n)) = stimulus.patches{di}.color(2);
            b(((di-1)*n+1):(di*n)) = stimulus.patches{di}.color(3);
        end
    end
end

drop = isnan(x);
if any(drop)
    x = x(~drop);
    y = y(~drop);
    r = r(~drop);
    g = g(~drop);
    b = b(~drop);
end

% randomly sort x/y/r/g/b so that overlapping patches render correctly
perm = randperm(length(x));
x = x(perm); y = y(perm); r = r(perm); g = g(perm); b = b(perm);

afPoints(x,y,stimulus.dotScale,[r' g' b']);

% draw all the dots at once
% mglPoints2c(x,y,stimulus.dotScale,r,g,b);

mglStencilSelect(0);

function afPoints(x,y,scale,c)

cFlag = size(c,1)==1;
% draw the dots one at a time with mglGluDisk
for di = 1:length(x)
    if cFlag
        mglGluDisk(x(di),y(di),scale,c);
    else
        mglGluDisk(x(di),y(di),scale,c(di,:));
    end
end

function drawTarget(task)

global stimulus

if stimulus.cue==1
    stimulus.patches{task.thistrial.target}.dots = updateDots(stimulus.patches{task.thistrial.target}.dots,1,false);
    color = stimulus.colors.mean;
else
    % if we we cued color set the coherence to zero so that there's no
    % direction information
    stimulus.patches{task.thistrial.target}.dots = updateDots(stimulus.patches{task.thistrial.target}.dots,0,false);
    color = ang2rgb(stimulus.dotColors(task.thistrial.target));
end
% compute the offset position
offX = stimulus.patches{task.thistrial.target}.xcenter - stimulus.patches{task.thistrial.target}.dots.maxX/2;
offY = stimulus.patches{task.thistrial.target}.ycenter - stimulus.patches{task.thistrial.target}.dots.maxY/2;

% draw the actual points
mglStencilSelect(1);
afPoints(stimulus.patches{task.thistrial.target}.dots.x + offX,stimulus.patches{task.thistrial.target}.dots.y + offY,stimulus.dotScale,color);
mglStencilSelect(0);

function drawPicker(task)

global stimulus

if stimulus.cue==1
    % When we cue spatial/direction we need to draw the color picker
    for ti = 1:length(stimulus.thetas)
        theta = stimulus.thetas(ti) + task.thistrial.cwOffset;
        mglGluPartialDisk_(0,0,1,1.25,180/pi*(theta-stimulus.theta_/2),180/pi*stimulus.theta_,stimulus.colorwheel.rgb(ti,:));
    end
    % Also draw a little marker to indicate the current rotation
    mglGluPartialDisk_(0,0,1,1.25,180/pi*(task.thistrial.respAngle+task.thistrial.cwOffset)-2.5,5,[0.75 0.75 0.75]);
else
    % Don't rotate the marker using cwOffset
    mglGluPartialDisk_(0,0,1,1.25,180/pi*(task.thistrial.respAngle)-2.5,5,[0.75 0.75 0.75]);
end

function drawResp(angle)

global stimulus
% Draw the chosen color as a backgroundc ircle
if stimulus.cue==1
    mglFillOval(0,0,stimulus.fixWidth*[1 1],ang2rgb(angle));
else
    mglGluPartialDisk_(0,0,1,1.25,180/pi*angle-2.5,5,[0.75 0.75 0.75]);
end

function mglGluPartialDisk_(x,y,isize,osize,sangle,sweep,color)
% just a wrapper around mglGluPartialDisk which converst from REAL angles
% to MGL angles. I absolutely hate this aspect of MGL which I assume is
% inherited from OpenGL...
sangle = 90-sangle; % this sets 0 to be vertical and all coordinates go clockwise
mglGluPartialDisk(x,y,isize,osize,sangle,sweep,color);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus

mglClearScreen();

if task.thistrial.dead && mglGetSecs(stimulus.live.deadTime)>1
    task = jumpSegment(task,inf);
end

% skip screen updates if you are already dead
if task.thistrial.dead
    if task.thistrial.dead
        mglTextSet([],32,stimulus.colors.red);
        mglTextDraw('Eye Movement Detected',[0 0]);
    end
    return
end

if (task.thistrial.thisseg==stimulus.seg.resp)
    if stimulus.powerwheel
        mInfo = mglGetMouse(myscreen.screenNumber);
        curPos = -mInfo.x/90;
        task.thistrial.respAngle = task.thistrial.respAngle + curPos-stimulus.live.trackingAngle;
        stimulus.live.trackingAngle = curPos;
    else
        mInfo = mglGetMouse(myscreen.screenNumber);
        degx = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
        degy = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;
        if stimulus.cue==1
            task.thistrial.respAngle = atan2(degy,degx) - task.thistrial.cwOffset;
        else
            task.thistrial.respAngle = atan2(degy,degx);
        end
    end
    task.thistrial.respAngle = mod(task.thistrial.respAngle,2*pi);
    
    stimulus.data.mouseTrack(task.trialnum,stimulus.data.mouseTick) = task.thistrial.respAngle;
    stimulus.data.mouseTick = stimulus.data.mouseTick + 1;
    
    % note that respAngle is stored in *real* angles -- so that it
    % corresponds correctly to the direction task. This means that when you
    % transform into visual space you need to flip into MGL angles, see
    % mglGluDiskAnnulus_ which does this step
end

switch task.thistrial.thisseg
        
    case stimulus.seg.iti
        drawStim(task,false);
        drawFix(task,stimulus.colors.white);
        
    case stimulus.seg.cue
        % fixation
        drawStim(task,false);
        drawCue(task);
        drawFix(task,stimulus.colors.white);
        
    case stimulus.seg.isi
        drawStim(task,false);
        drawFix(task,stimulus.colors.white);
        
    case stimulus.seg.stim
        drawStim(task,true);
        drawFix(task,stimulus.colors.white);
        
    case stimulus.seg.delay
        drawFix(task,stimulus.colors.white);
        
    case stimulus.seg.resp
        drawTarget(task);
        drawPicker(task);
        if stimulus.cue==1
            % only draw the chosen color at fixation if we're doing cued
            % direction
            drawResp(task.thistrial.respAngle);
        end
        drawFix(task,stimulus.colors.white);
        
    case stimulus.seg.feedback
        drawResp(task.thistrial.targetAngle);
        drawFix(task,stimulus.colors.white);
        
end

drawAllBorders(stimulus.patches,stimulus.targetWidth/2);

% do eye position tracking, but only during some segments
if any(task.thistrial.thisseg==[stimulus.seg.cue stimulus.seg.stim stimulus.seg.resp])
    % check eye pos
    if (~stimulus.noeye) && (stimulus.eyewindow>0)

        % mouse version for testing with no eyetracker
        if stimulus.mousedebug
            mInfo = mglGetMouse(myscreen.screenNumber);
            degx = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
            degy = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;

            pos = [degx, degy];
        else
            [pos,~] = mglEyelinkGetCurrentEyePos;
        end
        % compute distance
        dist = hypot(pos(1),pos(2));
    end

    % Eye movement detection code
    if (~stimulus.noeye) && (stimulus.eyewindow>0) && ~task.thistrial.dead
        if ~any(isnan(pos))
            if dist > stimulus.eyewindow && stimulus.live.eyeCount > stimulus.eyeFrames
                disp('Eye movement detected!!!!');
                stimulus.live.deadTime = mglGetSecs;
                task.thistrial.dead = 1;
                return
            elseif dist > stimulus.eyewindow
                stimulus.live.eyeCount = stimulus.live.eyeCount + 1;
            end
        end
    end
end

function [task, myscreen] = getResponseCallback(task, myscreen)

if task.thistrial.dead, return; end

if task.thistrial.gotResponse==0
    % jump to the feedback segment
    task = jumpSegment(task);
end

function drawFix(task,color)

global stimulus;

if task.thistrial.dead
    mglGluDisk(0,0,[1 1],stimulus.colors.red,60);
else
    mglFixationCross(stimulus.fixWidth,1,color);
end

function drawBorder(x,y,r,c)

for t = 0:pi/4:2*pi
    mglGluPartialDisk(x,y,r,r+0.05,(180/pi)*(t-pi/16),22.5,c);
end

function drawAllBorders(locations,r)

% draw the borders
for li = 1:length(locations)
    drawBorder(locations{li}.xcenter,locations{li}.ycenter,r,[0.05 0.05 0.05]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helpers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rgb = ang2rgb(ang)

global stimulus

a = stimulus.colorwheel.distanceLab*cos(ang)+stimulus.colorwheel.acenter;
b = stimulus.colorwheel.distanceLab*sin(ang)+stimulus.colorwheel.bcenter;

rgb = mglLab2rgb([stimulus.backgroundLab(1) a b],stimulus.calib);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create dots for horizontal motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initDots(dots)

dots.dir = 0;

area = dots.maxX*dots.maxY;

dots.n = round(area * dots.density);

% make a some points
% dots.n = 500*dots.density;
% make sure it's an even number
dots.n = dots.n + mod(dots.n,2);

% set half to white and half to black
dots.con = repmat([1 2],1,dots.n/2);

dots.x = rand(1,dots.n)*dots.maxX;
dots.y = rand(1,dots.n)*dots.maxY;

% Why replace dots? Because if you don't then peripheral overlapping dot
% patches will rival!!

dots.alive = randi(dots.maxAlive,1,dots.n); % set to random up to 200 ms

dots.xdisp = dots.x;
dots.ydisp = dots.y;

% set incoherent dots to 0
dots.coherency = 1;
dots.incoherent = rand(1,dots.n) > dots.coherency;
dots.incoherentn = sum(dots.incoherent);
dots.coherent = ~dots.incoherent;

% track time
dots.time = mglGetSecs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots for horizontal motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updateDots(dots,coherence,repick)

elapsed = mglGetSecs-dots.time;
dots.time = mglGetSecs;

dots.alive = dots.alive+1;
rIdx = dots.alive>dots.maxAlive;
replace = sum(rIdx);
dots.x(rIdx) = rand(1,replace)*dots.maxX;
dots.y(rIdx) = rand(1,replace)*dots.maxY;
dots.alive(rIdx) = 0;

% get the coherent and incoherent dots
if repick
    dots.incoherent = rand(1,dots.n) > coherence;
    dots.incoherentn = sum(dots.incoherent);
    dots.coherent = ~dots.incoherent;
    dots.coherency = coherence;
elseif dots.coherency ~= coherence
    cohDiff = coherence - dots.coherency;
    numDots = round(abs(cohDiff) * dots.n); % actual number of dots to flip
    if numDots > dots.n, numDots = dots.n; end
    if cohDiff > 0
        % we need to add more coherent dots
        flipDots = [zeros(1,numDots) ones(1,sum(dots.incoherent)-numDots)];
        dots.incoherent(dots.incoherent) = flipDots(randperm(length(flipDots)));
    else
        % we need to add more incoherent dots
        flipDots = [ones(1,numDots) zeros(1,sum(dots.coherent)-numDots)];
        dots.incoherent(dots.coherent) = flipDots(randperm(length(flipDots)));
    end
    dots.incoherentn = sum(dots.incoherent);
    dots.coherent = ~dots.incoherent;
    dots.coherency = sum(dots.coherent)/dots.n;
end
dots.coherentn = dots.n-dots.incoherentn;

vectorLength = dots.speed*elapsed;

% move coherent dots
dots.x(dots.coherent) = dots.x(dots.coherent) + vectorLength * cos(dots.dir);
dots.y(dots.coherent) = dots.y(dots.coherent) + vectorLength * sin(dots.dir);

dots.x(dots.incoherent) = dots.x(dots.incoherent) + vectorLength * cos(rand(1,dots.incoherentn)*2*pi);
dots.y(dots.incoherent) = dots.y(dots.incoherent) + vectorLength * sin(rand(1,dots.incoherentn)*2*pi);

offscreen = dots.x > dots.maxX;
dots.x(offscreen) = mod(dots.x(offscreen),dots.maxX);
offscreen = dots.x < 0;
dots.x(offscreen) = mod(dots.x(offscreen),dots.maxX);

offscreen = dots.y > dots.maxY;
dots.y(offscreen) = mod(dots.y(offscreen),dots.maxY);
offscreen = dots.y < 0;
dots.y(offscreen) = mod(dots.y(offscreen),dots.maxY);
