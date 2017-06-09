function [ myscreen ] = distjdg( varargin )
%
% DISTANCE JUDGMENTS 
%    Distance judgment task from Klein et al (2016). 
%
%    Usage: distjdg(varargin)
%    Authors: Akshay Jagadeesh
%    Date: 06/05/2017
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
getArgs(varargin,{'scan=0','plots=0','noeye=0','debug=0','training=0'});
stimulus.scan = scan;
stimulus.plots = plots;
stimulus.noeye = noeye;
stimulus.debug = debug;
stimulus.noimp = noimp;
stimulus.training = training;
clear localizer invisible scan noeye task test2

if stimulus.scan
    warning('Not setup for scanning');
end

%% Stimulus parameters 
stimulus.horizDist = 10; %each pair of gabors is located 10 degrees horizontally away from the center.
stimulus.spatFreq = 1;
stimulus.sd = pi/8;
stimulus.ecc = 6;

%% Open Old Stimfile
stimulus.counter = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/distjdg/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/distjdg/%s/1*mat',mglGetSID));

    if length(files) >= 1
        fname = files(end).name;
        
        s = load(sprintf('~/data/distjdg/%s/%s',mglGetSID,fname));
        % copy staircases and run numbers
        stimulus.counter = s.stimulus.counter + 1;
        stimulus.staircase = s.stimulus.staircase;
        clear s;
        disp(sprintf('(distjdg) Data file: %s loaded.',fname));
    end
end
disp(sprintf('(distjdg) This is run #%i',stimulus.counter));

%% Setup Screen
myscreen = initScreen('VPixx');

% set background to grey
myscreen.background = 0.5;

%% Staircase
if ~isfield(stimulus,'staircase')
    disp('(distjdg) WARNING: New staircase');
    initStair();
else
    resetStair();
end

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
    
stimulus.responseKeys = [1 2]; % left right

%% Colors
initGammaTable(myscreen);
stimulus.colors.rmed = 127.5;

% We're going to add an equal number of reserved colors to the top and
% bottom, to try to keep the center of the gamma table stable.
stimulus.colors.reservedBottom = [1 0 0; 0 0 0]; % fixation cross colors
stimulus.colors.reservedTop = [1 1 1; 0 1 0]; % correct/incorrect colors
stimulus.colors.black = 1/255; stimulus.colors.white = 254/255;
stimulus.colors.red = 0/255; stimulus.colors.green = 255/255;
stimulus.colors.nReserved = 2; % this is /2 the true number, because it's duplicated
stimulus.colors.nUnreserved = 256-(2*stimulus.colors.nReserved);

stimulus.colors.mrmax = stimulus.colors.nReserved - 1 + stimulus.colors.nUnreserved;
stimulus.colors.mrmin = stimulus.colors.nReserved;

%% Setup Task

%%%%%%%%%%%%% PHASE ONE %%%%%%%%%%%%%%%%%
%%%%% PRIOR + ESTIMATE OF THRESHOLD %%%%%

stimulus.curTrial(1) = 0;

task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;

task{1}{1}.seglen = [inf 0.750 0.050 0.067 0.083 1.000]; % Fixate, delay, cue, ISI, stim, response
%task{1}{1}.segmax = [inf 0.750 0.050 0.067 0.083 1.000];

stimulus.seg = {};
stimulus.seg{1}.ITI1 = 1;
stimulus.seg{1}.delay = 2;
stimulus.seg{1}.cue = 3;
stimulus.seg{1}.ISI = 4;
stimulus.seg{1}.stim = 5;
stimulus.seg{1}.resp = 6;

if stimulus.noeye==1
    task{1}{1}.seglen(1) = 0.5;
end

task{1}{1}.synchToVol = zeros(size(task{1}{1}.seglen));
task{1}{1}.getResponse = zeros(size(task{1}{1}.seglen));
task{1}{1}.getResponse(stimulus.seg{1}.resp) = 1;

% Number of trials
task{1}{1}.numTrials = 50;
task{1}{1}.random = 1;

if stimulus.scan
    task{1}{1}.synchToVol(stimulus.seg.ITI) = 1;
end

% Task trial parameters
task{1}{1}.parameter.horizDist = stimulus.horizDist;
task{1}{1}.parameter.gaborSF = stimulus.spatFreq; % Spatial frequency of Gabor
task{1}{1}.parameter.contrast = 0.45;

% Task variables to be calculated later
task{1}{1}.randVars.calculated.cuePresent = nan; % was cue shown on this trial?
task{1}{1}.randVars.calculated.cueSide = nan; % 1 for Left, 2 for Right.
task{1}{1}.randVars.calculated.leftDist = nan; % distance between the pair of Gabors on the left
task{1}{1}.randVars.calculated.rightDist = nan; % distance between the pair of Gabors on the right
task{1}{1}.randVars.calculated.rotation = nan; % rotation of the grating
%task{1}{1}.randVars.calculated.contrast = nan; % contrast of the grating
task{1}{1}.randVars.calculated.detected = 0; % did they see the grating?
task{1}{1}.randVars.calculated.dead = 0;

% Delete these all later
% task{1}{1}.randVars.calculated.startRespAngle = nan;
% task{1}{1}.randVars.calculated.respAngle = nan;
% task{1}{1}.randVars.calculated.angle = nan; % angle at which displayed, depends on attention mode
% task{1}{1}.randVars.calculated.target = nan; 
% task{1}{1}.randVars.calculated.perisaccInt = nan;

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
disp(sprintf('(distjdg) Starting run number: %i.',stimulus.counter));

%% Main Task Loop

setGammaTable(1);
mglClearScreen(0.5); mglFixationCross(1,1,stimulus.colors.white);

mglFlush
mglClearScreen(0.5); mglFixationCross(1,1,stimulus.colors.white);

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
    disp('(distjdg) Displaying plots');
    dispInfo(stimulus);
end

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startTrialCallback(task,myscreen)
%%
% if ~isempty(task.lasttrial)
%     if task.lasttrial.thisphase == 1
%         disp(sprintf('Last trial, target: %02.0f, stim: %02.0f',task.lasttrial.target*180/pi,task.lasttrial.angle*180/pi));
%     elseif task.lasttrial.thisphase == 2
%         disp(sprintf('Last trial, target: %02.0f, stim: %02.0f, resp: %02.0f',task.lasttrial.target*180/pi,task.lasttrial.angle*180/pi,task.lasttrial.respAngle*180/pi));
%     end
% end
global stimulus

task.thistrial.dead = 0;
task.thistrial.detected = 0;

if (~isempty(task.lasttrial)) && (task.lasttrial.detected~=1) && ~task.lasttrial.dead
    stimulus.staircase = doStaircase('update',stimulus.staircase,task.lasttrial.detected);
    if(task.lasttrial.detected == 0) % only print this during phase 1
        disp(sprintf('No response reported'));
    end
end

stimulus.live.gotResponse = 0;
stimulus.curTrial(task.thistrial.thisphase) = stimulus.curTrial(task.thistrial.thisphase) + 1;

% compute missing variables
task.thistrial.cueSide = ceil(rand*2); %1 for left, 2 for right.
task.thistrial.cuePresent = (rand > 0.5);

task.thistrial.standardDist = 3;
task.thistrial.otherDist = 5; %varies with staircase

if rand < 0.5
    task.thistrial.leftDist = task.thistrial.otherDist;
    task.thistrial.rightDist = task.thistrial.standardDist;
else
    task.thistrial.leftDist = task.thistrial.standardDist;
    task.thistrial.rightDist = task.thistrial.otherDist;
end

task.thistrial.rotation = rand*2*pi;
% task.thistrial.startRespAngle = rand*2*pi;
% switch stimulus.att
%     case 1 %Endo
%         task.thistrial.angle = randn*task.thistrial.priorSTD; %normally distributed around the prior
%     case 2 %Exo
%         task.thistrial.angle = rand*2*pi; % stim angle is random
%         task.thistrial.target = rand*2*pi; % cue angle is random
%         if task.thistrial.thisphase == 1 && ~stimulus.test2
%             task.thistrial.visible = (rand > 0.5);
%         end
%     case 3 %Sacc
%         task.thistrial.angle = rand*2*pi;
%         task.thistrial.target = rand*2*pi;
%     otherwise
%         disp('Invalid attention condition. Quitting...');
% end

% contrast from staircase in both phases
%[task.thistrial.contrast, stimulus.staircase] = doStaircase('testValue',stimulus.staircase);

% Reset mouse to center of screen at start of every trial
mglSetMousePosition(960,540,1);
myscreen.flushMode = 0;

disp(sprintf('(distjdg) Trial (%i): cuePresent: %d, leftDist: %d, rightDist: %d, rotation: %02.0f, contrast: %02.02f%%',...
    task.trialnum, task.thistrial.cuePresent, task.thistrial.leftDist, task.thistrial.rightDist, task.thistrial.rotation*180/pi,100*task.thistrial.contrast));
    
stimulus.live.eyeCount = 0;

%setGammaTable(task.thistrial.contrast);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)
%%

global stimulus

stimulus.live.triggerWaiting = 0;
if any(task.thistrial.thisseg==[stimulus.seg{task.thistrial.thisphase}.ITI1])
    stimulus.live.triggerWaiting = 1;
    stimulus.live.centered = 0;
    stimulus.live.triggerTime = 0;
    stimulus.live.lastTrigger = -1;
end

stimulus.live.eyeDead = 0;
stimulus.live.resp = 0;
stimulus.live.fixColor = stimulus.colors.white;
stimulus.live.fix = 1;
stimulus.live.stim = 0;
stimulus.live.cue = 0;

if task.thistrial.thisseg==stimulus.seg{task.thistrial.thisphase}.stim
    stimulus.live.stim = 1;
elseif task.thistrial.thisseg==stimulus.seg{task.thistrial.thisphase}.resp
    stimulus.live.resp = 1;
    mInfo = mglGetMouse(myscreen.screenNumber);
    stimulus.live.trackingAngle = -mInfo.x/90;
elseif task.thistrial.thisseg==stimulus.seg{task.thistrial.thisphase}.cue
    if task.thistrial.cuePresent == 1
        stimulus.live.cue = 1;
    end
end

for i = 1:2
    mglClearScreen(0.5);
    if stimulus.live.stim
        x = task.thistrial.horizDist;
        yl = task.thistrial.leftDist/2;
        yr = task.thistrial.rightDist/2;
        mglBltTexture(stimulus.live.grating,[-x yl],0,0,task.thistrial.rotation*180/pi);
        mglBltTexture(stimulus.live.grating,[-x -yl],0,0,task.thistrial.rotation*180/pi);
        mglBltTexture(stimulus.live.grating,[x yr],0,0,task.thistrial.rotation*180/pi);
        mglBltTexture(stimulus.live.grating,[x -yr],0,0,task.thistrial.rotation*180/pi);
    elseif stimulus.live.cue && i == 1
        if task.thistrial.cueSide == 1
            x = -task.thistrial.horizDist;
        else
            x = task.thistrial.horizDist;
        end
        y = 0;
        mglPolygon([x-.15, x-.15, x+.15, x+.15], [y-.15, y+.15, y+.15, y-.15], stimulus.colors.white);
%        mglFlush();
%        mglClearScreen(0.5);
    end
    
    % resp is updated in screenUpdate
    
    if stimulus.live.fix
        upFix(stimulus);
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
if ~stimulus.noeye && ~any(task.thistrial.thisseg==[stimulus.seg{task.thistrial.thisphase}.ITI1]) && ~stimulus.scan
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

if stimulus.live.resp && (task.thistrial.thisphase==2)
    mglClearScreen(0.5);
    if stimulus.live.fix, upFix(stimulus); end
    if stimulus.live.resp, mglFillOval(stimulus.live.respx,stimulus.live.respy,[0.5 0.5],stimulus.colors.white); end
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

function convertRespXY(task)
global stimulus

stimulus.live.respx = task.thistrial.ecc*cos(stimulus.live.angle+task.thistrial.target);
stimulus.live.respy = task.thistrial.ecc*sin(stimulus.live.angle+task.thistrial.target);

function upFix(stimulus)
%%
% for this experiment use a circle to indicate where participants can
% fixate inside of (rather than a cross which might arbitrarily enforce
% poisitioning
% mglGluAnnulus(0,0,1.5,1.55,stimulus.live.fixColor,64);
mglFixationCross(1,1,stimulus.live.fixColor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus

if task.thistrial.dead, return; end
    
%if isfield(task.thistrial,'whichButton') && (task.thistrial.whichButton==stimulus.responseKeys(1)) && isempty(task.thistrial.mouseButton)
    % subject didn't see anything
%    task = jumpSegment(task,inf);
%    task.thistrial.detected = -1; %-1 means they reported not seeing anything
%    disp(sprintf('Subject reported not seeing %02.02f%% contrast stimulus', task.thistrial.contrast*100));
%    return
%end

validResponse = any(task.thistrial.whichButton == stimulus.responseKeys);

if validResponse
    if stimulus.live.gotResponse==0
        task.thistrial.detected = 1;
        disp(sprintf('Subject pressed button %d', task.thistrial.whichButton));
        %stimulus.staircase = doStaircase('update',stimulus.staircase,task.thistrial.detected);
        stimulus.live.fix = 0;
    else
        disp(sprintf('Subject responded multiple times: %i',stimulus.live.gotResponse));
    end
    stimulus.live.gotResponse=stimulus.live.gotResponse+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              HELPER FUNCTIONS                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function initStair()
global stimulus

stimulus.staircase = doStaircase('init','upDown',...
            'initialThreshold',0.25,...
            'initialStepsize',0.25/3,...
            'minThreshold=0.0001','maxThreshold=0.4','stepRule','pest',...
            'nTrials=40','maxStepsize=0.2','minStepsize=0.0001');
        
function resetStair()

global stimulus

if doStaircase('stop',stimulus.staircase)
    disp('(distjdg) Staircase is being reset');
    stimulus.staircase(end+1) = doStaircase('init',stimulus.staircase(end));
    if stimulus.staircase(end).s.threshold>0.3
        disp('(distjdg) Bad staircase threshold: setting to 0.3');
        stimulus.staircase(end).s.threshold=0.3;
    elseif stimulus.staircase(end).s.threshold<0
        disp('(distjdg) Bad staircase threshold: setting to 0.05');
        stimulus.staircase(end).s.threshold=0.05;
    end
end

function [trials] = totalTrials()
%%

% Counts trials + estimates the threshold based on the last 500 trials

% get the files list
files = dir(fullfile(sprintf('~/data/distjdg/%s/17*stim*.mat',mglGetSID)));

trials = 0;

for fi = 1:length(files)
    load(fullfile(sprintf('~/data/distjdg/%s/%s',mglGetSID,files(fi).name)));
    
    e = getTaskParameters(myscreen,task);
    e = e{1}; % why?!
    trials = trials + e.nTrials;
end

%%%%%%%%%%%%%%%%%%%%%%%
%    dispInfo    %
%%%%%%%%%%%%%%%%%%%%%%%
function dispInfo(rstimulus)
%%

doStaircase('threshold',rstimulus.staircase,'dispFig=1','type=weibull');
% ctask = task; cscreen = myscreen; % save this incase we need them

% compute % correct for valid and invalid trials, display learning over
% time (including history from other runs)
% exp = getTaskParameters(task,myscreen);

% get the files list
files = dir(fullfile(sprintf('~/data/distjdg/%s/17*stim*.mat',mglGetSID)));

% load the files and pull out the data (long form)
%  rrun # counter #    local trial     real trial   angle     respAngle    
%     1       2             3              4           5           6
%  target    startRespAngle     contrast     detected      ecc    priorsd
%     7            8                9           10          11      12
%    rotation
%       13
count = 1; data = zeros(10000,13);

for fi = 1:length(files)
    load(fullfile(sprintf('~/data/distjdg/%s/%s',mglGetSID,files(fi).name)));
    
    e = getTaskParameters(myscreen,task);
    if e{1}(1).nTrials>1
        e = e{1}(2); % why?!
    
        run = stimulus.counter;

        if stimulus.att==2
            data(count:count+(e.nTrials-1),:) = [repmat(fi,e.nTrials,1) repmat(run,e.nTrials,1) (1:e.nTrials)' (count:count+(e.nTrials-1))' ...
                e.randVars.angle' e.randVars.respAngle' e.randVars.target' ...
                e.randVars.startRespAngle' e.randVars.contrast' e.randVars.detected' ...
                e.parameter.ecc' e.parameter.priorSTD' e.randVars.rotation'];
        else
            data(count:count+(e.nTrials-1),:) = [repmat(fi,e.nTrials,1) repmat(run,e.nTrials,1) (1:e.nTrials)' (count:count+(e.nTrials-1))' ...
                e.randVars.angle' e.randVars.respAngle' e.parameter.target' ...
                e.randVars.startRespAngle' e.randVars.contrast' e.randVars.detected' ...
                e.parameter.ecc' e.parameter.priorSTD' e.randVars.rotation'];
        end

        count = count+e.nTrials;
    end
end

data = data(1:(count-1),:);

if any(data(:,6)>pi), data(data(:,6)>pi,6) = data(data(:,6)>pi,6)-2*pi; end
%% Compute angle-target and respAngle-target plot
h = figure; hold on

low = [0 0.075 0.081 0.09 inf];

%for i = 1:4
    %subplot(4,1,i); hold on
    % data(:,6) = data(:,6)-data(:,5);

    % remove no-response trials
    data_ = data(~isnan(data(:,6)),:);
    %data_ = data_(logical((data_(:,9)>low(i)).*(data_(:,9)<low(i+1))),:);
    % find the trials where stimulus is - relative to the prior
    flip = data_(:,5)<0; flip = flip*1;
    flip(flip==1) = -1; flip(flip==0) = 1;
    % flip all the stimulus-target to be in the positive space
    data_(:,5:6) = data_(:,5:6).* repmat(flip,1,2);
    Y = data_(:,6);
    X = [ones(size(Y)) data_(:,5)];
    b = X\Y;
    c = [X Y];

    bci = bootci(1000,@(x) x(:,1:2)\x(:,3),c);

    % [p,s] = polyfit(data_(:,5),data_(:,6),1);


    plot(data_(:,5),data_(:,6),'*');
    % plot constant line
    plot([-1 1],[-1 1],'--r');
    % plot fit
    x = -1:1;
    plot(x,b(1)+b(2)*x,'--k');
    % compute SD of residuals? 
    % todo
    xlabel('Stimulus - Target (deg)');
    ylabel('Resp - Target (deg)');
    title(sprintf('Bias %01.2f [%01.2f %01.2f], slope %01.2f [%01.2f %01.2f], Steeper = resp away, shallower = resp toward',b(1),bci(1,1),bci(2,1),b(2),bci(1,2),bci(2,2)));

    axis([0 1 -1 1]);
    axis square

    set(gca,'XTick',0:.5:1,'XTickLabel',round((0:.5:1)*180/pi,2),'YTick',-1:.5:1,'YTickLabel',round((-1:.5:1)*180/pi,2));

    drawPublishAxis;
    % h = figure;
    % hist(data(:,5)-data(:,6));
%end

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
grating = 251/2*mglMakeGrating(sz,sz,stimulus.spatFreq,0) + 255/2;
gauss = mglMakeGaussian(sz,sz,sz/6,sz/6);
alphamask = repmat(grating,1,1,4);
alphamask(:,:,4) = gauss*255;

% we'll adjust the gamma table to control contrast
stimulus.live.grating = mglCreateTexture(alphamask); % high contrast

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sets the gamma table so that we can have
% finest possible control over the stimulus contrast.
%
% stimulus.colors.reservedColors should be set to the reserved colors (for cue colors, etc).
% maxContrast is the maximum contrast you want to be able to display.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setGammaTable(maxContrast)

global stimulus;

% set the bottom
gammaTable(1:size(stimulus.colors.reservedBottom,1),1:size(stimulus.colors.reservedBottom,2)) = stimulus.colors.reservedBottom;

% set the gamma table
if maxContrast == 1
    % create the rest of the gamma table
    cmax = 1;cmin = 0;
    luminanceVals = cmin:((cmax-cmin)/(stimulus.colors.nUnreserved-1)):cmax;

    % now get the linearized range
    redLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,luminanceVals,'linear');
    greenLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,luminanceVals,'linear');
    blueLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,luminanceVals,'linear');
elseif maxContrast > 0
    % create the rest of the gamma table
    cmax = 0.5+maxContrast/2;cmin = 0.5-maxContrast/2;
    luminanceVals = cmin:((cmax-cmin)/(stimulus.colors.nUnreserved-1)):cmax;

    % now get the linearized range
    redLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,luminanceVals,'linear');
    greenLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,luminanceVals,'linear');
    blueLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,luminanceVals,'linear');
else
    % if we are asked for 0 contrast then simply set all the values to gray
    redLinearized = repmat(interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,.5,'linear'),1,stimulus.colors.nUnreserved);
    greenLinearized = repmat(interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,.5,'linear'),1,stimulus.colors.nUnreserved);
    blueLinearized = repmat(interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,.5,'linear'),1,stimulus.colors.nUnreserved);
end

% add to the table!
gammaTable((stimulus.colors.mrmin:stimulus.colors.mrmax)+1,:)=[redLinearized;greenLinearized;blueLinearized]';

% set the top
gammaTable = [gammaTable; stimulus.colors.reservedTop];

if size(gammaTable,1)~=256
    disp('(setGammaTable) Failure: Incorrect number of colors in gamma table produced');
    keyboard
end

% set the gamma table
succ = mglSetGammaTable(gammaTable);

if ~succ
    warning('Gamma table set failure');
    keyboard
end

% remember what the current maximum contrast is that we can display
stimulus.curMaxContrast = maxContrast;


function initGammaTable(myscreen)
global stimulus
%% Gamma Table Initialization

% get gamma table
if ~isfield(myscreen,'gammaTable')
  stimulus.linearizedGammaTable = mglGetGammaTable;
  disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
  disp(sprintf('(cuecon:initGratings) No gamma table found in myscreen. Contrast displays like this'));
  disp(sprintf('         should be run with a valid calibration made by moncalib for this monitor.'));
  disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
else
  % check to make sure this matches the calibration file
  
  % get each channel table that should have been set by mglGetGammaTable
  redTable = myscreen.initScreenGammaTable.redTable(:);
  greenTable = myscreen.initScreenGammaTable.greenTable(:);
  blueTable = myscreen.initScreenGammaTable.blueTable(:);
  % get what the calibration structure says it should have been set to
  gammaTable = myscreen.gammaTable(:);
  % table values are only good to 10 bits
  redTable = round(redTable*1024)/1024;
  greenTable = round(greenTable*1024)/1024;
  blueTable = round(blueTable*1024)/1024;
  gammaTable = round(gammaTable*1024)/1024;
  % compare, ignoring nans
  if ~isequaln(mglGetGammaTable,myscreen.initScreenGammaTable) || ~isequaln(redTable,gammaTable) || ~isequaln(greenTable,gammaTable) || ~isequaln(blueTable,gammaTable)
    disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
    disp(sprintf('(curecon:initGrating) Gamma table does not match calibration'));
    disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
    keyboard
  end
end
stimulus.linearizedGammaTable = myscreen.initScreenGammaTable;
