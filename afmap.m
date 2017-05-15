function [ myscreen ] = afmap( varargin )
%ATTENTIONFIELDMAPPING 
%
%   Map the attention field in the scanner. This function works by having a
%   participant perform an asynchronous attention task at fixation or in a
%   quarterfield region. A pre-determined poisson process generates random
%   flashes of rotating gratings throughout the visual field at low or high
%   contrast.
%
%   The probe stimuli have three sizes 1x1 deg, 2x2 deg, or 4x4 deg, to
%   help estimate different RF sizes and are placed at 2 degree
%   increments. The probability of a probe stimulus turning on is 2.5% per
%   TR and the dead time is five seconds. Using poisson processes in this
%   way ensures a random distribution at every location that is
%   uncorrelated to all other locations. The code starts with a 10 s blank
%   and has another 10 s blank every three minutes. Probes are at 20% and
%   80% contrast, each probe lasts two TRs.
%   
%   The attention task involves performing orientation judgments on gabors
%   at a location cued continuously by a circular aperture. Fixation is
%   maintained at the center and monitored within a 1.5 deg window. Gabor
%   is at full contrast to differentiate from the probe stimuli (and to
%   minimize the effect of the probes on performance)
%
%   The goal of the task is to provide an independent pRF dataset that can
%   be used to build channel encoding models that cover space, contrast,
%   and attention. The end goal is to link frontal areas that control
%   attention with visual areas. 
%
%   Project goals:
%       (1) Show that we can map pRFs hierarchically throughout all of
%       cortex
%       (2) Build a visuospatial encoding model that incorporates attention
%       (3) Show that our model predicts tuning shifts for non-visuospatial
%       features due to spatial attention

%%

global stimulus

stimulus = struct;

%% Initialize Variables

% add arguments later
scan = 0;
plots = 0;
noeye = 0;
debug = 0;
getArgs(varargin,{'scan=0','plots=0','noeye=0','debug=0'});
stimulus.scan = scan;
stimulus.plots = plots;
stimulus.noeye = noeye;
stimulus.debug = debug;
clear localizer invisible scan noeye task test2

if stimulus.scan
    warning('Not setup for scanning');
end

%% Stimulus parameters

stimulus.stimX = 24; % max ecc in any direction
stimulus.stimY = 12;
stimulus.stimR = 2; % deg between each stimulus

if mod(stimulus.stimX,stimulus.stimR)==1 || mod(stimulus.stimY,stimulus.stimR)==1
    warning('Your stimulus size is not correctly setup');
end

stimulus.stimx = -stimulus.stimX:stimulus.stimR:stimulus.stimX;
stimulus.stimy = -stimulus.stimY:stimulus.stimR:stimulus.stimY;

stimulus.probeOn = .025;
stimulus.probeUp = 2;
stimulus.probeDown = 10;

stimulus.live.grid = zeros(length(stimulus.stimx),length(stimulus.stimy));
stimulus.live.gridCon = zeros(length(stimulus.stimx),length(stimulus.stimy));
stimulus.live.gridCons = cell(length(stimulus.stimx),length(stimulus.stimy));
stimulus.live.gridSize = zeros(length(stimulus.stimx),length(stimulus.stimy));
stimulus.live.gridSizes = cell(length(stimulus.stimx),length(stimulus.stimy));

stimulus.gratingContrasts = [0.1 0.2 1.0];
stimulus.gratingSizes = [1 2 4];

%% Open Old Stimfile
stimulus.counter = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/afmap/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/afmap/%s/1*mat',mglGetSID));

    if length(files) >= 1
        fname = files(end).name;
        
        s = load(sprintf('~/data/afmap/%s/%s',mglGetSID,fname));
        % copy staircases and run numbers
        stimulus.counter = s.stimulus.counter + 1;
        stimulus.staircase = s.stimulus.staircase;
        clear s;
        disp(sprintf('(posjdg) Data file: %s loaded.',fname));
    end
end
disp(sprintf('(posjdg) This is run #%i',stimulus.counter));

%% Setup Screen
myscreen = initScreen('VPixx');

% set background to grey
myscreen.background = 0.5;

%% Staircase
if ~isfield(stimulus,'staircase')
    disp('(posjdg) WARNING: New staircase');
    initStair();
else
    resetStair();
end

%% Plot and return
if stimulus.plots==2
    dispInfo(stimulus);
    return
end

%% Initialize Stimulus

myscreen.stimulusNames{1} = 'stimulus';

localInitStimulus();
    
stimulus.responseKeys = [1 2]; % left right

%% Colors
stimulus.colors.white = [1 1 1]; stimulus.colors.red = [1 0 0];
stimulus.colors.green = [0 1 0]; stimulus.colors.black = [0 0 0];
% % initGammaTable(myscreen);
% stimulus.colors.rmed = 127.5;
% 
% % We're going to add an equal number of reserved colors to the top and
% % bottom, to try to keep the center of the gamma table stable.
% stimulus.colors.reservedBottom = [1 0 0; 0 0 0]; % fixation cross colors
% stimulus.colors.reservedTop = [1 1 1; 0 1 0]; % correct/incorrect colors
% stimulus.colors.black = 1/255; stimulus.colors.white = 254/255;
% stimulus.colors.red = 0/255; stimulus.colors.green = 255/255;
% stimulus.colors.nReserved = 2; % this is /2 the true number, because it's duplicated
% stimulus.colors.nUnreserved = 256-(2*stimulus.colors.nReserved);
% 
% stimulus.colors.mrmax = stimulus.colors.nReserved - 1 + stimulus.colors.nUnreserved;
% stimulus.colors.mrmin = stimulus.colors.nReserved;

%% Setup Probe Task

task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;
% task waits for fixation on first segment
task{1}{1}.seglen = 0.450;

stimulus.seg.stim = 1;

task{1}{1}.synchToVol = 1;
task{1}{1}.getResponse = 0;

task{1}{1}.numTrials = Inf;

task{1}{1}.random = 0;

if stimulus.scan
    task{1}{1}.synchToVol = 1;
end

task{1}{1}.randVars.calculated.probesOn = nan;

%% Setup Attention Task

task{2}{1} = struct;
task{2}{1}.waitForBacktick = 0;
% task waits for fixation on first segment
task{2}{1}.segmin = [0.200 1 0.200 1 1];
taks{2}{1}.segmax = [0.200 1 0.200 1 1];

stimulus.seg.stim1 = 1;
stimulus.seg.delay1 = 2;
stimulus.seg.stim2 = 3;
stimulus.seg.delay2 = 4;
stimulus.seg.resp = 5;

task{2}{1}.synchToVol = [0 0 0 0 0];
task{2}{1}.getResponse = [0 0 0 0 1];

task{2}{1}.numTrials = Inf;

task{2}{1}.random = 0;

if stimulus.scan
    task{2}{1}.synchToVol = 1;
end

task{2}{1}.randVars.calculated.probesOn = nan;

%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,[],@screenUpdateCallback1,[],@startTrialCallback,[],[]);
    [task{2}{phaseNum}, myscreen] = initTask(task{2}{phaseNum},myscreen,@startSegmentCallback2,@screenUpdateCallback2,@getResponseCallback2,@startTrialCallback2,[],[]);
end

%% EYE CALIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

% let the user know
disp(sprintf('(posjdg) Starting run number: %i.',stimulus.counter));

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
    [task{2}, myscreen, phaseNum] = updateTask(task{2},myscreen,phaseNum);
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
    disp('(posjdg) Displaying plots');
    dispInfo(stimulus);
end

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startTrialCallback1(task,myscreen)
%%
global stimulus

% Design the state space
for x = 1:length(stimulus.stimx)
    for y = 1:length(stimulus.stimy)
        % increment
        if stimulus.live.grid(x,y) > 1
            stimulus.live.grid(x,y) = stimulus.live.grid(x,y)-1;
        elseif stimulus.live.grid(x,y) == 1
            % shut down grid location
            stimulus.live.grid(x,y) = 0;
            stimulus.live.gridCon(x,y) = 0;
            stimulus.live.gridSize(x,y) = 0;
        else
            if rand < stimulus.probeOn
                % turn on grid location
                stimulus.live.grid(x,y) = stimulus.probeDown + stimulus.probeUp;
                % pick attributes
                cOpts = stimulus.gratingContrasts;
                sOpts = stimulus.gratingSizes;
            end
        end
    end
end
    
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
    if task.thistrial.visible == 1
        stimulus.live.stim = 1;
    end
elseif task.thistrial.thisseg==stimulus.seg{task.thistrial.thisphase}.resp
    stimulus.live.resp = 1;
    stimulus.live.angle = task.thistrial.startRespAngle;
    stimulus.live.anyAngleAdj = false;
    mInfo = mglGetMouse(myscreen.screenNumber);
    stimulus.live.trackingAngle = -mInfo.x/90;
    convertRespXY(task);
elseif stimulus.att == 2 && task.thistrial.thisseg==stimulus.seg{task.thistrial.thisphase}.cue
    stimulus.live.cue = 1;
end

for i = 1:2
    mglClearScreen(0.5);
    if stimulus.live.stim
        x = task.thistrial.ecc * cos(task.thistrial.angle+task.thistrial.target);
        y = task.thistrial.ecc * sin(task.thistrial.angle+task.thistrial.target);
        mglBltTexture(stimulus.live.grating,[x y],0,0,task.thistrial.rotation*180/pi);
    elseif stimulus.live.cue && i == 1
        x = task.thistrial.ecc * cos(task.thistrial.target);
        y = task.thistrial.ecc * sin(task.thistrial.target);
        mglPolygon([x-.15, x-.15, x+.15, x+.15], [y-.15, y+.15, y+.15, y-.15], stimulus.colors.white);
        mglFlush();
        mglClearScreen(0.5);
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

if (task.thistrial.thisseg==stimulus.seg{task.thistrial.thisphase}.resp) && stimulus.powerwheel
    mInfo = mglGetMouse(myscreen.screenNumber);
    curPos = -mInfo.x/90;
    stimulus.live.angle = stimulus.live.angle + curPos-stimulus.live.trackingAngle;
    if abs(curPos-stimulus.live.trackingAngle)>0
        stimulus.live.anyAngleAdj = true;
    end
    stimulus.live.trackingAngle = curPos;
    convertRespXY(task);
elseif task.thistrial.thisseg==stimulus.seg{task.thistrial.thisphase}.resp
    keys = find(mglGetKeys);
    if any(keys==19)
        stimulus.live.angle = stimulus.live.angle+0.01;
        stimulus.live.anyAngleAdj = true;
    elseif any(keys==20)
        stimulus.live.angle = stimulus.live.angle-0.01;
        stimulus.live.anyAngleAdj = true;
    end
    convertRespXY(task);
end

if stimulus.live.resp && ((task.thistrial.thisphase==2) || stimulus.test2)
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
    
if isfield(task.thistrial,'whichButton') && (task.thistrial.whichButton==stimulus.responseKeys(1))
    % subject didn't see anything
    task = jumpSegment(task,inf);
end

if stimulus.powerwheel
    validResponse = task.thistrial.mouseButton == 1;
else
    validResponse = any(task.thistrial.whichButton == stimulus.responseKeys);
end

if validResponse
    if stimulus.live.gotResponse==0
        if (task.thistrial.thisphase==1 && (~stimulus.test2))
            % they saw it
            task.thistrial.detected = 1;
            stimulus.staircase = doStaircase('update',stimulus.staircase,task.thistrial.detected);
            disp(sprintf('Subject saw %01.2f%% contrast',task.thistrial.contrast*100));
            stimulus.live.fix = 0;
        elseif ((task.thistrial.thisphase==2) || stimulus.test2)
            % they saw it
%             if stimulus.live.anyAngleAdj == true
%                 task.thistrial.detected = 1;
%                 stimulus.staircase = doStaircase('update',stimulus.staircase,task.thistrial.detected);
%             end
            % they are actually reporting locations
            task.thistrial.respAngle = stimulus.live.angle;
            disp(sprintf('Subject reported %02.0f real %02.0f at %02.0f%% contrast',task.thistrial.respAngle*180/pi,task.thistrial.angle*180/pi,task.thistrial.contrast*100));
            stimulus.live.fix = 0;
            stimulus.live.resp = 0;
            task = jumpSegment(task,inf);
        end
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
            'initialStepsize',0.025,...
            'minThreshold=0.0001','maxThreshold=0.4','stepRule','pest',...
            'nTrials=40','maxStepsize=0.2','minStepsize=0.0001');
        
function resetStair()

global stimulus

if doStaircase('stop',stimulus.staircase)
    disp('(posjdg) Staircase is being reset');
    stimulus.staircase(end+1) = doStaircase('init',stimulus.staircase(end));
    if stimulus.staircase(end).s.threshold>0.3
        disp('(posjdg) Bad staircase threshold: setting to 0.3');
        stimulus.staircase(end).s.threshold=0.3;
    elseif stimulus.staircase(end).s.threshold<0
        disp('(posjdg) Bad staircase threshold: setting to 0.05');
        stimulus.staircase(end).s.threshold=0.05;
    end
end

function [trials] = totalTrials()
%%

% Counts trials + estimates the threshold based on the last 500 trials

% get the files list
files = dir(fullfile(sprintf('~/data/posjdg/%s/17*stim*.mat',mglGetSID)));

trials = 0;

for fi = 1:length(files)
    load(fullfile(sprintf('~/data/posjdg/%s/%s',mglGetSID,files(fi).name)));
    
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
files = dir(fullfile(sprintf('~/data/posjdg_%s/%s/17*stim*.mat',rstimulus.condition,mglGetSID)));

% load the files and pull out the data (long form)
%  rrun # counter #    local trial     real trial   angle     respAngle    
%     1       2             3              4           5           6
%  target    startRespAngle     contrast     detected      ecc    priorsd
%     7            8                9           10          11      12
%    rotation
%       13
count = 1; data = zeros(10000,13);

for fi = 1:length(files)
    load(fullfile(sprintf('~/data/posjdg_%s/%s/%s',rstimulus.condition,mglGetSID,files(fi).name)));
    
    e = getTaskParameters(myscreen,task);
    if e{1}(1).nTrials>1
        e = e{1}(2); % why?!
    
        run = stimulus.counter;

        data(count:count+(e.nTrials-1),:) = [repmat(fi,e.nTrials,1) repmat(run,e.nTrials,1) (1:e.nTrials)' (count:count+(e.nTrials-1))' ...
            e.randVars.angle' e.randVars.respAngle' e.parameter.target' ...
            e.randVars.startRespAngle' e.randVars.contrast' e.randVars.detected' ...
            e.parameter.ecc' e.parameter.priorSTD' e.randVars.rotation'];

        count = count+e.nTrials;
    end
end

data = data(1:(count-1),:);

if any(data(:,6)>pi), data(data(:,6)>pi,6) = data(data(:,6)>pi,6)-2*pi; end
%% Compute angle-target and respAngle-target plot
h = figure; hold on

low = [0 0.075 0.081 0.09 inf];

for i = 1:4
    subplot(4,1,i); hold on
    % data(:,6) = data(:,6)-data(:,5);

    % remove no-response trials
    data_ = data(~isnan(data(:,6)),:);
    data_ = data_(logical((data_(:,9)>low(i)).*(data_(:,9)<low(i+1))),:);
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

    axis([-1 1 -1 1]);
    axis square

    set(gca,'XTick',-1:.5:1,'XTickLabel',round((-1:.5:1)*180/pi,2),'YTick',-1:.5:1,'YTickLabel',round((-1:.5:1)*180/pi,2));

    drawPublishAxis;
    % h = figure;
    % hist(data(:,5)-data(:,6));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localInitStimulus()

global stimulus

for ci = 1:length(stimulus.gratingContrasts)
    for si = 1:length(stimulus.gratingSizes)
        sz = 1.5 * stimulus.gratingSizes(si);
        % use total degs / num to compute size
        grating = stimulus.gratingContrasts(ci) * 255/2 * mglMakeGrating(sz,sz,2,0) + 255/2;
        gauss = mglMakeGaussian(sz,sz,sz/6,sz/6);
        alphamask = repmat(grating,1,1,4);
        alphamask(:,:,4) = gauss*255;

        stimulus.live.grating(ci,si)  = mglCreateTexture(alphamask); % high contrast
    end
end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % sets the gamma table so that we can have
% % finest possible control over the stimulus contrast.
% %
% % stimulus.colors.reservedColors should be set to the reserved colors (for cue colors, etc).
% % maxContrast is the maximum contrast you want to be able to display.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function setGammaTable(maxContrast)
% 
% global stimulus;
% 
% % set the bottom
% gammaTable(1:size(stimulus.colors.reservedBottom,1),1:size(stimulus.colors.reservedBottom,2)) = stimulus.colors.reservedBottom;
% 
% % set the gamma table
% if maxContrast == 1
%     % create the rest of the gamma table
%     cmax = 1;cmin = 0;
%     luminanceVals = cmin:((cmax-cmin)/(stimulus.colors.nUnreserved-1)):cmax;
% 
%     % now get the linearized range
%     redLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,luminanceVals,'linear');
%     greenLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,luminanceVals,'linear');
%     blueLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,luminanceVals,'linear');
% elseif maxContrast > 0
%     % create the rest of the gamma table
%     cmax = 0.5+maxContrast/2;cmin = 0.5-maxContrast/2;
%     luminanceVals = cmin:((cmax-cmin)/(stimulus.colors.nUnreserved-1)):cmax;
% 
%     % now get the linearized range
%     redLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,luminanceVals,'linear');
%     greenLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,luminanceVals,'linear');
%     blueLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,luminanceVals,'linear');
% else
%     % if we are asked for 0 contrast then simply set all the values to gray
%     redLinearized = repmat(interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,.5,'linear'),1,stimulus.colors.nUnreserved);
%     greenLinearized = repmat(interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,.5,'linear'),1,stimulus.colors.nUnreserved);
%     blueLinearized = repmat(interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,.5,'linear'),1,stimulus.colors.nUnreserved);
% end
% 
% % add to the table!
% gammaTable((stimulus.colors.mrmin:stimulus.colors.mrmax)+1,:)=[redLinearized;greenLinearized;blueLinearized]';
% 
% % set the top
% gammaTable = [gammaTable; stimulus.colors.reservedTop];
% 
% if size(gammaTable,1)~=256
%     disp('(setGammaTable) Failure: Incorrect number of colors in gamma table produced');
%     keyboard
% end
% 
% % set the gamma table
% succ = mglSetGammaTable(gammaTable);
% 
% if ~succ
%     warning('Gamma table set failure');
%     keyboard
% end
% 
% % remember what the current maximum contrast is that we can display
% stimulus.curMaxContrast = maxContrast;
% 
% 
% function initGammaTable(myscreen)
% global stimulus
% %% Gamma Table Initialization
% 
% % get gamma table
% if ~isfield(myscreen,'gammaTable')
%   stimulus.linearizedGammaTable = mglGetGammaTable;
%   disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
%   disp(sprintf('(cuecon:initGratings) No gamma table found in myscreen. Contrast displays like this'));
%   disp(sprintf('         should be run with a valid calibration made by moncalib for this monitor.'));
%   disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
% else
%   % check to make sure this matches the calibration file
%   
%   % get each channel table that should have been set by mglGetGammaTable
%   redTable = myscreen.initScreenGammaTable.redTable(:);
%   greenTable = myscreen.initScreenGammaTable.greenTable(:);
%   blueTable = myscreen.initScreenGammaTable.blueTable(:);
%   % get what the calibration structure says it should have been set to
%   gammaTable = myscreen.gammaTable(:);
%   % table values are only good to 10 bits
%   redTable = round(redTable*1024)/1024;
%   greenTable = round(greenTable*1024)/1024;
%   blueTable = round(blueTable*1024)/1024;
%   gammaTable = round(gammaTable*1024)/1024;
%   % compare, ignoring nans
%   if ~isequaln(mglGetGammaTable,myscreen.initScreenGammaTable) || ~isequaln(redTable,gammaTable) || ~isequaln(greenTable,gammaTable) || ~isequaln(blueTable,gammaTable)
%     disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
%     disp(sprintf('(curecon:initGrating) Gamma table does not match calibration'));
%     disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
%     keyboard
%   end
% end
% stimulus.linearizedGammaTable = myscreen.initScreenGammaTable;
