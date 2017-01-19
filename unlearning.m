function [ myscreen ] = unlearning( varargin )
%UNLEARNING 
%
% Unconscious learning of spatial patterns. This experiment runs two
% simultaneous tasks. One task asks subjects to notice and press a key
% (spacebar) when any location on the polar grid turns orange. This task is
% not staircased, but is sufficiently hard that it requires attending to
% all of the grid locations.
%
% The second task is the real unconscious learning task. The goal is to
% learn delayed match-to-sample between two patterns along a particular
% polar angle. Each pattern is an nCk generated bitwise pattern, randomly
% chosen to match or non-match. The patterns are drawn with dividers
% between them so that they can be (theoretically) isolated with receptive
% field mapping to unique sets of responsive voxels.
%
% Failing to identify a spacebar press causes a 10-s timeout and a loud
% beep noise to strongly encourage participants to attend to this task and
% not the other task.
%
% Participants are given no instructions about the pattern matching task.
%
% Neither task is difficulty staircase.
%

global stimulus

%% Open Old Stimfile
stimulus.initStair = 1;
stimulus.counter = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/unlearning/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/unlearning/%s/1*mat',mglGetSID));

    if length(files) >= 1
        fname = files(end).name;
        
        s = load(sprintf('~/data/unlearning/%s/%s',mglGetSID,fname));
        % copy staircases and run numbers
        stimulus = s.stimulus;
        stimulus.counter = stimulus.counter + 1;

        clear s;
        disp(sprintf('(unlearn) Data file: %s loaded.',fname));
        
    end
end
disp(sprintf('(unlearn) This is run #%i',stimulus.counter));

%% EXPERIMENT PARAMETERS
if ~isfield(stimulus,'cur'), stimulus.cur = {}; end

stimulus.cur{end+1} = struct;
stimulus.cur{end}.N = 4;
stimulus.cur{end}.K = 2;
stimulus.cur{end}.angle = 30;
stimulus.cur{end}.num = 360/stimulus.cur{end}.angle;
stimulus.cur{end}.buffer = 5; % buffer is used to stencil over the wedges
stimulus.cur{end}.isize = 1;
stimulus.cur{end}.osize = 10;

if ~isfield(stimulus,'learn')
    stimulus.learn = randi(stimulus.cur_.num);
    disp('(unlearn) Re-setting learn position');
end
disp(sprintf('(unlearn) Subject %s is learning %i',mglGetSID,stimulus.learn));

stimulus.cur_ = stimulus.cur{end};

%% Initialize Variables

% add arguments later
scan = 0;
plots = 0;
noeye = 0;
getArgs(varargin,{'scan=0','plots=0','noeye=1'});
stimulus.scan = scan;
stimulus.plots = plots;
stimulus.noeye = noeye;
clear localizer invisible scan category noeye task

if stimulus.scan
    warning('Not setup for scanning');
end

%% Setup missing initial variables

if ~isfield(stimulus,'counter')
    stimulus.counter = 1; % This keeps track of what "run" we are on.
end

%% Setup Screen

if stimulus.scan
    myscreen = initScreen('fMRIprojFlex');
else
    myscreen = initScreen('VPixx');
end

% set background to grey
myscreen.background = 0.5;

%% Plot and return

if stimulus.plots==2
    dispInfo(stimulus);
    return
end

%% Initialize Stimulus

myscreen = initStimulus('stimulus',myscreen);
    
if stimulus.scan
    stimulus.responseKeys = [2 1]; % corresponds to NOMATCH, MATCH
else
    stimulus.responseKeys = [2 1]; % corresponds to  NOMATCH, MATCH
end

stimulus.colors.black = [0 0 0];
stimulus.colors.white = [1 1 1];
stimulus.colors.green = [0 1 0];
stimulus.colors.red = [1 0 0];

%% Generate stencils
% The stencil is a series of arcs 
mglStencilCreateBegin(1);
% Draw an annulus at every buffer location
for i = 0:(stimulus.cur_.num-1)
    mglGluPartialDisk(0,0,stimulus.cur_.isize,stimulus.cur_.osize,i*stimulus.cur_.angle+stimulus.cur_.buffer/2,stimulus.cur_.angle-stimulus.cur_.buffer/2,[1 1 1],60,2);
end
mglStencilCreateEnd;
mglClearScreen(0.5);

%% Setup Task
task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;

stimulus.curTrial = 0;

task{1}{1}.segmin = [inf .650 1.000 .650 0.500 1.500 0.500];
task{1}{1}.segmax = [inf .650 1.000 .650 0.500 1.500 6.000];

if stimulus.noeye==1
    task{1}{1}.segmin(1) = 0.5;
    task{1}{1}.segmax(1) = 0.5;
end

stimulus.seg.ITI1 = 1; % waits for user input (button press + held) and eye fixation (within 2 degrees)
stimulus.seg.stim1 = 2;
stimulus.seg.delay1 = 3;
stimulus.seg.stim2 = 4;
stimulus.seg.delay2 = 5;
stimulus.seg.resp = 6;
stimulus.seg.ITI2 = 7;

task{1}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
task{1}{1}.getResponse = zeros(size(task{1}{1}.segmin)); task{1}{1}.getResponse(stimulus.seg.resp)=1;
task{1}{1}.numTrials = 50;
task{1}{1}.random = 1;
task{1}{1}.parameter.match = [0 1];
task{1}{1}.parameter.count = [0 0 0 0 0 0 1 2];

if stimulus.scan
    task{1}{1}.synchToVol(stimulus.seg.ITI) = 1;
end

stimulus.patterns = unique(perms([ones(1,stimulus.cur_.K) zeros(1,stimulus.cur_.N-stimulus.cur_.K)]),'rows');
stimulus.npatterns = size(stimulus.patterns,1);
task{1}{1}.parameter.pattern1 = 1:stimulus.npatterns; % which test pattern to use for the first stim

%% Tracking

% these are variables that we want to track for later analysis.
task{1}{1}.randVars.calculated.correct = nan;
task{1}{1}.randVars.calculated.cumcount = nan;
task{1}{1}.randVars.calculated.pattern2 = nan;

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

%% Get Ready...
% clear screen    
% % % mglWaitSecs(1);
% % % mglFixationCross(0.1,0.1,stimulus.colors.white);
% % % if stimulus.scan        
% % %     mglTextDraw('DO NOT MOVE',[0 1.5]);
% % % end
% % % mglFlush

% let the user know
disp(sprintf('(unlearn) Starting run number: %i.',stimulus.counter));
% if stimulus.unattended
% % myscreen.flushMode = 0;
% end

%% Main Task Loop

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
count = 0;
while true
    mglTextDraw('Report your count: 0',[0 0]);
    mglFlush
    keys = mglGetKeys;
    if keys(19)
        count = count+1;
    end
    if keys(20)
        break;
    end
end
mglTextDraw(sprintf('True count: %i',task.thistrial.cumcount),[0 2]);
myscreen.flushMode = 1;

stimulus.count = count;
stimulus.truecount = task.thistrial.cumcount;
stimulus.off = abs(stimulus.count-stimulus.truecount);
disp(sprintf('(unlearn) Subject was off by: %i. True: %i Report: %i',stimulus.off,stimulus.truecount,stimulus.count));

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

if stimulus.plots
    disp('(unlearn) Displaying plots');
    dispInfo(stimulus);
end

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startTrialCallback(task,myscreen)
%%

global stimulus

stimulus.curTrial = stimulus.curTrial + 1;

myscreen.flushMode = 0;

% Setup the displays
% .rings holds N rings consisting of num segments
stimulus.live.rings1 = zeros(stimulus.cur_.N,stimulus.cur_.num);
stimulus.live.rings2 = stimulus.live.rings1; % copy
for i = 1:stimulus.cur_.num
    if stimulus.learn==i
        stimulus.live.rings1(:,i) = stimulus.patterns(task.thistrial.pattern1,:)';
        if task.thistrial.match==1
            % if match, set equal to current pattern
            stimulus.live.rings2(:,i) = stimulus.live.rings1(:,i);
            task.thistrial.pattern2 = task.thistrial.pattern1;
        else
            notpatterns = 1:stimulus.npatterns;
            notpatterns = notpatterns(notpatterns~=task.thistrial.pattern1);
            task.thistrial.pattern2 = notpatterns(randi(length(notpatterns)));
            stimulus.live.rings2(:,i) = stimulus.patterns(task.thistrial.pattern2,:)';
        end
    else
        fpattern = randi(stimulus.npatterns);
        stimulus.live.rings1(:,i) = stimulus.patterns(fpattern,:)';
        if randi(2)==1
            stimulus.live.rings2(:,i) = stimulus.live.rings1(:,i);
        else
            notpatterns = 1:stimulus.npatterns;
            notpatterns = notpatterns(notpatterns~=fpattern);
            stimulus.live.rings2(:,i) = stimulus.patterns(notpatterns(randi(length(notpatterns))),:)';
        end
    end
end

opts = {'Non-match','Match'};
disp(sprintf('(unlearn) %s trial. Pattern A: %i, pattern B: %i',opts{task.thistrial.match+1},task.thistrial.pattern1,task.thistrial.pattern2));

stimulus.live.eyeCount = 0;
stimulus.dead = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)
%%

global stimulus

stimulus.live.triggerWaiting = 0;
if any(task.thistrial.thisseg==[stimulus.seg.ITI1])
    stimulus.live.triggerWaiting = 1;
    stimulus.live.centered = 0;
    stimulus.live.triggerTime = 0;
    stimulus.live.lastTrigger = -1;
end

stimulus.live.resp = 0;
stimulus.live.fixColor = stimulus.colors.white;
stimulus.live.fix = 1;
stimulus.live.stim = 0;

if any(task.thistrial.thisseg==[stimulus.seg.stim1 stimulus.seg.stim2])
    stimulus.live.stim = 1;
elseif task.thistrial.thisseg==stimulus.seg.resp
    stimulus.live.fix = 0;
elseif task.thistrial.thisseg==stimulus.seg.ITI2
    stimulus.live.fix = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus
mglClearScreen(0.5);
% check eye pos
% % % if ~stimulus.noeye && task.thistrial.thisseg~=stimulus.seg.ITI && ~stimulus.scan
% % %     [pos,time] = mglEyelinkGetCurrentEyePos;
% % %     if ~any(isnan(pos))
% % %         dist = hypot(pos(1),pos(2));
% % %         if dist > stimulus.ring.inner && stimulus.live.eyeCount > 30
% % %             mglTextSet([],32,stimulus.colors.red);
% % %             disp('Eye movement detected!!!!');
% % %             mglTextDraw('Eye Movement Detected',[0 0]);
% % %             mglFlush
% % %             myscreen.flushMode = 1;
% % %             stimulus.dead = 1;
% % %             return
% % %         elseif dist > stimulus.ring.inner-1
% % %             stimulus.live.eyeCount = stimulus.live.eyeCount + 1;
% % %         end
% % %     end
% % % end

if stimulus.live.stim
    if task.thistrial.thisseg==stimulus.seg.stim1
        rings = stimulus.live.rings1;
    else
        rings = stimulus.live.rings2;
    end
    mglStencilSelect(1);
    % draw rings
    
    % revert stencil
    mglStencilSelect(0);
end

if stimulus.live.fix
%      cover
    mglFillOval(0,0,[1 1],0.5);
    upFix(stimulus);
end

% if stimulus.live.triggerWaiting
%     now = mglGetSecs;
%     % check eye position, if 
%     [pos,time] = mglEyelinkGetCurrentEyePos;
%     if ~any(isnan(pos))
%         dist = hypot(pos(1),pos(2));
%         wasCentered = stimulus.live.centered;
%         stimulus.live.centered = dist<2;
%         if wasCentered && stimulus.live.centered && stimulus.live.lastTrigger>0
%             stimulus.live.triggerTime = stimulus.live.triggerTime + now-stimulus.live.lastTrigger;
%         end
%         stimulus.live.lastTrigger = now;
%     end
%     if stimulus.live.triggerTime > 0.75 % not in ms dummy, wait 1.5 seconds (reasonable slow time)
%         disp('Eye position centered');
%         task = jumpSegment(task);
%     end
% end

function upGrating(stimulus,task)

% phaseNum = floor(length(stimulus.grating.phases)*rem(mglGetSecs(task.thistrial.trialstart)*stimulus.grating.tf,1)+1);
mglBltTexture(stimulus.tex(stimulus.live.phaseNum),[stimulus.live.x stimulus.live.y 6],0,0,stimulus.live.dir*180/pi+90); % we add 90 so that it's aligned with the motion
mglBltTexture(stimulus.mask,[stimulus.live.x stimulus.live.y 6 6],0,0,stimulus.live.dir*180/pi+90);

function upFix(stimulus)
%%
% for this experiment use a circle to indicate where participants can
% fixate inside of (rather than a cross which might arbitrarily enforce
% poisitioning
% mglGluAnnulus(0,0,1.5,1.55,stimulus.live.fixColor,64);
mglFixationCross(1,1,stimulus.live.fixColor);

function stimulus = upDots(stimulus,myscreen)

stimulus.dots{stimulus.dot} = updateDotsRadial(stimulus.dots{stimulus.dot},stimulus.live.coherence,myscreen,true);

mglPoints2(stimulus.dots{stimulus.dot}.x+stimulus.live.x,stimulus.dots{stimulus.dot}.y+stimulus.live.y,...
    stimulus.dots{stimulus.dot}.dotsize,stimulus.live.dotColor);

mglBltTexture(stimulus.mask,[stimulus.live.x stimulus.live.y 6 6],0,0,stimulus.live.dir*180/pi+90);

% function stimulus = upDotsInc(stimulus,myscreen)
% stimulus.idots = updateDotsRadial(stimulus.idots,0,myscreen,true);
% mglPoints2(stimulus.idots.x,stimulus.idots.y,...
%     stimulus.idots.dotsize,stimulus.live.dotColor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus

if stimulus.dead, return; end
responseText = {'Incorrect','Correct'};
sideText = {'Left','Right'};
fixColors = {stimulus.colors.red,stimulus.colors.green};
    
responses = [1 0 2];
if any(task.thistrial.whichButton == stimulus.responseKeys)
    if task.thistrial.gotResponse == 0
        task.thistrial.correct = task.thistrial.whichButton == stimulus.responseKeys(responses(task.thistrial.rotation+2));
        disp(sprintf('Subject pressed %i: %s %s',task.thistrial.whichButton,sideText{task.thistrial.whichButton},responseText{task.thistrial.correct+1}));
        stimulus.live.fixColor = fixColors{task.thistrial.correct+1};
        stimulus.staircase{stimulus.feature} = doStaircase('update',stimulus.staircase{stimulus.feature},task.thistrial.correct);
        stimulus.live.resp = 1;
        stimulus.live.dots = 0;
        stimulus.live.fix=1;
    else
        disp(sprintf('(unlearn) Subject responded multiple times: %i',task.thistrial.gotResponse+1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              HELPER FUNCTIONS                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
%    initStaircase     %
%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initStaircase(stimulus)
%%
stimulus.staircase = {};
stimulus.staircase{1} = doStaircase('init','upDown',...
        'initialThreshold',0.25,... % radians of rotation (
        'initialStepsize',0.05,...
        'minThreshold=0.001','maxThreshold=0.35','stepRule','pest',...
        'nTrials=50','maxStepsize=0.1','minStepsize=0.001');
stimulus.staircase{2} = stimulus.staircase{1};

function stimulus = resetStair(stimulus)

for i = 1:2
    if doStaircase('stop',stimulus.staircase{i})
        disp('(unlearn) Initializing new staircase...');
        stimulus.staircase{i}(end+1) = doStaircase('init',stimulus.staircase{i}(end));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
%    dispInfo    %
%%%%%%%%%%%%%%%%%%%%%%%
function dispInfo(stimulus)
%%

if ~stimulus.localizer && ~stimulus.staircasing
    disp(sprintf('Participant %s has earned $%2.2f',mglGetSID,stimulus.run.points/100));
end
% load the luminance table
% % % load(myscreen.calibFullFilename)
% % % luminance = interp1(calib.tableCorrected.outputValues,calib.tableCorrected.luminance,0:1/255:255);
if stimulus.staircasing
    %%
    notstaircase = stimulus.staircase;
    thresholds = zeros(size(stimulus.run.stimLengths));
    for i = 1:length(stimulus.staircase)
        out = doStaircase('threshold',notstaircase{i},'type','weibull','dispFig=0');
        thresholds(i) = out.threshold;
    end
    % reorganize into matrix
    stimCons = unique(stimulus.run.stimCon);
    stimCons = sort(stimCons);
    stimLengths = unique(stimulus.run.stimLengths);
    stimLengths = sort(stimLengths);
    datamat = nan(length(stimCons),length(stimLengths),5);
    for ci = 1:length(stimCons)
        for li = 1:length(stimLengths)
            idxs = logical((stimulus.run.stimLengths==stimLengths(li)) .* (stimulus.run.stimCon==stimCons(ci)));
            datamat(ci,li,1:sum(idxs)) = thresholds(idxs);
        end
    end
    if any(thresholds<0) || any(thresholds>1)
        % remove errant thresholds
        warning('should remove some thresholds...');
    end
    datamat(datamat>1) = NaN;
    datamat(datamat<=0) = NaN;
    %%
    datamu = nanmean(datamat,3);
    datamu(datamu==0) = NaN;
    datamu = round((1-datamu)*255);
    datasd = nanstd(datamat,[],3);
    datasd(datasd==0) = NaN;
    %% plot
    cmap = brewermap(length(stimCons)+1,'Purples');
    cmap = cmap(2:end,:);
    figure, hold on
    legs = {};
    for i = 1:length(stimCons)
        plot(stimLengths,datamu(i,:),'o','MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',[1 1 1],'MarkerSize',10);
        errbar(stimLengths,datamu(i,:),datasd(i,:),'-','Color',cmap(i,:));
        legs{end+1} = sprintf('Stimulus luminance: %i/255',stimCons(i));
    end
    a = axis;
    axis([50 100 0 a(4)]);
    legend(legs)
    xlabel('Stimulus length (ms)');
    ylabel('Mask contrast at just noticeable difference (% luminance)');
    drawPublishAxis
elseif stimulus.localizer
else
%     perf = zeros(size(stimulus.istaircase));
%     for i = 1:length(stimulus.istaircase)
%         perf(i) = mean(stimulus.istaircase(i).response);
%     end
%     figure
%     plot(1:length(perf),perf,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
%     set(gca,'XAxisTick',1:length(perf));
%     drawPublishAxis
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create dots for optic flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initDots(dots,~)

% maximum depth of points
dots.dir = 0;

area = (dots.maxX-dots.minX)*(dots.maxY-dots.minY);

dots.n = area * dots.density;

dots.x = rand(1,dots.n)*(dots.maxX-dots.minX)+dots.minX;
dots.y = rand(1,dots.n)*abs(dots.maxY-dots.minY)+dots.minY;

dots.x = dots.x;
dots.y = dots.y;

% set incoherent dots to 0
dots.coherency = 1;
dots.incoherent = rand(1,dots.n) > dots.coherency;
dots.incoherentn = sum(dots.incoherent);
dots.coherent = ~dots.incoherent;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots for Radial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updateDots(dots,coherence,myscreen,repick)

% get the coherent and incoherent dots
if coherence==1
    dots.coherent = true(size(dots.x));
    dots.incoherent = false(size(dots.x));
    dots.incoherentn = 0;
elseif repick
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

freq_factor = dots.speed/myscreen.framesPerSecond;

if dots.coherentn>0
    % move coherent dots
    % dots.dir is an angle, so compute x/y based on this
    dots.y(dots.coherent) = dots.y(dots.coherent) + (freq_factor * sin(dots.dir));
    dots.x(dots.coherent) = dots.x(dots.coherent) + (freq_factor * cos(dots.dir));
end

if dots.incoherentn>0
    % these are for flipping into the other quadrants
    rdir = rand(1,dots.incoherentn)*pi*2;
    dots.y(dots.incoherent) = dots.y(dots.incoherent) + (freq_factor * sin(rdir));
    dots.x(dots.incoherent) = dots.x(dots.incoherent) + (freq_factor * cos(rdir));
end

offscreen = dots.x > dots.maxX;
dots.x(offscreen) = dots.x(offscreen) - abs(dots.maxX - dots.minX);
offscreen = dots.x < dots.minX;
dots.x(offscreen) = dots.x(offscreen) + abs(dots.maxX - dots.minX);

offscreen = dots.y > dots.maxY;
dots.y(offscreen) = dots.y(offscreen) - abs(dots.maxY - dots.minY);
offscreen = dots.y < dots.minY;
dots.y(offscreen) = dots.y(offscreen) + abs(dots.maxY - dots.minY);


function setGT(myscreen,stimulus)

% multipliers relative
max = 95;
% load the calibration
load(myscreen.calibFullFilename);

% get the max value from the calibration file
localMax = calib.tableCorrected.luminance(end);

% get an approximate factor
factor = round(localMax/max);

disp(sprintf('Setting gamma table to use %02.2f%% of available space.',1/factor*100));
% if not 1, set the gammaTable accordingly

if factor > 1
    newTable = interp1(linspace(0,1,256),stimulus.linearizedGammaTable.redTable,linspace(0,1/factor,256));
    succ = mglSetGammaTable(repmat(newTable',1,3));

    if ~succ
        warning('Gamma table set failure');
        keyboard
    end
end

function stimulus = initGratings(stimulus,myscreen)

stimulus.maxIndex = 255;
disppercent(-inf,'Creating grating textures');

nContrasts = length(stimulus.contrast);
nPhases = length(stimulus.grating.phases);

gaussianWin = mglMakeGaussian(stimulus.grating.width,stimulus.grating.height,stimulus.grating.sdx,stimulus.grating.sdy);
if strcmp(stimulus.grating.windowType,'gabor')
  % a gaussian window
  win = stimulus.maxIndex-stimulus.maxIndex*gaussianWin;
else
  % a simple window
  win = stimulus.maxIndex-stimulus.maxIndex*(gaussianWin>exp(-1/2));
end
mask = ones(size(win,1),size(win,2),4)*1;%myscreen.grayIndex;
win = win / max(win(:)) * 255;
mask(:,:,4) = round(win);
stimulus.mask = mglCreateTexture(mask);

% make each one of he called for gratings
for iPhase = 1:nPhases
  for iContrast = 1:nContrasts
    disppercent(calcPercentDone(iPhase,nPhases,iContrast,nContrasts));
    % get the phase and contast
    thisPhase = (stimulus.grating.phase+stimulus.grating.phases(iPhase))*180/pi;
    thisContrast = stimulus.contrastOverride;
    % make the grating
    thisGrating = round(stimulus.maxIndex*((thisContrast*mglMakeGrating(stimulus.grating.width,nan,stimulus.grating.sf,0,thisPhase))+1)/2);
    % create the texture
%     thisGrating = repmat(thisGrating,453,1);
    
    stimulus.tex(iContrast,iPhase) = mglCreateTexture(thisGrating);
  end
end
disppercent(inf);
% stimulus.randMaskSize = [size(mask,1) size(mask,2)];
% stimulus.randMask = mglCreateTexture(floor(stimulus.maxIndex*rand(stimulus.randMaskSize)));


for iAngle = 1:length(stimulus.orientations)
  % get center of patch
  thisAngle = stimulus.orientations(iAngle);
  centerX = stimulus.x + stimulus.grating.radius*cos(pi*thisAngle/180);
  centerY = stimulus.y + stimulus.grating.radius*sin(pi*thisAngle/180);
  % now get top and bottom point of grating
  thisOrientation = thisAngle+90;
  radius = sqrt((stimulus.grating.width/2).^2 +(stimulus.grating.height/2).^2)-0.5;
  topX = centerX + radius*cos(pi*thisOrientation/180);
  topY = centerY + radius*sin(pi*thisOrientation/180);
  thisOrientation = thisOrientation+180;
  bottomX = centerX + radius*cos(pi*thisOrientation/180);
  bottomY = centerY + radius*sin(pi*thisOrientation/180);
  % place points
  stimulus.grating.refPoints.x{iAngle} = [topX bottomX];
  stimulus.grating.refPoints.y{iAngle} = [topY bottomY];
end

stimulus.waitForBacktickText = mglText('Hit backtick (`) key to start');
