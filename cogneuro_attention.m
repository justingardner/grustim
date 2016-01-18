
% cohCon
%
%      usage: myscreen=cohcon()
%         by: daniel birman
%       date: 11/10/14
%    purpose: contrast change detection with cued selective attention.
%
%        use: call flowAwe() to initialize. The first
%             run takes significantly longer due to loading stimuli.
%
%      flags: stimFileNum (-1/#) - Load a specific stimfile from a
%             subject's folder. Defaults to the last stimfile. Warning:
%             Only the first stimfile is saved with the image file data,
%             subsequent stimfiles only contain run data.
%             unattended (0/1) - If 1, runs a fixation task while the main
%             task just runs on idle in the background (no inputs)
%             plots (0/1) - Displays staircase plots (and estimated
%             psychophysic functions)
%             overrideTask (1/2) - Specifies the task to run: 1 =
%             coherence, 2 = contrast
%             projector (0/1) - Masks stimuli using the default projector
%             mask.
%             scan (0/1) - Scanner timing
%
%
%   TR .75 = 560 volumes (7:00 total)
%   TR 1.4 = 300 volumes (7:00 total)

function [myscreen] = cohcon(varargin)

global stimulus
clear fixStimulus
%% Initialize Variables

% add arguments later
stimFileNum = [];
plots = [];
overrideTask = [];
projector = [];
scan = [];
training = [];
nocatch = [];
getArgs(varargin,{'stimFileNum=-1','nocatch=0', ...
    'plots=1','overrideTask=0','scan=0','training=0'});
stimulus.scan = scan;
stimulus.plots = plots;
stimulus.training = training;
stimulus.nocatch = nocatch;

if stimulus.scan && ~mglGetParam('ignoreInitialVols')==16 && ~mglGetParam('ignoreInitialVols')==4
    warning('ignoreInitialVols is set to %i.',mglGetParam('ignoreInitialVols'));
    if ~strcmp('y',input('Is this correct? [y/n]'))
        mglSetParam('ignoreInitialVols',input('Please input the correct value (mux8 = 16, mux2 = 4): '));
    end
end

stimulus.counter = 1; % This keeps track of what "run" we are on.

%% Useful stimulus stuff

stimulus.pedestals.pedOpts = {'coherence','contrast'};

if stimulus.nocatch && stimulus.scan
    stimulus.pedestals.coherence = [1];
    stimulus.pedestals.contrast = [1];
else
    stimulus.pedestals.coherence = 1;
    stimulus.pedestals.contrast = 1;
end

stimulus.pedestals.initThresh.coherence = 0;
stimulus.pedestals.initThresh.contrast = 0;

stimulus.pedestals.catch.coherence = [0 exp(-3:.33:-.3)];
stimulus.pedestals.catch.contrast = [0 exp(-4:.33:-1.3)];

%% Setup Screen

if stimulus.scan
    myscreen = initScreen('fMRIprojFlex');
else
    myscreen = initScreen('VPixx');
end

%% Open Old Stimfile
stimulus.initStair = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/cohcon/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/cohcon/%s/1*mat',mglGetSID));

    if length(files) >= 1
        if stimFileNum == -1
            if length(files) > 1
                warning('Multiple stimfiles found, loading last one. If you wanted a different functionality use stimFileNum=#');
            end
            fname = files(end).name;
        else
            fname = files(stimFileNum).name;
        end
        s = load(sprintf('~/data/cohcon/%s/%s',mglGetSID,fname));
        stimulus.staircase = s.stimulus.staircase;
        stimulus.stairCatch = s.stimulus.stairCatch;
        stimulus.nocatchs.staircase = s.stimulus.nocatchs.staircase;
        stimulus.counter = s.stimulus.counter + 1;

        % load blocks too
        stimulus.runs = s.stimulus.runs;
        stimulus.runs.loaded = 1;

        clear s;
        stimulus.initStair = 0;
        disp(sprintf('(cohCon) Data file: %s loaded.',fname));
    end
end
disp(sprintf('(cohCon) This is run #%i',stimulus.counter));

%% Initialize Stimulus

myscreen = initStimulus('stimulus',myscreen);

if stimulus.scan
    stimulus.responseKeys = [1 2]; % corresponds to LEFT - RIGHT
else
    stimulus.responseKeys = [1 2]; % corresponds to LEFT - RIGHT
end

%% Colors
stimulus.colors.rmed = 127.75;

% We're going to add an equal number of reserved colors to the top and
% bottom, to try to keep the center of the gamma table stable.
stimulus.colors.reservedBottom = [0 0 0; 1 1 1]; % fixation cross colors
stimulus.colors.reservedTop = [1 0 0; 0 1 0]; % correct/incorrect colors
stimulus.colors.black = 0/255; stimulus.colors.white = 1/255;
stimulus.colors.red = 254/255; stimulus.colors.green = 255/255;
stimulus.colors.nReserved = 2; % this is /2 the true number, because it's duplicated
stimulus.colors.nUnreserved = 256-(2*stimulus.colors.nReserved);

stimulus.colors.mrmax = stimulus.colors.nReserved - 1 + stimulus.colors.nUnreserved;
stimulus.colors.mrmin = stimulus.colors.nReserved;

%% Everything else
stimulus.dots.xcenter = 0;
stimulus.dots.ycenter = 0;
stimulus.dots.dotsize = 4;
stimulus.dots.density = 1.4;
stimulus.dots.speed = 3.25;
stimulus.dots.centerOffset = 2;

stimulus.dotsR = stimulus.dots;
stimulus.dotsR.mult = 1;
stimulus.dotsL = stimulus.dots;
stimulus.dotsL.mult = -1;
stimulus = rmfield(stimulus,'dots');

stimulus.dotsR = initDotsRadial(stimulus.dotsR,myscreen);
stimulus.dotsL = initDotsRadial(stimulus.dotsL,myscreen);

if ~length(stimulus.pedestals.coherence)==length(stimulus.pedestals.contrast)
    disp('(cohcon) Contrast and coherence have different pedestal #, this can cause errors');
end

%% Gamma Table Initialization

% get gamma table
if ~isfield(myscreen,'gammaTable')
  stimulus.linearizedGammaTable = mglGetGammaTable;
  disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
  disp(sprintf('(cuecon:initGratings) No gamma table found in myscreen. Contrast displays like this'));
  disp(sprintf('         should be run with a valid calibration made by moncalib for this monitor.'));
  disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
end
stimulus.linearizedGammaTable = myscreen.initScreenGammaTable;
    
%% Mask
d = stimulus.dotsL;
stimulus.mask.x = repmat([d.minX+.125:.25:d.maxX-.125],1,40);
stimulus.mask.y = [d.minY+.125:.25:d.maxY-.125];
tmp = repmat(stimulus.mask.y,30,1);
stimulus.mask.y = transpose(tmp(:));

%% Character textures
mglTextSet('Helvetica',32,stimulus.colors.white,0,0,0,0,0,0,0);
stimulus.text.mTexW = mglText('M');
stimulus.text.cTexW = mglText('C');

mglTextSet('Helvetica',32,stimulus.colors.black,0,0,0,0,0,0,0);
stimulus.text.mTexK = mglText('M');
stimulus.text.cTexK = mglText('C');
%% Setup Task

% This is the contrast change detection task
task{1}{1}.waitForBacktick = 1;

stimulus.seg.ITI = 5; % the ITI is either 20s (first time) or 1s
stimulus.seg.stim = 1;
stimulus.seg.mask = 2;
stimulus.seg.ISI = 3;
stimulus.seg.resp = 4;
task{1}{1}.segmin = [2 0 0 1.5 6];
task{1}{1}.segmax = [2 0 0 1.5 15];

task{1}{1}.synchToVol = [0 0 0 0 0];
if stimulus.scan
    task{1}{1}.synchToVol(stimulus.seg.ITI) = 1;
end
task{1}{1}.getResponse = [0 0 0 1 0];
task{1}{1}.parameter.dirL = [1 2];
task{1}{1}.parameter.dirR = [1 2];
task{1}{1}.random = 1;
task{1}{1}.numTrials = 25;

if stimulus.scan
    task{1}{1}.numTrials = inf;
end

if stimulus.nocatch
    task{1}{1}.parameter.catch = -1;
end

%% Tracking

% these are variables that we want to track for later analysis.
task{1}{1}.randVars.calculated.task = nan; % Current task (calc per BLOCK)
task{1}{1}.randVars.calculated.correct = nan;
task{1}{1}.randVars.calculated.trialNum = nan;

stimulus.curTrial = 0;

%% Block setup
if isfield(stimulus,'runs') && isfield(stimulus.runs,'loaded')
    % We already have our blocks
    stimulus.runs = rmfield(stimulus.runs,'loaded'); % remove the load field, otherwise it gets saved across runs
    if stimulus.counter > length(stimulus.runs.taskList)
        stimulus.runs.taskList  = repmat(stimulus.runs.taskList,1,2);
    end
else
    % This is the first run, build up the blocks.
    stimulus.runs = struct;
    stimulus.runs.taskOpts = [1 2];
    stimulus.runs.taskOptsText = {'Attend Left','Attend Right'};
    stimulus.runs.taskList = stimulus.runs.taskOpts(randperm(2));
end

%% Task Override
if overrideTask > 0
    stimulus.runs.curTask = overrideTask;
else
    stimulus.runs.curTask = stimulus.runs.taskList(stimulus.counter);
end

%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);
end

%% Get Ready...
% clear screen    
mglWaitSecs(2);
mglClearScreen(0.5);
if stimulus.scan        
    mglTextDraw('DO NOT MOVE',[0 1.5]);
    mglTextDraw(stimulus.runs.taskOptsText{stimulus.runs.curTask},[0 0]);
else

    mglTextDraw(stimulus.runs.taskOptsText{stimulus.runs.curTask},[0 0]);
end
mglFlush

% let the user know
disp(sprintf('(cohCon) Starting run number: %i. Current task: %s',stimulus.counter,stimulus.runs.taskOptsText{stimulus.runs.curTask}));
% if stimulus.unattended
myscreen.flushMode = 1;
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
mglTextDraw('Run complete... please wait.',[0 0]);
mglFlush
myscreen.flushMode = 1;

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

if stimulus.plots
    disp('(cohCon) Displaying plots');
    dispStaircase(stimulus);
end

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startTrialCallback(task,myscreen)
%%

global stimulus

stimulus.curTrial = stimulus.curTrial + 1;

%  Set the current task
task.thistrial.task = stimulus.runs.curTask;


% Set the missing thistrial vars
task.thistrial.trialNum = stimulus.curTrial;

% set directions
flip = [-1 1];
stimulus.dotsL.dir = flip(task.thistrial.dirL);
stimulus.dotsR.dir = flip(task.thistrial.dirR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)
%%

myscreen.flushMode = 0;

global stimulus

stimulus.live.mt = 0;
switch task.thistrial.thisseg
    case stimulus.seg.ITI
        stimulus.live.dots = 0;
        stimulus.live.fixColor = stimulus.colors.black;
    case stimulus.seg.stim
        stimulus.live.dots = 1;
        stimulus.live.fixColor = stimulus.colors.black;
    case stimulus.seg.mask
        stimulus.live.dots = 0;
        stimulus.live.fixColor = stimulus.colors.black;
    case stimulus.seg.ISI
        stimulus.live.dots = 0;
        stimulus.live.fixColor = stimulus.colors.black;
    case stimulus.seg.resp
        stimulus.live.dots = 0;
        stimulus.live.fixColor = stimulus.colors.white;
end

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus

mglClearScreen(0.5);

if stimulus.live.dots==1, stimulus = upDots(task,stimulus,myscreen); end
upFix(task,stimulus);

function upFix(task,stimulus)
%%

if ~(task.thistrial.catch > 0) || stimulus.live.catchFix == 0
    mglFixationCross(1.5,1.5,stimulus.live.fixColor);
else
    if stimulus.runs.curTask==2
        if stimulus.live.fixColor==stimulus.colors.white
            mglBltTexture(stimulus.text.mTexW,[0 0]);
        else
            mglBltTexture(stimulus.text.mTexK,[0 0]);
        end
    else
        if stimulus.live.fixColor==stimulus.colors.white
            mglBltTexture(stimulus.text.cTexW,[0 0]);
        else
            mglBltTexture(stimulus.text.cTexK,[0 0]);
        end
    end
end

function stimulus = upDots(task,stimulus,myscreen)

% update the dots

tCoh = 1;
tCon = 1;

%% Old update code start here
stimulus.dotsL = updateDotsRadial(stimulus.dotsL,tCoh,myscreen,false);
stimulus.dotsR = updateDotsRadial(stimulus.dotsR,tCoh,myscreen,false);

% dotsR
% update +contrast

mglPoints2(stimulus.dotsR.xdisp(stimulus.dotsR.con==1),stimulus.dotsR.ydisp(stimulus.dotsR.con==1),...
    stimulus.dotsR.dotsize,[0 0 0]);
% update - contrast
mglPoints2(stimulus.dotsR.xdisp(stimulus.dotsR.con==2),stimulus.dotsR.ydisp(stimulus.dotsR.con==2),...
    stimulus.dotsR.dotsize,[1 1 1]);
% dotsL
% update +contrast
mglPoints2(stimulus.dotsL.xdisp(stimulus.dotsL.con==1),stimulus.dotsL.ydisp(stimulus.dotsL.con==1),...
    stimulus.dotsL.dotsize,[0 0 0]);
% update - contrast
mglPoints2(stimulus.dotsL.xdisp(stimulus.dotsL.con==2),stimulus.dotsL.ydisp(stimulus.dotsL.con==2),...
    stimulus.dotsL.dotsize,[1 1 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus

responseText = {'Incorrect','Correct'};
responsePos = {'Left','Right'};
fixColors = {stimulus.colors.red,stimulus.colors.green};

if any(task.thistrial.whichButton == stimulus.responseKeys)
    if task.thistrial.gotResponse == 0
        if task.thistrial.task==1
            cSide = task.thistrial.dirL;
        else
            cSide = task.thistrial.dirR;
        end
        task.thistrial.correct = task.thistrial.whichButton == stimulus.responseKeys(cSide);
        % Store whether this was correct
        stimulus.live.fixColor = fixColors{task.thistrial.correct+1};
        disp(sprintf('(cohCon) Response %s: %s',responseText{task.thistrial.correct+1},responsePos{find(stimulus.responseKeys==task.thistrial.whichButton)}));
    else
        disp(sprintf('(cohCon) Subject responded multiple times: %i',task.thistrial.gotResponse+1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              HELPER FUNCTIONS                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
%    getDeltaPed       %
%%%%%%%%%%%%%%%%%%%%%%%%

function [cohPed, conPed, stimulus] = getDeltaPed(task,stimulus)
%%
if stimulus.runs.curTask == 1
    % COHERENCE MAIN TASK    
    if ~stimulus.nocatch
        [cohPed, stimulus.staircase{1,curPedVal(task,1)}] = doStaircase('testValue',stimulus.staircase{1,curPedVal(task,1)});
    else
        [cohPed, stimulus.nocatchs.staircase{1,curPedVal(task,1)}] = doStaircase('testValue',stimulus.nocatchs.staircase{1,curPedVal(task,1)});
    end
    
    if task.thistrial.catch > 0
        [conPed, stimulus.stairCatch{2,curPedVal(task,2)}] = doStaircase('testValue',stimulus.stairCatch{2,curPedVal(task,2)});
    else
        conPed = stimulus.pedestals.catch.contrast(randi(length(stimulus.pedestals.catch.contrast)));
    end
else
    % CONTRAST MAIN TASK
    if ~stimulus.nocatch
        [conPed, stimulus.staircase{2,curPedVal(task,2)}] = doStaircase('testValue',stimulus.staircase{2,curPedVal(task,2)});
    else
        [conPed, stimulus.nocatchs.staircase{2,curPedVal(task,2)}] = doStaircase('testValue',stimulus.nocatchs.staircase{2,curPedVal(task,2)});
    end
    
    if task.thistrial.catch > 0
        [cohPed, stimulus.stairCatch{1,curPedVal(task,1)}] = doStaircase('testValue',stimulus.stairCatch{1,curPedVal(task,1)});
    else
        cohPed = stimulus.pedestals.catch.coherence(randi(length(stimulus.pedestals.catch.coherence)));
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%
%    dispStaircase    %
%%%%%%%%%%%%%%%%%%%%%%%
function dispStaircase(stimulus)

try
    taskOpts = {'catch - coherence','catch - contrast','coherence','contrast','Nocatch coherence','Nocatch contrast'};
    drawing = {'or' 'ob' '*r' '*b' '+r' '+b'};
    
    
    plotting = zeros(2,1);
    catchPlot = zeros(2,1);
    ci = zeros(2,1);
    nocatchplot = zeros(2,1);
    
    for task = 1:2
        for ped = 1:length(stimulus.pedestals.contrast)
            try
                each = [];
                for i = 1:length(stimulus.nocatchs.staircase{task,ped})
                    if stimulus.nocatchs.staircase{task,ped}(i).trialNum > 0
                        out = doStaircase('threshold',stimulus.nocatchs.staircase{task,ped}(i),'type','weibull'); % noise
                        each(end+1) = out.threshold;
                    end
                end
                nocatchplot(task,ped) = mean(each);
            catch
                nocatchplot(task,ped) = -1;
            end
            try
                each = [];
                for i = 1:length(stimulus.staircase{task,ped})
                    if stimulus.staircase{task,ped}(i).trialNum > 0
                        out = doStaircase('threshold',stimulus.staircase{task,ped}(i),'type','weibull'); % noise
                        each(end+1) = out.threshold;
                    end
                end
                plotting(task,ped) = mean(each);
                ci(task,ped) = std(each)/sqrt(length(each))*1.96;
            catch
                plotting(task,ped) = -1;
            end
            try
                outC = doStaircase('threshold',stimulus.stairCatch{task,ped},'type','weibull');
                catchPlot(task,ped) = outC.threshold;
            catch
                catchPlot(task,ped) = -1;
            end
        end
    end
    figure
    hold on
    plot(stimulus.pedestals.(taskOpts{3})(1),catchPlot(1,:),drawing{1});
    plot(stimulus.pedestals.(taskOpts{4})(1),catchPlot(2,:),drawing{2});
    plot(stimulus.pedestals.(taskOpts{3})(1),plotting(1,:),drawing{3});
    plot(stimulus.pedestals.(taskOpts{4})(1),plotting(2,:),drawing{4});
    plot(stimulus.pedestals.(taskOpts{3})(1),nocatchplot(1,:),drawing{5});
    plot(stimulus.pedestals.(taskOpts{4})(1),nocatchplot(2,:),drawing{6});
    legend(taskOpts);
    a = axis;
    axis([0 .7 a(3) a(4)]);
    hold off

catch
    disp('(cohCon) Figures were not generated successfully.');
end

%% checkStaircaseStop
function checkStaircaseStop()
global stimulus

for task = 1:2
    stimulus.stairCatch{task,1} = resetStair(stimulus.stairCatch{task,1});
end
% Check both staircases
for task = 1:2
    stimulus.staircase{task,1} = resetStair(stimulus.staircase{task,1});
end
% Check nocatch staircases

for task = 1:2
    if stimulus.scan
        for p = 1:size(stimulus.nocatchs.staircase,2)
            stimulus.nocatchs.staircase{task,p} = resetStair(stimulus.nocatchs.staircase{task,p});
        end
    else
        stimulus.nocatchs.staircase{task,1} = resetStair(stimulus.nocatchs.staircase{task,1});
    end
end

function s = resetStair(s)

if doStaircase('stop',s)
    % this is a bit of a pain... you can't pass an initialThreshold
    % argument do doStaircase('init',s, ...), it ignores everything and
    % resets using the calculated threshold. Because you can't override it
    [args, vals, ~] = getArgs(s(1).initArgs);
    threshPos = -1;
    stepPos = -1;
    for i = 1:length(args)
        switch args{i}
            case 'initialThreshold'
                threshPos = i;
            case 'initialStepsize'
                stepPos = i;
        end
    end
    out = doStaircase('threshold',s);
    in = input(sprintf('Resetting Staircase... Estimate is: %1.2f. Reset ([Y]/[C]ustom/[O]riginal): ',out.threshold),'s');
    switch in
        case 'Y'
            vals{threshPos} = out.threshold;
            vals{stepPos} = out.threshold / 3;
        case 'C'
            disp('Original values:');
            disp(sprintf('%s: %0.2f',args{threshPos},num2str(vals{threshPos})));
            val = str2double(input('New threshold value: ','s'));
            vals{threshPos} = val;
            vals{stepPos} = val / 3;
        case 'O'
    end
    if ~length(args) == 8
        disp('Args incorrect length...');
        keyboard
    end
    %             stimulus.staircase{task,ped}(end+1) = doStaircase('init',s,'initialThreshold',vals{threshPos},'initialStepsize',vals{stepPos});
    s(end+1) = doStaircase('init','upDown',args{1},vals{1},args{2},vals{2},args{3},vals{3},args{4},vals{4},args{5},vals{5},args{6},vals{6},args{7},vals{7},args{8},vals{8});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create dots for optic flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initDotsRadial(dots,~)

% maximum depth of points
dots.minX = 3.5;
dots.maxX = 11;
dots.minY = -5;
dots.maxY = 5;

dots.dir = 1;

% make a some points
dots.n = 500*dots.density;
% make sure it's an even number
dots.n = dots.n + mod(dots.n,2);

% set half to white and half to black
dots.con = repmat([1 2],1,dots.n/2);

dots.x = rand(1,dots.n)*(dots.maxX-dots.minX)+dots.minX;
dots.y = rand(1,dots.n)*abs(dots.maxY-dots.minY)+dots.minY;

dots.xdisp = dots.mult*dots.x;
dots.ydisp = dots.y;

% set incoherent dots to 0
dots.coherency = 1;
dots.incoherent = rand(1,dots.n) > dots.coherency;
dots.incoherentn = sum(dots.incoherent);
dots.coherent = ~dots.incoherent;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots for Radial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updateDotsRadial(dots,coherence,myscreen,repick)

% stuff to compute median speed
dots.oldx = dots.x;
dots.oldy = dots.y;

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

freq_factor = dots.speed/myscreen.framesPerSecond;

% move coherent dots
dots.x(dots.coherent) = dots.x(dots.coherent) + dots.dir*freq_factor;

% these are for flipping into the other quadrants
xmat = repmat([1 1 -1 -1],1,dots.incoherentn+4-mod(dots.incoherentn,4));
ymat = repmat([1 -1 1 -1],1,dots.incoherentn+4-mod(dots.incoherentn,4));
perms = randperm(dots.incoherentn);

% move incoherent dots
% get random vectors
dots.rX = rand(1,dots.incoherentn);
dots.rY = sqrt(1-(dots.rX.^2));
dots.rX = (dots.rX .* xmat(perms)) * freq_factor; % rescale to match the velocity
dots.rY = (dots.rY .* ymat(perms)) * freq_factor;
% dots.rX = (dots.rX .* xmat) * freq_factor .* ((1.75*rand(1,dots.incoherentn)).^2); % rescale to match the velocity
% dots.rY = (dots.rY .* ymat) * freq_factor .* ((1.75*rand(1,dots.incoherentn)).^2);
dots.x(dots.incoherent) = dots.x(dots.incoherent) + dots.rX;
dots.y(dots.incoherent) = dots.y(dots.incoherent) + dots.rY;

offscreen = dots.x > dots.maxX;
dots.x(offscreen) = dots.x(offscreen) - abs(dots.maxX - dots.minX);
offscreen = dots.x < dots.minX;
dots.x(offscreen) = dots.x(offscreen) + abs(dots.maxX - dots.minX);

offscreen = dots.y > dots.maxY;
dots.y(offscreen) = dots.y(offscreen) - abs(dots.maxY - dots.minY);
offscreen = dots.y < dots.minY;
dots.y(offscreen) = dots.y(offscreen) + abs(dots.maxY - dots.minY);

dots.xdisp = dots.mult*dots.x;
dots.ydisp = dots.y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sets the gamma table so that we can have
% finest possible control over the stimulus contrast.
%
% stimulus.colors.reservedColors should be set to the reserved colors (for cue colors, etc).
% maxContrast is the maximum contrast you want to be able to display.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setGammaTable_flowMax(maxContrast)

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
    redLinearized = repmat(.5,1,stimulus.colors.nUnreserved);
    greenLinearized = repmat(.5,1,stimulus.colors.nUnreserved);
    blueLinearized = repmat(.5,1,stimulus.colors.nUnreserved);
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
mglSetGammaTable(gammaTable);

% remember what the current maximum contrast is that we can display
stimulus.curMaxContrast = maxContrast;