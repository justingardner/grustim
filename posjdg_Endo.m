function [ myscreen ] = posjdg( varargin )
%POSITIONJUDGMENTS 
%
% Position judgment task with three attentional conditions. 

global stimulus

stimulus = struct;
%% Open Old Stimfile
stimulus.counter = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/posjdg/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/posjdg/%s/1*mat',mglGetSID));

    if length(files) >= 1
        fname = files(end).name;
        
        s = load(sprintf('~/data/posjdg/%s/%s',mglGetSID,fname));
        % copy staircases and run numbers
        stimulus.counter = s.stimulus.counter + 1;
        stimulus.staircase = s.stimulus.staircase;
        clear s;
        disp(sprintf('(posjdg) Data file: %s loaded.',fname));
    end
end
disp(sprintf('(posjdg) This is run #%i',stimulus.counter));


%% Initialize Variables

stimulus.attentionModes = {'Prior','Exo','Saccade'};
% add arguments later
scan = 0;
plots = 0;
noeye = 0;
debug = 0;
noimp = 0;
training = 0; attmode= 0;
getArgs(varargin,{'scan=0','attmode=1','plots=0','noeye=0','debug=0','training=0'});
stimulus.scan = scan;
stimulus.attentionMode = attmode;
stimulus.plots = plots;
stimulus.noeye = noeye;
stimulus.debug = debug;
stimulus.noimp = noimp;
stimulus.training = training;
clear localizer invisible scan noeye task

if stimulus.scan
    warning('Not setup for scanning');
end

disp(sprintf('Attention mode: %s',stimulus.attentionModes{stimulus.attentionMode}));

if any(stimulus.attentionMode==[2 3])
    keyboard
    disp('Not implemented');
end
%% Setup Screen

% if stimulus.scan
%     myscreen = initScreen('fMRIprojFlex');
% else
    myscreen = initScreen('VPixx');
% end

% set background to grey
myscreen.background = 0.5;

%% Staircase
if ~isfield(stimulus,'staircase')
    disp('(posjdg) WARNING: New staircase');
    stimulus.staircase = initStair();
else
    resetStair();
end

%% Setup missing initial variables

if ~isfield(stimulus,'counter')
    stimulus.counter = 1; % This keeps track of what "run" we are on.
end

%% Plot and return

%% Initialize Stimulus

myscreen.stimulusNames{1} = 'stimulus';

localInitStimulus();
    
if stimulus.scan
    stimulus.responseKeys = [2 1]; % corresponds to NOMATCH, MATCH
else
    stimulus.responseKeys = [2 1]; % corresponds to  NOMATCH, MATCH
end

stimulus.colors.black = [0 0 0];
stimulus.colors.white = [1 1 1];
stimulus.colors.green = [0 1 0];
stimulus.colors.red = [1 0 0];

stimulus.colors.valid = [255,143,143]/255;
stimulus.colors.impossible = [94,161,204]/255;
stimulus.colors.chance = [255,254,168]/255;

% % %% Generate stencils
% % % The stencil is a series of arcs 
% % mglStencilCreateBegin(1);
% % % Draw an annulus at every buffer location
% % for i = 0:(stimulus.cur_.num-1)
% %     partialDiskFuckOGL(0,0,stimulus.cur_.isize,stimulus.cur_.osize,i*stimulus.cur_.angle+stimulus.cur_.buffer/2,stimulus.cur_.angle-stimulus.cur_.buffer,[1 1 1],60,2);
% % end
% % mglStencilCreateEnd;
% % mglClearScreen(0.5);
% % myscreen.flushMode = 1;

if stimulus.plots==2
    dispInfo(stimulus);
    return
end

%% Setup Task
task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;

stimulus.curTrial = 0;

if stimulus.attentionMode==1
    % task waits for fixation on first segment
    task{1}{1}.segmin = [inf .200 0.200 1.000 0.250];
    task{1}{1}.segmax = [inf .200 0.800 1.000 0.750];

    stimulus.seg.ITI1 = 1; % waits for user input (button press + held) and eye fixation (within 2 degrees)
    stimulus.seg.stim = 2;
    stimulus.seg.delay = 3;
    stimulus.seg.resp = 4;
    stimulus.seg.ITI2 = 5;
else
    keyboard;
end

if stimulus.noeye==1
    task{1}{1}.segmin(1) = 0.5;
    task{1}{1}.segmax(1) = 0.5;
end

task{1}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
task{1}{1}.getResponse = zeros(size(task{1}{1}.segmin)); task{1}{1}.getResponse(stimulus.seg.resp)=1;
task{1}{1}.numTrials = 60;
task{1}{1}.random = 1;
task{1}{1}.parameter.ecc = 6; % eccentricity of display
if stimulus.attentionMode==1
    task{1}{1}.parameter.target = 6; % prior center (used as cue or saccade target in aM==2/3)
    task{1}{1}.parameter.priorSTD = 0.3; % radians
else
    keyboard
end
% task{1}{1}.parameter.flip = 1; % DO NOT USE flips from right angles to left angles

if stimulus.scan
    task{1}{1}.synchToVol(stimulus.seg.ITI) = 1;
end

%% Tracking

% these are variables that we want to track for later analysis.
task{1}{1}.randVars.calculated.correct = nan;
task{1}{1}.randVars.calculated.angle = nan; % angle at which displayed, depends on attention mode
task{1}{1}.randVars.calculated.rotation = nan; % rotation of the grating

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
disp(sprintf('(posjdg) Starting run number: %i.',stimulus.counter));

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
mglTextDraw('Please wait',[0 0]);
mglFlush
myscreen.flushMode = 1;

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

mglClearScreen(0.5);

if stimulus.plots
    disp('(posjdg) Displaying plots');
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

dopts = {'Right','Left'};
vopts = {'Valid','Impossible'};
disp(sprintf('(posjdg) %s:%s trial (%i), %s. Angle 1: %3.0f Angle 2: %3.0f Pattern A: %i, pattern B: %i',...
    opts{task.thistrial.match+1},dopts{task.thistrial.match+1},task.trialnum,vopts{task.thistrial.impossible+1},...
    task.thistrial.angle1*180/pi,task.thistrial.angle2*180/pi,...
    task.thistrial.pattern1,task.thistrial.pattern2));
    
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

stimulus.live.eyeDead =0 ;
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

for i = 1:2
    mglClearScreen(0.5);
    if stimulus.live.stim
        if task.thistrial.thisseg==stimulus.seg.stim1
            data = stimulus.live.data1;
            orient = stimulus.live.angles1;
        else
            data = stimulus.live.data2;
            orient = stimulus.live.angles2;
        end

        % draw background on debug
        if stimulus.debug
    %         partialDiskFuckOGL(0,0,stimulus.cur_.isize-0.5,stimulus.cur_.osize+3,(stimulus.learn-1)*stimulus.cur_.angle,stimulus.cur_.angle,[163 93 93]/255,6,2);
        end

        % draw rings
        upData(data,orient,stimulus);
        % revert stencil

    %     if stimulus.debug
    %         mglTextSet([],32,stimulus.colors.white);
    %         for si = 0:(stimulus.cur_.num-1)
    %             mglTextDraw(num2str(si+1),[(stimulus.cur_.osize+1)*cos(deg2rad(si*stimulus.cur_.angle+stimulus.cur_.angle/2)) (stimulus.cur_.osize+1)*sin(deg2rad(si*stimulus.cur_.angle+stimulus.cur_.angle/2))]);
    %         end
    %     end

        % draw V or I for valid/invalid trials
        if stimulus.debug
            text = {'V','I'};
            mglTextSet([],32,stimulus.colors.white);
            mglTextDraw(text{task.thistrial.impossible+1},[-7.5 -7.5]);
        end
    end
    
    if stimulus.live.fix
    %      cover
            mglFillOval(0,0,[1 1],0.5);
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
if stimulus.dead && mglGetSecs(task.thistrial.segStartSeconds)>1
    jumpSegment(task,inf); stimulus.dead=0;
end

% skip screen updates if you are already dead
if stimulus.dead
    if stimulus.dead && stimulus.live.eyeDead
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
if ~stimulus.noeye && ~any(task.thistrial.thisseg==[stimulus.seg.ITI1 stimulus.seg.ITI2 stimulus.seg.resp]) && ~stimulus.scan
    if ~any(isnan(pos))
        if dist > 2.5 && stimulus.live.eyeCount > 30
            disp('Eye movement detected!!!!');
            stimulus.dead = 1;
            stimulus.live.eyeDead=1;
            return
        elseif dist > 2.5
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

function upData(data,orient,stimulus)

xp = stimulus.cur_.pos([1:floor(length(stimulus.cur_.pos)/2) (length(stimulus.cur_.pos)-(floor(length(stimulus.cur_.pos)/2)-1)):length(stimulus.cur_.pos)]);
yp = fliplr(xp);

for x = 1:size(data,1)
    for y = 1:size(data,2)
        if data(x,y)>0
            mglBltTexture(stimulus.live.gratings{data(x,y)},[xp(y) yp(x)],0,0,orient(x,y)*180/pi);
        end
    end
end

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

if stimulus.dead, return; end
responseText = {'Incorrect','Correct'};
respText = {'-1','+5'};
sideText = {'Left','Right'};
matchText = {'Non-match','Match'};
fixColors = {stimulus.colors.red,stimulus.colors.green};
    
if any(task.thistrial.whichButton == stimulus.responseKeys)
    if task.thistrial.gotResponse == 0
        task.thistrial.correct = task.thistrial.whichButton == stimulus.responseKeys(task.thistrial.match+1);
        if ~task.thistrial.impossible && ~isfield(stimulus,'outThreshold')
            stimulus.staircase = doStaircase('update',stimulus.staircase,task.thistrial.correct);
        end
        disp(sprintf('Subject pressed %i/%s: %s %s',task.thistrial.whichButton,sideText{task.thistrial.whichButton},matchText{stimulus.responseKeys(task.thistrial.whichButton)},responseText{task.thistrial.correct+1}));
        stimulus.live.fixColor = fixColors{task.thistrial.correct+1};
        stimulus.live.resp = 1;
        stimulus.live.fix = 1;
        stimulus.live.respText = respText{task.thistrial.correct+1};
        for i = 1:2
            mglTextSet([],32,stimulus.live.fixColor);
            mglTextDraw(stimulus.live.respText,[0 0]);
            mglFlush
        end
    else
        disp(sprintf('(posjdg) Subject responded multiple times: %i',task.thistrial.gotResponse+1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              HELPER FUNCTIONS                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function staircase = initStair()

staircase = cell(1,3);
for i = 1:3
    staircase = doStaircase('init','upDown',...
                'initialThreshold',0.10,...
                'initialStepsize',0.025,...
                'minThreshold=0.0001','maxThreshold=0.4','stepRule','pest',...
                'nTrials=50','maxStepsize=0.2','minStepsize=0.0001');
end
        
function resetStair()
global stimulus

for i = 1:3
    if doStaircase('stop',stimulus.staircase{i})
        disp('(posjdg) Staircase is being reset');
        stimulus.staircase{i}(end+1) = doStaircase('init',stimulus.staircase{i}(end));
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
disp('dispInfo not implemented');

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
grating = 255/2*mglMakeGrating(sz,sz,6,0) + 255/2;
gauss = mglMakeGaussian(sz,sz,sz/6,sz/6);
alphamask = repmat(grating,1,1,4);
alphamask(:,:,4) = gauss*255;

% we'll adjust the gamma table to control contrast
stimulus.live.grating  = mglCreateTexture(alphamask); % high contrast
