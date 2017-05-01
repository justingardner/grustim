function [ myscreen ] = posjdg_Endo( varargin )
%POSITIONJUDGMENTS 
%
% Position judgment task with three attentional conditions. 

global stimulus

stimulus = struct;

%% Stimulus parameters 

% new prior for each run
stimulus.prior = rand*2*pi;
stimulus.sd = pi/4;
stimulus.ecc = 6;

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
clear localizer invisible scan noeye task

if stimulus.scan
    warning('Not setup for scanning');
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

myscreen.stimulusNames{1} = 'stimulus';

localInitStimulus();
    
if stimulus.scan
    stimulus.responseKeys = 3; % corresponds to NOMATCH, MATCH
else
    stimulus.responseKeys = 3; % corresponds to  NOMATCH, MATCH
end

%% Colors
initGammaTable(myscreen);
stimulus.colors.rmed = 127.5;

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

%% Setup Task

%%%%%%%%%%%%% PHASE ONE %%%%%%%%%%%%%%%%%
%%%%% PRIOR + ESTIMATE OF THRESHOLD %%%%%

stimulus.curTrial(1) = 0;

task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;

% task waits for fixation on first segment
task{1}{1}.segmin = [inf 0.000 .200 0.500];
task{1}{1}.segmax = [inf 2.000 .200 0.500];

stimulus.seg = {};

stimulus.seg{1}.ITI1 = 1; % waits for user input (button press + held) and eye fixation (within 2 degrees)
stimulus.seg{1}.delay = 2;
stimulus.seg{1}.stim = 3;
stimulus.seg{1}.resp = 4;

if stimulus.noeye==1
    task{1}{1}.segmin(1) = 0.5;
    task{1}{1}.segmax(1) = 0.5;
end

task{1}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
task{1}{1}.getResponse = zeros(size(task{1}{1}.segmin)); task{1}{1}.getResponse(stimulus.seg{1}.resp)=1;
task{1}{1}.numTrials = 30;
task{1}{1}.random = 1;
task{1}{1}.parameter.ecc = stimulus.ecc; % eccentricity of display
task{1}{1}.parameter.target = stimulus.prior; % prior center
task{1}{1}.parameter.priorSTD = stimulus.sd; % radians

if stimulus.scan
    task{1}{1}.synchToVol(stimulus.seg.ITI) = 1;
end

task{1}{1}.randVars.calculated.respAngle = nan;
task{1}{1}.randVars.calculated.angle = nan; % angle at which displayed, depends on attention mode
task{1}{1}.randVars.calculated.rotation = nan; % rotation of the grating
task{1}{1}.randVars.calculated.contrast = nan; % contrast of the grating

%%%%%%%%%%%%% PHASE TWO %%%%%%%%%%%%%%%%%
%%%%% PRIOR + ESTIMATE OF THRESHOLD %%%%%

stimulus.curTrial(2) = 0;

task{1}{2} = struct;
task{1}{2}.waitForBacktick = 1;

% task waits for fixation on first segment
task{1}{2}.segmin = [inf .200 0.000 inf];
task{1}{2}.segmax = [inf .200 2.000 inf];

stimulus.seg{2}.ITI1 = 1; % waits for user input (button press + held) and eye fixation (within 2 degrees)
stimulus.seg{2}.stim = 2;
stimulus.seg{2}.delay = 3;
stimulus.seg{2}.resp = 4;

if stimulus.noeye==1
    task{1}{1}.segmin(1) = 0.5;
    task{1}{1}.segmax(1) = 0.5;
end

task{1}{2}.synchToVol = zeros(size(task{1}{1}.segmin));
task{1}{2}.getResponse = zeros(size(task{1}{1}.segmin)); task{1}{1}.getResponse(stimulus.seg{2}.resp)=1;
task{1}{2}.numTrials = 100;
task{1}{2}.random = 1;
task{1}{2}.parameter.ecc = stimulus.ecc; % eccentricity of display
task{1}{2}.parameter.target = stimulus.prior; % prior center
task{1}{2}.parameter.priorSTD = stimulus.sd; % radians

if stimulus.scan
    task{1}{2}.synchToVol(stimulus.seg.ITI) = 1;
end

%% Tracking

% these are variables that we want to track for later analysis.
task{1}{2}.randVars.calculated.respAngle = nan;
task{1}{2}.randVars.calculated.angle = nan; % angle at which displayed, depends on attention mode
task{1}{2}.randVars.calculated.rotation = nan; % rotation of the grating
task{1}{2}.randVars.calculated.contrast = nan; % contrast of the grating

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
disp(sprintf('(posjdg) Starting run number: %i.',stimulus.counter));

%% Main Task Loop

setGammaTable(1);
mglClearScreen(0.5);
mglFlush
mglClearScreen(0.5);

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

if (task.thistrial.thisphase==1) && (~isempty(task.lasttrial))
    stimulus.staircase = doStaircase('update',stimulus.staircase,task.lasttrial.correct);
    disp(sprintf('Subject did not see %02.0f%% contrast',task.thistrial.contrast*100));
end

stimulus.curTrial(task.thistrial.thisphase) = stimulus.curTrial(task.thistrial.thisphase) + 1;

% compute missing variables
task.thistrial.angle = randn*task.thistrial.priorSTD+task.thistrial.target;
task.thistrial.rotation = rand*2*pi;

% contrast from staircase
[task.thistrial.contrast, stimulus.staircase] = doStaircase('testValue',stimulus.staircase);

myscreen.flushMode = 0;

disp(sprintf('(posjdg) Trial (%i): rotation: %02.0f, angle: %02.0f, contrast: %02.0f%%',...
    task.trialnum,task.thistrial.rotation,task.thistrial.angle,task.thistrial.contrast*100));
    
stimulus.live.eyeCount = 0;
stimulus.dead = 0;

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

stimulus.live.eyeDead =0 ;
stimulus.live.resp = 0;
stimulus.live.fixColor = stimulus.colors.white;
stimulus.live.fix = 1;
stimulus.live.stim = 0;

if task.thistrial.thisseg==stimulus.seg{task.thistrial.thisphase}.stim
    stimulus.live.stim = 1;
elseif task.thistrial.thisseg==stimulus.seg{task.thistrial.thisphase}.resp
    setGammaTable(1);
    stimulus.live.resp = 1;
    stimulus.live.angle = 0*2*pi;
    convertRespXY();
end
    
if stimulus.live.stim
    setGammaTable(task.thistrial.contrast);
end

for i = 1:2
    mglClearScreen(0.5);
    if stimulus.live.stim
        x = stimulus.ecc * acos(task.thistrial.angle);
        y = stimulus.ecc * asin(task.thistrial.angle);
        mglBltTexture(stimulus.live.grating,[x y],0,0,task.thistrial.rotation*180/pi);
    end
    
    % resp is updated in screenUpdate
    
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
if ~stimulus.noeye && ~any(task.thistrial.thisseg==[stimulus.seg{task.thistrial.thisphase}.ITI1 stimulus.seg{task.thistrial.thisphase}.resp]) && ~stimulus.scan
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

if (task.thistrial.thisseg==stimulus.seg{task.thistrial.thisphase}.resp) && 0
    
elseif task.thistrial.thisseg==stimulus.seg{task.thistrial.thisphase}.resp
    keys = find(mglGetKeys);
    if any(keys==19)
        stimulus.live.angle = stimulus.live.angle+0.02;
    elseif any(keys==20)
        stimulus.live.angle = stimulus.live.angle-0.02;
    end
    convertRespXY();
end

if stimulus.live.resp && task.thistrial.thisphase==2
    mglClearScreen(0.5);
    mglFillOval(stimulus.live.respx,stimulus.live.respy,[0.5 0.5],stimulus.colors.white);
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

function convertRespXY()
global stimulus

stimulus.live.respx = stimulus.ecc*cos(stimulus.live.angle);
stimulus.live.respy = stimulus.ecc*sin(stimulus.live.angle);

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
sideText = {'Left','Right'};
matchText = {'Non-match','Match'};
    
if any(task.thistrial.whichButton == stimulus.responseKeys)
    if task.thistrial.gotResponse == 0
        if task.thistrial.thisphase==1
            % they saw it, sot hey got it right
            task.thistrial.correct = 1;
            stimulus.staircase = doStaircase('update',stimulus.staircase,task.thistrial.correct);
            disp(sprintf('Subject saw %02.0f%% contrast',task.thistrial.contrast*100));
        elseif task.thistrial.thisphase==2
            % they are actually reporting locations
            task.thistrial.respAngle = stimulus.live.angle;
            disp(sprintf('Subject pressed %i/%s: %s %s',task.thistrial.whichButton,sideText{task.thistrial.whichButton},matchText{stimulus.responseKeys(task.thistrial.whichButton)},responseText{task.thistrial.correct+1}));
            stimulus.live.fix = 0;
        end
    else
        disp(sprintf('(posjdg) Subject responded multiple times: %i',task.thistrial.gotResponse+1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              HELPER FUNCTIONS                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function initStair()
global stimulus

stimulus.staircase = doStaircase('init','upDown',...
            'initialThreshold',0.5,...
            'initialStepsize',0.025,...
            'minThreshold=0.0001','maxThreshold=0.4','stepRule','pest',...
            'nTrials=25','maxStepsize=0.2','minStepsize=0.0001');
        
function resetStair()

global stimulus

if doStaircase('stop',stimulus.staircase)
    disp('(posjdg) Staircase is being reset');
    stimulus.staircase(end+1) = doStaircase('init',stimulus.staircase(end));
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