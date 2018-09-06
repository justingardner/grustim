function [ myscreen ] = afcom( varargin )
% ***scan info: 736 TRs (6 minute * 120 + 16)***
% *** set 'noeye=1' to turn of the eye tracker***
%
%Attention field and color mapping (afcom)
%
%   Map the effects of spatial attention and attention to specific colors 

%% Check version

if matlabVersionNumber<8.4
    error('rgb2lab and lab2rgb are not supported in this version of MATLAB -- upgrade to 2014b');
end

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

getArgs(varargin,{'scan=0','plots=0','noeye=0','powerwheel=0','eyewindow=1.5','practice=0','debug=0','replay=0','attend=1','run=0','build=0','mouse=0'});
stimulus.scan = scan;
stimulus.plots = plots;
stimulus.noeye = noeye;
stimulus.practice = practice;
stimulus.mousedebug = mouse;
stimulus.powerwheel = powerwheel;
stimulus.eyewindow = eyewindow;
stimulus.debug = debug;
stimulus.replay = replay;
stimulus.overrideRun = run;

clear localizer invisible scan noeye task test2 attend build eyewindow mouse practice powerwheel

if stimulus.scan
    warning('Disabling eyewindow');
    stimulus.eyewindow=0;
end

%% Replay mode
% if any(replay>0)
%     if ischar(replay)
%         % a file was called for, load it
%         loaded = load(replay);
%         stimulus = loaded.stimulus;
%         % check that this is actually an afmap file
%         if isempty(strfind(loaded.task{1}{1}.taskFilename,'afmap'))
%             disp(sprintf('File %s is not an afmap run.',replay));
%             return
%         end
%         % get the task parameters so that you can sync the replay correctly
%         disp('GET TASK PARAMETERS');
%         e = getTaskParameters(loaded.myscreen,loaded.task);
%         e1 = e{1};
%         % pull out the trial volumes
%         stimulus.tVolumes = e{1}.trialVolume;
%         disp('Stimulus volumes were found at:');
%         disp(stimulus.tVolumes);
%         
%         stimulus.replayFile = strcat(replay(1:(strfind(replay,'.mat')-1)),'_replay.mat');
%         stimulus.replay = true;
%     else
%         % do an entire subject
%         disp('Replay mode initiated');
%         folder = input('What folder would you like to replay? [Folder]: ');
%         files = dir(fullfile(folder,'*.mat'));
%         for fi = 1:length(files)
%             if isempty(strfind(files(fi).name,'replay')) && isempty(strfind(files(fi).name,'original'))
%                 afmap(sprintf('replay=%s/%s',folder,files(fi).name));
%             end
%         end
%         return;
%     end
% end

%% Open Old Stimfile
% if ~stimulus.replay
%     stimulus.counter = 1;
%     
%     if ~isempty(mglGetSID) && isdir(sprintf('~/data/afcom/%s',mglGetSID))
%         % Directory exists, check for a stimefile
%         files = dir(sprintf('~/data/afcom/%s/1*mat',mglGetSID));
%         
%         if length(files) >= 1
%             fname = files(end).name;
%             
%             s = load(sprintf('~/data/afcom/%s/%s',mglGetSID,fname));
%             % copy staircases and run numbers
%             stimulus.counter = s.stimulus.counter + 1;
%             stimulus.live.attend = mod(s.stimulus.live.attend+1,3);
%             if s.stimulus.attend ~= stimulus.attend
%                 error('(afcom) Cannot continue: stimfile parameters were generated with a different attention mode than you requested. You need to save the existing stimfiles into a backup folder');
%             end
%             clear s;
%             disp(sprintf('(afcom) Data file: %s loaded.',fname));
%         else
%             warning('(afcom) Unable to load previous data files. If this is *not* the first run there is something wrong.');
%         end
%     end
% end

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
    dispInfo(stimulus);
    return
end

%% Initialize Stimulus

if ~stimulus.replay
    myscreen.stimulusNames{1} = 'stimulus';
    
    stimulus.responseKeys = [1 2]; % left right
else
    localInitStimulus();
end

%% Colors
stimulus.colors.white = [1 1 1]; stimulus.colors.red = [1 0 0];
stimulus.colors.green = [0 1 0]; stimulus.colors.black = [0 0 0];

% get the lab space rgb values
stimulus.backgroundLab = rgb2lab([0.5 0.5 0.5]);

% setup color wheel
stimulus.colorwheel.acenter = stimulus.backgroundLab(2);
stimulus.colorwheel.bcenter = stimulus.backgroundLab(3);

% compute the ranges around 0 and pi for the colorwheel
theta_ = pi/64; % increment size
stimulus.colorwheel.thetaInc = theta_;
stimulus.colorwheel.thetaRange1 = 0 + (-pi/4:theta_:pi/4);
stimulus.colorwheel.thetaRange2 = pi + (-pi/4:theta_:pi/4);
stimulus.colorwheel.thetas = mod([stimulus.colorwheel.thetaRange1 stimulus.colorwheel.thetaRange2],2*pi);
% stimulus.colorwheel.dispThetas = 
stimulus.colorwheel.range1 = logical([ones(size(stimulus.colorwheel.thetaRange1)) zeros(size(stimulus.colorwheel.thetaRange2))]);
stimulus.colorwheel.range2 = ~stimulus.colorwheel.range1;

D = 70;
stimulus.colorwheel.distanceLab = D;

for ti = 1:length(stimulus.colorwheel.thetas)
    theta = stimulus.colorwheel.thetas(ti);
    
    a = D*cos(theta)+stimulus.colorwheel.acenter;
    b = D*sin(theta)+stimulus.colorwheel.bcenter;
    
    rgb = lab2rgb([stimulus.backgroundLab(1) a b],'ColorSpace','adobe-rgb-1998');
%     rgb = lab2rgb([stimulus.backgroundLab(1) a b]);
    stimulus.colorwheel.rgb(ti,:) = rgb;
end

%% Sizes
stimulus.fixWidth = 0.5;
stimulus.targetWidth = 1;

%% Setup patch colors
stimulus.patchcolors = stimulus.colorwheel.rgb(11:11:66,:);

%% Setup patches and stencils

% there will be 12 possible locations where we can show dots. We will use 6
% patches of dots 

patches = 6;
patchEcc = 5;

stimulus.dotScale = 2;
stimulus.dotDirs = [0.5 1.5]*pi;
 
dots = struct;

dots.density = 21;
dots.speed = 3;
dots.maxX = stimulus.targetWidth;
dots.maxY = stimulus.targetWidth;

for di = 1:patches
    stimulus.patches{di} = struct;
    
    % patch dots
    stimulus.patches{di}.dots = initDots(dots);
    
    % color
    stimulus.patches{di}.color = [1 1 1];
end

% setup the locations

locations = 12;

for li = 1:locations
    stimulus.locations{li} = struct;
    % patch location
    stimulus.locations{li}.theta = (li-1)/locations*2*pi;
    disp(stimulus.locations{li}.theta)
    stimulus.locations{li}.ecc = patchEcc;
    
    % x/y
    stimulus.locations{li}.xcenter = stimulus.locations{li}.ecc * cos(stimulus.locations{li}.theta);
    stimulus.locations{li}.ycenter = stimulus.locations{li}.ecc * sin(stimulus.locations{li}.theta);
end

% stencils
mglClearScreen(0);
mglStencilCreateBegin(1);
for li = 1:locations
    mglFillOval(stimulus.locations{li}.xcenter,stimulus.locations{li}.ycenter,[stimulus.targetWidth, stimulus.targetWidth],[1 1 1]);
end
mglFillOval(0,0,[stimulus.targetWidth stimulus.targetWidth],[1 1 1]);
mglFlush;
mglStencilCreateEnd;

%% Extra stuff
stimulus.live.trackingAngle = 0;

%% Create the cue patch

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

task{1}{1}.segmin = [0.5 0.5 1 2 2.5 0.5 2];
task{1}{1}.segmax = [0.5 0.5 1 2 2.5 0.5 8];

if stimulus.debug
    task{1}{1}.segmin = [1 2 1 2 1 5 0.5];
    task{1}{1}.segmax = [1 2 1 2 1 5 0.5];
end

task{1}{1}.waitForBacktick = 1;

task{1}{1}.getResponse = zeros(1,length(task{1}{1}.segmin));
task{1}{1}.getResponse(stimulus.seg.resp) = 1;

task{1}{1}.numTrials = Inf;

task{1}{1}.random = 1;

task{1}{1}.parameter.trialType = [1 2]; % 1 = spatial, 2 = feature
task{1}{1}.parameter.targets = patches/2; % number of elements in the target group
task{1}{1}.parameter.patches = patches; % number of elements in the target group

if ~stimulus.replay && stimulus.scan
    task{1}{1}.synchToVol = zeros(1,length(task{1}{1}.segmin));
    task{1}{1}.synchToVol(end) = 1;
end

% which patches are on
task{1}{1}.randVars.calculated.lOn1 = nan;
task{1}{1}.randVars.calculated.lOn2 = nan;
task{1}{1}.randVars.calculated.lOn3 = nan;
task{1}{1}.randVars.calculated.lOn4 = nan;
task{1}{1}.randVars.calculated.lOn5 = nan;
task{1}{1}.randVars.calculated.lOn6 = nan;

% attention targets
task{1}{1}.randVars.calculated.target1 = nan;
task{1}{1}.randVars.calculated.target2 = nan;
task{1}{1}.randVars.calculated.target3 = nan;

% feature target
task{1}{1}.randVars.calculated.dead = nan;
task{1}{1}.randVars.calculated.targetDir = nan;
task{1}{1}.randVars.calculated.targetAngle = nan;
task{1}{1}.randVars.calculated.respAngle = nan;
task{1}{1}.randVars.calculated.respDistance = nan;
task{1}{1}.randVars.calculated.targetNum = nan; % which of the patches is the target
task{1}{1}.randVars.calculated.cwOffset = nan; % colorwheel offset rotation

%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,@endTrialCallback,[]);
end

%% EYE CALIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~stimulus.replay
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


function [task, myscreen] = startTrialCallback(task,myscreen)
global stimulus

mglSetMousePosition(960,540,1);

if stimulus.powerwheel
    task.thistrial.respAngle = rand*2*pi;
end

% choose where the patches will be located
locationOpts = 1:length(stimulus.locations);
stimulus.live.locationsOn = locationOpts(randperm(6));
for i=1:6
    task.thistrial.(sprintf('lOn%i',i)) = stimulus.live.locationsOn(i);
end

% choose the targets
stimulus.live.targetNums = randsample(stimulus.live.locationsOn,task.thistrial.targets);
for i = 1:3
    task.thistrial.(sprintf('target%i',i)) = stimulus.live.targetNums(i);
end

% choose the true target
task.thistrial.targetNum = randsample(stimulus.live.targetNums,1);

% set the target colors randomly
colorOpts = [zeros(1,length(stimulus.patches)/2) ones(1,length(stimulus.patches)/2)];
colors = colorOpts(randperm(length(colorOpts)));

% set the colors of the patches
for di = 1:length(stimulus.patches)
    if colors(di)
        ctheta = randsample(stimulus.colorwheel.thetas(stimulus.colorwheel.range1),1);
    else
        ctheta = randsample(stimulus.colorwheel.thetas(stimulus.colorwheel.range2),1);
    end
    if stimulus.live.locationsOn(di)==task.thistrial.targetNum
        task.thistrial.targetAngle = ctheta;
    end
    stimulus.patches{di}.color = ang2rgb(ctheta);
end

if task.thistrial.trialType==1
    % if this is a spatial trial we need to set the directions randomly
    
    for di = 1:length(stimulus.patches)
        stimulus.patches{di}.dots.dir = randsample(stimulus.dotDirs,1);
    end
else
    % if this is a motion trial we need to set the directions to be the
    % same for the targets
    task.thistrial.targetDir = randsample(stimulus.dotDirs,1);
    offDir = mod(task.thistrial.targetDir+pi,2*pi);
    
    for di = 1:length(stimulus.patches)
        if any(stimulus.live.targetNums==di)
            stimulus.patches{di}.dots.dir = task.thistrial.targetDir;
        else
            stimulus.patches{di}.dots.dir = offDir;
        end
    end
end

% colorwheel random rotation
task.thistrial.cwOffset = rand*2*pi;

trialTypes = {'locations','colors'};
disp(sprintf('(afcom) Starting trial %i. Attending %s',task.trialnum,trialTypes{task.thistrial.trialType}));

% eye tracking 
task.thistrial.dead = 0;
stimulus.live.eyeCount=0;

function [task, myscreen] = endTrialCallback(task,myscreen)

% if isnan(task.thistrial.respDistance)
%     task.thistrial.respDistance = abs(mod(task.thistrial.respAngle,pi) - mod(task.thistrial.targetAngle,pi));
%     disp(sprintf('Recorded no-response - angle of %1.2f true %1.2f: %1.2f distance',task.thistrial.respAngle,task.thistrial.targetAngle,task.thistrial.respDistance));
% end


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

if task.thistrial.trialType==1
    % spatial - draw lines to attended locations
    
    for ti = 1:length(stimulus.live.targetNums)
        target = stimulus.live.targetNums(ti);
        
        % draw the line from fixWidth to 2*fixWidth
        x = stimulus.fixWidth * cos(stimulus.locations{target}.theta);
        y = stimulus.fixWidth * sin(stimulus.locations{target}.theta);
        mglLines2(x,y,2*x,2*y,1,[1 1 1]);
    end
else
    % feature - draw a half circle indicating motion direction
%     theta = task.thistrial.targetDir-pi/2;
%     
%     mglGluPartialDisk(0,0,stimulus.fixWidth,stimulus.fixWidth+0.05,180/pi*(theta-pi/2),180,[1 1 1]);

    % cue dots version
    stimulus.cueDots = updateDots(stimulus.cueDots,1,false);
    stimulus.cueDots.dir = task.thistrial.targetDir;
    
    mglStencilSelect(1);
    mglPoints2(stimulus.cueDots.x-stimulus.cueDots.maxX/2,stimulus.cueDots.y-stimulus.cueDots.maxY/2,stimulus.dotScale,[1 1 1]);
    mglStencilSelect(0);
    
end

function drawStim(~,grayscale)

global stimulus

mglStencilSelect(1);
% update and collapse x/y coordinates for drawing
n = stimulus.patches{1}.dots.n;

x = zeros(1,n*length(stimulus.patches));
y = x;

for di = 1:length(stimulus.patches)
    stimulus.patches{di}.dots = updateDots(stimulus.patches{di}.dots,1,false);
    
    offX = stimulus.locations{stimulus.live.locationsOn(di)}.xcenter - stimulus.patches{di}.dots.maxX/2;
    offY = stimulus.locations{stimulus.live.locationsOn(di)}.ycenter - stimulus.patches{di}.dots.maxY/2;
    
    if grayscale
        x(((di-1)*n+1):(di*n)) = offX + stimulus.patches{di}.dots.x;
        y(((di-1)*n+1):(di*n)) = offY + stimulus.patches{di}.dots.y;
    else
        mglPoints2(offX+stimulus.patches{di}.dots.x,offY+stimulus.patches{di}.dots.y,stimulus.dotScale,stimulus.patches{di}.color);
    end
end

% draw all the dots at once
if grayscale
    mglPoints2(x,y,stimulus.dotScale,[1 1 1]);
end
mglStencilSelect(0);


function drawFix(task,color)

global stimulus;

if task.thistrial.dead
    mglGluDisk(0,0,[1 1],stimulus.colors.red,60);
else
    mglFixationCross(stimulus.fixWidth,1,color);
end

function drawTarget(task)

global stimulus

mglStencilSelect(1);
mglFillRect(stimulus.locations{task.thistrial.targetNum}.xcenter,stimulus.locations{task.thistrial.targetNum}.ycenter,repmat(stimulus.targetWidth,1,2),stimulus.colors.white);
mglStencilSelect(0);

function drawPicker(task)

global stimulus

% Draw the color picker
for ti = 1:length(stimulus.colorwheel.thetas)
    theta = stimulus.colorwheel.thetas(ti);
    mglGluPartialDisk(0,0,1,1.25,180/pi*(theta-stimulus.colorwheel.thetaInc/2),180/pi*stimulus.colorwheel.thetaInc,stimulus.colorwheel.rgb(ti,:));
end
mglGluPartialDisk(0,0,1,1.25,180/pi*task.thistrial.respAngle-2.51,5,[0.75 0.75 0.75]);

function drawResp(angle)

global stimulus
% Draw the chosen color as a backgroundc ircle

mglFillOval(0,0,stimulus.fixWidth*[1 1],ang2rgb(angle));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus

mglClearScreen();


if (task.thistrial.thisseg==stimulus.seg.resp) && stimulus.powerwheel
    mInfo = mglGetMouse(myscreen.screenNumber);
    curPos = -mInfo.x/90;
    task.thistrial.respAngle = mod(task.thistrial.respAngle + curPos-stimulus.live.trackingAngle,2*pi);
    stimulus.live.trackingAngle = curPos;
elseif task.thistrial.thisseg==stimulus.seg.resp
    mInfo = mglGetMouse(myscreen.screenNumber);
    degx = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
    degy = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;
    task.thistrial.respAngle = mod(-atan2(degy,degx)+pi/2,2*pi);
end


switch task.thistrial.thisseg
    case stimulus.seg.cue
        % fixation
        drawStim(task,true);
        drawCue(task);
        drawFix(task,stimulus.colors.white);
        
    case stimulus.seg.isi
        drawStim(task,true);
        drawFix(task,stimulus.colors.white);
        
    case stimulus.seg.stim
        drawStim(task,false);
        drawFix(task,stimulus.colors.white);
        
    case stimulus.seg.delay
        drawFix(task,stimulus.colors.white);
        
    case stimulus.seg.resp
        drawTarget(task);
        drawPicker(task);
        drawResp(task.thistrial.respAngle);
        drawFix(task,stimulus.colors.white);
        
    case stimulus.seg.feedback
        drawResp(task.thistrial.targetAngle);
        drawFix(task,stimulus.colors.white);
        
    case stimulus.seg.iti
        drawStim(task,true);
        drawFix(task,stimulus.colors.white);
end

% if ~stimulus.replay
%     drawFix(myscreen);
% 
%     % check eye pos
%     if (~stimulus.noeye) && (stimulus.eyewindow>0)
% 
% 
%         % mouse version for testing with no eyetracker
%         if stimulus.mousedebug
%             mInfo = mglGetMouse(myscreen.screenNumber);
%             degx = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
%             degy = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;
% 
%             pos = [degx, degy];
%         else
%             [pos,~] = mglEyelinkGetCurrentEyePos;
%         end
% 
%         if stimulus.debug > 0
%             if stimulus.debug >= 10
%                 disp(sprintf('Mouse position: %1.1f %1.1f',pos(1),pos(2)));
%                 stimulus.debug = 1;
%             else
%                 stimulus.debug = stimulus.debug+1;
%             end
%         end
% 
% 
%         % compute distance
%         dist = hypot(pos(1),pos(2));
%     end
% 
%     % Eye movement detection code
%     if (~stimulus.noeye) && (stimulus.eyewindow>0) && ~task.thistrial.dead
%         if ~any(isnan(pos))
% 
%             if dist > stimulus.eyewindow && stimulus.live.eyeCount > 40
%                 disp('Eye movement detected!!!!');
%                 task.thistrial.dead = 1;
%                 return
%             elseif dist > stimulus.eyewindow
%                 stimulus.live.eyeCount = stimulus.live.eyeCount + 1;
%             end
%         end
%     end
% end

function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus

if task.thistrial.dead, return; end

if stimulus.powerwheel
    validResponse = task.thistrial.mouseButton == 1;
else
    validResponse = task.thistrial.whichButton == stimulus.responseKeys(1);
end

if validResponse
    if task.thistrial.gotResponse==0
        task.thistrial.respDistance = abs(mod(task.thistrial.respAngle,pi) - mod(task.thistrial.targetAngle,pi));
        disp(sprintf('Received response angle of %1.2f true %1.2f: %1.2f distance',task.thistrial.respAngle,task.thistrial.targetAngle,task.thistrial.respDistance));
    else
        disp(sprintf('Subject responded multiple times: %i',task.thistrial.gotResponse));
    end
    task.thistrial.gotResponse=task.thistrial.gotResponse+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helpers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function col = ang2rgb(ang)

global stimulus

col = interp1(stimulus.colorwheel.thetas',stimulus.colorwheel.rgb,ang);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create dots for horizontal motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initDots(dots)

dots.dir = 0;

area = dots.maxX*dots.maxY;

dots.n = area * dots.density;

% make a some points
% dots.n = 500*dots.density;
% make sure it's an even number
dots.n = dots.n + mod(dots.n,2);

% set half to white and half to black
dots.con = repmat([1 2],1,dots.n/2);

dots.x = rand(1,dots.n)*dots.maxX;
dots.y = rand(1,dots.n)*dots.maxY;

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
