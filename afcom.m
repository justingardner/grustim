function [ myscreen ] = afcom( varargin )
% ***scan info: 736 TRs (6 minute * 120 + 16)***
% *** set 'noeye=1' to turn of the eye tracker***
%
%Attention field and color mapping (afcom)
%
%   Map the effects of spatial attention and attention to specific colors 

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

getArgs(varargin,{'scan=0','plots=0','noeye=0','powerwheel=1','eyewindow=1.5','practice=0','debug=0','replay=0','run=0','build=0','mouse=0'});
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
    dispInfo(stimulus);
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
calib = load(fullfile(myscreen.calibFullFilename));
stimulus.calib = calib.calib;

%% Colors
if ~isfield(stimulus,'colors')
    stimulus.colors.white = [1 1 1]; stimulus.colors.red = [1 0 0];
    stimulus.colors.green = [0 1 0]; stimulus.colors.black = [0 0 0];
end

if ~isfield(stimulus,'colorwheel')
    % get the lab space rgb values
    stimulus.backgroundLab = rgb2lab([0.5 0.5 0.5]);

    % setup color wheel
    stimulus.colorwheel.acenter = stimulus.backgroundLab(2);
    stimulus.colorwheel.bcenter = stimulus.backgroundLab(3);

    % compute the ranges around 0 and pi for the colorwheel
    theta_ = pi/64; % increment size
    stimulus.colorwheel.thetaInc = theta_;
    stimulus.colorwheel.thetas = 0:theta_:2*pi;

    D = 60;
    stimulus.colorwheel.distanceLab = D;

    for ti = 1:length(stimulus.colorwheel.thetas)
        theta = stimulus.colorwheel.thetas(ti);

        a = D*cos(theta)+stimulus.colorwheel.acenter;
        b = D*sin(theta)+stimulus.colorwheel.bcenter;

        rgb = mglLab2rgb([stimulus.backgroundLab(1) a b],stimulus.calib);
    %     rgb = lab2rgb([stimulus.backgroundLab(1) a b]);
        stimulus.colorwheel.rgb(ti,:) = rgb;
    end

    % if any values are outside RGB space just cut them off
    stimulus.colorwheel.rgb(stimulus.colorwheel.rgb<0) = 0;
    stimulus.colorwheel.rgb(stimulus.colorwheel.rgb>1) = 1;
end

%% Sizes
stimulus.fixWidth = 0.5;
stimulus.targetWidth = 10;
stimulus.patchEcc = 8;

%% Setup patches and stencils

% there will be 12 possible locations where we can show dots. We will use 6
% patches of dots 

stimulus.dotScale = 0.6;
stimulus.cueScale = 0.15;

stimulus.dotDirs = [0.5 1.5 0.5 1.5]*pi;
 
dots = struct;

dots.density = 0.1;
dots.speed = 4;
dots.maxAlive = myscreen.framesPerSecond/6;
dots.maxX = stimulus.targetWidth;
dots.maxY = stimulus.targetWidth;

thetas = [0 0 pi pi];

for di = 1:4
    stimulus.patches{di} = struct;
    
    % patch dots
    stimulus.patches{di}.dots = initDots(dots);
    stimulus.patches{di}.dots.dir = stimulus.dotDirs(di);

    % color
    stimulus.patches{di}.color = [1 1 1];
    
    % location
    stimulus.patches{di}.theta = thetas(di);
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

task{1}{1}.segmin = [1 0.5 0.5 1 1.5 4 0.75];
task{1}{1}.segmax = [1 0.5 0.5 1 1.5 4 0.75];

if stimulus.debug
    task{1}{1}.segmin = [1 2 0.1 1.5 1 inf 0.5];
    task{1}{1}.segmax = [1 2 0.1 1.5 1 inf 0.5];
end

task{1}{1}.waitForBacktick = 1;

task{1}{1}.getResponse = zeros(1,length(task{1}{1}.segmin));
task{1}{1}.getResponse(stimulus.seg.resp) = 1;

if stimulus.scan==1
    task{1}{1}.numTrials = Inf;
else
    task{1}{1}.numTrials = 32;
end

task{1}{1}.random = 1;

task{1}{1}.parameter.trialType = [1 2]; % 1 = spatial, 2 = feature
task{1}{1}.parameter.target = [1 2 3 4]; % which patch is the target

if ~stimulus.replay && stimulus.scan
    task{1}{1}.synchToVol = zeros(1,length(task{1}{1}.segmin));
    task{1}{1}.synchToVol(end) = 1;
end

% feature target
task{1}{1}.randVars.calculated.dead = nan;
task{1}{1}.randVars.calculated.targetDir = nan;
task{1}{1}.randVars.calculated.targetAngle = nan; % angle of the target
task{1}{1}.randVars.calculated.distractorAngle = nan; % angle of the other thing you had to attend
task{1}{1}.randVars.claculated.distractor = nan;
task{1}{1}.randVars.calculated.angle1 = nan;
task{1}{1}.randVars.calculated.angle2 = nan;
task{1}{1}.randVars.calculated.angle3 = nan;
task{1}{1}.randVars.calculated.angle4 = nan;
task{1}{1}.randVars.calculated.respAngle = nan;
task{1}{1}.randVars.calculated.respDistance = nan;
task{1}{1}.randVars.calculated.distDistance = nan;
task{1}{1}.randVars.calculated.cwOffset = nan; % colorwheel offset rotation

%% Mouse movement storage data

stimulus.data.mouseTrack = zeros(50,200);
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

function dispInfo(stimulus)
%%
files = dir(fullfile('~/data/afcom/',mglGetSID,'*.mat'));

for fi = 1:length(files)
    load(fullfile('~/data/afcom/',mglGetSID,files(fi).name));
    exp = getTaskParameters(myscreen,task);
    e{fi} = exp{1};
end

stop = 1;

function [task, myscreen] = startTrialCallback(task,myscreen)
global stimulus

if stimulus.powerwheel
    mglSetMousePosition(myscreen.screenWidth/2,myscreen.screenHeight/2,1);
else
    mglSetMousePosition(myscreen.screenWidth/2,myscreen.screenHeight/2,2);
end

stimulus.cueDots.dir = stimulus.patches{task.thistrial.target}.dots.dir;

if task.thistrial.trialType==1
    distractors = [2 1 4 3];
else
    distractors = [3 4 1 2];
end
task.thistrial.distractor = distractors(task.thistrial.target);

% set the colors of the patches
for di = 1:length(stimulus.patches)
    ctheta = randsample(stimulus.colorwheel.thetas,1);
    
    if di==task.thistrial.target
        task.thistrial.targetAngle = ctheta;
    end
    if di==task.thistrial.distractor
        task.thistrial.distractorAngle = ctheta;
    end
    
    stimulus.patches{di}.color = ang2rgb(ctheta);
    
    task.thistrial.(sprintf('angle%i',di)) = ctheta;
end

% colorwheel random rotation
task.thistrial.cwOffset = rand*2*pi;
task.thistrial.respAngle = -task.thistrial.cwOffset;

trialTypes = {'locations','colors'};
disp(sprintf('(afcom) Starting trial %i. Attending %s',task.trialnum,trialTypes{task.thistrial.trialType}));

% eye tracking 
task.thistrial.dead = 0;
stimulus.live.eyeCount=0;

% mouse tracking
stimulus.data.mouseTick = 1;

function [task, myscreen] = endTrialCallback(task,myscreen)

if task.thistrial.dead, return; end

if isnan(task.thistrial.respDistance)
    task.thistrial.respDistance = circDist(task.thistrial.respAngle,task.thistrial.targetAngle);
    task.thistrial.distDistance = circDist(task.thistrial.respAngle,task.thistrial.distractorAngle);
    disp(sprintf('Recorded no-response - angle of %1.2f true %1.2f: %1.2f distance',task.thistrial.respAngle,task.thistrial.targetAngle,task.thistrial.respDistance));
end

function d = circDist(t1,t2)
if length(t1)>1 || length(t2)>1
    warning('Circle distance only works on numbers not arrays');
    d = [];
    return
end
% circular distance
t1 = mod(t1,2*pi);
t2 = mod(t2,2*pi);
if t1>t2
    d = t1-t2;
else
    d = t2-t1;
end
if d>pi
    d = 2*pi-d;
end

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
        
    % draw the line from fixWidth to 2*fixWidth
    x = stimulus.fixWidth * cos(stimulus.patches{task.thistrial.target}.theta);
    y = stimulus.fixWidth * sin(stimulus.patches{task.thistrial.target}.theta);
    mglLines2(x,y,2*x,2*y,1,[1 1 1]);
else
    % feature - draw a half circle indicating motion direction
%     theta = task.thistrial.targetDir-pi/2;
%     
%     mglGluPartialDisk(0,0,stimulus.fixWidth,stimulus.fixWidth+0.05,180/pi*(theta-pi/2),180,[1 1 1]);

    % cue dots version
    stimulus.cueDots = updateDots(stimulus.cueDots,1,false);
    
    mglStencilSelect(1);
    afPoints(stimulus.cueDots.x-stimulus.cueDots.maxX/2,stimulus.cueDots.y-stimulus.cueDots.maxY/2,stimulus.cueScale,[1 1 1]);
    mglStencilSelect(0);
end

function drawStim(~,grayscale)

global stimulus

mglStencilSelect(1);
% update and collapse x/y coordinates for drawing
n = stimulus.patches{1}.dots.n;

x = zeros(1,n*length(stimulus.patches));
y = x;
r = ones(1,n*length(stimulus.patches));
g = r;
b = r;

for di = 1:length(stimulus.patches)
    stimulus.patches{di}.dots = updateDots(stimulus.patches{di}.dots,1,false);
    
    offX = stimulus.patches{di}.xcenter - stimulus.patches{di}.dots.maxX/2;
    offY = stimulus.patches{di}.ycenter - stimulus.patches{di}.dots.maxY/2;
    
    x(((di-1)*n+1):(di*n)) = offX + stimulus.patches{di}.dots.x;
    y(((di-1)*n+1):(di*n)) = offY + stimulus.patches{di}.dots.y;
    if ~grayscale
        r(((di-1)*n+1):(di*n)) = stimulus.patches{di}.color(1);
        g(((di-1)*n+1):(di*n)) = stimulus.patches{di}.color(2);
        b(((di-1)*n+1):(di*n)) = stimulus.patches{di}.color(3);
    end
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
    drawBorder(locations{li}.xcenter,locations{li}.ycenter,r,[0.2 0.2 0.2]);
end

function drawTarget(task)

global stimulus

color = ang2rgb(task.thistrial.respAngle);

stimulus.patches{task.thistrial.target}.dots = updateDots(stimulus.patches{task.thistrial.target}.dots,1,false);
offX = stimulus.patches{task.thistrial.target}.xcenter - stimulus.patches{task.thistrial.target}.dots.maxX/2;
offY = stimulus.patches{task.thistrial.target}.ycenter - stimulus.patches{task.thistrial.target}.dots.maxY/2;

mglStencilSelect(1);
afPoints(stimulus.patches{task.thistrial.target}.dots.x + offX,stimulus.patches{task.thistrial.target}.dots.y + offY,stimulus.dotScale,color);
mglStencilSelect(0);

function drawPicker(task)

global stimulus

% Draw the color picker
for ti = 1:length(stimulus.colorwheel.thetas)
    theta = stimulus.colorwheel.thetas(ti) + task.thistrial.cwOffset;
    mglGluPartialDisk(0,0,1,1.25,180/pi*(theta-stimulus.colorwheel.thetaInc/2),180/pi*stimulus.colorwheel.thetaInc,stimulus.colorwheel.rgb(ti,:));
end
mglGluPartialDisk(0,0,1,1.25,180/pi*(task.thistrial.respAngle+task.thistrial.cwOffset)-2.5,5,[0.75 0.75 0.75]);

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

if (task.thistrial.thisseg==stimulus.seg.resp)
    if stimulus.powerwheel
        mInfo = mglGetMouse(myscreen.screenNumber);
        curPos = mInfo.x/90;
        task.thistrial.respAngle = mod(task.thistrial.respAngle + curPos-stimulus.live.trackingAngle,2*pi);
        stimulus.live.trackingAngle = curPos;
    else
        mInfo = mglGetMouse(myscreen.screenNumber);
        degx = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
        degy = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;
        task.thistrial.respAngle = mod(-atan2(degy,degx)+pi/2 - task.thistrial.cwOffset,2*pi);
    end
    
    stimulus.data.mouseTrack(task.trialnum,stimulus.data.mouseTick) = task.thistrial.respAngle;
    stimulus.data.mouseTick = stimulus.data.mouseTick + 1;
end

switch task.thistrial.thisseg
        
    case stimulus.seg.iti
        drawStim(1:4,true);
        drawAllBorders(stimulus.patches,stimulus.targetWidth/2);
        drawFix(task,stimulus.colors.white);
        
    case stimulus.seg.cue
        % fixation
        drawStim(1:4,true);
        drawCue(task);
        drawAllBorders(stimulus.patches,stimulus.targetWidth/2);
        drawFix(task,stimulus.colors.white);
        
    case stimulus.seg.isi
        drawStim(1:4,true);
        drawAllBorders(stimulus.patches,stimulus.targetWidth/2);
        drawFix(task,stimulus.colors.white);
        
    case stimulus.seg.stim
        drawStim(1:4,false);
        drawAllBorders(stimulus.patches,stimulus.targetWidth/2);
        drawFix(task,stimulus.colors.white);
        
    case stimulus.seg.delay
        drawAllBorders(stimulus.patches,stimulus.targetWidth/2);
        drawFix(task,stimulus.colors.white);
        
    case stimulus.seg.resp
        drawAllBorders(stimulus.patches,stimulus.targetWidth/2);
        drawTarget(task);
        drawPicker(task);
        drawResp(task.thistrial.respAngle);
        drawFix(task,stimulus.colors.white);
        
    case stimulus.seg.feedback
        drawAllBorders(stimulus.patches,stimulus.targetWidth/2);
        drawResp(task.thistrial.targetAngle);
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

if task.thistrial.gotResponse==0
    % jump to the feedback segment
    task = jumpSegment(task);
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
