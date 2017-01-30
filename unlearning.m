function [ myscreen ] = unlearning( varargin )
%UNLEARNING 
%
% Unconscious learning of spatial patterns. This experiment runs two
% simultaneous tasks. 
%

global stimulus

stimulus = struct;
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

% spc = 1.1;
% xs = 9;
% xpos = linspace(-floor(xs/2)*spc,floor(xs/2)*spc,xs);
% ys = 7;
% ypos = linspace(-floor(ys/2)*spc,floor(ys/2)*spc,ys);
% data = zeros(5,7);

stimulus.cur{end+1} = struct;
stimulus.cur{end}.spc = 1.1;
% rows/cols are per quadrant
stimulus.cur{end}.rows = 4;
stimulus.cur{end}.cols = 4;
stimulus.cur{end}.pos = linspace(-4*stimulus.cur{end}.spc,4*stimulus.cur{end}.spc,9);
stimulus.cur{end}.N = 8;
stimulus.cur{end}.K = 5;

stimulus.cur_ = stimulus.cur{end};

if ~isfield(stimulus,'learn')
    stimulus.learn = randi(8);
    disp('(unlearn) WARNING: New quadrant part chosen for learning');
end
disp(sprintf('(unlearn) Subject %s is learning quadrant part #%i',mglGetSID,stimulus.learn));

%% Initialize Variables

% add arguments later
scan = 0;
plots = 0;
noeye = 0;
debug = 0;
getArgs(varargin,{'scan=0','plots=0','noeye=1','debug=0'});
stimulus.scan = scan;
stimulus.plots = plots;
stimulus.noeye = noeye;
stimulus.debug = debug;
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

%% Setup Task
task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;

stimulus.curTrial = 0;

task{1}{1}.segmin = [inf .650 1.000 .650 0.500 1.500 0.500];
task{1}{1}.segmax = [inf .650 1.000 .650 0.500 1.500 1.500];

if stimulus.debug
    task{1}{1}.segmin = [inf 1.5 1.000 1.5 0.500 1.500 0.500];
    task{1}{1}.segmax = [inf 1.5 1.000 1.5 0.500 1.500 1.500];
end

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
task{1}{1}.parameter.impossible = [0 0 0 0 0 0 1 1 1 1];

if stimulus.scan
    task{1}{1}.synchToVol(stimulus.seg.ITI) = 1;
end

stimulus.patterns = unique(perms([ones(1,stimulus.cur_.K) zeros(1,stimulus.cur_.N-stimulus.cur_.K)]),'rows');
% select 5 patterns for permanent use
if ~isfield(stimulus,'patternopts')
    stimulus.patternopts = randperm(size(stimulus.patterns,1),5);
    stimulus.generalopts = setdiff(1:size(stimulus.patterns,1),stimulus.patternopts);
    disp('(unlearn) WARNING: New pattern options detected');
end
stimulus.npatterns = length(stimulus.patternopts);
task{1}{1}.parameter.pattern1 = 1:stimulus.npatterns; % which test pattern to use for the first stim

%% Tracking

% these are variables that we want to track for later analysis.
task{1}{1}.randVars.calculated.correct = nan;
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

task = buildData(task);

opts = {'Non-match','Match'};
% if  
dopts = {'Right','Left'};
vopts = {'Valid','Impossible'};
disp(sprintf('(unlearn) %s:%s trial, %s. Pattern A: %i, pattern B: %i',opts{task.thistrial.match+1},dopts{task.thistrial.match+1},vopts{task.thistrial.impossible+1},task.thistrial.pattern1,task.thistrial.pattern2));
    
stimulus.live.eyeCount = 0;
stimulus.dead = 0;

function task = buildData(task)

global stimulus
% Setup the displays
% .rings holds N rings consisting of num segments

gxpos = [1:2;3:4;5:6;7:8;7:8;5:6;3:4;1:2];
gypos = [1:4;1:4;1:4;1:4;5:8;5:8;5:8;5:8];

data2 = zeros(2*stimulus.cur_.rows,2*stimulus.cur_.cols); % controls the contrast 2 levels
data1 = zeros(2*stimulus.cur_.rows,2*stimulus.cur_.cols); % controls the contrast 1 levels
orient2 = zeros(2*stimulus.cur_.rows,2*stimulus.cur_.cols); % controls orientation 2
orient1 = zeros(2*stimulus.cur_.rows,2*stimulus.cur_.cols); % controls orientation 1


for group = 1:8
    % set contrast information
    if group==stimulus.learn
        % set same different pattern
        patB = stimulus.patterns(stimulus.patternopts(task.thistrial.pattern1),:);
        if task.thistrial.match==1
            patA = patB;
            task.thistrial.pattern2 = task.thistrial.pattern1;
        else
            notpatterns = 1:length(stimulus.patternopts);
            notpatterns = notpatterns(notpatterns~=task.thistrial.pattern1);
            task.thistrial.pattern2 = notpatterns(randi(length(notpatterns)));
            patA = stimulus.patterns(stimulus.patternopts(task.thistrial.pattern2),:);
        end
        for j = 1:stimulus.cur_.N
            if patA(j)==1 && round(rand)
                patA(j) = 2;
            end
            if patB(j)==1 && round(rand)
                patB(j) = 2;
            end
        end
        data2(gypos(group,:),gxpos(group,:)) = reshape(patB,4,2);
        data1(gypos(group,:),gxpos(group,:)) = reshape(patA,4,2);
    else
        for x = gxpos(group,:);
            for y = gypos(group,:);
                % add random stuff
                data2(y,x) = randi(3)-1;
                data1(y,x) = randi(3)-1;
            end
        end
    end
end

% set orientation information
mask1 = data1>0;

sum0=0; sum1=0;
while sum0<8 || sum1<8
    orient1 = round(rand(size(data1))).*mask1; % compute initial orientations
    sz = size(orient1); orient1 = orient1(:);

    pos0 = find(orient1==0);
    pos1 = find(orient1==1);
    sum0 = sum(pos0);
    sum1 = sum(pos1);
end

rand0 = randperm(length(pos0));
rand1 = randperm(length(pos1));

if task.thistrial.impossible
    % change a bunch randomly
    orient2(pos0(rand0(1:4))) = 1;
    orient2(pos1(rand1(1:4))) = 0;
elseif task.thistrial.match==1
    % increase verticals (ones)
    orient2(pos0(rand0(1:6))) = 1;
    orient2(pos1(rand1(1:2))) = 0;
else
    % increase horizontals (zeros) 
    orient2(pos0(rand0(1:2))) = 1;
    orient2(pos1(rand1(1:6))) = 0;
end

orient1 = reshape(orient1,sz);
orient2 = reshape(orient2,sz);
    
stimulus.live.data1 = data1;
stimulus.live.data2 = data2;
stimulus.live.orient1 = orient1;
stimulus.live.orient2 = orient2;

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
% if ~stimulus.noeye && ~stimulus.scan
%     [pos,~] = mglEyelinkGetCurrentEyePos;
%     if ~any(isnan(pos))
%         dist = hypot(pos(1),pos(2));
%         if dist > stimulus.ring.inner && stimulus.live.eyeCount > 30
%             mglTextSet([],32,stimulus.colors.red);
%             disp('Eye movement detected!!!!');
%             mglTextDraw('Eye Movement Detected',[0 0]);
%             mglFlush
%             myscreen.flushMode = 1;
%             stimulus.dead = 1;
%             return
%         elseif dist > stimulus.ring.inner-1
%             stimulus.live.eyeCount = stimulus.live.eyeCount + 1;
%         end
%     end
% end


if stimulus.live.stim
    if task.thistrial.thisseg==stimulus.seg.stim1
        data = stimulus.live.data1;
        orient = stimulus.live.orient1;
    else
        data = stimulus.live.data2;
        orient = stimulus.live.orient2;
    end
    
    % draw background on debug
    if stimulus.debug
        partialDiskFuckOGL(0,0,stimulus.cur_.isize-0.5,stimulus.cur_.osize+3,(stimulus.learn-1)*stimulus.cur_.angle,stimulus.cur_.angle,[163 93 93]/255,6,2);
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
end

if stimulus.live.fix
%      cover
    if stimulus.live.resp==1
        mglTextSet([],32,stimulus.live.fixColor);
        mglTextDraw(stimulus.live.respText,[0 0]);
    else
        mglFillOval(0,0,[1 1],0.5);
        upFix(stimulus);
    end
end

% if stimulus.live.triggerWaiting
%     now = mglGetSecs;
%     % check eye position, if 
%     [pos,~] = mglEyelinkGetCurrentEyePos;
%     if ~any(isnan(pos))
%         dist = hypot(pos(1),pos(2));
%         wasCentered = stimulus.live.centered;
%         stimulus.live.centered = dist<2;
%         if wasCentered && stimulus.live.centered && stimulus.live.lastTrigger>0
%             stimulus.live.triggerTime = stimulus.live.triggerTime + now-stimulus.live.lastTrigger;
%         end
%         stimulus.live.lastTrigger = now;
%     end
%     if stimulus.live.triggerTime > 0.5 % not in ms dummy, wait 1.5 seconds (reasonable slow time)
%         disp('Eye position centered');
%         task = jumpSegment(task);
%     end
% end

function upData(data,orient,stimulus)

xp = stimulus.cur_.pos([1:4 6:9]);
yp = fliplr(xp);

for x = 1:8
    for y = 1:8
        if data(x,y)>0
            mglBltTexture(stimulus.live.gratings{data(x,y)},[xp(y) yp(x)],0,0,orient(x,y)*90);
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
        disp(sprintf('Subject pressed %i/%s: %s %s',task.thistrial.whichButton,sideText{task.thistrial.whichButton},matchText{stimulus.responseKeys(task.thistrial.whichButton)},responseText{task.thistrial.correct+1}));
        stimulus.live.fixColor = fixColors{task.thistrial.correct+1};
        stimulus.live.resp = 1;
        stimulus.live.fix = 1;
        stimulus.live.respText = respText{task.thistrial.correct+1};
    else
        disp(sprintf('(unlearn) Subject responded multiple times: %i',task.thistrial.gotResponse+1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              HELPER FUNCTIONS                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
%    dispInfo    %
%%%%%%%%%%%%%%%%%%%%%%%
function dispInfo(task,myscreen,stimulus)
%%

% compute % correct for valid and invalid trials, display learning over
% time (including history from other runs)
% exp = getTaskParameters(task,myscreen);
disp('(unlearn) Display info not implemented yet');

function partialDiskFuckOGL(x,y,isize,osize,sangle,dist,col,slices,loops)
% fucking opengl! wtf!
mglGluPartialDisk(x,y,isize,osize,fuckopengl(sangle),-dist,col,slices,loops);

function deg = fuckopengl(deg)
% BECAUSE FUCK YOU OPENGL
deg = -deg+90;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localInitStimulus()

global stimulus

gratings = cell(1,stimulus.cur_.N);
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
    sz = 1;
    % use total degs / num to compute size
    grating = 255/2*mglMakeGrating(sz,sz,6,0) + 255/2;
    lgrating = (255*0.25)/2*mglMakeGrating(sz,sz,6,0) + 255/2;
%     gratings{gi} = mglCreateTexture(grating);
    gauss = mglMakeGaussian(sz,sz,sz/6,sz/6);
%     alphamask = zeros(size(gauss,1),size(gauss,2),4);
    alphamask = repmat(grating,1,1,4);
    alphamaskl = repmat(lgrating,1,1,4);
    alphamask(:,:,4) = gauss*255;
    alphamaskl(:,:,4) = gauss*255;
    gratings{1} = mglCreateTexture(alphamaskl);
    gratings{2} = mglCreateTexture(alphamask); % high contrast
%     mglBltTexture(gratings{gi,1},[crad 0],0,0,round(rand)*90);
% end
% 
stimulus.live.gratings = gratings;

% mglFlush

%% test rectangle
% mglClearScreen
% spc = 1.1;
% xs = 9;
% xpos = linspace(-floor(xs/2)*spc,floor(xs/2)*spc,xs);
% ys = 7;
% ypos = linspace(-floor(ys/2)*spc,floor(ys/2)*spc,ys);
% data = zeros(5,7);
% for x = 1:xs
%     for y = 1:ys
%         if (~xpos(x)==0) || (~ypos(y)==0)
%             mglBltTexture(gratings{randi(2)},[xpos(x),ypos(y)],0,0,round(rand)*90);
%         end
%     end
% end
% mglFlush
%% test polar
% add = 90;
% mglClearScreen(0.5)
% for out = 1:8
%     if mod(out,2) % odd
%         inc = 90 / (out+1);
%         pos = 0;
%         while pos<=360
%             [x,y] = pol2cart(deg2rad(pos),out);
%             if 1, mglBltTexture(gratings{randi(2)},[x,y],0,0,round(rand)*add+pos); end
%             pos = pos+inc;
%         end
%     else % even
%         inc = 90 / (out+2);
%         pos = inc/2;
%         while pos<=360
%             [x,y] = pol2cart(deg2rad(pos),out);
%             if 1, mglBltTexture(gratings{randi(2)},[x,y],0,0,round(rand)*add+pos); end
%             pos = pos+inc;
%         end
%     end
% end
% % mglBltTexture(gratings{2},[1 0],0,0,0);
% mglFlush