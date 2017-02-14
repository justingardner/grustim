function [ myscreen ] = unlearning( varargin )
%UNLEARNING 
%
% Unconscious learning of spatial patterns. This experiment runs two
% simultaneous tasks. 
%
% COLORS
%   Valid: FF8F8F, accent: B25252
%   Impossible: 5EA1CC, accent: 3F86B2
%   Neutral: FFFEA8

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
        stimulus.counter = s.stimulus.counter + 1;
        % copy patterns
        stimulus.patterns = s.stimulus.patterns;
        stimulus.patternopts = s.stimulus.patternopts;
        stimulus.holdoutopts = s.stimulus.holdoutopts;
        stimulus.npatterns = s.stimulus.npatterns;
        stimulus.learn = s.stimulus.learn;

        clear s;
        disp(sprintf('(unlearn) Data file: %s loaded.',fname));
        
    else
        if isfield(stimulus,'patternopts') || isfield(stimulus,'learn')
            disp('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX');
            disp(sprintf('WARNING: File was not loaded--but patterns are already set for this subject'));
            disp('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX');
        end
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
stimulus.cur{end}.spc = 2.2;
% rows/cols are per quadrant
stimulus.cur{end}.rows = 3;
stimulus.cur{end}.cols = 3;
stimulus.cur{end}.pos = linspace(-3*stimulus.cur{end}.spc,3*stimulus.cur{end}.spc,7);
stimulus.cur{end}.N = 36;
stimulus.cur{end}.K = 24;

stimulus.cur_ = stimulus.cur{end};

if ~isfield(stimulus,'learn')
% %     opts = [2 3 6 7];
% %     stimulus.learn = opts(randi(4));
    stimulus.learn = 1;
    disp('(unlearn) WARNING: New quadrant part chosen for learning');
end
disp(sprintf('(unlearn) Subject %s is learning quadrant part #%i',mglGetSID,stimulus.learn));

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
    dispInfo;
    return
end

%% Setup Task
task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;

stimulus.curTrial = 0;

task{1}{1}.segmin = [inf .500 0.750 .500 0.2 1.500 0.250];
task{1}{1}.segmax = [inf .500 0.750 .500 0.2 1.500 1.000];

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
task{1}{1}.numTrials = 60;
task{1}{1}.random = 1;
task{1}{1}.parameter.match = [0 1];
task{1}{1}.parameter.impossible = [0 0 0 0 0 0 1 1 1 1];
task{1}{1}.parameter.vertical1 = [10 12 14];
task{1}{1}.parameter.difficulty = [1 2];

if stimulus.scan
    task{1}{1}.synchToVol(stimulus.seg.ITI) = 1;
end

% select 6 patterns for permanent use
if ~isfield(stimulus,'patternopts')
    disp('(unlearn) WARNING: New pattern options are generating...');
    % caution: the patterns need to be controlled for having consistent
    % proportions of 1 / 0 at each location. This is key in making sure
    % that all spatial locations provide equal information. Fortunately
    % there are 15 111100 patterns, so every location can be randomized
    % from one of these to generate the full pattern. We will generate two
    % sets of patterns, a main pattern, and a holdout set.
    bases = unique(perms([1 1 1 1 0 0]),'rows');
    
    % caution--this can take a while, we want exactly 24 1s per row and 4
    % per column!
    patterns = zeros(36,6);
    holdout = zeros(36,6);
    while any(sum(patterns,1)~=24)
        for r = 1:36 % which row we are on
            patterns(r,:) = bases(randi(15),:);
        end
    end
    while any(sum(holdout,1)~=24)
        for r = 1:36 % which row we are on
            holdout(r,:) = bases(randi(15),:);
        end
    end
    
    stimulus.patterns = [patterns';holdout'];
    stimulus.patternopts = 1:6;
    stimulus.holdoutopts = 7:12;
end
stimulus.npatterns = length(stimulus.patternopts);
task{1}{1}.parameter.pattern1 = 1:stimulus.npatterns; % which test pattern to use for the first stim

%% Tracking

% these are variables that we want to track for later analysis.
task{1}{1}.randVars.calculated.correct = nan;
task{1}{1}.randVars.calculated.pattern2 = nan;
task{1}{1}.randVars.calculated.vertical2 = nan;

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
    dispInfo();
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
disp(sprintf('(unlearn) %s:%s trial (%i), %s. Pattern A: %i, pattern B: %i',opts{task.thistrial.match+1},dopts{task.thistrial.match+1},task.trialnum,vopts{task.thistrial.impossible+1},task.thistrial.pattern1,task.thistrial.pattern2));
    
stimulus.live.eyeCount = 0;
stimulus.dead = 0;

function task = buildData(task)

global stimulus
% Setup the displays
% .rings holds N rings consisting of num segments
gxpos = 1:6;
gypos = 1:6;

data2 = zeros(2*stimulus.cur_.rows,2*stimulus.cur_.cols); % controls the contrast 2 levels
data1 = zeros(2*stimulus.cur_.rows,2*stimulus.cur_.cols); % controls the contrast 1 levels
sz = size(data1);

gsize = 36;
gx = 6; gy = 6;

for group = 1:size(gxpos,1)
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
        % in patA and patB we're required to have 12=0, 12=1, and 12=2
        % here we set 12=2
        maskA = find(patA==1);
        maskA = maskA(randperm(length(maskA)));
        patA(maskA(1:12)) = 2;
        maskB = find(patB==1);
        maskB = maskB(randperm(length(maskB)));
        patB(maskB(1:12)) = 2;
        data2(gypos(group,:),gxpos(group,:)) = reshape(patB,gy,gx);
        data1(gypos(group,:),gxpos(group,:)) = reshape(patA,gy,gx);
    else
        disp('full pattern model has failed');
        % pick how many will be 0/1/2
        lgroup = randi(3,1,gsize)-1; count = 1;
        perm1 = lgroup(randperm(length(lgroup)));
        perm2 = lgroup(randperm(length(lgroup)));
        for x = gxpos(group,:);
            for y = gypos(group,:);
                % add random stuff
                data2(y,x) = perm2(count);
                data1(y,x) = perm1(count);
                count = count+1;
            end
        end
    end
end

% set orientation information
mask1 = data1>0;
mask2 = data2>0;

% set orient1 positions (depending on task.thistrial.vertical1)
mask1 = find(mask1(:));
mask1 = mask1(randperm(length(mask1))); % randomize

orient1 = zeros(1,gsize);
orient1(mask1(1:task.thistrial.vertical1)) = 1;

if task.thistrial.impossible
    % change a bunch randomly
    inc=0;
elseif task.thistrial.match==1
    % increase verticals (ones)
    inc=task.thistrial.difficulty;
else
    % increase horizontals (zeros) 
    inc=-task.thistrial.difficulty;
end
task.thistrial.vertical2 = task.thistrial.vertical1+inc;

pos = find(mask2(:));
pos = pos(randperm(length(pos)));
orient2 = zeros(1,gsize);
orient2(pos(1:task.thistrial.vertical2)) = 1;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus
mglClearScreen(0.5);

if stimulus.dead && mglGetSecs(task.thistrial.segStartSeconds)>1
    jumpSegment(task,inf); stimulus.dead=0;
end

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

if ~stimulus.noeye && ~any(task.thistrial.thisseg==[stimulus.seg.ITI1 stimulus.seg.ITI2 stimulus.seg.resp]) && ~stimulus.scan
    if ~any(isnan(pos))
        if dist > 3 && stimulus.live.eyeCount > 30
            disp('Eye movement detected!!!!');
            stimulus.dead = 1;
            stimulus.live.eyeDead=1;
            return
        elseif dist > 3
            stimulus.live.eyeCount = stimulus.live.eyeCount + 1;
        end
    end
end


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
    if stimulus.live.resp==1
        mglTextSet([],32,stimulus.live.fixColor);
        mglTextDraw(stimulus.live.respText,[0 0]);
    else
        mglFillOval(0,0,[1 1],0.5);
        upFix(stimulus);
    end
end

if ~stimulus.noeye && stimulus.live.triggerWaiting
    now = mglGetSecs;
    % check eye position, if 
    if ~any(isnan(pos))
        wasCentered = stimulus.live.centered;
        stimulus.live.centered = dist<3;
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
function dispInfo(task,myscreen)
%%
global stimulus
rstimulus = stimulus;

% ctask = task; cscreen = myscreen; % save this incase we need them

% compute % correct for valid and invalid trials, display learning over
% time (including history from other runs)
% exp = getTaskParameters(task,myscreen);
disp('(unlearn) Display info not implemented yet');

% get the files list
files = dir(fullfile(sprintf('~/data/unlearning/%s/*.mat',mglGetSID)));

% load the files and pull out the data (long form)
%  run #    local trial     real trial   impossible   match   vert1
%    1             2            3            4          5        6
%  vert2  pattern1    pattern2    response    correct
%     7      8           9           10         11
count = 1; data = zeros(10000,11);

for fi = 1:length(files)
    load(fullfile(sprintf('~/data/unlearning/%s/%s',mglGetSID,files(fi).name)));
    
    e = getTaskParameters(myscreen,task);
    e = e{1}; % why?!
    
    data(count:count+(e.nTrials-1),:) = [repmat(stimulus.counter,e.nTrials,1) (1:e.nTrials)' (count:count+(e.nTrials-1))' ...
        e.parameter.impossible' e.parameter.match' e.parameter.vertical1' ...
        e.randVars.vertical2' e.parameter.pattern1' e.randVars.pattern2' ...
        e.response' e.randVars.correct'];
    
    count = count+e.nTrials;
end

data = data(1:(count-1),:);

% separate data into impossible and valid
idata = data(data(:,4)==1,:);
vdata = data(data(:,4)==0,:);

% check statistics across sessions
uruns = unique(data(:,1));

vci_ = zeros(length(uruns),2);
ici_ = zeros(length(uruns),2);

for ri = 1:length(uruns)
    run = uruns(ri);
    % valid
    vdat = vdata(vdata(:,1)==run,11);
    vci = bootci(1000,@nanmean,vdat); vperf = mean(vci);
    vcis = sprintf('[%2.0f%% %2.0f%%]',vci(1)*100,vci(2)*100);
    vci_(ri,:) = vci;
    % impossible
    idat = idata(idata(:,1)==run,11);
    ici = bootci(1000,@nanmean,idat); iperf = mean(ici);
    icis = sprintf('[%2.0f%% %2.0f%%]',ici(1)*100,ici(2)*100);
    ici_(ri,:) = ici;
    
    
    disp(sprintf('Performance on run %i. Valid: %2.0f%% %s Impossible: %2.0f%% %s',run,100*vperf,vcis,100*iperf,icis));
end
%% write data
% header = {'Run','Trial'
% csvwriteh(fullfile(sprintf('~/data/unlearning/%s/%s_data.mat',mglGetSID,mglGetSID)),data,header);
save(fullfile(sprintf('~/data/unlearning/%s/%s_data.mat',mglGetSID,mglGetSID)),'data');

%% plot
h = figure; hold on;

offset = .05;
title(sprintf('Subj: %s performance and error over time',mglGetSID));
% valid
plot(uruns-offset,mean(vci_,2),'o','MarkerFaceColor',rstimulus.colors.valid','MarkerEdgeColor',[1 1 1],'MarkerSize',5);
errbar(uruns-offset,mean(vci_,2),vci_(:,2)-mean(vci_,2),'-','Color',rstimulus.colors.valid);

% invalid
plot(uruns+offset,mean(ici_,2),'o','MarkerFaceColor',rstimulus.colors.impossible','MarkerEdgeColor',[1 1 1],'MarkerSize',5);
errbar(uruns+offset,mean(ici_,2),ici_(:,2)-mean(ici_,2),'-','Color',rstimulus.colors.impossible);

legend({'Valid','Impossible'});
z = hline(0.5,'--k');
% set(z,'Color',stimulus.colors.chance);

xlabel('Run (#)');
ylabel('Performance (% correct)');

set(gca,'XTick',uruns);
set(gca,'YTick',[0 0.25 0.5 0.75 1],'YTickLabels',{'0','25%','50%','75%','100%'});

drawPublishAxis

savepdf(h,fullfile(sprintf('~/data/unlearning/%s/%s_performance.pdf',mglGetSID,mglGetSID)));

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
    sz = 2;
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