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
stimulus.cur{end}.N = 5;
stimulus.cur{end}.K = 3;
stimulus.cur{end}.angle = 30;
stimulus.cur{end}.num = 360/stimulus.cur{end}.angle;
stimulus.cur{end}.buffer = 8; % buffer is used to stencil over the wedges
stimulus.cur{end}.isize = 1.5;
stimulus.cur{end}.osize = 11;
stimulus.cur{end}.rotNum = 2; % number of wedges to rotate (2*30 = 60)

stimulus.cur_ = stimulus.cur{end};

if ~isfield(stimulus,'learn')
    stimulus.learn = randi(stimulus.cur_.num);
    disp('(unlearn) WARNING: New wedge chosen for learning');
end
disp(sprintf('(unlearn) Subject %s is learning wedge #%i',mglGetSID,stimulus.learn));

%% Initialize Variables

% add arguments later
scan = 0;
plots = 0;
noeye = 0;
debug = 0;
getArgs(varargin,{'scan=0','plots=0','noeye=1','debug=1'});
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
task{1}{1}.randVars.calculated.correlation = nan;

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

task = buildRings(task);

opts = {'Non-match','Match'};
% if  
dopts = {'Right','Left'};
disp(sprintf('(unlearn) %s:%s trial. Pattern A: %i, pattern B: %i',opts{task.thistrial.match+1},dopts{task.thistrial.match+1},task.thistrial.pattern1,task.thistrial.pattern2));
    
stimulus.live.eyeCount = 0;
stimulus.dead = 0;

function task = buildRings(task)

global stimulus
% Setup the displays
% .rings holds N rings consisting of num segments

stimulus.live.rings2 = zeros(stimulus.cur_.N,stimulus.cur_.num);
sz = size(stimulus.live.rings2);

patA = stimulus.patterns(stimulus.patternopts(task.thistrial.pattern1),:)';

if task.thistrial.match==1
    % track pattern
    patB = patA;
    task.thistrial.pattern2 = task.thistrial.pattern1;
else
    notpatterns = 1:length(stimulus.patternopts);
    notpatterns = notpatterns(notpatterns~=task.thistrial.pattern1);
    task.thistrial.pattern2 = notpatterns(randi(length(notpatterns)));
    patB = stimulus.patterns(stimulus.patternopts(task.thistrial.pattern2),:)';
end

% randomize pattern A 1/0.5
for j = 1:stimulus.cur_.N
    if patA(j)==1 && round(rand)
        patA(j) = 0.5;
    end
    if patB(j)==1 && round(rand)
        patB(j) = 0.5;
    end
end

sum1 = 0; sum5 =0; sum0 = 0;
while any([sum1 sum5 sum0]<=7)

    % fill in rings2
    for i = 1:stimulus.cur_.num
        stimulus.live.rings2(:,i) = (randi(3,stimulus.cur_.N,1)-1)/2;
    end

    stimulus.live.rings2 = stimulus.live.rings2(:);
    pos1 = stimulus.live.rings2==1;
    pos5 = stimulus.live.rings2==0.5;
    pos0 = stimulus.live.rings2==0;
    
    sum1 = sum(pos1); sum5 = sum(pos5); sum0 = sum(pos0);
end

stimulus.live.rings1 = stimulus.live.rings2; % copy
% at this point there are at minimum 7 positions with 1, 0.5 and 0 each
% we can safely increase decrease appropriately

pos1 = find(pos1); pos5 = find(pos5); pos0 = find(pos0);
rand1 = randperm(length(pos1)); rand5 = randperm(length(pos5)); rand0 = randperm(length(pos0));

% 
if task.thistrial.impossible
    % maintain overall contrast
    % drop 1->0.5: 3
    stimulus.live.rings1(rand1(1:3)) = 0.5;
    % drop 0.5->0: 3
    stimulus.live.rings1(rand5(1:3)) = 0;
    % inc 0.5->1: 3
    stimulus.live.rings1(rand5(4:6)) = 1;
    % inc 0->0.5: 3
    stimulus.live.rings1(rand0(1:3)) = 0.5;
elseif task.thistrial.match==1
    % increase overall contrast from 1->2, so decrease here
    % drop 1->0.5: 4
    stimulus.live.rings1(rand1(1:4)) = 0.5;
    % drop 0.5->0: 4
    stimulus.live.rings1(rand5(1:4)) = 0;
    % inc 0.5->1: 2
    stimulus.live.rings1(rand5(5:6)) = 1;
    % inc 0->0.5: 2
    stimulus.live.rings1(rand0(1:2)) = 0.5;
else
    % decrease overall contrast from 1->2, so increase here
    % drop 1->0.5: 2
    stimulus.live.rings1(rand1(1:2)) = 0.5;
    % drop 0.5->0: 2
    stimulus.live.rings1(rand5(1:2)) = 0;
    % inc 0.5->1: 4
    stimulus.live.rings1(rand5(3:6)) = 1;
    % inc 0->0.5: 4
    stimulus.live.rings1(rand0(1:2)) = 0.5;
end

% randomize orientation, cause why not?
for i = 1:length(stimulus.live.rings1)
    if rand>=0.5
        stimulus.live.rings1(i) = stimulus.live.rings1(i)*-1;
        stimulus.live.rings2(i) = stimulus.live.rings2(i)*-1;
    end
end

stimulus.live.rings1 = reshape(stimulus.live.rings1,sz);
stimulus.live.rings1(:,stimulus.learn) = patA;
stimulus.live.rings2 = reshape(stimulus.live.rings2,sz);
stimulus.live.rings2(:,stimulus.learn) = patB;

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
        rings = stimulus.live.rings1;
    else
        rings = stimulus.live.rings2;
    end
    
    % draw background on debug
    if stimulus.debug
        partialDiskFuckOGL(0,0,stimulus.cur_.isize-0.5,stimulus.cur_.osize+3,(stimulus.learn-1)*stimulus.cur_.angle,stimulus.cur_.angle,[163 93 93]/255,6,2);
    end
    
    % draw rings
    upRing(rings,stimulus);
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

function upRing(ring,stimulus)

for rn = 1:stimulus.cur_.N
    for si = 0:(stimulus.cur_.num-1)
        thetad = si*stimulus.cur_.angle+stimulus.cur_.angle/2;
        theta = deg2rad(thetad);
        xy = [stimulus.live.pos(rn)*cos(theta) stimulus.live.pos(rn)*sin(theta)];
        thetad = thetad + 90 * (ring(rn,si+1)<0);
        if ring(rn,si+1)==0.5 || ring(rn,si+1)==-0.5
            mglBltTexture(stimulus.live.gratings{rn,1},xy,0,0,thetad);
        elseif ring(rn,si+1)==1 || ring(rn,si+1)==-1
            mglBltTexture(stimulus.live.gratings{rn,2},xy,0,0,thetad);
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
function dispInfo(stimulus)
%%

% if ~stimulus.localizer && ~stimulus.staircasing
%     disp(sprintf('Participant %s has earned $%2.2f',mglGetSID,stimulus.run.points/100));
% end
% % load the luminance table
% % % % load(myscreen.calibFullFilename)
% % % % luminance = interp1(calib.tableCorrected.outputValues,calib.tableCorrected.luminance,0:1/255:255);
% if stimulus.staircasing
%     %%
%     notstaircase = stimulus.staircase;
%     thresholds = zeros(size(stimulus.run.stimLengths));
%     for i = 1:length(stimulus.staircase)
%         out = doStaircase('threshold',notstaircase{i},'type','weibull','dispFig=0');
%         thresholds(i) = out.threshold;
%     end
%     % reorganize into matrix
%     stimCons = unique(stimulus.run.stimCon);
%     stimCons = sort(stimCons);
%     stimLengths = unique(stimulus.run.stimLengths);
%     stimLengths = sort(stimLengths);
%     datamat = nan(length(stimCons),length(stimLengths),5);
%     for ci = 1:length(stimCons)
%         for li = 1:length(stimLengths)
%             idxs = logical((stimulus.run.stimLengths==stimLengths(li)) .* (stimulus.run.stimCon==stimCons(ci)));
%             datamat(ci,li,1:sum(idxs)) = thresholds(idxs);
%         end
%     end
%     if any(thresholds<0) || any(thresholds>1)
%         % remove errant thresholds
%         warning('should remove some thresholds...');
%     end
%     datamat(datamat>1) = NaN;
%     datamat(datamat<=0) = NaN;
%     %%
%     datamu = nanmean(datamat,3);
%     datamu(datamu==0) = NaN;
%     datamu = round((1-datamu)*255);
%     datasd = nanstd(datamat,[],3);
%     datasd(datasd==0) = NaN;
%     %% plot
%     cmap = brewermap(length(stimCons)+1,'Purples');
%     cmap = cmap(2:end,:);
%     figure, hold on
%     legs = {};
%     for i = 1:length(stimCons)
%         plot(stimLengths,datamu(i,:),'o','MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',[1 1 1],'MarkerSize',10);
%         errbar(stimLengths,datamu(i,:),datasd(i,:),'-','Color',cmap(i,:));
%         legs{end+1} = sprintf('Stimulus luminance: %i/255',stimCons(i));
%     end
%     a = axis;
%     axis([50 100 0 a(4)]);
%     legend(legs)
%     xlabel('Stimulus length (ms)');
%     ylabel('Mask contrast at just noticeable difference (% luminance)');
%     drawPublishAxis
% elseif stimulus.localizer
% else
% %     perf = zeros(size(stimulus.istaircase));
% %     for i = 1:length(stimulus.istaircase)
% %         perf(i) = mean(stimulus.istaircase(i).response);
% %     end
% %     figure
% %     plot(1:length(perf),perf,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
% %     set(gca,'XAxisTick',1:length(perf));
% %     drawPublishAxis
% end

function partialDiskFuckOGL(x,y,isize,osize,sangle,dist,col,slices,loops)
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
stimulus.live.pos = logspace(log10(stimulus.cur_.isize),log10(stimulus.cur_.osize),stimulus.cur_.N);

% mglClearScreen(0.5)
for gi = 1:stimulus.cur_.N
    % for each grating distance
    % calculate the center position to estimate the radius
    crad = stimulus.live.pos(gi);
    % get total degrees around circle
    degs = 2*pi*crad;
    sz = degs/stimulus.cur_.num;
    % use total degs / num to compute size
    grating = 255/2*mglMakeGrating(sz,sz,5/sqrt(crad),0) + 255/2;
    lgrating = (255*0.25)/2*mglMakeGrating(sz,sz,5/sqrt(crad),0) + 255/2;
%     gratings{gi} = mglCreateTexture(grating);
    gauss = mglMakeGaussian(sz,sz,sz/6,sz/6);
%     alphamask = zeros(size(gauss,1),size(gauss,2),4);
    alphamask = repmat(grating,1,1,4);
    alphamaskl = repmat(lgrating,1,1,4);
    alphamask(:,:,4) = gauss*255;
    alphamaskl(:,:,4) = gauss*255;
    gratings{gi,1} = mglCreateTexture(alphamaskl);
    gratings{gi,2} = mglCreateTexture(alphamask); % high contrast
%     mglBltTexture(gratings{gi,1},[crad 0],0,0,round(rand)*90);
end

stimulus.live.gratings = gratings;

% mglFlush
