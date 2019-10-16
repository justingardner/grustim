function [ myscreen ] = bananacat_search( varargin )
%% LOW vs. HIGH-LEVEL FB SEARCH 
%
%   Feature-based search for categorical images that have had their low
%   level features equalized by a texture matching algorithm. Unlimited
%   time search version.
%
%   This is a behavioral search task which uses eye saccades for reports.
%   The program loads the different categories of images in the folder
%   ~/proj/bananacat/images/*group*/out/*
%   It assumes the folders in images/group/ have the category names and
%   that the folders in out/ are named category1category2 with numbered
%   images inside (e.g. 11_*.png, 12_*.png, etc)
%
%   On each trial the search target is shown and then one target among some
%   number of distractors will appear. The observer must saccade to the
%   target and fixate it for 300 ms. In the noeye version (testing) the
%   code simply continues after 5 s. 
%
%   Performance is evaluated as RT relative to the number of distractors
%
%    Usage: bananacat(varargin)
%    Authors: Dan Birman
%    Date: 05/23/2018
%
%    Parameters:
%       ecc: the eccentricity of the targets
%       size: the diameter of the targets

global stimulus

% clear stimulus -- we will load the old stimulus information
if isfield(stimulus,'tex')
    tex = stimulus.tex;
    mask = stimulus.mask;
    stimulus = struct; 
    stimulus.tex = tex;
    stimulus.mask = mask;
else
    stimulus = struct;
end

%% Initialize Variables

% add arguments later
eccmax=0;sz=0;plots=0;debug=0;training=0;group='';scan=0;len=0;distractors=0;genTex=0;noeye=0;gray=0;
getArgs(varargin,{'eccmax=15','sz=5','distractors=10','plots=0','gray=0','noeye=0','debug=0','training=0','group=objects','scan=0','genTex=0'});
stimulus.scan = scan;
stimulus.eccmax = eccmax;
stimulus.sz = sz;
stimulus.dist=distractors;
stimulus.len = len;
stimulus.noeye = noeye;
stimulus.gray=gray;
stimulus.plots = plots;
stimulus.debug = debug;
stimulus.training = training;
stimulus.group = group;

clear plots debug training ecc sz group scan distractors len noeye

if stimulus.scan
    warning('Not setup for scanning');
end

%% Setup stimulus

% get category information
stimulus.folder = '~/proj/bananacat/images';
cats = dir(fullfile(stimulus.folder,stimulus.group));
cats = {cats.name};
cats = cats(cellfun(@(x) isempty(strfind(x,'.')),cats));
cats = cats(cellfun(@(x) isempty(strfind(x,'out')),cats));
stimulus.categories = cats;
stimulus.iNums = 1:4;

%% Open Old Stimfile
stimulus.counter = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/bananacat/%s',stimulus.condition,mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/bananacat/%s/1*mat',stimulus.condition,mglGetSID));

    if length(files) >= 1
        fname = files(end).name;
        
        s = load(sprintf('~/data/bananacat/%s/%s',stimulus.condition,mglGetSID,fname));
        % copy staircases and run numbers
        stimulus.counter = s.stimulus.counter + 1;
        stimulus.performance = s.stimulus.performance;
        
        clear s;
        disp(sprintf('(bananacat) Data file: %s loaded.',fname));
    end
end
disp(sprintf('(bananacat) This is run #%i',stimulus.counter));

%% Setup Screen
myscreen = initScreen('VPixx');

% set background to grey
myscreen.background = 0.5;

%% Staircase
if ~isfield(stimulus,'performance')
    disp('(bananacat) WARNING: New performance matrix');
    initPerformance();
end

%% Plot and return
if stimulus.plots==2
    dispInfo(stimulus);
    return
end

%% Initialize Stimulus

myscreen = initStimulus('stimulus',myscreen);

if genTex || ~isfield(stimulus,'tex')
    localInitStimulus(myscreen);
end
%     
% if stimulus.powerwheel
%     stimulus.responseKeys = 5; % no-detection key (when powerwheel is in use)
% else
%     stimulus.responseKeys = [1 2]; % 
% end

%% Colors
stimulus.colors.black = [0 0 0]; stimulus.colors.white = [1 1 1];
stimulus.colors.red = [1 0 0]; stimulus.colors.green = [0 1 0];
stimulus.colors.gray = [0.75 0.75 0.75];

%% Setup Task

stimulus.curTrial(1) = 0;

task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;

stimulus.seg.cue = 1;
stimulus.seg.isi = 2;
stimulus.seg.stim = 3;
stimulus.seg.iti = 4;

task{1}{1}.seglen = [0.250 0.750 inf 0.500];

if stimulus.noeye
    task{1}{1}.seglen(stimulus.seg.stim) = 5;
end

task{1}{1}.synchToVol = zeros(size(task{1}{1}.seglen));
task{1}{1}.getResponse = zeros(size(task{1}{1}.seglen));
task{1}{1}.numTrials = 20;
task{1}{1}.random = 1;

if stimulus.scan
    task{1}{1}.synchToVol(stimulus.seg.ITI) = 1;
end

% Task trial parameters
task{1}{1}.parameter.eccmax = stimulus.eccmax;
task{1}{1}.parameter.sz = stimulus.sz;
task{1}{1}.parameter.dist = stimulus.dist;
task{1}{1}.parameter.gray = stimulus.gray;

task{1}{1}.parameter.matchStyle = [0 1];
task{1}{1}.parameter.matchCat = [0 1];

task{1}{1}.randVars.calculated.targetStyle = nan;
task{1}{1}.randVars.calculated.targetCat = nan;
task{1}{1}.randVars.calculated.targetRot = nan;
task{1}{1}.randVars.calculated.tsi = nan;
task{1}{1}.randVars.calculated.tci = nan;

task{1}{1}.randVars.calculated.distractStyle = nan;
task{1}{1}.randVars.calculated.distractCat = nan;
% we don't keep track of the distractor rotation or numbers (too much info)

task{1}{1}.randVars.calculated.choice = nan;
task{1}{1}.randVars.calculated.correct = nan;
task{1}{1}.randVars.calculated.dead = false;

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
disp(sprintf('(bananacat) Starting run number: %i.',stimulus.counter));

%% Main Task Loop

mglClearScreen(0.5); mglFixationCross(1,1,stimulus.colors.white);
mglFlush
mglClearScreen(0.5); mglFixationCross(1,1,stimulus.colors.white);

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
mglFlush
myscreen.flushMode = 1;

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

if stimulus.plots
    disp('(bananacat) Displaying plots');
    dispInfo(stimulus);
end

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startTrialCallback(task,myscreen)
%%
global stimulus

task.thistrial.targetStyle = randi(length(stimulus.categories));
if task.thistrial.matchStyle
    task.thistrial.distractStyle = task.thistrial.targetStyle;
else
    opts = 1:length(stimulus.categories);
    opts = setdiff(opts,task.thistrial.targetStyle);
    task.thistrial.distractStyle = opts(randi(length(opts)));
end
task.thistrial.targetCat = randi(length(stimulus.categories));
if task.thistrial.matchCat
    task.thistrial.distractCat = task.thistrial.targetCat;
else
    opts = 1:length(stimulus.categories);
    opts = setdiff(opts,task.thistrial.targetCat);
    task.thistrial.distractCat = opts(randi(length(opts)));
end
task.thistrial.targetRot = rand*360;
task.thistrial.tsi = randsample(stimulus.iNums,1);
task.thistrial.tci = randsample(stimulus.iNums,1);

% get the textures
stimulus.live.targetTex = stimulus.tex{task.thistrial.targetStyle,task.thistrial.targetCat,task.thistrial.tsi,task.thistrial.tci,randi(5)};
for di = 1:stimulus.dist
    stimulus.live.distractTexs(di) = stimulus.tex{task.thistrial.distractStyle,task.thistrial.distractCat,randsample(stimulus.iNums,1),randsample(stimulus.iNums,1),randi(5)};
    stimulus.live.distractRots(di) = rand*360;
end

% generate random x/y positions within the max eccentricity distance
stimulus.live.x = rand(1,stimulus.dist+1) * task.thistrial.eccmax*2 - task.thistrial.eccmax;
stimulus.live.y = rand(1,stimulus.dist+1) * task.thistrial.eccmax*2 - task.thistrial.eccmax;

dists = hypot(stimulus.live.x,stimulus.live.y)>task.thistrial.eccmax;
while any(dists)
    stimulus.live.x(dists) = rand(1,sum(dists)) * task.thistrial.eccmax*2 - task.thistrial.eccmax;
    stimulus.live.y(dists) = rand(1,sum(dists)) * task.thistrial.eccmax*2 - task.thistrial.eccmax;
    dists = hypot(stimulus.live.x,stimulus.live.y)>task.thistrial.eccmax;
end

disp(sprintf('(bananacat) Trial (%i): %s%s hidden among %i %s%s, size %1.1f, len %i ms',...
    task.trialnum,stimulus.categories{task.thistrial.targetStyle},...
    stimulus.categories{task.thistrial.targetCat},...
    task.thistrial.dist,...
    stimulus.categories{task.thistrial.distractStyle},...
    stimulus.categories{task.thistrial.distractCat},...
    task.thistrial.sz));
    
stimulus.live.eyeCount = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)
%%

global stimulus

stimulus.live.triggerWaiting = 0;
if any(task.thistrial.thisseg==[stimulus.seg.iti])
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
stimulus.live.mask = 0;

if task.thistrial.thisseg==stimulus.seg.stim
    stimulus.live.stim = 1;
elseif task.thistrial.thisseg==stimulus.seg.cue
    stimulus.live.cue = 1;
    stimulus.live.fix = 0;
end

for i = 1:2    
    mglClearScreen(0.5);
    if stimulus.live.cue
        % place the cue in the center and skip the fixation cross
        mglBltTexture(stimulus.live.targetTex,[0 0],0,0,task.thistrial.targetRot);
    elseif stimulus.live.stim
        % get the rotation pattern 
        mglBltTexture([stimulus.live.targetTex stimulus.live.distractTexs],[stimulus.live.x' stimulus.live.y'],0,0,[task.thistrial.targetRot stimulus.live.distractRots]);
    elseif stimulus.live.resp
        % draw circles around the saccade targets
        for xi = 1:length(stimulus.live.x)
            if xi==1
                rot = task.thistrial.targetRot;
            else
                rot = stimulus.live.distractRots(xi-1);
            end
            drawCircle(stimulus.live.x(xi),stimulus.live.y(xi),stimulus.sz,rot);
        end
    end
    
    if stimulus.live.fix
        upFix(stimulus);
    end

    mglFlush
end

function drawCircle(x,y,rad,thetaoffset)
global stimulus

for i = 0:9 % use 10 pieces
    mglGluPartialDisk(x,y,rad/3-0.025,rad/3+0.025,i*36+thetaoffset,18,stimulus.colors.gray);
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

% if (task.thistrial.thisseg==stimulus.seg.resp) && stimulus.powerwheel
%     mInfo = mglGetMouse(myscreen.screenNumber);
%     curPos = -mInfo.x/90;
%     stimulus.live.angle = stimulus.live.angle + curPos-stimulus.live.trackingAngle;
%     if abs(curPos-stimulus.live.trackingAngle)>0
%         stimulus.live.anyAngleAdj = true;
%     end
%     stimulus.live.trackingAngle = curPos;
%     convertRespXY(task);
% elseif task.thistrial.thisseg==stimulus.seg{task.thistrial.thisphase}.resp % powerwheel==0
%     keys = find(mglGetKeys);
%     if any(keys==19)
%         stimulus.live.angle = stimulus.live.angle+0.01;
%         stimulus.live.anyAngleAdj = true;
%     elseif any(keys==20)
%         stimulus.live.angle = stimulus.live.angle-0.01;
%         stimulus.live.anyAngleAdj = true;
%     end
%     convertRespXY(task);
% end
% 
% if stimulus.live.resp && ((task.thistrial.thisphase==2) || stimulus.test2)
%     mglClearScreen(0.5);
%     if stimulus.live.fix, upFix(stimulus); end
%     if stimulus.live.resp, mglFillOval(stimulus.live.respx,stimulus.live.respy,[0.5 0.5],stimulus.colors.white); end
% end
% 
% % Trial trigger on eye fixation code  
% if ~stimulus.noeye && stimulus.live.triggerWaiting
%     now = mglGetSecs;
%     % check eye position, if 
%     if ~any(isnan(pos))
%         wasCentered = stimulus.live.centered;
%         stimulus.live.centered = dist<2.5;
%         if wasCentered && stimulus.live.centered && stimulus.live.lastTrigger>0
%             stimulus.live.triggerTime = stimulus.live.triggerTime + now-stimulus.live.lastTrigger;
%         end
%         stimulus.live.lastTrigger = now;
%     end
%     if stimulus.live.triggerTime > 0.5 % not in ms dummy, wait 1.5 seconds (reasonable slow time)
%         disp('Starting trial--eye centered and space pressed.');
%         task = jumpSegment(task);
%     end
% end

% function convertRespXY(task)
% global stimulus
% 
% stimulus.live.respx = task.thistrial.ecc*cos(stimulus.live.angle+task.thistrial.target);
% stimulus.live.respy = task.thistrial.ecc*sin(stimulus.live.angle+task.thistrial.target);

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
    
if isfield(task.thistrial,'whichButton') && (task.thistrial.whichButton==stimulus.responseKeys(1)) && isempty(task.thistrial.mouseButton)
    % subject didn't see anything
    task = jumpSegment(task,inf);
    task.thistrial.detected = -1; %-1 means they reported not seeing anything
    disp(sprintf('Subject reported not seeing %02.02f%% contrast stimulus', task.thistrial.contrast*100));
    return
end

if stimulus.powerwheel
    validResponse = task.thistrial.mouseButton == 1;
else
    validResponse = task.thistrial.whichButton == stimulus.responseKeys(4);
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
            task.thistrial.detected = 1;
            stimulus.staircase = doStaircase('update',stimulus.staircase,task.thistrial.detected);
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

function initPerformance()
global stimulus

stimulus.performance = [];

%%%%%%%%%%%%%%%%%%%%%%%
%    dispInfo    %
%%%%%%%%%%%%%%%%%%%%%%%
function dispInfo(rstimulus)
%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localInitStimulus(msc)

global stimulus

deg2pix = msc.screenWidth/msc.imageWidth; % how many pixels per degree

cropDim = 255;
bot = floor(cropDim/2);
top = floor(cropDim/2);

gauss = mglMakeGaussian(1,1,1/6,1/6,0,0,cropDim,cropDim);
% gauss = min(1,gauss ./ max(gauss(:)/5) );
% load every combination of categories for images 1:4
disppercent(-1/length(stimulus.categories));
for si = 1:length(stimulus.categories)
    for ci = 1:length(stimulus.categories)
        for sii = 1:length(stimulus.iNums)
            for cii = 1:length(stimulus.iNums)
                style = stimulus.categories{si};
                category = stimulus.categories{ci};
%                 disp(sprintf('%s%s',style,category));
                sNum = stimulus.iNums(sii);
                cNum = stimulus.iNums(cii);
                
                img = imread(fullfile(stimulus.folder,stimulus.group,'out',sprintf('%s%s',style,category),sprintf('%i%i_at_iteration_10.png',sNum,cNum)));
                
                if stimulus.gray
                    img = mean(img,3);
                    img = uint8(repmat(img,1,1,3));
                end
                
                % generate five random 256*256 crops
                xmin = bot;
                xmax = size(img,1)-top;
                ymin = bot;
                ymax = size(img,2)-top;
                
                for cropi = 1:5
                    xc = randi(xmax-xmin)+xmin;
                    yc = randi(ymax-ymin)+ymin;
                    crop = img(xc-bot:xc+top,yc-bot:yc+top,:);
%                     figure; imagesc(crop);
                    crop(:,:,4) = gauss*255;
                    
                    crop = imresize(crop,stimulus.sz*[deg2pix deg2pix]);
                    stimulus.tex{si,ci,sii,cii,cropi} = mglCreateTexture(permute(crop,[3 1 2]));
                end
                
%                 mglClearScreen(0.5);
%                 mglBltTexture(stimulus.tex{si,ci,sii,cii,1},[0 0],0,0,rand*360);
%                 mglFlush
%                 pause(0.25);
                
            end
        end
    end
    disppercent(si/length(stimulus.categories));
end
disppercent(inf/length(stimulus.categories));

for mi = 1:200
    if stimulus.gray
        mask = uint8(rand(cropDim,cropDim)*255);
        mask = repmat(mask,1,1,3);
    else
        mask = uint8(rand(cropDim,cropDim,3)*255);
    end
    mask(:,:,4) = gauss*255;
    
    stimulus.mask{mi} = mglCreateTexture(permute(mask,[3 1 2]));
%     mglClearScreen(0.5);
%     mglBltTexture(stimulus.mask{mi},[0 0],0,0,rand*360);
%     mglFlush
%     pause(0.25);
end