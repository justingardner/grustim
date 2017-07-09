function [ myscreen ] = afmapb( varargin )
%ATTENTIONFIELDMAPPING 
%
%   Map the attention field using behavior. This function works by having a
%   participant perform an asynchronous attention task at fixation or in a
%   quarterfield region. The entire background of the screen consists of
%   white noise. By comparing the white noise on trials that were correct
%   against trials that were incorrect we can estimate what regions of
%   space affect the task--thus obtaining an estimate of the attention
%   field (combined with the readout field?). 

%%

global stimulus

stimulus = struct;

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
clear localizer invisible scan noeye task test2

if stimulus.scan
    warning('Not setup for scanning');
end


%% Stimulus parameters

stimulus.gaussFWHM = 2;
stimulus.gaussSD = stimulus.gaussFWHM/(2*sqrt(2*log(2)));
stimulus.gaussX = 5;
stimulus.gaussY = 5;
stimulus.pixelSize = 8;

%% Open Old Stimfile
stimulus.counter = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/afmapb/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/afmapb/%s/1*mat',mglGetSID));

    if length(files) >= 1
        fname = files(end).name;
        
        s = load(sprintf('~/data/afmapb/%s/%s',mglGetSID,fname));
        % copy staircases and run numbers
        stimulus.counter = s.stimulus.counter + 1;
        stimulus.staircase = s.stimulus.staircase;
        clear s;
        disp(sprintf('(afmapb) Data file: %s loaded.',fname));
    end
end
disp(sprintf('(afmapb) This is run #%i',stimulus.counter));

%% Setup Screen
myscreen = initScreen('VPixx');

myscreen.stimulusNames{1} = 'stimulus';
% set background to grey
myscreen.background = 0.5;

%% Staircase
stimulus.useStair = true;
if ~isfield(stimulus,'staircase')
    disp('(afmapb) WARNING: New staircase');
    disp('(afmapb) Staircase is running');
    initStair();
elseif true %mod(stimulus.counter,5)==0
    disp('(afmapb) Staircase is running');
    resetStair();
else
%     stimulus.useStair = false;
%     stimulus.out = doStaircase('threshold',stimulus.staircase,'type=weibull','dispFig=0');
%     stimulus.contrast = stimulus.out.threshold;
%     while false
%         val = input(sprintf('Current value is %01.2f, [enter] or change: ',stimulus.contrast));
%         if ~isempty(val)
%             stimulus.contrast = val;
%         else
%             break;
%         end
%     end
%     disp(sprintf('(afmapb) Contrast is fixed at %01.2f',stimulus.contrast));
end

%% Plot and return
if stimulus.plots==2
    dispInfo(stimulus);
    return
end
%% White noise tracking
stimulus.wn.count = 1;
stimulus.wn.img = zeros(10000,myscreen.screenWidth/stimulus.pixelSize,myscreen.screenHeight/stimulus.pixelSize);
stimulus.wn.trials = cell(1,10000);


%% Initialize Stimulus
    
stimulus.responseKeys = [2 1]; % absent present

%% Colors
stimulus.colors.white = [1 1 1]; stimulus.colors.red = [1 0 0];
stimulus.colors.green = [0 1 0]; stimulus.colors.black = [0 0 0];
% % initGammaTable(myscreen);
% stimulus.colors.rmed = 127.5;
% 
% % We're going to add an equal number of reserved colors to the top and
% % bottom, to try to keep the center of the gamma table stable.
% stimulus.colors.reservedBottom = [1 0 0; 0 0 0]; % fixation cross colors
% stimulus.colors.reservedTop = [1 1 1; 0 1 0]; % correct/incorrect colors
% stimulus.colors.black = 1/255; stimulus.colors.white = 254/255;
% stimulus.colors.red = 0/255; stimulus.colors.green = 255/255;
% stimulus.colors.nReserved = 2; % this is /2 the true number, because it's duplicated
% stimulus.colors.nUnreserved = 256-(2*stimulus.colors.nReserved);
% 
% stimulus.colors.mrmax = stimulus.colors.nReserved - 1 + stimulus.colors.nUnreserved;
% stimulus.colors.mrmin = stimulus.colors.nReserved;

% stimulus.gaussX = 0;
% stimulus.gaussY = 15;
%% Setup Gaussian
[X,Y] = meshgrid((0.5:(myscreen.screenWidth/stimulus.pixelSize-0.5))-myscreen.screenWidth/(stimulus.pixelSize*2),(0.5:(myscreen.screenHeight/stimulus.pixelSize-0.5))-myscreen.screenHeight/(stimulus.pixelSize*2));
ppdw = myscreen.screenWidth/myscreen.imageWidth;
ppdh = myscreen.screenHeight/myscreen.imageHeight;
if ~(ppdw==ppdh)
    warning('PIXELS ARE NOT SQUARE');
end
stimulus.live.X = X*stimulus.pixelSize./ppdw;
stimulus.live.Y = Y*stimulus.pixelSize./ppdh;
% pre-compute distance from gaussian
stimulus.live.dist = normpdf(hypot(stimulus.live.X-stimulus.gaussX,stimulus.live.Y-stimulus.gaussY),0,stimulus.gaussSD)';
stimulus.live.dist = uint8(stimulus.live.dist ./ max(stimulus.live.dist(:)) * 255);
stimulus.live.dist = repmat(reshape(stimulus.live.dist,[1 size(stimulus.live.dist)]),3,1,1);
% imagesc(squeeze(stimulus.live.dist(1,:,:)));
%%
% imagesc(flipud(squeeze(stimulus.live.dist(1,:,:))'));
% axis equal
%% Setup Attention Task

stimulus.curTrial = 0;

task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;
% task waits for fixation on first segment
task{1}{1}.seglen = [inf 0.500 0.200 1];
% task{1}{1}.segmax = [inf 0.500 0.200 1];

if stimulus.noeye
    task{1}{1}.seglen(1) = 1;
end

stimulus.seg.ITI = 1;
stimulus.seg.delay1 = 2;
stimulus.seg.stim = 3;
stimulus.seg.resp = 4;

task{1}{1}.synchToVol = [0 0 0 0];
task{1}{1}.getResponse = [0 0 0 1];

task{1}{1}.numTrials = 60;

task{1}{1}.parameter.present = [0 1];

task{1}{1}.random = 1;

task{1}{1}.randVars.calculated.contrast = nan;
task{1}{1}.randVars.calculated.resp = nan;
task{1}{1}.randVars.calculated.correct = nan;
task{1}{1}.randVars.calculated.hit = false;
task{1}{1}.randVars.calculated.fa = false;
task{1}{1}.randVars.calculated.miss = false;
task{1}{1}.randVars.calculated.cr = false;
task{1}{1}.randVars.calculated.dead = false;

%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback1,@screenUpdateCallback1,@getResponseCallback1,@startTrialCallback1,[],[]);
end

%% EYE CALIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

% let the user know
disp(sprintf('(afmapb) Starting run number: %i.',stimulus.counter));

%% Main Task Loop

% setGammaTable(1);
mglClearScreen(0.5); mglFixationCross(1,1,stimulus.colors.black);

mglFlush
mglClearScreen(0.5); mglFixationCross(1,1,stimulus.colors.black);

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

stimulus.wn.img = stimulus.wn.img(1:(stimulus.wn.count-1),:,:);
stimulus.wn.trials = stimulus.wn.trials(1:(stimulus.wn.count-1));
stimulus.wn.count = stimulus.wn.count-1;

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

if stimulus.plots
    disp('(afmapb) Displaying plots');
    dispInfo(stimulus);
end

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback1(task,myscreen)
global stimulus

stimulus.live.triggerWaiting = 0;
if any(task.thistrial.thisseg==[stimulus.seg.ITI])
    stimulus.live.triggerWaiting = 1;
    stimulus.live.centered = 0;
    stimulus.live.triggerTime = 0;
    stimulus.live.lastTrigger = -1;
end

stimulus.live.eyeDead = 0;
stimulus.live.fixColor = stimulus.colors.white;

if task.thistrial.thisseg==stimulus.seg.ITI
    stimulus.live.fixColor = stimulus.colors.black;
elseif task.thistrial.thisseg==stimulus.seg.stim
    refreshWN(task,myscreen);
    mglBltTexture(stimulus.live.tex,[0 0 myscreen.imageWidth myscreen.imageHeight]);
    mglFlush
    mglClearScreen(0.5);
    mglFlush
end

function [task, myscreen] = startTrialCallback1(task,myscreen)
%%
global stimulus

stimulus.curTrial = stimulus.curTrial+1;

if stimulus.useStair
    [task.thistrial.contrast, stimulus.staircase] = doStaircase('testValue',stimulus.staircase);
else
    task.thistrial.contrast = stimulus.contrast;
end

stimulus.live.wnTimer = -1;

% refreshWN(task,myscreen); it will refresh on the first frame anyways, no
% need to do it here

disp(sprintf('(afmapb) Trial %i: %02.1f',stimulus.curTrial,task.thistrial.contrast*100));
    
stimulus.live.eyeCount = 0;

function refreshWN(task,myscreen)
global stimulus
wn = repmat(randi(256,1,myscreen.screenWidth/stimulus.pixelSize,myscreen.screenHeight/stimulus.pixelSize,'uint8')-1,3,1,1);
% save the white noise
stimulus.wn.img(stimulus.wn.count,:,:) = wn(1,:,:);
stimulus.wn.trials{stimulus.curTrial}(end+1) = stimulus.wn.count;
stimulus.wn.count = stimulus.wn.count+1;
% check whether we need to add the gaussian
if task.thistrial.present
    % add the gaussian
    wn = min(wn+task.thistrial.contrast*stimulus.live.dist,255);
end
wn(4,:,:) = 255;
if isfield(stimulus,'live') && isfield(stimulus.live,'wn')
    mglDeleteTexture(stimulus.live.tex);
end
stimulus.live.tex = mglCreateTexture(wn,[],0,{'GL_TEXTURE_MAG_FILTER','GL_NEAREST'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback1(task, myscreen)
%%
global stimulus

mglClearScreen(0.5);

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
if ~stimulus.noeye && ~any(task.thistrial.thisseg==[stimulus.seg.ITI]) && ~stimulus.scan
    if ~any(isnan(pos))
        if dist > 1.5 && stimulus.live.eyeCount > 20
            disp('Eye movement detected!!!!');
            task.thistrial.dead = 1;
            stimulus.live.eyeDead=1;
            return
        elseif dist > 1.5
            stimulus.live.eyeCount = stimulus.live.eyeCount + 1;
        end
    end
end

if (task.thistrial.thisseg==stimulus.seg.stim)
    1;
%     stimulus.live.wnTimer = stimulus.live.wnTimer+1;
%     if mod(stimulus.live.wnTimer,3)==0
%         refreshWN(task,myscreen);
%     end
%     mglBltTexture(stimulus.live.tex,[0 0 myscreen.imageWidth myscreen.imageHeight]);
    
elseif (task.thistrial.thisseg==stimulus.seg.ITI)
    for i = 1:8
        mglGluPartialDisk(stimulus.gaussX,stimulus.gaussY,stimulus.gaussFWHM/2-.01,stimulus.gaussFWHM/2+.01,(i-1)*360/8-11.25,360/16,stimulus.colors.white);
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

upFix(stimulus);

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
function [task, myscreen] = getResponseCallback1(task, myscreen)

global stimulus

colors = [stimulus.colors.red;stimulus.colors.green];
text = {'Incorrect','Correct'};
stext = {'Absent','Present'};
if any(task.thistrial.whichButton==stimulus.responseKeys)
    if task.thistrial.gotResponse==0
        task.thistrial.resp = stimulus.responseKeys(task.thistrial.whichButton)-1;
        
        task.thistrial.correct = task.thistrial.resp==task.thistrial.present;
        
        if task.thistrial.present
            if task.thistrial.correct
                task.thistrial.hit = true;
            else
                task.thistrial.miss = true;
            end
        else
            if task.thistrial.correct
                task.thistrial.cr = true;
            else
                task.thistrial.fa = true;
            end
        end
        
        stimulus.staircase = doStaircase('update',stimulus.staircase,task.thistrial.correct);
        
        stimulus.live.fixColor = colors(task.thistrial.correct+1,:);
        disp(sprintf('Subject responded %s: %s',stext{task.thistrial.whichButton},text{task.thistrial.correct+1}));
    else
        disp(sprintf('Subject responded multiple times: %i',task.thistrial.gotResponse));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              HELPER FUNCTIONS                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function initStair()
global stimulus

stimulus.staircase = doStaircase('init','upDown',...
            'initialThreshold',1,...
            'initialStepsize',0.1,...
            'minThreshold=0.0001','maxThreshold=1','stepRule','pest',...
            'nTrials=80','maxStepsize=0.3','minStepsize=0.0001');
        
function resetStair()

global stimulus

if doStaircase('stop',stimulus.staircase)
    disp('(afmapb) Staircase is being reset');
    stimulus.staircase(end+1) = doStaircase('init',stimulus.staircase(end));
    if stimulus.staircase(end).s.threshold>0.3
        disp('(afmapb) Bad staircase threshold: setting to 0.3');
        stimulus.staircase(end).s.threshold=0.3;
    elseif stimulus.staircase(end).s.threshold<0
        disp('(afmapb) Bad staircase threshold: setting to 0.05');
        stimulus.staircase(end).s.threshold=0.05;
    end
end

function [trials] = totalTrials()
%%

% Counts trials + estimates the threshold based on the last 500 trials

% get the files list
files = dir(fullfile(sprintf('~/data/afmapb/%s/17*stim*.mat',mglGetSID)));

trials = 0;

for fi = 1:length(files)
    load(fullfile(sprintf('~/data/afmapb/%s/%s',mglGetSID,files(fi).name)));
    
    e = getTaskParameters(myscreen,task);
    e = e{1}; % why?!
    trials = trials + e.nTrials;
end

%%%%%%%%%%%%%%%%%%%%%%%
%    dispInfo    %
%%%%%%%%%%%%%%%%%%%%%%%
function dispInfo(rstimulus)
%%

% ctask = task; cscreen = myscreen; % save this incase we need them

% compute % correct for valid and invalid trials, display learning over
% time (including history from other runs)
% exp = getTaskParameters(task,myscreen);

% get the files list
files = dir(fullfile(sprintf('~/data/afmapb/%s/17*stim*.mat',mglGetSID)));

% load the files and pull out the data (long form)
%  rrun # counter #    local trial     real trial   present     resp    
%     1       2             3              4          5         6
%  contrast    correct    hit    fa    miss    cr     dead    
%      7          8       9      10     11      12      13

% get the threshold estimate from the last run
out = doStaircase('threshold',rstimulus.staircase','type=weibull','dispFig=0');
threshold = out.threshold;

%%
wn = zeros(10000,240,135);
count = 1; data = zeros(10000,13);

%figure;
for fi = 1:length(files)
    load(fullfile(sprintf('~/data/afmapb/%s/%s',mglGetSID,files(fi).name)));
    
    e = getTaskParameters(myscreen,task);
    if e{1}.nTrials>1
        e = e{1}; % why?!
    
        run = stimulus.counter;

        data(count:count+(e.nTrials-1),:) = [repmat(fi,e.nTrials,1) repmat(run,e.nTrials,1) (1:e.nTrials)' (count:count+(e.nTrials-1))' ...
            e.parameter.present' e.randVars.resp' e.randVars.contrast' ...
            e.randVars.correct' e.randVars.hit' e.randVars.fa' ...
            e.randVars.miss' e.randVars.cr' e.randVars.dead'];

        for ti = 1:e.nTrials
            timg = squeeze(stimulus.wn.img(stimulus.wn.trials{ti},:,:));
            wn(count+(ti-1),:,:) = timg;
             %imagesc(timg);colormap('gray');
             %pause(.01);
        end
        count = count+e.nTrials;
    end
end

wn = wn(1:count,:,:);
data = data(1:count,:);

l = size(data,1);
disp(sprintf('Found %01.2f%% hits %01.2f%% fa %01.2f%% miss %01.2f%% cr',sum(data(:,9))/l*100,sum(data(:,10))/l*100,sum(data(:,11))/l*100,sum(data(:,12))/l*100));

%% test
% for i = 1:size(wn,1)
%     imagesc(squeeze(wn(i,:,:)));
%     pause(.01);
% end

x = stimulus.live.X;
y = stimulus.live.Y;
d = 1;

%% Split data by hit/fa/miss/cr
% img = struct;
% name = {'hit','miss','fa','cr'};
% idx = [9 10 11 12];
% h = figure;
% for ci = 1:4
%     subplot(2,2,ci); hold on
%     img.(name{ci}) = squeeze(mean(wn(logical(data(:,idx(ci))),:,:)));
%     img.(name{ci}) = img.(name{ci})/255;
%     imagesc(img.(name{ci})');
%     set(gca,'YDir','normal');    
%     set(gca,'YTick',1:50:270,'YTickLabel',round((-134.5:50:134.5)/(myscreen.screenHeight/myscreen.imageHeight/4)));
%     set(gca,'XTick',1:100:480,'XTickLabel',round((-230.5:100:230.5)/(myscreen.screenHeight/myscreen.imageHeight/4)));
%     title(name{ci});
%     colormap('gray');
%     caxis([0 1]);
%     axis equal
%     %rectangle('Position',[x y d d],'Curvature',[1 1],'EdgeColor','w');
% end

%% 

% set a threshold offset (how much +/- the threshold we are willing to
% allow data to come from)
t_var = inf;
t_idxs = (data(:,7)>(threshold-threshold*t_var)) .* (data(:,7)<(threshold+threshold*t_var));
disp(sprintf('Using %i of %i when accounting from threshold variability',sum(t_idxs),size(t_idxs,1)));
img = struct;
name = {'hit','miss','fa','cr'};
idx = [9 10 11 12];
h = figure; hold on
for ci = 1:4
    temp = squeeze(wn(logical(t_idxs .* data(:,idx(ci))),:,:));
    temp = temp/255;
    temp = temp - 0.5;
%     temp = reshape(squeeze(std(temp)),1,240,135);
    img.(name{ci}) = temp;
end
aimg = [img.hit ; img.fa ; -img.miss ; -img.cr];
% aimg = [img.fa ; - img.miss];

% % Bootstrap takes a VERY VERY VERY long time:
% am_ = squeeze(bootci(10,@nanmean,aimg));
% am  = squeeze(mean(am_));

am = squeeze(mean(aimg));
as = squeeze(std(aimg));

imagesc(stimulus.live.X(1,:),stimulus.live.Y(:,1)',flipud(am'));

colormap('gray');
colorbar
set(gca,'Clim',[min(am(:)) max(am(:))]);
% subplot(211)
% imagesc(pm'); colorbar; axis equal
% subplot(212)
% imagesc(mm'); colorbar; axis equal

set(gca,'YDir','normal');    
title('(Hits + FA) - (Miss + CR)');

x = rstimulus.gaussX - rstimulus.gaussFWHM/2;
y = rstimulus.gaussY - rstimulus.gaussFWHM/2;
d = rstimulus.gaussFWHM;
rectangle('Position',[x y d d],'Curvature',[1 1],'EdgeColor','w');

plot([-1 1],[0 0],'-w');
plot([0 0],[-1 1],'-w');

axis equal

drawPublishAxis

%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % sets the gamma table so that we can have
% % finest possible control over the stimulus contrast.
% %
% % stimulus.colors.reservedColors should be set to the reserved colors (for cue colors, etc).
% % maxContrast is the maximum contrast you want to be able to display.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function setGammaTable(maxContrast)
% 
% global stimulus;
% 
% % set the bottom
% gammaTable(1:size(stimulus.colors.reservedBottom,1),1:size(stimulus.colors.reservedBottom,2)) = stimulus.colors.reservedBottom;
% 
% % set the gamma table
% if maxContrast == 1
%     % create the rest of the gamma table
%     cmax = 1;cmin = 0;
%     luminanceVals = cmin:((cmax-cmin)/(stimulus.colors.nUnreserved-1)):cmax;
% 
%     % now get the linearized range
%     redLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,luminanceVals,'linear');
%     greenLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,luminanceVals,'linear');
%     blueLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,luminanceVals,'linear');
% elseif maxContrast > 0
%     % create the rest of the gamma table
%     cmax = 0.5+maxContrast/2;cmin = 0.5-maxContrast/2;
%     luminanceVals = cmin:((cmax-cmin)/(stimulus.colors.nUnreserved-1)):cmax;
% 
%     % now get the linearized range
%     redLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,luminanceVals,'linear');
%     greenLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,luminanceVals,'linear');
%     blueLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,luminanceVals,'linear');
% else
%     % if we are asked for 0 contrast then simply set all the values to gray
%     redLinearized = repmat(interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,.5,'linear'),1,stimulus.colors.nUnreserved);
%     greenLinearized = repmat(interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,.5,'linear'),1,stimulus.colors.nUnreserved);
%     blueLinearized = repmat(interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,.5,'linear'),1,stimulus.colors.nUnreserved);
% end
% 
% % add to the table!
% gammaTable((stimulus.colors.mrmin:stimulus.colors.mrmax)+1,:)=[redLinearized;greenLinearized;blueLinearized]';
% 
% % set the top
% gammaTable = [gammaTable; stimulus.colors.reservedTop];
% 
% if size(gammaTable,1)~=256
%     disp('(setGammaTable) Failure: Incorrect number of colors in gamma table produced');
%     keyboard
% end
% 
% % set the gamma table
% succ = mglSetGammaTable(gammaTable);
% 
% if ~succ
%     warning('Gamma table set failure');
%     keyboard
% end
% 
% % remember what the current maximum contrast is that we can display
% stimulus.curMaxContrast = maxContrast;
% 
% 
% function initGammaTable(myscreen)
% global stimulus
% %% Gamma Table Initialization
% 
% % get gamma table
% if ~isfield(myscreen,'gammaTable')
%   stimulus.linearizedGammaTable = mglGetGammaTable;
%   disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
%   disp(sprintf('(cuecon:initGratings) No gamma table found in myscreen. Contrast displays like this'));
%   disp(sprintf('         should be run with a valid calibration made by moncalib for this monitor.'));
%   disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
% else
%   % check to make sure this matches the calibration file
%   
%   % get each channel table that should have been set by mglGetGammaTable
%   redTable = myscreen.initScreenGammaTable.redTable(:);
%   greenTable = myscreen.initScreenGammaTable.greenTable(:);
%   blueTable = myscreen.initScreenGammaTable.blueTable(:);
%   % get what the calibration structure says it should have been set to
%   gammaTable = myscreen.gammaTable(:);
%   % table values are only good to 10 bits
%   redTable = round(redTable*1024)/1024;
%   greenTable = round(greenTable*1024)/1024;
%   blueTable = round(blueTable*1024)/1024;
%   gammaTable = round(gammaTable*1024)/1024;
%   % compare, ignoring nans
%   if ~isequaln(mglGetGammaTable,myscreen.initScreenGammaTable) || ~isequaln(redTable,gammaTable) || ~isequaln(greenTable,gammaTable) || ~isequaln(blueTable,gammaTable)
%     disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
%     disp(sprintf('(curecon:initGrating) Gamma table does not match calibration'));
%     disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
%     keyboard
%   end
% end
% stimulus.linearizedGammaTable = myscreen.initScreenGammaTable;

function ovals = osum(ivals)
% compute the "other" sum, i.e. vals[i] = sum(vals[~i])
ivals(ivals==0)=1;
for i = 1:length(ivals)
    ovals(i) = sum(ivals(setdiff(1:length(ivals),i)));
end