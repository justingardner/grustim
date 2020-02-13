function [ myscreen ] = spatobj( varargin )
% SPATOBJ (spatial attention object detection)
% *** set 'noeye=1' to turn off the eye tracker***
%
% This code is for a psychophysics experiment based on a recent paper
% [Lindsay & Miller, 2018] looking at whether a CNN architecture can
% replicate simple results from humans. 
%
% REFERENCES
% Lindsay, G. W., & Miller, K. D. (2018). How biological attention mechanisms improve task performance in a large-scale visual system model. eLife, 7, e38105.
% 


%% Over-write the stimulus struct

global stimulus 
stimulus = struct;
stimulus.imagesLoaded = false;

%% Initialize Variables
plots = 0;
noeye = 0;
eyewindow=0; 
rebuild = 0;

getArgs(varargin,{'plots=0','noeye=0','eyewindow=2','rebuild=1'});
stimulus.plots = plots;
stimulus.noeye = noeye;
stimulus.eyewindow = eyewindow;
stimulus.cats.focal = 1:20;
stimulus.cats.distributed = 1:20;
stimulus.rebuild = rebuild;
stimulus.categories = {'artichoke','bathtub','cabbage butterfly','computer','mortar','greenhouse','padlock','toaster','paintbrush','football helmet','horned beetle','home theatre','stone wall','baked goods','coffee','spider','banana','clock','ferris wheel','seashore'};

clear plots noeye eyewindow load

%% Load iamges

if stimulus.rebuild
    % load data images
    data = load('~/data/spatobj/stimuli_uint8.mat');
    data = data.data;
    
    metaData = zeros(size(data,1),3);
    metaData(:,1) = 1:size(metaData,1);
    temp = repmat(1:20,80,1);
    metaData(:,2) = temp(:);
    temp = repmat([zeros(1,40) ones(1,40)]',1,20);
    metaData(:,3) = temp(:);
    % meta data columns
    %   1        2               3
    % idx  target categ  distractors_only
    
    % load data csv
    
    % NOTE: for now don't load the headers or data, because the data
    % contain strings and csvreadh can't acutally handle that
    %[metaHeader, metaData] = csvreadh('~/data/spatobj/stimuli_meta.csv');
    %[exHeader, exData] = csvreadh('~/data/spatobj/exemplar_meta.csv');
    
    examples = load('~/data/spatobj/exemplars_uint8.mat');
    examples = examples.data;

    exData = zeros(size(examples,1),2);
    % exData columns
    %   1        2
    % idx    category
    exData(:,1) = 1:size(exData,1);
    temp = repmat(1:20,5,1);
    exData(:,2) = temp(:);
end

%% Open Old Stimfile

if ~stimulus.rebuild && ~isempty(mglGetSID) && isdir(sprintf('~/data/spatobj/%s',mglGetSID))
    % Directory exists, check for a stimfile
    files = dir(sprintf('~/data/spatobj/%s/1*mat',mglGetSID));

    if length(files) >= 1
        fname = files(end).name;
        s = load(sprintf('~/data/spatobj/%s/%s',mglGetSID,fname));
        % copy staircases and run numbers
        stimulus.remainCategory = s.stimulus.remainCategory;
        clear s;
        disp(sprintf('(spatobj) Data file: %s loaded.',fname));
    else
        warning('(spatobj) Unable to load previous data files. If this is *not* the first run there is something wrong.');
    end
end

disp(sprintf('(spatobj) %i categories remaining (total trials %i)',length(stimulus.remainCategory),length(stimulus.remainCategory)*30));

%% Setup Screen
myscreen = initScreen('test');
% myscreen = initScreen('VPixx');
% set background to black
myscreen.background = 0;
stimulus.eyeFrames = myscreen.framesPerSecond * 0.300; % eye movements occur when for 300 ms someone moves out of the fixation region

%% Reset the gamma table to be blank:
% should this even happen?
% gammaTable(:,1) = (0:1/255:1);
% gammaTable(:,2) = (0:1/255:1);
% gammaTable(:,3) = (0:1/255:1);
% mglSetGammaTable(gammaTable);
%% Sizes
stimulus.arrayWidth = round(9*myscreen.screenWidth/myscreen.imageWidth); % in pixels
disp(sprintf('Images will be resized to %1.2f pixels to match 9 degrees',stimulus.arrayWidth));

% resize all of the images
stop = 1;

dataTemp = zeros(size(data,1),stimulus.arrayWidth,stimulus.arrayWidth,size(data,4));
examplesTemp = zeros(size(examples,1),stimulus.arrayWidth,stimulus.arrayWidth,size(examples,4));

for i = 1:size(data,1)
    dataTemp(i,:,:,:) = imresize(squeeze(data(i,:,:,:)),[stimulus.arrayWidth, stimulus.arrayWidth]);
    examplesTemp(i,:,:,:) = imresize(squeeze(data(i,:,:,:)),[stimulus.arrayWidth, stimulus.arrayWidth]);
end

data = dataTemp;
examples = examplesTemp;

clear dataTemp examplesTemp

%% Plot and return
if stimulus.plots==2
    dispInfo;
    return
end

%% Initialize Stimulus
myscreen.stimulusNames{1} = 'stimulus';
stimulus.responseKeys = [1 2]; % present / absent

stimulus.colors.green = [0 1 0];
stimulus.colors.red = [1 0 0];
stimulus.colors.white = [1 1 1];
stimulus.colors.black = [0 0 0];

stimulus.live = struct;

%% Setup Examples Phase

exPhase = struct;
exPhase.seglen = [inf inf inf inf inf inf];
exPhase.waitForBackTick = 0;
exPhase.numTrials = 1;
exPhase.randVars.calculated.exampleNum = nan; % 1->5
exPhase.randVars.calculated.targetCategory = nan; % 1->20

%% Setup Task Phase

% task waits for fixation on first segment
stimulus.seg.iti = 1;
stimulus.seg.fix = 2;
stimulus.seg.stim = 3;
stimulus.seg.mask = 4;
stimulus.seg.resp = 5; % resp and stim overlap

taskPhase = struct;

taskPhase.segmin = [0.5 inf inf 1 1];
taskPhase.segmax = [1.5 inf inf 1 1];

if stimulus.noeye
    taskPhase.segmin(stimulus.seg.fix) = 0.5;
    taskPhase.segmax(stimulus.seg.fix) = 0.5;
end

taskPhase.waitForBacktick = 0;

taskPhase.getResponse = zeros(1,length(taskPhase.segmin));
taskPhase.getResponse(stimulus.seg.mask) = 1;
taskPhase.getResponse(stimulus.seg.resp) = 1;

taskPhase.numTrials = 40; % there are 20 images with the target and 20 without

taskPhase.random = 1;

taskPhase.parameter.imageNumber = 1:20;
taskPhase.parameter.targetPresent = [0 1];
taskPhase.randVars.calculated.duration = nan; % we will randomize duration from 50 ms to 150 ms;
taskPhase.randVars.calculated.targetCategory = nan;
taskPhase.randVars.calculated.dead = nan;
taskPhase.randVars.calculated.responsePresent = nan;
taskPhase.randVars.calculated.focal = nan;
taskPhase.randVars.calculated.correct = nan;
taskPhase.randVars.calculated.imgCat1 = nan;
taskPhase.randVars.calculated.imgCat2 = nan;
taskPhase.randVars.calculated.imgCat3 = nan;
taskPhase.randVars.calculated.imgCat4 = nan;

%% Set up all phases
for i = 1:length(stimulus.cats.focal)
    % add a phase for each focal category remaining
    task{1}{end+1} = exPhase;
    task{1}{end+1} = taskPhase;
    if i==1
        % if this is the first phase pair, set backtick to 1
        task{1}{1}.waitForBacktick = 1;
        task{1}{2}.waitForBacktick = 1;
    end
end

for i = 1:length(stimulus.cats.distributed)
    task{1}{end+1} = exPhase;
    task{1}{end+1} = taskPhase;
end

%% Full Setup
% Initialize task phases
for phaseNum = 1:2:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,[],@exampleUpdateCallback,[],@exampleTrialCallback,[],[]);
end
for phaseNum = 2:2:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);
end

%% EYE CALIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~stimulus.noeye
    myscreen = eyeCalibDisp(myscreen);
end

%% Put up the current target category
for i= 1:2
    mglClearScreen;
    mglFixationCross;
    mglFlush
end

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
mglClearScreen;
mglFlush
myscreen.flushMode = 1;

% track finished categories
r = input('Was the last block completed? [enter to confirm]','s');
if ~isempty(r)
    stimulus.doneCategories(end) = [];
end
stimulus.remainCategory = setdiff(stimulus.remainCategory,stimulus.doneCategories);

if ~stimulus.rebuild
    % remove extra fields
    stimulus = rmfield(stimulus,'images');
    stimulus = rmfield(stimulus,'distractors');
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

if stimulus.plots
    disp('(spatobj) Displaying plots');
    dispInfo(stimulus);
end

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

function dispInfo()
%%

function [task, myscreen] = startBlockCallback(task,myscreen)
global stimulus
% set the new category
stimulus.live.curCategory = randsample(stimulus.remainCategory,1);
if ~isfield(stimulus,'doneCategories')
    stimulus.doneCategories = [];
end
stimulus.doneCategories = [stimulus.doneCategories stimulus.live.curCategory];

disp(sprintf('Starting new block, category: %s',stimulus.categories{stimulus.live.curCategory}));




function [task, myscreen] = startTrialCallback(task,myscreen)
global stimulus

task.thistrial.category = stimulus.live.curCategory;


task.thistrial.duration = rand*0.75 + 0.25;
task.thistrial.seglen(stimulus.seg.stim) = task.thistrial.duration;

% create the texture for this trial
if task.thistrial.targetPresent
    img = permute(squeeze(stimulus.images(task.thistrial.category,task.thistrial.targetImg,:,:,:)),[3 1 2]);
    for i = 1:4
        task.thistrial.(sprintf('imgCat%i',i)) = stimulus.images_info(task.thistrial.category,task.thistrial.targetImg,i);
    end
else
    img = permute(squeeze(stimulus.distractors(task.thistrial.category,task.thistrial.targetImg,:,:,:)),[3 1 2]);
    for i = 1:4
        task.thistrial.(sprintf('imgCat%i',i)) = stimulus.distractors_info(task.thistrial.category,task.thistrial.targetImg,i);
    end
end
mask = img;
img(4,:,:) = 255;
stimulus.live.tex = mglCreateTexture(img);
mask = mask(:);
mask = mask(randperm(length(mask)));
mask = reshape(mask,3,stimulus.arrayWidth,stimulus.arrayWidth);
mask(4,:,:) = 255;
stimulus.live.mask = mglCreateTexture(mask);

stimulus.live.fixColor = stimulus.colors.white;

task.thistrial.dead = false;
stimulus.live.fixCount = 0;
stimulus.live.eyeCount = 0;

if task.thistrial.targetPresent
    disp(sprintf('(spatobj) There is a %s in the array.',stimulus.categories{task.thistrial.category}));
else
    disp(sprintf('(spatobj) No target'));
end

function [task, myscreen] = startSegmentCallback(task,myscreen)
%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus

mglClearScreen();
if task.thistrial.dead && mglGetSecs(stimulus.live.deadTime)>1
    task = jumpSegment(task,inf);
end

% skip screen updates if you are already dead
if task.thistrial.dead
    if task.thistrial.dead
        mglTextSet([],32,stimulus.colors.red);
        mglTextDraw('Eye Movement Detected',[0 0]);
    end
    return
end


if task.thistrial.thisseg==stimulus.seg.iti
    % if this is the first trial
    if mod(task.trialnum,30)==1
        task.thistrial.seglen(stimulus.seg.iti) = 3;
        mglTextDraw(sprintf('Look for: %s',stimulus.categories{stimulus.live.curCategory}),[0 1.5]);
    end
else
    if task.thistrial.thisseg==stimulus.seg.stim
        mglBltTexture(stimulus.live.tex,[0 0]);
    else
        mglBltTexture(stimulus.live.mask,[0 0]);
    end
end

mglGluDisk(0,0,0.5,stimulus.colors.black);
mglFixationCross(0.5,0.5,stimulus.live.fixColor);


% do eye position tracking, but only during some segments
if (~stimulus.noeye) && (stimulus.eyewindow>0) && any(task.thistrial.thisseg==[stimulus.seg.fix stimulus.seg.stim])
    % check eye pos
    [pos,~] = mglEyelinkGetCurrentEyePos;

    % compute distance
    dist = hypot(pos(1),pos(2));

    if task.thistrial.thisseg==stimulus.seg.fix
        if stimulus.live.fixCount > stimulus.eyeFrames
            task = jumpSegment(task);
        elseif ~any(isnan(pos))
            if dist < stimulus.eyewindow
                stimulus.live.fixCount = stimulus.live.fixCount + 1;
            else
                stimulus.live.fixCount = 0;
            end
        end
    else
        % Eye movement detection code
        if (~stimulus.noeye) && (stimulus.eyewindow>0) && ~task.thistrial.dead
            if ~any(isnan(pos))
                if dist > stimulus.eyewindow && stimulus.live.eyeCount > stimulus.eyeFrames
                    disp('Eye movement detected!!!!');
                    stimulus.live.deadTime = mglGetSecs;
                    task.thistrial.dead = 1;
                    return
                elseif dist > stimulus.eyewindow
                    stimulus.live.eyeCount = stimulus.live.eyeCount + 1;
                end
            end
        end
    end

end

function [task, myscreen] = getResponseCallback(task, myscreen)
global stimulus

responseText = {'Incorrect','Correct'};
targetText = {'absent','present'};

if any(task.thistrial.whichButton == stimulus.responseKeys)
    if task.thistrial.gotResponse == 0
        task.thistrial.responsePresent = task.thistrial.whichButton==1;
        task.thistrial.correct = task.thistrial.targetPresent==task.thistrial.responsePresent;

        disp(sprintf('(spatobj) Response %s, target was %s',responseText{task.thistrial.correct+1},targetText{task.thistrial.targetPresent+1}));
        if task.thistrial.correct
            stimulus.live.fixColor = stimulus.colors.green;
        else
            stimulus.live.fixColor = stimulus.colors.red;
        end
    end
end

if task.thistrial.thisseg<stimulus.seg.resp
    task = jumpSegment(task,stimulus.seg.resp);
end
