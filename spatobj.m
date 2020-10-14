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

global stimulus resizedData resizedExamples
stimulus = struct;
stimulus.imagesLoaded = false;

%% Initialize Variables
plots = 0;
noeye = 0;
eyewindow=0; 
rebuild = 0;
test =0;
practice = 0;
pause = 0;

getArgs(varargin,{'plots=0','noeye=0','eyewindow=2','rebuild=1','test=0','practice=0','pause=0'});
stimulus.plots = plots;
stimulus.noeye = noeye;
stimulus.eyewindow = eyewindow;
stimulus.cats.remaining = setdiff(0:20,1);
stimulus.rebuild = rebuild;
stimulus.practice = practice;
stimulus.pause = pause;

clear plots noeye eyewindow load

%% Open Old Stimfile

if ~isempty(mglGetSID) && isdir(sprintf('~/data/spatobj/%s',mglGetSID))
    % Directory exists, check for a stimfile
    files = dir(sprintf('~/data/spatobj/%s/*mat',mglGetSID));

    if length(files) >= 1
        fname = files(end).name;
        s = load(sprintf('~/data/spatobj/%s/%s',mglGetSID,fname));
        % copy staircases and run numbers
        if ~s.stimulus.practice
            stimulus.cats = s.stimulus.cats;
            disp(sprintf('(spatobj) Data file: %s loaded.',fname));
        end
        clear s;
    else
        warning('(spatobj) Unable to load previous data files. If this is *not* the first run there is something wrong.');
    end
end

disp(sprintf('(spatobj) %i categories remaining (total trials %i)',length(stimulus.cats.remaining),length(stimulus.cats.remaining)*80));

%% Check practic emode
if stimulus.practice
    disp(sprintf('(spatobj) PRACTICE MODE\n(spatobj)Loading only tree category'));
    stimulus.cats.remaining = 1;
end

%% Load iamges

if stimulus.rebuild
    % load data images
    data = load('~/data/spatobj/stimuli.mat');
    data = data.data;
    
    % meta data columns
    %   1        2               3
    % idx  target loc        target cat
    [metaHeader, metaData] = csvreadh('~/data/spatobj/stimuli_meta_ints.csv');
    metaData = metaData(1:size(metaData,1)-1,:);
    metaData(:,1) = metaData(:,1)+1; % add 1 for matlab indexing
    
    % NOTE: for now don't load the headers or data, because the data
    % contain strings and csvreadh can't acutally handle that
    
    examples = load('~/data/spatobj/exemplars.mat');
    stimulus.cats.categories = examples.cat_map;
    examples = examples.data;
    
    exData = zeros(size(examples,1),2);
    % exData columns
    %   1        2
    % idx    category
    exData(:,1) = 1:size(exData,1);
    temp = repmat(0:20,5,1);
    exData(:,2) = temp(:);
end

%% Rename some categories
orig = {'bakery','garden spider','long-horned beetle','analog clock','coffee','Ferris wheel'};
new = {'baked goods','spider','beetle','clock','coffee beans','ferris wheel'};

for oi = 1:length(orig)
    idx = cellfun(@(x) strcmp(x,orig{oi}),stimulus.cats.categories);
    if any(idx)
        disp(sprintf('(spatobj) Renaming %s to %s',orig{oi},new{oi}));
        stimulus.cats.categories{find(idx,1)} = new{oi};
    end
end

%% Setup Screen
% warning('using test screen');
if test
    myscreen = initScreen('test');
else
    myscreen = initScreen('VPixx');
end

% set background to black
myscreen.background = 0;
stimulus.eyeFrames = myscreen.framesPerSecond * 0.300; % eye movements occur when for 300 ms someone moves out of the fixation region

%% Reset the gamma table to be blank:
if ~test
    disp('(spatobj) Resetting gamma table to be blank');
    gammaTable(:,1) = (0:1/255:1);
    gammaTable(:,2) = (0:1/255:1);
    gammaTable(:,3) = (0:1/255:1);
    mglSetGammaTable(gammaTable);
end

%% Initialize Stimulus
myscreen.stimulusNames{1} = 'stimulus';
stimulus.responseKeys = [1 2]; % present / absent

stimulus.colors.green = [0 1 0];
stimulus.colors.red = [1 0 0];
stimulus.colors.white = [1 1 1];
stimulus.colors.black = [0 0 0];

stimulus.live = struct;
stimulus.live.last = struct;
stimulus.live.last.done = false;

%% Sizes
stimulus.arrayWidth = round(9*myscreen.screenWidth/myscreen.imageWidth); % in pixels
disp(sprintf('Images will be resized to %1.2f pixels to match 9 degrees',stimulus.arrayWidth));

if ~stimulus.rebuild && ~isempty(resizedData) && ~isempty(resizedExamples) && size(resizedData,2)==stimulus.arrayWidth
    disp('Resized images are being loaded -- if this is the first run for this subject, [esc] and clear all');
    data = resizedData;
    examples = resizedExamples;
else
    disp('Resizing images');
    dataTemp = zeros(size(data,1),stimulus.arrayWidth,stimulus.arrayWidth,size(data,4));
    examplesTemp = zeros(size(examples,1),stimulus.arrayWidth,stimulus.arrayWidth,size(examples,4));

    for i = 1:size(data,1)
        dataTemp(i,:,:,:) = imresize(squeeze(data(i,:,:,:)),[stimulus.arrayWidth, stimulus.arrayWidth]);
    end
    for i = 1:size(examples,1)
        examplesTemp(i,:,:,:) = imresize(squeeze(examples(i,:,:,:)),[stimulus.arrayWidth, stimulus.arrayWidth]);
    end

    data = dataTemp;
    examples = examplesTemp;
    
    resizedData = uint8(data);
    resizedExamples = uint8(examples);
end

    
stimulus.live.metaHeader = metaHeader;
stimulus.live.metaData = metaData;
stimulus.live.data = uint8(data);
stimulus.live.exHeader = {'index','category'};
stimulus.live.exData = exData;
stimulus.live.examples = uint8(examples);

clear metaHeader metaData exData examples data dataTemp examplesTemp

%% Plot and return
if stimulus.plots==2
    dispInfo;
    return
end

%% Setup Examples Phase

exPhase = struct;
exPhase.segmin = [inf inf inf inf inf inf inf];
exPhase.segmax = [inf inf inf inf inf inf inf];
exPhase.getResponse = ones(size(exPhase.segmin));

exPhase.waitForBacktick = 0;
exPhase.numTrials = 1;
exPhase.parameter.targetCategory = nan; % 0->19
exPhase.randVars.calculated.exampleNum = nan; % 1->5
exPhase.random = 0;

%% Setup Task Phase

% task waits for fixation on first segment
stimulus.seg.iti = 1;
stimulus.seg.fix = 2;
stimulus.seg.cue = 3;
stimulus.seg.isi = 4;
stimulus.seg.stim = 5;
stimulus.seg.mask = 6;
stimulus.seg.resp = 7; % resp and stim overlap

taskPhase = struct;

taskPhase.segmin = [0.25 inf 0.75 0.75  inf 1 1];
taskPhase.segmax = [0.25 inf 0.75 0.75  inf 1 1];

if stimulus.noeye
    taskPhase.segmin(stimulus.seg.fix) = 0.5;
    taskPhase.segmax(stimulus.seg.fix) = 0.5;
end

taskPhase.waitForBacktick = 0;

taskPhase.getResponse = zeros(1,length(taskPhase.segmin));
taskPhase.getResponse(stimulus.seg.stim) = 1;
taskPhase.getResponse(stimulus.seg.mask) = 1;
taskPhase.getResponse(stimulus.seg.resp) = 1;

taskPhase.numTrials = 80; % there are 20 images with the target and 20 without, for each of focal and distributed

taskPhase.random = 1;

taskPhase.parameter.imageNumber = 1:20;
taskPhase.parameter.targetPresent = [0 1];
taskPhase.parameter.targetCategory = nan;
taskPhase.parameter.targetPosition = nan;
taskPhase.parameter.focal = [0 1];
taskPhase.randVars.calculated.duration = nan; % we will randomize duration from 50 ms to 150 ms;
taskPhase.randVars.calculated.dead = nan;
taskPhase.randVars.calculated.responsePresent = nan;
taskPhase.randVars.calculated.correct = nan;
taskPhase.randVars.calculated.imgCat0 = nan;
taskPhase.randVars.calculated.imgCat1 = nan;
taskPhase.randVars.calculated.imgCat2 = nan;
taskPhase.randVars.calculated.imgCat3 = nan;

%% Set up all phases
catIdxOrder = randperm(length(stimulus.cats.remaining));

task{1} = {};

for i = 1:length(stimulus.cats.remaining)
    temp = taskPhase;
    temp.parameter.targetCategory = stimulus.cats.remaining(catIdxOrder(i));
    
    exTemp = exPhase;
    exTemp.parameter.targetCategory = temp.parameter.targetCategory;
    
    % now add the example phase and task phase
    task{1}{end+1} = exTemp;
    task{1}{end+1} = temp;
end

% set the first exPhase to wait for backtick
task{1}{1}.waitForBacktick = 1;

%% Full Setup
% Initialize task phases
for phaseNum = 1:2:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@exampleSegmentCallback,@exampleUpdateCallback,@exampleResponseCallback,@exampleTrialCallback,[],[]);
end
for phaseNum = 2:2:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,@endTrialCallback,[]);
end

%% EYE CALIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~stimulus.noeye
    myscreen = eyeCalibDisp(myscreen);
end

%% Main Task Loop

phaseNum = 1;
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

if stimulus.live.last.done
    %we just finished a category, save it
    stimulus.cats.remaining = setdiff(stimulus.cats.remaining,stimulus.live.last.category);
end

% remove extra fields from stimulus to prevent file size from being huge
stimulus = rmfield(stimulus,'live');

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

% if stimulus.practice
%     disp('(spatobj) Displaying plots');
%     dispInfo(stimulus);
% end

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

% function dispInfo(stimulus)
%%

function [task, myscreen] = exampleSegmentCallback(task,myscreen)
global stimulus

% if this is seg #1 do nothing, any other segment we get the image
if task.thistrial.thisseg > 1
    % build the current exemplar image
    imgIdx = task.thistrial.targetCategory * 5 + task.thistrial.thisseg-1;
    stimulus.live.exemplar = loadExemplar(imgIdx);
end

function [task, myscreen] = exampleTrialCallback(task,myscreen)
%%
global stimulus

if stimulus.live.last.done
    %we just finished a category, save it
    stimulus.cats.remaining = setdiff(stimulus.cats.remaining,stimulus.live.last.category);
end

stimulus.live.last = struct;
stimulus.live.last.done = false;

function [task, myscreen] = exampleUpdateCallback(task,myscreen)
%%
global stimulus
mglClearScreen();

if task.thistrial.thisseg == 1
    mglTextDraw('Search for: ',[0 1]);
    mglTextDraw(stimulus.cats.categories{task.thistrial.targetCategory+1},[0 -1]);
elseif task.thistrial.thisseg == 7  
    mglTextDraw('Ready?',[0 0]);
else
    mglBltTexture(stimulus.live.exemplar,[0 0]);
end

function [task, myscreen] = exampleResponseCallback(task,myscreen)

% when the user clicks, go to the next screen
task = jumpSegment(task);

function [task, myscreen] = endTrialCallback(task,myscreen)
global stimulus

if task.trialnum == task.numTrials
    stimulus.live.last.done = true;
end

function [task, myscreen] = startTrialCallback(task,myscreen)
%%
global stimulus

if task.trialnum == 1
    % set the last info, so this block will get saved
    stimulus.live.last = struct;
    stimulus.live.last.done = false;
    stimulus.live.last.category = task.thistrial.targetCategory;
end

% set the duration of this trial
if stimulus.pause
    task.thistrial.duration = inf;
else
    task.thistrial.duration = randsample([1 2 4 8 16 32],1);
end
stimulus.live.framesRemaining = task.thistrial.duration;

% stim segment is inf, so we have to jumpSegment when the frames are
% elapsed

% build the stimulus image for this trial

% first get images of this target category
cData = stimulus.live.metaData(stimulus.live.metaData(:,3)==task.thistrial.targetCategory,:);
if task.thistrial.targetPresent
    cData = cData(cData(:,2)>=0,:);
else
    cData = sel(cData,2,-1);
end

if task.thistrial.focal
    cData = cData(1:20,:);
    % we have to cheat here because Kai didn't include a
    % focal/distributed column, but it's organized in the right order so it
    % should be fine. If the images get re-generated into a different order
    % this would have to be changed.
else
    cData = cData(21:end,:);
end

cData = cData(task.thistrial.imageNumber,:);

if task.thistrial.focal
    if task.thistrial.targetPresent
        task.thistrial.targetPosition = cData(2); % this will be -1 when target is absent
    else
        % on focal trials with the target not present we have a problem
        % targetPosition is -1, but we still want to show a cue
        if cData(2)~=-1
            warning('serious error in data');
            keyboard
        end
        task.thistrial.targetPosition = floor(mod(cData(1),20)/5);
    end
else
    task.thistrial.targetPosition = cData(2);
end
task.thistrial.imgCat0 = cData(4);
task.thistrial.imgCat1 = cData(5);
task.thistrial.imgCat2 = cData(6);
task.thistrial.imgCat3 = cData(7);

imgIdx = cData(1);

[stimulus.live.tex, stimulus.live.mask] = loadData(imgIdx);

stimulus.live.fixColor = stimulus.colors.white;

task.thistrial.dead = false;
stimulus.live.fixCount = 0;
stimulus.live.eyeCount = 0;

if task.thistrial.targetPresent
    disp(sprintf('(spatobj) There is a %s in the array.',stimulus.cats.categories{task.thistrial.targetCategory+1}));
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

if task.thistrial.thisseg==stimulus.seg.stim
    if stimulus.live.framesRemaining > 0
        mglBltTexture(stimulus.live.tex,[0 0]);
        stimulus.live.framesRemaining = stimulus.live.framesRemaining - 1;
    else
        mglBltTexture(stimulus.live.mask,[0 0]);
    end 
else
    mglBltTexture(stimulus.live.mask,[0 0]);
end

mglGluDisk(0,0,0.5,stimulus.colors.black);
mglFixationCross(0.5,0.5,stimulus.live.fixColor);

% if this is the cue stimulus, color white in the directions to be attended
if task.thistrial.thisseg==stimulus.seg.cue
    if task.thistrial.focal
        switch task.thistrial.targetPosition
            case 0
                mglLines2(0,0,-0.25,0.25,1,stimulus.colors.white);
            case 1
                mglLines2(0,0,0.25,0.25,1,stimulus.colors.white);
            case 2
                mglLines2(0,0,-0.25,-0.25,1,stimulus.colors.white);
            case 3
                mglLines2(0,0,0.25,-0.25,1,stimulus.colors.white);
        end
    else
        mglLines2(-0.25,-0.25,0.25,0.25,1,stimulus.colors.white);
        mglLines2(-0.25,0.25,0.25,-0.25,1,stimulus.colors.white);
    end
end

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

function tex = loadExemplar(idx)
global stimulus
img = permute(imrotate(squeeze(stimulus.live.examples(idx,:,:,:)),-90),[3 1 2]);
img(4,:,:) = 255;
tex = mglCreateTexture(img);

function [tex, mask] = loadData(idx)  
global stimulus
img = permute(imrotate(squeeze(stimulus.live.data(idx,:,:,:)),-90),[3 1 2]);
mask = img;
img(4,:,:) = 255;
tex = mglCreateTexture(img);
mask = mask(:);
mask = mask(randperm(length(mask)));
mask = reshape(mask,3,stimulus.arrayWidth,stimulus.arrayWidth);
mask(4,:,:) = 255;
mask = mglCreateTexture(mask);