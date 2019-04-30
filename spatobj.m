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
% EXPERIMENT CALL:
% scanning
% afcom('scan=1');
% TESTING CALL:
% afcom('cue=#','noeye=1','powerwheel=0');

%% Over-write the stimulus struct

global stimulus 
stimulus = struct;
stimulus.imagesLoaded = false;

%% Initialize Variables
plots = 0;
noeye = 0;
eyewindow=0; 

getArgs(varargin,{'scan=0','plots=0','noeye=0','eyewindow=1.5'});
stimulus.plots = plots;
stimulus.noeye = noeye;
stimulus.eyewindow = eyewindow;
stimulus.remainCategory = 1:20;


clear plots noeye eyewindow

%% Open Old Stimfile

if ~isempty(mglGetSID) && isdir(sprintf('~/data/spatobj/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/spatobj/%s/1*mat',mglGetSID));

    if length(files) >= 1
        fname = files(end).name;

        s = load(sprintf('~/data/spatobj/%s/%s',mglGetSID,fname));
        % copy staircases and run numbers
        stimulus.images = s.stimulus.images;
        stimulus.categories = s.stimulus.categories;
        stimulus.images_info = s.stimulus.images_info;
        stimulus.distractors = s.stimulus.distractors;
        stimulus.distractors_info = s.stimulus.distractors_info;
        stimulus.imagesLoaded = s.stimulus.imagesLoaded;
        stimulus.remainCategory = s.stimulus.remainCategory;
        clear s;
        disp(sprintf('(spatobj) Data file: %s loaded.',fname));
    else
        warning('(spatobj) Unable to load previous data files. If this is *not* the first run there is something wrong.');
    end
end

disp(sprintf('(spatobj) %i categories remaining (total trials %i)',length(stimulus.remainCategory),length(stimulus.remainCategory)*30));

%% Load stimulus
if ~stimulus.imagesLoaded
    loadStimulus();
end

%% Setup Screen
myscreen = initScreen('VPixx');
% set background to black
myscreen.background = 0;

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

%% Sizes
stimulus.arrayWidth = myscreen.screenWidth/myscreen.imageWidth*7; % IN PIXELS! We'll resize images to this... 

%% Setup Task

% task waits for fixation on first segment
stimulus.seg.iti = 1;
stimulus.seg.fix = 2;
stimulus.seg.stim = 3;
stimulus.seg.mask = 4;
stimulus.seg.resp = 5; % resp and stim overlap

task{1}{1}.segmin = [0 inf inf 0.5 1];
task{1}{1}.segmax = [1 inf inf 0.5 1];

if stimulus.noeye
    task{1}{1}.segmin(stimulus.seg.fix) = 0;
    task{1}{1}.segmax(stimulus.seg.fix) = 0;
end

task{1}{1}.waitForBacktick = 1;

task{1}{1}.getResponse = zeros(1,length(task{1}{1}.segmin));
task{1}{1}.getResponse(stimulus.seg.stim) = 1;
task{1}{1}.getResponse(stimulus.seg.mask) = 1;
task{1}{1}.getResponse(stimulus.seg.resp) = 1;

task{1}{1}.numTrials = 30*20;

task{1}{1}.random = 1;

task{1}{1}.parameter.targetImg = 1:15;
task{1}{1}.parameter.targetPresent = [0 1];
task{1}{1}.randVars.calculated.duration = nan; % we will randomize duration from 50 ms to 150 ms;
task{1}{1}.randVars.calculated.targetCategory = nan;
task{1}{1}.randVars.calculated.dead = nan;
task{1}{1}.randVars.calculated.responsePresent = nan;
task{1}{1}.randVars.calculated.correct = nan;
task{1}{1}.randVars.calculated.imgCat1 = nan;
task{1}{1}.randVars.calculated.imgCat2 = nan;
task{1}{1}.randVars.calculated.imgCat3 = nan;
task{1}{1}.randVars.calculated.imgCat4 = nan;

task{1}{1}.synchToVol = zeros(1,length(task{1}{1}.segmin));

%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],@startBlockCallback);
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
stimulus.remainCategory = setdiff(stimulus.remainCategory,stimulus.live.curCategory);

function [task, myscreen] = startTrialCallback(task,myscreen)
global stimulus

task.thistrial.category = stimulus.live.curCategory;


task.thistrial.duration = rand*2 + 0.5;
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
mask = reshape(mask,3,224,224);
mask(4,:,:) = 255;
stimulus.live.mask = mglCreateTexture(mask);

stimulus.live.fixColor = stimulus.colors.white;

task.thistrial.dead = false;

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

switch task.thistrial.thisseg
    case stimulus.seg.iti
        % if this is the first trial
        if mod(task.trialnum,30)==1
            task.thistrial.seglen(stimulus.seg.iti) = 5;
            mglTextDraw(sprintf('Look for: %s',stimulus.categories{stimulus.live.curCategory}),[0 1.5]);
        end
    case stimulus.seg.stim
        mglBltTexture(stimulus.live.tex,[0 0]);
    case stimulus.seg.mask
        mglBltTexture(stimulus.live.mask,[0 0]);
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
        
function loadStimulus()

global stimulus

% load the file which tracks which images are shown to which participants
trackFile = '~/data/spatobj/tracking.mat';
if isfile(trackFile)
    track = load(trackFile);
    track = track.track;
else
    track = struct;
    track.todo = 1:95;
    track.done = {};
end
% mygroup is the number from 1:95 which corresponds to the images which
% will be shown to this participant, in theory if we run for all 95 groups
% we can match the dataset with the CNN dataset
mygroup = randsample(track.todo,1);
track.todo = setdiff(track.todo,mygroup);
save(trackFile,'track');

% for each category we will show all 15 examples once, and we will also
% show one random image from a different category. This means we need to
% load at least 15 extra random images from each category, although they
% won't all get used. 
aimages = zeros(20,15,224,224,3,'uint8');
aimages_info = zeros(20,15,4,'uint8');

adistractors = zeros(20,15,224,224,3,'uint8');
adistractors_info = zeros(20,15,4,'uint8');

%% distractors
% the first thing to do is load all of the info files, and concatenate them
% into one data array. We need to do this to search for distractor images
% which don't include each of the categories. 
ainfo = nan(20*15*20,7);
count = 1;

for cat = 0:19
    info = readNPY(sprintf('~/data/spatobj/images/dat%i_info.npy',cat));
    
    for ci = 1:15
        for gi = 1:95
            ainfo(count,:) = [cat squeeze(info(ci,gi,:))' ci gi];
            count = count+1;
        end
    end
end
clear info

% now we will go through each category (1:20) and pick 15 distractor images
% from the other categories which don't have this category in them
distractorList = zeros(20,15,2); % this will be the category # and group # of the distractor image

for cat = 0:19
    % remove all rows which include this category
    info = ainfo(~any(ainfo(:,1:5)==cat,2),:);
    % pick 15 random categories
    myDistCats = randsample(setdiff(0:19,cat),15);
    % for each category, pick a random group from the available options
    for di = 1:length(myDistCats)
        dcat = myDistCats(di);
        % check which groups are available
        cinfo = info(info(:,1)==dcat,:);
        % get 15 random images
        idxs = randsample(1:size(cinfo,1),1);
        % save these
        distractorList(cat+1,di,:) = cinfo(idxs,6:7);
    end
end

%% now load the actual images

disppercent(-1/20);
for cat = 0:19
    dat = readNPY(sprintf('~/data/spatobj/images/dat%i.npy',cat));
    info = readNPY(sprintf('~/data/spatobj/images/dat%i_info.npy',cat));
    
    % get all the images for my group
    images = squeeze(dat(:,mygroup,:,:,:));
    images_info = squeeze(info(:,mygroup,:));
    
    aimages(cat+1,:,:,:,:) = uint8(images);
    aimages_info(cat+1,:,:) = uint8(images_info);
    
    % go through the distractors and pull any images that are needed
    for dc = 0:19
        for di = 1:15
            if distractorList(dc+1,di,1)==cat
                adistractors(dc+1,di,:,:,:) = dat(distractorList(dc+1,di,1),distractorList(dc+1,di,2),:,:,:);
                adistractors_info(dc+1,di,:) = info(distractorList(dc+1,di,1),distractorList(dc+1,di,2),:);
            end
        end
    end
    
    disppercent(cat/20);
    
    clear dat
    clear info
end
disppercent(inf);

%% and save
stimulus.categories = {'paintbrush','wall clock','seashore','paddlewheel','padlock','garden spider','long-horned beetle','cabbage butterfly','toaster','greenhouse','bakery','stone wall','artichoke','modem','football helmet','stage','mortar','soup','dough','bathtub'};
stimulus.imagesLoaded = true;
stimulus.images = aimages;
stimulus.images_info = aimages_info;
stimulus.distractors = adistractors;
stimulus.distractors_info = adistractors_info;

%% test code (to check that images actually have their category)
% figure;`2221211221212111122112121222
% for cat = 1:20
%     % pick a random image
%     idx = randi(15);
%     % display
%     imagesc(squeeze(aimages(cat,idx,:,:,:)));
%     % tell me what this is
%     r = input(sprintf('Press [enter] if you see a %s',stimulus.categories{cat}),'s');
%     if ~isempty(r)
%         disp(aimages_info(cat,idx,:));
%         warning(sprintf('Category %i index %i',cat,idx));
%     end
% end