function [ myscreen ] = fbsear_pilot( varargin )
%
% FEATURE BASED SEARCH 
%    Feature based search in natural images using the COCO dataset. This
%    version runs the fMRI code. Each stimulus is 1.5 s. ITI is randomized
%    to average 2-7 s (avg 4.5). Subject is blocked each minute to attend
%    to a specific feature category (dog?, car?, people?, fixate). The
%    different tasks run for one minute before changing. Observers respond
%    on each trial whether the target was present or absent
%
%    Usage: fbsear(varargin)
%    Authors: Dan Birman
%    Date: 12/08/2017
%
%    Parameters:
%
%    You should use *736* volumes (60 trials * 6 secs * 2 + 16) you will
%    have to manually end each run by hitting ESC (if you use 736 a volume
%    is acquired after the screen goes dark, so don't do that)
%
%      On the first run call:
%       fbsear_pilot('run=1','genTex=1');
%
%      All subsequent runs:
%       fbsear_pilot('run=#');
%
%    If you close the screen inbetween runs you need to repeat the
%    'genTex=1' call, otherwise the textures will be incorrectly
%    identified.

global stimulus fixStimulus

if ~isstruct(stimulus), stimulus = struct; end
if ~isstruct(fixStimulus), fixStimulus = struct; end
% clear stimulus

% reset
f = fields(stimulus);
for fi = 1:length(f)
    if ~strcmp(f{fi},'fbsdata')
        stimulus = rmfield(stimulus,f{fi});
    end
end
% clear fix stimulus
% reset
f = fields(fixStimulus);
for fi = 1:length(f)
    fixStimulus = rmfield(fixStimulus,f{fi});
end

%% Initialize Variables

% add arguments later
plots = 0;
run = 0;
genTex = 0;
getArgs(varargin,{'plots=0','run=0','genTex=0'});
stimulus.plots = plots;
stimulus.run = run;
stimulus.generateTextures = genTex;

clear plots run genTex

if stimulus.run==0
    disp('Please set a run number');
    return;
end

if ~stimulus.generateTextures && ~isfield(stimulus,'fbsdata')
    disp('WARNING: You must generate textures');
    return
end

if ~stimulus.generateTextures
    disp('WARNING: If the screen was closed textures will not work');
end

%% Image categories
categories = {'person','car','personcar','null'};
attend_categories = {'attend_person','attend_car'};

%% Open Old Stimfile
stimulus.counter = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/fbsear_pilot/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/fbsear_pilot/%s/1*mat',mglGetSID));

    if length(files) >= 1
        fname = files(end).name;
        
        s = load(sprintf('~/data/fbsear_pilot/%s/%s',mglGetSID,fname));
        % copy staircases and run numbers
        stimulus.counter = s.stimulus.counter + 1;
        stimulus.runs = s.stimulus.runs;
        clear s;
        disp(sprintf('(fbsear_pilot) Data file: %s loaded.',fname));
    end
end
disp(sprintf('(fbsear_pilot) This is run #%i',stimulus.counter));

%% Setup Screen
myscreen = initScreen('VPixx');

if ~isfield(myscreen,'genTexTrack')
    myscreen.genTexTrack = rand*10000000;
end

% set background to grey
myscreen.background = 0.5;

%% Load images

if stimulus.generateTextures
    loadFBSearImages(categories,attend_categories);
end

%% Setup missing initial variables

if ~isfield(stimulus,'counter')
    stimulus.counter = 1; % This keeps track of what "run" we are on.
end

%% Initialize Stimulus

myscreen = initStimulus('stimulus',myscreen);
  
stimulus.responseKeys = [1 2]; % 

%% Colors

stimulus.colors.white = [1 1 1];
stimulus.colors.black = [0 0 0];
stimulus.colors.grey = [.7 .7 .7];
stimulus.colors.red = [1 0 0];
stimulus.colors.green = [0 1 0];
% initGammaTable(myscreen);
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

%% Setup runs

if ~isfield(stimulus,'runs')
    disp('WARNING WARNING WARNING WARNING WARNING WARNING');
    disp('WARNING WARNING WARNING WARNING WARNING WARNING');
    disp('WARNING WARNING WARNING WARNING WARNING WARNING');
    disp('WARNING WARNING WARNING WARNING WARNING WARNING');
    disp('WARNING WARNING WARNING WARNING WARNING WARNING');
    disp('WARNING: Building new runs!!!!!!!!!!!!! WARNING');
    disp('If you are in the middle of a scan session this is bad');
    disp('it means that the stim files were not saved/loaded correctly');
    disp('WARNING WARNING WARNING WARNING WARNING WARNING');
    disp('WARNING WARNING WARNING WARNING WARNING WARNING');
    disp('WARNING WARNING WARNING WARNING WARNING WARNING');
    disp('WARNING WARNING WARNING WARNING WARNING WARNING');
    disp('WARNING WARNING WARNING WARNING WARNING WARNING');
    stimulus.runs = struct;
    
    % build runs
    
    % Choose task types (fix/p/c/fix-cp/fix-cc)
    fixOpts = [1 0 0 1 1]; % whether to run fixation task or not
    maskOpts = [0 0 0 1 2]; % which mask set to use, boosted for people or boosted for cars
    attend = [0 1 2 0 0]; % 1 = attend person, 2 = attend car
    text = {'Fixate','Look for people','Look for cars','Fixate (contrast people)','Fixate (contrast cars)'};
    
    % Pick 45 images for each of the three runs
    idx45 = zeros(length(categories),45);
    for ci = 1:length(categories)
        idx_opts = randperm(length(stimulus.fbsdata.imgs.(categories{ci})));
        idx45(ci,:) = idx_opts(1:45);
    end
    % we now have the indexes for images from each set, copy out 15 for
    % each set
    img_idxs = zeros(length(categories),3,15);
    for ci = 1:length(categories)
        for rep = 1:3
            img_idxs(ci,rep,:) = idx45(ci,(rep-1)*15+1:(rep*15));
        end
    end
    
    stimulus.runs.imageIndexes = img_idxs;
    
    % Build run data
    for rep = 1:3
        % set the trial order info 
        trials = repmat(1:15,1,4);
        cats = [1*ones(1,15) 2*ones(1,15) 3*ones(1,15) 4*ones(1,15)];
        trialOrder = randperm(length(trials));
        
        repeatTrials = trials(trialOrder);
        repeatCats = cats(trialOrder);
        
        for run = 1:5
            runData = struct;

            runData.fixate = fixOpts(run);
            runData.useMask = maskOpts(run);
            runData.attend = attend(run);
            runData.text = text{run};
            
            runData.group = run;
            runData.repeat = rep;
            
            runData.trialOrder = repeatTrials;
            runData.catOrder = repeatCats;
            
            % Save
            stimulus.runs.runs{run,rep} = runData;
        end
    end
    
    % Copy so that we have at least 3x each run
    stimulus.runs.runOrder = [];
    stimulus.runs.runData = {};
    runNums = 1:5;
    for reps = 1:3
        % for each rep round, 
        runOrder = runNums(randperm(length(runNums)));
        for run = 1:5
            stimulus.runs.runData((reps-1)*5+run) = stimulus.runs.runs(runOrder(run),reps);
        end
        stimulus.runs.runOrder = [stimulus.runs.runOrder runOrder];
    end
end

% Choose run
stimulus.curRun = stimulus.runs.runData{stimulus.run};
disp('                                                                  ');
disp('`````````````````````````````````````````````````````````````````');
disp('`````````````````````````````````````````````````````````````````');
disp('```````````````RUN INFO: WRITE THIS DOWN`````````````````````````');
disp(sprintf('```` Current run type: %s',stimulus.curRun.text));
disp(sprintf('```` Group type: %i',stimulus.curRun.group));
disp(sprintf('```` This is repeat: %i',stimulus.curRun.repeat));
maskType = {'NO','YES'};
disp(sprintf('```` Using mask data set: %s',maskType{stimulus.curRun.useMask+1}));
disp('`````````````````````````````````````````````````````````````````');
disp('`````````````````````````````````````````````````````````````````');
%% Pull textures

% pull the textures for this run
for ci = 1:length(categories)
    if stimulus.curRun.useMask==0
        texs = stimulus.fbsdata.imgs.(categories{ci});
        texs = texs(squeeze(stimulus.runs.imageIndexes(ci,stimulus.curRun.repeat,:)));
        stimulus.curRun.stimulus{ci} = texs;
    else
        texs = stimulus.fbsdata.mimgs.(attend_categories{stimulus.curRun.useMask}).(categories{ci});
        texs = texs(squeeze(stimulus.runs.imageIndexes(ci,stimulus.curRun.repeat,:)));
        stimulus.curRun.stimulus{ci} = texs;
    end
end

%% Setup Task

%%%%%%%%%%%%% PHASE ONE %%%%%%%%%%%%%%%%%
%%%%% PRIOR + ESTIMATE OF THRESHOLD %%%%%

task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;
task{1}{1}.seglen = [1.9 3.9];
stimulus.seg.stim = 1;
stimulus.seg.iti = 2;

task{1}{1}.synchToVol = [1 1];

task{1}{1}.getResponse = [0 1];
task{1}{1}.numTrials = 60; % 1096 volumes for scan

task{1}{1}.random = 1;

% stimulus controls (which image gets shown)
% task{1}{1}.parameter.trial = 1:15;
% task{1}{1}.parameter.category = 1:4;

% Task variables to be calculated later
task{1}{1}.randVars.calculated.trial = nan;
task{1}{1}.randVars.calculated.category = nan;
task{1}{1}.randVars.calculated.task = nan; %0/1/2
task{1}{1}.randVars.calculated.group = nan;
task{1}{1}.randVars.calculated.repeat = nan;
task{1}{1}.randVars.calculated.taskPresent = nan;

%% Add fixation task

if stimulus.curRun.fixate
    
    if ~isfield(fixStimulus,'threshold') fixStimulus.threshold = 0.3; end
    if ~isfield(fixStimulus,'stairStepSize') fixStimulus.stairStepSize = 0.05; end

    [task{2}, myscreen] = fixStairInitTask(myscreen);
end

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
disp(sprintf('(fbsear_pilot) Starting run number: %i.',stimulus.counter));

%% Live

stimulus.live.correctButton = 2;
stimulus.live.lastCorrect = -1;
stimulus.live.changeTime = -1;

%% Main Task Loop

for i = 1:2
    mglClearScreen(0.5); 
    mglTextSet([],32,stimulus.colors.white);
    if ~stimulus.curRun.fixate
        mglFixationCross(1,1,stimulus.colors.white);
        if stimulus.curRun.attend==1
            mglTextDraw('Look for people',[0 2]);
        elseif stimulus.curRun.attend==2
            mglTextDraw('Look for cars',[0 2]);
        end
    else
        mglTextDraw('Fixate',[0 2]);
    end
    mglFlush
end
% myscreen.flushMode = 1;

phaseNum = 1;
% Again, only one phase.
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    % update the task
    [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);
    if stimulus.curRun.fixate
        [task{2}, myscreen] = updateTask(task{2},myscreen,1);
    end
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

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startTrialCallback(task,myscreen)
%%
global stimulus
stimulus.live.fixColor = [0 0 0];

task.thistrial.trial = stimulus.curRun.trialOrder(task.trialnum);
task.thistrial.category = stimulus.curRun.catOrder(task.trialnum);

correctButtons = {[1 2 1 2] [2 1 1 2]};
if stimulus.curRun.attend>0
    flip = [1 0];
    stimulus.live.correctButton= correctButtons{stimulus.curRun.attend}(task.thistrial.category);
    task.thistrial.taskpresent = flip(stimulus.live.correctButton);
end
task.thistrial.task = stimulus.curRun.attend; %0/1/2 - fixate/people/cars
task.thistrial.group = stimulus.curRun.group; % this is the ATTENTION group (1/2/3/4/5)
task.thistrial.repeat = stimulus.curRun.repeat;

disp(sprintf('Trial #%i: Correct button: %i',task.trialnum,stimulus.live.correctButton));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)
% pass
global stimulus
if task.thistrial.thisseg == stimulus.seg.stim
    stimulus.live.fixColor = stimulus.colors.black;
else
    stimulus.live.fixColor = stimulus.colors.grey;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus
mglClearScreen(0.5);
if task.thistrial.thisseg == stimulus.seg.stim
    % blt the current texture
    mglBltTexture(stimulus.curRun.stimulus{task.thistrial.category}{task.thistrial.trial},[0 0]);
end

if ~stimulus.curRun.fixate
    if (stimulus.live.changeTime>0) && ( mglGetSecs > (stimulus.live.changeTime+1) )
        stimulus.live.fixColor = stimulus.colors.black;
        stimulus.live.changeTime = -1;
    end
    upFix(stimulus);
end

function upFix(stimulus)
%%
% for this experiment use a circle to indicate where participants can
% fixate inside of (rather than a cross which might arbitrarily enforce
% poisitioning
mglGluDisk(0,0,[1 1],0.5,60);
% mglGluAnnulus(0,0,1.5,1.55,stimulus.live.fixColor,64);
mglFixationCross(1,3,stimulus.live.fixColor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus

if ~any(task.thistrial.whichButton==stimulus.responseKeys)
    return
end

if ~stimulus.curRun.fixate

    task.thistrial.gotResponse = task.thistrial.gotResponse + 1;

    if task.thistrial.gotResponse==1
        if task.thistrial.whichButton==stimulus.live.correctButton
            task.thistrial.correct = 1;
            stimulus.live.fixColor = stimulus.colors.green;
            disp('Correct');
        else
            task.thistrial.correct = 0;
            stimulus.live.fixColor = stimulus.colors.red;
            disp('Incorrect');
        end
        stimulus.live.changeTime = mglGetSecs;
    else
        disp(sprintf('Got %i responses',task.thistrial.gotResponse));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              HELPER FUNCTIONS                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function loadFBSearImages(cats,acats)
%% load and setup
global stimulus

stimulus.fbsdata = struct;

% get image data
load(fullfile('~/proj/fbsear/data/info_m.mat'));

%% prepare mgl textures (both contrast masked and non-masked)

stimulus.fbsdata.imgs = struct;
for ci = 1:length(cats)
    
    cat = cats{ci};
    
    img_texs = cell(size(info.imgs.(cat)));
    
    imgs = info.imgs.(cat);
    
    for ii = 1:length(imgs)
        img = double(imgs{ii});
        for rgb = 1:3
            img(:,:,rgb) = flipud(img(:,:,rgb));
        end
%         img = permute(img,[3 1 2]);
% %         img = reshape(img,size(img,3),size(img,1),size(img,2));
%         img(4,:,:) = ones(1,size(img,2),size(img,3));
        img_texs{ii} = mglCreateTexture(img);
    end
    
    stimulus.fbsdata.imgs.(cat) = img_texs;
end

stimulus.fbsdata.mimgs = struct;
for aci = 1:length(acats)
    acat = acats{aci};
    for ci = 1:length(cat)
        cat = cats{ci};
        mimg_texs = cell(size(info.mimgs.(acat).(cat)));
        mimgs = info.mimgs.(acat).(cat);
        for ii = 1:length(imgs)
            img = double(mimgs{ii});
            for rgb = 1:3
                img(:,:,rgb) = flipud(img(:,:,rgb));
            end
    %         img = permute(img,[3 1 2]);
    % %         img = reshape(img,size(img,3),size(img,1),size(img,2));
    %         img(4,:,:) = ones(1,size(img,2),size(img,3));
            mimg_texs{ii} = mglCreateTexture(img);
        end
        stimulus.fbsdata.mimgs.(acat).(cat) = mimg_texs;
    end
end