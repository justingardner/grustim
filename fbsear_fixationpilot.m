function [ myscreen ] = fbsear_fixationpilot( varargin )
%
% SEMANTIC MODELING TASK
%   Show a bunch of COCO images to try to estimate a semantic model of the
%   COCO dataset with the maximum number of categories. 
%
%    Usage: fbsear_fixationpilot(varargin)
%    Authors: Dan Birman
%    Date: 04/13/2017
%
%    Parameters:
%
%    You should use *836* volumes (270 images * 1.5 secs * 2 + 16) you will
%    have to manually end each run by hitting ESC (if you use 836+1 a volume
%    is acquired after the screen goes dark, so don't do that)
%
%      On the first run call:
%       fbsear_fixationpilot('run=1','genTex=1');
%
%      All subsequent runs:
%       fbsear_fixationpilot('run=#');
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
screen = '';
getArgs(varargin,{'plots=0','run=0','genTex=0','screen=fMRIProj32'});
stimulus.plots = plots;
stimulus.run = run;
stimulus.generateTextures = genTex;
stimulus.screen = screen;

if strcmp(screen,'fMRIProj32')    
    stimulus.scale = 1.35;
else
    stimulus.scale = 1;
end

clear plots run genTex screen

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

%% Open Old Stimfile
stimulus.counter = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/fbsear_fixationpilot/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/fbsear_fixationpilot/%s/1*mat',mglGetSID));

    if length(files) >= 1
        fname = files(end).name;
        
        s = load(sprintf('~/data/fbsear_fixationpilot/%s/%s',mglGetSID,fname));
        % copy staircases and run numbers
        stimulus.counter = s.stimulus.counter + 1;
        stimulus.runs = s.stimulus.runs;
        clear s;
        disp(sprintf('(fbsear_fixationpilot) Data file: %s loaded.',fname));
    end
end
disp(sprintf('(fbsear_fixationpilot) This is run #%i',stimulus.counter));

%% Setup Screen
myscreen = initScreen(stimulus.screen);

if ~isfield(myscreen,'genTexTrack')
    myscreen.genTexTrack = rand*10000000;
end

% set background to grey
myscreen.background = 0.5;

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

%% Load images

if stimulus.generateTextures
    loadSemImages();
end
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
    runNums = 1:4;
    
    % randomize the order of images
    order = randperm(stimulus.fbsdata.n);
    % copy into four runs
    idxs = 0:(1080/4):length(order);
    
    for ri = 1:4
        crun = struct;
        crun.idxs = order((idxs(ri)+1):idxs(ri+1));
        crun.text = sprintf('Run group %i',ri);
        crun.fixate = true;
        for reps = 1:3
            crun.repeat = reps;
            stimulus.runs.runs{ri,reps} = crun;
        end
    end
    
    % Copy so that we have at least 3x each run
    stimulus.runs.runOrder = [];
    stimulus.runs.runData = {};
    for reps = 1:3
        % for each rep round, 
        runOrder = runNums(randperm(length(runNums)));
        for run = 1:length(runNums)
            stimulus.runs.runData((reps-1)*length(runNums)+run) = stimulus.runs.runs(runOrder(run));
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
disp(sprintf('```` This is repeat: %i',stimulus.curRun.repeat));
disp('`````````````````````````````````````````````````````````````````');
disp('`````````````````````````````````````````````````````````````````');


%% Setup Task

%%%%%%%%%%%%% PHASE ONE %%%%%%%%%%%%%%%%%
%%%%% PRIOR + ESTIMATE OF THRESHOLD %%%%%

task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;
task{1}{1}.seglen = [1.4];
stimulus.seg.stim = 1;

task{1}{1}.synchToVol = [1 1];

% task{1}{1}.getResponse = [0 1];
task{1}{1}.numTrials = 1080/4; % 1096 volumes for scan

task{1}{1}.random = 1;

% stimulus controls (which image gets shown)
% task{1}{1}.parameter.trial = 1:15;
% task{1}{1}.parameter.category = 1:4;

% Task variables to be calculated later
task{1}{1}.randVars.calculated.trial = nan;
task{1}{1}.randVars.calculated.img = nan;
task{1}{1}.randVars.calculated.repeat = nan;

%% Add fixation task

if stimulus.curRun.fixate
    
    if ~isfield(fixStimulus,'threshold') fixStimulus.threshold = 0.3; end
    if ~isfield(fixStimulus,'stairStepSize') fixStimulus.stairStepSize = 0.05; end
    if ~isfield(fixStimulus,'diskSize') fixStimulus.diskSize = 0.5; end
    if ~isfield(fixStimulus,'fixWidth') fixStimulus.fixWidth = 0.5; end
    if ~isfield(fixStimulus,'fixLineWidth') fixStimulus.fixLineWidth = 2; end

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
disp(sprintf('(fbsear_fixationpilot) Starting run number: %i.',stimulus.counter));

%% Live
% 
% stimulus.live.correctButton = 2;
% stimulus.live.lastCorrect = -1;
% stimulus.live.changeTime = -1;

%% Main Task Loop

for i = 1:2
    mglClearScreen(0.5); 
    mglTextSet([],32,stimulus.colors.white);
    if ~stimulus.curRun.fixate
%         mglFixationCross(1,1,stimulus.colors.white);
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

disp(sprintf('Stimulus #%i: image %i.',task.trialnum,stimulus.curRun.idxs(task.trialnum)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)
% pass

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus
mglClearScreen(0.5);
    % blt the current texture
if stimulus.scale~=1
    pixPerDegW = myscreen.screenWidth/myscreen.imageWidth;
    pixPerDegH = myscreen.screenHeight/myscreen.imageHeight;
    cTex = stimulus.fbsdata.tex{stimulus.curRun.idxs(task.trialnum)};
    degWidth = stimulus.scale*cTex.imageWidth/pixPerDegW;
    degHeight = stimulus.scale*cTex.imageHeight/pixPerDegH;
    mglBltTexture(cTex,[0 0 degWidth degHeight]);
else
    mglBltTexture(stimulus.fbsdata.tex{stimulus.curRun.idxs(task.trialnum)},[0 0]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [task, myscreen] = getResponseCallback(task, myscreen)
% pass

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              HELPER FUNCTIONS                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function loadSemImages()
%% load and setup
global stimulus

stimulus.fbsdata = struct;

% get image data
load(fullfile('~/proj/fbsear/data/semantic.mat'));

%% prepare mgl textures (both contrast masked and non-masked)

stimulus.fbsdata.n = length(semdata.imgs);
stimulus.fbsdata.tex = cell(1,stimulus.fbsdata.n);

for ii = 1:length(semdata.imgs)
    img = double(semdata.imgs{ii});
    for rgb = 1:3
        img(:,:,rgb) = flipud(img(:,:,rgb));
    end
    
    stimulus.fbsdata.tex{ii} = mglCreateTexture(img);
end