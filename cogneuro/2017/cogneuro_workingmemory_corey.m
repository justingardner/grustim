
% cohCon
%
%      usage: myscreen=cohcon()
%         by: daniel birman
%       date: 11/10/14
%    purpose: contrast change detection with cued selective attention.
%
%        use: call flowAwe() to initialize. The first
%             run takes significantly longer due to loading stimuli.
%
%      flags: stimFileNum (-1/#) - Load a specific stimfile from a
%             subject's folder. Defaults to the last stimfile. Warning:
%             Only the first stimfile is saved with the image file data,
%             subsequent stimfiles only contain run data.
%             unattended (0/1) - If 1, runs a fixation task while the main
%             task just runs on idle in the background (no inputs)
%             plots (0/1) - Displays staircase plots (and estimated
%             psychophysic functions)
%             overrideTask (1/2) - Specifies the task to run: 1 =
%             coherence, 2 = contrast
%             projector (0/1) - Masks stimuli using the default projector
%             mask.
%             scan (0/1) - Scanner timing
%
%
%   TR .5 = 296 volumes (10 * 14 * 2 + 16)

function [myscreen] = cogneuro_workingmemory_corey(varargin)

global stimulus
clear fixStimulus
%% Initialize Variables

disp('*********************************************');
disp('** NEPR Into Cog Neuro Working Memory Task **');
disp('*********************************************');
% add arguments later
scan = [];
testing = [];
getArgs(varargin,{'testing=0','scan=1'});
stimulus.scan = scan;
stimulus.testing = testing;

if stimulus.scan && ~mglGetParam('ignoreInitialVols')==16 && ~mglGetParam('ignoreInitialVols')==4
    warning('ignoreInitialVols is set to %i.',mglGetParam('ignoreInitialVols'));
    if ~strcmp('y',input('Is this correct? [y/n]'))
        mglSetParam('ignoreInitialVols',input('Please input the correct value (mux8 = 16, mux2 = 4): '));
    end
end

stimulus.counter = 1; % This keeps track of what "run" we are on.


%% Colors
stimulus.colors.rmed = 127.5;

% We're going to add an equal number of reserved colors to the top and
% bottom, to try to keep the center of the gamma table stable.
stimulus.colors.reservedBottom = [0.8 0 0; 0 0 0]; % fixation cross colors
stimulus.colors.reservedTop = [0.8 0.8 0.8; 0 0.8 0]; % correct/incorrect colors
stimulus.colors.black = 1/255; stimulus.colors.white = 254/255;
stimulus.colors.red = 0/255; stimulus.colors.green = 255/255;
stimulus.colors.nReserved = 2; % this is /2 the true number, because it's duplicated
stimulus.colors.nUnreserved = 256-(2*stimulus.colors.nReserved);

stimulus.colors.mrmax = stimulus.colors.nReserved - 1 + stimulus.colors.nUnreserved;
stimulus.colors.mrmin = stimulus.colors.nReserved;


%% Stimulus (oriented gratings) properties
stimulus.ecc = 6;
stimulus.angleOpts = [30 120];
stimulus.shiftPos = 18;
stimulus.shiftRot = 5;

% set stimulus contrast, scanner absolute luminance is lower so the
% relative contrast needs to be higher to make up for this
if stimulus.scan
    stimulus.contrast = .1; % 4%
else
    stimulus.contrast = .5;
end

%% Setup Screen
% Settings for Oban/Psychophysics room launch
if stimulus.scan
    myscreen = initScreen('fMRIprojFlex');
else
    myscreen = initScreen('VPixx'); % this will default out if your computer isn't the psychophysics computer
end
myscreen.background = 0.5;

%% Open Old Stimfile
stimulus.counter = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/cogneuro_workingmemory_corey/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/cogneuro_workingmemory_corey/%s/1*mat',mglGetSID));

    if length(files) >= 1
        fname = files(end).name;
        
        s = load(sprintf('~/data/cogneuro_workingmemory_corey/%s/%s',mglGetSID,fname));
        % copy staircases and run numbers
        stimulus.counter = s.stimulus.counter + 1;
        clear s;
        disp(sprintf('(cn_corey) Data file: %s loaded.',fname));
    end
end
disp(sprintf('(cn_corey) This is run #%i',stimulus.counter));

%% Initialize Stimulus
% Setup for MGL
myscreen = initStimulus('stimulus',myscreen);
stimulus.responseKeys = [1 2]; % corresponds to LEFT - RIGHT

%% Setup Task

% Blocks the task so it doesn't start immediately
task{1}{1}.waitForBacktick = 1;

% Storing some helper info for callbacks later
stimulus.seg.wait = 1;
stimulus.seg.stim = 2;
stimulus.seg.mask = 3;
stimulus.seg.delay = 4;
stimulus.seg.resp = 5;
stimulus.seg.ITI = 6;

% Trial timing (see above for what each column corresponds to)
task{1}{1}.segmin = [0.500 1.000 0.250 12 2.500 2.500];
task{1}{1}.segmax = [0.500 1.000 0.250 12 2.500 7.500]; % 5 s average 

% When scanning we synchronize the stimulus to the scanner
task{1}{1}.synchToVol = zeros(1,length(task{1}{1}.segmin));
if stimulus.scan
    task{1}{1}.synchToVol(end) = 1;
end

% Get responses on the resp segment (4)
task{1}{1}.getResponse = zeros(length(task{1}{1}.segmin));
task{1}{1}.getResponse(stimulus.seg.resp) = 1;

% Parameters that control the direction of each oriented grating
task{1}{1}.parameter.axis = [0 45]; % 0 = on axis, 45 = off axis (45 degs off)
task{1}{1}.parameter.task = [1 2]; % 1 = circle judgment task, 2 = orientation judgment task
task{1}{1}.parameter.rotation = [1 2]; % counter/clock rotation
task{1}{1}.parameter.changeLoc = [1 2 3 4];
% This will randomize trials
task{1}{1}.random = 1;
% Outside the scanner fix the trial count
task{1}{1}.numTrials = 50;
% Inside the scanner, inf length so that we can stop the stimulus after the
% scanner stops running.
if stimulus.scan
    task{1}{1}.numTrials = inf;
end

%% Tracking
% these are variables that we want to track for later analysis.
task{1}{1}.randVars.calculated.correct = nan;
task{1}{1}.randVars.calculated.jitter = nan;

stimulus.curTrial = 0;

%% Generate stencil
mglStencilCreateBegin(1);
mglGluAnnulus(0,0,stimulus.ecc-1,stimulus.ecc+1,repmat(stimulus.colors.white,3,1),32);
mglStencilCreateEnd;

%% Gratings
sz = 4;
g = mglMakeGrating(sz*2,sz*2,1,0,0);
gauss = mglMakeGaussian(sz*2,sz*2,sz/6,sz/6);
% normalize
g = (g .* gauss + 1) / 2; % bounded 0-1
g =  (stimulus.colors.nUnreserved-stimulus.colors.nReserved)* g +stimulus.colors.nReserved + 1;
% save texture for faster rendering
stimulus.live.grating = mglCreateTexture(g);

%% Full Setup
% Initialize task (note phase == 1) (this is for MGL)
for phaseNum = 1:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);
end

%% Get Ready...
% clear screen
mglTextSet('Helvetica',32,stimulus.colors.white);
stimulus.live.Otext = mglText('O');
stimulus.live.Ptext = mglText('P');

mglClearScreen(0.5);
mglFixationCross(1,1,stimulus.colors.black);
mglFlush
mglClearScreen(0.5);
mglFixationCross(1,1,stimulus.colors.black);
mglFlush
setGammaTable_flowMax(stimulus.contrast);

% let the user know
disp(sprintf('(cn_corey) Starting run number: %i',stimulus.counter));

%% Main Task Loop

% Run the task, this is for MGL
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
mglTextDraw('Run complete... please wait.',[0 0]);
mglFlush
myscreen.flushMode = 1;

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS AND CALLBACKS FOLLOW %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startTrialCallback(task,myscreen)
%%

global stimulus

task.thistrial.jitter = randn*10;

% Increment the trial number
stimulus.curTrial = stimulus.curTrial + 1;
task.thistrial.trialNum = stimulus.curTrial;

task.thistrial.length = task.thistrial.seglen(4);
 
tasks = {'position','orientation'};
sides = {'left','right'};
disp(sprintf('Trial %i: task %s, rotation %s',stimulus.curTrial,tasks{task.thistrial.task},sides{task.thistrial.rotation}));

stimulus.live.positions = (0:90:270) + 45*task.thistrial.axis + task.thistrial.jitter;
stimulus.live.rotations = 90*ones(1,4); % all vertical gratings

% WN mask
wn = repmat(stimulus.colors.mrmin+randi(251,1,myscreen.screenWidth/4,myscreen.screenHeight/4,'uint8')-1,3,1,1);
wn(4,:,:) = 255;
stimulus.live.wn = mglCreateTexture(wn,[],0,{'GL_TEXTURE_MAG_FILTER','GL_NEAREST'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    runs at the start of each segment       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)
%%
myscreen.flushMode = 0;
global stimulus

stimulus.live.fix = 1;
stimulus.live.cue = 0;
stimulus.live.grate = 0;
stimulus.live.mask = 0;
stimulus.live.fixColor = stimulus.colors.white;

% All we're going to do is adjust what gets shown based on what segment we
% are in, the actual updating happens in the next function.
switch task.thistrial.thisseg
    case stimulus.seg.stim
        stimulus.live.grate = 1;
    case stimulus.seg.mask
        stimulus.live.mask = 1;
        stimulus.live.fix = 0;
    case stimulus.seg.resp
        stimulus.live.grate = 1;
        flip = [1 -1];
        if task.thistrial.task==1
            stimulus.live.positions(task.thistrial.changeLoc) = stimulus.live.positions(task.thistrial.changeLoc) + flip(task.thistrial.rotation)*stimulus.shiftPos;
        else
            stimulus.live.rotations = stimulus.live.rotations + flip(task.thistrial.rotation)*stimulus.shiftRot;
        end
    case stimulus.seg.ITI
        stimulus.live.fixColor = stimulus.colors.black;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus
mglClearScreen(0.5);

if stimulus.live.mask
    mglStencilSelect(1);
    mglBltTexture(stimulus.live.wn,[0 0 myscreen.imageWidth myscreen.imageWidth]);
    mglStencilSelect(0);
    if task.thistrial.task==1
        mglBltTexture(stimulus.live.Ptext,[0 0]);
    else
        mglBltTexture(stimulus.live.Otext,[0 0]);
    end
end

% actuall display the grating
if stimulus.live.grate, upGrate(task); end
% display the fixation cross (always)
if stimulus.live.fix, upFix(stimulus); end
% display the cue (this draws over the cross)
if stimulus.live.cue, upCue(task); end

function upCue(task)
%%
global stimulus
% some helper code for training
%atex = {'Attend Left', 'Attend Right'};
%mglTextDraw(atex{task.thistrial.attend},[0 1.5]);

% otherwise just draw the line over the fixation cross
if task.thistrial.attend==1
    % left
    mglLines2(0,0,-.5,0,1,stimulus.colors.black);
elseif task.thistrial.attend==2
    mglLines2(0,0,.5,0,1,stimulus.colors.black);
end

function upFix(stimulus)
%%
mglFixationCross(1,1,stimulus.live.fixColor);

function upGrate(task)
global stimulus
% distance is to the center of the grating, in degrees

if task.thistrial.thisseg==stimulus.seg.stim
    x = stimulus.ecc * cos(deg2rad(stimulus.live.positions));
    y = stimulus.ecc * sin(deg2rad(stimulus.live.positions));
    mglGluDisk(x,y,1,stimulus.colors.white,16);
elseif task.thistrial.thisseg==stimulus.seg.resp
    x = stimulus.ecc * cos(deg2rad(stimulus.live.positions));
    y = stimulus.ecc * sin(deg2rad(stimulus.live.positions));
    for i = 1:4
        mglBltTexture(stimulus.live.grating,[x(i) y(i)],0,0,stimulus.live.rotations(i));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus

% setup
responseText = {'Incorrect','Correct'};
responsePos = {'Left','Right'};
fixColors = {stimulus.colors.red, stimulus.colors.green};

% check that we got a response and it's a key we want
if any(task.thistrial.whichButton == stimulus.responseKeys)
    % check that we haven't already responded
    if task.thistrial.gotResponse < 2
        
        % check if subject was correct
        task.thistrial.correct = task.thistrial.whichButton == task.thistrial.rotation;

        stimulus.live.fixColor = fixColors{task.thistrial.correct+1};
        disp(sprintf('(cn_corey) Response %s: %s',responseText{task.thistrial.correct+1},responsePos{stimulus.responseKeys(task.thistrial.whichButton)}));
    else
        disp(sprintf('(cn_corey) Subject responded multiple times: %i',task.thistrial.gotResponse+1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              HELPER FUNCTIONS                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sets the gamma table so that we can have
% finest possible control over the stimulus contrast.
%
% stimulus.colors.reservedColors should be set to the reserved colors 
% (for cue colors, etc). maxContrast is the maximum contrast you want to 
% be able to display.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setGammaTable_flowMax(maxContrast)

% definitely don't mess with this function.

global stimulus;

if ~isfield(stimulus,'linearizedGammaTable')
    stimulus.linearizedGammaTable = mglGetGammaTable;
end

% set the bottom
gammaTable(1:size(stimulus.colors.reservedBottom,1),1:size(stimulus.colors.reservedBottom,2)) = stimulus.colors.reservedBottom;

% set the gamma table
if maxContrast == 1
    % create the rest of the gamma table
    cmax = 1;cmin = 0;
    luminanceVals = cmin:((cmax-cmin)/(stimulus.colors.nUnreserved-1)):cmax;

    % now get the linearized range
    redLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,luminanceVals,'linear');
    greenLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,luminanceVals,'linear');
    blueLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,luminanceVals,'linear');
elseif maxContrast > 0
    % create the rest of the gamma table
    cmax = 0.5+maxContrast/2;cmin = 0.5-maxContrast/2;
    luminanceVals = cmin:((cmax-cmin)/(stimulus.colors.nUnreserved-1)):cmax;

    % now get the linearized range
    redLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,luminanceVals,'linear');
    greenLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,luminanceVals,'linear');
    blueLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,luminanceVals,'linear');
else
    % if we are asked for 0 contrast then simply set all the values to gray
    redLinearized = repmat(.5,1,stimulus.colors.nUnreserved);
    greenLinearized = repmat(.5,1,stimulus.colors.nUnreserved);
    blueLinearized = repmat(.5,1,stimulus.colors.nUnreserved);
end

% add to the table!
gammaTable((stimulus.colors.mrmin:stimulus.colors.mrmax)+1,:)=[redLinearized;greenLinearized;blueLinearized]';

% set the top
gammaTable = [gammaTable; stimulus.colors.reservedTop];

if size(gammaTable,1)~=256
    disp('(setGammaTable) Failure: Incorrect number of colors in gamma table produced');
    keyboard
end

% set the gamma table
mglSetGammaTable(gammaTable);

% remember what the current maximum contrast is that we can display
stimulus.curMaxContrast = maxContrast;