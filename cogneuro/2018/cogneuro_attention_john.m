
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

function [myscreen] = cogneuro_attention_john(varargin)

global stimulus
clear fixStimulus
%% Initialize Variables

disp('****************************************');
disp('** NEPR Into Cog Neuro Attention Task **');
disp('****************************************');
% add arguments later
stimFileNum = [];
scan = [];
testing = [];
getArgs(varargin,{'stimFileNum=-1','testing=0','scan=0'});
stimulus.scan = scan;
stimulus.testing = testing;

if stimulus.scan && ~mglGetParam('ignoreInitialVols')==16 && ~mglGetParam('ignoreInitialVols')==6
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
stimulus.colors.reservedBottom = [0 0 0; .8 .8 .8]; % fixation cross colors
stimulus.colors.reservedTop = [.55 0 0; 0 .55 0]; % correct/incorrect colors
stimulus.colors.black = 0/255; stimulus.colors.white = 1/255;
stimulus.colors.red = 254/255; stimulus.colors.green = 255/255;
stimulus.colors.nReserved = 2; % this is /2 the true number, because it's duplicated
stimulus.colors.nUnreserved = 256-(2*stimulus.colors.nReserved);

stimulus.colors.mrmax = stimulus.colors.nReserved - 1 + stimulus.colors.nUnreserved;
stimulus.colors.mrmin = stimulus.colors.nReserved;


%% Stimulus (oriented gratings) properties
stimulus.pedestals.angle = 90;
stimulus.pedestals.initThresh.angle = 10.0; % angle that the staircase will start at

% set stimulus contrast, scanner absolute luminance is lower so the
% relative contrast needs to be higher to make up for this
if stimulus.scan
    stimulus.contrast = .5; % 4%
else
    stimulus.contrast = .5;
end

%% Setup Screen
% Settings for Oban/Psychophysics room launch
if stimulus.scan
    myscreen = initScreen('fMRIproj32');
else
    myscreen = initScreen('VPixx'); % this will default out if your computer isn't the psychophysics computer
end

%% Open Old Stimfile
% Do not modify this code
stimulus.initStair = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/cogneuro_attention/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/cogneuro_attention/%s/1*mat',mglGetSID));
    
    if length(files) >= 1
        if stimFileNum == -1
            if length(files) > 1
                warning('Multiple stimfiles found, loading last one. If you wanted a different functionality use stimFileNum=#');
            end
            fname = files(end).name;
        else
            fname = files(stimFileNum).name;
        end
        s = load(sprintf('~/data/cogneuro_attention/%s/%s',mglGetSID,fname));
        stimulus.staircase = s.stimulus.staircase;
        stimulus.counter = s.stimulus.counter + 1;
        
        clear s;
        stimulus.initStair = 0;
        disp(sprintf('(cogneuro_att) Data file: %s loaded.',fname));
    end
end
disp(sprintf('(cogneuro_att) This is run #%i',stimulus.counter));

%% Initialize Stimulus
% Setup for MGL
myscreen = initStimulus('stimulus',myscreen);
stimulus.responseKeys = [1 2]; % corresponds to LEFT - RIGHT

%% Setup Task

% Blocks the task so it doesn't start immediately
task{1}{1}.waitForBacktick = 1;

% Storing some helper info for callbacks later
stimulus.seg.cue = 1;
stimulus.seg.isi = 2;
stimulus.seg.stim = 3;
stimulus.seg.resp = 4;
stimulus.seg.ITI = 5;

% Trial timing (see above for what each column corresponds to)
task{1}{1}.segmin = [0.5 1 1 2 1];
task{1}{1}.segmax = [0.5 1 1 2 8];

% When scanning we synchronize the stimulus to the scanner
task{1}{1}.synchToVol = [0 0 0 0 0];
if stimulus.scan
    task{1}{1}.synchToVol(end) = 1;
end

% Get responses on the resp segment (4)
task{1}{1}.getResponse = [0 0 0 1 0];

% Parameters that control the direction of each oriented grating
task{1}{1}.parameter.dir1 = [-1 1];
task{1}{1}.parameter.dir2 = [-1 1];
% Parameter to control which side should be attended on this trial
task{1}{1}.parameter.attend = [1 2];
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
task{1}{1}.randVars.calculated.trialNum = nan;
task{1}{1}.randVars.calculated.dir = nan;
task{1}{1}.randVars.calculated.rot = nan;

stimulus.curTrial = 0;

%% Full Setup
% Initialize task (note phase == 1) (this is for MGL)
for phaseNum = 1:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init staircase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We will make the task consistently difficult by using a staircase, we
% initialize this here.
if stimulus.initStair
    disp(sprintf('(cogneuro_att) Initializing staircases from scratch...'));
    stimulus = initStaircase(stimulus);
else
    disp('(cogneuro_att) Re-using staircase from previous run...');
end

%% Get Ready...
% clear screen
mglTextSet('Helvetica',32,0);
mglWaitSecs(2);
setGammaTable_flowMax(stimulus.contrast);
mglClearScreen(0.5);
if stimulus.scan
    mglTextDraw('DO NOT MOVE',[0 1.5]);
end
mglFlush

% let the user know
disp(sprintf('(cogneuro_att) Starting run number: %i',stimulus.counter));
% block screen refresh
myscreen.flushMode = 1;

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

% Increment the trial number
stimulus.curTrial = stimulus.curTrial + 1;
task.thistrial.trialNum = stimulus.curTrial;

% Store the "correct" direction to attend (could be useful for analysis)
if task.thistrial.attend==1
    task.thistrial.dir = task.thistrial.dir1;
else
    task.thistrial.dir = task.thistrial.dir2;
end

% Get the amount of rotation for this trial from the staircase
[rot, stimulus] = getDeltaPed(stimulus);
% Store this info and tell the operator
task.thistrial.rot = rot;
aT = {'left','right'}; aR = {'right','','left'};
disp(sprintf('(cogneuro_att) Attending: %s, Respond: %s, Rot: %2.2f deg',aT{task.thistrial.attend},aR{task.thistrial.dir+2},rot));

% Now make our gratings (we could save these ahead of time and rotate, but
% they generate in only a few ms so whatever, screw it).
g = mglMakeGrating(10,10,0.5,stimulus.pedestals.angle + rot*task.thistrial.dir1,0);
gauss = mglMakeGaussian(10,10,1.25,1.25);
% normalize
g = (g .* gauss + 1) / 2; % bounded 0-1
g =  (stimulus.colors.nUnreserved-stimulus.colors.nReserved)* g +stimulus.colors.nReserved + 1;
% save texture for faster rendering
stimulus.live.tex1 = mglCreateTexture(g);

g = mglMakeGrating(10,10,0.5,stimulus.pedestals.angle + rot*task.thistrial.dir2,0);
g = (g .* gauss + 1) / 2; % bounded 0-1
g =  (stimulus.colors.nUnreserved-stimulus.colors.nReserved)* g +stimulus.colors.nReserved + 1;
stimulus.live.tex2 = mglCreateTexture(g);
%%%%%%%%%%%%%%%%%%%%%%%%
%    getDeltaPed       %
%%%%%%%%%%%%%%%%%%%%%%%%

function [rot, stimulus] = getDeltaPed(stimulus)
[rot, stimulus.staircase] = doStaircase('testValue',stimulus.staircase);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    runs at the start of each segment       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)
%%
myscreen.flushMode = 0;
global stimulus

stimulus.live.cue = 0;
stimulus.live.grate = 0;
stimulus.live.fixColor = 0;

% All we're going to do is adjust what gets shown based on what segment we
% are in, the actual updating happens in the next function.
switch task.thistrial.thisseg
    case stimulus.seg.cue
        stimulus.live.cue = 1;
    case stimulus.seg.stim
        stimulus.live.grate = 1;
    case stimulus.seg.resp
        stimulus.live.fixColor = 1/255;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus
mglClearScreen(0.5);
% actuall display the grating
if stimulus.live.grate, stimulus = upGrate(stimulus); end
% display the fixation cross (always)
upFix(stimulus);
% display the cue (this draws over the cross)
if stimulus.live.cue, upCue(task); end

function upCue(task)
%%

% some helper code for training
%atex = {'Attend Left', 'Attend Right'};
%mglTextDraw(atex{task.thistrial.attend},[0 1.5]);

% otherwise just draw the line over the fixation cross
if task.thistrial.attend==1
    % left
    mglLines2(0,0,-.75,0,1.5,1/255);
else
    mglLines2(0,0,.75,0,1.5,1/255);
end

function upFix(stimulus)
%%
mglFixationCross(1.5,1.5,stimulus.live.fixColor);


function stimulus = upGrate(stimulus)
% distance is to the center of the grating, in degrees
mglBltTexture(stimulus.live.tex1,[-5 0]);
mglBltTexture(stimulus.live.tex2,[5 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus

% setup
responseText = {'Incorrect','Correct'};
responsePos = {'Left','Right'};
fixColors = {254/255,1};

% check that we got a response and it's a key we want
if any(task.thistrial.whichButton == stimulus.responseKeys)
    % check that we haven't already responded
    if task.thistrial.gotResponse < 2
        % determine which button is correct
        if task.thistrial.dir==-1
            cButt = 2;
        else
            cButt = 1;
        end
        % check if subject was correct
        corr = task.thistrial.whichButton == cButt;
        % store this info
        task.thistrial.correct = corr;
        stimulus.staircase = doStaircase('update',stimulus.staircase,corr);
        stimulus.live.fixColor = fixColors{corr+1};
        disp(sprintf('(cogneuro_att) Response %s: %s',responseText{corr+1},responsePos{find(stimulus.responseKeys==task.thistrial.whichButton)}));
    else
        disp(sprintf('(cogneuro_att) Subject responded multiple times: %i',task.thistrial.gotResponse+1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              HELPER FUNCTIONS                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%
%    initStaircase     %
%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initStaircase(stimulus)
%%
stimulus.staircase = doStaircase('init','upDown',...
        'initialThreshold',stimulus.pedestals.initThresh.angle,...
        'initialStepsize',stimulus.pedestals.initThresh.angle/3,...
        'minThreshold=0.001','maxThreshold=20','stepRule','pest', ...
        'nTrials=50','maxStepsize=5','minStepsize=.001');
%%

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