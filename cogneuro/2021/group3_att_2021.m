% group3_att_2021 
%
%      usage: myscreen=group3_att_2021()
%         by: akshay jagadeesh
%       date: 04/05/2021
%    purpose: contrast change detection with cued selective attention.
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

function [myscreen] = group3_att_2021(varargin)

global stimulus
clear fixStimulus
%% Initialize Variables

disp('****************************************');
disp('** NEPR Into Cog Neuro Attention (Orientation Discrimination) Task **');
disp('****************************************');
% add arguments later
stimFileNum = [];
scan = [];
testing = [];
getArgs(varargin,{'stimFileNum=-1','testing=0','scan=0'});
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
stimulus.colors.reservedBottom = [0 0 0; .6 .6 .6]; % fixation cross colors
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
    stimulus.contrast = .04; % 4%
else
    stimulus.contrast = .15;
end

%% Setup Screen
% Settings for Oban/Psychophysics room launch
myscreen.background = 128/255;
if stimulus.scan
    myscreen.displayName = 'fMRIprojFlex';
else
    myscreen.displayName = 'VPixx2';
end
myscreen=initScreen(myscreen);
%% Open Old Stimfile
% Do not modify this code
stimulus.initStair = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/group3_att_2021/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/group3_att_2021/%s/21*mat',mglGetSID));
    
    if length(files) >= 1
        if stimFileNum == -1
            if length(files) > 1
                warning('Multiple stimfiles found, loading last one. If you wanted a different functionality use stimFileNum=#');
            end
            fname = files(end).name;
        else
            fname = files(stimFileNum).name;
        end
        s = load(sprintf('~/data/group3_att_2021/%s/%s',mglGetSID,fname));
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
stimulus.seg.stim = 1;
stimulus.seg.change = 2;

% Trial timing (see above for what each column corresponds to)
task{1}{1}.segmin = [1 2];
task{1}{1}.segmax = [1 2];

% When scanning we synchronize the stimulus to the scanner
task{1}{1}.synchToVol = zeros(task{1}{1}.segmin);
if stimulus.scan
    task{1}{1}.synchToVol(end) = 1;
    task{1}{1}.segmin(end) = max(0, task{1}{1}.segmin(end) - 0.100);
    task{1}{1}.segmax(end) = max(0, task{1}{1}.segmax(end) - 0.200);

end

% Get responses on the resp segment (4)
task{1}{1}.getResponse = [0 1];

% Parameters that control the direction of each oriented grating
stimulus.nGratings = 8;
stimulus.gratingRadius = 10;
for i = 1:stimulus.nGratings
    task{1}{1}.parameter.(sprintf('dir%i', i)) = [-1,1];
    task{1}{1}.parameter.(sprintf('change%i', i)) = [0,1];
end

% Parameter to control which side should be attended on this trial
% This will randomize trials
task{1}{1}.random = 1;
% Outside the scanner fix the trial count
task{1}{1}.numTrials = 80;
% Inside the scanner, inf length so that we can stop the stimulus after the
% scanner stops running.

%% Tracking
% these are variables that we want to track for later analysis.
task{1}{1}.randVars.calculated.attend = nan;
task{1}{1}.randVars.calculated.cuedChange = nan;
task{1}{1}.randVars.calculated.correct = nan;
task{1}{1}.randVars.calculated.trialNum = nan;
task{1}{1}.randVars.calculated.dir = nan;
task{1}{1}.randVars.calculated.rot = nan;

stimulus.curTrial = 0;
stimulus.fixWidth = 1;

%% Full Setup
for i = 1:2
  mglFixationCross(stimulus.fixWidth, 1, [1 1 1]);
  mglFlush;
end
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
task.thistrial.attend = 1+mod(task.trialnum, stimulus.nGratings); % Go through each grating one by one.
task.thistrial.dir = task.thistrial.(sprintf('dir%i', task.thistrial.attend));
task.thistrial.cuedChange = task.thistrial.(sprintf('change%i', task.thistrial.attend));

% Get the amount of rotation for this trial from the staircase
[rot, stimulus] = getDeltaPed(stimulus);
% Store this info and tell the operator
task.thistrial.rot = rot;
aT = {'left','right'}; aR = {'right','','left'};
disp(sprintf('(cogneuro_att) Attending: %i; Did cued change: %i',task.thistrial.attend, task.thistrial.cuedChange));

% Now make our gratings (we could save these ahead of time and rotate, but
gauss = mglMakeGaussian(5,5,1,1);
for i = 1:stimulus.nGratings
  g = mglMakeGrating(5,5,0.5, stimulus.pedestals.angle + rot * task.thistrial.(sprintf('dir%i', i)),0);
  g = (g .* gauss + 1) / 2; % bounded 0-1
  g =  (stimulus.colors.nUnreserved-stimulus.colors.nReserved)* g +stimulus.colors.nReserved + 1;
  stimulus.live.(sprintf('tex%i_1',i)) = mglCreateTexture(g);

  if task.thistrial.(sprintf('change%i', i))
    g = mglMakeGrating(5,5,0.5, stimulus.pedestals.angle + 2*rot * task.thistrial.(sprintf('dir%i', i)),0);
    g = (g .* gauss + 1) / 2; % bounded 0-1
    g =  (stimulus.colors.nUnreserved-stimulus.colors.nReserved)* g +stimulus.colors.nReserved + 1;
  end
  stimulus.live.(sprintf('tex%i_2',i)) = mglCreateTexture(g);
end

if task.thistrial.cuedChange == 1
    task.thistrial.correct = 0;
else
    task.thistrial.correct = 1;
end

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

stimulus.live.cue = 1;
stimulus.live.grate1 = 0;
stimulus.live.grate2 = 0;
stimulus.live.fixColor = 0;
stimulus.live.green = [0 1 0];
stimulus.live.red = [1 0 0];

% All we're going to do is adjust what gets shown based on what segment we
% are in, the actual updating happens in the next function.
switch task.thistrial.thisseg
    case stimulus.seg.stim
        stimulus.live.grate1 = 1;
    case stimulus.seg.change
        stimulus.live.grate2 = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus
mglClearScreen(0.5);
% display the fixation cross (always)
upFix(stimulus);

% display the grating
if stimulus.live.grate1
    upGrate1(stimulus);
elseif stimulus.live.grate2
    upGrate2(stimulus);
end
% display the cue (this draws over the cross)
if stimulus.live.cue
    upCue(task, stimulus); 
end


%%%%%%%%%%%%
function upCue(task, stimulus)

theta = (task.thistrial.attend-1)*(2*pi/stimulus.nGratings);
x=[1,.75,2,.75];
y=[0,.5,0,-.5];
mglPolygon(cos(theta)*x - sin(theta)*y, sin(theta)*x + cos(theta)*y, [0,0,0])

%%%%%%%%%%%%
function upFix(stimulus, fixColor)

if ieNotDefined('fixColor')
  fixColor = stimulus.live.fixColor;
end
mglFixationCross(1.5,1.5,fixColor);

%%%%%%%%%%%%
function upGrate1(stimulus)
rad = stimulus.gratingRadius;
for i = 1:stimulus.nGratings
  theta = (i-1)*(2*pi/stimulus.nGratings);
  mglBltTexture(stimulus.live.(sprintf('tex%i_1',i)), [rad*cos(theta),rad*sin(theta)]);
end

%%%%%%%%%%%%
function upGrate2(stimulus)
rad = stimulus.gratingRadius;
for i = 1:stimulus.nGratings
  theta = (i-1)*(2*pi/stimulus.nGratings);
  mglBltTexture(stimulus.live.(sprintf('tex%i_2',i)), [rad*cos(theta),rad*sin(theta)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus

% setup
responseText = {'Incorrect','Correct'};
responsePos = {'Left','Right'};
fixColors = {stimulus.live.red,stimulus.live.green};

% check that we got a response and it's a key we want
if any(task.thistrial.whichButton == stimulus.responseKeys)
    % check that we haven't already responded
    if task.thistrial.gotResponse < 2
        % determine which button is correct
        if task.thistrial.cuedChange == 1
            task.thistrial.correct = 1;
        else
            task.thistrial.correct = 0;
        end
        corr = task.thistrial.correct;
        stimulus.staircase = doStaircase('update',stimulus.staircase,corr);
        stimulus.live.fixColor = fixColors{corr+1};
        disp(sprintf('(cogneuro_att) Response %s',responseText{corr+1}));
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
