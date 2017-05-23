
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

function [myscreen] = cogneuro_attention_eshed(varargin)

global stimulus
clear fixStimulus
%% Initialize Variables

disp('****************************************');
disp('** NEPR Into Cog Neuro Attention Task **');
disp('****************************************');
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
stimulus.colors.reservedBottom = [1 0 0; 0 0 0]; % fixation cross colors
stimulus.colors.reservedTop = [1 1 1; 0 1 0]; % correct/incorrect colors
stimulus.colors.black = 1/255; stimulus.colors.white = 254/255;
stimulus.colors.red = 0/255; stimulus.colors.green = 255/255;
stimulus.colors.nReserved = 2; % this is /2 the true number, because it's duplicated
stimulus.colors.nUnreserved = 256-(2*stimulus.colors.nReserved);

stimulus.colors.mrmax = stimulus.colors.nReserved - 1 + stimulus.colors.nUnreserved;
stimulus.colors.mrmin = stimulus.colors.nReserved;


%% Stimulus (oriented gratings) properties
stimulus.angleOpts = [30 120];
stimulus.shiftOpts = [3 6];

% set stimulus contrast, scanner absolute luminance is lower so the
% relative contrast needs to be higher to make up for this
if stimulus.scan
    stimulus.contrast = .25; % 4%
else
    stimulus.contrast = 1;
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

if ~isempty(mglGetSID) && isdir(sprintf('~/data/cogneuro_attention_eshed/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/cogneuro_attention_eshed/%s/1*mat',mglGetSID));

    if length(files) >= 1
        fname = files(end).name;
        
        s = load(sprintf('~/data/cogneuro_attention_eshed/%s/%s',mglGetSID,fname));
        % copy staircases and run numbers
        stimulus.counter = s.stimulus.counter + 1;
        clear s;
        disp(sprintf('(cn_eshed) Data file: %s loaded.',fname));
    end
end
disp(sprintf('(cn_eshed) This is run #%i',stimulus.counter));

%% Initialize Stimulus
% Setup for MGL
myscreen = initStimulus('stimulus',myscreen);
stimulus.responseKeys = [1 2]; % corresponds to LEFT - RIGHT

%% Setup Task

% Blocks the task so it doesn't start immediately
task{1}{1}.waitForBacktick = 1;

% Storing some helper info for callbacks later
stimulus.seg.wait = 1;
stimulus.seg.cue = 2;
stimulus.seg.isi = 3;
stimulus.seg.stim = 4;
stimulus.seg.resp = 5;
stimulus.seg.ITI = 6;

% Trial timing (see above for what each column corresponds to)
task{1}{1}.segmin = [0.500 0.500 1.000 4.000 1.500 2.500];
task{1}{1}.segmax = [0.500 0.500 1.000 4.000 1.500 9.500]; % 6.5 s average 

% When scanning we synchronize the stimulus to the scanner
task{1}{1}.synchToVol = zeros(1,length(task{1}{1}.segmin));
if stimulus.scan
    task{1}{1}.synchToVol(end) = 1;
end

% Get responses on the resp segment (4)
task{1}{1}.getResponse = zeros(length(task{1}{1}.segmin));
task{1}{1}.getResponse(stimulus.seg.resp) = 1;

% Parameters that control the direction of each oriented grating
task{1}{1}.parameter.rotL = [0 1];
task{1}{1}.parameter.attend = [1 2];
% Parameter to control which side should be attended on this trial
task{1}{1}.parameter.change = [0 0 1 1 1 1 1 1 1 1]; % 80% of the time a change occurs
task{1}{1}.parameter.valid = [0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; % 85% valid cue
task{1}{1}.parameter.diff = [1 2];
task{1}{1}.parameter.dir = [-1 1];
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
task{1}{1}.randVars.calculated.rotR = nan;
task{1}{1}.randVars.calculated.changeSide = nan;

stimulus.curTrial = 0;

%% Gratings

sz = 8;
% % use total degs / num to compute size
% grating = 255/2 * mglMakeGrating(sz,sz,0.5,0) + 255/2;
% gauss = mglMakeGaussian(sz,sz,sz/6,sz/6);
% alphamask = repmat(grating,1,1,4);
% alphamask(:,:,4) = gauss*255;
% 
% stimulus.live.grating  = mglCreateTexture(alphamask); % high contrast

g = mglMakeGrating(sz*1.5,sz*1.5,0.5,0,0);
gauss = mglMakeGaussian(sz*1.5,sz*1.5,sz/6,sz/6);
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
mglTextSet('Helvetica',32,0);
mglClearScreen(0.5);
mglFixationCross(1,1,stimulus.colors.black);
mglFlush
mglClearScreen(0.5);
mglFixationCross(1,1,stimulus.colors.black);
mglFlush
setGammaTable_flowMax(stimulus.contrast);

% let the user know
disp(sprintf('(cn_eshed) Starting run number: %i',stimulus.counter));

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

task.thistrial.rotR = ~task.thistrial.rotL;

if task.thistrial.valid
    task.thistrial.changeSide = task.thistrial.attend;
else
    flip = [2 1];
    task.thistrial.changeSide = flip(task.thistrial.attend);
end
 
sides = {'left','right'};
if task.thistrial.change
    changeText = sides{task.thistrial.changeSide};
else
    changeText = 'none';
end
disp(sprintf('Trial %i: cue %s, change %s, difficulty +%i',stimulus.curTrial,sides{task.thistrial.attend},changeText,stimulus.shiftOpts(task.thistrial.diff)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    runs at the start of each segment       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)
%%
myscreen.flushMode = 0;
global stimulus

stimulus.live.cue = 0;
stimulus.live.grate = 0;
stimulus.live.fixColor = stimulus.colors.white;

stimulus.live.rotL = task.thistrial.rotL;
stimulus.live.rotR = task.thistrial.rotR;
% All we're going to do is adjust what gets shown based on what segment we
% are in, the actual updating happens in the next function.
switch task.thistrial.thisseg
    case stimulus.seg.cue
        stimulus.live.cue = 1;
    case stimulus.seg.stim
        stimulus.live.grate = 1;
        stimulus.live.rotLex = 0;
        stimulus.live.rotRex = 0;
    case stimulus.seg.resp
        stimulus.live.grate = 1;
        if task.thistrial.change && (task.thistrial.changeSide==1)
            stimulus.live.rotLex = stimulus.shiftOpts(task.thistrial.diff);
        elseif task.thistrial.change
            stimulus.live.rotRex = stimulus.shiftOpts(task.thistrial.diff);
        else
            % no change is occurring
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
% actuall display the grating
if stimulus.live.grate, upGrate(task); end
% display the fixation cross (always)
upFix(stimulus);
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
else
    mglLines2(0,0,.5,0,1,stimulus.colors.black);
end

function upFix(stimulus)
%%
mglFixationCross(1,1,stimulus.live.fixColor);

function upGrate(task)
global stimulus
% distance is to the center of the grating, in degrees
ctime = ceil(10*mglGetSecs-task.thistrial.segStartSeconds);

if mod(ctime,2)==1 % if time is a multiple of 100 ms display, otherwise no, 10 hz flip rate
    mglBltTexture(stimulus.live.grating,[-6 0],0,0,stimulus.angleOpts(stimulus.live.rotL+1)+task.thistrial.dir*stimulus.live.rotLex);
    mglBltTexture(stimulus.live.grating,[6 0],0,0,stimulus.angleOpts(stimulus.live.rotR+1)+task.thistrial.dir*stimulus.live.rotRex);
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
        task.thistrial.correct = task.thistrial.whichButton == task.thistrial.changeSide;

        stimulus.live.fixColor = fixColors{task.thistrial.correct+1};
        disp(sprintf('(cn_eshed) Response %s: %s',responseText{task.thistrial.correct+1},responsePos{stimulus.responseKeys(task.thistrial.whichButton)}));
    else
        disp(sprintf('(cn_eshed) Subject responded multiple times: %i',task.thistrial.gotResponse+1));
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