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

stimulus.gaussSz = 1:5;

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
if ~isfield(stimulus,'staircase')
    disp('(afmapb) WARNING: New staircase');
    initStair();
else
    resetStair();
end

%% Plot and return
if stimulus.plots==2
    dispInfo(stimulus);
    return
end

%% Initialize Stimulus
    
stimulus.responseKeys = [1 2]; % left right

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

%% Setup Probe Task

%% Setup Attention Task

stimulus.curTrial = 0;

task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;
% task waits for fixation on first segment
task{1}{1}.segmin = [1 0.500 0.200 1];
task{1}{1}.segmax = [1 0.500 0.200 1];

stimulus.seg.ITI = 1;
stimulus.seg.delay1 = 2;
stimulus.seg.stim = 3;
stimulus.seg.resp = 4;

task{1}{1}.synchToVol = [0 0 0 0];
task{1}{1}.getResponse = [0 0 0 1];

task{1}{1}.numTrials = Inf;

task{1}{1}.parameter.present = [0 1];

task{1}{1}.random = 1;

task{1}{1}.randVars.calculated.contrast = nan;
task{1}{1}.randVars.calculated.resp = nan;
task{1}{1}.randVars.calculated.correct = nan;

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

stimulus.live.fixColor = stimulus.colors.white;

if task.thistrial.thisseg==stimulus.seg.ITI
    stimulus.live.fixColor = stimulus.colors.black;
end

function [task, myscreen] = startTrialCallback1(task,myscreen)
%%
global stimulus

stimulus.curTrial = stimulus.curTrial+1;
[task.thistrial.contrast, stimulus.staircase] = doStaircase('testValue',stimulus.staircase);

stimulus.live.wn = mglCreateTexture(repmat(128+randi(127,1,myscreen.screenWidth,myscreen.screenHeight,'uint8'),4,1,1));

stimulus.gaussian = {};
sz = 5;
gauss = mglMakeGaussian(sz,sz,sz/6,sz/6);
alphamask = repmat(255*ones(size(gauss)),1,1,4);
alphamask(:,:,4) = gauss*255*task.thistrial.contrast;
stimulus.gaussian = mglCreateTexture(alphamask);

disp(sprintf('(afmapb) Trial %i: %02.1f',stimulus.curTrial,task.thistrial.contrast*100));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback1(task, myscreen)
%%
global stimulus

mglClearScreen(0.5);

% draw gratings for probe task
stimulus.live.aX = 5;
stimulus.live.aY = 5;

if (task.thistrial.thisseg==stimulus.seg.stim)
    mglBltTexture(stimulus.live.wn,[0 0]);
    
    if task.thistrial.present
        mglBltTexture(stimulus.gaussian,[stimulus.live.aX,stimulus.live.aY]);
    end
elseif (task.thistrial.thisseg==stimulus.seg.ITI)
    for i = 1:8
        mglGluPartialDisk(stimulus.live.aX,stimulus.live.aY,0.99,1.01,(i-1)*360/8-11.25,360/16,stimulus.colors.white);
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
        task.thistrial.resp = task.thistrial.whichButton-1;
        
        task.thistrial.correct = task.thistrial.resp==task.thistrial.present;
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
            'initialThreshold',0.5,...
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
%  rrun # counter #    local trial     real trial   angle     respAngle    
%     1       2             3              4           5           6
%  target    startRespAngle     contrast     detected      ecc    priorsd
%     7            8                9           10          11      12
%    rotation
%       13
% count = 1; data = zeros(10000,13);
% 
% for fi = 1:length(files)
%     load(fullfile(sprintf('~/data/afmapb_%s/%s/%s',rstimulus.condition,mglGetSID,files(fi).name)));
%     
%     e = getTaskParameters(myscreen,task);
%     if e{1}(1).nTrials>1
%         e = e{1}(2); % why?!
%     
%         run = stimulus.counter;
% 
%         data(count:count+(e.nTrials-1),:) = [repmat(fi,e.nTrials,1) repmat(run,e.nTrials,1) (1:e.nTrials)' (count:count+(e.nTrials-1))' ...
%             e.randVars.angle' e.randVars.respAngle' e.parameter.target' ...
%             e.randVars.startRespAngle' e.randVars.contrast' e.randVars.detected' ...
%             e.parameter.ecc' e.parameter.priorSTD' e.randVars.rotation'];
% 
%         count = count+e.nTrials;
%     end
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localInitStimulus()

global stimulus

for ci = 1:length(stimulus.gratingContrasts)
    for si = 1:length(stimulus.gratingSizes)
        sz = 1.5 * stimulus.gratingSizes(si);
        % use total degs / num to compute size
        grating = stimulus.gratingContrasts(ci) * 255/2 * mglMakeGrating(sz,sz,2,0) + 255/2;
        gauss = mglMakeGaussian(sz,sz,sz/6,sz/6);
        alphamask = repmat(grating,1,1,4);
        alphamask(:,:,4) = gauss*255;

        stimulus.live.grating(ci,si)  = mglCreateTexture(alphamask); % high contrast
    end
end
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