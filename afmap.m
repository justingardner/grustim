function [ myscreen ] = afmap( varargin )
%ATTENTIONFIELDMAPPING 
%
%   Map the attention field in the scanner. This function works by having a
%   participant perform an asynchronous attention task at fixation or in a
%   quarterfield region. A pre-determined poisson process generates random
%   flashes of rotating gratings throughout the visual field at low or high
%   contrast.
%
%   The probe stimuli have three sizes 1x1 deg, 2x2 deg, or 4x4 deg, to
%   help estimate different RF sizes and are placed at 2 degree
%   increments. The probability of a probe stimulus turning on is 2.5% per
%   TR and the dead time is five seconds. Using poisson processes in this
%   way ensures a random distribution at every location that is
%   uncorrelated to all other locations. The code starts with a 10 s blank
%   and has another 10 s blank every three minutes. Probes are at 20% and
%   80% contrast, each probe lasts two TRs.
%   
%   The attention task involves performing orientation judgments on gabors
%   at a location cued continuously by a circular aperture. Fixation is
%   maintained at the center and monitored within a 1.5 deg window. Gabor
%   is at full contrast to differentiate from the probe stimuli (and to
%   minimize the effect of the probes on performance)
%
%   The goal of the task is to provide an independent pRF dataset that can
%   be used to build channel encoding models that cover space, contrast,
%   and attention. The end goal is to link frontal areas that control
%   attention with visual areas. 
%
%   Project goals:
%       (1) Show that we can map pRFs hierarchically throughout all of
%       cortex
%       (2) Build a visuospatial encoding model that incorporates attention
%       (3) Show that our model predicts tuning shifts for non-visuospatial
%       features due to spatial attention

%%

global stimulus

stimulus = struct;

%% Initialize Variables

% add arguments later
scan = 0;
plots = 0;
noeye = 0;
debug = 0;
getArgs(varargin,{'scan=1','plots=0','noeye=0','debug=0'});
stimulus.scan = scan;
stimulus.plots = plots;
stimulus.noeye = noeye;
stimulus.debug = debug;
clear localizer invisible scan noeye task test2

if stimulus.scan
    warning('Not setup for scanning');
end

%% Stimulus parameters

stimulus.stimX = 25; % max ecc in any direction
stimulus.stimY = 13;
stimulus.stimR = 2; % deg between each stimulus

if mod(stimulus.stimX,stimulus.stimR)==1 || mod(stimulus.stimY,stimulus.stimR)==1
    warning('Your stimulus size is not correctly setup');
end

stimulus.stimx = -stimulus.stimX:stimulus.stimR:stimulus.stimX;
stimulus.stimy = -stimulus.stimY:stimulus.stimR:stimulus.stimY;

stimulus.probeOn = .002;
stimulus.live.probeOnGrid = zeros(length(stimulus.stimx),length(stimulus.stimy));
% stimulus.probeMaxLag = 30;
stimulus.probeUp = 4;
stimulus.probeDown = 12;

stimulus.gridCount = 1;
stimulus.grid.on = zeros(5000,length(stimulus.stimx),length(stimulus.stimy));
stimulus.grid.con = stimulus.grid.on;
stimulus.grid.sz = stimulus.grid.on;

stimulus.live.grid = zeros(length(stimulus.stimx),length(stimulus.stimy));
stimulus.live.gridCon = zeros(length(stimulus.stimx),length(stimulus.stimy));
stimulus.live.gridSize = zeros(length(stimulus.stimx),length(stimulus.stimy));
stimulus.live.gridPhase = zeros(length(stimulus.stimx),length(stimulus.stimy));

stimulus.gratingContrasts = [0.1 1.0 1.0];
stimulus.live.gridCons = zeros(1,length(stimulus.gratingContrasts)-1);
stimulus.gratingSizes = [0.5 1 2];
stimulus.live.gridSizes = zeros(1,length(stimulus.gratingSizes));

stimulus.live.rotations = zeros(length(stimulus.stimx),length(stimulus.stimy));

stimulus.live.attend = 0;

stimulus.blanks.rotTR = 20; % every 20 TRs (10s) rotate through blanks
stimulus.blanks.names = {'None','NW','NE','SE','SW','None','All'};
stimulus.blanks.xmin = [-inf 0 0 -inf 0 -inf];
stimulus.blanks.xmax = [0 inf inf 0 0 inf];
stimulus.blanks.ymin = [0 0 -inf -inf 0 -inf];
stimulus.blanks.ymax = [inf inf 0 0 0 inf];
stimulus.live.cBlank = 0;
stimulus.live.rotation = 1;
stimulus.live.blankTime = stimulus.blanks.rotTR;
stimulus.live.numBlanks = length(stimulus.blanks.ymax);

%% Open Old Stimfile
stimulus.counter = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/afmap/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/afmap/%s/1*mat',mglGetSID));

    if length(files) >= 1
        fname = files(end).name;
        
        s = load(sprintf('~/data/afmap/%s/%s',mglGetSID,fname));
        % copy staircases and run numbers
        stimulus.counter = s.stimulus.counter + 1;
        stimulus.staircase = s.stimulus.staircase;
        stimulus.live.attend = mod(s.stimulus.live.attend+1,3);
        clear s;
        disp(sprintf('(afmap) Data file: %s loaded.',fname));
    end
end
disp(sprintf('(afmap) This is run #%i',stimulus.counter));

%% Setup attention

stimulus.attendX = [0 5 5];
stimulus.attendY = [0 5 -5];
stimulus.live.aX = stimulus.attendX(stimulus.live.attend+1);
stimulus.live.aY = stimulus.attendY(stimulus.live.attend+1);

%% Setup Screen
myscreen = initScreen('VPixx');

% set background to grey
myscreen.background = 0.5;

%% Staircase
if ~isfield(stimulus,'staircase')
    disp('(afmap) WARNING: New staircase');
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

myscreen.stimulusNames{1} = 'stimulus';

if ~isfield(stimulus.live,'grating')
    localInitStimulus();
end
    
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

task{1}{1} = struct;
task{1}{1}.waitForBacktick = 1;
% task waits for fixation on first segment
if stimulus.scan
    task{1}{1}.seglen = 0.450;
else
    task{1}{1}.seglen = 0.500;
end

stimulus.seg.stim = 1;

if stimulus.scan
    task{1}{1}.synchToVol = 1;
end
task{1}{1}.getResponse = 0;

task{1}{1}.numTrials = Inf;

task{1}{1}.random = 0;

if stimulus.scan
    task{1}{1}.synchToVol = 1;
end

task{1}{1}.randVars.calculated.probesOn = nan;

%% Setup Attention Task

stimulus.curTrial = 0;

global fixStimulus

fixStimulus.diskSize = 0.45;
fixStimulus.fixWidth = 0.4;
fixStimulus.fixLineWidth = 1;
[task{2}, myscreen] = fixStairInitTask(myscreen);
% task{2}{1} = struct;
% task{2}{1}.waitForBacktick = 0;
% % task waits for fixation on first segment
% task{2}{1}.segmin = [0.500 1 0.200 1];
% task{2}{1}.segmax = [2.500 1 0.200 1];
% 
% stimulus.seg.ITI = 1;
% stimulus.seg.delay1 = 2;
% stimulus.seg.stim = 3;
% stimulus.seg.resp = 4;
% 
% task{2}{1}.synchToVol = [0 0 0 0];
% task{2}{1}.getResponse = [0 0 0 1];
% 
% task{2}{1}.numTrials = Inf;
% 
% task{2}{1}.parameter.rotation = [-1 1];
% 
% task{2}{1}.random = 1;
% 
% if stimulus.scan
%     task{2}{1}.synchToVol = 1;
% end
% 
% task{2}{1}.randVars.calculated.resp = nan;
% task{2}{1}.randVars.calculated.correct = nan;

%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback1,@screenUpdateCallback1,[],@startTrialCallback1,[],[]);
%     [task{2}{phaseNum}, myscreen] = initTask(task{2}{phaseNum},myscreen,@startSegmentCallback2,@screenUpdateCallback2,@getResponseCallback2,@startTrialCallback2,[],[]);
end

%% EYE CALIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

% let the user know
disp(sprintf('(afmap) Starting run number: %i.',stimulus.counter));

%% Main Task Loop

% setGammaTable(1);
mglClearScreen(0.5); %mglFixationCross(1,1,stimulus.colors.white);

mglFlush
mglClearScreen(0.5); %mglFixationCross(1,1,stimulus.colors.white);

phaseNum = 1;
% Again, only one phase.
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    % update the task
    [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);
    [task{2}, myscreen, phaseNum] = updateTask(task{2},myscreen,phaseNum);
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

stimulus.grid.on = stimulus.grid.on(1:(stimulus.gridCount-1),:,:);
stimulus.grid.con = stimulus.grid.con(1:(stimulus.gridCount-1),:,:);
stimulus.grid.sz = stimulus.grid.sz(1:(stimulus.gridCount-1),:,:);

if stimulus.plots
    disp('(afmap) Displaying plots');
    dispInfo(stimulus);
end

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback1(task,myscreen)
% pass


function [task, myscreen] = startTrialCallback1(task,myscreen)
%%
global stimulus

% Deal with blanks
stimulus.live.blankTime = stimulus.live.blankTime - 1;
if stimulus.live.blankTime == 0
    stimulus.live.blankTime = stimulus.blanks.rotTR;
    stimulus.live.cBlank = stimulus.live.cBlank + 1;
    if stimulus.live.cBlank > stimulus.live.numBlanks
        stimulus.live.cBlank = 0;
        stimulus.live.rotation = stimulus.live.rotation + 1;
        disp(sprintf('Rotation %i starting',stimulus.live.rotation));
    end
end

if stimulus.live.rotation > 6
    disp(sprintf('All rotations complete'));
    return
end

% Design the state space
for x = 1:length(stimulus.stimx)
    for y = 1:length(stimulus.stimy)
        % increment
        if stimulus.live.grid(x,y) > 1
            stimulus.live.grid(x,y) = stimulus.live.grid(x,y)-1;
            stimulus.live.gridPhase(x,y) = ~stimulus.live.gridPhase(x,y);
        elseif stimulus.live.grid(x,y) == 1
            % shut down grid location
            stimulus.live.grid(x,y) = 0;
            stimulus.live.gridCon(x,y) = 0;
            stimulus.live.gridSize(x,y) = 0;
            stimulus.live.rotations(x,y) = 0;
            stimulus.live.probeOnGrid(x,y) = 0;
        else
            xp = stimulus.stimx(x); yp = stimulus.stimy(y);
            if (stimulus.live.cBlank>0)
                nocheck = (xp<stimulus.blanks.xmax(stimulus.live.cBlank) && xp > stimulus.blanks.xmin(stimulus.live.cBlank)) && (yp<stimulus.blanks.ymax(stimulus.live.cBlank) && yp>stimulus.blanks.ymin(stimulus.live.cBlank));
            else
                nocheck = false;
            end
            % if nocheck is true we skip this position (i.e. don't turn on in
            % blank locations)
            
            if ~nocheck
                probeOn = stimulus.live.probeOnGrid(x,y);
                stimulus.live.probeOnGrid(x,y) = min(1,stimulus.live.probeOnGrid(x,y)+stimulus.probeOn);
                if rand < probeOn
                    % turn on grid location
                    stimulus.live.grid(x,y) = stimulus.probeDown + stimulus.probeUp;
                    % pick attributes
                    conOpts = osum(stimulus.live.gridCons);
                    conOpts = cumsum(conOpts)/sum(conOpts);
                    r = rand;
                    conChoice = find(r<conOpts,1);
                    stimulus.live.gridCons(conChoice) = stimulus.live.gridCons(conChoice) + 1;
                    stimulus.live.gridCon(x,y) = conChoice;

                    sizeOpts = osum(stimulus.live.gridSizes);
                    sizeOpts = cumsum(sizeOpts)/sum(sizeOpts);
                    r = rand;
                    sizeChoice = find(r<sizeOpts,1);
                    stimulus.live.gridSizes(sizeChoice) = stimulus.live.gridSizes(sizeChoice) + 1;
                    stimulus.live.gridSize(x,y) = sizeChoice;

                    stimulus.live.rotations(x,y) = rand*2*pi;
                    
                    stimulus.live.gridPhase(x,y) = ~stimulus.live.gridPhase(x,y);
                end
            else
                
            end
        end
    end
end

stimulus.grid.on(stimulus.gridCount,:,:) = stimulus.live.grid;
stimulus.grid.con(stimulus.gridCount,:,:) = stimulus.live.gridCon;
stimulus.grid.sz(stimulus.gridCount,:,:) = stimulus.live.gridSize;
stimulus.gridCount = stimulus.gridCount + 1;

task.thistrial.probesOn = sum(stimulus.live.grid(:)>stimulus.probeDown);
    

disp(sprintf('(afmap) Probes: %i, current blank: %s',task.thistrial.probesOn,stimulus.blanks.names{stimulus.live.cBlank+1}));
% REFRESH THE SCREEN

mglClearScreen();

% draw gratings for probe task

for xi = 1:length(stimulus.stimx)
    for yi = 1:length(stimulus.stimy)
        if stimulus.live.grid(xi,yi) > stimulus.probeDown
            x = stimulus.stimx(xi);
            y = stimulus.stimy(yi);
            con = stimulus.live.gridCon(xi,yi);
            sz = stimulus.live.gridSize(xi,yi);
            ph = stimulus.live.gridPhase(xi,yi)+1;
            
            mglBltTexture(stimulus.grating(con,sz,ph),[x y],0,0,stimulus.live.rotations(xi,yi)*180/pi);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback1(task, myscreen)
%%
global stimulus



function [task, myscreen] = startTrialCallback2(task,myscreen)
%%
global stimulus

stimulus.curTrial = stimulus.curTrial + 1;

[rotation, stimulus.staircase] = doStaircase('testValue',stimulus.staircase);
task.thistrial.rotation = rotation*task.thistrial.rotation;

rotText = {'-','+'};
disp(sprintf('Trial %i: %s%01.2f',stimulus.curTrial,rotText{1+(task.thistrial.rotation>0)},abs(task.thistrial.rotation)*180/pi));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback2(task, myscreen)
%%

global stimulus

stimulus.live.fix = 1;
stimulus.live.stim = 0;
stimulus.live.fixColor = stimulus.colors.white;

if task.thistrial.thisseg==stimulus.seg.ITI
    stimulus.live.fixColor = stimulus.colors.black;
elseif task.thistrial.thisseg == stimulus.seg.stim
    stimulus.live.rotation = task.thistrial.rotation*180/pi;
    stimulus.live.stim = 1;
% elseif task.thistrial.thisseg == stimulus.seg.stim2
%     stimulus.live.rotation = task.thistrial.rotation*180/pi;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback2(task, myscreen)
%%
global stimulus

% mglClearScreen();

if stimulus.live.fix
    upFix(stimulus);
end

if stimulus.live.stim
    upStim(stimulus);
end

upAttend(stimulus);

function upAttend(stimulus)
%%
for i = 1:8
    mglGluPartialDisk(stimulus.live.aX,stimulus.live.aY,0.99,1.01,(i-1)*360/8-11.25,360/16,stimulus.colors.white);
end

function upFix(stimulus)
%%
% for this experiment use a circle to indicate where participants can
% fixate inside of (rather than a cross which might arbitrarily enforce
% poisitioning
% mglGluAnnulus(0,0,1.5,1.55,stimulus.live.fixColor,64);
mglFixationCross(1,1,stimulus.live.fixColor);

function upStim(stimulus)

mglBltTexture(stimulus.grating(3,1),[stimulus.live.aX,stimulus.live.aY],0,0,stimulus.live.rotation+90);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [task, myscreen] = getResponseCallback2(task, myscreen)

global stimulus

colors = [stimulus.colors.red;stimulus.colors.green];
text = {'Incorrect','Correct'};
stext = {'Left','Right'};
if any(task.thistrial.whichButton==stimulus.responseKeys)
    if task.thistrial.gotResponse==0
        task.thistrial.correct = (task.thistrial.whichButton==1 && task.thistrial.rotation>0) || (task.thistrial.whichButton==2 && task.thistrial.rotation<0);
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
            'initialThreshold',0.15,...
            'initialStepsize',0.03,...
            'minThreshold=0.0001','maxThreshold=0.4','stepRule','pest',...
            'nTrials=80','maxStepsize=0.2','minStepsize=0.0001');
        
function resetStair()

global stimulus

if doStaircase('stop',stimulus.staircase)
    disp('(afmap) Staircase is being reset');
    stimulus.staircase(end+1) = doStaircase('init',stimulus.staircase(end));
    if stimulus.staircase(end).s.threshold>0.3
        disp('(afmap) Bad staircase threshold: setting to 0.3');
        stimulus.staircase(end).s.threshold=0.3;
    elseif stimulus.staircase(end).s.threshold<0
        disp('(afmap) Bad staircase threshold: setting to 0.05');
        stimulus.staircase(end).s.threshold=0.05;
    end
end

function [trials] = totalTrials()
%%

% Counts trials + estimates the threshold based on the last 500 trials

% get the files list
files = dir(fullfile(sprintf('~/data/afmap/%s/17*stim*.mat',mglGetSID)));

trials = 0;

for fi = 1:length(files)
    load(fullfile(sprintf('~/data/afmap/%s/%s',mglGetSID,files(fi).name)));
    
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
files = dir(fullfile(sprintf('~/data/afmap/%s/17*stim*.mat',mglGetSID)));

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
%     load(fullfile(sprintf('~/data/afmap_%s/%s/%s',rstimulus.condition,mglGetSID,files(fi).name)));
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
        sz = 4 * stimulus.gratingSizes(si); % make it twice as a big, so that the FWHM can be equal to the size
        % use total degs / num to compute size
        for phase = 1:2
            grating = stimulus.gratingContrasts(ci) * 255/2 * mglMakeGrating(sz,sz,4/stimulus.gratingSizes(si),0,(phase-1)*180) + 255/2;
            gauss = mglMakeGaussian(sz,sz,stimulus.gratingSizes(si)/(2*sqrt(2*log(2)))/2,stimulus.gratingSizes(si)/(2*sqrt(2*log(2)))/2);
            alphamask = repmat(grating,1,1,4);
            alphamask(:,:,4) = gauss*255;

            stimulus.grating(ci,si,phase)  = mglCreateTexture(alphamask); % high contrast
        end
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