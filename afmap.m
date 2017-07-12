function [ myscreen ] = afmap( varargin )
%ATTENTIONFIELDMAPPING 
%
%   Map the attention field in the scanner. This function works by having a
%   participant perform an asynchronous attention task at fixation or in a
%   quarterfield region. A pre-determined poisson process generates random
%   flashes of rotating gratings throughout the visual field at low or high
%   contrast.
%
%   The probe stimuli have three sizes 0.5x0.5, 1x1 or 2x2 deg, to
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
replay = 0;
attend = 0;
getArgs(varargin,{'scan=1','plots=0','noeye=0','debug=0','replay=0','attend=1'});
stimulus.scan = scan;
stimulus.plots = plots;
stimulus.noeye = noeye;
stimulus.debug = debug;
stimulus.replay = replay;
stimulus.attend = attend; % controls attention location: 0 = fixation, 1 = 5,5, 2 = 5,-5
clear localizer invisible scan noeye task test2 attend

if stimulus.scan
    warning('Not setup for scanning');
end

%% Replay mode
if any(replay>0)
    if isstr(replay)
        % a file was called for, load it
        loaded = load(replay);
        stimulus =  loaded.stimulus;
        stimulus.replayFile = strcat(replay(1:(strfind(replay,'.mat')-1)),'_replay.mat');
        stimulus.replay = true;
    else
        % do an entire subject
        disp('Replay mode initiated');
        folder = input('What folder would you like to replay? [Folder]: ');
        files = dir(fullfile(folder,'*.mat'));
        for fi = 1:length(files)
            if isempty(strfind(files(fi).name,'replay')) && isempty(strfind(files(fi).name,'original'))
                afmap(sprintf('replay=%s/%s',folder,files(fi).name));
            end
        end
        return;
    end
end

%% Open Old Stimfile
if ~stimulus.replay
    stimulus.counter = 1;

    if ~isempty(mglGetSID) && isdir(sprintf('~/data/afmap/%s',mglGetSID))
        % Directory exists, check for a stimefile
        files = dir(sprintf('~/data/afmap/%s/1*mat',mglGetSID));

        if length(files) >= 1
            fname = files(end).name;

            s = load(sprintf('~/data/afmap/%s/%s',mglGetSID,fname));
            % copy staircases and run numbers
            stimulus.counter = s.stimulus.counter + 1;
            stimulus.build = s.stimulus.build;
            stimulus.builds = s.stimulus.builds;
            stimulus.attention = s.stimulus.attention;
            stimulus.staircase = s.stimulus.staircase;
            stimulus.live.attend = mod(s.stimulus.live.attend+1,3);
            clear s;
            disp(sprintf('(afmap) Data file: %s loaded.',fname));
        end
    end
    disp(sprintf('(afmap) This is run #%i',stimulus.counter));
end

%% Stimulus parameters
if ~stimulus.replay
    stimulus.stimX = 25; % max ecc in any direction
    stimulus.stimY = 13;
    stimulus.stimR = 2; % deg between each stimulus

    if mod(stimulus.stimX,stimulus.stimR)==0 || mod(stimulus.stimY,stimulus.stimR)==0
        warning('Your stimulus size is not correctly setup');
    end

    stimulus.stimx = -stimulus.stimX:stimulus.stimR:stimulus.stimX;
    stimulus.stimy = -stimulus.stimY:stimulus.stimR:stimulus.stimY;
    
    if any(stimulus.stimx==0) || any(stimulus.stimy==0)
        warning('Your stimulus overlaps the fixation cross and the attention task!! You totally fd up!');
    end

    % how many times each probe location turns on per contrast and size
    % condition
    stimulus.probeOn = 2;
    % how long a probe stays up for (in TR, 4 = 2.0s)
    stimulus.probeUp = 4; % this must be EVEN!!
    % how long a probe is guaranteed to stay down (in TR, 12 = 6.0s)
    stimulus.probeDown = 12;

    % stimulus.live will hold what actually gets displayed on the screen
    stimulus.live.con = zeros(length(stimulus.stimx),length(stimulus.stimy));
    stimulus.live.sz = zeros(length(stimulus.stimx),length(stimulus.stimy));
    stimulus.live.ph = zeros(length(stimulus.stimx),length(stimulus.stimy));
    stimulus.live.theta = zeros(length(stimulus.stimx),length(stimulus.stimy));

    % gratingContrasts and gratingsizes control the possible sizes 
    stimulus.gratingContrasts = [0.1 1.0];
    stimulus.gratingSizes = [0.5 1 2];

    % when we are doing the attention task
    stimulus.live.attend = 0;

    % blank options
    stimulus.blanks.none.range = [0 0 0 0];
    stimulus.blanks.all.range = [-inf inf -inf inf];
    stimulus.blanks.NW.range = [-inf 0 0 inf];
    stimulus.blanks.NE.range = [0 inf 0 inf];
    stimulus.blanks.SE.range = [0 inf -inf 0];
    stimulus.blanks.SW.range = [-inf 0 -inf 0];
    stimulus.blanks.opts = {'NW','NE','SE','SW'};
end

%% Attention stimulus
if ~stimulus.replay && ~isfield(stimulus,'attention') 
    stimulus.attention = struct;
    stimulus.attention.attendX = [0 5 5];
    stimulus.attention.attendY = [0 5 -5];
    stimulus.attention.rotate = length(stimulus.attention.attendX);
    stimulus.attention.curAttend = 1;
end

%% Build stimulus
if ~stimulus.replay && ~isfield(stimulus,'build')
    stimulus.build = struct;
    
    stimulus.build.curBuild = 0; % will be incremented later
    
    stimulus.build.uniques = 4; % how many unique patterns to generate
    stimulus.build.rotate = 3; % how many patterns to rotate through (set to 4 or 5 for 2x repeat runs)
        
    stimulus.build.cycles = 6;
    stimulus.build.cycleLength = 120;
    stimulus.build.availableTRs = stimulus.build.cycles*stimulus.build.cycleLength; % how long the task should run for
    
    stimulus.build.conditions = length(stimulus.gratingContrasts)*length(stimulus.gratingSizes);
    stimulus.build.conditionsRep = stimulus.build.conditions * stimulus.probeOn; % total # of displays per location
    
    if stimulus.build.conditionsRep*(stimulus.probeUp+stimulus.probeDown) > (stimulus.build.availableTRs*2/3) % ~1/3 of the time the screen will be blank
        warning('The code has insufficient TRs available for the stimulus program you requested');
        keyboard
    end
    
    for bi = 1:stimulus.build.uniques
        build = struct; % initialize
        
        build.con = zeros(stimulus.build.availableTRs,length(stimulus.stimx),length(stimulus.stimy));
        build.sz = zeros(stimulus.build.availableTRs,length(stimulus.stimx),length(stimulus.stimy));
        build.ph = zeros(stimulus.build.availableTRs,length(stimulus.stimx),length(stimulus.stimy));
        build.theta = zeros(stimulus.build.availableTRs,length(stimulus.stimx),length(stimulus.stimy));
        
        % build the blackout positions for this run
        % each cycle goes ALL, {NW/NE/SE/SW}, NONE in 10 s increments (60 s
        % total)
        % get the rotations for each blackout and put them in place
        blackout = cell(stimulus.build.cycles,6);
        for c = 1:stimulus.build.cycles
            order = randperm(4);
            blackout{c,1} = 'none';
            for r = 1:4
                blackout{c,r+1} = stimulus.blanks.opts{order(r)};
            end
            blackout{c,6} = 'all';
        end
        build.blackout = blackout;
        
        % now compute the X/Y ranges that are subject to the blackout
        bminx = zeros(stimulus.build.cycles,6);
        bmaxx = zeros(stimulus.build.cycles,6);
        bminy = zeros(stimulus.build.cycles,6);
        bmaxy = zeros(stimulus.build.cycles,6);
        for c = 1:stimulus.build.cycles
            for d = 1:6
                bminx(c,d) = stimulus.blanks.(blackout{c,d}).range(1);
                bmaxx(c,d) = stimulus.blanks.(blackout{c,d}).range(2);
                bminy(c,d) = stimulus.blanks.(blackout{c,d}).range(3);
                bmaxy(c,d) = stimulus.blanks.(blackout{c,d}).range(4);
            end
        end
        
        disppercent(-1/length(stimulus.stimx));
        for x = 1:length(stimulus.stimx)
            for y = 1:length(stimulus.stimy)
                
                % determine which cycles will get which of the conditions
                % on which of their repeats
                %    1     2       3       4
                %  cycle  con     size    repeat
                redo = true;
                while redo
                    conditionTiming = zeros(stimulus.build.conditionsRep,4);
                    count = 1;
                    for i = 1:length(stimulus.gratingContrasts)
                        for j = 1:length(stimulus.gratingSizes)
                            for k = 1:stimulus.probeOn
                                conditionTiming(count,:) = [randi(stimulus.build.cycles) i j k];
                                count = count + 1;
                            end
                        end
                    end
                    redo = false;
                    for c = 1:stimulus.build.cycles
                        cycleCount = sum(conditionTiming(:,1)==c);
                        if (cycleCount==0) || (cycleCount >= (2 / stimulus.build.cycles * size(conditionTiming,1)))
                            % if we accidentally generated a set of cycles
                            % where there is a poor distribution of probes,
                            % we just repeat it and do it again. 
                            redo = true;
                        end
                    end
                end
                
                % pre-compute the blackout positions
                sx = stimulus.stimx(x);
                sy = stimulus.stimy(y);
                bpos = logical((sx>bminx).*(sx<bmaxx).*(sy>bminy).*(sy<bmaxy));
                for c = 1:stimulus.build.cycles
                    % get all the indexes
                    cStart = (c-1) * stimulus.build.cycleLength + 1;
                    cEnd = c * stimulus.build.cycleLength;
                    indexes = cStart:cEnd;
                    % get the blackout indexes for this cycle
                    for d = 1:6
                        if bpos(c,d)
                            % blackout
                            rmindexes = cStart-1 + (((d-1)*20+1):(d*20));
                            indexes = setdiff(indexes,rmindexes);
                        end
                    end
                    
                    % for each cycle, build a stimulus display timeseries,
                    % this means picking the actual timing of the various
                    % events within the time block
                    events = conditionTiming(conditionTiming(:,1)==c,:);
                    eventTimes = randsample(indexes,size(events,1));
                    attempts = 0;
                    while any(diff(eventTimes)<(stimulus.probeUp+stimulus.probeDown))
                        attempts = attempts + 1;
                        eventTimes = randsample(indexes,size(events,1));
%                         if attempts>10000
%                             disp(sprintf('(afmap) Attempted %i times to create a functional cycle and failed--stimulus properties are fd up',attempts));
%                             keyboard
%                         end
                    end
                    if attempts>100000
                        warning(sprintf('It took a whole lot of attempts: %i, to build that cycle',attempts));
                    end
                    for ei = 1:size(events,1)
                        eidxs = eventTimes(ei):(eventTimes(ei)+stimulus.probeUp-1);
                        build.con(eidxs,x,y) = events(ei,2);
                        build.sz(eidxs,x,y) = events(ei,3);
                        build.ph(eidxs,x,y) = repmat([1 2],1,length(eidxs)/2);
                        build.theta(eidxs,x,y) = rand*2*pi;
                    end
                    
                end
            end
            disppercent(x/length(stimulus.stimx));
        end
        disppercent(inf/length(stimulus.stimx));
        % test code
% %         keyboard
% %         figure
% %         colormap('gray');
% %         caxis([0 1]);
% % %         build.con = build.con>0;
% %         tb = build.con/max(build.con(:));
% %         for i = 1:720
% %             imagesc(squeeze(tb(i,:,:)));
% %             pause(.01);
% %         end
        % test code end
        disp(sprintf('(afmap) Pre-build of build %i has finished, saving.',bi));
        stimulus.builds{bi} = build;
    end
    disp(sprintf('(afmap) Pre-build complete. Created %i unique builds which will rotate every %i runs.',stimulus.build.uniques,stimulus.build.rotate));
end

%% Get the current build number

if ~stimulus.replay
    stimulus.build.curBuild = stimulus.build.curBuild + 1;
    if stimulus.build.curBuild > stimulus.build.rotate
        stimulus.build.curBuild = 1;
        stimulus.attention.curAttend = stimulus.attention.curAttend + 1;
        if stimulus.attention.curAttend > stimulus.attention.rotate
            stimulus.attention.curAttend = 1;
        end
    end
    stimulus.attention.curAttendX = stimulus.attention.attendX(stimulus.attention.curAttend);
    stimulus.attention.curAttendY = stimulus.attention.attendY(stimulus.attention.curAttend);
    disp(sprintf('(afmap) Build %i selected',stimulus.build.curBuild));
    disp(sprintf('(afmap) Attending X: %i Y: %i selected',stimulus.attention.curAttendX,stimulus.attention.curAttendY));
end

%% Load the current build
if ~stimulus.replay
    % stimulus.grid is used to track the grid for this run -- we will
    % update this every time the screen changes and save the time it
    % occurred
    stimulus.grid.buildNumber = stimulus.build.curBuild;
    stimulus.grid.t = zeros(1,stimulus.build.availableTRs);
    stimulus.grid.con = stimulus.builds{stimulus.build.curBuild}.con;
    stimulus.grid.sz = stimulus.builds{stimulus.build.curBuild}.sz;
    stimulus.grid.ph = stimulus.builds{stimulus.build.curBuild}.ph;
    stimulus.grid.theta = stimulus.builds{stimulus.build.curBuild}.theta;
    disp(sprintf('(afmap) Build %i loaded from pre-build',stimulus.build.curBuild));
end

%% Setup attention
if ~stimulus.replay
    stimulus.attendX = [0 5 5];
    stimulus.attendY = [0 5 -5];
    stimulus.live.aX = stimulus.attendX(stimulus.live.attend+1);
    stimulus.live.aY = stimulus.attendY(stimulus.live.attend+1);
end

%% Setup Screen
if stimulus.replay
    myscreen = initScreen('replayScreen');
else
    myscreen = initScreen('VPixx');
end

% set background to grey
myscreen.background = 0.5;

%% Staircase
if ~stimulus.replay
    if ~isfield(stimulus,'staircase')
        disp('(afmap) WARNING: New staircase');
        initStair();
    else
        resetStair();
    end
end

%% Plot and return
if stimulus.plots==2
    dispInfo(stimulus);
    return
end

%% Initialize Stimulus

if ~stimulus.replay
    myscreen.stimulusNames{1} = 'stimulus';

    if ~isfield(stimulus.live,'grating')
        localInitStimulus();
    end

    stimulus.responseKeys = [1 2]; % left right
end

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
% task waits for fixation on first segment
if stimulus.replay
    task{1}{1}.waitForBacktick = 0;
    task{1}{1}.seglen = repmat(0.050,1,120);
else
    task{1}{1}.waitForBacktick = 1;
    task{1}{1}.seglen = repmat(0.500,1,120);
end

stimulus.seg.stim = 1;

task{1}{1}.getResponse = 0;

if stimulus.replay
    task{1}{1}.numTrials = size(stimulus.grid.con,1);
else
    task{1}{1}.numTrials = stimulus.build.cycles;
end

task{1}{1}.random = 0;

if ~stimulus.replay && stimulus.scan
    task{1}{1}.synchToVol = zeros(1,length(task{1}{1}.seglen));
    task{1}{1}.synchToVol(end) = 1;
end

task{1}{1}.randVars.calculated.probesOn = nan;

%% Setup Attention Task

stimulus.curTrial = 0;

if ~stimulus.replay
    global fixStimulus %#ok<TLEV>

    fixStimulus.diskSize = 0.75;
    fixStimulus.fixWidth = 0.75;
    fixStimulus.fixLineWidth = 1;
    fixStimulus.stimTime = 0.35;
    fixStimulus.interTime = 1.4;
    fixStimulus.stairUsePest = 1;
    fixStimulus.pos = [stimulus.attention.curAttendX stimulus.attention.curAttendY];
    [task{2}, myscreen] = gruFixStairInitTask(myscreen);
    
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
end

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
if ~stimulus.replay
    myscreen = eyeCalibDisp(myscreen);

    % let the user know
    disp(sprintf('(afmap) Starting run number: %i.',stimulus.counter));
end

%% Main Task Loop

% setGammaTable(1);
if stimulus.replay
    mglClearScreen(0);
    mglFlush;
    mglClearScreen(0);
else
    mglClearScreen(0.5); %mglFixationCross(1,1,stimulus.colors.white);
    if stimulus.attention.curAttendX>0 || stimulus.attention.curAttendY > 0
        mglFixationCross(1,1,stimulus.colors.black);
    end
    mglFlush
    mglClearScreen(0.5); %mglFixationCross(1,1,stimulus.colors.white);
    if stimulus.attention.curAttendX>0 || stimulus.attention.curAttendY > 0
        mglFixationCross(1,1,stimulus.colors.black);
    end
end

phaseNum = 1;
% Again, only one phase.
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    % update the task
    [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);
    if length(task)>1
        [task{2}, myscreen, phaseNum] = updateTask(task{2},myscreen,phaseNum);
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

% save stimulus
if stimulus.replay
    s = load(replay);

    screenWidth = myscreen.screenWidth; screenHeight = myscreen.screenHeight;
    imageWidth = myscreen.imageWidth; imageHeight = myscreen.imageHeight;
    [pRFstim.x, pRFstim.y] = ndgrid(-imageWidth/2:imageWidth/(screenWidth-1):imageWidth/2, -imageHeight/2:imageHeight/(screenHeight-1):imageHeight/2);

    pRFstim.t = 1:size(stimulus.frames,3);
    
    if size(stimulus.frames,3)~=720
        warning('Number of frames does not match the expected length (840)--padding end with blank screen');
        input('Press [enter] to confirm: ');
        stimulus.frames(:,:,end:720) = 0;
    end
    
    pRFstim.im = stimulus.frames;
    
    s.pRFStimImage = pRFstim;
    save(replay,'-struct','s');
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

if ~stimulus.replay && stimulus.plots
    disp('(afmap) Displaying plots');
    dispInfo(stimulus);
end

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback1(task,myscreen)
%%
global stimulus

stimulus.curTrial = stimulus.curTrial + 1;

stimulus.live.con = squeeze(stimulus.grid.con(stimulus.curTrial,:,:));
stimulus.live.sz = squeeze(stimulus.grid.sz(stimulus.curTrial,:,:));
stimulus.live.ph = squeeze(stimulus.grid.ph(stimulus.curTrial,:,:));
stimulus.live.theta = squeeze(stimulus.grid.theta(stimulus.curTrial,:,:));
stimulus.grid.t(stimulus.curTrial) = mglGetSecs;

if stimulus.replay
    mglClearScreen(0);
    myscreen.flushMode = 1;
else
    mglClearScreen();
    myscreen.flushMode = 1;
end

% draw gratings for probe task

for xi = 1:length(stimulus.stimx)
    for yi = 1:length(stimulus.stimy)
        if stimulus.live.con(xi,yi)>0
            x = stimulus.stimx(xi);
            y = stimulus.stimy(yi);
            con = stimulus.live.con(xi,yi);
            sz = stimulus.live.sz(xi,yi);
            ph = stimulus.live.ph(xi,yi);
            theta = stimulus.live.theta(xi,yi);

            if stimulus.replay
                % just draw a circle
                % /2 because the FWHM defines a diameter of 1/2/3 degree
                mglBltTexture(stimulus.gaussian(con,sz,ph),[x y],0,0,0);
    %                 mglFillOval(x,y,repmat(stimulus.gratingSizes(sz)/(2*sqrt(2*log(2)))*2,1,2),stimulus.gratingContrasts(con)*[1 1 1]);
            else
                mglBltTexture(stimulus.grating(con,sz,ph),[x y],0,0,theta*180/pi);
            end
        end
    end
end

if stimulus.replay
    mglFlush % the screen will blank after the frame, but whatever
    frame = mglFrameGrab;
    if ~isfield(stimulus,'frames')
        stimulus.frames = zeros(myscreen.screenWidth,myscreen.screenHeight,stimulus.gridCount);
    end
    stimulus.frames(:,:,stimulus.curTrial) = frame(:,:,1);
end

% disp(sprintf('(afmap) Starting trial %01.0f',stimulus.curTrial));
disppercent(stimulus.curTrial/120,'(afmap) Running: ');

function [task, myscreen] = startTrialCallback1(task,myscreen)
global stimulus

disp(sprintf('(afmap) Starting cycle %01.0f',(stimulus.curTrial/120)+1));
disppercent(-1/120,'(afmap) Running: ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback1(task, myscreen)
%%
global stimulus

if stimulus.attention.curAttendX>0 || stimulus.attention.curAttendY > 0
    mglFixationCross(1,1,stimulus.colors.black);
end


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
            
            % make the grating
            stimulus.grating(ci,si,phase)  = mglCreateTexture(alphamask); % high contrast
            % make a gaussian (for when we display, make sure to use the
            % actual contrast for this setting or we'll fuck up later)
            gData = stimulus.gratingContrasts(ci)*255*ones(size(grating,1),size(grating,2),4);
            gData(:,:,4) = gauss*255;
            stimulus.gaussian(ci,si,phase) = mglCreateTexture(gData);
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