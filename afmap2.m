function [ myscreen ] = afmap2( varargin )
% ***scan info: 736 TRs (6 minute * 120 + 16)***
% *** set 'noeye=1' to turn of the eye tracker***
%
%ATTENTIONFIELDMAPPING (afmap2)
%
%   Map the attention field in the scanner. This function works by having a
%   participant perform an asynchronous attention task at fixation or in a
%   quarterfield region. A random process generates flashes of rotated 
%   gratings throughout the visual field at low or high contrast.
%   
%   The attention task involves performing the standard luminance decrement
%   task at an off-fixation location. The task is sped up to be more
%   continuous than it usually is and is staircased.
%
%   If you accidentally crash a run or need to stop the scanner you can
%   specify a run by setting the flag, 'run=#'.
%
%   During scanning use 'noeye=1'

%%

global stimulus fixStimulus

stimulus = struct;
fixStimulus = struct;

%% Initialize Variables

% add arguments later
scan = 0;
plots = 0;
noeye = 0;
debug = 0;
replay = 0;
attend = 0; run = 0; build = 0; eyewindow=0; mouse=0;
getArgs(varargin,{'scan=1','plots=0','noeye=0','eyewindow=1.5','debug=0','replay=0','attend=1','run=0','build=0','mouse=0'});
stimulus.scan = scan;
stimulus.plots = plots;
stimulus.noeye = noeye;
stimulus.mousedebug = mouse;
stimulus.eyewindow = eyewindow;
stimulus.debug = debug;
stimulus.replay = replay;
stimulus.overrideRun = run;
stimulus.attend = attend; % controls whether attention mode runs
stimulus.buildOverride = build;
if ~stimulus.attend
    warning('*****ATTENTION MODE IS DISABLED*****');
end
clear localizer invisible scan noeye task test2 attend build eyewindow mouse

%% Replay mode
if any(replay>0)
    if ischar(replay)
        % a file was called for, load it
        loaded = load(replay);
        stimulus = loaded.stimulus;
        % check that this is actually an afmap file
        if isempty(strfind(loaded.task{1}{1}.taskFilename,'afmap'))
            disp(sprintf('File %s is not an afmap run.',replay));
            return
        end
        % get the task parameters so that you can sync the replay correctly
        disp('GET TASK PARAMETERS');
        e = getTaskParameters(loaded.myscreen,loaded.task);
        e1 = e{1};
        % pull out the trial volumes
        stimulus.tVolumes = e{1}.trialVolume;
        disp('Stimulus volumes were found at:');
        disp(stimulus.tVolumes);
        
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
    stimulus.curRun = 1;

    if ~isempty(mglGetSID) && isdir(sprintf('~/data/afmap2/%s',mglGetSID))
        % Directory exists, check for a stimefile
        files = dir(sprintf('~/data/afmap2/%s/1*mat',mglGetSID));

        if length(files) >= 1
            fname = files(end).name;
            
            s = load(sprintf('~/data/afmap2/%s/%s',mglGetSID,fname));
            % copy staircases and run numbers
            stimulus.counter = s.stimulus.counter + 1;
            stimulus.curRun = s.stimulus.curRun + 1;
            stimulus.build = s.stimulus.build;
            stimulus.builds = s.stimulus.builds;
            stimulus.attention = s.stimulus.attention;
            stimulus.staircases = s.stimulus.staircases;
            stimulus.order = s.stimulus.order;
            stimulus.live.attend = mod(s.stimulus.live.attend+1,3);
            if s.stimulus.attend ~= stimulus.attend
                error('(afmap2) Cannot continue: stimfile parameters were generated with a different attention mode than you requested. You need to save the existing stimfiles into a backup folder');
            end
            clear s;
            disp(sprintf('(afmap2) Data file: %s loaded.',fname));
        else
            warning('(afmap2) Unable to load previous data files. If this is *not* the first run there is something wrong.');
        end
    end
end

%% Override run
if ~stimulus.replay && stimulus.overrideRun>0
    stimulus.curRun = stimulus.overrideRun;
end

%% Display run info
if ~stimulus.replay
    disp('*************************');
    disp(sprintf('(afmap2) This is scan #%i',stimulus.counter));
    disp(sprintf('(afmap2) This is run #%i',stimulus.curRun));
    disp('*************************');
end

%% Stimulus parameters
if ~stimulus.replay
    % we'll interpolate between these along an arc
    stimulus.thetaMin = 0;
    stimulus.thetaMax = 2*pi;
    stimulus.minEcc = 1;
    
    stimulus.maxOnScreen = 36;

    % how 
    stimulus.probeOn = 2;
    % how long a probe stays up for (in TR, 4 = 2.0s)
    stimulus.probeUp = 4; % this must be EVEN!!
    % how long a probe is guaranteed to stay down (in TR, 12 = 6.0s)
    stimulus.probeDown = 12;

    % stimulus.live will hold what actually gets displayed on the screen
% %     stimulus.live.con = zeros(length(stimulus.stimx),length(stimulus.stimy));
% %     stimulus.live.sz = zeros(length(stimulus.stimx),length(stimulus.stimy));
% %     stimulus.live.ph = zeros(length(stimulus.stimx),length(stimulus.stimy));
% %     stimulus.live.theta = zeros(length(stimulus.stimx),length(stimulus.stimy));

    % gratingContrasts and gratingsizes control the possible sizes 
    stimulus.gratingContrasts = [0.1 1.0];
    % design eccs
    stimulus.designEccs = logspace(0,log10(6),8);
    stimulus.drawEccs = logspace(0,log10(28),8);
    % the ratios are approximately the sigma / ecc ratio for V1, V4, and
    % higher regions (MT/LO/VO/TO)
    stimulus.gratingRatios = [.15 .27 0.5];

    % when we are doing the attention task
    stimulus.live.attend = 0;

    % blank options
    stimulus.blanks.none.range = [2*pi 0];
    stimulus.blanks.all.range = [0 2*pi];
    stimulus.blanks.NW.range = [pi/2 pi];
    stimulus.blanks.NE.range = [0 pi/2];
    stimulus.blanks.SE.range = [3/2*pi 2*pi];
    stimulus.blanks.SW.range = [pi 3/2*pi];
    stimulus.blanks.opts = {'NW','NE','SE','SW'};
end

%% Attention stimulus
if ~stimulus.replay && ~isfield(stimulus,'attention') 
    stimulus.attention = struct;
    stimulus.attention.attendX = [0 5 0];
    stimulus.attention.attendY = [0 -5 0];

    if stimulus.attend
        stimulus.attention.rotate = length(stimulus.attention.attendX);
    else
        stimulus.attention.rotate = 1;
    end
    stimulus.attention.curAttend = 1;
end

%% Build stimulus
if ~stimulus.replay && ~isfield(stimulus,'build')
    stimulus.build = struct;
    
    stimulus.build.curBuild = 0; % will be incremented later
    
    stimulus.build.uniques = 3; % how many unique patterns to generate
    stimulus.build.rotate = 3; % how many patterns to rotate through (set to 4 or 5 for 2x repeat runs)
        
    stimulus.build.cycles = 6;
    stimulus.build.cycleLength = 120;
    stimulus.build.availableTRs = stimulus.build.cycles*stimulus.build.cycleLength; % how long the task should run for
    
    stimulus.build.conditions = length(stimulus.gratingContrasts)*length(stimulus.gratingRatios);
    stimulus.build.conditionsRep = stimulus.build.conditions * stimulus.probeOn; % total # of displays per location
    
    if stimulus.build.conditionsRep*(stimulus.probeUp+stimulus.probeDown) > (stimulus.build.availableTRs*2/3) % ~1/3 of the time the screen will be blank
        warning('The code has insufficient TRs available for the stimulus program you requested');
        keyboard
    end
    
    disp(sprintf('Building %i unique probe sequences. Each sequence consists of %i cycles of %i TRs each.',stimulus.build.uniques,stimulus.build.cycles,stimulus.build.cycleLength));
    
    for bi = 1:stimulus.build.uniques
        build = struct; % initialize
        
        build.t = zeros(stimulus.build.availableTRs,stimulus.maxOnScreen);
        build.ecc = zeros(stimulus.build.availableTRs,stimulus.maxOnScreen);
        build.con = zeros(stimulus.build.availableTRs,stimulus.maxOnScreen);
        build.sz = zeros(stimulus.build.availableTRs,stimulus.maxOnScreen);
        build.ph = zeros(stimulus.build.availableTRs,stimulus.maxOnScreen);
        build.theta = zeros(stimulus.build.availableTRs,stimulus.maxOnScreen);
        
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
        bmint = zeros(stimulus.build.cycles,6);
        bmaxt = zeros(stimulus.build.cycles,6);
        for c = 1:stimulus.build.cycles
            for d = 1:6
                bmint(c,d) = stimulus.blanks.(blackout{c,d}).range(1);
                bmaxt(c,d) = stimulus.blanks.(blackout{c,d}).range(2);
            end
        end
        
        mult = 3;
        onScreenNum = [36 27 27 27 27 0]*mult;
        
        % polar angles where eccentricity was measured
        mesAngles =  [-90 -60 -45 -30 -15 0 15 30 45 60 75 90];
        % eccentricity measured at the above angles
        mesEcc = ([15 17 20 28 38 35 30 26 24 22 21 21]-1)*1;

        % linear interpolate (in radial coordinates) to make smoother
        angles = -90:90;
        ecc = interp1(mesAngles,mesEcc,angles,'linear');
        angles = angles * pi/180;
        
        angles = [angles pi-angles];
        ecc = [ecc ecc];
        
        angles = mod(angles,2*pi);

%         figure; hold on
%         for i = 1:length(angles)
%             x = ecc(i) * cos(angles(i));
%             y = ecc(i) * sin(angles(i));
%             plot(x,y,'*');
%         end
        
        disppercent(-1/stimulus.build.availableTRs);
        
        availableIdxs = ones(stimulus.build.availableTRs,ceil(stimulus.maxOnScreen*mult/3));
        
        for cycle = 1:stimulus.build.cycles
            for innerCycle = 1:6
                mint = bmint(cycle,innerCycle);
                maxt = bmaxt(cycle,innerCycle);
                
                % get the TRs for the current (inner) cycle
                TRs = (cycle-1)*120 + (((innerCycle-1)*20+1):innerCycle*20);
                
                % get how many onScreen we should generate (random, 50-100%
                % of max)
                maxOn = onScreenNum(innerCycle);
                turnOn = round(rand*maxOn*.5 + maxOn*.5);
                
                % get a random index for each one
                idxs = randsample(TRs,turnOn,true);
                
                % get parameters, t (angle), ecc, size, etc
                t = rand(size(idxs))*2*pi;
                while any((t>mint).*(t<maxt))
                    ridxs = logical((t>mint).*(t<maxt));
                    t(ridxs) = rand(1,sum(ridxs))*2*pi;
                end
                
                % for ecc, flip across and then round to get maxEcc values
                t2 = t;
%                 t2 = mod(t2,2*pi);
                
                eccMax = zeros(size(t2));
                
                for ti = 1:length(t2)
                    [~,mini] = min(abs(angles-t2(ti)));
                    eccMax(ti) = ecc(mini);
                end
                
                necc = stimulus.minEcc + rand(size(eccMax)).*(eccMax-stimulus.minEcc);
                
                con = randi(2,1,turnOn);
                sz = randi(3,1,turnOn);
                initph = rand(1,turnOn)>.5;
                theta = rand(1,turnOn)*2*pi;
                
                % sort indexes
                [sidxs,i] = sort(idxs);
                
                for oi = 1:turnOn
                    cStartidx = sidxs(oi);
                    ct = t(i(oi));
                    cecc = necc(i(oi));
                    ccon = con(i(oi));
                    csz = sz(i(oi));
                    cinitph = initph(i(oi));
                    cph = [cinitph ~cinitph cinitph ~cinitph]+1;
                    ctheta = theta(i(oi));
                    
                    % first obtain the indexes that this stimulus will be
                    % up for
                    lidxs = cStartidx:(cStartidx+stimulus.probeUp-1);
                    
                    % now check where we can index this stimulus in the
                    % index tracker
                    storeIdx = find(availableIdxs(cStartidx,:),1);
                    % block these indexes in the tracker
                    availableIdxs(lidxs,storeIdx) = 0;
                    
                    % now block off the appropriate indexes in all the
                    % build section
                    build.t(lidxs,storeIdx) = ct; % theta does not change
                    build.ecc(lidxs,storeIdx) = cecc;
                    build.con(lidxs,storeIdx) = ccon; % contrast does not change
                    build.sz(lidxs,storeIdx) = csz;
                    build.ph(lidxs,storeIdx) = cph;
                    build.theta(lidxs,storeIdx) = ctheta;
                end
            end
        end
        disppercent(inf/stimulus.build.availableTRs);
        % test code
        
        %%
        % test code end
        disp(sprintf('(afmap2) Pre-build of build %i has finished (will be saved with stimfile).',bi));
        stimulus.builds{bi} = build;
    end
    disp(sprintf('(afmap2) Pre-build complete. Created %i unique builds which will rotate every %i runs.',stimulus.build.uniques,stimulus.build.rotate));
end

%% Build the order
if ~stimulus.replay && ~isfield(stimulus,'order')
    stimulus.order = struct;
    
    stimulus.order.BaseAtt = zeros(1,stimulus.build.rotate*stimulus.attention.rotate);
    for ai = 1:stimulus.attention.rotate
        stimulus.order.BaseAtt((ai-1)*stimulus.build.rotate+1:ai*stimulus.build.rotate) = ai*ones(1,stimulus.build.rotate);
    end
    stimulus.order.BaseBui = repmat(1:stimulus.build.rotate,1,stimulus.attention.rotate);
    stimulus.order.n = length(stimulus.order.BaseAtt);
    stimulus.order.curOrder = [];
    stimulus.order.doneMat = zeros(stimulus.attention.rotate,stimulus.build.rotate);
end

%% Staircase
if ~stimulus.replay
    if ~isfield(stimulus,'staircases')
        disp('(afmap2) WARNING: New staircase');
        initStair();
    else
        resetStair();
    end
end

%% Get the current build number
if ~stimulus.replay
    % need to set: stimulus.build.curBuild and stimulus.attention.curAttend
    
    if stimulus.curRun > length(stimulus.order.curOrder)
        nOrder = randperm(stimulus.order.n);
        stimulus.order.curOrder = [stimulus.order.curOrder nOrder];
    end
    
    orderIdx = stimulus.order.curOrder(stimulus.curRun);
    stimulus.build.curBuild = stimulus.order.BaseBui(orderIdx);
    % override build
    stop = 1;
    
    if stimulus.buildOverride
        stimulus.build.curBuild = stimulus.buildOverride;
    end
    
    stimulus.attention.curAttend = stimulus.order.BaseAtt(orderIdx);
    
    stimulus.attention.curAttendX = stimulus.attention.attendX(stimulus.attention.curAttend);
    stimulus.attention.curAttendY = stimulus.attention.attendY(stimulus.attention.curAttend);
    
    stimulus.staircase = stimulus.staircases{stimulus.attention.curAttend};
end

%% Display completion information
if ~stimulus.replay

    strs = {};
    initStrs = {'Build\t','#\t'};
    for bi = 1:stimulus.build.rotate
        if bi==1
            strs{end+1} = sprintf('\t');
        elseif bi>(length(initStrs)+1)
            strs{end+1} = sprintf('\t');
        else
            strs{end+1} = sprintf(initStrs{bi-1});
        end
        for ai = 1:stimulus.attention.rotate
            strs{end} = sprintf('%s%i\t',strs{end},stimulus.order.doneMat(ai,bi));
        end
    end
    if stimulus.attend
        disp('******************************');
        disp(sprintf('\tAttention condition'));
    %     disp(aconds);
        for bi = 1:stimulus.build.rotate
            disp(sprintf('%s',strs{bi}));
        end
        disp('******************************');
    end
end

%% Display build info
if ~stimulus.replay
    disp('******************************');
    disp(sprintf('(afmap2) Build %i selected',stimulus.build.curBuild));
    disp(sprintf('(afmap2) Attending X: %i Y: %i selected',stimulus.attention.curAttendX,stimulus.attention.curAttendY));
    disp('******************************');
end

%% If tVolumes are not % 120 = 1, then we need to adjust the build
if stimulus.replay
    m120 = mod(stimulus.tVolumes,120);
    
    if any(m120>1)
        disp('WARNING: Volume timing failed. Build output must be adjusted to contain additional blanks.');
        
        % get the indexes for each build (i.e. 1:120, 121:240, etc)
        def = [1:120:720];
        bIndexes = zeros(6,120); nIndexes = zeros(6,120);
        for bi = 1:length(def)
            bIndexes(bi,:) = def(bi):(def(bi)+119);
            nIndexes(bi,:) = stimulus.tVolumes(bi):(stimulus.tVolumes(bi)+119);
        end
        
        % save orig
        stimulus.grid.preTRfix = stimulus.grid;
        
        nMax = max(nIndexes(:));
        
        orig_sz = size(stimulus.grid.con);
        
        ncon = zeros(nMax,orig_sz(2),orig_sz(3));
        nsz = ncon;
        nph = ncon;
        ntheta = ncon;
        for bi = 1:size(bIndexes,1)
            disp(sprintf('Build %i needs to be offset by %i from %i to %i',bi,m120(bi),def(bi),stimulus.tVolumes(bi)));
            
            nt(nIndexes(bi,:),:) = stimulus.grid.t(bIndexes(bi,:),:);
            necc(nIndexes(bi,:),:) = stimulus.grid.ecc(bIndexes(bi,:),:);
            ncon(nIndexes(bi,:),:) = stimulus.grid.con(bIndexes(bi,:),:);
            nsz(nIndexes(bi,:),:) = stimulus.grid.sz(bIndexes(bi,:),:);
            nph(nIndexes(bi,:),:) = stimulus.grid.ph(bIndexes(bi,:),:);
            ntheta(nIndexes(bi,:),:) = stimulus.grid.theta(bIndexes(bi,:),:);
        end
        
        % replace
        stimulus.grid.t = nt(1:orig_sz(1),:);
        stimulus.grid.ecc = necc(1:orig_sz(1),:);
        stimulus.grid.con = ncon(1:orig_sz(1),:);
        stimulus.grid.sz = nsz(1:orig_sz(1),:);
        stimulus.grid.ph = nph(1:orig_sz(1),:);
        stimulus.grid.theta = ntheta(1:orig_sz(1),:);
    end
end

%% Load the current build
if ~stimulus.replay
    % stimulus.grid is used to track the grid for this run -- we will
    % update this every time the screen changes and save the time it
    % occurred
    stimulus.grid.buildNumber = stimulus.build.curBuild;
    stimulus.grid.time = zeros(1,stimulus.build.availableTRs);
    stimulus.grid.t = stimulus.builds{stimulus.build.curBuild}.t;
    stimulus.grid.ecc = stimulus.builds{stimulus.build.curBuild}.ecc;
    stimulus.grid.con = stimulus.builds{stimulus.build.curBuild}.con;
    stimulus.grid.sz = stimulus.builds{stimulus.build.curBuild}.sz;
    stimulus.grid.ph = stimulus.builds{stimulus.build.curBuild}.ph;
    stimulus.grid.theta = stimulus.builds{stimulus.build.curBuild}.theta;
    disp(sprintf('(afmap2) Build %i loaded from pre-build',stimulus.build.curBuild));
end

%% Setup attention
if ~stimulus.replay
    stimulus.attendX = stimulus.attention.attendX;
    stimulus.attendY = stimulus.attention.attendY;
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
else
    localInitStimulus();
end

%% Colors
stimulus.colors.white = [1 1 1]; stimulus.colors.red = [1 0 0];
stimulus.colors.green = [0 1 0]; stimulus.colors.black = [0 0 0];

%% Setup Probe Task

task{1}{1} = struct;
% task waits for fixation on first segment
if stimulus.replay
    task{1}{1}.waitForBacktick = 0;
    task{1}{1}.seglen = repmat(0.050,1,120);
else
    task{1}{1}.waitForBacktick = 1;
    task{1}{1}.seglen = repmat(0.500,1,120);
    if stimulus.scan
        task{1}{1}.seglen(end) = 0.050; % make the last segment short so that it will synchtovol and hopefully align the runs correctly
    end
end

stimulus.seg.stim = 1;

task{1}{1}.getResponse = 0;

task{1}{1}.numTrials = stimulus.build.cycles;

task{1}{1}.random = 0;

if ~stimulus.replay && stimulus.scan
    task{1}{1}.synchToVol = zeros(1,length(task{1}{1}.seglen));
    task{1}{1}.synchToVol(end) = 1;
end

task{1}{1}.randVars.calculated.probesOn = nan;

%% Setup Attention Task

stimulus.curStep = 0;

if ~stimulus.replay
    fixStimulus.diskSize = 1;
    fixStimulus.fixWidth = 1.25;
    fixStimulus.fixLineWidth = 4;
    fixStimulus.stairStepSize = 0.02;
    fixStimulus.stimTime = 0.25;
    fixStimulus.interTime = 0.35;
    fixStimulus.stairUsePest = 1;
    fixStimulus.responseTime = 1;
    fixStimulus.staircase = stimulus.staircase;
    fixStimulus.pos = [stimulus.attention.curAttendX stimulus.attention.curAttendY];
    [task{2}, myscreen] = gruFixStairInitTask_afmap(myscreen);
end

%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback1,@screenUpdateCallback1,[],@startTrialCallback1,[],[]);
end

%% EYE CALIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~stimulus.replay
    myscreen = eyeCalibDisp(myscreen);

    % let the user know
    disp(sprintf('(afmap2) Starting run number: %i.',stimulus.counter));
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
        mglFixationCross(1,3,stimulus.colors.black);
    end
    mglFlush
    mglClearScreen(0.5); %mglFixationCross(1,1,stimulus.colors.white);
    if stimulus.attention.curAttendX>0 || stimulus.attention.curAttendY > 0
        mglFixationCross(1,3,stimulus.colors.black);
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
        warning('Number of frames does not match the expected length (720)--padding end with blank screen');
        input('Press [enter] to confirm: ');
        stimulus.frames(:,:,end:720) = 0;
    end
    
    pRFstim.im = stimulus.frames;
    
    s.pRFStimImage = pRFstim;
    save(replay,'-struct','s');
end

% save info
if ~stimulus.replay
    % save which run we just finished
    stimulus.order.doneMat(stimulus.attention.curAttend,stimulus.build.curBuild) = stimulus.order.doneMat(stimulus.attention.curAttend,stimulus.build.curBuild)+1;
    % save staircase
    stimulus.staircase = fixStimulus.staircase;
    stimulus.staircases{stimulus.attention.curAttend} = stimulus.staircase;
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

if ~stimulus.replay && stimulus.plots
    disp('(afmap2) Displaying plots');
    dispInfo(stimulus);
end

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%


function [task, myscreen] = startTrialCallback1(task,myscreen)
global stimulus

disp(sprintf('(afmap2) Starting cycle %01.0f',(stimulus.curStep/120)+1));
stimulus.live.dead = 0;
stimulus.live.eyeCount=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback1(task,myscreen)
%%
global stimulus fixStimulus

stimulus.curStep = stimulus.curStep + 1;

% check online eye
if (stimulus.live.dead==2) && all(fixStimulus.thisColor==[0 1 1])
    % we can reset the colors now
    stimulus.live.dead= 0;
    stimulus.live.eyeCount = 0;
elseif (stimulus.live.dead==1) && all(fixStimulus.thisColor==[1 1 0])
    % set the timer
    stimulus.live.dead = 2;
end

stimulus.live.t = squeeze(stimulus.grid.t(stimulus.curStep,:));
stimulus.live.ecc = squeeze(stimulus.grid.ecc(stimulus.curStep,:));
stimulus.live.con = squeeze(stimulus.grid.con(stimulus.curStep,:));
stimulus.live.sz = squeeze(stimulus.grid.sz(stimulus.curStep,:));
stimulus.live.ph = squeeze(stimulus.grid.ph(stimulus.curStep,:));
stimulus.live.theta = squeeze(stimulus.grid.theta(stimulus.curStep,:));
stimulus.grid.t(stimulus.curStep) = mglGetSecs;

if stimulus.replay
    mglClearScreen(0);
    drawGratings();
    mglFlush
    mglClearScreen(0);
    drawGratings();
    myscreen.flushMode = 1;
else
    mglClearScreen();
    drawGratings();
    if stimulus.attention.curAttendX>0 || stimulus.attention.curAttendY > 0
        mglFixationCross(1,3,stimulus.colors.black);
    end
    drawFix(myscreen);
    mglFlush;
    mglClearScreen();
    drawGratings();
    if stimulus.attention.curAttendX>0 || stimulus.attention.curAttendY > 0
        mglFixationCross(1,3,stimulus.colors.black);
    end
    drawFix(myscreen);
end

% draw gratings for probe task

if stimulus.replay
    mglFlush % the screen will blank after the frame, but whatever
    frame = mglFrameGrab;
    if ~isfield(stimulus,'frames')
        stimulus.frames = zeros(myscreen.screenWidth,myscreen.screenHeight,stimulus.build.availableTRs);
    end
    stimulus.frames(:,:,stimulus.curStep) = frame(:,:,1);
end

function drawFix(myscreen)

global stimulus fixStimulus;

if stimulus.live.dead
    mglGluDisk(fixStimulus.pos(1),fixStimulus.pos(2),fixStimulus.diskSize*[1 1],stimulus.colors.red,60);
else
    if fixStimulus.trainingMode,mglClearScreen;end

    if ~isempty(fixStimulus.displayText)
      mglBltTexture(fixStimulus.displayText,fixStimulus.displayTextLoc);
    end
    mglGluDisk(fixStimulus.pos(1),fixStimulus.pos(2),fixStimulus.diskSize*[1 1],myscreen.background,60);

    mglFixationCross(fixStimulus.fixWidth,fixStimulus.fixLineWidth,fixStimulus.thisColor,fixStimulus.pos);
end

% draw an annulus around the cross that stays up permanently (helps to
% avoid losing the cross due to adaptation)
for i = 0:45:(360-1)
    mglGluPartialDisk(fixStimulus.pos(1),fixStimulus.pos(2),fixStimulus.diskSize,fixStimulus.diskSize*1.1,i,22.5,stimulus.colors.white);
end

function drawGratings

global stimulus

live = find(stimulus.live.con>0);

for si = 1:length(live)
    cidx = live(si);
    ct = stimulus.live.t(cidx);
    cecc = stimulus.live.ecc(cidx);
    x = cecc*cos(ct);
    y = cecc*sin(ct);
    con = stimulus.live.con(cidx);
    sz = stimulus.live.sz(cidx);
    ph = stimulus.live.ph(cidx);
    theta = stimulus.live.theta(cidx);
    
    [~,eccIdx] = min(abs(cecc-stimulus.drawEccs));
    
    if stimulus.replay
        % just draw a circle
        % /2 because the FWHM defines a diameter of 1/2/3 degree
        mglBltTexture(stimulus.gaussian(con,eccIdx,sz,ph),[x y],0,0,0);
    %                 mglFillOval(x,y,repmat(stimulus.gratingSizes(sz)/(2*sqrt(2*log(2)))*2,1,2),stimulus.gratingContrasts(con)*[1 1 1]);
    else
        % organized: contrast, ecc, size, ph
        mglBltTexture(stimulus.grating(con,eccIdx,sz,ph),[x y],0,0,theta*180/pi);
    end
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback1(task, myscreen)
%%
global stimulus

drawFix(myscreen);

% check eye pos
if (~stimulus.noeye) && (stimulus.eyewindow>0)
    
    % mouse version for testing with no eyetracker
    if stimulus.mousedebug
        mInfo = mglGetMouse(myscreen.screenNumber);
        degx = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
        degy = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;

        pos = [degx, degy];
    else
        [pos,~] = mglEyelinkGetCurrentEyePos; 
    end
    
    if stimulus.debug > 0
        if stimulus.debug >= 10
            disp(sprintf('Mouse position: %1.1f %1.1f',pos(1),pos(2)));
            stimulus.debug = 1;
        else
            stimulus.debug = stimulus.debug+1;
        end
    end
        
    
    % compute distance
    dist = hypot(pos(1),pos(2));
end

% Eye movement detection code
if (~stimulus.noeye) && (stimulus.eyewindow>0) && ~stimulus.live.dead
    if ~any(isnan(pos))
        if dist > stimulus.eyewindow && stimulus.live.eyeCount > 30
            disp('Eye movement detected!!!!');
            stimulus.live.dead = 1;
            return
        elseif dist > stimulus.eyewindow
            stimulus.live.eyeCount = stimulus.live.eyeCount + 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              HELPER FUNCTIONS                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function initStair()
global stimulus

for ai = 1:stimulus.attention.rotate
    stimulus.staircases{ai} = doStaircase('init','upDown',...
                'initialThreshold',0.40,...
                'initialStepsize',0.03,...
                'minThreshold=0.0001','maxThreshold=0.4','stepRule','pest',...
                'nTrials=80','maxStepsize=0.2','minStepsize=0.0001');
end

function resetStair()

global stimulus

for ai = 1:stimulus.attention.rotate
    s = stimulus.staircases{ai};
    if doStaircase('stop',s)
        disp('(afmap2) Staircase is being reset');
        s(end+1) = doStaircase('init',s(end));
        if s(end).s.threshold>1
            disp('(afmap2) Bad staircase threshold: setting to 1');
            s(end).s.threshold=1;
        elseif s(end).s.threshold<0
            disp('(afmap2) Bad staircase threshold: setting to 0.05');
            s(end).s.threshold=0.05;
        end
        stimulus.staircases{ai} = s;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localInitStimulus()
%%
global stimulus

fwhm_sd = 2*sqrt(2*log(2));

for ci = 1:length(stimulus.gratingContrasts)
    for si = 1:length(stimulus.designEccs)
        for ri = 1:length(stimulus.gratingRatios)
            sz = stimulus.gratingRatios(ri) * stimulus.designEccs(si);
            % use total degs / num to compute size
            for phase = 1:2
                grating = stimulus.gratingContrasts(ci) * 255/2 * mglMakeGrating(sz*4,sz*4,2/sz,0,(phase-1)*180) + 255/2;
                gauss = mglMakeGaussian(sz*4,sz*4,sz/fwhm_sd,sz/fwhm_sd);
                alphamask = repmat(grating,1,1,4);
                alphamask(:,:,4) = gauss*255;

                % make the grating
                stimulus.grating(ci,si,ri,phase) = mglCreateTexture(alphamask); % high contrast
                % make a gaussian (for when we display, make sure to use the
                % actual contrast for this setting or we'll fuck up later)
                gData = stimulus.gratingContrasts(ci)*255*ones(size(grating,1),size(grating,2),4);
                gData(:,:,4) = gauss*255;
                stimulus.gaussian(ci,si,ri,phase) = mglCreateTexture(gData);
            end
        end
    end
end

%% testing

% mglClearScreen;
% for si = 1:7
%     for ri = 1:3
%         mglBltTexture(stimulus.gaussian(2,si,ri,1),[randn*10 randn*10]);
%     end
% end
% mglFlush