function [ myscreen ] = afmap( varargin )
%
% TODO
% * fix rotation order
% * fix fixation cross task (way too hard at 5,5)
%
%
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
attend = 0; run = 0; build = 0;
getArgs(varargin,{'scan=1','plots=0','noeye=0','debug=0','replay=0','attend=1','run=0','build=0'});
stimulus.scan = scan;
stimulus.plots = plots;
stimulus.noeye = noeye;
stimulus.debug = debug;
stimulus.replay = replay;
stimulus.overrideRun = run;
stimulus.attend = attend; % controls whether attention mode runs
stimulus.buildOverride = build;
if ~stimulus.attend
    warning('*****ATTENTION MODE IS DISABLED*****');
end
clear localizer invisible scan noeye task test2 attend build

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

    if ~isempty(mglGetSID) && isdir(sprintf('~/data/afmap/%s',mglGetSID))
        % Directory exists, check for a stimefile
        files = dir(sprintf('~/data/afmap/%s/1*mat',mglGetSID));

        if length(files) >= 1
            fname = files(end).name;

            s = load(sprintf('~/data/afmap/%s/%s',mglGetSID,fname));
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
                error('Cannot continue: stimfile parameters were generated with a different attention mode than you requested. You need to save the existing stimfiles into a backup folder');
            end
            clear s;
            disp(sprintf('(afmap) Data file: %s loaded.',fname));
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
    disp(sprintf('(afmap) This is scan #%i',stimulus.counter));
    disp(sprintf('(afmap) This is run #%i',stimulus.curRun));
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
    stimulus.attention.attendX = [0 4 4];
    stimulus.attention.attendY = [0 4 -4];

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
%         figure
%         colormap('gray');
%         caxis([0 1]);
% %         build.con = build.con>0;
%         tb = build.con/max(build.con(:));
%         for i = 1:720
%             imagesc(squeeze(tb(i,:,:)));
%             pause(.01);
%         end
        %%
        % test code end
        disp(sprintf('(afmap) Pre-build of build %i has finished (will be saved with stimfile).',bi));
        stimulus.builds{bi} = build;
    end
    disp(sprintf('(afmap) Pre-build complete. Created %i unique builds which will rotate every %i runs.',stimulus.build.uniques,stimulus.build.rotate));
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
        disp('(afmap) WARNING: New staircase');
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
    disp(sprintf('(afmap) Build %i selected',stimulus.build.curBuild));
    disp(sprintf('(afmap) Attending X: %i Y: %i selected',stimulus.attention.curAttendX,stimulus.attention.curAttendY));
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
    if stimulus.scan
        task{1}{1}.seglen(end) = 0.050; % make the last segment short so that it will synchtovol and hopefully align the runs
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

stimulus.curTrial = 0;

if ~stimulus.replay
    global fixStimulus %#ok<TLEV>

    fixStimulus.diskSize = 0.75;
    fixStimulus.fixWidth = 1;
    fixStimulus.fixLineWidth = 3;
    fixStimulus.stimTime = 0.35;
    fixStimulus.interTime = 1.4;
    fixStimulus.stairUsePest = 1;
    fixStimulus.responseTime = 2;
    fixStimulus.staircase = stimulus.staircase;
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
    disp('(afmap) Displaying plots');
    dispInfo(stimulus);
end

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%


function [task, myscreen] = startTrialCallback1(task,myscreen)
global stimulus

disp(sprintf('(afmap) Starting cycle %01.0f',(stimulus.curTrial/120)+1));
% disppercent(-1/120,'(afmap) Running: ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback1(task,myscreen)
%%
global stimulus

stimulus.curTrial = stimulus.curTrial + 1;

stimulus.live.t = squeeze(stimulus.grid.t(stimulus.curTrial,:));
stimulus.live.ecc = squeeze(stimulus.grid.ecc(stimulus.curTrial,:));
stimulus.live.con = squeeze(stimulus.grid.con(stimulus.curTrial,:));
stimulus.live.sz = squeeze(stimulus.grid.sz(stimulus.curTrial,:));
stimulus.live.ph = squeeze(stimulus.grid.ph(stimulus.curTrial,:));
stimulus.live.theta = squeeze(stimulus.grid.theta(stimulus.curTrial,:));
stimulus.grid.t(stimulus.curTrial) = mglGetSecs;

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
    stimulus.frames(:,:,stimulus.curTrial) = frame(:,:,1);
end

% disp(sprintf('(afmap) Starting trial %01.0f',stimulus.curTrial));
% disppercent(stimulus.curTrial/120,'(afmap) Running: ');

function drawFix(myscreen)

global fixStimulus;

if fixStimulus.trainingMode,mglClearScreen;end

if ~isempty(fixStimulus.displayText)
  mglBltTexture(fixStimulus.displayText,fixStimulus.displayTextLoc);
end
mglGluDisk(fixStimulus.pos(1),fixStimulus.pos(2),fixStimulus.diskSize*[1 1],myscreen.background,60);

mglFixationCross(fixStimulus.fixWidth,fixStimulus.fixLineWidth,fixStimulus.thisColor,fixStimulus.pos);

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
% global stimulus



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
        disp('(afmap) Staircase is being reset');
        s(end+1) = doStaircase('init',s(end));
        if s(end).s.threshold>1
            disp('(afmap) Bad staircase threshold: setting to 1');
            s(end).s.threshold=1;
        elseif s(end).s.threshold<0
            disp('(afmap) Bad staircase threshold: setting to 0.05');
            s(end).s.threshold=0.05;
        end
        stimulus.staircases{ai} = s;
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

%%
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
