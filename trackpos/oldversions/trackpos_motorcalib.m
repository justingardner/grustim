%        $Id: $
%      usage: trackpos motor calibration
%         by: Josh Ryu
%       date: 12/13/2021
%    purpose: 

% todo: make noiseLum more flexible

function myscreen = trackpos_motorcalib(varargin)
 

%getArgs(varargin,{'subjectID=s999','centerX=10','centerY=0','diameter=16'}); getArgs(varargin,{'subjectID=-1'});
% set up screen
myscreen = struct();
if isempty(mglGetSID)
    myscreen.subjectID  = -1;
else
    myscreen.subjectID  = mglGetSID;
end

%myscreen.displayName = 'debug'; myscreen.screenNumber = 1; 
%myscreen.screenWidth = 860; myscreen.screenHeight = 600; 
%myscreen.hideCursor = 1;
myscreen                = initScreen(myscreen);

%% parameters
% Experimenter parameters
%todo:  check these throughout the code!!
exp.debug               = 0; % debug code; no eyetracking`
exp.noeye               = 1; % 1 if no eyetracking (mouse for eye); 0 if there is eye tracking `
exp.showmouse           = 0; 
exp.grabframe           = 0; % to save frames of the task
exp.fixateCenter        = 0; % fixate center
exp.eyemousedebug       = 0; % debug eyetracker with mouse
exp.downsample_timeRes  = 1;
exp.usejoystick         = 0;

design_block = struct();
design_block.steady_thresh  = floor(2 * myscreen.framesPerSecond);%if cursor is 
design_block.waitsecs       = 2;
design_block.backLum        = 90; %160;%90;  % background luminance; units: luminance 
design_block.stimLum        = 255;
design_block.stimStd        = [0.4, 0.4]; %, 1.4, 0.4, 1.4];
design_block.time           = [300, 10]; % trial maximum time.
design_block.maxtrials      = [5, 300]; % maxmimum number of trials

randomize_order             = true;

%% Set up tasks 

selected_packages = {'stillblob', 'stillblob'};
design = load_motor_calib_packages(myscreen, selected_packages, randomize_order);
task = {};
task = configureExperiment(task, myscreen, design, design_block);

% count trials
numTrials        = 0; 
totaltime        = 0;
for phaseNum = 1:length(task{1})
    numTrials = numTrials + task{1}{phaseNum}.numTrials;
    totaltime = totaltime + task{1}{phaseNum}.numTrials * sum(task{1}{phaseNum}.segmax);
end
disp(['Approx task duration = ' num2str(totaltime/60/60) ' hours']);


%% initialize
% intiailize task
disp(' Initializing Task....')

for phaseN = 1:length(task{1})
    [task{1}{phaseN} myscreen] = initTask(task{1}{phaseN},myscreen,...
        @startSegmentCallback,@screenUpdateCallback,@responseCallback,@initTrialCallback,@endTrialCallback,[]);
end

% initialize stimulus
global stimulus; stimulus = struct;

disp(' Initializing Stimulus....') 

myscreen = initStimulus('stimulus',myscreen); % what does this do???
stimulus.exp = exp;
stimulus.design = design;
stimulus.stimStd = 0.4;
stimulus.stimLum = 122;
stimulus = trackposInitStimulus(stimulus,myscreen); %

if stimulus.exp.grabframe
    global frame
    frame = {};
end

%% Eye calibration
if ~exp.noeye && ~exp.debug
    disp(' Calibrating Eye ....')
    myscreen = eyeCalibDisp(myscreen);
    
    % let the user know
    disp(sprintf('(trackpos) Starting Run...'));
end

if stimulus.exp.usejoystick
global joy; joy = vrjoystick(1); % use simulink 3d animation to load joystick object
if isempty(joy)
    stimulus.exp.usejoystick = 0;
    exp = stimulus.exp;
    disp(' FAILED TO FIND JOYSTICK! MAKE SURE SIMULINK 3D ANIMATION PACKAGE IS INSTALLED AND THE JOYSTICK IS PROPERLY CONNECTED');
    disp(' USING MOUSE FOR TRACKING');
else
    joy_params              = struct();
    joy_params.maxv         = 0.2;
    joy_params.deadzone     = 0.02;
    joy_params.sensitivity  = 2;
    joy_params.poly_order   = 1.2;
    stimulus.joy_params     = joy_params;
end
end


%% run the task
disp(' Running Task....'); stimulus.t0 = mglGetSecs; % 

% let the experimentee know too...
% mglClearScreen(task{1}{1}.parameter.backLum/255);
mglTextDraw('task (trackpos_motorcalib) starting... ', [0 0.5])
mglTextDraw('Once the red cursor appears, move your mouse to the brightest point on the screen',[0 -0.5]);
mglTextDraw('When you are ready, press backtick to go to start task',[0 -1.5]);
mglFlush

if ~exp.showmouse, mglDisplayCursor(0);, end %hide cursor

phaseNum = 1;
stimulus.started = false; 
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);     % update the task
    myscreen = tickScreen(myscreen,task);     % flip screen
    if phaseNum > 1 && ~stimulus.started
        mglTextDraw('When you are ready, press backtick to go to next block',[0 -1.5]);
        mglFlush
    end
end

%% End task
disp(' End Task....')

myscreen = endTask(myscreen,task);
mglClose 
endScreen(myscreen); mglDisplayCursor(1) %show cursor

if stimulus.exp.grabframe
    savemat = '/Users/jryu/proj/grustim/trackpos/trackpos_motorcalib.mat';
    save(savemat, 'frame')
end

end

%% Initialize trials 
function [task myscreen] = initTrialCallback(task, myscreen)
    % before backtick..
    global stimulus 
    
    condNum = task.thistrial.thisphase; % each phase is a condition
    n = task.trialnum; % each phase is a condition
    stimulus.currdesign = stimulus.design(condNum);
    disp(['(trackpos_motorcalib) running stimulus ', stimulus.currdesign.name])
    
    % set stimulus position
    stimulus.position   = stimulus.currdesign.start_pos(n,:);
    stimulus.stimLum    = task.thistrial.stimLum;
    stimulus.backLum    = task.thistrial.backLum;
    stimulus.stimStd    = task.thistrial.stimStd;
    
    stimulus = trackposInitStimulus(stimulus,myscreen); %centerX,Y, diameter called by getArgs.
   
    %{ 
    mglClearScreen(0);
    mglBltTexture(stimulus.gaussian,stimulus.position);
    mglFlush
    stimulus        = updateTarget(stimulus,myscreen,task); % update position.
    pause(1/myscreen.framesPerSecond)
    end
    %}
    
    %% frame counter.
    task.thistrial.framecount = 0;
    
    if stimulus.exp.grabframe
        global frame
        %save('/Users/jryu/Dropbox/Stanford/Current/FYP/FYP talk/figures/trackposTask_nonoise_160back.mat', 'frame','-v7.3')
        %save('/Users/jryu/Dropbox/Stanford/Current/FYP/FYP talk/figures/trackposTask_noise_160back.mat', 'frame','-v7.3')
        frame = {};
        frame{task.segmax(1)*myscreen.framesPerSecond} = [];
    end

    if mod(task.trialnum,ceil(task.numTrials/20)) == 1
        disp(['(trackpos) '  num2str(task.trialnum/task.numTrials) ...
            '% finished: Trial ' num2str(task.trialnum) ' / ' num2str(task.numTrials)]);
    end
end

%% Start segment
function [task myscreen] = startSegmentCallback(task, myscreen)
% S1: Stimulus (10s)
global stimulus 

stimulus.started = true;

% log start time.
tic
stimulus.trial_starttime = datetime('now'); 

if ~stimulus.exp.usejoystick
    % set mouse position to the middle of the screen
    x_img = 0; y_img = 0; 
    % x_img = stimulus.position(1);  y_img = stimulus.position(2);
    x_screen = x_img*myscreen.screenWidth/myscreen.imageWidth + myscreen.screenWidth/2;
    y_screen = y_img*myscreen.screenHeight/myscreen.imageHeight + myscreen.screenHeight/2;
    mglSetMousePosition(ceil(x_screen),floor(y_screen), myscreen.screenNumber); % correct for screen resolution???
end

if ~stimulus.exp.showmouse, mglDisplayCursor(0);, end %hide cursor

% number of frames for which the stimulus is stead
stimulus.cursor_steady = 0;

end

function [task myscreen] = endTrialCallback(task,myscreen)
    global stimulus
    mglClearScreen(stimulus.backLum/255);
    mglFlush
    
    stimulus.started = false;
end

%% screen update
function [task myscreen] = screenUpdateCallback(task, myscreen)
% S1: Stimulus (10s)

%% Update Screen
global stimulus % call stimulus. Takes ~0.000013 s.
if stimulus.exp.usejoystick, global joy;,end
% stimulus.timedebug(9,task.thistrial.framecount+1) = mglGetSecs(stimulus.t0); % takes ~0.0087425s
mglClearScreen(stimulus.backLum/255);

if (task.thistrial.thisseg== 1)
    task.thistrial.framecount = task.thistrial.framecount + 1;
    
    if stimulus.exp.fixateCenter == 1 % fixation below others.
        mglGluAnnulus(0,0,0.2,0.3,stimulus.fixColor.response,60,1);
        mglGluDisk(0,0,0.1,rand(1,3),60,1);
    end
    
    % update stimulus position
    stimulus.velocity = stimulus.currdesign.update(...
        stimulus.position, stimulus.velocity, stimulus.currdesign, task.trialnum, myscreen, stimulus);
    stimulus.position = stimulus.position + stimulus.velocity;
    mglBltTexture(stimulus.gaussian,stimulus.position); % draw stimulus

    % **&display mouse position
    % (datetime('now')-stimulus.trial_starttime) < task.thistrial.waitsecs
    if toc < task.thistrial.waitsecs
        % set pointer to middle (green, can't move yet)
        
        if ~stimulus.exp.usejoystick
            mglSetMousePosition(ceil(myscreen.screenWidth/2),floor(myscreen.screenHeight/2), myscreen.screenNumber);
        end
        
        mglGluDisk(0, 0, 0.1, [0 1 0])
        
        % ***record stimulus position and mouse position  
        task.thistrial.trackResp(task.thistrial.framecount,:) = [nan, nan];
        
        if stimulus.exp.usejoystick
            [vx, vy] = joy2vel(joy, stimulus.joy_params, myscreen);
            if abs(vx) + abs(vy) > 0
                tic; % new tic 
                disp(['Delaying start, please reset joystick'])
            end
        end
    else 
        if ~stimulus.exp.usejoystick
            mInfo = mglGetMouse(myscreen.screenNumber);
            stimulus.pointer(1) = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
            stimulus.pointer(2) = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;
        else
            [vx, vy] = joy2vel(joy, stimulus.joy_params, myscreen);
            stimulus = update_pointer(stimulus, [vx, vy], myscreen);
            task.thistrial.trackJoy(task.thistrial.framecount,:)  = axis(joy);
        end
        
        mglGluDisk(stimulus.pointer(1), stimulus.pointer(2), 0.1, [1 0 0])
            
        task.thistrial.trackResp(task.thistrial.framecount,:) = stimulus.pointer;
        
        % jumpsegment if the mouse is still for a second    
        if sum(abs([stimulus.pointer] - stimulus.position)) < stimulus.stimStd/2
            stimulus.cursor_steady = stimulus.cursor_steady + 1;
        end
        
    end

    % ***record stimulus position and mouse position  
    task.thistrial.trackStim(task.thistrial.framecount,:) = stimulus.position;
    task.thistrial.trackTime(task.thistrial.framecount)   = mglGetSecs(stimulus.t0);
    
    % if stimulus is steady for long
    if stimulus.cursor_steady > task.thistrial.steady_thresh
        task = jumpSegment(task); 
    end
end

%% eye tracking
% track for task
% *** track for fixation???
if (~stimulus.exp.noeye) && any(task.thistrial.thisseg==[1])
    % mouse version for testing with no eyetracker
    if stimulus.exp.eyemousedebug
        mInfo = mglGetMouse(myscreen.screenNumber);
        degx = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
        degy = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;

        pos = [degx, degy];
    else  % check eye pos
        [pos,postime] = mglEyelinkGetCurrentEyePos; % is this in image coordinates?
    end
        
    task.thistrial.trackEye(task.thistrial.framecount,:)  = pos;
    task.thistrial.trackEyeTime(task.thistrial.framecount) = postime;
end



if stimulus.exp.grabframe && (task.thistrial.thisseg== 1)
    global frame; frame{task.thistrial.framecount} = mglFrameGrab;
end

end

%% Get response; do nothing. 
function [task myscreen] = responseCallback(task, myscreen)

global stimulus
 
end

function [stimx, stimy] = convertNearestPixel(myscreen,x_img,y_img)
    x_screen = x_img*myscreen.screenWidth/myscreen.imageWidth + myscreen.screenWidth/2;
    y_screen = y_img*myscreen.screenHeight/myscreen.imageHeight + myscreen.screenHeight/2;
    
    stimx = (ceil(x_screen)-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
    stimy = (floor(y_screen)-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight; % what is imagewidth???
end

%% stimuli to use
function design = load_motor_calib_packages(myscreen, selected_packages, randomize_order)
% selected_packages: cell of strings indicating which packages to use

design(length(selected_packages)) = struct('name',[], 'n',[], 'start_pos', [], 'vel', [], 'update', []);

%% possible predefined stimuli

% still target

if any(cellfun(@(a) strcmp(a, 'stillblob'),selected_packages))
    stillblob = struct();
    stillblob.name = 'stillblob';
    
    maxr    = min(myscreen.imageWidth/2,myscreen.imageHeight/2) - 5; % few degrees in 
    % maxr = 9.66; % for gru setup
    % radius  = [maxr/5, maxr/4, maxr/3]; % start from the edge 
    radius  = [1, 2, 3, 4, 5, 8]; % start from the edge 
    
    da = pi/4;
    angle  = 0:da:(2*pi-da);
    
    start_pos = [];
    for r = radius
        for a = angle
            start_pos = [start_pos; ...
                r * cos(a), r * sin(a)];
        end
    end
    
    stillblob.n = size(start_pos,1);
    
    if randomize_order
        start_pos = start_pos(randperm(stillblob.n),:);
    end    
    
    stillblob.start_pos     = start_pos;
    stillblob.vel           = zeros(stillblob.n,2);
    stillblob.update        = @(pos,vel,sb, n, m, s) update_stillblob(pos, vel, sb, n,m,s);
end

% linear target movement
if any(cellfun(@(a) strcmp(a, 'linearblob'), selected_packages))
    linearblob = struct();
    linearblob.name = 'linearblob';
    
    maxr    = min(myscreen.imageWidth/2,myscreen.imageHeight/2) - 5; % few degrees in 
    % maxr = 9.66
    % radius  = [maxr/5, maxr/4, maxr/3]; % start from the edge 
    radius  = [2.5]; % start from the edge 
    da      = pi/4;
    angle   = 0:da:(2*pi-da);
    speed   = [0.7, 1.3]; % degrees per second
    
    start_pos = [];
    vel = [];
    for s0 = speed
    for r = radius
    for a = angle
        s = s0 / myscreen.framesPerSecond; % linear motion`
        start_pos = [start_pos; r * cos(a), r * sin(a)];
        
        if mod(a,pi/2) == 0 % on cardianl 
            % go cardinal direction. 
            % e.g. pos left => vel up/down; 
            % e.g. pos up => vel left/right
            vel  = [vel; s * cos(a-pi/2), s * sin(a-pi/2)];
            
            % go the other cardinal direction
            start_pos = [start_pos; r * cos(a), r * sin(a)];
            vel  = [vel; s * cos(a+pi/2), s * sin(a+pi/2)];
            
        else % off cardinal
            % go the same direction
            vel = [vel; s * cos(a), s * sin(a)];
            
            % move cardinally
            start_pos = [start_pos; r * cos(a), r * sin(a)];
            vel = [vel; s * sign(cos(a)), 0];
            
            % move cardinally
            start_pos = [start_pos; r * cos(a), r * sin(a)];
            vel = [vel; 0, s * sign(sin(a))];
        end
    end
    end
    end
    
    linearblob.n = size(start_pos,1);
    
    if randomize_order
        perm = randperm(linearblob.n);
        start_pos   = start_pos(perm,:);
        vel         = vel(perm,:);
    end    
    
    linearblob.start_pos = start_pos;
    linearblob.vel = vel;
    linearblob.update = @(p,v, lb, n, m, s)  update_linearblob(p, v, lb, n, m,s);
end
  
%% add all the stimuli
for n = 1:length(selected_packages)
    eval(['design(', num2str(n), ') = ', selected_packages{n} ';'])
end

end

function stimulus = update_pointer(stimulus, vel, myscreen)
    pos = stimulus.pointer;
    [horz_out, vert_out] = check_oob(pos + vel, myscreen, stimulus);
    stimulus.pointer(1) = stimulus.pointer(1) + (1-horz_out)*vel(1);
    stimulus.pointer(2) = stimulus.pointer(2) + (1-vert_out)*vel(2);
end

function newvel = update_stillblob(pos, vel, stillblob, n, myscreen, stimulus)
    newvel = [0,0]; % do not move.
end

function newvel = update_linearblob(pos, vel, linearblob, n, myscreen, stimulus)
    newvel = linearblob.vel(n,:);
    
    % check out of bounds
    [horz_out, vert_out] = check_oob(pos + newvel, myscreen, stimulus);
    if any([horz_out, vert_out])
        newvel = [0,0]; % stop moving
    end
end



%% configure experiment
function task = configureExperiment(task, myscreen, design, design_block)
    offset = 0; 
    ncond = length(design);
    
    repeat_stims = 1;
    
    % all blocks have an training period and a tracking period
    for condNum = 1:ncond
        phaseNum = offset+ condNum;
        
        t = design_block.time(condNum);
        task{1}{phaseNum}.segmin           = [t]; % 10s max per trial
        task{1}{phaseNum}.segmax           = [t]; 
        
        maxN = design_block.maxtrials(condNum);
        task{1}{phaseNum}.numTrials        = min(maxN, repeat_stims * design(condNum).n);
        task{1}{phaseNum}.getResponse      = [0]; %segment to get response.
        task{1}{phaseNum}.synchToVol       = [0]; %segment to wait for backtick
        task{1}{phaseNum}.waitForBacktick  = 1; %wait for backtick before starting task phase
        
        % parameters
        task{1}{phaseNum}.parameter.steady_thresh   = design_block.steady_thresh;
        task{1}{phaseNum}.parameter.waitsecs        = design_block.waitsecs;
        task{1}{phaseNum}.parameter.backLum         = design_block.backLum; %160;%90;  % background luminance; units: luminance 
        task{1}{phaseNum}.parameter.stimLum         = design_block.stimLum;
        task{1}{phaseNum}.parameter.stimStd         = design_block.stimStd(condNum);

        % calculated parameters
        maxframes = ceil(task{1}{phaseNum}.segmax(1)*myscreen.framesPerSecond)+10; %
        task{1}{phaseNum}.randVars.calculated.trackStim    = nan(maxframes,2);
        task{1}{phaseNum}.randVars.calculated.trackResp    = nan(maxframes,2);
        task{1}{phaseNum}.randVars.calculated.trackEye     = nan(maxframes,2);
        task{1}{phaseNum}.randVars.calculated.trackJoy     = nan(maxframes,4);
        task{1}{phaseNum}.randVars.calculated.trackTime    = nan(1,maxframes);
        task{1}{phaseNum}.randVars.calculated.trackEyeTime = nan(1,maxframes); % for referencing edf file
    end
end
    