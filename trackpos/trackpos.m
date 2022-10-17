%        $Id: $
%      usage: trackpos
%         by: Josh Ryu
%       date: 05/05/2021
%    purpose: 

% stimulus struct has object structs: 
% target, pointer struct and background,
% fixation and cell of other objects
% each object struct has fields: position and img. (and other properties
% native and updated by the dynamics modules
% tracking variables. (all other variables to keep track of should be
% defined in the task dynamics)
% task.thistrial.framecount
% task.thistrial.trackResp(task.thistrial.framecount,:)
% task.thistrial.trackStim(task.thistrial.framecount,:)
% task.thistrial.trackTime
% task.thistrial.trackEye(task.thistrial.framecount,:)
% task.thistrial.trackEyeTime(task.thistrial.framecount)


% todo: make noiseLum more flexible
% todo: currently only segment 1 is the tracking period (calculated
% variables stored only here including eye tracking, stimulus and response
% etc... make this more flexible? framecount is a bit awkward too.
% indicate tracking period or 
% todo: only load background once
% todo: move cursor into the task stimuli 
% todo: move things to startsegment. clear screen at inittrial
%

function myscreen = trackpos(varargin)
%getArgs(varargin,{'subjectID=s999','centerX=10','centerY=0','diameter=16'}); getArgs(varargin,{'subjectID=-1'});
% set up screen
myscreen = struct();
if isempty(mglGetSID)
    myscreen.subjectID  = -1;
else
    myscreen.subjectID  = mglGetSID;
end

%myscreen.screenWidth = 860; myscreen.screenHeight = 600; 
%myscreen.hideCursor = 1;
myscreen.displayName        = 'vpixx';
myscreen.calibType          = 'Specify particular calibration';
myscreen.calibFilename      = '0001_dn0a221834_221005.mat';
myscreen.calibFullFilename  = '/Users/gru/proj/mgl/task/displays/0001_dn0a221834_221005';
myscreen.saveData           = 1; % save stimfile to data directory
myscreen.datadir            = '/Users/gru/data/';

myscreen                = initScreen(myscreen);

% set to argb2101010 pixel format
mglMetalSetViewColorPixelFormat(4);

%% experiment parameters
% Experimenter parameters`
%todo:  check these throughout the code!!
exp.debug               = 0; % debug code
exp.trackEye            = 0; % 1 if no eyetracking (mouse for eye); 0 if there is eye tracking `
exp.showMouse           = 0; % show mouse during everything

exp.fixateCenter        = 1; % fixate center
exp.controlMethod       = 'mouse'; %todo: 1: mouse; 2: eye; 3:joystick

exp.downsample_timeRes  = 1; % downsample temporal resolution of background noise the by this factor.
exp.phasescrambleOn     = 0; % load background, if specified by task

exp.grabframe           = 0; % capture frames; todo: change. specify save directory

global stimulus; stimulus = struct;
stimulus.exp = exp;

task = {}; 

%% specify task design
% sb      = stillblob(myscreen, 'backLum', 0, 'maxtrialtime',30,'steady_thresh_frame',50);
% phase1  = sb.configureExperiment(stimulus,task,myscreen);

% no noise run
cps = {};
stimStepStdList     = [1];
stimStdList         = [1]; %[0.5, 1 ,2];
stimLums            = [0.05, 0.1, 0.2, 0.4]; %[0.1, 0.2, 0.5]; 
backLum             = 0.4;

ntrial_learn        = 3;
ntrials             = 12; % trials per condition


for stimStepStd = stimStepStdList
    % 3 learning phase
    cps{end+1} = brownian(myscreen, 'maxtrials', ntrial_learn, 'noiseLum', 0, 'backLum', backLum, ...
        'stimLum', 1, 'stimColor', 'k', 'stimStd', [1], 'stimStepStd', stimStepStd, ...
        'pointLum',1, 'pointColor', 'r','pointStd', 0, 'pointStepStd', 0, ...
        'bgfile', []);
    
    for stimStd = stimStdList
        ntrials_phase = length(stimLums) * ntrials;
        cps{end+1} = brownian(myscreen, 'maxtrials', ntrials_phase, 'noiseLum', 0, 'backLum', backLum, ...
            'stimLum', stimLums, 'stimColor','k', 'stimStd', stimStd, 'stimStepStd', stimStepStd, ...
            'pointLum', 1, 'pointColor', 'r','pointStd', 0, 'pointStepStd', 0, ...
            'bgfile', []);
        
        cps{end+1} = brownian(myscreen, 'maxtrials', 18, 'noiseLum', 0, 'backLum', backLum, ...
            'stimLum', stimLums, 'stimColor', 'r', 'stimStd', stimStd, 'stimStepStd', stimStepStd, ...
            'pointLum', 1, 'pointColor', 'r','pointStd', 0, 'pointStepStd', 0, ...
            'bgfile', []);
        
        cps{end+1} = brownian(myscreen, 'maxtrials', 18, 'noiseLum', 0, 'backLum', backLum, ...
            'stimLum', stimLums, 'stimColor', 'b', 'stimStd', stimStd, 'stimStepStd', stimStepStd, ...
            'pointLum', 1, 'pointColor', 'r','pointStd', 0, 'pointStepStd', 0, ...
            'bgfile', []);
    end
end
  
stimulus.task = cps;

% cp2     = brownian(myscreen, 'noiseLum', 0, 'stimLum', 96, 'stimStd', stimsize, 'stimStepStd', stimStepStd);
% cp2     = brownian(myscreen, 'noiseLum', 0, 'stimLum', 96, 'stimStd', stimsize, 'stimStepStd', stimStepStd);
% cp2     = brownian(myscreen, 'noiseLum', 0, 'stimLum', 96, 'stimStd', stimsize, 'stimStepStd', stimStepStd);
% 
% stimsize = 1; stimStepStd = 1; 
% cp3     = brownian(myscreen, 'noiseLum', 32, 'stimLum', 96, 'stimStd', 0,  'pointLum', 96, 'pointStd', stimsize);
% cp4     = brownian(myscreen, 'noiseLum', 32, 'stimLum', 64, 'stimStd', stimsize);
% cp5     = brownian(myscreen, 'noiseLum', 32, 'stimLum', 64, 'stimStd', 0,  'pointLum', 64, 'pointStd', stimsize);
% cp6     = brownian(myscreen, 'noiseLum', 32, 'stimLum', 0, 'stimStd', 0);
%stimulus.task = {circ1,circ2};
%stimulus.task = {cp2, cp3, cp4, cp5}; %, cp3, cp4, cp5, cp6};
% 
% % run the task at multiple radii.
% stimulus.task = {};
% thetaStd0_r1 = pi/3; % angular velocity at r=1; 
% for r = [2,4,6,8,10]
%     thetaStd0 = thetaStd0_r1/r; % keep linear velocity constant
%     circ   = circular(myscreen, 'noiseLum', 32, 'stimLum', 80, 'stimStd', stimsize, ...
%                       'r', [r], 'thetaStd0', thetaStd0, ...
%                       'maxtrials',5);
%     stimulus.task{end+1} = circ;                   
% end

%% configure task
task{1} = cell(length(stimulus.task),1);
for ts = 1:length(stimulus.task)
    task{1}{ts} = stimulus.task{ts}.configureExperiment(task,myscreen,stimulus);
end
task{1} = add_calculated_params(task{1}, myscreen);
totaldur = approximate_total_task_dur(task);

%% initialize
% intiailize task
disp(' Initializing Task....')

for phaseN = 1:length(task{1})
    [task{1}{phaseN} myscreen] = initTask(task{1}{phaseN},myscreen,...
        @startSegmentCallback,@screenUpdateCallback,@responseCallback,@initTrialCallback);
end

myscreen = initStimulus('stimulus',myscreen); % save the stimulus into stimfile

%% load background
% todo: make stimulus-specific background?
% todo: looping over 1 phase right now. do we need more?
for phaseN = 1:length(task{1})
    currtask = stimulus.task{phaseN};
    if ~isempty(currtask.bgfile) && exp.phasescrambleOn
        disp(['Loading phase scrambled background noise phase ' num2str(phaseN) '...'])

        savefile = currtask.bgfile;
        % savefile            = '/Users/joshua/data/trackpos_2afc/trackpos.mat'; % just use noise 1 and permute
        if ~exist(savefile,'file')
            disp(' THE SPECIFIED BACKGROUND FILE DOES NOT EXIST')
            someinput = input('press any button to continue');
            continue
        end

        load(savefile,'backgroundnoise_rgb');

        % create all background textures and then load them later
        if stimulus.exp.debug 
            nnn = 100;
        else
            nnn = size(backgroundnoise_rgb,4);
        end
        for idx = 1:nnn %too big?? memory?
            stimulus.backnoise{phaseN}{idx} = mglCreateTexture(backgroundnoise_rgb(:,:,:,idx));
        end

        clearvars('backgroundnoise_rgb')
        disp('Finished loading phase scrambled background noise...')
    end
end

%% grabframe
if stimulus.exp.grabframe
    global frame
    frame = {};
end

%% Eye calibration and check joystick

if strcmp(exp.controlMethod,'eye') && ~stimulus.exp.trackEye
    disp(' Need to track eye..  setting  exp.trackEye to true')
    stimulus.exp.trackEye = true;
    exp.trackEye = true;
end

if  strcmp(exp.controlMethod,'eye')
    [positions_target, positions_eye] = calibrateEye(myscreen, stimulus, true);
    eyecalib = struct();
    eyecalib.target{1}      = positions_target;
    eyecalib.eye{1}         = positions_eye;
    stimulus.eyecalib       = eyecalib;
elseif stimulus.exp.trackEye
    disp(' Calibrating Eye ....')
    myscreen = eyeCalibDisp(myscreen); % calibrate eye every time.
end

if ~strcmp(stimulus.exp.controlMethod, 'mouse') && ...
        ~strcmp(stimulus.exp.controlMethod, 'eye') &&...
        ~strcmp(stimulus.exp.controlMethod, 'joystick')
    disp('control method not found! using mouse for tracking')
    stimulus.exp.controlMethod = 'mouse';
end

if strcmp(stimulus.exp.controlMethod, 'joystick')
    stimulus = calibrateJoy(myscreen, stimulus);
end
exp = stimulus.exp;

%% run the task
disp(' Running Task....'); stimulus.t0 = mglGetSecs; % 

% let the experimentee know too...
% mglClearScreen(task{1}{1}.parameter.backLum/255);
%mglFlush;

mglBltTexture(mglText('task starting... '), [0 0.5]);
mglBltTexture(mglText('Track the target with the red pointer'),[0 -0.5]);
mglBltTexture(mglText('When you are ready, press backtick to go to next trial'),[0 -1.5]);
mglFlush;

if ~exp.showMouse, mglDisplayCursor(0);, end %hide cursor

phaseNum = 1;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    [task{1}, myscreen, newphaseNum] = updateTask(task{1},myscreen,phaseNum); % update the task
    if newphaseNum ~= phaseNum
        mglClearScreen(0.5);
        
        if strcmp(exp.controlMethod,'eye')
            [positions_target, positions_eye] = calibrateEye(myscreen, stimulus, true);
            stimulus.eyecalib.target{newphaseNum}    = positions_target;
            stimulus.eyecalib.eye{newphaseNum}       = positions_eye;
        end
        
        mglTextDraw('Press backtick to go to next trial',[0 0]);
        mglFlush;
    end
    phaseNum = newphaseNum;
    myscreen = tickScreen(myscreen,task);     % flip screen
end

%% End task
disp(' End Task....')

if isfield(stimulus,'backnoise')
    for idx = 1:length(stimulus.backnoise)
        mglDeleteTexture(stimulus.backnoise{idx});
    end
end

myscreen = endTask(myscreen,task);
mglClose 
endScreen(myscreen); mglDisplayCursor(1) %show cursor

if stimulus.exp.grabframe
    save('/Users/joshryu/Dropbox/GardnerLab/Stanford/Current/FYP/FYP talk/trackposTask.mat', 'frame')
end

end

%% Initialize trials 
function [task myscreen] = initTrialCallback(task, myscreen)
    global stimulus    
    phaseNum = task.thistrial.thisphase; % phaseNum?
       
    % save seed for generating random textures
    rng('shuffle','twister'); 
    s=rng; 
    task.thistrial.randomSeed = s.Seed;
    rng(task.thistrial.randomSeed,'twister');
    
    % noise: permute background and save
    if stimulus.exp.phasescrambleOn == 1 && ~isempty(stimulus.task{phaseNum}.bgfile)
        nframes = myscreen.framesPerSecond*task.segmax(1) + 20;%/downsample_timeRes; 
        task.thistrial.bgpermute(1:nframes) = randi(length(stimulus.backnoise{1}),nframes,1);
    end
    
    % task status report
    if mod(task.trialnum,ceil(task.numTrials/20)) == 1
        disp(['(trackpos) '  num2str(task.trialnum/task.numTrials) ...
            '% finished: Trial ' num2str(task.trialnum) ' / ' num2str(task.numTrials)]);
    end  

    
    % these should be updated by the dynamics function 
    stimulus.target = struct();
    stimulus.pointer = struct();
    stimulus.otherObjs = {};
    
    % task initTrial
    stimulus = stimulus.task{phaseNum}.initTrial(task, myscreen, stimulus);
end

%% Start segment
function [task myscreen] = startSegmentCallback(task, myscreen)
    global stimulus 
    phaseNum = task.thistrial.thisphase;
   
    stimulus = stimulus.task{phaseNum}.startSegment(task, myscreen, stimulus);
    
    if stimulus.task{phaseNum}.doTrack
        task.thistrial.framecount = 1;

        % set mouse to the middle
        % we extract only relative position (velocity) from the mouse positions
        if strcmp(stimulus.exp.controlMethod, 'mouse')
            x_screen = myscreen.screenWidth/2;
            y_screen = myscreen.screenHeight/2;
            mglSetMousePosition(ceil(x_screen),floor(y_screen), myscreen.screenNumber);

            mInfo = mglGetMouse(myscreen.screenNumber);
            [x,y] = screen2deg(mInfo.x, mInfo.y, myscreen);
            stimulus.pointer.mouse0 = [x,y];
        end

    else
        task.thistrial.framecount = [];
    end

    if stimulus.exp.grabframe
        global frame
        %save('/Users/jryu/Dropbox/Stanford/Current/FYP/FYP talk/figures/trackposTask_nonoise_160back.mat', 'frame','-v7.3')
        %save('/Users/jryu/Dropbox/Stanford/Current/FYP/FYP talk/figures/trackposTask_noise_160back.mat', 'frame','-v7.3')
        frame = {};
        frame{task.segmax(1)*myscreen.framesPerSecond} = [];
    end

end

%% screen update
function [task, myscreen] = screenUpdateCallback(task, myscreen)
    global stimulus
    phaseNum = task.thistrial.thisphase;
    
    % set background luminance
    if task.thistrial.backLum > 1
        mglClearScreen(task.thistrial.backLum/255);
    else
        mglClearScreen(task.thistrial.backLum);
    end

    % eye tracking
    if stimulus.exp.trackEye && stimulus.task{phaseNum}.doTrack
        % mouse version for testing with no eyetracker
        [pos,postime] = mglEyelinkGetCurrentEyePos; % is this in image coordinates?
        task.thistrial.trackEye(task.thistrial.framecount,:)    = pos;
        task.thistrial.trackEyeTime(task.thistrial.framecount)  = postime;
    end
    
    % move cursor        
    if stimulus.task{phaseNum}.movecursor 
        % **&display mouse position
        if strcmp(stimulus.exp.controlMethod, 'mouse')
            % extract how much the mouse moved mouse position
            mInfo = mglGetMouse(myscreen.screenNumber);
            [x,y] = screen2deg(mInfo.x, mInfo.y, myscreen);
            vx = x - stimulus.pointer.mouse0(1);
            vy = y - stimulus.pointer.mouse0(2);

            if norm([vx, vy]) > 1e-2 && stimulus.exp.debug
                disp(['mouse vel: ' num2str(vx) ','  num2str(vy)]) % todo: delete this line after checking
                disp(['pointer pos: ' num2str(stimulus.pointer.position)]) % todo: delete this line after checking
            end

            % reset mouse position
            if mInfo.x < myscreen.screenWidth * 0.1 || mInfo.x > myscreen.screenWidth *0.9 ...
                    || mInfo.y <  myscreen.screenHeight*0.1 || mInfo.y > myscreen.screenHeight*0.9
                mglSetMousePosition(ceil(myscreen.screenWidth/2),...
                    floor(myscreen.screenHeight/2), myscreen.screenNumber);      
                stimulus.pointer.mouse0 = [0,0];
            else
                stimulus.pointer.mouse0 = [x,y];
            end
            
            [x,y] = update_pointer(stimulus.pointer, [vx, vy], myscreen);
            stimulus.pointer.position = [x,y];
        elseif strcmp(stimulus.exp.controlMethod, 'joystick')
            [vx, vy] = joy2vel(stimulus.joy, stimulus.joy_params, myscreen);
            [x,y] = update_pointer(stimulus.pointer, [vx, vy], myscreen);
            stimulus.pointer.position = [x,y];
        elseif strcmp(stimulus.exp.controlMethod, 'eye')
            % cheat a bit take previous frame position...
            stimulus.pointer.position = ...
                task.thistrial.trackEye(task.thistrial.framecount,:); 
            
            if false && stimulus.exp.debug
                disp(['(trackpos) pointer pos: ' num2str(stimulus.pointer.position)]) % todo: delete this line after checking
                disp(['(trackpos) target pos: ' num2str(stimulus.target.position)]) % todo: delete this line after checking
            end
        end
    end

    % update task
    stimulus  = stimulus.task{phaseNum}.update(task, myscreen, stimulus);
    
    % if we are in the tracking period,  save tracking variables
    if stimulus.task{phaseNum}.doTrack
        % update framecount
        task.thistrial.trackResp(task.thistrial.framecount,:) = stimulus.pointer.position;
        task.thistrial.trackStim(task.thistrial.framecount,:) = stimulus.target.position;
        task.thistrial.trackTime(task.thistrial.framecount)   = mglGetSecs(stimulus.t0); 
        task.thistrial.framecount = task.thistrial.framecount + 1;
    end

    % display background
    if ~isempty(stimulus.task{phaseNum}.bgfile) && task.thistrial.noiseLum > 0 && stimulus.exp.phasescrambleOn
        mglBltTexture(...
            stimulus.backnoise{phaseN}{task.thistrial.bgpermute(task.thistrial.framecount)}, ...
            [0,0,myscreen.imageWidth, myscreen.imageHeight]);
    end
    
    % display target
    if ~isempty(stimulus.target)
        mglBltTexture(stimulus.target.img, stimulus.target.position);
    end
    
    % display pointer
    if ~isempty(stimulus.pointer)
        if isempty(stimulus.pointer.img)
            mglGluDisk(stimulus.pointer.position(1), stimulus.pointer.position(2), 0.2, [1,0,0],60,1); 
        else
            mglBltTexture(stimulus.pointer.img, stimulus.pointer.position);
        end
    end
    
    % display other objects
    for ij = 1:length(stimulus.otherObjs)
        mglBltTexture(stimulus.otherObjs{ij}.img, stimulus.otherObjs{ij}.position);
    end
    
    % display fixation
    if stimulus.exp.fixateCenter == 1 && stimulus.task{phaseNum}.displayFix % fixation below others.
        mglGluAnnulus(0,0,0.2,0.4,[1 1 1],60,1);
        mglGluDisk(0, 0, 0.1,rand(1,3),60,1); 

    end
    
    if ~stimulus.exp.showMouse, mglDisplayCursor(0);, end %hide cursor

    % grabframe
    if stimulus.exp.grabframe && any(task.thistrial.thisseg==stimulus.task{phaseNum}.doTrack)
        global frame; frame{task.thistrial.framecount} = mglFrameGrab;
    end
end

%% Get response; do nothing. 
function [task myscreen] = responseCallback(task, myscreen)

global stimulus
 
end

%% for position updates
function [x,y] = update_pointer(obj, vel, myscreen)
    pos = obj.position;
    
    [horz_out, vert_out] = check_oob(pos + vel, myscreen, obj.std);
    x = pos(1) + (1-horz_out)*vel(1);
    y = pos(2) + (1-vert_out)*vel(2);
end

function [x,y] = screen2deg(x_screen, y_screen, myscreen)
    x = (x_screen-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
    y = (y_screen-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;
end

function [x_screen,y_screen] = deg2screen(x, y, myscreen)
    x_screen = x*myscreen.screenWidth/myscreen.imageWidth + myscreen.screenWidth/2;
    y_screen = y*myscreen.screenHeight/myscreen.imageHeight + myscreen.screenHeight/2;
end

%% utility
function totaldur = approximate_total_task_dur(task)
    % count trials
    numTrials        = 0; 
    totaltime        = 0;
    for phaseNum = 1:length(task{1})
        numTrials = numTrials + task{1}{phaseNum}.numTrials;
        totaltime = totaltime + task{1}{phaseNum}.numTrials * sum(task{1}{phaseNum}.segmax);
    end
    totaldur = totaltime/60/60;
    disp(['Approx task duration = ' num2str(totaldur) ' hours']);
end

function task = add_calculated_params(task, myscreen)
    for phaseNum = 1:length(task)
        maxframes = ceil(task{phaseNum}.segmax(1)*myscreen.framesPerSecond) + 20;
        task{phaseNum}.randVars.calculated.randomSeed   = nan;
        task{phaseNum}.randVars.calculated.bgpermute    = nan(1,maxframes); % nframes x 1 for the background
        task{phaseNum}.randVars.calculated.initStim     = [nan nan];
        task{phaseNum}.randVars.calculated.trackStim    = nan(maxframes,2);
        task{phaseNum}.randVars.calculated.trackResp    = nan(maxframes,2);
        task{phaseNum}.randVars.calculated.trackEye     = nan(maxframes,2);
        task{phaseNum}.randVars.calculated.trackTime    = nan(1,maxframes);
        task{phaseNum}.randVars.calculated.trackEyeTime = nan(1,maxframes); % for referencing edf file
    end
end