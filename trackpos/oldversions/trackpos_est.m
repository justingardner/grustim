%
%        $Id: $
%      usage: trackpos_est
%         by: Josh Ryu
%       date: 05/06/2021
%    purpose: estimate position of blob
%
% staircase not implemented yet
% consider: adding another noise period after the stimulus

% S1: random period of fixation (random ~0.5s)
% S2: stimulus period (stimdur s)
% S3: repsonse period + feedback (arbitrary)
% S4: feedback (1s)

% fixed position difference or samples from continuous distributions?

function myscreen = trackpos_est(varargin)

%% set up screen and experiment
% set input arguments
if isempty(mglGetSID)
    myscreen.subjectID  = -1;
else
    myscreen.subjectID  = mglGetSID;
    myscreen.saveData = 1;
end

myscreen.saveData       = 1; % save stimfile to data directory
myscreen = initScreen(myscreen);

% Experimenter parameters
exp                 = struct();
exp.noeye           = true;
exp.eyemousedebug   = false;
exp.showmouse       = false;
exp.grabframe       = false; 
exp.debug           = false;
exp.phasescrambleOn = true;
exp.backprecompute  = true;
exp.feedback        = true; 
exp.estim_horiz     = true; % do hoiztonal estimation
exp.estim_verti     = false; % do vertical estimation
exp.colorfix        = false;

%% task parameters
% S1: random period of fixation (random ~0.5s)
% S2: stimulus period (stimdur s)
% S3: repsonse period + feedback (arbitrary)
% S4: feedback (1s)

% stimulus and background
task{1}{1}.random               = 1;
task{1}{1}.parameter.backLum    = 90; %32;  % background luminance; units: fraction of full luminance 
task{1}{1}.parameter.noiseLum   = 32; % noise luminance, if there is one.

teststimLum                     = task{1}{1}.parameter.noiseLum*[1]; %SNR
teststimDur                     = [2/60 5/60 10/60 15/60]; %[2/60 5/60 10/60 15/60]; %frames/hz
posDiff                         = linspace(0,0.5,11); % in degs; minimum and maximum offset from fixation
% test trials in various ranges of min/max pos
% minPos
% maxPos

task{1}{1}.parameter.stimright  = [0, 1]; % 1 if stimulus is to the right of fixation
task{1}{1}.parameter.posDiff    = posDiff; % forst fixed values
task{1}{1}.parameter.stimLum 	= teststimLum;

% note: seglen are changed later
if ~exp.feedback 
    task{1}{1}.segmin           = [0.4 nan inf];
    task{1}{1}.segmax           = [0.8 nan inf]; 
    task{1}{1}.getResponse      = [0 0 1]; %segment to get response.
else
    task{1}{1}.segmin           = [0.4 nan inf 1];
    task{1}{1}.segmax           = [0.8 nan inf 1]; 
    task{1}{1}.getResponse      = [0 0 1 0]; %segment to get response.
    task{1}{1}.synchToVol       = [0 0 0 1]; %segmet to wait for backtick

end

%ntrials x l/r x stimDur x stimLum x posDiff
task{1}{1}.numTrials        = 10*2*length(teststimDur) * length(teststimLum)*length(posDiff); 
taskdur = (nansum(task{1}{1}.segmax(isfinite(task{1}{1}.segmax))) + 1 + max(teststimDur)) * ...
    task{1}{1}.numTrials/60/60; % approximate duration in hours
disp(['Approx task duration = ' num2str(taskdur) ' hours']);

task{1}{1}.waitForBacktick  = 1; %wait for backtick before starting each trial 

% calculated variables
maxframes = ceil((task{1}{1}.segmax(1)+max(teststimDur))...
    *myscreen.framesPerSecond)+10; % with some additional overflow

task{1}{1}.randVars.calculated.stimdur      = nan; % save stimdur as a calculated.

task{1}{1}.randVars.calculated.pos_estim    = nan(1,2); % nframes x 2 for position estimate
task{1}{1}.randVars.calculated.bgpermute    = nan(1,maxframes); % nframes x 1 for the background
task{1}{1}.randVars.calculated.stimON       = nan(1,maxframes); % nframes x 1 for the stimulus

task{1}{1}.randVars.calculated.trackTime    = nan(1,maxframes);
task{1}{1}.randVars.calculated.trackEye     = nan(maxframes,2);
task{1}{1}.randVars.calculated.trackEyeTime = nan(1,maxframes); % for referencing edf file

%% task blocks. 

% change stimulus duration
for trialN = 1:task{1}{1}.numTrials
    fixdur      = rand*(task{1}{1}.segmax(1) - task{1}{1}.segmin(1)) + task{1}{1}.segmin(1);
    stimdur     = teststimDur(randi(length(teststimDur)));
    if ~exp.feedback
        task{1}{1}.seglenPrecompute.seglen{trialN} = [fixdur stimdur inf];
    else
        task{1}{1}.seglenPrecompute.seglen{trialN} = [fixdur stimdur inf 1];
    end
end

% initializing task...
disp(' Initializing Task....')

for phaseN = 1:length(task{1})
    [task{1}{phaseN} myscreen] = initTask(task{1}{phaseN},myscreen,...
        @startSegmentCallback,@screenUpdateCallback,@responseCallback,@initTrialCallback);
end

%% initialize stimulus
global stimulus;
stimulus = [];

stimulus.exp = exp;
stimulus.teststimDur        = teststimDur; % not saved in the task.

stimulus = myInitStimulus(stimulus,myscreen,task);
myscreen = initStimulus('stimulus',myscreen); % what does this do???

if stimulus.exp.phasescrambleOn == 1;
    disp('Loading phase scrambled background noise...')

    tic
    if stimulus.exp.backprecompute == 1;
        savefile = '/Users/gru/data/trackpos/trackpos.mat';
        % savefile            = '/Users/joshua/data/trackpos_2afc/trackpos.mat'; % just use noise 1 and permute
        if ~exist(savefile,'file')
            error('need background file')
        end
        
        load(savefile,'backgroundnoise_rgb');

        if isfield(stimulus,'backnoise')
            for idx = 1:length(stimulus.backnoise)
                mglDeleteTexture(stimulus.backnoise{idx});
            end
        end

        % create all background textures and then load them later
        if stimulus.exp.debug 
            nnn = 200;
        else
            nnn = size(backgroundnoise_rgb,4);
        end
        for idx = 1:nnn %too big?? memory?
            stimulus.backnoise{idx} = mglCreateTexture(backgroundnoise_rgb(:,:,:,idx));
        end

        clearvars('backgroundnoise_rgb')
    end
    toc
end

%% Eye calibration
if ~stimulus.exp.noeye && ~ stimulus.exp.debug
    disp(' Calibrating Eye ....')
    myscreen = eyeCalibDisp(myscreen);
    
    % let the experimenter know
    disp(sprintf('(trackpos) Starting Run...'));
end

%% run the task
if stimulus.exp.grabframe, global frame; frame = {};, end
stimulus.t0 = mglGetSecs; % 

mglDisplayCursor(0); %hide cursor
mglClearScreen(task{1}{1}.parameter.backLum/255);
mglTextDraw('task (trackpos_est) starting... ', [0 0.5])
mglTextDraw('After the stimulus is presented, point the cursor to the center of the stimulus. Press 1 when you are done',[0 -0.5]);
mglFlush

if ~exp.showmouse, mglDisplayCursor(0);, end %hide cursor

phaseNum = 1;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);     % update the task
    myscreen = tickScreen(myscreen,task);     % flip screen
end

myscreen = endTask(myscreen,task);
mglClose 
endScreen(myscreen); mglDisplayCursor(1) %show cursor

if stimulus.exp.grabframe
    save('/Users/jryu/data/trackpos_est/taskSetup/frame.mat', 'frame')
end
end

%% Initialize trials; set staircuse
function [task myscreen] = initTrialCallback(task, myscreen)
    global stimulus
    
    %hide cursor
    if ~stimulus.exp.showmouse, mglDisplayCursor(0);, end 
    
    % print trial number every 5%. 
    if mod(task.trialnum,ceil(task.numTrials/20)) == 1
        disp(['Trial ' num2str(task.trialnum) ' / ' num2str(task.numTrials)]);
    end
    
    stimulus.stimLum    = task.thistrial.stimLum;
    stimulus.backLum    = task.thistrial.backLum;
    stimulus.noiseLum   = task.thistrial.noiseLum;
    
    task.thistrial.framecount = 0;
    
    task.thistrial.stimdur    = task.thistrial.seglen(2);
    
    stimulus = myInitStimulus(stimulus,myscreen,task); %centerX,Y, diameter called by getArgs.
    
    if stimulus.exp.phasescrambleOn == 1 && stimulus.exp.backprecompute == 1;
        nframes = ceil(sum(task.thistrial.seglen(1:2))*myscreen.framesPerSecond)+10; % with some additional overflow
        task.thistrial.bgpermute(1:nframes) = randi(length(stimulus.backnoise),nframes,1);
    end
        
    if stimulus.exp.grabframe
        global frame
        %save('/Users/jryu/Dropbox/Stanford/Current/FYP/FYP talk/figures/trackpos_2afc_32back_500ms.mat', 'frame','-v7.3')
        frame = {};
        frame{sum(task.segmax)*myscreen.framesPerSecond} = [];
    end

end

%% Start segment
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus

% center the mouse
if task.thistrial.thisseg == 3
    x_img = 0;  y_img = 0;
    x_screen = x_img*myscreen.screenWidth/myscreen.imageWidth + myscreen.screenWidth/2;
    y_screen = y_img*myscreen.screenHeight/myscreen.imageHeight + myscreen.screenHeight/2;
    mglSetMousePosition(x_screen,y_screen, myscreen.screenNumber);
end

end

%% screen update
function [task myscreen] = screenUpdateCallback(task, myscreen)
% S1: random period of fixation (~0.5s)
% S2: stimulus period (stimdur s)
% S3: repsonse period + feedback (1s + arbitrary)

global stimulus % call stimulus

% set background luminance
if task.thistrial.backLum > 1
    mglClearScreen(task.thistrial.backLum/255);
else
    mglClearScreen(task.thistrial.backLum);
end

task.thistrial.framecount = task.thistrial.framecount + 1;
task.thistrial.stimON(task.thistrial.framecount) = 0; %count stimulus

% add fixation
if stimulus.exp.colorfix
    % changing fixation colors
    % mglGluAnnulus(0,0,0.2,0.3,stimulus.fixColor,60,1);
    mglGluDisk(0,0,0.1,rand(1,3),60,1);
else
%     % blue fixation dot ??
%     mglGluDisk(0, 0, 0.1, [0 0 1],60,1); 
    
    if isfield(stimulus,'pointer') && ~isempty(stimulus.pointer)
        mglBltTexture(stimulus.pointer.img, [0,0]);
    else
        mglGluDisk(0, 0, 0.1, [1 0 0]) % small red dot (cursor)
    end
end

% inject noise, track time
if any(task.thistrial.thisseg == [1,2]) 
    if stimulus.exp.phasescrambleOn == 1 
        idx = task.thistrial.bgpermute(task.thistrial.framecount);
        mglBltTexture(stimulus.backnoise{idx},...
            [0 0 myscreen.imageWidth myscreen.imageHeight])
    end
    
    task.thistrial.trackTime(task.thistrial.framecount) = mglGetSecs(stimulus.t0);
end

% draw blob or response feedback
stim_pos = (2*task.thistrial.stimright-1)*task.thistrial.posDiff;
if task.thistrial.thisseg == 2 % stimulus
    task.thistrial.stimON(task.thistrial.framecount) = 1;
    mglBltTexture(stimulus.gaussian,[stim_pos 0]);
elseif task.thistrial.thisseg == 3 %response period; 
    % show cursor
    mInfo = mglGetMouse(myscreen.screenNumber);
    degx = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
    degy = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;
    
    if ~stimulus.exp.estim_horiz, degx = 0; end
    if ~stimulus.exp.estim_verti, degy = 0; end
    
    if isfield(stimulus,'pointer') && ~isempty(stimulus.pointer)
        mglBltTexture(stimulus.pointer.img, [degx, degy]);
    else
        mglGluDisk(0, 0, 0.1, [1 0 0]) % small red dot (cursor)
    end
    mglGluDisk(degx, degy, 0.1, [1 0 0]);
    
    task.thistrial.stimON(task.thistrial.framecount) = 1; %count stimulus
elseif task.thistrial.thisseg == 4 % feedback period
    mglGluDisk(stim_pos, 0, 0.1, [1 0 0]) ;    % draw center of blob
end

% track eye
if (~stimulus.exp.noeye) && any(task.thistrial.thisseg==[1,2]) 
    % mouse version for testing with no eyetracker
    if stimulus.exp.eyemousedebug
        mInfo = mglGetMouse(myscreen.screenNumber);
        degx = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
        degy = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;

        pos = [degx, degy];
    else  % check eye pos
        [pos,postime] = mglEyelinkGetCurrentEyePos; % is this in image coordinates?
    end
        
    task.thistrial.trackEye(task.thistrial.framecount,:)    = pos;
    task.thistrial.trackEyeTime(task.thistrial.framecount)  = postime;
end

if stimulus.exp.grabframe
    global frame; frame{task.thistrial.framecount} = mglFrameGrab;
end

end

%% Get response 
function [task myscreen] = responseCallback(task, myscreen)

global stimulus

% if the button 1 is pressed, record position
if task.thistrial.whichButton == 1
    % record position of the mouse
    mInfo = mglGetMouse(myscreen.screenNumber);
    degx = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
    degy = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;
    
    if ~stimulus.exp.estim_horiz, degx = 0; end
    if ~stimulus.exp.estim_verti, degy = 0; end
    task.thistrial.pos_estim = [degx, degy]; % save position estimates
    
    stim_pos = (2*task.thistrial.stimright-1)*task.thistrial.posDiff;

    disp(['Stimulus Hor. Pos.: ' num2str(stim_pos) '; ' ...
        'Response: ' num2str(degx)  ...
        '; Error: ' num2str(degx - stim_pos)]);
    
    % go to next segment
    task = jumpSegment(task);
end

end

%% Initialize stimulus, initialize blob
function stimulus = myInitStimulus(stimulus,myscreen,task)  
    % set standard deviation of stimulus
    
    % stimulus size
    if ~isfield(stimulus,'stimStd'), stimulus.stimStd = 1;,end %unit: imageX, in deg. 
    %GardnerLab: stimstd = 1; CSNL stimStd = 0.4.. (why...?)
    stimulus.patchsize = min(6*stimulus.stimStd,min(myscreen.imageWidth,myscreen.imageHeight));
    
    %stimulus initial position. uniform distribution across the screen
    x_img = min(3*stimulus.stimStd,1/3*myscreen.imageWidth)*(2*rand(1)-1); 
    y_img = min(3*stimulus.stimStd,1/3*myscreen.imageWidth)*(2*rand(1)-1);
    stimulus.position = [x_img, y_img];
        
    % stimulus luminance
    if ~isfield(stimulus,'stimLum'), stimulus.stimLum = 122;,end %unit: luminance
            
    % background noise
    if ~isfield(stimulus,'noiseLum'), stimulus.noiseLum = 122;,end; % unit: luminance
    
    % background luminance
    if ~isfield(stimulus,'backLum'), stimulus.backLum = 32;,end; % unit: luminance

    % initialize stimulus
    if isfield(stimulus,'gaussian'), mglDeleteTexture(stimulus.gaussian);, end 
    gaussian    =  mglMakeGaussian(stimulus.patchsize,stimulus.patchsize,...
        stimulus.stimStd,stimulus.stimStd)*(stimulus.stimLum);
    gaussian_rgb           = 255*ones(4,size(gaussian,2),size(gaussian,1),'uint8');
    gaussian_rgb(4,:,:)    = round(gaussian');
    gaussian_rgb           = uint8(gaussian_rgb);

    stimulus.gaussian = mglCreateTexture(gaussian_rgb);
    
    % fixation cross
    stimulus.fixColors.response = [1 1 1];
    stimulus.fixColors.correct  = [0 1 0];
    stimulus.fixColors.incorrect = [1 0 0];
end

