%
%        $Id: $
%      usage: trackpos_2afc
%         by: Josh Ryu
%       date: 04/13/2021
%    purpose: 2 choice task for position of one blob against fixation.
%
% staircase not implemented yet
% consider: adding another noise period after the stimulus

% S1: random period of fixation (random ~0.5s)
% S2: stimulus period (stimdur s)
% S3: repsonse period + feedback (1s + arbitrary)


function myscreen = trackpos_2c(varargin)

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
exp.noeye           = false;
exp.eyemousedebug   = false;
exp.grabframe       = false; 
exp.staircase       = false; 
exp.debug           = false;

%% task parameters
% S1: random period of fixation (~0.5s)
% S2: stimulus period (stimdur s)
% S3: repsonse period + feedback (1s + arbitrary)

% stimulus and background
task{1}{1}.parameter.backLum    = 90; %32;  % background luminance; units: fraction of full luminance 
task{1}{1}.parameter.noiseLum   = 32; % noise luminance, if there is one.
task{1}{1}.parameter.stimLum    = 255 - task{1}{1}.parameter.backLum;  % stimulus luminance (out of 255)
% teststimLum                     = linspace(task{1}{1}.parameter.stimLum, task{1}{1}.parameter.noiseLum,3);
teststimLum                     = task{1}{1}.parameter.noiseLum*[0.5, 1, 1.5, 2]; %SNR
teststimDur                     = [2/60 5/60 10/60 15/60]; %frames/hz
posDiff                         = [0.05, 0.1, 0.15, 0.2];

task{1}{1}.parameter.stimright  = [0, 1]; % 1 if stimulus is to the right of fixation
task{1}{1}.parameter.posDiff    = posDiff; % forst fixed values
task{1}{1}.parameter.stimLum 	= teststimLum;

task{1}{1}.segmin           = [0.4 nan 1];
task{1}{1}.segmax           = [0.8 nan 1]; 
%ntrials x l/r x stimDur x stimLum x posDiff
task{1}{1}.numTrials        = 22*2*length(teststimDur) * length(teststimLum)*length(posDiff); 
taskdur = 2 * task{1}{1}.numTrials/60/60; % approximate duration in hours
disp(['Task duration = ' num2str(taskdur) ' hours']);

task{1}{1}.getResponse      = [0 0 1]; %segment to get response.
task{1}{1}.synchToVol       = [0 0 1]; % segmet to wait for backtick
task{1}{1}.waitForBacktick  = 1; %wait for backtick before starting each trial 

% Run fixed intervals (1) or staircase (0)
if exp.staircase
    task{1}{1}.runfixedint  = 0; % staircase
else
    task{1}{1}.runfixedint  = 1; % fixed interval
end
stimulus.phaseScramble          = 1;
stimulus.backprecompute         = 1;

% calculated variables
task{1}{1}.randVars.calculated.subjcorrect  = nan; 

maxframes = ceil((task{1}{1}.segmax(1)+max(teststimDur))...
    *myscreen.framesPerSecond)+10; % with some additional overflow

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
    task{1}{1}.seglenPrecompute.seglen{trialN} = [fixdur stimdur 1];
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

% for staircase
if stimulus.exp.staircase
    stimulus.stairUp        = 1;
    stimulus.stairDown      = 2;
    stimulus.stairStepSize  = 0.01;
    stimulus.stairUseLevitt = 0;
    stimulus.stairUsePest   = 1;    % use PEST
    stimulus.stairRep       = 50;   % repeat staircase every [stairRep] trials
    stimulus.stairN         = 0;    % keeps track of how many staircases they did
    stimulus.threshold      = 1;    % starting threshold?
else % run fixed intervals (set the values here) 
end

stimulus.teststimDur        = teststimDur;

stimulus.phasescrambleOn = 1;
stimulus.backprecompute = 1;

stimulus = myInitStimulus(stimulus,myscreen,task);
myscreen = initStimulus('stimulus',myscreen); % what does this do???

if stimulus.phasescrambleOn == 1;
    disp('Loading phase scrambled background noise...')

    tic
    if stimulus.backprecompute == 1;
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
mglTextDraw('task (trackpos_2afc) starting... ', [0 0.5])
mglTextDraw('Press 1 if the second stimulus is to the left of the first stimulus; and 2 otherwise',[0 -0.5]);
mglFlush

phaseNum = 1;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);     % update the task
    myscreen = tickScreen(myscreen,task);     % flip screen
end

myscreen = endTask(myscreen,task);
mglClose 
endScreen(myscreen); mglDisplayCursor(1) %show cursor

if stimulus.exp.grabframe
    save('/Users/jryu/data/trackpos_2afc/taskSetup/frame.mat', 'frame')
    % save('/Users/jryu/data/trackpos_2afc/taskSetup/frame.mat', 'frame')
    % save('/Users/joshryu/Dropbox/GardnerLab/data/trackpos_2afc/taskSetup/frame.mat', 'frame')
end
end

%% Initialize trials; set staircuse
function [task myscreen] = initTrialCallback(task, myscreen)
    global stimulus
    
    % print trial number every 5%. 
    if mod(task.trialnum,ceil(task.numTrials/20)) == 1
        disp(['Trial = ' num2str(task.trialnum) ' / ' num2str(task.numTrials)]);
    end
    
    stimulus.stimLum    = task.thistrial.stimLum;
    stimulus.backLum    = task.thistrial.backLum;
    stimulus.noiseLum   = task.thistrial.noiseLum;
    
    task.thistrial.framecount = 0;
    
    task.thistrial.stimdur    = task.thistrial.seglen(2);
    
    if task.runfixedint == 0 % run staircase
        if task.trialnum == 1 %does the task 
            stimulus.stairN = 0;
        end

        %initialize staircase  
        if mod(stimulus.stairN, stimulus.stairRep) == 0
            stimulus = initStaircase(stimulus);
        end

        [stimulus.threshold, stimulus.staircase] = doStaircase('testValue',stimulus.staircase); %left threshold
    %     [stimulus.threshold(2) stimulus.staircase(2)] = doStaircase('testValue',stimulus.staircase(2)); %right threshold
    else % run fixed interval
        
    end
    
    stimulus = myInitStimulus(stimulus,myscreen,task); %centerX,Y, diameter called by getArgs.
    
    if stimulus.phasescrambleOn == 1 && stimulus.backprecompute == 1;
        nframes = ceil(sum(task.thistrial.seglen(1:end-1))*myscreen.framesPerSecond)+10; % with some additional overflow
        task.thistrial.bgpermute(1:nframes) = randi(length(stimulus.backnoise),nframes,1);
    end
        
    if stimulus.exp.grabframe
        global frame
        %save('/Users/jryu/Dropbox/Stanford/Current/FYP/FYP talk/figures/trackpos_2afc_32back_500ms.mat', 'frame','-v7.3')
        frame = {};
        frame{sum(task.segmax)*myscreen.framesPerSecond} = [];
    end

end

function stimulus = initStaircase(stimulus)
% set up left and right staircase
if stimulus.stairUseLevitt
    stimulus.staircase = doStaircase('init','upDown','nup',stimulus.stairUp,...
        'ndown',stimulus.stairDown,'initialThreshold',stimulus.threshold,...
        'initialStepsize',stimulus.stairStepSize,'testType=levitt','minThreshold',stimulus.stairStepSize,'maxThreshold',45); % maximum has to be 45. 
%     stimulus.staircase(2) = doStaircase('init','upDown','nup',stimulus.stairUp,...
%         'ndown',stimulus.stairDown,'initialThreshold',stimulus.threshold(2),...
%         'initialStepsize',stimulus.stairStepSize,'testType=levitt','minThreshold',stimulus.stairStepSize,'maxThreshold',45); % maximum has to be 45. 
elseif stimulus.stairUsePest
    stimulus.staircase = doStaircase('init','upDown','nup',stimulus.stairUp,...
        'ndown',stimulus.stairDown,'initialThreshold',stimulus.threshold,...
        'initialStepsize',stimulus.stairStepSize,'testType=pest','minThreshold',stimulus.stairStepSize,'maxThreshold',45);
%     stimulus.staircase(2) = doStaircase('init','upDown','nup',stimulus.stairUp,...
%         'ndown',stimulus.stairDown,'initialThreshold',stimulus.threshold(2),...
%         'initialStepsize',stimulus.stairStepSize,'testType=pest','minThreshold',stimulus.stairStepSize,'maxThreshold',45);
else
    stimulus.staircase = doStaircase('init','upDown','nup',stimulus.stairUp,...
        'ndown',stimulus.stairDown,'initialThreshold',stimulus.threshold,...
        'initialStepsize',stimulus.stairStepSize,'minThreshold',0,'maxThreshold',45);
%     stimulus.staircase(2) = doStaircase('init','upDown','nup',stimulus.stairUp,...
%         'ndown',stimulus.stairDown,'initialThreshold',stimulus.threshold(2),...
%         'initialStepsize',stimulus.stairStepSize,'minThreshold',0,'maxThreshold',45);
end

end

%% Start segment
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus

if task.thistrial.thisseg == 3 %response feedback
    stimulus.fixColor = stimulus.fixColors.response;
end

end

%% screen update
function [task myscreen] = screenUpdateCallback(task, myscreen)
% S1: random period of fixation (~0.5s)
% S2: stimulus period (stimdur s)
% S3: repsonse period + feedback (1s + arbitrary)

global stimulus % call stimulus

mglClearScreen(stimulus.backLum/255);

task.thistrial.framecount = task.thistrial.framecount + 1;
task.thistrial.stimON(task.thistrial.framecount) = 0; %count stimulus

% add fixation
mglGluDisk(0, 0, 0.1, [1 0 0]) % small red dot (cursor)

% draw blob or response feedback
if task.thistrial.thisseg == 2 % stimulus
    task.thistrial.stimON(task.thistrial.framecount) = 1;
    pos = (2*task.thistrial.stimright-1)*task.thistrial.posDiff;
    mglBltTexture(stimulus.gaussian,[pos 0]);
elseif task.thistrial.thisseg == 3 %response feedback
    % no fixation cross until response.
    if any(stimulus.fixColor ~= stimulus.fixColors.response)
        mglGluAnnulus(0,0,0.2,0.3,stimulus.fixColor,60,1);
    end
end

% inject noise, track time
if any([1,2] == task.thistrial.thisseg) 
    if stimulus.phasescrambleOn == 1 
        idx = task.thistrial.bgpermute(task.thistrial.framecount);
        mglBltTexture(stimulus.backnoise{idx},...
            [0 0 myscreen.imageWidth myscreen.imageHeight])
    end
    
    task.thistrial.trackTime(task.thistrial.framecount) = mglGetSecs(stimulus.t0);
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

% record responses. correct/incorrect
if any(task.thistrial.whichButton == [1 2])
    respIsRight = (task.thistrial.whichButton == 2);
    correct = (task.thistrial.stimright == respIsRight); % correct if first is right and response is 2.
    task.thistrial.subjcorrect = correct;
    
    if task.runfixedint == 0 
        stimulus.stairN = stimulus.stairN+1; %count how many times 
    end
else
    stimIs1 = task.thistrial.stimright; 
    respIs1 = nan; correct = nan;
end

% change color of fixation for feedback.  
if isnan(correct)
    stimulus.fixColor = stimulus.fixColors.response;
elseif correct
    stimulus.fixColor = stimulus.fixColors.correct;
else % incorrect
    stimulus.fixColor = stimulus.fixColors.incorrect;
end

% change fixation color (doesn't this go back to the screenupate function?)
mglGluAnnulus(0,0,0.2,0.3,stimulus.fixColor,60,1); 

% Output response to the screen. 
if task.thistrial.whichButton == 1, respSide = 'left';
elseif task.thistrial.whichButton == 2, respSide = 'right'; end

if correct == 0, corrString = 'incorrect';
elseif correct == 1, corrString = 'correct';
else, corrString = 'no response'; end

%stimulus.staircase(stimulus.leftcorrect+1) = doStaircase('update',stimulus.staircase(stimulus.leftcorrect+1),correctIncorrect,abs(task.thistrial.dirDiff));
%[stimulus.threshold(stimulus.leftcorrect+1), stimulus.staircase(stimulus.leftcorrect+1)] = doStaircase('testValue',stimulus.staircase(stimulus.leftcorrect+1));

if task.runfixedint == 0
    stimulus.staircase = doStaircase('update',stimulus.staircase,correct,abs(task.thistrial.posDiff));
    % [stimulus.threshold, stimulus.staircase] =
    % doStaircase('testValue',stimulus.staircase); % I do this in the
    % beginning of segment
end

posdiff = (2*task.thistrial.stimright-1)*task.thistrial.posDiff;
disp(['Position difference: ' num2str(posdiff) '; ' ...
    'Response: ' respSide '; ' corrString])

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

