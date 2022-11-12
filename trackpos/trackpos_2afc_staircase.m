%        $Id: $
%      usage: trackpos_2afc
%         by: Josh Ryu
%       date: 04/13/2021
%    purpose: 2 choice task for position of one blob against fixation.

% S1 (w): wait time while the main task calls this subtask.
% S2 (ITI): random period of fixation (random ~0.5s)
% S3 (s): stimulus period (stimdur s)
% S4 (md): mask ITI
% S5 (m): mask
% S6 (r): repsonse period (inf)
% S7 (fb): feedback (1s)

% todo: keep track of which segment is which number in a stimulus struct?
% => change code accordingly

function [task, myscreen] = trackpos_2afc_staircase(varargin)
%% set up screen and experiment
% set input arguments
if isempty(mglGetSID)
    myscreen.subjectID  = -1;
else
    myscreen.subjectID  = mglGetSID;
    myscreen.saveData = 1;
end

% /Users/JRyu/github/mgl/task/displays, 0001_dn0a22167c_220912.mat
myscreen.displayName    = 'vpixx';
myscreen.calibType      = 'Specify particular calibration';
myscreen.calibFilename  = '/Users/gru/proj/mgl/task/displays/0001_dn0a221834_221005.mat';%'0001_dn0a221834_221005.mat'; % '0001_dn0a22167c_220912.mat';
myscreen.saveData       = 1; % save stimfile to data directory
myscreen.datadir        = '/Users/gru/data/';

myscreen = initScreen(myscreen);

% set to argb2101010 pixel format
mglMetalSetViewColorPixelFormat(4);

% Experimenter parameters
exp                 = struct();
exp.debug           = false;
exp.noeye           = true;
exp.eyemousedebug   = false;
exp.showmouse       = false;
exp.phasescrambleOn = false;
exp.backprecompute  = false;
exp.feedback        = true; 
exp.feedback_center = true;  % feedback about the exact center
exp.estim_horiz     = true;  % do hoiztonal estimation
exp.estim_verti     = false; % do vertical estimation
exp.colorfix        = false; % colored fixation
exp.block_design    = false; % in each block, present all combinations of parameters
exp.noise_mask      = '/Users/gru/proj/grustim/trackpos/noise/white0.mat'; % in each block, present all combinations of parameters

%% task parameters
% stimulus and background
params            = struct();
params.backLum    = 0.4; %32;  % background luminance; units: fraction of full luminance 
params.noiseLum   = 0; % noise luminance, if there is one.

% main task parameters
params.stimLum      = [0.05, 0.1, 0.2, 0.4]; % [0.1,0.2,0.5] % [16,32,48,96]
if exp.debug
    params.stimDur      = [2/60, 15/60]; %[2/60 5/60 10/60 15/60]; %frames/hz
else
    params.stimDur      = [2/60, 3/60, 4/60, 6/60, 10/60, 15/60]; %[2/60 5/60 10/60 15/60]; %frames/hz
end

params.stimStd          = [1]; % [1,1.5]
params.stimColor        = 'k';

% mask parameters
params.maskDur          = 3/60; % mask duration
params.mask_TOff2MOn    = 5/60; % stimulus offset to mask onset (Neisser 1967)
params.maskLum          = 0.05;

% staircase parameters
trialpercond        = 40;
if exp.debug, trialpercond = 2; end

params.initThreshold    = 0.1;
params.initThresholdSd  = 0.1;

% count conditions
nconditions             = length(params.stimDur) * length(params.stimStd) * length(params.stimLum);
params.trialpercond     = trialpercond ;
params.numTrials        = trialpercond * nconditions;

disp(['Number of conditions = ' num2str(nconditions)]);

% stimulus and background
task{1}.random               = 1;
task{1}.parameter.backLum    = params.backLum; 
task{1}.parameter.noiseLum   = params.noiseLum;
task{1}.parameter.stimright  = [0,1];
task{1}.parameter.stimLum 	 = params.stimLum;
task{1}.parameter.stimStd 	 = params.stimStd;
task{1}.parameter.stimColor	 = params.stimColor;
task{1}.parameter.stimDur	 = params.stimDur;

% note: seglen are changed later
% note: rewrite this
if mglIsFile(exp.noise_mask)
    maskDur = params.maskDur;
    maskDel = params.mask_TOff2MOn;
    task{1}.segmin           = [0.1 0.4 inf maskDel maskDur inf 1];
    task{1}.segmax           = [0.1 0.8 inf maskDel maskDur inf 1]; 
    task{1}.getResponse      = [0 0 0 0 0 1 0]; %segment to get response.
else
    % basically skip mask period
    task{1}.segmin           = [0.1 0.4 inf inf 1];
    task{1}.segmax           = [0.1 0.8 inf inf 1]; 
    task{1}.getResponse      = [0 0 0 1 0]; %segment to get response.
end


if exp.block_design
    task{1}.numBlocks        = ceil(trialpercond / 2); % dont count stimright as condition %with some overflow
else
    task{1}.numTrials        = params.numTrials; % dont count stimright as condition %with some overflow
end

taskdur = (0.5 + max(params.stimDur) + 1 + 1) * nconditions * trialpercond/60/60; % approximate duration in hours
disp(['Approx task duration = ' num2str(taskdur) ' hours']);

task{1}.waitForBacktick  = 1;

% calculated variables
maxframes = ceil((task{1}.segmax(2)+max(params.stimDur))...
    *myscreen.framesPerSecond)+10; % with some additional overflow

task{1}.randVars.calculated.subjcorrect  = nan; 
task{1}.randVars.calculated.stimDur0     = nan; 
task{1}.randVars.calculated.posDiff      = nan; % set by the staircase 

task{1}.randVars.calculated.bgpermute    = nan(1,maxframes); % nframes x 1 for the background
task{1}.randVars.calculated.stimON       = nan(1,maxframes); % nframes x 1 for the stimulus

task{1}.randVars.calculated.trackTime    = nan(1,maxframes);
task{1}.randVars.calculated.trackEye     = nan(maxframes,2);
task{1}.randVars.calculated.trackEyeTime = nan(1,maxframes); % for referencing edf file

%% construct staircase 
tpnames    = ["backLum","noiseLum","stimLum","stimDur", "stimStd", "stimColor", "staircase"];
tparams    = cell(1,length(tpnames));
for backLum     = params.backLum
for noiseLum    = params.noiseLum
for stimLum     = params.stimLum
for stimDur     = params.stimDur
for stimStd     = params.stimStd
for stimColor	= params.stimColor
    staircase = doStaircase('init','quest',...
        ['initialThreshold=' num2str(params.initThreshold)],['initialThresholdSd=' num2str(params.initThresholdSd)],...
        'nTrials',trialpercond);
    tparams{1} = [tparams{1}; backLum];
    tparams{2} = [tparams{2}; noiseLum];
    tparams{3} = [tparams{3}; stimLum];
    tparams{4} = [tparams{4}; stimDur];
    tparams{5} = [tparams{5}; stimStd];
    tparams{6} = [tparams{6}; stimColor];
    tparams{7} = [tparams{7}; {staircase}];
end
end
end
end
end
end

%% stimulus struct.
global stimulus
stimulus = [];

stimulus.exp            = exp;
stimulus.params         = params; % not saved in the task.
stimulus.target         = trackposInitStimulus(stimulus,myscreen);

stimulus.fixColors.response = [1 1 1];
stimulus.fixColors.stim     = [0 1 0]; % green
stimulus.fixColors.est      = [1 0 0]; % red
stimulus.fixColors.afc      = [1 1 1]; % blue

stimulus.t0 = mglGetSecs; % keeps track of trackTime

stimulus.staircaseTable = table(tparams{:},'VariableNames', tpnames); % save staircase

myscreen = initStimulus('stimulus',myscreen); % what does this do???

% phase scrambled noise deprecated -- get rid of
if stimulus.exp.phasescrambleOn == 1
    disp('Loading phase scrambled background noise...')

    tic
    if stimulus.exp.backprecompute == 1
        savefile = '/Users/gru/proj/grustim/trackpos/trackpos.mat';
        %savefile = '/Users/jryu/data/trackpos/trackpos.mat'; 
        % savefile = '/Users/joshua/data/trackpos_2afc/trackpos.mat'; % just use noise 1 and permute
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

% noise mask
if mglIsFile(stimulus.exp.noise_mask)
    stimulus.noise_mask = load(stimulus.exp.noise_mask);
end

%% Eye calibration
if ~stimulus.exp.noeye && ~ stimulus.exp.debug
    disp(' Calibrating Eye ....')
    myscreen = eyeCalibDisp(myscreen);
    
    % let the experimenter know
    disp(sprintf('(trackpos) Starting Run...'));
end

%% initializing task...
disp(' Initializing Task....')

for phaseN = 1:length(task)
    [task{phaseN}, myscreen] = initTask(task{phaseN},myscreen,...
        @startSegmentCallback,@screenUpdateCallback,@responseCallback,@initTrialCallback);
end

%% run task
phaseNum = 1;phaseNum2=1;phaseNum3=1;
while (phaseNum <= length(task)) && ~myscreen.userHitEsc
    [task, myscreen, phaseNum]   = updateTask(task,myscreen,phaseNum);     % update the task
    myscreen = tickScreen(myscreen,task);     % flip screen
end

myscreen = endTask(myscreen,task);
mglClose; endScreen(myscreen); mglDisplayCursor(1) %show cursor

% save staircase
staircase = stimulus.staircaseTable;
save([myscreen.stimfile(1:end-4),'_staircase.mat'], 'staircase');

function [task, myscreen] = initTrialCallback(task, myscreen)
    global stimulus

    backLum     = task.thistrial.backLum;
    noiseLum    = task.thistrial.noiseLum;
    stimLum     = task.thistrial.stimLum;
    stimDur     = task.thistrial.stimDur;
    stimStd     = task.thistrial.stimStd;
    stimColor   = task.thistrial.stimColor;
    idx         = findCondIdx(stimulus.staircaseTable,backLum,noiseLum,stimLum,stimDur,stimStd,stimColor);
    
    stimulus.staircaseIdx = idx;
    [s, stimulus.staircaseTable.staircase{idx}] = doStaircase('testValue',stimulus.staircaseTable.staircase{idx});
    
    task.thistrial.posDiff = s;
    
    task.thistrial.seglen(3) = task.thistrial.stimDur;

    % noise mask
    if mglIsFile(stimulus.exp.noise_mask)     
        % seglen vs segmax?
        nframes = myscreen.framesPerSecond*task.thistrial.seglen(4) + 20; %/downsample_timeRes; 
        stimulus.noise_mask_trial = randi(size(stimulus.noise_mask.backgroundnoise_rgb,4),nframes,1); % sample with replacement
        for idx = 1:nframes
            midx    = stimulus.noise_mask_trial(idx);
            maskimg = stimulus.noise_mask.backgroundnoise_rgb(:,:,:,midx);
            maskimg(4,:,:) = stimulus.params.maskLum * maskimg(4,:,:);
            % maskimg(4,:,:) = 0 * maskimg(4,:,:);
            maskimg = permute(maskimg,[2,3,1]);
            stimulus.noise_mask_texture{idx} = mglMetalCreateTexture(maskimg);
        end
    end
    

%% Start segment
function [task, myscreen] = startSegmentCallback(task, myscreen)
global stimulus

% set flushMode based on noiseLum
if stimulus.exp.phasescrambleOn == 1 && stimulus.exp.backprecompute == 1 && task.thistrial.noiseLum >0
    myscreen.flushMode = 0;
else
    myscreen.flushMode = 1; %1
end

if task.thistrial.thisseg == 1
elseif task.thistrial.thisseg == 2
    if ~stimulus.exp.showmouse, mglDisplayCursor(0);, end 

    % start the task.
    stimulus.lum        = task.thistrial.stimLum;
    stimulus.std        = task.thistrial.stimStd;
    stimulus.color      = task.thistrial.stimColor;
    stimulus.backLum    = task.thistrial.backLum;
    stimulus.noiseLum   = task.thistrial.noiseLum;
    
    task.thistrial.framecount = 0;
    task.thistrial.stimDur0 = task.thistrial.seglen(3);
    
    stimulus.target = trackposInitStimulus(stimulus,myscreen); %centerX,Y, diameter called by getArgs.
    
    if stimulus.exp.phasescrambleOn == 1 && stimulus.exp.backprecompute == 1&& stimulus.noiseLum
        nframes = length(task.thistrial.bgpermute);
        task.thistrial.bgpermute(1:nframes) = randi(length(stimulus.backnoise),nframes,1);
    end
elseif task.thistrial.thisseg == 3
elseif task.thistrial.thisseg == 4
    % mask
    myscreen.flushMode = 0; % refresh every frame
    task.thistrial.framecount = 0; % restart framecount
    % actual stimulus length is updated by updateTask; this should be almost the same as Dur0
    task.thistrial.stimDur = task.thistrial.seglen(3); 
    %disp(['Segment duration error1: ', num2str(task.thistrial.stimDur - task.thistrial.stimDur0)]);
end

% blt screen once before screenUpdates loops
if task.thistrial.thisseg > 1
    [task, myscreen] = screenUpdateCallback(task, myscreen);
    mglFlush;
    if task.thistrial.thisseg == 3
        stimulus.start = mglGetSecs;
    elseif task.thistrial.thisseg == 4
        stimulus.length = mglGetSecs - stimulus.start;
        disp(['Segment duration error: ', num2str(stimulus.length - task.thistrial.stimDur)]);
    end
end

%% screen update
function [task, myscreen] = screenUpdateCallback(task, myscreen)

global stimulus % call stimulus

if task.thistrial.thisseg== 1
%% don't do anything
    
else
    
%% do the task
% set background luminance
if task.thistrial.backLum > 1
    mglClearScreen(stimulus.backLum/255);
else
    mglClearScreen(stimulus.backLum);
end

task.thistrial.framecount = task.thistrial.framecount + 1;
task.thistrial.stimON(task.thistrial.framecount) = 0; %count stimulus

% inject noise, track time, add fixation
if any(task.thistrial.thisseg == [2, 3]) 
    if stimulus.exp.phasescrambleOn == 1 
        idx = task.thistrial.bgpermute(task.thistrial.framecount);
        mglBltTexture(stimulus.backnoise{idx},...
            [0 0 myscreen.imageWidth myscreen.imageHeight])
    end
    
    task.thistrial.trackTime(task.thistrial.framecount) = mglGetSecs(stimulus.t0);
    if stimulus.exp.colorfix
        % changing fixation colors
        % mglGluAnnulus(0,0,0.2,0.3,stimulus.fixColor,60,1);
        mglGluDisk(0,0,0.1,rand(1,3),60,1);
    else
        mglGluDisk(0, 0, 0.1, stimulus.fixColors.stim,60,1); 
    end
elseif task.thistrial.thisseg == 5
    mglMetalBltTexture(stimulus.noise_mask_texture{task.thistrial.framecount},[0 0]);
elseif any(task.thistrial.thisseg == [6,7])
    % fixation indicating estimation task
    mglGluDisk(0, 0, 0.1, stimulus.fixColors.afc,60,1); 
end

% draw blob or response feedback
if task.thistrial.thisseg == 3 % stimulus
    stim_pos = (2*task.thistrial.stimright-1)*task.thistrial.posDiff;
    task.thistrial.stimON(task.thistrial.framecount) = 1;
    mglBltTexture(stimulus.target.img,[stim_pos 0]);
elseif task.thistrial.thisseg == 7 %feedback period
    % no fixation cross until response.
    mglGluAnnulus(0,0,0.2,0.5,stimulus.currfixColor,60,1);
    
    % feedback about presented position
    if stimulus.exp.feedback_center
        stim_pos = (2*task.thistrial.stimright-1)*task.thistrial.posDiff;
        mglGluDisk(stim_pos, 0, 0.1, [0 0 1]) ;    % draw center of blob
    end
end

% track eye
if (~stimulus.exp.noeye) && any(task.thistrial.thisseg==[2,3]) 
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

end
%% Get response 
function [task, myscreen] = responseCallback(task, myscreen)

global stimulus

% record responses. correct/incorrect
if any(task.thistrial.whichButton == [1 2])
    respIsRight = (task.thistrial.whichButton == 2);
    correct = (task.thistrial.stimright == respIsRight); % correct if first is right and response is 2.
    task.thistrial.subjcorrect = correct;
    
else % if they pressed other keys, record as nan, but still go on.
    correct = nan;
end

% change color of fixation for feedback.  
if isnan(correct)
    stimulus.currfixColor = [1 1 1]; % white
elseif correct
    stimulus.currfixColor = [0 1 0]; % green
else % incorrect
    stimulus.currfixColor = [1 0 0]; % red
end

% Output response to the screen. 
if task.thistrial.whichButton == 1, respSide = 'left';
elseif task.thistrial.whichButton == 2, respSide = 'right'; 
else respSide = 'missed'; 
end

if correct == 0, corrString = 'incorrect';
elseif correct == 1, corrString = 'correct';
else corrString = 'missed'; 
end

posdiff = (2*task.thistrial.stimright-1)*task.thistrial.posDiff;
disp(['(subtask_2c) Position difference: ' num2str(posdiff) '; ' ...
      'Response: ' respSide '; ' corrString])

idx = stimulus.staircaseIdx;
stimulus.staircaseTable.staircase{idx} = ...
    doStaircase('update',stimulus.staircaseTable.staircase{idx},correct);

task = jumpSegment(task); % go to next segment