%        $Id: $
%      usage: trackpos_2afc
%         by: Josh Ryu
%       date: 04/13/2021
%    purpose: 2 choice task for position of one blob against fixation.

% S1: wait time while the main task calls this subtask.
% S2: random period of fixation (random ~0.5s)
% S3: stimulus period (stimdur s)
% S4: delay to mask
% S5: mask
% S6: repsonse period (inf)
% S7: feedback (1s)

function [task, myscreen] = trackpos_sub_2afc(myscreen, params, exp)
% stimulus and background
task{1}.random               = 1;
task{1}.parameter.backLum    = params.backLum; 
task{1}.parameter.noiseLum   = params.noiseLum;
task{1}.parameter.stimright  = [0,1];
task{1}.parameter.stimLum 	 = params.stimLum;
task{1}.parameter.stimStd 	 = params.stimStd;
task{1}.parameter.stimColor	 = params.stimColor;
task{1}.parameter.stimDur	 = params.stimDur;

% todo: change for no feedback and no mask
if mglIsFile(exp.noise_mask)
    maskDur = params.maskDur;
    maskDel = params.mask_TOff2MOn;
    task{1}.segmin           = [inf 0.4 inf maskDel maskDur inf 1];
    task{1}.segmax           = [inf 0.8 inf maskDel maskDur inf 1]; 
    task{1}.getResponse      = [0 0 0 0 0 1 0]; %segment to get response.
else
    % basically skip mask period
    task{1}.segmin           = [0.1 0.4 inf inf 1];
    task{1}.segmax           = [0.1 0.8 inf inf 1]; 
    task{1}.getResponse      = [0 0 0 1 0]; %segment to get response.
end

global stimulus
% if exist(stimulus.exp.staircase, 'file') 
%     task{1}.randVars.calculated.posDiff         = nan; % save posDiff as random variable
% else
if stimulus.exp.staircase
    task{1}.randVars.calculated.posDiff = nan;
else
    task{1}.parameter.posDiff                   = params.posDiff; 
end

% if exp.block_design  % with some overflow.
%     task{1}.numBlocks        = params.trialpercond; % dont count stimright as condition %with some overflow
% else
task{1}.numTrials        = params.numTrials; % dont count stimright as condition %with some overflow

% calculated variables
maxframes = ceil((task{1}.segmax(2)+max(params.stimDur))...
    *myscreen.framesPerSecond)+10; % with some additional overflow

task{1}.randVars.calculated.subjcorrect  = nan; 
task{1}.randVars.calculated.stimDur0     = nan; 

task{1}.randVars.calculated.bgpermute    = nan(1,maxframes); % nframes x 1 for the background
task{1}.randVars.calculated.stimON       = nan(1,maxframes); % nframes x 1 for the stimulus

task{1}.randVars.calculated.trackTime    = nan(1,maxframes);
task{1}.randVars.calculated.trackEye     = nan(maxframes,2);
task{1}.randVars.calculated.trackEyeTime = nan(1,maxframes); % for referencing edf file

%% task blocks. 
% initializing task...
disp(' Initializing Task....')

for phaseN = 1:length(task)
    [task{phaseN}, myscreen] = initTask(task{phaseN},myscreen,...
        @startSegmentCallback,@screenUpdateCallback,@responseCallback,@initTrialCallback);
end


%% Start segment
function [task, myscreen] = initTrialCallback(task, myscreen)

global stimulus

backLum     = task.thistrial.backLum;
noiseLum    = task.thistrial.noiseLum;
stimLum     = task.thistrial.stimLum;
stimDur     = task.thistrial.stimDur; 
stimStd     = task.thistrial.stimStd;
stimColor   = task.thistrial.stimColor;

% if exist(stimulus.exp.staircase, 'file')
%     idx         = findCondIdx(stimulus.staircaseTable,backLum,noiseLum,stimLum,stimDur,stimStd,stimColor);
%         
%     trial_idx   = stimulus.staircaseTable.trial_idx_2afc(idx);
%     posDiff     = stimulus.staircaseTable.posDiffs_trial_2afc{idx}(trial_idx);
%     
%     task.thistrial.posDiff = posDiff;
%     stimulus.staircaseTable.trial_idx_2afc(idx) = stimulus.staircaseTable.trial_idx_2afc(idx) + 1;
%     
%     task.thistrial.seglen(3) = task.thistrial.stimDur;
% else
if stimulus.exp.staircase
    idx         = findCondIdx(stimulus.staircaseTable,backLum,noiseLum,stimLum,stimDur,stimStd,stimColor);
    stimulus.staircaseIdx = idx;
    [s, stimulus.staircaseTable.staircase{idx}] = doStaircase('gettestvalue',stimulus.staircaseTable.staircase{idx});
    task.thistrial.posDiff = s;
    task.thistrial.seglen(3) = task.thistrial.stimDur;
end

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

function [task, myscreen] = startSegmentCallback(task, myscreen)
global stimulus

% set flushMode based on noiseLum
if stimulus.exp.phasescrambleOn == 1 && stimulus.exp.backprecompute == 1 && task.thistrial.noiseLum >0
    myscreen.flushMode = 0;
else
    myscreen.flushMode = 1;
end

if task.thistrial.thisseg == 1
    % if we come back to beginning, means the previous trial finished. 
    % wait for the next task 
    myscreen.flushMode = 0; % update screen
    stimulus.currtask = 'done';
    
elseif task.thistrial.thisseg == 2
    if ~stimulus.exp.noeye, myscreen.flushMode = 0; end
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
    
    if stimulus.exp.phasescrambleOn == 1 && stimulus.exp.backprecompute == 1&& stimulus.noiseLum;
        nframes = length(task.thistrial.bgpermute);
        task.thistrial.bgpermute(1:nframes) = randi(length(stimulus.backnoise),nframes,1);
    end
elseif task.thistrial.thisseg == 3
    if ~stimulus.exp.noeye, myscreen.flushMode = 0; end
elseif task.thistrial.thisseg == 4
    task.thistrial.framecount = 0; % restart framecount
    task.thistrial.stimDur = task.thistrial.seglen(3); 
    % stimulus length recorded by updateTask; make sure is same as Dur0 % disp(['Segment duration error1: ', num2str(task.thistrial.stimDur - task.thistrial.stimDur0)])
elseif task.thistrial.thisseg == 5
    myscreen.flushMode = 0; % refresh every frame
end

% blt screen once before screenUpdates loops
if task.thistrial.thisseg > 1
    [task, myscreen] = screenUpdateCallback(task, myscreen);
    mglFlush;
    if task.thistrial.thisseg == 3
        stimulus.start = mglGetSecs;
    elseif task.thistrial.thisseg == 4
        stimulus.length = mglGetSecs - stimulus.start;
        disp(['Segment duration error: ', num2str(stimulus.length - task.thistrial.stimDur)])
    end
end

%% screen update
function [task, myscreen] = screenUpdateCallback(task, myscreen)

global stimulus % call stimulus

if task.thistrial.thisseg== 1
    %% waiting for task to start
    if strcmp(stimulus.currtask,'2c')
        stimulus.currtask = 'running 2afc';
        task = jumpSegment(task); 
    end

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
            % mglGluDisk(0,0,0.1,0.5+0.5*rand(1,3),60,1);
            mglMetalDots([0;0;0], [0.5+0.5*rand(3,1);1], [0.1;0.1], 1, 1);
        else
            % mglGluDisk(0, 0, 0.1, stimulus.fixColors.stim,60,1); 
            mglMetalDots([0;0;0], [stimulus.fixColors.stim';1], [0.1;0.1], 1, 1);
        end
        
    elseif task.thistrial.thisseg == 5
        mglMetalBltTexture(stimulus.noise_mask_texture{task.thistrial.framecount},[0 0]);
    elseif any(task.thistrial.thisseg == [6,7])
        % fixation indicating estimation task
        % mglGluDisk(0, 0, 0.1, stimulus.fixColors.afc,60,1); 
        mglMetalDots([0;0;0], [stimulus.fixColors.afc';1], [0.1;0.1], 1, 1);
    end

    
    % draw blob or response feedback
    if task.thistrial.thisseg == 3 % stimulus
        stim_pos = (2*task.thistrial.stimright-1)*task.thistrial.posDiff;
        task.thistrial.stimON(task.thistrial.framecount) = 1;
        mglBltTexture(stimulus.target.img,[stim_pos 0]);
    elseif task.thistrial.thisseg == 7 %feedback period
        % no fixation cross until response.
        % mglGluAnnulus(0,0,0.2,0.5,stimulus.currfixColor,60,1);
        mglMetalArcs([0;0;0], [stimulus.currfixColor'; 1], [0.3;0.5],[0;2*pi], 1);
        
        % feedback about presented position
        if stimulus.exp.feedback_center
            stim_pos = (2*task.thistrial.stimright-1)*task.thistrial.posDiff;
            %mglGluDisk(stim_pos, 0, 0.1, [0 0 1]) ;    % draw center of blob
            mglMetalDots([stim_pos;0;0], [1;1;1;1], [0.1;0.1], 1, 1);
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

if stimulus.exp.staircase && ~isnan(correct)
    idx = stimulus.staircaseIdx;
    stimulus.staircaseTable.staircase{idx} = ...
        doStaircase('update',stimulus.staircaseTable.staircase{idx},correct);
end

task = jumpSegment(task); % go to next segment