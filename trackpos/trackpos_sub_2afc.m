%        $Id: $
%      usage: trackpos_2afc
%         by: Josh Ryu
%       date: 04/13/2021
%    purpose: 2 choice task for position of one blob against fixation.

% S1: wait time while the main task calls this subtask.
% S2: cue (1s)
% S3: random period of fixation (random ~0.5s)
% S4: stimulus period (stimdur s)
% S5: delay to mask
% S6: mask
% S7: repsonse period (inf)
% S8: feedback (1s)

function [task, myscreen] = trackpos_sub_2afc(myscreen, params, exp)
% stimulus and background
task{1}.random               = 1;

fieldnames = fields(params.task);
for fidx = 1:length(fieldnames) % check that all parameters needed are there
    fieldname = fieldnames{fidx};
    task{1}.parameter.(fieldname) = params.task.(fieldname);
end
task{1}.parameter.stimright  = [0,1]; % or polar angle?

% todo: change for no feedback and no mask
if exp.feedback, fb = 0.1;, else, fb = 0;,end

if mglIsFile(exp.noise_mask)
                             % wait, cue, fix, stim, del, mask, resp, feedback 
    task{1}.segmin           = [inf   1   0.5  inf   inf   inf  inf    fb];
    task{1}.segmax           = [inf   1   0.5  inf   inf   inf  inf    fb]; 
    task{1}.getResponse      = [0 0 0 0 0 0 1 0]; %segment to get response.
else
                             % wait, cue, fix, stim, resp, feedback 
    task{1}.segmin           = [inf 1 0.4 inf inf fb];
    task{1}.segmax           = [inf 1 0.8 inf inf fb]; 
    task{1}.getResponse      = [0 0 0 0 1 0]; %segment to get response.
end

global stimulus

% presentation schedule
if strcmp(param.presSched, 'staircase')
    % set up staircase
    task{1}.randVars.calculated.posDiff = nan;

    [paramnames, paramvals]     = countconditions(params.task);
    condition_combs             = allcomb(paramvals{:});
    allparamnames               = fields(params.task);
    tparams                     = cell(1,length(allparamnames)+1);
    tpnames                     = allparamnames;
    tpnames{end+1}              = 'staircase';
    
    if isfield(params.staircase, 'staircase_init')
        if isfile(params.staircase_init)
            disp(["(trackpos_sub_2afc) Initializing staircase with file " params.staircase.staircase_init])
            a = load(params.staircase.staircase_init);
            saved_staircase = a.staircase;
        else
            disp('staircase initialization file not found')
        end
    end

    %todo: check this.
    for combidx = 1:size(condition_combs,1)
        thistrial = struct();
        for aparamidx = 1:length(allparamnames)
            pname = allparamnames{aparamidx};
            paramidx = find(strcmp(pname,paramnames));
            if isempty(paramidx)
                thistrial.(pname) = params.task.(pname);
                tparams{aparamidx} = [tparams{aparamidx}; params.task.(pname)];
            else
                thistrial.(pname) = condition_combs(combidx,paramidx);
                tparams{aparamidx} = [tparams{aparamidx}; condition_combs(combidx,paramidx)];
            end
        end

        if isfield(params.staircase, 'staircase_init') && isfile(params.staircase.staircase_init)
            saved_idx           = findCondIdx(saved_staircase, thistrial);
            disp(['saved_idx = ' num2str(saved_idx)])
            if isempty(saved_idx)
                disp("(trackpos_sub_2afc) WARNING: the staircase file provided does not have a matching condition")
                init_thresh = params.staircase.initThreshold;
            else
                saved_idx_staircase = saved_staircase.staircase{saved_idx};
                init_thresh = saved_idx_staircase.s.pThreshold;
            end
        else
            init_thresh = params.staircase.initThreshold;
        end

        staircase = doStaircase('init','quest','nTrials',params.trialpercond,...
            ['initialThreshold=' num2str(init_thresh)], ...
            ['initialThresholdSd=' num2str(params.staircase.initThresholdSd)],...
            'verbose=0');
        tparams{end} = [tparams{end}; {staircase}];
    end
    
    task{1}.private.staircaseTable = table(tparams{:},'VariableNames', tpnames); % save staircase to stimulus

elseif strcmp(param.presSched, 'fixed')
    task{1}.parameter.posDiff            = params.task.posDiff; 
else
    task{1}.parameter.posDiff            = params.task.posDiff; 
end

% privates

task{1}.private.presSched           = param.presSched;
if isfield(params.staircase, 'threshstd_thresh') 
    task{1}.private.threshstd_thresh    = params.staircase.threshstd_thresh;
end

% trial numbers
if exp.block_design  % with some overflow.
    task{1}.numBlocks        = params.trialpercond; % dont count stimright as condition %with some overflow
else
    task{1}.numTrials        = params.numTrials; % dont count stimright as condition %with some overflow
end

% calculated variables
trialdur = task{1}.segmax(2) + task{1}.segmax(3)+max(params.task.stimDur);
% if mglIsFile(exp.noise_mask)
%     trialdur = trialdur + max(params.mask_TOff2MOn) + max(params.maskDur);
% end
maxframes = ceil(trialdur)*myscreen.framesPerSecond+10; % with some additional overflow

task{1}.randVars.calculated.subjcorrect  = nan; 
task{1}.randVars.calculated.stimDur0     = nan; 

if ~isfield(params.task, 'polarAngle') && ~isfield(params.task, 'displAngle')
    task{1}.randVars.calculated.polarAngle   = 0;
    task{1}.randVars.calculated.displAngle   = 0;
end

task{1}.randVars.calculated.bgpermute    = nan(1,maxframes); % nframes x 1 for the background
task{1}.randVars.calculated.stimON       = nan(1,maxframes); % nframes x 1 for the stimulus

task{1}.randVars.calculated.segTime      = nan(1,length(task{1}.segmin));

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
global stimulus;


function [task, myscreen] = initTrial_post_wait(task,myscreen)
global stimulus;

% start the task.
stimulus.lum        = task.thistrial.stimLum;
stimulus.std        = task.thistrial.stimStd;
stimulus.color      = task.thistrial.stimColor;
stimulus.backLum    = task.thistrial.backLum;
stimulus.noiseLum   = task.thistrial.noiseLum;

task.thistrial.framecount   = 0;
task.thistrial.seglen(4)    = task.thistrial.stimDur;
task.thistrial.stimDur0     = task.thistrial.seglen(4);

stimulus.reference  = struct();
stimulus.target     = trackposInitStimulus(stimulus,myscreen); %centerX,Y, diameter called by getArgs.

if stimulus.exp.phasescrambleOn == 1 && stimulus.exp.backprecompute == 1&& stimulus.noiseLum;
    nframes = length(task.thistrial.bgpermute);
    task.thistrial.bgpermute(1:nframes) = randi(length(stimulus.backnoise),nframes,1);
end

if strcmp(task.private.presSched, 'staircase')
    idx         = findCondIdx(task.private.staircaseTable,task.thistrial);
    stimulus.staircaseIdx = idx;
    [s, task.private.staircaseTable.staircase{idx}] = ...
        doStaircase('gettestvalue',task.private.staircaseTable.staircase{idx});

     % check for out of bounds
    if s > (myscreen.imageWidth/2 - stimStd - task.thistrial.pointerOffset)
        s = max(0.1, (myscreen.imageWidth/2 - 3 * stimStd - task.thistrial.pointerOffset));
    end

    task.private.staircaseTable.staircase{idx}.lastTestValue = s;

    task.thistrial.posDiff = s;
    
    if isfield(task.private, 'threshstd_thresh')
        % todo: check quest sd
        if QuestSd(task.private.staircaseTable.staircase{idx}) < task.private.threshstd_thresh
            stimulus.skip = true;
        else
            stimulus.skip = false;
        end
    end
end

% noise mask
if mglIsFile(stimulus.exp.noise_mask)     
    task.thistrial.seglen(5) = task.thistrial.mask_TOff2MOn;
    task.thistrial.seglen(6) = task.thistrial.maskDur;
    maskLum = task.thistrial.maskLum;

    nframes = myscreen.framesPerSecond*task.thistrial.seglen(6) + 20; %/downsample_timeRes; 
    stimulus.noise_mask_trial = randi(size(stimulus.noise_mask.backgroundnoise_rgb,4),nframes,1); % sample with replacement

    % delete texture
    if isfield(stimulus,'noise_mask_texture') 
       for idx = 1:length(stimulus.noise_mask_texture)
           mglDeleteTexture(stimulus.noise_mask_texture{idx});
       end
    end

    % create texture
    for idx = 1:nframes
        midx    = stimulus.noise_mask_trial(idx);
        maskimg = stimulus.noise_mask.backgroundnoise_rgb(:,:,:,midx);
        maskimg(4,:,:) = maskLum * maskimg(4,:,:);
        % maskimg(4,:,:) = 0 * maskimg(4,:,:);
        maskimg = permute(maskimg,[2,3,1]);
        stimulus.noise_mask_texture{idx} = mglMetalCreateTexture(maskimg);
    end
end

% angleSets
if isfield(task.thistrial, 'angleSet') 
    [polarAngle, displAngle] = angleSet(task.thistrial.angleSet);
    task.thistrial.polarAngle = polarAngle;
    task.thistrial.displAngle = displAngle;
end    

if isfield(task.thistrial, 'polarAngle') && isfield(task.thistrial, 'displAngle')
    pa = task.thistrial.polarAngle;
    da = task.thistrial.displAngle;
    r0 = task.thistrial.pointerOffset;
    dr = (2*task.thistrial.stimright-1)*task.thistrial.posDiff;
    
    x0 = r0 * cos(pa);
    y0 = r0 * sin(pa);

    x = x0 + dr*cos(da);
    y = y0 + dr*sin(da);
    
    stimulus.target.position = [x,y];
    stimulus.reference.position = [x0, y0];
else           
    stimulus.target.position = [stim_pos 0];
    stimulus.reference.position = [task.thistrial.pointerOffset, 0];
end



function [task, myscreen] = startSegmentCallback(task, myscreen)
global stimulus

% set flushMode based on noiseLum
if stimulus.exp.phasescrambleOn == 1 && stimulus.exp.backprecompute == 1 && task.thistrial.noiseLum >0
    myscreen.flushMode = 0;
else
    myscreen.flushMode = 0; % default 0 frames
end
 
if task.thistrial.thisseg == 1
    % if we come back to beginning, means the previous trial finished. 
    % wait for the next task 
    myscreen.flushMode = 0; % update screen
    stimulus.currtask   = 'done';
    stimulus.skip       = false;
elseif task.thistrial.thisseg == 2
    if stimulus.exp.trackEye, myscreen.flushMode = 0; end
    if ~stimulus.exp.showmouse, mglDisplayCursor(0);, end 
    
    [task, myscreen] = initTrial_post_wait(task,myscreen);    
    
elseif task.thistrial.thisseg == 4
    if stimulus.exp.trackEye, myscreen.flushMode = 0; end

elseif task.thistrial.thisseg == 5
    task.thistrial.framecount = 0; % restart framecount
    task.thistrial.stimDur = task.thistrial.seglen(4); 
    % stimulus length recorded by updateTask; make sure is same as Dur0 % disp(['Segment duration error1: ', num2str(task.thistrial.stimDur - task.thistrial.stimDur0)])

elseif task.thistrial.thisseg == 6
    task.thistrial.framecount = 0; % restart framecount
    myscreen.flushMode = 0; % refresh every frame
end

% blt screen once before screenUpdates loops
if task.thistrial.thisseg > 1 && task.thistrial.seglen(task.thistrial.thisseg) > 0
    [task, myscreen] = screenUpdateCallback(task, myscreen);
    mglFlush;
    if task.thistrial.thisseg == 4
        stimulus.start = mglGetSecs;
    elseif task.thistrial.thisseg == 5
        stimulus.length = mglGetSecs - stimulus.start;
        disp(['Segment duration error: ', num2str(stimulus.length - task.thistrial.stimDur)])
    end
end

% task.thistrial.segTime(task.thistrial.thisseg) = mglGetSecs;
% disp(['segment' num2str(task.thistrial.thisseg)])
% task.thistrial.seglen
% disp('continuing')

%% screen update
function [task, myscreen] = screenUpdateCallback(task, myscreen)

global stimulus % call stimulus

if task.thistrial.thisseg== 1 % waiting for task to start
    if strcmp(stimulus.currtask,'2afc') % start task
        stimulus.currtask = 'running 2afc';
        task = jumpSegment(task); 
        stimulus.skip = false;
    end
elseif stimulus.skip 
    task = jumpSegment(task);
else
    %% do the task
    % set background luminance
    if task.thistrial.backLum > 1
        mglClearScreen(task.thistrial.backLum/255);
    else
        mglClearScreen(task.thistrial.backLum);
    end
    
    task.thistrial.framecount = task.thistrial.framecount + 1;
    task.thistrial.stimON(task.thistrial.framecount) = 0; %count stimulus
    
    % inject noise, track time
    if any(task.thistrial.thisseg == [2, 3, 4]) 
        if stimulus.exp.phasescrambleOn == 1 
            idx = task.thistrial.bgpermute(task.thistrial.framecount);
            mglBltTexture(stimulus.backnoise{idx},...
                [0 0 myscreen.imageWidth myscreen.imageHeight])
        end
        task.thistrial.trackTime(task.thistrial.framecount) = mglGetSecs(stimulus.t0);
    end

    % draw blob, mask
    if task.thistrial.thisseg == 4 % stimulus
        task.thistrial.stimON(task.thistrial.framecount) = 1;
        mglMetalBltTexture(stimulus.target.img, stimulus.target.position);
        
    elseif task.thistrial.thisseg == 6 % mask
         mglMetalBltTexture(stimulus.noise_mask_texture{task.thistrial.framecount},...
             stimulus.reference.position);
    end

    % reference/fixation helper
    if task.thistrial.thisseg == 2
        mglMetalArcs([stimulus.reference.position, 0], [stimulus.fixColors.stim'; 0.3], [0.4;0.7],[0;2*pi], 1);
        mglMetalArcs([0;0;0], [stimulus.fixColors.afc'; 0.3], [0.2;0.4],[0;2*pi], 1);
    end

    % add fixation
    if any(task.thistrial.thisseg == [2, 3,4])
        if stimulus.exp.colorfix
            % changing fixation colors
            % mglGluAnnulus(0,0,0.2,0.3,stimulus.fixColor,60,1);
            % mglGluDisk(0,0,0.1,0.5+0.5*rand(1,3),60,1);
            mglMetalDots([0;0;0], [0.5+0.5*rand(3,1);1], [stimulus.pointerR; stimulus.pointerR], 1, 1);
        else
            % mglGluDisk(0, 0, 0.1, stimulus.fixColors.stim,60,1); 
            % white fixation
            mglMetalDots([0;0;0], [stimulus.fixColors.afc';1], [stimulus.pointerR; stimulus.pointerR], 1, 1);
        end
        
        % reference
        mglMetalDots([stimulus.reference.position, 0], [stimulus.fixColors.stim';1], ...
            [stimulus.pointerR; stimulus.pointerR], 1, 1);
    elseif any(task.thistrial.thisseg == [5,6,7,8])
        % afc response fixation 
        mglMetalDots([0;0;0], [stimulus.fixColors.afc';1], [0.1;0.1], 1, 1);

        % reference
        mglMetalDots([stimulus.reference.position, 0], [stimulus.fixColors.stim';1], ...
            [stimulus.pointerR; stimulus.pointerR], 1, 1);
    end
    
    % response direction arrow
    if isfield(stimulus.exp, 'respDirArrow') && stimulus.exp.respDirArrow
        if isfield(task.thistrial, 'polarAngle') && isfield(task.thistrial, 'displAngle')
            pa = task.thistrial.polarAngle;
            da = task.thistrial.displAngle;
        else
            pa = 0;
            da = 0;
        end
        r0 = task.thistrial.pointerOffset;

        dr = 1; % length of arrow
        arm_ratio = 1/3;
        arm_angle = pi/6;
        
        x0 = r0 * cos(pa);
        y0 = r0 * sin(pa);
        
        mglMetalArrow(x0,y0,da,dr,arm_ratio, arm_angle, 0.1, stimulus.fixColors.stim);
    end

    % feedback
    if stimulus.exp.feedback && task.thistrial.thisseg == 8
        % no fixation cross until response.
        % mglGluAnnulus(0,0,0.2,0.5,stimulus.currfixColor,60,1);
        mglMetalArcs([stimulus.reference.position, 0], [stimulus.currfixColor'; 1], [0.3;0.5],[0;2*pi], 1);
        
        % feedback about presented position
        if isfield(params, 'feedback_center') && params.feedback_center
            mglMetalDots([stimulus.target.position,0], [stimulus.fixColors.fb';1], [stimulus.pointerR; stimulus.pointerR], 1, 1);
        end
    end
        
    
    % track eye
    if stimulus.exp.trackEye && any(task.thistrial.thisseg==[2,3,4]) 
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

if strcmp(task.private.presSched, 'staircase') && ~isnan(correct)
    idx = stimulus.staircaseIdx;
    task.private.staircaseTable.staircase{idx} = ...
        doStaircase('update',task.private.staircaseTable.staircase{idx},correct);
end

task = jumpSegment(task); % go to next segment


%% predefined angleSet
function [polarAngle, displAngle] = angleSet(angleSet)
    % set polar angle based on remainder
    polarAngle = angleCode2angle(mod(angleSet, 100), 4);
        
    if floor(angleSet/100) == 0 % tangential motion
        displAngle = mod(angleCode2angle(mod(angleSet, 100), 4)+pi/2,2*pi);
    end
    
    
function angle = angleCode2angle(code, div)
    % code: natural number
    angle = mod((code-1) * 2*pi/div, 2*pi);
end

    
