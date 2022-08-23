%        $Id: $
%      usage: trackpos_2afc
%         by: Josh Ryu
%       date: 04/13/2021
%    purpose: 2 choice task for position of one blob against fixation.

% S1: wait time while the main task calls this subtask.
% S2: random period of fixation (random ~0.5s)
% S3: stimulus period (stimdur s)
% S4: repsonse period (inf)
% S5: feedback (1s)

function [task, myscreen] = trackpos_sub_2c(myscreen, params, exp)
% stimulus and background
task{1}.random               = 1;
task{1}.parameter.backLum    = params.backLum; 
task{1}.parameter.noiseLum   = params.noiseLum;
task{1}.parameter.stimright  = [0,1];
task{1}.parameter.posDiff    = params.posDiff; 
task{1}.parameter.stimLum 	 = params.stimLum;
task{1}.parameter.stimStd 	 = params.stimStd;

% note: seglen are changed later
if ~exp.feedback 
    task{1}.segmin           = [inf 0.4 nan inf];
    task{1}.segmax           = [inf 0.8 nan inf]; 
    task{1}.getResponse      = [0 0 0 1]; %segment to get response.
else
    task{1}.segmin           = [inf 0.4 nan inf 1];
    task{1}.segmax           = [inf 0.8 nan inf 1]; 
    task{1}.getResponse      = [0 0 0 1 0]; %segment to get response.
end

task{1}.numTrials        = (params.numTrials)/2 + 100; % (with some overflow)

% calculated variables
maxframes = ceil((task{1}.segmax(2)+max(params.stimDur))...
    *myscreen.framesPerSecond)+10; % with some additional overflow

task{1}.randVars.calculated.subjcorrect  = nan; 
task{1}.randVars.calculated.stimDur      = nan; 

task{1}.randVars.calculated.bgpermute    = nan(1,maxframes); % nframes x 1 for the background
task{1}.randVars.calculated.stimON       = nan(1,maxframes); % nframes x 1 for the stimulus

task{1}.randVars.calculated.trackTime    = nan(1,maxframes);
task{1}.randVars.calculated.trackEye     = nan(maxframes,2);
task{1}.randVars.calculated.trackEyeTime = nan(1,maxframes); % for referencing edf file

%% task blocks. 

% change stimulus duration
for trialN = 1:task{1}.numTrials
    fixdur      = rand*(task{1}.segmax(2) - task{1}.segmin(2)) + task{1}.segmin(2);
    stimdur     = params.stimDur(randi(length(params.stimDur)));
    if ~exp.feedback
        task{1}.seglenPrecompute.seglen{trialN} = [inf fixdur stimdur inf];
    else
        task{1}.seglenPrecompute.seglen{trialN} = [inf fixdur stimdur inf 1];
    end
end

% initializing task...
disp(' Initializing Task....')

for phaseN = 1:length(task)
    [task{phaseN}, myscreen] = initTask(task{phaseN},myscreen,...
        @startSegmentCallback,@screenUpdateCallback,@responseCallback,[]);
end


%% Start segment
function [task, myscreen] = startSegmentCallback(task, myscreen)
global stimulus

if task.thistrial.thisseg == 1
    % if we come back to beginning, means the previous trial finished. 
    % wait for the next task 
    stimulus.currtask = 'done';
    
elseif task.thistrial.thisseg == 2
    if ~stimulus.exp.showmouse, mglDisplayCursor(0);, end 

    % start the task.
    stimulus.lum        = task.thistrial.stimLum;
    stimulus.std        = task.thistrial.stimStd;
    stimulus.backLum    = task.thistrial.backLum;
    stimulus.noiseLum   = task.thistrial.noiseLum;
    
    task.thistrial.framecount = 0;
    task.thistrial.stimdur    = task.thistrial.seglen(3);
    
    stimulus.target = trackposInitStimulus(stimulus,myscreen); %centerX,Y, diameter called by getArgs.
    
    if stimulus.exp.phasescrambleOn == 1 && stimulus.exp.backprecompute == 1;
        nframes = length(task.thistrial.bgpermute);
        task.thistrial.bgpermute(1:nframes) = randi(length(stimulus.backnoise),nframes,1);
    end

end

%% screen update
function [task, myscreen] = screenUpdateCallback(task, myscreen)
% S1: wait time while the main task calls this subtask.
% S2: random period of fixation (random ~0.5s)
% S3: stimulus period (stimdur s)
% S4: repsonse period (inf)
% S5: feedback (1s)

global stimulus % call stimulus

if task.thistrial.thisseg== 1
%% waiting for task to start
if strcmp(stimulus.currtask,'2c')
    stimulus.currtask = 'running 2c';
    task = jumpSegment(task); 
end

else
%% do the task
mglClearScreen(stimulus.backLum/255);

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
    
elseif any(task.thistrial.thisseg == [4,5])
    % fixation indicating estimation task
    mglGluDisk(0, 0, 0.1, stimulus.fixColors.afc,60,1); 
end

% draw blob or response feedback
if task.thistrial.thisseg == 3 % stimulus
    stim_pos = (2*task.thistrial.stimright-1)*task.thistrial.posDiff;
    task.thistrial.stimON(task.thistrial.framecount) = 1;
    mglBltTexture(stimulus.target.img,[stim_pos 0]);
    
elseif task.thistrial.thisseg == 5 %feedback period
    % no fixation cross until response.
    mglGluAnnulus(0,0,0.2,0.3,stimulus.currfixColor,60,1);
    
    % feedback about presented position
%     stim_pos = (2*task.thistrial.stimright-1)*task.thistrial.posDiff;
%     mglGluDisk(stim_pos, 0, 0.1, [1 0 0]) ;    % draw center of blob

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

task = jumpSegment(task); % go to next segment