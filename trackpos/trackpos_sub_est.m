%
%        $Id: $
%      usage: trackpos_est
%         by: Josh Ryu
%       date: 05/06/2021
%    purpose: estimate position of blob
%
% staircase not implemented yet
% consider: adding another noise period after the stimulus

% S1: wait time while the main task calls this subtask.
% S2: random period of fixation (random ~0.5s)
% S3: stimulus period (stimdur s)
% S4: repsonse period (inf)
% S5: feedback (1s)

% fixed position difference or samples from continuous distributions?

function [task, myscreen] = trackpos_sub_est(myscreen, params, exp)
%% task parameters
% stimulus and background
task{1}.random               = 1;
task{1}.parameter.backLum    = params.backLum; 
task{1}.parameter.noiseLum   = params.noiseLum;
task{1}.parameter.stimright  = [0,1];
task{1}.parameter.posDiff    = params.posDiff; 
task{1}.parameter.stimLum 	 = params.stimLum;

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

task{1}.randVars.calculated.stimDur      = nan; 
task{1}.randVars.calculated.pos_estim    = nan(1,2); % nframes x 2 for position estimate
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
    % start the task.
    if ~stimulus.exp.showmouse, mglDisplayCursor(0);, end 

    
    stimulus.stimLum    = task.thistrial.stimLum;
    stimulus.backLum    = task.thistrial.backLum;
    stimulus.noiseLum   = task.thistrial.noiseLum;

    task.thistrial.framecount = 0;
    task.thistrial.stimDur    = task.thistrial.seglen(3); % sstimulus period

    stimulus = trackposInitStimulus(stimulus,myscreen); %centerX,Y, diameter called by getArgs.

    if stimulus.exp.phasescrambleOn == 1 && stimulus.exp.backprecompute == 1;
        nframes = length(task.thistrial.bgpermute);
        task.thistrial.bgpermute(1:nframes) = randi(length(stimulus.backnoise),nframes,1);
    end
    
elseif task.thistrial.thisseg == 4
    % response period
    % center the mouse 
    x_img = 0;  y_img = 0;
    x_screen = x_img*myscreen.screenWidth/myscreen.imageWidth + myscreen.screenWidth/2;
    y_screen = y_img*myscreen.screenHeight/myscreen.imageHeight + myscreen.screenHeight/2;
    mglSetMousePosition(x_screen,y_screen, myscreen.screenNumber);
end

%% screen update
function [task, myscreen] = screenUpdateCallback(task, myscreen)
% S1: wait time while the main task calls this subtask.
% S2: random period of fixation (random ~0.5s)
% S3: stimulus period (stimdur s)
% S4: repsonse period (inf)
% S5: feedback (1s)

global stimulus  % call stimulus

if task.thistrial.thisseg== 1
%% waiting for task to start
if strcmp(stimulus.currtask,'est')
    stimulus.currtask = 'running est';
    task = jumpSegment(task); 
end

else
%% do the task
mglClearScreen(stimulus.backLum/255);

task.thistrial.framecount = task.thistrial.framecount + 1;
task.thistrial.stimON(task.thistrial.framecount) = 0; %count stimulus

% inject noise, track time, add fixation
if any(task.thistrial.thisseg == [2,3]) 
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
    mglGluDisk(0, 0, 0.1, stimulus.fixColors.est,60,1); 
end

% draw blob or response feedback
if task.thistrial.thisseg == 3 % stimulus
    stim_pos = (2*task.thistrial.stimright-1)*task.thistrial.posDiff;
    task.thistrial.stimON(task.thistrial.framecount) = 1;
    mglBltTexture(stimulus.gaussian,[stim_pos 0]);
    
elseif task.thistrial.thisseg == 4 %response period; 
    % show cursor
    mInfo = mglGetMouse(myscreen.screenNumber);
    degx = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
    degy = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;
    
    if ~stimulus.exp.estim_horiz, degx = 0; end
    if ~stimulus.exp.estim_verti, degy = 0; end
    mglGluDisk(degx, degy, 0.1, [1 0 0]);
    
elseif task.thistrial.thisseg == 5 % feedback period
    stim_pos = (2*task.thistrial.stimright-1)*task.thistrial.posDiff;
    mglGluDisk(stim_pos, 0, 0.1, [1 0 0]) ;    % draw center of blob
    
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

% if the button 3 is pressed, record position
if task.thistrial.whichButton == 3
    % record position of the mouse
    mInfo = mglGetMouse(myscreen.screenNumber);
    degx = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
    degy = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;
    
    if ~stimulus.exp.estim_horiz, degx = 0; end
    if ~stimulus.exp.estim_verti, degy = 0; end
    task.thistrial.pos_estim = [degx, degy]; % save position estimates
    
    stim_pos = task.thistrial.posDiff;

    disp(['(subtask_est) Stimulus Hor. Pos.: ' num2str(stim_pos) '; ' ...
            'Response: ' num2str(degx)  ...
            '; Error: ' num2str(degx - stim_pos)]);
    
    % go to next segment
    task = jumpSegment(task);
end
