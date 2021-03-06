% attentionComp.m
%
%      usage: myscreen = attentionComp()
%         by: minyoung lee
%       date: 07/01/2016
%    purpose: biased competition experiment
%
%       task: Vertical/CW from vertical

function myscreen = attentionComp(varargin)

clear global stimulus
global stimulus

sensory = 0;
debug = 0;
% evaluate the arguments
getArgs(varargin);
stimulus.sensory = sensory;

% Setup parameters for stimulus
stimulus.width = 6; % gabor
stimulus.contrastLow = .2;
stimulus.contrastHigh = .8;
stimulus.sf = 1.8;
stimulus.orientation = 90;
% angles  (1->4 clockwise from top)
% loc1: 60  loc2: 30  loc3: 240 loc4: 210
stimulus.angles = [65 25 245 205] * pi /180;
stimulus.locations = 1:length(stimulus.angles);
stimulus.radius = 6; % eccentricity
stimulus.upperLoc = [1 2];
stimulus.lowerLoc = [3 4];

% fixation cross and cue
stimulus.fixWidth = 1.2;
stimulus.fixColor = [1 1 1];
stimulus.cueColor = [1 1 1];

stimulus.interval = [2 4];

%set parameters of staircase
stimulus.initialThreshold = [3 3]; %for low and high separately
stimulus.stepsize = .05;
stimulus.minStepsize = 0.005;
stimulus.minThreshold = 0;
stimulus.maxThreshold = 45;

% initalize the screen
myscreen.background = 0.5;
myscreen = initScreen(myscreen);

% if this is sensory expt, initiate the fixation task
if stimulus.sensory
global fixStimulus
fixStimulus.fixWidth = stimulus.fixWidth;
% setting the fixation task to task 2
[task{2} myscreen] = fixStairInitTask(myscreen);
end

%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%%

% task{1}.waitForBacktick = 0;
task{1}{1}.waitForBacktick = 1;
% trial: 4.2s + ITI(1~14s)
task{1}{1}.segmin = [1 0.7 0.3 0.7 1.5 1];
task{1}{1}.segmax = [1 0.7 0.3 0.7 1.5 14];
if ~stimulus.sensory
    task{1}{1}.getResponse = [0 0 0 0 1 0];
else
    task{1}{1}.getResponse = [0 0 0 0 0 0];
end
task{1}{1}.synchToVol = [0 0 0 0 0 1];
task{1}{1}.random = 1;
if debug
    task{1}{1}.waitForBacktick = 0;
    task{1}{1}.synchToVol(end) = 0;
    task{1}{1}.segmax(end) = 1;
end

% parameters & randomization
if ~stimulus.sensory
task{1}{1}.parameter.targetContrast = [1 2]; %1 low /2 high
task{1}{1}.randVars.uniform.targetLoc = stimulus.locations;
% task{1}.randVars.uniform.targetContrast = [1 2]; %low/high
else
    task{1}{1}.randVars.block.unattendedStim2 = [-1 1 2];
    task{1}{1}.randVars.block.unattendedStimPos2 = [1 2];    
end
task{1}{1}.randVars.block.unattendedStim = [-1 1 2]; %low alone/ high alone/ both
task{1}{1}.randVars.block.unattendedStimPos = [1 2]; %1:low-high 2:high-low

task{1}{1}.randVars.block.orientation1 = [1 2];
task{1}{1}.randVars.block.orientation2 = [1 2];

% calculated
if ~stimulus.sensory
task{1}{1}.randVars.calculated.attendContrast = nan;
task{1}{1}.randVars.calculated.attendAway = nan;
task{1}{1}.randVars.calculated.attendedPair= {nan};
task{1}{1}.randVars.calculated.unattendedPair = {nan};
task{1}{1}.randVars.calculated.correct = nan;
task{1}{1}.randVars.calculated.unattendedQuad = nan;
task{1}{1}.randVars.calculated.attendedQuad = nan;
end

% initialize the task
for phaseNum = 1:length(task{1})
    if ~stimulus.sensory
  [task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
    else
        [task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback);
    end
end

% init the stimulus
myscreen = initStimulus('stimulus',myscreen);

% to initialize the stimulus for your experiment.
stimulus = initGabor(stimulus,myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
  % update the task
  [task{1} myscreen phaseNum] = updateTask(task{1},myscreen,phaseNum);
  if stimulus.sensory
    % update the fixation task
    [task{2} myscreen] = updateTask(task{2},myscreen,1);
  end
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)
global stimulus
if task.thistrial.thisseg == 1
    % decide which stim interval will contain tilted garbor at each of 4 locations 
    task.thistrial.whichInterval = round(rand(1,length(stimulus.locations))) + 1;%(round(rand(1,length(stimulus.locations))) + 1)*2;%%
    stimulus.cueColor = [1 1 1];
    
    % orientation (vertical/horizontal) for loc 1-4
    if task.thistrial.orientation1 == 1 % vertical in loc 1 & horizontal in loc 2
        task.thistrial.orient(1) = 0;
        task.thistrial.orient(2) = -90;
    else % hor in loc 1 & ver in loc 2
        task.thistrial.orient(1) = -90;
        task.thistrial.orient(2) = 0;
    end
    if task.thistrial.orientation2 == 1
        task.thistrial.orient(3) = 0;
        task.thistrial.orient(4) = -90;
    else
        task.thistrial.orient(3) = -90;
        task.thistrial.orient(4) = 0;
    end
    
    % Unattended stimuli configs...
    if task.thistrial.unattendedStimPos == 1
        low(1) = 1; high(1) =2;
    else
        low(1) = 2; high(1) = 1;
    end
    
    if ~stimulus.sensory   
    if any(task.thistrial.targetLoc == stimulus.upperLoc)
        attendedPair = stimulus.upperLoc;
        unattendedPair = stimulus.lowerLoc;
        task.thistrial.attendedQuad = 1;
        task.thistrial.unattendedQuad = 2;
    else
        attendedPair = stimulus.lowerLoc;
        unattendedPair = stimulus.upperLoc;
        task.thistrial.attendedQuad = 2;
        task.thistrial.unattendedQuad = 1;
    end
        
    % Attended side    
    if task.thistrial.targetContrast == 1 % Attend to LOW
        task.thistrial.tex{task.thistrial.targetLoc} = stimulus.tex(1);
        task.thistrial.tex{attendedPair(task.thistrial.targetLoc ~= attendedPair)} = stimulus.tex(2);
        task.thistrial.contrast{task.thistrial.targetLoc} = 'low';
        task.thistrial.contrast{attendedPair(task.thistrial.targetLoc ~= attendedPair)} = 'high';
        task.thistrial.attendContrast = 1;
    else   % Attend to HIGH
        task.thistrial.tex{task.thistrial.targetLoc} = stimulus.tex(2);
        task.thistrial.tex{attendedPair(task.thistrial.targetLoc ~= attendedPair)} = stimulus.tex(1);
        task.thistrial.contrast{task.thistrial.targetLoc} = 'high';
        task.thistrial.contrast{attendedPair(task.thistrial.targetLoc ~= attendedPair)} = 'low';
        task.thistrial.attendContrast = 2;
    end
    
    % Unattended stimuli configs...
     task.thistrial.tex{unattendedPair(low)} = stimulus.tex(1);
     task.thistrial.tex{unattendedPair(high)} = stimulus.tex(2);
     task.thistrial.contrast{unattendedPair(low)} = 'low';
     task.thistrial.contrast{unattendedPair(high)} = 'high';
    if task.thistrial.unattendedStim == -1
        task.thistrial.tex{unattendedPair(high)} = nan;
        task.thistrial.contrast{unattendedPair(high)} = nan;
        task.thistrial.attendAway = -1;
    elseif task.thistrial.unattendedStim == 1
        task.thistrial.tex{unattendedPair(low)} = nan;
        task.thistrial.contrast{unattendedPair(low)} = nan;
        task.thistrial.attendAway = 1;
    else
        task.thistrial.attendAway = 2;
    end

    % save
    task.thistrial.attendedPair = attendedPair;
    task.thistrial.unattendedPair = unattendedPair;
    
    disp(sprintf('(attentionComp) Trial %i: Location %i ITI %0.2f s', task.trialnum, task.thistrial.targetLoc, task.thistrial.seglen(end)));
    
    
    elseif stimulus.sensory
        if task.thistrial.unattendedStimPos2 == 1
            low(2) = 1; high(2) = 2;
        else
            low(2) = 2; high(2) = 1;
        end
        Quad = [1 2; 3 4];
        for q = 1:2
            task.thistrial.tex{Quad(q,low(q))} = stimulus.tex(1);
            task.thistrial.tex{Quad(q,high(q))} = stimulus.tex(2);
            
            task.thistrial.contrast{Quad(q,low(q))} = 'low';
            task.thistrial.contrast{Quad(q,high(q))} = 'high';
            if q == 1
                if task.thistrial.unattendedStim == -1
                    task.thistrial.tex{Quad(q,high(q))} = nan;
                    task.thistrial.contrast{Quad(q,high(q))} = nan;
                elseif task.thistrial.unattendedStim == 1
                    task.thistrial.tex{Quad(q,low(q))} = nan;
                    task.thistrial.contrast{Quad(q,low(q))} = nan;
                end
            else
                if task.thistrial.unattendedStim2 == -1
                    task.thistrial.tex{Quad(q,high(q))} = nan;
                    task.thistrial.contrast{Quad(q,high(q))} = nan;
                elseif task.thistrial.unattendedStim2 == 1
                    task.thistrial.tex{Quad(q,low(q))} = nan;
                    task.thistrial.contrast{Quad(q,low(q))} = nan;
                end
            end
                    
                    
        end

        disp(sprintf('(attentionComp:Sensory) Trial %i: %s %s %s %s, ITI %0.2f s', task.trialnum, task.thistrial.contrast{1}, ...
            task.thistrial.contrast{2}, task.thistrial.contrast{3}, task.thistrial.contrast{4}, task.thistrial.seglen(end)));
        
    end

elseif any(task.thistrial.thisseg == stimulus.interval)

    for loc = 1:length(stimulus.locations)
        if ~isnan(task.thistrial.contrast{loc})
            if ~stimulus.sensory
            % which staircase? (target low or high)
            if strcmp(task.thistrial.contrast(loc), 'low')
                stairNum = 1;
            else
                stairNum = 2; 
            end
            if task.thistrial.thisseg == task.thistrial.whichInterval(loc)*2 % Segment 2 or 4
                task.thistrial.thisRotation(loc) = task.thistrial.orient(loc)-stimulus.stair{stairNum}.threshold;
            else
                task.thistrial.thisRotation(loc) = task.thistrial.orient(loc);
            end
            else
                if task.thistrial.thisseg == task.thistrial.whichInterval(loc)*2
                    task.thistrial.thisRotation(loc) = task.thistrial.orient(loc)-.5;
                else
                    task.thistrial.thisRotation(loc) = task.thistrial.orient(loc);
                end
            end
        else
            task.thistrial.thisRotation(loc) = nan;
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)
global stimulus
mglClearScreen;

% draw stim, 2nd/4th seg
if any(task.thistrial.thisseg == stimulus.interval)  
    for i = 1:length(stimulus.locations)
        if ~isnan(task.thistrial.contrast{i})
            mglBltTexture(task.thistrial.tex{i}, [stimulus.x(i) stimulus.y(i)], 0,0, task.thistrial.thisRotation(i));
        end
    end
end

if ~stimulus.sensory
% draw the cue
if task.thistrial.thisseg ~= 6    
    mglLines2(stimulus.cueLines.x0(task.thistrial.targetLoc), stimulus.cueLines.y0(task.thistrial.targetLoc),...
        stimulus.cueLines.x1(task.thistrial.targetLoc), stimulus.cueLines.y1(task.thistrial.targetLoc), 1.5, stimulus.cueColor);    
end

%draw fixation cross
mglFixationCross(stimulus.fixWidth,1.5,stimulus.fixColor);
end

% for i = 1:length(stimulus.locations)
%      mglLines2(stimulus.refLines.x1{i},stimulus.refLines.y1{i},stimulus.refLines.x2{i},stimulus.refLines.y2{i},1,[1 1 1]);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)

global stimulus

if ~task.thistrial.gotResponse
if task.thistrial.whichInterval(task.thistrial.targetLoc) == task.thistrial.whichButton
    % correct
    task.thistrial.correct = true;
    % feedback
    stimulus.cueColor = [0 1 0];
    %update staircase
    stimulus.stair{task.thistrial.targetContrast} = upDownStaircase(stimulus.stair{task.thistrial.targetContrast},1);
    task.thistrial.correct = 1;
    disp(sprintf('(attentionComp) Trial %i: %s contrast %0.4f correct', task.trialnum, task.thistrial.contrast{task.thistrial.targetLoc},stimulus.stair{task.thistrial.targetContrast}.threshold));
else
    % incorrect
    task.thistrial.correct = false;
    stimulus.cueColor = [1 0 0];
    %update staircase
    stimulus.stair{task.thistrial.targetContrast} = upDownStaircase(stimulus.stair{task.thistrial.targetContrast},0);
    task.thistrial.correct = 0;
    disp(sprintf('(attentionComp) Trial %i: %s contrast %0.4f incorrect', task.trialnum, task.thistrial.contrast{task.thistrial.targetLoc}, stimulus.stair{task.thistrial.targetContrast}.threshold));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initGabor(stimulus,myscreen)

% compute the grating
grating = mglMakeGrating(stimulus.width, stimulus.width, stimulus.sf, stimulus.orientation,0);
gaussian = mglMakeGaussian(stimulus.width, stimulus.width,stimulus.width/10, stimulus.width/10); %sd?????

gaborLow(:,:,1) = 255*(grating+1)/2;
gaborLow(:,:,2) = 255*(grating+1)/2;
gaborLow(:,:,3) = 255*(grating+1)/2;
gaborLow(:,:,4) =255*(stimulus.contrastLow*gaussian);

gaborHigh(:,:,1) = 255*(grating+1)/2;
gaborHigh(:,:,2) = 255*(grating+1)/2;
gaborHigh(:,:,3) = 255*(grating+1)/2;
gaborHigh(:,:,4) =255*(stimulus.contrastHigh*gaussian);

%create texture
stimulus.tex(1) = mglCreateTexture(gaborLow);
stimulus.tex(2) = mglCreateTexture(gaborHigh);

if stimulus.sensory
    for i = 1:length(stimulus.angles)
        %stim centers
      [stimulus.x(i), stimulus.y(i)] = pol2cart(stimulus.angles(i),stimulus.radius);
    end
        
else
%cue lines / stimulus centers
for i = 1:length(stimulus.angles)
      % starting point of lines
      [stimulus.cueLines.x0(i), stimulus.cueLines.y0(i)] = pol2cart(stimulus.angles(i), 0.05);
      % ending point
      [stimulus.cueLines.x1(i), stimulus.cueLines.y1(i)] = pol2cart(stimulus.angles(i), stimulus.fixWidth*0.6);
      
      %stim centers
      [stimulus.x(i), stimulus.y(i)] = pol2cart(stimulus.angles(i),stimulus.radius);
      
%         % circular reference lines
%     stimulus.refLines.x1{i} = [];
%     stimulus.refLines.y1{i} = [];
%     stimulus.refLines.x2{i} = [];
%     stimulus.refLines.y2{i} = [];
%   d = 0:1:360;
%   for dIndex = 1:length(d)-1
%     stimulus.refLines.x1{i} = [stimulus.refLines.x1{i} stimulus.x(i)+stimulus.width/2*cos(d2r(d(dIndex)))];
%     stimulus.refLines.y1{i} = [stimulus.refLines.y1{i} stimulus.y(i)+stimulus.width/2*sin(d2r(d(dIndex)))];
%     stimulus.refLines.x2{i} = [stimulus.refLines.x2{i} stimulus.x(i)+stimulus.width/2*cos(d2r(d(dIndex+1)))];
%     stimulus.refLines.y2{i} = [stimulus.refLines.y2{i} stimulus.y(i)+stimulus.width/2*sin(d2r(d(dIndex+1)))];
%   end

end

% see if there was a previous staircase
if ~isempty(mglGetSID) && isdir(sprintf('~/data/attentionComp/%s', mglGetSID))
    files = dir(sprintf('~/data/attentionComp/%s/1*', mglGetSID));
    if ~isempty(files)
        s = load(sprintf('~/data/attentionComp/%s/%s', mglGetSID,files(end).name));
    end
else
    s = getLastStimfile(myscreen);
end

if exist('s') && isfield(s, 'stimulus') && isfield(s.stimulus,'stair')
    if ~strcmp(s.myscreen.computerShortname, myscreen.computerShortname)
        % set starting thershold to staircase value
        stimulus.initialThreshold(1) = s.stimulus.stair{1}.threshold;
        stimulus.initialThreshold(2) = s.stimulus.stair{2}.threshold;
        
        % now create staircases for each orientation
        for i = 1:2
     % init a 2 down 1 up staircase
         stimulus.stair{i} = upDownStaircase(1,2,stimulus.initialThreshold(i),[stimulus.stepsize, stimulus.minStepsize, stimulus.stepsize], 'pest');
         stimulus.stair{i}.minThreshold = stimulus.minThreshold;
         stimulus.stair{i}.maxThreshold = stimulus.maxThreshold;
        end

        % display what we are doing
        disp(sprintf('(attentionComp) Setting starting thershold based on last stimfile to: %f %f',stimulus.initialThreshold(1), stimulus.initialThreshold(2)));
    else
        % set starting thershold to staircase value
%         stimulus.initialThreshold(1) = s.stimulus.stair{1}.threshold;
%         stimulus.initialThreshold(2) = s.stimulus.stair{2}.threshold;
        stimulus.stair{1} = s.stimulus.stair{1};
        stimulus.stair{2} = s.stimulus.stair{2};
        % display what we are doing
        disp(sprintf('(attentionComp) Setting starting thershold based on last stimfile to: %f %f',stimulus.initialThreshold(1), stimulus.initialThreshold(2)));

    end
else
        % now create staircases for each orientation
        for i = 1:2
     % init a 2 down 1 up staircase
         stimulus.stair{i} = upDownStaircase(1,2,stimulus.initialThreshold(i),[stimulus.stepsize, stimulus.minStepsize, stimulus.stepsize], 'pest');
         stimulus.stair{i}.minThreshold = stimulus.minThreshold;
         stimulus.stair{i}.maxThreshold = stimulus.maxThreshold;
        end
end
clear s;
end
