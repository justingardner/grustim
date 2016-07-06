function myscreen = attentionComp()

%task
%localizer

% Contrasts
% Low/High 12.5/50%

% Orientations
% Vertical/CW to vertical => staircased

global stimulus
stimulus.width = 1.2; % gabor
stimulus.contrastLow = .125;
stimulus.contrastHigh = .5;
stimulus.sf = 4.8;
stimulus.orientation = 90;
% angles  (1->4 clockwise from top)`
% loc1: pi/3 (60)  loc2: pi/6 (30)  loc3: pi*4/3 (240) loc4: pi*7/6 (210)
stimulus.angles = [pi/3 pi/6 pi*4/3 pi*7/6];
stimulus.locations = 1:length(stimulus.angles);
stimulus.radius = 6;
stimulus.upperLoc = [1 2];
stimulus.lowerLoc = [3 4];

stimulus.fixWidth = 1;
stimulus.fixColor = [1 1 1];
stimulus.cueColor = [1 1 1];

stimulus.interval = [2 4];

%set parameters of staircase
stimulus.initialThreshold = 10; 
stimulus.stepsize = .25;
stimulus.minThreshold = 0;
stimulus.maxThreshold = 45;
% stimulus.lowCon = stair;
% stimulus.highCon = stair;

% initalize the screen
myscreen.background = 0.5;
myscreen = initScreen(myscreen);

%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%%

task{1}.waitForBacktick = 1;
task{1}.segmin = [1 0.7 0.3 0.7 1.5 1.5];
task{1}.segmax = [1 0.7 0.3 0.7 1.5 5];
task{1}.getResponse = [0 0 0 0 1 0];
% task{1}.synchToVol = [0 0 0 0 0 1];
task{1}.random = 1;

% parameters

task{1}.randVars.uniform.targetLoc = stimulus.locations;
task{1}.randVars.uniform.targetContrast = [1 2]; %low/high
task{1}.randVars.block.unattendedStim = [-1 1 2]; %low alone/ high alone/ both
task{1}.randVars.block.unattendedStimPos = [1 2];

% initialize the task
for phaseNum = 1:length(task)
  [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end

% init the stimulus
myscreen = initStimulus('stimulus',myscreen);

% to initialize the stimulus for your experiment.
stimulus = myInitStimulus(stimulus,myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task)) && ~myscreen.userHitEsc
  % update the task
  [task myscreen phaseNum] = updateTask(task,myscreen,phaseNum);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)
global stimulus
if task.thistrial.thisseg == 1
    task.thistrial.whichInterval = round(rand(1,length(stimulus.locations))) + 1;%(round(rand(1,length(stimulus.locations))) + 1)*2;%%
    stimulus.cueColor = [1 1 1];
    
    if any(task.thistrial.targetLoc == stimulus.upperLoc)
        attendedPair = stimulus.upperLoc;
        unattendedPair = stimulus.lowerLoc;
    else
        attendedPair = stimulus.lowerLoc;
        unattendedPair = stimulus.upperLoc;
    end
        
    if task.thistrial.targetContrast == 1
        task.thistrial.tex{task.thistrial.targetLoc} = stimulus.tex(1);
        task.thistrial.tex{attendedPair(task.thistrial.targetLoc ~= attendedPair)} = stimulus.tex(2);
        task.thistrial.contrast{task.thistrial.targetLoc} = 'low';
        task.thistrial.contrast{attendedPair(task.thistrial.targetLoc ~= attendedPair)} = 'high';
    else
        task.thistrial.tex{task.thistrial.targetLoc} = stimulus.tex(2);
        task.thistrial.tex{attendedPair(task.thistrial.targetLoc ~= attendedPair)} = stimulus.tex(1);
        task.thistrial.contrast{task.thistrial.targetLoc} = 'high';
        task.thistrial.contrast{attendedPair(task.thistrial.targetLoc ~= attendedPair)} = 'low';     
    end
        
    if task.thistrial.unattendedStimPos == 1
        low = 1; high =2;
    else
        low = 2; high = 1;
    end
    
     task.thistrial.tex{unattendedPair(low)} = stimulus.tex(1);
     task.thistrial.tex{unattendedPair(high)} = stimulus.tex(2);
     task.thistrial.contrast{unattendedPair(low)} = 'low';
     task.thistrial.contrast{unattendedPair(high)} = 'high';
    if task.thistrial.unattendedStim == -1
        task.thistrial.tex{unattendedPair(high)} = nan;
        task.thistrial.contrast{unattendedPair(high)} = nan;
    elseif task.thistrial.unattendedStim == 1
        task.thistrial.tex{unattendedPair(low)} = nan;
        task.thistrial.contrast{unattendedPair(low)} = nan;
    end

    task.thistrial.attendedPair = attendedPair;
    task.thistrial.unattendedPair = unattendedPair;

elseif any(task.thistrial.thisseg == stimulus.interval)

    for loc = 1:length(stimulus.locations)
        if ~isnan(task.thistrial.contrast{loc})
            if strcmp(task.thistrial.contrast(loc), 'low')
                stairNum = 1;
            else
                stairNum = 2; 
            end
            if task.thistrial.thisseg == task.thistrial.whichInterval(loc)*2
                task.thistrial.thisRotation(loc) = -stimulus.stair{stairNum}.threshold;
            else
                task.thistrial.thisRotation(loc) = 0;
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

%draw fixation cross
mglFixationCross(stimulus.fixWidth,1.5,stimulus.fixColor);

% draw stim, 2nd/4th seg
if any(task.thistrial.thisseg == stimulus.interval)  
    for i = 1:length(stimulus.locations)
        if ~isnan(task.thistrial.contrast{i})
            mglBltTexture(task.thistrial.tex{i}, [stimulus.x(i) stimulus.y(i)], 0,0, task.thistrial.thisRotation(i));
        end
    end
end

if task.thistrial.thisseg ~= 6 % draw the cue    
    mglLines2(stimulus.cueLines.x0(task.thistrial.targetLoc), stimulus.cueLines.y0(task.thistrial.targetLoc),...
        stimulus.cueLines.x1(task.thistrial.targetLoc), stimulus.cueLines.y1(task.thistrial.targetLoc), 1.5, stimulus.cueColor);    
end



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
    disp(sprintf('(ExptName) target: %s contrast %0.4f correct', task.thistrial.contrast{task.thistrial.targetLoc},stimulus.stair{task.thistrial.targetContrast}.threshold));
else
    % incorrect
    task.thistrial.correct = false;
    stimulus.cueColor = [1 0 0];
    %update staircase
    stimulus.stair{task.thistrial.targetContrast} = upDownStaircase(stimulus.stair{task.thistrial.targetContrast},0);
    disp(sprintf('(ExptName) target: %s contrast %0.4f incorrect', task.thistrial.contrast{task.thistrial.targetLoc}, stimulus.stair{task.thistrial.targetContrast}.threshold));
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen)

% compute the grating
grating = mglMakeGrating(stimulus.width, stimulus.width, stimulus.sf, stimulus.orientation,0);
gaussian = mglMakeGaussian(stimulus.width, stimulus.width, stimulus.width/8, stimulus.width/8); %sd?????

gaborLow = 255*(stimulus.contrastLow*grating.*gaussian+1)/2;
gaborHigh = 255*(stimulus.contrastHigh*grating.*gaussian+1)/2;

%create texture
stimulus.tex(1) = mglCreateTexture(gaborLow);
stimulus.tex(2) = mglCreateTexture(gaborHigh);

%stim loc
% angles  (1->4 clockwise)
% loc1: pi/3 (60)  loc2: pi/6 (30)  loc3: pi*4/3 (240) loc4: pi*7/6 (210)
%cue lines / stimulus centers
for i = 1:length(stimulus.angles)
      % starting point of lines
      [stimulus.cueLines.x0(i), stimulus.cueLines.y0(i)] = pol2cart(stimulus.angles(i), 0.05);
      % ending point
      [stimulus.cueLines.x1(i), stimulus.cueLines.y1(i)] = pol2cart(stimulus.angles(i), 0.6);
      
      %stim centers
      [stimulus.x(i), stimulus.y(i)] = pol2cart(stimulus.angles(i),stimulus.radius);
end

% now create staircases for each orientation
for i = 1:2
  % init a 2 down 1 up staircase
  stimulus.stair{i} = upDownStaircase(1,2,stimulus.initialThreshold,stimulus.stepsize);
  stimulus.stair{i}.minThreshold = stimulus.minThreshold;
  stimulus.stair{i}.maxThreshold = stimulus.maxThreshold;
  
end

