% attentionCompLoc.m
%
%      usage: myscreen = attentionCompLoc()
%         by: minyoung lee
%       date: 07/07/2016
%    purpose: localizer for biased competition experiment

% .5 TR
% 10 cycle * 48 vol/cycle + 8(mux)*2 = 496 vols

function myscreen = attentionCompLoc()

% initalize the screen
myscreen.background = 0.5;
myscreen = initScreen(myscreen);
% myscreen = initScreen('fMRIprojFlex');
% myscreen.background = 0.5;

% set the first task to be the fixation staircase task
clear global fixStimulus;
[task{1} myscreen] = fixStairInitTask(myscreen);


TR = .5;
task{2}{1}.seglen = [6 6 6 6-TR/2];
task{2}{1}.synchToVol = [0 0 0 1];
task{2}{1}.numTrials = 10;
task{2}{1}.waitForBacktick = 1;

% initialize our task
for phaseNum = 1:length(task{2})
  [task{2}{phaseNum} myscreen] = initTask(task{2}{phaseNum},myscreen,@startSegmentCallback,@updateScreenCallback);
end

% init the stimulus
clear global stimulus;
global stimulus;

myscreen = initStimulus('stimulus',myscreen);
stimulus.width = 6; % gabor
stimulus.sf = 1.8;
stimulus.orientation = 90;
stimulus.contrast = 1;
% angles  (1->4 clockwise from top)
% loc1: 60  loc2: 30  loc3: 240 loc4: 210
stimulus.angles = [65 25 245 205] * pi /180;
% stimulus.angles = [67.5 22.5 247.5 202.5] * pi /180;
stimulus.locations = 1:length(stimulus.angles);
stimulus.radius = 6;

stimulus.flickerRate = 2; % 2 Hz
stimulus.flickerDur = 1/stimulus.flickerRate;
stimulus.flickerNFrame = round(myscreen.framesPerSecond/stimulus.flickerRate);

% to initialize the stimulus for your experiment.
stimulus = initGabor(stimulus,myscreen);
stimulus.cycleTime = mglGetSecs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task{2})) && ~myscreen.userHitEsc
  % update the gabors
  [task{2} myscreen phaseNum] = updateTask(task{2},myscreen,phaseNum);
  % update the fixation task
  [task{1} myscreen] = updateTask(task{1},myscreen,1);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;

if task.thistrial.thisseg == 1
  
  % show the amount of time the last cycle took
  if task.trialnum > 1
    disp(sprintf('%i: %f last cycle length',task.trialnum,mglGetSecs(stimulus.cycleTime)));
    stimulus.cycleTime = mglGetSecs;
  end
  
end
stimulus.flickerIndex = 0;
stimulus.gaborIndex = 0;
stimulus.gaborNum = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task, myscreen)
mglClearScreen
global stimulus

% if (mglGetSecs-task{tnum}.thistrial.segstart)
if mod(stimulus.flickerIndex, stimulus.flickerNFrame) == 0 
    stimulus.gaborIndex = stimulus.gaborIndex + 1;
    if mod(stimulus.gaborIndex,2)==1
        stimulus.gaborNum = 1;
    else
        stimulus.gaborNum = 2;
    end
end
stimulus.flickerIndex = stimulus.flickerIndex+1;
    
for i = 1:2
    mglBltTexture(stimulus.tex(stimulus.gaborNum), [stimulus.x(stimulus.seq(i,task.thistrial.thisseg)) stimulus.y(stimulus.seq(i,task.thistrial.thisseg))], 0,0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initGabor(stimulus,myscreen)

% compute the grating
grating{1} = mglMakeGrating(stimulus.width, stimulus.width, stimulus.sf, stimulus.orientation,0);
grating{2} = mglMakeGrating(stimulus.width, stimulus.width, stimulus.sf, stimulus.orientation,180);
gaussian = mglMakeGaussian(stimulus.width, stimulus.width,stimulus.width/10, stimulus.width/10); %sd?????

for n = 1:2
gabor{n}(:,:,1) = 255*(grating{n}+1)/2;
gabor{n}(:,:,2) = 255*(grating{n}+1)/2;
gabor{n}(:,:,3) = 255*(grating{n}+1)/2;
gabor{n}(:,:,4) =255*(stimulus.contrast*gaussian);

%create texture
stimulus.tex(n) = mglCreateTexture(gabor{n});
end

% stimulus centers
for i = 1:length(stimulus.angles)
      %stim centers
      [stimulus.x(i), stimulus.y(i)] = pol2cart(stimulus.angles(i),stimulus.radius);      
end

stimulus.seq = [1 2 3 4; 2 3 4 1];
