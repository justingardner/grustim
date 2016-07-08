% attentionCompLoc.m
%
%      usage: myscreen = attentionCompLoc()
%         by: minyoung lee
%       date: 07/07/2016
%    purpose: localizer for biased competition experiment
%
function myscreen = attentionCompLoc()

% initalize the screen
myscreen.background = 0.5;
myscreen = initScreen(myscreen);
% myscreen = initScreen('fMRIprojFlex');
% myscreen.background = 0.5;

% set the first task to be the fixation staircase task
clear global fixStimulus;
[task{1} myscreen] = fixStairInitTask(myscreen);

task{2}{1}.seglen = [6 6 6 6];
% task{2}{1}.synchToVol = [0 0 0 0];
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task, myscreen)
mglClearScreen
global stimulus
for i = 1:2
    mglBltTexture(stimulus.tex, [stimulus.x(stimulus.seq(i,task.thistrial.thisseg)) stimulus.y(stimulus.seq(i,task.thistrial.thisseg))], 0,0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initGabor(stimulus,myscreen)

% compute the grating
grating = mglMakeGrating(stimulus.width, stimulus.width, stimulus.sf, stimulus.orientation,0);
gaussian = mglMakeGaussian(stimulus.width, stimulus.width,stimulus.width/10, stimulus.width/10); %sd?????

gabor(:,:,1) = 255*(grating+1)/2;
gabor(:,:,2) = 255*(grating+1)/2;
gabor(:,:,3) = 255*(grating+1)/2;
gabor(:,:,4) =255*(stimulus.contrast*gaussian);

%create texture
stimulus.tex = mglCreateTexture(gabor);

% stimulus centers
for i = 1:length(stimulus.angles)
      %stim centers
      [stimulus.x(i), stimulus.y(i)] = pol2cart(stimulus.angles(i),stimulus.radius);      
end

stimulus.seq = [1 2 3 4; 2 3 4 1];
