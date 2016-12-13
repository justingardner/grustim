% testLoc.m
%
%      usage: myscreen=testLoc()
%         by: minyoung lee
%       date: 11/02/16
%    purpose: generates stimuli for sensory uncertainty (location/contrast) experiment

function myscreen = testLoc(varargin)

clear global stimulus
global stimulus

scan = 0;
getArgs(varargin,{'scan=0'});
stimulus.scan = scan;

% stimulus parameters
stimulus.width = 12;
stimulus.sf = 1.8;

% initialize the screen
myscreen.background = 0.5;
myscreen = initScreen(myscreen);

% set the first task to be the fixation staircase task
[task{1} myscreen] = fixStairInitTask(myscreen);

%set params for flicker
stimulus.flickerRate = 2; % 2 Hz
stimulus.flickerDur = 1/stimulus.flickerRate;
stimulus.flickerNFrame = round(myscreen.framesPerSecond/stimulus.flickerRate);

%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%%
task{2}{1}.waitForBacktick = 1;
task{2}{1}.segmin = [2 1];
task{2}{1}.segmax = [2 5];

if stimulus.scan
    task{2}{1}.synchToVol = [0 1];
end

% parameters & randomization
task{2}{1}.parameter.contrast = [.25 .5 1];
task{2}{1}.parameter.location = [-14 -8 8 14];
task{2}{1}.randVars.block.orientation = 0:45:135;
task{2}{1}.random = 1;

% initialize the task
for phaseNum = 1:length(task{2})
  [task{2}{phaseNum} myscreen] = initTask(task{2}{phaseNum},myscreen,@startSegmentCallback,@updateScreenCallback);
end

% init the stimulus
myscreen = initStimulus('stimulus',myscreen);

% to initialize the stimulus for your experiment.
stimulus = initGabor(stimulus,myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task{2})) && ~myscreen.userHitEsc
  % update the gabor
  [task{2} myscreen phaseNum] = updateTask(task{2},myscreen,phaseNum);
  % update the fixation task
  [task{1} myscreen] = updateTask(task{1},myscreen,1);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);
clear global stimulus;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)
global stimulus;
if task.thistrial.thisseg == 1
    stimulus.contrast = task.thistrial.contrast;
    stimulus.location = task.thistrial.location;
    stimulus.orientation = task.thistrial.orientation;
    stimulus.flickerIndex = 0;
    stimulus.gaborIndex = 0;
    stimulus.gaborNum = 1;
    
    stimulus = updateGabor(stimulus,myscreen);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task, myscreen)
global stimulus
mglClearScreen;
if task.thistrial.thisseg == 2
    
    
% show a flickering gabor patch
if mod(stimulus.flickerIndex, stimulus.flickerNFrame) == 0 
    stimulus.gaborIndex = stimulus.gaborIndex + 1;
    if mod(stimulus.gaborIndex,2)==1
        stimulus.gaborNum = 1;
    else
        stimulus.gaborNum = 2;
    end
end
stimulus.flickerIndex = stimulus.flickerIndex+1;
    
mglBltTexture(stimulus.tex(stimulus.gaborNum), [stimulus.location, 0], 0,0, stimulus.orientation);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%`%%%%%%%%%%%%
function stimulus = initGabor(stimulus,myscreen)
% compute the grating
stimulus.grating1 = mglMakeGrating(stimulus.width,stimulus.width,stimulus.sf, 0,0);
stimulus.grating2 = mglMakeGrating(stimulus.width,stimulus.width,stimulus.sf, 0,180);
stimulus.gaussian = mglMakeGaussian(stimulus.width,stimulus.width,stimulus.width/8, stimulus.width/8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to update the stimulus and draw it to screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = updateGabor(stimulus,myscreen)
% set contrast & create texture
stimulus.gabor1 = (255*(stimulus.contrast*stimulus.grating1.*stimulus.gaussian+1)/2);
stimulus.gabor2 = (255*(stimulus.contrast*stimulus.grating2.*stimulus.gaussian+1)/2);
stimulus.tex(1) = mglCreateTexture(stimulus.gabor1);
stimulus.tex(2) = mglCreateTexture(stimulus.gabor2);
