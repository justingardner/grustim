function myscreen = testLoc(varargin)

clear global stimulus
global stimulus

scan = 0;
getArgs(varargin,{'scan=0'});

stimulus.scan = scan;
stimulus.width = 12;
stimulus.sf = 1.8;

% initialize the screen
myscreen.background = 0.5;
myscreen = initScreen(myscreen);

% set the first task to be the fixation staircase task
[task{1} myscreen] = fixStairInitTask(myscreen);

%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%%
task{2}{1}.waitForBacktick = 1;
task{2}{1}.segmin = [1 2 1];
task{2}{1}.segmax = [1 2 3];

if stimulus.scan
    task{2}{1}.synchToVol = [0 0 1];
end

% parameters & randomization
task{2}{1}.parameter.contrast = [.2 .8];
task{2}{1}.parameter.location = [-10 -8 8 10];
task{2}{1}.randVars.block.orientation = [0 45 90 135];
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

    stimulus.contrast = task.thistrial.contrast;
    stimulus.location = task.thistrial.location;
    stimulus.orientation = task.thistrial.orientation;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task, myscreen)

global stimulus
mglClearScreen;
if task.thistrial.thisseg == 2
stimulus = updateGabor(stimulus,myscreen);

% mglBltTexture(stimulus.tex, [stimulus.location, 0], 0,0, stimulus.orientation);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%`%%%%%%%%%%%%
function stimulus = initGabor(stimulus,myscreen)
% compute the grating
stimulus.grating = mglMakeGrating(stimulus.width,stimulus.width,stimulus.sf, 0,0);
stimulus.gaussian = mglMakeGaussian(stimulus.width,stimulus.width,stimulus.width/10, stimulus.width/10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to update the stimulus and dra`w it to screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = updateGabor(stimulus,myscreen)
% create texture
stimulus.gabor = (255*(stimulus.contrast*stimulus.grating.*stimulus.gaussian+1)/2);
stimulus.tex = mglCreateTexture(stimulus.gabor);
% 
mglBltTexture(stimulus.tex, [stimulus.location, 0], 0,0, stimulus.orientation);

