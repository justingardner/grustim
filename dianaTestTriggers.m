%    Diana.m
%
%       
%       usage: Diana
%       by: josh wilson
%       date: October 2022
%       
%       Test stimulus for DIANA imaging.
%      




function myscreen = Diana(varargin)

% check arguments
getArgs(varargin);

% initilaize the screen
myscreen = initScreen();% mglVisualAngleCoordinates(57,[16 12]); mglClearScreen;
task{1}.waitForBacktick = 1;

% init the stimulus
global stimulus;

myscreen = initStimulus('stimulus',myscreen)
    
% set where the quad will be
stimulus.xPosMin = -5;
stimulus.xPosMax = 5;
stimulus.yPosMin = -5;
stimulus.yPosMax = 5;
% compute location of quad
stimulus.quadX = [stimulus.xPosMin stimulus.xPosMin stimulus.xPosMax stimulus.xPosMax]';
stimulus.quadY = [stimulus.yPosMin stimulus.yPosMax stimulus.yPosMax stimulus.yPosMin]';

% set task parameters
task{1}.segmin = [0.1 0.5];
task{1}.segmax = [0.1 0.5];
task{1}.getResponse = [0 0];

task{1}.numTrials = inf;

% synch to vol
task{1}.synchToVol = [1 0];


% initialize the task
for phaseNum = 1:length(task)
    [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback);
end


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

if task.thistrial.thisseg == 1
    disp(sprintf('(dianaTestTriggers) Starting trial %i',task.trialnum));
    mglClearScreen(0);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%resp%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen, checkboardBltTex)

global stimulus
    

if task.thistrial.thisseg == 2

    if isodd(myscreen.tick)
      mglQuad(stimulus.quadX,stimulus.quadY,[1;1;1]);
    else
      mglQuad(stimulus.quadX,stimulus.quadY,[0;0;0]);
    end
        
else
    mglClearScreen(0);
end





%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = getResponseCallback(task,myscreen)

