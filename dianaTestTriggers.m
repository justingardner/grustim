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
myscreen = initScreen(); mglVisualAngleCoordinates(57,[16 12]); mglClearScreen; task{1}.waitForBacktick = 1;

% init the stimulus
global stimulus;

myscreen = initStimulus('stimulus',myscreen)
    
    % Gabor stimulus %
    gaussian = mglMakeGaussian(3,3,5,5); 
    stimulus = mglCreateTexture(gaussian);

    % checkerboard stimulus %
    %grating1 = mglMakeGrating(10,10,1.5,45,0);
    %grating2 = mglMakeGrating(10,10,1.5,135,0);
    %checkerboard = 255*(sign(grating1/2+grating2/2))/2;
    %stimulus = mglCreateTexture(checkerboard);

% set task parameters
task{1}.segmin = [1 1 1]; task{1}.segmax = [1 1 1]; task{1}.getResponse = [0 0 0]; %length of segments and response segment

task{1}.numTrials = 10;

    % synch to vol
    task{1}.synchToVol = [1 0 0];
    task{1}.segmin = task{1}.segmin - [0 0 0];
    task{1}.segmax = task{1}.segmax - [0 0 0];


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

mglClearScreen(0);
mglFlush;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%resp%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen, checkboardBltTex)

if task.thistrial.thisseg == 2

    global stimulus
    
    mglClearScreen(0);
    checkerboardBltTex = mglBltTexture(stimulus,[0 0]);

end

mglFlush;



%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = getResponseCallback(task,myscreen)

