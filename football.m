% football.m
%
%        $Id: taskTemplate.m 833 2010-06-28 07:41:54Z justin $
%      usage: taskTemplate
%         by: justin gardner
%       date: 09/07/06
%  copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%    purpose: 
%
function myscreen = football

% check arguments
if ~any(nargin == [0])
  help taskTemplate
  return
end

% initalize the screen
myscreen = initScreen;

task{1}.waitForBacktick = 1;
% fix: the task defined here has two segments, one that
% is 3 seconds long followed by another that is 
% 6-9 seconds (randomized in steps of 1.5 seconds)
% change this to what you want for your trial
task{1}.segmin = [1 1 1];
task{1}.segmax = [1 1 1];
% fix: enter the parameter of your choice
task{1}.parameter.myParameter = [0 30 90];
task{1}.random = 1;

% initialize the task
for phaseNum = 1:length(task)
  [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@stimStartSegmentCallback,@stimDrawStimulusCallback);
end

% init the stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);

% fix: you will change the funciton myInitStimulus
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = stimStartSegmentCallback(task, myscreen)

global stimulus;

if task.thistrial.thisseg == 1
  mglClearScreen(0.5);
  % make stimulus locations
  for i = 1:stimulus.nStimuli
    s = rand*4+2;
    stimulus.loc(i,1:4) = [rand*stimulus.deviceWidth-stimulus.deviceWidth/2 rand*stimulus.deviceHeight-stimulus.deviceHeight/2 s s];
    
  end

  stimulus.targetLoc = stimulus.loc(stimulus.nStimuli,:);
elseif task.thistrial.thisseg == 3
  stimulus.targetLoc = stimulus.loc(stimulus.nStimuli-1,:);
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = stimDrawStimulusCallback(task, myscreen)

global stimulus


fixAccuracy = 1;
mglClearScreen;
for i = 1:stimulus.nStimuli
  mglBltTexture(stimulus.tex,stimulus.loc(i,:));
end

if ((myscreen.eyetracker.eyepos(1) > (stimulus.targetLoc(1)-stimulus.targetLoc(3)/2-fixAccuracy/2)) &&...
    (myscreen.eyetracker.eyepos(1) < (stimulus.targetLoc(1)+stimulus.targetLoc(3)/2+fixAccuracy/2)) &&...
    (myscreen.eyetracker.eyepos(2) > (stimulus.targetLoc(2)-stimulus.targetLoc(4)/2-fixAccuracy/2)) &&...
    (myscreen.eyetracker.eyepos(2) < (stimulus.targetLoc(2)+stimulus.targetLoc(4)/2+fixAccuracy/2)))
  % display stimuli
  mglBltTexture(stimulus.greentex,stimulus.targetLoc);
else
  mglBltTexture(stimulus.redtex,stimulus.targetLoc);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen)

% load the stimulus
im1 = imread('ball.tif');
im2(:,:,1) = double(flipud(im1(:,:,1)));
im2(:,:,2) = im2(:,:,1);
im2(:,:,3) = im2(:,:,1);
im2(:,:,4) = im1(:,:,2);
stimulus.tex = mglCreateTexture(im2);

% make red one
im2(:,:,1) = double(flipud(im1(:,:,1)));
im2(:,:,2) = round(0.2*im2(:,:,1));
im2(:,:,3) = round(0.2*im2(:,:,1));
im2(:,:,4) = im1(:,:,2);
stimulus.redtex = mglCreateTexture(im2);

% make green one
im2(:,:,2) = double(flipud(im1(:,:,1)));
im2(:,:,1) = round(0.2*im2(:,:,2));
im2(:,:,3) = round(0.2*im2(:,:,2));
im2(:,:,4) = im1(:,:,2);
stimulus.greentex = mglCreateTexture(im2);

% make a texture
stimulus.deviceWidth = mglGetParam('deviceWidth');
stimulus.deviceHeight = mglGetParam('deviceHeight');

% number of stimuli
stimulus.nStimuli = 20;