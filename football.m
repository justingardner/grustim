% football.m
%
%        $Id: taskTemplate.m 833 2010-06-28 07:41:54Z justin $
%      usage: taskTemplate
%         by: justin gardner
%       date: 09/07/06
%  copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%    purpose: 
%
function football

% check arguments
if ~any(nargin == [0])
  help taskTemplate
  return
end

% initalize the screen
myscreen = initScreen(myscreen);

task{1}.waitForBacktick = 0;
% Seg 1 is first fixation. Seg 3 the target changes
% Seg 4 the target moves. seg 5 is response
task{1}.segmin = [1 1 0 1 1];
task{1}.segmax = [1 1 1 1 1];
task{1}.getResponse = [0 0 1 1 1];
% fix: enter the parameter of your choice
task{1}.randVars.uniform.saccadeContingent = [0 1];
task{1}.randVars.uniform.presentAbsent = [0 1];
task{1}.random = 1;
task{1}.numTrials = 100;

% initialize the task
for phaseNum = 1:length(task)
  [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@drawStimulusCallback,@responseCallback);
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

if task{1}.trialnum > 10
  e = getTaskParameters(myscreen,task);
  saccadeContingent = e.randVars.saccadeContingent;
  presentAbsent = e.randVars.presentAbsent;
  response = (e.response == 1);
  % compute saccadeContingent hits and falseAlarams
  hits = sum(presentAbsent(find(saccadeContingent)) & response(find(saccadeContingent)));
  hitsn = sum(presentAbsent(find(saccadeContingent)));
  falseAlarms = sum(~presentAbsent(find(saccadeContingent)) & response(find(saccadeContingent)));
  falseAlarmsN = sum(~presentAbsent(find(saccadeContingent)));
  disp(sprintf('(football) SaccadeContingent hits: %0.2f%% (%i/%i) false alarms: %0.2f%% (%i/%i)',100*hits/hitsn,hits,hitsn,100*falseAlarms/falseAlarmsN,falseAlarms,falseAlarmsN));

  hits = sum(presentAbsent(find(~saccadeContingent)) & response(find(~saccadeContingent)));
  hitsn = sum(presentAbsent(find(~saccadeContingent)));
  falseAlarms = sum(~presentAbsent(find(~saccadeContingent)) & response(find(~saccadeContingent)));
  falseAlarmsN = sum(~presentAbsent(find(~saccadeContingent)));
  disp(sprintf('(football) Not saccadeContingent hits: %0.2f%% (%i/%i) false alarms: %0.2f%% (%i/%i)',100*hits/hitsn,hits,hitsn,100*falseAlarms/falseAlarmsN,falseAlarms,falseAlarmsN));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;

if task.thistrial.thisseg == 1
  if task.trialnum > 5
    stimulus.displayEyeContingent = 0;
  end
  % start trial without the extra target
  stimulus.addTarget = 0;

  % delete old text texture
  if ~isempty(stimulus.textTex)
    mglDeleteTexture(stimulus.textTex);
    stimulus.textTex = [];
  end

  % clear screen
  mglClearScreen(0.5);
  
  % make stimulus locations
  for i = 1:stimulus.nStimuli
    s = rand*4+2;
    stimulus.loc(i,1:4) = [rand*stimulus.deviceWidth-stimulus.deviceWidth/2 rand*stimulus.deviceHeight-stimulus.deviceHeight/2 s s];
  end

  % get target locations
  stimulus.targetLoc1 = stimulus.loc(stimulus.nStimuli,:);
  stimulus.targetLoc2 = stimulus.loc(stimulus.nStimuli-1,:);
  stimulus.targetLoc = stimulus.targetLoc1;
  
  % remember some targetLocs
  stimulus.locs{task.trialnum}.targetLoc1 = stimulus.targetLoc1; 
  stimulus.locs{task.trialnum}.targetLoc2 = stimulus.targetLoc2; 
  stimulus.locs{task.trialnum}.newTarget = stimulus.loc(1,:); 
  stimulus.locs{task.trialnum}.locs = stimulus.loc;
elseif task.thistrial.thisseg == 3
  stimulus.targetLoc = stimulus.targetLoc2;
end  


%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)

global stimulus

if task.thistrial.gotResponse < 1
  if task.thistrial.presentAbsent && (task.thistrial.whichButton == 1)
    stimulus.type(task.trialnum) = 'H';
    stimulus.correct(task.trialnum) = 1;
  elseif task.thistrial.presentAbsent && (task.thistrial.whichButton == 2)
    stimulus.type(task.trialnum) = 'M';
    stimulus.correct(task.trialnum) = 0;
  elseif ~task.thistrial.presentAbsent && (task.thistrial.whichButton == 2)
    stimulus.type(task.trialnum) = 'C';
    stimulus.correct(task.trialnum) = 1;
  elseif ~task.thistrial.presentAbsent && (task.thistrial.whichButton == 1)
    stimulus.type(task.trialnum) = 'F';
    stimulus.correct(task.trialnum) = 0;
  end
end

dispText = sprintf('%s: saccadeContingent: %i',stimulus.type(task.trialnum),task.thistrial.saccadeContingent);
disp(sprintf('(footbal) %s',dispText));

stimulus.textTex = mglText(dispText);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = drawStimulusCallback(task, myscreen)

global stimulus


fixAccuracy = 5;
mglClearScreen;

if task.thistrial.thisseg == 5
  if ~isempty(stimulus.textTex)
    if stimulus.correct(task.trialnum)
      mglBltTexture(stimulus.greentex,[0 0]);
    else
      mglBltTexture(stimulus.redtex,[0 0]);
    end
%    mglBltTexture(stimulus.textTex,[0 0]);
  end
  return
end

for i = 2:stimulus.nStimuli
  mglBltTexture(stimulus.tex,stimulus.loc(i,:));
end

% check if we have fixed target 1
if ((myscreen.eyetracker.eyepos(1) > (stimulus.targetLoc1(1)-stimulus.targetLoc1(3)/2-fixAccuracy/2)) &&...
    (myscreen.eyetracker.eyepos(1) < (stimulus.targetLoc1(1)+stimulus.targetLoc1(3)/2+fixAccuracy/2)) &&...
    (myscreen.eyetracker.eyepos(2) > (stimulus.targetLoc1(2)-stimulus.targetLoc1(4)/2-fixAccuracy/2)) &&...
    (myscreen.eyetracker.eyepos(2) < (stimulus.targetLoc1(2)+stimulus.targetLoc1(4)/2+fixAccuracy/2)))
  stimulus.targetFixed1 = 1;
else
  stimulus.targetFixed1 = 0;
end

% check if we have fixed target 2
if ((myscreen.eyetracker.eyepos(1) > (stimulus.targetLoc2(1)-stimulus.targetLoc2(3)/2-fixAccuracy/2)) &&...
    (myscreen.eyetracker.eyepos(1) < (stimulus.targetLoc2(1)+stimulus.targetLoc2(3)/2+fixAccuracy/2)) &&...
    (myscreen.eyetracker.eyepos(2) > (stimulus.targetLoc2(2)-stimulus.targetLoc2(4)/2-fixAccuracy/2)) &&...
    (myscreen.eyetracker.eyepos(2) < (stimulus.targetLoc2(2)+stimulus.targetLoc2(4)/2+fixAccuracy/2)))
  stimulus.targetFixed2 = 1;
else
  stimulus.targetFixed2 = 0;
end

% display the target
mglBltTexture(stimulus.redtex,stimulus.targetLoc);

if stimulus.displayEyeContingent
  % display in green if we have saccaded to target
  if stimulus.targetFixed1
    mglBltTexture(stimulus.greentex,stimulus.targetLoc1);
  end

  % display in green if we are looking at second target
  if (task.thistrial.thisseg >= 2) && stimulus.targetFixed2
    mglBltTexture(stimulus.greentex,stimulus.targetLoc2);
  end

  % display eye positiona
  mglGluDisk(myscreen.eyetracker.eyepos(1),myscreen.eyetracker.eyepos(2), 0.1,[0.3 05 0.8]);
end

% if we are supposed to add a target
if task.thistrial.presentAbsent
  if task.thistrial.saccadeContingent
    if (task.thistrial.thisseg == 3) && (stimulus.targetFixed1 == 0)
      stimulus.addTarget = 1;
    end
  elseif task.thistrial.thisseg == 4
    stimulus.addTarget = 1;
  end
end
% add the extra target if necessary
if stimulus.addTarget
  mglBltTexture(stimulus.tex,stimulus.loc(1,:));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen)

% load the stimulus
im1 = imread('~/proj/grustim/images/changeBlindness/ball.tif');
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
stimulus.deviceWidth = mglGetParam('deviceWidth')*3/4;
stimulus.deviceHeight = mglGetParam('deviceHeight')*3/4;

% number of stimuli
stimulus.nStimuli = 10;

stimulus.textTex = [];

stimulus.displayEyeContingent = 1;