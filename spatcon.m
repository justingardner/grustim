% od
%
%      usage: myscreen=od(stimulus)
%         by: justin gardner
%       date: 04/15/06
%    purpose: orientation discrimination task
%
%
%
function myscreen = spatcon

% check arguments
if ~any(nargin == [0])
  help orientationDiscrimination
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other screen parameters
myscreen.autoCloseScreen = 1;
myscreen.allowpause = 0;
myscreen.eatkeys = 0;
myscreen.displayname = 'projector';
myscreen.background = 'gray';

myscreen = initScreen(myscreen);

global MGL;
clear global stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);

% parameters
stimulus.grating.orientations = [90];
stimulus.grating.contrasts = [0:0.1:1];
stimulus.grating.sf = 2;
stimulus.grating.tf = 1;
stimulus.grating.width = 7;
stimulus.grating.height = 7;
stimulus.grating.phase = 0;
stimulus.grating.windowType = 'gabor'; % should be gabor or thresh
stimulus.grating.sdx = stimulus.grating.width/6;
stimulus.grating.sdy = stimulus.grating.width/6;

stimulus.fix.width = 0.5;
stimulus.fix.linewidth = 1;
stimulus.fix.disksize = [2 2];
stimulus.fix.startColor = 1;
stimulus.fix.correctColor = [0 1 0];
stimulus.fix.incorrectColor = [1 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up segments of trials
task{1}.parameter.type = [0 1];
task{1}.random = 1;
task{1}.seglen = [3];
task{1}.synchToVol = [0 0 0];
task{1}.getResponse = [0 0 1];
task{1}.waitForBacktick = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialze tasks and stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stimulus = initGratings(stimulus,task,myscreen);

% initialze tasks
[task{1} myscreen] = initTask(task{1},myscreen,@startSegmentCallback,@trialStimulusCallback,@trialResponseCallback);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the tasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set which phase is active
tnum = 1;

while (tnum <= length(task)) && ~myscreen.userHitEsc
  % updatethe task
  [task myscreen tnum] = updateTask(task,myscreen,tnum);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     STIMULUS Gratings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initGratings(stimulus,myscreen,task)

disppercent(-inf,'Creating grating textures');
% make each one of he called for gratings
for i = 1:length(stimulus.grating.contrasts)
  disppercent((i-1)/length(stimulus.grating.contrasts));
  % get the orientation we want to create
  thisOrientation = stimulus.grating.orientations(1);
  % make the grating
  thisGrating = 255*(mglMakeGrating(stimulus.grating.width,stimulus.grating.height,...
				    stimulus.grating.sf,thisOrientation,...
				    stimulus.grating.phase)+1)/2;
  
  thisGaussian = mglMakeGaussian(stimulus.grating.width,stimulus.grating.height,...
				 stimulus.grating.sdx,stimulus.grating.sdy);

  % create an rgb/a matrix
  thisStimulus(:,:,1) = thisGrating;
  thisStimulus(:,:,2) = thisGrating;
  thisStimulus(:,:,3) = thisGrating;

  % make the gaussian window
  if strcmp(stimulus.grating.windowType,'gabor')
    % create an rgb/a matrix
    mask = stimulus.grating.contrasts(i)*255*thisGaussian;
  else
    mask = stimulus.grating.contrasts(i)*255*(thisGaussian>exp(-1/2));
  end

  thisStimulus(:,:,4) = mask;
  thatStimulus(:,:,4) = mask;

  % create the texture
  stimulus.tex(i) = mglCreateTexture(thisStimulus);
end
disppercent(inf);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to set stimulus parameters at
% the beginning of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task,myscreen)

global stimulus;
stimulus.fix.color = stimulus.fix.startColor;
task.thistrial.r = ceil(length(stimulus.tex)*rand(8));
task.thistrial.thisloc = ceil(8*rand(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to display stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = trialStimulusCallback(task,myscreen)

global stimulus;
mglClearScreen;

angles = (0:45:359)-45;
radius = 10;

for i = 1:length(angles)
  x = radius*cos(pi*angles(i)/180);
  y = radius*sin(pi*angles(i)/180);

  if 1
    if any(i == 2)
      contrastNum = 2;
    else
      contrastNum = 10;
    end    
  end

  if 0  
    if i ~= task.thistrial.thisloc
      if task.thistrial.type == 1
	contrastNum = task.thistrial.r(i);
      else
	contrastNum = task.thistrial.r(1);
      end
    else
      contrastNum = 4;
    end    
  end

  mglBltTexture(stimulus.tex(contrastNum),[x y]);
end

% draw fixation
mglFillOval(0,0,stimulus.fix.disksize,myscreen.background);
mglFixationCross(stimulus.fix.width,stimulus.fix.linewidth,stimulus.fix.color);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to set stimulus parameters at
% the beginning of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = trialResponseCallback(task,myscreen)

global stimulus;

