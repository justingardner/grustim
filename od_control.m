% od
%
%      usage: myscreen=od(stimulus)
%         by: justin gardner
%       date: 04/15/06
%    purpose: orientation discrimination task
%
%
%
function myscreen = od_control1

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
stimulus.grating.orientations = [90-22 90 90+22];
stimulus.grating.orientationDiffs = [-3 -1.5 -0.75 0 0.75 1.5 3];
%stimulus.grating.orientationDiffs = [0.0];
%stimulus.grating.orientationDiffs = [-15 0 15];
stimulus.grating.sf = 2;
stimulus.grating.tf = 1;
stimulus.grating.width = 20;
stimulus.grating.height = 20;
%stimulus.grating.width = 7;
%stimulus.grating.height = 7;
stimulus.grating.phase = 0;
stimulus.grating.windowType = 'thresh'; % should be gabor or thresh
stimulus.grating.sdx = stimulus.grating.width/2;
stimulus.grating.sdy = stimulus.grating.width/2;
stimulus.grating.contrast = 50;

stimulus.cue.width = stimulus.grating.width/2+1.5;
stimulus.cue.linewidth = 10;
stimulus.cue.startColor = 1;
stimulus.cue.cueColor = 0;
stimulus.cue.disksize = [stimulus.grating.width+.5 stimulus.grating.width+.5];
stimulus.cue.origin = [0 0];
stimulus.cue.tilt = stimulus.grating.orientations;

stimulus.fix.width = 0.5;
stimulus.fix.linewidth = 1;
stimulus.fix.disksize = [2 2];
stimulus.fix.startColor = 1;
stimulus.fix.correctColor = [0 1 0];
stimulus.fix.incorrectColor = [1 0 0];
stimulus.timer = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up segments of trials
task{1}.parameter.orientNum = 1:length(stimulus.grating.orientations);
task{1}.parameter.orientDiffNum = 1:length(stimulus.grating.orientationDiffs);
task{1}.random = 1;
task{1}.segmin = [0.8*3 .2 0.8*3];
task{1}.segmax = [0.8*9 .2 0.8*3];
task{1}.segquant = [0.8 0 0];
task{1}.synchToVol = [0 0 1];
%task{1}.segmin = [0.26 1 1.26*4];
%task{1}.segmax = [0.26 1 1.26*7];
%task{1}.synchToVol = [0 0 0];
task{1}.getResponse = [0 0 1];
task{1}.waitForBacktick = 1;
%task{1}.numBlocks = 25;
task{1}.writeTrace{2}.tracenum = [1 2];
task{1}.writeTrace{2}.tracevar{1} = 'orientNum';
task{1}.writeTrace{2}.tracevar{2} = 'orientDiffNum';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialze tasks and stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stimulus = initGratings(stimulus,task,myscreen);

% initialze tasks
task{1} = initTask(task{1},myscreen,@startSegmentCallback,@trialStimulusCallback,@trialResponseCallback);

task{1}.responses = [];
task{1}.orientationDiff = [];
task{1}.orientation = [];
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

%task{1}.pright = task{1}.right./(task{1}.right+task{1}.left)

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
function stimulus = initGratings(stimulus,task,myscreen)

disppercent(-inf,'Creating grating textures');
% make each one of he called for gratings
for i = 1:length(stimulus.grating.orientations)
  for j = 1:length(stimulus.grating.orientationDiffs)
    disppercent((i-1)/length(stimulus.grating.orientations)+j/length(stimulus.grating.orientationDiffs));
    % get the orientation we want to create
    thisOrientation = stimulus.grating.orientations(i)+stimulus.grating.orientationDiffs(j);
    % make the grating
    thisGrating = 255*(makeGrating(stimulus.grating.width,stimulus.grating.height,...
				   stimulus.grating.sf,thisOrientation,...
				   stimulus.grating.phase)+1)/2;

    thatGrating = 255*(makeGrating(stimulus.grating.width,stimulus.grating.height,...
				   stimulus.grating.sf,thisOrientation,...
				   stimulus.grating.phase+180)+1)/2;

    thisGaussian = makeGaussian(stimulus.grating.width,stimulus.grating.height,...
				stimulus.grating.sdx,stimulus.grating.sdy);

    % create an rgb/a matrix
    thisStimulus(:,:,1) = thisGrating;
    thisStimulus(:,:,2) = thisGrating;
    thisStimulus(:,:,3) = thisGrating;

    thatStimulus(:,:,1) = thatGrating;
    thatStimulus(:,:,2) = thatGrating;
    thatStimulus(:,:,3) = thatGrating;

    % make the gaussian window
    if strcmp(stimulus.grating.windowType,'gabor')
      % create an rgb/a matrix
      mask = stimulus.grating.contrast*255*thisGaussian/100;
    else
      mask = stimulus.grating.contrast*255*(thisGaussian>exp(-1/2))/100;
    end

      thisStimulus(:,:,4) = mask;
      thatStimulus(:,:,4) = mask;

    % create the texture
    stimulus.thistex(i,j) = mglCreateTexture(thisStimulus);
    stimulus.thattex(i,j) = mglCreateTexture(thatStimulus);
  end

  stimulus.cue.x(i) = cos(pi*(stimulus.cue.tilt(i))/180)*stimulus.cue.width;
  stimulus.cue.y(i) = sin(pi*(stimulus.cue.tilt(i))/180)*stimulus.cue.width;
  
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
if (task.thistrial.thisseg == 2)
    stimulus.grating.time0 = mglGetSecs;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to display stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = trialStimulusCallback(task,myscreen)

global stimulus;
mglClearScreen;

if task.thistrial.thisseg < 3
	mglLines2(-stimulus.cue.x(task.thistrial.orientNum), -stimulus.cue.y(task.thistrial.orientNum), stimulus.cue.x(task.thistrial.orientNum), stimulus.cue.y(task.thistrial.orientNum), stimulus.cue.linewidth, stimulus.cue.cueColor);
	mglFillOval(0,0,stimulus.cue.disksize,myscreen.background);
end

if task.thistrial.thisseg == 2
  if (mod(mglGetSecs(stimulus.grating.time0),1/stimulus.grating.tf) < 1/(2*stimulus.grating.tf))
    mglBltTexture(stimulus.thistex(task.thistrial.orientNum,task.thistrial.orientDiffNum),[0 0]);
  else
    mglBltTexture(stimulus.thattex(task.thistrial.orientNum,task.thistrial.orientDiffNum),[0 0]);
  end
end

mglFillOval(0,0,stimulus.fix.disksize,myscreen.background);
mglFixationCross(stimulus.fix.width,stimulus.fix.linewidth,stimulus.fix.color);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to set stimulus parameters at
% the beginning of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = trialResponseCallback(task,myscreen)

global stimulus;

if task.thistrial.thisseg == 3
%response will be either 1 for button 1 or 2 for button 2
response = 2*task.thistrial.buttonState(1)+task.thistrial.buttonState(2);

myscreen = writeTrace(response,myscreen.stimtrace+2,myscreen,1);

% update right/left counts
task.responses(end+1) = response;
task.orientationDiff(end+1) = stimulus.grating.orientationDiffs(task.thistrial.orientDiffNum);
task.orientation(end+1) = stimulus.grating.orientations(task.thistrial.orientNum);

% start by assuming incorrect
correct = 0;

% correct
if (((stimulus.grating.orientationDiffs(task.thistrial.orientDiffNum) < 0) && (response==1)) || ...
    ((stimulus.grating.orientationDiffs(task.thistrial.orientDiffNum) > 0) && (response==2)))
  correct = 1;
% ambiguous (i.e. stimulus is not left or right
elseif (stimulus.grating.orientationDiffs(task.thistrial.orientDiffNum) == 0)
  correct = 0.5;
end

%disp(sprintf('orientation: (%i) %0.2f response: %i correct: %f', task.thistrial.orientDiffNum,stimulus.grating.orientationDiffs(task.thistrial.orientDiffNum),response,correct));
disp(sprintf('orientation: %0.2f deg -- response: %i -- correct: %0.1f', stimulus.grating.orientationDiffs(task.thistrial.orientDiffNum),response,correct));

if rand <= correct
  stimulus.fix.color = stimulus.fix.correctColor;
else
  stimulus.fix.color = stimulus.fix.incorrectColor;
end
end



