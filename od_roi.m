% od
%
%      usage: myscreen=od(stimulus)
%         by: justin gardner
%       date: 04/15/06
%    purpose: orientation discrimination task
%
%
%
function myscreen = od_roi

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

stimulus.annulusIn = 2;
stimulus.annulusOut = 20;

stimulus.grating{1}.orientations = [90];
stimulus.grating{1}.orientationDiffs = [0];
stimulus.grating{1}.sf = 2;
stimulus.grating{1}.tf = 5;
stimulus.grating{1}.width = 35;
stimulus.grating{1}.height = 35;
stimulus.grating{1}.phase = 0;
stimulus.grating{1}.windowType = 'thresh'; % should be gabor or thresh
stimulus.grating{1}.sdx = stimulus.grating{1}.width/2;
stimulus.grating{1}.sdy = stimulus.grating{1}.width/2;
stimulus.grating{1}.contrast = 50;

stimulus.grating{2} = stimulus.grating{1};
stimulus.grating{2}.width = stimulus.annulusOut;
stimulus.grating{2}.height = stimulus.annulusOut;
stimulus.grating{2}.sdx = stimulus.grating{2}.width/2;
stimulus.grating{2}.sdy = stimulus.grating{2}.width/2;

stimulus.grating{3} = stimulus.grating{1};
stimulus.grating{3}.width = stimulus.annulusIn;
stimulus.grating{3}.height = stimulus.annulusIn;
stimulus.grating{3}.sdx = stimulus.grating{3}.width/2;
stimulus.grating{3}.sdy = stimulus.grating{3}.width/2;

%stimulus.dots.size = [0.3 0.3];
%stimulus.dots.eccentricity = stimulus.grating.width/2+0.5;
%stimulus.dots.color = 1;

stimulus.fix.width = 0.5;
stimulus.fix.linewidth = 1;
stimulus.fix.disksize = [2 2];
stimulus.fix.startColor = 1;
stimulus.fix.correctColor = [0 1 0];
stimulus.fix.incorrectColor = [1 0 0];
stimulus.timer = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up baseline task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%triallen = 1;

% set the first task to be the fixation staircase task
global fixStimulus;
fixStimulus.diskSize = 0;
[task{1} myscreen] = fixStairInitTask(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up segments of trials
task{2}{1}.parameter.orientNum = 1:length(stimulus.grating{1}.orientations);
task{2}{1}.parameter.orientDiffNum = 1:length(stimulus.grating{1}.orientationDiffs);
task{2}{1}.random = 1;
task{2}{1}.seglen = [15 15];
task{2}{1}.random = 1;
task{2}{1}.waitForBacktick = 1;
task{2}{1}.numTrials = 11;
task{2}{1}.writeTrace{1}.tracenum = [1 2];
task{2}{1}.writeTrace{1}.tracevar{1} = 'orientNum';
task{2}{1}.writeTrace{1}.tracevar{2} = 'orientDiffNum';
task{2}{1}.timeInVols = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialze tasks and stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stimulus = generateGrating(stimulus);

for tasknum = 1:length(task{2})
  task{2}{tasknum} = initTask(task{2}{tasknum},myscreen,@startSegmentCallback,@trialStimulusCallback);
end
% initialze tasks
%task{1} = initTask(task{1},myscreen,@startSegmentCallback,@trialStimulusCallback,@trialResponseCallback);

%task{1}.responses = [];
%task{1}.orientationDiff = [];
%task{1}.orientation = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the tasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set which phase is active
tnum = 1;

while (tnum <= length(task{2})) && ~myscreen.userHitEsc
  % updatethe task
  [task{2} myscreen tnum] = updateTask(task{2},myscreen,tnum);
  % update the fixation task
  [task{1} myscreen] = updateTask(task{1},myscreen,1);
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
%if (task.thistrial.thisseg == 2)
    stimulus.timer = mglGetSecs;
%end
task.thistrial.thisphase = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to display stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = trialStimulusCallback(task,myscreen)

global stimulus;
mglClearScreen;

if task.thistrial.thisseg == 1
  if (mod(mglGetSecs(stimulus.timer),1/stimulus.grating{1}.tf) < 1/(2*stimulus.grating{1}.tf))
      mglBltTexture(stimulus.grating{1}.thistex(task.thistrial.orientNum,task.thistrial.orientDiffNum),[0 0]);
      mglFillOval(0,0,[stimulus.annulusOut stimulus.annulusOut],myscreen.background);
      mglBltTexture(stimulus.grating{3}.thistex(task.thistrial.orientNum,task.thistrial.orientDiffNum),[0 0]);
  else
      mglBltTexture(stimulus.grating{1}.thattex(task.thistrial.orientNum,task.thistrial.orientDiffNum),[0 0]);
      mglFillOval(0,0,[stimulus.annulusOut stimulus.annulusOut],myscreen.background);
      mglBltTexture(stimulus.grating{3}.thattex(task.thistrial.orientNum,task.thistrial.orientDiffNum),[0 0]);
  end
elseif task.thistrial.thisseg == 2
  if (mod(mglGetSecs(stimulus.timer),1/stimulus.grating{1}.tf) < 1/(2*stimulus.grating{1}.tf))
      mglBltTexture(stimulus.grating{2}.thistex(task.thistrial.orientNum,task.thistrial.orientDiffNum),[0 0]);
      mglFillOval(0,0,[stimulus.annulusIn stimulus.annulusIn],myscreen.background);
  else
      mglBltTexture(stimulus.grating{2}.thattex(task.thistrial.orientNum,task.thistrial.orientDiffNum),[0 0]);
      mglFillOval(0,0,[stimulus.annulusIn stimulus.annulusIn],myscreen.background);
  end
end

% draw fixation
%mglFillOval(0,0,stimulus.fix.disksize,myscreen.background);
%mglFixationCross(stimulus.fix.width,stimulus.fix.linewidth,stimulus.fix.color);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     STIMULUS Gratings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = generateGrating(stimulus)

for grat = 1:length(stimulus.grating)
	thisStimulus = [];
	thatStimulus = [];

	disppercent(-inf,'Creating grating textures');
	% make each one of he called for gratings
	for i = 1:length(stimulus.grating{grat}.orientations)
	  for j = 1:length(stimulus.grating{grat}.orientationDiffs)
	    disppercent((i-1)/length(stimulus.grating{grat}.orientations)+j/length(stimulus.grating{grat}.orientationDiffs));
	    % get the orientation we want to create
	    thisOrientation = stimulus.grating{grat}.orientations(i)+stimulus.grating{grat}.orientationDiffs(j);
	    % make the grating
	    thisGrating = 255*(makeGrating(stimulus.grating{grat}.width,stimulus.grating{grat}.height,...
					   stimulus.grating{grat}.sf,thisOrientation,...
					   stimulus.grating{grat}.phase)+1)/2;
	
	    thatGrating = 255*(makeGrating(stimulus.grating{grat}.width,stimulus.grating{grat}.height,...
					   stimulus.grating{grat}.sf,thisOrientation,...
					   stimulus.grating{grat}.phase+180)+1)/2;
	
	    thisGaussian = makeGaussian(stimulus.grating{grat}.width,stimulus.grating{grat}.height,...
					stimulus.grating{grat}.sdx,stimulus.grating{grat}.sdy);

	    % create an rgb/a matrix
	    thisStimulus(:,:,1) = thisGrating;
	    thisStimulus(:,:,2) = thisGrating;
	    thisStimulus(:,:,3) = thisGrating;
	
	    thatStimulus(:,:,1) = thatGrating;
	    thatStimulus(:,:,2) = thatGrating;
	    thatStimulus(:,:,3) = thatGrating;
	
	    % make the gaussian window
	    if strcmp(stimulus.grating{grat}.windowType,'gabor')
	      % create an rgb/a matrix
	      mask = stimulus.grating{grat}.contrast*255*thisGaussian/100;
	    else
	      mask = stimulus.grating{grat}.contrast*255*(thisGaussian>exp(-1/2))/100;
	    end
	
	      thisStimulus(:,:,4) = mask;
	      thatStimulus(:,:,4) = mask;
	
	    % create the texture
	    stimulus.grating{grat}.thistex(i,j) = mglCreateTexture(thisStimulus);
	    stimulus.grating{grat}.thattex(i,j) = mglCreateTexture(thatStimulus);
	  end
	  % get location of dots
	%  stimulus.dots.x(i) = cos(pi*(stimulus.grating.orientations(i))/180)*stimulus.dots.eccentricity;
	%  stimulus.dots.y(i) = sin(pi*(stimulus.grating.orientations(i))/180)*stimulus.dots.eccentricity;
	end
	disppercent(inf);
end
	


