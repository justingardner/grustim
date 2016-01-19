% memorycontrast - fmri experiment, contrast response and 'remembered contrast' response
%
%      usage: [  ] = memorycontrast( subject, tasktype)
%         by: denis schluppeck / yue xing
%       date: 2009-06-12
%
%     inputs: subject   e.g. 'ds'
%             tasktype ==1: with comparison task; ==0: without
%
%  copyright: based on code from JG.
%             (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%
%
%    purpose: in this part of the study, we will look at remembered or
%    imagined stimuli  of different contrasts. need to check grating code.
%    for imagined contrast stimulus there are 2 levels: [0.3 0.7]
%
%             the task for the subject is to remember the contrast of the
%             stimulus (as given by the number during the cue condition).
%             during the comparison interval, after the delay period, the
%             subject responds by pressing 1/2 if the the comparison
%             stimulus had lower or higher contrast...  this version of the
%             code has response feedback in it.
%
%    accepts  'debug', 'ori', 'sf', 'ori' and 'sd'  with evalargs notation
%
%        e.g: memorycontrast( 'ab', 1)
%
%             memorycontrast( 'ab', 1, 'ori=[0 90]', 'sf=0.5')
%

function [  ]=cogneuro_memory( tasktype, varargin )

% check arguments
if nargin < 1
    help memorycontrast
  return
end

% evaluate the input arguments
eval(evalargs(varargin));

% setup default arguments
if ieNotDefined('debug'),debug=1;end
if ieNotDefined('sf'),sf=0.75;end
if ieNotDefined('ori'),ori=[45 135];end
if ieNotDefined('sd'),sd=0; end % if sd=0, make a circular NOT gaussian patch
if ieNotDefined('scan'),scan=0; end % if sd=0, make a circular NOT gaussian patch

% scanning params
if ieNotDefined('TR'),TR=0.5;end

% which task
if ~exist('tasktype','var')
    tasktype = 1; %with comparison task
    %tasktype =0; %no comparision task
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other screen parameters
myscreen.autoCloseScreen = 1;
myscreen.saveData = 1;
myscreen.allowpause = 0;
myscreen.eatkeys = 1;
myscreen.displayname = 'projector';
myscreen.background = 'gray';
TR = TR; % TR = 1000ms for scanner / for psychophysics set to 1.0
myscreen.tasktype = tasktype;

if scan
    myscreen = initScreen('fMRIproj32');
else
    myscreen = initScreen('VPixx');
end

myscreen.keyboard.backtick = mglCharToKeycode({'5'}); %TR = 1500ms
myscreen.keyboard.nums = mglCharToKeycode({'1' '2' '3' '4'    '6' '7' '8' '9' '0'});

% init the stimulus
global stimulus
myscreen = initStimulus('stimulus',myscreen);
%gratings
stimulus.ori = 0; % will be set in experiment
stimulus.visible =0; % controls visibility in segments
stimulus.radius = 5;%radius of patch that we want
stimulus.sd = sd;%standard deviation of grating-needed to be changed according to the radius of patch
stimulus.sf = sf;
stimulus.contrasts = 0.3;%[0.3 0.7];%[0.1 0.5 0.9];index into which contrast, changes with task
stimulus.stepsize = [0 0.4];
stimulus.cpcontrasts = [0.1 0.5 0.9]; % comparison contrasts, if we want the fixed
stimulus.outcontrasts =stimulus.contrasts+stimulus.stepsize;
stimulus.totalcontrasts = [stimulus.contrasts+stimulus.stepsize,stimulus.cpcontrasts];
stimulus.nsteps = 30; % phase steps for one complete cycle
stimulus.curPh = 0; % start phase of grating at 0
stimulus.curContrast = 0;
stimulus.textvisible = 0; %whether numbers are on the screen (changes during task)
stimulus.thisText =[];
stimulus.fontColor = [1, 0.25, 0.25];
stimulus.fontSize = 32;
stimulus.textpos = [0.5 0.5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up baseline task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
triallen = 1;

global fixStimulus
fixStimulus.responseTime = 2;
fixStimulus.stimTime = 0.6;
fixStimulus.interTime = 0.5;
fixStimulus.diskSize = 0.5;
fixStimulus.width = 0.5; %1.0;
fixStimulus.linewidth = 2;
fixStimulus.thisColor = [1 0 0];
fixStimulus.origin = [0 0];

% call initialization

% set the first task to be the fixation staircase task
% [task{1} myscreen] = fixStairInitTask(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set our task to have two phases.
% one starts out with moving grating
task{1}{1}.waitForBacktick = 0;
task{1}{1}.seglen = 0;
task{1}{1}.numBlocks = 1;
task{1}{1}.parameter.contrast = 0;
task{1}{1}.parameter.onsegs =1;
task{1}{1}.parameter.cuenum = 0;
task{1}{1}.parameter.cpcontrast = 0;
task{1}{1}.randVars.uniform.ori = [0]; % flip direction?

% the second phase

% task{2}{2}.segmin = [1 1 3]*TR;%[1 0.5]*TR;
% task{2}{2}.segmax = [1 1 5]*TR;%[1 1]*TR;
% task{2}{2}.segquant = [0 0 1]*TR;%[0 0]*TR;

task{1}{2}.parameter.contrast = stimulus.outcontrasts; % two contrast levels for now
task{1}{2}.randVars.uniform.ori = ori; % flip direction? suggested orientation
task{1}{2}.parameter.onsegs = [2 4];
task{1}{2}.parameter.cuenum = [1 2];
task{1}{2}.parameter.cpcontrast = stimulus.cpcontrasts;

switch tasktype

    case{1},
        % make event related
        % the timing here is chosen to be consisten with the Harrison &
        % Tong study on WM. Though the delay period interval is fixed,
        % this allows us to use their version of classification on the time
        % period following the randomized stimulus.
        %
        % for now we only have two contrast levels: stm1 and stim2
                           %pre st1 st2   cue re cp resp iti
        %                  1  2 3   4 5   6 7 8 9 10
        task{1}{2}.segmin = [1 1 0.2 1 0.2 1 5.6 1 1 3]*TR;
        task{1}{2}.segmax = [1 1 0.2 1 0.2 1 5.6 1 1 5]*TR;
        task{1}{2}.segquant = [0 0 0 0 0 0 0 0 0 1]*TR;
        task{1}{2}.getResponse = [0 0 0 0 0 0 0 0 1 0];

    otherwise
        task{1}{2}.segmin = [1 1 0.2 1 0.2 1 6 3]*TR;
        task{1}{2}.segmax = [1 1 0.2 1 0.2 1 6 5]*TR;
        task{1}{2}.segquant = [0 0 0 0 0 0 0 1];
end

%                    pre stim1 d1 stim2 d2 stim3 d4 cue ret comp_st response iti
%                    1 2  3  4  5  6 7   8 9 10 11 12
% task{2}{2}.segmin = [1 1 0.2 1 0.2 1 0.2 1 6 1 1 3]*TR;
% task{2}{2}.segmax = [1 1 0.2 1 0.2 1 0.2 1 6 1 1 5]*TR;
% task{2}{2}.segquant = [0 0 0 0 0 0 0 0 0 0 0 1]*TR;
% task{2}{2}.getResponse = [0 0 0 0 0 0 0 0 0 0 1 0]*TR;

task{1}{2}.random = 1;
task{1}{2}.numTrials = 40; % number of trials to go through


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialze tasks and stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stimulus = myInitStimulus(stimulus,task,myscreen);
stimulus = initGratings(stimulus,myscreen,task);


%initialize tasks
for phaseNum = 1:length(task{1})
  [task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@trialStimulusCallback,@fixResponseCallback);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
  % update the dots
  [task{1} myscreen phaseNum] = updateTask(task{1},myscreen,phaseNum);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     STIMULUS Gratings start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initGratings(stimulus,myscreen,task)

% actually a square patch of dots that get stenciled
% so calculate width and height

    %mglMakeGrating(width,height,sf,angle,phase,<xDeg2pix>,<yDeg2pix>)

if myscreen.tasktype ==0
    ncontrasts = length(stimulus.outcontrasts);%how many contrasts
else
    ncontrasts = length(stimulus.totalcontrasts);
end

% set parameters for text

mglTextSet('Helvetica',stimulus.fontSize,stimulus.fontColor,0,0,0);%show the number for each corresponding coherence levels

 for j=1:ncontrasts
  %fprintf('now creating stimuli for contrast %.2f\n',stimulus.contrasts(j));%this means the sequence of showing contrast is same for each trial

  for i=1:stimulus.nsteps
    stimulus.grating = mglMakeGrating(stimulus.radius.*2,stimulus.radius.*2,stimulus.sf,0,i*360/stimulus.nsteps);
     if stimulus.sd > 0
      stimulus.gaussian = mglMakeGaussian(stimulus.radius.*2,stimulus.radius.*2,stimulus.sd, stimulus.sd);
    else
      stimulus.gaussian = mglMakeGaussian(stimulus.radius.*2,stimulus.radius.*2,stimulus.radius/3, stimulus.radius/3);
      stimulus.gaussian = stimulus.gaussian > 0.001;
    end

    if ismac
      stimulus.grating = 255*(stimulus.grating+1)/2;
      stimulus.g(:,:,1)=stimulus.grating;
      stimulus.g(:,:,2)=stimulus.grating;
      stimulus.g(:,:,3)=stimulus.grating;
      curContrast = stimulus.totalcontrasts(j);%curContrast always get the current contrast: 0.2(1) or 0.9(2)
      stimulus.g(:,:,4)=255*stimulus.gaussian*curContrast;
    else
      stimulus.g=stimulus.grating*128.*stimulus.gaussian+127;
    end
    % now create the texture
    stimulus.tex(i,j) = mglCreateTexture(stimulus.g);%stimulus.tex has 2 changing parameters: i:phase j:contrast
  end

  % for each contrast level, make some text to go along.

  thisText = mglText(sprintf('%d',j));%thisText: sequence is always 1=>2
  stimulus.numbers(j) = thisText;%stimulus.numbers is always 1=>2
end


% display each texture to back buffer, to make sure
% everything is cached properly. this may not be necessary... but doesn't
% cost much time at this point.
for i = 1:stimulus.nsteps
  for j=1:ncontrasts
    mglBltTexture(stimulus.tex(i,j),[0 0]);
  end
end


% assign appropriate update function
stimulus.updateFunction = 'updateGratings';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to update grating and draw them to screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = updateGratings(stimulus,myscreen)
global fixStimulus
% change the image # (phase step in animation...) + draw.
% for now we just draw the image.

% figure out the current phase that's appropriate for the current frame:

% keyboard
% myscreen.time gives current timestamp.
% if we start with frame 1
stimulus.curPh = mod(stimulus.curPh+1, stimulus.nsteps-1)+1;

thisPhase = stimulus.curPh;
if myscreen.tasktype==0
    thisContrast = find(stimulus.outcontrasts == stimulus.curContrast);%localize curContrast in stimulus.contrasts
%   thisContrast is the number of the current-showing contrast
else
    thisContrast =find(stimulus.totalcontrasts == stimulus.curContrast);
end

mglClearScreen;

if stimulus.visible
    %calculate next phase step to display
    mglBltTexture(stimulus.tex(thisPhase, thisContrast),[0 0],0,0,stimulus.ori);%draw the grating with phase and contrast
end


if stimulus.textvisible
    mglBltTexture(stimulus.thisText,stimulus.textpos,0,0,0);
end

if stimulus.textvisible && stimulus.visible
   mglBltTexture(stimulus.numbers(thisContrast), stimulus.textpos,0,0,0);
end


mglFixationCross(fixStimulus.width,fixStimulus.linewidth,fixStimulus.thisColor);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)
% switching on /off visibility of stimuli, text, etc.

global stimulus fixStimulus

% for now we only have two contrast levels: stm1 and stim2
                   %pre  st1  st2   cue re cp resp iti
%                   1    2 3  4 5   6   7  8    9   10
%
% this timing sequence looks complicated, but we want to display two
% stimuli at the start st1 and st2, followed by a numeric cue.. so segments
% 2/3 deal with the first, 4/5 with the second grating, segm

if any(task.thistrial.thisseg == [1 9 10])
    % reset fixation color.
    fixStimulus.thisColor = [1 0 0];
end

if any(task.thistrial.thisseg == task.thistrial.onsegs)
    stimulus.visible = 1;
    stimulus.thisText =[];
    stimulus.ori = task.thistrial.ori;
    stimulus.curContrast = task.thistrial.contrast;
    stimulus.textvisible =1;
    fprintf('thisseg: %d, ori: %d, contrast: %.2f\n',task.thistrial.thisseg, stimulus.ori, stimulus.curContrast);

elseif ( task.thistrial.thisseg == 2 && task.thistrial.onsegs == 4 ) ...
        || ( task.thistrial.thisseg == 4 && task.thistrial.onsegs == 2 )
    % current segement is 2, but the stimulus subject will have to remember
    % is in 4 or vice versa... decide on the other stimulus contrast.

    stimulus.visible = 1;
    stimulus.thisText =[];
    stimulus.ori= task.thistrial.ori;
    stimulus.textvisible =1;
    % if the current contrast is bigger than mid (0.5)... should fix this
    % and remove "magic number" here and in next branching statement...
    if task.thistrial.contrast > 0.5
        stimulus.curContrast = stimulus.contrasts; %0.3
        fprintf('thisseg: %d, ori: %d, contrast: %.2f\n',task.thistrial.thisseg,stimulus.ori, stimulus.curContrast);
    else
        stimulus.curContrast = stimulus.contrasts + stimulus.stepsize(find(stimulus.stepsize~=0));
        fprintf('thisseg: %d,ori: %d, contrast: %.2f\n',task.thistrial.thisseg,stimulus.ori, stimulus.curContrast);
    end

elseif task.thistrial.thisseg ==6
    % cue / text only
     stimulus.thisText=mglText(sprintf('%d',task.thistrial.cuenum));
	 stimulus.textvisible =1;
	 stimulus.visible =0;

elseif task.thistrial.thisseg == 8
  % comparison stimulus.

  stimulus.textvisible =0;
  if myscreen.tasktype ==1
    stimulus.visible =1;
    stimulus.ori = task.thistrial.ori;
    stimulus.curContrast = task.thistrial.cpcontrast;
    fprintf('thisseg: %d,ori: %d, contrast: %.2f\n',task.thistrial.thisseg,stimulus.ori, stimulus.curContrast);
  else
    stimulus.visible =0;
  end
else
  stimulus.visible = 0;
  stimulus.textvisible =0;
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to display stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = trialStimulusCallback(task,myscreen)

global stimulus;
mglClearScreen;
stimulus = feval(stimulus.updateFunction,stimulus,myscreen);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stimulus task myscreen] = myInitStimulus(stimulus,task,myscreen)

if isfield(stimulus,'gratings')
   stimulus = initGratings(stimulus,myscreen,task); % normal gratings
    stimulus.updateFunction = @updateGratings;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that change the fixation after comparison between the response
% and onsegs(signals on)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = fixResponseCallback(task,myscreen)

global  stimulus fixStimulus


% a button has been pressed, set the responded flag
stimulus.responded = 1;

% collect observer response and display it
%task.getResponse(task.thistrial.thisseg==11);
disp(sprintf('task.thistrial.whichButton %d', task.thistrial.whichButton));

if task.thistrial.cuenum == 1 %1 is always the smaller contrast, 0.3 in the example case
    if task.thistrial.whichButton ==1 && task.thistrial.cpcontrast > (stimulus.contrasts+min(stimulus.stepsize))
        response = 0;
    elseif task.thistrial.whichButton ==2 && task.thistrial.cpcontrast <(stimulus.contrasts+min(stimulus.stepsize))
        response =0;
    else
        response =1;
    end
else
    if task.thistrial.whichButton ==1 && task.thistrial.cpcontrast < (stimulus.contrasts+min(stimulus.stepsize))
    response = 1;
elseif task.thistrial.whichButton ==2 && task.thistrial.cpcontrast >(stimulus.contrasts+min(stimulus.stepsize))
    response =1;
    else
    % incorrect
    response =0;
    end
end



if response == 1
  disp('correct response');
  % set to correct fixation color
  fixStimulus.thisColor = [0 1 0];
else
  disp('incorrect response');
  % set to incorrect fixation color
  fixStimulus.thisColor = [1 0 0];
end

fprintf('fixStimulus.thisColor %.2f, %.2f, %.2f \n', fixStimulus.thisColor);


end