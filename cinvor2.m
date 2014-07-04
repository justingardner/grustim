% cinvor.m
%
%        $Id: cinvor2.m 835 2010-06-29 04:04:08Z justin $
%      usage: cinvor2
%         by: justin gardner, modified by taosheng liu
%       date: 09/24/2013
%  copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%    purpose: program to generate stimuli for Contrast-INVariant-ORientation experiment
%
function myscreen = cinvor2

% check arguments
if ~any(nargin == [0])
  help cinvor2
  return
end

% create

% stimulus parameters
global stimulus;
stimulus.contrast = [0.2 0.8];
stimulus.orientation = linspace(0,157.5,8); %[0:180/12:179];
stimulus.sf = 0.7;
stimulus.phaseFreq = 5;
stimulus.size = 10; 
stimulus.eccentricity = 8;

% initalize the screen
myscreen.background = 'gray';
myscreen = initScreen(myscreen);

% init stimulus
stimulus = myInitStimulus(stimulus,myscreen);
% and register with myscreen
myscreen = initStimulus('stimulus',myscreen);

inmagnet = 1;
framePeriod = 1.28;
stimDur=4;
trialDur=7*framePeriod;

task{1}{1}.waitForBacktick = inmagnet;
% first phase is a single trial that is just run to be discarded
task{1}{1}.seglen = 4*framePeriod-0.1*inmagnet;  %TSL: shouldn't this be fixed 4 TRs? (numtrial vs. numblock)
task{1}{1}.numTrials = 1;
if inmagnet
    task{1}{1}.synchToVol = 1; 
end

task{1}{2}.seglen = [stimDur trialDur-stimDur-0.1*inmagnet];
task{1}{2}.parameter.contrast = repmat(stimulus.contrast,2,1);
task{1}{2}.parameter.orientation = repmat(stimulus.orientation,2,1);
% task{1}{2}.numTrials=2*length(stimulus.contrast)*length(stimulus.orientation);
task{1}{2}.numBlocks=3;
task{1}{2}.random = 1;
if inmagnet
    task{1}{2}.synchToVol =[0 1];
end

% initialize the main task
for taskNum = 1:length(task)
  for phaseNum = 1:length(task{taskNum})
    [task{taskNum}{phaseNum} myscreen] = initTask(task{taskNum}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback);
    % remember taskID
%     stimulus.taskID(taskNum,phaseNum) = task{taskNum}{phaseNum}.taskID; %TSL: what's taskID?
  end
end

% init the fixation task
global fixStimulus
fixStimulus.verbose = 0;
[task{2} myscreen] = fixStairInitTask(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1; %[1 1 1 1];
% do loop until all locations have finished
while phaseNum <= length(task{1}) && ~myscreen.userHitEsc
  % clear screen
  mglClearScreen; %TSL: clear screen here?
  % update the task
  [task{1} myscreen phaseNum(1)] = updateTask(task{1},myscreen,phaseNum(1));
  % update the fixtation task
  [task{2} myscreen] = updateTask(task{2},myscreen,1);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

% delete textures / etc
endStimulus;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;

% init phase to random
if task.thistrial.thisphase==2 && task.thistrial.thisseg == 1

  % display trial info
%   disp(sprintf('Trial: %i OriLeft: %0.1f ContrLeft: %0.2f OriRight: %0.1f ContrRight: %0.2f  (%s)',task.trialnum,task.thistrial.orientation(1),task.thistrial.contrast(1),task.thistrial.orientation(2),task.thistrial.contrast(2),num2str(task.thistrial.seglen,'%0.2f ')));
  fprintf('Trial: %i OriLeft: %0.1f ContrLeft: %0.2f OriRight: %0.1f ContrRight: %0.2f  (%s)\n',task.trialnum,task.thistrial.orientation(1),task.thistrial.contrast(1),task.thistrial.orientation(2),task.thistrial.contrast(2),num2str(task.thistrial.seglen,'%0.2f '));
  stimulus.phase = (randi(16,1,2)-1)/16*360;  %16 possible discreate phases uniformly spaced between 0 and 360
%   stimulus.phaseClock = mglGetSecs-rand(1,2)*stimulus.phaseUpdateTime;
  stimulus.phaseClockInit = mglGetSecs-rand(1,2)*stimulus.phaseUpdateTime;
  stimulus.phaseClock=stimulus.phaseClockInit;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)
global stimulus

if task.thistrial.thisphase==2 && task.thistrial.thisseg == 1 %task.thistrial.stimseg
    % get x and y for this position
    x = stimulus.x;
    y = stimulus.y;
    
    % get contrast num
    %   contrastNum = find(stimulus.contrast == task.thistrial.contrast);
    [temp, contrastNum] = ismember(task.thistrial.contrast, stimulus.contrast);
    
    for i=1:2
        % and choose a phase
        if mglGetSecs(stimulus.phaseClock(i)) > stimulus.phaseUpdateTime
            % reset the phase
            stimulus.phase(i) = (randi(16)-1)/16*360;
            % and restart phase clock (reset to actual time remaining till enxt phase change)
            %     stimulus.phaseClock(task.thistrial.loc) = 2*mglGetSecs-(stimulus.phaseClock(task.thistrial.loc)+stimulus.phaseUpdateTime);
            stimulus.phaseClock(i) = stimulus.phaseClock(i)+stimulus.phaseUpdateTime;
        end
        %   phase = stimulus.phase(task.thistrial.loc);
        
        % calculate x/y shift needed for phase shift
        xPhaseShift = (1/stimulus.sf)*((stimulus.phase(i)-180)/360)*cos(d2r(task.thistrial.orientation));
        yPhaseShift = (1/stimulus.sf)*((stimulus.phase(i)-180)/360)*sin(d2r(task.thistrial.orientation));
        
        mglBltTexture(stimulus.tex1dsin(contrastNum(i)),[x(i)+xPhaseShift(i) y(i)+yPhaseShift(i) stimulus.texsize stimulus.texsize],0,0,task.thistrial.orientation(i));
        % blt window
        mglBltTexture(stimulus.tex2dGaussWin,[x(i) y(i) stimulus.masksize stimulus.masksize],0,0,task.thistrial.orientation(i)); %TSL: is it better to blt twice, instead of once?
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen)

stimulus.x=[-stimulus.eccentricity, stimulus.eccentricity];
stimulus.y= [0 0];
stimulus.gray = myscreen.grayIndex;
% convert to pixels
texSizePixels = round(mglGetParam('xDeviceToPixels')*stimulus.size);

% make the sine wave grating longer than needed, so that we can move it underneath the gaussian window
stimulus.texsize = stimulus.size+(1/stimulus.sf);
numCycles = stimulus.sf*stimulus.texsize;
for iContrast = 1:length(stimulus.contrast)
  stimulus.tex1dsin(iContrast) = mglCreateTexture(255*(stimulus.contrast(iContrast)*sin(0:numCycles*2*pi/(1.5*texSizePixels-1):numCycles*2*pi)+1)/2);
end

% make a gaussian window
stimulus.masksize = stimulus.texsize+(1/stimulus.sf);
[win xWin yWin] = mglMakeGaussian(stimulus.masksize,stimulus.masksize,stimulus.size/7,stimulus.size/7);   %TSL: why generate Gaussian window but don't use it?

% now threshold out at the right size (comment out for gabor)
win = sqrt((xWin.^2)+(yWin.^2)) < (stimulus.size/2);
m(:,:,4) = 255*(1-win);
m(:,:,1) = stimulus.gray;
m(:,:,2) = stimulus.gray;
m(:,:,3) = stimulus.gray;
stimulus.tex2dGaussWin = mglCreateTexture(m);

% clear screen
mglClearScreen(stimulus.gray);mglFlush;
mglClearScreen(stimulus.gray);mglFlush;

% get the time after which the phase needs to be updated
stimulus.phaseUpdateTime = 1/stimulus.phaseFreq;

return
% test code - display a stimulus with these parameters
phase = 270;
orientation = 45;
contrast = 0.8;
x = 4;y = 4;

% calculate x/y shift needed for phase shift
xPhaseShift = (1/stimulus.sf)*((phase-180)/360)*cos(d2r(orientation));
yPhaseShift = (1/stimulus.sf)*((phase-180)/360)*sin(d2r(orientation));

% get matching contrastNum
contrastNum = find(stimulus.contrast==contrast);

% blt grating texture 
mglBltTexture(stimulus.tex1dsin(contrastNum),[x+xPhaseShift y+yPhaseShift stimulus.texsize stimulus.texsize],0,0,orientation);

% blt window
mglBltTexture(stimulus.tex2dGaussWin,[x y stimulus.masksize stimulus.masksize],0,0,orientation);

% flush
mglFlush;
 
keyboard

%%%%%%%%%%%%%%%%%%%%%
%    endStimulus    %
%%%%%%%%%%%%%%%%%%%%%
function endStimulus %TSL: why are you doing this?

global stimulus
if isfield(stimulus,'tex1dsin') 
  for iTex = 1:length(stimulus.tex1dsin)
    mglDeleteTexture(stimulus.tex1dsin(iTex));
  end
  stimulus = rmfield(stimulus,'tex1dsin');
end

if isfield(stimulus,'tex2dGaussWin') 
  for iTex = 1:length(stimulus.tex2dGaussWin)
    mglDeleteTexture(stimulus.tex2dGaussWin(iTex));
  end
  stimulus = rmfield(stimulus,'tex2dGaussWin');
end
