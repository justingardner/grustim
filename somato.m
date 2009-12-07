% somato - present somatosensory stimuli (cue or via NI)
%
%      usage: [  ] = somato( stimType )
%         by: denis schluppeck
%       date: 2007-05-08
%        $Id:
%     inputs: direction [inc/dec]
%    outputs: 
%
%    purpose: somatosensory stimulus presentation. while national
%    instruments interface is being built for piezo-electric device,
%    just cue visually as to what's happening...
%
% 
function [  ]=somato( direction, displayname )

% check arguments
if ~( any(nargin == [1 2]) && any(strcmp(direction,{ 'dec', 'inc', 'er'})) )
  help somato
  return
end

if nargin < 2
  myscreen.displayname = 'projector';
else
  if ischar(displayname)
    myscreen.displayname = displayname;
  else
    disp('(UHOH) displayname not set properly; reset to empty string')
    myscreen.displayname = '';
  end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other screen parameters
myscreen.autoCloseScreen = 0; % 0
myscreen.allowpause = 0;
myscreen.eatkeys = 0;
myscreen.background = 'black';

myscreen = initScreen(myscreen);
myscreen.keyboard.backtick = mglCharToKeycode({'5'});
myscreen.keyboard.nums= mglCharToKeycode({'1' '2' '3' '4' '6' '7' '8' '9'});
myscreen.saveData = 1;
myscreen.TR = 2.4;


global stimulus
myscreen = initStimulus('stimulus',myscreen);
stimulus.fixSize = 0.1;
stimulus.fixColor = [1 1 1];
stimulus.targetSize = 0.1;
stimulus.targetColor = [1 0 0];
stimulus.digits.showtext = 1; % show text?
stimulus.digits.location = [0 0]; % where to put the text

% stimulus that shows an icon moving (for subject to mimic)
stimulus.digits.driver.freq = 2; % Hz
stimulus.digits.driver.rmax = 0; % rmax
stimulus.digits.driver.loc = [0 -4]; % 

stimulus.digits.usepiezo = 0; % use piezo device (later)

stimulus = initDigits(stimulus,myscreen);

% ----------------------------------------------------
% add horizontal and vertical offset by changing the MODELVIEW matrix
% this is to make stimuli display correctly at Nott'm scanner setup
stimulus.xoffset = 0;
stimulus.yoffset = 0;
mglTransform('GL_MODELVIEW','glTranslate',stimulus.xoffset,stimulus.yoffset,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up baseline task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
triallen = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prescan
task{1}{1}.numBlocks = 1;
task{1}{1}.seglen = 0; 
task{1}{1}.getResponse = 1;
task{1}{1}.waitForBacktick = 0;

% somatotopy
% task{1}{2}.numBlocks = 1;
if strcmp(direction, 'inc')
  task{1}{2}.parameter.digit = 1:5;
elseif strcmp(direction, 'dec')
  task{1}{2}.parameter.digit = 5:-1:1;
elseif strcmp(direction, 'er')
  task{1}{2}.parameter.digit = [2 3 4];
else
  disp(sprintf('(UHOH) direction needs to be "inc" or "dec"; %s is not valid', direction));
  mglClose;
  return;
end

% Trial definition
% ------------------------------------------------------------------
% [ cue ] [ drive ]  [ x ] [ x ]
% ------------------------------------------------------------------

if any(strcmp(direction, {'inc', 'dec'}))
  % block design
  task{1}{2}.random = 0;
  % task{1}{2}.seglen = [.4 4.2 .1 .1]; % make 24/5s = 4.8s
  task{1}{2}.seglen = [1/6 1+3/4 1/24 1/24]*myscreen.TR; % make 24/5s = 4.8s
  task{1}{2}.waitForBacktick = 1;
  task{1}{2}.numBlocks = 10;
else 
  % event related
  task{1}{2}.random = 1;
  % task{1}{2}.seglen = [.4 4.2 .1 .1 ]; % make 24/5s = 4.8s
  task{1}{2}.segmin = [1/6+3/4 1 1/24 1/24 2 ]*myscreen.TR;
  task{1}{2}.segmax = [1/6+3/4 1 1/24 1/24 4 ]*myscreen.TR;
  task{1}{2}.segquant = [0     0 0    0    1 ]*myscreen.TR;
  task{1}{2}.waitForBacktick = 1;
end

% task{1}{2}.timeInVols = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialze tasks and stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialze task
% pre-scan
task{1}{1} = initTask(task{1}{1},myscreen,@pre_startSegmentCallback,@pre_trialStimulusCallback);
% saccadotopy - phase encoded version
task{1}{2} = initTask(task{1}{2},myscreen,@startSegmentCallback,@trialStimulusCallback,@trialResponseCallback);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the tasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set which phase is active
tnum = 1;

while (tnum <= length(task{1})) && ~myscreen.userHitEsc
  % updatethe task
  [task{1} myscreen tnum] = updateTask(task{1},myscreen,tnum);
  
  % update the dots (if the flag is set in that segment)
  stimulus = updateDigits(stimulus, myscreen);

  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     STIMULUS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
% pre-scan callback functions... before experiment 
% starts...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = pre_startSegmentCallback(task,myscreen)

global stimulus;
mglClearScreen;

% draw a red marker for now
mglGluDisk(0, 0, stimulus.fixSize ,  stimulus.fixColor, 72, 2);
stimulus.digits.showdriver  = 1;  
  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to display stimulus during
% pre-scan period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = pre_trialStimulusCallback(task,myscreen)

global stimulus;
mglClearScreen;

% draw a red marker for now
mglGluDisk(0, 0, stimulus.fixSize ,  stimulus.fixColor, 72, 2);
stimulus.digits.showdriver = 0;   
  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to set stimulus parameters at
% the beginning of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task,myscreen)

global stimulus;
% set the stimulus parameters
if task.thistrial.thisseg == 1
  
  %if task.blocknum == 1
  %  global MGL;
  %  mglPlaySound(find(strcmp(MGL.soundNames,'Submarine')));
  %  disp('block 1 seg 1')
  %end
  
  % reset the phase to 0, so driver starts at 0
  stimulus.digits.driver.currentPh = 0;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to display stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = trialStimulusCallback(task,myscreen)

global stimulus;
mglClearScreen;

% values for targets +- jitter
%theta = task.thistrial.theta + task.thistrial.private.thetaJitter;
%r = task.thistrial.r + task.thistrial.private.rJitter;
%x = sind(theta).*r;
%y = cosd(theta).*r;
y = stimulus.digits.driver.loc(2);

% do some display logic here:
if task.thistrial.thisseg == 1 
  % show brief cue
  mglGluDisk(0, y/2, stimulus.fixSize*3, stimulus.fixColor, 24, 2);
  mglBltTexture( stimulus.digits.tex( task.thistrial.digit ), stimulus.digits.location);
  stimulus.digits.showdriver = 0;

elseif task.thistrial.thisseg == 2
  
  % display text and sinusoidal driver stimulus
  mglGluDisk(0, y/2, stimulus.fixSize*3, stimulus.targetColor, 24, 2);
  mglBltTexture( stimulus.digits.tex( task.thistrial.digit ), stimulus.digits.location);
  stimulus.digits.showdriver = 0;

elseif task.thistrial.thisseg == 3 

  mglGluDisk(0, y/2, stimulus.fixSize*3, stimulus.fixColor, 24, 2);
  stimulus.digits.showdriver = 0;

elseif task.thistrial.thisseg == 4 

  mglGluDisk(0, y/2, stimulus.fixSize*3, stimulus.fixColor, 24, 2);
  stimulus.digits.showdriver = 0;

elseif task.thistrial.thisseg > 4 

  mglGluDisk(0, y/2, stimulus.fixSize*3, stimulus.fixColor, 24, 2);
  stimulus.digits.showdriver = 0;

  
end

% now update the dots, by calling update function
% stimulus = feval(stimulus.updateFunction,stimulus,myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to display stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = trialResponseCallback(task,myscreen)

global stimulus;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the digit stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function stimulus = initDots(stimulus,myscreen)
function stimulus = initDigits(stimulus,myscreen)

global MGL

% convert the passed in parameters to real units
if ~isfield(stimulus,'digits') || ~isfield(stimulus.digits,'fontsize'), stimulus.digits.fontsize = 48;,end
if ~isfield(stimulus,'digits') || ~isfield(stimulus.digits,'showdriver'), stimulus.digits.showdriver = 0;,end


mglClearScreen;
% set up digits
digits = {'1 - thumb', '2 - index', '3 - middle', '4 - ring', '5 - little'};
mglTextSet('Helvetica Neue',stimulus.digits.fontsize,[1 1 1],0,0,0,1,1,0,0);
for iDigit = 1:length(digits)
  mglClearScreen;
  stimulus.digits.tex(iDigit) = mglText(digits{iDigit});
  mglWaitSecs(0.5);
  mglBltTexture(stimulus.digits.tex(iDigit),[0 0]);
  mglFlush;
end
mglClearScreen;

% set up driver stimulus
stimulus.digits.driver.currentPh = 0; % phase progresses at 
% how many cycles/s: X Hz
% how many rad/cycle: 2pi
% how many rad/s: 2*pi / X Hz
stimulus.digits.driver.phaseStep = stimulus.digits.driver.freq.*2*pi./(MGL.frameRate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to update dot positions and draw them to screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = updateDigits(stimulus,myscreen)

% provide a sinusoidal driver in lieu of somatosensory stimuli for the moment
stimulus.digits.driver.currentPh = stimulus.digits.driver.currentPh+stimulus.digits.driver.phaseStep;

% phase wrap 0 < phase < 2*pi
if stimulus.digits.driver.currentPh > 2*pi
  stimulus.digits.driver.currentPh = stimulus.digits.driver.currentPh -2*pi;
elseif stimulus.digits.driver.currentPh < 0
  stimulus.digits.driver.currentPh = stimulus.digits.driver.currentPh +2*pi;
end

% now show
if stimulus.digits.showtext == 1 && stimulus.digits.showdriver    
  xdelta = stimulus.digits.driver.rmax .* sin(stimulus.digits.driver.currentPh);
  x0 = stimulus.digits.driver.loc(1); 
  x = xdelta + x0;
  y = stimulus.digits.driver.loc(2);
  mglGluDisk(0, y/2, stimulus.fixSize*3, [1 0 0], 24, 2);
  mglGluDisk(x, y, stimulus.fixSize, stimulus.fixColor, 24, 2);
else 
  % an interval in which we need to do nothing
end




