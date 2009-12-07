% screensize - determining placement and size for visual stimuli
%
%      usage: [ sizes ] = screensize(  )
%         by: denis schluppeck
%       date: 2007-03-14
%        $Id: screensize.m,v 1.2 2007/03/19 18:00:07 ds Exp $:
%     inputs: 
%    outputs: sizes.r [visible radius]
%             sizes.x0y0 [xy center]
%
%    purpose: utility to measure visible area of visual stimuli, e.g.
%             the radius and center and xy extent of a circular region
%             that can been seen easily from within the MRI 
%             scanner (3T and 7T)
%
%             assumes distances, display details, etc. have been set up 
%             correctly in task/initScreen
%
%        e.g: screensize
%
function [ sizes ]=screensize(  )

% check arguments
if ~any(nargin == [0])
  help screensize
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other screen parameters
myscreen.autoCloseScreen = 1;
myscreen.allowpause = 0;
myscreen.eatkeys = 0;
%myscreen.displayname = 'projector';
%myscreen.background = 'black';

myscreen = initScreen(myscreen);

global stimulus
myscreen = initStimulus('stimulus',myscreen);
% a simple circle that has radius R and is centered at x0y0
stimulus.circle.r = 5; 
stimulus.circle.rdelta = 0.5; % radius increment 
stimulus.circle.x0y0 = [0 0];
stimulus.circle.xydelta = [.5 .5]; % xy increment

stimulus.updateFunction = @updateCircle;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up baseline task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
triallen = 1;

% set the first task to be the fixation staircase task
% [task{1} myscreen] = fixStairInitTask(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change xy
% task{1}{1}.numBlocks = 1;
task{1}{1}.seglen = 1;
task{1}{1}.getResponse = 1;
task{1}{1}.waitForBacktick = 0;

% change r
task{1}{2}.numBlocks = 1;
task{1}{2}.parameter.dir = [0;0];
task{1}{2}.parameter.coherence = [0;0];
task{1}{2}.random = 0;
task{1}{2}.seglen = 1;
task{1}{2}.waitForBacktick = 0;
task{1}{2}.timeInVols = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialze tasks and stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% stimulus = 

% initialze tasks
for tasknum = 1:length(task{1})
  task{1}{tasknum} = initTask(task{1}{tasknum},myscreen,@startSegmentCallback,@trialStimulusCallback,@trialResponseCallback);
end

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
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);
% report measurements
sizes.r = stimulus.circle.r;
sizes.x0y0 = stimulus.circle.x0y0;
mydisp(sprintf('-----------------------------\nfinal screen size\nR: %.1f deg, xy0: %.1f, %.1f\n-----------------------------\n', stimulus.circle.r, stimulus.circle.x0y0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     STIMULUS Dots start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to update dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = updateCircle(stimulus,myscreen)

% plot the circle
x0y0 = stimulus.circle.x0y0;
r = stimulus.circle.r;
mglGluDisk(x0y0(1), x0y0(2), r ,  [0.1 0.6 0], 72, 2);
% and a little black fixation marker on top
mglGluDisk(x0y0(1), x0y0(2), 0.2 ,  [0 0 0], 72, 2);



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
% set the stimulus parameters
% set direction of dots


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to display stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = trialStimulusCallback(task,myscreen)

global stimulus;
mglClearScreen;

% now update the dots, by calling update function
stimulus = feval(stimulus.updateFunction,stimulus,myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to display stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = trialResponseCallback(task,myscreen)

global stimulus;

% if the next button press would bring us out of display
max_xy = myscreen.displaySize/2;
xy = stimulus.circle.x0y0;
outofrange = ( xy(1) < -max_xy(1) ) || ( xy(1) > +max_xy(1)) || ...
    ( xy(2) < -max_xy(2) ) || ( xy(2) > +max_xy(2) );

toosmall = (stimulus.circle.r < 2);
toobig = (stimulus.circle.r > max_xy(1));

if ~outofrange
  % buttons 1-4 move in x and y
  m = [-1 0; +1 0; 0 -1; 0 +1; zeros(6,2)]; % 10 buttons
  stimulus.circle.x0y0 = stimulus.circle.x0y0 + ...
      stimulus.circle.xydelta.*m(find(task.thistrial.buttonState),:) ;
  
  % buttons 5/6 change r - add check for min R, max R
  if task.thistrial.buttonState(5)
    stimulus.circle.r = stimulus.circle.r + stimulus.circle.rdelta;
  elseif task.thistrial.buttonState(6)
    stimulus.circle.r = stimulus.circle.r - stimulus.circle.rdelta;
  end
else
  stimulus.circle.x0y0 = [0 0];
end

if toosmall,  stimulus.circle.r = 2; end
if toobig,  stimulus.circle.r = max_xy(1); end

mydisp(sprintf('r:%.1f, xy: %+.1f,%+.1f\n',stimulus.circle.r, stimulus.circle.x0y0));

% now update the circle, by calling update function
stimulus = feval(stimulus.updateFunction,stimulus,myscreen);



