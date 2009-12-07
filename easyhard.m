% easyhard.m
%
%        $Id: easyhard.m,v 1.10 2007/02/15 22:33:33 justin Exp $
%      usage: testExperiment
%         by: denis schluppeck
%       date: 2007-01-10
%    purpose: easy hard program for experiments extending rees findings.
%             based on testExperiment and code written for dotsDualTask by justin.
%
%      tasks: (one) sang-hun task at fixation - letters and numbers came up in a rapid 
%                   sequence. subject has to indicate whether a letter  or number was the
%                   repeated item. responses (1) for number, (2) for letter
%
%             (two) task on the dots. there are two intervals. one contains motion up/left, 
%                   the other up/right. subject has to indicate interval '1', '2'.
%                   the difficulty of the task is adapted by changing the level of _coherence_.
%
%             parameters to set for fixStimulus:
%                .colors = {[1 0 0],[1 0.5 0]}; % increase
%                .letters = ['a':'b' '0':'1'];  % increase
%                .repeatLoc = 2;    % increase to 4 for sang-hun's parameters
%                .nItems = 5;       % how many items to show
%                .stimTime = 0.5;   % how long each letter/number is shown
%                .responseTime = 0.75;  % response interval after nItems have been shown
% TODO
% 
%   * green and blue colors don't show for rendered text on my machine (macintel??)
%   * fix traces to save data
%   * figure out exact timing for dual task + responses
%   * DONE dots - task: direction of motion discrimination [interval]
%   * ++ TrialResponseCallback ++ trial response callbacks get called
%     for each task; fix the code so each trialresponse callback deals
%     with the right number of responses
%
%   * staircase OK check that performance is recorded correctly

function [myscreen, task] = easyhard

global debugFlag
debugFlag = 0; 

% check arguments
if ~any(nargin == [0])
  help easyhard
  return
end

% initalize the screen
myscreen.allowpause = 1;
myscreen = initScreen(myscreen);
myscreen.sync = -1;

% decide which task to do. two task structures. both task response
% trialback functions will be updated (if a button is pressed twice,
% both trialback functions will be called twice!). determine which
% response gets priority by setting the order: 1, 2 (or NaN if you
% never want to consider the response)
sanghunTask_response = 1;
dotTask_response = 2;

% feedback off? set parameters... [ clean this up in the code! ]
% myscreen.colorClockLen = 0 ; % ... for dots
% task{1}{1}.fixColorClockLen = 0; % ... for sanghun task

% set up the stencils - draw on the same stencil?
myscreen.stencil.dots = 1;
myscreen.stencil.fixation = 1;

% ------------------------------------------------------------
% set some fixation details:
% ------------------------------------------------------------
myscreen.fix.loc = [0 0];
myscreen.fix.size = [.5 .5];
myscreen.fix.diskSize = [1 1]; % black background disk
myscreen.fix.color = [1 0 0];
myscreen.fix.linewidth = 2;
myscreen.fix.disksz = [1 1];
myscreen.background = [0 0 0];

% ------------------------------------------------------------
% a couple of parameters for feedback - which can be passed 
% to initDots
% ------------------------------------------------------------
myscreen.colorClockLen = 20;

% ------------------------------------------------------------
% set the first task to be the sang-hun task
% ------------------------------------------------------------
task = {cell(1)}; % initialize to empty
task{1}{1}.waitForBacktick = 1;
task{1}{1}.fontSize= 36;
task{1}{1}.fontName = 'Helvetica'; 
task{1}{1}.colors = {[1 0 0],[1 0.5 0],[0 1 0]}; 
task{1}{1}.letters = ['a':'d' '0':'3']; 
task{1}{1}.repeatLoc = 3; 
task{1}{1}.fixColor = [1 1 1]; 
task{1}{1}.fixColorClockLen = 30;
% timing
task{1}{1}.nItems = 10; 
task{1}{1}.stimTime = 0.5; 
task{1}{1}.responseTime = 1.5; 
% the duration of a trial here is:    nItems*stimTime + 1*responseTime
task{1}{1}.numBlocks = 50;

task{1}{1}.whichResponse = sanghunTask_response;  % which response to consider

[task myscreen] = fixSanghunInitTask(myscreen, task);

% change .seglen
% change .getResponse
% to fit with the other task?
% if you change the # of segments, check that the drawStimulusCallback
% function does the right thing.

% ------------------------------------------------------------
% second task is the task on the dots
% ------------------------------------------------------------

% set our task to have two phases. 
% one starts out with dots moving for incohrently for x  seconds
task{2}{1}.waitForBacktick = 1;
task{2}{1}.seglen = 0 ; 
task{2}{1}.numBlocks = 1;
task{2}{1}.parameter.dir = 0;
task{2}{1}.parameter.coherence = 0;

 % task{2}{1}.randVars.uniform.sigint = [1 2];

task{2}{1}.getResponse = [0];
task{2}{2}.waitForBacktick = 0;

% task on the dots - trial structure is as follows
% [pre] [interval 1] [iti] [interval2] [resp] [post]
%
% during the two intervals the coherence increases, and two directions
% of motion are presented the staircase adjusts the coherence level to
% threshold for 1up/3down for the two pre-set directions
task{2}{2}.segmin =     [1  1.0 .5 1.0 inf];
task{2}{2}.segmax =     [2  1.0 .5 1.0 inf];
task{2}{2}.getResponse = [0 0 0 0 0 0]; % response interval will be set temporarily after response to sanghun task
task{2}{2}.whichResponse = dotTask_response;  % which response to consider when multiple buttons are pressed

task{2}{2}.parameter.sigint = [1 2];
task{2}{2}.directions = {[+135 +45], [+45 +135]}; % 135 is target
task{2}{2}.parameter.coherence = [1]; % controlled by staircase
task{2}{2}.random = 0; % explicitly randomize at the beginning of each trial
% coherence threshold for a fixed difference in direction of motion
task{2}{2}.threshold = 0.5; % starting value 
task{2}{2}.stepsize = 0.05; % adjust per subject
task{2}{2}.stair = upDownStaircase(1,3, task{2}{2}.threshold, task{2}{2}.stepsize,1); % levitt rule

% and set to remember the values for dir and sigint
task{2}{2}.writetrace{1}.tracenum = 1;
task{2}{2}.writetrace{1}.tracevar{1} = 'dir';
task{2}{2}.writetrace{1}.usenum = 1;
task{2}{2}.writetrace{1}.tracenum = 2;
task{2}{2}.writetrace{1}.tracevar{1} = 'sigint';
task{2}{2}.writetrace{1}.usenum = 1;

task{2}{2}.correct = 0;
task{2}{2}.incorrect = 0;
task{2}{2}.missed = 0;

% initialize our task
% NB! there are two different @stimStartSegmentCallback for task{2}{1} !
task{2}{1} = initTask(task{2}{1},myscreen,@preStimStartSegmentCallback,@stimDrawStimulusCallback, @stimTrialResponseCallback, @endTrialCallback);
task{2}{2} = initTask(task{2}{2},myscreen,@stimStartSegmentCallback,@stimDrawStimulusCallback, @stimTrialResponseCallback, @endTrialCallback);


% init the stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);
stimulus = initDots(stimulus,myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tnum = 1;
while (tnum <= length(task{2})) && ~myscreen.userHitEsc

  % update the dots
  [task{2} myscreen tnum] = updateTask(task{2},myscreen,tnum);

  % update the fixation task
  [task{1} myscreen] = updateTask(task{1},myscreen,1);
  
  % if a response is detected in task{1}, then 'getResponse' field for
  % this segment should be set to 1. after the next response is
  % received (appropriate for the task{2} trial response callback),
  % then this is set back to zero, to make it quiet for the rest of
  % the time.
  
  if task{1}{1}.thistrial.gotResponse > 0
    % .. fast forward to last segment, so that dot task
    % feedback workd
    task{2}{tnum}.thistrial.thisseg = length(task{2}{tnum}.segmin);
    % also. set to accept responses here!
    task{2}{tnum}.getResponse(task{2}{tnum}.thistrial.thisseg) = 1;
    
  end
   
  % check for synching pulse
  if myscreen.sync
    % << fast forward by setting the seglen of this trial to 0
    task{2}{tnum}.thistrial.seglen(task{2}{tnum}.thistrial.thisseg) = 0;
    myscreen.sync = 0;
  end
  
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each TRIAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = stimStartTrialCallback(task, myscreen)

global stimulus;
% write trace for this... << FIX
myscreen.thistrial.gotResponse = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each SEGMENT 
% but only for the 0 coherence case for task{2}{1}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = preStimStartSegmentCallback(task, myscreen)

global stimulus;
% simply show 0 coherence dots.
stimulus.dots.coherence = 0;
stimulus.dots.color = [1 1 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each SEGMENT
% for task{2}{2}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = stimStartSegmentCallback(task, myscreen)

global stimulus;

% trials look as follows.
% [pre] [interval 1] [iti] [interval2] [post]

if (task.thistrial.thisseg == 1)
  mydisp(sprintf('(dottask sigint) %i\n',task.thistrial.sigint))
  
  % determine threshold coherence
  task.thistrial.threshold = task.stair.threshold;
  
  % pre interval, coherence is 0
  stimulus.dots.coherence = 0;
  % at the beginning to each trial, set color
  stimulus.dots.color = [1 1 1];
  stimulus.dots.colorClock = 0;
elseif (task.thistrial.thisseg == 2) 
  % 2ifc, interval 1
  stimulus.dots.coherence = task.stair.threshold;
  stimulus.dots.dir = task.directions{ task.thistrial.sigint }(1);
elseif (task.thistrial.thisseg == 3) 
  % time between interval1 and interval2
  stimulus.dots.coherence = 0;
  stimulus.dots.color = [1 1 1]; 
elseif (task.thistrial.thisseg == 4)  
  % 2ifc, interval 2
  % set the coherence to what the staircase has produced
  stimulus.dots.coherence = task.stair.threshold; 
  stimulus.dots.dir = task.directions{ task.thistrial.sigint }(2);
  stimulus.dots.color = [1 1 1];
elseif ( task.thistrial.thisseg > 4 )  
  % post interval, coherence is 0
  stimulus.dots.coherence = 0;
  % color could be white (or green/red depending on feedback
  % options)
  if myscreen.colorClockLen == 0
    % reset color to white
    stimulus.dots.color = [1 1 1]; 
    % otherwise let the dot drawing callback take care of color.
  end
  
end

% stimulus.dots.dir = task.thistrial.dir;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each FRAME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = stimDrawStimulusCallback(task, myscreen)

global stimulus
mglClearScreen;

% if we are in the response interval, and the color clock is running, then count down
if (stimulus.dots.colorClock >= 1 && task.thistrial.thisseg == 5) 
  stimulus.dots.colorClock = stimulus.dots.colorClock-1;
  if (stimulus.dots.colorClock <= 0)
    stimulus.dots.color = [1 1 1];
  end
end


stimulus = updateDots(stimulus,myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the end of each TRIAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = endTrialCallback(task, myscreen)

global stimulus

% syncing between the two tasks is done by the dominant task: sanghun
% letter. when we get to the end of a trial there, a syncing pulse is
% set and the seglen of the current dot trial is set to 0, hence
% quitting it.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initDots(stimulus,myscreen)

% convert the passed in parameters to real units
if ~isfield(stimulus,'dots') || ~isfield(stimulus.dots,'rmax'), stimulus.dots.rmax = 15;,end
if ~isfield(stimulus.dots,'xcenter'), stimulus.dots.xcenter = 0;,end
if ~isfield(stimulus.dots,'ycenter'), stimulus.dots.ycenter = 0;,end
if ~isfield(stimulus.dots,'dotsize'), stimulus.dots.dotsize = 3;,end
if ~isfield(stimulus.dots,'density'), stimulus.dots.density = 3;,end
if ~isfield(stimulus.dots,'coherence'), stimulus.dots.coherence = 1;,end
if ~isfield(stimulus.dots,'speed'), stimulus.dots.speed = 5;,end
if ~isfield(stimulus.dots,'dir'), stimulus.dots.dir = 0;,end
if ~isfield(stimulus.dots,'color'), stimulus.dots.color = [1 .5 .5]; end % pink
if ~isfield(stimulus.dots,'colorClock'), stimulus.dots.colorClock = 0; end % clock not set
% if ~isfield(stimulus.dots,'colorClockLen'), stimulus.dots.colorClockLen = 20; end % used to be 20.
stimulus.dots.colorClockLen = myscreen.colorClockLen;

% actually a square patch of dots that get stenciled
% so calculate width and height
stimulus.dots.width = stimulus.dots.rmax*2;
stimulus.dots.height = stimulus.dots.rmax*2;

% get the number of dots
stimulus.dots.n = round(stimulus.dots.width*stimulus.dots.height*stimulus.dots.density);

% get max and min points for dots
stimulus.dots.xmin = -stimulus.dots.width/2;
stimulus.dots.xmax = stimulus.dots.width/2;
stimulus.dots.ymin = -stimulus.dots.height/2;
stimulus.dots.ymax = stimulus.dots.height/2;

% get initial position
stimulus.dots.x = rand(1,stimulus.dots.n)*stimulus.dots.width;
stimulus.dots.y = rand(1,stimulus.dots.n)*stimulus.dots.height;

% get the step size
stimulus.dots.stepsize = stimulus.dots.speed/myscreen.framesPerSecond;

% create stencil
mglClearScreen;

mglStencilCreateBegin(myscreen.stencil.dots);
mglClearScreen;
% and draw that oval
% mglGluDisk(0,0,[1 1].*myscreen.fix.diskSize,[1 1 1],60);
mglGluDisk(stimulus.dots.xcenter,stimulus.dots.ycenter,[stimulus.dots.rmax stimulus.dots.rmax]/2,[1 1 1],60);
mglStencilCreateEnd;
mglClearScreen;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to update dot positions and draw them to screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = updateDots(stimulus,myscreen)


% get the dots step
stimulus.dots.xstep = cos(pi*stimulus.dots.dir/180)*stimulus.dots.stepsize;
stimulus.dots.ystep = sin(pi*stimulus.dots.dir/180)*stimulus.dots.stepsize;

% pick a random set of dots
stimulus.dots.coherent = rand(1,stimulus.dots.n) < stimulus.dots.coherence;

% now move those dots in the right direction
stimulus.dots.x(stimulus.dots.coherent) = stimulus.dots.x(stimulus.dots.coherent)+stimulus.dots.xstep;
stimulus.dots.y(stimulus.dots.coherent) = stimulus.dots.y(stimulus.dots.coherent)+stimulus.dots.ystep;

% randomwalk rule
thisdir = rand(1,sum(~stimulus.dots.coherent))*2*pi;
stimulus.dots.x(~stimulus.dots.coherent) = stimulus.dots.x(~stimulus.dots.coherent)+cos(thisdir)*stimulus.dots.stepsize;
stimulus.dots.y(~stimulus.dots.coherent) = stimulus.dots.y(~stimulus.dots.coherent)+sin(thisdir)*stimulus.dots.stepsize;

% movshon noise
%stimulus.dots.x(~stimulus.dots.coherent) = rand(1,sum(~stimulus.dots.coherent))*stimulus.dots.width;
%stimulus.dots.y(~stimulus.dots.coherent) = rand(1,sum(~stimulus.dots.coherent))*stimulus.dots.height;

% make sure we haven't gone off the patch
% do the dots separately for left and right hand side
stimulus.dots.x(stimulus.dots.x < stimulus.dots.xmin) = stimulus.dots.x(stimulus.dots.x < stimulus.dots.xmin)+stimulus.dots.width;
stimulus.dots.x(stimulus.dots.x > stimulus.dots.xmax) = stimulus.dots.x(stimulus.dots.x > stimulus.dots.xmax)-stimulus.dots.width;
stimulus.dots.y(stimulus.dots.y < stimulus.dots.ymin) = stimulus.dots.y(stimulus.dots.y < stimulus.dots.ymin)+stimulus.dots.height;
stimulus.dots.y(stimulus.dots.y > stimulus.dots.ymax) = stimulus.dots.y(stimulus.dots.y > stimulus.dots.ymax)-stimulus.dots.height;

% draw the dots
mglStencilSelect(myscreen.stencil.dots);
mglPoints2(stimulus.dots.x,stimulus.dots.y,stimulus.dots.dotsize, stimulus.dots.color);
mglStencilSelect(0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to get the response in the dot task...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = stimTrialResponseCallback(task,myscreen)

global stimulus

% update number of responses happens in task code.

% the response to the dot task is either the first or second,
% depending on which task is prioritized
% mydisp(sprintf('resp: %i', task.thistrial.gotResponse))

if ( task.thistrial.gotResponse == task.whichResponse-1 ) 
    
  % now that we have got the response, reset the getResponse field for this segment back to zero
  task.getResponse(task.thistrial.thisseg) = 0;
  
  % update staircase
  if isfield(task,'stair')
    task.stair = upDownStaircase(task.stair,task.thistrial.buttonState(task.thistrial.sigint));
    task.threshold = task.stair.threshold;
  end
  % keep reaction time
  task.thistrial.reactionTime = mglGetSecs-task.thistrial.segStartSeconds;
  
  if (task.thistrial.buttonState(task.thistrial.sigint))
    disp(sprintf('correct threshold=%0.2f',task.threshold));
    stimulus.dots.colorClock = stimulus.dots.colorClockLen;
    if stimulus.dots.colorClock > 0
      % it's possible that colorClockLen is set to 0 - so never turn green or red!
      stimulus.dots.color = [0 1 0];
    end  
    %myscreen = writetrace(task.threshold+1,task.taskTracenum,myscreen);
    task.correct = task.correct+1;
  else
    disp(sprintf('incorrect threshold=%0.2f',task.threshold));
    stimulus.dots.colorClock = stimulus.dots.colorClockLen;
    if stimulus.dots.colorClock > 0
      stimulus.dots.color = [1 0 0];
    end
    %myscreen = writetrace(-1,task.taskTracenum,myscreen);
    task.incorrect = task.incorrect+1;
  end  
end


