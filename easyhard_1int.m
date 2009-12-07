% easyhard_1int.m
%
%        $Id: easyhard_1int.m,v 1.9 2007/06/26 13:38:59 ds Exp $
%      usage: easyhard_1int( whichExperiment, fixed_coh, stair_coh, stair_dir )
%         by: denis schluppeck
%       date: 2007-01-10
%    purpose: easy hard program for experiments extending rees findings.
%             based on testExperiment and code written for dotsDualTask by justin.
%
%      tasks: (one) sang-hun task at fixation - letters and numbers came up in a rapid 
%                   sequence. subject has to indicate whether a letter  or number was the
%                   repeated item. responses (1) for number, (2) for letter
%
%             (two) task on the dots. there is one interval. 
%                   - get coherence threshold for large delta-dir (2-task)
%                   - get coherence threshold for small delta-dir 
%                   - get dir threshold for fixed coherence (2-task)
%
%           e.g.: easyhard_1int( 4, 0.5, [], [18 5 1 2])  % scanner 
%                 % dir staircase (start at 18, step 5, 1up/2down) at coh=0.5
%
%                        (priority task)     |    2nd task ( if any )
%   ------------------------------------------------------------------------
%   whichExperiment: 1 - sanghun task; coherence thresh for delta-dir=120 
%                        args: [], stair_coh, []
%
%                    2 - coherence thresh for delta-dir=120; sh task
%                        args: [], stair_coh, []
%
%                    3 - delta-dir thresh for fixed coherence; sh task
%                        args: fixed_coh, [], stair_dir
%
%            scanner 4 - staircase dir (fine), fixed coherence; sh task
%                        args: fixed_coh, [], stair_dir
%
%            scanner 5 - sanghun task; fixed coherence, dir (passive)
%                        args: fixed_coh, [], []

function [myscreen, task] = easyhard_1int( whichExperiment, fixed_coh, stair_coh, stair_dir)

global debugFlag
debugFlag = 0; 

% check arguments + inputs
if ~any(nargin == [1 2 3 4])
  help easyhard_1int
  return
end

if nargin ==  1; error('please specify which experiment'); end
switch whichExperiment,
 case {1, 2},
  %  1 - sanghun task; coherence thresh for delta-dir=120 
  %  args: [], stair_coh, []
  %  2 - coherence thresh for delta-dir=120; sh task
  %  args: [], stair_coh, []
  if  numel(stair_coh) ~= 4 , error('(UHOH) check inputs'); end
 case {3, 4, 6},
  %  3 - delta-dir thresh for fixed coherence; sh task
  %      args: fixed_coh, [], stair_dir
  %  4 [scanner] - staircase dir (fine), fixed coherence; sh task
  %      args: fixed_coh, [], stair_dir
  if (any( numel(fixed_coh)==[0]) ) || ( numel(stair_dir) ~= 4 ), error('(UHOH) check inputs'); end
 case 5,
  %  5 [scanner] - sanghun task; fixed coherence, dir (passive)
  if ( any(numel(fixed_coh)==[0])), error('(UHOH) check inputs'); end
 otherwise
  error('(UHOH) this experiment does not exist!')
end
  
% initalize the screen
myscreen = initScreen;
myscreen.sync = -1;

myscreen.TR = 1; % TR in seconds
% ------------------------------------------------------------------------
§%  which experiment to run.
% ------------------------------------------------------------------------
%
%       *_response variables
%
% decide which task to do. two task structures. both task response
% trialback functions will be updated (if a button is pressed twice,
% both trialback functions will be called twice!). determine which
% response gets priority by setting the order: 1, 2 (or NaN if you
% never want to consider the response)

% fixed_coherence = 0.5; % for whichExperiment 3,4,5
% start_delta_dir = 10; % for whichExperiment 3,4: starting point for
                      % threshold staircase
% start_coherence = 0.5; % for whichExperiment 3,4: starting point for
                      % threshold staircase
		      
task{1}{1}.itiMinTime = 0;
task{1}{1}.itiMaxTime = 0;
% display dots in isi interval
myscreen.isiDots = 1;

if whichExperiment == 1
  %  1 - sanghun task; coherence thresh for delta-dir=120 
  exptype = 'coherence_thresholds';
  sanghunTask_response = 1;
  dotTask_response = 2;
  myscreen.fixed_coherence = [];
  % feedback
  sanghun_colorClockLen = 20;
  dot_colorClockLen = 20;
  
elseif whichExperiment == 2
  %  2 - coherence thresh for delta-dir=120; sh task
  exptype = 'coherence_thresholds';  
  sanghunTask_response = 2;
  dotTask_response = 1;
  myscreen.fixed_coherence = [];
  % feedback
  saghun_colorClockLen = 20;
  dot_colorClockLen = 20;
  
elseif whichExperiment == 3
  %  3 - delta-dir thresh for fixed coherence; sh task
  exptype = 'fine_direction_thresholds';  
  sanghunTask_response = 2;
  dotTask_response = 1;
  myscreen.fixed_coherence = fixed_coh;
  % feedback
  sanghun_colorClockLen = 20;
  dot_colorClockLen = 20;
  
elseif whichExperiment == 4
  %  scanner 4 - staircase dir (fine), fixed coherence; sh task (passive)
  exptype = 'fine_direction_thresholds';  
  sanghunTask_response = nan;
  dotTask_response = 1;
  myscreen.fixed_coherence = fixed_coh;
  % no feedback
  sanghun_colorClockLen = 0;
  dot_colorClockLen = 20;
  task{1}{1}.itiMinTime = 1*myscreen.TR;
  task{1}{1}.itiMaxTime = 3*myscreen.TR;

elseif whichExperiment == 5
  %  scanner 5 - sanghun task; fixed coherence, dir (passive)
  exptype = 'fixed_params';  
  sanghunTask_response = 1;
  dotTask_response = nan;
  myscreen.fixed_coherence = fixed_coh;
  % no feedback
  sanghun_colorClockLen = 20;
  dot_colorClockLen = 0;
  task{1}{1}.itiMinTime = 1*myscreen.TR;
  task{1}{1}.itiMaxTime = 3*myscreen.TR;
elseif whichExperiment == 6
  %  scanner 6 - staircase dir (fine), fixed coherence; sh task (passive) with no dots in isi)
  exptype = 'fine_direction_thresholds';  
  sanghunTask_response = nan;
  dotTask_response = 1;
  myscreen.fixed_coherence = fixed_coh;
  % no feedback
  sanghun_colorClockLen = 0;
  dot_colorClockLen = 20;
  task{1}{1}.itiMinTime = 1*myscreen.TR;
  task{1}{1}.itiMaxTime = 3*myscreen.TR;
  myscreen.isiDots = 0;
else
 error('(UHOH) whichExperiment not specified')
end

% make visible for stimulus callbacks
myscreen.exptype = exptype; 
myscreen.whichExperiment = whichExperiment;
myscreen.sigpresent = -1;

% feedback off? set parameters... [ clean this up in the code! ]
% myscreen.colorClockLen = 0 ; % ... for dots
% task{1}{1}.fixColorClockLen = 0; % ... for sanghun task
myscreen.colorClockLen = dot_colorClockLen;

% set up the stencils - draw on the same stencil?
myscreen.stencil.dots = 1;
myscreen.stencil.fixation = 1;

% ------------------------------------------------------------
% set some fixation details:
% ------------------------------------------------------------
myscreen.fix.loc = [0 0];
myscreen.fix.size = [.5 .5];
myscreen.fix.diskSize = [2.5 2.5]; % black background disk
myscreen.fix.color = [1 0 0];
myscreen.fix.linewidth = 2;
myscreen.fix.disksz = [1 1];
myscreen.background = [0 0 0];

% ------------------------------------------------------------
% set the first task to be the sang-hun task (DOMINANT TASK)
% ------------------------------------------------------------
task{1}{1}.waitForBacktick = 1;
task{1}{1}.fontSize= 36;
task{1}{1}.fontName = 'Helvetica'; 
task{1}{1}.colors = {[1 0 0],[1 0.5 0],[0 1 0]}; 
task{1}{1}.letters = ['a':'d' '0':'3']; 
task{1}{1}.repeatLoc = 3; 
task{1}{1}.fixColor = [1 1 1]; 
task{1}{1}.fixColorClockLen = sanghun_colorClockLen;

% timing
task{1}{1}.nItems = 10; 
task{1}{1}.stimTime = myscreen.TR/2;  % approx. 500ms 
task{1}{1}.responseTime = myscreen.TR; 
% the duration of a trial here is:    nItems*stimTime + 1*responseTime
task{1}{1}.numBlocks = 50;         % <<<< change this

task{1}{1}.whichResponse = sanghunTask_response;  % which response to consider

[task myscreen] = fixSanghunInitTask(myscreen, task);

% ------------------------------------------------------------
% second task is the task on the dots (SUBSIDIARY TASK)
% ------------------------------------------------------------

% set our task to have two phases. 
% one starts out with dots moving for incohrently for x  seconds
%task{2}{1}.waitForBacktick = 1;
%task{2}{1}.seglen = 10*myscreen.TR; 
%task{2}{1}.numBlocks = 1;
%task{2}{1}.parameter.dir = 0;
%task{2}{1}.parameter.coherence = 0;
%task{2}{1}.getResponse = [0];

% second phase - the actual task
task{2}{1}.waitForBacktick = 1;

% task on the dots - trial structure is as follows
% [pre] [interval 1] [post which can accept a response]
%
task{2}{1}.segmin =      [inf 0  1  inf].*myscreen.TR;
task{2}{1}.segmax =      [inf 3  1  inf].*myscreen.TR;
task{2}{1}.segquant =    [inf 1  0  0 ].*myscreen.TR;
task{2}{1}.getResponse = [inf 0  0  0 ]; % response interval will be set temporarily after response to sanghun task
task{2}{1}.whichResponse = dotTask_response;  % which response to consider when multiple buttons are pressed
task{2}{1}.randVars.uniform{1}.sigpresent = [ 0 1 ];
task{2}{1}.directions = [+150 +30]; % 145 is target
task{2}{1}.randVars.uniform{2}.pedestaldir = task{2}{1}.directions;
task{2}{1}.randVars.uniform{2}.blocklen_ = 10;
% ------------------------------------------------------------
% settings for different tasks
% ------------------------------------------------------------

if strcmp(myscreen.exptype, 'coherence_thresholds')
  task{2}{1}.parameter.coherence = [1]; % controlled by staircase...
  % coherence threshold for a fixed difference in direction of motion
  task{2}{1}.threshold = stair_coh(1); % starting value 
  task{2}{1}.stepsize = stair_coh(2); % adjust per subject
  task{2}{1}.stair{1} = upDownStaircase(stair_coh(3),stair_coh(4), task{2}{1}.threshold, ...
				     task{2}{1}.stepsize,1); % levitt rule

elseif strcmp(myscreen.exptype, 'fine_direction_thresholds')
  % direction threshold for a fixed coherence
  task{2}{1}.parameter.coherence = myscreen.fixed_coherence; % determined by pilot experiment
  task{2}{1}.threshold = stair_dir(1); % starting value, direction difference (degrees)
  task{2}{1}.stepsize = stair_dir(2); % adjust per subject
  for i = 1:length(myscreen.fixed_coherence)
    task{2}{1}.stair{i} = upDownStaircase(stair_dir(3),stair_dir(4), task{2}{1}.threshold, ...
				       task{2}{1}.stepsize,1); % levitt rule
  end
elseif strcmp(myscreen.exptype, 'fixed_params')
  task{2}{1}.parameter.coherence = myscreen.fixed_coherence; % determined by pilot experiment
  task{2}{1}.threshold = 0; % starting value, direction difference (degrees)
  task{2}{1}.stepsize = 0; % adjust per subject
  % dummy staircase - does not really get used
  for i = 1:length(myscreen.fixed_coherence)
    task{2}{1}.stair{i} = upDownStaircase(1,2, task{2}{1}.threshold, ...
				     task{2}{1}.stepsize,1); % levitt rule
  end
else
  error('(UHOH) mistyped exptype?');
end

% pass the directions to myscreen, so dot stimuli update functions
% know about the directions 
myscreen.directions = unique( task{2}{1}.directions ); 
myscreen.pedestaldir = nan; % only initially
  
task{2}{1}.correct = 0;
task{2}{1}.incorrect = 0;
task{2}{1}.missed = 0;

% initialize our task
% NB! there are two different @stimStartSegmentCallback for task{2}{1} !
[task{2}{1} myscreen] = initTask(task{2}{1},myscreen,@preStimStartSegmentCallback,@stimDrawStimulusCallback, @stimTrialResponseCallback, @endTrialCallback);
[task{2}{1} myscreen] = initTask(task{2}{1},myscreen,@stimStartSegmentCallback,@stimDrawStimulusCallback, @stimTrialResponseCallback, @endTrialCallback);


% init the stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);
stimulus.dots.rmax = 10; % stimulus size. xcenter, ycenter are 0
stimulus.dots.speed = 8;

stimulus.startSubsidiary = 0;
stimulus.endSubsidiary = 0;

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

  
  % update the fixation task
  [task{1} myscreen] = updateTask(task{1},myscreen,1);
  
  % update the dots - SUBSIDIARY TASK
  [task{2} myscreen tnum] = updateTask(task{2},myscreen,tnum);

  
  % if a response is detected in task{1}, then 'getResponse' field for
  % this segment should be set to 1. after the next response is
  % received (appropriate for the task{2} trial response callback),
  % then this is set back to zero, to make it quiet for the rest of
  % the time.
  
  if task{1}{1}.thistrial.gotResponse > 0
    % also. set to accept responses here!
    task{2}{tnum}.getResponse(task{2}{tnum}.thistrial.thisseg) = 1;
    % mydisp(sprintf('thistrial.thisseg [task 1 2] %i %i\n',task{1}{1}.thistrial.thisseg, task{2}{1}.thistrial.thisseg ));
  end
   
  % check for synching pulse
  %if myscreen.sync
  %  % << fast forward by setting the seglen of this trial to 0
  %  task{2}{tnum}.thistrial.seglen(task{2}{tnum}.thistrial.thisseg) = 0;
  %  myscreen.sync = 0;
  %end
  
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
% but only for the 0 coherence case for task  (OBSOLETE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = preStimStartSegmentCallback(task, myscreen)

global stimulus;
% simply show 0 coherence dots.
stimulus.dots.coherence = 0;
stimulus.dots.color = [1 1 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each SEGMENT
% for task{2}{1}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = stimStartSegmentCallback(task, myscreen)

global stimulus;

% trials look as follows.
% [pre] [interval 1] [resp] [post]

stairnum = find(task.thistrial.coherence == task.parameter.coherence);

if (task.thistrial.thisseg == 1)
  disp('(XXXX) Seg 1 for subsidiary task to start')

elseif (task.thistrial.thisseg == 2)
  disp('(XXXX) Seg 2 subsidiary task started')
  
  mydisp(sprintf('(dottask sigpresent) %i\n',task.thistrial.sigpresent))
  
  % set the get response not to pass responses
  % since we reset this when we get the singal
  % from the sanghun task
  task.getResponse(:) = 0;
  
  % determine threshold coherence
  task.thistrial.threshold = task.stair{stairnum}.threshold;
  
  % pre interval, coherence is 0
  stimulus.dots.coherence = 0;
  % at the beginning to each trial, set color
  stimulus.dots.color = [1 1 1];
  stimulus.dots.colorClock = 0;
  
  myscreen.pedestaldir = nan; % line hints for discrim.

  
elseif (task.thistrial.thisseg == 3) 
  
  myscreen.sigpresent = task.thistrial.sigpresent;
  myscreen.pedestaldir = nan;
  
  % interval 1
  if strcmp(myscreen.exptype, 'fine_direction_thresholds')  
    % if the main task is fine d.o.m. discrim
    % randomly pick one of the two main directions as the pedestal
%    pedestaldir = task.directions( randsample(2,1) );
    pedestaldir = task.thistrial.pedestaldir;
    % write out a trace for old times sake.
    myscreen = writeTrace(pedestaldir,myscreen.stimtrace,myscreen,1);
    myscreen.pedestaldir = pedestaldir;
    if task.thistrial.sigpresent 
      stimulus.dots.dir = pedestaldir - task.thistrial.threshold/2; % more clockwise
    else
      stimulus.dots.dir = pedestaldir + task.thistrial.threshold/2;
    end
    disp(sprintf('stimulus.dots.dir: %0.1f %0.1f',stimulus.dots.dir,task.thistrial.threshold));
    % and fixed coherence
    stimulus.dots.coherence = task.thistrial.coherence;%myscreen.fixed_coherence;
    disp(sprintf('USING coherence of %0.2f',stimulus.dots.coherence));
  elseif  strcmp(myscreen.exptype, 'coherence_thresholds')  
    % sanghun task is the main task
    % select dir 1 if sigpresent, dir 2 if not.
    stimulus.dots.dir = task.directions(2-task.thistrial.sigpresent);
    stimulus.dots.coherence = task.stair{1}.threshold;
  elseif  strcmp(myscreen.exptype, 'fixed_params')  
    pedestaldir = task.thistrial.pedestaldir;
    myscreen.pedestaldir = pedestaldir;
    stimulus.dots.dir = task.directions(2-task.thistrial.sigpresent);
    % and fixed coherence
    stimulus.dots.coherence = task.thistrial.coherence;
  else
    error('(UHOH) mistyped exptype??')
  end
  
elseif ( task.thistrial.thisseg == 4 )
  % response interval
  stimulus.dots.coherence = 0;
  
elseif (task.thistrial.thisseg == 5 )
  
  % post interval, coherence is 0
  stimulus.dots.coherence = 0;
  % color could be white (or green/red depending on feedback
  % options)
  if myscreen.colorClockLen == 0
    % reset color to white
    stimulus.dots.color = [1 1 1]; 
    % otherwise let the dot drawing callback take care of color.
  end
  myscreen.pedestaldir = nan; % line hints for discrim.

end

% stimulus.dots.dir = task.thistrial.dir;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each FRAME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = stimDrawStimulusCallback(task, myscreen)

global stimulus

% if we want the subsidiary task to START
if(stimulus.startSubsidiary == 1)
 stimulus.startSubsidiary = 0;
 task = jumpSegment(task);
end

% if we want the subsidiary task to STOP
if(stimulus.endSubsidiary == 1)
 stimulus.endSubsidiary = 0;
 task = jumpSegment(task,inf);
end

mglClearScreen;

% if we are in the response interval, and the color clock is running, then count down
if stimulus.dots.colorClock >= 1
  stimulus.dots.colorClock = stimulus.dots.colorClock-1;
  if (stimulus.dots.colorClock <= 0)
    stimulus.dots.color = [1 1 1];
  end
end

if myscreen.isiDots || (task.thistrial.thisseg == 2)
  stimulus = updateDots(stimulus,myscreen);
end

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
mglGluDisk(stimulus.dots.xcenter,stimulus.dots.ycenter,[stimulus.dots.rmax stimulus.dots.rmax],[1 1 1],60);
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

mglStencilSelect(0);

if any(strcmp(myscreen.exptype, {'fine_direction_thresholds','fixed_params'}))
  % then we are doing a fine discrimination task
  % draw the marks for direction discrim
  % two markers - one each at dir1 and dir2
  r0 = 1.1*stimulus.dots.rmax;
  r1 = 1.2*r0;

  
  if isnan( myscreen.pedestaldir ) % hint lines are white
    % by default
    thetas = (myscreen.directions(:)');
    mark_x0 = [r0].*cosd(thetas);
    mark_y0 = [r0].*sind(thetas);
    mark_x1 = [+r1].*cosd(thetas);
    mark_y1 = [+r1].*sind(thetas);
    mglLines2(mark_x0, mark_y0, mark_x1, mark_y1, 8, [1 1 1]); % white
  else
    % label yellow where the task happens
    thetas = (myscreen.pedestaldir);
    mark_x0 = [r0].*cosd(thetas);
    mark_y0 = [r0].*sind(thetas);
    mark_x1 = [+r1].*cosd(thetas);
    mark_y1 = [+r1].*sind(thetas);
    mglLines2(mark_x0, mark_y0, mark_x1, mark_y1, 8, [1 1 0]); %
    
    % the other one
    thetas = (setdiff(myscreen.directions(:), myscreen.pedestaldir)');
    mark_x0 = [r0].*cosd(thetas);
    mark_y0 = [r0].*sind(thetas);
    mark_x1 = [+r1].*cosd(thetas);
    mark_y1 = [+r1].*sind(thetas);
    mglLines2(mark_x0, mark_y0, mark_x1, mark_y1, 8, [1 1 1]); %
  end
  
elseif strcmp(myscreen.exptype, 'coherence_thresholds')
  % then we are doing a coarse discrimination task
  % draw the marks for direction discrim
  % one markers - one at 90 (UP)
  r0 = 1.1*stimulus.dots.rmax;
  r1 = 1.2*r0;
  thetas = 90;
  mark_x0 = [r0].*cosd(thetas);
  mark_y0 = [r0].*sind(thetas);
  mark_x1 = [+r1].*cosd(thetas);
  mark_y1 = [+r1].*sind(thetas);
  mglLines2(mark_x0, mark_y0, mark_x1, mark_y1, 8, [1 1 0]);
   
end

mglStencilSelect(myscreen.stencil.dots);

% draw the dots
mglPoints2(stimulus.dots.x,stimulus.dots.y,stimulus.dots.dotsize, ...
	   stimulus.dots.color);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to get the response in the dot task...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = stimTrialResponseCallback(task,myscreen)

global stimulus

% update number of responses happens in task code.

% the response to the dot task is either the first or second,
% depending on which task is prioritized - check which one:
if ( task.thistrial.gotResponse == task.whichResponse-1 ) 
  % get the relevant staircase number
  stairnum = find(task.thistrial.coherence == task.parameter.coherence);
  disp(sprintf('Using staircase: %i',stairnum));
  if strcmp(myscreen.exptype, 'fine_direction_thresholds')%  &&  myscreen.pedestaldir < 90
    % then response contingencies are reversed:
    % button one means: in the more ccw direction 
    % button two means: in the more cw direction
    % -- this changes if myscreen.pedestaldir goes < 90
    correctButton = 1+task.thistrial.sigpresent;
  elseif strcmp(myscreen.exptype, 'coherence_thresholds')
    correctButton = 2-task.thistrial.sigpresent;
  else
    correctButton = 1+task.thistrial.sigpresent;
  end
    
  % update staircase
  if isfield(task,'stair')
    task.stair{stairnum} = upDownStaircase(task.stair{stairnum},task.thistrial.buttonState(correctButton));
    task.threshold = task.stair{stairnum}.threshold;
  end
  % keep reaction time
  task.thistrial.reactionTime = mglGetSecs-task.thistrial.segStartSeconds;

    
  if (task.thistrial.buttonState( correctButton ))
    disp(sprintf('correct threshold=%0.2f',task.threshold));
    if strcmp(myscreen.exptype,'fine_direction_thresholds')
      % change fixation cross clor
      global fixStimulus;
      fixStimulus.fixColor = [0 1 0];
      fixStimulus.fixColorClock = stimulus.dots.colorClockLen;
    else
      stimulus.dots.colorClock = stimulus.dots.colorClockLen;
      if stimulus.dots.colorClock > 0
	% it's possible that colorClockLen is set to 0 - so never turn green or red!
	stimulus.dots.color = [0 1 0];
      end  
    end
    %myscreen = writetrace(task.threshold+1,task.taskTracenum,myscreen);
    task.correct = task.correct+1;
  else
    disp(sprintf('incorrect threshold=%0.2f',task.threshold));
    if strcmp(myscreen.exptype,'fine_direction_thresholds')
      % change fixation cross clor
      global fixStimulus;
      fixStimulus.fixColor = [1 0 0];
      fixStimulus.fixColorClock = stimulus.dots.colorClockLen;
    else
      stimulus.dots.colorClock = stimulus.dots.colorClockLen;
      if stimulus.dots.colorClock > 0
	stimulus.dots.color = [1 0 0];
      end
    end      
    %myscreen = writetrace(-1,task.taskTracenum,myscreen);
    task.incorrect = task.incorrect+1;
  end  
end


