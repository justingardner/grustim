% stimprog.m
%
%      usage: [stimulus,task,myscreen]=stimprog(type,stimulus)
%         by: justin gardner
%       date: 01/13/06
%    purpose: stimulus program
%             type can be grating, dots or gratingloc
%       e.g.: stimprog('grating')
%
function [stimulus, task, myscreen] = stimprog(type,stimulus)

% check arguments
if ~any(nargin == [1 2])
  help stimprog
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other screen parameters
myscreen.autoclosescreen = 1;
myscreen.allowpause = 0;
myscreen.displayname = 'projector';

if ~isempty(strfind(type,'dots'))
  myscreen.background = 'black';
else
  myscreen.background = 'gray';
end  

myscreen = initScreen(myscreen);

% check to see if stimulus is defined
if exist('stimulus') ~= 1
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % set up the dot stimulus
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ~isempty(strfind(type,'dots'))
    stimulus.dots.xcenter = 0;
    stimulus.dots.ycenter = 0;
    stimulus.dots.width = myscreen.imageWidth;
    stimulus.dots.height = myscreen.imageHeight;
    stimulus.dots.density = 3;
    stimulus.dots.opticflowspeed = 5/myscreen.framesPerSecond;
    stimulus.dots.dotsize = 3;
    stimulus.dots.rmax = 10;
    stimulus.dots.lifespan = 5;
    %stimulus.dots.dotmotion = 'opticflow';
    stimulus.dots.dotmotion = 'linear2';
    %stimulus.dots.dottype = 'randomposition';
    %stimulus.dots.dottype = 'randomdirection';
    stimulus.dots.dottype = 'randomwalk';
    stimulus.dots.dotsame = 0;
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % set up the grating stimulus
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ~isempty(strfind(type,'grating'))
    stimulus.grating.sf = 1.5;
    stimulus.grating.updateFreq = [];
    stimulus.grating.phase = [2*pi*(1:16)/16];
    stimulus.grating.phase = [0 pi];
    stimulus.grating.width = myscreen.imageHeight;
    stimulus.grating.height = myscreen.imageHeight;
    stimulus.grating.orientation = 15:60:180;
%    stimulus.grating.orientation = 22.5:22.5:180;
%    stimulus.grating.orientation = [0 180];
%    stimulus.grating.orientation = 36:36:180;
    stimulus.grating.double = 1;
    stimulus.grating.maskangle = 15;
    stimulus.grating.fun = 'square';
    stimulus.grating.checker = 0;
    %stimulus.grating.gaborWidth = 2;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up baseline task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
triallen = 30;
if ~isempty(strfind(type,'localizer'))
  if ~isempty(strfind(type,'dots'))
    % set up incoherent dots
    task{1}.parameter.direction = [0;0];
    task{1}.parameter.coherence = [0 1;1 0];
    task{1}.parameter.speed = [2;2];
    task{1}.parameter.randomizeDirection = 1;
    task{1}.randomize = 0;
    % set up the segments
    task{1}.seglen = 0.5*ones(1,round(6*2));
    % set to run this task for one block of one trial and end
    task{1}.numblocks = inf;
    for i = 1:length(task{1}.seglen)
      % set up the traces to be displayed
      task{1}.writetrace{i}.tracenum = 1;
      task{1}.writetrace{i}.tracevar{1} = 'coherence';
    end
  else
    % set up phase and orientation randomized stimulus
    task{1}.parameter.contrast = [[0 1];[1 0]];
    task{1}.parameter.randomizeOrientation = 1;
    % set up the segments
    task{1}.seglen = repmat(.250,1,48);
    for i = 1:length(task{1}.seglen)
      % set up the traces to be displayed
      task{1}.writetrace{i}.tracenum = 1;
      task{1}.writetrace{i}.tracevar{1} = 'contrast';
    end
  end
  % 0 is no randomization, 1 is complete and 2 is sequential
  task{1}.random = 0;
  % wait for backtick before starting
  task{1}.waitForBacktick = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this part is for 8 directions/orientations task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
  if ~isempty(strfind(type,'dots'))
    % set up incoherent dots
    task{1}.parameter.direction = [0;0];
    task{1}.parameter.coherence = [0;0];
    task{1}.parameter.speed = [2;2];
    task{1}.parameter.randomizeDirection = 0;
    % set up the segments
    task{1}.seglen = 12;
    % set to run this task for one block of one trial and end
    task{1}.numblocks = 1;
  else
    % set up phase and orientation randomized stimulus
    task{1}.parameter.orientation = [stimulus.grating.orientation ; stimulus.grating.orientation];
    task{1}.parameter.contrast = [1;1];
    task{1}.parameter.randomizeOrientation = 0;
    % set up the segments
    task{1}.seglen = repmat(0.250,1,2);
    % set to run for 6 seconds
    task{1}.numtrials = 12;%triallen*2;
  end
  % 0 is no randomization, 1 is complete and 2 is sequential
  task{1}.random = 1;
  % this parameter is just for display on the output
  task{1}.parameter.traceoutput = -1;
  for i = 1:length(task{1}.seglen)
    % set up the traces to be displayed
    task{1}.writetrace{i}.tracenum = [1 2];
    task{1}.writetrace{i}.tracevar{1} = 'traceoutput';
    task{1}.writetrace{i}.tracevar{2} = 'traceoutput';
  end
  % wait for backtick before starting
  task{1}.waitForBacktick = 1;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % set up task
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % set up parameters over which we will randomize
  if ~isempty(strfind(type,'dots'))
%    dirs8 = [30 60 120 150 210 240 300 330];
    dir3= [15 135 255];
    % set up dots stimului parameters
    task{2}.parameter.direction = [ dir3 ; dir3];
    task{2}.parameter.coherence = [1;1];
    task{2}.parameter.speed = [2;2];
    task{2}.parameter.randomizeDirection = 0;
    % set up segments of trials
    task{2}.seglen = [0.250 (triallen-0.250)];
    % set up which traces to write out
    task{2}.writetrace{1}.tracenum = [1 2];
    task{2}.writetrace{1}.tracevar{1} = 'direction';
    task{2}.writetrace{1}.tracevar{2} = 'direction';
    task{2}.writetrace{1}.tracerow = [1 2];
    task{2}.writetrace{1}.usenum = [1 1];
  else
    % set up grating stimului parameters
    task{2}.parameter.orientation = [stimulus.grating.orientation ; stimulus.grating.orientation];
    task{2}.parameter.contrast = [1;1];
    task{2}.parameter.randomizeOrientation = 0;
    % set up segments of trials
    task{2}.seglen = repmat(0.250,1,triallen/0.250);
    % set up which traces to write out
    task{2}.writetrace{1}.tracenum = [1 2];
    task{2}.writetrace{1}.tracevar{1} = 'orientation';
    task{2}.writetrace{1}.tracevar{2} = 'orientation';
    task{2}.writetrace{1}.tracerow = [1 2];
    task{2}.writetrace{1}.usenum = [1 1];
  end
  % 0 is no randomization, 1 is complete and 2 is sequential
  task{2}.random = 0;
  % set up whether to keep time in frame ticks or ms.
  task{2}.timeInTicks = 0;
  % set up fixation color on each segment of task
  task{2}.fixcolor = {};
  % set up whether to get response on each segment of trial
  task{2}.getresponse = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialze tasks and stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(stimulus,'init') || (stimulus.init == 0)
  [stimulus myscreen] = initstimulus(stimulus,myscreen);
end
% make into a cell array of tasks
if isstruct(task)
  temp{1} = task;
  task = temp;
  clear temp;
end

% initialze tasks
for tasknum = 1:length(task)
  task{tasknum} = inittask(task{tasknum},myscreen);
end

% set which task is active
tnum = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyecalibdisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the tasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = myscreentic(myscreen);
while (tnum <= length(task))
  % set the screen background color
  if (myscreen.background ~= 0)
    Screen('FillRect', myscreen.w, myscreen.background);
  end
  % update the task
  [stimulus task myscreen tnum] = updatetask(stimulus,task,myscreen,tnum);
  % display fixation cross
  myscreen = fixdispstair(myscreen);
  % flip screen
  myscreen = tickscreen(myscreen);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finish up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (myscreen.autoclosescreen)
  Screen('closeall');
else
  Screen('flip',myscreen.w);
end
myscreen = endscreen(myscreen);
% save stim data
savestimdata(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of program, usually gets here from error thrown
% from tickscreen, which captures the user hitting ESC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
catch
  if (myscreen.autoclosescreen),Screen('closeall');,end
  % see if the last error was just the error
  % thrown by ending the task
  err = lasterror;
  % if not rethrow the error
  if (isempty(strfind(err.message,'taskend')))
    rethrow(err);
  else
    mydisp(sprintf('End task\n'));
    % otherwise we are done
    % package up the output variables
    myscreen.task = task;
    return
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     functions to run trials/task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     STIMULUS DOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init the dots stimuli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initdots_linear(myscreen,task)

% get the passed in parameters
dots.width = task.width;
dots.height = task.height;
dots.xcenter = task.xcenter+myscreen.imageWidth/2-dots.width/2;
dots.ycenter = task.ycenter+myscreen.imageHeight/2-dots.height/2;
dots.density = task.density;
dots.dotsize = task.dotsize;
dots.coherence = -1;
dots.type = task.dottype;
dots.same = task.dotsame;
dots.dotmotion = task.dotmotion;

% get the number of dots
dots.n = round(dots.width*dots.height*dots.density);

% get max and min points for dots
dots.xmin = 0;
dots.xmax = task.width;
dots.ymin = 0;
dots.ymax = task.height;

% set direction of dots
dots.dir = 0;
dots.speed = 0;

% get initial position
dots.x = rand(1,dots.n)*dots.width;
dots.y = rand(1,dots.n)*dots.height;

% get the step size
dots.stepsize = dots.speed/myscreen.framesPerSecond;
dots.xstep = cos(dots.dir)*dots.stepsize;
dots.ystep = sin(dots.dir)*dots.stepsize;

% just for random direction, assign random directions to everybody
dots.rtheta = rand(1,dots.n)*2*pi;
dots.rxstep = cos(dots.rtheta)*dots.stepsize;
dots.rystep = sin(dots.rtheta)*dots.stepsize;

% if dots have a lifteime
dots.lifespan = task.lifespan;
dots.age = rand(1,dots.n)*dots.lifespan;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init the dots stimuli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initdots_linear2(myscreen,task)

% get the passed in parameters
dots.width = task.width;
dots.height = task.height;
dots.xcenter = task.xcenter+myscreen.imageWidth/2-dots.width/2;
dots.ycenter = task.ycenter+myscreen.imageHeight/2-dots.height/2;
dots.density = task.density;
dots.dotsize = task.dotsize;
dots.coherence = [-1 ; -1];
dots.type = task.dottype;
dots.same = task.dotsame;
dots.dotmotion = task.dotmotion;

% get the number of dots
dots.n = round(dots.width*dots.height*dots.density);

% get max and min points for dots
dots.xmin = 0;
dots.xmax = task.width;
dots.xhalf = (dots.xmax-dots.xmin)/2;
dots.ymin = 0;
dots.ymax = task.height;
dots.yhalf = (dots.ymax-dots.ymin)/2;

% set direction of dots
dots.dir = [0 0];
dots.speed = 0;

% get initial position
dots.x = rand(1,dots.n)*dots.width;
dots.y = rand(1,dots.n)*dots.height;

% get the step size
dots.stepsize = dots.speed/myscreen.framesPerSecond;
dots.xstep = cos(dots.dir)*dots.stepsize;
dots.ystep = sin(dots.dir)*dots.stepsize;

% just for random direction, assign random directions to everybody
dots.rtheta = rand(1,dots.n)*2*pi;
dots.rxstep = cos(dots.rtheta)*dots.stepsize;
dots.rystep = sin(dots.rtheta)*dots.stepsize;

% if dots have a lifteime
dots.lifespan = task.lifespan;
dots.age = rand(1,dots.n)*dots.lifespan;

% if we are using dots that appear on one half or the other screen
dots.lefthalf = dots.x < dots.xhalf;
dots.righthalf = dots.x >= dots.xhalf;

% mask for upper portion of screen
pointlist1(1,:) = [0 1];
pointlist1(2,:) = [-myscreen.imageWidth/2 myscreen.imageWidth/2];
pointlist1(3,:) = [myscreen.imageWidth/2 myscreen.imageWidth/2];
pointlist1(4,:) = [0 1];

% mask for lower portion of screen
pointlist2(1,:) = [0 -1];
pointlist2(2,:) = [-myscreen.imageWidth/2 -myscreen.imageWidth/2];
pointlist2(3,:) = [myscreen.imageWidth/2 -myscreen.imageWidth/2];
pointlist2(4,:) = [0 -1];

% convert masks from degrees to points
for i = 1:size(pointlist1,1)
  dots.mask.pointlist1(i,1) = pointlist1(i,1)*myscreen.xdeg2pix+myscreen.centerx;
  dots.mask.pointlist1(i,2) = pointlist1(i,2)*myscreen.ydeg2pix+myscreen.centery;
end
for i = 1:size(pointlist2,1)
  dots.mask.pointlist2(i,1) = pointlist2(i,1)*myscreen.xdeg2pix+myscreen.centerx;
  dots.mask.pointlist2(i,2) = pointlist2(i,2)*myscreen.ydeg2pix+myscreen.centery;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init the dots stimuli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initdots_opticflow(myscreen,task)

% set type
dots.width = task.width;
dots.height = task.height;
dots.xcenter = task.xcenter+myscreen.imageWidth/2-dots.width/2;
dots.ycenter = task.ycenter+myscreen.imageHeight/2-dots.height/2;
dots.density = task.density;
dots.dotsize = task.dotsize;
dots.coherence = -1;
dots.type = task.dottype;
dots.same = task.dotsame;
dots.dotmotion = task.dotmotion;

% focal length to projection plane
% projection plane is defined to be 
% 1 unit wide and high, so with 
% this focal length, we are looking at
% a view of the world with a 90 deg fov
dots.f = .5;

% translation and rotation matrices
dots.T = [0 0 task.opticflowspeed];
dots.R = [0 0 0];

% maximum depth of points
dots.maxZ = 10-dots.f;dots.minZ = dots.f;
dots.maxX = 10;
dots.maxY = 10*myscreen.imageHeight/myscreen.imageWidth;

% make a brick of points
dots.n = round(myscreen.imageWidth*myscreen.imageHeight*dots.density);

% initial position of dots
dots.X = 2*dots.maxX*rand(1,dots.n)-dots.maxX;
dots.Y = 2*dots.maxY*rand(1,dots.n)-dots.maxY;
dots.Z = (dots.maxZ-dots.minZ)*rand(1,dots.n)+dots.minZ;

% get projection on to plane
dots.xproj = dots.f*dots.X./dots.Z;
dots.yproj = dots.f*dots.Y./dots.Z;

% generate a random transformation matrix for each incoherent point
dots.randT = rand(3,dots.n)-0.5;
% and normalize the transformation to have the same length
% (i.e. speed) as the real transformation matrix
dots.randT = sqrt(sum(dots.T.^2))*dots.randT./([1 1 1]'*sqrt(sum(dots.randT.^2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init the dots stimuli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initdots_radial(myscreen,task)

% get the passed in parameters
dots.rmin = 0;
dots.rmax = task.rmax;
dots.xcenter = task.xcenter+myscreen.imageWidth/2;
dots.ycenter = task.ycenter+myscreen.imageHeight/2;
dots.density = task.density;
dots.dotsize = task.dotsize;
dots.coherence = -1;
dots.type = task.dottype;
dots.same = task.dotsame;
dots.dotmotion = task.dotmotion;

% get the number of dots
dots.n = round(pi*(dots.rmax^2)*dots.density);

% get initial position
dots.r = rand(1,dots.n)*dots.rmax;
dots.theta = rand(1,dots.n)*2*pi;

% get the step size
dots.stepsize = task.speed/myscreen.framesPerSecond;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots for opticflow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updatedots_randomdirection_opticflow(dots,coherence,samerule,myscreen)

% get the coherent and incoherent dots
if (~samerule || (dots.coherence ~= coherence))
  dots.incoherent = rand(1,dots.n) > coherence;
  dots.incoherentn = sum(dots.incoherent);
  dots.coherent = ~dots.incoherent;
  dots.coherence = coherence;
end

% update relative position of dots in 3-space to observer
dots.X(dots.coherent) = dots.X(dots.coherent)-dots.T(1);
dots.Y(dots.coherent) = dots.Y(dots.coherent)-dots.T(2);
dots.Z(dots.coherent) = dots.Z(dots.coherent)-dots.T(3);

% now move the incoherent points according to the random trasnformation
dots.X(dots.incoherent) = dots.X(dots.incoherent)-dots.randT(1,dots.incoherent);
dots.Y(dots.incoherent) = dots.Y(dots.incoherent)-dots.randT(2,dots.incoherent);
dots.Z(dots.incoherent) = dots.Z(dots.incoherent)-dots.randT(3,dots.incoherent);

% get all points that have fallen off the screen
offscreen = dots.Z<dots.minZ;

% and put them at the furthest distance
dots.Z(offscreen) = dots.maxZ;

% get all points that have fallen out of view
offscreen = dots.Z>dots.maxZ;
% and move them to the front plane
dots.Z(offscreen) = dots.minZ;

% put points fallen off the X edge back
offscreen = dots.X < -dots.maxX;
dots.X(offscreen) = dots.X(offscreen)+2*dots.maxX;
offscreen = dots.X > dots.maxX;
dots.X(offscreen) = dots.X(offscreen)-2*dots.maxX;

% put points fallen off the Y edge back
offscreen = dots.Y < -dots.maxY;
dots.Y(offscreen) = dots.Y(offscreen)+2*dots.maxY;
offscreen = dots.Y > dots.maxY;
dots.Y(offscreen) = dots.Y(offscreen)-2*dots.maxY;

% project on to screen
dots.xproj = dots.f*dots.X./dots.Z;
dots.yproj = dots.f*dots.Y./dots.Z;

% get actual screen coordinates
dots.x = myscreen.xdeg2pix*(dots.xproj*myscreen.imageWidth+myscreen.imageWidth/2);
dots.y = myscreen.ydeg2pix*(dots.yproj*myscreen.imageHeight+myscreen.imageHeight/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots for opticflow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updatedots_randomwalk_opticflow(dots,coherence,samerule,myscreen)

% get the coherent and incoherent dots
if (~samerule || (dots.coherence ~= coherence))
  dots.incoherent = rand(1,dots.n) > coherence;
  dots.incoherentn = sum(dots.incoherent);
  dots.coherent = ~dots.incoherent;
  dots.coherence = coherence;
end

% generate a random transformation matrix for each incoherent point
dots.randT = rand(3,dots.incoherentn)-0.5;
% and normalize the transformation to have the same length
% (i.e. speed) as the real transformation matrix
dots.randT = sqrt(sum(dots.T.^2))*dots.randT./([1 1 1]'*sqrt(sum(dots.randT.^2)));


% update relative position of dots in 3-space to observer
dots.X(dots.coherent) = dots.X(dots.coherent)-dots.T(1);
dots.Y(dots.coherent) = dots.Y(dots.coherent)-dots.T(2);
dots.Z(dots.coherent) = dots.Z(dots.coherent)-dots.T(3);

% now move the incoherent points according to the random trasnformation
dots.X(dots.incoherent) = dots.X(dots.incoherent)-dots.randT(1,:);
dots.Y(dots.incoherent) = dots.Y(dots.incoherent)-dots.randT(2,:);
dots.Z(dots.incoherent) = dots.Z(dots.incoherent)-dots.randT(3,:);

% get all points that have fallen off the screen
offscreen = dots.Z<dots.minZ;

% and put them at the furthest distance
dots.Z(offscreen) = dots.maxZ;

% get all points that have fallen out of view
offscreen = dots.Z>dots.maxZ;
% and move them to the front plane
dots.Z(offscreen) = dots.minZ;

% put points fallen off the X edge back
offscreen = dots.X < -dots.maxX;
dots.X(offscreen) = dots.X(offscreen)+2*dots.maxX;
offscreen = dots.X > dots.maxX;
dots.X(offscreen) = dots.X(offscreen)-2*dots.maxX;

% put points fallen off the Y edge back
offscreen = dots.Y < -dots.maxY;
dots.Y(offscreen) = dots.Y(offscreen)+2*dots.maxY;
offscreen = dots.Y > dots.maxY;
dots.Y(offscreen) = dots.Y(offscreen)-2*dots.maxY;

% project on to screen
dots.xproj = dots.f*dots.X./dots.Z;
dots.yproj = dots.f*dots.Y./dots.Z;

% get actual screen coordinates
dots.x = myscreen.xdeg2pix*(dots.xproj*myscreen.imageWidth+myscreen.imageWidth/2);
dots.y = myscreen.ydeg2pix*(dots.yproj*myscreen.imageHeight+myscreen.imageHeight/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots for opticflow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updatedots_randomposition_opticflow(dots,coherence,samerule,myscreen)

% get the coherent and incoherent dots
if (~samerule || (dots.coherence ~= coherence))
  dots.incoherent = rand(1,dots.n) > coherence;
  dots.incoherentn = sum(dots.incoherent);
  dots.coherent = ~dots.incoherent;
  dots.coherence = coherence;
end

% update relative position of dots in 3-space to observer
dots.X(dots.coherent) = dots.X(dots.coherent)-dots.T(1);
dots.Y(dots.coherent) = dots.Y(dots.coherent)-dots.T(2);
dots.Z(dots.coherent) = dots.Z(dots.coherent)-dots.T(3);

% replot the other ones
dots.X(dots.incoherent) = 2*dots.maxX*rand(1,dots.incoherentn)-dots.maxX;
dots.Y(dots.incoherent) = 2*dots.maxY*rand(1,dots.incoherentn)-dots.maxY;
dots.Z(dots.incoherent) = (dots.maxZ-dots.minZ)*rand(1,dots.incoherentn)+dots.minZ;
% get all points that have fallen off the screen
offscreen = dots.Z<dots.minZ;

% and put them at the furthest distance
dots.Z(offscreen) = dots.maxZ;

% get all points that have fallen out of view
offscreen = dots.Z>dots.maxZ;
% and move them to the front plane
dots.Z(offscreen) = dots.minZ;

% put points fallen off the X edge back
offscreen = dots.X < -dots.maxX;
dots.X(offscreen) = dots.X(offscreen)+2*dots.maxX;
offscreen = dots.X > dots.maxX;
dots.X(offscreen) = dots.X(offscreen)-2*dots.maxX;

% put points fallen off the Y edge back
offscreen = dots.Y < -dots.maxY;
dots.Y(offscreen) = dots.Y(offscreen)+2*dots.maxY;
offscreen = dots.Y > dots.maxY;
dots.Y(offscreen) = dots.Y(offscreen)-2*dots.maxY;

% project on to screen
dots.xproj = dots.f*dots.X./dots.Z;
dots.yproj = dots.f*dots.Y./dots.Z;

% get actual screen coordinates
dots.x = myscreen.xdeg2pix*(dots.xproj*myscreen.imageWidth+myscreen.imageWidth/2);
dots.y = myscreen.ydeg2pix*(dots.yproj*myscreen.imageHeight+myscreen.imageHeight/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots randomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updatedots_randomposition_linear(dots,coherence,samerule,myscreen)

if (~samerule || (dots.coherence ~= coherence))
  % pick a random set of dots
  dots.coherent = rand(1,dots.n) < coherence;
  dots.coherence = coherence;
end

% now move those dots in the right direction
dots.x(dots.coherent) = dots.x(dots.coherent)+dots.xstep;
dots.y(dots.coherent) = dots.y(dots.coherent)+dots.ystep;

% other dots, get randomly replotted
dots.x(~dots.coherent) = rand(1,sum(~dots.coherent))*dots.width;
dots.y(~dots.coherent) = rand(1,sum(~dots.coherent))*dots.height;

% make sure we haven't gone off the patch
dots.x(dots.x > dots.xmax) = dots.x(dots.x > dots.xmax)-dots.width;
dots.x(dots.x < dots.xmin) = dots.x(dots.x < dots.xmin)+dots.width;
dots.y(dots.y > dots.ymax) = dots.y(dots.y > dots.ymax)-dots.height;
dots.y(dots.y < dots.ymin) = dots.y(dots.y < dots.ymin)+dots.height;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots randomposition, radial motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updatedots_randomposition_radial(dots,coherence,samerule,myscreen)

if (~samerule || (dots.coherence ~= coherence))
  % pick a random set of dots
  dots.coherent = rand(1,dots.n) < coherence;
  dots.coherence = coherence;
end

% now move those dots in the right direction
dots.r(dots.coherent) = dots.r(dots.coherent)+dots.stepsize;

% other dots, get randomly replotted
dots.r(~dots.coherent) = rand(1,sum(~dots.coherent))*dots.rmax;
dots.theta(~dots.coherent) = rand(1,sum(~dots.coherent))*2*pi;

% make sure we haven't gone off the patch
dots.r(dots.r > dots.rmax) = dots.r(dots.r > dots.rmax)-dots.rmax;
dots.r(dots.r < dots.rmin) = dots.r(dots.r < dots.rmin)+dots.rmax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots randomwalk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updatedots_randomwalk_linear2(dots,coherence,samerule,myscreen)

if (~samerule || ~isequal(dots.coherence,coherence))
  % pick a random set of dots
  dots.coherent(dots.lefthalf) = rand(1,sum(dots.lefthalf)) < coherence(1);
  dots.coherent(dots.righthalf) = rand(1,sum(dots.righthalf)) < coherence(2);
  dots.coherence = coherence;
end

% now move those dots in the right direction
dots.x(dots.coherent&dots.lefthalf) = dots.x(dots.coherent&dots.lefthalf)+dots.xstep(1);
dots.y(dots.coherent&dots.lefthalf) = dots.y(dots.coherent&dots.lefthalf)+dots.ystep(1);
dots.x(dots.coherent&dots.righthalf) = dots.x(dots.coherent&dots.righthalf)+dots.xstep(2);
dots.y(dots.coherent&dots.righthalf) = dots.y(dots.coherent&dots.righthalf)+dots.ystep(2);

% other dots, get moved in a random directon
dots.x(~dots.coherent) = dots.x(~dots.coherent)+cos(rand(1,sum(~dots.coherent))*2*pi)*dots.stepsize;
dots.y(~dots.coherent) = dots.y(~dots.coherent)+sin(rand(1,sum(~dots.coherent))*2*pi)*dots.stepsize;

% make sure we haven't gone off the patch
% do the dots separately for left and right hand side
dots.x(dots.lefthalf & (dots.x < dots.xmin)) = dots.x(dots.lefthalf & (dots.x < dots.xmin))+dots.width;
dots.x(dots.lefthalf & (dots.x > dots.xhalf)) = dots.x(dots.lefthalf & (dots.x > dots.xhalf))-dots.width/2;
dots.x(dots.righthalf & (dots.x < dots.xhalf)) = dots.x(dots.righthalf & (dots.x < dots.xhalf))+dots.width/2;
dots.x(dots.righthalf & (dots.x > dots.xmax)) = dots.x(dots.righthalf & (dots.x > dots.xmax))-dots.width/2;
dots.y(dots.y > dots.ymax) = dots.y(dots.y > dots.ymax)-dots.height;
dots.y(dots.y < dots.ymin) = dots.y(dots.y < dots.ymin)+dots.height;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots randomwalk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updatedots_randomwalk_linear(dots,coherence,samerule,myscreen)

if (~samerule || (dots.coherence ~= coherence))
  % pick a random set of dots
  dots.coherent = rand(1,dots.n) < coherence;
  dots.coherence = coherence;
end

% now move those dots in the right direction
dots.x(dots.coherent) = dots.x(dots.coherent)+dots.xstep;
dots.y(dots.coherent) = dots.y(dots.coherent)+dots.ystep;

% other dots, get moved in a random directon
dots.x(~dots.coherent) = dots.x(~dots.coherent)+cos(rand(1,sum(~dots.coherent))*2*pi)*dots.stepsize;
dots.y(~dots.coherent) = dots.y(~dots.coherent)+sin(rand(1,sum(~dots.coherent))*2*pi)*dots.stepsize;

% make sure we haven't gone off the patch
dots.x(dots.x > dots.xmax) = dots.x(dots.x > dots.xmax)-dots.width;
dots.x(dots.x < dots.xmin) = dots.x(dots.x < dots.xmin)+dots.width;
dots.y(dots.y > dots.ymax) = dots.y(dots.y > dots.ymax)-dots.height;
dots.y(dots.y < dots.ymin) = dots.y(dots.y < dots.ymin)+dots.height;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots randomdirection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updatedots_randomdirection_linear(dots,coherence,samerule,myscreen)

if (~samerule || (dots.coherence ~= coherence))
  % pick a random set of dots
  dots.coherent = rand(1,dots.n) < coherence;
  dots.coherence = coherence;
end

% now move those dots in the right direction
dots.x(dots.coherent) = dots.x(dots.coherent)+dots.xstep;
dots.y(dots.coherent) = dots.y(dots.coherent)+dots.ystep;

% other dots, get moved in a random directon
dots.x(~dots.coherent) = dots.x(~dots.coherent)+dots.rxstep(~dots.coherent);
dots.y(~dots.coherent) = dots.y(~dots.coherent)+dots.rystep(~dots.coherent);

% make sure we haven't gone off the patch
dots.x(dots.x > dots.xmax) = dots.x(dots.x > dots.xmax)-dots.width;
dots.x(dots.x < dots.xmin) = dots.x(dots.x < dots.xmin)+dots.width;
dots.y(dots.y > dots.ymax) = dots.y(dots.y > dots.ymax)-dots.height;
dots.y(dots.y < dots.ymin) = dots.y(dots.y < dots.ymin)+dots.height;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots randomdirection with lifetime
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updatedots_randomdirection_lifetime_linear(dots,coherence,samerule,myscreen)

if (~samerule || (dots.coherence ~= coherence))
  % pick a random set of dots
  dots.coherent = rand(1,dots.n) < coherence;
  dots.coherence = coherence;
end

% now move those dots in the right direction
dots.x(dots.coherent) = dots.x(dots.coherent)+dots.xstep;
dots.y(dots.coherent) = dots.y(dots.coherent)+dots.ystep;

% other dots, get moved in a random directon
dots.x(~dots.coherent) = dots.x(~dots.coherent)+dots.rxstep(~dots.coherent);
dots.y(~dots.coherent) = dots.y(~dots.coherent)+dots.rystep(~dots.coherent);

% dots that have expired, get a new position
dots.age = dots.age+1;
deaddots = dots.age > dots.lifespan;
dots.x(deaddots) = rand(1,sum(deaddots))*dots.width;
dots.y(deaddots) = rand(1,sum(deaddots))*dots.height;
dots.age(deaddots) = 0;

% make sure we haven't gone off the patch
dots.x(dots.x > dots.xmax) = dots.x(dots.x > dots.xmax)-dots.width;
dots.x(dots.x < dots.xmin) = dots.x(dots.x < dots.xmin)+dots.width;
dots.y(dots.y > dots.ymax) = dots.y(dots.y > dots.ymax)-dots.height;
dots.y(dots.y < dots.ymin) = dots.y(dots.y < dots.ymin)+dots.height;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display the dots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = dispdots_linear(myscreen,dots)

% plot the points
Screen('glPoint',myscreen.w,255,(dots.xcenter+dots.x)*myscreen.xdeg2pix,(dots.ycenter+dots.y)*myscreen.ydeg2pix,dots.dotsize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display the dots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = dispdots_linear2(myscreen,dots)

% plot the points
Screen('glPoint',myscreen.w,255,(dots.xcenter+dots.x)*myscreen.xdeg2pix,(dots.ycenter+dots.y)*myscreen.ydeg2pix,dots.dotsize);

% draw the mask
Screen('FillPoly',myscreen.w,[0 0 0],dots.mask.pointlist1);
Screen('FillPoly',myscreen.w,[0 0 0],dots.mask.pointlist2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display the dots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = dispdots_opticflow(myscreen,dots)

% plot the points
Screen('glPoint',myscreen.w,255,dots.x,dots.y,dots.dotsize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display the dots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = dispdots_radial(myscreen,dots)

% get cartesian position
dots.x = dots.r.*cos(dots.theta);
dots.y = dots.r.*sin(dots.theta);

% plot the points
Screen('glPoint',myscreen.w,255,(dots.xcenter+dots.x)*myscreen.xdeg2pix,(dots.ycenter+dots.y)*myscreen.ydeg2pix,dots.dotsize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the dot dir and speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = setdirspeed(dots,dir,speed,myscreen)

% set direction of dots
dots.dir = dir;
dots.speed = speed;

% get the step size
for i = 1:length(dots.speed)
  dots.stepsize = dots.speed(i)/myscreen.framesPerSecond;
  dots.xstep(i) = cos(dots.dir(i))*dots.stepsize;
  dots.ystep(i) = sin(dots.dir(i))*dots.stepsize;
end

% just for random direction, assign random directions to everybody
dots.rtheta = rand(1,dots.n)*2*pi;
dots.rxstep = cos(dots.rtheta)*dots.stepsize;
dots.rystep = sin(dots.rtheta)*dots.stepsize;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     STIMULUS DOTS END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     STIMULUS Gratings start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init gratings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initGrating(stimulus,myscreen)

% get the dimensions
stimulus.grating.pixwidth = round(stimulus.grating.width*myscreen.xdeg2pix);
stimulus.grating.pixheight = round(stimulus.grating.height*myscreen.ydeg2pix);

% make sure pixwidth is odd
if iseven(stimulus.grating.pixwidth)
  stimulus.grating.pixwidth = stimulus.grating.pixwidth-1;
end

% get grid over which we should calculate
xhalfwidth = round(stimulus.grating.width/2);
yhalfwidth = round(stimulus.grating.height/2);
% get x and y
x = -xhalfwidth:2*xhalfwidth/(stimulus.grating.pixwidth-1):xhalfwidth;
y = -yhalfwidth:2*yhalfwidth/(stimulus.grating.pixheight-1):yhalfwidth;
% find xcenter
xleft = x(x<0);stimulus.grating.pixwidthleft = length(xleft);
xright = x(x>0);stimulus.grating.pixwidthright = length(xright);
% turn into grid
[xMesh,yMesh] = meshgrid(x,y);
[xLeftMesh,yLeftMesh] = meshgrid(xleft,y);
[xRightMesh,yRightMesh] = meshgrid(xright,y);

% now compute the mask
gratingMask = ((xMesh.^2)/xhalfwidth^2 + (yMesh.^2)/yhalfwidth^2) < 1;
gratingMaskLeft = ((xLeftMesh.^2)/xhalfwidth^2 + (yLeftMesh.^2)/yhalfwidth^2) < 1;
gratingMaskRight = ((xRightMesh.^2)/xhalfwidth^2 + (yRightMesh.^2)/yhalfwidth^2) < 1;
% create gaussian window
if isfield(stimulus.grating,'gaborWidth') && isfield(stimulus.grating,'gaborHeight')
  gratingMask = exp(-((xMesh.^2)/stimulus.grating.gaborWidth^2 + (yMesh.^2)/stimulus.grating.gaborHeight^2));
  gratingMaskLeft = exp(-((xLeftMesh.^2)/stimulus.grating.gaborWidth^2 + (yLeftMesh.^2)/stimulus.grating.gaborHeight^2));
  gratingMaskRight = exp(-((xMeshRight.^2)/stimulus.grating.gaborWidth^2 + (yMeshRight.^2)/stimulus.grating.gaborHeight^2));
end
% cut out appropriate angle
gratingMaskLeft(abs(r2d(atan(yLeftMesh./xLeftMesh)))>(90-stimulus.grating.maskangle)) = 0;
gratingMaskRight(abs(r2d(atan(yRightMesh./xRightMesh)))>(90-stimulus.grating.maskangle)) = 0;

% some defaults
if ~isfield(stimulus.grating,'checker')
  stimulus.grating.checker = 0;
end

% compute each asked for grating
disppercent(-inf,'Generating gratings');
for orientNum = 1:length(stimulus.grating.orientation)
  disppercent(orientNum/length(stimulus.grating.orientation));
  for phaseNum = 1:length(stimulus.grating.phase)
    % get this phase
    phase = stimulus.grating.phase(phaseNum);
    % get orientation
    angle = d2r(stimulus.grating.orientation(orientNum));
    % get spatial frequency
    f=stimulus.grating.sf*2*pi; 
    a=cos(angle)*f;
    b=sin(angle)*f;
    % also get orthogonal grating
    aorth=cos(angle+pi/2)*f;
    borth=sin(angle+pi/2)*f;
    % compute grating
    if stimulus.grating.double
      % create gratings
      mleft = eval(sprintf('%s(a*xLeftMesh+b*yLeftMesh+phase).*gratingMaskLeft',stimulus.grating.fun));
      mright = eval(sprintf('%s(a*xRightMesh+b*yRightMesh+phase).*gratingMaskRight',stimulus.grating.fun));
      % if checkerboard grating get the orthogonal grating
      % and multiple by it
      if (stimulus.grating.checker)
	% create gratings
	mleftOrth = eval(sprintf('%s(aorth*xLeftMesh+borth*yLeftMesh+phase).*gratingMaskLeft',stimulus.grating.fun));
	mrightOrth = eval(sprintf('%s(aorth*xRightMesh+borth*yRightMesh+phase).*gratingMaskRight',stimulus.grating.fun));
	mleft = mleft.*mleftOrth;
	mright = mright.*mrightOrth;
      end
      % now convert them to textures
      stimulus.grating.ltex(orientNum,phaseNum) = Screen('MakeTexture', myscreen.screenNumber, gammaCorrect(myscreen.grayIndex+myscreen.inc*mleft,myscreen));
      stimulus.grating.rtex(orientNum,phaseNum) = Screen('MakeTexture', myscreen.screenNumber, gammaCorrect(myscreen.grayIndex+myscreen.inc*mright,myscreen));
    else
      m = eval(sprintf('%s(a*xMesh+b*yMesh+phase).*gratingMask',stimulus.grating.fun));
      % if checkerboard grating get the orthogonal grating
      % and multiple by it
      if (stimulus.grating.checker)
	morth = eval(sprintf('%s(aorth*xMesh+borth*yMesh+phase).*gratingMask',stimulus.grating.fun));
	m = m.*morth;
      end
      % convert it to a texture
      stimulus.grating.tex(orientNum,phaseNum) = Screen('MakeTexture', myscreen.screenNumber, gammaCorrect(myscreen.grayIndex+myscreen.inc*m,myscreen));
    end
  end
end
disppercent(inf);

% get source and destination locations
if stimulus.grating.double
  % grating on left
  stimulus.grating.rdestrect = SetRect(myscreen.centerx,myscreen.centery-stimulus.grating.pixheight/2,myscreen.centerx+stimulus.grating.pixwidthleft-1,myscreen.centery+stimulus.grating.pixheight/2-1);
  % grating on right
  stimulus.grating.ldestrect = SetRect(myscreen.centerx-stimulus.grating.pixwidthright+1,myscreen.centery-stimulus.grating.pixheight/2+1,myscreen.centerx,myscreen.centery+stimulus.grating.pixheight/2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to update grating stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = updateGrating(stimulus,myscreen)

if (stimulus.grating.double)
  if stimulus.grating.contrast(1) > 0
    Screen('DrawTexture', myscreen.w, stimulus.grating.ltex(stimulus.grating.orientNum(1),stimulus.grating.phaseNum(1)),[],stimulus.grating.ldestrect);
  end
  if stimulus.grating.contrast(2) > 0
    Screen('DrawTexture', myscreen.w, stimulus.grating.rtex(stimulus.grating.orientNum(2),stimulus.grating.phaseNum(2)),[],stimulus.grating.rdestrect);
  end
else
  if stimulus.grating.contrast > 0
    Screen('DrawTexture', myscreen.w, stimulus.grating.tex(stimulus.grating.orientNum,stimulus.grating.phaseNum));
  end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the orientation and phase to display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grating = setGratingParameters(grating,orient,phase,contrast,myscreen)

if (grating.double)
  for i = 1:length(orient)
    grating.orientNum(i) = find(grating.orientation == orient(i));
    grating.phaseNum(i) = find(grating.phase == phase(i));
    grating.contrast(i) = contrast(i);
  end
else
  grating.orientNum = find(grating.orientation == orient);
  grating.phaseNum = find(grating.phase == phase);
  grating.contrast = contrast;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     STIMULUS Gratings end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to handle observer response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial = trialresponse(buttons,trial)

% do nothing if we already got response
if (trial.gotresponse)
  return
end
trial.gotresponse = 1;
% get the reaction time
trial.reactiontime = GetSecs-trial.segstart;

% this should check subject response
% see if user is correct or not
%if (buttons(trial.sig) && ~buttons(1+mod(trial.sig,2)))
%  trial.correct = 1;
%  % set the color of the fixation spot
% trial.fixcolor{trial.thisseg} = [0 255 0];
%else
%  trial.correct = 0;
%  % set the color of the fixation spot
%  trial.fixcolor{trial.thisseg} = [255 0 0];
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to set stimulus parameters at
% the beginning of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = setStimulusParameters(stimulus,task,myscreen)

if isfield(stimulus,'dots')
  % set random direction
  if (task.thistrial.randomizeDirection)
    task.thistrial.direction = 360*rand(1,2);
  end
  % set the parameters
  stimulus.dots = setdirspeed(stimulus.dots,d2r(task.thistrial.direction)',task.thistrial.speed',myscreen);
end

if isfield(stimulus,'grating')
  % contrast reversing phase change
  if mod(task.thistrial.thisseg,2)
    task.thistrial.phase(1:2) = stimulus.grating.phase(end);
  else
    task.thistrial.phase(1:2) = stimulus.grating.phase(end/2);
  end
  % get a random phase
  %thisphase(1) = stimulus.grating.phase(ceil(rand*length(stimulus.grating.phase)));
  %thisphase(2) = stimulus.grating.phase(ceil(rand*length(stimulus.grating.phase)));
  % get a random orientation if called for
  if (task.thistrial.randomizeOrientation)
    randorient(1) = stimulus.grating.orientation(ceil(rand(1)*length(stimulus.grating.orientation)));
    randorient(2) = stimulus.grating.orientation(ceil(rand(1)*length(stimulus.grating.orientation)));
    stimulus.grating = setGratingParameters(stimulus.grating,randorient,task.thistrial.phase,task.thistrial.contrast,myscreen);
  else
    % set the orientation and phase for this segment
    stimulus.grating = setGratingParameters(stimulus.grating,task.thistrial.orientation,task.thistrial.phase,task.thistrial.contrast,myscreen);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to display stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stimulus task] = trialstimulus(stimulus,task,myscreen)


% display dots
if isfield(stimulus,'dots')
  % and update the stimulus
  stimulus.dots = eval(sprintf('updatedots_%s_%s(stimulus.dots,task.thistrial.coherence,stimulus.dots.same,myscreen);',stimulus.dots.type,stimulus.dots.dotmotion));
  eval(sprintf('dispdots_%s(myscreen,stimulus.dots);',stimulus.dots.dotmotion));
end

% display grating
if isfield(stimulus,'grating')
  % display the grating only on odd segmens
  %if mod(task.thistrial.thisseg,2)
  updateGrating(stimulus,myscreen);
  %end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end trial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function task = endtrial(stimulus,task,myscreen)

% this is called at end of trial, can be
% used for updating staircase

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stimulus myscreen] = initstimulus(stimulus,myscreen)

% init dots
if isfield(stimulus,'dots')
  % check whether stimulus dots are appropriate
  if (isempty(which(sprintf('initdots_%s',stimulus.dots.dotmotion))) | ...
      isempty(which(sprintf('updatedots_%s_%s',stimulus.dots.dottype,stimulus.dots.dotmotion))))
    disp(sprintf('UHOH: dot type %s with %s not implemented',stimulus.dots.dotmotion,stimulus.dots.dottype));
    return
  end

  % init the dots 
  eval(sprintf('stimulus.dots = initdots_%s(myscreen,stimulus.dots);',stimulus.dots.dotmotion));
  
  % tick the dots once, to get them started
  dots = eval(sprintf('updatedots_%s_%s(stimulus.dots,[0 0],%i,myscreen);',stimulus.dots.type,stimulus.dots.dotmotion,stimulus.dots.same));
end

% init gratings
if isfield(stimulus,'grating')
  stimulus = initGrating(stimulus,myscreen);
end

% remember that we have been initialized
stimulus.init = 1;
