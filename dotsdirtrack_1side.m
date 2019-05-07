%
%        $Id: $
%      usage: dotsdirtrck_1side
%         by: Josh Ryu
%       date: 02/12/2019
%    purpose: 

function myscreen = dotsdirtrack_1side(varargin)

% set input arguments
%getArgs(varargin,{'subjectID=s999','centerX=10','centerY=0','diameter=16'});
getArgs(varargin,{'subjectID=-1','centerX=10','centerY=0','diameter=16'});

% set up screen
myscreen.subjectID = subjectID;
myscreen.saveData = 1;
%myscreen.displayName = 'screen1';
%myscreen.displayName = 'testVpixx';
%myscreen.displayName = 'test'; myscreen.screenNumber = 1;
myscreen = initScreen(myscreen);

%% parameters

% Task design
% S1: Stimulus (30s)
% S2: Fixation (3s)
task{1}{1}.segmin = [30 3]; 
task{1}{1}.segmax = [30 3]; 
task{1}{1}.numTrials = 10;
task{1}{1}.getResponse = [1 0]; %segment to get response.
task{1}{1}.waitForBacktick = 1; %wait for backtick before starting each trial 

% random task parameters
task{1}{1}.parameter.respSide = [0 1]; %Is the response side left (1) or right (0)?
task{1}{1}.randVars.uniform.initDir = [0:45:359]; %initial directions
task{1}{1}.parameter.dirNoise = 1; %wait for backtick before starting each trial 

coherence = [1 0.8 0.6 0.4 0.2 0.1]; %1; %[1 0.8 0.6 0.4 0.2 0.1]; %[1 0.8 0.6 0.4 0.2 0.1];

% calculated parameters
task{1}{1}.randVars.calculated.respAngle = nan; %response angle
task{1}{1}.randVars.calculated.stimDir = nan;
task{1}{1}.randVars.calculated.trackTick = nan;
task{1}{1}.randVars.calculated.trackStim = nan(ceil(task{1}{1}.segmax(1)*myscreen.framesPerSecond),1);
task{1}{1}.randVars.calculated.trackResp = nan(ceil(task{1}{1}.segmax(1)*myscreen.framesPerSecond),1);

%% initialize
for phaseN = 1:length(coherence)
    task{1}{phaseN} = task{1}{1};
    task{1}{phaseN}.parameter.coherence = coherence(phaseN);
    [task{1}{phaseN} myscreen] = initTask(task{1}{phaseN},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback,@initTrialCallback);
end

% initialize stimulus
global stimulus;
stimulus = [];

myscreen = initStimulus('stimulus',myscreen); % what does this do???
stimulus = myInitStimulus(stimulus,myscreen,task,centerX,centerY,diameter); %centerX,Y, diameter called by getArgs.
stimulus.powerwheel = 0; %1; % powerwheel (1)  or mouse (0)

stimulus.grabframe = 0; %save frames into matrices
if stimulus.grabframe
    global frame
    frame = {};
end

phaseNum = 1;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);     % update the task
    myscreen = tickScreen(myscreen,task);     % flip screen
end

myscreen = endTask(myscreen,task);
mglClose
mglDisplayCursor(1)

if stimulus.grabframe
    save('/Users/joshryu/Dropbox/GardnerLab/data/FYP/dotsdircont/frame.mat', 'frame')
end

end

%% Initialize trials 
function [task myscreen] = initTrialCallback(task, myscreen)
    global stimulus    
end

%% Start segment
function [task myscreen] = startSegmentCallback(task, myscreen)
% S1: Stimulus (30s)
% S2: Fixation (3s)
global stimulus 

%change stimulus accordingly
if task.thistrial.thisseg == 1
    stimulus.fixColor = stimulus.fixColors.response;    
    
    %set the initial response to the stimulus direction
    task.thistrial.respAngle = task.thistrial.initDir; 
    task.thistrial.stimDir = task.thistrial.initDir; 
    
    % set mouse position to the stimulus direction. 
    if stimulus.powerwheel
        %each x movement (+1 unit) counts as -1/2 degree turn (note the
        %screensize when doing this)
        x_screen = ceil(myscreen.screenWidth/2-task.thistrial.respAngle*2);
        mglSetMousePosition(x_screen,floor(myscreen.screenHeight/2), myscreen.screenNumber);
        
        %correct the angle by the set mouse position
        task.thistrial.respAngle = mod(-(x_screen-myscreen.screenWidth/2)/2,360); %(in degrees)
    else
        theta = mod(task.thistrial.respAngle/360*2*pi,2*pi);   
        x_img = 5*cos(theta); y_img = 5*sin(theta);
        x_screen = x_img*myscreen.screenWidth/myscreen.imageWidth + myscreen.screenWidth/2;
        y_screen = y_img*myscreen.screenHeight/myscreen.imageHeight + myscreen.screenHeight/2;
        mglSetMousePosition(ceil(x_screen),floor(y_screen), myscreen.screenNumber); %identify main screen?
        
        % note that the mouse position is an approximation of the initial angle
        % it is the top left grid point from the circle of radius 5, with
        % the given angle.
        % thus we correct the angle here: 
        distx = (ceil(x_screen)-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
        disty = (floor(y_screen)-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight; % what is imagewidth???
        task.thistrial.respAngle = atan2(disty,distx)/(2*pi)*360;
    end
    
    mglDisplayCursor(0) %hide cursor
    
    task.thistrial.trackTick = 0;
    
else %intertrial interval or feedback
    stimulus.fixColor = stimulus.fixColors.interStim;
end    

end

%% screen update
function [task myscreen] = screenUpdateCallback(task, myscreen)
% S1: Stimulus (30s)
% S2: Fixation (3s)

mglClearScreen % clear screen

global stimulus % call stimulus

if task.thistrial.thisseg == 1 % [stimulus period] draw dots   
    % get variables for this task
    coherence = task.thistrial.coherence; %task.thistrial.coherence; 
    direction = task.thistrial.stimDir;

    % choose stencil with L&R holes
    mglStencilSelect(1);
    
    if task.thistrial.respSide == 1
    % draw the left patch
        stimulus.dots{1}.dir = direction;
        stimulus.dots{1} = updateDots(stimulus.dots{1},coherence,myscreen);
    else
    % draw the right patch
        stimulus.dots{2}.dir = direction;
        stimulus.dots{2} = updateDots(stimulus.dots{2},coherence,myscreen);
    end
    % return to unstenciled drawing
    mglStencilSelect(0);
    
    [task myscreen] = drawCenterCue(task,myscreen,1);

    % draw response
    [task, myscreen] = getTurnResponse(task, myscreen);
    
    % record direction and response
    task.thistrial.trackTick = task.thistrial.trackTick +1;
    task.thistrial.trackResp(task.thistrial.trackTick) = task.thistrial.respAngle;
    task.thistrial.trackStim(task.thistrial.trackTick) = task.thistrial.stimDir;
    
    % update direction. Noise is distributed normally in angle space. 
    task.thistrial.stimDir = mod(task.thistrial.stimDir + normrnd(0,task.thistrial.dirNoise),360);
else
    % draw fixation
    mglGluAnnulus(0,0,0.5,0.75,stimulus.fixColor,60,1);
end

% draw circle for stimulus patches 
mglGluAnnulus(-stimulus.centerX,stimulus.centerY,...
    stimulus.circleSize(1)/2,stimulus.circleSize(1)/2+0.1,[1 1 1],100,1)
mglGluAnnulus(stimulus.centerX,stimulus.centerY,...
    stimulus.circleSize(1)/2,stimulus.circleSize(1)/2+0.1,[1 1 1],100,1)

if stimulus.grabframe
    global frame; frame{task.thistrial.thisseg} = mglFrameGrab;
end

end

function [task myscreen] = drawCenterCue(task,myscreen,isStim)
global stimulus         
    mglGluAnnulus(0,0,0.5,0.75,stimulus.fixColor,60,1);
end

function [task myscreen]  = getTurnResponse(task, myscreen)
global stimulus % call stimulus        
    if stimulus.powerwheel
        mInfo = mglGetMouse(myscreen.screenNumber); %each movement by x=1 moves the  cursor 1/2 degrees
        nextrespangle = -(mInfo.x-myscreen.screenWidth/2)/2*2*pi/360; %(in radians)
        nextrespangle_rad = mod(nextrespangle,2*pi);
        nextrespangle_deg = mod(nextrespangle_rad/(2*pi)*360,360);
        
        if mInfo.x < 0 || mInfo.x > myscreen.screenWidth
            mglSetMousePosition(ceil(myscreen.screenWidth/2 - nextrespangle_deg*2),floor(myscreen.screenHeight/2), myscreen.screenNumber);
            mglDisplayCursor(0)% hide cursor
        end
    else
        % mglSetMousePosition(ceil(x_screen),floor(y_screen)); %identify main screen?
        mInfo = mglGetMouse(myscreen.screenNumber);
        distx = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
        disty = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight; % what is imagewidth???
        nextrespangle = atan2(disty,distx);
        nextrespangle_rad = mod(nextrespangle,2*pi);
        nextrespangle_deg = mod(nextrespangle_rad/(2*pi)*360,360);
    end
       
    % note that respAngle is stored in *real* angles -- so that it
    % corresponds correctly to the direction task. This means that when you
    % transform into visual space you need to flip into MGL angles, see
    % mglGluDiskAnnulus_ which does this step
    
    task.thistrial.respAngle = nextrespangle_deg;
    theta = mod(task.thistrial.respAngle/360*2*pi,2*pi);   
    
    mglLines2(0.9*cos(theta), 0.9*sin(theta), 1.4*cos(theta), 1.4*sin(theta),5,stimulus.fixColors.response)
end

%% Get response 
function [task myscreen] = responseCallback(task, myscreen)

global stimulus
 
end

%% Initialize stimulus
function stimulus = myInitStimulus(stimulus,myscreen,task,centerX,centerY,diameter)  
  % stimulus field is square. 
  stimulus.width = 0.5*floor(myscreen.imageHeight/0.5); %half of the screen (for L/R side)
  stimulus.height = stimulus.width;
  
  % make stencils
  stimulus.stencilAngle = 10;
  stencilRadius = max(stimulus.width,stimulus.height)+1;
  
  fixDiskSize = 1; 
  distFromEdge = 0.5;
  
  if ~isempty(diameter), circleSize = diameter;   % select circle size
  else circleSize = (myscreen.imageWidth/2) - fixDiskSize - distFromEdge; end  
  circleSize = [circleSize circleSize];
  stimulus.circleSize = circleSize;
  
  if ~isempty(centerX), stimulus.centerX = centerX; %select X center (pass on to stimulus struct)
  else stimulus.centerX = fixDiskSize+circleSize/2; end
  if ~isempty(centerY), stimulus.centerY = centerY; % select Y center
  else stimulus.centerY = 0; end
  
  mglStencilCreateBegin(1);
  mglGluDisk(stimulus.centerX,stimulus.centerY,circleSize(1)/2,[1 1 1],128);
  mglGluDisk(-stimulus.centerX,stimulus.centerY,circleSize(1)/2,[1 1 1],128);
  mglStencilCreateEnd
  
  % create patches of dots
  stimulus.dots = {};
  dots.rmax = circleSize(1)/2;
  dots.xcenter = -stimulus.centerX; %left patch
  dots.ycenter = stimulus.centerY;
  stimulus.dots{1} = initDots(myscreen,dots);
  dots.xcenter = stimulus.centerX; % right patch
  dots.ycenter = stimulus.centerY;
  stimulus.dots{2} = initDots(myscreen,dots);
  mglClearScreen(0);mglFlush;

  % fixation cross
  stimulus.fixColors.stim = [0 1 1];
  stimulus.fixColors.interStim = [1 1 1];
  stimulus.fixColors.response = [1 1 0];
end

%% initialize dots
function dots = initDots(myscreen,dots)
    if ieNotDefined('dots'),dots = [];end

    % convert the passed in parameters to real units
    if ~isfield(dots,'type'), dots.type = 'Linear';,end
    if ~isfield(dots,'rmax'), dots.rmax = max(myscreen.imageWidth,myscreen.imageHeight)/2;,end %radius to fill the entire screen unit of imagewidth?
    if ~isfield(dots,'xcenter'), dots.xcenter = 0;,end %define centers 
    if ~isfield(dots,'ycenter'), dots.ycenter = 0;,end
    if ~isfield(dots,'dotsize'), dots.dotsize = 4;,end
    if ~isfield(dots,'density'), dots.density = 5;,end
    if ~isfield(dots,'coherence'), dots.coherence = 1;,end
    if ~isfield(dots,'speed'), dots.speed = 6;,end
    if ~isfield(dots,'dir'), dots.dir = 0;,end

    % define a square patch
    dots.width = dots.rmax*2;
    dots.height = dots.rmax*2;
    dots.xmin = -dots.width/2;
    dots.xmax = dots.width/2;
    dots.ymin = -dots.height/2;
    dots.ymax = dots.height/2;
    
    % number of dots 
    dots.n = round(dots.width*dots.height*dots.density);

    % initialize positions, uniform distribution. 
    dots.x = rand(1,dots.n)*dots.width;
    dots.y = rand(1,dots.n)*dots.height;

    % get the step size
    dots.stepsize = dots.speed/myscreen.framesPerSecond;
end

%% Update/plot dots. 
function dots = updateDots(dots,coherence,myscreen)   
    % convert steps to direction gradients
    dots.xstep = cos(pi*dots.dir/180)*dots.stepsize;
    dots.ystep = sin(pi*dots.dir/180)*dots.stepsize;

    % pick a random set of coherent dots
    dots.coherent = rand(1,dots.n) < coherence;

    % Move coherent dots
    dots.x(dots.coherent) = dots.x(dots.coherent)+dots.xstep;
    dots.y(dots.coherent) = dots.y(dots.coherent)+dots.ystep;

    % Move incoherent dots, randomwalk rule
    thisdir = rand(1,sum(~dots.coherent))*2*pi;
    dots.x(~dots.coherent) = dots.x(~dots.coherent)+cos(thisdir)*dots.stepsize;
    dots.y(~dots.coherent) = dots.y(~dots.coherent)+sin(thisdir)*dots.stepsize;
    
    % bring the dots back into the patch. 
    dots.x(dots.x < dots.xmin) = dots.x(dots.x < dots.xmin)+dots.width;
    dots.x(dots.x > dots.xmax) = dots.x(dots.x > dots.xmax)-dots.width;
    dots.y(dots.y < dots.ymin) = dots.y(dots.y < dots.ymin)+dots.height;
    dots.y(dots.y > dots.ymax) = dots.y(dots.y > dots.ymax)-dots.height;
    
    % draw the dots
    mglPoints2(dots.x+dots.xcenter,dots.y+dots.ycenter,dots.dotsize,[1 1 1]);
end
