%
%        $Id: $
%      usage: dotsdircont
%         by: Josh Ryu
%       date: 11/15/2018
%    purpose: Compares motion direction in left or right patches of dots. 
% Cues to one side, two side or both sides

function myscreen = dotsdircont(varargin)

% set input arguments
%getArgs(varargin,{'subjectID=s999','centerX=10','centerY=0','diameter=16'});
getArgs(varargin,{'subjectID=-1','centerX=10','centerY=0','diameter=16'});

% set up screen
myscreen.subjectID = subjectID;
myscreen.saveData = 1;
myscreen.displayName = 'screen1';
myscreen = initScreen(myscreen);

% Go straight to task.
% S1: stimulus cue period (1.5s)
% S2: stimulus period (0.5s)
% S3: repsonse period (infs)
% S4: feedback period (1.5s)
% S5: random period of fixation (1~3s)
task{1}{1}.segmin = [1.5 0.5 inf 1.5 1];
task{1}{1}.segmax = [1.5 0.5 inf 1.5 3];
%task{1}{1}.numBlocks = 1;
task{1}{1}.numTrials = 300;
task{1}{1}.getResponse = [0 0 0 0 1 0]; %segment to get response.
task{1}{1}.waitForBacktick = 1; %wait for backtick before starting each trial 
%task{1}{1}.random = 1; %randomize order of parameter presentation. 

%task parameters
%dirDiff = [0, 1, 5, 10]; dirDiff = [dirDiff -dirDiff(2:end)];`
%task{1}{1}.parameter.dirDiff = dirDiff;
%task{1}{1}.randVars.calculated.direction = nan; %"non-crucial" variables, block randomized. 
%task{1}{1}.randVars.calculated.coherence = 1;
% task{1}{1}.randVars.calculated.dirDiff = nan; having same randVars as
% parameters creates an infinite loop in calculatedRandomization.m??
% task{1}{1}.randVars.calculated.correctIncorrect = nan; %store values calculated during the task. 
task{1}{1}.randVars.calculated.leftDir = nan;
task{1}{1}.randVars.calculated.righttDir = nan;
% task{1}{1}.randVars.calculated.cohPres = nan; %??? does it not save the parameters of interest automatically? 
task{1}{1}.randVars.calculated.respAngle = nan; %response angle
task{1}{1}.randVars.calculated.respStable = nan;
%task{1}{1}.randVars.calculated.mouseStart = struct(); %how do we save not preallocated randVars???

task{1}{1}.randVars.uniform.direction = [0:1:359];
task{1}{1}.randVars.uniform.dirDiff = [-45:0.25:45];
task{1}{1}.randVars.uniform.respAngle = [0:1:359];

task{1}{1}.parameter.distAttention = [0 1]; % cue both sides?
task{1}{1}.parameter.respSide = [0 1]; %Is the response side left? (1) or right (0)
%task{1}{1}.directions = [0:1:359]; %direction to be presented on the left side. 
%task{1}{1}.dirDiff = [-45:0.25:45]; %direction to be presented on the right, relative to the left. 
task{1}{1}.parameter.coherence = [1 0.8 0.6 0.4];

[task{1}{1} myscreen] = initTask(task{1}{1},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback,@initTrialCallback);

% add a trace for the mouse position ??
[task{1}{1} myscreen] = addTraces(task{1}{1},myscreen,'mouseTrack');

% initialize stimulus
global stimulus;
stimulus = [];

myscreen = initStimulus('stimulus',myscreen); % what does this do???
stimulus = myInitStimulus(stimulus,myscreen,task,centerX,centerY,diameter); %centerX,Y, diameter called by getArgs.
stimulus.powerwheel = 0; % powerwheel (1)  or mouse (0)

phaseNum = 1;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);     % update the task
    myscreen = tickScreen(myscreen,task);     % flip screen
end

myscreen = endTask(myscreen,task);
mglClose
end

%% Initialize trials 
function [task myscreen] = initTrialCallback(task, myscreen)
    global stimulus    
end

%% Start segment
function [task myscreen] = startSegmentCallback(task, myscreen)
% S1: stimulus cue period (1.5s)
% S2: stimulus period (0.5s)
% S3: repsonse period (infs)
% S4: feedback period (1.5s)
% S5: random period of fixation (1~3s)
global stimulus 

%change stimulus accordingly
if any(task.thistrial.thisseg == [1 2])
    stimulus.fixColor = stimulus.fixColors.stim;
elseif task.thistrial.thisseg == 4
    stimulus.fixColor = stimulus.fixColors.response;    
    % set mouse position to the middle. 
    if stimulus.powerwheel
        theta = mod(task.thistrial.respAngle/360*2*pi,2*pi);   
        x = 5*cos(theta);
        mglSetMousePosition(x,myscreen.screenHeight/2);
    else
        theta = mod(task.thistrial.respAngle/360*2*pi,2*pi);   
        x_img = 5*cos(theta); y_img = 5*sin(theta);
        x_screen = x_img*myscreen.screenWidth/myscreen.imageWidth + myscreen.screenWidth/2;
        y_screen = y_img*myscreen.screenHeight/myscreen.imageHeight + myscreen.screenHeight/2;
        mglSetMousePosition(x_screen,y_screen); %identify main screen?
    end
else %intertrial interval or feedback
    stimulus.fixColor = stimulus.fixColors.interStim;
end    

end

%% screen update
function [task myscreen] = screenUpdateCallback(task, myscreen)
% S1: stimulus cue period (1.5s)
% S2: stimulus period (0.5s)
% S3: repsonse period (infs)
% S4: feedback period (1.5s)
% S5: random period of fixation (1~3s)

mglClearScreen % clear screen

global stimulus % call stimulus

if task.thistrial.thisseg == 1 %[cue period] draw cue  
    [task myscreen] = drawCenterCue(task,myscreen,1);

elseif task.thistrial.thisseg == 2 % [stimulus period] draw dots   
    % get variables for this task
    coherence = task.thistrial.coherence; %task.thistrial.coherence; 
    direction = task.thistrial.direction;
    dirDiff = task.thistrial.dirDiff;

    % choose stencil with L/R holes
    mglStencilSelect(1);
    
    % draw the left patch
    stimulus.dots{1}.dir = direction;
    stimulus.dots{1} = updateDots(stimulus.dots{1},coherence,myscreen);
  
    % draw the right patch
    stimulus.dots{2}.dir = mod(direction+dirDiff,360);
    stimulus.dots{2} = updateDots(stimulus.dots{2},coherence,myscreen);
  
    % return to unstenciled drawing
    mglStencilSelect(0);
    
    [task myscreen] = drawCenterCue(task,myscreen,1);
elseif (task.thistrial.thisseg==3) %[response segment] move the response bar. 
    [task, myscreen] = getTurnResponse(task, myscreen);
    [task myscreen] = drawCenterCue(task,myscreen,0);
elseif (task.thistrial.thisseg== 4) %[feedback segment] move the response bar. 
    theta = mod(task.thistrial.respAngle/360*2*pi,2*pi);
    mglLines2(0.9*cos(theta), 0.9*sin(theta), 1.4*cos(theta), 1.4*sin(theta),5,stimulus.fixColors.response)
    
    if task.thistrial.respSide
        theta = task.thistrial.direction*2*pi/360;
    else
        theta = (task.thistrial.direction+task.thistrial.dirDiff)*2*pi/360;
    end
    
    mglLines2(0.9*cos(theta), 0.9*sin(theta), 1.4*cos(theta), 1.4*sin(theta),5,[1 0 0 ])
    
    % draw fixation
    [task myscreen] = drawCenterCue(task,myscreen,0);
else
    % draw fixation
    mglGluAnnulus(0,0,0.5,0.75,stimulus.fixColor,60,1);
end

% draw circle for stimulus patches
mglGluAnnulus(-stimulus.centerX,stimulus.centerY,...
    stimulus.circleSize(1)/2,stimulus.circleSize(1)/2+0.1,[1 1 1],100,1)
mglGluAnnulus(stimulus.centerX,stimulus.centerY,...
    stimulus.circleSize(1)/2,stimulus.circleSize(1)/2+0.1,[1 1 1],100,1)
end

function [task myscreen] = drawCenterCue(task,myscreen,isStim)
global stimulus     


if isStim
    colors = [stimulus.fixColors.interStim', stimulus.fixColors.stim'];
    if task.thistrial.distAttention
%         mglLines2(0.9, 0, 1.4, 0, 10,stimulus.fixColors.stim)
%         mglLines2(-0.9, 0, -1.4, 0, 10,stimulus.fixColors.stim)

        mglGluAnnulus(0,0,0.5,0.75,stimulus.fixColor,60,1);
    else
    %         mglLines2(0.9*(1-2*task.thistrial.respSide), 0, ...
    %             1.4*(1-2*task.thistrial.respSide), 0,...
    %             10,stimulus.fixColors.stim)
        startAngles = [0; 180];

        if ~(task.thistrial.respSide) %if the stimulus is on the right side.
            colors = [colors(:,2), colors(:,1)];
        end

        mglGluPartialDisk([0;0],[0;0],[0.5;0.5],[0.75;0.75],startAngles,[180;180],colors,[60;60],[1;1]);
    end
else %for response cue
    colors = [stimulus.fixColors.interStim', stimulus.fixColors.response'];

%         mglLines2(0.9*(1-2*task.thistrial.respSide), 0, ...
%             1.4*(1-2*task.thistrial.respSide), 0,...
%             10,stimulus.fixColors.stim)
    startAngles = [0; 180];

    if ~(task.thistrial.respSide) %if the stimulus is on the right side.
        colors = [colors(:,2), colors(:,1)];
    end

    mglGluPartialDisk([0;0],[0;0],[0.5;0.5],[0.75;0.75],startAngles,[180;180],colors,[60;60],[1;1]);
end

end

function [task myscreen]  = getTurnResponse(task, myscreen)
global stimulus % call stimulus        
    if stimulus.powerwheel
        mInfo = mglGetMouse(myscreen.screenNumber);
        nextrespangle = -mInfo.x/90*2*pi/360; %(in radians)
        %what happens if we move offscreen???
    else
        mInfo = mglGetMouse(myscreen.screenNumber);
        distx = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
        disty = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight; % what is imagewidth???
        nextrespangle = atan2(disty,distx);
    end
    
    nextrespangle_rad = mod(nextrespangle,2*pi);
    nextrespangle_deg = mod(nextrespangle_rad/(2*pi)*360,360);
    
    %task.thistrial.mouseTrack(task.trialnum,stimulus.data.mouseTick) = task.thistrial.respAngle-stimulus.live.mouseStart;
    %task.thistrial.mouseTick = stimulus.data.mouseTick + 1;
    
    % note that respAngle is stored in *real* angles -- so that it
    % corresponds correctly to the direction task. This means that when you
    % transform into visual space you need to flip into MGL angles, see
    % mglGluDiskAnnulus_ which does this step

    if (task.thistrial.respAngle - nextrespangle_deg) < 1e-10
        task.thistrial.respStable = task.thistrial.respStable + 1;
    else
        task.thistrial.respStable = 0;
    end
    
    task.thistrial.respAngle = nextrespangle_deg;
    theta = mod(task.thistrial.respAngle/360*2*pi,2*pi);   
    
    % stimulus.data.mouseTrack(task.trialnum,stimulus.data.mouseTick) = task.thistrial.respAngle-stimulus.live.mouseStart;
    mglLines2(0.9*cos(theta), 0.9*sin(theta), 1.4*cos(theta), 1.4*sin(theta),5,stimulus.fixColors.response)
    
    if task.thistrial.respStable == 1.5*myscreen.framesPerSecond %if response stable for a second. 
        % save variables
        task.thistrial.leftDir = task.thistrial.direction; 
        task.thistrial.righttDir = task.thistrial.direction+task.thistrial.dirDiff; 

        if task.thistrial.distAttention, att = 'distributed'; else att = 'focal'; end
        if task.thistrial.respSide, respSide = 'left'; else respSide = 'right'; end

        disp(['Coherence: ' num2str(task.thistrial.coherence) '; ' ...
            'Directions: ' num2str(task.thistrial.leftDir) ' (l) vs ' num2str(task.thistrial.righttDir) ' (r); ' ...
            'Difference: ' num2str(task.thistrial.dirDiff) '; '... 
            att ' attention condition, response to ' respSide ' side at angle : ' ...
            num2str(task.thistrial.respAngle) ''])

        % end the segment when subject responses. 
        task = jumpSegment(task);
    end

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
    if ~isfield(dots,'rmax'), dots.rmax = max(myscreen.imageWidth,myscreen.imageHeight)/2;,end %radius to fill the entire screen
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
