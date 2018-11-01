%
%        $Id: $
%      usage: dotsdiraflr
%         by: Josh Ryu
%       date: 10/22/2018
%    purpose: Compare motion direction in left or right patches of dots. 
%
% Uses staircase. 
% staircase termination might overlap. 

function myscreen = dotsdiraflr(varargin)

% set input arguments
%getArgs(varargin,{'subjectID=s999','centerX=10','centerY=0','diameter=16'});
getArgs(varargin,{'centerX=10','centerY=0','diameter=16'});

% set up screen
myscreen.subjectID = subjectID;
myscreen.saveData = 1;
myscreen.displayName = 'screen1';
myscreen = initScreen(myscreen);

% Go straight to task.
% S1: random period of fixation (1~11s)
% S2: stimulus period (1s)
% S3: repsonse period (5s)
task{1}{1}.segmin = [1 0.5 2];
task{1}{1}.segmax = [3 0.5 2];
%task{1}{1}.numBlocks = 1;
task{1}{1}.numTrials = 1000;
task{1}{1}.getResponse = [0 0 1]; %segment to get response.
task{1}{1}.waitForBacktick = 1; %wait for backtick before starting each trial 
%task{1}{1}.random = 1; %randomize order of parameter presentation. 

%task parameters
%coherence = [1 0.8 0.6];
%dirDiff = [0, 1, 5, 10]; dirDiff = [dirDiff -dirDiff(2:end)];
%task{1}{1}.parameter.coherence = coherence;
%task{1}{1}.parameter.dirDiff = dirDiff;
task{1}{1}.randVars.calculated.direction = nan; %"non-crucial" variables, block randomized. 
task{1}{1}.randVars.calculated.coherence = 1;
task{1}{1}.randVars.calculated.dirDiff = nan;
task{1}{1}.randVars.calculated.correctIncorrect = nan; %store values calculated during the task. 
task{1}{1}.randVars.calculated.leftDir = nan;
task{1}{1}.randVars.calculated.righttDir = nan;
task{1}{1}.randVars.calculated.cohPres = nan; %??? does it not save the parameters of interest automatically? 

%initialize task 
[task{1}{1} myscreen] = initTask(task{1}{1},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback,@initTrialCallback);
% task{1}{1} = task; % do we need to do this ??? 

% initialize stimulus
global stimulus;
stimulus = [];
stimulus.stairUp = 1;
stimulus.stairDown = 2;
stimulus.stairStepSize = 2;
stimulus.stairUseLevitt = 0;
stimulus.stairUsePest = 1;
stimulus.stairRep = 200; %repeat staircase every [stairRep] trials
stimulus.stairN = 0; %keeps track of how many staircases they did
stimulus.threshold(1) = 10;
stimulus.threshold(2) = 10;

%stimulus = initStaircase(stimulus); %initialize staircase. 

myscreen = initStimulus('stimulus',myscreen); % what does this do???
stimulus = myInitStimulus(stimulus,myscreen,task,centerX,centerY,diameter); %centerX,Y, diameter called by getArgs.
directions = [0:1:359]; %direction to be presented on the left side. 
stimulus.directions = directions;

phaseNum = 1;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);     % update the task
    myscreen = tickScreen(myscreen,task);     % flip screen
end

myscreen = endTask(myscreen,task);
mglClose
end

function stimulus = initStaircase(stimulus)
% set up left and right staircase
if stimulus.stairUseLevitt
    stimulus.staircase(1) = doStaircase('init','upDown','nup',stimulus.stairUp,...
        'ndown',stimulus.stairDown,'initialThreshold',stimulus.threshold(1),...
        'initialStepsize',stimulus.stairStepSize,'testType=levitt','minThreshold',0.1,'maxThreshold',45); % maximum has to be 45. 
    stimulus.staircase(2) = doStaircase('init','upDown','nup',stimulus.stairUp,...
        'ndown',stimulus.stairDown,'initialThreshold',stimulus.threshold(2),...
        'initialStepsize',stimulus.stairStepSize,'testType=levitt','minThreshold',0.1,'maxThreshold',45); % maximum has to be 45. 
elseif stimulus.stairUsePest
    stimulus.staircase(1) = doStaircase('init','upDown','nup',stimulus.stairUp,...
        'ndown',stimulus.stairDown,'initialThreshold',stimulus.threshold(1),...
        'initialStepsize',stimulus.stairStepSize,'testType=pest','minThreshold',0,'maxThreshold',45);
    stimulus.staircase(2) = doStaircase('init','upDown','nup',stimulus.stairUp,...
        'ndown',stimulus.stairDown,'initialThreshold',stimulus.threshold(2),...
        'initialStepsize',stimulus.stairStepSize,'testType=pest','minThreshold',0,'maxThreshold',45);
else
    stimulus.staircase(1) = doStaircase('init','upDown','nup',stimulus.stairUp,...
        'ndown',stimulus.stairDown,'initialThreshold',stimulus.threshold(1),...
        'initialStepsize',stimulus.stairStepSize,'minThreshold',0,'maxThreshold',45);
    stimulus.staircase(2) = doStaircase('init','upDown','nup',stimulus.stairUp,...
        'ndown',stimulus.stairDown,'initialThreshold',stimulus.threshold(2),...
        'initialStepsize',stimulus.stairStepSize,'minThreshold',0,'maxThreshold',45);
end

end


%% Initialize trials 
function [task myscreen] = initTrialCallback(task, myscreen)
    global stimulus
    
    %initialize staircase. 
    if mod(stimulus.stairN, stimulus.stairRep) == 0
        stimulus = initStaircase(stimulus);
    end
    
    [stimulus.threshold(1) stimulus.staircase(1)] = doStaircase('testValue',stimulus.staircase(1)); %left threshold
    [stimulus.threshold(2) stimulus.staircase(2)] = doStaircase('testValue',stimulus.staircase(2)); %right threshold
end

%% Start segment
function [task myscreen] = startSegmentCallback(task, myscreen)
% S1: random period of fixation (1~11s)
% S2: stimulus period (1s)
% S3: repsonse period (5s)
% S4: feedback period (1s)
global stimulus 

%change stimulus accordingly
if task.thistrial.thisseg == 2
    % decide if left side is correct
    stimulus.leftcorrect = (rand(1)<0.5); %1 if correct side is on the left. 0 if left. 
    
    % alternate between the two threshold values
    stimulus.whichthreshold = (rand(1)<0.5)+1;
    
    % choose orientation of left stimulus
    task.thistrial.direction = stimulus.directions(randperm(length(stimulus.directions),1));
    
    % choose orientation of right stimulus, relative to the left stimulus. 
    task.thistrial.dirDiff = (1-2*(stimulus.leftcorrect))*stimulus.threshold(stimulus.whichthreshold);
    
    %fixation color
    stimulus.fixColor = stimulus.fixColors.stim;
        
elseif task.thistrial.thisseg == 3
    stimulus.fixColor = stimulus.fixColors.response;
else %any(task.thistrial.thisseg == [1,3])
    stimulus.fixColor = stimulus.fixColors.interStim;
end    

end

%% screen update
function [task myscreen] = screenUpdateCallback(task, myscreen)
mglClearScreen % clear screen

global stimulus % call stimulus

% draw dots
if task.thistrial.thisseg == 2   
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
end

% draw fixation
mglFixationCross(1,2,stimulus.fixColor);
end

%% Get response 
function [task myscreen] = responseCallback(task, myscreen)

global stimulus

% record responses. correct/incorrect
if any(task.thistrial.whichButton == [1 2])
    resIsLeft = (task.thistrial.whichButton == 1); %1 if the subject chose left
    correctIncorrect = (stimulus.leftcorrect == resIsLeft); %1 if correct
    stimulus.stairN = stimulus.stairN+1; %count how many times 
else
    stimIsLeft = nan; resIsLeft = nan; correctIncorrect = nan;
end

% change color of fixation for feedback.  
if isnan(correctIncorrect)
    stimulus.fixColor = [0 0 1];
elseif correctIncorrect
    stimulus.fixColor = [0 1 0];
else
    stimulus.fixColor = [1 0 0];
end

% save variables
task.thistrial.correctIncorrect = correctIncorrect; 
task.thistrial.leftDir = task.thistrial.direction; 
task.thistrial.righttDir = task.thistrial.direction+task.thistrial.dirDiff; 
task.thistrial.randVars.calculated.cohPres = task.thistrial.coherence; %??

% Output response to the screen. 
if task.thistrial.whichButton == 1, respSide = 'left';
elseif task.thistrial.whichButton == 2, respSide = 'right'; end
if correctIncorrect == 0, corrString = 'incorrect';
elseif correctIncorrect == 1, corrString = 'correct';
else, corrString = 'no response'; end

stimulus.staircase(stimulus.leftcorrect+1) = doStaircase('update',stimulus.staircase(stimulus.leftcorrect+1),correctIncorrect,abs(task.thistrial.dirDiff));
[stimulus.threshold(stimulus.leftcorrect+1), stimulus.staircase(stimulus.leftcorrect+1)] = doStaircase('testValue',stimulus.staircase(stimulus.leftcorrect+1));

disp(['Coherence: ' num2str(task.thistrial.coherence) '; ' ...
    'Directions: ' num2str(task.thistrial.leftDir) ' (l) vs ' num2str(task.thistrial.righttDir) ' (r); ' ...
    'Difference: ' num2str(task.thistrial.dirDiff) '; '... 
    'Response: ' respSide '; ' corrString])

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
    if ~isfield(dots,'speed'), dots.speed = 10;,end
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
