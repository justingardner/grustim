%
%        $Id: $
%      usage: dotscat
%         by: justin gardner
%       date: 06/20/2012
%    purpose: direction category experiment 
%
function myscreen = dotscat(varargin)

% set input arguments
getArgs(varargin,{'subjectID=s9999','centerX=10','centerY=0','diameter=16'});

global stimulus
stimulus = [];

% init screen
myscreen.subjectID = subjectID;
myscreen.saveData = 1;
myscreen = initScreen(myscreen);

% set our task to have two phases
% category boundary phase
task{1}{1}.waitForBacktick = 1;
task{1}{1}.seglen = 4;
task{1}{1}.numBlocks = 1;

% dot categorization task
directions = [60:60:419];
task{1}{2}.segmin = [1 1 6 1 2];
task{1}{2}.segmax = [11 1 6 1 2];
task{1}{2}.synchToVol = [1 0 0 0 0 ];
task{1}{2}.getResponse = [0 0 0 0 1];
task{1}{2}.waitForBacktick = 1;
task{1}{2}.parameter.direction1 = [directions;directions];
task{1}{2}.parameter.direction2 = [directions;directions];
task{1}{2}.random = 1;
task{1}{2}.randVars.calculated.correctIncorrect = nan;
task{1}{2}.randVars.calculated.matchNonmatch = nan;
task{1}{2}.randVars.calculated.oddDirection11 = nan;
task{1}{2}.randVars.calculated.oddDirection12 = nan;
task{1}{2}.randVars.calculated.oddDirection21 = nan;
task{1}{2}.randVars.calculated.oddDirection22 = nan;
task{1}{2}.randVars.calculated.evenDirection11= nan;
task{1}{2}.randVars.calculated.evenDirection12 = nan;
task{1}{2}.randVars.calculated.evenDirection21 = nan;
task{1}{2}.randVars.calculated.evenDirection22 = nan;
task{1}{2}.numBlocks = 1
task{1}{2}.numTrials = 2000

% init task
for phaseNum = 1:length(task{1})
[task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback,@initTrialCallback);
end

% register stimulus in myscreen
global stimulus;
myscreen = initStimulus('stimulus',myscreen);

% do our initialization which creates dots
stimulus = myInitStimulus(stimulus,myscreen,task,centerX,centerY,diameter);

% create the direction categories
stimulus.category1 = directions(directions < 210);
stimulus.category2 = directions(directions >= 210);

% which stimulus is the one to do the discrimination on (1 is left, 2 is right)
stimulus.categorySide = 2;
stimulus.notCategorySide = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
  % update the task
  [task{1} myscreen phaseNum] = updateTask(task{1},myscreen,phaseNum);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   initTrialCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = initTrialCallback(task, myscreen)

% get a random direction sequence for block design
global stimulus;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;

% set the fixation color
if any(task.thistrial.thisseg == [2 4])
  stimulus.fixColor = stimulus.fixColors.stim;
elseif task.thistrial.thisseg == 5
  stimulus.fixColor = stimulus.fixColors.response;
else
  stimulus.fixColor = stimulus.fixColors.interStim;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus;

% clear the screen
mglClearScreen;
    
% draw category boundary
if task.thistrial.thisphase == 1
    mglLines2(stimulus.dots{2}.xcenter+((stimulus.dots{2}.rmax)*(cos(7*(pi)/6))), stimulus.dots{2}.ycenter+ ((stimulus.dots{2}.rmax)*(sin(7*(pi)/6))), stimulus.dots{2}.xcenter+((stimulus.dots{2}.rmax)*(cos((pi)/6))), stimulus.dots{2}.ycenter+((stimulus.dots{2}.rmax)*(sin((pi)/6))), 4, [1 1 1]);
end

if any(task.thistrial.thisseg == [2 4])
  % use the different directions for each segment
  if task.thistrial.thisseg == 2
    direction = task.thistrial.direction1;
  else
    direction = task.thistrial.direction2;
  end
  
  % get coherence for this phase
  coherence = 1;

  % select the stencil to make the dot patterns round patches
  mglStencilSelect(1);
  
  %if stimulus.categorySide == 1
    % draw the left patch
    stimulus.dots{1}.dir = direction(1);
    stimulus.dots{1} = updateDots(stimulus.dots{1},coherence,myscreen);
  %else
    % draw the right patch
    stimulus.dots{2}.dir = direction(2);
    stimulus.dots{2} = updateDots(stimulus.dots{2},coherence,myscreen);
  %end
  % return to unstenciled drawing
  mglStencilSelect(0);
  
  
end
 
% fixation cross
mglFixationCross(1,2,stimulus.fixColor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets subject  response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task, myscreen)

global stimulus;
stimulus.fixColor = [0 0 1];

% get the directions that were presented
dir1 = task.thistrial.direction1(stimulus.categorySide);
dir2 = task.thistrial.direction2(stimulus.categorySide);
dir3 = task.thistrial.direction1(stimulus.notCategorySide);
dir4 = task.thistrial.direction2(stimulus.notCategorySide);

% get the "category" for each direction. 
cat1 = 1+any(dir1 == stimulus.category2);
cat2 = 1+any(dir2 == stimulus.category2);

% match or nonmatch
catMatch = (cat1==cat2);
if catMatch,catMatchStr = 'match';else,catMatchStr = 'nonmatch';end

% see if the subject got it right
if (catMatch && (task.thistrial.whichButton==1)) || (~catMatch && (task.thistrial.whichButton==2))
  correctIncorrect = true;
  correctIncorrectStr = 'correct';
  stimulus.fixColor = [0 1 0];
else
  correctIncorrect = false;
  correctIncorrectStr = 'incorrect';
  stimulus.fixColor = [1 0 0];
end

% display 
disp(sprintf('(dotscat) %i: %0.2f vs %0.2f (cat:%i vs cat:%i) %s - %s (%i)',task.trialnum,dir1,dir2,cat1,cat2,catMatchStr,correctIncorrectStr,task.thistrial.whichButton));

% keep the correct incorrect in a calculated variable
task.thistrial.correctIncorrect = correctIncorrect;

% store matchNonmatch as a calculated variable (0=nonmatch, 1=category1
% match, 2=category2 match)
if any(dir1 == stimulus.category2) & any(dir2 == stimulus.category2);
task.thistrial.matchNonmatch = 2
end
if any(dir1==stimulus.category1)& any(dir2==stimulus.category1);
    task.thistrial.matchNonmatch = 1
end
if ~catMatch;
    task.thistrial.matchNonmatch=0;
end


% sort trials by odd and even and store odd even as a calculated variable
% (note: randVars direction label is hardcoded, but linked to categorySide
% or nonCategorySide. Pay close attention to this if changing side that
% categorization takes place on)
if isodd(task.trialnum);
task.thistrial.oddDirection12 = dir1;
task.thistrial.oddDirection11 = dir3;
task.thistrial.oddDirection22 = dir2;
task.thistrial.oddDirection21 = dir4;
task.thistrial.evenDirection12 = 0;
task.thistrial.evenDirection11 = 0;
task.thistrial.evenDirection22 = 0;
task.thistrial.evenDirection21 = 0;
end

if iseven(task.trialnum);
task.thistrial.evenDirection12 = dir1;
task.thistrial.evenDirection11 = dir3;
task.thistrial.evenDirection22 = dir2;
task.thistrial.evenDirection21 = dir4;
task.thistrial.oddDirection12 = 0;
task.thistrial.oddDirection11 = 0;
task.thistrial.oddDirection22 = 0;
task.thistrial.oddDirection21 = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen,task,centerX,centerY,diameter)

% use circular aperature or not
stimulus.circularAperture = 1;

% use single side
stimulus.singleAperture = 0;

% size of stimulus
stimulus.width = 0.5*floor(myscreen.imageHeight/0.5);
stimulus.height = stimulus.width;

% stencil locations
stimulus.stencilAngle = 10;
stencilRadius = max(stimulus.width,stimulus.height)+1;

% screen places
screenRight = myscreen.imageWidth/2;
screenTop = myscreen.imageHeight/2;
screenBottom = -myscreen.imageWidth/2;
screenLeft = -myscreen.imageHeight/2;

% make right cutout
x = [0 cos(d2r(90-stimulus.stencilAngle))*stencilRadius screenRight ...
     screenRight cos(d2r(270+stimulus.stencilAngle))*stencilRadius 0];
y = [0 sin(d2r(90-stimulus.stencilAngle))*stencilRadius screenTop ...
     screenBottom sin(d2r(270+stimulus.stencilAngle))*stencilRadius 0];

% or size of circular cutot
if stimulus.circularAperture
  fixDiskSize = 1;
  distFromEdge = 0.5;
  if ~isempty(diameter)
    circleSize = diameter;
  else
    circleSize = (myscreen.imageWidth/2) - fixDiskSize - distFromEdge;
  end
  if ~isempty(centerX)
    stimulus.centerX = centerX;
  else
    stimulus.centerX = fixDiskSize+circleSize/2;
  end
  if ~isempty(centerY)
    stimulus.centerY = centerY;
  else
    stimulus.centerY = 0;
  end
  circleSize = [circleSize circleSize];
else
  stimulus.centerX = 0;
end

% right stencil
mglStencilCreateBegin(1);
if stimulus.circularAperture
  mglFillOval(stimulus.centerX,stimulus.centerY,circleSize,[1 1 1]);
  mglGluDisk(stimulus.centerX,stimulus.centerY,circleSize(1)/2,[1 1 1],128);
else
  % for the single grating
  if stimulus.singleAperture
    mglGluDisk(0,0,myscreen.imageWidth,[1 1 1],128);
  else
    mglPolygon(x,y,[1 1 1]);
  end
end
%mglStencilCreateEnd;

% make left stencil
%mglClearScreen;
%mglStencilCreateBegin(2);
if stimulus.circularAperture
  mglGluDisk(-stimulus.centerX,stimulus.centerY,circleSize(1)/2,[1 1 1],128);
else
  if ~stimulus.singleAperture
    mglPolygon(-x,y,[1 1 1]);
  end
end
mglStencilCreateEnd;

% create dots stimulus here, easier to do once we have stencil info
stimulus.dots = {};
dots.rmax = circleSize(1)/2;
dots.xcenter = -stimulus.centerX;
dots.ycenter = stimulus.centerY;
stimulus.dots{1} = initDots(myscreen,dots);
dots.xcenter = stimulus.centerX;
dots.ycenter = stimulus.centerY;
stimulus.dots{2} = initDots(myscreen,dots);
mglClearScreen(0);mglFlush;
mglClearScreen(0);mglFlush;

% set fixation colors
stimulus.fixColors.stim = [0 1 1];
stimulus.fixColors.interStim = [1 1 1];
stimulus.fixColors.response = [1 1 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initDots(myscreen,dots)

% init dots structure
if ieNotDefined('dots'),dots = [];end

% convert the passed in parameters to real units
if ~isfield(dots,'type'), dots.type = 'Linear';,end
if ~isfield(dots,'rmax'), dots.rmax = max(myscreen.imageWidth,myscreen.imageHeight)/2;,end
if ~isfield(dots,'xcenter'), dots.xcenter = 0;,end
if ~isfield(dots,'ycenter'), dots.ycenter = 0;,end
if ~isfield(dots,'dotsize'), dots.dotsize = 4;,end
if ~isfield(dots,'density'), dots.density = 5;,end
if ~isfield(dots,'coherence'), dots.coherence = 1;,end
if ~isfield(dots,'speed'), dots.speed = 6;,end
if ~isfield(dots,'dir'), dots.dir = 0;,end

% init the dots
dots = feval(sprintf('initDots%s',dots.type),dots,myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create dots for linear2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initDotsLinear(dots,myscreen)

% actually a square patch of dots that get stenciled
% so calculate width and height
dots.width = dots.rmax*2;
dots.height = dots.rmax*2;

% get the number of dots
dots.n = round(dots.width*dots.height*dots.density);

% get max and min points for dots
dots.xmin = -dots.width/2;
dots.xmax = dots.width/2;
dots.ymin = -dots.height/2;
dots.ymax = dots.height/2;

% get initial position
dots.x = rand(1,dots.n)*dots.width;
dots.y = rand(1,dots.n)*dots.height;

% get the step size
dots.stepsize = dots.speed/myscreen.framesPerSecond;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updateDots(dots,coherence,myscreen)

% update the dots
dots = feval(sprintf('updateDots%s',dots.type),dots,coherence,myscreen);

% draw the dots
mglPoints2(dots.x+dots.xcenter,dots.y+dots.ycenter,dots.dotsize,[1 1 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots for linear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updateDotsLinear(dots,coherence,myscreen)

% get the dots step
dots.xstep = cos(pi*dots.dir/180)*dots.stepsize;
dots.ystep = sin(pi*dots.dir/180)*dots.stepsize;

% pick a random set of dots
dots.coherent = rand(1,dots.n) < coherence;

% now move those dots in the right direction
dots.x(dots.coherent) = dots.x(dots.coherent)+dots.xstep;
dots.y(dots.coherent) = dots.y(dots.coherent)+dots.ystep;

% randomwalk rule
thisdir = rand(1,sum(~dots.coherent))*2*pi;
dots.x(~dots.coherent) = dots.x(~dots.coherent)+cos(thisdir)*dots.stepsize;
dots.y(~dots.coherent) = dots.y(~dots.coherent)+sin(thisdir)*dots.stepsize;

% movshon noise
%dots.x(~dots.coherent) = rand(1,sum(~dots.coherent))*dots.width;
%dots.y(~dots.coherent) = rand(1,sum(~dots.coherent))*dots.height;

% make sure we haven't gone off the patch
% do the dots separately for left and right hand side
dots.x(dots.x < dots.xmin) = dots.x(dots.x < dots.xmin)+dots.width;
dots.x(dots.x > dots.xmax) = dots.x(dots.x > dots.xmax)-dots.width;
dots.y(dots.y < dots.ymin) = dots.y(dots.y < dots.ymin)+dots.height;
dots.y(dots.y > dots.ymax) = dots.y(dots.y > dots.ymax)-dots.height;



