%
%      usage: vfAttend
%         by: justin gardner & josh wilson
%       date: 10/2020
%   purpose: Contrast discrimination task with spatial attentional component.
%   Seperate staircases for each quadrant and each attentional cue
%	(8 total - 4 quadrants; distributed and side cues).




function myscreen = vfAttend(varargin);

% check arguments
if ~any(nargin == [0 1])
help taskTemplateContrast10bit
return
end

% get arguments
getArgs(varargin,{'quickThreshold=0','testXvals=[8 1 5 -4 -6 -8 2 7]','testYvals=[11 5 -6 -9 5 -6 10 7]','gradientThresh=0','randThresh=1','Eye=1','eyewindow=5'},'verbose=1');

% set up screen
screenParams = mglGetScreenParams;
stimulus.displayDistance = screenParams{1}.displayDistance*.01;
myscreen = initScreen;

% init the stimulus
clear global stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);

% these are the reserved colors, if you need them later
% you can display them by setting your color to the appropriate
% index in stimulus.colors.reservedColor e.g. to get the
% second color, in this case white, you would do
% mglClearScreen(stimulus.colors.reservedColor(2));
stimulus.colors.reservedColors = [0 0 0; 1 1 1; 0 .65 0; 1 0 0];

% eye tracking
stimulus.eye = Eye
stimulus.eyewindow = eyewindow
stimulus.eyeFrames = myscreen.framesPerSecond * 0.300; % eye movements occur when for 300 ms someone moves out of the fixation region

% grating parameters
stimulus.grating.radius = 2.25;
stimulus.grating.n = 4;
stimulus.grating.sf = 2;
stimulus.grating.tf = 2;
stimulus.grating.width = 6;
stimulus.grating.height = 6;
stimulus.grating.windowType = 'gabor'; % should be gabor or thresh
stimulus.grating.sdx = stimulus.grating.width/7;
stimulus.grating.sdy = stimulus.grating.width/7;

%staircase
stimulus.useStaircase = 1;
stimulus.initialThreshold = 7;
stimulus.initialStepsize = 1;
stimulus.minThreshold = 1;
stimulus.maxThreshold = 5;
stimulus.minStepsize = 0.05;
stimulus.maxStepsize = .75;
stimulus.restartStaircase = 1

stimulus = initStair(stimulus);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up task and initialize parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

task{1}{1}.seglen = [1 1 .5 .5 .75 1.5]; %% empty, cue, empty, stimulus, empty, response
if stimulus.eye == 1; task{1}{1}.seglen = [1 1 inf .5 .75 1.5]; end; %% wait for fixation to begin eye tracking trials
task{1}{1}.getResponse = [0 0 0 0 0 1];
task{1}{1}.parameter.quadrant = [1:4];
task{1}{1}.parameter.contrast = [.5]; %base contrast
task{1}{1}.randVars.uniform.jitterA = [-.1:.01:0];  %contrast jitter for each 4 quadrants
task{1}{1}.randVars.uniform.jitterB = [-.1:.01:0];
task{1}{1}.randVars.uniform.jitterC = [-.1:.01:0];
task{1}{1}.randVars.uniform.jitterD = [-.1:.01:0];
task{1}{1}.random = 1;
task{1}{1}.location.x = 11; %localize stimuli, masks, and attention/response cues
task{1}{1}.location.y = 8;
task{1}{1}.randVars.whichDiff = [0 1]; %% which stimuli in the quadrant will be the comparison (0 top, 1 bottom)
task{1}{1}.parameter.attend = [1 2]; %% how many segments the subject attends to = 2 segments = 1 side
task{1}{1}.parameter.contDiff = [2]; %% base contrast is divided by this value

% empty arrays to store stuff for analysis
task{1}{1}.parameter.calculated.att = nan;
task{1}{1}.randVars.calculated.resp = nan;
task{1}{1}.randVars.calculated.correct = nan
task{1}{1}.randVars.calculated.staircor = {nan;nan;nan;nan;nan;nan;nan;nan}
task{1}{1}.randVars.calculated.rt = nan;
task{1}{1}.randVars.calculated.bottom = nan;
task{1}{1}.parameter.calculated.quad = nan;
task{1}{1}.parameter.calculated.stairDiff = nan;
taskPhase.randVars.calculated.dead = nan;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialze tasks and stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stimulus = initGratings(stimulus,myscreen,task);

% initialze tasks
for phaseNum = 1:length(task{1})
[task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@updateScreenCallback,@responseCallback);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find rough starting threshold value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if quickThreshold
doThreshold(task,stimulus,randThresh,testXvals,testYvals)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if stimulus.eye
myscreen = eyeCalibDisp(myscreen);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the tasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set which phase is active
tnum = 1;

phaseNum = 1;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
% update the task
[task{1} myscreen phaseNum] = updateTask(task{1},myscreen,phaseNum);
% flip screen
myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);


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
% set the maximum contrast we can display
disp(sprintf('(taskTemplateContrast10bit:startSegmentCallback) Displaying contrast of %f',task.thistrial.contrast));

if task.thistrial.thisseg ==1
setGammaTableForMaxContrast(1);
end
stimulus.live.fixCount = 0;
stimulus.live.eyeCount = 0;
task.thistrial.dead = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to display stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task,myscreen)
global stimulus;

% clear screen to gray
mglClearScreen(stimulus.colors.grayColor);
if task.thistrial.dead && mglGetSecs(stimulus.live.deadTime)>1
task = jumpSegment(task,3);
end
% skip screen updates if you are already dead
if task.thistrial.dead
if task.thistrial.dead
mglTextSet([],32,stimulus.colors.red);
mglTextDraw('Eye Movement Detected',[0 0]);
end
return
end


%% display gratings %%
if task.thistrial.thisseg == 4


% get the contrast index
jitter(1:2) = task.thistrial.jitterA; %make indexable jitter string (same quadrants have same jitter)
jitter(3:4) = task.thistrial.jitterB;
jitter(5:6) = task.thistrial.jitterC;
jitter(7:8) = task.thistrial.jitterD;


% choose proper staircase and get comparison contrast
stairNum = pickStair(task); stimulus.stairNum = stairNum;
[testValue, stimulus.stair(stairNum)] = doStaircase('testValue', stimulus.stair(stairNum));
task.thistrial.stairDiff = testValue;

for index = 1:8;
if (2*task.thistrial.quadrant-task.thistrial.whichDiff) == index;
contrastIndex(index) = getContrastIndex((task.thistrial.contrast+jitter(index))/task.thistrial.stairDiff);
else;
contrastIndex(index) = getContrastIndex(task.thistrial.contrast+jitter(index));
end;
end

% blt texture in 8 locations - quadrants numbered clockwise from top right
% update ^ that's wrong; figure out specific numbers later
mglBltTexture(stimulus.tex(contrastIndex(1)),[task.location.x task.location.y-stimulus.grating.height/2 stimulus.grating.height]);
mglBltTexture(stimulus.tex(contrastIndex(2)),[task.location.x task.location.y+stimulus.grating.height/2 stimulus.grating.height]);
mglBltTexture(stimulus.tex(contrastIndex(3)),[task.location.x -task.location.y-stimulus.grating.height/2 stimulus.grating.height]);
mglBltTexture(stimulus.tex(contrastIndex(4)),[task.location.x -task.location.y+stimulus.grating.height/2 stimulus.grating.height]);
mglBltTexture(stimulus.tex(contrastIndex(5)),[-task.location.x -task.location.y-stimulus.grating.height/2 stimulus.grating.height]);
mglBltTexture(stimulus.tex(contrastIndex(6)),[-task.location.x -task.location.y+stimulus.grating.height/2 stimulus.grating.height]);
mglBltTexture(stimulus.tex(contrastIndex(7)),[-task.location.x task.location.y-stimulus.grating.height/2 stimulus.grating.height]);
mglBltTexture(stimulus.tex(contrastIndex(8)),[-task.location.x task.location.y+stimulus.grating.height/2 stimulus.grating.height]);


% and mask with the gaussian
mglBltTexture(stimulus.mask,[task.location.x task.location.y-stimulus.grating.height/2]);
mglBltTexture(stimulus.mask,[task.location.x task.location.y+stimulus.grating.height/2]);
mglBltTexture(stimulus.mask,[task.location.x -task.location.y-stimulus.grating.height/2]);
mglBltTexture(stimulus.mask,[task.location.x -task.location.y+stimulus.grating.height/2]);
mglBltTexture(stimulus.mask,[-task.location.x task.location.y-stimulus.grating.height/2]);
mglBltTexture(stimulus.mask,[-task.location.x task.location.y+stimulus.grating.height/2]);
mglBltTexture(stimulus.mask,[-task.location.x -task.location.y-stimulus.grating.height/2]);
mglBltTexture(stimulus.mask,[-task.location.x -task.location.y+stimulus.grating.height/2]);
end


if task.thistrial.thisseg == 6 %%report segment and identification
if task.thistrial.quadrant == 1; x = 1; y=1;elseif task.thistrial.quadrant == 2; x = 1; y=-1; elseif task.thistrial.quadrant == 3; x = -1; y=-1; elseif task.thistrial.quadrant == 4; x = -1; y=1; end;
mglFillOval(x, y, [.5 .5],  [0 0 0]);
end



if task.thistrial.thisseg == 2 %%  cue attention segment

if task.thistrial.attend == 1
if task.thistrial.quadrant == 1; x = 1; y=1;elseif task.thistrial.quadrant == 2; x = 1; y=-1; elseif task.thistrial.quadrant == 3; x = -1; y=-1; elseif task.thistrial.quadrant == 4; x = -1; y=1; end;
end

if task.thistrial.attend == 2
if task.thistrial.quadrant == 1 || task.thistrial.quadrant == 2; x = 1; y = 1; a = 1; b = -1; elseif task.thistrial.quadrant == 3 || task.thistrial.quadrant == 4; x = -1; y=1; a = -1; b = -1; end;
mglFillOval(a, b, [.5 .5],  [0 0 0]);
mglFillOval(-a, -b, [.5 .5],  [0 0 0]);
mglFillOval(-x, -y, [.5 .5],  [0 0 0]);
end

mglFillOval(x, y, [.5 .5],  [0 0 0]);

end

% fixation cross and cue boxes
mglFixationCross(2,2,[0 0 0]);
x0 = task.location.x-stimulus.grating.width/2; y0 = task.location.y-stimulus.grating.height; x1 = task.location.x+stimulus.grating.width/2; y1 = task.location.y+stimulus.grating.height;
mglLines2(x0, y0, x0, y1, 2, [0 0 0]);mglLines2(x0, y0, x1, y0, 2, [0 0 0]);mglLines2(x1, y0, x1, y1, 2, [0 0 0]);mglLines2(x0, y1, x1, y1, 2, [0 0 0]);
mglLines2(-x0, y0, -x0, y1, 2, [0 0 0]);mglLines2(-x0, y0, -x1, y0, 2, [0 0 0]);mglLines2(-x1, y0, -x1, y1, 2, [0 0 0]);mglLines2(-x0, y1, -x1, y1, 2, [0 0 0]);
mglLines2(x0, -y0, x0, -y1, 2, [0 0 0]);mglLines2(x0, -y0, x1, -y0, 2, [0 0 0]);mglLines2(x1, -y0, x1, -y1, 2, [0 0 0]);mglLines2(x0, -y1, x1, -y1, 2, [0 0 0]);
mglLines2(-x0, -y0, -x0, -y1, 2, [0 0 0]);mglLines2(-x0, -y0, -x1, -y0, 2, [0 0 0]);mglLines2(-x1, -y0, -x1, -y1, 2, [0 0 0]);mglLines2(-x0, -y1, -x1, -y1, 2, [0 0 0]);


% do eye position tracking, but only during some segments
if (stimulus.eye) && (stimulus.eyewindow>0) && task.thistrial.thisseg==3
% check eye pos
[pos,~] = mglEyelinkGetCurrentEyePos;

% compute distance
dist = hypot(pos(1),pos(2));

if task.thistrial.thisseg==3
if stimulus.live.fixCount > stimulus.eyeFrames
task = jumpSegment(task);
elseif ~any(isnan(pos))
if dist < stimulus.eyewindow
stimulus.live.fixCount = stimulus.live.fixCount + 1;
else
stimulus.live.fixCount = 0;
end
end
else
% Eye movement detection code
if (stimulus.eye) && (stimulus.eyewindow>0) && ~task.thistrial.dead
if ~any(isnan(pos))
if dist > stimulus.eyewindow && stimulus.live.eyeCount > stimulus.eyeFrames
disp('Eye movement detected!!!!');
stimulus.live.deadTime = mglGetSecs;
task.thistrial.dead = 1;
return
elseif dist > stimulus.eyewindow
stimulus.live.eyeCount = stimulus.live.eyeCount + 1;
end
end
end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)
global stimulus
% here, we just check whether this is the first time we got a response
if ~task.thistrial.gotResponse
if task.thistrial.whichDiff == 0 && task.thistrial.whichButton == 1
task.thistrial.correct == 1
elseif task.thistrial.whichDiff == 0 && task.thistrial.whichButton == 2
task.thistrial.correct == 0
elseif task.thistrial.whichDiff == 1 && task.thistrial.whichButton == 1
task.thistrial.correct == 1
elseif task.thistrial.whichDiff == 1 && task.thistrial.whichButton == 2
task.thistrial.correct == 0
end

task.thistrial.rt = task.thistrial.reactionTime;

if (task.thistrial.whichDiff == 0 && task.thistrial.whichButton == 1) || (task.thistrial.whichDiff == 1 && task.thistrial.whichButton == 2) %% press 1 for top, 2 for bottom
task.thistrial.correct = 1;
else task.thistrial.correct = 0; end;

end


task.randVars.calculated.resp = [task.randVars.calculated.resp task.thistrial.whichButton];
task.randVars.calculated.bottom = [task.randVars.calculated.bottom task.thistrial.whichDiff];
task.randVars.calculated.correct = [task.randVars.calculated.correct task.thistrial.correct];
task.randVars.calculated.staircor{stimulus.stairNum} = [task.randVars.calculated.staircor{stimulus.stairNum} task.thistrial.correct];
task.parameter.calculated.quad = [task.parameter.calculated.quad task.thistrial.quadrant] ;
task.parameter.calculated.att = [task.parameter.calculated.att task.thistrial.attend];
task.parameter.calculated.stairDiff = [task.parameter.calculated.stairDiff task.thistrial.stairDiff];

disp(sprintf('contrast %f',task.thistrial.stairDiff));
stimulus.stair(stimulus.stairNum) = doStaircase('update', stimulus.stair(stimulus.stairNum), task.thistrial.correct);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     STIMULUS Gratings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initGratings(stimulus,myscreen,task)

% set maximum color index (for 24 bit color we have 8 bits per channel, so 255)
maxIndex = 255;

% get gamma table
if ~isfield(myscreen,'gammaTable')
stimulus.linearizedGammaTable = mglGetGammaTable;
disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
disp(sprintf('(taskTemplateContrast10bit:initGratings) No gamma table found in myscreen. Contrast'));
disp(sprintf('         displays like this should be run with a valid calibration made by moncalib'));
disp(sprintf('         for this monitor.'));
disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
end
stimulus.linearizedGammaTable = myscreen.initScreenGammaTable;

disppercent(-inf,'Creating grating textures');

% calculate some colors information
%  number of reserved colors
stimulus.colors.nReservedColors = size(stimulus.colors.reservedColors,1);

% number of colors possible for gratings, make sure that we
% have an odd number
stimulus.colors.nGratingColors = maxIndex+1-stimulus.colors.nReservedColors;
if iseven(stimulus.colors.nGratingColors)
stimulus.colors.nGratingColors = stimulus.colors.nGratingColors-1;
end
% min, mid and max index of gratings colors (index values are 0 based)
stimulus.colors.minGratingIndex = maxIndex+1-stimulus.colors.nGratingColors;
stimulus.colors.midGratingIndex = stimulus.colors.minGratingIndex+floor(stimulus.colors.nGratingColors/2);
stimulus.colors.maxGratingIndex = maxIndex;
% number of contrasts we can display (not including 0 contrast)
stimulus.colors.nDisplayContrasts = floor(stimulus.colors.nGratingColors/2);

% get the color value for gray (i.e. the number between 0 and 1 that corresponds to the midGratingIndex)
stimulus.colors.grayColor = stimulus.colors.midGratingIndex/maxIndex;

% set the reserved colors - this gives a convenient value between 0 and 1 to use the reserved colors with
for i = 1:stimulus.colors.nReservedColors
stimulus.colors.reservedColor(i) = (i-1)/maxIndex;
end

% make the window through with the gratings will be displayed
gaussianWin = mglMakeGaussian(stimulus.grating.width,stimulus.grating.height,stimulus.grating.sdx,stimulus.grating.sdy);
if strcmp(stimulus.grating.windowType,'gabor')
% a gaussian window
win = maxIndex-maxIndex*gaussianWin;
else
% a simple window
win = maxIndex-maxIndex*(gaussianWin>exp(-1/2));
end
mask = ones(size(win,1),size(win,2),4)*stimulus.colors.midGratingIndex;
mask(:,:,4) = win;
stimulus.mask = mglCreateTexture(mask);

% make all the 1D gratings. We compute all possible contrast values given the
% range of indexes available to us. The 1st texture is gray the nth texture is full
% contrast for the current gamma setting
for iContrast = 0:stimulus.colors.nDisplayContrasts
disppercent(iContrast/stimulus.colors.nDisplayContrasts);
if myscreen.userHitEsc,mglClose;keyboard,end
% make the grating
thisGrating = round(iContrast*mglMakeGrating(stimulus.grating.width,nan,stimulus.grating.sf,0,0)+stimulus.colors.midGratingIndex);
% create the texture
stimulus.tex(iContrast+1) = mglCreateTexture(thisGrating);
end
disppercent(inf);



%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getContrastIndex    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function contrastIndex = getContrastIndex(desiredContrast,verbose)

if nargin < 2,verbose = 0;end

global stimulus;
if desiredContrast < 0, desiredContrast = 0;end

% now find closest matching contrast we can display with this gamma table
contrastIndex = min(round(stimulus.colors.nDisplayContrasts*desiredContrast/stimulus.currentMaxContrast),stimulus.colors.nDisplayContrasts);

% display the desired and actual contrast values if verbose is set
if verbose
actualContrast = stimulus.currentMaxContrast*(contrastIndex/stimulus.colors.nDisplayContrasts);
disp(sprintf('(getContrastIndex) Desired contrast: %0.4f Actual contrast: %0.4f Difference: %0.4f',desiredContrast,actualContrast,desiredContrast-actualContrast));
end

% out of range check
if round(stimulus.colors.nDisplayContrasts*desiredContrast/stimulus.currentMaxContrast)>stimulus.colors.nDisplayContrasts
disp(sprintf('(getContrastIndex) Desired contrast (%0.9f) out of range max contrast : %0.9f',desiredContrast,stimulus.currentMaxContrast));
keyboard
end

% 1 based indexes (0th index is gray, nDisplayContrasts+1 is full contrast)
contrastIndex = contrastIndex+1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sets the gamma table so that we can have
% finest possible control over the stimulus contrast.
%
% stimulus.reservedColors should be set to the reserved colors (for cue colors, etc).
% maxContrast is the maximum contrast you want to be able to display.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setGammaTableForMaxContrast(maxContrast)

global stimulus;
% if you just want to show gray, that's ok, but to make the
% code work properly we act as if you want to display a range of contrasts
if maxContrast <= 0,maxContrast = 0.01;end

% set the reserved colors
gammaTable(1:size(stimulus.colors.reservedColors,1),1:size(stimulus.colors.reservedColors,2))=stimulus.colors.reservedColors;

% set the gamma table
if maxContrast > 0
% create the rest of the gamma table
cmax = 0.5+maxContrast/2;cmin = 0.5-maxContrast/2;
luminanceVals = cmin:((cmax-cmin)/(stimulus.colors.nGratingColors-1)):cmax;

% now get the linearized range
redLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,luminanceVals,'linear');
greenLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,luminanceVals,'linear');
blueLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,luminanceVals,'linear');

% add these values to the table
gammaTable((stimulus.colors.minGratingIndex:stimulus.colors.maxGratingIndex)+1,:)=[redLinearized;greenLinearized;blueLinearized]';
else
% if we are asked for 0 contrast then simply set all the values to gray
gammaTable((stimulus.colors.minGratingIndex:stimulus.colors.maxGratingIndex)+1,1)=interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,0.5,'linear');
gammaTable((stimulus.colors.minGratingIndex:stimulus.colors.maxGratingIndex)+1,2)=interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,0.5,'linear');
gammaTable((stimulus.colors.minGratingIndex:stimulus.colors.maxGratingIndex)+1,3)=interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,0.5,'linear');
end

% set the gamma table
mglSetGammaTable(gammaTable);

% remember what the current maximum contrast is that we can display
stimulus.currentMaxContrast = maxContrast;



function stimulus = initStair(stimulus)

% init the staircase
stimulus.stair(1) = doStaircase('init','upDown', 'nup=1','ndown=2','initialThreshold', stimulus.initialThreshold, 'initialStepsize',stimulus.initialStepsize,'minStepsize',stimulus.minStepsize,'maxStepsize',stimulus.maxStepsize,'minThreshold',stimulus.minThreshold,'maxThreshold', stimulus.maxThreshold,'stepRule=pest','dispFig=1');
stimulus.stair(2) = doStaircase('init','upDown', 'nup=1','ndown=2','initialThreshold', stimulus.initialThreshold, 'initialStepsize',stimulus.initialStepsize,'minStepsize',stimulus.minStepsize,'maxStepsize',stimulus.maxStepsize,'minThreshold',stimulus.minThreshold,'maxThreshold', stimulus.maxThreshold,'stepRule=pest','dispFig=1');
stimulus.stair(3) = doStaircase('init','upDown', 'nup=1','ndown=2','initialThreshold', stimulus.initialThreshold, 'initialStepsize',stimulus.initialStepsize,'minStepsize',stimulus.minStepsize,'maxStepsize',stimulus.maxStepsize,'minThreshold',stimulus.minThreshold,'maxThreshold', stimulus.maxThreshold,'stepRule=pest','dispFig=1');
stimulus.stair(4) = doStaircase('init','upDown', 'nup=1','ndown=2','initialThreshold', stimulus.initialThreshold, 'initialStepsize',stimulus.initialStepsize,'minStepsize',stimulus.minStepsize,'maxStepsize',stimulus.maxStepsize,'minThreshold',stimulus.minThreshold,'maxThreshold', stimulus.maxThreshold,'stepRule=pest','dispFig=1');
stimulus.stair(5) = doStaircase('init','upDown', 'nup=1','ndown=2','initialThreshold', stimulus.initialThreshold, 'initialStepsize',stimulus.initialStepsize,'minStepsize',stimulus.minStepsize,'maxStepsize',stimulus.maxStepsize,'minThreshold',stimulus.minThreshold,'maxThreshold', stimulus.maxThreshold,'stepRule=pest','dispFig=1');
stimulus.stair(6) = doStaircase('init','upDown', 'nup=1','ndown=2','initialThreshold', stimulus.initialThreshold, 'initialStepsize',stimulus.initialStepsize,'minStepsize',stimulus.minStepsize,'maxStepsize',stimulus.maxStepsize,'minThreshold',stimulus.minThreshold,'maxThreshold', stimulus.maxThreshold,'stepRule=pest','dispFig=1');
stimulus.stair(7) = doStaircase('init','upDown', 'nup=1','ndown=2','initialThreshold', stimulus.initialThreshold, 'initialStepsize',stimulus.initialStepsize,'minStepsize',stimulus.minStepsize,'maxStepsize',stimulus.maxStepsize,'minThreshold',stimulus.minThreshold,'maxThreshold', stimulus.maxThreshold,'stepRule=pest','dispFig=1');
stimulus.stair(8) = doStaircase('init','upDown', 'nup=1','ndown=2','initialThreshold', stimulus.initialThreshold, 'initialStepsize',stimulus.initialStepsize,'minStepsize',stimulus.minStepsize,'maxStepsize',stimulus.maxStepsize,'minThreshold',stimulus.minThreshold,'maxThreshold', stimulus.maxThreshold,'stepRule=pest','dispFig=1');


function stairNum = pickStair(task)
if task.thistrial.quadrant == 1 && task.thistrial.attend == 1; stairNum = 1;
elseif task.thistrial.quadrant == 1 && task.thistrial.attend == 2; stairNum = 2;
elseif task.thistrial.quadrant == 2 && task.thistrial.attend == 1; stairNum = 3;
elseif task.thistrial.quadrant == 2 && task.thistrial.attend == 2; stairNum = 4;
elseif task.thistrial.quadrant == 3 && task.thistrial.attend == 1; stairNum = 5;
elseif task.thistrial.quadrant == 3 && task.thistrial.attend == 2; stairNum = 6;
elseif task.thistrial.quadrant == 4 && task.thistrial.attend == 1; stairNum = 7;
elseif task.thistrial.quadrant == 4 && task.thistrial.attend == 2; stairNum = 8;end;




function doThreshold(task,stimulus,randThresh,testXvals,testYvals)
setGammaTableForMaxContrast(1);
comparison = [4 1 6 2 8 7 4 1]; %will cycle through in order; 1 is no difference
if length(comparison) ~= 8; sprintf('!!! comparison array should be 8 values !!!'), keyboard; end

%initiate gaussians for comparison
for index = 1:8;
contrastIndex(index) = getContrastIndex(task{1}{1}.parameter.contrast/(index));
end

%%%% random presentation %%%%
if randThresh == 1
mglClearScreen(stimulus.colors.grayColor); mglFlush;
keyboard
%stimuli here
mglClearScreen(stimulus.colors.grayColor); mglBltTexture(stimulus.tex(contrastIndex(1)),[testXvals(1) testYvals(1)-stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.mask,[testXvals(1) testYvals(1)-stimulus.grating.height/2]);mglBltTexture(stimulus.tex(contrastIndex(comparison(1))),[testXvals(1) testYvals(1)+stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.mask,[testXvals(1) testYvals(1)+stimulus.grating.height/2]); mglFlush;pause(.5);mglClearScreen(stimulus.colors.grayColor); mglFlush;keyboard
mglClearScreen(stimulus.colors.grayColor); mglBltTexture(stimulus.tex(contrastIndex(comparison(2))),[testXvals(2) testYvals(2)-stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.mask,[testXvals(2) testYvals(2)-stimulus.grating.height/2]);mglBltTexture(stimulus.tex(contrastIndex(1)),[testXvals(2) testYvals(2)+stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.mask,[testXvals(2) testYvals(2)+stimulus.grating.height/2]); mglFlush;pause(.5);mglClearScreen(stimulus.colors.grayColor); mglFlush;keyboard
mglClearScreen(stimulus.colors.grayColor); mglBltTexture(stimulus.tex(contrastIndex(1)),[testXvals(3) testYvals(3)-stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.mask,[testXvals(3) testYvals(3)-stimulus.grating.height/2]);mglBltTexture(stimulus.tex(contrastIndex(comparison(3))),[testXvals(3) testYvals(3)+stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.mask,[testXvals(3) testYvals(3)+stimulus.grating.height/2]); mglFlush;pause(.5);mglClearScreen(stimulus.colors.grayColor); mglFlush; keyboard
mglClearScreen(stimulus.colors.grayColor); mglBltTexture(stimulus.tex(contrastIndex(comparison(4))),[testXvals(4) testYvals(4)-stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.mask,[testXvals(4) testYvals(4)-stimulus.grating.height/2]);mglBltTexture(stimulus.tex(contrastIndex(1)),[testXvals(4) testYvals(4)+stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.mask,[testXvals(4) testYvals(4)+stimulus.grating.height/2]); mglFlush;pause(.5); mglClearScreen(stimulus.colors.grayColor); mglFlush;keyboard
mglClearScreen(stimulus.colors.grayColor); mglBltTexture(stimulus.tex(contrastIndex(1)),[testXvals(5) testYvals(5)-stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.mask,[testXvals(5) testYvals(5)-stimulus.grating.height/2]);mglBltTexture(stimulus.tex(contrastIndex(comparison(5))),[testXvals(5) testYvals(5)+stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.mask,[testXvals(5) testYvals(5)+stimulus.grating.height/2]); mglFlush;pause(.5);mglClearScreen(stimulus.colors.grayColor); mglFlush;keyboard
mglClearScreen(stimulus.colors.grayColor); mglBltTexture(stimulus.tex(contrastIndex(1)),[testXvals(6) testYvals(6)-stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.mask,[testXvals(6) testYvals(6)-stimulus.grating.height/2]);mglBltTexture(stimulus.tex(contrastIndex(comparison(6))),[testXvals(6) testYvals(6)+stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.mask,[testXvals(6) testYvals(6)+stimulus.grating.height/2]); mglFlush;pause(.5);mglClearScreen(stimulus.colors.grayColor); mglFlush;keyboard
mglClearScreen(stimulus.colors.grayColor); mglBltTexture(stimulus.tex(contrastIndex(1)),[testXvals(7) testYvals(7)-stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.mask,[testXvals(7) testYvals(7)-stimulus.grating.height/2]);mglBltTexture(stimulus.tex(contrastIndex(comparison(7))),[testXvals(7) testYvals(7)+stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.mask,[testXvals(7) testYvals(7)+stimulus.grating.height/2]); mglFlush;pause(.5);mglClearScreen(stimulus.colors.grayColor); mglFlush;keyboard
mglClearScreen(stimulus.colors.grayColor); mglBltTexture(stimulus.tex(contrastIndex(comparison(8))),[testXvals(8) testYvals(8)-stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.mask,[testXvals(8) testYvals(8)-stimulus.grating.height/2]);mglBltTexture(stimulus.tex(contrastIndex(1)),[testXvals(8) testYvals(8)+stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.mask,[testXvals(8) testYvals(8)+stimulus.grating.height/2]); mglFlush;pause(.5);mglClearScreen(stimulus.colors.grayColor); mglFlush;keyboard
end

%%%% gradient presentation %%%%%
if gradientThresh==1
gradientX=[-18:(36/7):18] %x and y coordinates of the the gradient comparisons
gradientY = -10
mglClearScreen(stimulus.colors.grayColor);
mglBltTexture(stimulus.tex(contrastIndex(1)),[gradientX(1) gradientY-stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.tex(contrastIndex(comparison(1))),[gradientX(1) gradientY+stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.mask,[gradientX(1) gradientY-stimulus.grating.height/2]); mglBltTexture(stimulus.mask,[gradientX(1) gradientY+stimulus.grating.height/2]);
mglBltTexture(stimulus.tex(contrastIndex(1)),[gradientX(2) gradientY-stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.tex(contrastIndex(comparison(2))),[gradientX(2) gradientY+stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.mask,[gradientX(2) gradientY-stimulus.grating.height/2]); mglBltTexture(stimulus.mask,[gradientX(2) gradientY+stimulus.grating.height/2]);
mglBltTexture(stimulus.tex(contrastIndex(1)),[gradientX(3) gradientY-stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.tex(contrastIndex(comparison(3))),[gradientX(3) gradientY+stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.mask,[gradientX(3) gradientY-stimulus.grating.height/2]); mglBltTexture(stimulus.mask,[gradientX(3) gradientY+stimulus.grating.height/2]);
mglBltTexture(stimulus.tex(contrastIndex(1)),[gradientX(4) gradientY-stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.tex(contrastIndex(comparison(4))),[gradientX(4) gradientY+stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.mask,[gradientX(4) gradientY-stimulus.grating.height/2]); mglBltTexture(stimulus.mask,[gradientX(4) gradientY+stimulus.grating.height/2]);
mglBltTexture(stimulus.tex(contrastIndex(1)),[gradientX(5) gradientY-stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.tex(contrastIndex(comparison(5))),[gradientX(5) gradientY+stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.mask,[gradientX(5) gradientY-stimulus.grating.height/2]); mglBltTexture(stimulus.mask,[gradientX(5) gradientY+stimulus.grating.height/2]);
mglBltTexture(stimulus.tex(contrastIndex(1)),[gradientX(6) gradientY-stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.tex(contrastIndex(comparison(6))),[gradientX(6) gradientY+stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.mask,[gradientX(6) gradientY-stimulus.grating.height/2]); mglBltTexture(stimulus.mask,[gradientX(6) gradientY+stimulus.grating.height/2]);
mglBltTexture(stimulus.tex(contrastIndex(1)),[gradientX(7) gradientY-stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.tex(contrastIndex(comparison(7))),[gradientX(7) gradientY+stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.mask,[gradientX(7) gradientY-stimulus.grating.height/2]); mglBltTexture(stimulus.mask,[gradientX(7) gradientY+stimulus.grating.height/2]);
mglBltTexture(stimulus.tex(contrastIndex(1)),[gradientX(8) gradientY-stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.tex(contrastIndex(comparison(8))),[gradientX(8) gradientY+stimulus.grating.height/2 stimulus.grating.height]);mglBltTexture(stimulus.mask,[gradientX(8) gradientY-stimulus.grating.height/2]); mglBltTexture(stimulus.mask,[gradientX(8) gradientY+stimulus.grating.height/2]);
mglFlush
keyboard
end