
%
%      usage: westheimer
%         by: justin gardner
%       date: 11/09/2023
%    purpose: display stimuli for Westheimer sensitatization paradigm as 
%             created in:
%
%             Norcia AM, Yakovleva A, Jehangir N, Goldberg JL (2022) Preferential Loss of Contrast
%             Decrement Responses in Human Glaucoma. Investig Ophthalmol Vis Sci 63:16.
%
%             Norcia AM, Yakovleva A, Hung B, Goldberg JL (2020) Dynamics of Contrast Decrement
%             and Increment Responses in Human Visual Cortex. Transl Vis Sci Technol 9:6.
%
%
function myscreen = westheimerTask(varargin)
getArgs(varargin,{'block=[]'});

if isempty(block)
    warning('Missing input: set "block" to "on", "off", or "onoff"');
    return
end

global stimulus
stimulus.block = block;
% initalize stimulus variable with settings
if strcmp(block, 'on') || strcmp(block, 'onoff')
    stimulus.increment = init('sawtoothProfile','increment');
    stimulus.backgroundLuminance = stimulus.increment.backgroundLuminance;
end
if strcmp(block, 'off') || strcmp(block, 'onoff')
    stimulus.decrement = init('sawtoothProfile','decrement');
    stimulus.backgroundLuminance = stimulus.decrement.backgroundLuminance;
end

%displayName = 'test_small_screen';
%displayName = '';
%myscreen = initScreen(displayName);
myscreen = initScreen;
% mglClearScreen(stimulus.backgroundLuminance);

% set the first task to be the fixation staircase task 
global fixStimulus
fixStimulus.diskSize = 0.0;
fixStimulus.fixWidth = 1;
fixStimulus.fixLineWidth = 0.3;
fixStimulus.threshold = 1;
[task{1} myscreen] = fixStairInitTaskLR(myscreen);

% if strcmp(block, 'on') || strcmp(block, 'onoff')
%     update(0,stimulus.increment,1);
% elseif strcmp(block, 'off')
%     update(0,stimulus.decrement,1);
% end
% mglFlush;

% set our task to have two phases. 
% one starts out with dots moving for incohrently for 10 seconds
task{2}{2}.waitForBacktick = 0;
task{2}{2}.seglen = [12 11.5];
task{2}{2}.numBlocks = 11;
task{2}{2}.synchToVol = [0 1];

task{2}{1}.waitForBacktick = 0;
task{2}{1}.seglen = 0.5;
task{2}{1}.numBlocks = 1;
task{2}{1}.numTrials = 1;
task{2}{1}.synchToVol = 1;

% initialize our task
% for phaseNum = 1:length(task{2})
%   [task{2}{phaseNum} myscreen] = initTask(task{2}{phaseNum},myscreen,@startSegmentCallback,@updateScreenCallback);
% end
[task{2}{2}, myscreen] = initTask(task{2}{2},myscreen,@startSegmentCallback,@updateScreenCallback);
[task{2}{1}, myscreen] = initTask(task{2}{1},myscreen,@startSegmentCallback,@updateScreenCallbackPre);

% init the stimulus
% global stimulus;
myscreen = initStimulus('stimulus',myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task{2})) && ~myscreen.userHitEsc
  % update the dots
  [task{2} myscreen phaseNum] = updateTask(task{2},myscreen,phaseNum);
  % update the fixation task
  [task{1} myscreen] = updateTask(task{1},myscreen,1);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);
clear global stimulus;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;

if (task.thistrial.thisseg == 1)
  stimulus.f = 0;
    if strcmp(stimulus.block, 'on') || strcmp(stimulus.block, 'onoff')
        stimulus.show = 'on';
    elseif strcmp(stimulus.block, 'off')
        stimulus.show = 'off';
    end
else
    if strcmp(stimulus.block, 'onoff')
        stimulus.show = 'off';
    else % on or off block then static
        stimulus.show = 'static';
    end
    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task, myscreen)

global stimulus
mglClearScreen(stimulus.backgroundLuminance);
if strcmp(stimulus.show, 'on')
    update(stimulus.f,stimulus.increment,1);
elseif strcmp(stimulus.show, 'off')
    update(stimulus.f,stimulus.decrement,1);
    
elseif strcmp(stimulus.show, 'static')
  if strcmp(stimulus.block, 'on')
    update(stimulus.f,stimulus.increment,0);
  else
    update(stimulus.f,stimulus.decrement,0);
  end
end

stimulus.f = stimulus.f + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallbackPre(task, myscreen)

global stimulus
mglClearScreen(stimulus.backgroundLuminance);
if strcmp(stimulus.show, 'on')
    update(stimulus.f,stimulus.increment,1);
elseif strcmp(stimulus.show, 'off')
    update(stimulus.f,stimulus.decrement,1);
    
elseif strcmp(stimulus.show, 'static')
  if strcmp(stimulus.block, 'on')
    update(stimulus.f,stimulus.increment,0);
  else
    update(stimulus.f,stimulus.decrement,0);
  end
end

% stimulus.f = stimulus.f + 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init westheimer stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = init(varargin)

% init settings variable
s = [];

s.dispContrastMin = 0.1;
s.dispContrastMax = 1;

% get arguments
getArgs(varargin,{'backgroundLuminance=0.0','pedestalLuminance=0.5','probeWeberContrast=0.4','sawtoothFrequency=3','elementRadius=0.25','pedestalProbeSizeRatio=5','frameRate=60','elementGridSpacing',2/sqrt(3),'sawtoothProfile=increment'});

% Luminance values (normalized monitor units)
s.backgroundLuminanceOrig = backgroundLuminance;
s.pedestalLuminanceOrig = pedestalLuminance;
% Converted Luminance values (using limited contrast range)
s.backgroundLuminance = (s.dispContrastMax - s.dispContrastMin) * s.backgroundLuminanceOrig + s.dispContrastMin;
s.pedestalLuminance = (s.dispContrastMax - s.dispContrastMin) * s.pedestalLuminanceOrig + s.dispContrastMin;
% modulation depth of sawtooth, as Weber contrast 
s.probeWeberContrast = probeWeberContrast;
% frequency for sawtooth modulation in hertex
s.sawtoothFrequency = sawtoothFrequency; 
s.frameRate = frameRate;
% radius of the base element in degrees (this will later get M-scaled)
s.elementRadius = elementRadius; 
% ratio of size of pedestal to probe
s.pedestalProbeSizeRatio = pedestalProbeSizeRatio;
% grid spacing of elements (as a ratio of the element size)
s.elementGridSpacing = elementGridSpacing;

% display
dispHeader('westheimer stimulus settings')
disp(sprintf('(westheimer:init) backgroundLuminance=%0.2f pedestalLuminance=%0.2f',s.backgroundLuminanceOrig,s.pedestalLuminanceOrig));
disp(sprintf('(westheimer:init) backgroundLuminanceCorrected=%0.2f pedestalLuminanceCorrected=%0.2f',s.backgroundLuminance, s.pedestalLuminance));
disp(sprintf('(westheimer:init) probeWeberContrast=%0.2f',s.probeWeberContrast));
disp(sprintf('(westheimer:init) sawtoothFrequency=%0.2f Hz (frameRate=%i Hz)',s.sawtoothFrequency,s.frameRate));
disp(sprintf('(westheimer:init) elementRadius=%0.2f pedestalProbeSizeRatio=%0.2f elementGridSpacing=%f',s.elementRadius,s.pedestalProbeSizeRatio,s.elementGridSpacing));

if strcmp(lower(sawtoothProfile),'increment')
  s.increment = true;
  disp(sprintf('(westheimer:init) sawtoothProfile=increment'));
else
  s.increment = false;
  disp(sprintf('(westheimer:init) sawtoothProfile=decrement'));
end

% derived sizes
s.probeRadius = s.elementRadius / s.pedestalProbeSizeRatio;
s.elementWidth = sqrt(3)*s.elementRadius;
s.elementHeight = 3*s.elementRadius/2;
% derived luminances
s.incrementProbeLuminanceOrig = s.pedestalLuminanceOrig + s.probeWeberContrast * s.pedestalLuminanceOrig;
s.decrementProbeLuminanceOrig = s.pedestalLuminanceOrig - s.probeWeberContrast * s.pedestalLuminanceOrig;
% converted luminances
s.incrementProbeLuminance = (s.dispContrastMax - s.dispContrastMin) * s.incrementProbeLuminanceOrig + s.dispContrastMin;
s.decrementProbeLuminance = (s.dispContrastMax - s.dispContrastMin) * s.decrementProbeLuminanceOrig + s.dispContrastMin;
% derived times
s.framesPerCycle = s.frameRate / s.sawtoothFrequency;

% Loop through rows and columns to plot hexagons
numRows = 8;
numCols = 8;
s.pedestalHexagonX = [];
s.pedestalHexagonY = [];
s.probeHexagonX = [];
s.probeHexagonY = [];
for yPos = -numRows:numRows
  % offset odd rows
  if isodd(yPos)
    % offset
    xOffset = (s.elementGridSpacing*s.elementWidth)/2;
    % add an extra starting element
    xStart = -numCols-1;
  else
    % offset
    xOffset = 0;
    % add an extra starting element
    xStart = -numCols;
  end

  for xPos = xStart:numCols
    % Calculate the center of the hexagon
    xCenter = xPos * s.elementWidth * s.elementGridSpacing + xOffset;
    yCenter = yPos * s.elementHeight * s.elementGridSpacing;
        
    if ~((xPos == 0) && (yPos == 0))
      % get pedestal hexagon
      [xPolygon yPolygon] = getHexagon(xCenter, yCenter, s.elementRadius);
      % split into two quads - this is just so that we can use the
      % mglQuad function which allows sending multiple quads and works
      % faster than mglPolygon
      s.pedestalHexagonX(1:4,end+1:end+2) = [xPolygon(1:4)' xPolygon(4:7)'];
      s.pedestalHexagonY(1:4,end+1:end+2) = [yPolygon(1:4)' yPolygon(4:7)'];

      % get probe hexagon
      [xPolygon yPolygon] = getHexagon(xCenter, yCenter, s.probeRadius);
      s.probeHexagonX(1:4,end+1:end+2) = [xPolygon(1:4)' xPolygon(4:7)'];
      s.probeHexagonY(1:4,end+1:end+2) = [yPolygon(1:4)' yPolygon(4:7)'];
    end
  end
end
s.nPedestals = size(s.pedestalHexagonX,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display westheimer stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function probeLuminance = update(frameNum,s,dynamic)

% calculate probeLuminance
sawtoothValue = ((s.framesPerCycle - mod(frameNum,s.framesPerCycle))/s.framesPerCycle);

if ~dynamic
    probeLuminance = s.pedestalLuminance;
    
else % dynamic

    if s.increment
        % increments
%         probeLuminance = (s.incrementProbeLuminance - s.pedestalLuminance)*sawtoothValue + s.pedestalLuminance;
        probeLuminanceOrig = (s.incrementProbeLuminanceOrig - s.pedestalLuminanceOrig)*sawtoothValue + s.pedestalLuminanceOrig;
        probeLuminance = (s.dispContrastMax - s.dispContrastMin) * probeLuminanceOrig + s.dispContrastMin;

    else
            % decrements
%         probeLuminance = (s.decrementProbeLuminance - s.pedestalLuminance)*sawtoothValue + s.pedestalLuminance;
        probeLuminanceOrig = (s.decrementProbeLuminanceOrig - s.pedestalLuminanceOrig)*sawtoothValue + s.pedestalLuminanceOrig;
        probeLuminance = (s.dispContrastMax - s.dispContrastMin) * probeLuminanceOrig + s.dispContrastMin;
    end
end
% display the pedestals
mglQuad(s.pedestalHexagonX(:,:),s.pedestalHexagonY(:,:),repmat(s.pedestalLuminance,1,s.nPedestals),1);

% display the probes
mglQuad(s.probeHexagonX(:,:),s.probeHexagonY(:,:),repmat(probeLuminance,1,s.nPedestals),1);

%%%%%%%%%%%%%%%%%%
% get hexagon coordinates
%%%%%%%%%%%%%%%%%%
function [x y] = getHexagon(xCenter,yCenter,elementRadius)

% start with angles of each vertex
alpha = pi/6:pi/3:2*pi+pi/6;

% get coordinates of hexagonal centers
x = elementRadius * cos(alpha) + xCenter;
y = elementRadius * sin(alpha) + yCenter;

% convert to polar coordinates
[theta, rho] = cart2pol(x,y);

% apply cortical mangiciation to eccentricity
% where k is cortical mangification factor
k = 0.7;
rho = rho .* (1 + (rho/k));

% convert back to cartesian coordinates
[x,y] = pol2cart(theta, rho);
