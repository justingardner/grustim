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
function westheimer(varargin)

getArgs(varargin,{'sawtoothProfile=increment'});

% initalize stimulus variable with settings
global stimulus;
stimulus.westheimer = init('sawtoothProfile',sawtoothProfile);

% open mgl screen
initScreen;
mglClearScreen(stimulus.westheimer.backgroundLuminance);mglFlush;

% loop for numSecs, displaying westheime stimulus.
numSecs = 5;
probeLuminancePlot = [];
for iFrame = 0:stimulus.westheimer.frameRate*numSecs
  % clear screen to background luminance
  mglClearScreen(stimulus.westheimer.backgroundLuminance);
  % update stimulus
  probeLuminance = update(iFrame,stimulus.westheimer);
  % flush
  mglFlush;
  % keep luminance to plot what we are doing
  probeLuminancePlot(end+1) = probeLuminance;
end
mglClearScreen(stimulus.westheimer.backgroundLuminance);mglFlush;

mlrSmartfig('westheimerSawtooth','reuse');
clf;
plot(probeLuminancePlot);
title(sprintf('Sawtooth freq %0.1f (Monitor frame rate: %i Hz)',stimulus.westheimer.sawtoothFrequency,stimulus.westheimer.frameRate))
xlabel('Frame Number');
ylabel('Luminance Value (normalized luminance)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init westheimer stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = init(varargin)

% init settings variable
s = [];

% get arguments
getArgs(varargin,{'backgroundLuminance=0.25','pedestalLuminance=0.5','probeWeberContrast=0.2','sawtoothFrequency=3','elementRadius=0.3','pedestalProbeSizeRatio=5','frameRate=60','elementGridSpacing',2/sqrt(3),'sawtoothProfile=increment'});

% Luminance values (normalized monitor units)
s.backgroundLuminance = backgroundLuminance;
s.pedestalLuminance = pedestalLuminance;
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
disp(sprintf('(westheimer:init) backgroundLuminance=%0.2f pedestalLuminance=%0.2f',s.backgroundLuminance,s.pedestalLuminance));
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
s.incrementProbeLuminance = s.pedestalLuminance + s.probeWeberContrast * s.pedestalLuminance;
s.decrementProbeLuminance = s.pedestalLuminance - s.probeWeberContrast * s.pedestalLuminance;
% derived times
s.framesPerCycle = s.frameRate / s.sawtoothFrequency;

% Loop through rows and columns to plot hexagons
numRows = 6;
numCols = 6;
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
s.nPedestals = size(s.pedestalHexagonX,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display westheimer stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function probeLuminance = update(frameNum,s)

% calculate probeLuminance
sawtoothValue = ((s.framesPerCycle - mod(frameNum,s.framesPerCycle))/s.framesPerCycle);

if s.increment
  % increments
  probeLuminance = (s.incrementProbeLuminance - s.pedestalLuminance)*sawtoothValue + s.pedestalLuminance;
else
  % decrements
  probeLuminance = (s.decrementProbeLuminance - s.pedestalLuminance)*sawtoothValue + s.pedestalLuminance;
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