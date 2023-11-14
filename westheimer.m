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

getArgs(varargin,{'sawtoothProfile=increment'})

% Luminance values (normalized monitor units)
backgroundLuminance = 0.25;
pedestalLuminance = 0.5;
% modulation depth of sawtooth, as Weber contrast 
probeWeberContrast = 0.2;
% frequency for sawtooth modulation in hertex
sawtoothFrequency = 3; 
frameRate = 60;
% radius of the base element in degrees (this will later get M-scaled)
elementRadius = 0.8; 
% ratio of size of pedestal to probe
pedestalProbeSizeRatio = 5;
% grid spacing of elements (as a ratio of the element size)
elementGridSpacing = 2/sqrt(3);

% derived sizes
probeRadius = elementRadius / pedestalProbeSizeRatio;
elementWidth = sqrt(3)*elementRadius;
elementHeight = 3*elementRadius/2;
% derived luminances
incrementProbeLuminance = pedestalLuminance + probeWeberContrast * pedestalLuminance;
decrementProbeLuminance = pedestalLuminance - probeWeberContrast * pedestalLuminance;
% derived times
framesPerCycle = frameRate / sawtoothFrequency;

% Loop through rows and columns to plot hexagons
numRows = 2;
numCols = 2;
pedestalHexagon = zeros(numCols*2+1,numRows*2+1,2,7);
probeHexagon = zeros(numCols*2+1,numRows*2+1,2,7);
for yPos = -numRows:numRows
    for xPos = -numCols:numCols
        % Calculate the center of the hexagon
        xCenter = xPos * elementWidth * elementGridSpacing;
        yCenter = yPos * elementHeight * elementGridSpacing;
        
        % offset every other row
        if isodd(yPos)
          xCenter = xCenter + (elementGridSpacing*elementWidth)/2;
        end

        % get index values
        xIndex = xPos + numCols + 1;
        yIndex = yPos + numRows + 1;
        
        % get pedestal hexagon
        [pedestalHexagon(xIndex,yIndex,1,:),pedestalHexagon(xIndex,yIndex,2,:)] = getHexagon(xCenter, yCenter, elementRadius);

        % get probe hexagon
        [probeHexagon(xIndex,yIndex,1,:),probeHexagon(xIndex,yIndex,2,:)] = getHexagon(xCenter, yCenter, probeRadius);
    end
end

mglOpen;
mglVisualAngleCoordinates(57,[40 30]);
mglClearScreen(backgroundLuminance);
mglFlush;

numSecs = 5;
probeLuminancePlot = [];
for iFrame = 0:frameRate*numSecs
  % clera screen to background luminance
  mglClearScreen(backgroundLuminance);
  % calculate probeLuminance
  sawtoothValue = ((framesPerCycle - mod(iFrame,framesPerCycle))/framesPerCycle);
  % increments
  probeLuminance = (incrementProbeLuminance - pedestalLuminance)*sawtoothValue + pedestalLuminance;
  % decrements
  %probeLuminance = (decrementProbeLuminance - pedestalLuminance)*sawtoothValue + pedestalLuminance;
  % cycle over all hexagons
  for yPos = -numRows:numRows
    for xPos = -numCols:numCols
      % get index values
      xIndex = xPos + numCols + 1;
      yIndex = yPos + numRows + 1;
      % draw pedestal hexagon
      mglPolygon(pedestalHexagon(xIndex,yIndex,1,:),pedestalHexagon(xIndex,yIndex,2,:),pedestalLuminance);
      % draw probe hexagon
      mglPolygon(probeHexagon(xIndex,yIndex,1,:),probeHexagon(xIndex,yIndex,2,:),probeLuminance);
    end
  end
  mglFlush;
  probeLuminancePlot(end+1) = probeLuminance;
end

mlrSmartfig('westheimerSawtooth','reuse');
clf;
plot(probeLuminancePlot);
title(sprintf('Sawtooth freq %0.1f (Monitor frame rate: %i Hz)',sawtoothFrequency,frameRate))
xlabel('Frame Number');
ylabel('Luminance Value (normalized luminance)')

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