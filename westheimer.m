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
elementRadius = 0.3; 
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
numRows = 6;
numCols = 6;
pedestalHexagonX = [];
pedestalHexagonY = [];
probeHexagonX = [];
probeHexagonY = [];
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
        [xPolygon yPolygon] = getHexagon(xCenter, yCenter, elementRadius);
        % split into two quads - this is just so that we can use the
        % mglQuad function which allows sending multiple quads and works
        % faster than mglPolygon
        pedestalHexagonX(1:4,end+1:end+2) = [xPolygon(1:4)' xPolygon(4:7)'];
        pedestalHexagonY(1:4,end+1:end+2) = [yPolygon(1:4)' yPolygon(4:7)'];

        % get probe hexagon
        [xPolygon yPolygon] = getHexagon(xCenter, yCenter, probeRadius);
        probeHexagonX(1:4,end+1:end+2) = [xPolygon(1:4)' xPolygon(4:7)'];
        probeHexagonY(1:4,end+1:end+2) = [yPolygon(1:4)' yPolygon(4:7)'];        
    end
end
nPedestals = size(pedestalHexagonX,2);

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
   % display
  mglQuad(pedestalHexagonX(:,:),pedestalHexagonY(:,:),repmat(pedestalLuminance,1,nPedestals),1);
  mglQuad(probeHexagonX(:,:),probeHexagonY(:,:),repmat(probeLuminance,1,nPedestals),1);
  % flush
  mglFlush;
  % keep luminance to plot what we are doing
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