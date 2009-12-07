% drawWedge - draw a wedge with alternating black and white pattern
%
%      usage: [  ] = drawWedge(startAngle, sweepAngle, nSlices, nRings, segPhase)
%         by: denis schluppeck
%       date: 2007-03-20
%        $Id: drawWedge.m,v 1.1 2007/03/20 20:04:44 ds Exp $:
%     inputs: 
%    outputs: 
%
%    purpose: 
%
%        e.g:
%
% function [ ] =drawWedge(startAngle, sweepAngle, nSlices, nRings, segPhase)

rmin = 1;
rmax = 10;
nRings = 4;
nSlices = 3;
nElements = nSlices * nRings;
r = linspace(rmin, rmax, nRings+1);
isize = repmat( r(1:end-1),nSlices,1);
osize = repmat( r(2:end),nSlices,1);
x = zeros(1,nElements);
y = zeros(1,nElements);

theta0 = 0;
thetaSweep = 45;

theta = linspace(theta0, thetaSweep, nSlices+1)
startAngles = repmat(theta(1:end-1)', 1, nRings);
sweepAngles = repmat(diff(theta)', 1, nRings);;



segPhase = 0;
colors = [1 1 1; 0 0 0]';
allcolors = [];
for iRing = 1:nRings
  colors = fliplr(colors); % alternate for each slice
  allcolors = [allcolors repmat(colors, 1, nSlices)]
end
allcolors = allcolors(:, 1:nElements);

mglOpen
mglVisualAngleCoordinates(57,[40 30]);
mglGluPartialDisk(x, y, isize, osize, startAngles, sweepAngles, allcolors, 60, 2);
mglFlush
mglWaitSecs(0.2);
mglGluPartialDisk(x, y, isize, osize, startAngles, sweepAngles, fliplr(allcolors), 60, 2);
mglFlush

