function gruProjStencil(myscreen)
%MGLPROJSTENCIL Generate a mask stencil for the GE projecter @ Stanford
%   Builds a stencil based on the input list, which corresponds to the
%   visual eccentricity at angles [90, -90, 0, 15, 30, 45, 60, 75, -15, -30, -45, -60, -75]
% use by adding mglStencilSelect(1); and then mglStencilSelect(0); to turn
% off

% This approximates the full screen (needs to be adjusted):

% botangs1 = [-15 -30];
% botangs2 = [-45 -60 -75];
% transEccs = [10./cos(deg2rad(abs(botangs1))) 10./cos(deg2rad(90+botangs2))];
% eccArray = [15,15,12.5,27,24,21,18,15, transEccs];
% eccArray = eccArray;

% This is the actual recordings we made for what you can see in the 32
% channel coil:

% eccArray = [5.5, 7.5, 17, 18, 19, 11, 8, 6,14 , 15, 11, 8, 7.5];

% polar angles where eccentricity was measured
angles = [90, -90, 0, 15, 30, 45, 60, 75, -15, -30, -45, -60, -75];

% width and height of screen in degrees when measured
mesWidth = 51.5;
mesHeight = 29.8;

% eccentricity measured at the above angles
eccArray = [11, 14, 21, 17, 15, 13, 12, 11, 25, 23, 16,15, 14]*32/21;

% check match of eccArray and angles
if length(eccArray) ~= length(angles)
  error('Eccentricity array has incorrect length.');
end

% get current size of display (note that we do not
% use imageWidth and imageHeight in myscreen since these values
% may be croppped/scaled, so we get directly from mglGetParam)
curWidth = mglGetParam('deviceWidth');
curHeight = mglGetParam('deviceHeight');

% rescale eccArray to current display size
for iAngle = 1:length(angles)
  % get the angle
  thisAngle = pi*angles(iAngle)/180;
  % get ratio of original measurement to this screen
  mesLen = sqrt((cos(thisAngle)*mesWidth)^2 + (sin(thisAngle)*mesHeight)^2);
  curLen = sqrt((cos(thisAngle)*curWidth)^2 + (sin(thisAngle)*curHeight)^2);
  % scale the eccArray accordingly
  eccArray(iAngle) = eccArray(iAngle)*curLen/mesLen;
end

% make the other side as a symmetric copy of the side measured
angles2 = [180-angles(3:8) angles(9:end)-90];
angles = [angles angles2];
angles = mod(angles-90,360);
eccArray = [eccArray eccArray(3:8) fliplr(eccArray(9:end))];

% we will draw partial disks at every 3 deg angle, and using the linear
% interpolated position based on the eccArray

xs = [];
ys = [];
eccs = [];
sAs = [];

for curBlock = 0:15:345
  for withinB = 0:0.5:14.5
    % all gluPartial disks will start in the center of the screen
    % note that we use the shiftOrigin values of myscreen
    % to make sure that we keep the stencil in the same location
    % even if we move the origin
    xs = [xs -myscreen.shiftOrigin(1)];
    ys = [ys -myscreen.shiftOrigin(2)];
    % add on the current start angle to the list of start angles
    cSA = curBlock+withinB;
    sAs = [sAs cSA];
    % linearly interpolate the eccentricity array for this set of angles
    origEcc = eccArray(angles==curBlock);
    finalEcc = eccArray(angles==mod(curBlock+15,360));
    interpEcc = (finalEcc-origEcc) * withinB / 15 + origEcc;
    %         disp(sprintf('%i: %.03f',curBlock,interpEcc));
    % add add to the eccentricity list
    eccs = [eccs interpEcc];
  end
end

% draw the partial disks for one side
mglGluPartialDisk(xs,ys,zeros(size(xs)),eccs,sAs,ones(size(xs))*2,ones(3,size(xs,2)));
% and the other
mglGluPartialDisk(xs,ys,zeros(size(xs)),eccs,sAs,-ones(size(xs))*2,ones(3,size(xs,2)));




