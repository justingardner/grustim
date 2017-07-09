% gruProjStencil.m
%
%        $Id:$ 
%      usage: gruProjStencil(myscreen)
%         by: Dan Birman
%       date: 08/13/15
%    purpose: Generate a mask stencil for the GE projecter @ Stanford
%             Builds a stencil based on the hardcoded list of 
%             eccentricities (mesEcc) that could be seen at the far
%             extremes of each polar angle (mesAngles) in the magnet
%             with the projection screen.
%
%             This is the screenMaskFunction we use in mglEditScreenParams
%             and gets called by initScreen
%             
%             It adjusts to compensate for flipping, shifting or scaling
%             the screen.
%
function gruProjStencil(myscreen)

% check input arguments
if nargin ~= 1
  help gruProjStencil
  return
end

% polar angles where eccentricity was measured
mesAngles =  [-15 0 15 30 45 60 75 90];
% eccentricity measured at the above angles
mesEcc = ([38 35 30 26 24 22 21 21]-1)*1;

% check match of eccentricity and angles
if length(mesEcc) ~= length(mesAngles)
  error('Eccentricity array has incorrect length.');
end

% linear interpolate (in radial coordinates) to make smoother
angles = -15:90;
ecc = interp1(mesAngles,mesEcc,angles,'linear');

% width and height of screen in degrees when measured
mesWidth = 74.9675043592654;
mesHeight = 44.4029738290369;

% get current size of display (note that we do not
% use imageWidth and imageHeight in myscreen since these values
% may be croppped/scaled, so we get directly from mglGetParam)
curWidth = mglGetParam('deviceWidth');
curHeight = mglGetParam('deviceHeight');

% convert to cartesian and adjust for difference in screen size
x = (curWidth/mesWidth) * cos(pi*angles/180).*ecc;
y = (curHeight/mesHeight) * sin(pi*angles/180).*ecc;

% now flip it left/right to get the other side and add bottom corners
% note that the curvature is just at the top so this uses the
% bottom corneres to make the bottom square

% dan's fix: adjusts the corners based on the relative size, so that the
% stencil displays outside the screen if necessary
lx = -curWidth/2*(curWidth/mesWidth);
rx = curWidth/2*(curWidth/mesWidth);
ly = -curHeight/2*(curHeight/mesHeight);
ry = -curHeight/2*(curHeight/mesHeight);
x = [lx -x fliplr(x) rx];
y = [ly y fliplr(y) ry];
x = [-x fliplr(x)];
y = [y fliplr(y)];

% account for origin shift
x = x-myscreen.shiftOrigin(1);
y = y-myscreen.shiftOrigin(2);

% flip if necessary
if myscreen.flipHV(1),x = -x;end
if myscreen.flipHV(2),y = -y;end

% draw the shape in white - shifting to offset any shiftOrigin
% so that the stencil stays locked to screen even if we shift
% the origin to deal with subject issues
mglPolygon(x,y,1);



