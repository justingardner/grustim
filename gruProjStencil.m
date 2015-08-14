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
mesEcc = [25 21 17 15 13 12 11 11]*32/21;

% check match of eccentricity and angles
if length(mesEcc) ~= length(mesAngles)
  error('Eccentricity array has incorrect length.');
end

% linear interpolate (in radial coordinates) to make smoother
angles = -15:90;
ecc = interp1(mesAngles,mesEcc,angles,'linear');

% width and height of screen in degrees when measured
mesWidth = 51.5;
mesHeight = 29.8;

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
x = [-curWidth/2 -x fliplr(x) curWidth/2];
y = [-curHeight/2 y fliplr(y) -curHeight/2];

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



