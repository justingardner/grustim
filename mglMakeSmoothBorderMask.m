% mglMakeSmoothBorderMask.m
%
%      usage: mglMakeSmoothBorderMask(width,height,diameter,borderWidth,sigmaFromPeak,<xDeg2pix>,<yDeg2pix>)
%         by: Georgios Keliris (modified mglMakeGrating)
%       date: 12/27/12
%    purpose: create a radail 2D grating. You should start mgl
%             and use mglVisualAngleCoordinates before using.
%
%             width and height are in degrees of visual angle
%             sf is in cycles/degrees
%             phase is in degrees 
%
%             xDeg2pix and yDeg2pix are optional arguments that specify the
%             number of pixels per visual angle in the x and y dimension, respectively.
%             If not specified, these values are derived from the open mgl screen (make
%             sure you set mglVisualAngleCoordinates).
%       e.g.:
%
% mglOpen;
% mglVisualAngleCoordinates(57,[52 32.7]);
% g = mglMakeGrating(10,10,1.5,90,0);
% m = mglMakeSmoothBorderMask(10,10,8,1.5,3);
% g = g.*m;
% g = 255*(g+1)/2;
% tex = mglCreateTexture(g);
% mglClearScreen(0.5);
% mglBltTexture(tex,[0 0]);
% mglFlush;

function m = mglMakeSmoothBorderMask(width,height,diameter,borderWidth,sigmaFromPeak,xDeg2pix,yDeg2pix)

% check arguments
m = [];
if ~any(nargin == [3 4 5 6 7])
  help mglMakeSmoothBorderMask
  return
end

if ieNotDefined('borderWidth'),borderWidth = 0.5;end
if ieNotDefined('sigmaFromPeak'),sigmaFromPeak = 3;end

% defaults for xDeg2pix
if ieNotDefined('xDeg2pix')
  if isempty(mglGetParam('xDeviceToPixels'))
    disp(sprintf('(mglMakeSmoothBorderMask) mgl is not initialized'));
    return
  end
  xDeg2pix = mglGetParam('xDeviceToPixels');
end

% defaults for yDeg2pix
if ieNotDefined('yDeg2pix')
  if isempty(mglGetParam('yDeviceToPixels'))
    disp(sprintf('(mglMakeSmoothBorderMask) mgl is not initialized'));
    return
  end
  yDeg2pix = mglGetParam('yDeviceToPixels');
end

% get size in pixels
widthPixels = round(width*xDeg2pix);
heightPixels = round(height*yDeg2pix);
widthPixels = widthPixels + mod(widthPixels+1,2);
heightPixels = heightPixels + mod(heightPixels+1,2);

if diameter <= 0
    m=zeros(widthPixels,heightPixels);
    return
end

if 2*borderWidth > diameter
    disp(sprintf('(mglMakeSmootBorderMask) borderWidth cannot be bigger than half diameter (radius)'));
    return
end


% calculate function in 1D
max_R = round(sqrt((0.5*width)^2+(0.5*height)^2)+1);
resolution = 2*xDeg2pix;
s = zeros(1,floor(diameter/2*resolution)+1); 
xs = 0:1/2/xDeg2pix:diameter/2;
xs_ext = diameter/2+1:1/2/xDeg2pix:max_R;
sigma=borderWidth/sigmaFromPeak;
x0=diameter/2-borderWidth;
points=length(s)+1-[1:resolution*borderWidth];
s(points)=1/(sqrt(2*pi)*sigma) * exp(-(xs(points)-x0).^2./(2*sigma^2));
s(s==0)=max(s);
s=(s-min(s))/(max(s)-min(s));

% 2D
% get a grid of x and y coordinates that has
% the correct number of pixels
x = -width/2:width/(widthPixels-1):width/2;
y = -height/2:height/(heightPixels-1):height/2;
[xMesh,yMesh] = meshgrid(x,y);
[thMesh,rMesh] = cart2pol(xMesh,yMesh);

% compute mask
xs = [xs xs_ext]; s = [s zeros(size(xs_ext))];
m = interp1(xs,s,rMesh,'linear');