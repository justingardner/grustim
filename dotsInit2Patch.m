% dotsInit2Patch.m
%
%        $Id:$ 
%      usage: dotsInit2Patch (dotsInit)
%         by: austin kuo (justin gardner)
%       date: 11/15/19 (07/08/17)
%    purpose: Create and update a data structure for random dots
%             that can be used to display different kinds of dot
%             movement, coherence and contrast
%             
%             On init, you can set the following properties
%             xCenter: x center in degrees of dots patch
%             yCenter: y center in degrees of dots patch
%             width: width in degrees
%             speed: speed of dots
%             dir: direction of dots (in degrees for linear)
%             coherence: coherence of dots (from 0 to 1)
%             type: type of dots. linear (default) is 2D linear
%                   motion. opticFlow is an optic flow field (note 
%                   that direction for opticFlow goes from -1 to 1 (in to out)
%             dotSize: size in pixels of dots
%             density: dots per deg^2
%             contrast: contrast of dots
%             framesPerSecond: This defaults to 60, but should be set
%               if you want the speeds to be correct
%             mask: set to 1 for a circular mask
%             drawType: Can be plot or mgl (for drawing in a figure or 
%                using mgl to draw
%
%            Note that after init, you should not change the
%            structure parameters yourself (as their may be 
%            some calculated fields that need to also be changed)
%            Instead, use the callbacks to change the stimulus
%
%            dots = dots.setSpeed(dots, speed) - speeds is in deg/s
%            dots = dots.dir(dots, dir) - dir is in deg or -1 - 1 for flow
%            dots = dots.setContrast(dots,contrast) - contrast is 0 - 1
%            dots = dots.setCenter(dots,x,y) - sets the x,y center in deg
%
%            To draw the dots, you update (which changes their position and
%            needs to be called at framesPerSecond) and draw
%  
%            dots = dots.update(dots);
%            dots = dots.draw(dots);
%
% e.g.
%
% % init the dots (using plot to draw them - typically, mgl is default)
% dots = dotsInit('drawType=plot');
%
% % set direction, speed and contrast
% dots = dots.setSpeed(dots,4);
% dots = dots.setDir(dots,45);
% dots = dots.setContrast(dots,0.5);
%
% % dispplay 50 frames
% for i = 1:50
%   % update the dots
%   dots = dots.update(dots);
%
%   % draw the dots
%   dots = dots.draw(dots);
% end
%
%
function dots = dotsInit2Patch(varargin)

% parse arguments
contrast = 1;type='linear';dir = 1;
getArgs(varargin,{'xCenter=0','yCenter=0','width=5','height=[]','speed=1','dir=1','coherence=1','type=linear','dotSize=4','density=5','contrast=1','framesPerSecond=60','mask=1','drawType=mgl'});

% set parameters of dots
dots.width = width;
if ~isempty(height)
  % set the height
  dots.height = height;
else
  % if height is empty, then make a square patch
  dots.height = width;
end
dots.speed = speed;
dots.dir = dir;
dots.coherence = coherence;
dots.type = lower(type);
dots.dotSize = dotSize;
dots.density = density;
dots.contrast = contrast;
dots.framesPerSecond = framesPerSecond;
dots.setCenter = @setCenter;
dots = dots.setCenter(dots,xCenter,yCenter);
dots.setCoherence = @setCoherence;
dots = dots.setCoherence(dots,coherence); % probably unnecessary line

% mask
dots.mask = mask;
if dots.mask
  dots.applyMask = @applyMask;
else
  dots.applyMask = @noMask;
end
  
% dot type specific functions
if strcmp(dots.type,'linear')
  dots.setSpeed = @setSpeedLinear;
  dots.setDir = @setDirLinear;
  dots.update = @updateDotsLinear;
  % init the dots
  dots = initDotsLinear(dots);
elseif strcmp(dots.type,'opticflow')
  dots.setSpeed = @setSpeedOpticflow;
  dots.setDir = @setDirOpticflow;
  dots.update = @updateDotsOpticflow;
  % init the dots
  dots = initDotsOpticflow(dots);
else
  disp(sprintf('(dotsInit) Unknown dot type: %s',dots.type));
  keyboard
end

% set draw function
if strcmp(lower(drawType),'plot')
  dots.draw = @plotDotsDraw;
elseif strcmp(lower(drawType),'mgl')
  dots.draw = @mglDotsDraw;
else
  disp(sprintf('(dotsInit) Unknown draw type: %s',dots.drawType));
  keyboard
end

% set color of each dot
dots.blackOrWhite  = 2*(rand(1,dots.n) > 0.5)-1;
dots.black = dots.blackOrWhite == -1;
dots.white = dots.blackOrWhite == 1;

% set contrast
dots.setContrast = @setContrast;
dots = dots.setContrast(dots,dots.contrast);

% update them once (to make sure all fields get filled in)
dots = dots.update(dots);

%%%%%%%%%%%%%%%%%%%
%    setCenter    %
%%%%%%%%%%%%%%%%%%%
function dots = setCenter(dots,xCenter,yCenter)

% set x,y center of dots
dots.xCenter = xCenter;
dots.yCenter = yCenter;

%%%%%%%%%%%%%%%%%%%%%
%    setContrast    %
%%%%%%%%%%%%%%%%%%%%%
function dots = setContrast(dots,contrast)

% set the contrast
dots.contrast = contrast;
dots.blackColor = repmat(0.5-dots.contrast/2,1,3);
dots.whiteColor = repmat(0.5+dots.contrast/2,1,3);

%%%%%%%%%%%%%%%%%%%%%%
%    setCoherence    %
%%%%%%%%%%%%%%%%%%%%%%
function dots = setCoherence(dots,coherence)

% set the coherence
dots.coherence = coherence;

%%%%%%%%%%%%%%%%%%%%%
%    mglDotsDraw    %
%%%%%%%%%%%%%%%%%%%%%
function dots = mglDotsDraw(dots)

% get dots passed through mask
[blackX blackY whiteX whiteY] = dots.applyMask(dots);

% draw black dots
mglPoints2(blackX+dots.xCenter,blackY+dots.yCenter,dots.dotSize,dots.blackColor);
% draw white dots
mglPoints2(whiteX+dots.xCenter,whiteY+dots.yCenter,dots.dotSize,dots.whiteColor);

%%%%%%%%%%%%%%%%
%    noMask    %
%%%%%%%%%%%%%%%%
function [blackX blackY whiteX whiteY] = noMask(dots)

blackX = dots.x(dots.black);
blackY = dots.y(dots.black);
whiteX = dots.x(dots.white);
whiteY = dots.y(dots.white);

%%%%%%%%%%%%%%%%%%%
%    applyMask    %
%%%%%%%%%%%%%%%%%%%
function [blackX blackY whiteX whiteY] = applyMask(dots)

goodDots = (dots.x.^2 + dots.y.^2) <= (dots.width/2)^2;
blackX = dots.x(dots.black & goodDots);
blackY = dots.y(dots.black & goodDots);
whiteX = dots.x(dots.white & goodDots);
whiteY = dots.y(dots.white & goodDots);

%%%%%%%%%%%%%%%%%%%%%%
%    plotDotsDraw    %
%%%%%%%%%%%%%%%%%%%%%%
function dots = plotDotsDraw(dots)

% get dots passed through mask
[blackX blackY whiteX whiteY] = dots.applyMask(dots);

% clear figure
cla;
% set background to gray
set(gca,'Color',[0.5 0.5 0.5]);

% draw black dots
plot(blackX,blackY,'.','Color',dots.blackColor);hold on

% draw white dots
plot(whiteX,whiteY,'.','Color',dots.whiteColor);

% set axis
axis([dots.xCenter-dots.width/2 dots.xCenter+dots.width/2 dots.yCenter-dots.height/2 dots.yCenter+dots.height/2]);
axis square

% and draw
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the dots speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = setSpeedLinear(dots,speed)

% get the step size
dots.speed = speed;
dots.stepsize = speed/dots.framesPerSecond;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the dots direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = setDirLinear(dots,direction)

% get the direction
dots.dir = direction;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots for linear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updateDotsLinear(dots,coherence)

% get the dots step
dots.xstep = cos(pi*dots.dir/180)*dots.stepsize;
dots.ystep = sin(pi*dots.dir/180)*dots.stepsize;

% pick a random set of dots
dots.coherent = rand(1,dots.n) < dots.coherence;

% now move those dots in the right direction
dots.x(dots.coherent) = dots.x(dots.coherent)+dots.xstep;
dots.y(dots.coherent) = dots.y(dots.coherent)+dots.ystep;

% randomwalk rule
thisdir = rand(1,sum(~dots.coherent))*2*pi;
dots.x(~dots.coherent) = dots.x(~dots.coherent)+cos(thisdir)*dots.stepsize;
dots.y(~dots.coherent) = dots.y(~dots.coherent)+sin(thisdir)*dots.stepsize;

% movshon noise
%dots.x(~dots.coherent) = rand(1,sum(~dots.coherent))*dots.width;
%dots.y(~dots.coherent) = rand(1,sum(~dots.coherent))*dots.height;

% make sure we haven't gone off the patch
% do the dots separately for left and right hand side
dots.x(dots.x < dots.xmin) = dots.x(dots.x < dots.xmin)+dots.width;
dots.x(dots.x > dots.xmax) = dots.x(dots.x > dots.xmax)-dots.width;
dots.y(dots.y < dots.ymin) = dots.y(dots.y < dots.ymin)+dots.height;
dots.y(dots.y > dots.ymax) = dots.y(dots.y > dots.ymax)-dots.height;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots for opticflow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updateDotsOpticflow(dots)

% get the coherent and incoherent dots
dots.isIncoherent = rand(1,dots.n) > dots.coherence;
dots.incoherentn = sum(dots.isIncoherent);
dots.isCoherent = ~dots.isIncoherent;

% generate a random transformation matrix for each incoherent point
dots.randT = rand(3,dots.incoherentn)-0.5;
% and normalize the transformation to have the same length
% (i.e. speed) as the real transformation matrix
dots.randT = sqrt(sum(dots.T.^2))*dots.randT./([1 1 1]'*sqrt(sum(dots.randT.^2)));

% update relative position of dots in 3-space to observer
dots.X(dots.isCoherent) = dots.X(dots.isCoherent)-dots.T(1);
dots.Y(dots.isCoherent) = dots.Y(dots.isCoherent)-dots.T(2);
dots.Z(dots.isCoherent) = dots.Z(dots.isCoherent)-dots.T(3);

% now move the incoherent points according to the random trasnformation
dots.X(dots.isIncoherent) = dots.X(dots.isIncoherent)-dots.randT(1,:);
dots.Y(dots.isIncoherent) = dots.Y(dots.isIncoherent)-dots.randT(2,:);
dots.Z(dots.isIncoherent) = dots.Z(dots.isIncoherent)-dots.randT(3,:);

% get all points that have fallen off the screen
offscreen = dots.Z<dots.minZ;

% and put them at the furthest distance
dots.Z(offscreen) = dots.maxZ;

% get all points that have fallen out of view
offscreen = dots.Z>dots.maxZ;
% and move them to the front plane
dots.Z(offscreen) = dots.minZ;

% put points fallen off the X edge back
offscreen = dots.X < -dots.maxX;
dots.X(offscreen) = dots.X(offscreen)+2*dots.maxX;
offscreen = dots.X > dots.maxX;
dots.X(offscreen) = dots.X(offscreen)-2*dots.maxX;

% put points fallen off the Y edge back
offscreen = dots.Y < -dots.maxY;
dots.Y(offscreen) = dots.Y(offscreen)+2*dots.maxY;
offscreen = dots.Y > dots.maxY;
dots.Y(offscreen) = dots.Y(offscreen)-2*dots.maxY;

% project on to screen
dots.xproj = dots.f*dots.X./dots.Z;
dots.yproj = dots.f*dots.Y./dots.Z;

% stuff to compute median speed
dots.oldx = dots.x;
dots.oldy = dots.y;

% put into screen coordinates
dots.x = dots.xproj*dots.width+dots.xCenter;
dots.y = dots.yproj*dots.height+dots.yCenter;

%medianSpeed = median(sqrt((dots.oldx-dots.x).^2+(dots.oldy-dots.y).^2)*myscreen.framesPerSecond);
%minSpeed = min(sqrt((dots.oldx-dots.x).^2+(dots.oldy-dots.y).^2)*myscreen.framesPerSecond);
%disp(sprintf('min: %f median: %f',minSpeed,medianSpeed));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create dots for linear2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initDotsLinear(dots)

% get the number of dots
dots.n = round(dots.width*dots.height*dots.density);

% get max and min points for dots
dots.xmin = -dots.width/2;
dots.xmax = dots.width/2;
dots.ymin = -dots.height/2;
dots.ymax = dots.height/2;

% get initial position
dots.x = rand(1,dots.n)*dots.width;
dots.y = rand(1,dots.n)*dots.height;

% get the step size
dots.stepsize = dots.speed/dots.framesPerSecond;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the dots speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = setSpeedOpticflow(dots,speed)

% get the step size
dots.speed = speed;
dots.T = [0 0 dots.speed/dots.framesPerSecond];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the dots direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = setDirOpticflow(dots,direction)

% get the step size
dots.T = [0 0 direction*dots.speed/dots.framesPerSecond];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create dots for optic flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initDotsOpticflow(dots)

% focal length to projection plane
% projection plane is defined to be 
% 1 unit wide and high, so with 
% this focal length, we are looking at
% a view of the world with a 90 deg fov
dots.f = .5;

% translation and rotation matrices
dots = setSpeedOpticflow(dots,dots.speed);
dots.R = [0 0 0];

% maximum depth of points
dots.maxZ = 10;dots.minZ = dots.f;
dots.maxX = 10;
dots.maxY = 10;

% make a brick of points
dots.n = round(dots.width*dots.height*dots.density);

% initial position of dots
dots.X = 2*dots.maxX*rand(1,dots.n)-dots.maxX;
dots.Y = 2*dots.maxY*rand(1,dots.n)-dots.maxY;
dots.Z = (dots.maxZ-dots.minZ)*rand(1,dots.n)+dots.minZ;

% get projection on to plane
dots.xproj = dots.f*dots.X./dots.Z;
dots.yproj = dots.f*dots.Y./dots.Z;

% put into screen coordinates
dots.x = dots.xproj*dots.width+dots.xCenter;
dots.y = dots.yproj*dots.height+dots.yCenter;

% set incoherent dots to 0
dots.isIncoherent = rand(1,dots.n) > dots.coherence;
dots.incoherentn = sum(dots.isIncoherent);
dots.isCoherent = ~dots.isIncoherent;

% init random translation
dots.randT = zeros(3,dots.incoherentn);
