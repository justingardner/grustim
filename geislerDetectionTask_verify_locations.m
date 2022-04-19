function geislerDetectionTask_verify_locations

mglClose        % close MGL if it's open
clear all, close all, clc

myscreen.screenNumber = 2;
myscreen.saveData = 0;

mglSetParam('abortedStimfilesDir', '~/proj/data/geislerDetectionTask/aborted',1);
myscreen = initScreen(myscreen);  

load('geislerDetectionTask_pinkFilter.mat');
stimulus.pink_filter = pink_filter;

%%%%% set stimulus parameter
stimulus.responsekeys = [44,48];   % space bar
stimulus.noise.size = 15;   % visual angle
stimulus.noise.contrasts = [.2];

stimulus.gabor.size = .5;    % visual angle
stimulus.gabor.tilt = 315;
stimulus.gabor.cycle = 6;
nLoc = 25;   % 25 for the real experiment
stimulus.gabor.contrasts = .05;

%% define locations
nLayer = floor((nLoc-1)/8);

radius_va = linspace(0, stimulus.noise.size/2+1, nLayer+2);     % radius in visual angle
radius_va = radius_va(2:end-1);

% theta
theta = linspace(0, 2*pi, 9);
theta(end) = [];

% determine locations - in visual angle
locations = [0, 0];
cTheta = 0;     % current theta
cLayer = 1;     % current layer
for cLoc = 1:nLoc-1
    cTheta = cTheta + 1;
    x_pos = radius_va(cLayer) * cos(theta(cTheta));
    y_pos = radius_va(cLayer) * sin(theta(cTheta));    
    locations = [locations; [x_pos, y_pos]];
    
    if cTheta == 8, cTheta = 0; end
    if mod(cLoc,8) == 0, cLayer = cLayer+1; end
end
locations_va = locations; 

% convert visual angle of the locations to pixels
displaySize = max(size(stimulus.pink_filter));
gabor_locations = visualAngleToPixels(locations_va, displaySize);

%% create noise
w = myscreen.screenWidth;
h = myscreen.screenHeight;
sz = max(w,h);

% make the odd size of the image
if mod(sz,2)==0, sz = sz-1; end

% fft on white noise
white = randn(sz,sz);
fwhite = fftshift(fft2(white));
phase = angle(fwhite);
    
% create new magnitude
new_mag = fwhite .* stimulus.pink_filter;
new_Fourier = new_mag .* (cos(phase) + sqrt(-1)*sin(phase));
im = ifft2(ifftshift(new_Fourier));

N = length(im(:));
m_im = mean(im(:));
coeff = sqrt((N*stimulus.noise.contrasts^2) / sum((im(:)-m_im).^2));
noise = coeff .* im;

%% create gabor
grating = mglMakeGrating(stimulus.gabor.size, stimulus.gabor.size, ...
    stimulus.gabor.cycle, stimulus.gabor.tilt, 0);
contrast = stimulus.gabor.contrasts;
grating = grating * contrast;
gaussian = mglMakeGaussian(stimulus.gabor.size, stimulus.gabor.size, 1, 1);
gabor = (grating.*gaussian);

%% superimpose gabor on the noise image - on all 25 locations
% make circular gabor patch
radius_va = stimulus.gabor.size/2;
radius_cm = 2 .* mglGetParam('devicePhysicalDistance') * tan(radius_va /2 * (pi/180));
radius_pix = radius_cm * min([mglGetParam('xDeviceToPixels'), mglGetParam('yDeviceToPixels')]);
radius = round(radius_pix);
[stencil_x, stencil_y] = meshgrid(-ceil(size(gabor,1)/2)+1:ceil(size(gabor,1)/2)-1, ...
    -ceil(size(gabor,2)/2)+1:ceil(size(gabor,2)/2)-1);
stencil = (sqrt(stencil_x.^2 + stencil_y.^2) <= radius);
gabor_circle = stencil' .* gabor;

gabor_position = zeros(size(noise,1), size(noise,2));
for loc = 1:nLoc    
    % determine the location to display
    x_lims = [gabor_locations(loc,1)-ceil(size(gabor_circle,1)/2)+1, gabor_locations(loc,1)+ceil(size(gabor_circle,1)/2)-1];
    y_lims = [gabor_locations(loc,2)-ceil(size(gabor_circle,2)/2)+1, gabor_locations(loc,2)+ceil(size(gabor_circle,2)/2)-1];
        
    gabor_position(x_lims(1):x_lims(2),y_lims(1):y_lims(2)) = gabor_circle;
end

% add gabor to the noise
final_im = noise + gabor_position;

% scale it to [0 255], for both stimulus images
final_im =  255 .* ((final_im + 1) ./ 2);

%% Make stencil
mglClearScreen(.5);
mglStencilCreateBegin(1);
mglVisualAngleCoordinates(myscreen.displayDistance,myscreen.displaySize);
mglFillOval(0, 0, [stimulus.noise.size, stimulus.noise.size]);
mglStencilCreateEnd;
mglClearScreen(.5);

%% display the image and white circle
mglClearScreen(.5)
tex = mglCreateTexture(final_im);
mglStencilSelect(1);
mglBltTexture(tex,[0 0])
mglStencilSelect(0);
mglFillOval(0,0,[.2 .2],0)

for loc = [2, 10, 18];
% loc = 18;
    sz = size(stimulus.pink_filter);
    target_location = gabor_locations(loc,:);
    target_location = pixelsToVisualAngle(target_location,sz);
    mglGluAnnulus(target_location(1), target_location(2), .35, .4, ...
        [1 1 1], 120, 2)
end
mglFlush;






