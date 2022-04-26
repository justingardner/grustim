function geislerDetectionTask_verify_locations

mglClose        % close MGL if it's open
clear all, close all, clc

mglSetSID('test')
myscreen.screenNumber = 2;
% myscreen.displayName = 'bowmore';
myscreen.saveData = 0;

mglSetParam('abortedStimfilesDir', '~/proj/data/geislerDetectionTask/aborted',1);
myscreen = initScreen(myscreen);  

load('geislerDetectionTask_pinkFilter.mat');
stimulus.pink_filter = pink_filter;

%%%%% set stimulus parameter
stimulus.responsekeys = [44,48];   % space bar
stimulus.noise.size = 15;   % visual angle
stimulus.noise.contrasts = [.2];

stimulus.gabor.size = 1;    % visual angle
stimulus.gabor.tilt = 315;
stimulus.gabor.cycle = 6;
nLoc = 25;   % 25 for the real experiment
stimulus.gabor.contrasts = .5;


%% create noise
noise_with_buffer = stimulus.noise.size + 3;    % visual angle
noise_frame_pixel = visualAngleToPixels(noise_with_buffer, ...
    [myscreen.screenWidth, myscreen.screenHeight]);
% make the size of the image an odd number
if mod(noise_frame_pixel,2)==0, noise_frame_pixel = noise_frame_pixel+1; end    

filter_sz = size(stimulus.pink_filter);
pink_filter = stimulus.pink_filter(...
    floor(filter_sz(1)/2)+1-(noise_frame_pixel-1)/2:floor(filter_sz(1)/2)+1+(noise_frame_pixel-1)/2, ...
    floor(filter_sz(2)/2)+1-(noise_frame_pixel-1)/2:floor(filter_sz(2)/2)+1+(noise_frame_pixel-1)/2);
    
% fft on white noise
white = randn(noise_frame_pixel, noise_frame_pixel);
fwhite = fftshift(fft2(white));
phase = angle(fwhite);
    
% create new magnitude
new_mag = fwhite .* pink_filter;
new_Fourier = new_mag .* (cos(phase) + sqrt(-1)*sin(phase));
im = ifft2(ifftshift(new_Fourier));

N = length(im(:));
m_im = mean(im(:));
coeff = sqrt((N*stimulus.noise.contrasts^2) / sum((im(:)-m_im).^2));
noise = coeff .* im;

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
gabor_locations_va = locations; 
gabor_locations_onNoise = visualAngleToPixels(gabor_locations_va, noise_frame_pixel);

% gabor_locations_onNoise(end-1,:) = [273 80];


%% create gabor
grating = mglMakeGrating(stimulus.gabor.size, stimulus.gabor.size, ...
    stimulus.gabor.cycle, stimulus.gabor.tilt, 0);
contrast = stimulus.gabor.contrasts;
grating = grating * contrast;
gaussian = mglMakeGaussian(stimulus.gabor.size, stimulus.gabor.size, 1, 1);
gabor = (grating.*gaussian);

%% superimpose gabor on the noise image - on all 25 locations
% make circular gabor patch
radius_pixel = visualAngleToPixels(stimulus.gabor.size/2);
gabor_sz = size(gabor);
[x y] = meshgrid(-(gabor_sz(1)-1)/2:(gabor_sz(1)-1)/2, ...
    -(gabor_sz(2)-1)/2:(gabor_sz(2)-1)/2);
stencil = sqrt(x.^2 + y.^2) <= radius_pixel;
gabor_circle = stencil' .* gabor;

gabor_position = zeros(size(noise,1), size(noise,2));
for loc = 1:nLoc
    pos = gabor_locations_onNoise(loc,:);
    gabor_position([pos(1)-(gabor_sz(1)-1)/2:pos(1)+(gabor_sz(1)-1)/2], ...
        [pos(2)-(gabor_sz(2)-1)/2:pos(2)+(gabor_sz(2)-1)/2]) = gabor_circle;
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

%% present screen
mglClearScreen(.5);
tex = mglCreateTexture(final_im);
mglStencilSelect(1)
mglBltTexture(tex,[0,0])
mglStencilSelect(0)

for loc = 1:nLoc;
    mglGluAnnulus(gabor_locations_va(loc,1), gabor_locations_va(loc,2), ...
            stimulus.gabor.size/2-.07, stimulus.gabor.size/2, ...
            [1 1 1], 120, 2)
end
mglFlush








