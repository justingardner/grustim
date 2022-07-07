function geisler_verify_locations

mglClose        % close MGL if it's open
clear all, close all, clc

mglSetSID('test')
myscreen.screenNumber = 2;
% myscreen.displayName = 'bowmore';
myscreen.saveData = 0;
testing = 'Search';  % Detection or Search
disp(['Verify locations of the ' testing ' task'])

mglSetParam('abortedStimfilesDir', '~/proj/data/geislerDetectionTask/aborted',1);
myscreen = initScreen(myscreen);

%%%%% set stimulus parameter
stimulus.responsekeys = [44,48];   % space bar
stimulus.noise.size = 15;   % visual angle
stimulus.noise.size_pix = visualAngleToPixels(stimulus.noise.size, ...
    [myscreen.screenWidth, myscreen.screenHeight]);
stimulus.noise.contrasts = [.2];

stimulus.gabor.size = 1;    % visual angle
stimulus.gabor.tilt = 315;
stimulus.gabor.cycle = 6;
stimulus.gabor.contrasts = .5;


%% create noise
noise_with_buffer = stimulus.noise.size + 3;    % visual angle
noise_frame_pixel = visualAngleToPixels(noise_with_buffer, ...
    [myscreen.screenWidth, myscreen.screenHeight]);
% make the size of the image an odd number
if mod(noise_frame_pixel,2)==0, noise_frame_pixel = noise_frame_pixel+1; end

% pink filter is already saved as a file. just load the variable
if exist('geislerDetectionTask_pinkFilter.mat') == 2
    load('geislerDetectionTask_pinkFilter.mat')
end

% pink filter has not created before or needs to re-created
if exist('geislerDetectionTask_pinkFilter.mat') == 0 || size(pink_filter,1) < double(noise_frame_pixel)
    clear pink_filter
    fprintf('[geisler] Creating pink filter... \n')
    pink_filter = createPinkFilter(myscreen);
    save('geislerDetectionTask_pinkFilter.mat', 'pink_filter')
    fprintf('[geisler] Process done! \n')
else
    fprintf('[geisler] Pink filter size is compatible with the current monitor setup \n')
end
stimulus.pink_filter = pink_filter;

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
if strcmp(testing, 'Detection')
    nLoc = 25;   % 25 for the real experiment
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
    
elseif strcmp(testing, 'Search')
    h_dist = 1.5;     % horizontal distance betweeen any two gabors
    v_dist = sqrt(h_dist^2-(h_dist/2)^2);
    y_lim = stimulus.noise.size - 1;    
    
    trigrid = [];
    y_current = 0;
    xx = 0;
    displacement = 0;
    while y_current < y_lim
        if displacement == 0
            xx = [0:h_dist:y_lim]';
            yy = ones(length(xx), 1)*y_current;
            displacement = 1;
        else
            xx = [h_dist/2:h_dist:y_lim]';
            yy = ones(length(xx), 1)*y_current;
            displacement = 0;
        end
        trigrid = [trigrid;[xx, yy]];
        y_current = y_current+v_dist;
    end

    trigrid = trigrid - repmat(max(trigrid)./2,size(trigrid,1),1);
    inside = sqrt(trigrid(:,1).^2+trigrid(:,2).^2) <= y_lim/2;    
    
    locations = trigrid(inside,:);
    nLoc = size(locations,1);
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
final_im = final_im;

%% Make stencil
fprintf('[geisler] Creating stencil \n')
mglClearScreen(.5);
mglFlush();

mglStencilCreateBegin(1);
mglVisualAngleCoordinates(myscreen.displayDistance,myscreen.displaySize);
mglFillOval(0, 0, [stimulus.noise.size, stimulus.noise.size]);
% mglFillOval(0, 0, [30 30]);
mglStencilCreateEnd;
fprintf('[geisler] Process done! \n')

%% present screen
fprintf('[geisler] Generating gabors at each locations \n')
mglClearScreen(.5);
tex = mglCreateTexture(final_im);
mglStencilSelect(1);
mglBltTexture(tex,[0,0]);
mglStencilSelect(0);

for loc = 1:nLoc
    mglGluAnnulus(gabor_locations_va(loc,2), gabor_locations_va(loc,1), ...
        stimulus.gabor.size/2-.07, stimulus.gabor.size/2, ...
        [1 1 1], 120, 2);
    mglTextSet([],20,1);    
    mglTextDraw(num2str(loc), ...
        [gabor_locations_va(loc,2) gabor_locations_va(loc,1)]);
end
fprintf('[geisler] Process done! \n')

mglFlush();


%%%%%% helper function
function pink_filter = createPinkFilter(myscreen)
w = myscreen.screenWidth;
h = myscreen.screenHeight;
sz = max(w,h);

% make the odd size of the image
if mod(sz,2)==0, sz = sz-1; end

% make pink filter
last_freq = ceil(sz/2);
pink_filter = zeros(sz,sz);
[x y] = meshgrid(-ceil(sz/2)+1:ceil(sz/2)-1, -ceil(sz/2)+1:ceil(sz/2)-1);
index = sqrt(x.^2 + y.^2);
for f = 1:last_freq
    pink_filter(index > f-1 & index < f+1) = 1/f;
end

