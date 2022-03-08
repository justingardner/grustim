% geislerDetectionTask_jy.m
% 
%      usage: geislerDetectionTask_jy.m
%         by: jiwon yeon
%       date: 
%  copyright: (c) 2022 Jiwon Yeon
%    purpose: Replicating Najemnik&Geisler's 2005 study
%
%             The study has two experiments: Detection and Search tasks.
%             The current code only replicates the detection task.
%

% making one line change

function geislerDetectionTask_jy
clear all, close all
myscreen.screenNumber = 2;
myscreen = initScreen(myscreen);

% load pink_filter
if exist([cd '/geislerDetectionTask_pinkFilter.mat']) ~= 0
    load('geislerDetectionTask_pinkFilter.mat');
    task.noise.pink_filter = pink_filter;
else
    task = createPinkFilter(myscreen);
end

% create noise
% % % % % task.noise.luminance = 4;
task.noise.size = 15;   % visual angle

task.noise.contrasts = [0, .05, .10, .20];
% task.noise.contrasts = 0;

index = randsample(length(task.noise.contrasts),1);
task.thistrial.noise_contrast = task.noise.contrasts(index);
task = createPinkNoise(myscreen, task);

% create a gabor patch
task.gabor.size = 2;    % visual angle
task.gabor.tilt = 315;
task.gabor.cycle = 6;

% task.gabor.contrasts = [.2, .1, .075, .05];
task.gabor.contrasts = 1;

index = randsample(length(task.gabor.contrasts),1);
task.thistrial.gabor_contrast = task.gabor.contrasts(index);
task = createGabor(task);

% define locations
task.nLoc = 9;   % 25 for the real experiment
task = defineLocations(task);
task.thistrial.gabor_location = task.locations(randsample(1:task.nLoc,1),:);

% combine noise and gabor
task = combinedStimulus(task);

% make stencils
% noise stencil
mglClearScreen(.5);
mglStencilCreateBegin(1);
mglVisualAngleCoordinates(myscreen.displayDistance,myscreen.displaySize);
mglFillOval(0, 0, [task.noise.size, task.noise.size]);
mglStencilCreateEnd;

% show stimulus
mglClearScreen(task.thistrial.bg_color);
mglStencilSelect(1);
tex_noise = mglCreateTexture(task.thistrial.im);
mglBltTexture(tex_noise,[0 0]);
mglStencilSelect(0);

mglFlush



%%%%%%%%%% helper functions
function task = createPinkFilter(myscreen)
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
task.noise.pink_filter = pink_filter;

function task = createPinkNoise(myscreen, task)
w = myscreen.screenWidth;
h = myscreen.screenHeight;
sz = max(w,h);

% make the odd size of the image
if mod(sz,2)==0, sz = sz-1; end

% fft on white noise
white = randn(sz,sz);
fwhite = fftshift(fft2(white));
phase = angle(fwhite);

% % % % % % change luminance
% % % % % task.noise.pink_filter(ceil(size(task.noise.pink_filter,1)), ...
% % % % %     ceil(size(task.noise.pink_filter,2))) = task.noise.luminance;
% % % % % % % % % % %  final image has imaginary values when change the DC

% create new magnitude
new_mag = fwhite .* task.noise.pink_filter;
new_Fourier = new_mag .* (cos(phase) + sqrt(-1)*sin(phase));
im = ifft2(ifftshift(new_Fourier));

% change contrast
contrast = task.thistrial.noise_contrast;
N = length(im(:));
m_im = mean(im(:));
coeff = sqrt((N*contrast^2) / sum((im(:)-m_im).^2));
task.noise.im = coeff .* im;

function task = createGabor(task)
grating = mglMakeGrating(task.gabor.size, task.gabor.size, task.gabor.cycle, ...
    task.gabor.tilt, 0);
contrast = task.thistrial.gabor_contrast;
grating = grating * contrast;
gaussian = mglMakeGaussian(task.gabor.size, task.gabor.size, 1, 1);
task.gabor.im = (grating.*gaussian)/2;

function task = defineLocations(task)
% determine how many layers to have
% maximum number of locations per layer is 8
nLoc = task.nLoc;
if nLoc > 8 && mod(nLoc,8) ~= 0 
    if mod(nLoc,8) < 4
        nLoc = floor(nLoc/8) * 8 + 1;
    else 
        nLoc = ceil(nLoc/8) * 8 + 1;
    end
    task.nLoc = nLoc;
    nLayer = floor(nLoc/8);
elseif nLoc <= 8
    nLayer = 1;
end

% create x and y positions
[x y] = meshgrid(-task.noise.size:.1:task.noise.size, ...
    -task.noise.size:.1:task.noise.size);

% radius
rim = task.noise.size - 2;  % make some margin 
radius = round(linspace(0,rim,nLayer+2),1);
radius = radius(2:end-1);

% theta
if nLoc < 8
    theta = linspace(0, 2*pi, nLoc+1);
else
    theta = linspace(0, 2*pi, 9);
end
theta(end) = [];

% determine locations
locations = [0, 0];
cTheta = 0; 
cRadius = 1;
for cLoc = 1:nLoc-1
    cTheta = cTheta + 1;
    x_pos = radius(cRadius) * cos(theta(cTheta));
    y_pos = radius(cRadius) * sin(theta(cTheta));    
    locations = [locations; [x_pos, y_pos]];
    
    if cTheta == 8, cTheta = 0; end
    if mod(cLoc,8) == 0, cRadius = cRadius+1; end
end
task.locations = locations;

function task = combinedStimulus(task)
noise = task.noise.im;
gabor = task.gabor.im;
c_location = task.thistrial.gabor_location;   % gabor's center

% make circular gabor patch
radius = floor(min(size(gabor))/2);
[stencil_x, stencil_y] = meshgrid(-floor(size(gabor,1)/2):floor(size(gabor,1)/2), ...
    -floor(size(gabor,2)/2):floor(size(gabor,2)/2));
stencil = (sqrt(stencil_x.^2 + stencil_y.^2) <= radius);
gabor_circle = stencil' .* gabor;

% determine the location to display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% how to change the location in to visual angle?  %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

location = round([c_location(1) + floor(size(noise,1)/2), c_location(2) + floor(size(noise,2)/2)]);
x_loc = [location(1)-floor(size(gabor_circle,1)/2), location(1)+floor(size(gabor_circle,1)/2)];
y_loc = [location(2)-floor(size(gabor_circle,2)/2), location(2)+floor(size(gabor_circle,2)/2)];
gabor_position = zeros(size(noise,1), size(noise,2));
gabor_position(x_loc(1):x_loc(2), y_loc(1):y_loc(2)) = gabor_circle;

% add gabor to the noise
final_im = noise + gabor_position;

% scale it to [0 255]
final_im =  255 * (final_im - min(final_im(:))) ./ (max(final_im(:)) - min(final_im(:)));
task.thistrial.im = final_im;

% decide background color
if task.thistrial.noise_contrast == 0
    bg_color = final_im(1,1);    
else    
    bg_color = mean(final_im(:));
end
task.thistrial.bg_color = bg_color;


