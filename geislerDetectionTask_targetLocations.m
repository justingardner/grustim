% geislerDetectionTask_targetLocations.m
% 
%      usage: geislerDetectionTask_targetLocations(nTarget).m
%         by: jiwon yeon
%       date: 
%  copyright: (c) 2022 Jiwon Yeon
%    purpose: Showing target locations for Najemnik&Geisler's 2005 study
%             The original study used 25 targets for the detection task
%


function geislerDetectionTask_targetLocations

clear all, close all
myscreen.screenNumber = 2;
myscreen = initScreen(myscreen);

% load pink_filter
load('geislerDetectionTask_pinkFilter.mat');
task.noise.pink_filter = pink_filter;

% create noise
task.noise.size = 15;   % visual angle
task.noise.contrasts = .05;

task.thistrial.noise_contrast = task.noise.contrasts;
task = createPinkNoise(myscreen, task);

% define locations
nTarget = 25;
task.gabor.nLoc = nTarget;   % 25 for the real experiment
task = defineLocations(task);

task = combinedStimulus(task);

% make stencils
% noise stencil
mglClearScreen(.5);
mglStencilCreateBegin(1);
mglVisualAngleCoordinates(myscreen.displayDistance,myscreen.displaySize);
mglFillOval(0, 0, [task.noise.size, task.noise.size]);
mglStencilCreateEnd;

% show stimulus
mglClearScreen(.5);
mglStencilSelect(1);
tex_noise = mglCreateTexture(task.final_im);
mglBltTexture(tex_noise,[0 0]);
mglStencilSelect(0);

mglFlush


%%%%%%%%%% helper functions
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


function task = defineLocations(task)
% determine how many layers to have
% maximum number of locations per layer is 8
nLoc = task.gabor.nLoc;
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

radius_va = linspace(0, task.noise.size/2+1, nLayer+2);     % radius in visual angle
radius_va = radius_va(2:end-1);

% theta
if nLoc < 8
    theta = linspace(0, 2*pi, nLoc+1);
else
    theta = linspace(0, 2*pi, 9);
end
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
clear locations

% convert the locations from visual angle to pixel
% get parameters
pix_noise_height = size(task.noise.im,1);
pix_noise_width = size(task.noise.im,2);
x_nPixPerCm = mglGetParam('xDeviceToPixels');
y_nPixPerCm = mglGetParam('yDeviceToPixels');
screen_distance = mglGetParam('devicePhysicalDistance');

% convert angle to cm
locations_cm = 2 .* screen_distance .* tan(locations_va ./2 .* (pi/180));

% convert cm to pixel
locations_pix = round([locations_cm(:,1) .* x_nPixPerCm, ...
    locations_cm(:,2) .* y_nPixPerCm]);

% centering
locations = [locations_pix(:,1) + ceil(pix_noise_width/2), ...
    locations_pix(:,2) + ceil(pix_noise_height/2)];

task.gabor.locations = locations;


function task = combinedStimulus(task)
noise = task.noise.im;
location = task.gabor.locations;

% dot_positions = zeros(size(noise,1), size(noise,2));
for loc = 1:task.gabor.nLoc;
    noise(location(loc,1)-5:location(loc,1)+5, ...
        location(loc,2)-5:location(loc,2)+5) = 1;
end

% add gabor to the noise
final_im = noise;

% scale it to [0 255]
final_im =  255 * ((final_im + 1) ./ 2);
task.final_im = final_im';





