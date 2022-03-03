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
myscreen.screenNumber = 2;
myscreen = initScreen(myscreen);
mglVisualAngleCoordinates(57,[51 38]);

% build pink noise background
noise = createPinkNoise(myscreen);

% create a gabor patch
task.gabor.size = 2;    % visual angle
task.gabor.tilt = 45;
task.gabor.cycle = 6;
gabor = createGabor(task);

% define locations

% make stencils
mglStencilCreateBegin(1);
mglFillOval(0,0,[15,15]);
mglStencilCreateEnd;
mglClearScreen(.5);

mglStencilCreateBegin(2);
mglFillOval(loc(1), loc(2), [task.gabor.size, task.gabor.size]);

% 
mglClose
mglOpen(2);
mglClearScreen(.5);
mglVisualAngleCoordinates(57,myscreen.displaySize);
texture = mglCreateTexture(gabor);
mglBltTexture(texture, [0 0 5 5], 0, 0, 0)
mglFlush


mglFlush;
mglClose


%%%%%%%%%% helper functions
function norm_im = createPinkNoise(myscreen)

w = myscreen.screenWidth;
h = myscreen.screenHeight;
sz = max(w,h);

% make the odd size of the image
if mod(sz,2)~=0, sz = sz-1; end

% fft on white noise
white = randn(sz,sz);
fwhite = fftshift(fft2(white));
phase = angle(f_white);

% make pink filter
pink_filter = zeros(sz,sz);
[x y] = meshgrid(-ceil(sz/2)+1:ceil(sz/2)-1, -ceil(sz/2)+1:ceil(sz/2)-1);
last_freq = ceil(sz/2);
for f = 1:last_freq
    index = (sqrt(x.^2 + y.^2) > f-1) & (sqrt(x.^2 + y.^2) < f+1);
    pink_filter(index) = 1/f;
end

% create new magnitude
new_mag = fwhite .* pink_filter;
new_Fourier = new_mag .* (cos(phase) + sqrt(-1)*sin(phase));
im = ifft2(ifftshift(new_Fourier));
norm_im = (im-min(im(:))) / (max(im(:)) - min(im(:))) .* 255; 

function gabor = createGabor(task)
grating = mglMakeGrating(task.gabor.size, task.gabor.size, task.gabor.cycle, ...
    task.gabor.tilt, 0);
gaussian = mglMakeGaussian(task.gabor.size, task.gabor.size, 1, 1);
gabor = 255*(grating.*gaussian+1)/2;




