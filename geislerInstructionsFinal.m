%  geislerInstructionFinal.m
% 
%      usage: geislerInstructionFinal.m
%       date: Aug 2022
%    purpose: Provide instructions for participants before starting
%             geislerDetectionTask_staircase_jy.m

function geislerInstructionsFinal
%% POTENTIAL BUGS
%
%  (HOLD) the white stencils both show up upon opening program (hold)
%  (HOLD) I'm not sure if stimulus.gabor.size = 1 is the correct stimulus
%  size. I tried editing it to .4 to match gluAnnulus size but then it
%  seemed way too small? regardless, should be easy to fix!

%% set up global variables 
% initalize 
mglClose        % close MGL if it's open
clear all, close all, clc

mglSetSID('test')
myscreen.screenNumber = 1;
myscreen.saveData = 0;
myscreen.keyboard.left = 44;
myscreen.keyboard.right = 48;
myscreen = initScreen(myscreen);

%%%%% set stimulus parameter
stimulus.responsekeyleft = 44;  
stimulus.responsekeyright = 48;  
stimulus.responsekeys = [44,48];   % '<,' & '>.'
stimulus.noise.size = 15;   % visual angle
stimulus.noise.contrasts = .2;

stimulus.gabor.size = 1;    % visual angle
stimulus.gabor.tilt = 315;
stimulus.gabor.cycle = 6;
stimulus.gabor.contrasts = .25;

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

% noise x2
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
noise2 = coeff .* im;
 %% define locations

nLoc = 1;

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

gabor_locations_va = [0,0]; 
%gabor_locations_va = locations; // changed
gabor_locations_onNoise = visualAngleToPixels(gabor_locations_va, noise_frame_pixel);


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
noise_only = noise2;


% scale it to [0 255], for both stimulus images
final_im =  255 .* ((final_im + 1) ./ 2);
final_im = final_im;

noise_only = 255 .* ((noise_only + 1) ./ 2);
 

%% create stencil - double (stencil 1)

% reset
mglOpen(1);
mglClearScreen(.5);
mglFlush();

% create stencil  

mglStencilCreateBegin(1);
mglVisualAngleCoordinates(myscreen.displayDistance,myscreen.displaySize);
x = [-stimulus.noise.size/2, stimulus.noise.size/2 + 1];
y = [0, 0];
mglFillOval(x, y, [stimulus.noise.size, stimulus.noise.size]);
mglStencilCreateEnd;


%% create stencil - single (stencil 2)

% create stencil 2

% reset
mglOpen(1);
mglClearScreen(.5);
mglFlush();

% create stencil 2

mglStencilCreateBegin(2);
mglVisualAngleCoordinates(myscreen.displayDistance,myscreen.displaySize);
mglFillOval(0, 0, [stimulus.noise.size, stimulus.noise.size]);
mglStencilCreateEnd;

%% change backtick to space bar 
backTickCharacter = ' ';
backTickCode = mglCharToKeycode({' '});
myscreen.keyboard.backtick = backTickCode;
%% present screen with text instructions

mglClearScreen(.5);

mglTextDraw('Thank you for participating in our research experiment.' ,[0,4]);
mglTextDraw('Before beginning the experiment, we will first give you brief instructions and a' ,[0,0]);
mglTextDraw('chance to test your understanding of the task.' ,[0,-2]);
mglTextDraw('Press the SPACE BAR when you are ready to continue to the next instruction screen.' ,[0,-14]);


mglFlush; 

while 1
    k = mglGetKeys;
    if k(myscreen.keyboard.backtick)==1, break; end
end

mglTextDraw(['In this experiment, we will present two screens in quick succession. One screen will have a target while the other will not.'] ,[0,4]);
mglTextDraw(['Your task is to respond with which screen had the target.'] ,[0,2]);
mglTextDraw(['While you are performing the experiment, we will concurrently collect your eye movement data.'] ,[0,-2]);
mglTextDraw(['With the collected data, we aim to build a model that can predict eye movements. '] ,[0,-4]);
mglTextDraw(['Press the SPACE BAR when you are ready to continue to the next instruction screen.'] ,[0,-14]);

mglFlush;

while 1
    k = mglGetKeys;
    if k(myscreen.keyboard.backtick)==1, break; end
end

mglTextDraw(['When a trial begins, you will see a black dot like the one below at the center of the screen'] ,[0,14]);
mglFillOval(0,0,[.2 .2],0)
mglTextDraw('Please focus your eyes at the center of the screen for the entire duration of the experiment.' ,[0,12]);
mglTextDraw(['Press the SPACE BAR when you are ready to continue to the next instruction screen.'] ,[0,-14]);
mglFlush;

while 1
    k = mglGetKeys;
    if k(myscreen.keyboard.backtick)==1, break; end
end

mglTextDraw(['In each trial, two stimuli will flash on the screen in sequence.'], [0, 14]);
mglTextDraw(['A white dot will be briefly presented between screens.'], [0, 12]);
mglTextDraw(['Example stimuli are shown on the screen.'] ,[0,10]);

% PRESENT EXAMPLE STIMULI

tex = mglCreateTexture(final_im);
tex2 = mglCreateTexture(noise_only);
mglStencilSelect(1);
mglBltTexture([tex tex2],[-stimulus.noise.size/2,0; stimulus.noise.size/2 + 1,0]);
mglStencilSelect(0);

%mglFlush();

%%%

mglTextDraw(['As you can see, the LEFT stimulus has a target at the center while the RIGHT stimulus does not.'] ,[0,-10]);
mglTextDraw('Press the SPACE BAR when you are ready to continue to the next instruction screen.' ,[0,-14]);
mglFlush;

while 1
    k = mglGetKeys;
    if k(myscreen.keyboard.backtick)==1, break; end
end

mglTextDraw('After presenting the two stimuli, you will be asked to identify which stimulus contained the target.' ,[0, 10]);
mglTextDraw('Which screen showed the target?',[0,0])
mglTextDraw('1(<)  or  2(>)',[0, -2]);

% PRESENT EXAMPLE STIMULI

%tex = mglCreateTexture(final_im);
%tex2 = mglCreateTexture(noise_only);
%mglStencilSelect(1);
%mglBltTexture([tex tex2],[-stimulus.noise.size/2,0; stimulus.noise.size/2 + 1,0]);
%mglStencilSelect(0);

%%%

mglTextDraw('As the message indicates, press the COMMA button (<) to indicate that the FIRST image contained.' ,[0,-8]);
mglTextDraw(['the target, and the PERIOD button (>) to indicate that the SECOND image contained the target.'] ,[0,-10]);
mglTextDraw(['Press the SPACE BAR when you are ready to continue to the next instruction screen.'] ,[0,-14]);
mglFlush;

while 1
    k = mglGetKeys;
    if k(myscreen.keyboard.backtick)==1, break; end
end

mglTextDraw(['If your answer is correct, a green dot will appear.'] ,[0,4]);
mglFillOval(0,0,[.2 .2],[0 1 0])
mglTextDraw(['If your answer is incorrect, a red dot will appear. '] ,[0,-4]);
%mglFillOval(0,0,[.2 .2],[1 0 0])
mglTextDraw(['Then the dot will turn white to notify you to prepare for the next trial.'] ,[0,-6]);
%mglFillOval(0,-8,[.2 .2],1)
mglTextDraw(['Press the SPACE BAR when you are ready to continue to the next instruction screen.'] ,[0,-14]);
mglFlush;

while 1
    k = mglGetKeys;
    if k(myscreen.keyboard.backtick)==1, break; end
end

mglTextDraw(['Now we will give you three example trials.'] ,[0, 0]);
mglTextDraw(['Press the SPACE BAR when you are ready to continue to the next instruction screen.'] ,[0, -14]);
mglFlush;

while 1    
    k = mglGetKeys;
    if k(myscreen.keyboard.backtick)==1, break; end
end

%% Three easy trials with high contrast value and long display time
 
mglTextDraw('THREE EXAMPLE TRIALS' ,[0,0]);
mglFlush;
mglWaitSecs(1);

% PRESENT EXAMPLE STIMULI
% trial 1 set up

mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],0)
mglFlush();

% trial 1 actual start 
while 1
    keycode = mglGetKeys;
    if any(keycode==1), break; end
end
%%
% insert

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
noise2 = coeff .* im;

grating = mglMakeGrating(stimulus.gabor.size, stimulus.gabor.size, ...
    stimulus.gabor.cycle, stimulus.gabor.tilt, 0);
contrast = 0.5;
grating = grating * contrast;
gaussian = mglMakeGaussian(stimulus.gabor.size, stimulus.gabor.size, 1, 1);
gabor = (grating.*gaussian);

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
noise_only = noise2;


% scale it to [0 255], for both stimulus images
final_im =  255 .* ((final_im + 1) ./ 2);
final_im = final_im;

noise_only = 255 .* ((noise_only + 1) ./ 2);
%%
tex = mglCreateTexture(final_im);
tex2 = mglCreateTexture(noise_only);
mglStencilSelect(2);
mglBltTexture(tex,[0,0]);
mglFlush();
mglWaitSecs(1);
mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],1)
mglFlush();
mglWaitSecs(0.5); 
mglBltTexture(tex2,[0,0]);
mglStencilSelect(0);
mglFlush();
mglWaitSecs(1); 
mglClearScreen(.5)
mglTextSet([],32,1);
mglFillOval(0,0,[.2 .2],1)
mglTextDraw('Which screen showed the target?',[0,-8])
mglTextDraw(['1(<)  or  2(>)'],[0 -10]);
mglFlush();

% get response 
%myscreen.keyboard.left = mglCharToKeycode({','});
stimulus.responsekeys = [44,48]; 

while 1
    k = mglGetKeys;
    if k(44)==1, mglFillOval(0,0,[.2 .2],[0 1 0]), break; end
    if k(48)==1, mglFillOval(0,0,[.2 .2],[1 0 0]), break; end
end
mglFlush();
mglWaitSecs(0.75);
mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],0) 
mglFlush();

% trial 2
while 1
    keycode = mglGetKeys;
    if any(keycode==1), break; end
end
%%
%%
% insert

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
noise2 = coeff .* im;

grating = mglMakeGrating(stimulus.gabor.size, stimulus.gabor.size, ...
    stimulus.gabor.cycle, stimulus.gabor.tilt, 0);
contrast = 0.5;
grating = grating * contrast;
gaussian = mglMakeGaussian(stimulus.gabor.size, stimulus.gabor.size, 1, 1);
gabor = (grating.*gaussian);

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
noise_only = noise2;


% scale it to [0 255], for both stimulus images
final_im =  255 .* ((final_im + 1) ./ 2);
final_im = final_im;

noise_only = 255 .* ((noise_only + 1) ./ 2);
%%
tex = mglCreateTexture(final_im);
tex2 = mglCreateTexture(noise_only);
mglStencilSelect(2);
mglBltTexture(tex2,[0,0]);
mglFlush();
mglWaitSecs(1);

mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],1)
mglFlush();
mglWaitSecs(0.5); 

mglBltTexture(tex,[0,0]);
mglStencilSelect(0);
mglFlush();
mglWaitSecs(1);
mglClearScreen(.5)
mglTextSet([],32,1);
mglFillOval(0,0,[.2 .2],1)
mglTextDraw(['Which screen showed the target?'],[0,-8])
mglTextDraw('1(<)  or  2(>)',[0 -10]);
mglFlush();

% get response 
while 1
    k = mglGetKeys;
    if k(44)==1, mglFillOval(0,0,[.2 .2],[1 0 0]), break; end
    if k(48)==1, mglFillOval(0,0,[.2 .2],[0 1 0]), break; end
end

mglFlush();
mglWaitSecs(0.75);
mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],0)
mglFlush();

% trial 3

while 1
    keycode = mglGetKeys;
    if any(keycode==1), break; end
end
%%
%%
% insert

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
noise2 = coeff .* im;

grating = mglMakeGrating(stimulus.gabor.size, stimulus.gabor.size, ...
    stimulus.gabor.cycle, stimulus.gabor.tilt, 0);
contrast = 0.5;
grating = grating * contrast;
gaussian = mglMakeGaussian(stimulus.gabor.size, stimulus.gabor.size, 1, 1);
gabor = (grating.*gaussian);

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
noise_only = noise2;


% scale it to [0 255], for both stimulus images
final_im =  255 .* ((final_im + 1) ./ 2);
final_im = final_im;

noise_only = 255 .* ((noise_only + 1) ./ 2);
%%
tex = mglCreateTexture(final_im);
tex2 = mglCreateTexture(noise_only);
mglStencilSelect(2);
mglBltTexture(tex,[0,0]);
mglFlush();
mglWaitSecs(1);

mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],1)
mglFlush();
mglWaitSecs(0.5); 

mglBltTexture(tex2,[0,0]);
mglStencilSelect(0);
mglFlush();
mglWaitSecs(1);
mglClearScreen(.5)
mglTextSet([],32,1);
mglFillOval(0,0,[.2 .2],1)
mglTextDraw(['Which screen showed the target?'],[0,-8])
mglTextDraw('1(<)  or  2(>)',[0 -10]);
mglFlush();

% get response 
while 1
    k = mglGetKeys;
    if k(44)==1, mglFillOval(0,0,[.2 .2],[0 1 0]), break; end
    if k(48)==1, mglFillOval(0,0,[.2 .2],[1 0 0]), break; end
end

mglFlush();
mglWaitSecs(0.75);


%%%
mglTextDraw(['Great work so far! The real experiment will be more difficult.'], [0, 4]);
mglTextDraw(['The target will be less clear and presented for shorter.'], [0, 2]);

mglTextDraw('Now you will do 10 practice trials that will increase in difficulty.', [0, -4]);
mglTextDraw(['Press the SPACE BAR when you are ready to continue to the next instruction screen.'] ,[0,-14]);
mglFlush;

while 1
    k = mglGetKeys;
    if k(myscreen.keyboard.backtick)==1, break; end
end

%% 10 practice trials while decreasing time length and contrast

mglTextDraw('TEN EXAMPLE TRIALS WITH INCREASING DIFFICULTY' ,[0,0]);
mglFlush;
mglWaitSecs(1);

% PRESENT EXAMPLE STIMULI
stimulus.responsekeys = [44,48]; 

% trial 1 set up

mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],0)
mglFlush();

% trial 1 real
while 1
    keycode = mglGetKeys;
    if any(keycode==1), break; end
end
%%
%%
% insert

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
noise2 = coeff .* im;

grating = mglMakeGrating(stimulus.gabor.size, stimulus.gabor.size, ...
    stimulus.gabor.cycle, stimulus.gabor.tilt, 0);
contrast = 0.5;
grating = grating * contrast;
gaussian = mglMakeGaussian(stimulus.gabor.size, stimulus.gabor.size, 1, 1);
gabor = (grating.*gaussian);

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
noise_only = noise2;


% scale it to [0 255], for both stimulus images
final_im =  255 .* ((final_im + 1) ./ 2);
final_im = final_im;

noise_only = 255 .* ((noise_only + 1) ./ 2);
%%
tex = mglCreateTexture(final_im);
tex2 = mglCreateTexture(noise_only);
mglStencilSelect(2);
mglBltTexture(tex2,[0,0]);
mglFlush();
mglWaitSecs(.75);

mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],1)
mglFlush();
mglWaitSecs(0.5); 

mglBltTexture(tex,[0,0]);
mglStencilSelect(0);
mglFlush();
mglWaitSecs(0.25);
mglClearScreen(.5)
mglTextSet([],32,1);
mglFillOval(0,0,[.2 .2],1)
mglTextDraw(['Which screen showed the target?'],[0,-8])
mglTextDraw(['1(<)  or  2(>)'],[0 -10]);
mglFlush();

% get response 
while 1
    k = mglGetKeys;
    if k(44)==1, mglFillOval(0,0,[.2 .2],[1 0 0]),break; end
    if k(48)==1, mglFillOval(0,0,[.2 .2],[0 1 0]), break; end
end

mglFlush();
mglWaitSecs(0.75);
mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],1)
mglFlush();


% trial 2 - 0.75 second 
while 1
    keycode = mglGetKeys;
    if any(keycode==1), break; end
end
%%
%%
% insert

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
noise2 = coeff .* im;

grating = mglMakeGrating(stimulus.gabor.size, stimulus.gabor.size, ...
    stimulus.gabor.cycle, stimulus.gabor.tilt, 0);
contrast = 0.5;
grating = grating * contrast;
gaussian = mglMakeGaussian(stimulus.gabor.size, stimulus.gabor.size, 1, 1);
gabor = (grating.*gaussian);

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
noise_only = noise2;


% scale it to [0 255], for both stimulus images
final_im =  255 .* ((final_im + 1) ./ 2);
final_im = final_im;

noise_only = 255 .* ((noise_only + 1) ./ 2);
%%
tex = mglCreateTexture(final_im);
tex2 = mglCreateTexture(noise_only);
mglStencilSelect(2);
mglBltTexture(tex2,[0,0]);
mglFlush();
mglWaitSecs(0.75);

mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],1)
mglFlush();
mglWaitSecs(0.5); 

mglBltTexture(tex,[0,0]);
mglStencilSelect(0);
mglFlush();
mglWaitSecs(0.75);
mglClearScreen(.5)
mglTextSet([],32,1);
mglFillOval(0,0,[.2 .2],1)
mglTextDraw(['Which screen showed the target?'],[0,-8])
mglTextDraw(['1(<)  or  2(>)'],[0 -10]);
mglFlush();

% get response 
while 1
    k = mglGetKeys;
    if k(44)==1, mglFillOval(0,0,[.2 .2],[1 0 0]),break; end
    if k(48)==1, mglFillOval(0,0,[.2 .2],[0 1 0]), break; end
end

mglFlush();
mglWaitSecs(0.75);
mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],1)
mglFlush();

% trial 3 - 0.5 second 
while 1
    keycode = mglGetKeys;
    if any(keycode==1), break; end
end
%%
%%
% insert

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
noise2 = coeff .* im;

grating = mglMakeGrating(stimulus.gabor.size, stimulus.gabor.size, ...
    stimulus.gabor.cycle, stimulus.gabor.tilt, 0);
contrast = 0.5;
grating = grating * contrast;
gaussian = mglMakeGaussian(stimulus.gabor.size, stimulus.gabor.size, 1, 1);
gabor = (grating.*gaussian);

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
noise_only = noise2;


% scale it to [0 255], for both stimulus images
final_im =  255 .* ((final_im + 1) ./ 2);
final_im = final_im;

noise_only = 255 .* ((noise_only + 1) ./ 2);
%%
tex = mglCreateTexture(final_im);
tex2 = mglCreateTexture(noise_only);
mglStencilSelect(2);
mglBltTexture(tex,[0,0]);
mglFlush();
mglWaitSecs(0.75);

mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],1)
mglFlush();
mglWaitSecs(0.5); 

mglBltTexture(tex2,[0,0]);
mglStencilSelect(0);
mglFlush();
mglWaitSecs(0.75);
mglClearScreen(.5)
mglTextSet([],32,1);
mglFillOval(0,0,[.2 .2],1)
mglTextDraw(['Which screen showed the target?'],[0,-8])
mglTextDraw(['1(<)  or  2(>)'],[0 -10]);
mglFlush();

% get response 
while 1
    k = mglGetKeys;
    if k(44)==1, mglFillOval(0,0,[.2 .2],[0 1 0]),break; end
    if k(48)==1, mglFillOval(0,0,[.2 .2],[1 0 0]), break; end
end

mglFlush();
mglWaitSecs(0.75);
mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],1)
mglFlush();

% trial 4 - 0.5 second and 0.3 contrast
while 1
    keycode = mglGetKeys;
    if any(keycode==1), break; end
end
%%
%%
% insert

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
noise2 = coeff .* im;

grating = mglMakeGrating(stimulus.gabor.size, stimulus.gabor.size, ...
    stimulus.gabor.cycle, stimulus.gabor.tilt, 0);
contrast = 0.3;
grating = grating * contrast;
gaussian = mglMakeGaussian(stimulus.gabor.size, stimulus.gabor.size, 1, 1);
gabor = (grating.*gaussian);

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
noise_only = noise2;


% scale it to [0 255], for both stimulus images
final_im =  255 .* ((final_im + 1) ./ 2);
final_im = final_im;

noise_only = 255 .* ((noise_only + 1) ./ 2);
%%
tex = mglCreateTexture(final_im);
tex2 = mglCreateTexture(noise_only);
mglStencilSelect(2);
mglBltTexture(tex2,[0,0]);
mglFlush();
mglWaitSecs(0.5);

mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],1)
mglFlush();
mglWaitSecs(0.5); 

mglBltTexture(tex,[0,0]);
mglStencilSelect(0);
mglFlush();
mglWaitSecs(0.5);
mglClearScreen(.5)
mglTextSet([],32,1);
mglFillOval(0,0,[.2 .2],1)
mglTextDraw(['Which screen showed the target?'],[0,-8])
mglTextDraw(['1(<)  or  2(>)'],[0 -10]);
mglFlush();

% get response 
while 1
    k = mglGetKeys;
    if k(44)==1, mglFillOval(0,0,[.2 .2],[1 0 0]),break; end
    if k(48)==1, mglFillOval(0,0,[.2 .2],[0 1 0]), break; end
end

mglFlush();
mglWaitSecs(0.75);
mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],1)
mglFlush();

% trial 5 - 0.25 second 
while 1
    keycode = mglGetKeys;
    if any(keycode==1), break; end
end
%%
%%
% insert

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
noise2 = coeff .* im;

grating = mglMakeGrating(stimulus.gabor.size, stimulus.gabor.size, ...
    stimulus.gabor.cycle, stimulus.gabor.tilt, 0);
contrast = 0.3;
grating = grating * contrast;
gaussian = mglMakeGaussian(stimulus.gabor.size, stimulus.gabor.size, 1, 1);
gabor = (grating.*gaussian);

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
noise_only = noise2;


% scale it to [0 255], for both stimulus images
final_im =  255 .* ((final_im + 1) ./ 2);
final_im = final_im;

noise_only = 255 .* ((noise_only + 1) ./ 2);
%%
tex = mglCreateTexture(final_im);
tex2 = mglCreateTexture(noise_only);
mglStencilSelect(2);
mglBltTexture(tex,[0,0]);
mglFlush();
mglWaitSecs(0.25);

mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],1)
mglFlush();
mglWaitSecs(0.5); 

mglBltTexture(tex2,[0,0]);
mglStencilSelect(0);
mglFlush();
mglWaitSecs(0.25);
mglClearScreen(.5)
mglTextSet([],32,1);
mglFillOval(0,0,[.2 .2],1)
mglTextDraw(['Which screen showed the target?'],[0,-8])
mglTextDraw(['1(<)  or  2(>)'],[0 -10]);
mglFlush();

% get response 
while 1
    k = mglGetKeys;
    if k(44)==1, mglFillOval(0,0,[.2 .2],[0 1 0]),break; end
    if k(48)==1, mglFillOval(0,0,[.2 .2],[1 0 0]), break; end
end

mglFlush();
mglWaitSecs(0.75);
mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],1)
mglFlush();

stimulus.gabor.contrasts = .2;

% trial 6 - 0.25 second, lower contrast .2
while 1
    keycode = mglGetKeys;
    if any(keycode==1), break; end
end
%%
%%
% insert

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
noise2 = coeff .* im;

grating = mglMakeGrating(stimulus.gabor.size, stimulus.gabor.size, ...
    stimulus.gabor.cycle, stimulus.gabor.tilt, 0);
contrast = 0.2;
grating = grating * contrast;
gaussian = mglMakeGaussian(stimulus.gabor.size, stimulus.gabor.size, 1, 1);
gabor = (grating.*gaussian);

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
noise_only = noise2;


% scale it to [0 255], for both stimulus images
final_im =  255 .* ((final_im + 1) ./ 2);
final_im = final_im;

noise_only = 255 .* ((noise_only + 1) ./ 2);
%%
tex = mglCreateTexture(final_im);
tex2 = mglCreateTexture(noise_only);
mglStencilSelect(2);
mglBltTexture(tex2,[0,0]);
mglFlush();
mglWaitSecs(0.25);

mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],1)
mglFlush();
mglWaitSecs(0.5); 

mglBltTexture(tex,[0,0]);
mglStencilSelect(0);
mglFlush();
mglWaitSecs(0.25);
mglClearScreen(.5)
mglTextSet([],32,1);
mglFillOval(0,0,[.2 .2],1)
mglTextDraw(['Which screen showed the target?'],[0,-8])
mglTextDraw(['1(<)  or  2(>)'],[0 -10]);
mglFlush();

% get response 
while 1
    k = mglGetKeys;
    if k(44)==1, mglFillOval(0,0,[.2 .2],[1 0 0]),break; end
    if k(48)==1, mglFillOval(0,0,[.2 .2],[0 1 0]), break; end
end

mglFlush();
mglWaitSecs(0.75);
mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],1)
mglFlush();

stimulus.gabor.contrasts = .1;

% trial 7 - 0.25 second, lower contrast .1
while 1
    keycode = mglGetKeys;
    if any(keycode==1), break; end
end
%%
%%
% insert

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
noise2 = coeff .* im;

grating = mglMakeGrating(stimulus.gabor.size, stimulus.gabor.size, ...
    stimulus.gabor.cycle, stimulus.gabor.tilt, 0);
contrast = 0.2;
grating = grating * contrast;
gaussian = mglMakeGaussian(stimulus.gabor.size, stimulus.gabor.size, 1, 1);
gabor = (grating.*gaussian);

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
noise_only = noise2;


% scale it to [0 255], for both stimulus images
final_im =  255 .* ((final_im + 1) ./ 2);
final_im = final_im;

noise_only = 255 .* ((noise_only + 1) ./ 2);
%%
tex = mglCreateTexture(final_im);
tex2 = mglCreateTexture(noise_only);
mglStencilSelect(2);
mglBltTexture(tex2,[0,0]);
mglFlush();
mglWaitSecs(0.25);

mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],1)
mglFlush();
mglWaitSecs(0.5); 

mglBltTexture(tex,[0,0]);
mglStencilSelect(0);
mglFlush();
mglWaitSecs(0.25);
mglClearScreen(.5)
mglTextSet([],32,1);
mglFillOval(0,0,[.2 .2],1)
mglTextDraw('Which screen showed the target?',[0,-8])
mglTextDraw('1(<)  or  2(>)',[0 -10]);
mglFlush();

% get response 
while 1
    k = mglGetKeys;
    if k(44)==1, mglFillOval(0,0,[.2 .2],[1 0 0]),break; end
    if k(48)==1, mglFillOval(0,0,[.2 .2],[0 1 0]), break; end
end

mglFlush();
mglWaitSecs(0.75);
mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],1)
mglFlush();

stimulus.gabor.contrasts = .1;

% trial 8 - 0.25 second, lower contrast .15
while 1
    keycode = mglGetKeys;
    if any(keycode==1), break; end
end
%%
%%
% insert

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
noise2 = coeff .* im;

grating = mglMakeGrating(stimulus.gabor.size, stimulus.gabor.size, ...
    stimulus.gabor.cycle, stimulus.gabor.tilt, 0);
contrast = 0.15;
grating = grating * contrast;
gaussian = mglMakeGaussian(stimulus.gabor.size, stimulus.gabor.size, 1, 1);
gabor = (grating.*gaussian);

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
noise_only = noise2;


% scale it to [0 255], for both stimulus images
final_im =  255 .* ((final_im + 1) ./ 2);
final_im = final_im;

noise_only = 255 .* ((noise_only + 1) ./ 2);
%%
tex = mglCreateTexture(final_im);
tex2 = mglCreateTexture(noise_only);
mglStencilSelect(2);
mglBltTexture(tex,[0,0]);
mglFlush();
mglWaitSecs(0.25);

mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],1)
mglFlush();
mglWaitSecs(0.5); 

mglBltTexture(tex2,[0,0]);
mglStencilSelect(0);
mglFlush();
mglWaitSecs(0.25);
mglClearScreen(.5)
mglTextSet([],32,1);
mglFillOval(0,0,[.2 .2],1)
mglTextDraw(['Which screen showed the target?'],[0,-8])
mglTextDraw(['1(<)  or  2(>)'],[0 -10]);
mglFlush();

% get response 
while 1
    k = mglGetKeys;
    if k(44)==1, mglFillOval(0,0,[.2 .2],[0 1 0]),break; end
    if k(48)==1, mglFillOval(0,0,[.2 .2],[1 0 0]), break; end
end

mglFlush();
mglWaitSecs(0.75);
mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],1)
mglFlush();

stimulus.gabor.contrasts = .05;

% trial 9 - 0.25 second, lower contrast .15
while 1
    keycode = mglGetKeys;
    if any(keycode==1), break; end
end
%%
%%
% insert

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
noise2 = coeff .* im;

grating = mglMakeGrating(stimulus.gabor.size, stimulus.gabor.size, ...
    stimulus.gabor.cycle, stimulus.gabor.tilt, 0);
contrast = 0.15;
grating = grating * contrast;
gaussian = mglMakeGaussian(stimulus.gabor.size, stimulus.gabor.size, 1, 1);
gabor = (grating.*gaussian);

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
noise_only = noise2;


% scale it to [0 255], for both stimulus images
final_im =  255 .* ((final_im + 1) ./ 2);
final_im = final_im;

noise_only = 255 .* ((noise_only + 1) ./ 2);
%%
tex = mglCreateTexture(final_im);
tex2 = mglCreateTexture(noise_only);
mglStencilSelect(2);
mglBltTexture(tex2,[0,0]);
mglFlush();
mglWaitSecs(0.25);

mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],1)
mglFlush();
mglWaitSecs(0.5); 

mglBltTexture(tex,[0,0]);
mglStencilSelect(0);
mglFlush();
mglWaitSecs(0.25);
mglClearScreen(.5)
mglTextSet([],32,1);
mglFillOval(0,0,[.2 .2],1)
mglTextDraw(['Which screen showed the target?'],[0,-8])
mglTextDraw(['1(<)  or  2(>)'],[0 -10]);
mglFlush();

% get response 
while 1
    k = mglGetKeys;
    if k(44)==1, mglFillOval(0,0,[.2 .2],[1 0 0]),break; end
    if k(48)==1, mglFillOval(0,0,[.2 .2],[0 1 0]), break; end
end

mglFlush();
mglWaitSecs(0.75);
mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],1)
mglFlush();

stimulus.gabor.contrasts = .05;

% trial 10, 0.25 seconds, 0.1 contrast

while 1
    keycode = mglGetKeys;
    if any(keycode==1), break; end
end
%%
%%
% insert

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
noise2 = coeff .* im;

grating = mglMakeGrating(stimulus.gabor.size, stimulus.gabor.size, ...
    stimulus.gabor.cycle, stimulus.gabor.tilt, 0);
contrast = 0.1;
grating = grating * contrast;
gaussian = mglMakeGaussian(stimulus.gabor.size, stimulus.gabor.size, 1, 1);
gabor = (grating.*gaussian);

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
noise_only = noise2;


% scale it to [0 255], for both stimulus images
final_im =  255 .* ((final_im + 1) ./ 2);
final_im = final_im;

noise_only = 255 .* ((noise_only + 1) ./ 2);
%%
tex = mglCreateTexture(final_im);
tex2 = mglCreateTexture(noise_only);
mglStencilSelect(2);
mglBltTexture(tex,[0,0]);
mglFlush();
mglWaitSecs(.25);

mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],1)
mglFlush();
mglWaitSecs(0.5); 

mglBltTexture(tex2,[0,0]);
mglStencilSelect(0);
mglFlush();
mglWaitSecs(.25);
mglClearScreen(.5)
mglTextSet([],32,1);
mglFillOval(0,0,[.2 .2],1)
mglTextDraw(['Which screen showed the target?'],[0,-8])
mglTextDraw(['1(<)  or  2(>)'],[0 -10]);
mglFlush();

% get response 
while 1
    k = mglGetKeys;
    if k(44)==1, mglFillOval(0,0,[.2 .2],[0 1 0]),break; end
    if k(48)==1, mglFillOval(0,0,[.2 .2],[1 0 0]), break; end
end

mglFlush();
mglWaitSecs(0.75);


%%%
mglTextDraw(['The target has only appeared at the center so far. However, in the real'] ,[0,6]);
mglTextDraw(['task the target will appear in other positions inside the circle.'] ,[0,4]);
mglTextDraw(['We will let you know in advance where the target will appear in a new location.'] ,[0,0]);
mglTextDraw(['Please keep your eyes at the center even if the target appears somewhere else.'] ,[0,-4]);
mglTextDraw(['For your final practice, we will present 5 trials where the target appears in new positions.'] ,[0,-6]);
mglTextDraw(['Press the SPACE BAR when you are ready to continue to the next instruction screen.'] ,[0,-14]);
mglFlush;

while 1
    k = mglGetKeys;
    if k(myscreen.keyboard.backtick)==1, break; end
end
%% five practice at different locations
% screen 1 - use black dot to indicate loc
% screen 2 and 3 - target trials

mglTextDraw(['FIVE EXAMPLE TRIALS AT DIFFERENT LOCATIONS'] ,[0,0]);
mglFlush;
mglWaitSecs(1);

% PRESENT EXAMPLE STIMULI
stimulus.responsekeys = [44,48]; 

% create new magnitude
new_mag = fwhite .* pink_filter;
new_Fourier = new_mag .* (cos(phase) + sqrt(-1)*sin(phase));
im = ifft2(ifftshift(new_Fourier));

N = length(im(:));
m_im = mean(im(:));
coeff = sqrt((N*stimulus.noise.contrasts^2) / sum((im(:)-m_im).^2));
noiseStable = coeff .* im;

white = randn(noise_frame_pixel, noise_frame_pixel);
fwhite = fftshift(fft2(white));
phase = angle(fwhite);

noiseStable = 255 .* ((noiseStable + 1) ./ 2);
%%
texStable = mglCreateTexture(noiseStable);

% trial 1 loc 

mglClearScreen(.5)
mglTextSet([],32,1);
mglTextDraw(['Target location for this task'], [0, 12])
mglStencilSelect(2);
mglBltTexture(texStable,[0,0]);
mglGluAnnulus(0, -3.5, .35, .4, ...
    [1 1 1], 120, 2)
mglFlush;
mglWaitSecs(2);

% trial 1 set up - loc 1 0.25 seconds

mglClearScreen(.5)
mglStencilSelect(0);
mglFillOval(0,0,[.2 .2],0)
mglFlush();

% actual trial 1
while 1
    keycode = mglGetKeys;
    if any(keycode==1), break; end
end
%%
%%
% insert

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
noise2 = coeff .* im;

grating = mglMakeGrating(stimulus.gabor.size, stimulus.gabor.size, ...
    stimulus.gabor.cycle, stimulus.gabor.tilt, 0);
contrast = 0.15;
grating = grating * contrast;
gaussian = mglMakeGaussian(stimulus.gabor.size, stimulus.gabor.size, 1, 1);
gabor = (grating.*gaussian);

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
noise_only = noise2;


% scale it to [0 255], for both stimulus images
final_im =  255 .* ((final_im + 1) ./ 2);
final_im = final_im;

noise_only = 255 .* ((noise_only + 1) ./ 2);
%%
tex = mglCreateTexture(final_im);
tex2 = mglCreateTexture(noise_only);
mglStencilSelect(2);
mglBltTexture(tex2,[0,0]);
mglFlush();
mglWaitSecs(0.25);

mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],1)
mglFlush();
mglWaitSecs(0.5); 

mglBltTexture(tex,[0,-3.5]);
mglStencilSelect(0);
mglFlush();
mglWaitSecs(0.25);
mglClearScreen(.5)
mglTextSet([],32,1);
mglFillOval(0,0,[.2 .2],1)
mglTextDraw(['Which screen showed the target?'],[0,-8])
mglTextDraw(['1(<)  or  2(>)'],[0 -10]);
mglFlush();

% get response 
while 1
    k = mglGetKeys;
    if k(44)==1, mglFillOval(0,0,[.2 .2],[1 0 0]),break; end
    if k(48)==1, mglFillOval(0,0,[.2 .2],[0 1 0]), break; end
end

mglFlush();
mglWaitSecs(0.75);

% trial 2 loc 

mglClearScreen(.5)
mglTextSet([],32,1);
mglTextDraw(['Target location for this task'], [0, 12])
mglStencilSelect(2);
mglBltTexture(texStable,[0,0]);
mglGluAnnulus(2, 1.75, .35, .4, ...
    [1 1 1], 120, 2)
mglFlush;
mglWaitSecs(2);
mglClearScreen(.5)
mglStencilSelect(0);
mglFillOval(0,0,[.2 .2],1)
mglFlush();

% trial 2 - loc 2 0.25 seconds

while 1
    keycode = mglGetKeys;
    if any(keycode==1), break; end
end
%%
% insert

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
noise2 = coeff .* im;

grating = mglMakeGrating(stimulus.gabor.size, stimulus.gabor.size, ...
    stimulus.gabor.cycle, stimulus.gabor.tilt, 0);
contrast = 0.15;
grating = grating * contrast;
gaussian = mglMakeGaussian(stimulus.gabor.size, stimulus.gabor.size, 1, 1);
gabor = (grating.*gaussian);

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
noise_only = noise2;


% scale it to [0 255], for both stimulus images
final_im =  255 .* ((final_im + 1) ./ 2);
final_im = final_im;

noise_only = 255 .* ((noise_only + 1) ./ 2);
%%
tex = mglCreateTexture(final_im);
tex2 = mglCreateTexture(noise_only);
mglStencilSelect(2);
mglBltTexture(tex,[2,1.75]);
mglFlush();
mglWaitSecs(0.25);

mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],1)
mglFlush();
mglWaitSecs(0.5); 

mglBltTexture(tex2,[0,0]);
mglStencilSelect(0);
mglFlush();
mglWaitSecs(0.25);
mglClearScreen(.5)
mglTextSet([],32,1);
mglFillOval(0,0,[.2 .2],1)
mglTextDraw(['Which screen showed the target?'],[0,-8])
mglTextDraw(['1(<)  or  2(>)'],[0 -10]);
mglFlush();

% get response 
while 1
    k = mglGetKeys;
    if k(44)==1, mglFillOval(0,0,[.2 .2],[0 1 0]),break; end
    if k(48)==1, mglFillOval(0,0,[.2 .2],[1 0 0]), break; end
end

mglFlush();
mglWaitSecs(0.75);

% trial 3 loc 

mglClearScreen(.5)
mglTextSet([],32,1);
mglTextDraw(['Target location for this task'], [0, 12])
mglStencilSelect(2);
mglBltTexture(texStable,[0,0]);
mglGluAnnulus(.22, -.22, .35, .4, ...
    [1 1 1], 120, 2)
mglFlush;
mglWaitSecs(2);

mglClearScreen(.5)
mglStencilSelect(0);
mglFillOval(0,0,[.2 .2],1)
mglFlush();
 


% trial 3 - loc 3 0.25 seconds

while 1
    keycode = mglGetKeys;
    if any(keycode==1), break; end
end
%%
% insert

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
noise2 = coeff .* im;

grating = mglMakeGrating(stimulus.gabor.size, stimulus.gabor.size, ...
    stimulus.gabor.cycle, stimulus.gabor.tilt, 0);
contrast = 0.15;
grating = grating * contrast;
gaussian = mglMakeGaussian(stimulus.gabor.size, stimulus.gabor.size, 1, 1);
gabor = (grating.*gaussian);

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
noise_only = noise2;


% scale it to [0 255], for both stimulus images
final_im =  255 .* ((final_im + 1) ./ 2);
final_im = final_im;

noise_only = 255 .* ((noise_only + 1) ./ 2);
%%
tex = mglCreateTexture(final_im);
tex2 = mglCreateTexture(noise_only);
mglStencilSelect(2);
mglBltTexture(tex,[.22,-.22]);
mglFlush();
mglWaitSecs(0.25);

mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],1)
mglFlush();
mglWaitSecs(0.5); 

mglBltTexture(tex2,[0,0]);
mglStencilSelect(0);
mglFlush();
mglWaitSecs(0.25);
mglClearScreen(.5)
mglTextSet([],32,1);
mglFillOval(0,0,[.2 .2],1)
mglTextDraw(['Which screen showed the target?'],[0,-8])
mglTextDraw(['1(<)  or  2(>)'],[0 -10]);
mglFlush();

% get response 
while 1
    k = mglGetKeys;
    if k(44)==1, mglFillOval(0,0,[.2 .2],[0 1 0]),break; end
    if k(48)==1, mglFillOval(0,0,[.2 .2],[1 0 0]), break; end
end

mglFlush();
mglWaitSecs(0.75);

% trial 4 loc 

mglClearScreen(.5)
mglTextSet([],32,1);
mglTextDraw(['Target location for this task'], [0, 12])
mglStencilSelect(2);
mglBltTexture(texStable,[0,0]);
mglGluAnnulus(-2, 2, .35, .4, ...
    [1 1 1], 120, 2)
mglFlush;
mglWaitSecs(2);

mglClearScreen(.5)
mglStencilSelect(0);
mglFillOval(0,0,[.2 .2],1)
mglFlush();

% trial 4 - loc 4 0.25 seconds

while 1
    keycode = mglGetKeys;
    if any(keycode==1), break; end
end
%%
% insert

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
noise2 = coeff .* im;

grating = mglMakeGrating(stimulus.gabor.size, stimulus.gabor.size, ...
    stimulus.gabor.cycle, stimulus.gabor.tilt, 0);
contrast = 0.15;
grating = grating * contrast;
gaussian = mglMakeGaussian(stimulus.gabor.size, stimulus.gabor.size, 1, 1);
gabor = (grating.*gaussian);

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
noise_only = noise2;


% scale it to [0 255], for both stimulus images
final_im =  255 .* ((final_im + 1) ./ 2);
final_im = final_im;

noise_only = 255 .* ((noise_only + 1) ./ 2);
%%
tex = mglCreateTexture(final_im);
tex2 = mglCreateTexture(noise_only);
mglStencilSelect(2);
mglBltTexture(tex,[-2,2]);
mglFlush();
mglWaitSecs(0.25);

mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],1)
mglFlush();
mglWaitSecs(0.5); 

mglBltTexture(tex2,[0,0]);
mglStencilSelect(0);
mglFlush();
mglWaitSecs(0.25);
mglClearScreen(.5)
mglTextSet([],32,1);
mglFillOval(0,0,[.2 .2],1)
mglTextDraw('Which screen showed the target?',[0,-8])
mglTextDraw('1(<)  or  2(>)',[0 -10]);
mglFlush();

% get response 
while 1
    k = mglGetKeys;
    if k(44)==1, mglFillOval(0,0,[.2 .2],[0 1 0]),break; end
    if k(48)==1, mglFillOval(0,0,[.2 .2],[1 0 0]), break; end
end

mglFlush();
mglWaitSecs(0.75);

% trial 5 loc 

mglClearScreen(.5)
mglTextSet([],32,1);
mglTextDraw(['Target location for this task'], [0, 12])
mglStencilSelect(2);
mglBltTexture(texStable,[0,0]);
mglGluAnnulus(1, -1, .35, .4, ...
    [1 1 1], 120, 2)
mglFlush;
mglWaitSecs(2);

mglClearScreen(.5)
mglStencilSelect(0);
mglFillOval(0,0,[.2 .2],1)
mglFlush();
 
% trial 5 - loc 5 0.25 seconds

while 1
    keycode = mglGetKeys;
    if any(keycode==1), break; end
end
%%
% insert

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
noise2 = coeff .* im;

grating = mglMakeGrating(stimulus.gabor.size, stimulus.gabor.size, ...
    stimulus.gabor.cycle, stimulus.gabor.tilt, 0);
contrast = 0.15;
grating = grating * contrast;
gaussian = mglMakeGaussian(stimulus.gabor.size, stimulus.gabor.size, 1, 1);
gabor = (grating.*gaussian);

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
noise_only = noise2;


% scale it to [0 255], for both stimulus images
final_im =  255 .* ((final_im + 1) ./ 2);
final_im = final_im;

noise_only = 255 .* ((noise_only + 1) ./ 2);
%%
tex = mglCreateTexture(final_im);
tex2 = mglCreateTexture(noise_only);
mglStencilSelect(2);
mglBltTexture(tex,[1,-1]);
mglFlush();
mglWaitSecs(0.25);

mglClearScreen(.5)
mglFillOval(0,0,[.2 .2],1)
mglFlush();
mglWaitSecs(0.5); 

mglBltTexture(tex2,[0,0]);
mglStencilSelect(0);
mglFlush();
mglWaitSecs(0.25);
mglClearScreen(.5)
mglTextSet([],32,1);
mglFillOval(0,0,[.2 .2],1)
mglTextDraw('Which screen showed the target?',[0,-8])
mglTextDraw('1(<)  or  2(>)',[0 -10]);
mglFlush();

% get response 
while 1
    k = mglGetKeys;
    if k(44)==1, mglFillOval(0,0,[.2 .2],[0 1 0]),break; end
    if k(48)==1, mglFillOval(0,0,[.2 .2],[1 0 0]), break; end
end
 
mglFlush;
mglWaitSecs(0.75);

%

mglTextDraw('Now we will move on to the real experiment. Please call the experimenter.' ,[0,0]);
mglTextDraw('If you have any questions, please ask the experimenter now.' ,[0,-2]);
mglFlush;
%% END INSTRUCTIONS.

%% present screen with double visual display sample code [OLD]

%mglClearScreen(.5);
%tex = mglCreateTexture(final_im);
%tex2 = mglCreateTexture(noise_only);
%mglStencilSelect(1);
%mglBltTexture([tex tex2],[-stimulus.noise.size/2,0; stimulus.noise.size/2 + 1,0]);
%mglStencilSelect(0);

%mglFlush();
%% helper function for pink noise
%function pink_filter = createPinkFilter(myscreen)
%w = myscreen.screenWidth;
%h = myscreen.screenHeight;
%sz = max(w,h);
        
% make the odd size of the image
%if mod(sz,2)==0, sz = sz-1; end

% make pink filter
%last_freq = ceil(sz/2);
%pink_filter = zeros(sz,sz);
%[x, y] = meshgrid(-ceil(sz/2)+1:ceil(sz/2)-1, -ceil(sz/2)+1:ceil(sz/2)-1);
%index = sqrt(x.^2 + y.^2);
%for f = 1:last_freq
    %pink_filter(index > f-1 & index < f+1) = 1/f;
%end

function createPinkNoise(myscreen)
global stimulus
noise_with_buffer = stimulus.noise.size + 3;    % visual angle
noise_frame_pixel = visualAngleToPixels(noise_with_buffer, ...
    [myscreen.screenWidth, myscreen.screenHeight]);
% make the size of the image an odd number
if mod(noise_frame_pixel,2)==0, noise_frame_pixel = noise_frame_pixel+1; end    

filter_sz = size(stimulus.pink_filter);
pink_filter = stimulus.pink_filter(...
    floor(filter_sz(1)/2)+1-(noise_frame_pixel-1)/2:floor(filter_sz(1)/2)+1+(noise_frame_pixel-1)/2, ...
    floor(filter_sz(2)/2)+1-(noise_frame_pixel-1)/2:floor(filter_sz(2)/2)+1+(noise_frame_pixel-1)/2);

for images = 1:2    % create two noise images
    % fft on white noise
    white = randn(noise_frame_pixel, noise_frame_pixel);
    fwhite = fftshift(fft2(white));
    phase = angle(fwhite);
    
    % create new magnitude
    new_mag = fwhite .* pink_filter;
    new_Fourier = new_mag .* (cos(phase) + sqrt(-1)*sin(phase));
    im = ifft2(ifftshift(new_Fourier));
    
    % change contrast
    contrast = stimulus.noise_contrast;
    N = length(im(:));
    m_im = mean(im(:));
    coeff = sqrt((N*contrast^2) / sum((im(:)-m_im).^2));
    stimulus.noise.im{images} = coeff .* im;
end
