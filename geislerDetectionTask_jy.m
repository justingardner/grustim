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

function geislerDetectionTask_jy
mglClose        % close MGL if it's open
clear all, close all, clc
global stimulus

myscreen.screenNumber = 2;
myscreen.responseKeys = {'50'}; % respond only with the space bar
myscreen = initScreen(myscreen);

% load pink_filter
if exist([cd '/geislerDetectionTask_pinkFilter.mat']) ~= 0
    load('geislerDetectionTask_pinkFilter.mat');
    stimulus.pink_filter = pink_filter;
else
    createPinkFilter(myscreen);
end

%%%%% define task timings and responses
task{1}.waitForBacktick = 1;
task{1}.segmin = [inf, .25, .5, .25, inf];  
task{1}.segmax = [inf, .25, .5, .25, inf];  
%  = button press - stim1(250ms) - int(500ms) - stim2(250ms) - response
task{1}.getResponse = [1 0 0 0 1];
stimulus.nBlocks = 2;
stimulus.TrialsPerBlock = 2;
task{1}.numTrials = stimulus.nBlocks * stimulus.TrialsPerBlock;

%%%%% set stimulus parameter
stimulus.responseKeys = [50];   % space bar
stimulus.noise.size = 15;   % visual angle
stimulus.noise.contrasts = [0, .05, .10, .20];
stimulus.gabor.size = 2;    % visual angle
stimulus.gabor.tilt = 315;
stimulus.gabor.cycle = 6;
stimulus.gabor.nLoc = 25;   % 25 for the real experiment
stimulus.gabor.contrasts = [.2, .1, .075, .05];
stimulus.contrast_combinations = [1: ...
    length(stimulus.noise.contrasts) * length(stimulus.gabor.contrasts)];
stimulus.nPossibleContrasts = length(stimulus.contrast_combinations);
defineLocations;


%%%%% things to be randomized 
task{1}.randVars.uniform.whichseg = [2 4];    % at which segment to present the stimulus
task{1}.randVars.calculated.noise_contrast = nan;  
task{1}.randVars.calculated.gabor_contrast = nan;
task{1}.randVars.calculated.gabor_location = [nan, nan];    % start with a random position

%%%%% initialize stimulus
myscreen = initStimulus('stimulus', myscreen);

%%%%% make stencil
mglClearScreen(.5);
mglStencilCreateBegin(1);
mglVisualAngleCoordinates(myscreen.displayDistance,myscreen.displaySize);
mglFillOval(0, 0, [stimulus.noise.size, stimulus.noise.size]);
mglStencilCreateEnd;
mglClearScreen(.5);

%%%%% initialize task
[task{1} myscreen] = initTask(task{1},myscreen,...
    @startSegmentCallback, @updateScreenCallback, @getResponseCallback, ...
    @startTrialCallback);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while task{1}.trialnum <= task{1}.numTrials ...
        || ~myscreen.userHitEsc
    % update the task
    [task myscreen] = updateTask(task,myscreen,1);
    % flip the screen
    myscreen = tickScreen(myscreen, task);
end
% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% startTrialCallback
%   prepare the noise and stimulus images to present
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startTrialCallback(task, myscreen)
global stimulus

% for every new block, update stimulus 
if mod(task.trialnum, stimulus.TrialsPerBlock) == 1
    % decide on contrasts
    stimulus.currentContrasts = randsample(stimulus.contrast_combinations,1);
    stimulus.contrast_combinations(stimulus.contrast_combinations == stimulus.currentContrasts) = [];
    
    % decide on location to present the gabor
    stimulus.current_gabor_location = randsample(stimulus.gabor.nLoc,1);
end
% decide on contrasts
index = reshape(1:stimulus.nPossibleContrasts, ...
    length(stimulus.noise.contrasts), length(stimulus.gabor.contrasts));
[noise_contrast, gabor_contrast] = find(index == stimulus.currentContrasts);
task.thistrial.noise_contrast = stimulus.noise.contrasts(noise_contrast);
task.thistrial.gabor_contrast = stimulus.gabor.contrasts(gabor_contrast);

% decide on location
task.thistrial.gabor_location = stimulus.gabor_locations(stimulus.current_gabor_location,:);

% generate noise images
task.thistrial.noise_contrast
createPinkNoise(myscreen, task);

% generate gabor
createGabor(task);

% combine noise and gabor
combinedStimulus(task);

% create texture in advance
stimulus.tex_target = mglCreateTexture(stimulus.final_im{1});
stimulus.tex_nontarget = mglCreateTexture(stimulus.final_im{2});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)
global stimulus

if task.thistrial.thisseg == 1
    % show a fixation cross and wait for the button press
    mglClearScreen(.5)
    mglFixationCross(1,1,255,[0 0])    
    
elseif task.thistrial.thisseg == 3
    % present a screen with a black dot
    mglClearScreen(.5)
    mglFillOval(0,0,[.2 .2],0)
    
elseif task.thistrial.thisseg == 5
    % present a screen with a white dot and wait for the response
    mglClearScreen(.5)
    mglFillOval(0,0,[.2 .2],255)

elseif task.thistrial.thisseg == task.thistrial.whichseg
    % present noise with gabor stimulus
    mglClearScreen(stimulus.bg_color{1});
    mglStencilSelect(1);
    mglBltTexture(stimulus.tex_target,[0 0])
    mglStencilSelect(0);
    
else
    % present noise only screen
    mglClearScreen(stimulus.bg_color{2});
    mglStencilSelect(1);    
    mglBltTexture(stimulus.tex_nontarget,[0 0]);
    mglStencilSelect(0);
    
end
mglFlush


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task, myscreen)
%%%%% this function is left empty since there's no component to be updated
%%%%% by framewise


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to get responses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = getResponseCallback(task, myscreen)
global stimulus

if task.thistrial.thisseg == 1
%     while 
%         k = mglGetKeys;
%         if (
%         continue;
%     end
    
    % remove the button responses to get response for the last segment
%     task.thistrial.whichButton = [];

    % move to the next segment
    task = jumpSegment(task);
    
elseif task.thistrial.thisseg == 5
    % get the response 
    while ~any(stimulus.responseKeys == task.thistrial.whichButton)
        continue;
    end
    
    % start a new trial
    task = jumpSegment(task,inf);
end


%%%%%%%%%% helper functions
function createPinkFilter(myscreen)
global stimulus
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
stimulus.pink_filter = pink_filter;

function createPinkNoise(myscreen, task)
global stimulus
w = myscreen.screenWidth;
h = myscreen.screenHeight;
sz = max(w,h);

% make the odd size of the image
if mod(sz,2)==0, sz = sz-1; end

for images = 1:2    % create two noise images
    % fft on white noise
    white = randn(sz,sz);
    fwhite = fftshift(fft2(white));
    phase = angle(fwhite);
    
    % create new magnitude
    new_mag = fwhite .* stimulus.pink_filter;
    new_Fourier = new_mag .* (cos(phase) + sqrt(-1)*sin(phase));
    im = ifft2(ifftshift(new_Fourier));
    
    % change contrast
    contrast = task.thistrial.noise_contrast;
    N = length(im(:));
    m_im = mean(im(:));
    coeff = sqrt((N*contrast^2) / sum((im(:)-m_im).^2));
    stimulus.noise.im{images} = coeff .* im;
end

function createGabor(task)
global stimulus
grating = mglMakeGrating(stimulus.gabor.size, stimulus.gabor.size, ...
    stimulus.gabor.cycle, stimulus.gabor.tilt, 0);
contrast = task.thistrial.gabor_contrast;
grating = grating * contrast;
gaussian = mglMakeGaussian(stimulus.gabor.size, stimulus.gabor.size, 1, 1);
stimulus.gabor.im = (grating.*gaussian)/2;

function defineLocations
global stimulus
% determine how many layers to have
% maximum number of locations per layer is 8
nLoc = stimulus.gabor.nLoc;
if nLoc > 8 && mod(nLoc,8) ~= 0 
    if mod(nLoc,8) < 4
        nLoc = floor(nLoc/8) * 8 + 1;
    else 
        nLoc = ceil(nLoc/8) * 8 + 1;
    end
    stimulus.gabor.nLoc = nLoc;
    nLayer = floor(nLoc/8);
elseif nLoc <= 8
    nLayer = 1;
end

radius_va = linspace(0, stimulus.noise.size/2+1, nLayer+2);     % radius in visual angle
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
sz = max(size(stimulus.pink_filter));
pix_noise_height = sz;
pix_noise_width = sz;
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

stimulus.gabor_locations = locations;

function combinedStimulus(task)
global stimulus
noise = stimulus.noise.im{1};
gabor = stimulus.gabor.im;
location = task.thistrial.gabor_location;   % gabor's center

% make circular gabor patch
radius_va = stimulus.gabor.size/2;
radius_cm = 2 .* mglGetParam('devicePhysicalDistance') * tan(radius_va /2 * (pi/180));
radius_pix = radius_cm * min([mglGetParam('xDeviceToPixels'), mglGetParam('yDeviceToPixels')]);
radius = round(radius_pix);
[stencil_x, stencil_y] = meshgrid(-ceil(size(gabor,1)/2)+1:ceil(size(gabor,1)/2)-1, ...
    -ceil(size(gabor,2)/2)+1:ceil(size(gabor,2)/2)-1);
stencil = (sqrt(stencil_x.^2 + stencil_y.^2) <= radius);
gabor_circle = stencil' .* gabor;

% determine the location to display
x_lims = [location(1)-ceil(size(gabor_circle,1)/2)+1, location(1)+ceil(size(gabor_circle,1)/2)-1];
y_lims = [location(2)-ceil(size(gabor_circle,2)/2)+1, location(2)+ceil(size(gabor_circle,2)/2)-1];

gabor_position = zeros(size(noise,1), size(noise,2));
gabor_position(x_lims(1):x_lims(2),y_lims(1):y_lims(2)) = gabor_circle;

% add gabor to the noise
final_im = noise + gabor_position;

% scale it to [0 255], for both stimulus images
final_im =  255 .* ((final_im + 1) ./ 2);
stimulus.final_im{1} = final_im';
stimulus.final_im{2} = 255 .* (stimulus.noise.im{2} + 1 ./ 2);

% decide background color
for image = 1:2
    if task.thistrial.noise_contrast == 0
        bg_color = stimulus.final_im{image}(1,1);
    else
        bg_color = mean(stimulus.final_im{image});
    end
    stimulus.bg_color{image} = bg_color;
end


