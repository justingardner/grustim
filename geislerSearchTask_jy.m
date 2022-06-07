% geislerSearchTask_jy.m
% 
%      usage: geislerSearchTask_jy.m
%         by: jiwon yeon
%       date: 
%  copyright: (c) 2022 Jiwon Yeon
%    purpose: Replicating Najemnik&Geisler's 2005 study
%
%             The search task. To properly run this task, subject must have
%             completed the detection task (to set the target's contrast 
%             level) A stimulus screen would be presented briefly and then 
%             subjects have to indicate with the mouse where the target 
%             appeared.
%

function geislerSearchTask_jy
mglClose        % close MGL if it's open
clear all, close all, clc
global stimulus

myscreen.screenNumber = 2;
myscreen.saveData = 0;
myscreen.eyetracker = 0;
myscreen.datadir = '~/proj/jiwon/data/geisler/geislerSearchTask';
mglSetParam('abortedStimfilesDir', '~/proj/data/geislerSearchTask/aborted',1);

myscreen.keyboard.nums = [50]; % respond only with the space bar
myscreen = initScreen(myscreen);

% load pink_filter
if exist([cd '/geislerDetectionTask_pinkFilter.mat']) ~= 0
    load('geislerDetectionTask_pinkFilter.mat');
    stimulus.pink_filter = pink_filter;
else
    createPinkFilter(myscreen);
end
stimulus.nBlocks = 2;   % 6 blocks
stimulus.cBlock = 0;    % current block
stimulus.TrialsPerBlock = 2;    % 32 trials per block

%%%%% set stimulus parameter
stimulus.responsekeys = [50];   % space bar
stimulus.noise.size = 15;   % visual angle
stimulus.noise.contrasts = .2;   % two levels of noise contrasts

stimulus.gabor.size = 1;    % visual angle
stimulus.gabor.tilt = 315;
stimulus.gabor.cycle = 6;
% 6 levels of target contrasts, that computed from the dprime
% d' = [3, 3.5, 4, 5, 6, 7];
stimulus.gabor.contrasts = 1;
defineLocations

%%%%% define task timings and responses
task{1}{1}.waitForBacktick = 0;
task{1}{1}.segmin = [inf, .1, inf, inf ];  
task{1}{1}.segmax = [inf, .5, inf, inf];
        %  fixation-int-search-response (need feedback?)
task{1}{1}.getResponse = [1 0 1 0 ];
task{1}{1}.numTrials = stimulus.nBlocks * stimulus.TrialsPerBlock;

%%%%% things to be randomized or to be saved
task{1}{1}.randVars.block.gabor_location = 1:stimulus.gabor.nLoc;
task{1}{1}.randVars.block.gabor_contrast = stimulus.gabor.contrasts; 
task{1}{1}.randVars.calculated.detection_rt = nan;
task{1}{1}.randVars.calculated.decision_rt = nan;
task{1}{1}.randVars.calculated.mousePos = [nan nan];
task{1}{1}.randVars.calculated.response_offset = [nan nan];

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
% Eye calibration (optional)
if ~myscreen.eyetracker
    disp(' Calibrating Eye ....')
    myscreen = eyeCalibDisp(myscreen);
end

% hide cursor 
mglDisplayCursor(0)

% notify the starting of the task
mglClearScreen(.5);
mglTextSet([],32,1);
mglTextDraw('Starting the experiment',[0,0]);
mglFlush;
while 1
    k = mglGetKeys;
    if k(myscreen.keyboard.backtick)==1, break; end
end

% main loop
while (task{1}{1}.trialnum <= task{1}{1}.numTrials) && ~myscreen.userHitEsc
    % update the task
    [task myscreen] = updateTask(task,myscreen,1);
    % flip the screen
    myscreen = tickScreen(myscreen, task);
end

% ends the task
mglClearScreen(0.5);
mglTextSet([],32,1);
mglTextDraw('Experiment ends',[0, .7]);
mglTextDraw('Please wait..', [0, -.7]);
mglFlush

% if we got here, we are at the end of the experiment
mglWaitSecs(3);
myscreen = endTask(myscreen,task);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% startTrialCallback
%   prepare the noise and stimulus images to present
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startTrialCallback(task, myscreen)
global stimulus

% for every new block, update stimulus 
if mod(task.trialnum, stimulus.TrialsPerBlock) == 1
    % update current block
    stimulus.cBlock = stimulus.cBlock + 1;
end

% generate noise images
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
    % show how many blocks are left 
    if mod(task.trialnum, stimulus.TrialsPerBlock) == 1
        mglClearScreen(.5)
        mglTextSet([],32,1);
        mglTextDraw(sprintf('Starting block %d out of %d blocks', ...
            stimulus.cBlock, stimulus.nBlocks),[0,0])
        mglFlush
        mglWaitSecs(2)
    end
    
    % show a fixation cross and wait for the button press
    mglClearScreen(.5)
    mglFillOval(0,0,[.2 .2],0)    
    mglFlush
    
elseif task.thistrial.thisseg == 2 
    % present an empty screen
    mglClearScreen(.5)    
    mglFlush
    
elseif task.thistrial.thisseg == 3
    % present a stimulus
    mglClearScreen(stimulus.bg_color);
    mglStencilSelect(1);
    mglBltTexture(stimulus.tex_target,[0 0])
    mglStencilSelect(0);    
    mglFlush
    
    disp(['current target:' task.thistrial.gabor_location])
    
    target_location = stimulus.gabor_locations_va(task.thistrial.gabor_location,:);    
    mglGluAnnulus(target_location(2), target_location(1), .35, .4, ...
        [1 1 1], 120, 2)
    mglTextSet([],20,1);    
    mglTextDraw(num2str(task.thistrial.gabor_location), ...
        [target_location(task.thistrial.gabor_location,2) target_location(task.thistrial.gabor_location,1)]);
    
elseif task.thistrial.thisseg == 4
    % show mouse cursor at the initial location
    mglSetMousePosition(myscreen.screenWidth/2, myscreen.screenHeight/2, ...
        myscreen.screenNumber)
    mglDisplayCursor(1)

    % decision prompt
    mglClearScreen(stimulus.bg_color);
    mglStencilSelect(1);
    mglBltTexture(stimulus.tex_nontarget,[0 0])
    mglStencilSelect(0);
    mglTextDraw('Click on the screen where the target appeared', [0,10])   
    mglFlush
    
    % start response time recording
    stimulus.t0 = mglGetSecs;
    
    % start recording mouse positions
    mInfo = mglGetMouse(myscreen.screenNumber);
    mousePos(1,:) = [mInfo.x, mInfo.y];
    
    % keep recording until responding
    while 1
        mInfo = mglGetMouse(myscreen.screenNumber);
        mousePos(end+1,:) = [mInfo.x, mInfo.y];
        if mInfo.buttons
            task.thistrial.decision_rt = mglGetSecs(stimulus.t0);
            mglDisplayCursor(0)
            break            
        end        
    end
    
    % save response info
    task.thistrial.mousePos = mousePos;     % is this in va or pixels?
    task.thistrial.response_offset = task.thistrial.gabor_location - [mInfo.x, mInfo.y];
    
    task = jumpSegment(task);    
% elseif task.thistrial.thisseg == 5
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task, myscreen)
%%% screen doesn't have to be updated
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to get responses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = getResponseCallback(task, myscreen)
global stimulus

if task.thistrial.thisseg == 1
    % waiting for the subject to start the trial
    while 1
        keycode = mglGetKeys;
        if any(keycode(stimulus.responsekeys)==1)
            break
        end
    end
    % move to the next segment
    task = jumpSegment(task);
    
elseif task.thistrial.thisseg == 3
    % get detection response
    while 1
        keycode = mglGetKeys;
        if any(keycode(stimulus.responsekeys)==1)
            task.thistrial.detection_rt = task.thistrial.reactionTime;
            break
        end
    end
    % move to the next segment
    task = jumpSegment(task);
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

function createGabor(task)
global stimulus
grating = mglMakeGrating(stimulus.gabor.size, stimulus.gabor.size, ...
    stimulus.gabor.cycle, stimulus.gabor.tilt, 0);
contrast = task.thistrial.gabor_contrast;
grating = grating * contrast;
gaussian = mglMakeGaussian(stimulus.gabor.size, stimulus.gabor.size, 1, 1);
stimulus.gabor.im = (grating.*gaussian);

function defineLocations
global stimulus
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
locations_va = trigrid(inside,:);

% convert visual angle of the locations to pixels
displaySize = max(size(stimulus.pink_filter));
locations = visualAngleToPixels(locations_va, displaySize);

stimulus.gabor_locations_va = locations_va;
stimulus.gabor_locations = locations;   % in pixels
stimulus.gabor.nLoc = size(locations,1);

function combinedStimulus(task)
global stimulus
noise = stimulus.noise.im;
gabor = stimulus.gabor.im;
location = stimulus.gabor_locations(task.thistrial.gabor_location,:);   % gabor's center

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
x_lims = [location(1)-floor(size(gabor_circle,1)/2), location(1)+floor(size(gabor_circle,1)/2)];
y_lims = [location(2)-floor(size(gabor_circle,2)/2), location(2)+floor(size(gabor_circle,2)/2)];

gabor_position = zeros(size(noise,1), size(noise,2));
gabor_position(x_lims(1):x_lims(2),y_lims(1):y_lims(2)) = gabor_circle;

% add gabor to the noise
final_im = noise + gabor_position;

% scale it to [0 255], for both stimulus images
final_im =  255 .* ((final_im + 1) ./ 2);
stimulus.final_im{1} = final_im;
stimulus.final_im{2} = 255 .* ((noise + 1) ./ 2);

% decide background color
bg_color = mean(stimulus.final_im{2}(:));
stimulus.bg_color = bg_color;


