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

%%%%
%%%% need to change how to decide contrast and noise levels
%%%%

function geislerSearchTask_jy
mglClose        % close MGL if it's open
clear all, close all, clc
global stimulus

myscreen.screenNumber = 2;
myscreen.saveData = 1;
myscreen.datadir = '~/proj/data/geislerSearchTask';
myscreen.eyetracker = 0;
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

%%%%% define task timings and responses
task{1}.waitForBacktick = 1;
task{1}.segmin = [inf, .1, inf, inf];  
task{1}.segmax = [inf, .5, inf, inf];
%  fixation-int-search-response

task{1}.getResponse = [1 0 1 0];
stimulus.nBlocks = 1;
stimulus.cBlock = 0;    % current block
stimulus.TrialsPerBlock = 2;
task{1}.numTrials = stimulus.nBlocks * stimulus.TrialsPerBlock;

%%%%% set stimulus parameter
stimulus.responsekeys = [50];   % space bar
stimulus.noise.size = 15;   % visual angle
stimulus.noise.contrasts = [.05, .2];   % two levels of noise contrasts
stimulus.noise.contrasts = .05;

stimulus.gabor.size = 1;    % visual angle
stimulus.gabor.tilt = 315;
stimulus.gabor.cycle = 6;
stimulus.gabor.nLoc = 85;   % 85 for the real experiment

% 6 levels of target contrasts, that computed from the dprime
% d' = [3, 3.5, 4, 5, 6, 7];
stimulus.gabor.contrasts = [1, .5, .25, .1, .075, .05];
stimulus.gabor.contrasts = 1;

stimulus.contrast_combinations = [1: ...
    length(stimulus.noise.contrasts) * length(stimulus.gabor.contrasts)];
stimulus.nPossibleContrasts = length(stimulus.contrast_combinations);
defineLocations;

%%%%% things to be randomized or to be saved
task{1}.randVars.calculated.noise_contrast = nan;  
task{1}.randVars.calculated.gabor_contrast = nan;
task{1}.randVars.calculated.gabor_location = [nan, nan]; 
task{1}.randVars.calculated.mousePos = [nan nan];
task{1}.randVars.calculated.detection_rt = nan;
task{1}.randVars.calculated.decision_rt = nan;
task{1}.randVars.calculated.response_offset = [nan nan];

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

while (task{1}.trialnum <= task{1}.numTrials) && ~myscreen.userHitEsc
    % update the task
    [task myscreen] = updateTask(task,myscreen,1);
    % flip the screen
    myscreen = tickScreen(myscreen, task);
end

% task ended
mglClearScreen(0.5);
mglTextSet([],32,1);
% get count
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

    % decide on contrasts
    stimulus.currentContrasts = randsample(stimulus.contrast_combinations,1);
    stimulus.contrast_combinations(stimulus.contrast_combinations == stimulus.currentContrasts) = [];
end

% decide on the location to present the gabor
stimulus.current_gabor_location = randsample(stimulus.gabor.nLoc,1);
task.thistrial.gabor_location = stimulus.gabor_locations(stimulus.current_gabor_location,:);

% decide on contrasts
index = reshape(1:stimulus.nPossibleContrasts, ...
    length(stimulus.noise.contrasts), length(stimulus.gabor.contrasts));
[noise_contrast, gabor_contrast] = find(index == stimulus.currentContrasts);
task.thistrial.noise_contrast = stimulus.noise.contrasts(noise_contrast);
task.thistrial.gabor_contrast = stimulus.gabor.contrasts(gabor_contrast);

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
    % present noise with gabor stimulus
    mglClearScreen(stimulus.bg_color);
    mglStencilSelect(1);
    mglBltTexture(stimulus.tex_target,[0 0])
    mglStencilSelect(0);    
    mglFlush
    
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
    mglTextDraw('Click on the location where the target appeared', [0,10])   
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
    task.thistrial.mousePos = mousePos;
    task.thistrial.response_offset = task.thistrial.gabor_location - [mInfo.x, mInfo.y];
    
    task = jumpSegment(task);    
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

% change contrast
contrast = task.thistrial.noise_contrast;
N = length(im(:));
m_im = mean(im(:));
coeff = sqrt((N*contrast^2) / sum((im(:)-m_im).^2));
stimulus.noise.im = coeff .* im;

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

% convert visual angle of the locations to pixels
displaySize = max(size(stimulus.pink_filter));
locations = visualAngleToPixels(locations_va, displaySize);

stimulus.gabor_locations = locations;

function combinedStimulus(task)
global stimulus
noise = stimulus.noise.im;
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
x_lims = [location(1)-floor(size(gabor_circle,1)/2), location(1)+floor(size(gabor_circle,1)/2)];
y_lims = [location(2)-floor(size(gabor_circle,2)/2), location(2)+floor(size(gabor_circle,2)/2)];

gabor_position = zeros(size(noise,1), size(noise,2));
gabor_position(x_lims(1):x_lims(2),y_lims(1):y_lims(2)) = gabor_circle;

% add gabor to the noise
final_im = noise + gabor_position;

% scale it to [0 255], for both stimulus images
final_im =  255 .* ((final_im + 1) ./ 2);
stimulus.final_im{1} = final_im';
stimulus.final_im{2} = (255 .* ((noise + 1)./2))';

% decide background color
bg_color = mean(stimulus.final_im{2}(:));
stimulus.bg_color = bg_color;


