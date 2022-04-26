% geislerDetectionTask_staircase_jy.m
% 
%      usage: geislerDetectionTask_staircase_jy.m
%         by: jiwon yeon
%       date: 
%  copyright: (c) 2022 Jiwon Yeon
%    purpose: Get contrast values for Najemnik&Geisler's 2005 experiment
%
%             The task presented two stimulus screens in a trial,
%             one with the task target (grating) and the other without.
%             Subjects have to indicate which screen contained the target
%             grating stimulus. 
%             The target locations are limited to left, right, up, and down
%             from the center. 
%             location throughout a block.
%

function geislerDetectionTask_staircase_jy
mglClose        % close MGL if it's open
clear all, close all, clc
global stimulus

testingLoc = input('Testing location?: ');
mglSetSID('test')
myscreen.eyetracker = 1;
myscreen.screenNumber = 2;
myscreen.saveData = 1;
myscreen.datadir = '~/proj/jiwon/data/geisler';
if ~exist(myscreen.datadir), mkdir(myscreen.datadir); end
myscreen.eyetracker = 1;
mglSetParam('abortedStimfilesDir', '~/proj/jiwon/data/geisler/aborted',1);

myscreen.keyboard.nums = [44,48]; % ',<' for 1, '.>' for 2
myscreen = initScreen(myscreen);  

% load pink_filter
if exist('~/proj/grustim/geislerDetectionTask_pinkFilter.mat') ~= 0
    load('~/proj/grustim/geislerDetectionTask_pinkFilter.mat');
    stimulus.pink_filter = pink_filter;
else
    createPinkFilter(myscreen);
end

%%%%% define task timings and responses
task{1}.waitForBacktick = 1;
task{1}.seglen = [inf, 1, .5, 1, inf, .5];  
%  fixation-stim1-int-stim2-response-feedback

task{1}.getResponse = [1 0 0 0 1 0];
stimulus.nBlocks = 1;
stimulus.cBlock = 1;    % current block
% nContrasts = 7;
stimulus.TrialsPerBlock = 40;
task{1}.numTrials = stimulus.nBlocks * stimulus.TrialsPerBlock;
stimulus.gabor.nLoc = 25;

%%%%% set stimulus parameter
stimulus.responsekeys = [44,48];   % '<,' & '>.'
stimulus.noise.size = 15;   % visual angle

stimulus.gabor.size = .5;    % visual angle
stimulus.gabor.tilt = 315;   % 315 degree
stimulus.gabor.cycle = 6;

%%%%% parameters used for the constant stimuli method %%%
% contrast_minmax = [.05 .2]; %[.2, .1, .075, .05];
% gabor_contrasts = logspace(contrast_minmax(1), contrast_minmax(2), nContrasts);
% gabor_contrasts = log10(gabor_contrasts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % define4Locations;   % four locations: right,left,up,& down
defineLocations;    % locations with predefined numbers
stimulus.gaborLoc_thisblock = testingLoc;    % one location per block
if stimulus.gaborLoc_thisblock < 10
    init_threshold = .13;
else
    init_threshold = .2;
end

% stimulus.locations_left = repmat(1:size(stimulus.gabor_locations,1), 1, ...
%     stimulus.nBlocks/size(stimulus.gabor_locations,1));

%%%%% initialize staircase
% stimulus.stair = doStaircase('init','fixed',['fixedVals=' num2str(gabor_contrasts)], ...
%     ['nTrials=' num2str(task{1}.numTrials)]);
stimulus.stair = doStaircase('init','upDown','nup=1','ndown=2',...
    ['initialThreshold=' num2str(init_threshold)], ...
    'initialStepsize=.05', ...
    'nTrials=40');


%%%%% things to be randomized or to be saved
task{1}.waitForBacktick = 0;
task{1}.parameter.noise_contrast = [.2];
task{1}.random = 1;

task{1}.randVars.uniform.whichseg = [2 4];    % at which segment to present the stimulus
task{1}.randVars.gabor_location = nan;
task{1}.randVars.calculated.gabor_contrast = nan;
task{1}.randVars.calculated.correct = nan;
task{1}.randVars.calculated.rt = nan;

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
% % Eye calibration (optional)
if isfield(myscreen,'eyetracker')
    disp(' Calibrating Eye ....')
    myscreen = eyeCalibDisp(myscreen);
end

mglClearScreen(.5);
mglTextSet([],32,1);
mglTextDraw('Starting experiment',[0,0]);
mglFlush;
mglWaitSecs(2)

while (task{1}.trialnum <= task{1}.numTrials) && ~myscreen.userHitEsc    
    % update the task
    [task myscreen] = updateTask(task,myscreen,1);
    % flip the screen
    myscreen = tickScreen(myscreen, task);
end

% task ended
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

% for every new block, give a short break
if (mod(task.trialnum, stimulus.TrialsPerBlock) == 1) && (task.trialnum ~= 1)
    mglClearScreen(.5)
    mglTextSet([],32,1);
    mglTextDraw(['Take a short break'],[0,.7])
    mglTextDraw(['Press any keys when you are ready'], [0,-.7])
    mglFlush    
    
    % update current block
    stimulus.cBlock = stimulus.cBlock + 1;
    
    % Listen keys
    while 1
        k= mglGetKeys;
        if (any(k)), break; end
    end
end

% % decide on the gabor location
% if mod(task.trialnum, stimulus.TrialsPerBlock)==1
%     index = randsample(1:length(stimulus.locations_left),1);
%     stimulus.gaborLoc_thisblock = stimulus.locations_left(index);
%     stimulus.locations_left(index) = [];
% end
task.thistrial.gabor_location = stimulus.gaborLoc_thisblock;

% generate noise images
createPinkNoise(myscreen, task);

% generate gabor
[gabor_contrast, stimulus.stair] = doStaircase('testValue', stimulus.stair);
task.thistrial.gabor_contrast = gabor_contrast;
createGabor(task);

% convert gabor locations to pixels
if task.trialnum == 1
    % convert visual angle of the locations to pixels
    displaySize = size(stimulus.noise.im{1});
    stimulus.gabor_locations_pix = visualAngleToPixels(stimulus.gabor_locations_deg, displaySize);
end

% combine noise and gabor
combinedStimulus(task);

% create texture in advance
stimulus.tex_target = mglCreateTexture(stimulus.final_im{1});
stimulus.tex_nontarget = mglCreateTexture(stimulus.final_im{2});

disp(sprintf('Trial: %i Contrast: %0.2f',task.trialnum,task.thistrial.gabor_contrast));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)
global stimulus

if task.thistrial.thisseg == 1
    if mod(task.trialnum, stimulus.TrialsPerBlock) == 1
        % show how many blocks are left 
        mglClearScreen(.5)
        mglTextSet([],32,1);
        mglTextDraw(sprintf('Starting block %d out of %d blocks', ...
            stimulus.cBlock, stimulus.nBlocks),[0,0])
        mglFlush
        mglWaitSecs(2)
        
        % show where the target will appear
        mglClearScreen(.5)
        mglTextSet([],32,1);
        mglTextDraw(['Target location for this block'], [0, 10])
        sz = size(stimulus.pink_filter,1);
        target_location = stimulus.gabor_locations_deg(task.thistrial.gabor_location,:);
        mglGluAnnulus(target_location(1), target_location(2), .35, .4, ...
            [1 1 1], 120, 2)
        mglFillOval(0,0,[.2 .2],0)
        mglFlush;
        mglWaitSecs(3);
    end
        
    % show a fixation cross and wait for the button press
    mglClearScreen(.5)
    mglFillOval(0,0,[.2 .2],0)
    
elseif task.thistrial.thisseg == 3
    % present a screen with a black dot
    mglClearScreen(.5)
    mglFillOval(0,0,[.2 .2],1)
    
elseif task.thistrial.thisseg == 5
    % present a screen with a white dot and wait for the response
    mglClearScreen(.5)
    mglTextSet([],32,1);
    mglTextDraw(['Which screen showed the target?'],[0,1])
    mglTextDraw(['1(<)  or  2(>)'],[0 -1]);
    
elseif task.thistrial.thisseg == 6
    % feedback
    mglClearScreen(.5)
    mglFillOval(0,0,[.2 .2],stimulus.feedback_color);
    
    % show the target location
    sz = size(stimulus.pink_filter,1);
    target_location = stimulus.gabor_locations_deg(task.thistrial.gabor_location,:);    
    mglGluAnnulus(target_location(1), target_location(2), .35, .4, ...
        [1 1 1], 120, 2)
    
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
    % waiting for the subject to start the trial
    while 1
        keycode = mglGetKeys;
        if any(keycode(stimulus.responsekeys)==1)
            break
        end
    end
    % move to the next segment
    task = jumpSegment(task);
    
elseif task.thistrial.thisseg == 5
    % get the response 
    while 1
        keycode = mglGetKeys;
        if any(keycode(stimulus.responsekeys)==1)
            % save whether the response was correct
            if (keycode(stimulus.responsekeys(1)) == 1 && task.thistrial.whichseg == 2) ...
                    || (keycode(stimulus.responsekeys(2)) == 1 && task.thistrial.whichseg == 4)
                % correct trial
                task.thistrial.correct = 1;
                task.thistrial.rt = task.thistrial.reactionTime;
                stimulus.feedback_color = [0 1 0];  % green                
            else
                % incorrect trial
                task.thistrial.correct = 0;                
                task.thistrial.rt = task.thistrial.reactionTime;
                stimulus.feedback_color = [1 0 0];  % red
            end
            
            % update staircase
            stimulus.stair = doStaircase('update', stimulus.stair, ...
                task.thistrial.correct, task.thistrial.gabor_contrast);            
            break
        end
    end
    
    % start a new trial
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
stimulus.gabor.im = grating.*gaussian;

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
stimulus.gabor_locations_deg = locations;

function define4Locations
global stimulus
% for the staircase, use only 4 locations - up, down, left, and right
nLoc = 4;
nLayer = 1;
radius_va = linspace(0, stimulus.noise.size/2, nLayer+2);     % radius in visual angle
radius_va = radius_va(2:end-1);

theta = linspace(0, 2*pi, nLoc+1);
theta(end) = [];

% determine locations - in visual angle
locations = [];
cTheta = 0;     % current theta
cLayer = 1;     % current layer
for cLoc = 1:nLoc
    cTheta = cTheta + 1;
    x_pos = radius_va(cLayer) * cos(theta(cTheta));
    y_pos = radius_va(cLayer) * sin(theta(cTheta));    
    locations = [locations; [x_pos, y_pos]];
    
    if cTheta == 8, cTheta = 0; end
    if mod(cLoc,8) == 0, cLayer = cLayer+1; end
end
stimulus.gabor_locations_deg = locations;

function combinedStimulus(task)
global stimulus
noise = stimulus.noise.im{1};
gabor = stimulus.gabor.im;
location_num = task.thistrial.gabor_location;   
location = stimulus.gabor_locations_pix(location_num,:);

% make circular gabor patch
radius_pixel = visualAngleToPixels(stimulus.gabor.size/2);
gabor_sz = size(gabor);
[stencil_x, stencil_y] = meshgrid(-(gabor_sz(1)-1)/2:(gabor_sz(1)-1)/2, ...
    -(gabor_sz(2)-1)/2:(gabor_sz(2)-1)/2);
stencil = (sqrt(stencil_x.^2 + stencil_y.^2) <= radius_pixel);
gabor_circle = stencil' .* gabor;

% determine the location to display
gabor_position = zeros(size(noise,1), size(noise,2));
gabor_position([location(1)-(gabor_sz(1)-1)/2:location(1)+(gabor_sz(1)-1)/2], ...
    [location(2)-(gabor_sz(2)-1)/2:location(2)+(gabor_sz(2)-1)/2]) = gabor_circle;

% add gabor to the noise
final_im = noise + gabor_position;

% clip to max and min
final_im(final_im > 1) = 1;
final_im(final_im < -1) = -1;

% scale it to [0 255], for both stimulus images
final_im =  255 .* ((final_im + 1) ./ 2);
stimulus.final_im{1} = final_im';
stimulus.final_im{2} = 255 .* ((stimulus.noise.im{2} + 1) ./ 2);

% decide background color
for image = 1:2
    if task.thistrial.noise_contrast == 0
        bg_color = stimulus.final_im{image}(1,1);
    else
        bg_color = mean(stimulus.final_im{image}(:));
    end
    stimulus.bg_color{image} = bg_color;
end


