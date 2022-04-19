% geislerDetectionTask_jy.m
% 
%      usage: geislerDetectionTask_jy.m
%         by: jiwon yeon
%       date: 
%  copyright: (c) 2022 Jiwon Yeon
%    purpose: Replicating Najemnik&Geisler's 2005 study
%
%             The detection task presented two stimulus screens in a trial,
%             one with the task target (grating) and the other without.
%             Subjects have to indicate which screen contained the target
%             grating stimulus. The grating target presented at the same
%             location throughout a block.
%

function geislerDetectionTask_jy
mglClose        % close MGL if it's open
clear all, close all, clc
global stimulus

myscreen.screenNumber = 2;
myscreen.saveData = 1;
myscreen.datadir = '~/proj/data/geislerDetectionTask';
myscreen.eyetracker = 1;
mglSetParam('abortedStimfilesDir', '~/proj/data/geislerDetectionTask/aborted',1);

myscreen.keyboard.nums = [44,48]; % ',<' for 1, '.>' for 2
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
task{1}.seglen = [inf, .25, .5, .25, inf, .6];  
task{1}.getResponse = [1 0 0 0 1 0];
%  fixation-stim1-int-stim2-response-feedback

%%%%% set stimulus parameter
stimulus.responsekeys = [44,48];   % space bar
stimulus.noise.size = 15;   % visual angle
stimulus.noise.contrasts = [.2];

stimulus.gabor.size = .5;    % visual angle
stimulus.gabor.tilt = 315;
stimulus.gabor.cycle = 6;
stimulus.gabor.nLoc = 25;   % 25 for the real experiment
contrasts = logspace(.1, .2, 4);
stimulus.gabor.contrasts = [.05, log10(contrasts)];%[.2, .1, .075, .05];

defineLocations;
stimulus.locations_left = 1:stimulus.gabor.nLoc;

%%%%% parameters regarding to the experiment duration
stimulus.nBlocks = stimulus.gabor.nLoc;
stimulus.cBlock = 1;    % current block
stimulus.TrialsPerBlock = length(stimulus.gabor.contrasts) * 5;
task{1}.numTrials = stimulus.nBlocks * stimulus.TrialsPerBlock;

%%%%% things to be randomized or to be saved
task{1}.random = 1;
task{1}.parameter.gabor_contrast = repmat(stimulus.gabor.contrasts, 1, ...
    stimulus.TrialsPerBlock/length(stimulus.gabor.contrasts));
task{1}.parameter.noise_contrast = stimulus.noise.contrasts;

task{1}.randVars.uniform.whichseg = [2 4];    % at which segment to present the stimulus
task{1}.randVars.calculated.gabor_location = nan;    % start with a random position
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
% if ~myscreen.eyetracker
%     disp(' Calibrating Eye ....')
%     myscreen = eyeCalibDisp(myscreen);
% end

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

% for every new block,
% give a short break, if it is not the first block
if (task.trialnum ~= 1) && (mod(task.trialnum, stimulus.TrialsPerBlock) == 1)
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

if mod(task.trialnum, stimulus.TrialsPerBlock) == 1
    % decide on the target location    
    index = randsample(1:length(stimulus.locations_left),1);
    stimulus.gaborLoc_thisblock = stimulus.locations_left(index);
    stimulus.locations_left(index) = [];
end
task.thistrial.gabor_location = stimulus.gaborLoc_thisblock;

% generate noise images
createPinkNoise(myscreen, task);

% generate gabor
createGabor(task);

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
        % present how many blocks are left
        mglClearScreen(.5)
        mglTextSet([],32,1);
        mglTextDraw(sprintf('Starting block %d out of %d blocks', ...
            stimulus.cBlock, stimulus.nBlocks),[0,0])
        mglFlush
        mglWaitSecs(2)
        
        % show where the target will appear
        mglClearScreen(.5)
        mglTextSet([],32,1);
        mglTextDraw(['Target location for this block'], [0, 12])
        sz = size(stimulus.pink_filter);
        target_location = stimulus.gabor_locations(task.thistrial.gabor_location,:);
        target_location = pixelsToVisualAngle(target_location,sz);
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

if task.thistrial.thisseg == 1 || task.thistrial.thisseg==2
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
                task.thistrial.correct = 1;
                task.thistrial.rt = task.thistrial.reactionTime;
                stimulus.feedback_color = [0 1 0];  % green
            else
                task.thistrial.correct = 0;                
                task.thistrial.rt = task.thistrial.reactionTime;
                stimulus.feedback_color = [1 0 0];  % red
            end            
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
noise = stimulus.noise.im{1};
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
x_lims = [location(1)-ceil(size(gabor_circle,1)/2)+1, location(1)+ceil(size(gabor_circle,1)/2)-1];
y_lims = [location(2)-ceil(size(gabor_circle,2)/2)+1, location(2)+ceil(size(gabor_circle,2)/2)-1];

gabor_position = zeros(size(noise,1), size(noise,2));
gabor_position(x_lims(1):x_lims(2),y_lims(1):y_lims(2)) = gabor_circle;

% add gabor to the noise
final_im = noise + gabor_position;

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


