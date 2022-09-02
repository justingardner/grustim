%  geislerInstructionFinal.m
% 
%      usage: geislerInstructionFinal.m
%         by: Lila Shroff, Jiwon Yeon
%       date: Aug 2022      
%    purpose: Provide instructions for participants before starting
%             geislerDetectionTask_staircase_jy.m

function geislerInstructionsFinal_jy

%% POTENTIAL BUGS
%
%  (HOLD) the white stencils both show up upon opening program (hold)

%% set up global variables 
% initalize 
mglClose        % close MGL if it's open
clear all, close all, clc

mglSetSID('test')
myscreen.displayName = 'dell-wuTsai';
% myscreen.displayName = 'home';
myscreen.saveData = 0;
myscreen.keyboard.left = 44;
myscreen.keyboard.right = 48;
myscreen = initScreen(myscreen);

%%%%% set stimulus parameter
global stimulus
stimulus.responsekeyleft = 44;  
stimulus.responsekeyright = 48;  
stimulus.responsekeys = [44,48];   % '<,' & '>.'
stimulus.noise.size = 15;   % visual angle
stimulus.noise.contrasts = .2;

stimulus.gabor.size = 1;    % visual angle
stimulus.gabor.tilt = 315;
stimulus.gabor.cycle = 6;
stimulus.gabor.contrasts = .5;
stimulus.gabor.nLoc = 25;

stimulus.noise.buffer = stimulus.noise.size + 3;     % visual angle
stimulus.noise.frame_pixel = visualAngleToPixels(stimulus.noise.buffer, ...
    [myscreen.screenWidth, myscreen.screenHeight]);

% define locations
define_locations(myscreen)

% create noise images
create_noise(myscreen)

% create gabor
create_gabor(myscreen)

% combine noise and gabor
stimulus.gabor.present_loc = 1;     % initial target position
combine_noise_gabor(myscreen)

% creating stencils
mglClearScreen(.5);
mglFlush();

% left stencil
mglStencilCreateBegin(1);
mglVisualAngleCoordinates(myscreen.displayDistance,myscreen.displaySize);
x = -stimulus.noise.size/2-1;
y = 0;
mglFillOval(x, y, [stimulus.noise.size, stimulus.noise.size]);
mglStencilCreateEnd;
keyboard;

% right stencil
mglStencilCreateBegin(2);
mglVisualAngleCoordinates(myscreen.displayDistance,myscreen.displaySize);
x = stimulus.noise.size/2+1;
y = 0;
mglFillOval(x, y, [stimulus.noise.size, stimulus.noise.size]);
mglStencilCreateEnd;

% center stencil
mglStencilCreateBegin(3);
mglVisualAngleCoordinates(myscreen.displayDistance,myscreen.displaySize);
mglFillOval(0, 0, [stimulus.noise.size, stimulus.noise.size]);
mglStencilCreateEnd;

%% change backtick to space bar 
backTickCharacter = ' ';
backTickCode = mglCharToKeycode({' '});
myscreen.keyboard.backtick = backTickCode;

%% present screen with text instructions
% page 1: welcome
mglClearScreen(.5);
mglTextSet('Helvetica', 32, [1 1 1]);
mglTextDraw('Thank you for participating in our research experiment.' ,[0,3]);
mglTextDraw('Before beginning the experiment, we will first give you brief instructions and a' ,[0,0]);
mglTextDraw('chance to test your understanding of the task.' ,[0,-1]);
mglTextDraw('Press the SPACE BAR when you are ready to continue to the next instruction screen.' ,[0,-10]);
mglFlush; 
mglWaitSecs(.5);
while 1
    k = mglGetKeys;
    if k(myscreen.keyboard.backtick)==1, break; end
end

% page 2: Written explanation about the task
mglTextDraw(['In this experiment, we will present two screens in quick succession.'] ,[0,2]);
mglTextDraw(['One screen will have a target while the other will not.'],[0,1]);
mglTextDraw(['Your task is to respond with which screen had the target.'] ,[0,0]);
mglTextDraw(['While you are performing the experiment, we will concurrently collect your eye movement data.'] ,[0,-2]);
mglTextDraw(['With the collected data, we aim to build a model that can predict eye movements. '] ,[0,-3]);
mglTextDraw(['Press the SPACE BAR when you are ready to continue to the next instruction screen.'] ,[0,-10]);
mglFlush;
mglWaitSecs(.5);
while 1
    k = mglGetKeys;
    if k(myscreen.keyboard.backtick)==1, break; end
end

% page 3: Black fixation dot
mglTextDraw(['When a trial begins, you will see a black dot like the one below at the center of the screen'] ,[0,10]);
mglTextDraw('The trial will be started when you hit a response key (which will be introduced shortly).' ,[0,9]);
mglTextDraw('Please focus your eyes at the center of the screen for the entire duration of the experiment.' ,[0,7]);
mglTextDraw(['Press the SPACE BAR when you are ready to continue to the next instruction screen.'] ,[0,-10]);
mglFillOval(0,0,[.2 .2],0);
mglFlush;
mglWaitSecs(.5);
while 1
    k = mglGetKeys;
    if k(myscreen.keyboard.backtick)==1, break; end
end

% page 4 PRESENT EXAMPLE STIMULI
mglTextDraw(['In each trial, two stimuli will flash on the screen in sequence.'], [0, 12]);
mglTextDraw(['A white dot will be briefly presented between screens.'], [0, 11]);
mglTextDraw(['Example stimuli are shown on the screen.'] ,[0,10]);
mglTextDraw(['As you can see, the LEFT stimulus has a target at the center while the RIGHT stimulus does not.'] ,[0,9]);
mglTextDraw('Press the SPACE BAR when you are ready to continue to the next instruction screen.' ,[0,-10]);

tex = mglCreateTexture(stimulus.final_im{1});   % image with target
mglStencilSelect(1);    
mglBltTexture(tex, [-stimulus.noise.size/2-1, 0]);
mglStencilSelect(0);    

tex = mglCreateTexture(stimulus.final_im{2});   % noise only image
mglStencilSelect(2);
mglBltTexture(tex, [stimulus.noise.size/2+1, 0]);
mglStencilSelect(0);
mglFlush;
mglWaitSecs(.5);
while 1
    k = mglGetKeys;
    if k(myscreen.keyboard.backtick)==1, break; end
end

% page 5: making responses
mglTextDraw('After presenting the two stimuli, you will be asked to identify which stimulus contained the target.' ,[0, 10]);
mglTextDraw('As the message indicates, press the COMMA button (<) to indicate that the FIRST image contained.' ,[0,9]);
mglTextDraw(['Otherwise, press the PERIOD button (>) to indicate that the SECOND image contained the target.'] ,[0,8]);
mglTextDraw('Which screen showed the target?',[0,1]);
mglTextDraw('1(<)  or  2(>)',[0, -2]);
mglTextDraw(['Press the SPACE BAR when you are ready to continue to the next instruction screen.'] ,[0,-10]);
mglFlush;
mglWaitSecs(.5);
while 1
    k = mglGetKeys;
    if k(myscreen.keyboard.backtick)==1, break; end
end

% page 6: feedback
mglTextDraw(['If your answer is correct, a green dot will appear.'] ,[0,10]);
mglTextDraw(['If your answer is incorrect, a red dot will appear. '] ,[0,9]);
mglTextDraw(['Then the dot will turn white to notify you to prepare for the next trial.'] ,[0,8]);
mglFillOval(0,0,[.2 .2],[0 1 0]);
mglTextDraw(['Press the SPACE BAR when you are ready to continue to the next instruction screen.'] ,[0,-10]);
mglFlush;
mglWaitSecs(.5);
while 1
    k = mglGetKeys;
    if k(myscreen.keyboard.backtick)==1, break; end
end

% page 7: Three example trials
mglTextDraw(['Now we will give you three example trials.'] ,[0, 1]);
mglTextDraw(['Remember that each trial will begin when you press one of the response buttons.'] ,[0, 0]);
mglTextDraw(['Press the SPACE BAR when you are ready to start the three practice trials.'] ,[0, -5]);
mglFlush;
mglWaitSecs(.5);
while 1    
    k = mglGetKeys;
    if k(myscreen.keyboard.backtick)==1, break; end
end

%% Three easy trials with high contrast value and long display time
mglTextDraw('THREE EXAMPLE TRIALS' ,[0,0]);
mglFlush;
mglWaitSecs(1);

% PRESENT EXAMPLE STIMULI
stimulus.gabor.contrasts = .4;  % adjust gabor a little bit
thistrial = 1;
while thistrial <= 3
    % decide which screen would show the target
    whichseg = randsample(2,1);

    % create new stimuli
    create_noise(myscreen)
    create_gabor(myscreen)
    combine_noise_gabor(myscreen)
    
    % present a black fixation dot
    mglClearScreen(.5);
    mglFillOval(0,0,[.2 .2],0);
    mglFlush();
    while 1
        keycode = mglGetKeys;
        if any(keycode==1), break; end
    end

    for thisseg = 1:2
        if thisseg == whichseg
            % target
            tex = mglCreateTexture(stimulus.final_im{1});
        else
            % noise
            tex = mglCreateTexture(stimulus.final_im{2});
        end
        
        % present the stimulus image
        mglClearScreen(.5);
        mglStencilSelect(3);
        mglBltTexture(tex,[0,0]);
        mglStencilSelect(0);
        mglFlush();
        mglWaitSecs(.7);
        
        if thisseg == 1
            % white dot
            mglClearScreen(.5);
            mglFillOval(0,0,[.2 .2],1);
            mglFlush();
            mglWaitSecs(.5);
        end
    end

    % get response
    mglClearScreen(.5);
    mglTextDraw('Which screen showed the target?',[0,1]);
    mglTextDraw('1(<)  or  2(>)',[0, -2]);
    mglFlush();
    while 1
        k = mglGetKeys;
        if any(k(stimulus.responsekeys)==1)     % get into the loop only when valid key is pressed 
            % correct 
            if (whichseg == 1 && k(stimulus.responsekeys(1))==1) || ...
                    (whichseg == 2 &&  k(stimulus.responsekeys(2))==1)
                feedback_color = [0 1 0];
                break;

            % incorrect
            else 
                feedback_color = [1 0 0];
                break;
            end
        end
    end
    
    % show feedback
    mglClearScreen(.5);
    mglFillOval(0,0,[.2 .2], feedback_color);
    mglFlush();
    mglWaitSecs(.7);

    % continue to the next trial
    thistrial = thistrial+1;    
end


%% Instruction for 10 additional trials
mglTextDraw(['Great work so far!'], [0, 2]);
mglTextDraw(['The real experiment will be more difficult. The target will be less clear and presented for shorter.'], [0, 1]);
mglTextDraw('Now you will do 10 practice trials that will increase in difficulty.', [0, -1]);
mglTextDraw('Remember to keep fixating your eyes at the center even though the black dot disappears.', [0, -2]);

mglTextDraw(['Press the SPACE BAR when you are ready to do 10 practice trials.'] ,[0,-8]);
mglFlush;
mglWaitSecs(.5);
while 1
    k = mglGetKeys;
    if k(myscreen.keyboard.backtick)==1, break; end
end

%% 10 practice trials while decreasing time length and contrast
mglTextDraw('TEN EXAMPLE TRIALS WITH INCREASING DIFFICULTY' ,[0,0]);
mglFlush;
mglWaitSecs(1);

stimulus.gabor.contrasts = .4;      % initial contrast level
t_stimulus = .7;    % initial stimulus presentation time
thistrial = 1;
while thistrial <= 10
    % decide which screen would show the target
    whichseg = randsample(2,1);

    % create new stimuli
    create_noise(myscreen)
    create_gabor(myscreen)
    combine_noise_gabor(myscreen)
    
    % present a black fixation dot
    mglClearScreen(.5);
    mglFillOval(0,0,[.2 .2],0);
    mglFlush();
    while 1
        keycode = mglGetKeys;
        if any(keycode==1), break; end
    end

    for thisseg = 1:2
        if thisseg == whichseg
            % target
            tex = mglCreateTexture(stimulus.final_im{1});
        else
            % noise
            tex = mglCreateTexture(stimulus.final_im{2});
        end
        
        % present the stimulus image
        mglClearScreen(.5);
        mglStencilSelect(3);
        mglBltTexture(tex,[0,0]);
        mglStencilSelect(0);
        mglFlush();
        mglWaitSecs(t_stimulus);

        if thisseg == 1
            % white dot
            mglClearScreen(.5);
            mglFillOval(0,0,[.2 .2],1);
            mglFlush();
            mglWaitSecs(.5);
        end
    end

    % get response
    mglClearScreen(.5);
    mglTextDraw('Which screen showed the target?',[0,1]);
    mglTextDraw('1(<)  or  2(>)',[0, -2]);
    mglFlush();
    while 1
        k = mglGetKeys;
        if any(k(stimulus.responsekeys)==1)     % get into the loop only when valid key is pressed 
            % correct 
            if (whichseg == 1 && k(stimulus.responsekeys(1))==1) || ...
                    (whichseg == 2 &&  k(stimulus.responsekeys(2))==1)
                feedback_color = [0 1 0];
                break;

            % incorrect
            else 
                feedback_color = [1 0 0];
                break;
            end
        end
    end

    % show feedback
    mglClearScreen(.5);
    mglFillOval(0,0,[.2 .2], feedback_color);
    mglFlush();
    mglWaitSecs(.7);
    
    % adjust contrast level and presentation time
    if stimulus.gabor.contrasts > .2
        stimulus.gabor.contrasts = stimulus.gabor.contrasts - .05;
    end
    if t_stimulus > .25
        t_stimulus = t_stimulus - .05;
    end

    % continue to the next trial
    thistrial = thistrial+1;    
end

%% Last five trials where target locations change
mglTextDraw(['Nice job!'] ,[0,4]);
mglTextDraw(['While previous practice trials showed the target only at the center,'], [0 3]);
mglTextDraw(['in the real experiment, the target will appear in other positions.'] ,[0,2]);
mglTextDraw(['We will let you know in advance where the target will appear in a new location.'] ,[0,0]);
mglTextDraw(['Please keep your eyes at the center even if the target appears somewhere else.'] ,[0,-1]);
mglTextDraw(['For your final practice, we will present 5 trials where the target appears in new positions.'] ,[0,-3]);
mglTextDraw(['Press the SPACE BAR when you are ready to do the last 5 practice trials.'] ,[0,-10]);
mglFlush;
mglWaitSecs(.5);
while 1
    k = mglGetKeys;
    if k(myscreen.keyboard.backtick)==1, break; end
end


thistrial = 1;
stimulus.gabor.contrasts = .45;     % adjust target to be more visible
while thistrial <= 5
    % decide which screen would show the target
    whichseg = randsample(2,1);

    % decide which location the target would appear
    stimulus.gabor.present_loc = randsample(stimulus.gabor.nLoc,1);

    % create new stimuli
    create_noise(myscreen)
    create_gabor(myscreen)
    combine_noise_gabor(myscreen)
    
    % present the target location with the black fixation dot
    mglClearScreen(.5);
    mglTextDraw(['Target location for this task'], [0, 10]);
    mglFillOval(0,0,[.2 .2],0);
    mglGluAnnulus(stimulus.gabor_locations_deg(stimulus.gabor.present_loc,1), ...
        stimulus.gabor_locations_deg(stimulus.gabor.present_loc,2), ...
        .8, 1, [1 1 1], 120, 2);
    mglFlush;
    mglWaitSecs(2);

    % present a black fixation dot
    mglClearScreen(.5);
    mglFillOval(0,0,[.2 .2],0);
    mglFlush();
    while 1
        keycode = mglGetKeys;
        if any(keycode==1), break; end
    end

    for thisseg = 1:2
        if thisseg == whichseg
            % target
            tex = mglCreateTexture(stimulus.final_im{1});
        else
            % noise
            tex = mglCreateTexture(stimulus.final_im{2});
        end
        
        % present the stimulus image
        mglClearScreen(.5);
        mglStencilSelect(3);
        mglBltTexture(tex,[0,0]);
        mglStencilSelect(0);
        mglFlush();
        mglWaitSecs(t_stimulus);        

        if thisseg == 1
            % white dot
            mglClearScreen(.5);
            mglFillOval(0,0,[.2 .2],1);
            mglFlush();
            mglWaitSecs(.5);
        end
    end

    % get response
    mglClearScreen(.5);
    mglTextDraw('Which screen showed the target?',[0,1]);
    mglTextDraw('1(<)  or  2(>)',[0, -2]);
    mglFlush();
    while 1
        k = mglGetKeys;
        if any(k(stimulus.responsekeys)==1)     % get into the loop only when valid key is pressed 
            % correct 
            if (whichseg == 1 && k(stimulus.responsekeys(1))==1) || ...
                    (whichseg == 2 &&  k(stimulus.responsekeys(2))==1)
                feedback_color = [0 1 0];
                break;

            % incorrect
            else 
                feedback_color = [1 0 0];
                break;
            end
        end
    end

    % show feedback
    mglClearScreen(.5);
    mglFillOval(0,0,[.2 .2], feedback_color);
    mglFlush();
    mglWaitSecs(.7);
  
    % continue to the next trial
    thistrial = thistrial+1;    
end


%% END INSTRUCTIONS.
mglTextDraw('Now we will move on to the real experiment. Please call the experimenter.' ,[0,0]);
mglTextDraw('If you have any questions, please ask the experimenter now.' ,[0,-2]);
mglFlush;

mglWaitSecs(5);
mglClose;


%%%%%% Helper functions %%%%%%
%% create_noise
function create_noise(myscreen)
global stimulus
% make the size of the image an odd number
frame_pixel = stimulus.noise.frame_pixel;
if mod(frame_pixel,2)==0, frame_pixel = frame_pixel+1; end

% pink filter is already saved as a file. just load the variable
if exist('geislerDetectionTask_pinkFilter.mat') == 2
    load('geislerDetectionTask_pinkFilter.mat')
end

% pink filter has not created before or needs to re-created
if exist('geislerDetectionTask_pinkFilter.mat') == 0 || size(pink_filter,1) < double(frame_pixel)
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
    floor(filter_sz(1)/2)+1-(frame_pixel-1)/2:floor(filter_sz(1)/2)+1+(frame_pixel-1)/2, ...
    floor(filter_sz(2)/2)+1-(frame_pixel-1)/2:floor(filter_sz(2)/2)+1+(frame_pixel-1)/2);

% create two noise images
for i_num = 1:2
    % fft on white noise
    white = randn(frame_pixel, frame_pixel);
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
    
    stimulus.noise.im{i_num} = noise;
end

%% create_gabor
function create_gabor(myscreen)
global stimulus
grating = mglMakeGrating(stimulus.gabor.size, stimulus.gabor.size, ...
    stimulus.gabor.cycle, stimulus.gabor.tilt, 0);
contrast = stimulus.gabor.contrasts;
grating = grating * contrast;
gaussian = mglMakeGaussian(stimulus.gabor.size, stimulus.gabor.size, 1, 1);
stimulus.gabor.im = (grating.*gaussian);

%% combine_noise_gabor
function combine_noise_gabor(myscreen)
global stimulus
% make circular gabor patch
radius_pixel = visualAngleToPixels(stimulus.gabor.size/2);
gabor_sz = size(stimulus.gabor.im);
[x y] = meshgrid(-(gabor_sz(1)-1)/2:(gabor_sz(1)-1)/2, ...
    -(gabor_sz(2)-1)/2:(gabor_sz(2)-1)/2);
stencil = sqrt(x.^2 + y.^2) <= radius_pixel;
gabor_circle = stencil' .* stimulus.gabor.im;

pos = stimulus.gabor_locations_onNoise(stimulus.gabor.present_loc,:);
gabor_position = zeros(size(stimulus.noise.im{1},1), size(stimulus.noise.im{1},2));
gabor_position([pos(2)-(gabor_sz(1)-1)/2:pos(2)+(gabor_sz(1)-1)/2], ...
    [pos(1)-(gabor_sz(2)-1)/2:pos(1)+(gabor_sz(2)-1)/2]) = gabor_circle;
gabor_position = flipud(gabor_position);

% add gabor to the noise
final_im = stimulus.noise.im{1} + gabor_position;

% scale it to [0 255], for both stimulus images
final_im =  255 .* ((final_im + 1) ./ 2);
final_im = final_im;
noise_only = 255 .* ((stimulus.noise.im{2} + 1) ./ 2);

% save it to stimulus structure
stimulus.final_im{1} = final_im;
stimulus.final_im{2} = noise_only;


%% define_locations
function define_locations(myscreen)
global stimulus
nLoc = stimulus.gabor.nLoc;
nLayer = floor(nLoc/8);
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
stimulus.gabor_locations_deg = locations;
gabor_locations_onNoise = visualAngleToPixels(locations, stimulus.noise.frame_pixel);
stimulus.gabor_locations_onNoise = gabor_locations_onNoise;
