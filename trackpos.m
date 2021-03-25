%        $Id: $
%      usage: trackpos
%         by: Josh Ryu
%       date: 05/15/2019
%    purpose: 

function myscreen = trackpos(varargin)
 

%getArgs(varargin,{'subjectID=s999','centerX=10','centerY=0','diameter=16'}); getArgs(varargin,{'subjectID=-1'});
% set up screen
if isempty(mglGetSID)
    myscreen.subjectID  = -1;
else
    myscreen.subjectID  = mglGetSID;
end

%myscreen.displayName = 'debug'; myscreen.screenNumber = 1; 
%myscreen.screenWidth = 860; myscreen.screenHeight = 600; 
%myscreen.hideCursor = 1;
myscreen                = initScreen(myscreen);

%% parameters
global stimulus; stimulus = struct;

% Experimenter parameters
noeye           = 0; % 1 if no eyetracking (mouse for eye); 0 if there is eye tracking `
showmouse       = 0; 
grabframe       = 0; 
whitenoiseOn    = 0; % 1: white noise; 2: 
fixateCenter    = 1;
phasescrambleOn = 1;
backprecompute  = 1;

% Task design (might be changed later, so check this)
% S1: Stimulus (30s) 
% S2: Fixation (3s)
task{1}{1}.segmin = [30 3]; % fixation for shorter bc of the segment start takes time.
task{1}{1}.segmax = [30 3]; 
task{1}{1}.numTrials = 5; % changed later depending on the condition
task{1}{1}.getResponse = [1 0]; %segment to get response.
task{1}{1}.waitForBacktick = 0; %wait for backtick before starting each trial 

% task parameters for adaptation conditions
if whitenoiseOn == 1 || phasescrambleOn == 1
    task{1}{1}.parameter.phasescrambleOn    = 1;
    task{1}{1}.parameter.backLum            = 96; %160;%90;  % background luminance; units: luminance 
    task{1}{1}.parameter.noiseLum           = 32;
else 
    task{1}{1}.parameter.backLum = 32;  % background luminance; units: fraction of full luminance 
end
task{1}{1}.parameter.stimLum = 255 - task{1}{1}.parameter.backLum;  % stimulus luminance (out of 255)
% task{1}{1}.parameter.stimStep = 6; % stimulus velocity (standard deviation) in deg/sec

% The pilot test has three main parts:
% 1. motor gain estimation through testing different stimulus speed
% 2. adaptation trials so that the subjects learn priors
% 3. test different uncertainty values.

% calculated parameters
task{1}{1}.randVars.calculated.randomSeed = nan;
task{1}{1}.randVars.calculated.initStim = [nan nan];
task{1}{1}.randVars.calculated.trackStim = nan(ceil(task{1}{1}.segmax(1)*myscreen.framesPerSecond),2);
task{1}{1}.randVars.calculated.trackResp = nan(ceil(task{1}{1}.segmax(1)*myscreen.framesPerSecond),2);
task{1}{1}.randVars.calculated.trackEye  = nan(ceil(task{1}{1}.segmax(1)*myscreen.framesPerSecond),2);
task{1}{1}.randVars.calculated.trackTime = nan(1,ceil(task{1}{1}.segmax(1)*myscreen.framesPerSecond));
task{1}{1}.randVars.calculated.trackEyeTime = nan(1,ceil(task{1}{1}.segmax(1)*myscreen.framesPerSecond)); % for referencing edf file

%% Set up tasks 

%motor gain estimation
%{
% The pilot test has three main parts:
% 1. motor gain estimation through testing different stimulus speed at full
luminance
% 2. adaptation trials so that the subjects learn priors
% 3. test different uncertainty values.

phaseN = 1; 
task{1}{phaseN} = task{1}{1};
task{1}{phaseN}.parameter.backLum = 255*0.25;  % background luminance; units: luminance 
task{1}{phaseN}.parameter.stimLum = 255 - task{1}{phaseN}.parameter.backLum;  % stimulus luminance (out of 255)
task{1}{phaseN}.numTrials = 20;
task{1}{phaseN}.parameter = rmfield(task{1}{phaseN}.parameter,'stimStep');
task{1}{phaseN}.randVars.uniform.stimStep = [2,6,8,12,16];
% [task{1}{phaseN} myscreen] = addTraces(task{1}{phaseN},myscreen,'trackStimX','trackStimY','trackRespX','trackRespY');

% changing stimulus luminance
phaseN = 1; %code phase 2 first cuz its vanilla.
task{1}{phaseN} = task{1}{1};
% [task{1}{phaseN} myscreen] = addTraces(task{1}{phaseN},myscreen,'trackStimX','trackStimY','trackRespX','trackRespY');

teststimLum = [40,30,25,20,15,10]; % main task 
mainphaseN = 2; %starting phase number for the main task 
for phaseN = mainphaseN:(mainphaseN-1+length(teststimLum)) %adaptation trials
    task{1}{phaseN} = task{1}{mainphaseN-1};
    task{1}{phaseN}.parameter.stimLum = teststimLum(phaseN-mainphaseN+1);
    % [task{1}{phaseN} myscreen] = addTraces(task{1}{phaseN},myscreen,'trackStimX','trackStimY','trackRespX','trackRespY');
end

% changing stimulus speed
teststimSteps = [5,10,15,20,25];
teststimLum   = 20;

for idx = 1:length(teststimSteps)
    % adaptation task
    task{1}{2*idx-1} = task{1}{1};  
    task{1}{2*idx-1}.parameter.stimStep = teststimSteps(idx);

    task{1}{2*idx}   = task{1}{2*idx-1};
    task{1}{2*idx}.parameter.stimLum    = teststimLum;
end
%}

% change stimulus speed and luminance; cross conditions.
teststimSteps = [1.5]; %[0.75,1.5,2.25]; %[5,15,25];
teststimLum   = linspace(task{1}{1}.parameter.stimLum, task{1}{1}.parameter.noiseLum,3);


for stepIdx = 1:length(teststimSteps)
    idx = (stepIdx-1)*(length(teststimLum)+1) +1;
    % adaptation condition, at full luminance
    task{1}{idx}                                = task{1}{1};  
    task{1}{idx}.parameter.stimStep             = teststimSteps(stepIdx);
    task{1}{idx}.parameter.phasescrambleOn      = 0;
    task{1}{idx}.parameter.noiseLum             = 0;    
    task{1}{idx}.segmin         = [5 3]; %fixation time constrained by the texture loading
    task{1}{idx}.segmax         = [5 3]; 
    task{1}{idx}.numTrials      = 1;
    
    for lumIdx = 1:length(teststimLum)
        % adaptation task
        idx = (stepIdx-1)*(length(teststimLum)+1) + lumIdx+1; %cycle through the luminances first.
        
        task{1}{idx}            = task{1}{1};  
        
        task{1}{idx}.parameter.stimStep         = teststimSteps(stepIdx);
        task{1}{idx}.parameter.stimLum          = teststimLum(lumIdx);
        task{1}{idx}.parameter.phasescrambleOn  = 1;
        task{1}{idx}.segmin     = [5 3]; %fixation time constrained a bit by the texture loading
        task{1}{idx}.segmax     = [5 3]; 
        task{1}{idx}.numTrials  = 1;
        task{1}{idx}.parameter.noiseLum           = 32;
    end
end

%% initialize
% intiailize task
disp(' Initializing Task....')

for phaseN = 1:length(task{1})
    [task{1}{phaseN} myscreen] = initTask(task{1}{phaseN},myscreen,...
        @startSegmentCallback,@screenUpdateCallback,@responseCallback,@initTrialCallback);
end

% initialize stimulus
disp(' Initializing Stimulus....') 

myscreen = initStimulus('stimulus',myscreen); % what does this do???
stimulus.whitenoiseOn       = whitenoiseOn;
stimulus.noeye              = noeye;
stimulus.showmouse          = showmouse;
stimulus.fixateCenter       = fixateCenter;
stimulus.backprecompute     = backprecompute;

stimulus.grabframe = grabframe; %save frames into matrices
if stimulus.grabframe
    global frame
    frame = {};
end

%% Eye calibration
if ~stimulus.noeye
    disp(' Calibrating Eye ....')
    myscreen = eyeCalibDisp(myscreen);
    
    % let the user know
    disp(sprintf('(trackpos) Starting Run...'));
end

%% run the task
disp(' Running Task....'); stimulus.t0 = mglGetSecs; % 

% let the experimentee know too...
mglClearScreen(task{1}{1}.parameter.backLum/255);
mglTextDraw('task (trackpos) starting... ', [0 0.5])
mglTextDraw('Track the brightest point of the screen with the red mouse cursor',[0 -0.5]);
mglFlush

if ~stimulus.showmouse, mglDisplayCursor(0);, end %hide cursor

phaseNum = 1;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);     % update the task
    myscreen = tickScreen(myscreen,task);     % flip screen
end

%% End task
disp(' End Task....')

if isfield(stimulus,'backnoise')
    for idx = 1:length(stimulus.backnoise)
        mglDeleteTexture(stimulus.backnoise{idx});
    end
end

myscreen = endTask(myscreen,task);
mglClose 
endScreen(myscreen); mglDisplayCursor(1) %show cursor

if stimulus.grabframe
    save('/Users/joshryu/Dropbox/GardnerLab/Stanford/Current/FYP/FYP talk/trackposTask.mat', 'frame')
end

end

%% Initialize trials 
function [task myscreen] = initTrialCallback(task, myscreen)
    global stimulus    
       
    stimulus.stepStd    = task.thistrial.stimStep/sqrt(myscreen.framesPerSecond); % convert to std per frame; i.e. variance of shz gaussian RVs.
    stimulus.stimLum    = task.thistrial.stimLum;
    
    stimulus.phasescrambleOn    = task.thistrial.phasescrambleOn;
    stimulus.backLum            = task.thistrial.backLum;
    stimulus.noiseLum           = task.thistrial.noiseLum;
    
    stimulus = myInitStimulus(stimulus,myscreen,task); %centerX,Y, diameter called by getArgs.
    
    % save seed for generating random textures
    rng('shuffle','twister'); 
    s=rng; 
    task.thistrial.randomSeed = s.Seed;
    rng(task.thistrial.randomSeed,'twister');
    
    %% noise
    % background noise
    downsample_timeRes  = 1; stimulus.downsample_timeRes = downsample_timeRes;
    downsample_spatRes  = 10;
        
    nframes             = myscreen.framesPerSecond*task.segmax(1);%/downsample_timeRes; 

    if stimulus.phasescrambleOn == 1 && stimulus.noiseLum > 1
        disp('Loading phase scrambled background noise')
        tic
        if stimulus.backprecompute == 1;
            savefile            = '/Users/joshua/data/trackpos/trackpos.mat';
            if ~exist(savefile,'file')
                savefile            = '/Users/jryu/proj/grustim/trackpos.mat';
            end
            backsetidx          = randi(50); %1; %randi(63); % choose a random background set
            load(savefile,['backgroundnoise_rgb' num2str(backsetidx)]);
          
            if isfield(stimulus,'backnoise')
                for idx = 1:nframes/downsample_timeRes +1 
                    mglDeleteTexture(stimulus.backnoise{idx});
                end
            end

            % permute time to add more randomness
            perm = randperm(size(eval(['backgroundnoise_rgb' num2str(backsetidx)]),4)); 
            
            for idx = 1:nframes/downsample_timeRes + 1
                stimulus.backnoise{idx} = mglCreateTexture(eval(['backgroundnoise_rgb' num2str(backsetidx) '(:,:,:,perm(idx))'])); %0.02s
            end
            
            clearvars(['backgroundnoise_rgb' num2str(backsetidx)])
            toc
        else
            if ~isfield(stimulus,'gaussianFFT')
                xsize_deg          = round(myscreen.imageWidth/downsample_spatRes);
                ysize_deg          = round(myscreen.imageHeight/downsample_spatRes);

                backgaussian = mglMakeGaussian(xsize_deg,ysize_deg,...
                    stimulus.stimStd/downsample_spatRes,stimulus.stimStd/downsample_spatRes)*255;
                stimulus.gaussianFFT = getHalfFourier(backgaussian);
            end 

            if isfield(stimulus,'backnoise')
                for idx = 1:nframes/downsample_timeRes
                    mglDeleteTexture(stimulus.backnoise{idx});
                end
            end

            for idx = 1:nframes/downsample_timeRes + 1
                back                        = stimulus.gaussianFFT; %0.02s
                back.phase                  = rand(size(back.mag))*2*pi; % scramble phase % 0.02s
                backgroundnoise             = round(reconstructFromHalfFourier(back));   %0.04s
                if idx == 1 % to save time; only allocate memory once.
                    backgroundnoise_rgb         = 255*ones(4,size(backgroundnoise,2),size(backgroundnoise,1),'uint8'); %0.1165 s
                end % otherwise just overwrite.
                backgroundnoise_rgb(4,:,:)  = backgroundnoise'/max(max(backgroundnoise))*stimulus.noiseLum;  % normalize contrast %0.025s
                backgroundnoise_rgb         = uint8(backgroundnoise_rgb); %0.02s
                stimulus.backnoise{idx}     = mglCreateTexture(backgroundnoise_rgb); %0.02s
            end
            toc
            %mglClearScreen(0);mglBltTexture(stimulus.backnoise{1},[0 0 myscreen.imageWidth myscreen.imageHeight]);mglFlush
        end
    elseif stimulus.whitenoiseOn == 1
        tic
        text_xsize          = round(myscreen.screenWidth/downsample_spatRes);
        text_ysize          = round(myscreen.screenHeight/downsample_spatRes);

        backgroundnoise         = round(rand(text_xsize,text_ysize,nframes)*stimulus.noiseLum);  
        backgroundnoise_rgb     = 255*ones(4,size(backgroundnoise,1),size(backgroundnoise,2),size(backgroundnoise,3),'uint8');
        backgroundnoise_rgb(4,:,:,:)    = backgroundnoise; % change alpha. full rgb.
        backgroundnoise_rgb             = uint8(backgroundnoise_rgb); 
        
        if isfield(stimulus,'backnoise')
            for idx = 1:nframes/downsample_timeRes + 1
                mglDeleteTexture(stimulus.backnoise{idx});
                stimulus.backnoise{idx} = mglCreateTexture(backgroundnoise_rgb(:,:,:,idx));    
            end
        else
            for idx = 1:nframes/downsample_timeRes + 1
                stimulus.backnoise{idx} = mglCreateTexture(backgroundnoise_rgb(:,:,:,idx));    
            end
        end
        toc
    end
end

%% Start segment
function [task myscreen] = startSegmentCallback(task, myscreen)
% S1: Stimulus (30s)
% S2: Fixation (3s)
global stimulus 

if task.thistrial.thisseg == 1       
    %% set stimulus
    % after the noise bc we want to recreate the texture with the seed.
    stimulus = myInitStimulus(stimulus,myscreen,task); %centerX,Y, diameter called by getArgs.
    
    % set mouse position to the stimulus direction. 
    x_img = stimulus.position(1);  y_img = stimulus.position(2);
    x_screen = x_img*myscreen.screenWidth/myscreen.imageWidth + myscreen.screenWidth/2;
    y_screen = y_img*myscreen.screenHeight/myscreen.imageHeight + myscreen.screenHeight/2;
    mglSetMousePosition(ceil(x_screen),floor(y_screen), myscreen.screenNumber); % correct for screen resolution???
    if ~stimulus.showmouse, mglDisplayCursor(0);, end %hide cursor
    
    task.thistrial.initStim = stimulus.position;% [stimx, stimy];
    
    %{ 
    see the stimulus
    for idx = 1:length(stimulus.backnoise)
        mglClearScreen(0);
        mglBltTexture(stimulus.backnoise{idx},[0 0 myscreen.imageWidth myscreen.imageHeight])
        mglBltTexture(stimulus.gaussian,stimulus.position);
        mglFlush
        stimulus        = updateTarget(stimulus,myscreen,task); % update position.
        pause(1/myscreen.framesPerSecond)
    end
    %}   
    %% frame counter.
    task.thistrial.framecount = 0;
    
    if stimulus.grabframe
        global frame
        %save('/Users/jryu/Dropbox/Stanford/Current/FYP/FYP talk/figures/trackposTask_nonoise_160back.mat', 'frame','-v7.3')
        %save('/Users/jryu/Dropbox/Stanford/Current/FYP/FYP talk/figures/trackposTask_noise_160back.mat', 'frame','-v7.3')
        frame = {};
        frame{task.segmax(1)*myscreen.framesPerSecond} = [];
    end
    
    % time debugging
    % stimulus.timedebug = nan(15,ceil(task.segmax(1)*myscreen.framesPerSecond)); %nanmean(diff(stimulus.timedebug),2)
else %intertrial interval
    % *** fixation cross.
end    

end

%% screen update
function [task myscreen] = screenUpdateCallback(task, myscreen)
% S1: Stimulus (30s)
% S2: Fixation (3s)

%% Update Screen

global stimulus % call stimulus. Takes ~0.000013 s.
% stimulus.timedebug(9,task.thistrial.framecount+1) = mglGetSecs(stimulus.t0); % takes ~0.0087425s
mglClearScreen(stimulus.backLum/255);

if (task.thistrial.thisseg== 1)
    task.thistrial.framecount = task.thistrial.framecount + 1;
    
    % stimulus.timedebug(10,task.thistrial.framecount) = mglGetSecs(stimulus.t0); % takes ~ 6.264918 e-5 s
    % stimulus.timedebug(1,task.thistrial.framecount+1) = mglGetSecs(stimulus.t0);

    stimulus        = updateTarget(stimulus,myscreen,task); % update position.
    
    % stimulus.timedebug(2,task.thistrial.framecount+1) = mglGetSecs(stimulus.t0); % takes ~0.00012s 
    
    if stimulus.phasescrambleOn == 1 && mod(task.thistrial.framecount, stimulus.downsample_timeRes) == 1
        mglBltTexture(stimulus.backnoise{ceil(task.thistrial.framecount/stimulus.downsample_timeRes)},...
            [0 0 myscreen.imageWidth myscreen.imageHeight])
    end
    mglBltTexture(stimulus.gaussian,stimulus.position);
    
    % stimulus.timedebug(3,task.thistrial.framecount+1) = mglGetSecs(stimulus.t0); % takes ~0.000495s

    % **&display mouse position
    mInfo = mglGetMouse(myscreen.screenNumber);
    mimg_x = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
    mimg_y = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;
    mglGluDisk(mimg_x, mimg_y, 0.1, [1 0 0])
    % stimulus.timedebug(4,task.thistrial.framecount+1) = mglGetSecs(stimulus.t0); %takes ~9.2678e-5 s

    % ***record stimulus position and mouse position  
    task.thistrial.trackStim(task.thistrial.framecount,:) = stimulus.position;
    task.thistrial.trackResp(task.thistrial.framecount,:) = [mimg_x, mimg_y];
    task.thistrial.trackTime(task.thistrial.framecount)   = mglGetSecs(stimulus.t0);
    
    % stimulus.timedebug(5,task.thistrial.framecount+1) = mglGetSecs(stimulus.t0); %takes ~7.67331e-5 s

    if stimulus.fixateCenter == 1
        mglGluAnnulus(0,0,0.2,0.3,stimulus.fixColor,60,1);
        mglGluDisk(0,0,0.1,rand(1,3),60,1);
    end
    
elseif (task.thistrial.thisseg == 2)
    if stimulus.fixateCenter == 1 % stop the flashing
        rng(task.thistrial.randomSeed,'twister');
        mglGluDisk(0,0,0.1,rand(1,3),60,1);
    end
    
    mglGluAnnulus(0,0,0.2,0.3,stimulus.fixColor,60,1);
end

% fixation cross for all tasks. 


% stimulus.timedebug(6,task.thistrial.framecount+1) = mglGetSecs(stimulus.t0); % takes ~2.48e-5 s
% mglFlush
% stimulus.timedebug(7,task.thistrial.framecount+1) = mglGetSecs(stimulus.t0); % takes ~0.007 s ..i.e. compensating for the other steps

%% eye tracking
% track for task
% *** track for fixation???
if (~stimulus.noeye) && any(task.thistrial.thisseg==[1])
    % mouse version for testing with no eyetracker
    if stimulus.eyemousedebug
        mInfo = mglGetMouse(myscreen.screenNumber);
        degx = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
        degy = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;

        pos = [degx, degy];
    else  % check eye pos
        [pos,postime] = mglEyelinkGetCurrentEyePos; % is this in image coordinates?
    end
        
    task.thistrial.trackEye(task.thistrial.framecount,:)  = pos;
    task.thistrial.trackEyeTime(task.thistrial.framecount) = postime;
end

% stimulus.timedebug(8,task.thistrial.framecount+1) = mglGetSecs(stimulus.t0); % takes ~2.86661e-5 s

if stimulus.grabframe && (task.thistrial.thisseg== 1)
    global frame; frame{task.thistrial.framecount} = mglFrameGrab;
end

end

%% Get response; do nothing. 
function [task myscreen] = responseCallback(task, myscreen)

global stimulus
 
end

%% Initialize stimulus
function stimulus = myInitStimulus(stimulus,myscreen,task)  
    % set standard deviation of stimulus
    
    % stimulus size
    if ~isfield(stimulus,'stimStd'), stimulus.stimStd = 0.4;,end %unit: imageX, in deg. 
    stimulus.patchsize = min(6*stimulus.stimStd,min(myscreen.imageWidth,myscreen.imageHeight));
    
    %stimulus initial position. uniform distribution across the screen
    x_img = min(3*stimulus.stimStd,1/3*myscreen.imageWidth)*(2*rand(1)-1); 
    y_img = min(3*stimulus.stimStd,1/3*myscreen.imageWidth)*(2*rand(1)-1);
    stimulus.position = [x_img, y_img];
    
    % stimulus speed
    if ~isfield(stimulus,'stepStd'), stimulus.stepStd = 3/myscreen.framesPerSecond;,end %unit: cm/s to deg/frame
    % this might change based on effective sampling rate.
    
    % stimulus luminance
    if ~isfield(stimulus,'stimLum'), stimulus.stimLum = 122;,end %unit: luminance
            
    % background noise
    if ~isfield(stimulus,'noiseLum'), stimulus.noiseLum = 122;,end; % unit: luminance
    
    % background luminance
    if ~isfield(stimulus,'backLum'), stimulus.backLum = 32;,end; % unit: luminance

    % initialize stimulus
    if isfield(stimulus,'gaussian'), mglDeleteTexture(stimulus.gaussian);, end 
    gaussian    =  mglMakeGaussian(stimulus.patchsize,stimulus.patchsize,...
        stimulus.stimStd,stimulus.stimStd)*(stimulus.stimLum);
    gaussian_rgb           = 255*ones(4,size(gaussian,2),size(gaussian,1),'uint8');
    gaussian_rgb(4,:,:)    = round(gaussian');
    gaussian_rgb           = uint8(gaussian_rgb);

    stimulus.gaussian = mglCreateTexture(gaussian_rgb);
    
    % fixation cross
    stimulus.fixColor = [1 1 1];
end

function stimulus = updateTarget(stimulus,myscreen,task)
    % convert steps to direction gradients
    xstep = normrnd(0,stimulus.stepStd);
    ystep = normrnd(0,stimulus.stepStd);
    
    %[xstep, ystep] = convertNearestPixel(myscreen,xstep,ystep);
    
    stimulus.position = stimulus.position + [xstep, ystep];
    
    % if thre circle goes out of the image 
    if stimulus.position(1)+3*stimulus.stimStd > myscreen.imageWidth/2 ...
            || stimulus.position(1)-3*stimulus.stimStd < -myscreen.imageWidth/2,
        stimulus.position(1) = stimulus.position(1) - 2*xstep;
    end
        
    if stimulus.position(2)+3*stimulus.stimStd > myscreen.imageHeight/2 ...
            || stimulus.position(2)-3*stimulus.stimStd < -myscreen.imageHeight/2,
        stimulus.position(2) = stimulus.position(2) - 2*ystep;
    end
    
    % if eye fixation is enforced
    if stimulus.fixateCenter == 1
        if norm(stimulus.position) < 1 % center of the circle is never in the middle
            stimulus.position(1) = stimulus.position(1) - 2*xstep;
            stimulus.position(2) = stimulus.position(2) - 2*ystep;
        end
    end
    
    %disp(['Stimulus position (x,y): ' num2str(stimulus.position)])
end

function [stimx, stimy] = convertNearestPixel(myscreen,x_img,y_img)
    x_screen = x_img*myscreen.screenWidth/myscreen.imageWidth + myscreen.screenWidth/2;
    y_screen = y_img*myscreen.screenHeight/myscreen.imageHeight + myscreen.screenHeight/2;
    
    stimx = (ceil(x_screen)-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
    stimy = (floor(y_screen)-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight; % what is imagewidth???
end


function precomputBackground(myscreen,stimulus)
% generate a lot of noise
% takes about 12 minutes to generate on the stimulus test computer. 
savefile = '/Users/joshua/data/trackpos/trackpos.mat';

%% noise
% background noise
downsample_timeRes  = 1;
downsample_spatRes  = 10;

ntrials             = 50; %7*9; %7 trials * 9 conditions
nframes             = myscreen.framesPerSecond*30; % 30s; %/downsample_timeRes; 

for idx = 1:ntrials
    tic
    if ~isfield(stimulus,'gaussianFFT')
        xsize_deg          = round(myscreen.imageWidth/downsample_spatRes);
        ysize_deg          = round(myscreen.imageHeight/downsample_spatRes);

        backgaussian = mglMakeGaussian(xsize_deg,ysize_deg,...
            stimulus.stimStd/downsample_spatRes,stimulus.stimStd/downsample_spatRes)*255;
        stimulus.gaussianFFT = getHalfFourier(backgaussian);
    end 

    for idx2 = 1:nframes
        back                        = stimulus.gaussianFFT; %0.02s
        back.phase                  = rand(size(back.mag))*2*pi; % scramble phase % 0.02s
        backgroundnoise             = round(reconstructFromHalfFourier(back));   %0.04s
        if idx2 == 1 % to save time; only allocate memory once.
            backgroundnoise_rgb         = 255*ones(4,size(backgroundnoise,2),size(backgroundnoise,1),nframes,'uint8'); %0.1165 s
        end % otherwise just overwrite.
        backgroundnoise_rgb(4,:,:,idx2)  = backgroundnoise'/max(max(backgroundnoise))*stimulus.noiseLum;  % normalize contrast %0.025s
    end
    % backgroundnoise_rgb = uint8(backgroundnoise_rgb);
    eval(['backgroundnoise_rgb' num2str(idx) '=uint8(backgroundnoise_rgb);']); %0.02s 
    toc
end
clearvars('backgroundnoise_rgb')
save(savefile, 'backgroundnoise_rgb*','-v7.3')

end
