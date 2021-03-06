%        $Id: $
%      usage: trackpos
%         by: Josh Ryu
%       date: 05/05/2021
%    purpose: 

% todo: make noiseLum more flexible

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
%todo:  check these throughout the code!!
exp.noeye               = 0; % 1 if no eyetracking (mouse for eye); 0 if there is eye tracking `
exp.showmouse           = 0; 
exp.grabframe           = 0; 
exp.whitenoiseOn        = 0; % 1: white noise; 2: 
exp.fixateCenter        = 1;
exp.phasescrambleOn     = 1;
exp.backprecompute      = 1;
exp.eyemousedebug       = 0; % debug eyetracker with mouse
exp.debug               = 0; % debug code
exp.downsample_timeRes  = 1;


% Task design (might be changed later, so check this)
% S1: Stimulus (30s) 
% S2: Fixation (3s)
task{1}{1}.segmin           = [30 0.1]; % fixation for shorter bc of the segment start takes time.
task{1}{1}.segmax           = [30 0.1]; 
task{1}{1}.numTrials        = 5; % changed later depending on the condition
task{1}{1}.getResponse      = [0 0]; %segment to get response.
task{1}{1}.synchToVol       = [0 1]; %segmet to wait for backtick
task{1}{1}.waitForBacktick  = 1; %wait for backtick before starting task
task{1}{1}.random           = 1;

% task parameters for adaptation conditions
task{1}{1}.random               = 1;
if exp.whitenoiseOn == 1 || exp.phasescrambleOn == 1
    task{1}{1}.parameter.phasescrambleOn    = 1;
    task{1}{1}.parameter.backLum            = 90; %160;%90;  % background luminance; units: luminance 
    task{1}{1}.parameter.noiseLum           = 32;
else 
    task{1}{1}.parameter.backLum = 90;  % background luminance; units: fraction of full luminance 
end
% task{1}{1}.parameter.stimLum = 255 - task{1}{1}.parameter.backLum;  % stimulus luminance (out of 255)


% calculated parameters
maxframes = ceil(task{1}{1}.segmax(1)*myscreen.framesPerSecond) + 20;
task{1}{1}.randVars.calculated.randomSeed   = nan;
task{1}{1}.randVars.calculated.bgpermute    = nan(1,maxframes); % nframes x 1 for the background
task{1}{1}.randVars.calculated.perm         = nan(maxframes,1);
task{1}{1}.randVars.calculated.initStim     = [nan nan];
task{1}{1}.randVars.calculated.trackStim    = nan(maxframes,2);
task{1}{1}.randVars.calculated.trackResp    = nan(maxframes,2);
task{1}{1}.randVars.calculated.trackEye     = nan(maxframes,2);
task{1}{1}.randVars.calculated.trackTime    = nan(1,maxframes);
task{1}{1}.randVars.calculated.trackEyeTime = nan(1,maxframes); % for referencing edf file

%% Set up tasks 

% change stimulus speed and luminance; cross conditions.
teststimSteps = [1]; %[0.75, 1.5, 2.25];
teststimLum   = task{1}{1}.parameter.noiseLum*[0.5, 1, 1.5, 2]; %SNR
% teststimLum   = linspace(task{1}{1}.parameter.stimLum, task{1}{1}.parameter.noiseLum,3);

for stepIdx = 1:length(teststimSteps)
    idx = (stepIdx-1)*(2) + 1;
    % adaptation condition, at full luminance
    task{1}{idx}                                = task{1}{1};  
    task{1}{idx}.parameter.stimStep             = teststimSteps(stepIdx);
    task{1}{idx}.parameter.phasescrambleOn      = 0;
    task{1}{idx}.parameter.noiseLum             = 0;  
    task{1}{idx}.parameter.stimLum              = 32;
    task{1}{idx}.segmin         = [30 3]; %fixation time constrained by the texture loading
    task{1}{idx}.segmax         = [30 3]; 
    task{1}{idx}.numTrials      = 5;
    
    % adaptation task
    idx = (stepIdx-1)*(2) + 2; %cycle through the luminances first.

    task{1}{idx}                            = task{1}{1};  
    task{1}{idx}.parameter.phasescrambleOn  = 1;
    task{1}{idx}.parameter.noiseLum         = 32;
    task{1}{idx}.parameter.stimStep         = teststimSteps(stepIdx);
    task{1}{idx}.parameter.stimLum          = teststimLum;
    task{1}{idx}.segmin                     = [30 3]; %fixation time constrained a bit by the texture loading
    task{1}{idx}.segmax                     = [30 3]; 
    task{1}{idx}.numTrials                  = 10 * length(teststimLum);
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
stimulus.exp = exp;

if stimulus.exp.phasescrambleOn == 1;
    disp('Loading phase scrambled background noise...')

    tic
    if stimulus.exp.backprecompute == 1;
        savefile = '/Users/gru/data/trackpos/trackpos.mat';
        % savefile            = '/Users/joshua/data/trackpos_2afc/trackpos.mat'; % just use noise 1 and permute
        if ~exist(savefile,'file')
            error('need background file')
        end
        
        load(savefile,'backgroundnoise_rgb');

        if isfield(stimulus,'backnoise')
            for idx = 1:length(stimulus.backnoise)
                mglDeleteTexture(stimulus.backnoise{idx});
            end
        end

        % create all background textures and then load them later
        if stimulus.exp.debug 
            nnn = 200;
        else
            nnn = size(backgroundnoise_rgb,4);
        end
        for idx = 1:nnn %too big?? memory?
            stimulus.backnoise{idx} = mglCreateTexture(backgroundnoise_rgb(:,:,:,idx));
        end

        clearvars('backgroundnoise_rgb')
    end
    toc
end

if stimulus.exp.grabframe
    global frame
    frame = {};
end

%% Eye calibration
if ~exp.noeye
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

if ~exp.showmouse, mglDisplayCursor(0);, end %hide cursor

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

if stimulus.exp.grabframe
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
    if stimulus.exp.phasescrambleOn == 1 && stimulus.exp.backprecompute == 1
        nframes = myscreen.framesPerSecond*task.segmax(1) + 20;%/downsample_timeRes; 
        task.thistrial.bgpermute(1:nframes) = randi(length(stimulus.backnoise),nframes,1);
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
    if ~stimulus.exp.showmouse, mglDisplayCursor(0);, end %hide cursor
    
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
    
    if stimulus.exp.grabframe
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
    % inject noise
    if task.thistrial.phasescrambleOn == 1 && mod(task.thistrial.framecount, stimulus.exp.downsample_timeRes) == 0
        idx = task.thistrial.bgpermute(task.thistrial.framecount);
        mglBltTexture(stimulus.backnoise{idx},...
            [0 0 myscreen.imageWidth myscreen.imageHeight])
    end
    
    % draw stimulus
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

    if stimulus.exp.fixateCenter == 1
        mglGluAnnulus(0,0,0.2,0.3,stimulus.fixColor,60,1);
        mglGluDisk(0,0,0.1,rand(1,3),60,1);
    end
    
elseif (task.thistrial.thisseg == 2)
    if stimulus.exp.fixateCenter == 1 % stop the flashing
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
if (~stimulus.exp.noeye) && any(task.thistrial.thisseg==[1])
    % mouse version for testing with no eyetracker
    if stimulus.exp.eyemousedebug
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

if stimulus.exp.grabframe && (task.thistrial.thisseg== 1)
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
    if ~isfield(stimulus,'stimStd'), stimulus.stimStd = 1;,end %unit: imageX, in deg. 
    %GardnerLab: stimstd = 2; CSNL stimStd = 0.4.. (why...?)
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
    if stimulus.exp.fixateCenter == 1
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
savefile = '/Users/gru/data/trackpos/trackpos.mat';

%% noise
% background noise
downsample_timeRes  = 1;
downsample_spatRes  = 10;

%ntrials             = 1; %7*9; %7 trials * 9 conditions
nframes             = 10000; %myscreen.framesPerSecond*30; % 30s; %/downsample_timeRes; 

%for idx = 1:ntrials
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
backgroundnoise_rgb = uint8(backgroundnoise_rgb);
%eval(['backgroundnoise_rgb =uint8(backgroundnoise_rgb);']); %0.02s 
toc
%end
save(savefile, 'backgroundnoise_rgb*','-v7.3')

end
