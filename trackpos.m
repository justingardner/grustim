%
%        $Id: $
%      usage: trackpos
%         by: Josh Ryu
%       date: 05/15/2019
%    purpose: 

function myscreen = trackpos(varargin)

% set input arguments
%getArgs(varargin,{'subjectID=s999','centerX=10','centerY=0','diameter=16'});
getArgs(varargin,{'subjectID=-1'});

% set up screen
myscreen.subjectID = subjectID;
myscreen.saveData = 1;
myscreen = initScreen(myscreen);

%% parameters

% Experimenter parameters
eyewindow       = 0; %whats this?
noeye           = 1; % 1 if no eyetracking (mouse for eye); 0 if there is eye tracking 
eyemousedebug   = 0; % do i need this? 
grabframe       = 0; 
whitenoiseOn    = 0;

% Task design
% S1: Stimulus (30s)
% S2: Fixation (3s)
task{1}{1}.segmin = [30 3]; 
task{1}{1}.segmax = [30 3]; 
task{1}{1}.numTrials = 10;
task{1}{1}.getResponse = [1 0]; %segment to get response.
task{1}{1}.waitForBacktick = 1; %wait for backtick before starting each trial 

% task parameters
task{1}{1}.parameter.backLum = 0.5;  % background luminance; units: fraction of full luminance 
task{1}{1}.parameter.stimLum = 122;  % stimulus luminance (out of 255)
task{1}{1}.parameter.stimStep = 4;   % stimulus velocity in cm/sec

% The pilot test has three main parts:
% 1. motor gain estimation through testing different stimulus speed
% 2. adaptation trials so that the subjects learn priors
% 3. test different uncertainty values.

% calculated parameters
task{1}{1}.randVars.calculated.initStim = [nan nan];
task{1}{1}.randVars.calculated.trackStim = nan(ceil(task{1}{1}.segmax(1)*myscreen.framesPerSecond),2);
task{1}{1}.randVars.calculated.trackResp = nan(ceil(task{1}{1}.segmax(1)*myscreen.framesPerSecond),2);
task{1}{1}.randVars.calculated.trackEye  = nan(ceil(task{1}{1}.segmax(1)*myscreen.framesPerSecond),2);
task{1}{1}.randVars.calculated.trackTime = nan(ceil(task{1}{1}.segmax(1)*myscreen.framesPerSecond),1);

%% initialize stuff 
% initialize task
disp(' Initializing Task....')

% The pilot test has three main parts:
% 1. motor gain estimation through testing different stimulus speed
% 2. adaptation trials so that the subjects learn priors
% 3. test different uncertainty values.
phaseN = 1; %motor gain estimation
task{1}{phaseN} = task{1}{1};
task{1}{phaseN}.numTrials = 30;
task{1}{phaseN}.parameter = rmfield(task{1}{phaseN}.parameter,'stimStep');
task{1}{phaseN}.randVars.uniform.stimStep = [2,4,6,8,10,12];

phaseN = 2; %adaptation trialsgegetg
task{1}{phaseN} = task{1}{1};

teststimLum = [120:-20:20]; % main task 
for phaseN = 3:(2+length(teststimLum)) %adaptation trials
    task{1}{phaseN} = task{1}{1};
    task{1}{phaseN}.parameter.stimLum = teststimLum(phaseN-2);
end

for phaseN = 1:(2+length(teststimLum))
    [task{1}{phaseN} myscreen] = initTask(task{1}{phaseN},myscreen,...
        @startSegmentCallback,@screenUpdateCallback,@responseCallback,@initTrialCallback);
end

% initialize stimulus
disp(' Initializing Stimulus....')
global stimulus; stimulus = struct;

myscreen = initStimulus('stimulus',myscreen); % what does this do???
stimulus.whitenoiseOn = whitenoiseOn;
stimulus.eyewindow  = eyewindow;
stimulus.noeye      = noeye;
stimulus.eyemousedebug = eyemousedebug;

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
disp(' Running Task....')

phaseNum = 1;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);     % update the task
    myscreen = tickScreen(myscreen,task);     % flip screen
end


disp(' End Task....')

myscreen = endTask(myscreen,task);
mglClose
mglDisplayCursor(1)

if stimulus.grabframe
    save('/Users/joshryu/Dropbox/GardnerLab/data/FYP/trackpos/frame_nored.mat', 'frame')
end

end

%% Initialize trials 
function [task myscreen] = initTrialCallback(task, myscreen)
    global stimulus    
       
    stimulus.stimStep = task.thistrial.stimStep;
    stimulus.stimLum  = task.thistrial.stimLum;
    stimulus.backLum  = task.thistrial.backLum;
    
    stimulus = myInitStimulus(stimulus,myscreen,task); %centerX,Y, diameter called by getArgs.

    % mglDisplayCursor(0) %hide cursor

end

%% Start segment
function [task myscreen] = startSegmentCallback(task, myscreen)
% S1: Stimulus (30s)
% S2: Fixation (3s)
global stimulus 

%change stimulus accordingly
if task.thistrial.thisseg == 1    
    % task{1}{1}.randVars.calculated.initStim = [nan nan];
    
    %% stimulus
    % initialize stimulus
    stimulus.gaussian = mglMakeGaussian(stimulus.circlerad,stimulus.circlerad,...
        stimulus.stimStd,stimulus.stimStd)*(stimulus.stimLum) + 255*stimulus.backLum;
    
    %stimulus initial position. uniform distribution across the screen
    x_img = myscreen.imageWidth/2*(rand(1)-0.5); 
    y_img = myscreen.imageHeight/2*(rand(1)-0.5);
    
    x_screen = x_img*myscreen.screenWidth/myscreen.imageWidth + myscreen.screenWidth/2;
    y_screen = y_img*myscreen.screenHeight/myscreen.imageHeight + myscreen.screenHeight/2;
    
    % set mouse position to the stimulus direction. 
    mglSetMousePosition(ceil(x_screen),floor(y_screen), myscreen.screenNumber); % correct for screen resolution???

    % convert screen coordinates to image coordinates.
    stimx = (ceil(x_screen)-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
    stimy = (floor(y_screen)-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight; % what is imagewidth???
    
    task.thistrial.initStim = [stimx, stimy];
    stimulus.position = [stimx, stimy];
    
    %% noise
    if stimulus.whitenoiseOn == 1
        nosietime = size(stimulus.backgroundnoise,4);
        repeatnoise = ceil((task.segmax(1)*myscreen.framesPerSecond)/nosietime);
        task.thistrial.noiseperm = [];
        for rp = 1:repeatnoise
            task.thistrial.noiseperm = [task.thistrial.noiseperm randperm(nosietime)]; % change to random sampling with replacement???
        end
    end
    %% frame counter.
    task.thistrial.framecount = 0;
    
else %intertrial interval
    % *** fixation cross.
    
end    

end

%% screen update
function [task myscreen] = screenUpdateCallback(task, myscreen)
% S1: Stimulus (30s)
% S2: Fixation (3s)

%% Update Screen
global stimulus % call stimulus
mglClearScreen(stimulus.backLum);

if (task.thistrial.thisseg== 1)
    task.thistrial.framecount = task.thistrial.framecount + 1;

    % myscreen = initScreen(myscreen);mglClearScreen

    % background noise
    if stimulus.whitenoiseOn == 1
        backgroundnoise = stimulus.backgroundnoise(:,:,:,task.thistrial.noiseperm(task.thistrial.framecount));
        texture         = mglCreateTexture(backgroundnoise);
    
    %{
    backgroundnoise             = round(rand(myscreen.screenWidth/2,myscreen.screenHeight/2)*stimulus.noiseLum);  
    backgroundnoise_rgb         = 255*ones(4,size(backgroundnoise,1),size(backgroundnoise,2),size(backgroundnoise,3),'uint8');
    backgroundnoise_rgb(4,:,:,:)= backgroundnoise;
    backgroundnoise_rgb         = uint8(backgroundnoise_rgb);
    
    texture         = mglCreateTexture(backgroundnoise_rgb);
    
    
    mglClearScreen;
    mglBltTexture(texture,[0,0,myscreen.imageWidth,myscreen.imageHeight]);
    mglFlush
    %}
    end

    % stimulus
    stimulus        = updateTarget(stimulus,myscreen,task); % update position.
    texture         = mglCreateTexture(stimulus.gaussian);
    mglBltTexture(texture,stimulus.position);

    % **&display mouse position
    mInfo = mglGetMouse(myscreen.screenNumber);

    mimg_x = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
    mimg_y = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;

    mglGluDisk(mimg_x, mimg_y, 0.1, [1 0 0])

    % ***record stimulus position and mouse position  
    task.thistrial.trackStim(task.thistrial.framecount,:) = stimulus.position;
    task.thistrial.trackResp(task.thistrial.framecount,:) = [mimg_x, mimg_y];
    task.thistrial.trackTime(task.thistrial.framecount)   = rem(now,1)*24*60*60; %serial date time in seconds 

elseif (task.thistrial.thisseg == 2)
    % draw fixation;
    mglGluAnnulus(0,0,0.5,0.75,stimulus.fixColor,60,1);
end

mglFlush

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
        [pos,~] = mglEyelinkGetCurrentEyePos; % is this in image coordinates?
    end
        
    task.thistrial.trackEye(task.thistrial.framecount,:)  = pos;
end

if stimulus.grabframe
    global frame; frame{task.thistrial.thisseg} = mglFrameGrab;
end

end

%% Get response; do nothing. 
function [task myscreen] = responseCallback(task, myscreen)

global stimulus
 
end

%% Initialize stimulus
function stimulus = myInitStimulus(stimulus,myscreen,task)  
    tic
    % set standard deviation of stimulus
    
    % stimulus size
    if ~isfield(stimulus,'stimStd'), stimulus.stimStd = 1;,end %unit: imageX, in cm. 
    stimulus.circlerad = 10*stimulus.stimStd;
    stimulus.position = [0 0]; 
    
    % stimulus speed
    if ~isfield(stimulus,'stepStd'), stimulus.stepStd = 4/myscreen.framesPerSecond;,end %unit: cm/s
    % this might change based on effective sampling rate.
    
    % stimulus luminance
    if ~isfield(stimulus,'stimLum'), stimulus.stimLum = 122;,end %unit: luminance
            
    % background noise
    if ~isfield(stimulus,'noiseLum'), stimulus.noiseLum = 122;,end; % unit: luminance
    
    % background luminance
    if ~isfield(stimulus,'backLum'), stimulus.backLum = 0.5;,end; % unit: percentage full luminance

    if stimulus.whitenoiseOn == 1;
        load('trackpos.mat') % *** still very slow....
        stimulus.gaussian = gaussian_rgb;
        stimulus.backgroundnoise = backgroundnoise_rgb; %round(rand(myscreen.screenHeight-1,myscreen.screenWidth-1,50)*stimulus.noiseLum);
    else
%         load('trackpos.mat','gaussian_rgb') % *** still very slow....
%         stimulus.gaussian = gaussian_rgb;
        stimulus.gaussian = mglMakeGaussian(stimulus.circlerad,stimulus.circlerad,...
            stimulus.stimStd,stimulus.stimStd)*(stimulus.stimLum-stimulus.backLum) + 255*stimulus.backLum; 
    end    
    
    % fixation cross
    stimulus.fixColor = [1 1 1];
    toc
end

function stimulus = updateTarget(stimulus,myscreen,task)
    % convert steps to direction gradients
    xstep = normrnd(0,stimulus.stepStd);
    ystep = normrnd(0,stimulus.stepStd);
    
    %[xstep, ystep] = convertNearestPixel(myscreen,xstep,ystep);
    
    stimulus.position = stimulus.position + [xstep, ystep];
    
    % if thre circle goes out of the image 
    if stimulus.position(1)+stimulus.circlerad > myscreen.imageWidth/2 ...
            || stimulus.position(1)-stimulus.circlerad < -myscreen.imageWidth/2,
        stimulus.position(1) = stimulus.position(1) - 2*xstep;
    end
        
    if stimulus.position(2)+stimulus.circlerad > myscreen.imageHeight/2 ...
            || stimulus.position(2)-stimulus.circlerad < -myscreen.imageHeight/2,
        stimulus.position(2) = stimulus.position(2) - 2*ystep;
    end
    
    %`disp(['Stimulus position (x,y): ' num2str(stimulus.position)])
end

function [stimx, stimy] = convertNearestPixel(myscreen,x_img,y_img)
    x_screen = x_img*myscreen.screenWidth/myscreen.imageWidth + myscreen.screenWidth/2;
    y_screen = y_img*myscreen.screenHeight/myscreen.imageHeight + myscreen.screenHeight/2;
    
    stimx = (ceil(x_screen)-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
    stimy = (floor(y_screen)-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight; % what is imagewidth???
end

%% only used for generating noise beforehand... the task does not run this
% Takes about 5 minutes to run. 
% Make sure the stimulus parameters are consistent with the version of the
% task you are running. 
function backgroundnoise = generatebackgroundnoise(myscreen,stimulus)
    % run with debugger after initializations, before loading the
    % background in myInitStimulus
    % initialize background noise before trial;
    % each frame of background noise takes about 0.05s to generate;
    % i.e. 3s to generate 1s of background noise, and 3 minutes to generate
    % 30s of noise.
    tic 
    
    nframes = 60*10; 
    downsample_screenRes        = 2;
    
    backgroundnoise             = round(rand(myscreen.screenWidth/downsample_screenRes,myscreen.screenHeight/downsample_screenRes,nframes)*stimulus.noiseLum);  
    backgroundnoise_rgb         = 255*ones(4,size(backgroundnoise,1),size(backgroundnoise,2),size(backgroundnoise,3),'uint8');
    backgroundnoise_rgb(4,:,:,:)= backgroundnoise;
    backgroundnoise_rgb         = uint8(backgroundnoise_rgb);
    
    gaussian               = mglMakeGaussian(stimulus.circlerad,stimulus.circlerad,stimulus.stimStd,stimulus.stimStd)*stimulus.stimLum; 
    gaussian_rgb           = 255*ones(4,size(gaussian,1),size(gaussian,2),'uint8');
    gaussian_rgb(4,:,:)    = round(gaussian);
    gaussian_rgb           = uint8(gaussian_rgb);

    save('trackpos.mat', 'backgroundnoise_rgb','gaussian_rgb','stimulus','-v7.3');
    toc

end
