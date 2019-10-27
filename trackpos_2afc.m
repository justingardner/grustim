%
%        $Id: $
%      usage: trackpos_2afc
%         by: Josh Ryu
%       date: 10/27/2019
%    purpose: 2 alternative forced choice task for position of two blobs
%

% Uses staircase procedures. 
% S1: stimulus period (0.5s)
% S2: fixation period (0.5s)
% S3: stimulus 2 period (0.5s)
% S4: repsonse period (1.5s)
% S5: random period of fixation (1s)

function myscreen = trackpos_2afc(varargin)

%% set up screen and experiment
% set input arguments
if isempty(mglGetSID)
    myscreen.subjectID  = -1;
else
    myscreen.subjectID  = mglGetSID;
    myscreen.saveData = 1;
end

myscreen.displayName = 'debug'; myscreen.screenNumber = 1; 
myscreen.screenWidth = 860; myscreen.screenHeight = 600; 
myscreen.hideCursor = 1;
myscreen = initScreen(myscreen);

% Experimenter parameters
grabframe       = 1; 

%%
% Go straight to task.
% S1: stimulus period (0.5s)
% S3: fixation period (0.5s)
% S3: stimulus 2 period (0.5s)
% S4: repsonse period (2s)
% S5: random period of fixation (1~3s)
task{1}{1}.segmin = [0.5 0.5 0.5 1.5 1];
task{1}{1}.segmax = [0.5 0.5 0.5 1.5 2];
task{1}{1}.numTrials = 150;
task{1}{1}.getResponse = [0 0 0 1 0]; %segment to get response.
task{1}{1}.waitForBacktick = 1; %wait for backtick before starting each trial 

% Run fixed intervals (1) or staircase (0)
task{1}{1}.runfixedint  = 0; % run staircase
task{1}{1}.blankrun     = 1;

% stimulus and background
task{1}{1}.parameter.backLum    = 32; %32;  % background luminance; units: fraction of full luminance 
task{1}{1}.parameter.noiseLum   = 0; % noise luminance, if there is one.
task{1}{1}.parameter.stimLum    = 255 - task{1}{1}.parameter.backLum;  % stimulus luminance (out of 255)
teststimLum   = linspace(task{1}{1}.parameter.stimLum, task{1}{1}.parameter.noiseLum,3);
task{1}{1}.parameter.stimLum    = teststimLum(2); % choose medium condition.
teststimDur   = [0.3 0.2 0.1];

%task parameters
task{1}{1}.randVars.calculated.posDiff      = nan;
task{1}{1}.randVars.calculated.firstcorr    = nan; %store values calculated during the task. 
task{1}{1}.randVars.calculated.firstpos     = nan;
task{1}{1}.randVars.calculated.secondpos    = nan;
task{1}{1}.randVars.calculated.subjcorrect  = nan; %store values calculated during the task. 
% task{1}{1}.randVars.calculated.lumCurr      = nan; %??? does it not save the parameters of interest automatically? 

%% task blocks. change stimulus duration.
for durIdx = 1:length(teststimDur)
    phaseN = durIdx;
    stimDur = teststimDur(durIdx);
    % adaptation condition, at full luminance
    task{1}{phaseN}             = task{1}{1};  
    task{1}{phaseN}.segmin      = [stimDur 0.5 stimDur 1.5 1];
    task{1}{phaseN}.segmax      = [stimDur 0.5 stimDur 1.5 2];    
end

% add a blank run with 10 trials
if task{1}{1}.blankrun == 1
    phaseN = length(teststimDur)+1;
    task{1}{phaseN} = task{1}{1};
    task{1}{phaseN}.numTrials = 10;
    task{1} = task{1}([phaseN, 1:phaseN-1]);
end

% initializing task...
disp(' Initializing Task....')

for phaseN = 1:length(task{1})
    [task{1}{phaseN} myscreen] = initTask(task{1}{phaseN},myscreen,...
        @startSegmentCallback,@screenUpdateCallback,@responseCallback,@initTrialCallback);
end

%% initialize stimulus
global stimulus;
stimulus = [];

if task{1}{1}.runfixedint == 0
    stimulus.stairUp        = 1;
    stimulus.stairDown      = 2;
    stimulus.stairStepSize  = 0.05;
    stimulus.stairUseLevitt = 0;
    stimulus.stairUsePest   = 1;    % use PEST
    stimulus.stairRep       = 50;   % repeat staircase every [stairRep] trials
    stimulus.stairN         = 0;    % keeps track of how many staircases they did
    stimulus.threshold      = 1;    % starting threshold?
else % run fixed intervals (set the values here) 
end

stimulus = myInitStimulus(stimulus,myscreen,task);
myscreen = initStimulus('stimulus',myscreen); % what does this do???

%% run the task
stimulus.grabframe = grabframe; %save frames into matrices for outputting task image
if stimulus.grabframe, global frame; frame = {};, end

mglDisplayCursor(0); %hide cursor
mglClearScreen(task{1}{1}.parameter.backLum/255);
mglTextDraw('task (trackpos_2afc) starting... ', [0 0.5])
mglTextDraw('Press 1 if the first stimulus is to the left of the second stimulus',[0 -0.5]);
mglFlush

phaseNum = 1;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);     % update the task
    myscreen = tickScreen(myscreen,task);     % flip screen
end

myscreen = endTask(myscreen,task);
mglClose 
endScreen(myscreen); mglDisplayCursor(1) %show cursor

if stimulus.grabframe
    save('/Users/jryu/data/trackpos_2afc/taskSetup/frame.mat', 'frame')
    % save('/Users/joshryu/Dropbox/GardnerLab/data/trackpos_2afc/taskSetup/frame.mat', 'frame')
end
end

%% Initialize trials; set staircuse
function [task myscreen] = initTrialCallback(task, myscreen)
    global stimulus
    
    stimulus.stimLum    = task.thistrial.stimLum;
    stimulus.backLum    = task.thistrial.backLum;
    stimulus.noiseLum   = task.thistrial.noiseLum;
    
    task.thistrial.framecount = 0;
    
    if task.runfixedint == 0
        if task.trialnum == 1 %does the task 
            stimulus.stairN = 0;
        end

        %initialize staircase  
        if mod(stimulus.stairN, stimulus.stairRep) == 0
            stimulus = initStaircase(stimulus);
        end

        [stimulus.threshold, stimulus.staircase] = doStaircase('testValue',stimulus.staircase); %left threshold
    %     [stimulus.threshold(2) stimulus.staircase(2)] = doStaircase('testValue',stimulus.staircase(2)); %right threshold
    else % run fixed interval
        
    end
    
    if stimulus.grabframe
        global frame
        %save('/Users/jryu/Dropbox/Stanford/Current/FYP/FYP talk/figures/trackpos_2afc_32back_500ms.mat', 'frame','-v7.3')
        %save('/Users/jryu/Dropbox/Stanford/Current/FYP/FYP talk/figures/trackpos_2afc_32back_100ms.mat', 'frame','-v7.3')
        frame = {};
        frame{sum(task.segmax)*myscreen.framesPerSecond} = [];
    end

end

function stimulus = initStaircase(stimulus)
% set up left and right staircase
if stimulus.stairUseLevitt
    stimulus.staircase = doStaircase('init','upDown','nup',stimulus.stairUp,...
        'ndown',stimulus.stairDown,'initialThreshold',stimulus.threshold,...
        'initialStepsize',stimulus.stairStepSize,'testType=levitt','minThreshold',stimulus.stairStepSize,'maxThreshold',45); % maximum has to be 45. 
%     stimulus.staircase(2) = doStaircase('init','upDown','nup',stimulus.stairUp,...
%         'ndown',stimulus.stairDown,'initialThreshold',stimulus.threshold(2),...
%         'initialStepsize',stimulus.stairStepSize,'testType=levitt','minThreshold',stimulus.stairStepSize,'maxThreshold',45); % maximum has to be 45. 
elseif stimulus.stairUsePest
    stimulus.staircase = doStaircase('init','upDown','nup',stimulus.stairUp,...
        'ndown',stimulus.stairDown,'initialThreshold',stimulus.threshold,...
        'initialStepsize',stimulus.stairStepSize,'testType=pest','minThreshold',stimulus.stairStepSize,'maxThreshold',45);
%     stimulus.staircase(2) = doStaircase('init','upDown','nup',stimulus.stairUp,...
%         'ndown',stimulus.stairDown,'initialThreshold',stimulus.threshold(2),...
%         'initialStepsize',stimulus.stairStepSize,'testType=pest','minThreshold',stimulus.stairStepSize,'maxThreshold',45);
else
    stimulus.staircase = doStaircase('init','upDown','nup',stimulus.stairUp,...
        'ndown',stimulus.stairDown,'initialThreshold',stimulus.threshold,...
        'initialStepsize',stimulus.stairStepSize,'minThreshold',0,'maxThreshold',45);
%     stimulus.staircase(2) = doStaircase('init','upDown','nup',stimulus.stairUp,...
%         'ndown',stimulus.stairDown,'initialThreshold',stimulus.threshold(2),...
%         'initialStepsize',stimulus.stairStepSize,'minThreshold',0,'maxThreshold',45);
end

end

%% Start segment
function [task myscreen] = startSegmentCallback(task, myscreen)
% S1: stimulus period (0.5s)
% S2: repsonse period (2s)
% S3: random period of fixation (1~3s)
global stimulus 

%change stimulus accordingly
%find stimulus position for left and right
if task.thistrial.thisseg == 1
    % after the noise bc we want to recreate the texture with the seed.
    stimulus = myInitStimulus(stimulus,myscreen,task); %centerX,Y, diameter called by getArgs.
        
    %stimulus initial position. uniform distribution somewhere middle of screen
    x_img = (2*stimulus.stimStd)*(2*rand(1)-1); 
    
    task.thistrial.firstpos     = x_img;
    task.thistrial.firstcorr    = (rand(1)<0.5);   %1 if correct side is on the left. 0 if left. 
    task.thistrial.posDiff      = (2*task.thistrial.firstcorr-1)*stimulus.threshold;
    task.thistrial.secondpos    = x_img + task.thistrial.posDiff; %x_img is higher (more right) in the second if left is the correct one. 
        
elseif task.thistrial.thisseg == 4
    stimulus.fixColor = stimulus.fixColors.response;
else
   % no fixation (it becomes a cue)
end    

end

%% screen update
function [task myscreen] = screenUpdateCallback(task, myscreen)
global stimulus % call stimulus

mglClearScreen(stimulus.backLum/255);

task.thistrial.framecount = task.thistrial.framecount + 1;

% draw blob
if task.thistrial.thisseg == 1   
    
    mglBltTexture(stimulus.gaussian,[task.thistrial.firstpos 0]);

elseif task.thistrial.thisseg == 3
    mglBltTexture(stimulus.gaussian,[task.thistrial.secondpos 0]);
elseif task.thistrial.thisseg == 2
    % no fixation.
elseif task.thistrial.thisseg == 4 %response feedback
    mglGluAnnulus(0,0,0.5,0.75,stimulus.fixColor,60,1);
end

if stimulus.grabframe
    global frame; frame{task.thistrial.framecount} = mglFrameGrab;
end

end

%% Get response 
function [task myscreen] = responseCallback(task, myscreen)

global stimulus

% record responses. correct/incorrect
if any(task.thistrial.whichButton == [1 2])
    respIs1 = (task.thistrial.whichButton == 1); %1 if the subject chose left
    correct = (task.thistrial.firstcorr == respIs1); % correct if first is correct and response is 1.
    task.thistrial.subjcorrect = correct;
    
    if task.runfixedint == 0
        stimulus.stairN = stimulus.stairN+1; %count how many times 
    end
else
    stimIs1 = task.thistrial.firstcorr; respIs1 = nan; correct = nan;
end

% change color of fixation for feedback.  
if isnan(correct)
    stimulus.fixColor = [1 1 1];
elseif correct
    stimulus.fixColor = [0 1 0];
else
    stimulus.fixColor = [1 0 0];
end

% Output response to the screen. 
if task.thistrial.whichButton == 1, respSide = 'seg1';
elseif task.thistrial.whichButton == 2, respSide = 'seg2'; end
if correct == 0, corrString = 'incorrect';
elseif correct == 1, corrString = 'correct';
else, corrString = 'no response'; end

%stimulus.staircase(stimulus.leftcorrect+1) = doStaircase('update',stimulus.staircase(stimulus.leftcorrect+1),correctIncorrect,abs(task.thistrial.dirDiff));
%[stimulus.threshold(stimulus.leftcorrect+1), stimulus.staircase(stimulus.leftcorrect+1)] = doStaircase('testValue',stimulus.staircase(stimulus.leftcorrect+1));

if task.runfixedint == 0
    stimulus.staircase = doStaircase('update',stimulus.staircase,correct,abs(task.thistrial.posDiff));
    % [stimulus.threshold, stimulus.staircase] =
    % doStaircase('testValue',stimulus.staircase); % I do this in the
    % beginning of segment
end

disp(['Position difference: ' num2str(task.thistrial.posDiff) '; ' ...
    'Position: ' num2str(task.thistrial.firstpos) ' (seg1) vs ' num2str(task.thistrial.secondpos) ' (seg2); ' ...
    'Response: ' respSide '; ' corrString])

end

%% Initialize stimulus, initialize blob
function stimulus = myInitStimulus(stimulus,myscreen,task)  
    % set standard deviation of stimulus
    
    % stimulus size
    if ~isfield(stimulus,'stimStd'), stimulus.stimStd = 2;,end %unit: imageX, in deg. 
    stimulus.patchsize = min(6*stimulus.stimStd,min(myscreen.imageWidth,myscreen.imageHeight));
    
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
    stimulus.fixColors.response = [1 1 1];

end

