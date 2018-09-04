function [ myscreen ] = afcom( varargin )
% ***scan info: 736 TRs (6 minute * 120 + 16)***
% *** set 'noeye=1' to turn of the eye tracker***
%
%Attention field and color mapping (afcom)
%
%   Map the effects of spatial attention and attention to specific colors 
%%

global stimulus fixStimulus

stimulus = struct;
fixStimulus = struct;

%% Initialize Variables

% add arguments later
scan = 0;
plots = 0;
noeye = 0;
debug = 0;
replay = 0;
powerwheel = 0;
run = 0; 
eyewindow=0; 
mouse=0; 
practice=0; 

getArgs(varargin,{'scan=0','plots=0','noeye=0','powerwheel=0','eyewindow=1.5','practice=0','debug=0','replay=0','attend=1','run=0','build=0','mouse=0'});
stimulus.scan = scan;
stimulus.plots = plots;
stimulus.noeye = noeye;
stimulus.practice = practice;
stimulus.mousedebug = mouse;
stimulus.powerwheel = powerwheel;
stimulus.eyewindow = eyewindow;
stimulus.debug = debug;
stimulus.replay = replay;
stimulus.overrideRun = run;

clear localizer invisible scan noeye task test2 attend build eyewindow mouse practice powerwheel

if stimulus.scan
    warning('Disabling eyewindow');
    stimulus.eyewindow=0;
end

%% Replay mode
% if any(replay>0)
%     if ischar(replay)
%         % a file was called for, load it
%         loaded = load(replay);
%         stimulus = loaded.stimulus;
%         % check that this is actually an afmap file
%         if isempty(strfind(loaded.task{1}{1}.taskFilename,'afmap'))
%             disp(sprintf('File %s is not an afmap run.',replay));
%             return
%         end
%         % get the task parameters so that you can sync the replay correctly
%         disp('GET TASK PARAMETERS');
%         e = getTaskParameters(loaded.myscreen,loaded.task);
%         e1 = e{1};
%         % pull out the trial volumes
%         stimulus.tVolumes = e{1}.trialVolume;
%         disp('Stimulus volumes were found at:');
%         disp(stimulus.tVolumes);
%         
%         stimulus.replayFile = strcat(replay(1:(strfind(replay,'.mat')-1)),'_replay.mat');
%         stimulus.replay = true;
%     else
%         % do an entire subject
%         disp('Replay mode initiated');
%         folder = input('What folder would you like to replay? [Folder]: ');
%         files = dir(fullfile(folder,'*.mat'));
%         for fi = 1:length(files)
%             if isempty(strfind(files(fi).name,'replay')) && isempty(strfind(files(fi).name,'original'))
%                 afmap(sprintf('replay=%s/%s',folder,files(fi).name));
%             end
%         end
%         return;
%     end
% end

%% Open Old Stimfile
if ~stimulus.replay
    stimulus.counter = 1;
    
    if ~isempty(mglGetSID) && isdir(sprintf('~/data/afcom/%s',mglGetSID))
        % Directory exists, check for a stimefile
        files = dir(sprintf('~/data/afcom/%s/1*mat',mglGetSID));
        
        if length(files) >= 1
            fname = files(end).name;
            
            s = load(sprintf('~/data/afcom/%s/%s',mglGetSID,fname));
            % copy staircases and run numbers
            stimulus.counter = s.stimulus.counter + 1;
            stimulus.live.attend = mod(s.stimulus.live.attend+1,3);
            if s.stimulus.attend ~= stimulus.attend
                error('(afcom) Cannot continue: stimfile parameters were generated with a different attention mode than you requested. You need to save the existing stimfiles into a backup folder');
            end
            clear s;
            disp(sprintf('(afcom) Data file: %s loaded.',fname));
        else
            warning('(afcom) Unable to load previous data files. If this is *not* the first run there is something wrong.');
        end
    end
end

%% Display run info
if ~stimulus.replay
    disp('*************************');
    disp(sprintf('(afcom) This is scan #%i',stimulus.counter));
    disp('*************************');
end

%% Setup Screen
if stimulus.replay
    myscreen = initScreen('replayScreen');
else
    myscreen = initScreen('VPixx');
end
% set background to grey
myscreen.background = 0.5;


%% Plot and return
if stimulus.plots==2
    dispInfo(stimulus);
    return
end

%% Initialize Stimulus

if ~stimulus.replay
    myscreen.stimulusNames{1} = 'stimulus';
    
    stimulus.responseKeys = [1 2]; % left right
else
    localInitStimulus();
end

%% Colors
stimulus.colors.white = [1 1 1]; stimulus.colors.red = [1 0 0];
stimulus.colors.green = [0 1 0]; stimulus.colors.black = [0 0 0];

% setup color wheel
stimulus.colorwheel.stops = [255 0 0 ; 255 255 0 ; 0 255 0 ; 0 255 255 ; 0 0 255 ; 255 0 255];
stimulus.colorwheel.thetas = 0:(2*pi/6):(2*pi-.05);

%% Sizes
stimulus.fixWidth = 0.5;
stimulus.probeWidth = 0.5;

%% Setup Probe Task

task{1}{1} = struct;
% task waits for fixation on first segment
stimulus.seg.cue = 1;
stimulus.seg.stim = 2;
stimulus.seg.delay = 3;
stimulus.seg.resp = 4;
stimulus.seg.iti = 5;

task{1}{1}.segmin = [1 2 2 3 2];
task{1}{1}.segmax = [1 2 2 3 8];

task{1}{1}.waitForBacktick = 1;

task{1}{1}.getResponse = zeros(1,length(task{1}{1}.segmin));
task{1}{1}.getResponse(stimulus.seg.resp) = 1;

task{1}{1}.numTrials = Inf;

task{1}{1}.random = 1;

task{1}{1}.parameter.attend = [0 1 2]; % 0 = both, 1 = left, 2 = right
task{1}{1}.parameter.color = [0 1 2]; % 0 = both, 1 = red, 2 = green
stimulus.targetColors = {'black','red','green'};
task{1}{1}.parameter.on = 6; % number of elements in the target group

if ~stimulus.replay && stimulus.scan
    task{1}{1}.synchToVol = zeros(1,length(task{1}{1}.segmin));
    task{1}{1}.synchToVol(end) = 1;
end

task{1}{1}.randVars.calculated.dead = nan;
task{1}{1}.randVars.calculated.targetAngle = nan;
task{1}{1}.randVars.calculated.respAngle = nan;
task{1}{1}.randVars.calculated.respDistance = nan;
task{1}{1}.randVars.calculated.target = nan; % which of the 12 is the target

%% Setup grid

stimulus.gridX = -10:4:10;
stimulus.gridY = -9:6:9;

stimulus.grid = zeros(length(stimulus.gridY),length(stimulus.gridX),50); % 6 minutes should end up being ~ 24 trials, 50 is a safe upper bound


stimulus.gridCount = 1;
stimulus.live.grid = stimulus.grid(:,:,1);
stimulus.live.right = stimulus.live.grid;
stimulus.live.right(:,stimulus.gridX>0) = 1;
stimulus.live.right = logical(stimulus.live.right);
stimulus.live.left = ~stimulus.live.right;

stimulus.live.leftIdxs =  find(stimulus.live.left(:));
stimulus.live.rightIdxs = find(stimulus.live.right(:));
stimulus.live.allIdxs = 1:(size(stimulus.grid,1)*size(stimulus.grid,2));

stimulus.live.trackingAngle = 0;

%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,@endTrialCallback,[]);
end

%% EYE CALIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~stimulus.replay
    myscreen = eyeCalibDisp(myscreen);
    
    % let the user know
    disp(sprintf('(afcom) Starting run number: %i.',stimulus.counter));
end

%% Main Task Loop

% setGammaTable(1);
mglClearScreen(0.5);
mglFlush;
mglClearScreen(0.5);
mglFlush;

phaseNum = 1;
% Again, only one phase.
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    % update the task
    [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);
    % flip screen
    myscreen = tickScreen(myscreen,task);
end

% task ended
mglClearScreen(0.5);
mglTextSet([],32,stimulus.colors.white);
% get count
mglTextDraw('Please wait',[0 0]);
mglFlush
myscreen.flushMode = 1;

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

if ~stimulus.replay && stimulus.plots
    disp('(afcom) Displaying plots');
    dispInfo(stimulus);
end

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%


function [task, myscreen] = startTrialCallback(task,myscreen)
global stimulus

sideOpts = {'all','left','right'};
colorOpts = {'','red','green'};

% pick the locations for the targets
targetLocs = task.thistrial.on;

sides = {'left','right'};

for side = 1:2
    pos = zeros(size(stimulus.live.grid(stimulus.live.(sides{side}))));
    
    % turn on some portion of these
    perm = randperm(length(pos));
    idxs = perm(1:targetLocs);
    
    pos(idxs(1:(targetLocs/2))) = 1;
    pos(idxs((targetLocs/2+1):end)) = 2;

    stimulus.live.grid(stimulus.live.(sides{side})) = pos;
end

stimulus.live.targets = stimulus.live.grid(:);
if task.thistrial.attend>0
    if task.thistrial.attend==1
        stimulus.live.sideIdxs = stimulus.live.leftIdxs;
    else
        stimulus.live.sideIdxs = stimulus.live.rightIdxs;
    end
else
   stimulus.live.sideIdxs = stimulus.live.allIdxs; 
end
if task.thistrial.color>0
    stimulus.live.targetIdxs = find(stimulus.live.targets==task.thistrial.color);
else
    stimulus.live.targetIdxs = find(stimulus.live.targets>0);
end

% pick the target
task.thistrial.target = randsample(intersect(stimulus.live.sideIdxs,stimulus.live.targetIdxs),1);
task.thistrial.targetAngle = 0;
task.thistrial.respAngle = rand*2*pi;

stimulus.live.target = zeros(size(stimulus.live.grid));
stimulus.live.target(task.thistrial.target) = 1;

disp(sprintf('(afcom) Starting trial %i. Attending %s %s',task.trialnum,sideOpts{task.thistrial.attend+1},colorOpts{task.thistrial.color+1}));

task.thistrial.dead = 0;
stimulus.live.eyeCount=0;

function [task, myscreen] = endTrialCallback(task,myscreen)

global stimulus

stimulus.grid(:,:,stimulus.gridCount) = stimulus.live.grid;
stimulus.gridCount = stimulus.gridCount + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task,myscreen)
%%
global stimulus

disp(task.thistrial.thisseg)

if task.thistrial.thisseg == stimulus.seg.cue
        % fixation
        mglClearScreen(0.5);
        drawFix(task);
        mglFlush;
        mglClearScreen(0.5);
        drawFix(task);
        mglFlush;
elseif task.thistrial.thisseg == stimulus.seg.stim
        mglClearScreen(0.5);
        drawStim();
        drawFix(task,stimulus.colors.black);
        mglFlush;
        mglClearScreen(0.5);
        drawStim();
        drawFix(task,stimulus.colors.black);
        mglFlush;
elseif task.thistrial.thisseg == stimulus.seg.delay
        mglClearScreen(0.5);
        drawFix(task,stimulus.colors.black);
        mglFlush;
        mglClearScreen(0.5);
        drawFix(task,stimulus.colors.black);
        mglFlush;
        
elseif task.thistrial.thisseg == stimulus.seg.resp
        % drawing happens in updateScreen
        
elseif task.thistrial.thisseg == stimulus.seg.iti
        mglClearScreen(0.5);
        drawFix(task,stimulus.colors.black);
        mglFlush;
        mglClearScreen(0.5);
        drawFix(task,stimulus.colors.black);
        mglFlush;
end

function drawStim()

global stimulus

colors = {'red','green'};

for yi = 1:size(stimulus.live.grid,1)
    for xi = 1:size(stimulus.live.grid,2)
        if stimulus.live.grid(yi,xi)>0
            mglFillRect(stimulus.gridX(xi),stimulus.gridY(yi),repmat(stimulus.probeWidth,1,2),stimulus.colors.(colors{stimulus.live.grid(yi,xi)}));
        end
    end
end

function drawFix(task,color)

global stimulus;

if task.thistrial.dead
    mglGluDisk(0,0,[1 1],stimulus.colors.red,60);
elseif task.thistrial.thisseg==stimulus.seg.cue   
    if task.thistrial.attend==0
        % both
        leftColor = stimulus.colors.(stimulus.targetColors{task.thistrial.color+1});
        rightColor = leftColor;
    elseif task.thistrial.attend==1
        leftColor = stimulus.colors.(stimulus.targetColors{task.thistrial.color+1});
        rightColor = [];
    else
        leftColor = [];
        rightColor = stimulus.colors.(stimulus.targetColors{task.thistrial.color+1});
    end
    
    mglLines2(0,-stimulus.fixWidth/2,0,stimulus.fixWidth/2,1,stimulus.colors.black);
    % left
    if ~isempty(leftColor)
        mglLines2(-stimulus.fixWidth/2,0,0,0,1,leftColor);
    end
    % right
    if ~isempty(rightColor)
        mglLines2(0,0,stimulus.fixWidth/2,0,1,rightColor);
    end
else
    mglFixationCross(stimulus.fixWidth,1,color);
end

function col = angle2rgb(ang)

global stimulus

col = interp1(stimulus.colorwheel.thetas',stimulus.colorwheel.stops,ang);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus

if (task.thistrial.thisseg==stimulus.seg.resp)
    for yi = 1:size(stimulus.live.grid,1)
        for xi = 1:size(stimulus.live.grid,2)
            if stimulus.live.target(yi,xi)>0
                mglFillRect(stimulus.gridX(xi),stimulus.gridY(yi),repmat(stimulus.probeWidth,1,2),stimulus.colors.white);
            end
        end
    end

    % Draw the color picker
    for theta = 0:10:359
        mglGluPartialDisk(0,0,1,1.25,theta-5,10,angle2rgb(pi/180*theta)/255);
    end
    
    cColor = angle2rgb(task.thistrial.respAngle)/255;

    % Draw the current angle
    mglGluPartialDisk(0,0,1,1.25,180/pi*(task.thistrial.respAngle-pi/16),22.5,cColor);

    % Draw the current color
    mglFillOval(0,0,repmat(stimulus.fixWidth,1,2),cColor);

    % Fixation on top
    drawFix(task,stimulus.colors.black);

end

if (task.thistrial.thisseg==stimulus.seg.resp) && stimulus.powerwheel
    mInfo = mglGetMouse(myscreen.screenNumber);
    curPos = -mInfo.x/90;
    stimulus.live.angle = stimulus.live.angle + curPos-stimulus.live.trackingAngle;
    stimulus.live.trackingAngle = curPos;
elseif task.thistrial.thisseg==stimulus.seg.resp
    mInfo = mglGetMouse(myscreen.screenNumber);
    degx = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
    degy = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;
    task.thistrial.respAngle = -atan2(degy,degx)+pi/2;
end



% if ~stimulus.replay
%     drawFix(myscreen);
% 
%     % check eye pos
%     if (~stimulus.noeye) && (stimulus.eyewindow>0)
% 
% 
%         % mouse version for testing with no eyetracker
%         if stimulus.mousedebug
%             mInfo = mglGetMouse(myscreen.screenNumber);
%             degx = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
%             degy = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;
% 
%             pos = [degx, degy];
%         else
%             [pos,~] = mglEyelinkGetCurrentEyePos;
%         end
% 
%         if stimulus.debug > 0
%             if stimulus.debug >= 10
%                 disp(sprintf('Mouse position: %1.1f %1.1f',pos(1),pos(2)));
%                 stimulus.debug = 1;
%             else
%                 stimulus.debug = stimulus.debug+1;
%             end
%         end
% 
% 
%         % compute distance
%         dist = hypot(pos(1),pos(2));
%     end
% 
%     % Eye movement detection code
%     if (~stimulus.noeye) && (stimulus.eyewindow>0) && ~task.thistrial.dead
%         if ~any(isnan(pos))
% 
%             if dist > stimulus.eyewindow && stimulus.live.eyeCount > 40
%                 disp('Eye movement detected!!!!');
%                 task.thistrial.dead = 1;
%                 return
%             elseif dist > stimulus.eyewindow
%                 stimulus.live.eyeCount = stimulus.live.eyeCount + 1;
%             end
%         end
%     end
% end

function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus

if task.thistrial.dead, return; end

if stimulus.powerwheel
    validResponse = task.thistrial.mouseButton == 1;
else
    validResponse = task.thistrial.whichButton == stimulus.responseKeys(1);
end

if validResponse
    if stimulus.live.gotResponse==0
        disp(sprintf('Received response angle of %i',stimulus.live.angle));
    else
        disp(sprintf('Subject responded multiple times: %i',stimulus.live.gotResponse));
    end
    stimulus.live.gotResponse=stimulus.live.gotResponse+1;
end