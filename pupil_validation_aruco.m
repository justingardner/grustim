%%% Drawing some dots at different timing points and locations with sound
%%% cue. Observers are supposed to look at the dots according to the sound
%%% cue.
%%%
%%% Jiwon Yeon, April 2023

function pupil_validation_aruco

global stimulus device
mglClose;   % close any MGL windows if opened

% specify the participant
id = input('ID? : ', 's');
mglSetSID('test');

% using eyelink
eyetracker = input('Eyelink?(Y/N) : ', 's');

% connect pupil labs
ip = '10.35.127.131';
device = pyrunfile("connect_pupilLabs.py","device", ip = ip);
if isempty(device.phone_id)
    error("Pupil lab device is not connected.")
end

% task repetition
task_repetition = 3;

% initialize myscreen
myscreen.displayName = 'dell_wuTsai';
myscreen.saveData = 1;
myscreen.datadir = ['Users/jyeon/Documents/GitHub/pupillabs_validation_exp/data'];
if ~exist(myscreen.datadir), mkdir(myscreen.datadir); end
mglSetParam('abortedStimfilesDir', [myscreen.datadir '/aborted'],1);
myscreen.keyboard.nums = [37,50];  % Enter, Spacebar
myscreen = initScreen(myscreen);

% set sitmulus
% target positions
% width = mglGetParam('screenWidth');
% height = mglGetParam('screenHeight');
% positions = ...
%     [[width/2, height/2]; [40, 40]; [width-40, 40]; ...
%     [40, height/2];  [width-40, height/2]; ...
%     [40, height-40]; [width-40, height-40]; ...
%     [width/2-500, height/2]; [width/2-350, height/2-350]; [width/2, height/2-500]; [width/2+350, height/2-350]; ...
%     [width/2+500, height/2]; [width/2+350, height/2+350]; [width/2, height/2+500]; [width/2-350, height/2+350]; ...
%     [width/2-200, height/2]; [width/2-140, height/2-140]; [width/2, height/2-200]; [width/2+140, height/2-140]; ...
%     [width/2+200, height/2]; [width/2+140, height/2+140]; [width/2, height/2+200]; [width/2-140, height/2+140]];

% determine target locations (angle)
disp_size = floor([mglGetParam('deviceWidth'), mglGetParam('deviceHeight')]);
disp_dist = mglGetParam('devicePhysicalDistance');
screen_limit_ang = floor(rad2deg(atan(floor(disp_size/2)/disp_dist)));  % screen limit in angle
screen_limit_ang = screen_limit_ang-1;
xlim = screen_limit_ang(1);
ylim = screen_limit_ang(2);

positions = [[-xlim, ylim]; [0, ylim]; [xlim, ylim]; [xlim, 0]; ...
    [xlim, -ylim]; [0, -ylim]; [-xlim, -ylim]; [-xlim, 0]];
positions = [positions; floor(positions./2)];
positions = [[0, 0]; positions];
stimulus.positions = positions;
task{1}.parameter.positions = positions;

% % sound_cue
% Fs = 10000;
% t = linspace(0,.1,Fs);
% w = 2*pi*10000;
% s = sin(t*w);
% stimulus.sound = s;

% call image marker
imagepath = [pwd '/images/aruco/'];
aruco5 = imread(fullfile(imagepath, 'aruco5.png'));
aruco5 = squeeze(aruco5(:,:,1));
marker_startend = mglCreateTexture(aruco5);
stimulus.marker_startend = marker_startend;

aruco61 = imread(fullfile(imagepath, 'aruco61.png'));
aruco61 = squeeze(aruco61(:,:,1));
marker_target = mglCreateTexture(aruco61);
stimulus.marker_target = marker_target;

% other setups
task{1}.segmin = [nan inf];     % first interval: ITI, second: response
task{1}.segmax = [nan inf];
task{1}.segdur{1} = .8:.2:2.6;
task{1}.getResponse = [0 1];
task{1}.randVars.block.target_location = 1:length(stimulus.positions);
task{1}.randVars.calculated.timing_flip = [nan nan];
task{1}.randVars.calculated.rt = nan;
task{1}.numTrials = length(stimulus.positions)*task_repetition;
task{1}.parameter.marker_size = 1.5;        % visual angle?

% initialize the task
[task{1} myscreen] = initTask(task{1},myscreen,...
    @startSegmentCallback, @updateScreenCallback, @getResponseCallback);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eye calibration (optional)
if strcmpi(eyetracker, 'y')
    disp(' Calibrating Eye ....')
    myscreen = eyeCalibDisp(myscreen);
end

% wait for pupil lab start recording
device.recording_start()
mglWaitSecs(2);

% start the experiment
mglClearScreen(0.5);
mglTextSet('Arial',32,1);
mglBltTexture(stimulus.marker_startend,[0, 0, task{1}.parameter.marker_size, ...
    task{1}.parameter.marker_size])
txt = ['Targets would appear as a marker that you''re seeing above.'];
mglTextDraw(txt,[0, -2]);
txt = ['Press Spacebar when you fixate your eyes at the target, then it will move on to the next target.'];
mglTextDraw(txt,[0, -3.5]);
txt = ['Press Enter to start the session.'];
mglTextDraw(txt,[0, -7]);
mglFlush;
disp('Press Enter to start the experiment')

while 1
    k = mglGetKeys;
    if k(myscreen.keyboard.nums(1))==1, break; end
end
myscreen.timing_experimentStart = mglGetSecs;
device.send_event("experiment start")

% % show the target locations
% for d = 1:size(positions,1)
%     mglGluDisk(positions(d,1), positions(d,2), [1 1], [.5 .5 .5]);
% end
% mglGluDisk(positions(1,1), positions(1,2), [1 1], [.7 .7 .7]);
% mglFlush();
% mglWaitSecs(1);

% task
phaseNum = 1;
while (task{1}.numTrials >= task{1}.trialnum) ...
        && ~myscreen.userHitEsc
    % update the task
    [task myscreen] = updateTask(task,myscreen,phaseNum);
    % flip screen
    myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen.timing_experimentEnd = mglGetSecs;
device.send_event("experiment end")
mglBltTexture(stimulus.marker_startend,[0, 0, ...
    task{1}.parameter.marker_size, ...
    task{1}.parameter.marker_size]);
mglTextDraw(['Test ended!'], [0, -3]);
mglFlush;
mglWaitSecs(3);
device.recording_stop_and_save()
device.close()
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)
global stimulus device

% if it is the first segment, do nothing
if task.thistrial.thisseg == 1
    task.thistrial.timing_flip(1) = mglGetSecs;
    mglClearScreen(.5);
    myscreen.flushMode=1;

    % if it is the second segment, present the target
elseif task.thistrial.thisseg == 2
    t = task.thistrial.target_location;
    position = stimulus.positions(t,:);
    mglBltTexture(stimulus.marker_target, [position(1), position(2), ...
        task.thistrial.marker_size, task.thistrial.marker_size]);
    task.thistrial.timing_flip(2) = mglGetSecs;
    device.send_event('target')
    myscreen.flushMode = 1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task, myscreen)


%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = getResponseCallback(task,myscreen)
global device 
if task.thistrial.buttonState(2) == 1
    % log event
    task.thistrial.rt = mglGetSecs;
    device.send_event('response')

    % move to the next segment
    task = jumpSegment(task);
end
