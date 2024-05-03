%%% Drawing some dots at different timing points and locations with sound
%%% cue. Observers are supposed to look at the dots according to the sound
%%% cue.
%%%
%%% updated the calibration window
%%%
%%% Jiwon Yeon, Aug 2023

function validation_drawDots
mglClose;   % close any MGL windows if opened

% % specify the participant
% id = input('ID? : ', 's');
% saving_name = [id '_' datestr(now,'mmddyy-HHMMSS')];
% mglSetSID('test')
% 
% % using eyelink
% eyetracker = input('Eyelink?(Y/N) : ', 's');
% 
% % connect pupil labs
% global device
% ip = '10.35.121.121';
% device = pyrunfile("connect_pupilLabs.py","device", ip = ip);
% if isempty(device.phone_id)
%     error("Pupil lab device is not connected.")
% end

% task repetition
task_repetition = 10;

% initialize myscreen
myscreen.displayName = 'vpixx';
myscreen.saveData = 1;
myscreen.datadir = [pwd '/data'];
if ~exist(myscreen.datadir), mkdir(myscreen.datadir); end
mglSetParam('abortedStimfilesDir', [pwd '/data/aborted'],1);
myscreen.keyboard.nums = [37,50];  % Enter, Spacebar
myscreen = initScreen(myscreen);

% set sitmulus
% target positions
global stimulus

% determine target locations (angle)
stimulus.deviceRect = floor(mglGetParam('deviceRect'));  % screen limit in angle
xlim = stimulus.deviceRect(3)-15;
ylim = stimulus.deviceRect(4)-10;

positions = [[-xlim, ylim]; [0, ylim]; [xlim, ylim]; [xlim, 0]; ...
    [xlim, -ylim]; [0, -ylim]; [-xlim, -ylim]; [-xlim, 0]];
positions = [positions; floor(positions./2)];
positions = [[0, 0]; positions];
stimulus.positions = positions;
task{1}.parameter.positions = positions;

% sound_cue
Fs = 10000;
t = linspace(0,.1,Fs);
w = 2*pi*10000;
s = sin(t*w);
stimulus.sound = s;
myscreen = initStimulus('stimulus', myscreen);

% other setups
task{1}.seglen = inf;
task{1}.getResponse = 1;
task{1}.randVars.block.target_location = 1:length(stimulus.positions);
task{1}.randVars.calculated.rt = nan;
task{1}.numTrials = length(stimulus.positions)*task_repetition;

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
mglClearScreen(0);
mglTextDraw('Press Enter to start',[0,0]);
mglFlush;
disp('Press Enter to start the experiment')
while 1
    k = mglGetKeys;
    if k(myscreen.keyboard.nums(1))==1, break; end
end
device.send_event("experiment start")

% show the target locations
for d = 1:size(positions,1)
    mglGluDisk(positions(d,1), positions(d,2), [.5 .5], [.5 .5 .5]);
end
mglGluDisk(positions(1,1), positions(1,2), [.5 .5], [1, 1, 1]);
mglFlush();
mglWaitSecs(2);

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
device.send_event("experiment end")
mglTextDraw(['Test ended!'], [0, 0]);
mglFlush;
mglWaitSecs(3);
device.recording_stop_and_save()
device.close()
myscreen = endTask(myscreen,task);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)
global stimulus device

% present target
positions = stimulus.positions;

% current target
t = task.thistrial.target_location;

% show target
for d = 1:size(positions,1)
    mglGluDisk(positions(d,1), positions(d,2), [.5 .5], [.5 .5 .5]);
end
mglGluDisk(positions(t,1), positions(t,2), [.5 .5], [1 0 0]);

% play sound
soundsc(stimulus.sound(1:length(stimulus.sound)/10))

% log event
device.send_event('target')

% flush everything
myscreen.flushMode = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task, myscreen)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = getResponseCallback(task,myscreen)
global device

if task.thistrial.buttonState(2) == 1
    % log event
    device.send_event('response')

    % move to the next segment
    task = jumpSegment(task);
end
end
