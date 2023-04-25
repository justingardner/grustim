%%% Drawing some dots at different timing points and locations with sound
%%% cue. Observers are supposed to look at the dots according to the sound
%%% cue. 
%%% 
%%% Jiwon Yeon, April 2023

function validation_drawingDots
mglClose;   % close any MGL windows if opened

% specify the participant
id = input('ID? : ', 's');
saving_name = [id '_' datestr(now,'mmddyy-HHMMSS')];
mglSetSID('test')

% using eyelink
eyetracker = input('Eyelink?(Y/N) : ', 's');

% task repetition
task_repetition = 3;

% initialize myscreen
myscreen.displayName = 'vpixx';
myscreen.saveData = 1;
myscreen.datadir = '/Users/gru/proj/jiwon';
if ~exist(myscreen.datadir), mkdir(myscreen.datadir); end
mglSetParam('abortedStimfilesDir', '/Users/gru/proj/jiwon/aborted',1);
myscreen.keyboard.nums = [37,50];  % Enter, Spacebar
myscreen = initScreen(myscreen);  

% set sitmulus
% target positions
global stimulus
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

% sound_cue
Fs = 10000;
t = linspace(0,.1,Fs);
w = 2*pi*10000;
s = sin(t*w);
stimulus.sound = s;

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

% start the experiment 
mglClearScreen(0);
mglTextSet('Arial',32,1);
mglTextDraw('Press Enter to start',[0,0]);
mglFlush;
disp('Press Enter to start the experiment')
while 1
    k = mglGetKeys;
    if k(myscreen.keyboard.nums(1))==1, break; end
end

% show the target locations
for d = 1:size(positions,1)
    mglGluDisk(positions(d,1), positions(d,2), [1 1], [.5 .5 .5]);
end
mglGluDisk(positions(1,1), positions(1,2), [1 1], [1 0 0]);
mglFlush();
mglWaitSecs(1);

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
mglTextDraw(['Test ended!'], [0, 0])
mglFlush;
mglWaitSecs(3);
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)
global stimulus

% present target
positions = stimulus.positions;
% current target
t = task.thistrial.target_location;

% show target 
for d = 1:size(positions,1)
    mglGluDisk(positions(d,1), positions(d,2), [1 1], [.5 .5 .5]);
end
mglGluDisk(positions(t,1), positions(t,2), [1 1], [1 0 0]);

% play sound
soundsc(stimulus.sound(1:length(stimulus.sound)/10))

% flush everything
myscreen.flushMode = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task, myscreen)


%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = getResponseCallback(task,myscreen)

if task.thistrial.buttonState(2) == 1
    % move to the next segment
    task = jumpSegment(task);
end












