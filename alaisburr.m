% alaisburr.m
%
%      usage: myscreen = alaisburr()
%         by: minyoung lee
%       date: 
%    purpose: replication of Alais & Burr, 2004

function myscreen = alaisburr(varargin)
 
clear global stimulus
mglEatKeys('12`');
global stimulus
 
% get arguments
width = 32; visual = 0; auditory = 0; bimodal = 0; disp = 0;
getArgs(varargin,{'width=32','visual=0','auditory=0','bimodal=0','disp=0'},'verbose=1');

if sum([visual,auditory,bimodal]) > 1
    warning('(alaisburr) More than one task type detected.');
    return
elseif sum([visual,auditory,bimodal]) == 0
    warning('(alaisburr) Task type unspecified. Running visual task...')
    return
end
if visual
    stimulus.task = 1;
elseif auditory
    stimulus.task = 2;
else
    stimulus.task = 3;
end
stimulus.disp = disp;
%%%%%%%%%%%%%%%%%%%%%%%%%
stimulus.width = width;
stimulus.stimDur = .015; % 15ms
stimulus.gaussainDur = .015; % 15ms
stimulus.clickDur = 0.0015; % 1.5ms
stimulus.ISI = .500; % 500ms
stimulus.contrast = .1; % 10% contrast
stimulus.interval = [2 4];
% fixation cross
stimulus.fixWidth = 1;
stimulus.fixColor = [1 1 1];

% initalize the screen
myscreen.background = 0;  %black
myscreen = initScreen(myscreen);

%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%% 
task{1}{1}.waitForBacktick = 1;
% trial: Fixation + Stim1 (.015s/.0015s) + ISI (.5s) + Stim2 (.015s/.0015s) + Resp + ITI
task{1}{1}.segmin = [1 stimulus.stimDur stimulus.ISI stimulus.stimDur 1.5 1];
task{1}{1}.segmax = [1 stimulus.stimDur stimulus.ISI stimulus.stimDur 1.5 1];
task{1}{1}.getResponse = [0 0 0 0 1 0];
task{1}{1}.numBlocks = 8;
% parameters & randomization
task{1}{1}.parameter.centerWhich = [1 2]; % centered in which interval
task{1}{1}.random = 1;
task{1}{1}.parameter.posDiff = [-15 -10 -5 -2.5 -1.25 0 1.25 2.5 5 10 15]; 

task{1}{1}.randVars.calculated.resp = nan;
task{1}{1}.randVars.calculated.correct = nan;
task{1}{1}.randVars.calculated.diff = nan;
task{1}{1}.randVars.calculated.rt = nan;
task{1}{1}.randVars.calculated.centerInt = nan;
 
% initialize the task
for phaseNum = 1:length(task{1})
  [task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end
 
% init the stimulus
myscreen = initStimulus('stimulus',myscreen);

% to initialize the stimulus for your experiment.
stimulus = initGaussian(stimulus,myscreen);
stimulus = initClick(stimulus,myscreen);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
stimulus.endflag = 0;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc && ~stimulus.endflag
  % update the task
  [task{1} myscreen phaseNum] = updateTask(task{1},myscreen,phaseNum);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end
 
% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

if stimulus.disp
dispPsychometric(task{1}{1});
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)
global stimulus
if task.thistrial.thisseg == 1
    stimulus.fixColor = [1 1 1];
    task.thistrial.jitter = rand - 0.5; %random jittering between -0.5 and 0.5 deg
    % horizontal position of first, second stim
    if ~task.thistrial.posDiff
        task.thistrial.xpos = [task.thistrial.jitter, task.thistrial.jitter];
    else
        if task.thistrial.centerWhich == 1
            task.thistrial.xpos = [task.thistrial.jitter, task.thistrial.posDiff + task.thistrial.jitter];
            task.thistrial.centerint = 1;
        else
            task.thistrial.xpos = [task.thistrial.jitter - task.thistrial.posDiff, task.thistrial.jitter];
            task.thistrial.centerint = 2;
        end
    end

    task.thistrial.diff = task.thistrial.posDiff;
    
    if stimulus.task ~= 1 %auditory or bimodal condition
        for int = 1:2
            stimulus.sound(int) = createITD(stimulus,task.thistrial.xpos(int));
        end
    end

end
if task.thistrial.thisseg == 6
    if exist('task.thistrial.reactionTime', 'var')
        task.thistrial.rt = task.thistrial.reactionTime;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)
global stimulus
mglClearScreen(0);

if stimulus.task ~= 2 %visual or bimodal condition
    if task.thistrial.thisseg == stimulus.interval(1)
        mglBltTexture(stimulus.tex, [task.thistrial.xpos(1), 0]);
    elseif task.thistrial.thisseg == stimulus.interval(2)
        mglBltTexture(stimulus.tex, [task.thistrial.xpos(2), 0]);
    end
end
if stimulus.task ~= 1 %auditory or bimodal condition
    if task.thistrial.thisseg == stimulus.interval(1)
        mglPlaySound(stimulus.sound(1));
    elseif task.thistrial.thisseg == stimulus.interval(2)
        mglPlaySound(stimulus.sound(2));
    end
end

%draw fixation cross
if task.thistrial.thisseg == 5 || task.thistrial.thisseg == 6
    mglFixationCross(stimulus.fixWidth,1.5,stimulus.fixColor*.5);
else
    mglFixationCross(stimulus.fixWidth,1.5,stimulus.fixColor);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)
global stimulus
% here, we just check whether this is the first time we got a response
if ~task.thistrial.gotResponse
    % which one seemed more to the LEFT
    % centerWhich (1/2) first or second one more eccentric
    if (task.thistrial.posDiff > 0 && task.thistrial.whichButton == 1) || ...
            (task.thistrial.posDiff < 0 && task.thistrial.whichButton == 2)
        % correct
        task.thistrial.correct = 1;
        % feeback
        stimulus.fixColor = [0 1 0];
        disp(sprintf('(alaisburr) Trial %i: %0.4f resp %i correct', ...
            task.trialnum, task.thistrial.posDiff, task.thistrial.whichButton))
    else
        % incorrect
        task.thistrial.correct = 0;
        stimulus.fixColor = [1 0 0];
        disp(sprintf('(alaisburr) Trial %i: %0.4f resp %i incorrect', ...
            task.trialnum, task.thistrial.posDiff, task.thistrial.whichButton))
    end
        
    task.thistrial.resp = task.thistrial.whichButton;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initGaussian(stimulus,myscreen)

% compute the guassian
gauss = mglMakeGaussian(stimulus.width,stimulus.width, stimulus.width/8,stimulus.width/8);

gaussian = zeros(size(gauss,1), size(gauss,2), 4);
for i = 1:3
    gaussian(:,:,i) = 255*ones(size(gauss,1), size(gauss,2));
end
    gaussian(:,:,4) = 255*gauss*stimulus.contrast;
    
%create texture
stimulus.tex = mglCreateTexture(gaussian);

 %stim centers
% [stimulus.x, stimulus.y] = pol2cart(0*pi/180,stimulus.eccentricity);

function stimulus = initClick(stimulus,myscreen)
% sampling frequency (samples per sec)
duration = 0.0015; % 1.5ms
fs = 8192;
t = 0:2*pi/fs:2*pi*duration;
% frequency of signal in hz
hz = 440;
amplitude = 0.5;
stimulus.wav = amplitude * sin(2*pi*hz*t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to update the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = createITD(stimulus,theta)
% theta: angle to sound source (-90 to left and +90 to right)

%%% interaural time difference
% radius of head in meter (assuming spherical)
r = 0.09;
% speed of sound at 23 deg celcius (m/s)
c = 345;
fs = 8192;
% for low frequency sound
td = (2*sind(theta)) * r / c;
td_a = 0:1/fs:abs(td);
clear waveform  s
if td > 0
    waveform(1,:) = [stimulus.wav, zeros(1,length(td_a))];
    waveform(2,:) = [zeros(1,length(td_a)), stimulus.wav];
elseif td < 0
    waveform(2,:) = [stimulus.wav, zeros(1,length(td_a))];
    waveform(1,:) = [zeros(1,length(td_a)), stimulus.wav];
else
    waveform(1,:) = stimulus.wav;
    waveform(2,:) = stimulus.wav;
end

s = mglInstallSound(waveform);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display psychometric functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispPsychometric(task)
posDiff = task.parameter.posDiff;
n = zeros(1,length(posDiff)); k = zeros(1,length(posDiff));
% percent interval 2 (probe "Left")
for i = 1:length(posDiff)
    resp{i} = task.randVars.resp(task.randVars.diff == posDiff(i));
    n(i) = sum(resp{i} == 1 | resp{i}==2);
    k(i) = sum(resp{i} == 2);
end
percent = k./n;

figure;
h = plot(posDiff, percent, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
ylabel('Proportion of trials probe seen "left"');
xlabel('Displacement of probe (degs)');
axis([-20 20 0 1]); box off;

