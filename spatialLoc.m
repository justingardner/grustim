% spatialLoc.m
%
%      usage: myscreen = spatialLoc()
%         by: minyoung lee
%       date: 10/14/16
%    purpose: spatial localization task
%
function myscreen = spatialLoc()

clear global stimulus
global stimulus

stimulus.sdLow = 1.2;
stimulus.sdHigh = .4;
stimulus.width = 6;
stimulus.contrast = 1;
stimulus.eccentricity = 6;
stimulus.meanXpos = [-6 6];

stimulus.stimDur = .015; % 15ms ---> refresh rate???
stimulus.ISI = .1; % 100ms

stimulus.interval = [2 4];
stimulus.string = {'Low','High'};

% fixation cross
stimulus.fixWidth = 1;
stimulus.fixColor = [1 1 1];

% get screen size in visual angle

% initalize the screen
myscreen.background = 0;%0.5;
myscreen = initScreen(myscreen);

%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%%

task{1}{1}.waitForBacktick = 0;
% task{1}{1}.waitForBacktick = 1;
% trial: Fixation + Stim1 (.015s) + ISI (.1s) + Stim2 (.015s) + Resp + ITI
task{1}{1}.segmin = [1 stimulus.stimDur stimulus.ISI stimulus.stimDur 1.5 1];
task{1}{1}.segmax = [1 stimulus.stimDur stimulus.ISI stimulus.stimDur 1.5 1];
task{1}{1}.getResponse = [0 0 0 0 1 0];

% parameters & randomization
task{1}{1}.parameter.reliability = [1 2]; % Low High
task{1}{1}.parameter.whichHemifield = [1 2]; % Left Right (-ecc +ecc)
task{1}{1}.parameter.posDiff = [-1.5 -.75:.25:.75 1.5];
task{1}{1}.random = 1;

task{1}{1}.randVars.calculated.resp = nan;
task{1}{1}.randVars.calculated.correct = nan;

task{1}{1}.randVars.calculated.rel = nan;
task{1}{1}.randVars.calculated.diff = nan;
task{1}{1}.randVars.calculated.hemi = nan;

% initialize the task
for phaseNum = 1:length(task{1})
  [task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end

% init the stimulus
myscreen = initStimulus('stimulus',myscreen);

% to initialize the stimulus for your experiment.
stimulus = initGaussian(stimulus,myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
  % update the task
  [task{1} myscreen phaseNum] = updateTask(task{1},myscreen,phaseNum);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

dispPsychometric(task{1}{1});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)
global stimulus
if task.thistrial.thisseg == 1
    stimulus.fixColor = [1 1 1];
    task.thistrial.xpos = [0 0];
    
    if task.thistrial.posDiff ~= 0 % <0: 1st left, 2nd right  // >0: 1st right, 2nd left
        task.thistrial.xpos = [stimulus.meanXpos(task.thistrial.whichHemifield) + task.thistrial.posDiff/2, ...
            stimulus.meanXpos(task.thistrial.whichHemifield) - task.thistrial.posDiff/2];
    else % diff = 0
        task.thistrial.xpos = [stimulus.meanXpos(task.thistrial.whichHemifield), ...
            stimulus.meanXpos(task.thistrial.whichHemifield)];
    end
    
    task.thistrial.tex = stimulus.tex(task.thistrial.reliability);
    
    task.thistrial.rel = task.thistrial.reliability;
    task.thistrial.diff = task.thistrial.posDiff;
    task.thistrial.hemi = task.thistrial.whichHemifield;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)
global stimulus
mglClearScreen;

if task.thistrial.thisseg == stimulus.interval(1)
    mglBltTexture(task.thistrial.tex, [task.thistrial.xpos(1), 0]);
elseif task.thistrial.thisseg == stimulus.interval(2)
    mglBltTexture(task.thistrial.tex, [task.thistrial.xpos(2), 0]);
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
if ~task.thistrial.gotResponse
    % posDiff<0: 1st left, 2nd right  // >0: 1st right, 2nd left
    if (task.thistrial.posDiff < 0 && task.thistrial.whichButton == 2) || (task.thistrial.posDiff > 0 && task.thistrial.whichButton == 1) 
        % correct
        task.thistrial.correct = 1;
        % feeback
        stimulus.fixColor = [0 1 0];
        disp(sprintf('(spatialdiscr) Trial %i: %s Reliability %0.4f resp %i correct', ...
            task.trialnum, stimulus.string{task.thistrial.reliability}, task.thistrial.posDiff, task.thistrial.whichButton))
    else
        % incorrect
        task.thistrial.correct = 0;
        stimulus.fixColor = [1 0 0];
        disp(sprintf('(spatialdiscr) Trial %i: %s Reliability %0.4f resp %i incorrect', ...
            task.trialnum, stimulus.string{task.thistrial.reliability}, task.thistrial.posDiff, task.thistrial.whichButton))
    end
    task.thistrial.resp = task.thistrial.whichButton;
else
    disp(sprintf('(spatialdiscr) Trial %i: %s Reliability %0.4f no resp', ...
            task.trialnum, stimulus.string{task.thistrial.reliability}, task.thistrial.posDiff))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initGaussian(stimulus,myscreen)

% compute the guassian

gaussLow = mglMakeGaussian(stimulus.width, stimulus.width, stimulus.sdLow, stimulus.sdLow);
gaussHigh = mglMakeGaussian(stimulus.width, stimulus.width, stimulus.sdHigh, stimulus.sdHigh);

for i = 1:4
    gaussianLow(:,:,i) = 255*gaussLow*stimulus.contrast;
    gaussianHigh(:,:,i) = 255*gaussHigh*stimulus.contrast;
end

%create texture
stimulus.tex(1) = mglCreateTexture(gaussianLow);
stimulus.tex(2) = mglCreateTexture(gaussianHigh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display psychometric functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispPsychometric(task)
posDiff = [-1.5 -.75:.25:.75 1.5];
%percent Interval 1
%high rel
highRel = zeros(1,length(posDiff));
for i = 1:length(posDiff)
    eval(sprintf('hResp.d%i = task.randVars.resp(task.randVars.rel == 2 & task.randVars.diff == posDiff(i));',i));
    eval(sprintf('highRel(i) = sum(hResp.d%i == 1) / sum(hResp.d%i == 1 | hResp.d%i == 2)*100;', i,i,i));
    eval(sprintf('highN(i) = sum(hResp.d%i == 1 | hResp.d%i == 2);', i,i));
end

%low rel
lowRel = zeros(1,length(posDiff));
for i = 1:length(posDiff)
    eval(sprintf('lResp.d%i = task.randVars.resp(task.randVars.rel == 1 & task.randVars.diff == posDiff(i));',i));
    eval(sprintf('lowRel(i) = sum(lResp.d%i == 1) / sum(lResp.d%i == 1 | lResp.d%i == 2)*100;', i,i,i));
    eval(sprintf('lowN(i) = sum(lResp.d%i == 1 | lResp.d%i == 2);', i,i));
end

figure;
subplot(1,2,1)
h1 = plot(posDiff, highRel, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
title(sprintf('High reliabilty spatial psychometric function (N=%i)', sum(highN)));
ylabel('Percent choices Interval 1 (%)');
xlabel('Postition difference of targets between interval 1 and 2 (deg)');
axis([-2 2 0 100]); box off;
xlabh = get(gca,'xLabel');
set(xlabh,'Position', get(xlabh, 'Position') + [2,0,0]);

subplot(1,2,2)
h2 = plot(posDiff, lowRel, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
title(sprintf('Low reliabilty spatial psychometric function (N=%i)', sum(lowN)));
% ylabel('Percent choices Interval 1 (%)');
% xlabel('Postition difference of targets between interval 1 and 2 (deg)');
axis([-2 2 0 100]); box off;



