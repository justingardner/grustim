% posdisc.m
%
%      usage: myscreen = posdisc()
%         by: minyoung lee
%       date: 10/14/16
%    purpose: spatial localization task
%
function myscreen = posdisc(stimType,eccNum)

% gabor VS gaussian
if ~exist('stimType','var')
    stimType = 'gaussian';
end
if ~exist('eccNum','var')
    eccNum = 0;
end

clear global stimulus
mglEatKeys('12`');
global stimulus
stimulus.stimType = stimType;
stimulus.width = 10;
stimulus.sf = 1.8;
stimulus.eccNum = eccNum;
stimulus.eccList = [7.5 10 12.5 15];
if stimulus.eccNum ~= 0
    stimulus.eccentricity = stimulus.eccList(stimulus.eccNum);
end

stimulus.stimDur = .015; % 15ms
stimulus.ISI = .1; % 100ms

stimulus.interval = [2 4];
% stimulus.string = {'Low','High'};

% fixation cross
stimulus.fixWidth = 1;
stimulus.fixColor = [1 1 1];

stimulus.n = 0;

% initalize the screen
myscreen.background = 0.5;
myscreen = initScreen(myscreen);

%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%%

% task{1}{1}.waitForBacktick = 0;
task{1}{1}.waitForBacktick = 1;
% trial: Fixation + Stim1 (.015s) + ISI (.1s) + Stim2 (.015s) + Resp + ITI
task{1}{1}.segmin = [1 stimulus.stimDur stimulus.ISI stimulus.stimDur 1.5 1];
task{1}{1}.segmax = [1 stimulus.stimDur stimulus.ISI stimulus.stimDur 1.5 1];
task{1}{1}.getResponse = [0 0 0 0 1 0];

task{1}{1}.numBlocks = 8;

% parameters & randomization
% task{1}{1}.parameter.reliability = [1 2]; % Low High
% task{1}{1}.parameter.whichHemifield = [1 2]; % Left Right
task{1}{1}.parameter.contrast = [0.0625 0.125 0.25 0.5 1];
if stimulus.eccNum == 0
task{1}{1}.parameter.eccentricity = [7.5 10 12.5 15];
end
task{1}{1}.parameter.posDiff = [-5 -2.5 -1.25 -.5 -.25 -.125 0 .125 .25 .5 1.25 2.5 5];

task{1}{1}.random = 1;

task{1}{1}.randVars.calculated.resp = nan;
task{1}{1}.randVars.calculated.correct = nan;

task{1}{1}.randVars.calculated.rel = nan;
task{1}{1}.randVars.calculated.diff = nan;
% task{1}{1}.randVars.calculated.hemi = nan;
task{1}{1}.randVars.calculated.rt = nan;
task{1}{1}.randVars.calculated.con = nan;
task{1}{1}.randVars.calculated.ecc = nan;

% initialize the task
for phaseNum = 1:length(task{1})
  [task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end

% init the stimulus
myscreen = initStimulus('stimulus',myscreen);

% to initialize the stimulus for your experiment.
% stimulus = initGaussian(stimulus,myscreen);

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

dispPsychometric(task{1}{1}, stimulus);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)
global stimulus
if task.thistrial.thisseg == 1
    stimulus.n = stimulus.n+1;

        if stimulus.n>1 && isnan(task.randVars.correct(stimulus.n-1))
        disp(sprintf('(posdisc) Trial %i: contrast %0.4f ecc %0.1f diff %0.4f No resp', ...
        stimulus.n-1, stimulus.contrast, stimulus.eccentricity, task.randVars.diff(stimulus.n-1)))
        end
    
    stimulus.fixColor = [1 1 1];
    task.thistrial.xpos = [0 0];
    
    stimulus.contrast = task.thistrial.contrast;
    if stimulus.eccNum == 0
        stimulus.eccentricity = task.thistrial.eccentricity;
        task.thistrial.ecc = task.thistrial.eccentricity;
    else
        task.thistrial.ecc = stimulus.eccentricity;
    end
    task.thistrial.diff = task.thistrial.posDiff;
    task.thistrial.con = task.thistrial.contrast;

%     task.thistrial.hemi = task.thistrial.whichHemifield;
    stimulus = initGaussian(stimulus,myscreen);
    
    if task.thistrial.posDiff ~= 0 % <0: 1st left, 2nd right  // >0: 1st right, 2nd left
        task.thistrial.xpos = [stimulus.x + task.thistrial.posDiff/2, ...
            stimulus.x - task.thistrial.posDiff/2];
%         task.thistrial.xpos = [stimulus.meanXpos(task.thistrial.whichHemifield) + task.thistrial.posDiff/2, ...
%             stimulus.meanXpos(task.thistrial.whichHemifield) - task.thistrial.posDiff/2];
    else % diff = 0
        task.thistrial.xpos = [stimulus.x, ...
            stimulus.x];
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
mglClearScreen(0.5);

if task.thistrial.thisseg == stimulus.interval(1)
    mglBltTexture(stimulus.tex, [task.thistrial.xpos(1), stimulus.y]);
elseif task.thistrial.thisseg == stimulus.interval(2)
    mglBltTexture(stimulus.tex, [task.thistrial.xpos(2), stimulus.y]);
end

%draw fixation cross
if task.thistrial.thisseg == 5 || task.thistrial.thisseg == 6
    mglFixationCross(stimulus.fixWidth,1.5,stimulus.fixColor*.75);
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
    % posDiff<0: 1st left, 2nd right  // >0: 1st right, 2nd left
    if (task.thistrial.posDiff < 0 && task.thistrial.whichButton == 2) || (task.thistrial.posDiff > 0 && task.thistrial.whichButton == 1) 
        % correct
        task.thistrial.correct = 1;
        % feeback
        stimulus.fixColor = [0 1 0];
        disp(sprintf('(posdisc) Trial %i: contrast %0.4f ecc %0.1f diff %0.4f resp %i correct', ...
            task.trialnum, stimulus.contrast, stimulus.eccentricity, task.thistrial.posDiff, task.thistrial.whichButton))
    else
        % incorrect
        task.thistrial.correct = 0;
        stimulus.fixColor = [1 0 0];
        disp(sprintf('(posdisc) Trial %i: contrast %0.4f ecc %0.1f diff %0.4f resp %i incorrect', ...
            task.trialnum, stimulus.contrast, stimulus.eccentricity, task.thistrial.posDiff, task.thistrial.whichButton))
    end 
    task.thistrial.resp = task.thistrial.whichButton;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initGaussian(stimulus,myscreen)

if strcmpi(stimulus.stimType,'gaussian')
% compute the guassian
gauss = mglMakeGaussian(stimulus.width,stimulus.width,stimulus.width/8, stimulus.width/8);
% gaussHigh = mglMakeGaussian(stimulus.width,stimulus.width,stimulus.width/8, stimulus.width/8);
gaussian = zeros(size(gauss,1), size(gauss,2), 4);
for i = 1:3
    gaussian(:,:,i) = 255*ones(size(gauss,1), size(gauss,2));
end
    gaussian(:,:,4) = 255*gauss*stimulus.contrast;
    
%create texture
stimulus.tex = mglCreateTexture(gaussian);

elseif strcmpi(stimulus.stimType, 'gabor')
    grating = mglMakeGrating(stimulus.width,stimulus.width,stimulus.sf, 0,0);
    gaussian = mglMakeGaussian(stimulus.width,stimulus.width,stimulus.width/8, stimulus.width/8);
    gabor = (255*(stimulus.contrast*grating.*gaussian+1)/2);
    stimulus.tex = mglCreateTexture(gabor);
end

 %stim centers
[stimulus.x, stimulus.y] = pol2cart(0*pi/180,stimulus.eccentricity);
% [stimulus.x, stimulus.y] = pol2cart(30*pi/180,stimulus.eccentricity);
% [stimulus.x, stimulus.y] = pol2cart(330*pi/180,stimulus.eccentricity);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display psychometric functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispPsychometric(task,stimulus)
posDiff = task.parameter.posDiff;
contrast = task.parameter.contrast;
if stimulus.eccNum == 0
eccentricity = task.parameter.eccentricity;
%percent Interval 1
for c = 1:length(contrast)
    for ecc = 1:length(eccentricity)
    for i = 1:length(posDiff)
        resp{c}{ecc}{i} = task.randVars.resp(task.randVars.con == contrast(c) & ...
            task.randVars.ecc == eccentricity(ecc) & task.randVars.diff == posDiff(i));
        n{c}{ecc}(i) = sum(resp{c}{ecc}{i} == 1 | resp{c}{ecc}{i} == 2);
        k{c}{ecc}(i) = sum(resp{c}{ecc}{i} == 1);
    end
    percent{c}{ecc} = k{c}{ecc}./n{c}{ecc};
    end
end

figure;
for c = 1:length(contrast)
    subplot(2,3,c)
    for ecc = 1:length(eccentricity)
        plot(posDiff, percent{c}{ecc}*100, getcolor(ecc,getsymbol(ecc)), 'MarkerSize', 7, 'MarkerFaceColor', getcolor(ecc));
        hold on;
    end
    title(sprintf('Contrast=%0.4f', contrast(c)));
    legend('ecc=7.5','ecc=10','ecc=12.5','ecc=15','Location','SouthEast');
    box off;
    ylabel('Percent choices Interval 1 (%)');
    xlabel('Postition difference of targets between interval 1 and 2 (deg)');
    yaxis(0,100);
end
else
    ecc = 1;
    %percent Interval 1
for c = 1:length(contrast)
    
    for i = 1:length(posDiff)
        resp{c}{ecc}{i} = task.randVars.resp(task.randVars.con == contrast(c) & ...
           task.randVars.diff == posDiff(i));
        n{c}{ecc}(i) = sum(resp{c}{ecc}{i} == 1 | resp{c}{ecc}{i} == 2);
        k{c}{ecc}(i) = sum(resp{c}{ecc}{i} == 1);
    end
    percent{c}{ecc} = k{c}{ecc}./n{c}{ecc};

end
figure;
for c = 1:length(contrast)
    subplot(2,3,c)
        plot(posDiff, percent{c}{ecc}*100, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
        hold on;

    title(sprintf('Contrast=%0.4f', contrast(c)));
    legend(sprintf('ecc=%0.1f',stimulus.eccentricity),'Location','SouthEast');
    box off;
    ylabel('Percent choices Interval 1 (%)');
    xlabel('Postition difference of targets between interval 1 and 2 (deg)');
    yaxis(0,100);
end

end
