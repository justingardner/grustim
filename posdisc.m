% posdisc.m
%
%      usage: myscreen = posdisc()
%         by: minyoung lee
%       date: 10/14/16
%    purpose: spatial localization task
%       e.g.: myscreen = posdisc()

%             contrast staircase at a fixed eccentricity:
%             myscreen = posdisc('gabor','staircase=1')

%      flags: 
%      usage: posdisc('gabor=1', 'stairCon=1')
%             posdisc('gabor=1', 'fixCon=1', 'contrast=0.125', 'ecc=10')
%             posdisc('gabor=1', 'fixEcc=1', 'contrast=0.50', 'ecc=20')
%             posdisc('gabor=1', 'constant=1', 'ecc=20')
%
function myscreen = posdisc(varargin)

clear global stimulus
mglEatKeys('12`');
global stimulus

% get arguments
% stim type
gabor=0; gaussian=0;
% task type
stairCon=0; fixCon=0; fixEcc=0; constant=0;
ecc=[]; contrast=[]; plots=0; laststim=0;
stairDist=1;
getArgs(varargin,{'gabor=0','gaussian=0','stairCon=0','constant=0','fixCon=0','fixEcc=0','plots=0','ecc=[]', 'contrast=[]', 'laststim=0'},'verbose=1');
stimType = 'gabor';
if gabor
    stimType='gabor';
elseif gaussian
    stimType='gaussian';
elseif ~gabor && ~gaussian
    warning('(posdisc) Target stimulus undefined. Setting to default: gabor.');
end
if sum([stairCon, fixCon, fixEcc, constant]) > 1
    warning('(posdisc) More than one task type detected.');
    return
elseif sum([stairCon, fixCon, fixEcc, constant]) == 0
    warning('(posdisc) Must specify task type.')
    return
end
    
if stairCon 
    if ~isempty(contrast)
     warning('(posdisc) Staircasing contrast. Overriding contrast argument...');
     contrast=[];
    end
    % Not staircasing position difference
    stairDist = 0;
    ecc=[];
    taskType=1;
    disp('(posdisc) stairCon=1');
else
    if ~isempty(ecc) && ~any(ecc == [10, 15, 17.5, 20])
        warning('(posdisc) Eccentricity must be one from [10 15 17.5 20]');
        return
    elseif isempty(ecc)
        warning('(posdisc) Must specify ecc');
    end
    
    if fixCon 
        if ~isempty(contrast) && ~any(contrast == [0.125 0.25 0.5]) 
            warning('(posdisc) Contrast value must be one from [0.125 0.25 0.5]');
            return
        elseif isempty(contrast)
            warning('(posdisc) Must specify contrast');
            return
        end
        taskType=2;
        disp('(posdisc) fixCon=1');
    elseif fixEcc 
        if ~isempty(contrast) && ~any(contrast == [0.125 0.25 0.5]) 
            warning('(posdisc) Contrast value must be one from [0.125 0.25 0.5]');
            return
        elseif isempty(contrast)
            warning('(posdisc) Must specify contrast');
            return
        end
        taskType=3;
        disp('(posdisc) fixEcc=1');
    elseif constant
        stairDist = 0;
        contrast = [];
        taskType=4;
        disp('(posdisc) constant=1');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%
%%  %%
%%%%%%%%%%%%%%%%%%%%%%%%%
stimulus.stimType = stimType;
stimulus.eccentricity = ecc;
stimulus.stairCon = stairCon;
stimulus.stairDist = stairDist;
stimulus.contrast = contrast;
stimulus.taskType = taskType;

stimulus.width = 10;
stimulus.sf = 1.8;
stimulus.stimDur = .015; % 15ms
stimulus.ISI = .1; % 100ms
stimulus.interval = [2 4];
stimulus.n = 0; %count n trials
% fixation cross
stimulus.fixWidth = 1;
stimulus.fixColor = [1 1 1];
if stimulus.eccentricity == 20
    stimulus.fixOrigin = [-5 0];
else
    stimulus.fixOrigin = [0 0];
end
% [7.5 10 12.5 15]
% stimulus.eccList = [10 12.5 15 17.5 20] + stimulus.fixOrigin(1); 
if ~isempty(stimulus.eccentricity)
    stimulus.eccentricity_org = stimulus.eccentricity;
    stimulus.eccentricity = stimulus.eccentricity + stimulus.fixOrigin(1);
end

% set up staircase
% if stimulus.stairDist
%     stimulus.stair.dist.initialThreshold = 1.25;
%     stimulus.stair.dist.minStepsize = 0.005;
%     stimulus.stair.dist.minThreshold = 0;
%     stimulus.stair.dist.maxThreshold = 5;
% end
if stimulus.stairCon
    
%     stimulus.stair.con.initialThreshold = .25;
%     stimulus.stair.con.minStepsize = 0.005;
%     stimulus.stair.con.minThreshold = 0.01;
%     stimulus.stair.con.maxThreshold = 100;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    stimulus.posDiff = 0.5; % fixing position difference
    stimulus.eccentricity = 15 + stimulus.fixOrigin(1); % fixing eccentricity
end

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

% task{1}{1}.numBlocks = 8;
if stimulus.taskType==4
    task{1}{1}.numBlocks = 4;
else
    task{1}{1}.numTrials = 110;
end
% parameters & randomization
% task{1}{1}.parameter.whichHemifield = [1 2]; % Left Right
task{1}{1}.parameter.whichInterval = [1 2]; % more eccentric (to the right) in which interval (1 or 2)
task{1}{1}.random = 1;

if stimulus.taskType == 4 % constant stimuli
    task{1}{1}.parameter.contrast = [0.125 0.5];
    task{1}{1}.parameter.posDiff = [0 .125 .25 .5 1.25 2.5 5];
end
task{1}{1}.randVars.calculated.resp = nan;
task{1}{1}.randVars.calculated.correct = nan;

% task{1}{1}.randVars.calculated.rel = nan;
task{1}{1}.randVars.calculated.diff = nan;
% task{1}{1}.randVars.calculated.hemi = nan;
task{1}{1}.randVars.calculated.rt = nan;
task{1}{1}.randVars.calculated.con = nan;
task{1}{1}.randVars.calculated.ecc = nan;
task{1}{1}.randVars.calculated.whichint = nan;

% initialize the task
for phaseNum = 1:length(task{1})
  [task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end

% init the stimulus
myscreen = initStimulus('stimulus',myscreen);
% init the staircase
stimulus = initStair(stimulus);
% to initialize the stimulus for your experiment.
% stimulus = initGaussian(stimulus,myscreen);

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

%dispPsychometric(task{1}{1}, stimulus);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)
global stimulus
if task.thistrial.thisseg == 1
    if (stimulus.stairCon && (stimulus.stair.con.s.reversaln==20)) || ...
            (stimulus.stairDist && (stimulus.stair.dist.s.reversaln==20))
        stimulus.endflag=1;
    end
        
    stimulus.n = stimulus.n+1;

        if stimulus.n>1 && isnan(task.randVars.correct(stimulus.n-1))
        disp(sprintf('(posdisc) Trial %i: contrast %0.4f ecc %0.1f diff %0.4f No resp', ...
        stimulus.n-1, stimulus.contrast, stimulus.eccentricity-stimulus.fixOrigin(1), task.randVars.diff(stimulus.n-1)))
        end
    
    stimulus.fixColor = [1 1 1];
    task.thistrial.xpos = [0 0];
    
    task.thistrial.ecc = stimulus.eccentricity-stimulus.fixOrigin(1);
    
    switch stimulus.taskType
        case 1 % stair contrast
            task.thistrial.diff = stimulus.posDiff;
            task.thistrial.posDiff = stimulus.posDiff;
            
            stimulus.contrast = stimulus.stair.con.s.threshold;
            task.thistrial.con = stimulus.stair.con.s.threshold;
        case 4 % constant stimuli
            task.thistrial.diff = task.thistrial.posDiff;
            
            stimulus.contrast = task.thistrial.contrast;
            task.thistrial.con = task.thistrial.contrast;
        otherwise % fixed contrast/eccentricity
            task.thistrial.posDiff = stimulus.stair.dist.s.threshold;
            task.thistrial.diff = stimulus.stair.dist.s.threshold;
            
            task.thistrial.con=stimulus.contrast;
            
    end

%     task.thistrial.hemi = task.thistrial.whichHemifield;
    stimulus = initTarget(stimulus,myscreen);
    
    if task.thistrial.whichInterval == 1 % first one more eccentric
        task.thistrial.xpos = [stimulus.x + task.thistrial.posDiff/2, ...
            stimulus.x - task.thistrial.posDiff/2];
        task.thistrial.whichint = 1;
    else % second one more eccentric
        task.thistrial.xpos = [stimulus.x - task.thistrial.posDiff/2, ...
            stimulus.x + task.thistrial.posDiff/2];
        task.thistrial.whichint = 2;
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
    mglFixationCross(stimulus.fixWidth,1.5,stimulus.fixColor*.75,stimulus.fixOrigin);
else
    mglFixationCross(stimulus.fixWidth,1.5,stimulus.fixColor,stimulus.fixOrigin);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)
global stimulus
% here, we just check whether this is the first time we got a response
if ~task.thistrial.gotResponse
    if task.thistrial.posDiff && (task.thistrial.whichInterval == task.thistrial.whichButton)
        % correct
        task.thistrial.correct = 1;
        % feeback
        stimulus.fixColor = [0 1 0];
        
        disp(sprintf('(posdisc) Trial %i: contrast %0.4f ecc %0.1f diff %0.4f resp %i correct', ...
            task.trialnum, stimulus.contrast, stimulus.eccentricity-stimulus.fixOrigin(1), task.thistrial.posDiff, task.thistrial.whichButton))
    else
        % incorrect
        task.thistrial.correct = 0;
        stimulus.fixColor = [1 0 0];
        
        disp(sprintf('(posdisc) Trial %i: contrast %0.4f ecc %0.1f diff %0.4f resp %i incorrect', ...
            task.trialnum, stimulus.contrast, stimulus.eccentricity-stimulus.fixOrigin(1), task.thistrial.posDiff, task.thistrial.whichButton))
    end 
    task.thistrial.resp = task.thistrial.whichButton;
    
    if stimulus.stairCon
        [testValue, stimulus.stair.con] = doStaircase('getTestValue',stimulus.stair.con);
         stimulus.stair.con = doStaircase('update', stimulus.stair.con, task.thistrial.correct);
    elseif stimulus.stairDist
        [testValue, stimulus.stair.dist] = doStaircase('getTestValue',stimulus.stair.dist);
         stimulus.stair.dist = doStaircase('update', stimulus.stair.dist, task.thistrial.correct);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initTarget(stimulus,myscreen)

switch stimulus.stimType
    case 'gaussian'

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

    case 'gabor'
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
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initStair(stimulus)

% init Stair
if stimulus.stairCon
stimulus.stair.con = doStaircase('init','upDown','nup=1','ndown=3',...
    'initialThreshold=.50', 'initialStepsize=0.05', ...
    'minStepsize=0.005','maxStepsize=0.1','minThreshold=0.01','maxThreshold=1', ...
    'nTrials=100', 'dispFig=1', 'stepRule=Pest');
elseif stimulus.stairDist
    stimulus.stair.dist = doStaircase('init','upDown','nup=1','ndown=3',...
    'initialThreshold=5','initialStepsize=0.5', ...
    'minStepsize=0.005','maxStepsize=1','minThreshold=0','maxThreshold=5', ...
    'nTrials=100', 'dispFig=1', 'stepRule=Pest');
end

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
