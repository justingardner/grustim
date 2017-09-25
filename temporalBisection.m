function myscreen = temporalBisection(varargin)

clear global stimulus
mglEatKeys('12`');
global stimulus

% get arguments
high = 0; low = 0;
getArgs(varargin,{'high=0','low=0'},'verbose=1');

if high
	stimulus.tone.hz = 1700;
	stimulus.tone.duration = .01;
	stimulus.tone.amplitude = 0.6;
elseif low
	stimulus.tone.hz = 200;
	stimulus.tone.duration = .08;
	stimulus.tone.amplitude = 0.2;
else
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stimulus.disk.diameter = 2;
stimulus.disk.duration = 1/60;% .025;%1/60; % one(or two) frame 
stimulus.disk.contrast = .20;

stimulus.tone.samplesPerSecond = 44100;
% stimulus.tone.hz = 1700;
% stimulus.tone.duration = .01;

stimulus.standardOnset1 = 0.06; % 60ms
stimulus.midPoint= 0.46; % MIDPOINT ONSET 60ms + 400ms
stimulus.standardOnset3 = 0.86; % 60ms + 800ms
stimulus.earlyOnset1 = 0;
stimulus.earlyOnset3 = 0.8; % 800ms
stimulus.lateOnset1 = 0.12; % 120ms
stimulus.lateOnset3 = 0.92; % 120 + 800

stimulus.interval = [3];
% fixation cross
stimulus.fixWidth = 1;
stimulus.fixColor = [1 1 1];

stimulus.trialDur = 0.8 + 0.24 + 0.1;

stimulus.initDelay = 0.025;
stimulus.initDelaySd = 0.01;

% initalize the screen
myscreen.background = 0;
myscreen = initScreen;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%% 
task{1}{1}.waitForBacktick = 1;
task{1}{1}.segmin = [1 0.5 stimulus.trialDur 1.5 1];
task{1}{1}.segmax = [1 0.5 stimulus.trialDur 1.5 1];
task{1}{1}.getResponse = [0 0 0 1 0];

% parameters & randomization
task{1}{1}.randVars.uniform.closeTo = [1 3];
task{1}{1}.parameter.condition = {'vision','auditory','noOffset','posOffset','negOffset'};

% task{1}{1}.parameter.delay = [-0.1 0 0.1];

task{1}{1}.numTrials = 150;
task{1}{1}.random = 1;

task{1}{1}.randVars.calculated.resp = nan;
task{1}{1}.randVars.calculated.correct = nan;
task{1}{1}.randVars.calculated.probeDelay = nan;
task{1}{1}.randVars.calculated.noise = nan;
task{1}{1}.randVars.calculated.rt = nan;
task{1}{1}.randVars.calculated.condNum = nan;
task{1}{1}.randVars.calculated.hz = nan;
 
% initialize the task
for phaseNum = 1:length(task{1})
  [task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end
 
% init the stimulus
myscreen = initStimulus('stimulus',myscreen);

% init the staircase
stimulus = initStair(stimulus);

% to initialize the stimulus for your experiment.
% stimulus = initDisc(stimulus,myscreen);
% stimulus = initClick(stimulus,task);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% if stimulus.disp
dispPsychometric(task{1}{1});
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)
global stimulus
if task.thistrial.thisseg == 1
	stimulus.fixColor = [1 1 1];
	
	task.thistrial.condNum = find(strcmp(char(task.thistrial.condition),{'vision','auditory','noOffset','posOffset','negOffset'}));

	% get Test Value
	[testValue, stimulus.stair{task.thistrial.condNum}] = doStaircase('testValue', stimulus.stair{task.thistrial.condNum});
	if testValue > 0.04 - stimulus.tone.duration
		testValue = 0.04 - stimulus.tone.duration;
	end
	task.thistrial.noise = 0.080 * randn(1); % random number from a gaussian distribution with a std of 80ms
	while (stimulus.midPoint + (testValue + task.thistrial.noise) <= stimulus.lateOnset1+stimulus.tone.duration) ||...
	 (stimulus.midPoint + (testValue + task.thistrial.noise) >= stimulus.earlyOnset3-stimulus.tone.duration)
		task.thistrial.noise = 0.080 * randn(1);
	end
	if task.thistrial.closeTo == 1
		sign = -1;
	else
		sign = 1;
	end
		
	task.thistrial.probeDelay = sign * (testValue + task.thistrial.noise);

	if ~strcmp(char(task.thistrial.condition), 'vision')
		stimulus = initClick(stimulus,task);
		task.thistrial.hz = stimulus.thisHz;
	end

elseif task.thistrial.thisseg == 3
	stimulus.t0 = mglGetSecs;

elseif task.thistrial.thisseg == 5
	stimulus.fixColor = [1 1 1]*0.3;

elseif task.thistrial.thisseg == 4
	stimulus.fixColor = [1 1 0];%[1 1 1] * 0.2;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)
global stimulus
mglClearScreen(0);
mglFixationCross(stimulus.fixWidth,1.5,stimulus.fixColor);


if task.thistrial.thisseg == 3
	
	t0 = stimulus.t0;

	thisOnset2 = stimulus.midPoint + task.thistrial.probeDelay;

 switch char(task.thistrial.condition)
	case 'vision'
		thisOnset1 = stimulus.standardOnset1;
		thisOnset3 = stimulus.standardOnset3;

		if mglGetSecs(t0) >= thisOnset1 && mglGetSecs(t0) <= thisOnset1 + stimulus.disk.duration
			mglFillOval(0,0,[stimulus.disk.diameter stimulus.disk.diameter], [1 1 1] * stimulus.disk.contrast);
		elseif mglGetSecs(t0) > thisOnset1+stimulus.disk.duration && mglGetSecs(t0) < thisOnset2
			mglClearScreen;
		elseif mglGetSecs(t0) >= thisOnset2 && mglGetSecs(t0) <= thisOnset2 + stimulus.disk.duration
			mglFillOval(0,0,[stimulus.disk.diameter stimulus.disk.diameter],[1 1 1] * stimulus.disk.contrast);
		elseif mglGetSecs(t0) > thisOnset2 + stimulus.disk.duration && mglGetSecs(t0) < thisOnset3
			mglClearScreen;
		elseif mglGetSecs(t0) >= thisOnset3 && mglGetSecs(t0) <= thisOnset3 + stimulus.disk.duration
			mglFillOval(0,0,[stimulus.disk.diameter stimulus.disk.diameter],[1 1 1] * stimulus.disk.contrast);
		else
			mglClearScreen;
		end
	case 'auditory'
		mglClearScreen;
		mglPlaySound(stimulus.wav);
	case 'noOffset'
		thisOnset1 = stimulus.standardOnset1;
		thisOnset3 = stimulus.standardOnset3;

		mglPlaySound(stimulus.wav);

		if mglGetSecs(t0) >= thisOnset1 && mglGetSecs(t0) <= thisOnset1 + stimulus.disk.duration
			mglFillOval(0,0,[stimulus.disk.diameter stimulus.disk.diameter], [1 1 1] * stimulus.disk.contrast);
		elseif mglGetSecs(t0) > thisOnset1+stimulus.disk.duration && mglGetSecs(t0) < thisOnset2
			mglClearScreen;
		elseif mglGetSecs(t0) >= thisOnset2 && mglGetSecs(t0) <= thisOnset2 + stimulus.disk.duration
			mglFillOval(0,0,[stimulus.disk.diameter stimulus.disk.diameter],[1 1 1] * stimulus.disk.contrast);
		elseif mglGetSecs(t0) > thisOnset2 + stimulus.disk.duration && mglGetSecs(t0) < thisOnset3
			mglClearScreen;
		elseif mglGetSecs(t0) >= thisOnset3 && mglGetSecs(t0) <= thisOnset3 + stimulus.disk.duration
			mglFillOval(0,0,[stimulus.disk.diameter stimulus.disk.diameter],[1 1 1] * stimulus.disk.contrast);
		else
			mglClearScreen;
		end

	case 'posOffset' % auditory early, vision late
		thisOnset1 = stimulus.lateOnset1;
		thisOnset3 = stimulus.lateOnset3;

		mglPlaySound(stimulus.wav);

		if mglGetSecs(t0) >= thisOnset1 && mglGetSecs(t0) <= thisOnset1 + stimulus.disk.duration
			mglFillOval(0,0,[stimulus.disk.diameter stimulus.disk.diameter], [1 1 1] * stimulus.disk.contrast);
		elseif mglGetSecs(t0) > thisOnset1+stimulus.disk.duration && mglGetSecs(t0) < thisOnset2
			mglClearScreen;
		elseif mglGetSecs(t0) >= thisOnset2 && mglGetSecs(t0) <= thisOnset2 + stimulus.disk.duration
			mglFillOval(0,0,[stimulus.disk.diameter stimulus.disk.diameter],[1 1 1] * stimulus.disk.contrast);
		elseif mglGetSecs(t0) > thisOnset2 + stimulus.disk.duration && mglGetSecs(t0) < thisOnset3
			mglClearScreen;
		elseif mglGetSecs(t0) >= thisOnset3 && mglGetSecs(t0) <= thisOnset3 + stimulus.disk.duration
			mglFillOval(0,0,[stimulus.disk.diameter stimulus.disk.diameter],[1 1 1] * stimulus.disk.contrast);
		else
			mglClearScreen;
		end

	case 'negOffset' % auditory late, vision early
		thisOnset1 = stimulus.earlyOnset1;
		thisOnset3 = stimulus.earlyOnset3;

		mglPlaySound(stimulus.wav);

		if mglGetSecs(t0) >= thisOnset1 && mglGetSecs(t0) <= thisOnset1 + stimulus.disk.duration
			mglFillOval(0,0,[stimulus.disk.diameter stimulus.disk.diameter], [1 1 1] * stimulus.disk.contrast);
		elseif mglGetSecs(t0) > thisOnset1+stimulus.disk.duration && mglGetSecs(t0) < thisOnset2
			mglClearScreen;
		elseif mglGetSecs(t0) >= thisOnset2 && mglGetSecs(t0) <= thisOnset2 + stimulus.disk.duration
			mglFillOval(0,0,[stimulus.disk.diameter stimulus.disk.diameter],[1 1 1] * stimulus.disk.contrast);
		elseif mglGetSecs(t0) > thisOnset2 + stimulus.disk.duration && mglGetSecs(t0) < thisOnset3
			mglClearScreen;
		elseif mglGetSecs(t0) >= thisOnset3 && mglGetSecs(t0) <= thisOnset3 + stimulus.disk.duration
			mglFillOval(0,0,[stimulus.disk.diameter stimulus.disk.diameter],[1 1 1] * stimulus.disk.contrast);
		else
			mglClearScreen;
		end
 end
			
	% if mglGetSecs(t0) >= .06 && mglGetSecs(t0) <= .06 + stimulus.disk.duration
	% 	mglFillOval(0,0,[stimulus.disk.diameter stimulus.disk.diameter],[1 1 1] * stimulus.disk.contrast);
	% elseif mglGetSecs(t0) > .06+stimulus.disk.duration && mglGetSecs(t0) < .46
	% 	mglClearScreen;
	% elseif mglGetSecs(t0) >= .46 && mglGetSecs(t0) <= .46 + stimulus.disk.duration
	% 	mglFillOval(0,0,[stimulus.disk.diameter stimulus.disk.diameter],[1 1 1] * stimulus.disk.contrast);
	% elseif mglGetSecs(t0) > .46+stimulus.disk.duration && mglGetSecs(t0) < .86
	% 	mglClearScreen;
	% elseif mglGetSecs(t0) >= .86 && mglGetSecs(t0) <= .86 + stimulus.disk.duration
	% 	mglFillOval(0,0,[stimulus.disk.diameter stimulus.disk.diameter],[1 1 1] * stimulus.disk.contrast);
	% else 
	% 	mglClearScreen;
	% end

% else
% 	mglFixationCross(stimulus.fixWidth,1.5,stimulus.fixColor);
elseif task.thistrial.thisseg == 2
	mglClearScreen;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)
global stimulus
% here, we just check whether this is the first time we got a response
if ~task.thistrial.gotResponse

	if (task.thistrial.probeDelay < 0 && task.thistrial.whichButton == 1) || ... % closer to the first
		(task.thistrial.probeDelay > 0 && task.thistrial.whichButton == 2)  % closer to the third 
		% correct
		task.thistrial.correct = 1;
		% if any(task.thistrial.condNum == [1 2 3])	
		% % feedback
		% stimulus.fixColor = [0 1 0];
		% end
		disp(sprintf('(temporalBisection) %i:%s delay %0.3f resp %i correct', ...
            task.trialnum, char(task.thistrial.condition), task.thistrial.probeDelay, task.thistrial.whichButton))

	else
		% incorrect
		task.thistrial.correct = 0;
		% if any(task.thistrial.condNum == [1 2 3])
		% stimulus.fixColor = [1 0 0];
		% end
		disp(sprintf('(temporalBisection) %i:%s delay %0.3f resp %i incorrect', ...
            task.trialnum, char(task.thistrial.condition), task.thistrial.probeDelay, task.thistrial.whichButton))

	end
	stimulus.fixColor = [0 0 1];
	stimulus.stair{task.thistrial.condNum} = doStaircase('update', stimulus.stair{task.thistrial.condNum}, task.thistrial.correct, ...
		abs(task.thistrial.probeDelay));
	 
	task.thistrial.resp = task.thistrial.whichButton;
    task.thistrial.rt = task.thistrial.reactionTime;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initStair(stimulus)
	condNames = {'vision','auditory','noOffset','posOffset','negOffset'};
for cond = 1:5
	stimulus.stair{cond} = doStaircase('init','quest', 'initialThreshold', stimulus.initDelay, 'initialThresholdSd', stimulus.initDelaySd, ...
		'pThreshold', 0.75,'dispFig=1','subplotRows=5','subplotCols=1','subplotNum',cond,'subplotName',condNames{cond});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stimulus = initClick(stimulus,task);
fs = stimulus.tone.samplesPerSecond - 1;
t = 0 : 1/fs : stimulus.tone.duration;
stimulus.thisHz = stimulus.tone.hz + (stimulus.tone.hz/10)*randn;

wav = stimulus.tone.amplitude*sin(2*pi*stimulus.thisHz*t);
if isodd(length(t))
	mu = t((1+length(t))/2);
else
	mu = t(round((1+length(t))/2));
end
sigma = stimulus.tone.duration;
envelope = normpdf(t,mu,sigma);
%normalize
envelope = envelope / max(envelope);
waveform = wav .* envelope;

trialDur = 0:1/fs:stimulus.trialDur;
thiswav = zeros(1,length(trialDur));

thisOnset2 = stimulus.midPoint + task.thistrial.probeDelay;
switch char(task.thistrial.condition)
	case {'auditory','noOffset'}
		thisOnset1 = stimulus.standardOnset1;
		thisOnset3 = stimulus.standardOnset3;
	case 'posOffset'
		thisOnset1 = stimulus.earlyOnset1;
		thisOnset3 = stimulus.earlyOnset3;
	case 'negOffset'
		thisOnset1 = stimulus.lateOnset1;
		thisOnset3 = stimulus.lateOnset3;
end

if ~thisOnset1
	onset1 = 0;
else
	onset1 = length(0:1/fs:thisOnset1);
end
thiswav(onset1+1:onset1+length(waveform)) = waveform;

onset2 = length(0:1/fs:thisOnset2);
thiswav(onset2+1:onset2+length(waveform)) = waveform;

onset3 = length(0:1/fs:thisOnset3);
thiswav(onset3+1:onset3+length(waveform)) = waveform;

stimulus.wav = mglInstallSound(thiswav, stimulus.tone.samplesPerSecond);

% trialDur = 0:1/fs:stimulus.trialDur;
% wav = zeros(1,length(trialDur));
% start = length(0:1/fs:.06);
% wav(start:start+length(waveform)-1) = waveform;
% start = length(0:1/fs:.46);
% wav(start:start+length(waveform)-1) = waveform;
% start = length(0:1/fs:.86);
% wav(start:start+length(waveform)-1) = waveform;

% stimulus.wav = mglInstallSound(wav, stimulus.tone.samplesPerSecond);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display psychometric functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispPsychometric(task)
uni = figure; bi = figure;
condNames = {'vision','auditory','noOffset','posOffset','negOffset'};
% condition = {};
% for b = 1:length(task.block)
%   condition = [condition, task.block(b).parameter.condition];
% end
condition = task.randVars.condNum;

for cond = 1:5
  thisTrials = find(condition == cond);%find(strcmp(condition, condNames{cond}));
  thisDelay = task.randVars.probeDelay(thisTrials);
  thisResp = task.randVars.resp(thisTrials);
  isLate = (thisResp == 2);

  binCenter = -0.175:.05:0.175; space = max(diff(binCenter));
  for b = 1:length(binCenter)
    switch b
      case 1
        binnedTrial{b} = find(thisDelay<= min(binCenter)+space/2); % -200 ~ -150 ms
      case length(binCenter)
        binnedTrial{b} = find(thisDelay > max(binCenter)-space/2);
      otherwise
        binnedTrial{b} = find((thisDelay > min(binCenter) + space*(b-2)) & (thisDelay <= min(binCenter) + space*(b-1)));
    end
  end
  
  nTotalBin = cellfun(@(x) length(x), binnedTrial);
  binnedTrial(nTotalBin==0) = [];
  binCenter(nTotalBin==0) = [];
  nTotalBin(nTotalBin==0) = [];

  nLateBin = cellfun(@(x) sum(isLate(x)), binnedTrial);
  pLateBin = nLateBin./nTotalBin;

  % xs = -0.2:.001:0.2;
  % fit{cond} = fitCG(binCenter, nTotalBin, nLateBin);

  % mu = fit{cond}.fitparams(1);
  % sigma = fit{cond}.fitparams(2);
  % lambda = fit{cond}.fitparams(3);
  % cgfunc = fit{cond}.cgfunc;
  % fitVal = arrayfun(@(x) cgfunc(mu,sigma,lambda,x), xs);
  % pfitValBin = arrayfun(@(x) cgfunc(mu,sigma,lambda,x), binCenter);
  % error = sqrt(pfitValBin.*(1-pfitValBin)./nTotalBin);

  switch cond
    case {1,2} 
      figure(uni);
      hold on;
      myerrorbar(binCenter, pLateBin,'Symbol',getsymbol(cond),'Color',getcolor(cond), 'MarkerSize',6);
      % myerrorbar(binCenter+(rand-0.5)*0.002, pLateBin,'yError',error, 'Symbol',getsymbol(cond),'Color',getcolor(cond), 'MarkerSize',6);
      % plot(xs, fitVal, 'Color',getcolor(cond),'LineWidth',0.8);
    otherwise
    	figure(bi);
     hold on;
      myerrorbar(binCenter, pLateBin,'Symbol',getsymbol(cond),'Color',getcolor(cond), 'MarkerSize',6);
      % myerrorbar(binCenter+(rand-0.5)*0.002, pLateBin,'yError',error, 'Symbol',getsymbol(cond),'Color',getcolor(cond), 'MarkerSize',6);
      % plot(xs, fitVal, 'Color',getcolor(cond),'LineWidth',0.8);
  end

end

figure(bi); mylegend({'No offset','Pos offset','Neg offset'}, {{getcolor(3)},{ getcolor(4)},{ getcolor(5)}});
figure(uni); mylegend({'vision','auditory'},{{getcolor(1)},{getcolor(2)}});
