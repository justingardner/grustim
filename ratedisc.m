function myscreen = ratedisc(varargin)

clear global stimulus
mglEatKeys('12`');
global stimulus

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stimulus.auditory.hz = 220;
stimulus.auditory.samplesPerSecond = 44100;
stimulus.auditory.ampRatio = sqrt(10)^((50.1-41.1)/10);
stimulus.auditory.highAmp = 0.2;
stimulus.auditory.lowAmp = stimulus.auditory.highAmp/stimulus.auditory.ampRatio;
stimulus.auditory.noiseRatio = sqrt(10)^((54.9-50.1)/10); % to high amp
stimulus.auditory.noiseAmp = 0.15;%stimulus.auditory.highAmp * stimulus.auditory.noiseRatio;
stimulus.auditory.duration = 0.01;

stimulus.visual.width = 10;
stimulus.visual.lumRatio = 66.5/14.5;
stimulus.visual.highLumRange = [0.2 0.3];
stimulus.visual.lowLumRange = stimulus.visual.highLumRange/stimulus.visual.lumRatio;
stimulus.visual.noiseRatio = sqrt(10)^(1/10); % to high amp
stimulus.visual.gaussianLum = 0.4;%mean(stimulus.visual.highLumRange) / stimulus.visual.noiseRatio;
stimulus.visual.duration = 0.01; % 10ms
stimulus.visual.backgroundLoc = [0 7.5];

stimulus.interval.short = 0.06;
stimulus.interval.long = 0.120;

stimulus.rate.boundary = 10.5; % events/s
stimulus.rate.range = 7:14;
stimulus.rate.pairLow = 7:12;
stimulus.rate.delta = 2;

stimulus.eventDur = 1;

% fixation cross
stimulus.fixWidth = 1;
stimulus.fixColor = [0 0 0];

% initalize the screen
myscreen.background = 0;
myscreen = initScreen;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%% 
task{1}{1}.waitForBacktick = 1;
task{1}{1}.segmin = [0.5 0.25 stimulus.eventDur 2 1];
task{1}{1}.segmax = [0.5 0.75 stimulus.eventDur 2 1];
task{1}{1}.getResponse = [0 0 0 1 0];

% parameters & randomization
task{1}{1}.parameter.condition = {'visual','auditory','noConflict','posConflict','negConflict'};
task{1}{1}.parameter.rate = stimulus.rate.range;
task{1}{1}.parameter.visualRel = [0,1];
task{1}{1}.parameter.auditoryRel = [0 1];
% negConflict -> visual = auditory -2;
% posConflict -> visual = auditory + 2;
% task{1}{1}.numTrials = 140;
task{1}{1}.random = 1;

task{1}{1}.randVars.calculated.resp = nan;
task{1}{1}.randVars.calculated.correct = nan;
task{1}{1}.randVars.calculated.actualrate = nan;
task{1}{1}.randVars.calculated.visualrate = nan;
task{1}{1}.randVars.calculated.auditoryrate = nan;
task{1}{1}.randVars.calculated.noise = nan;
task{1}{1}.randVars.calculated.rt = nan;
task{1}{1}.randVars.calculated.condNum = nan;
task{1}{1}.randVars.calculated.luminance = nan;

% initialize the task
for phaseNum = 1:length(task{1})
  [task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end
 

% init the stimulus
myscreen = initStimulus('stimulus',myscreen);

stimulus = initAuditory(stimulus);
stimulus = initVisual(stimulus);

% to initialize the stimulus for your experiment.
% stimulus = initDisc(stimulus,myscreen);
% stimulus = initClick(stimulus,task);

mglStencilCreateBegin(1);
mglFillRect(0,7.5,[10 10]);
mglFillOval(0,0,[2 2]);
mglStencilCreateEnd;
mglClearScreen;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)
global stimulus
if task.thistrial.thisseg == 1

	stimulus.fixColor = [1 1 1];

	task.thistrial.condNum = find(strcmp(char(task.thistrial.condition),{'visual','auditory','noConflict','posConflict','negConflict'}));

	switch char(task.thistrial.condition)
	case {'visual','auditory','noConflict'}
		task.thistrial.actualrate = task.thistrial.rate;
		task.thistrial.auditoryrate = task.thistrial.actualrate;
		task.thistrial.visualrate = task.thistrial.actualrate; 
		% stimulus = generateSequence(stimulus, task.thistrial.rate);
	case 'posConflict'
		task.thistrial.auditoryrate = stimulus.rate.pairLow(randi(length(stimulus.rate.pairLow)));
		task.thistrial.visualrate = task.thistrial.auditoryrate + stimulus.rate.delta;
	case 'negConflict'
		task.thistrial.visualrate = stimulus.rate.pairLow(randi(length(stimulus.rate.pairLow)));
		task.thistrial.auditoryrate = task.thistrial.visualrate + stimulus.rate.delta;
		% negConflict -> visual = auditory -2;
		% posConflict -> visual = auditory + 2;
	end
	stimulus = updateAuditory(stimulus,task);
	% create sequence for visual separately
	stimulus = generateSequence(stimulus,task.thistrial.visualrate);
	stimulus.visual.eventOn = []; stimulus.visual.intOn = [];
	stimulus.visual.eventOn(1) = stimulus.interval.delay;
	for i = 1:length(stimulus.interval.sequence)
		stimulus.visual.intOn(i) = stimulus.visual.eventOn(i) + stimulus.visual.duration;
		stimulus.visual.eventOn(i+1) = stimulus.visual.intOn(i) + stimulus.interval.sequence(i);
	end

	if task.thistrial.visualRel == 1
		stimulus.thistex = stimulus.visual.hightex;
	else
		stimulus.thistex = stimulus.visual.lowtex;
	end
		
elseif task.thistrial.thisseg == 3
	stimulus.t0 = mglGetSecs;

elseif task.thistrial.thisseg == 4
	stimulus.fixColor = [1 1 1]*0.3;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)
global stimulus
mglClearScreen(0.5);
mglFixationCross(stimulus.fixWidth,1.5,stimulus.fixColor);

if task.thistrial.thisseg == 3
	t0 = stimulus.t0;

 switch char(task.thistrial.condition)
	case 'visual'
		mglStencilSelect(1);
				mglBltTexture(stimulus.visual.noisetex, stimulus.visual.backgroundLoc);
		
		for i = 1:length(stimulus.visual.eventOn)
			
			while mglGetSecs(t0) >= stimulus.visual.eventOn(i) && mglGetSecs(t0) < stimulus.visual.eventOn(i) + stimulus.visual.duration
				mglStencilSelect(1);
				mglBltTexture(stimulus.visual.noisetex, stimulus.visual.backgroundLoc);
				mglBltTexture(stimulus.thistex,stimulus.visual.backgroundLoc);
			% else
			% 	mglStencilSelect(1);
			% 	mglBltTexture(stimulus.visual.noisetex, stimulus.visual.backgroundLoc);
			end
		end
 
	case 'auditory'
		% mglFixationCross(stimulus.fixWidth,1.5,stimulus.fixColor);
		mglPlaySound(stimulus.wav);
	otherwise
		mglPlaySound(stimulus.wav);
		mglStencilSelect(1);
				mglBltTexture(stimulus.visual.noisetex, stimulus.visual.backgroundLoc);

		for i = 1:length(stimulus.visual.eventOn)
			while mglGetSecs(t0) >= stimulus.visual.eventOn(i) && mglGetSecs(t0) < stimulus.visual.eventOn(i) + stimulus.visual.duration
				mglStencilSelect(1);
				mglBltTexture(stimulus.visual.noisetex, stimulus.visual.backgroundLoc);
				mglBltTexture(stimulus.thistex,stimulus.visual.backgroundLoc);
			% else
			% 	mglStencilSelect(1);
			% 	mglBltTexture(stimulus.visual.noisetex, stimulus.visual.backgroundLoc);
			end
		end

	end
	
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

	switch char(task.thistrial.condition)
	case {'visual','auditory','noConflict'}

	if (task.thistrial.rate > stimulus.rate.boundary && task.thistrial.whichButton == 2) || ... % closer to the first
		  (task.thistrial.rate < stimulus.rate.boundary && task.thistrial.whichButton == 1)% closer to the third 
		% correct
		task.thistrial.correct = 1;
		% feedback
		stimulus.fixColor = [0 1 0];
		disp(sprintf('(ratedisc) %i:%s rate %i resp %i correct', ...
            task.trialnum, char(task.thistrial.condition), task.thistrial.rate, task.thistrial.whichButton))

	else
		% incorrect
		task.thistrial.correct = 0;
		stimulus.fixColor = [1 0 0];
		disp(sprintf('(ratedisc) %i:%s rate %i resp %i incorrect', ...
            task.trialnum, char(task.thistrial.condition), task.thistrial.rate, task.thistrial.whichButton))

	end
	otherwise
		% random feedback
		task.thistrial.correct = randi([0 1]);
		if task.thistrial.correct
			stimulus.fixColor = [0 1 0];
		else
			stimulus.fixColor = [1 0 0];
		end
		disp(sprintf('(ratedisc) %i:%s rate %i resp %i randomfeedback %i', ...
            task.trialnum, char(task.thistrial.condition), task.thistrial.rate, task.thistrial.whichButton, task.thistrial.correct))

	end

	task.thistrial.resp = task.thistrial.whichButton;
    task.thistrial.rt = task.thistrial.reactionTime;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stimulus = initAuditory(stimulus)
fs = stimulus.auditory.samplesPerSecond - 1;
t = 0 : 1/fs : stimulus.auditory.duration;
event = sin(2*pi*stimulus.auditory.hz*t);
stimulus.auditory.lowEvent = event * stimulus.auditory.lowAmp;
stimulus.auditory.highEvent = event * stimulus.auditory.highAmp;
totalDur = 0:1/fs:stimulus.eventDur;
stimulus.auditory.noise = stimulus.auditory.noiseAmp * randn(1,length(totalDur));

fc = 2000; % cutoff frequency
% 5th order Butterworth filter
% [b,a] = butter(5, fc/(stimulus.tone.samplesPerSecond/2));
a = [1	-4.07876493416512	6.72527084144657	-5.59474636818042	2.34559680959441	-0.396133028715511];
b = [3.82287493728255e-05	0.000191143746864127	0.000382287493728255	0.000382287493728255	0.000191143746864127	3.82287493728255e-05];

stimulus.auditory.noise = filter(b,a,stimulus.auditory.noise);

stimulus.auditory.event = event;

function stimulus = updateAuditory(stimulus,task)
fs = stimulus.auditory.samplesPerSecond - 1;
totalDur = 0:1/fs:stimulus.eventDur;
beepDur = 0 : 1/fs : stimulus.auditory.duration;
thiswav = stimulus.auditory.noise;%zeros(1,length(totalDur));
stimulus = generateSequence(stimulus, task.thistrial.auditoryrate);
if task.thistrial.auditoryRel == 1
	thisbeep = stimulus.auditory.highEvent;
else
	thisbeep = stimulus.auditory.lowEvent;
end
onset = length(0:1/fs:stimulus.interval.delay) + 1;
% first event
thiswav(onset:onset+length(beepDur)-1) = thiswav(onset:onset+length(beepDur)-1)+thisbeep;
nextonset = onset+length(beepDur);
% second - last
for i = 1:length(stimulus.interval.sequence)
	thisIntLength = length(0:1/fs:stimulus.interval.sequence(i));
	nextonset = nextonset + thisIntLength;
	if nextonset+length(beepDur)-1 > length(thiswav)
		mglClose;
		dbquit
	end
	thiswav(nextonset:nextonset+length(beepDur)-1) = thiswav(nextonset:nextonset+length(beepDur)-1)+thisbeep;
	nextonset = nextonset+length(beepDur);
end
stimulus.wav = mglInstallSound(thiswav, stimulus.auditory.samplesPerSecond);
stimulus.thiswav = thiswav;
% trialDur = 0:1/fs:stimulus.trialDur;
% wav = zeros(1,length(trialDur));
% start = length(0:1/fs:.06);
% wav(start:start+length(waveform)-1) = waveform;
% start = length(0:1/fs:.46);
% wav(start:start+length(waveform)-1) = waveform;
% start = length(0:1/fs:.86);
% wav(start:start+length(waveform)-1) = waveform;

% stimulus.wav = mglInstallSound(wav, stimulus.tone.samplesPerSecond);

function stimulus = initVisual(stimulus);
xDeg2pix = mglGetParam('xDeviceToPixels');
yDeg2pix = mglGetParam('yDeviceToPixels');
width = 10; height = 10;
% get size in pixels
widthPixels = round(width*xDeg2pix);
heightPixels = round(height*yDeg2pix);
widthPixels = widthPixels + mod(widthPixels+1,2);
heightPixels = heightPixels + mod(heightPixels+1,2);

stimulus.visual.noisetex = mglCreateTexture(round(stimulus.visual.gaussianLum*rand(widthPixels,heightPixels)*255));
rect = zeros(widthPixels, heightPixels, 4);
for i = 1:3
	rect(:,:,i) = ones(widthPixels ,heightPixels)*255;
end
rect(:,:,4) = stimulus.visual.highLumRange(2) * 255;
stimulus.visual.hightex = mglCreateTexture(rect);
rect(:,:,4) = stimulus.visual.lowLumRange(1) * 255;
stimulus.visual.lowtex = mglCreateTexture(rect);


function stimulus = generateSequence(stimulus, rate)
% rate = task.thistrial.rate;
stimDur = rate * stimulus.auditory.duration;
blankDur = stimulus.eventDur - stimDur;
nInterval = rate - 1;
switch rate
case {7,8,9}
	nshort = 2*(rate-7);
otherwise
	nshort = 2*(rate-8) + 1;
end
nlong = nInterval - nshort;
intervalSeq = [repmat(stimulus.interval.short, [1 nshort]), repmat(stimulus.interval.long, [1 nlong])];
stimulus.interval.sequence = intervalSeq(randperm(nInterval));
stimulus.interval.nInterval = nInterval;
stimulus.interval.blankDur = blankDur;
stimulus.interval.delay = (blankDur - sum(stimulus.interval.sequence)) * rand;

% % short:long
% 7 ->6 0:6 120*6 = 720
% 8 ->7 2:5 60*2 + 120*5 = 720

% 9 ->8 4:4 60*4 + 120*4 = 720

% 10->9 4:5 60*5 + 120*4 = 780
% 11->10 6:4 60*7 + 120*3 = 780
% 12->11 8:3 60*9 + 120*2 =780
% 13->12 10:2 60*11 + 120*1 = 780
% 14 ->13 12:1 60*13 + 12*0 = 780


% % short:long
% 7 ->6 0:6 120*6 = 720
% 8 ->7 2:5 60*2 + 120*5 = 720

% 9 ->8 4:4 60*4 + 120*4 = 720
% % 9 ->8 2:6 60*2 + 120*6 = 840

% 10->9 4:5 60*4 + 120*5 = 840
% 11->10 6:4 60*6 + 120*4 = 840
% 12->11 8:3 60*8 + 120*3 = 840
% 13->12 10:2 60*10 + 120*2 = 840
% 14 ->13 12:1 60*12 + 12*1 = 840
% 15 ->14 14:0 60*14 + 120*0 = 840


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
