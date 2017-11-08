function myscreen = flashIllusion(varargin)

clear global stimulus
mglEatKeys('1234`');
global stimulus

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stimulus.disc.width = 2;
stimulus.disc.eccentricity = 5;
stimulus.disc.color = [1 1 1];
stimulus.disc.dur = 1/60;%0.01;
stimulus.disc.isi = 0.05;

stimulus.beep.hz = 440;
stimulus.beep.samplesPerSecond = 44100;
stimulus.beep.dur = 1/60;%0.01;
stimulus.beep.isi = 0.057;

stimulus.stimsegDur = 0.75;

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
task{1}{1}.segmin = [0.5 0.25 stimulus.stimsegDur 2 1];
task{1}{1}.segmax = [0.5 0.75 stimulus.stimsegDur 2 1];
task{1}{1}.getResponse = [0 0 0 1 0];

% parameters & randomization
task{1}{1}.parameter.numFlash = [1:4];
task{1}{1}.parameter.numBeep = [0:4];
task{1}{1}.numBlocks = 5;
task{1}{1}.random = 1;

task{1}{1}.randVars.calculated.resp = nan;
task{1}{1}.randVars.calculated.rt = nan;
task{1}{1}.randVars.calculated.nbeep = nan;
task{1}{1}.randVars.calculated.nflash = nan;

% initialize the task
for phaseNum = 1:length(task{1})
  [task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end
 
% init the stimulus
myscreen = initStimulus('stimulus',myscreen);

stimulus = initAuditory(stimulus);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)
global stimulus
if task.thistrial.thisseg == 1

	stimulus.fixColor = [1 1 1];
	if task.thistrial.numBeep
		stimulus = updateAuditory(stimulus,task);
	end
	% create sequence for visual separately
	stimulus.disc.on = [];
	for i = 1:task.thistrial.numFlash
		stimulus.disc.on(i) = (stimulus.disc.dur+stimulus.disc.isi)*(i-1);
	end

	task.thistrial.nflash = task.thistrial.numFlash;
	task.thistrial.nbeep = task.thistrial.numBeep;
		
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
mglClearScreen(0);
mglFixationCross(stimulus.fixWidth,1.5,stimulus.fixColor);

if task.thistrial.thisseg == 3
	if task.thistrial.numBeep
		mglPlaySound(stimulus.wav);
	end
	t0 = stimulus.t0;
	for i = 1:task.thistrial.numFlash
		while mglGetSecs(t0) >= stimulus.disc.on(i) && mglGetSecs(t0) < stimulus.disc.on(i) + stimulus.disc.dur
			mglFillOval(stimulus.disc.eccentricity,0,[stimulus.disc.width,stimulus.disc.width], stimulus.disc.color);
		% else
			
			% mglFixationCross(stimulus.fixWidth,1.5,stimulus.fixColor);
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

	stimulus.fixColor = [0 1 1];

	task.thistrial.resp = task.thistrial.whichButton;
    task.thistrial.rt = task.thistrial.reactionTime;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stimulus = initAuditory(stimulus)
fs = stimulus.beep.samplesPerSecond - 1;
t = 0 : 1/fs : stimulus.beep.dur;
stimulus.beep.beep = sin(2*pi*stimulus.beep.hz*t);

function stimulus = updateAuditory(stimulus,task)
fs = stimulus.beep.samplesPerSecond - 1;
totalDur = length(0:1/fs:stimulus.stimsegDur);
beepDur = length(0 : 1/fs : stimulus.beep.dur);
isi = length(0:1/fs:stimulus.beep.isi);

thiswav = zeros(1,totalDur);
for b = 1:task.thistrial.numBeep
	thisStart = (beepDur+isi)*(b-1) + 1;
	thiswav(thisStart:thisStart+beepDur-1) = stimulus.beep.beep;
end
stimulus.wav = mglInstallSound(thiswav, stimulus.beep.samplesPerSecond);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display psychometric functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotResults(task)
% single flash + 1~4 beeps
% y = perceived flash


%number of flashes * actual number of flashes (no sound)
% y = perceived flashes

% number of flashes  + 1 beep (catch trial)


figure;
subplot(2,1,1);
singleFlash = (task.randVars.nflash==1);
for n = 1:4
	nBeeps{n} = (task.randVars.nbeep == n);
	perceived{n} = task.randVars.resp(singleFlash && nBeeps{n});
end


task.randVars.resp

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
