function myscreen = spjdg(varargin)

clear global stimulus
mglEatKeys('12`sf-=cd');
global stimulus

% get arguments
high = 0; low = 0; bimodal=0; visual=0; auditory=0; tenbit = 1; practice = 0; auditoryTrain = 0; visualTrain=0; 
visualMan=0; auditoryMan=0; spencer=[]; noiseHigh=0; noiseLow=0; unimodal=0; easy=0;
getArgs(varargin,{'high=0','low=0','bimodal=0','visual=0','auditory=0','tenbit=1','practice=0',...
  'auditoryTrain=0','visualTrain=0','visualMan=0','auditoryMan=0','spencer=[]','noiseHigh=0','noiseLow=0','unimodal=0','easy=0'},'verbose=1');
if isempty(spencer) || ~any(spencer==[0 1])
    error('(spjdg) spencer?');
end
stimulus.gaussian.contrast = 0.05;
stimulus.gaussian.diameterHighRel = 6;
stimulus.gaussian.diameterLowRel = 36;
% % if ~noiseHigh && ~noiseLow
% stimulus.noiseHighRel = 0.150;
% stimulus.noiseLowRel = 0.55;
% % % else
% %   stimulus.noiseHighRel = noiseHigh;
% %   stimulus.noiseLowRel = noiseLow;
% % end`

if ~noiseHigh && ~noiseLow
stimulus.noiseHighRel = 0.15;
stimulus.noiseLowRel = 0.65;
else
  stimulus.noiseHighRel = noiseHigh;
  stimulus.noiseLowRel = noiseLow;
end

if visualMan
  stimulus.noiseContrast = 0;
  stimulus.gaussian.diameter = stimulus.gaussian.diameterHighRel;
  stimulus.offset = 10;
  stimulus.visrel = 1;
elseif easy
    stimulus.noiseContrast = 0;
    stimulus.gaussian.diameterLowRel = 10;
    stimulus.gaussian.diameter = stimulus.gaussian.diameterLowRel;
    stimulus.gaussian.contrast = 0.5;
    stimulus.offset = 10;
  stimulus.visrel = 2;
else

if high
    stimulus.gaussian.diameter = stimulus.gaussian.diameterHighRel;
    stimulus.noiseContrast = stimulus.noiseHighRel;
elseif low
    stimulus.gaussian.diameter = stimulus.gaussian.diameterLowRel;
    stimulus.noiseContrast = stimulus.noiseLowRel;

else 
    if auditoryTrain || auditory || auditoryMan
        tenbit =0;
        stimulus.gaussian.diameter = nan;
        stimulus.gaussian.contrast = nan;
      stimulus.noiseContrast = stimulus.noiseHighRel;
    % else
    %     return
    end
    if auditoryMan
      stimulus.offset = 10;
      stimulus.visrel = 1;
    end
end

end
stimulus.high = high;
stimulus.low = low;
stimulus.tenbit = tenbit;
stimulus.practice = practice;
stimulus.auditoryTrain = auditoryTrain;
stimulus.visualTrain = visualTrain;
stimulus.bimodal = bimodal;
stimulus.visual = visual;
stimulus.auditory = auditory;
stimulus.visualMan = visualMan;
stimulus.auditoryMan = auditoryMan;
stimulus.spencer= spencer;
stimulus.unimodal=unimodal;
stimulus.easy = easy;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%f%%%%%%%%%%%%%%%%%%%%%%
stimulus.gaussian.duration = .015;%.015;% .025;%1/60; % one(or two) frame 
    if ~stimulus.tenbit
        stimulus.gaussian.contrast = .1;
    end
stimulus.colors.reservedColors = [1 1 1; 0.2 0.2 0.2; 0 1 0; 1 0 0; 0 1 1; 0 0 1; 0 0.75 0.75; 0.75 0.75 0.75];

stimulus.tone.samplesPerSecond = 44100;
stimulus.tone.duration = .0015;

screenParams = mglGetScreenParams;
stimulus.fps = screenParams{1}.framesPerSecond;
stimulus.flickerRate = 2;
stimulus.screenupdate = (stimulus.flickerRate/stimulus.fps);

% stimulus.pos1 = -17.5;
% stimulus.pos3 = 17.5;
stimulus.midPoint = 0; %(stimulus.pos1 + stimulus.pos3)/2;
stimulus.pos1 = stimulus.midPoint;

% stimulus.pos1 = -7.5;
% stimulus.pos3 = stimulus.pos1+30;
% stimulus.midPoint = (stimulus.pos1 + stimulus.pos3)/2; %7.5 deg from center
stimulus.delta=2.5;

% fixation cross
stimulus.fixWidth = 2.5;
stimulus.fixColor = [1 1 1];
stimulus.fixYpos = -4;
if ~(visualMan || auditoryMan || auditoryTrain || visualTrain || easy)
stimulus.coordshift = 15;
else
stimulus.coordshift = 0;
end

stimulus.initOffset = log10(5 + stimulus.coordshift);
stimulus.initOffsetSd = log10(2);

if stimulus.auditoryTrain || stimulus.visualTrain
    stimulus.initialThreshold = 10;
    stimulus.initialStepsize = 2.5;
    stimulus.minThreshold = 0;
    stimulus.maxThreshold = 15;
    stimulus.minStepsize = 0.75;
    stimulus.maxStepsize = 5;
end

stimulus.transientDur = 1;
stimulus.nRefresh = stimulus.transientDur*stimulus.fps;
if ~isinteger(stimulus.nRefresh)
  stimulus.nRefresh = round(stimulus.nRefresh);
end
still = round(stimulus.nRefresh/2);

transpi = linspace(0,round(255*stimulus.noiseHighRel),still);
stimulus.transpi.highRel = transpi;
stimulus.transpi.highRel(still+1:stimulus.nRefresh+1) = max(stimulus.transpi.highRel);

transpi = linspace(0,round(255*stimulus.noiseLowRel),still);
stimulus.transpi.lowRel = transpi;
stimulus.transpi.lowRel(still+1:stimulus.nRefresh+1) = max(stimulus.transpi.lowRel);

transpd = linspace(round(255*stimulus.noiseHighRel),0,still);
stimulus.transpd.highRel = transpd;
stimulus.transpd.highRel(still+1:stimulus.nRefresh+1) = min(stimulus.transpd.highRel);

transpd = linspace(round(255*stimulus.noiseLowRel),0,still);
stimulus.transpd.lowRel = transpd;
stimulus.transpd.lowRel(still+1:stimulus.nRefresh+1) = min(stimulus.transpd.lowRel);


% initalize the screen
myscreen.keyboard.nums = [124,125,127,126,2,4,50,37,28,25, 9, 3]; % 1 2 up down s f space enter - =(+) (c=red) (d=green)  **124 left  125 right   19 1, 20 2
myscreen = initScreen(myscreen);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%% 
task{1}{1}.waitForBacktick = 1;
stimDur = stimulus.gaussian.duration;
if ~(stimulus.visualMan || stimulus.auditoryMan || stimulus.easy)
task{1}{1}.segmin = [1 0.5 stimDur 0.5 stimDur inf 0.5];
task{1}{1}.getResponse = [0 0 0 0 0 1 0];
elseif stimulus.visualMan
task{1}{1}.segmin = [2 0.5 stimDur*15 0.5*8 stimDur*15 inf inf];
task{1}{1}.getResponse = [0 0 0 0 0 1 1];
elseif stimulus.auditoryMan
  task{1}{1}.segmin = [2 0.5 stimDur 0.5*4 stimDur inf inf];
task{1}{1}.getResponse = [0 0 0 0 0 1 1];
elseif stimulus.easy
    task{1}{1}.segmin = [1 0.5 1 0.5 1 inf 0.5];
task{1}{1}.getResponse = [0 0 0 0 0 1 0];
    

end

task{1}{1}.segmax = task{1}{1}.segmin;


% parameters & randomization
if stimulus.practice
    task{1}{1}.parameter.condition = {'vision','auditory','noOffset'};
    task{1}{1}.numTrials = 10*length(task{1}{1}.parameter.condition);
    task{1}{1}.randVars.uniform.closeTo = [1 2];
elseif stimulus.auditoryTrain || stimulus.auditoryMan
%     task{1}{1}.numTrials = 30;
    task{1}{1}.parameter.condition = {'auditory'};
%     task{1}{1}.parameter.offset = [7.5 10];
    task{1}{1}.parameter.closeTo = [1 2];
    if stimulus.auditoryMan
    task{1}{1}.randVars.uniform.visRel = stimulus.visrel;
    end
    % figure;
    % xlabel('trial'); ylabel('correct'); 
    % xaxis([0 30]); yaxis([-1 2]);
elseif stimulus.visualTrain || stimulus.visualMan || stimulus.easy
  task{1}{1}.parameter.condition = {'vision'};
  task{1}{1}.parameter.closeTo = [1 2];
  if stimulus.visualMan || stimulus.easy
    task{1}{1}.randVars.uniform.visRel = stimulus.visrel;
   end

else

nVisRel = 2;
condNames = {'vision','auditory','posOffset','negOffset'};
cond2 = {'vision','posOffset','negOffset'};
nCond = length(condNames);
nRepeat = 10;
stimulus.nTrialTotal = nVisRel*(nCond-1)*nRepeat + nRepeat;
visRel = [repmat(1,1,(nCond)*nRepeat), repmat(0,1,(nCond-1)*nRepeat)];
condition = [repmat(condNames,1,nRepeat), repmat(cond2,1,nRepeat)];  %,'auditory','posOffset','negOffset'
randSeq = randperm(length(condition));%randperm(nVisRel*nCond*nRepeat);

task{1}{1}.numTrials = length(condition);%nVisRel*nCond*nRepeat;
stimulus.nTrialTotal = task{1}{1}.numTrials;
 %  task{1}{1}.parameter.visRel = [1,0]; % high, low
 %    task{1}{1}.parameter.condition = {'vision','auditory','noOffset','posOffset','negOffset'};
  % task{1}{1}.numTrials = 10*length(task{1}{1}.parameter.condition);
% using my own random sequence
  task{1}{1}.randVars.visRel = visRel(randSeq);
  task{1}{1}.randVars.condition = condition(randSeq);

  task{1}{1}.randVars.uniform.closeTo = [1 2];


% if stimulus.high || stimulus.auditory
%   task{1}{1}.randVars.uniform.visRel = 1;
% elseif stimulus.low
%   task{1}{1}.randVars.uniform.visRel = 0;
% end
%   task{1}{1}.parameter.offset = [-7.5 -2.5 -1.25 0 1.25 2.5 7.5];
%   if stimulus.visual
% 	task{1}{1}.parameter.condition = {'vision'};
%     task{1}{1}.parameter.closeTo = [1 2];
%     task{1}{1}.numBlocks = 6;
%   elseif stimulus.auditory
%      task{1}{1}.parameter.condition = {'auditory'};
%      task{1}{1}.parameter.closeTo = [1 2];
%     task{1}{1}.numBlocks = 6;
%   else
%   task{1}{1}.parameter.condition = {'noOffset','posOffset','negOffset'};
%   task{1}{1}.parameter.closeTo = [1 2];
%    task{1}{1}.numBlocks = 2;
%    % task{1}{1}.numTrials = 90;
%   end
%   stimulus.nTrialTotal = length(task{1}{1}.parameter.offset) * length(task{1}{1}.parameter.condition) * length(task{1}{1}.parameter.closeTo) * task{1}{1}.numBlocks;

end

task{1}{1}.random = 1;

task{1}{1}.randVars.calculated.probeOffset = nan; % from midpoint
task{1}{1}.randVars.calculated.resp = nan;
task{1}{1}.randVars.calculated.correct = nan;
task{1}{1}.randVars.calculated.noise = nan;
task{1}{1}.randVars.calculated.rt = nan;
task{1}{1}.randVars.calculated.condNum = nan;
task{1}{1}.randVars.calculated.tr = nan;
task{1}{1}.randVars.calculated.relNum = nan;
task{1}{1}.randVars.calculated.respRight = nan;
task{1}{1}.randVars.calculated.direction = nan;

if stimulus.visualTrain || stimulus.auditoryTrain || stimulus.visualMan || stimulus.auditoryMan || stimulus.easy
  task{1}{1}.randVars.calculated.visRel = nan;
end
if stimulus.visualMan || stimulus.auditoryMan
  task{1}{1}.randVars.calculated.stimDur = nan;
  task{1}{1}.randVars.calculated.isi = nan;
  task{1}{1}.randVars.calculated.noiseContrast = nan;
  task{1}{1}.randVars.calculated.curOffset = nan;
end


% initialize the task
for phaseNum = 1:length(task{1})
  [task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end
 
% init the stimulus
myscreen = initStimulus('stimulus',myscreen);

% init the staircase
stimulus = initStair(stimulus);
% to initialize the stimulus for your experiment.
stimulus = initGaussian(stimulus,myscreen);
stimulus = initClick(stimulus,task);
% stimulus = createNoise(stimulus,myscreen);

mglWaitSecs(1);
mglClearScreen(stimulus.colors.black);
mglTextSet([],85,stimulus.colors.white);
mglTextDraw('READY',[0 0]);
mglFlush;
mglClearScreen(stimulus.colors.black);
mglTextDraw('READY',[0 0]);
mglFlush;
mglWaitSecs(1);

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
dispPsychometric(task{1}{1},stimulus);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)
global stimulus
if task.thistrial.thisseg == 1
%   if strcmp(char(task.thistrial.condition),'auditory') && task.thistrial.visRel == 0
%     disp('(spjdg) Skipping auditory low rel trial');
%     task = jumpSegment(task,inf);
%   end
  stimulus.f = 0;
  task.thistrial.r = 0;
	stimulus.fixColor = stimulus.colors.blue;
	task.thistrial.condNum = find(strcmp(char(task.thistrial.condition),{'vision','auditory','noOffset','posOffset','negOffset'}));
  if (stimulus.visualTrain && stimulus.high) || stimulus.auditoryTrain %|| (stimulus.visual && stimulus.high) || stimulus.auditory
    task.thistrial.visRel = 1;
  elseif (stimulus.visualTrain && stimulus.low) %||(stimulus.visual && stimulus.low)
    task.thistrial.visRel = 0;
  elseif stimulus.visualMan || stimulus.auditoryMan || stimulus.easy
    task.thistrial.visRel = stimulus.visrel;
  end
  if ~stimulus.visualMan && ~stimulus.easy
   if task.thistrial.visRel == 1
          task.thistrial.thisCreateNoiseCon = stimulus.noiseHighRel;
          task.thistrial.thisTranspd = stimulus.transpd.highRel;
          task.thistrial.thisTranspi = stimulus.transpi.highRel;
          task.thistrial.relNum = 1;
    else
          task.thistrial.thisCreateNoiseCon = stimulus.noiseLowRel;
          task.thistrial.thisTranspd = stimulus.transpd.lowRel;
          task.thistrial.thisTranspi = stimulus.transpi.lowRel;
          task.thistrial.relNum = 2;
    end
  

  elseif stimulus.visualMan || stimulus.easy
      task.thistrial.thisCreateNoiseCon = stimulus.noiseContrast;
      task.thistrial.relNum = task.thistrial.visRel;
      % task.thistrial.thisTranspd = stimulus.transpd.highRel;
      % task.thistrial.thisTranspi = stimulus.transpi.highRel;

  % elseif stimulus.high || stimulus.auditory || stimulus.auditoryTrain 
  %         task.thistrial.thisCreateNoiseCon = stimulus.noiseHighRel;
  %         task.thistrial.thisTranspd = stimulus.transpd.highRel;
  %         task.thistrial.thisTranspi = stimulus.transpi.highRel;
  %         task.thistrial.relNum = 1;
  % else
  %         task.thistrial.thisCreateNoiseCon = stimulus.noiseLowRel;
  %         task.thistrial.thisTranspd = stimulus.transpd.lowRel;
  %         task.thistrial.thisTranspi = stimulus.transpi.lowRel;
  %         task.thistrial.relNum = 2;
  end

    % if stimulus.auditoryTrain || stimulus.visualTrain || stimulus.visualMan || stimulus.auditoryMan
      if task.thistrial.closeTo == 1
        sign = -1;
      else
        sign = 1;
      end
      task.thistrial.direction = sign;
    % end
    if stimulus.auditoryTrain || stimulus.visualTrain
    
    	[testValue, stimulus.stair] = doStaircase('testValue', stimulus.stair);
%         testValue = task.thistrial.offset;
        task.thistrial.noise = 0;
        % if stimulus.visualTrain && stimulus.low
        %   task.thistrial.visRel = 0;
        % else
        %   task.thistrial.visRel = 1;
        % end
        task.thistrial.probeOffset = sign * testValue;
    elseif stimulus.auditoryMan || stimulus.visualMan || stimulus.easy
      task.thistrial.probeOffset = stimulus.offset * sign;
    else
      % task.thistrial.probeOffset = task.thistrial.offset;
        %% get Test Value
        [testValue, stimulus.stair{task.thistrial.relNum}{task.thistrial.condNum}] = doStaircase('testValue', stimulus.stair{task.thistrial.relNum}{task.thistrial.condNum});
 
        testValue = 10^(testValue);
        if testValue > 10 + stimulus.coordshift
          testValue = 10 + stimulus.coordshift;
        elseif testValue < -10 + stimulus.coordshift
          testValue = -10 + stimulus.coordshift;
        end
        task.thistrial.noise = 5 * randn(1); % random number from a gaussian distribution with a std of 4 deg
        while (testValue + task.thistrial.noise < -10+stimulus.coordshift) || (testValue + task.thistrial.noise > 10+stimulus.coordshift)
            task.thistrial.noise = 5 * randn(1);
        end
        task.thistrial.probeOffset = testValue + task.thistrial.noise;

    end


 	switch char(task.thistrial.condition)

		case 'vision'
      task.thistrial.xposV = [stimulus.pos1 stimulus.midPoint+task.thistrial.probeOffset];
      task.thistrial.xposA = [nan nan nan];
    case 'auditory'
      task.thistrial.xposV = [nan nan nan];
      task.thistrial.xposA = [stimulus.pos1 stimulus.midPoint+task.thistrial.probeOffset ];
    case 'noOffset'
      task.thistrial.xposV = [stimulus.pos1 stimulus.midPoint+task.thistrial.probeOffset ];
      task.thistrial.xposA = [stimulus.pos1 stimulus.midPoint+task.thistrial.probeOffset ];
    case 'posOffset' % 1st/3rd tone shifted LEFTward & disc shifted RIGHTward ===> target tone shift Rightward, target disc shift Leftward
      task.thistrial.xposV = [stimulus.pos1 stimulus.midPoint+task.thistrial.probeOffset-stimulus.delta ];
      task.thistrial.xposA = [stimulus.pos1 stimulus.midPoint+task.thistrial.probeOffset+stimulus.delta ];
    case 'negOffset' % 1st/3rd tone shifted LEFTward & disc shifted RIGHTward ===>  vice versa
      task.thistrial.xposV = [stimulus.pos1 stimulus.midPoint+task.thistrial.probeOffset+stimulus.delta ];
      task.thistrial.xposA = [stimulus.pos1 stimulus.midPoint+task.thistrial.probeOffset-stimulus.delta ];
	end
  if ~(stimulus.auditoryTrain || stimulus.visualTrain || stimulus.visualMan || stimulus.auditoryMan || stimulus.easy)
    task.thistrial.xposV(2) = task.thistrial.xposV(2) - stimulus.coordshift;
    task.thistrial.xposA(2) = task.thistrial.xposA(2) - stimulus.coordshift;
  end

	if task.thistrial.condNum ~= 1 % if not vision condition
		% stimulus = initClick(stimulus,task);
        for int = 1:2
            stimulus.sound(int) = createITD(stimulus,task.thistrial.xposA(int));
        end
  end

  if task.trialnum == 1
    stimulus = createNoise(stimulus,task,myscreen);
  end
  task.thistrial.thisNoiseTex = stimulus.noisetex;
  task.thistrial.thisnoise = stimulus.noise;
elseif any(task.thistrial.thisseg == [3 5])
    stimulus.hasPlayed = 0;
elseif task.thistrial.thisseg == 6%8
  
	stimulus.fixColor = stimulus.colors.lightgrey;


if ~(stimulus.auditoryTrain || stimulus.visualTrain || stimulus.visualMan || stimulus.auditoryMan || stimulus.easy)


 if task.trialnum < stimulus.nTrialTotal

  if task.randVars.visRel(task.trialnum+1) == 1
    task.thistrial.thisCreateNoiseCon = stimulus.noiseHighRel;
    task.thistrial.thisTranspi = stimulus.transpi.highRel;
  else
    task.thistrial.thisCreateNoiseCon = stimulus.noiseLowRel;
    task.thistrial.thisTranspi = stimulus.transpi.lowRel;
  end

  % if stimulus.high || stimulus.auditory || stimulus.auditoryTrain%task.randVars.visRel(task.trialnum+1) == 1
  %   task.thistrial.thisCreateNoiseCon = stimulus.noiseHighRel;
  %   task.thistrial.thisTranspi = stimulus.transpi.highRel;
  % else
  %   task.thistrial.thisCreateNoiseCon = stimulus.noiseLowRel;
  %   task.thistrial.thisTranspi = stimulus.transpi.lowRel;
  % end
 end


end

  stimulus = createNoise(stimulus,task,myscreen);
  task.thistrial.nextNoiseTex = stimulus.noisetex;
  task.thistrial.nextnoise = stimulus.noise;

elseif task.thistrial.thisseg == 7
  
  if stimulus.visualMan || stimulus.auditoryMan
  task.thistrial.thisCreateNoiseCon = stimulus.noiseContrast;
   stimulus = createNoise(stimulus,task,myscreen);
  task.thistrial.nextNoiseTex = stimulus.noisetex;
  task.thistrial.nextnoise = stimulus.noise;
  end

else
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)
global stimulus
mglClearScreen(stimulus.colors.black);
% mglBltTexture(stimulus.noisetex,[0 0 200 200]);
% mglFixationCross(stimulus.fixWidth,1.5,stimulus.fixColor);
% stimulus.f = stimulus.f+stimulus.screenupdate;
if any(task.thistrial.thisseg == [3 5])
  switch char(task.thistrial.condition)
	case 'vision'
    if task.thistrial.visRel == 1 %stimulus.high
		  mglBltTexture(stimulus.tex.highRel, [task.thistrial.xposV(floor(task.thistrial.thisseg/2)), 10]);
    else
      mglBltTexture(stimulus.tex.lowRel, [task.thistrial.xposV(floor(task.thistrial.thisseg/2)), 10]);
    end
    mglBltTexture(task.thistrial.thisNoiseTex,[0 0 70 70]);
    mglFillOval(0,stimulus.fixYpos,[stimulus.fixWidth,stimulus.fixWidth],stimulus.fixColor);
    if stimulus.easy
        mglTextSet([],85,stimulus.colors.white);
        if task.thistrial.thisseg == 3
            text = 'FLASH 1';
        else
            text = 'FLASH 2';
        end
        mglTextDraw(text,[task.thistrial.xposV(floor(task.thistrial.thisseg/2)) 12.5]);
    end

	case 'auditory'
        if ~stimulus.hasPlayed
		mglPlaySound(stimulus.sound(floor(task.thistrial.thisseg/2)));
        stimulus.hasPlayed = 1;
        end
    mglBltTexture(task.thistrial.thisNoiseTex,[0 0 70 70]);
    mglFillOval(0,stimulus.fixYpos,[stimulus.fixWidth,stimulus.fixWidth],stimulus.fixColor);
	otherwise
		if ~stimulus.hasPlayed
		mglPlaySound(stimulus.sound(floor(task.thistrial.thisseg/2)));
        stimulus.hasPlayed = 1;
        end
		if task.thistrial.visRel == 1
      mglBltTexture(stimulus.tex.highRel, [task.thistrial.xposV(floor(task.thistrial.thisseg/2)), 10]);
    else
      mglBltTexture(stimulus.tex.lowRel, [task.thistrial.xposV(floor(task.thistrial.thisseg/2)), 10]);
    end
    mglBltTexture(task.thistrial.thisNoiseTex,[0 0 70 70]);
    mglFillOval(0,stimulus.fixYpos,[stimulus.fixWidth,stimulus.fixWidth],stimulus.fixColor);
        
  end		

elseif task.thistrial.thisseg == 7 %9
  if ~(stimulus.visualMan || stimulus.auditoryMan ||stimulus.easy)
  if task.thistrial.r < 60
  task.thistrial.r = task.thistrial.r + 1;
  end
  [stimulus task] = updateNoise(stimulus,task,myscreen);
  mglBltTexture(task.thistrial.thisNoiseTex,[0 0 70 70]);
  mglBltTexture(task.thistrial.nextNoiseTex,[0 0 70 70]);
  mglFillOval(0,stimulus.fixYpos,[stimulus.fixWidth,stimulus.fixWidth],stimulus.fixColor);
  task.thistrial.tr = task.thistrial.r;
  mglDeleteTexture(task.thistrial.nextNoiseTex);
  mglDeleteTexture(task.thistrial.thisNoiseTex);
  else
    mglBltTexture(task.thistrial.thisNoiseTex,[0 0 70 70]);
    mglFillOval(0,stimulus.fixYpos,[stimulus.fixWidth,stimulus.fixWidth],stimulus.fixColor);

  end
elseif task.thistrial.thisseg == 6
     mglBltTexture(task.thistrial.thisNoiseTex,[0 0 70 70]);
   mglFillOval(0,stimulus.fixYpos,[stimulus.fixWidth,stimulus.fixWidth],stimulus.fixColor);
    mglTextSet([],85,stimulus.colors.white);
    mglTextDraw('RESPOND',[0 -1.5]);

else
   mglBltTexture(task.thistrial.thisNoiseTex,[0 0 70 70]);
   mglFillOval(0,stimulus.fixYpos,[stimulus.fixWidth,stimulus.fixWidth],stimulus.fixColor);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)
global stimulus
if task.thistrial.thisseg == 7 && (stimulus.visualMan || stimulus.auditoryMan)
  if isodd(task.thistrial.gotResponse)
    stimulus.fixColor = stimulus.colors.lightblue;
  else
    stimulus.fixColor = stimulus.colors.cyan;
  end
  if task.thistrial.whichButton == 7 %space bar
    task = jumpSegment(task);
  elseif task.thistrial.whichButton == 6  % f key make it faster
    if stimulus.visualMan
    task.segmin(3:5) = task.segmin(3:5) * 0.8;
    elseif stimulus.auditoryMan
      task.segmin(4) = task.segmin(4) * 0.8;
    end

    % task.segmin(3) = task.segmin(3) * 0.8;
    % task.segmin(5) = task.segmin(5) * 0.8;
    % task.segmin(4) = task.segmin(4) * 0.9;
    if task.segmin(3) < stimulus.gaussian.duration
      task.segmin(3) = stimulus.gaussian.duration;
      task.segmin(5) = stimulus.gaussian.duration;
    end
    if task.segmin(4) < 0.5;
      task.segmin(4) = 0.5;
    end
    task.segmax = task.segmin;
    disp(sprintf('(spjdg) %i: stimDur %0.4f isi %0.4f', task.trialnum,task.segmin(3),task.segmin(4)));
  elseif task.thistrial.whichButton == 5 % s key slower
    % task{1}{1}.segmin = [2.5 0.5 stimDur*120 (0.5-stimDur)*120 stimDur*120 inf inf];
    task.segmin(3:5) = task.segmin(3:5) * 1.2;
    % task.segmin(3) = task.segmin(3) * 1.2;
    % task.segmin(5) = task.segmin(5) * 1.2;
    % task.segmin(4) = task.segmin(4) * 1.1;
    task.segmax = task.segmin;
    disp(sprintf('(spjdg) %i: stimDur %0.4f isi %0.4f', task.trialnum,task.segmin(3),task.segmin(4)));

  elseif task.thistrial.whichButton == 3 % up -> increase noise contrast (harder)
    stimulus.noiseContrast = stimulus.noiseContrast + 0.025;
    if stimulus.noiseContrast > 0.8
      stimulus.noiseContrast = 0.8;
    end
    disp(sprintf('(spjdg) %i: noiseCon %0.3f', task.trialnum,stimulus.noiseContrast));


  elseif task.thistrial.whichButton == 4 % down -> decrease noise contrast (easier)
    stimulus.noiseContrast = stimulus.noiseContrast - 0.025;
    if stimulus.noiseContrast < 0;
      stimulus.noiseContrast = 0;
    end
    disp(sprintf('(spjdg) %i: noiseCon %0.3f', task.trialnum,stimulus.noiseContrast));

  elseif task.thistrial.whichButton == 9 % "-" key -> decrease offset
    stimulus.offset = stimulus.offset - 0.5;
    if stimulus.offset < 0
      stimulus.offset = 0;
    end
    disp(sprintf('(spjdg) %i: offset %0.2f', task.trialnum,stimulus.offset));
  elseif task.thistrial.whichButton == 10 % "= (+)" key -> increase offset
    stimulus.offset = stimulus.offset + 0.5;
    disp(sprintf('(spjdg) %i: offset %0.2f', task.trialnum,stimulus.offset));
  elseif task.thistrial.whichButton == 8 % RETURN --> thres set
    if stimulus.visualMan && ~isfield(stimulus,'manual')
      stimulus.manual.highVision.noiseContrast = stimulus.noiseContrast;
      stimulus.manual.highVision.gaussianContrast = stimulus.gaussian.contrast;
      stimulus.manual.highVision.gaussianDiameter = stimulus.gaussian.diameter;
      stimulus.manual.highVision.guassianDur = task.segmin(3);
      stimulus.manual.highVision.isi = task.segmin(4);
      stimulus.manual.highVision.offset = stimulus.offset;

      % No longer doing high rel cond
      stimulus.gaussian.diameter = stimulus.gaussian.diameterLowRel;
      stimulus.visrel = 0;
      stimulus.offset = 10;

      % reset segment durs
      task.segmin = [2 0.5 stimulus.gaussian.duration*75 0.5*10 stimulus.gaussian.duration*75 inf inf];
      task.segmax = task.segmin;
      disp(sprintf('(spjdg) manual mode: high vision done. starting low vision'));

    elseif stimulus.visualMan && isfield(stimulus,'manual') && isfield(stimulus.manual,'highVision')
      stimulus.manual.lowVision.noiseContrast = stimulus.noiseContrast;
      stimulus.manual.lowVision.gaussianContrast = stimulus.gaussian.contrast;
      stimulus.manual.lowVision.gaussianDiameter = stimulus.gaussian.diameter;
      stimulus.manual.lowVision.guassianDur = task.segmin(3);
      stimulus.manual.lowVision.isi = task.segmin(4);
      stimulus.manual.lowVision.offset = stimulus.offset;
      disp(sprintf('(spjdg) manual mode: low vision done'));


    elseif stimulus.auditoryMan
      stimulus.manual.audition.noiseContrast = stimulus.noiseContrast;
      stimulus.manual.audition.guassianDur = task.segmin(3);
      stimulus.manual.audition.isi = task.segmin(4);
      stimulus.manual.auditory.offset = stimulus.offset;

      disp('(spjdg) manual mode: auditory done');
    end



  end

elseif task.thistrial.thisseg ==6
% here, we just check whether this is the first time we got a response
if ~task.thistrial.gotResponse
  
	if (task.thistrial.probeOffset-stimulus.coordshift < 0 && (task.thistrial.whichButton == 1 || (stimulus.spencer==1&&task.thistrial.whichButton==11) )) || ... % closer to the first
		(task.thistrial.probeOffset-stimulus.coordshift > 0 && (task.thistrial.whichButton == 2 || (stimulus.spencer==1&&task.thistrial.whichButton==12) ))  % closer to the third 
		% correct
		task.thistrial.correct = 1;
		% if any(task.thistrial.condNum == [1 2 3])	
		% % feedback
		if stimulus.practice || stimulus.auditoryTrain || stimulus.visualTrain || stimulus.visualMan || stimulus.auditoryMan
			stimulus.fixColor = stimulus.colors.green;
        end
		disp(sprintf('(spjdg) %i:%s offset %0.3f resp %i correct', ...
            task.trialnum, char(task.thistrial.condition), task.thistrial.probeOffset-stimulus.coordshift, task.thistrial.whichButton));
        if stimulus.easy
            mglTextSet([],85,stimulus.colors.white);
        mglTextDraw('CORRECT',[0 -1.5]);
        end

	else
		% incorrect
		task.thistrial.correct = 0;
		% if any(task.thistrial.condNum == [1 2 3])
		if stimulus.practice || stimulus.auditoryTrain || stimulus.visualTrain || stimulus.visualMan || stimulus.auditoryMan
			stimulus.fixColor = stimulus.colors.red;
		end
		disp(sprintf('(spjdg) %i:%s offset %0.3f resp %i incorrect', ...
            task.trialnum, char(task.thistrial.condition), task.thistrial.probeOffset-stimulus.coordshift, task.thistrial.whichButton));
        if stimulus.easy
            mglTextSet([],85,stimulus.colors.white);
        mglTextDraw('INCORRECT',[0 -1.5]);
        end
	end

	if ~(stimulus.practice || stimulus.auditoryTrain || stimulus.visualTrain || stimulus.visualMan || stimulus.auditoryMan || stimulus.easy)
		stimulus.fixColor = stimulus.colors.cyan;

    if task.thistrial.whichButton == 1
      task.thistrial.respRight = 0;
    elseif task.thistrial.whichButton == 2
      task.thistrial.respRight = 1;
    end

    stimulus.stair{task.thistrial.relNum}{task.thistrial.condNum} = doStaircase('update', stimulus.stair{task.thistrial.relNum}{task.thistrial.condNum}, ...
      task.thistrial.respRight, log10(task.thistrial.probeOffset));

  elseif stimulus.auditoryTrain || stimulus.visualTrain
    stimulus.stair = doStaircase('update', stimulus.stair, task.thistrial.correct);
  end


    % if ~(stimulus.auditoryTrain||stimulus.visualTrain)
        
    %     if strcmp(char(task.thistrial.condition),'vision') || strcmp(char(task.thistrial.condition),'auditory')
    %         stimulus.stair{task.thistrial.relNum}{task.thistrial.condNum} = doStaircase('update', stimulus.stair{task.thistrial.relNum}{task.thistrial.condNum}, task.thistrial.correct, ...
    %         abs(task.thistrial.probeOffset));        
    %     else
    %       randomnumbers = [1 1 1 1 1 1 1 1 1 0];
    %       thisrandom = randomnumbers(randi(length(randomnumbers),1));
    %        stimulus.stair{task.thistrial.relNum}{task.thistrial.condNum} = doStaircase('update', stimulus.stair{task.thistrial.relNum}{task.thistrial.condNum}, thisrandom, ...
    %         abs(task.thistrial.probeOffset));
    %     end
    % else
    %     stimulus.stair = doStaircase('update', stimulus.stair, task.thistrial.correct);
    % end

	task.thistrial.resp = task.thistrial.whichButton;
    task.thistrial.rt = task.thistrial.reactionTime;
    task = jumpSegment(task);

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initStair(stimulus)
if stimulus.auditoryTrain || stimulus.visualTrain
	 stimulus.stair = doStaircase('init','upDown', 'nup=1','ndown=2',...
        'initialThreshold', stimulus.initialThreshold, 'initialStepsize',stimulus.initialStepsize, ...
    'minStepsize',stimulus.minStepsize,'maxStepsize',stimulus.maxStepsize,'minThreshold',stimulus.minThreshold,'maxThreshold', stimulus.maxThreshold,...
     'stepRule=pest','dispFig=1');
elseif ~(stimulus.visualMan || stimulus.auditoryMan ||stimulus.easy)

	condNames = {'vision','auditory','noOffset','posOffset','negOffset'};
  relNames = {'high','low'};
for rel = 1:2
for cond = 1:5
% 	stimulus.stair{rel}{cond} = doStaircase('init','quest', 'initialThreshold', stimulus.initOffset, 'initialThresholdSd', stimulus.initOffsetSd, ...
% 		'pThreshold', 0.75,'dispFig=1','subplotRows=5','subplotCols=2','subplotNum',2*(cond-1)+rel,'subplotName',[relNames{rel},' ', condNames{cond}]);
  stimulus.stair{rel}{cond} = doStaircase('init','quest', 'tGuess', stimulus.initOffset, 'tGuessSd', stimulus.initOffsetSd, ...
    'pThreshold', 0.75,'dispFig=1','subplotRows=5','subplotCols=2','subplotNum',2*(cond-1)+rel,'subplotName',[relNames{rel},' ', condNames{cond}]);

  % stimulus.stair{rel}{cond} = doStaircase('init','upDown', 'nup=1','ndown=2',...
  %       'initialThreshold', stimulus.initialThreshold, 'initialStepsize',stimulus.initialStepsize, ...
  %   'minStepsize',stimulus.minStepsize,'maxStepsize',stimulus.maxStepsize,'minThreshold',stimulus.minThreshold,'maxThreshold', stimulus.maxThreshold,...
  %   'dispFig=1','subplotRows=5','subplotCols=2','subplotNum',2*(cond-1)+rel,'subplotName',[relNames{rel},' ', condNames{cond}]);
end
end
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initGaussian(stimulus,myscreen)
global stimulus;

if stimulus.tenbit
% set maximum color index (for 24 bit color we have 8 bits per channel, so 255)
maxIndex = 255;

% get gamma table
if ~isfield(myscreen,'gammaTable')
  stimulus.linearizedGammaTable = mglGetGammaTable;
  disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
  disp(sprintf('(alaisburr:initGratings) No gamma table found in myscreen. Contrast'));
  disp(sprintf('         displays like this should be run with a valid calibration made by moncalib'));
  disp(sprintf('         for this monitor.'));
  disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
end
stimulus.linearizedGammaTable = myscreen.initScreenGammaTable;

% disppercent(-inf,'Creating gaussian textures');

% calculate some colors information
%  number of reserved colors
stimulus.colors.nReservedColors = size(stimulus.colors.reservedColors,1);
% number of colors possible for gratings, make sure that we 
% have an odd number
stimulus.colors.nGaussianColors = maxIndex+1-stimulus.colors.nReservedColors;
% if iseven(stimulus.colors.nGaussianColors)
%   stimulus.colors.nGaussianColors = stimulus.colors.nGaussianColors-1;
% end

% min,mid,max index of gaussian colors
stimulus.colors.minGaussianIndex = maxIndex+1 - stimulus.colors.nGaussianColors;
stimulus.colors.midGaussianIndex = stimulus.colors.minGaussianIndex + floor(stimulus.colors.nGaussianColors/2);
stimulus.colors.maxGaussianIndex = maxIndex;
% number of contrasts we can display (not including 0 contrast)
stimulus.colors.nDisplayContrasts = floor(stimulus.colors.nGaussianColors-1);

% set the reserved colors - this gives a convenient value between 0 and 1 to use the reserved colors with
for i = 1:stimulus.colors.nReservedColors
  stimulus.colors.reservedColor(i) = (i-1)/maxIndex;
end

% setGammaTableForMaxContrast(stimulus.gaussian.contrast);
maxCon = 0.8;%stimulus.noiseLowRel;%max([stimulus.gaussian.contrast, stimulus.noiseContrast]);
setGammaTableForMaxContrast(maxCon);
contrastIndex = getContrastIndex(stimulus.gaussian.contrast,1);

% make all the 1D gaussians. We compute all possible contrast values given the
% range of indexes available to us. The 1st texture is black the nth texture is full
% contrast for the current gamma setting

% high Rel gaussian
gaussian = mglMakeGaussian(stimulus.gaussian.diameterHighRel, stimulus.gaussian.diameterHighRel, stimulus.gaussian.diameterHighRel/7, stimulus.gaussian.diameterHighRel/7);
iContrast = contrastIndex-1;
thisGaussian = round(iContrast*gaussian + stimulus.colors.minGaussianIndex);
% thisGaussian = zeros(size(gaussian,1), size(gaussian,2),4);
% gauss10 = round(iContrast*gaussian + stimulus.colors.minGaussianIndex);
% for i = 1:3
%     thisGaussian(:,:,i) = gauss10;
% end
% thisGaussian(:,:,4) =  round(gaussian * 255 * stimulus.gaussian.contrast);
stimulus.tex.highRel = mglCreateTexture(thisGaussian);

% low rel gaussian
gaussian = mglMakeGaussian(stimulus.gaussian.diameterLowRel, stimulus.gaussian.diameterLowRel, stimulus.gaussian.diameterLowRel/7, stimulus.gaussian.diameterLowRel/7);
thisGaussian = round(iContrast*gaussian + stimulus.colors.minGaussianIndex);
stimulus.tex.lowRel = mglCreateTexture(thisGaussian);


% get the color value for black (i.e. the number between 0 and 1 that corresponds to the minGaussianIndex)
stimulus.colors.black = stimulus.colors.minGaussianIndex/maxIndex;
stimulus.colors.grey =  stimulus.colors.midGaussianIndex/maxIndex;
% get the color values (i.e. reserved color)
stimulus.colors.white = stimulus.colors.reservedColor(1);
stimulus.colors.darkgrey = stimulus.colors.reservedColor(2);
stimulus.colors.green = stimulus.colors.reservedColor(3);
stimulus.colors.red = stimulus.colors.reservedColor(4);
stimulus.colors.cyan = stimulus.colors.reservedColor(5);
stimulus.colors.blue = stimulus.colors.reservedColor(6);
stimulus.colors.lightblue = stimulus.colors.reservedColor(7);
stimulus.colors.lightgrey = stimulus.colors.reservedColor(8);


% % compute the guassian
% gauss = mglMakeGaussian(stimulus.width,stimulus.width, stimulus.width/8,stimulus.width/8);

% gaussian = zeros(size(gauss,1), size(gauss,2), 4);
% for i = 1:3
%     gaussian(:,:,i) = 255*ones(size(gauss,1), size(gauss,2));
% end
%     gaussian(:,:,4) = 255*gauss*stimulus.contrast;
    
% %create texture
% stimulus.tex = mglCreateTexture(gaussian);

else
  %high rel gauss
	gauss = mglMakeGaussian(stimulus.gaussian.diameterHighRel, stimulus.gaussian.diameterHighRel, stimulus.gaussian.diameterHighRel/7, stimulus.gaussian.diameterHighRel/7);
	gaussian = zeros(size(gauss,1), size(gauss,2), 4);
	for i = 1:3
    	gaussian(:,:,i) = 255*ones(size(gauss,1), size(gauss,2));
	end
    gaussian(:,:,4) = round(255*gauss*stimulus.gaussian.contrast);
    stimulus.tex.highRel = mglCreateTexture(gaussian);
  %low rel gauss
  gauss = mglMakeGaussian(stimulus.gaussian.diameterLowRel, stimulus.gaussian.diameterLowRel, stimulus.gaussian.diameterLowRel/7, stimulus.gaussian.diameterLowRel/7);
  gaussian = zeros(size(gauss,1), size(gauss,2), 4);
  for i = 1:3
      gaussian(:,:,i) = 255*ones(size(gauss,1), size(gauss,2));
  end
    gaussian(:,:,4) = round(255*gauss*stimulus.gaussian.contrast);
    stimulus.tex.lowRel = mglCreateTexture(gaussian);

    % get the color value for black (i.e. the number between 0 and 1 that corresponds to the minGaussianIndex)
stimulus.colors.black = 0;
% get the color values (i.e. reserved color)
stimulus.colors.white = 1;
stimulus.colors.darkgrey = 0.2;
stimulus.colors.green = [0 1 0];
stimulus.colors.red = [1 0 0];
stimulus.colors.cyan = [0 1 1];
stimulus.colors.blue = [0 0 1];
stimulus.colors.lightblue =[0 0.75 0.75];
stimulus.colors.lightgrey =[0.75 0.75 0.75];

end   
	
%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getContrastIndex    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function contrastIndex = getContrastIndex(desiredContrast,verbose)

if nargin < 2,verbose = 0;end

global stimulus;
if desiredContrast < 0, desiredContrast = 0;end

% now find closest matching contrast we can display with this gamma table
contrastIndex = min(round(stimulus.colors.nDisplayContrasts*desiredContrast/stimulus.currentMaxContrast),stimulus.colors.nDisplayContrasts);

% display the desired and actual contrast values if verbose is set
if verbose
  actualContrast = stimulus.currentMaxContrast*(contrastIndex/stimulus.colors.nDisplayContrasts);
  disp(sprintf('(getContrastIndex) Desired contrast: %0.4f Actual contrast: %0.4f Difference: %0.4f',desiredContrast,actualContrast,desiredContrast-actualContrast));
end

% out of range check
if round(stimulus.colors.nDisplayContrasts*desiredContrast/stimulus.currentMaxContrast)>stimulus.colors.nDisplayContrasts
 disp(sprintf('(getContrastIndex) Desired contrast (%0.9f) out of range max contrast : %0.9f',desiredContrast,stimulus.currentMaxContrast));
 keyboard
end

% 1 based indexes (0th index is gray, nDisplayContrasts+1 is full contrast)
contrastIndex = contrastIndex+1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sets the gamma table so that we can have
% finest possible control over the stimulus contrast.
%
% stimulus.reservedColors should be set to the reserved colors (for cue colors, etc).
% maxContrast is the maximum contrast you want to be able to display.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setGammaTableForMaxContrast(maxContrast)

global stimulus;
% if you just want to show gray, that's ok, but to make the
% code work properly we act as if you want to display a range of contrasts
if maxContrast <= 0,maxContrast = 0.01;end

% set the reserved colors
gammaTable(1:size(stimulus.colors.reservedColors,1),1:size(stimulus.colors.reservedColors,2))=stimulus.colors.reservedColors;

% set the gamma table
if maxContrast > 0
  % create the rest of the gamma table
%   cmax = 0.5+maxContrast/2;cmin = 0.5-maxContrast/2;
  cmin = 0;
  cmax = maxContrast;
  luminanceVals = cmin:((cmax-cmin)/(stimulus.colors.nGaussianColors-1)):cmax;

  % replace NaN in gamma tables with zero
  stimulus.linearizedGammaTable.redTable(isnan(stimulus.linearizedGammaTable.redTable)) = 0;
  stimulus.linearizedGammaTable.greenTable(isnan(stimulus.linearizedGammaTable.greenTable)) = 0;
  stimulus.linearizedGammaTable.blueTable(isnan(stimulus.linearizedGammaTable.blueTable)) = 0;

  % now get the linearized range
  redLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,luminanceVals,'linear');
  greenLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,luminanceVals,'linear');
  blueLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,luminanceVals,'linear');
  
  % add these values to the table
  gammaTable((stimulus.colors.minGaussianIndex:stimulus.colors.maxGaussianIndex)+1,:)=[redLinearized;greenLinearized;blueLinearized]';
else
  % if we are asked for 0 contrast then simply set all the values to BLACK
  gammaTable((stimulus.colors.minGaussianIndex:stimulus.colors.maxGaussianIndex)+1,1)=interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,0,'linear');
  gammaTable((stimulus.colors.minGaussianIndex:stimulus.colors.maxGaussianIndex)+1,2)=interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,0,'linear');
  gammaTable((stimulus.colors.minGaussianIndex:stimulus.colors.maxGaussianIndex)+1,3)=interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,0,'linear');
end

% set the gamma table
mglSetGammaTable(gammaTable);

% remember what the current maximum contrast is that we can display
stimulus.currentMaxContrast = maxContrast;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = createNoise(stimulus,task,myscreen)
xDeg2pix = mglGetParam('xDeviceToPixels');
yDeg2pix = mglGetParam('yDeviceToPixels');
width = 2.5; height = 2.5;
% get size in pixels
widthPixels = round(width*xDeg2pix);
heightPixels = round(height*yDeg2pix);
widthPixels = widthPixels + mod(widthPixels+1,2);
heightPixels = heightPixels + mod(heightPixels+1,2);

if ~stimulus.tenbit
  stimulus.noise = zeros(widthPixels,heightPixels,4);  
  stimulus.noise(:,:,4) =  round(ones(widthPixels,heightPixels) * 255 * task.thistrial.thisCreateNoiseCon);
  n = round(rand(widthPixels,heightPixels) * 255);
  for i = 1:3
    stimulus.noise(:,:,i) = n;
  end
  stimulus.noisetex = mglCreateTexture(stimulus.noise,[],0,{'GL_TEXTURE_MAG_FILTER','GL_NEAREST'});
 
% stimulus.noisetex = mglCreateTexture(round(stimulus.noiseContrast*rand(widthPixels,heightPixels)*255));

else
    
contrastIndex = getContrastIndex(task.thistrial.thisCreateNoiseCon,0);
iContrast = contrastIndex-1;
stimulus.noise = zeros(widthPixels,heightPixels,4);  
stimulus.noise(:,:,4) =  round(ones(widthPixels,heightPixels) * 255 * task.thistrial.thisCreateNoiseCon);
n = round(rand(widthPixels,heightPixels) * (255-stimulus.colors.minGaussianIndex) + stimulus.colors.minGaussianIndex);
  for i = 1:3
    stimulus.noise(:,:,i) = n;
  end
  stimulus.noisetex = mglCreateTexture(stimulus.noise,[],0,{'GL_TEXTURE_MAG_FILTER','GL_NEAREST'});

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stimulus task] = updateNoise(stimulus,task,myscreen)
task.thistrial.thisnoise(:,:,4) = round(task.thistrial.thisTranspd(task.thistrial.r));
task.thistrial.nextnoise(:,:,4) = round(task.thistrial.thisTranspi(task.thistrial.r));
task.thistrial.thisNoiseTex = mglCreateTexture(task.thistrial.thisnoise,[],0,{'GL_TEXTURE_MAG_FILTER','GL_NEAREST'});
task.thistrial.nextNoiseTex = mglCreateTexture(task.thistrial.nextnoise,[],0,{'GL_TEXTURE_MAG_FILTER','GL_NEAREST'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sound
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initClick(stimulus,myscreen)
% sampling frequency (samples per sec)
fs = stimulus.tone.samplesPerSecond-1;
t = 0:1/fs:stimulus.tone.duration;

wav = 0.5 * randn(1,length(t));

fc = 2000; % cutoff frequency
% 5th order Butterworth filter
% [b,a] = butter(5, fc/(stimulus.tone.samplesPerSecond/2));
a = [1	-4.07876493416512	6.72527084144657	-5.59474636818042	2.34559680959441	-0.396133028715511];
b = [3.82287493728255e-05	0.000191143746864127	0.000382287493728255	0.000382287493728255	0.000191143746864127	3.82287493728255e-05];
wav = randn(1,length(t));
wavFiltered = filter(b,a,wav);
stimulus.wav = wavFiltered;

% stimulus.wav = wav;

% stimulus.thisHz = stimulus.tone.hz + (stimulus.tone.hz/10)*randn;

% wav = sin(2*pi*stimulus.thisHz*t);
% if isodd(length(t))
% 	mu = t((1+length(t))/2);
% else
% 	mu = t(round((1+length(t))/2));
% end
% sigma = stimulus.tone.duration;
% envelope = normpdf(t,mu,sigma);
% %normalize
% envelope = envelope / max(envelope);
% waveform = wav .* envelope;
% stimulus.wav = waveform;

% frequency of signal in hz
% stimulus.wav = sin(2*pi*stimulus.tone.hz*t);
stimulus.sound(1) = mglInstallSound(stimulus.wav, stimulus.tone.samplesPerSecond);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to update the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = createITD(stimulus,theta)
% theta: angle to sound source (-90 to left and +90 to right)

len = length(stimulus.wav);
% load impulse response
elev = 0;
if mod(theta, 5) == 0
    if theta >= 0
        temp = readhrtf(elev,theta,'L');
        IR_L = temp(1,:);
        IR_R = temp(2,:);
    else
        temp = readhrtf(elev,-theta,'L');
        IR_L = temp(2,:);
        IR_R = temp(1,:);
    end
    % fft
    yl = fft(IR_L, len);
    yr = fft(IR_R, len);
else
    if theta >= 0
        temp{1} = readhrtf(elev, floor(theta/5) * 5, 'L');
        temp{2} = readhrtf(elev, ceil(theta/5) * 5, 'L');
        for i = 1:2
            IR_L{i} = temp{i}(1,:);
            IR_R{i} = temp{i}(2,:);
        end
    else
        theta = -theta;
        temp{1} = readhrtf(elev, floor(theta/5) * 5, 'L');
        temp{2} = readhrtf(elev, ceil(theta/5) * 5, 'L');
        for i = 1:2
            IR_L{i} = temp{i}(2,:);
            IR_R{i} = temp{i}(1,:);
        end
    end
    for i = 1:2
        YL{i} = fft(IR_L{i},len);
        YR{i} = fft(IR_R{i},len);
    end
    
    yl = ((ceil(theta/5)*5 - theta)/5 * YL{1}) + ((theta - floor(theta/5)*5)/5 * YL{2});
    yr = ((ceil(theta/5)*5 - theta)/5 * YR{1}) + ((theta - floor(theta/5)*5)/5 * YR{2});
end
yw = fft(stimulus.wav,len);
if size(yw,1) ~= size(yl,1);
    yl = yl';
    yr = yr';
end

leftFunc = yw .* yl;
rightFunc = yw .* yr;

leftwav = ifft(leftFunc, len);
rightwav = ifft(rightFunc,len);
clear waveform
if ~stimulus.spencer
if size(leftwav,1) > size(leftwav,2)
    waveform(1,:) = leftwav';
    waveform(2,:) = rightwav';
else
    waveform(1,:) = leftwav;
    waveform(2,:) = rightwav;
end
else
  if size(leftwav,1) > size(leftwav,2)
    waveform(1,:) = rightwav';
    waveform(2,:) = leftwav';
else
    waveform(1,:) = rightwav;
    waveform(2,:) = leftwav;
end
  
end
s = mglInstallSound(waveform, stimulus.tone.samplesPerSecond);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display psychometric functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispPsychometric(task,stimulus)
uni = figure; 
if ~(stimulus.auditoryTrain || stimulus.visualTrain)
bi = figure;
condNames = {'vision','auditory','noOffset','posOffset','negOffset'};
nCon = 5;
else
	nCon = 2;
	
end

% condition = {};
% for b = 1:length(task.block)
%   condition = [condition, task.block(b).parameter.condition];
% end
condition = task.randVars.condNum;

for cond = 1:nCon
  thisTrials = find(condition == cond);%find(strcmp(condition, condNames{cond}));
  thisOffset = task.randVars.probeOffset(thisTrials);
  thisResp = task.randVars.resp(thisTrials);
  isThird = (thisResp == 2);

	binCenter = -14.5:1:14.5; space = max(diff(binCenter));
  for b = 1:length(binCenter)
    switch b
      case 1
        binnedTrial{b} = find(thisOffset<= min(binCenter)+space/2); % -200 ~ -150 ms
      case length(binCenter)
        binnedTrial{b} = find(thisOffset > max(binCenter)-space/2);
      otherwise
        binnedTrial{b} = find((thisOffset > min(binCenter) + space*(b-2)) & (thisOffset <= min(binCenter) + space*(b-1)));
    end
  end
  % binCenter = -12.5:1:12.5;%-0.175:.05:0.175;
  nTotalBin = cellfun(@(x) length(x), binnedTrial);
  binnedTrial(nTotalBin==0) = [];
  binCenter(nTotalBin==0) = [];
  nTotalBin(nTotalBin==0) = [];

  nThirdBin = cellfun(@(x) sum(isThird(x)), binnedTrial);
  pThirdBin = nThirdBin./nTotalBin;

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
      myerrorbar(binCenter, pThirdBin,'Symbol',getsymbol(cond),'Color',getcolor(cond), 'MarkerSize',6);
      % myerrorbar(binCenter+(rand-0.5)*0.002, pLateBin,'yError',error, 'Symbol',getsymbol(cond),'Color',getcolor(cond), 'MarkerSize',6);
      % plot(xs, fitVal, 'Color',getcolor(cond),'LineWidth',0.8);
    otherwise
    	figure(bi);
     hold on;
      myerrorbar(binCenter, pThirdBin,'Symbol',getsymbol(cond),'Color',getcolor(cond), 'MarkerSize',6);
      % myerrorbar(binCenter+(rand-0.5)*0.002, pLateBin,'yError',error, 'Symbol',getsymbol(cond),'Color',getcolor(cond), 'MarkerSize',6);
      % plot(xs, fitVal, 'Color',getcolor(cond),'LineWidth',0.8);
  end

end

if ~(stimulus.auditoryTrain || stimulus.visualTrain)
figure(bi); mylegend({'No offset','Pos offset','Neg offset'}, {{getcolor(3)},{ getcolor(4)},{ getcolor(5)}});
end
figure(uni); mylegend({'vision','auditory'},{{getcolor(1)},{getcolor(2)}});

% 
% function dispPerformance(task)
% hold on;
% if task.thistrial.correct
% plot(task.trialnum, task.thistrial.correct, 'go', 'markerFaceColor',[0 1 0]);
% else
%     plot(task.trialnum, task.thistrial.correct, 'ro', 'markerFaceColor',[1 0 0]);
% end
% drawnow;
% if any(task.trialnum == [10 20 30])
%     disp(sprintf('(spatialBisection) Percent Correct of Last 10 Trials: %0.2f%', sum(task.randVars.correct(task.trialnum-9:task.trialnum))/task.trialnum*100));
% end