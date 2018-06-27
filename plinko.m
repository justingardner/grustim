% plinko.m
%
%      usage: plinko(varargin)
%         by: minyoung lee
%       date: 06/2018
%    purpose: run plinko task (from Gerstenberg & Siegel & Tenenbaum, 2018) with mgl
%			  To run the audio-visual prediction task (Expt 2a)
%			  plinko('prediction=1')
%			  
%			  To run the audio-visual inference task (Expt 2b)
%			  plinko('inference=1')
%
%			  To run the audio-visual inference with occlusion task (Expt 3)
%			  plinko('occlusion=1')

function myscreen = plinko(varargin)
clear global
% check arguments
prediction = 0; inference = 0; occlusion = 0;
getArgs(varargin,{'prediction=0','inference=0','occlusion=0'},'verbose=1');

if ~(prediction || inference || occlusion)
	display('(plinko) task type not defined. running prediction task...');
	prediction = 1;
end

% world = [2,6,9,10,12,13,14,15,16,20,26,27,32,38,42,43,45,47,48,50,56,59,62,63,70,74,76,79,80,82,84,86,90,91,92,94,99,104,105,109,111,112,113,115,118];
% hole1 = [2,6,9,10,12,13,14,16,20,26,27,32,38,42,43,45,47,48,50,56,59,62,70,74,80,82,84,86,90,91,92,94,99,104,105,109,111,112,113,115,118];
world = [2,6,9,10,12,13,14,16,20,26,27,32,38,42,45,47,48,50,56,59,62,70,74,80,82,84,86,90,91,92,94,99,104,105,109,111,112,113,115,118];
% hole3 = [2,6,9,10,12,13,14,16,20,26,27,32,38,42,45,47,48,50,56,59,62,63,70,74,76,80,82,84,86,90,91,92,94,99,104,105,109,111,112,113,115,118];
stimulus.prediction = prediction;
stimulus.inference = inference;
stimulus.occlusion = occlusion;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the screen
myscreen = initScreen;
myscreen.keyboard.nums = [19,20,21,50]; % 1,2,3,space

% set up task
task{1}{1}.waitForBacktick = 1;
task{1}{1}.segmin = [2 1 4 4 3 inf];
task{1}{1}.segmax = [2 1 4 4 3 inf];
task{1}{1}.getResponse = [0 0 0 0 0 1];
task{1}{1}.random = 1;
task{1}{1}.numBlocks = 1;

task{1}{1}.parameter.world = world;
if stimulus.prediction
	task{1}{1}.parameter.hole = [1,2,3];
end

task{1}{1}.randVars.calculated.resp = nan;
task{1}{1}.randVars.calculated.rt = nan;
task{1}{1}.randVars.calculated.clickPosPix = nan;
task{1}{1}.randVars.calculated.clickPosDeg = nan;
task{1}{1}.randVars.calculated.whichHole = nan;
task{1}{1}.randVars.calculated.thisWorld = nan;


% initialize the task
for phaseNum = 1:length(task{1})
  [task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end

% init the stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);

stimulus = initBox(stimulus);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% myscreen = eyeCalibDisp(myscreen);

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

global stimulus;

if task.thistrial.thisseg == 1
	initWorld(task);
	initSound(task);
	if stimulus.prediction
		mglDisplayCursor(0);
	end

	task.thistiral.thisWorld = task.thistrial.world;
	if stimulus.prediction
		task.thistrial.whichHole = task.thistrial.hole;
	end
elseif any(task.thistrial.thisseg == [3 4])
	stimulus.tic = mglGetSecs;
elseif (task.thistrial.thisseg == 6)
	% response
	if stimulus.prediction
		mglDisplayCursor(1);
		stimulus.xposout = nan;
		stimulus.rawx = nan;
	end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus

mglClearScreen;
if (task.thistrial.thisseg == 1)
  mglFixationCross(1,1,[1 1 1]);
elseif (task.thistrial.thisseg == 2)
	% show occluded box
	if stimulus.prediction
		mglBltTexture(stimulus.tex.hole(task.thistrial.hole), [0 0 stimulus.widthDeg stimulus.heightDeg],0,0,180);
	elseif stimulus.inference || stimulus.occlusion
		mglBltTexture(stimulus.tex.cover, [0 0 stimulus.widthDeg stimulus.heightDeg],0,0,180);
	end
elseif any(task.thistrial.thisseg == [3 4])
	% show occluded box
	if stimulus.prediction
		mglBltTexture(stimulus.tex.hole(task.thistrial.hole), [0 0 stimulus.widthDeg stimulus.heightDeg],0,0,180);
	elseif stimulus.inference || stimulus.occlusion
		mglBltTexture(stimulus.tex.cover, [0 0 stimulus.widthDeg stimulus.heightDeg],0,0,180);
	end
	% play sound
	mglPlaySound(stimulus.sound);
	stimulus.toc = mglGetSecs(stimulus.tic);
	if stimulus.toc > (stimulus.soundDur + 2);
		task = jumpSegment(task);
	end
elseif (task.thistrial.thisseg == 5)
	% show the obstacles
	mglBltTexture(stimulus.tex.thisWorld, [0 0 stimulus.widthDeg stimulus.heightDeg],0,0,180);
elseif (task.thistrial.thisseg == 6)
	% show the obstacles
	mglBltTexture(stimulus.tex.thisWorld, [0 0 stimulus.widthDeg stimulus.heightDeg],0,0,180);
	if stimulus.prediction
		mglTextSet([],40,[0 0 0]);
    	mglTextDraw('Where did the ball land?',[0 0]);
		% get mouse response
		click = mglGetMouseEvent();
    	if exist('click','var') && click.clickState && (click.x >= 0 && click.x <= myscreen.screenWidth)
    	    stimulus.xposout = (click.x - myscreen.screenWidth/2) / stimulus.xDeg2pix;
    	end
    	if ~isnan(stimulus.xposout)
    		mglFillOval(stimulus.xposout, -stimulus.heightDeg/2 + 1.5, [1 1],[1 0 0]);
    	end
    else
    	mglTextSet([],40,[0 0 0]);
    	mglTextDraw('In which hole was the ball dropped?',[0 0]);
	end

end  


%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)

global stimulus

if stimulus.inference || stimulus.occlusion
% here, we just check whether this is the first time we got a response
% this trial and display what the subject's response was and the reaction time
	if task.thistrial.gotResponse < 1
		task.thistrial.resp = task.thistrial.whichButton;
    	task.thistrial.rt = task.thistrial.reactionTime;
    	task = jumpSegment(task);
  		disp(sprintf('Subject response: %i Reaction time: %0.2fs',task.thistrial.whichButton,task.thistrial.reactionTime));
	end
elseif stimulus.prediction
	if task.thistrial.whichButton == 4
        task.thistrial.resp = stimulus.xposout;
        task.thistrial.rt = task.thistrial.reactionTime;
        task.thistrial.clickPosPix = stimulus.xposout * stimulus.xDeg2pix;
		task.thistrial.clickPosDeg = stimulus.xposout;
        task = jumpSegment(task);
        disp(sprintf('Subject response: %i Reaction time: %0.2fs',task.thistrial.resp,task.thistrial.reactionTime));

    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load covered boxes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initBox(stimulus)
if stimulus.prediction
	for h = 1:3
		clear img
		img = double(imread(sprintf('~/proj/plinko_fmri/figures/stimuli/normal/misc/cover_hole%i.jpg',h)));
		stimulus.tex.hole(h) = mglCreateTexture(img);
	end
elseif stimulus.inference
	img = double(imread('~/proj/plinko_fmri/figures/stimuli/normal/misc/cover.jpg'));
	stimulus.tex.cover = mglCreateTexture(img);
elseif stimulus.occlusion
	img = double(imread('~/proj/plinko_fmri/figures/stimuli/normal/misc/cover_split.jpg'));
	stimulus.tex.cover = mglCreateTexture(img);
end

stimulus.xDeg2pix = mglGetParam('xDeviceToPixels');
stimulus.yDeg2pix = mglGetParam('yDeviceToPixels');
% get size in degrees
% stimulus.height = size(img,1) / stimulus.xDeg2pix;
% stimulus.width = size(img,2) / stimulus.xDeg2pix;
stimulus.heightDeg = 20;
stimulus.widthDeg = 24;
stimulus.heightPix = stimulus.heightDeg * stimulus.yDeg2pix;
stimulus.widthPix = stimulus.widthDeg * stimulus.xDeg2pix;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load world with obstacles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initWorld(task)
global stimulus
if stimulus.prediction
	world = sprintf('~/proj/plinko_fmri/figures/stimuli/normal/initial/hole%i_world%i.jpg',task.thistrial.hole,task.thistrial.world);
elseif stimulus.inference
	world = sprintf('~/proj/plinko_fmri/figures/stimuli/normal/final/final_world%i.jpg',task.thistrial.world);
elseif stimulus.occlusion
	world = sprintf('~/proj/plinko_fmri/figures/stimuli/normal/occluder/occluder_world%i.jpg',task.thistrial.world);
end
img = double(imread(world));
stimulus.tex.thisWorld = mglCreateTexture(img);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load wav
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initSound(task)
global stimulus
if stimulus.inference || stimulus.occlusion
	wavfile = sprintf('~/proj/plinko_fmri/audio/inference/world%i.wav', task.thistrial.world);
elseif stimulus.prediction
	wavfile = sprintf('~/proj/plinko_fmri/audio/prediction/world%ihole%i.wav', task.thistrial.world,task.thistrial.hole);
end
[y,fs] = audioread(wavfile);
stimulus.sound = mglInstallSound(y',fs);
stimulus.soundDur = length(y) / fs;


