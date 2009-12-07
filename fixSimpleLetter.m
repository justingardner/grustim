% fixStairInitTask.m
%
%        $Id: fixSimpleLetter.m,v 1.2 2007/01/17 21:11:10 ds Exp $
%      usage: [fixTask myscreen] = fixStairInitTask(myscreen)
%         by: justin gardner
%       date: 09/07/06
%  copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%    purpose: 
%
function [task myscreen] = fixSimpleLetter(myscreen)

global debugFlag

% check arguments
if ~any(nargin == [1])
  help fixSanghunInitTask
  return
end

% create the stimulus for the experiment, use defaults if they are
% not already set
myscreen = initStimulus('fixStimulus',myscreen);
global fixStimulus;

% this means that the repeat will happen in the
% uniform distribution over the last repeatLoc places

if ~isfield(fixStimulus,'fontSize'), fixStimulus.fontSize= 36; end
if ~isfield(fixStimulus,'fontName'), fixStimulus.fontName = 'Helvetica'; end
if ~isfield(fixStimulus,'colors'), fixStimulus.colors = {[1 1 1]*0.7; [1 1 1]*0.5; [1 1 1]*0.25; [1 1 1]*0.1}; end
if ~isfield(fixStimulus,'letters'), fixStimulus.letters = ['a':'d' '0':'3']; end
if ~isfield(fixStimulus,'repeatLoc'), fixStimulus.repeatLoc = 4; end
if ~isfield(fixStimulus,'fixColor'), fixStimulus.fix.color = [1 1 1]; end
if ~isfield(fixStimulus,'colorClockLen'), fixStimulus.fix.colorClockLen = 30; end
% timing
if ~isfield(fixStimulus,'nItems'), fixStimulus.nItems = 16; end
if ~isfield(fixStimulus,'stimTime'), fixStimulus.stimTime = 0.25; end
if ~isfield(fixStimulus,'responseTime'), fixStimulus.responseTime = 0.75; end

fixStimulus.colors = {[1 1 1]*0.7; [1 1 1]*0.5; [1 1 1]*0.25; [1 1 1]*0.1};

% set some default options for the text
mglTextSet('',0,0,0,0,0,0,0,0,0);

% create all the stimuli
fixStimulus.n = 0;
% disppercent(-inf,'Rendering text');
disp('Rendering text')
for colorNum = 1:length(fixStimulus.colors)
  % disppercent((colorNum-1)/length(fixStimulus.colors));
  for letterNum = 1:length(fixStimulus.letters)
    mglTextSet(fixStimulus.fontName,fixStimulus.fontSize,fixStimulus.colors{colorNum});
    fixStimulus.n = fixStimulus.n+1;
    fixStimulus.textures(fixStimulus.n) = mglText(fixStimulus.letters(letterNum));
    mydisp(sprintf('.%i',fixStimulus.n));
    % and keep information about what it is
    fixStimulus.thisLetter(fixStimulus.n) = fixStimulus.letters(letterNum);
    fixStimulus.thisColor(fixStimulus.n) = colorNum;
    fixStimulus.isDigit(fixStimulus.n) = any(fixStimulus.letters(letterNum)=='0':'9');
  end
end
mydisp(sprintf('.\nDone\n'));
% disppercent(inf);

% waitForKeypress;

nItems = fixStimulus.n; % number of letters/numbers to display

% create a fixation task
task{1}{1}.seglen = [ones(1,nItems)*fixStimulus.stimTime [1].*fixStimulus.responseTime];
task{1}{1}.getResponse = [zeros(1, nItems) 0]; % << reset to get response in last interval!
task{1}{1}.random = 0;
task{1}{1}.traceNum = 1;
task{1}{1}.writeSegmentsTrace = 0;
task{1}{1}.waitForBacktick = 0;

% init the task
task{1}{1} = initTask(task{1}{1},myscreen,...
		   @fixStartSegmentCallback,...
		   @fixDrawStimulusCallback,...
		   @fixTrialResponseCallback);

% keep the correct and incorrect counts
task{1}{1}.correct = 0;
task{1}{1}.incorrect = 0;



if debugFlag
  mglClearScreen;
  mglBltTexture([fixStimulus.textures(:)], [linspace(-10,10,fixStimulus.n)' zeros(fixStimulus.n,1)]);
  mglFlush;
  pause
end


% keep info about which textures are digits and which ones are letters
fixStimulus.digitNums = find(fixStimulus.isDigit);
fixStimulus.numDigits = length(fixStimulus.digitNums);
fixStimulus.letterNums = find(~fixStimulus.isDigit);
fixStimulus.numLetters = length(fixStimulus.letterNums);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fixation task callback: TrialCallback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = fixTrialStartCallback(task,myscreen)

global fixStimulus

% set trace
myscreen = writetrace(0,task.traceNum,myscreen);
% set color clock 
fixStimulus.fix.colorClock = 0;
% set fixation color to white
fixStimulus.fix.color = [1 1 1];
% create a sequence that will be displayed
% first figure out how many letters we will have
numSegs = length(task.seglen)-1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fixation task callback: fixStartSegmentCallback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = fixStartSegmentCallback(task,myscreen)

global fixStimulus;

% first segment, init the trial
if (task.thistrial.thisseg < length(task.seglen))
  %seqNum = fixStimulus.thisSequence(task.thistrial.thisseg);
  %mydisp(sprintf('%c %i, ',fixStimulus.thisLetter(seqNum),fixStimulus.thisColor(seqNum)));
  fixStimulus.iLetter = task.thistrial.thisseg;
  
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fixation task callback: fixTrialStimulusCallback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = fixDrawStimulusCallback(task,myscreen)

global fixStimulus;

% draw letters one after the other.

mglStencilSelect(myscreen.stencil.fixation);
% mglClearScreen;
mglGluDisk(0,0,myscreen.fix.diskSize,myscreen.background,60);
if (task.thistrial.thisseg >= 1) && (task.thistrial.thisseg < length(task.seglen))
  % draw the correct text
  mglBltTexture(fixStimulus.textures(fixStimulus.iLetter),myscreen.fix.loc);
else % last segment... this is the response interval
     % along x=0
     mglLines2(-myscreen.fix.size(1)+myscreen.fix.loc(1),myscreen.fix.loc(2),...
	       +myscreen.fix.size(1)+myscreen.fix.loc(1),myscreen.fix.loc(2),...
	       myscreen.fix.linewidth,fixStimulus.fix.color);
     % along y=0
     mglLines2(myscreen.fix.loc(1),-myscreen.fix.size(2)+myscreen.fix.loc(2),...
	       myscreen.fix.loc(1),+myscreen.fix.size(2)+myscreen.fix.loc(2),...
	       myscreen.fix.linewidth,fixStimulus.fix.color);
end
mglStencilSelect(0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fixation task callback: fixTrialResponseCallback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = fixTrialResponseCallback(task,myscreen)

global fixStimulus;

% update number of responses
task.thistrial.gotResponse = task.thistrial.gotResponse+1;

% only the first response
if (task.thistrial.gotResponse == 1)
  mydisp(sprintf('buttonstate(1): %i, buttonstate(2): %i\n',task.thistrial.buttonState(1),task.thistrial.buttonState(1)));
end



