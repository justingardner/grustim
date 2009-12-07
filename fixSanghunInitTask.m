% fixStairInitTask.m
%
%        $Id: fixSanghunInitTask.m,v 1.8 2007/06/26 13:39:00 ds Exp $
%      usage: [fixTask myscreen] = fixStairInitTask(myscreen)
%         by: justin gardner
%       date: 09/07/06
%  copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%    purpose: 
%
function [task myscreen] = fixSanghunInitTask(myscreen, task)

global debugFlag

% check arguments
if ~any(nargin == [1 2])
  help fixSanghunInitTask
  return
end

if nargin < 2
  task = {cell(1)};
else % check that task is of size {1}{1}
  if prod(size(task)) ~= 1, error('(UHOH) input "task" has wrong size'), end 
end

% create the stimulus for the experiment, use defaults if they are
% not already set
global fixStimulus;
myscreen = initStimulus('fixStimulus',myscreen);

% this means that the repeat will happen in the
% uniform distribution over the last repeatLoc places

% ------------------------------------------------------------
% set defaults
% ------------------------------------------------------------
fixStimulus.fontSize= 36;
fixStimulus.fontName = 'Helvetica'; 
fixStimulus.colors = {[1 0 0],[1 0.5 0]}; 
fixStimulus.letters = ['a':'b' '0':'1']; 
fixStimulus.repeatLoc = 2; 
fixStimulus.fixColor = [1 1 1]; 
fixStimulus.fixColorClockLen = 0; % 30 
% timing
fixStimulus.nItems = 5; 
fixStimulus.stimTime = 0.5; 
fixStimulus.responseTime = 0.75; 
fixStimulus.itiMinTime = 0;
fixStimulus.itiMaxTime = 0;

% ------------------------------------------------------------
% then reset fields according to what was passed in
% ------------------------------------------------------------
allFields = {'fontSize', 'fontName', 'colors', 'letters', 'repeatLoc','fixColor','fixColorClockLen','nItems','stimTime','responseTime','itiMinTime','itiMaxTime'};
for iField = allFields
  % check whether info was passed in in 'task' and (re)assign to fixStimulus
   if isfield(task{1}{1},iField{1})
     fixStimulus.(iField{1}) = task{1}{1}.(iField{1}); 
     if debugFlag, mydisp(sprintf('reassigned %s. \n',iField{1})), end
   else
     mydisp(sprintf('++ using default vals for field %s of task.\n', iField{1}));
   end
end


nItems = fixStimulus.nItems; % number of letters/numbers to display

% create a fixation task - with defaults if not supplied
% seglens depend on parameters stimTime and responseTime
if ~isfield(task{1}{1},'seglen')
  task{1}{1}.segmin = [ones(1,nItems)*fixStimulus.stimTime fixStimulus.responseTime fixStimulus.itiMinTime]; 
  task{1}{1}.segmax = [ones(1,nItems)*fixStimulus.stimTime fixStimulus.responseTime fixStimulus.itiMaxTime]; 
end

% response interval at the end, after all letters have been
% displayed. even if the task is switched off, accept responses (and
% deal with this in the task response callback function)

if ~isfield(task{1}{1},'getResponse')
  task{1}{1}.getResponse = [zeros(1, nItems) 1 0]; 
end

if ~isfield(task{1}{1},'waitForBacktick')
  task{1}{1}.waitForBacktick = 0; 
  disp('[set waitForBackTick to 0]')
end

% init the task
[task{1}{1} myscreen] = initTask(task{1}{1},myscreen,...
				 @fixStartSegmentCallback,...
				 @fixDrawStimulusCallback,...
				 @fixTrialResponseCallback,...
				 @fixTrialStartCallback,...
				 @fixEndTrialCallback);

% keep the correct and incorrect counts
task{1}{1}.correct = 0;
task{1}{1}.incorrect = 0;
task{1}{1}.nTrials = 0;

% set some default options for the text
mglTextSet('',0,0,0,0,0,0,0,0,0);

% create all the stimuli
fixStimulus.n = 0;
disppercent(-inf,'Rendering text');
for colorNum = 1:length(fixStimulus.colors)
  disppercent((colorNum-1)/length(fixStimulus.colors));
  for letterNum = 1:length(fixStimulus.letters)
    mglTextSet(fixStimulus.fontName,fixStimulus.fontSize,fixStimulus.colors{colorNum});
    fixStimulus.n = fixStimulus.n+1;
    fixStimulus.textures(fixStimulus.n) = mglText(fixStimulus.letters(letterNum));
    
    % and keep information about what it is
    fixStimulus.thisLetter(fixStimulus.n) = fixStimulus.letters(letterNum);
    fixStimulus.thisColor(fixStimulus.n) = colorNum;
    fixStimulus.isDigit(fixStimulus.n) = any(fixStimulus.letters(letterNum)=='0':'9');
  end
end

if debugFlag
  mglClearScreen;
  mglBltTexture([fixStimulus.textures(:)], [linspace(-10,10,fixStimulus.n)' zeros(fixStimulus.n,1)]);
  mglFlush;
  pause
  mglClearScreen; mglFlush;
end

% keep info about which textures are digits and which ones are letters
fixStimulus.digitNums = find(fixStimulus.isDigit);
fixStimulus.numDigits = length(fixStimulus.digitNums);
fixStimulus.letterNums = find(~fixStimulus.isDigit);
fixStimulus.numLetters = length(fixStimulus.letterNums);

disppercent(inf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fixation task callback: fixStartTrialCallback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = fixTrialStartCallback(task,myscreen)

global fixStimulus stimulus

myscreen.thistrial.gotResponse = 0;;

% increase the counter of no. of trials
task.nTrials = task.nTrials+1;

% set trace
% myscreen = writetrace(0,task.traceNum,myscreen);

% set color clock 
fixStimulus.fixColorClock = 0;
% set fixation color to white
fixStimulus.fixColor = [1 1 1];
% create a sequence that will be displayed
% first figure out how many letters we will have
numSegs = task.numsegs-2;
% generate a random permutations of letters and numbers
fixStimulus.thisLetters = fixStimulus.letterNums(randperm(fixStimulus.numLetters));
fixStimulus.thisDigits = fixStimulus.digitNums(randperm(fixStimulus.numDigits));
% now make these lists only as long as half the segments 
% we need, since that is how many we will show.
fixStimulus.thisLetters = fixStimulus.thisLetters(1:round(numSegs/2));
fixStimulus.thisDigits = fixStimulus.thisDigits(1:round(numSegs/2));
% put them into a list of all the letters/numbers we want for the
% first and second part of the sequence. We do this so that we can
% have the end of the sequence always have an equal number of
% digits and letters so that there won't be a cue about the task
% from looking at the last four letters/numbers
% first we split the letters into ones that will appear in the
% first or second part. If repeatLoc results in an odd number of
% positions at the end, we will have to have one more letter or
% one more digit. We decide that randomly by adding a 0.5 bias
% this does nothing if there is an even numer of spaces. Note,
% if this doesn't make any sense and you won't to avoid this,
% just make repeatLoc an even number so that you can have the
% same number of digits and letters at the end
letterBias = (rand>0.5)*0.5;
digitBias = (1-2*letterBias)*0.5;
fixStimulus.thisLettersFirst = fixStimulus.thisLetters(1:floor(letterBias+(numSegs-fixStimulus.repeatLoc)/2));
fixStimulus.thisLettersSecond = fixStimulus.thisLetters((floor(letterBias+(numSegs-fixStimulus.repeatLoc)/2)+1):length(fixStimulus.thisLetters));
% then we do the same for the digits. 
fixStimulus.thisDigitsFirst = fixStimulus.thisDigits(1:floor(digitBias+(numSegs-fixStimulus.repeatLoc)/2));
fixStimulus.thisDigitsSecond = fixStimulus.thisDigits((floor(digitBias+(numSegs-fixStimulus.repeatLoc)/2)+1):length(fixStimulus.thisDigits));
% now we construct the combined sequence and randomize both
fixStimulus.thisSequenceFirst = [fixStimulus.thisLettersFirst fixStimulus.thisDigitsFirst];
fixStimulus.thisSequenceSecond = [fixStimulus.thisLettersSecond fixStimulus.thisDigitsSecond];
% and randomize these two
fixStimulus.thisSequenceFirst = fixStimulus.thisSequenceFirst(randperm(length(fixStimulus.thisSequenceFirst)));
fixStimulus.thisSequenceSecond = fixStimulus.thisSequenceSecond(randperm(length(fixStimulus.thisSequenceSecond)));
% and put them together
fixStimulus.thisSequence = [fixStimulus.thisSequenceFirst fixStimulus.thisSequenceSecond];
% and cut it down to size, this will only be necessary if numSegs
% is not an even number
fixStimulus.thisSequence = fixStimulus.thisSequence(((length(fixStimulus.thisSequence)-numSegs)+1):length(fixStimulus.thisSequence));
% now get a random position at the end of the sequence to pull the
% repeat from
fixStimulus.thisRepeatLoc = numSegs-ceil(fixStimulus.repeatLoc*rand(1))+1;
% and get a position from the rest of the segment that does not include
% this location or its neighbors to choose the repeated location from
possibleLocs = setdiff(1:(numSegs-fixStimulus.repeatLoc),fixStimulus.thisRepeatLoc+(-1:1));
fixStimulus.thisRepeatLoc(2) = possibleLocs(ceil(rand(1)*length(possibleLocs)));
% get what the repeated character is
fixStimulus.thisRepeat = fixStimulus.thisSequence(fixStimulus.thisRepeatLoc(1));
% place the repeat in the repeat location
fixStimulus.thisSequence(fixStimulus.thisRepeatLoc(2)) = fixStimulus.thisRepeat;
% record what the target is
fixStimulus.thisTarget = fixStimulus.isDigit(fixStimulus.thisRepeat);
% and print out information
disp(sprintf('target=%c %i at pos %i %i',fixStimulus.thisLetter(fixStimulus.thisRepeat),fixStimulus.thisColor(fixStimulus.thisRepeat),fixStimulus.thisRepeatLoc(1),fixStimulus.thisRepeatLoc(2)));

% start the subsidiary task: dots
stimulus.startSubsidiary = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fixation task callback: fixStartSegmentCallback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = fixStartSegmentCallback(task,myscreen)

global fixStimulus;

% first segment, init the trial
if (task.thistrial.thisseg < (task.numsegs-1))
  seqNum = fixStimulus.thisSequence(task.thistrial.thisseg);
  mydisp(sprintf('%c %i, ',fixStimulus.thisLetter(seqNum),fixStimulus.thisColor(seqNum)));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fixation task callback: fixTrialStimulusCallback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = fixDrawStimulusCallback(task,myscreen)

global fixStimulus;

% update the color (this keeps the display of correct/incorrect
% as color of dots for as many frames as
% fixStimulus.dots.colorClockLen is set to.
fixStimulus.fixColorClock = fixStimulus.fixColorClock-1;
if (fixStimulus.fixColorClock <= 0)
  fixStimulus.fixColor = [1 1 1];
end


% draw oval behind fixation
% mglFillOval(myscreen.fix.loc(1), myscreen.fix.loc(2), myscreen.fix.disksz, myscreen.background);

% mglClearScreen;
mglGluDisk(0,0,myscreen.fix.diskSize,myscreen.background,60);
mglStencilSelect(myscreen.stencil.fixation)

if (task.thistrial.thisseg >= 1) && (task.thistrial.thisseg < (task.numsegs-1))
  % draw the correct text
  mglBltTexture(fixStimulus.textures(fixStimulus.thisSequence(task.thistrial.thisseg)),myscreen.fix.loc);
else
  % along x=0
  mglLines2(-myscreen.fix.size(1)+myscreen.fix.loc(1),myscreen.fix.loc(2),...
	    +myscreen.fix.size(1)+myscreen.fix.loc(1),myscreen.fix.loc(2),...
	    myscreen.fix.linewidth,fixStimulus.fixColor);
  % along y=0
  mglLines2(myscreen.fix.loc(1),-myscreen.fix.size(2)+myscreen.fix.loc(2),...
	    myscreen.fix.loc(1),+myscreen.fix.size(2)+myscreen.fix.loc(2),...
	    myscreen.fix.linewidth,fixStimulus.fixColor);
end
mglStencilSelect(0)

%global stimulus;
%if task.getResponse(task.thistrial.thisseg)
%  if ~stimulus.dots.colorClock
%    stimulus.dots.color = [0 1 1];
%  end
%else
%  if ~stimulus.dots.colorClock
%    stimulus.dots.color = [1 1 1];
%  end
%end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fixation task callback: fixTrialResponseCallback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = fixTrialResponseCallback(task,myscreen)

global fixStimulus;

% update number of responses happens in task code

if (task.thistrial.gotResponse == task.whichResponse-1)
  mydisp(sprintf('[sanghun task] response interval: %i\n', task.thistrial.gotResponse))
  if (task.thistrial.buttonState(1) & fixStimulus.thisTarget) | (task.thistrial.buttonState(2) & ~fixStimulus.thisTarget)
    fixStimulus.fixColor = [0 1 0];
    fixStimulus.fixColorClock = fixStimulus.fixColorClockLen;
    % myscreen = writeTrace(1,task.traceNum,myscreen);
    task.correct = task.correct+1;
    mydisp(sprintf('correct\n'));
  else
    fixStimulus.fixColor = [1 0 0];
    fixStimulus.fixColorClock = fixStimulus.fixColorClockLen;
    % myscreen = writeTrace(-1,task.traceNum,myscreen);
    task.incorrect = task.incorrect+1;
    mydisp(sprintf('incorrect\n'));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the end of each TRIAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = fixEndTrialCallback(task, myscreen)

global stimulus

% this enables syncing between the two tasks. whichever task gets to
% the end first sets the sync pulse

% set synching pulse
% myscreen.sync = 1;
stimulus.endSubsidiary = 1;

disp('stimulus.endSubsidiary = 1')
disp('end of letter trial')



