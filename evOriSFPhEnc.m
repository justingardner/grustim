% evOriSF.m
%
%      usage: evOriSF()
%         by: eli merriam
%       date: 07/26/22
%    purpose: 
%
function retval = evOriSFPhEnc(varargin)

% check arguments
if ~any(nargin == [0])
  help evOriSF
  return
end

% evaluate the input arguments
getArgs(varargin, [], 'verbose=0');

% set default parameters
if ieNotDefined('atScanner'),atScanner = 0;end
if ieNotDefined('recompITI'),recompITI = 0;end

mglSetSID(-1);

% initalize the screen
myscreen.background = 'gray';
myscreen.autoCloseScreen = 0;
myscreen.allowpause = 1;
myscreen.saveData = 0;
% myscreen.displayName = '3tb';
% myscreen.displayName = 'test';
% myscreen.displayName = 'fMRIprojFlex';
myscreen.displayName = 'fMRIproj_akuo2';
myscreen = initScreen(myscreen);

global stimulus;
myscreen = initStimulus('stimulus',myscreen);
phasedur = 0.25;
nseg = 64;
% orientation
orientation = [0 90];
stimulus.orientation = orientation;
task{1}{1}.parameter.orientation = orientation;
% ascending, descending sf conditions
stimulus.nsfs = 16;
stimulus.sf = 2.^linspace(-3,2,stimulus.nsfs);
stimulus.sfdirection = [1 -1];
task{1}{1}.parameter.sfdirection = stimulus.sfdirection;
% size
stimulus.height = 10;
stimulus.width = 10;

task{1}{1}.random = 0;
task{1}{1}.numTrials = Inf;
task{1}{1}.collectEyeData = true;
task{1}{1}.waitForBacktick = 1;
task{1}{1}.segmin = [repmat(phasedur, 1, nseg) nan];
task{1}{1}.segmax = [repmat(phasedur, 1, nseg) nan];
% duration of the ISI's
task{1}{1}.segdur{nseg+1} = [4] - 0.1;
if recompITI
  n = 10000;
  probs = 1+exprnd(2,n,1);
  task{1}{1}.segprob{nseg+1} = hist(probs, task{1}{1}.segdur{nseg+1})/n;
else
  % hard code
  % probs = [0.4585 0.2862 0.1356 0.0652 0.0309 0.0119 0.0066 0.0028 0.0013 0.001];
  probs = [1];
  task{1}{1}.segprob{nseg+1} = probs; 
end
% sync to scanner
task{1}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
if atScanner
  task{1}{1}.synchToVol(end) = 1;
end

global fixStimulus
fixStimulus.diskSize = 0.5;
fixStimulus.fixWidth = 0.8;
fixStimulus.fixLineWidth = 3;
[task{2} myscreen] = fixStairInitTask(myscreen);

% initialize the task
for phaseNum = 1:length(task{1})
  [task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end

% do our initialization which creates the gratings
stimulus = myInitStimulus(stimulus,myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mglClearScreen(); mglFlush;
phaseNum = 1;
while (phaseNum <= length(task)) && ~myscreen.userHitEsc
  % update the ori task
  [task{1} myscreen] = updateTask(task{1}, myscreen, 1);
  % update the fixation task
  [task{2} myscreen phaseNum] = updateTask(task{2},myscreen, 1);
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

% randomize the current phase of the stimulus
newPhase = ceil(rand(1)*stimulus.numPhases);
while stimulus.phaseNum == newPhase;
  newPhase = ceil(rand(1)*stimulus.numPhases);
end
stimulus.phaseNum = newPhase;

% set the current segment orientation
if task.thistrial.orientation == 0
    stimulus.oriInd = 1;
elseif task.thistrial.orientation == 90
    stimulus.oriInd = 2;
end

% determine what sf we are using based on segnum and sfdirection
if task.thistrial.sfdirection == 1
    stimulus.sfInd = ceil(task.thistrial.thisseg/4);
elseif task.thistrial.sfdirection == -1
    stimulus.sfInd = (stimulus.nsfs+1) - ceil(task.thistrial.thisseg/4);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus;

% clear the screen
mglClearScreen;

% draw the texture
if task.thistrial.thisseg<65
  mglBltTexture(stimulus.tex{stimulus.oriInd,stimulus.sfInd,stimulus.phaseNum}, [0 0 stimulus.height stimulus.height], 0, 0, 0);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen,task)

stimulus.pixRes = min(myscreen.screenHeight/myscreen.imageHeight, myscreen.screenWidth/myscreen.imageWidth);

% which phases (of the grating) we will have
stimulus.phaseNum = 1;
stimulus.numPhases = 16;
stimulus.phases = 0:(360-0)/stimulus.numPhases:360;

if isfield(stimulus, 'tex')
  disp(sprintf('(evORISF) Attention: Using precomputed stimulus textures!!'));
else  
  disppercent(-inf,'Creating the stimulus textures');
  for iOri=1:length(stimulus.orientation)
    for iSF=1:length(stimulus.sf)
      for iPhase=1:length(stimulus.phases)
          
        % make a grating  but now scale it
        grating = mglMakeGrating(stimulus.width, stimulus.height, stimulus.sf(iSF), stimulus.orientation(iOri), stimulus.phases(iPhase), stimulus.pixRes, stimulus.pixRes);

        % make a circular aperture
        grating = grating .*  mkDisc(size(grating), (length(grating)/2)-2, (size(grating)+1)/2, 1);

        % scale to range of display
        grating = 255*(grating+1)/2;

        % create a texture
        stimulus.tex{iOri, iSF, iPhase} = mglCreateTexture(grating, [], 1);
      end
    end
    disppercent(iOri/length(stimulus.orientation));
  end
  disppercent(inf);
end
