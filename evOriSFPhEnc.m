% evOriSFPhEnc.m
%
%      usage: evOriSFPhEnc(varargin)
%         by: austin kuo
%       date: 01/25/23
%    purpose: display phase encoded spatial frequency stimuli with specified orientation and low-high/high-low SF steps
%
%       args: varargin - ori=0: horizontal grating
%                        ori=90: vertical grating
%                        sfdir=1: low to high SFs
%                        sfdir=-1: high to low SFs

function retval = evOriSFPhEnc(varargin)

% % check arguments
if ~any(nargin == [2:4])
  help evOriSFPhEnc
  return
end

% evaluate the input arguments
getArgs(varargin, [], 'verbose=0');

% set default parameters
if ieNotDefined('atScanner'),atScanner = 0;end
if ieNotDefined('saveParam'),saveParam = 0;end

if ieNotDefined('ori') || ieNotDefined('sfdir')
    error('Specify a grating orientation and SF direction (e.g. ''ori=0'', ''sfdir=1'')')
end

% initalize the screen
myscreen.background = 'gray';
myscreen.autoCloseScreen = 0;
myscreen.allowpause = 1;
myscreen.saveData = saveParam;
myscreen.displayName = 'fMRIproj_akuo2';
myscreen = initScreen(myscreen);

global stimulus;
myscreen = initStimulus('stimulus',myscreen);
phasedur = 0.25;
nseg = 4;

% orientation
if ori == 0
    orientation = 0;
    stimulus.oriInd = 1;
elseif ori == 90
    orientation = 90;
    stimulus.oriInd = 2;
else
    error('Specify ''ori'' as either a horizontal or vertical orientation (0 or 90)')
end
stimulus.orientations = [0 90];
task{1}{1}.parameter.orientation = orientation;

% ascending, descending sf conditions
stimulus.nsfs = 24;
stimulus.reverseIdx = fliplr(1:stimulus.nsfs);
if sfdir == 1
    stimulus.sfdirection = 1;
    stimulus.sf = 2.^linspace(-2,3,stimulus.nsfs);
    
elseif sfdir == -1
    stimulus.sfdirection = -1;
    stimulus.sf = fliplr(2.^linspace(-2,3,stimulus.nsfs));
    
else
    error('Specify ''sfdir'' as either ''1'' or ''-1'' to indicate SF direction')
end
stimulus.sf = [stimulus.sf(stimulus.nsfs/2 + 1:end) stimulus.sf(1:stimulus.nsfs/2)]; % start on half-cycle
task{1}{1}.parameter.sf = stimulus.sf;

% size
stimulus.height = 20;
stimulus.width = 20;

task{1}{1}.random = 0;
task{1}{1}.numTrials = Inf;
task{1}{1}.collectEyeData = true;
task{1}{1}.waitForBacktick = 1;
task{1}{1}.seglen = [repmat(0.25, 1, nseg)];

% sync to scanner
task{1}{1}.synchToVol = zeros(size(task{1}{1}.seglen));
if atScanner
  task{1}{1}.fudgeLastVolume = 1;
  task{1}{1}.seglen(end) = task{1}{1}.seglen(end)-0.1;
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus;

% clear the screen
mglClearScreen;

sfInd = find(task.parameter.sf==task.thistrial.sf);
if stimulus.sfdirection == -1
    sfInd = stimulus.reverseIdx(sfInd);
end

% draw the texture
if task.thistrial.thisseg<5 
    mglBltTexture(stimulus.tex{stimulus.oriInd,sfInd,stimulus.phaseNum}, [0 0 stimulus.height stimulus.height], 0, 0, 0);
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
stimulus.phases = stimulus.phases(1:end-1);

if isfield(stimulus, 'tex')
    fprintf('(evOriSFPhEnc) Attention: Using precomputed stimulus textures!!');
else
    disppercent(-inf,'Creating the stimulus textures');
    for iOri=1:length(stimulus.orientations)
        for iSF=1:length(stimulus.sf)
            for iPhase=1:length(stimulus.phases)
                
                % make a grating  but now scale it
                grating = mglMakeGrating(stimulus.width, stimulus.height, stimulus.sf(iSF), stimulus.orientations(iOri), stimulus.phases(iPhase), stimulus.pixRes, stimulus.pixRes);
                
                % make a circular aperture
                grating = grating .*  mkDisc(size(grating), (length(grating)/2)-2, (size(grating)+1)/2, 1);
                
                % scale to range of display
                grating = 255*(grating+1)/2;
                
                % create a texture
                stimulus.tex{iOri, iSF, iPhase} = mglCreateTexture(grating, [], 1);
                
            end
        end
        disppercent(iOri/length(stimulus.orientations));
    end
    disppercent(inf);
end
