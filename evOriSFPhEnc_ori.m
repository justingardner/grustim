% evOriSFPhEnc.m
%
%      usage: evOriSFPhEnc(varargin)
%         by: austin kuo
%       date: 03/08/23
%    purpose: display phase encoded orientation stimuli with specified spatial frequency
%
%       args: varargin - oridir=1: clockwise
%                        oridir=-1: counter clockwise
%                        sf=0.5: low SF
%                        sf=1.5: high SF

function retval = evOriSFPhEnc_ori(varargin)

% % check arguments
if ~any(nargin == [2:5])
  help evOriSFPhEnc_ori
  return
end

% evaluate the input arguments
getArgs(varargin, [], 'verbose=0');

% set default parameters
if ieNotDefined('atScanner'),atScanner = 0;end
if ieNotDefined('saveParam'),saveParam = 0;end
if ieNotDefined('screenParam')
    myscreen.displayName = 'fMRIproj_akuo2';
else
    myscreen.displayName = screenParam;
end

if ieNotDefined('sfq') || ieNotDefined('oridir')
    error('Specify a grating SF and orientation direction (e.g. ''sfq=0.5'', ''oridir=1'')')
end

% initalize the screen
myscreen.background = 'gray';
myscreen.autoCloseScreen = 0;
myscreen.allowpause = 1;
myscreen.saveData = saveParam;
myscreen = initScreen(myscreen);

global stimulus;
myscreen = initStimulus('stimulus',myscreen);
phasedur = 0.25;
nseg = 4;

% spatial frequency
stimulus.sf = [0.5 4];
if sfq == stimulus.sf(1)
    stimulus.sfInd = 1;
elseif sfq == stimulus.sf(2)
    stimulus.sfInd = 2;
else
    error('Specify ''sfq'' as either (%0.1f or %0.1f)',stimulus.sf(1),stimulus.sf(2))
end
task{1}{1}.parameter.sf = stimulus.sf(stimulus.sfInd);

% ccw, cw ori conditions
stimulus.noris = 24;
stimulus.reverseIdx = fliplr(1:stimulus.noris);
if oridir == 1
    stimulus.oridirection = 1;
    stimulus.orientations = linspace(0,180 - 180/stimulus.noris,stimulus.noris);
    
elseif oridir == -1
    stimulus.oridirection = -1;
    stimulus.orientations = fliplr(linspace(0,180 - 180/stimulus.noris,stimulus.noris));
    
else
    error('Specify ''oridir'' as either ''1'' or ''-1'' to indicate orientation direction')
end
stimulus.orientations = [stimulus.orientations(stimulus.noris/2 + 1:end) stimulus.orientations(1:stimulus.noris/2)]; % start on half-cycle
task{1}{1}.parameter.orientations = stimulus.orientations;

% size
stimulus.height = 20;
stimulus.width = 20;
stimulus.aperOuterHeight = 16; % minor axis diameter
stimulus.aperOuterWidth = 19; % major axis diameter
stimulus.outerHeightRatio = stimulus.aperOuterHeight/stimulus.height;
stimulus.outerWidthRatio = stimulus.aperOuterWidth/stimulus.width;

stimulus.aperInnerHeight = 9; % minor axis diameter
stimulus.aperInnerWidth = 9; % major axis diameter
stimulus.innerHeightRatio = stimulus.aperInnerHeight/stimulus.height;
stimulus.innerWidthRatio = stimulus.aperInnerWidth/stimulus.width;

task{1}{1}.random = 0;
task{1}{1}.numTrials = Inf;
task{1}{1}.collectEyeData = true;
task{1}{1}.waitForBacktick = 1;
task{1}{1}.seglen = [0.25 0.25 0.25 0.15 0.1]; % [repmat(0.25, 1, nseg)];

% sync to scanner
task{1}{1}.synchToVol = zeros(size(task{1}{1}.seglen));
if atScanner
  task{1}{1}.fudgeLastVolume = 1;
  task{1}{1}.seglen(end) = task{1}{1}.seglen(end)-0.06;
  task{1}{1}.synchToVol(end) = 1;
end

global fixStimulus
fixStimulus.diskSize = 5;
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

oriInd = find(task.parameter.orientations==task.thistrial.orientations);
if stimulus.oridirection == -1
    oriInd = stimulus.reverseIdx(oriInd);
end

% draw the texture
if task.thistrial.thisseg<5
    mglBltTexture(stimulus.tex{oriInd,stimulus.sfInd,stimulus.phaseNum}, [0 0 stimulus.height stimulus.height], 0, 0, 0);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen,task)

stimulus.pixRes = min(myscreen.screenHeight/myscreen.imageHeight, myscreen.screenWidth/myscreen.imageWidth);

% which phases (of the grating) we will have
stimulus.phaseNum = 1;
stimulus.numPhases = 8;
stimulus.phases = 0:(360-0)/stimulus.numPhases:360;
stimulus.phases = stimulus.phases(1:end-1);

% test - remove later
% stimulus.numPhases = 2;
% stimulus.phases = [0 0];

if isfield(stimulus, 'tex')
    fprintf('(evOriSFPhEnc) Attention: Using precomputed stimulus textures!!');
else
    disppercent(-inf,'Creating the stimulus textures');
    for iOri=1:length(stimulus.orientations)
        for iSF=1:length(stimulus.sf)
            for iPhase=1:length(stimulus.phases)
                
                % make a grating  but now scale it
                grating = mglMakeGrating(stimulus.width, stimulus.height, stimulus.sf(iSF), stimulus.orientations(iOri), stimulus.phases(iPhase), stimulus.pixRes, stimulus.pixRes);
                
                % make a elliptical (circular) aperture
                % grating = grating .*  mkDisc(size(grating), (length(grating)/2)-2, (size(grating)+1)/2, 1);
                stimSize = size(grating,1);
                apertureOuter = circ([stimSize/2*stimulus.outerHeightRatio stimSize/2*stimulus.outerWidthRatio],[stimSize stimSize],[ceil(stimSize/2) ceil(stimSize/2)]);
                apertureInner = ~circ([stimSize/2*stimulus.innerHeightRatio stimSize/2*stimulus.innerWidthRatio],[stimSize stimSize],[ceil(stimSize/2) ceil(stimSize/2)]);
                aperture = and(apertureOuter,apertureInner);
                grating = grating .* aperture;
                
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
