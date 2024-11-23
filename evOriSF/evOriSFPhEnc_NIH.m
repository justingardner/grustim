% evOriSFPhEnc.m
%
%      usage: evOriSFPhEnc(varargin)
%         by: austin kuo
%       date: 03/08/23
%    purpose: display phase encoded orientation stimuli with specified spatial frequency
%
%       args: varargin - oridir=1: clockwise
%                        oridir=-1: counter clockwise
%                        sfq=0.5: spatial frequency = 0.5 cpd

function retval = evOriSFPhEnc_NIH(varargin)

% % check arguments
if ~any(nargin == [2:5])
  help evOriSFPhEnc_NIH
  return
end

% evaluate the input arguments
getArgs(varargin, [], 'verbose=0');

% set default parameters
if ieNotDefined('atScanner'),atScanner = 0;end
if ieNotDefined('saveParam'),saveParam = 0;end
if ieNotDefined('displayName')
    myscreen.displayName = 'test';
else
    myscreen.displayName = displayName;
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

% spatial frequency
stimulus.sf = sfq;
task{1}{1}.parameter.sf = stimulus.sf;

% ccw, cw ori conditions
stimulus.segdur = 0.25;
if ieNotDefined('noris')
    stimulus.noris = 8;
else
    stimulus.noris = noris;
end
stimulus.trialdur = 24/stimulus.noris;
stimulus.nseg = stimulus.trialdur/stimulus.segdur;
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

% gradient timings
% the total time of the gradient is stimulus.nGradSteps * 1/frameUpdateFreq * 2 (e.g., 1/60); keep this value < 14 (or 0.467 seconds) so that we don't run out of single seg boundary if each step update is 1/60
if ieNotDefined('gradsteps')
    stimulus.nGradSteps = 6;
else
    stimulus.nGradSteps = gradsteps;
end
stimulus.gradientUpdateFrameHz = 60;
stimulus.gradientTimesForward = linspace(1/stimulus.gradientUpdateFrameHz,(stimulus.nGradSteps+1)/stimulus.gradientUpdateFrameHz,stimulus.nGradSteps+1)';
stimulus.gradientTimesBackward = stimulus.segdur - linspace(1/stimulus.gradientUpdateFrameHz,(stimulus.nGradSteps+1)/stimulus.gradientUpdateFrameHz,stimulus.nGradSteps+1)';

% size
stimulus.height = 20;
stimulus.width = 20;
stimulus.aperOuterHeight = 10; % minor axis diameter
stimulus.aperOuterWidth = 10; % major axis diameter
stimulus.outerHeightRatio = stimulus.aperOuterHeight/stimulus.height;
stimulus.outerWidthRatio = stimulus.aperOuterWidth/stimulus.width;

stimulus.aperInnerHeight = 6; % minor axis diameter
stimulus.aperInnerWidth = 6; % major axis diameter
stimulus.innerHeightRatio = stimulus.aperInnerHeight/stimulus.height;
stimulus.innerWidthRatio = stimulus.aperInnerWidth/stimulus.width;

task{1}{1}.random = 0;
task{1}{1}.numTrials = Inf;
task{1}{1}.collectEyeData = true;
task{1}{1}.waitForBacktick = 1;
task{1}{1}.seglen = stimulus.segdur*ones(1,4*stimulus.trialdur); % [repmat(0.25, 1, nseg)];

% sync to scanner
task{1}{1}.synchToVol = zeros(size(task{1}{1}.seglen));
if atScanner
  task{1}{1}.fudgeLastVolume = 1;
  task{1}{1}.seglen(end) = task{1}{1}.seglen(end)-0.06;
  task{1}{1}.synchToVol(end) = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%% define fixStimulus 
%%%%%%%%%%%%%%%%%%%%%%%%

clear global fixStimulus
global fixStimulus
fixStimulus.fixTask = 'rsvp'; % either fixStair or rsvp
if strcmpi(fixStimulus.fixTask,'fixStair')
    fixStimulus.diskSize = 1.5; % radius of fixation disk
    if fixStimulus.diskSize*2 > stimulus.aperInnerHeight
        warning('The disk size for the fixation cross is larger than your inner aperture size.')
    end
    fixStimulus.fixWidth = 1; % arm length = half of this value
    fixStimulus.fixLineWidth = 0.2;
    [task{2} myscreen] = fixStairInitTaskMetal(myscreen);
elseif strcmpi(fixStimulus.fixTask,'rsvp')
    fixStimulus.diskSize = 1.5;
    fixStimulus.rate = 2; % letters/numbers per second (default=2)
    fixStimulus.RSVPSize = 22; % size of text (default=22)
    [task{2} myscreen] = fixRSVP(myscreen);
end

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

% calculate hit rate and false alarm rate
nTrials = task{2}{1}.trialnum;
thisScanTotalPossibleHits = sum(fixStimulus.targetSeq(1:nTrials));
thisScanTotalPossibleFA = nTrials - thisScanTotalPossibleHits;
hitRate = sum(task{2}{1}.randVars.hit) / thisScanTotalPossibleHits;
faRate = sum(task{2}{1}.randVars.falseAlarm) / thisScanTotalPossibleFA;

% avoid infinite or undefined z-scores
hitRate = max(min(hitRate, 0.99),0.01);
faRate = max(min(faRate, 0.99),0.01);

zHit = norminv(hitRate);
zFA = norminv(faRate);
dprime = zHit - zFA;
if dprime < 2
    encouragingMessage = sprintf('(This score relies on the ratio of hit rate to false alarm rate. Try to aim for a d'' of 2!)');
else
    encouragingMessage = sprintf('d'' > 2 - great job, keep it up!');
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

mglTextSet('Helvetica',32,[1 1 1]);
mglTextDraw(sprintf('Hit rate: %0.2f%%',hitRate*100),[0 1]);
mglTextDraw(sprintf('False alarm rate: %0.2f%%',faRate*100),[0 0]);
mglTextDraw(sprintf('d'' = %0.2f',dprime),[0 -1])
mglTextDraw(encouragingMessage,[0 -2])
mglFlush;

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

if task.thistrial.thisseg == stimulus.nseg || task.thistrial.thisseg == 1
    stimulus.segStart = mglGetSecs;
    myscreen = writeTrace(stimulus.segStart,myscreen.stimtrace+1,myscreen,1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus;

% clear the screen
mglClearScreen;

oriInd = find(task.parameter.orientations==task.thistrial.orientations);
% if stimulus.oridirection == -1
%     oriInd = stimulus.reverseIdx(oriInd);
% end

% draw the texture
if isfield(stimulus.tex,'gradient')
    if task.thistrial.thisseg < stimulus.nseg && task.thistrial.thisseg > 1
        %mglBltTexture([stimulus.tex.grating{stimulus.phaseNum} stimulus.tex.aperture], [0 0 stimulus.height stimulus.height; 0 0 stimulus.height*2 stimulus.height*2], 0, 0, [stimulus.orientations(oriInd) 0]);
        mglBltTexture(stimulus.tex.grating{stimulus.phaseNum}, [0 0 stimulus.height stimulus.height], 0, 0, stimulus.orientations(oriInd));
        mglBltTexture(stimulus.tex.aperture, [0 0 stimulus.height*2 stimulus.height*2], 0, 0, 0);
    elseif task.thistrial.thisseg == 1 % if we are starting the trial, make the contrast ramp up
        frameTime = mglGetSecs(stimulus.segStart);
        [~,gIdx] = pdist2(stimulus.gradientTimesForward,frameTime,'euclidean','Smallest',1);
        if gIdx <= stimulus.nGradSteps
            mglBltTexture(stimulus.tex.gradient{gIdx,stimulus.phaseNum}, [0 0 stimulus.height stimulus.height], 0, 0, stimulus.orientations(oriInd));
            mglBltTexture(stimulus.tex.aperture, [0 0 stimulus.height*2 stimulus.height*2], 0, 0, 0);
        else
            mglBltTexture(stimulus.tex.grating{stimulus.phaseNum}, [0 0 stimulus.height stimulus.height], 0, 0, stimulus.orientations(oriInd));
            mglBltTexture(stimulus.tex.aperture, [0 0 stimulus.height*2 stimulus.height*2], 0, 0, 0);
        end
    elseif task.thistrial.thisseg == stimulus.nseg % if we are ending the trial, make the contrast ramp down
        frameTime = mglGetSecs(stimulus.segStart);
        [~,gIdx] = pdist2(stimulus.gradientTimesBackward,frameTime,'euclidean','Smallest',1);
        if gIdx <= stimulus.nGradSteps
            mglBltTexture(stimulus.tex.gradient{gIdx,stimulus.phaseNum}, [0 0 stimulus.height stimulus.height], 0, 0, stimulus.orientations(oriInd));
            mglBltTexture(stimulus.tex.aperture, [0 0 stimulus.height*2 stimulus.height*2], 0, 0, 0);
        else
            mglBltTexture(stimulus.tex.grating{stimulus.phaseNum}, [0 0 stimulus.height stimulus.height], 0, 0, stimulus.orientations(oriInd));
            mglBltTexture(stimulus.tex.aperture, [0 0 stimulus.height*2 stimulus.height*2], 0, 0, 0);
        end
    end
else
    mglBltTexture(stimulus.tex.grating{stimulus.phaseNum}, [0 0 stimulus.height stimulus.height], 0, 0, stimulus.orientations(oriInd));
    mglBltTexture(stimulus.tex.aperture, [0 0 stimulus.height*2 stimulus.height*2], 0, 0, 0);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the grating stimulus %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen,task)

stimulus.pixRes = min(myscreen.screenHeight/myscreen.imageHeight, myscreen.screenWidth/myscreen.imageWidth);

% which phases (of the grating) we will have
stimulus.phaseNum = 1;
stimulus.numPhases = 8;
stimulus.phases = 0:(360-0)/stimulus.numPhases:360;
stimulus.phases = stimulus.phases(1:end-1);

disppercent(-inf,'Creating the stimulus textures');

% create the grating texture (just need phase and sf, mglBltTexture will handle the rotations)
for iPhase=1:length(stimulus.phases)

    % make a grating  but now scale it
    grating = mglMakeGrating(stimulus.width, stimulus.height, stimulus.sf, 0, stimulus.phases(iPhase), stimulus.pixRes, stimulus.pixRes);

    % scale to range of display
    grating = 255*(grating+1)/2;

    % create a texture
    stimulus.tex.grating{iPhase} = mglCreateTexture(grating, [], 1);

    disppercent(iPhase/length(stimulus.phases));
end

% create the gradient textures
gradientVals = linspace(0,1,stimulus.nGradSteps+1);
gradientVals = gradientVals(1:end-1);
for iGradient = 1:length(gradientVals)
    for iPhase=1:length(stimulus.phases)

        % make a grating  but now scale it
        grating = gradientVals(iGradient) * mglMakeGrating(stimulus.width, stimulus.height, stimulus.sf, 0, stimulus.phases(iPhase), stimulus.pixRes, stimulus.pixRes);

        % scale to range of display
        grating = 255*(grating+1)/2;

        % create a texture
        stimulus.tex.gradient{iGradient,iPhase} = mglCreateTexture(grating, [], 1);
    end
end

disppercent(inf);

% create the aperture texture - make a elliptical (circular) aperture
% grating = grating .*  mkDisc(size(grating), (length(grating)/2)-2, (size(grating)+1)/2, 1);
stimSize = size(grating,1);
apertureOuter = circStim([stimSize/2*stimulus.outerHeightRatio stimSize/2*stimulus.outerWidthRatio],[stimSize stimSize],[ceil(stimSize/2) ceil(stimSize/2)]);
apertureInner = ~circStim([stimSize/2*stimulus.innerHeightRatio stimSize/2*stimulus.innerWidthRatio],[stimSize stimSize],[ceil(stimSize/2) ceil(stimSize/2)]);
apertureAlpha = ~and(apertureOuter,apertureInner);
apertureAlpha = padarray(apertureAlpha,[(stimSize-1)/2 (stimSize-1)/2],1,'both'); % pad the array so that we block out the grating edges

% make rgb matrix (n x m x 3)
apertureRGB = 0.5 * ones(size(apertureAlpha,1),size(apertureAlpha,2),3);

% tack on the alpha values
aperture = 255 * cat(3,apertureRGB,apertureAlpha);
% keyboard
% create a texture
stimulus.tex.aperture = mglCreateTexture(aperture, [], 1);

