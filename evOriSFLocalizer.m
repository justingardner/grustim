% evOriSF.m
%
%      usage: evOriSF()
%         by: eli merriam
%       date: 07/26/22
%    purpose:
%
function retval = evOriSFLocalizer(varargin)

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

% initalize the screen
myscreen.background = 'gray';
myscreen.autoCloseScreen = 0;
myscreen.allowpause = 1;
myscreen.saveData = 0;
% myscreen.displayName = '3tb';
% myscreen.displayName = 'laptop';
myscreen = initScreen(myscreen);

global stimulus;
myscreen = initStimulus('stimulus',myscreen);
% orientation
orientation = linspace(0, 90, 2);
% orientation = orientation(1:end-1);
stimulus.orientation = orientation;
task{1}{1}.parameter.orientation = stimulus.orientation;
% spatial frequency
stimulus.sf = [0.25 0.5 1 2 4];
task{1}{1}.parameter.sf = stimulus.sf;
% ring/anti-ring
stimulus.ring = [1 2];
task{1}{1}.parameter.ring = stimulus.ring;
% location
task{1}{1}.parameter.location = 0;
% size
stimulus.height = 5.5;
stimulus.width = 5.5;


task{1}{1}.random = 1;
task{1}{1}.numTrials = Inf;
task{1}{1}.collectEyeData = true;
task{1}{1}.waitForBacktick = 1;
task{1}{1}.segmin = [repmat(0.25, 1, 12) nan];
task{1}{1}.segmax = [repmat(0.25, 1, 12) nan];
% duration of the ISI's
task{1}{1}.segdur{13} = [1.5:1.5:4.5]-0.1;
if recompITI
    n = 10000;
    probs = 1+exprnd(2,n,1);
    task{1}{1}.segprob{13} = hist(probs, task{1}{1}.segdur{13})/n;
else
    % hard code
    % probs = [0.4585 0.2862 0.1356 0.0652 0.0309 0.0119 0.0066 0.0028 0.0013 0.001];
    probs = [0.5 0.25 0.25];
    task{1}{1}.segprob{13} = probs;
end
% sync to scanner
task{1}{1}.synchToVol = zeros(size(task{1}{1}.segmin));
if atScanner
    task{1}{1}.synchToVol(end) = 1;
end


global fixStimulus
fixStimulus.diskSize = 0.25;
fixStimulus.fixWidth = 0.4;
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
oriInd = find(task.parameter.orientation==task.thistrial.orientation);
ringInd = find(task.parameter.ring==task.thistrial.ring);

% draw the texture
if task.thistrial.thisseg<13
    if ringInd == 2
        mglBltTexture(stimulus.tex{oriInd,sfInd,stimulus.phaseNum,ringInd}, [task.thistrial.location 0 myscreen.imageWidth myscreen.imageHeight], 0, 0, 0);
    else
        
        mglBltTexture(stimulus.tex{oriInd,sfInd,stimulus.phaseNum,ringInd}, [task.thistrial.location 0 stimulus.height stimulus.height], 0, 0, 0);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen,task)

stimulus.pixRes = min(myscreen.screenHeight/myscreen.imageHeight, myscreen.screenWidth/myscreen.imageWidth);

% which phases we will have
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
                
                for iRing = 1:length(stimulus.ring) % for ring/antiring
                    if stimulus.ring(iRing) == 1 % make a ring
                        
                        % make a ring
                        grating = grating .*  (mkDisc(size(grating), (length(grating)/2)-2, (size(grating)+1)/2, 0, [1,0]) + ... % outer
                                               mkDisc(size(grating), (length(grating)/2)-2 - 50, (size(grating)+1)/2, 0, [0,1]) - 1); % inner
                        
                        % scale to range of display
                        grating = 255*(grating+1)/2;
                        
                        % create a texture
                        stimulus.tex{iOri, iSF, iPhase, iRing} = mglCreateTexture(grating, [], 1);
                        
                    elseif stimulus.ring(iRing) == 2 % make an anti-ring
                        
                        gratingFull = mglMakeGrating(myscreen.imageWidth, myscreen.imageHeight, stimulus.sf(iSF), stimulus.orientation(iOri), stimulus.phases(iPhase), stimulus.pixRes, stimulus.pixRes);
                        
                        sizeDiff = size(gratingFull) - size(grating);
                        
                        % make an anti-ring
                        rings = (mkDisc(size(grating), (length(grating)/2)-2, (size(grating)+1)/2, 0, [0,1]) + ... % outer
                                 mkDisc(size(grating), (length(grating)/2)-2 - 50, (size(grating)+1)/2, 0, [1,0])); % inner
                             
                        padrings = padarray(rings,[sizeDiff(1)/2, sizeDiff(2)/2],1,'both');
                        
                        grating = gratingFull .*  padrings;
                        
                        % scale to range of display
                        grating = 255*(grating+1)/2;
                        
                        % create a texture
                        stimulus.tex{iOri, iSF, iPhase, iRing} = mglCreateTexture(grating, [], 1);
                        
                    end
                end
            end
        end
        disppercent(iOri/length(stimulus.orientation));
    end
    disppercent(inf);
end
