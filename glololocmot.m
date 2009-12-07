% glolo.m
%
%      usage: [myscreen,task,stimulus]=glolo(stimulus)
%         by: justin gardner
%       date: 01/27/06
%    purpose: glolo program (shows multiple directions)
%             
%       e.g.: myscreen = glolo
%
function [myscreen, task, stimulus] = glololocmot(randstate)

% check arguments
if ~any(nargin == [0 1 2])
  help glolo
  return
end

if ieNotDefined('randstate'),randstate = [];end

global debugGrating;
debugGrating = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other screen parameters
myscreen.autoCloseScreen = 1;
myscreen.allowpause = 1;
myscreen.eatkeys = 0;
myscreen.displayname = 'projector';
myscreen.background = 'gray';
myscreen.displayname = 'monitor';
myscreen = initScreen(myscreen,randstate);

clear global stimulus;
global stimulus;
myscreen = initStimulus('stimulus',myscreen);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up two patches, one for the right and left
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spatial and temporal frequencies
stimulus.patches.sf = 0.8;
stimulus.patches.tf = 6;
stimulus.patches.flickertf = 0;%6;
stimulus.patches.flickerMaxContrast = 0.5;
stimulus.patches.maxContrast = 0.5;
% how many frames apart each row will be shown for global motion
stimulus.patches.toffset = 8;%10;%6
% number of frames stimulus will be on. i.e. stimulus
% will have non-zero contrast for this many frames
stimulus.patches.pdur = 25;%25;
% width and height of patches in degrees
stimulus.patches.width = 3;
stimulus.patches.height = 2;
% orientation of patches
stimulus.patches.primaryOrientation = 90;
stimulus.patches.displayOrientations = [0 180];
% sinusoidal funcion used
stimulus.patches.fun = 'sin';
% spatial phase at which the flickering grating will be displayed
stimulus.patches.flickerPhase = 0;
% size of spatial envelope
stimulus.patches.gaussianWidth = stimulus.patches.width/4;
stimulus.patches.gaussianHeight = stimulus.patches.height/4;
% number of rows and columns
stimulus.patches.numrows = 9;
stimulus.patches.numcols = 4;
% spacing between patches in degrees
stimulus.patches.xspacing = 0.0;
stimulus.patches.yspacing = 0.0;%0.5;
% offset of left and right patches from fixation
stimulus.patches.xoffset = -2;
stimulus.patches.yoffset = 0;%1.5;

[stimulus myscreen] = initPatches(stimulus,myscreen);

% set the first task to be the fixation staircase task
[task{1} myscreen] = fixStairInitTask(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up baseline task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
triallen = 15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up patches stimului parameters
task{2}{1}.numBlocks = 1;
task{2}{1}.parameter.dirnum = [1;1];
task{2}{1}.parameterCode.localdir = -1;
task{2}{1}.parameterCode.globaldir = inf;
task{2}{1}.parameterCode.orientation = inf;
task{2}{1}.random = 0;
task{2}{1}.seglen = [0.1 11.9];
task{2}{1}.timeInVols = 0;
task{2}{1}.waitForBacktick = 1;

task{2}{2}.parameter.dirnum = [1:6;1:6];
task{2}{2}.parameterCode.localdir =  [-1 -1 -1 1 1 1];
task{2}{2}.parameterCode.globaldir = [0 -1 1 0 -1 1];
task{2}{2}.parameterCode.orientation = ones(1,max(task{2}{2}.parameter.dirnum(:)));
task{2}{2}.random = 1;
% set up segments of trials
%rampDownTime = stimulus.patches.pdur/myscreen.framesPerSecond;
%task{2}{2}.segmin = [4.8 3.5];
%task{2}{2}.segmax = [4.8 9.5];
task{2}{2}.segmin = [12 5.9];
task{2}{2}.segmax = [12 5.9];
task{2}{2}.synchToVol = [0 1];
task{2}{2}.timeInVols = 0;
task{2}{2}.waitForBacktick = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialze tasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialze tasks
for phasenum = 1:length(task{2})
  [task{2}{phasenum} myscreen] = initTask(task{2}{phasenum},myscreen,@startSegmentCallback,@trialStimulusCallback,@trialResponseCallback,[],@endTrialCallback);
end

% set which task is active
pnum = 1;

% start eating all keys that myscreen will need
mglEatKeys(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the tasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set which phase is active
tnum = 1;

while (tnum <= length(task{2})) && ~myscreen.userHitEsc
  % updatethe task
  [task{2} myscreen tnum] = updateTask(task{2},myscreen,tnum);
  % display fixation cross
  [task{1} myscreen] = updateTask(task{1},myscreen,1);
  % flip screen
  [myscreen task] = tickScreen(myscreen,task);
end

% stop eating keys
mglEatKeys([]);

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     STIMULUS Gratings start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init gratings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stimulus myscreen] = initPatches(stimulus,myscreen)

global MGL;
global debugGrating;

% get the dimensions
stimulus.patches.pixwidth = round(stimulus.patches.width*MGL.xDeviceToPixels);
stimulus.patches.pixheight = round(stimulus.patches.height*MGL.yDeviceToPixels);

% make sure pixwidth is odd
if iseven(stimulus.patches.pixwidth)
  stimulus.patches.pixwidth = stimulus.patches.pixwidth-1;
end

% get grid over which we should calculate
xhalfwidth = round(stimulus.patches.width/2);
yhalfwidth = round(stimulus.patches.height/2);

% get x and y
x = -xhalfwidth:2*xhalfwidth/(stimulus.patches.pixwidth-1):xhalfwidth;
y = -yhalfwidth:2*yhalfwidth/(stimulus.patches.pixheight-1):yhalfwidth;
% turn into grid
[xMesh,yMesh] = meshgrid(x,y);

% now compute the mask
gratingMask = ((xMesh.^2)/xhalfwidth^2 + (yMesh.^2)/yhalfwidth^2) < 1;

% create gaussian window
if isfield(stimulus.patches,'gaussianWidth') && isfield(stimulus.patches,'gaussianHeight')
  gratingMask = exp(-((xMesh.^2)/stimulus.patches.gaussianWidth^2 + (yMesh.^2)/stimulus.patches.gaussianHeight^2));
end

% this is the old way of computing phases, when
% we wanted phase to be independent of contrast
% now we always show the same phases during a
% pdur length stimulus presentations
if 0
  % compute which phases we need
  % find minimum phase difference from frame to frame
  minphase = stimulus.patches.tf/myscreen.framesPerSecond;
  % now find what multiple of that phase will give you 
  % an exact integer...or give up after n multiples
  % this will give us all the correct phases needed
  % for displaying the stimulus starting at cosine
  % phase.
  maxphases = 100;
  phasesneeded = maxphases;
  for i = 1:maxphases
    if nearlyequal(round(i*minphase),i*minphase,10)
      phasesneeded = i;
      break;
    end
  end
  % now calculate the phases we need
  for i = 0:phasesneeded-1
    stimulus.patches.phase(i+1) = mod(i*minphase*2*pi,2*pi);
  end
end

% compute what contrasts we need, make sure that we 
% have pdur number of non zero contrasts
stimulus.patches.contrast = (sin((0:2*pi/((stimulus.patches.pdur+2)-1):2*pi)+3*pi/2)+1)/2;
% remove starting and ending 0 contrast
stimulus.patches.contrast = stimulus.patches.contrast(2:end-1);

%new way of computing phases, we only need the the
% phases that will be shown during each pdur. Making
% sure to get 0 phase (i.e. sin phase if we have
% set stimulus.patches.fun to sin) in the middle of pdur
% first get the center of the pdur
midpdur = round(length(stimulus.patches.contrast)/2);
% now just compute phase in phase increments away from
% the center
phasestep = 2*pi*stimulus.patches.tf/myscreen.framesPerSecond;
stimulus.patches.phase = phasestep*((1:length(stimulus.patches.contrast))-midpdur);

% now compute flicker contrast
% get necessary temporal frequency advance for flickertf
phasestep = 2*pi*stimulus.patches.flickertf/myscreen.framesPerSecond;
flickerphases = phasestep*((1:length(stimulus.patches.contrast))-midpdur);
% we want to insure maximum flicker contrast at center of pdur
% so use cos phase
stimulus.patches.flickerContrast = stimulus.patches.flickerMaxContrast * cos(flickerphases);
% and multily by envelope
stimulus.patches.flickerContrast = stimulus.patches.flickerContrast.*stimulus.patches.contrast;

% now add one contrast and phase for 0 contrast
stimulus.patches.flickerContrast(stimulus.patches.pdur+1) = 0;
stimulus.patches.contrast(stimulus.patches.pdur+1) = 0;
stimulus.patches.phase(stimulus.patches.pdur+1) = 0;
zerocontrast = stimulus.patches.pdur+1;

% compute each asked for grating
disppercent(-inf,'Generating gratings');
for contrastNum = 1:length(stimulus.patches.contrast)
  disppercent(contrastNum/length(stimulus.patches.contrast));
  % get this phase and contrast
  phase = stimulus.patches.phase(contrastNum);
  contrast = stimulus.patches.maxContrast * stimulus.patches.contrast(contrastNum);
  % get orientation
  angle = d2r(stimulus.patches.primaryOrientation);
  % get spatial frequency
  f=stimulus.patches.sf*2*pi; 
  a=cos(angle)*f;
  b=sin(angle)*f;
  % compute grating
  m = eval(sprintf('%s(a*xMesh+b*yMesh+phase).*gratingMask',stimulus.patches.fun));
  % multiply by contrast
  m = contrast*m;
  % save the stimulus
  if ~debugGrating
    stimulus.patches.stim{contrastNum,1} = mglCreateTexture(myscreen.grayIndex+myscreen.inc*m);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%
  % for debugging, displays a cross section of drifting grating
  %%%%%%%%%%%%%%%%%%%%%%%%
  if debugGrating
    if debugGrating == 1
      mglClose;
      clf;
      debugGrating = 2;
    end
    subplot(2,2,1);
    cla
    % get a slice of the stimulus
    vertslice = myscreen.grayIndex+myscreen.inc*m(:,ceil(size(m,2)/2));
    % display it
    plot(vertslice);
    hold on;hline(128);yaxis(0,255);
    vline(size(m,1)/2);
    % set info in title
    title(sprintf('drifting: contrast = %0.4f (%i)',contrast,contrastNum));
    subplot(2,2,3);
    plot(vertslice,getcolor(contrastNum,'-'));
    hold on;hline(128);yaxis(0,255);
    vline(size(m,1)/2);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%
  % end debugging display
  %%%%%%%%%%%%%%%%%%%%%%%%
  % make the phase reversing grating
  % compute grating
  m = eval(sprintf('%s(a*xMesh+b*yMesh+stimulus.patches.flickerPhase).*gratingMask',stimulus.patches.fun));
  % multiply by contrast
  m = stimulus.patches.flickerContrast(contrastNum)*m;
  % store in the second position
  if ~debugGrating
    stimulus.patches.stim{contrastNum,2} = mglCreateTexture(myscreen.grayIndex+myscreen.inc*m);
    stimulus.patches.stim{contrastNum,3} = mglCreateTexture(myscreen.grayIndex-myscreen.inc*m);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%
  % for debugging, displays a cross section of flicker grating
  %%%%%%%%%%%%%%%%%%%%%%%%
  if debugGrating
    subplot(2,2,2);
    cla
    % get a slice of the stimulusn
    vertslice = myscreen.grayIndex+myscreen.inc*m(:,ceil(size(m,2)/2));
    % display it
    plot(vertslice);
    hold on;hline(128);yaxis(0,255);
    vline(size(m,1)/2);
    % set info in title
    title(sprintf('flickering: contrast = %0.4f',contrast));
    subplot(2,2,4);
    plot(vertslice,getcolor(contrastNum,'-'));
    hold on;hline(128);yaxis(0,255);
    vline(size(m,1)/2);
    % and pause, so it doesn't go by too fast
    input('next');
  end
  %%%%%%%%%%%%%%%%%%%%%%%%
  % end debugging display
  %%%%%%%%%%%%%%%%%%%%%%%%
end
disppercent(inf);

if debugGrating,keyboard,end

% set initial values
stimulus.patches.contrastNum = 1;

% set the contrast sequence, make sequences for
% local motion in both directions
stimulus.patches.wrap = 1;
if stimulus.patches.wrap
  zeroContrastSequence = zerocontrast*ones(1,stimulus.patches.toffset*(stimulus.patches.numrows)-stimulus.patches.pdur);
  stimulus.patches.contrastSeq{1} = [zeroContrastSequence stimulus.patches.pdur:-1:1];
  stimulus.patches.contrastSeq{2} = [zeroContrastSequence 1:stimulus.patches.pdur];
else
  zeroContrastSequence = zerocontrast*ones(1,stimulus.patches.toffset*(stimulus.patches.numrows-1));
  stimulus.patches.contrastSeq{1} = [zeroContrastSequence stimulus.patches.pdur:-1:1];
  stimulus.patches.contrastSeq{2} = [zeroContrastSequence 1:stimulus.patches.pdur];
end

% get the last sequence number that is gray
stimulus.patches.lastZeroContrastSeqNum = length(zeroContrastSequence);

% set the contrast sequence length
stimulus.patches.contrastSeqLen = length(stimulus.patches.contrastSeq{1});
stimulus.patches.stimPeriod = stimulus.patches.contrastSeqLen/myscreen.framesPerSecond;

stimulus.patches.colPixHeight = size(m,1);

yoffset = stimulus.patches.yoffset + ((stimulus.patches.numrows*stimulus.patches.height+stimulus.patches.yspacing*(stimulus.patches.numrows-1))/2);

% set left and right patch rect
stimulus.patches.rowcolsize = [stimulus.patches.numcols stimulus.patches.numrows];
for rownum = 1:stimulus.patches.numrows
  for colnum = 1:stimulus.patches.numcols
    % calculate spacing between patches
    xspacer = (colnum-1)*(round(stimulus.patches.xspacing)+stimulus.patches.width);
    yspacer = (rownum-1)*(round(stimulus.patches.yspacing)+stimulus.patches.height);
    % set locations for left and right rows
    stimulus.patches.leftdest(sub2ind(stimulus.patches.rowcolsize,colnum,rownum),:) = [stimulus.patches.xoffset-xspacer yoffset-yspacer];
    stimulus.patches.rightdest(sub2ind(stimulus.patches.rowcolsize,colnum,rownum),:) = [-stimulus.patches.xoffset+xspacer yoffset-yspacer];
  end
end

% blt all the textures to the back buffer once, to get them cued up
for j = 1:size(stimulus.patches.stim,2)
  for i = 1:size(stimulus.patches.stim,1)
    mglBltTexture(stimulus.patches.stim{i,j},[0 0]);
  end
end

% contrast sequence number defaults to gray forever
stimulus.patches.contrastSeqNum(1:2,1:stimulus.patches.numcols,1:stimulus.patches.numrows) = -inf;
stimulus.patches.desiredContrastSeqNum = stimulus.patches.contrastSeqNum;
stimulus.patches.thislocaldir(1:2,1:stimulus.patches.numcols,1:stimulus.patches.numrows) = 1;
stimulus.patches.desiredLocaldir(1:2,1:stimulus.patches.numcols,1:stimulus.patches.numrows) = 1;
stimulus.patches.orientation = [inf inf];
stimulus.patches.motiontype(1:2,1:stimulus.patches.numcols,1:stimulus.patches.numrows) = 1;
stimulus.patches.desiredMotiontype = stimulus.patches.motiontype;
stimulus.patches.size = size(stimulus.patches.contrastSeqNum);
stimulus.patches.motiontypeCountdown = -inf;
stimulus.patches.orientations(1:2,1:stimulus.patches.numcols,1:stimulus.patches.numrows) = stimulus.patches.displayOrientations(1);
stimulus.patches.desiredOrientations(1:2,1:stimulus.patches.numcols,1:stimulus.patches.numrows) = stimulus.patches.displayOrientations(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to set the patch local/global dir
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function myscreen = setPatchDir(myscreen)

global stimulus;
% set starting frame number
stimulus.patches.framenum = 0;

for side = 1:2
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % first set up the global motion
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % get the contrast sequence number for this global direction
  % first get default direction (since opposite direction will
  % just reverse this array

  % first compute a random offset so that the global motion won't start
  % in the same place. 
  randOffset = floor(rand*stimulus.patches.numrows)*stimulus.patches.toffset;

  % also find out the current phase of the stimulus so that we can seemlessly
  % switch to the new stimulus w/out jumping phase
  toffsetPhase = min(stimulus.patches.contrastSeqNum(side,:))-1;
  if isinf(toffsetPhase),toffsetPhase = 0;end

  % now go through and set the contrast seq numbers to create a global pattern
  for colnum = 1:stimulus.patches.numcols
    % set global pattern
%    stimulus.patches.desiredContrastSeqNum(side,colnum,:) = (1:stimulus.patches.toffset:stimulus.patches.numrows*stimulus.patches.toffset)+toffsetPhase;
    stimulus.patches.desiredContrastSeqNum(side,colnum,:) = (1:stimulus.patches.toffset:stimulus.patches.numrows*stimulus.patches.toffset)+toffsetPhase+stimulus.patches.toffset*floor(rand*stimulus.patches.numrows);
    % add add the random offset so that it doesn't start at the same place
    stimulus.patches.desiredContrastSeqNum(side,colnum,:) = mod(stimulus.patches.desiredContrastSeqNum(side,colnum,:)+randOffset-1,length(stimulus.patches.contrastSeq{side}))+1;
  end

  % if the global direction is opposite, reverse it
  if stimulus.patches.globaldir(side) == -1
    % set the contrast sequence
    stimulus.patches.desiredContrastSeqNum(side,:,:) = flipdim(stimulus.patches.desiredContrastSeqNum(side,:,:),3);

  % if gloabl direction is set to 0, then set the patches to all come on at
  % the same time (old glolo)
  elseif stimulus.patches.globaldir(side) == 0
    % set the sequence numbers to all start in middle
    for colnum = 1:stimulus.patches.numcols
      stimulus.patches.desiredContrastSeqNum(side,colnum,:) = stimulus.patches.toffset*(round(stimulus.patches.numrows/2)-1)+1;
    end
  
  % if global direction is set to infinity that means to randomize when targets come
  % on and off, do that by randomly permuting the locations
  elseif stimulus.patches.globaldir(side) == inf
    for colnum = 1:stimulus.patches.numcols
      randNums(colnum,:) = randperm(stimulus.patches.numrows)-1;
    end
    % now randomize across cols
    for rownum = 1:stimulus.patches.numrows
      randNums(:,rownum) = randNums(randperm(stimulus.patches.numcols),rownum);
    end
    % set the sequence numbers to be random
    stimulus.patches.desiredContrastSeqNum(side,:,:) = randNums*(stimulus.patches.toffset)+toffsetPhase+1;
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % now deal with the orientations
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if isinf(stimulus.patches.orientation(side))
    stimulus.patches.desiredOrientations(side,1:stimulus.patches.numcols,1:stimulus.patches.numrows) = stimulus.patches.displayOrientations(ceil(rand(stimulus.patches.numcols,stimulus.patches.numrows)*length(stimulus.patches.displayOrientations)));
  else
    % set all the patches to the orientation called for
    stimulus.patches.desiredOrientations(side,1:stimulus.patches.numcols,1:stimulus.patches.numrows) = ones(stimulus.patches.numcols,stimulus.patches.numrows)*stimulus.patches.orientation(side);
  end  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % now deal with the local motion
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if stimulus.patches.localdir(side)==0
    % randomly add 1, since this is the static case, and the randomization
    % will choose which phase the flickering or static grating will be in
%    stimulus.patches.desiredMotiontype(side,:,:) = 2 + (rand(stimulus.patches.numrows,stimulus.patches.numcols)>0.5);
    stimulus.patches.desiredMotiontype(side,:) = stimulus.patches.motiontype(side,:);
    stimulus.patches.desiredMotiontype(side,stimulus.patches.desiredMotiontype(side,:) == 0) = 2;
    stimulus.patches.desiredMotiontype(side,stimulus.patches.desiredMotiontype(side,:) == 1) = 2;
  else
    % motiontype 1 is for a moving grating
    % motiontype 2 and 3 are for flickering grating
    stimulus.patches.desiredMotiontype(side,:) = 1;
  end
  stimulus.patches.desiredLocaldir(side,:) = ceil((-stimulus.patches.localdir(side)+2)/2);
end

% don't change local motion right away, wait 
stimulus.patches.motiontypeCountdown = stimulus.patches.pdur;

% offset the contrast sequence nat the very beginning now, when
% all the patches start off with infinity.
if all(isinf(stimulus.patches.contrastSeqNum))
  stimulus.patches.contrastSeqNum = stimulus.patches.desiredContrastSeqNum-stimulus.patches.pdur+stimulus.patches.toffset;
  stimulus.patches.desiredContrastSeqNum = stimulus.patches.contrastSeqNum;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     STIMULUS Gratings end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to handle observer response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = trialResponseCallback(task,myscreen)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to set stimulus parameters at
% the beginning of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task,myscreen)

global stimulus;

if (task.thistrial.thisseg == 1)
  % set local and global directions
  stimulus.patches.localdir = task.parameterCode.localdir(task.thistrial.dirnum);
  stimulus.patches.globaldir = task.parameterCode.globaldir(task.thistrial.dirnum);
  stimulus.patches.orientation = task.parameterCode.orientation(task.thistrial.dirnum);
  % have both do the same thing
  if 0
    stimulus.patches.localdir(2) = stimulus.patches.localdir(1);
    stimulus.patches.globaldir(2) = stimulus.patches.globaldir(1);
    stimulus.patches.orientation(2) = stimulus.patches.orientation(1);
  end
  % display settings
  disp(sprintf('%i:orient=%i local=%i global=%i | %i:orient=%i local=%i global=%i',task.thistrial.dirnum(1),stimulus.patches.orientations(1),stimulus.patches.localdir(1),stimulus.patches.globaldir(1),task.thistrial.dirnum(2),stimulus.patches.orientation(2),stimulus.patches.localdir(2),stimulus.patches.globaldir(2)));
  % actually set the patch directions
  myscreen = setPatchDir(myscreen);
elseif (task.thistrial.thisseg == 2)
  stimulus.patches.localdir(:) = 1;
  stimulus.patches.globaldir(:) = inf;
  stimulus.patches.orientation(:) = inf;
  myscreen = setPatchDir(myscreen);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to display stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = trialStimulusCallback(task,myscreen)

global stimulus;

% update contrast Seqnum
% increment patch contrast sequence number
stimulus.patches.contrastSeqNum(stimulus.patches.contrastSeqNum>0) = ...
    mod(stimulus.patches.contrastSeqNum(stimulus.patches.contrastSeqNum>0),stimulus.patches.contrastSeqLen)+1;
stimulus.patches.contrastSeqNum(stimulus.patches.contrastSeqNum<=0) = ...
    stimulus.patches.contrastSeqNum(stimulus.patches.contrastSeqNum<=0)+1;

%increment desired patch contrast sequence number
stimulus.patches.desiredContrastSeqNum(stimulus.patches.desiredContrastSeqNum>0) = ...
    mod(stimulus.patches.desiredContrastSeqNum(stimulus.patches.desiredContrastSeqNum>0),stimulus.patches.contrastSeqLen)+1;
stimulus.patches.desiredContrastSeqNum(stimulus.patches.desiredContrastSeqNum<=0) = ...
    stimulus.patches.desiredContrastSeqNum(stimulus.patches.desiredContrastSeqNum<=0)+1;

% check for any patches that are gray
grayPatches = find(stimulus.patches.contrastSeqNum(:) <= stimulus.patches.lastZeroContrastSeqNum);
desiredGrayPatches = find(stimulus.patches.desiredContrastSeqNum(:) <= stimulus.patches.lastZeroContrastSeqNum);

% get the gray values these point to
grayVals = stimulus.patches.contrastSeqNum(grayPatches);
desiredGrayVals = stimulus.patches.desiredContrastSeqNum(desiredGrayPatches);

% find all the static patches that are gray and randomly flip their phase
grayStaticPatches = find(((stimulus.patches.motiontype(:) == 2) | (stimulus.patches.motiontype(:) == 3)) & (stimulus.patches.contrastSeqNum(:) <= stimulus.patches.lastZeroContrastSeqNum));
stimulus.patches.motiontype(grayStaticPatches) = 2 + (rand(1,length(grayStaticPatches))>0.5);

% now see if there is a mismatch, between what we have and what we want
if any(stimulus.patches.contrastSeqNum(desiredGrayPatches) ~= desiredGrayVals)
  % get the patches we can change, that is, ones that are gray in both the actual
  % being displayed sequence and the desired sequence
  desiredChangePatches = intersect(grayPatches,desiredGrayPatches);
  [desiredChangeVals sortindex] = sort(stimulus.patches.desiredContrastSeqNum(desiredChangePatches));
  desiredChangePatches = desiredChangePatches(sortindex);
    
  % and now get the patches we can really change, that is, we look for
  % the desired phase, and if it is available then we assign it to 
  % a change patch and all the phases that are assigned to change patches
  % get randomly assigned to the remaining patches (that way we are
  % not changing the total number of phases present but just their
  % distribution on the screen).
  sortGrayVals = sort(grayVals);
  i = 1;nonChangeVals = [];changePatches = [];changeVals = [];
  for j = 1:length(sortGrayVals)
    % if there is a match
    if (i <= length(desiredChangeVals)) && (desiredChangeVals(i) == sortGrayVals(j))
      changePatches(end+1) = desiredChangePatches(i);
      changeVals(end+1) = desiredChangeVals(i);
      i = i+1;
    else
      nonChangeVals(end+1) = sortGrayVals(j);
    end
  end
  nonChangePatches = setdiff(grayPatches,changePatches);
  nonChangeVals = nonChangeVals(randperm(length(nonChangeVals)));

  % and reset to make the changePatches match the desired phase
  if ~isempty(changePatches)
    stimulus.patches.contrastSeqNum(changePatches) = changeVals;
  end
  if ~isempty(nonChangePatches)
    stimulus.patches.contrastSeqNum(nonChangePatches) = nonChangeVals;
  end
end

% now change the motion type of any gray patches
stimulus.patches.motiontypeCountdown = stimulus.patches.motiontypeCountdown-1;
if stimulus.patches.motiontypeCountdown <= 0
  % if you are gray then the orientation can local direction is
  % allowed to change to the desired orienation and localdir
  stimulus.patches.orientations(grayPatches) = stimulus.patches.desiredOrientations(grayPatches);
  stimulus.patches.thislocaldir(grayPatches) = stimulus.patches.desiredLocaldir(grayPatches);
end
% for interstimulus interval, we force the contrastSeqNum to stop
% incrementing at 1, so that all stimuli go to gray
%if any(task.thistrial.thisseg==[1 3])
%  stimulus.patches.contrastSeqNum(stimulus.patches.contrastSeqNum<=stimulus.patches.lastZeroContrastSeqNum)=1;
%end


% randomize the orientation of any gray patches that also have their
% orientation set to inf, assume for now that this is not done
% independently on left and right, i.e. they will always be
% at a random orientaiton on both sides
if isinf(stimulus.patches.orientation(1))
  stimulus.patches.orientations(grayPatches) = stimulus.patches.displayOrientations(ceil(rand(1,length(grayPatches))*length(stimulus.patches.displayOrientations)));
end

% get this contrast seq, so we can make some that negative values
% mean to stay at gray
thisContrastSeqNum = stimulus.patches.contrastSeqNum;
thisContrastSeqNum(thisContrastSeqNum<=0) = 1;

% grab the appropriate stimuli
for rownum = 1:stimulus.patches.numrows
  for colnum = 1:stimulus.patches.numcols
    % get the correct patch
    leftPatches(sub2ind(stimulus.patches.rowcolsize,colnum,rownum)) = stimulus.patches.stim{stimulus.patches.contrastSeq{stimulus.patches.thislocaldir(1,colnum,rownum)}(thisContrastSeqNum(1,colnum,rownum)),stimulus.patches.motiontype(1,colnum,rownum)};
    rightPatches(sub2ind(stimulus.patches.rowcolsize,colnum,rownum)) = stimulus.patches.stim{stimulus.patches.contrastSeq{stimulus.patches.thislocaldir(2,colnum,rownum)}(thisContrastSeqNum(2,colnum,rownum)),stimulus.patches.motiontype(2,colnum,rownum)};
  end
end

% and blit them all together to the screen
%mglBltTexture(leftPatches,stimulus.patches.leftdest,1,-1);
%mglBltTexture(rightPatches,stimulus.patches.rightdest,-1,-1);
mglBltTexture(leftPatches,stimulus.patches.leftdest,1,-1,stimulus.patches.orientations(1,:,:));
mglBltTexture(rightPatches,stimulus.patches.rightdest,-1,-1,stimulus.patches.orientations(2,:,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end trial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = endTrialCallback(task,myscreen)

% this is called at end of trial, can be
% used for updating staircase
