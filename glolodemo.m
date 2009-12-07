% glolo.m
%
%      usage: [stimulus,task,myscreen]=glolo(stimulus)
%         by: justin gardner
%       date: 01/27/06
%    purpose: glolo program
%             
%       e.g.: stimulus = glolo
%
%             once it is run once, you can reuse the stimulus
%             as long as the psychtoolbox screen has not been
%             closed
%
%             glolo(stimulus);
%
function [stimulus, task, myscreen] = glolodemo

% check arguments
if ~any(nargin == [0 1 2])
  help glolo
  return
end

global debugGrating;
debugGrating = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other screen parameters
myscreen.autoCloseScreen = 1;
myscreen.allowpause = 0;
myscreen.eatkeys = 0;
myscreen.displayname = 'projector';
myscreen.background = 'gray';

myscreen = initScreen(myscreen);

clear global stimulus;
global stimulus;
myscreen = initStimulus('stimulus',myscreen);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up two patches, one for the right and left
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spatial and temporal frequencies
stimulus.patches.sf = 0.8;
stimulus.patches.tf = 6;
stimulus.patches.flickertf = 6;
% how many frames apart each row will be shown for global motion
stimulus.patches.toffset = 6;
% number of frames stimulus will be on. i.e. stimulus
% will have non-zero contrast for this many frames
stimulus.patches.pdur = 25;
% width and height of patches in degrees
stimulus.patches.width = 3;
stimulus.patches.height = 2;
% orientation of patches
stimulus.patches.orientation = 90;
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
stimulus.patches.yspacing = 0.0;
% offset of left and right patches from fixation
stimulus.patches.xoffset = -0.5;
stimulus.patches.yoffset = 0;
% calculate stim period
stimulus.patches.stimPeriod = (stimulus.patches.pdur + stimulus.patches.toffset*(stimulus.patches.numrows-1))/myscreen.framesPerSecond;

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
task{2}{1}.parameterCode.localdir = [0];
task{2}{1}.parameterCode.globaldir = [0];
task{2}{1}.random = 0;
task{2}{1}.seglen = 0;
task{2}{1}.timeInVols = 0;
task{2}{1}.waitForBacktick = 1;

task{2}{2}.parameter.dirnum = [1:4;1:4];
task{2}{2}.parameterCode.localdir =  [-1 1 -1 1];
task{2}{2}.parameterCode.globaldir = [0 0 1 -1];
%task{2}{2}.parameterCode.localdir =  [-1  1];
%task{2}{2}.parameterCode.globaldir = [1 -1];
%task{2}.parameter.dirnum = [1;1];
%task{2}.parameterCode.localdir =  [1]; 
%task{2}.parameterCode.globaldir = [0];
task{2}{2}.random = 0;
% set up segments of trials
task{2}{2}.segmin = [1];
task{2}{2}.segmax = [1];
task{2}{2}.timeInVols = 1;
task{2}{2}.waitForBacktick = 0;
% set up which traces to write out
task{2}{2}.writetrace{1}.tracenum = [1 2];
task{2}{2}.writetrace{1}.tracevar{1} = 'dirnum';
task{2}{2}.writetrace{1}.tracevar{2} = 'dirnum';
task{2}{2}.writetrace{1}.tracerow = [1 2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialze tasks and stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[stimulus task myscreen] = initPatches(stimulus,task,myscreen);

% initialze tasks
for phasenum = 1:length(task{2})
  task{2}{phasenum} = initTask(task{2}{phasenum},myscreen,@startSegmentCallback,@trialStimulusCallback,@trialResponseCallback,[],@endTrialCallback);
end

% set which task is active
pnum = 1;

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
  myscreen = tickScreen(myscreen,task);
end

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
function [stimulus task myscreen] = initPatches(stimulus,task,myscreen)

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
stimulus.patches.flickerContrast = cos(flickerphases);
% and multily by envelope
stimulus.patches.flickerContrast = stimulus.patches.flickerContrast.*stimulus.patches.contrast;

% calculate spacing between patches
xspacer = zeros(stimulus.patches.pixheight,round(MGL.xDeviceToPixels*stimulus.patches.xspacing));
if isodd(size(xspacer,2))
  xspacer = zeros(stimulus.patches.pixheight,round(MGL.xDeviceToPixels*stimulus.patches.xspacing)+1);
end

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
  contrast = stimulus.patches.contrast(contrastNum);
  % get orientation
  angle = d2r(stimulus.patches.orientation);
  % get spatial frequency
  f=stimulus.patches.sf*2*pi; 
  a=cos(angle)*f;
  b=sin(angle)*f;
  % compute grating
  m = eval(sprintf('%s(a*xMesh+b*yMesh+phase).*gratingMask',stimulus.patches.fun));
  % multiply by contrast
  m = contrast*m;
  mrow = m;
  % now chunk into a complete row
  for columnNum = 2:stimulus.patches.numcols
    mrow = [mrow xspacer m];
  end
  % save the stimulus
  if ~debugGrating
    stimulus.patches.stimrow{contrastNum,1} = mglCreateTexture(myscreen.grayIndex+myscreen.inc*mrow);
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
    vertslice = myscreen.grayIndex+myscreen.inc*mrow(:,ceil(size(m,2)/2));
    % display it
    plot(vertslice);
    hold on;hline(128);yaxis(0,255);
    vline(size(mrow,1)/2);
    % set info in title
    title(sprintf('drifting: contrast = %0.4f (%i)',contrast,contrastNum));
    subplot(2,2,3);
    plot(vertslice,getcolor(contrastNum,'-'));
    hold on;hline(128);yaxis(0,255);
    vline(size(mrow,1)/2);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%
  % end debugging display
  %%%%%%%%%%%%%%%%%%%%%%%%
  % make the phase reversing grating
  % compute grating
  m = eval(sprintf('%s(a*xMesh+b*yMesh+stimulus.patches.flickerPhase).*gratingMask',stimulus.patches.fun));
  % multiply by contrast
  m = stimulus.patches.flickerContrast(contrastNum)*m;
  mrow = m;
  % now chunk into a complete row
  for columnNum = 2:stimulus.patches.numcols
    mrow = [mrow xspacer m];
  end
  % store in the second position
  if ~debugGrating
    stimulus.patches.stimrow{contrastNum,2} = mglCreateTexture(myscreen.grayIndex+myscreen.inc*mrow);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%
  % for debugging, displays a cross section of flicker grating
  %%%%%%%%%%%%%%%%%%%%%%%%
  if debugGrating
    subplot(2,2,2);
    cla
    % get a slice of the stimulusn
    vertslice = myscreen.grayIndex+myscreen.inc*mrow(:,ceil(size(m,2)/2));
    % display it
    plot(vertslice);
    hold on;hline(128);yaxis(0,255);
    vline(size(mrow,1)/2);
    % set info in title
    title(sprintf('flickering: contrast = %0.4f',contrast));
    subplot(2,2,4);
    plot(vertslice,getcolor(contrastNum,'-'));
    hold on;hline(128);yaxis(0,255);
    vline(size(mrow,1)/2);
    % and pause, so it doesn't go by too fast
    input('next');
  end
  %%%%%%%%%%%%%%%%%%%%%%%%
  % end debugging display
  %%%%%%%%%%%%%%%%%%%%%%%%
end
disppercent(inf);

if debugGrating,keyboard,end

% get size of row patch
stimulus.patches.rowPixwidth = size(mrow,2);

% set initial values
stimulus.patches.contrastNum = 1;

% set the contrast sequence, make sequences for
% local motion in both directions
stimulus.patches.contrastSeq{1} = [zerocontrast*ones(1,stimulus.patches.toffset*(stimulus.patches.numrows-1)) stimulus.patches.pdur:-1:1];
stimulus.patches.contrastSeq{2} = [zerocontrast*ones(1,stimulus.patches.toffset*(stimulus.patches.numrows-1)) 1:stimulus.patches.pdur];

% set the contrst sequence length
stimulus.patches.contrastSeqLen = length(stimulus.patches.contrastSeq{1});

stimulus.patches.colPixHeight = size(m,1);

yoffset = stimulus.patches.yoffset + ((stimulus.patches.numrows*stimulus.patches.height+stimulus.patches.yspacing*(stimulus.patches.numrows-1))/2);

% set left and right patch rect
for i = 1:stimulus.patches.numrows
  % calculate spacing between patches
  yspacer = (i-1)*(round(stimulus.patches.yspacing)+stimulus.patches.height);
  % set locations for left and right rows
  stimulus.patches.leftdest{i} = [stimulus.patches.xoffset yoffset-yspacer];
  stimulus.patches.rightdest{i} = [-stimulus.patches.xoffset yoffset-yspacer];
end

% blt all the textures to the back buffer once, to get them cued up
for j = 1:size(stimulus.patches.stimrow,2)
  for i = 1:size(stimulus.patches.stimrow,1)
    mglBltTexture(stimulus.patches.stimrow{i,j},[0 0]);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to set the patch local/global dir
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setPatchDir()

global stimulus;
% set starting frame number
stimulus.patches.framenum = 0;

for j = 1:2
  % get the contrast sequence number for this global direction
  % first get default direction (since opposite direction will
  % just reverse this array
  stimulus.patches.contrastSeqNum{j} = 1:stimulus.patches.toffset:stimulus.patches.numrows*stimulus.patches.toffset;
  % if the global direction is opposite, reverse it
  if stimulus.patches.globaldir(j) == -1
    % set the contrast sequence
    stimulus.patches.contrastSeqNum{j} = fliplr(stimulus.patches.contrastSeqNum{j});
  elseif stimulus.patches.globaldir(j) == 0
    % set the sequence numbers to all start in middle
    stimulus.patches.contrastSeqNum{j}(:) = stimulus.patches.toffset*(round(stimulus.patches.numrows/2)-1)+1;
  end
  
  % motion type of 1 means drifting grating motion
  if stimulus.patches.localdir(j)==0
    stimulus.patches.motiontype(j) = 2;
  else
    stimulus.patches.motiontype(j) = 1;
  end
  stimulus.patches.thislocaldir(j) = ceil((-stimulus.patches.localdir(j)+2)/2);
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
  disp(sprintf('local=%i global=%i | local=%i global=%i',task.parameterCode.localdir(task.thistrial.dirnum(1)),task.parameterCode.globaldir(task.thistrial.dirnum(1)),task.parameterCode.localdir(task.thistrial.dirnum(2)),task.parameterCode.globaldir(task.thistrial.dirnum(2))));

  % set local and global directions
  stimulus.patches.localdir = task.parameterCode.localdir(task.thistrial.dirnum);
  stimulus.patches.globaldir = task.parameterCode.globaldir(task.thistrial.dirnum);
else
  stimulus.patches.localdir(:) = 0;
  stimulus.patches.globaldir(:) = 0;
end

setPatchDir;

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to display stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = trialStimulusCallback(task,myscreen)

global stimulus;

% set contrast
stimulus.patches.contrastSeqNum{1} = mod(stimulus.patches.contrastSeqNum{1},stimulus.patches.contrastSeqLen)+1;
stimulus.patches.contrastSeqNum{2} = mod(stimulus.patches.contrastSeqNum{2},stimulus.patches.contrastSeqLen)+1;

if task.thistrial.thisseg == 1
% draw the left stimulus
for rownum = 1:stimulus.patches.numrows
  % this awful looking thing selects the correct row texture and
  % displays it in the right location.
  mglBltTexture(stimulus.patches.stimrow{stimulus.patches.contrastSeq{stimulus.patches.thislocaldir(1)}(stimulus.patches.contrastSeqNum{1}(rownum)),stimulus.patches.motiontype(1)},stimulus.patches.leftdest{rownum},1,-1);
  mglBltTexture(stimulus.patches.stimrow{stimulus.patches.contrastSeq{stimulus.patches.thislocaldir(2)}(stimulus.patches.contrastSeqNum{2}(rownum)),stimulus.patches.motiontype(2)},stimulus.patches.rightdest{rownum},-1,-1);
end
else
  mglClearScreen;
end

% display the icons in corner of screen
% for debugging purposes
if 0
  ldir(1) = task.parameterCode.localdir(task.thistrial.dirnum(1));
  gdir(1) = task.parameterCode.globaldir(task.thistrial.dirnum(1));
  Screen('FillRect',myscreen.w,(myscreen.inc-1)*ldir(1)+myscreen.grayIndex,SetRect(0,0,5,5));
  Screen('FillRect',myscreen.w,(myscreen.inc-1)*(-ldir(1))+myscreen.grayIndex,SetRect(0,5,5,10));
  Screen('FillRect',myscreen.w,(myscreen.inc-1)*gdir(1)+myscreen.grayIndex,SetRect(5,0,10,5));
  Screen('FillRect',myscreen.w,(myscreen.inc-1)*(-gdir(1))+myscreen.grayIndex,SetRect(5,5,10,10));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end trial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = endTrialCallback(task,myscreen)

% this is called at end of trial, can be
% used for updating staircase

