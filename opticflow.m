% opticflow.m
%
%      usage: opticflow
%        $Id: opticflow.m,v 1.12 2008/03/16 15:42:36 eli Exp $
%         by: justin gardner
%       date: 01/17/07
%    purpose: to display optic flow dots
%      usage: function myscreen = opticflow(runNum,posNum)
%       
%
function myscreen = opticflow(runNum,posNum)

% check arguments
if ~any(nargin == [2])
    help opticflow
    return
end

% initalize the screen
myscreen.autoCloseScreen = 0;
myscreen.allowpause = 1;
myscreen.displayname = 'projector';
myscreen.background = 'gray';
myscreen = initScreen(myscreen);
myscreen.waitForBacktick =1;

% init stimulus
clear global stimulus;
global stimulus;
stimulus.seqPos = 0;
stimulus.runNum = runNum;
myscreen = initStimulus('stimulus',myscreen);

% to initialize the stimulus for your experiment.
stimulus = myInitStimulus(stimulus,myscreen);

global fixStimulus
fixStimulus.diskSize = .5;
switch(posNum)
  case 1
    %fixStimulus.pos = [-2*stimulus.size -stimulus.size];
     fixStimulus.pos = [-6 -6-stimulus.yoffset]; 
  case 2
    %fixStimulus.pos = [-2*stimulus.size stimulus.size];
    fixStimulus.pos = [-6 6-stimulus.yoffset]; 
  case 3
    %fixStimulus.pos = [2*stimulus.size stimulus.size];
     fixStimulus.pos = [6 6-stimulus.yoffset];
  case 4
    %fixStimulus.pos = [2*stimulus.size -stimulus.size];
    fixStimulus.pos = [6 -6-stimulus.yoffset];
  case 5
    fixStimulus.pos = [0 0-stimulus.yoffset];

end
[task{1} myscreen] = fixStairInitTask(myscreen);

task{2}{1}.waitForBacktick = 1;
% fix: the task defined here has two segments, one that
% is 3 seconds long followed by another that is 
% 6-9 seconds (randomized in steps of 1.5 seconds)
% change this to what you want for your trial
task{2}{1}.segmin = [0.4 0.4 0.4 0.3];
task{2}{1}.segmax = [0.4 0.4 0.4 0.3];
task{2}{1}.synchToVol = [0 0 0 1];

% fix: enter the parameter of your choice
task{2}{1}.parameter.myParameter = [0 30 90];
task{2}{1}.random = 1;
% and set to save the parameter setting in a trace
% so that you can later decode what trial was presented
% when.
% fix: You will need to change this to save out the
% parameters that you are actually want to keep track of.
% task{2}.writeTrace{2}.tracenum = 1;
% task{2}.writeTrace{2}.tracevar{2} = 'myParameter';
% task{2}.writeTrace{2}.usenum = 1;

% initialize the task
for phaseNum = 1:length(task{2})
    task{2}{phaseNum} = initTask(task{2}{phaseNum},myscreen,@stimStartSegmentCallback,@stimDrawStimulusCallback);
end


% fix: you will change the funciton myInitStimulus
% to initialize the stimulus for your experiment.
stimulus = myInitStimulus(stimulus,myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);
myscreen.background = 0;
mglClearScreen(0);mglFlush;mglClearScreen(0);mglFlush;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task)) && ~myscreen.userHitEsc
    % update the task
    [task{2} myscreen tnum] = updateTask(task{2}, myscreen, phaseNum);
    % display fixation cross
    [task{1} myscreen] = updateTask(task{1},myscreen,1);
    % flip screen
    myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = stimStartSegmentCallback(task, myscreen)

global stimulus;
%newSpeed = [0 0 rand(1)*0.5-0.25];

% set which mseq index we are on
if task.thistrial.thisseg == 1
    stimulus.seqPos = stimulus.seqPos+1;
    stimulus.drawdots = stimulus.seq(stimulus.seqPos,:);
end

for i = 1:stimulus.n
    % translation and rotation matrices
    %stimulus.dots{i}.T = [rand(1)*0.5-0.25 rand(1)*0.5-0.25 rand(1)*0.5-0.25];
    stimulus.dots{i}.T = [rand(1)*0.01-0.005 rand(1)*0.01-0.005 (round(rand(1))*2-1)*(rand(1)+.2)];
    stimulus.dots{i}.T = 0.3*stimulus.dots{i}.T./norm(stimulus.dots{i}.T);
    stimulus.dots{i}.Tarray = repmat(stimulus.dots{i}.T,stimulus.dots{i}.n,1)';
    %    stimulus.dots{i}.T = [0 0 rand(1)*1.5-0.75];
    %stimulus.dots{i}.T = newSpeed;
    stimulus.dots{i}.rotmatrix = makerotmatrix(rand(1)*4-2);
    % pick a random direction
    stimulus.dots{i}.direction = rand*360;
end

% read from the msequence
% stimulus.drawdots = round(rand(stimulus.n,1));

% fix: do anything that needs to be done at the beginning
% of a segment (like for example setting the stimulus correctly
% according to the parameters etc).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = stimDrawStimulusCallback(task, myscreen)

global stimulus

% fix: display your stimulus here, for this code we just display 
% a fixation cross that changes color depending on the segment
% we are on.

mglClearScreen;
for i = 1:stimulus.n
    if stimulus.drawdots(i)==1
        mglPoints2(stimulus.dots{i}.x,stimulus.dots{i}.y,3);
        stimulus.dots{i} = updateDotsLinear(stimulus.dots{i},1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen)

stimulus.init = 1;
stimulus.size = 2.5;
stimulus.radius = stimulus.size/2;
stimulus.spacing = .5;
stimulus.nx = 10;
stimulus.ny = 6;
stimulus.yoffset = 2;

stimulus.xcenter = [-stimulus.spacing/2-stimulus.radius-(stimulus.nx/2-1)*(stimulus.size+stimulus.spacing):stimulus.size+stimulus.spacing:stimulus.spacing/2+stimulus.radius+(stimulus.nx/2-1)*(stimulus.size+stimulus.spacing)];
stimulus.ycenter = -stimulus.size+[-stimulus.spacing/2-stimulus.radius-(stimulus.ny/2-1)*(stimulus.size+stimulus.spacing):stimulus.size+stimulus.spacing:stimulus.spacing/2+stimulus.radius+(stimulus.ny/2-1)*(stimulus.size+stimulus.spacing)] - stimulus.yoffset;

[stimulus.xcenter stimulus.ycenter] = meshgrid(stimulus.xcenter,stimulus.ycenter);
stimulus.xcenter = stimulus.xcenter(:);
stimulus.ycenter = stimulus.ycenter(:);
stimulus.n = length(stimulus.xcenter);
for i = 1:stimulus.n
    stimulus.dots{i} = initDotsLinear(myscreen,stimulus.xcenter(i),stimulus.ycenter(i),stimulus.radius);
end


%stimulus.p = 6*8;
stimulus.useMseq=0;
if stimulus.useMseq==1                  % calculate the appropriate msequence
    stimulus.whichSeq = 1;
    stimulus.seq = [];
    stimulus.Delta = 10;
    stimulus.ms = mseq2(7,4,0,stimulus.whichSeq);
    stimulus.ms(find(stimulus.ms~=1))=0;
    shift = 0;
    disp(sprintf('(opticflow) Using an msequence stimulus: numPixels=%i, Delta=%i', stimulus.n, stimulus.Delta));
    for i=1:stimulus.n
        stimulus.seq(:,i) = circshift(stimulus.ms,shift);
        shift =  shift + stimulus.Delta;
        fprintf('.');
    end
    disp(sprintf(' '));

else                                  % use a random stimulus sequence
    if stimulus.runNum<3
        stimulus.prob=.1;   % 0.03125;  % .125  .0625    0.03125
        disp(sprintf('(opticflow) Using random stimulus from a uniform distribution: prob=%i', stimulus.prob));
        stimulus.seqLength=256*2;
        rand('state', sum([1974 12 01 6.30 1 1]));
        stimulus.seq = (stimulus.prob>rand(stimulus.seqLength,stimulus.n));
        stimulus.overlap=10;
        stimulus.breakpoints=[0 256];
        stimulus.seq = circshift(stimulus.seq, stimulus.overlap-stimulus.breakpoints(stimulus.runNum));
    else % localizer
        blockLength=14;
        foo = zeros(stimulus.n/2, blockLength);
        foo(1:(stimulus.n/2), 1:(blockLength/2)) = 1;
        foo((stimulus.n/2)+1:stimulus.n, (blockLength/2)+1:blockLength) = 1;
        foo = repmat(foo, 1, 11);
        %foo([19 20 23 24 25 26 29 30],:) = 0;
        foo([25 26 29 30 31 32 35 36],:) = 0;
        stimulus.seq=foo';
    end
end

    %stimulus.seq = circshift(stimulus.seq, stimulus.Delta -  (length(stimulus.seq)/3)stimulus.runNum;
    % for i=1:24
    %  plot(flipud(ifft(fft(stimulus.seq(:,5)) .* conj(fft(stimulus.seq(:,i)))))/length(stimulus.seq)); hold on
    % end

    stimulus.updateFreq = 0.5;
    stimulus.seqPos = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initDotsLinear(myscreen,xcenter,ycenter,radius)

% dot parameters
dots.dotdensity = 6;
dots.speed = 10;
dots.direction = 45;
dots.xcenter = xcenter;
dots.ycenter = ycenter;
dots.radius = radius;

% set type
dots.type = 'linear';

% get the number of dots
dots.n = round(dots.dotdensity*((2*dots.radius)^2));

% get initial position
dots.x = rand(1,dots.n)*dots.radius*2+dots.xcenter-dots.radius;
dots.y = rand(1,dots.n)*dots.radius*2+dots.ycenter-dots.radius;

% get borders
dots.max = [dots.xcenter+dots.radius dots.ycenter+dots.radius];
dots.min = [dots.xcenter-dots.radius dots.ycenter-dots.radius];
dots.span = [2*dots.radius 2*dots.radius];

dots = updateDotsLinear(dots,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to update the dots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updateDotsLinear(dots,coherence)

global MGL;
% get the step size
dots.stepsize = [cos(d2r(dots.direction)) sin(d2r(dots.direction))]*dots.speed/MGL.frameRate;

% update position of dots
dots.x = dots.x+dots.stepsize(1);
dots.y = dots.y+dots.stepsize(2);

% remap offscreen dots
offscreen = dots.x > dots.max(1);
dots.x(offscreen) = dots.x(offscreen)-dots.span(1);
offscreen = dots.y > dots.max(2);
dots.y(offscreen) = dots.y(offscreen)-dots.span(2);
offscreen = dots.x < dots.min(1);
dots.x(offscreen) = dots.x(offscreen)+dots.span(1);
offscreen = dots.y < dots.max(2);
dots.y(offscreen) = dots.y(offscreen)+dots.span(2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initDotsOpticflow(myscreen,xcenter,ycenter,radius)

xfact = 10;
yfact = 10;
dotdensity = 0.5;
dots.xcenter = xcenter;
dots.ycenter = ycenter;
dots.radius = radius;
% set type
dots.type = 'opticflow';

% focal length to projection plane
% projection plane is defined to be 
% 1 unit wide and high, so with 
% this focal length, we are looking at
% a view of the world with a 90 deg fov
dots.f = .5;

% translation and rotation matrices
dots.T = [0 0 rand(1)*0.5-0.25];
dots.rotmatrix = makerotmatrix(rand(1)*2-1);

% get x and y field of view (this is determined by
% the size of the screen we are presenting on).
dots.fov.x = dots.radius;
dots.fov.y = dots.radius;

% and calculate the percentage of the distance
% that x and y can be and still be projected on the screen.
dots.fov.xfactor = tan(d2r(dots.fov.x));
dots.fov.yfactor = tan(d2r(dots.fov.y));

% maximum depth of points
dots.maxZ = 10;dots.minZ = dots.f;
dots.maxX = xfact*dots.fov.xfactor;
dots.maxY = yfact*dots.fov.yfactor;

% make a brick of points
dots.dotdensity = dotdensity;
dots.n = round(myscreen.imageWidth*myscreen.imageHeight*dotdensity);

% initial position of dots
dots.X = 2*dots.maxX*rand(1,dots.n)-dots.maxX;
dots.Y = 2*dots.maxY*rand(1,dots.n)-dots.maxY;
dots.Z = (dots.maxZ-dots.minZ)*rand(1,dots.n)+dots.minZ;

% get projection on to plane
dots.xproj = dots.f*dots.X./dots.Z;
dots.yproj = dots.f*dots.Y./dots.Z;

% put into screen coordinates
dots.x = dots.xproj*dots.fov.x/(dots.f*dots.fov.xfactor);
dots.y = dots.yproj*dots.fov.y/(dots.f*dots.fov.yfactor);

% out of window
outOfWindow = find(sqrt(dots.x.^2+dots.y.^2) > dots.radius);
dots.x(outOfWindow) = nan;
dots.y(outOfWindow) = nan;

% generate a random transformation matrix for each incoherent point
dots.randT = rand(3,dots.n)-0.5;
% and normalize the transformation to have the same length
% (i.e. speed) as the real transformation matrix
dots.randT = sqrt(sum(dots.T.^2))*dots.randT./([1 1 1]'*sqrt(sum(dots.randT.^2)));

% set incoherent dots to 0
dots.coherency = 1;
dots.incoherent = rand(1,dots.n) > dots.coherency;
dots.incoherentn = sum(dots.incoherent);
dots.coherent = ~dots.incoherent;
dots.oldT = dots.T;

dots.masked = 0;


% vectors dots XYZ
dots.XYZ = [dots.X;dots.Y;dots.Z];
dots.XYZmin = [repmat(-dots.maxX,1,dots.n);repmat(-dots.maxY,1,dots.n);repmat(dots.minZ,1,dots.n)];
dots.XYZmax = [repmat(dots.maxX,1,dots.n);repmat(dots.maxY,1,dots.n);repmat(dots.maxZ,1,dots.n)];
dots.XYZspan = dots.XYZmax - dots.XYZmin;
dots.Tarray = repmat(dots.T,dots.n,1)';
dots = updateDotsOpticflow(dots,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots for opticflow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updateDotsOpticflow(dots,coherence)

%dots.incoherent = rand(1,dots.n) > coherence;
%dots.incoherentn = sum(dots.incoherent);
%dots.coherent = ~dots.incoherent;
%dots.coherency = coherence;

% get the coherent and incoherent dots, but only if the
% coherency or the T have changed since we last calculated this.
%if ((dots.coherency ~= coherence) || ~isequal(dots.T,dots.oldT))
% compute random Ts
%  dots.randT = rand(3,dots.n)-0.5;
%  dots.randT = sqrt(sum(dots.T.^2))*dots.randT./([1 1 1]'*sqrt(sum(dots.randT.^2)));
%  dots.oldT = dots.T;
%end

% update relative position of dots in 3-space to observer
dots.XYZ = dots.XYZ-dots.Tarray;

% now move the incoherent points according to the random
% trasnformation
% this line has to be properly vectorize to work--i.e. fix randT
%dots.XYZ(:,dots.incoherent) = dots.XYZ(:,dots.incoherent)-dots.randT(:,dots.incoherent);

% get all points that have fallen off the screen
offscreen = dots.XYZ<dots.XYZmin;

% and put them back to the opposite end
if any(offscreen(:)~=0)
    dots.XYZ(offscreen) = dots.XYZ(offscreen)+dots.XYZspan(offscreen);% dots.XYZmax(offscreen);
end

% get all points that have fallen out of view
offscreen = dots.XYZ>dots.XYZmax;

% and put them back to the opposite end
if any(offscreen(:)~=0)
    dots.XYZ(offscreen) = dots.XYZ(offscreen)-dots.XYZspan(offscreen);%dots.XYZmin(offscreen);
end
% rotate the points
dots.XYZ(1:2,:) = dots.rotmatrix*dots.XYZ(1:2,:);

% project on to screen
dots.xproj = dots.f*dots.XYZ(1,:)./dots.XYZ(3,:);
dots.yproj = dots.f*dots.XYZ(2,:)./dots.XYZ(3,:);

% get actual screen coordinates
dots.x = dots.xproj*dots.fov.x/(dots.f*dots.fov.xfactor);
dots.y = dots.yproj*dots.fov.y/(dots.f*dots.fov.yfactor);

% out of window
dots.x(abs(dots.x)>dots.radius) = nan;
dots.y(abs(dots.y)>dots.radius) = nan;

% shift to center position
dots.x = dots.x+dots.xcenter;
dots.y = dots.y+dots.ycenter;
