% taskTemplateMovie.m
%
%        $Id: taskTemplate.m 835 2010-06-29 04:04:08Z justin $
%      usage: taskTemplate
%         by: justin gardner
%       date: 1/20/2013
%  copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%    purpose: example program to show how to use the task structure
%             and display a movie clip with a fixation cross over it.
%
% McGurk stimuli from Experiment 2 in Basu Mallick D, Magnotti JF, Beauchamp MS. Variability and stability in the McGurk effect: contributions of participants, stimuli, time, and response type. Psychonomic Bulletin and Review. (2015) 22:1299–1307. DOI 10.3758/s13423-015-0817-4

function myscreen = mcgurk(varargin)
mglEatKeys('123`');
global stimulus

% get arguments
visualOnly=0; auditoryOnly=0; av=0; avcong=0;avincong=0;
getArgs(varargin,{'visualOnly=0','auditoryOnly=0','avcong=0','avincong=0','av=0'},'verbose=1');

% % check arguments
% if ~any(nargin == [0])
%   help taskTemplate
%   return
% end
stimulus.visualOnly = visualOnly;
stimulus.auditoryOnly = auditoryOnly;
stimulus.av = av;
stimulus.avcong = avcong;
stimulus.avincong = avincong;

% initalize the screen
mglSetParam('movieMode',1);
myscreen = initScreen;

% load the movies
frameRate = 30;
stimulus = myInitStimulus(stimulus,myscreen,frameRate);

% set up task
task{1}.waitForBacktick = 0;
task{1}.seglen = [60/frameRate inf 0.5];
task{1}.getResponse = [0 1 0];

if visualOnly || auditoryOnly || avcong
task{1}.parameter.movieNum = 1:stimulus.congruent.numMovies;
task{1}.numBlocks = 1;
elseif avincong
  task{1}.parameter.movieNum = 1:stimulus.conflict.numMovies;
  task{1}.numBlocks = 3;
elseif av
  task{1}.parameter.condition = {'da','ba','ga','conflict','conflict','conflict'};
  task{1}.parameter.movieNum = 1:8;
  task{1}.numBlocks = 1;
end
task{1}.random = 1;
task{1}.waitForBacktick = 1;

task{1}.randVars.calculated.label = nan;
task{1}.randVars.calculated.thisCond = nan;
task{1}.randVars.calculated.rt = nan;
task{1}.randVars.calculated.resp = nan;

% initialize the task
for phaseNum = 1:length(task)
  [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end

% init the stimulus
myscreen = initStimulus('stimulus',myscreen);

mglWaitSecs(1);
mglClearScreen(0);
mglTextSet([],50,[1 1 1]);
mglTextDraw('READY',[0 0]);
mglFlush;
mglClearScreen(0);
mglTextDraw('READY',[0 0]);
mglFlush;
mglWaitSecs(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task)) && ~myscreen.userHitEsc
  % update the task
  [task myscreen phaseNum] = updateTask(task,myscreen,phaseNum);
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
% clear screen to transparent
mglClearScreen([0 0 0 0]);

% on first segment
if (task.thistrial.thisseg == 1)
  if stimulus.visualOnly
    vol=mglVolume([0,0]);
  else
    vol=mglVolume([0.5,0.5]);
  end

  % show movie
  startMovie(task);
  if stimulus.auditoryOnly
    mglFillRect(0,0,[1024 728],[0 0 0]);
    mglFixationCross(2,2,[0 1 1]);
  end
 
  mglFlush;

elseif task.thistrial.thisseg == 2
  % stop movie 
  stopMovie;
  % display fixation cross in white
  % mglFixationCross(1,1,[1 1 1]);
  mglTextSet([],72,[1 1 1]);
  mglTextDraw('Ba',[-5 0]);
  mglTextDraw('Da', [0 0]);
  mglTextDraw('Ga', [5 0]);
  mglFlush;
else
  mglFixationCross(2,2,[1 1 1]);
  mglFlush;

end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus
if task.thistrial.thisseg == 1
  stepMovie;
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)

global stimulus

% here, we just check whether this is the first time we got a response
% this trial and display what the subject's response was and the reaction time
if task.thistrial.gotResponse < 1
  disp(sprintf('Subject response: %i Reaction time: %0.2fs',task.thistrial.whichButton,task.thistrial.reactionTime));

  task.thistrial.rt = task.thistrial.reactionTime;
  task.thistrial.resp = task.thistrial.whichButton;

  task = jumpSegment(task);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen,frameRate)

% clear stimulus;

% remember the desired frame rate
stimulus.frameRate = frameRate;

% look for directory where movies live
moviePath = fullfile(fileparts(which('mcgurk')),'mcgurk_stim');
if ~isdir(moviePath)
  disp(sprintf('(taskTemplateMovie) Could not find movie dir %s',movieDir));
  mglClose;
  keyboard
end

movieDir_da = dir([moviePath,'/congruent_da']);
movieDir_ba = dir([moviePath,'/congruent_ba']);
movieDir_ga = dir([moviePath,'/congruent_ga']);
movieDir_conflict = dir([moviePath,'/conflict']);

% load names of files in movie directory
% movieDir = dir(moviePath);
movieNames_congruent = {};
movieNames_conflict = {};
stimulus.congruent.label = {};
stimulus.conflict.label = {};
stimulus.congruent.frameTimes ={};
stimulus.conflict.frameTimes={};

% see if any names end with the right extensions and keep a list of those

for i = 1:length(movieDir_da)
  [path name ext] = fileparts(movieDir_da(i).name);
  if any(strcmp({'.mp4','.mov'},ext))
    movieNames_congruent{end+1} = fullfile([moviePath,'/congruent_da'],movieDir_da(i).name);
    stimulus.congruent.label{end+1} = 'da';
  end
end
for i = 1:length(movieDir_ba)
  [path name ext] = fileparts(movieDir_ba(i).name);
  if any(strcmp({'.mp4','.mov'},ext))
    movieNames_congruent{end+1} = fullfile([moviePath,'/congruent_ba'],movieDir_ba(i).name);
    stimulus.congruent.label{end+1} = 'ba';
  end
end
for i = 1:length(movieDir_ga)
  [path name ext] = fileparts(movieDir_ga(i).name);
  if any(strcmp({'.mp4','.mov'},ext))
    movieNames_congruent{end+1} = fullfile([moviePath,'/congruent_ga'],movieDir_ga(i).name);
    stimulus.congruent.label{end+1} = 'ga';
  end
end

% check to make sure we found some movie names
if length(movieNames_congruent) < 1
  disp(sprintf('(taskTemplateMovie) Could not find any movies in directory %s\n',moviePath));
  mglClose;
  keyboard
end

for i = 1:length(movieDir_conflict)
  [path name ext] = fileparts(movieDir_conflict(i).name);
  if any(strcmp({'.mp4','.mov'},ext))
    movieNames_conflict{end+1} = fullfile([moviePath,'/conflict'],movieDir_conflict(i).name);
    stimulus.conflict.label{end+1} = 'conflict';
  end
end
% check to make sure we found some movie names
if length(movieNames_conflict) < 1
  disp(sprintf('(taskTemplateMovie) Could not find any movies in directory %s\n',moviePath));
  mglClose;
  keyboard
end

if ~stimulus.avincong
% go through each name
for iMovie = 1:length(movieNames_congruent)

  % display what we are doing
  % mglClearScreen;
  disp(sprintf('Loading cong movies (%i/%i)',iMovie,length(movieNames_congruent)));
  % mglFlush;

  % load movie
  m = mglMovie(movieNames_congruent{iMovie});
  if ~isempty(m)
    stimulus.congruent.m(iMovie) = m;
  else
    mglClose;
    keyboard
  end

  % if first one, show and hide
  if iMovie == 1
    % get the current position and move offscreen
    originalPos = mglMovie(stimulus.congruent.m(iMovie),'getPosition');
    offscreenPos = originalPos;
    offscreenPos(1) = 100000;
    mglMovie(stimulus.congruent.m(1),'move',offscreenPos);
    % display and hide
    mglMovie(stimulus.congruent.m(1),'show');
    mglMovie(stimulus.congruent.m(1),'hide');
    % move back to original position
    mglMovie(stimulus.congruent.m(1),'move',originalPos);
  end

  % count the number of frames, remembering the current time for each frame
  % so that we can use that to step the movie
  mglMovie(stimulus.congruent.m(iMovie),'gotoBeginning');
  frameNum = 1;
  stimulus.congruent.frameTimes{iMovie}{frameNum} = mglMovie(stimulus.congruent.m(iMovie),'getCurrentTime');
  mglMovie(stimulus.congruent.m(iMovie),'stepForward');
  % the movie is over when we find a frame with the same time stamp as the last one
  while ~strcmp(stimulus.congruent.frameTimes{iMovie}{end},mglMovie(stimulus.congruent.m(iMovie),'getCurrentTime'))
    % step forward again
    frameNum = frameNum+1;
    stimulus.congruent.frameTimes{iMovie}{frameNum} = mglMovie(stimulus.congruent.m(iMovie),'getCurrentTime');
    mglMovie(stimulus.congruent.m(iMovie),'stepForward');
  end
  stimulus.congruent.numFrames(iMovie) = frameNum;

  % go back to beginning to get it ready to play again
  mglMovie(stimulus.congruent.m(iMovie),'setCurrentTime',stimulus.congruent.frameTimes{iMovie}{1});
  
  % display number of frames counted
  disp(sprintf('(taskTemplateMovie) Counted %i frames in cong movie %s',stimulus.congruent.numFrames(iMovie),movieNames_congruent{iMovie}));
end

stimulus.congruent.numMovies = length(movieNames_congruent);
end

if stimulus.av || stimulus.avincong
% go through each name
for iMovie = 1:length(movieNames_conflict)

  % display what we are doing
  % mglClearScreen;
  disp(sprintf('Loading conf movies (%i/%i)',iMovie,length(movieNames_conflict)));
  % mglFlush;

  % load movie
  m = mglMovie(movieNames_conflict{iMovie});
  if ~isempty(m)
    stimulus.conflict.m(iMovie) = m;
  else
    mglClose;
    keyboard
  end

  % count the number of frames, remembering the current time for each frame
  % so that we can use that to step the movie
  mglMovie(stimulus.conflict.m(iMovie),'gotoBeginning');
  frameNum = 1;
  stimulus.conflict.frameTimes{iMovie}{frameNum} = mglMovie(stimulus.conflict.m(iMovie),'getCurrentTime');
  mglMovie(stimulus.conflict.m(iMovie),'stepForward');
  % the movie is over when we find a frame with the same time stamp as the last one
  while ~strcmp(stimulus.conflict.frameTimes{iMovie}{end},mglMovie(stimulus.conflict.m(iMovie),'getCurrentTime'))
    % step forward again
    frameNum = frameNum+1;
    stimulus.conflict.frameTimes{iMovie}{frameNum} = mglMovie(stimulus.conflict.m(iMovie),'getCurrentTime');
    mglMovie(stimulus.conflict.m(iMovie),'stepForward');
  end
  stimulus.conflict.numFrames(iMovie) = frameNum;

  % go back to beginning to get it ready to play again
  mglMovie(stimulus.conflict.m(iMovie),'setCurrentTime',stimulus.conflict.frameTimes{iMovie}{1});
  
  % display number of frames counted
  disp(sprintf('(taskTemplateMovie) Counted %i frames in conf movie %s',stimulus.conflict.numFrames(iMovie),movieNames_conflict{iMovie}));
end

stimulus.conflict.numMovies = length(movieNames_conflict);
end
% % display what we are doing
% mglClearScreen;mglFlush;
% mglClearScreen;mglFlush;


%%%%%%%%%%%%%%%%%%%%
%    startMovie    %
%%%%%%%%%%%%%%%%%%%%
function startMovie(task)

global stimulus;
if stimulus.av
  if strcmp(char(task.thistrial.condition),'da')
    stimulus.thismovie = stimulus.congruent.m(task.thistrial.movieNum);
    stimulus.thisNumFrames = stimulus.congruent.numFrames(task.thistrial.movieNum);
    stimulus.thisFrameTimes = stimulus.congruent.frameTimes{task.thistrial.movieNum};
  elseif strcmp(char(task.thistrial.condition),'ba')
    stimulus.thismovie = stimulus.congruent.m(8 + task.thistrial.movieNum);
    stimulus.thisNumFrames = stimulus.congruent.numFrames(8 + task.thistrial.movieNum);
    stimulus.thisFrameTimes = stimulus.congruent.frameTimes{8 + task.thistrial.movieNum};
  elseif strcmp(char(task.thistrial.condition),'ga')
    stimulus.thismovie = stimulus.congruent.m(16 + task.thistrial.movieNum);
    stimulus.thisNumFrames = stimulus.congruent.numFrames(16 + task.thistrial.movieNum);
    stimulus.thisFrameTimes = stimulus.congruent.frameTimes{16 + task.thistrial.movieNum};
  elseif strcmp(char(task.thistrial.condition),'conflict')
    stimulus.thismovie = stimulus.conflict.m(task.thistrial.movieNum);
    stimulus.thisNumFrames = stimulus.conflict.numFrames(task.thistrial.movieNum);
    stimulus.thisFrameTimes = stimulus.conflict.frameTimes{task.thistrial.movieNum};
  end
elseif stimulus.avincong
  stimulus.thismovie = stimulus.conflict.m(task.thistrial.movieNum);
  stimulus.thisNumFrames = stimulus.conflict.numFrames(task.thistrial.movieNum);
  stimulus.thisFrameTimes = stimulus.conflict.frameTimes{task.thistrial.movieNum};
else
  stimulus.thismovie = stimulus.congruent.m(task.thistrial.movieNum);
  stimulus.thisNumFrames = stimulus.congruent.numFrames(task.thistrial.movieNum);
  stimulus.thisFrameTimes = stimulus.congruent.frameTimes{task.thistrial.movieNum};
end
% % set the current movie number
% stimulus.movieNum = movieNum;

% show the movie
% mglMovie(stimulus.m(stimulus.movieNum),'show');
mglMovie(stimulus.thismovie,'show');

% debugging - display which frame is being drawn
mydisp(sprintf('Displaying frame: '));

% remember the start time
stimulus.movieStartTime = mglGetSecs;

%%%%%%%%%%%%%%%%%%%
%    stepMovie    %
%%%%%%%%%%%%%%%%%%%
function stepMovie

global stimulus;

% figure out which frame to show
frameNum = ceil(mglGetSecs(stimulus.movieStartTime)*stimulus.frameRate);
frameNum = min(frameNum,stimulus.thisNumFrames);

% set the current time to that frame
mglMovie(stimulus.thismovie,'setCurrentTime',stimulus.thisFrameTimes{frameNum});

% debugging - display which frame is being drawn
mydisp(sprintf('%i ',frameNum));

%%%%%%%%%%%%%%%%%%%
%    stopMovie    %
%%%%%%%%%%%%%%%%%%%
function stopMovie

global stimulus;

% hide the movie
mglMovie(stimulus.thismovie,'hide');

% show how much time was taken
disp(sprintf('\nElapsed time: %f movie frameNum: %i',mglGetSecs(stimulus.movieStartTime),find(strcmp(stimulus.thisFrameTimes,mglMovie(stimulus.thismovie,'getCurrentTime')))));

% go back to beginning to get it ready to play again
mglMovie(stimulus.thismovie,'setCurrentTime',stimulus.thisFrameTimes{1});

