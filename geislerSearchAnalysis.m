% WRITTEN BY: 
% (1) Yehia Elkersh
% (2) Josh Wilson (loading stimfiles)

% DESCRIPTION: 
% This scripts runs the analysis for the search task in the Najemnik & Gesiler 2005 Nature paper. 
% As of July 29, 2021, the only analysis being done is plotting the eyes traces for a particular trial

% NOTES:
% (1) This needs to be double-checked, but I think there is a parameter in d called reactionTime, which corresponds to the exact time when the subject pressed the
% response button (here the number 2). I am also pretty sure that the x and y positions if the eye are indexed by time, so to figure out where the subject was
% fixating when they pressed the response button (in order to determine if they were fixating close to the target or not), you can get the
% reactionTime of the trial and use it to index into the x and y position of their eye and then compare it with the position of the target on that trial




function geislerSearchAnalysis
% default return argument
fit = [];

% default to working on current directory
if nargin < 1, stimfileNames = [];end

% parse arguments
% getArgs(varargin);
% get filenames and path
[e.path stimfileNames] = getStimfileNames(stimfileNames);
if isempty(e.path),return,end

% check for valid files as we go through
% nFiles will contain how many vaild files we have
e.nFiles = 0;
e.visualStaircase = {};
e.auditoryStaircase = {};

% cycle through all files
for iFile = 1:length(stimfileNames)
  % display what is happening
  dispHeader(sprintf('(geislerSearchAnalysis) Analyzing file (%i/%i): %s',iFile,length(stimfileNames),stimfileNames{iFile}));
  
  % load and parse the stimfile
  % d stand for task (d)ata
  d = loadStimfile(fullfile(e.path,stimfileNames{iFile}));
  % t stands for eye (t)race data
  t = getTaskEyeTraces(fullfile(e.path,stimfileNames{iFile}));

  % valid file, so keep its information
  if ~isempty(d)
      % update count
      e.nFiles = e.nFiles + 1;
      % and list of filenames
      e.filenames{e.nFiles} = stimfileNames{iFile};
      % and the data
      e.d{e.nFiles} = d;
  end
end

% if no valid files found return
if e.nFiles == 0
  disp(sprintf('(geislerSearchAnalysis) No files found'));
  return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots traces of the eye movements on the first trial
x = t.eye.xPos(1, :);
y = t.eye.yPos(1, :);
plot(x,y)

% STOP HERE WHEN DEBUGGING
k=2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    loadStimfile    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = loadStimfile(stimfileName)

% default to empty
d = [];

% add mat extension
stimfileName = setext(stimfileName,'mat');

% see if it exists
if ~isfile(stimfileName)
  disp(sprintf('(geislerSearchAnalysis:loadStimfile) Could not find stimfile: %s',stimfileName));
  return
end

% load the file
s = load(stimfileName);
if ~isfield(s,'myscreen') || ~isfield(s,'task')
  disp(sprintf('(geislerSearchAnalysis:loadStimfile) No myscreen or task in stimfile: %s',stimfileName));
  return
end

% check task filename
taskFilename = s.task{1}.taskFilename;
if isempty(strfind(lower(taskFilename),'search')) & isempty(strfind(lower(taskFilename),'search'))
  disp(sprintf('(geislerSearchAnalysis:loadStimfile) Incorrect task in stimfile: %s',taskFilename));
  return
end

% parse into parameters
d = getTaskParameters(s.myscreen,s.task);
%d = d{1};

% print what we found
disp(sprintf('(geislerSearchAnalysis:loadStimfile) Found task: %s (%i trials) SID: %s Date: %s',taskFilename,d.nTrials,s.myscreen.SID,s.myscreen.starttime));

% get the variables
d.myscreen = s.myscreen;
d.task = s.task;
%d.stimulus = s.stimulus;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getStimfileNames    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stimfilePath stimfileNames] = getStimfileNames(stimfileNames)

%  check if we are examining a single mat file
if isfile(setext(stimfileNames,'mat'))
  % make sure the extension is .mat
  stimfileNames = setext(stimfileNames,'mat');
  % first check to see if it has a path
  if ~isempty(fileparts(stimfileNames))
    stimfilePath = fileparts(stimfileNames);
    stimfileNames = getLastDir(stimfileNames);
  else
    stimfilePath = pwd;
  end
else
  % not a single file, so go look for the path
  if isempty(stimfileNames)
    % get current directory
    stimfilePath = pwd;
    % check if it is a directory
  elseif isdir(stimfileNames)
    % then use that as the path
    stimfilePath = stimfileNames;
    % see if it is a subject ID
  elseif ((length(stimfileNames)>1) && (isequal(lower(stimfileNames(1)),'s'))) || isnumeric(stimfileNames)
    % turn a numeric into a string SID
    if isnumeric(stimfileNames)
      stimfileNames = sprintf('s%03i',stimfileNames);
    end
    % get the path for this subject
    stimfilePath = fullfile('~/data/search',stimfileNames);
  else
    disp(sprintf('(geislerSearchAnalysis) Could not find %s',stimfileNames));
    stimfilePath = '';
    stimfileNames = '';
    return
  end
  % get everything in the directory that is a mat file
  matfiles = dir(fullfile(stimfilePath,'*.mat'));
  stimfileNames = {matfiles(:).name};
end

% make sure we are returning a cell array
stimfileNames = cellArray(stimfileNames);
