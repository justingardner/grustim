% DESCRIPTION: 
% This script runs the analysis for the detection task in the Najemnik & Gesiler 2005 Nature paper. read geislerDetectionAnalysisMultipleLocs for
% a full description. 
% This script is designed to only analyze data for one location and thus is compatible with geislerDetectionTaskOneLoc only

function geislerDetectionAnalysis

% default return argument
fit = [];

% default to working on current directory
if nargin < 1, stimfileNames = [];end

% parse arguments
%getArgs(varargin);
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
  dispHeader(sprintf('(geislerDetectionAnalysis) Analyzing file (%i/%i): %s',iFile,length(stimfileNames),stimfileNames{iFile}));
  
  % load and parse the stimfile
  d = loadStimfile(fullfile(e.path,stimfileNames{iFile}));

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
  disp(sprintf('(geislerDetectionAnalysis) No files found'));
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setting the location variable
location = d.task{1}.location;
d.location.x = location(1);
d.location.y = location(2);

% Each condition (i.e location) ran for 544 trials
trialNums = [1:544];
% The contrast values for the first 200 trials
d.contrast = d.parameter.contrast(trialNums);
% A unique set (no duplicates) of the contrasts used in the first 200 trials
d.uniquecontrast = unique(d.contrast);
% The responses for the first 200 trials
correct = d.task{1}.response.correct(trialNums);
    
% We are going to iterate through the unique contrasts
for i = 1:length(d.uniquecontrast)
    % find() will return the indeces of the the matrix that has the contrast values for the first 200 trials for which the contrast value
    % matches the contrast value we are currently iterating on. Those indeces are also trial numbers
    whichTrials = find(d.contrast == d.uniquecontrast(i));
    % The number of trials that had the contrast value we are iterating on (should be equal for each contrast if we randomized correctly)
    nTrials = length(whichTrials);
    % correct(whichTrials) returns the values of the response array for the trials that had the contrast value we are iterationg on
    % Summing them and dividing by the number of trials gives us the percent correct for that contrast value
    d.correctBinned(i) = sum(correct(whichTrials))/nTrials;
    % This just saves the number of trials for that contrast value into a variable
    d.nTrials(i) = nTrials;
end


% Setiing all performance values below 0.5 to 0.5 becuase they are theoretically at chance performance
for i = 1: length(d.correctBinned)
    if d.correctBinned(i) < 0.5
        d.correctBinned(i) = 0.5;
    end
end
% Scaling so that all values are between 0 and 1 (important for fitting the cumalitve gaussian)
d.correctBinned =  ( 2 * d.correctBinned ) - 1;
% fit a cumulative gaussian to data
d.fit = fitCumulativeGaussian(d.uniquecontrast,d.correctBinned);
% find the threshold contrast (target contrast at 82.01% 0r 82.02% performance)
idx = find(d.fit.fitY > 0.820 & d.fit.fitY < 0.823);
d.thresholdContrast = d.fit.fitX(idx);


% STOP HERE WHEN DEBUGGING
k=2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTES:             
% To draw a scatter plot
scatter(d.cond(1).uniquecontrast, d.cond(1).correctBinned)

% To fit a curve
plot(d.fit(1).fitX, d.fit(1).fitY)

% Error
errors = [];
n = 32;
for i=1:length(d.cond(1).correctBinned)
    p = d.cond(1).correctBinned(i);
    error = sqrt( (p*(1-p)) / n );
    errors = [errors error];
end

errorbar(d.cond(1).uniquecontrast, d.cond(1).correctBinned, errors);




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
  disp(sprintf('(geislerDetectionAnalysis:loadStimfile) Could not find stimfile: %s',stimfileName));
  return
end

% load the file
s = load(stimfileName);
if ~isfield(s,'myscreen') || ~isfield(s,'task')
  disp(sprintf('(geislerDetectionAnalysis:loadStimfile) No myscreen or task in stimfile: %s',stimfileName));
  return
end

% Create task filename
taskFilename = s.task{1}.taskFilename;
% Check task filename (Commented out, so code now assumes you are using the correct file)
%if isempty(strfind(lower(taskFilename),'search')) & isempty(strfind(lower(taskFilename),'search'))
%   disp(sprintf('(geislerDetectionAnalysis:loadStimfile) Incorrect task in stimfile: %s',taskFilename));
%   return
%end

% parse into parameters
d = getTaskParameters(s.myscreen,s.task);
%d = d{1};

% print what we found
disp(sprintf('(geislerDetectionAnalysis:loadStimfile) Found task: %s (%i trials) SID: %s Date: %s',taskFilename,d.nTrials,s.myscreen.SID,s.myscreen.starttime));

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
    disp(sprintf('(geislerDetectionAnalysis) Could not find %s',stimfileNames));
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
