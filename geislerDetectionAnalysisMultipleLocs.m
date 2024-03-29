% WRITTEN BY: 
% (1) Yehia Elkersh
% (2) Josh Wilson (loading stimfiles)

% DESCRIPTION: 
% This scripts runs the analysis for the detection task in the Najemnik & Gesiler 2005 Nature paper. It is expecting a stimfile where the task was run
% at multiple locations in the same experimet, and it also expects an array of the locations used in the experiment (in the order that they were used) 
% and turns each location into a "condition" (here called 'cond') and runs the analysis on each condition.

% NOTES:
% (1) As of July 15, 2021 this script is compatible with the geislerDetectionTaskMultipleLocs file. For instance, it can only analyze two locations, 
% with 544 trials at each location. Since the aforementioned task file needs to undergo some major changes (such as being broken up into phases), 
% this script will need to be adjusted accordingly.
% (2) In order to have access to the data in the command line (e.g. to plot the psychometric curve), you should stop the script at the line 'k = 2' 
% because when the script ends, the data is no longer available


function geislerDetectionAnalysis
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
% Analysis   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A condition(cond) refers to a particular (x,y) location, so cond{1} will contain the data related to location 1 
% (which is the first location in which the target was presented)
d.nCond = [1:length(d.task{1}.locations)];

% creates a condition for each location
for i=1:length(d.nCond)
    d.cond(i).location = d.task{1}.locations(i,:);
    d.cond(i).x = d.cond(i).location(1);
    d.cond(i).y = d.cond(i).location(2);
end

% Each condition (i.e location) ran for 544 trials
d.condTrialNums{1} = [1:544];
d.condTrialNums{2} = [545:1088];

% Commented for the first iteration
for iCond = 1:length(d.nCond)
    % Indeces for the trials in the first condition: [1, 2, ... , 544]
    trialNums = d.condTrialNums{iCond}; 
    % The contrast values for the first 544 trials
    d.cond(iCond).contrast = d.parameter.contrast(trialNums);
    % A unique set (no duplicates) of the contrasts used in the first 544 trials
    d.cond(iCond).uniquecontrast = unique(d.cond(iCond).contrast);
    % The responses for the first 544 trials
    correct = d.task{1}.response.correct(trialNums);
    
    % We are going to iterate through the unique contrasts
    for iVal = 1:length(d.cond(iCond).uniquecontrast)
        % find() will return the indeces of the the matrix that has the contrast values for the first 544 trials for which the contrast value
        % matches the contrast value we are currently iterating on. Those indeces are also trial numbers.
        whichTrials = find(d.cond(iCond).contrast == d.cond(iCond).uniquecontrast(iVal));
        % nTrials is the number of trials that had the contrast value we are iterating on (should be equal for each contrast if we randomized correctly)
        nTrials = length(whichTrials);
        % correct(whichTrials) returns the values of the response array for the trials that had the contrast value we are iterationg on
        % Summing them and dividing by the number of trials gives us the percent correct for that contrast value
        d.cond(iCond).correctBinned(iVal) = sum(correct(whichTrials))/nTrials;
        % This just saves the number of trials for that contrast value into a variable
        d.cond(iCond).nTrials(iVal) = nTrials;
    end
end

for iCond = 1:length(d.nCond)
    % Setiing all performance values below 0.5 to 0.5 becuase they are theoretically at chance performance
    for i = 1: length(d.cond(iCond).correctBinned)
        if d.cond(iCond).correctBinned(i) < 0.5
            d.cond(iCond).correctBinned(i) = 0.5;
        end
    end
    % Scaling so that all values are between 0 and 1 (important for fitting the cumalitve gaussian)
    d.cond(iCond).correctBinned =  ( 2 * d.cond(iCond).correctBinned ) - 1;
    % fit a cumulative gaussian to data
    d.fit(iCond) = fitCumulativeGaussian(d.cond(iCond).uniquecontrast,d.cond(iCond).correctBinned);
    % find the threshold contrast (target contrast at 82.02% performance)
    idx = find(d.fit(iCond).fitY > 0.8201 & d.fit(iCond).fitY < 0.8203);
    d.cond(iCond).thresholdContrast = d.fit(iCond).fitX(idx);
end

% STOP HERE WHEN DEBUGGING
k=2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTES:             
% To draw a scatter plot
% scatter(d.cond(1).uniquecontrast, d.cond(1).correctBinned)

% To fit a curve
% plot(d.fit(1).fitX, d.fit(1).fitY)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loadStimfile    
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
% getStimfileNames    
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
