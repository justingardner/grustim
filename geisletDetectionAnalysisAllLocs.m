
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Making dataMatrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make a matrix where the rows represent data points (x, y, mean, std, thresholdContrast) and the columns represent files
dataMatrix = [];


% setting the location variable
for iFile = 1:e.nFiles
    location = e.d{iFile}.task{1}.location;
    e.d{iFile}.location.x = location(1);
    e.d{iFile}.location.y = location(2);

    % Each condition (i.e location) ran for 175 trials 
    trialNums = [1:175];
    % The contrast values for the first 200 trials
    e.d{iFile}.contrast = e.d{iFile}.parameter.contrast(trialNums);
    % A unique set (no duplicates) of the contrasts used in the first 200 trials
    e.d{iFile}.uniquecontrast = unique(e.d{iFile}.contrast);
    % The responses for the first 200 trials
    correct = e.d{iFile}.task{1}.response.correct(trialNums);
    
    % We are going to iterate through the unique contrasts
    for i = 1:length(e.d{iFile}.uniquecontrast)
        % find() will return the indeces of the the matrix that has the contrast values for the first 200 trials for which the contrast value
        % matches the contrast value we are currently iterating on. Those indeces are also trial numbers
        whichTrials = find(e.d{iFile}.contrast == e.d{iFile}.uniquecontrast(i));
        % The number of trials that had the contrast value we are iterating on (should be equal for each contrast if we randomized correctly)
        nTrials = length(whichTrials);
        % correct(whichTrials) returns the values of the response array for the trials that had the contrast value we are iterationg on
        % Summing them and dividing by the number of trials gives us the percent correct for that contrast value
        e.d{iFile}.correctBinned(i) = sum(correct(whichTrials))/nTrials;
        % This just saves the number of trials for that contrast value into a variable
        e.d{iFile}.nTrials(i) = nTrials;
    end

    % Setiing all performance values below 0.5 to 0.5 becuase they are theoretically at chance performance
    for i = 1: length(e.d{iFile}.correctBinned)
        if e.d{iFile}.correctBinned(i) < 0.5
            e.d{iFile}.correctBinned(i) = 0.5;
        end
    end
    % Scaling so that all values are between 0 and 1 (important for fitting the cumalitve gaussian)
    e.d{iFile}.correctBinned =  ( 2 * e.d{iFile}.correctBinned ) - 1;
    % Fit a cumulative gaussian to data
    e.d{iFile}.fit = fitCumulativeGaussian(e.d{iFile}.uniquecontrast,e.d{iFile}.correctBinned);
    % Find the threshold contrast (target contrast at 82% performance)
    % For the cummalitve gaussian (the function we used to fit the data), the y-values are the intergrals of a Normal(mu, sigma) distribution from -inf to the x-values
    mu = e.d{iFile}.fit.mean;
    sigma = e.d{iFile}.fit.std;
    thresholdContrast = norminv(0.82,mu,sigma);
    e.d{iFile}.thresholdContrast = thresholdContrast;

    % Fill out dataMatrix
    dataMatrix(1, iFile) = e.d{iFile}.location.x;
    dataMatrix(2, iFile) = e.d{iFile}.location.y;
    dataMatrix(3, iFile) = e.d{iFile}.fit.mean;
    dataMatrix(4, iFile) = e.d{iFile}.fit.std;
    dataMatrix(5, iFile) = e.d{iFile}.thresholdContrast;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Graphing 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Graph the psychometric cruves fir the four different eccentricities along each location axis
    eccen = num2str(sqrt( (dataMatrix(1, iFile))^2 + (dataMatrix(2, iFile))^2 ));
   
    plot(e.d{iFile}.fit.fitX, e.d{iFile}.fit.fitY, 'DisplayName', eccen)
    
    hold on
    title('Axis 0')
    
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Graph 25 psychometric curves            
    % To draw a scatter plot
    figure(iFile)
    scatter(e.d{iFile}.uniquecontrast, e.d{iFile}.correctBinned)
    hold on 

    % To fit a curve with titles
    plot(e.d{iFile}.fit.fitX, e.d{iFile}.fit.fitY)
    titleStr = sprintf('X location: %.3d. Y Location: %.3d // mean = %.3d, std = %.3d',e.d{iFile}.location.x,e.d{iFile}.location.y,e.d{iFile}.fit.mean,e.d{iFile}.fit.std)
    title(titleStr)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
end
legend show;
legend('Location', 'southeast')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating Statistics on dataMatrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setting data points to their indeces in the dataMatrix so code is more legible
X = 1;
Y = 2;
Mean = 3;
Std = 4;
threshCon = 5;
    

threshCons0 = []; stds0 = [];
threshCons225 = []; stds225 = [];
threshCons45 = []; stds45 = [];
threshCons675 = []; stds675 = [];

% Seperate by eccentricity (eccen)
for iCol=1:e.nFiles
    
    eccen = sqrt( (dataMatrix(X, iCol))^2 + (dataMatrix(Y, iCol))^2 );
    
    if  eccen == 0
        threshCons0 = [threshCons0 dataMatrix(threshCon, iCol)];
        stds0 = [stds0 dataMatrix(Std, iCol)];
    end
    
    if eccen > 1 & eccen < 3
        threshCons225 = [threshCons225 dataMatrix(threshCon, iCol)];
        stds225 = [stds225 dataMatrix(Std, iCol)];
    end
    
    if eccen > 4 & eccen < 5
        threshCons45 = [threshCons45 dataMatrix(threshCon, iCol)];
        stds45 = [stds45 dataMatrix(Std, iCol)];
    end
    
    if eccen > 6
        threshCons675 = [threshCons675 dataMatrix(threshCon, iCol)];
        stds675 = [stds675 dataMatrix(Std, iCol)];
    end
end

AvgthreshCon0 = sum(threshCons0) / length(threshCons0)
AvgthreshCon225 = sum(threshCons225) / length(threshCons225)
AvgthreshCon45 = sum(threshCons45) / length(threshCons45)
AvgthreshCon675 = sum(threshCons675) / length(threshCons675)

    
% STOP HERE WHEN DEBUGGING
k=2


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

