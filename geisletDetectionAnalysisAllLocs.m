% WRITTEN BY: 
% (1) Yehia Elkersh
% (2) Josh Wilson (loading stimfiles)

% DESCRIPTION: 
% This script runs the analysis for the detection task in the Najemnik & Gesiler 2005 Nature paper. For a full description of the 
% analysis, read the header comments in geislerDetectionAnalysisMultipleLocs.
% This script is designed to analyze data for all 25 locations. You should 'cd' into a directory that has 25 stimfiles, where each stimfile has data
% from running the task with one location only (i.e. from running the task file geislerDetectionTaskOneLoc)


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
% Making dataMatrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make a matrix where the rows represent data points (x, y, mean, std, thresholdContrast) and the columns represent different files
dataMatrix = [];

% Here, we are looping over all the files (iFile is a file number), which
% is equivalent to loopong over all locations because each file contained
% one location
for iFile = 1:e.nFiles
    % extracting the location from the file
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
    % Normin(x, mu, sigma) gives the integral from -inf to x, so it technically gives y-value of the function we used to fit the data (i.e. psychometric function)
    % for a particular x-value
    mu = e.d{iFile}.fit.mean;
    sigma = e.d{iFile}.fit.std;
    thresholdContrast = norminv(0.82,mu,sigma);
    e.d{iFile}.thresholdContrast = thresholdContrast;

    % Filling out dataMatrix
    % NOTE: each column corresponds to a file and each row corresponds to a data point (i.e. the first row is the x location from all the files)
    dataMatrix(1, iFile) = e.d{iFile}.location.x;
    dataMatrix(2, iFile) = e.d{iFile}.location.y;
    dataMatrix(3, iFile) = e.d{iFile}.fit.mean;
    dataMatrix(4, iFile) = e.d{iFile}.fit.std;
    dataMatrix(5, iFile) = e.d{iFile}.thresholdContrast;
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (2)
    % Graph 25 psychometric curves (one for each location) in 25 separate figures           
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
    
    
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (2)
    % Graph the psychometric cruves for the four different eccentricities along each location axis
    
    % Eccentricity is calculated in order to deterime coolor for the legend
    eccen = round(sqrt( (dataMatrix(1, iFile))^2 + (dataMatrix(2, iFile))^2 ), 2);
    if  eccen == 0
        color = 'red';
    elseif eccen > 1 & eccen < 3
        color = 'green';
    elseif eccen > 4 & eccen < 5
        color = 'blue';
    elseif eccen > 6
        color = 'magenta';
    end
    % This is used for titling the Axis of each subplots
    Titles = [0 45 90 135 180 225 270 315];
    
    x = e.d{iFile}.location.x;
    y = e.d{iFile}.location.y;
    
    % if this file is the foveal position, plot its psychomteric curver and scatter plot in every subplot
    if x == 0 & y == 0
        for thisSubplot = 1:8
            subplot(2,4,thisSubplot)  
            scatter(e.d{iFile}.uniquecontrast, e.d{iFile}.correctBinned, color)
            hold on 
            plot(e.d{iFile}.fit.fitX, e.d{iFile}.fit.fitY, color)
            errors = [];
            n = 25;
            for i=1:length(e.d{iFile}.correctBinned)
                p = e.d{iFile}.correctBinned(i);
                error = sqrt( (p*(1-p)) / n );
                errors = [errors error];
            end
            errorbar(e.d{iFile}.uniquecontrast, e.d{iFile}.correctBinned, errors, 'LineStyle', 'none', 'Color', color);
            hold on
        end
        
    % Divide the subplots based on the Axis from which each file was recorded
    elseif x > 0 & y == 0 
        thisSubplot = 1;
    elseif x > 0 & y > 0
        thisSubplot = 2;
    elseif x == 0 & y > 0
        thisSubplot = 3;
    elseif x < 0 & y > 0
        thisSubplot = 4;
    elseif x < 0 & y == 0
        thisSubplot = 5;
    elseif x < 0 & y < 0
        thisSubplot = 6;
    elseif x == 0 & y < 0
        thisSubplot = 7;
    elseif x > 0 & y < 0
        thisSubplot = 8;
    end
    
    % Plot scatter and psycometric curve in respective subplots
    subplot(2,4,thisSubplot)  
    title(['Axis: ' num2str(Titles(thisSubplot))])
    scatter(e.d{iFile}.uniquecontrast, e.d{iFile}.correctBinned, color)
    hold on 
    plot(e.d{iFile}.fit.fitX, e.d{iFile}.fit.fitY, color)
    % Calculate error
    errors = [];
    n = 25;
    for i=1:length(e.d{iFile}.correctBinned)
        p = e.d{iFile}.correctBinned(i);
        error = sqrt( (p*(1-p)) / n );
        errors = [errors error];
    end
    errorbar(e.d{iFile}.uniquecontrast, e.d{iFile}.correctBinned, errors, 'LineStyle', 'none', 'Color', color);
    hold on
    
    % Hack in order to have a legend for the graph based on the eccentricity of each curve
    h = zeros(4, 1);
    h(1) = plot(NaN,NaN,'red');
    h(2) = plot(NaN,NaN,'blue');
    h(3) = plot(NaN,NaN,'green');
    h(4) = plot(NaN,NaN,'magenta');
    legend(h, ['0' char(176)], ['2.25' char(176)], ['4.5' char(176)], ['6.75' char(176)], 'Location', 'southeast');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Building the Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1)
% Creating a relationship between the locations (Xs and Ys) and psychometric functions. The goal is to be able to draw a psychometric
% function for any location around the circle, which will be dones by interpolating between the 25 locations we recorded from

% Setting data points to their indeces in the dataMatrix (just makes code more legible)
X = 1;
Y = 2;
Mean = 3;
Std = 4;
threshCon = 5;
% Intializing vectors that will be used for interpolarion
Xs = [];
Ys = [];
Means = [];
Stds = [];
% Populating vectors that will be used for interpolarion
for iCol=1:e.nFiles
    Xs = [Xs dataMatrix(X, iCol)];
    Ys = [Ys dataMatrix(Y, iCol)];
    Means = [Means dataMatrix(Mean, iCol)];
    Stds = [Stds dataMatrix(Std, iCol)];
end
% Convert row vectors to column vectors for the scatterInterpolant() function
Xs = Xs';
Ys = Ys';
Means = Means';
Stds = Stds';
% Interpolating the function linearly
Fmean = scatteredInterpolant(Xs, Ys, Means, 'linear');
Fstd = scatteredInterpolant(Xs, Ys, Stds, 'linear');


% (2)
% Setting up fixed variables
LOCATIONS = 4;
c = 0.2;
alpha = 0.0525;
en = 0.25; 
dPrimeE = sqrt( (c^2) / (alpha*en) );
XmodelStd = 1 / dPrimeE;
priori = 1/LOCATIONS;
priorj = 1/LOCATIONS;

MAP = 1;
criterion = 0.7;

saccades = [0 0;];
fovea = [0 0];

Xidx = 1;
Yidx = 2; 
DPrimeIidx = 3; 
Presentidx = 4; 
XModelidx = 5;
Widx = 6;
Gidx = 7;
Postidx = 8;

% Each location is a rown and each column is a data entry (i.e. the first row has all the information about location 1)
visMatrix = [1.44 1.44 0 -0.5 normrnd(0,XmodelStd) 0 0 0;
             1.43 1.43 0 -0.5 normrnd(0,XmodelStd) 0 0 0;
             1.42 1.42 0 -0.5 normrnd(0,XmodelStd) 0 0 0;
             3.18 -3.18 0 -0.5 normrnd(0,XmodelStd) 0 0 0;];

% Adding 0.5 to the location where the target is present
targetPresent = randi(LOCATIONS);
visMatrix(targetPresent,Presentidx) = 0.5;

T = [];

for Tidx = 1:4
    for Locidx = 1: LOCATIONS
        relativeX = visMatrix(Locidx, Xidx) - fovea(1);
        relativeY = visMatrix(Locidx, Yidx) - fovea(2);
        mu = Fmean(relativeX, relativeY);
        sigma = Fstd(relativeX, relativeY);
        fc = normcdf(c,mu,sigma);
        % This is the dPrime calculated from psychometric data
        dPrimeF = sqrt(2) * norminv(fc);
        beta = ( (c^2) - ((dPrimeF^2)*alpha*en) ) / ( dPrimeF^2 );
        dPrimeI = sqrt( (c^2) / beta );
        visMatrix(Locidx, DPrimeIidx) = dPrimeI;
        Xmodel = visMatrix(Locidx, XModelidx);
        N = normrnd(0,(1/dPrimeI));
        Present = visMatrix(Locidx, Presentidx);
        W = Xmodel + N + Present;
        visMatrix(Locidx, Widx) = W;
        
        % (1) Calculating g
        % (1.1) add the d's for that location across time
        dPrimeIAgg{Locidx}{Tidx} = dPrimeI;
        dPrimeISquareSum = 0;
        for i=1:Tidx
            dPrimeISquareSum = dPrimeISquareSum + (dPrimeIAgg{Locidx}{i})^2;
        end
        
        g= ( dPrimeE * (dPrimeI)^2 ) / ( dPrimeE + dPrimeISquareSum );
        visMatrix(Locidx, Gidx) = g;
        T(:,:,Tidx) = visMatrix;
        
        posterior = calculatePosterior(Locidx, Tidx, T);
        visMatrix(Locidx, Postidx) = posterior;
        
        T(:,:,Tidx) = visMatrix;
    end
    
    [maxPost, rowLocidx] = max(visMatrix(:,Postidx));
    nextX = visMatrix(rowLocidx, Xidx);
    nextY = visMatrix(rowLocidx, Yidx);
    saccade = [nextX nextY];
    saccades(end+1, :) = saccade;
    fovea = saccade;
    
    if maxPost > criterion
        break
    end
end

    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate the Posterior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function posterior = calculatePosterior(Locidx, Tidx, T)
   % i here stands for the location for which you are trying to find the posterior. It is reference by Locidx. 
   % j is an andex for all other locations
   % Tidx stands for time point or fixation. Let us say I am at the 3rd time point (Tidx = 3), then my priori is the prior for location i at time point 3.
   % My giSum is the sum of the g's at location i across the 3 times points that have passed so far

   % (1) Summation over time points for location 
   gxWiSum = 0;
   for tidx = 1:Tidx
       tVisMatrix = T(:,:,tidx);
       giT = tVisMatrix(Locidx, Gidx);
       wiT = tVisMatrix(Locidx, Widx);
       gxWi = giT*wiT;
       gxWiSum = gxWiSum + gxWi;
   end
   numerator = priori * exp(gxWiSum);
   
   %(2) the sume of the priors of all locatiexp(ons at the time point given by Tidx
   sumOverLocs = 0;
   for jLoc = 1:LOCATIONS
       Locidx = jLoc;
       gxWjSum = 0;
       for tidx = 1:Tidx
            tVisMatrix = T(:,:,tidx);
            gjT = tVisMatrix(Locidx, Gidx);
            wjT = tVisMatrix(Locidx, Widx);
            gxWj = gjT*wjT;
            gxWjSum = gxWjSum + gxWj;
       end
       
       oneLoc = priorj * exp(gxWjSum);
       sumOverLocs = sumOverLocs + oneLoc;
   end
  
   posterior = numerator / sumOverLocs;
   
end

k=2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating Statistics on dataMatrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
% Setting data points to their indeces in the dataMatrix (just makes code more legible)
X = 1;
Y = 2;
Mean = 3;
Std = 4;
threshCon = 5;
%}


%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1)
% Making a 3D graph for the linear interpolation of the means and stds of the 25 psychometric curves as a function of (x,y) (i.e location)

Xs = [];
Ys = [];
Means = [];
Stds = [];

for iCol=1:e.nFiles
    Xs = [Xs dataMatrix(X, iCol)];
    Ys = [Ys dataMatrix(Y, iCol)];
    Means = [Means dataMatrix(Mean, iCol)];
    Stds = [Stds dataMatrix(Std, iCol)];
end

% Convert row vectors to column vectors for the scatterInterpolant() function
Xs = Xs';
Ys = Ys';
Means = Means';
Stds = Stds';

% Interpolating the function linearly
Fm = scatteredInterpolant(Xs, Ys, Means, 'linear');
Fs = scatteredInterpolant(Xs, Ys, Stds, 'linear');

% Creating a meshgrid for plotting
[X, Y] = meshgrid(Xs, Ys);
Vm = Fm(X,Y);
Vs = Fs (X,Y);

% Plotting
figure(1)
mesh(X,Y,Vm);
hold on
plot3(X, Y, Vm);
title('Mean')
xlabel('X'); ylabel('Y');

figure(2)
mesh(X,Y,Vs);
hold on
plot3(X, Y, Vs, 'o');
title('Std')
xlabel('X'); ylabel('Y');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}


%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2)
% Making graphs of the average thresholdContrast, std, and mean for each eccentricty 

% initialize arrays that will hold data points for each eccentricity (regardless of location)
threshCons0 = []; stds0 = []; means0 = [];
threshCons225 = []; stds225 = []; means225 = [];
threshCons45 = []; stds45 = []; means45 = [];
threshCons675 = []; stds675 = []; means675 = [];

% iterate through dataMatrix and pull out the thresholdContrasts, stds, and means based on eccentricity
for iCol=1:e.nFiles
    
    eccen = sqrt( (dataMatrix(X, iCol))^2 + (dataMatrix(Y, iCol))^2 );
    
    if  eccen == 0
        threshCons0 = [threshCons0 dataMatrix(threshCon, iCol)];
        stds0 = [stds0 dataMatrix(Std, iCol)];
        means0 = [means0 dataMatrix(Mean, iCol)];
    end
    
    if eccen > 1 & eccen < 3
        threshCons225 = [threshCons225 dataMatrix(threshCon, iCol)];
        stds225 = [stds225 dataMatrix(Std, iCol)];
        means225 = [means225 dataMatrix(Mean, iCol)];
    end
    
    if eccen > 4 & eccen < 5
        threshCons45 = [threshCons45 dataMatrix(threshCon, iCol)];
        stds45 = [stds45 dataMatrix(Std, iCol)];
        means45 = [means45 dataMatrix(Mean, iCol)];
    end
    
    if eccen > 6
        threshCons675 = [threshCons675 dataMatrix(threshCon, iCol)];
        stds675 = [stds675 dataMatrix(Std, iCol)];
        means675 = [means675 dataMatrix(Mean, iCol)];
    end
end

% Average thresholdContrast per eccentricity
AvgthreshCon0 = sum(threshCons0) / length(threshCons0);
AvgthreshCon225 = sum(threshCons225) / length(threshCons225);
AvgthreshCon45 = sum(threshCons45) / length(threshCons45);
AvgthreshCon675 = sum(threshCons675) / length(threshCons675);
% Average standard deviations of the cumalitive gaussian fits per eccentricity
Avgstd0 = sum(stds0) / length(stds0);
Avgstd225 = sum(stds225) / length(stds225);
Avgstd45 = sum(stds45) / length(stds45);
Avgstd675 = sum(stds675) / length(stds675);
% Average means of the cumalitive gaussian fits per eccentricity
Avgmean0 = sum(means0) / length(means0)
Avgmean225 = sum(means225) / length(means225)
Avgmean45 = sum(means45) / length(means45)
Avgmean675 = sum(means675) / length(means675)

% Creating the X-axis ange for plotting (the range of contrasts)
x = 0:0.01:0.3;
% Use the function normcdft (what we used to fit the data) to plot the psychometric functions (using the parameters mu and sigma)
plot(normcdf(x, Avgmean0, Avgstd0), 'DisplayName', '0')
hold on
plot(normcdf(x, Avgmean225, Avgstd225),'DisplayName', '2.25')
hold on
plot(normcdf(x, Avgmean45, Avgstd45),'DisplayName', '4.5')
hold on
plot(normcdf(x, Avgmean675, Avgstd675),'DisplayName', '6.75')
hold on
    
legend show ;
legend('Location', 'southeast');
title('Avg Psychometric Curve per Eccentricity')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}




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

end

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
end

end

