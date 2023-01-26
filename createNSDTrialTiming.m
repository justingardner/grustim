% createNSDTrialTiming.m - convert trial structure of inserting blanks found in NSD experiment presentation to trial-by-trial timescale
%
% timing: 75 trials - trials 1-3: blank
%                     trials 72-75: blank
%                     5 intermediate blanks (between 9-14 non-blank trials between each intermediate blank)

function [trialTiming] = createNSDTrialTiming()

mydir = '~/Documents/MATLAB/data/nsddesign';

filePattern = fullfile(mydir, '*.tsv');
theFiles = dir(filePattern);

for i = 1:length(theFiles)
    tdfread(strcat(theFiles(i).folder,'/',theFiles(i).name));
    trialTiming(:,i) = convertNSDTrialStruct(x0);
end

    %% %%%%%%%%%%%%%% %%
    % helper functions %
    %%%%%%%%%%%%%%%%%%%%

    function [trialVec] = convertNSDTrialStruct(NSD_vec)
        
        removeSupersampling = [0;NSD_vec(3:3:end,:)];
        trialVec = removeSupersampling > 0;
        
    end

end

