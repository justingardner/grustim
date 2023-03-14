% evOriSFPhEnc_oriWrapper.m - wrapper script to feed evOriSFPhEnc the proper stim parameters

clear all

% easy access parameters for scanning/testing
saveParam = 'saveParam=0';
scanner = 'atScanner=0';
mglSetSID(-1);
setRepeats = 6; % number of repeats of all conditions

% loose variables - MAKE SURE THESE MATCH PARAMETERS IN evOriSFPhEnc
sf = [0.5]; % 1];
dir = [1 -1];

for setIdx = 1:setRepeats
    runCounter = 1;
    for i = 1:length(sf)
        for j = 1:length(dir)
            sfParam = strcat('sf=',num2str(sf(i)));
            dirParam = strcat('oridir=',num2str(dir(j)));
            evOriSFPhEnc_ori(sfParam,dirParam,scanner,saveParam);
            
            nextRunMsg = sprintf('\nRun %d of %d (current set: %d of %d) completed - hit return to continue...\n',...
                         runCounter,length(sf)*length(dir),setIdx,setRepeats);
            input(nextRunMsg);
            runCounter = runCounter + 1;
        end
    end
end