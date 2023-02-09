% evOriSFPhEncWrapper.m - wrapper script to feed evOriSFPhEnc the proper stim parameters

clear all

% easy access parameters for scanning/testing
saveParam = 'saveParam=1';
scanner = 'atScanner=1';
mglSetSID('s601');
setRepeats = 3; % number of repeats of all conditions

% loose variables - MAKE SURE THESE MATCH PARAMETERS IN evOriSFPhEnc
ori = [0 90];
dir = [1 -1];

for setIdx = 1:setRepeats
    runCounter = 1;
    for i = 1:length(ori)
        for j = 1:length(dir)
            oriParam = strcat('ori=',num2str(ori(i)));
            dirParam = strcat('sfdir=',num2str(dir(j)));
            evOriSFPhEnc(oriParam,dirParam,scanner,saveParam);
            
            nextRunMsg = sprintf('\nRun %d of %d (current set: %d of %d) completed - hit return to continue...\n',...
                         runCounter,length(ori)*length(dir),setIdx,setRepeats);
            input(nextRunMsg);
            runCounter = runCounter + 1;
        end
    end
end