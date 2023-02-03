% evOriSFSetISIWrapper.m - wrapper script to feed evOriSFSetISI the proper stim/blank structure

clear

trialtiming = createNSDTrialTiming();

% add on some safety trials so it won't error out right after 75 trials
trialtiming(end+1:end+10,:) = 0;
mglSetSID(-1);
for i = 1:size(trialtiming,2)
    evOriSFSetISI(trialtiming(:,i),'atScanner=0');
    nextRunMsg = sprintf('\nRun %d completed - hit return to continue...\n',i);
    input(nextRunMsg);
end