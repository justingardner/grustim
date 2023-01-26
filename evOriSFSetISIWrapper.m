% evOriSFSetISIWrapper.m - wrapper script to feed evOriSFSetISI the proper stim/blank structure

trialtiming = createNSDTrialTiming();

for i = 1:size(trialtiming,2)
    evOriSFSetISI(trialtiming(:,i));
end