%% Stanford Retinotopy

% A script for a 'standard' retinotopy at the CNI scanner.

%% Wedges 1
mglRetinotopy('displayName=fMRIproj16','wedges=1','direction=1','numCycles=10','stimulusPeriod=24','stepsPerCycle',17,'doEyeCalib=-1','initialHalfCycle=1');
%%
mglRetinotopy('displayName=fMRIproj16','wedges=1','direction=-1','numCycles=10','stimulusPeriod=24','stepsPerCycle',17,'doEyeCalib=-1','initialHalfCycle=1');

%% Rings
mglRetinotopy('displayName=fMRIproj16','rings=1','direction=1','numCycles=10','stimulusPeriod=24','stepsPerCycle',17,'doEyeCalib=-1','initialHalfCycle=1');
%%
mglRetinotopy('displayName=fMRIproj16','rings=1','direction=-1','numCycles=10','stimulusPeriod=24','stepsPerCycle',17,'doEyeCalib=-1','initialHalfCycle=1');

%% Wedges 2
mglRetinotopy('displayName=fMRIproj16','wedges=1','direction=1','numCycles=10','stimulusPeriod=24','stepsPerCycle',17,'doEyeCalib=-1','initialHalfCycle=1');
%%
mglRetinotopy('displayName=fMRIproj16','wedges=1','direction=-1','numCycles=10','stimulusPeriod=24','stepsPerCycle',17,'doEyeCalib=-1','initialHalfCycle=1');

%% Bars

mglRetinotopy('displayName=fMRIproj16','bars=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',17,'blanks=3','doEyeCalib=-1');
%%
mglRetinotopy('displayName=fMRIproj16','bars=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',17,'blanks=3','doEyeCalib=-1');
%%
mglRetinotopy('displayName=fMRIproj16','bars=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',17,'blanks=3','doEyeCalib=-1');
%%
mglRetinotopy('displayName=fMRIproj16','bars=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',17,'blanks=3','doEyeCalib=-1');