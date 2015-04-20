%% Intro

% This script runs the coherentContrast experiment
% This revision is for March 2015, data collection for Dan's FYP project.


%% Retinotopy + MTLoc

% The 16ch head coil should be installed and you should have already
% collected an in-plane anatomy.

% % Retinotopy WEDGES 1:
% mglRetinotopy('displayName=fMRIproj32','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',17,'blanks=3','doEyeCalib=0','wedges=1','direction=-1');
% mglRetinotopy('displayName=fMRIproj32','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',17,'blanks=3','doEyeCalib=0','wedges=1','direction=1');
% 
% % Retinotopy RINGS:
% mglRetinotopy('displayName=fMRIproj32','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',17,'blanks=3','doEyeCalib=0','rings=1','direction=-1');
% mglRetinotopy('displayName=fMRIproj32','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',17,'blanks=3','doEyeCalib=0','rings=1','direction=1');
% 
% % Retinotopy WEDGES 2:
% mglRetinotopy('displayName=fMRIproj32','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',17,'blanks=3','doEyeCalib=0','wedges=1','direction=-1');
% mglRetinotopy('displayName=fMRIproj32','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',17,'blanks=3','doEyeCalib=0','wedges=1','direction=1');
% 
% % Retinotopy BARS x4:
% mglRetinotopy('displayName=fMRIproj32','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',17,'blanks=3','doEyeCalib=0','bars=1');

%% Unattended Task

% The 32ch coil should be installed.

% unattended x4-6
% cohcon('unattended=1','plots=0','scan=1');

% After running the unattended task collect an in-plane anatomical OR the
% canonical anatomical;

%% MT Localizer

% The 32ch coil should be installed.
disp(sprintf('Running MT Localizer.'));
while ~strcmp(input('Have you prepped the mux8 script? [y] ','s'),'y')
end

% mtLoc x2
% mtloc('0%',1.4);
mtloc('0%',.75);

%% Full Experiment

% The 32ch coil should be installed.

% exp script
cohcon('plots=0','scan=1');

% After running the experiment collect an in-plane anatomical OR the
% canonical anatomical;