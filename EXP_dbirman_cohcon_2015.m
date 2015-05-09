%% Intro

% This script runs the coherentContrast experiment
% This revision is for May 2015, data collection for Dan's FYP project.

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

%% MT Localizer

% The 32ch coil should be installed.
% disp(sprintf('Running MT Localizer.'));
% while ~strcmp(input('Have you prepped the mux8 script? [y] ','s'),'y')
% end
% 
% % mtLoc x2
% % mtloc('0%',1.4);
% mtloc('0%',.75);

%% Full Experiment

% The 32ch coil should be installed.

% Hi Steeve. Did you run the in plane anatomy yet? Do that. Also run the
% CAL sequence!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IN PLANE ANATOMY FIRST!! %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Then, run these first:
% You should be using the 2.5mm mux8 sequence. Total volumes = 840 (7:00).
%%
cohcon('plots=0','scan=1','nocatch=1');
%%
cohcon('plots=0','scan=1','nocatch=1');
%%
cohcon('plots=0','scan=1','nocatch=1');
%%
cohcon('plots=0','scan=1','nocatch=1');

% Then run four of these:
%%
cohcon('plots=0','scan=1');
%%
cohcon('plots=0','scan=1');
%%
cohcon('plots=0','scan=1');
%%
cohcon('plots=0','scan=1');

% Then run four more of those first ones:
%%
cohcon('plots=0','scan=1','nocatch=1');
%%
cohcon('plots=0','scan=1','nocatch=1');
%%
cohcon('plots=0','scan=1','nocatch=1');
%%
cohcon('plots=0','scan=1','nocatch=1');

% Then two more of these:
%%
cohcon('plots=0','scan=1');
%%
cohcon('plots=0','scan=1');

% Now just run as many of these as possible until I collapse and die.
% Please give me breaks once in a while.
%%
cohcon('plots=0','scan=1','nocatch=1');