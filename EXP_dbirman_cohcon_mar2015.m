%% Intro

% This script runs the coherentContrast experiment
% This revision is for March 2015, data collection for Dan's FYP project.

proj = 'projector=0';

%% Retinotopy + MTLoc

% The 16ch head coil should be installed and you should have already
% collected an in-plane anatomy.

% Retinotopy RINGS 1:
mglRetinotopy('displayName=fMRI_proj','doEyeCalib=0','rings=1','direction=-1');
mglRetinotopy('displayName=fMRI_proj','doEyeCalib=0','rings=1','direction=1');

% Retinotopy WEDGES 2:
mglRetinotopy('displayName=fMRI_proj','doEyeCalib=0','wedges=1','direction=-1');
mglRetinotopy('displayName=fMRI_proj','doEyeCalib=0','wedges=1','direction=1');
mglRetinotopy('displayName=fMRI_proj','doEyeCalib=0','wedges=1','direction=-1');
mglRetinotopy('displayName=fMRI_proj','doEyeCalib=0','wedges=1','direction=1');

% Retinotopy BARS 3:
mglRetinotopy('displayName=fMRI_proj','doEyeCalib=0','bars=1');
mglRetinotopy('displayName=fMRI_proj','doEyeCalib=0','bars=1');
mglRetinotopy('displayName=fMRI_proj','doEyeCalib=0','bars=1');
mglRetinotopy('displayName=fMRI_proj','doEyeCalib=0','bars=1');

%% Unattended Task

% The 32ch coil should be installed.

% mtLoc
coherentContrast('unattended=1','plots=0','mtloc=1',proj);

% unattended
coherentContrast('unattended=1','plots=0',proj);
coherentContrast('unattended=1','plots=0',proj);
coherentContrast('unattended=1','plots=0',proj);
coherentContrast('unattended=1','plots=0',proj);

% After running the unattended task collect an in-plane anatomical OR the
% canonical anatomical (time dependent);

%% Full Experiment

% The 32ch coil should be installed.

% exp script
coherentContrast(proj);

% After running the experiment collect an in-plane anatomical OR the
% canonical anatomical (time dependent);