%% Paperwork, payment form, and dominant eye

%% Set Subject Info
sid = 300;

%% Initialize Commands
% cd ~/proj/gru
% startup
% 
% open berlin_dms

%% Set SID
mglSetSID(sid);

%% Explanation (no eye tracker)
% Tell subject how experiment works:
%  - Fixate when the cross is present
%  - Memorize direction of first dots
%  - Respond counteclockwise / clockwise
%  - Subject does nothing--you tell them what they would do
berlin_dms('feature=2');

%% First practice run
% Tell subject to press 1 or 2 for counter/clockwise
berlin_dms('feature=2');

%% Run Grating Block
% First eyetracker run:
%  - Important: don't move head once calibrated
%  - Look at red dots, make sure they see the white center

aaac
aaac


a
c


% press ESC to accept all calibration

%% Run Motion Block
% Second stimulus:
%  - Same exact task, but using gratings that don't move
berlin_dms('feature=1','trigger=1');

%% Reset Data + Re-open MATLAB
% Copy subjects data folders to "training" folder

%% Run 2x motion
berlin_dms('feature=1','trigger=1');
%% Run 2x grating (then repeat motion)
berlin_dms('feature=2','trigger=1');

%% End - cleanup: run berlin_copy.m
berlin_copy;
mglClose