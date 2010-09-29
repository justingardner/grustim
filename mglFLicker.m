% mglFlicker.m
%
%        $Id:$ 
%      usage: mglFlicker()
%         by: justin gardner
%       date: 09/09/10
%    purpose: 
%
function retval = mglFlicker(varargin)

% check arguments
if ~any(nargin == [0 1 2 3 4 ])
  help mglFlicker
  return
end

freq = [];screenNumber = [];
getArgs(varargin,{'freq=[]','screenNumber=[]'});

% get the frame time
mglOpen(screenNumber);
frameTime = 1/mglGetParam('frameRate');

% handle specification of frequency
if isempty(freq)
  waitFrames = 1;
else
  waitFrames = round((1/freq)/frameTime);
end

% display settings
disp(sprintf('(mglFlicker) Testing frequency %f: %i frames',1/(waitFrames*frameTime),waitFrames));
disp(sprintf('(mglFlicker) Hit any key to end'));

% clear keyboard buffer
mglGetKeyEvent(0,1);

% init variables
currentColor= 1;
mglClearScreen(currentColor);
mglFlush;

% flicker screen
while isempty(mglGetKeyEvent)
  for i = 1:waitFrames
    mglClearScreen(double(currentColor));
    mglFlush;
  end
  currentColor = ~currentColor;
end
mglClose;

