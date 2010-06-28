% gruDispText.m
%
%        $Id:$ 
%      usage: gruDispText()
%         by: justin gardner
%       date: 06/28/10
%    purpose: 
%
function retval = gruDispText()

% check arguments
if ~any(nargin == [0])
  help gruDispText
  return
end

mglDispText([],'Press button once for YES twice for NO');
mglClose;
