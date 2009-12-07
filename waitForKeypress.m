% waitForKeypress - wait for key press using mgl functions
%
%      usage: [ key ] = waitForKeypress
%         by: denis schluppeck
%       date: 2007-01-16
%        $Id: waitForKeypress.m,v 1.1 2007/01/16 20:58:29 ds Exp $:
%     inputs: 
%    outputs: key - first key pressed [optional]
%
%    purpose: simple wrapper for mglGetKeys
%
%        e.g: waitForKeypress
%
function [  ]=waitForKeypress(  )

while (1); k=mglGetKeys; if (any(k)),break;end;end
keycode = find(k);
key = mglKeycodeToChar(keycode);


