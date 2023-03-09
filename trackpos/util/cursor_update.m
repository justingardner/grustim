function [dx, dy, mousestate] = cursor_update(myscreen,mousestate)
% calculate how much the mouse moved
% if the mouse went too far from center, then set it to the center again

% This funtion should be used when you only want to look at mouse
% velocities but can't reset the mouse every frame because of the mouse
% polling rate being smaller than the screen refresh rate

    thresh = 0.1;

    mInfo = mglGetMouse(myscreen.screenNumber);
    [x,y] = screen2deg(mInfo.x, mInfo.y, myscreen);

    dx = x - mousestate(1);
    dy = y - mousestate(2);
        
    % reset mouse position
    if mInfo.x < myscreen.screenWidth * thresh || mInfo.x > myscreen.screenWidth * (1-thresh) ...
            || mInfo.y <  myscreen.screenHeight * thresh || mInfo.y > myscreen.screenHeight *(1-thresh)
        
        [x_screen,y_screen] = deg2screen(0, 0, myscreen);
        mglSetMousePosition(ceil(x_screen),floor(y_screen), myscreen.screenNumber);
        
        mousestate = [0,0];
    else
        mousestate = [x,y];
    end

end

