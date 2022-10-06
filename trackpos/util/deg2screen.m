function [x_screen,y_screen] = deg2screen(x, y, myscreen)
    x_screen = x*myscreen.screenWidth/myscreen.imageWidth + myscreen.screenWidth/2;
    y_screen = y*myscreen.screenHeight/myscreen.imageHeight + myscreen.screenHeight/2;
end