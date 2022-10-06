function [x,y] = screen2deg(x_screen, y_screen, myscreen)
    x = (x_screen-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
    y = (y_screen-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;
end
