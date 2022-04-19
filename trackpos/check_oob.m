function [horz_out, vert_out] = check_oob(pos, myscreen, stimulus)    
    stimstd = stimulus.stimStd;
    xBound = myscreen.imageWidth/2;
    yBound = myscreen.imageHeight/2;
    
    horz_out = false;
    vert_out = false;
    
    % if thre circle goes out of the image 
    if pos(1)+3*stimstd > xBound || pos(1)-3*stimstd < (-1 * xBound),
        horz_out = true;
    end
        
    if pos(2)+3*stimstd > yBound || pos(2)-3*stimstd < (-1 * yBound),
        vert_out = true;
    end
end