function stimulus = joy2vel(joy, stimulus, myscreen, task)
    % convert the axes of the joystick to 
    
    a = axis(joy);
    vx = ;
    vy = ;
    
    stimulus.pointer(1) = stimulus.pointer(1) + vx / myscreen.framesPerSecond;
    stimulus.pointer(2) = stimulus.pointer(2) + vy /  myscreen.framesPerSecond;

end

function
end
