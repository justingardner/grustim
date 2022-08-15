function stimulus = calibrateJoy(myscreen, stimulus)
    %CALIBRATEJOY Summary of this function goes here
    %   Detailed explanation goes here
    
    disp(' (calibrateJoy) calibrating joystick ...');

    stimulus.joy = vrjoystick(1); % use simulink 3d animation to load joystick object
    if isempty(stimulus.joy)
        stimulus.exp.useJoystick = 0;
        exp = stimulus.exp;
        disp(' FAILED TO FIND JOYSTICK! MAKE SURE SIMULINK 3D ANIMATION PACKAGE IS INSTALLED AND THE JOYSTICK IS PROPERLY CONNECTED');
        disp(' USING MOUSE FOR TRACKING ...');
    else
        joy_params              = struct();
        joy_params.maxv         = 0.2;
        joy_params.deadzone     = 0.02;
        joy_params.sensitivity  = 1;
        joy_params.poly_order   = 0.8;
        joy_params.theta_range  = pi/6;
        joy_params.xyratio      = 1.1;
        joy_params.bias         = [0,0,0];
    end
    
    stimulus.pointer = [1,1];
    
    mglClearScreen;
    mglTextDraw('Calibrating joystick...  Please let go of the joystick.',[0 0]);
    mglFlush;
    myscreen = tickScreen(myscreen,[]); 
    a = axis(stimulus.joy);
    n=0; % need to be still for at least 100 frames
    while n < 100
        if all(a == axis(stimulus.joy))
            n = n+1;
            joy_params.bias = a(1:3);
        else
            n=0; 
            mglClearScreen;
            mglTextDraw('Please let go of the joystick....',[0 0]);
            a = axis(stimulus.joy);
        end
        myscreen = tickScreen(myscreen,[]); 
    end
    
    mglClearScreen;
    mglTextDraw('Press <space> to continue.',[0 0]);
    mglFlush;
    [keyCodes keyTimes] = mglGetKeyEvent([],1);
    while ~any(keyCodes==myscreen.keyboard.space)
        if any(keyCodes == myscreen.keyboard.esc)
          mglClearScreen; myscreen = tickScreen(myscreen,[]);
          return
        end
        [keyCodes keyTimes] = mglGetKeyEvent(0.1,1);
    end
    
    screencleartime = 5; % secs; clear screen every x seconds
    mglClearScreen
    tic
    while ~myscreen.userHitEsc  
        % move cursor        
        [vx, vy] = joy2vel(stimulus.joy, joy_params, myscreen);
        stimulus = update_pointer(stimulus, [vx, vy], myscreen);

        % display cursor
        mglGluDisk(stimulus.pointer(1), stimulus.pointer(2), 0.2, [1 0 0])
    
        % update parameters
        jp = joy_params;
        buttonstate = button(stimulus.joy);
        if buttonstate(3) == 1, jp.poly_order = max(jp.poly_order*0.98,0.05) ;,end
        if buttonstate(5) == 1, jp.poly_order = jp.poly_order*1.02 ;,end
        if buttonstate(4) == 1, jp.sensitivity = max(jp.sensitivity - 0.1,0) ;,end
        if buttonstate(6) == 1, jp.sensitivity = jp.sensitivity + 0.1 ;,end
        if buttonstate(7) == 1, jp.deadzone = jp.deadzone*0.95;,end
        if buttonstate(8) == 1, jp.deadzone = jp.deadzone*1.05;,end
        if buttonstate(9) == 1, jp.xyratio = jp.xyratio*0.99;,end
        if buttonstate(10) == 1, jp.xyratio = jp.xyratio*1.01;,end
        if buttonstate(11) == 1, jp.maxv = jp.maxv*0.95 ;,end
        if buttonstate(12) == 1, jp.maxv = jp.maxv*1.05 ;,end

        joy_params = jp;
        
        % clear screen
        if toc > screencleartime
            mglClearScreen
            tic
            joy_params;
        end
        
        myscreen = tickScreen(myscreen,[]);     % flip screen
    end

    stimulus.joy_params     = joy_params;
    mglClearScreen; mglFlush;
end

%% joystick update
function stimulus = update_pointer(stimulus, vel, myscreen)
    pos = stimulus.pointer;
    
    [horz_out, vert_out] = check_oob(pos + vel, myscreen, stimulus);
    stimulus.pointer(1) = stimulus.pointer(1) + (1-horz_out)*vel(1);
    stimulus.pointer(2) = stimulus.pointer(2) + (1-vert_out)*vel(2);
end