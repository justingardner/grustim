%        $Id: $
%      usage: 
%         by: Josh Ryu
%       date: 03/07/2023
%    purpose: 

% Tracking task for circular motion around an circle/ellipse
% Run in mglmetal 

function stimulus = calibrateWheel(myscreen, stimulus)

    disp(' (calibrateWheel) calibrating powerwheel ...');

    screencleartime = 5; % secs; clear screen every x seconds

    order_init          = 1;
    ecc_r               = 5;
    pointer_r           = 0.2;
    
    maxorder                = 4;
    wheel_params            = ornstein_uhlenbeck(maxorder, 0);
    wheel_params.maxorder   = order_init;
    state                   = zeros(maxorder,1);

    step    = 0.01;
    dt      = 1/myscreen.framesPerSecond;

    % set mouse position
    [x_screen,y_screen] = deg2screen(0, 0, myscreen);
    mglSetMousePosition(ceil(x_screen),floor(y_screen), myscreen.screenNumber);

    % show instructions
    mglClearScreen;
    mglTextDraw('Calibrating powerwheel...',[3 0]);
    mglTextDraw('Press 0/1/2/3/4 to adjust the order of dynamics (0:pos, 1:vel, etc..) ...',[2 0]);
    mglTextDraw('Press up/down to adjust decay constant...',[1 0]);
    mglTextDraw('Press left/right to adjust integration constant...',[0 0]);
    mglTextDraw('Press "<" or ">" to adjust step for the adjustments...',[-1 0]);
    mglTextDraw('Press <space> to continue...',[-3 0]);
     
    mglFlush;
    while ~any(keyCodes==myscreen.keyboard.space)
        if any(keyCodes == myscreen.keyboard.esc)
          mglClearScreen; myscreen = tickScreen(myscreen,[]);
          return
        end
        [keyCodes keyTimes] = mglGetKeyEvent(0.1,1);
    end

    % run adjustment until esc
    mglClearScreen
    tic
    while ~myscreen.userHitEsc  
        % display trajectory ring
        mglMetalArcs([0;0;0], [0; 0; 1; 1], [ecc_r-0.1; ecc_r+0.1], [0;2*pi], 1);

        % update parameters
        keystate = mglGetKeys;
        if any(keystate([30,19:21])) %keys: 0,1,2,3
            orderarray  = [0,1,2,3];
            orders_pressed = orderarray(keystate([29,18:20]));
            newmaxorder = orders_pressed(1);
            if newmaxorder ~= wheel_params.maxorder
                wheel_params.maxorder = newmaxorder;
            end
        end
        
        if keystate(127) % up key => less decay
            wheel_params.(['invtaus_decay', num2str(wheel_params.maxorder)]) = ...
                update_param(wheel_params.(['invtaus_decay', num2str(wheel_params.maxorder)]), 0, 1/dt, -1* step);
        elseif keystate(126) % down key => more decay
            wheel_params.(['invtaus_decay', num2str(wheel_params.maxorder)]) = ...
                update_param(wheel_params.(['invtaus_decay', num2str(wheel_params.maxorder)]), 0, 1/dt, step);
        end

        if keystate(124) % left key => instantaneous integration
            wheel_params.(['taus_int', num2str(wheel_params.maxorder)]) = ...
                update_param(wheel_params.(['taus_int', num2str(wheel_params.maxorder)]), dt, 10, -1*step);
        elseif keystate(125) % right key => longer integration
            wheel_params.(['taus_int', num2str(wheel_params.maxorder)]) = ...
                update_param(wheel_params.(['taus_int', num2str(wheel_params.maxorder)]), dt, 10, step);
        end

        if keystate(44) % < decrease step size
            step = update_param_scale(step, minval, maxval, 0.95);
        elseif keystate(48) %> increase step size
            step = update_param_scale(step, minval, maxval, 1.05);
        end

        % move cursor        
        mInfo = mglGetMouse(myscreen.screenNumber);
        [ux,uy] = screen2deg(mInfo.x, mInfo.y, myscreen);
                
        state = ou_update_state(state, ux, wheel_params, dt);
        [x_screen,y_screen] = deg2screen(0, 0, myscreen);
        mglSetMousePosition(ceil(x_screen),floor(y_screen), myscreen.screenNumber);

        % display cursor
        theta = stimulus.pointer.state(1);
        mglMetalArcs([ecc_r * cos(theta); ecc_r * sin(theta); 0], ...
            [0; 0; 1; 1], [ecc_r-0.1; ecc_r+0.1], [0;2*pi], 1);
        
        % clear screen
        if toc > screencleartime
            mglClearScreen
            tic
            joy_params;
        end
        
        myscreen = tickScreen(myscreen,[]);     % flip screen
    end

    stimulus.wheel_params = wheel_params;
    mglClearScreen; mglFlush;

    disp(' (calibrateWheel) finished calibration.');
end

function pval = update_param(pval, minval, maxval, step)
    pval = pval + step;
    pval = min(maxval, pval);
    pval = max(minval, pval);
end

function pval = update_param_scale(pval, minval, maxval, step)
    pval = pval*step;
    pval = min(maxval, pval);
    pval = max(minval, pval);
end