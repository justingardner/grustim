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
    use_numeric_keyboard = true;
    
    order_init          = 0;
    ecc_r               = 5;
    pointer_r           = 0.2;    
    
    maxorder                = 4;
    wheel_params            = ornstein_uhlenbeck(maxorder, 0);
    wheel_params.maxorder   = order_init;
    state                   = zeros(maxorder,1);

    step    = 0.01;
    dt      = 1/myscreen.framesPerSecond;
    mouse0  = [0,0];
    


    % set mouse position
    [x_screen,y_screen] = deg2screen(0, 0, myscreen);
    mglSetMousePosition(ceil(x_screen),floor(y_screen), myscreen.screenNumber);

    % show instructions
    mglClearScreen;
    mglTextDraw('Calibrating powerwheel...',[0 3]);
    mglTextDraw('Press 0/1/2/3/4 to adjust the order of dynamics (0:pos, 1:vel, etc..) ...',[0 2]);
    mglTextDraw('Press up/down to adjust decay constant...',[0 1]);
    mglTextDraw('Press left/right to adjust integration constant...',[0 0]);
    mglTextDraw('Press "<" or ">" to adjust step for the adjustments...',[0 -1]);
    mglTextDraw('Press <space> to continue...',[0 -3]);
     
    mglFlush;
    [keyCodes keyTimes] = mglGetKeyEvent(0.1,1);
    while ~any(keyCodes==myscreen.keyboard.space)
        if any(keyCodes == myscreen.keyboard.esc)
          mglClearScreen; myscreen = tickScreen(myscreen,[]);
          return
        end
        [keyCodes keyTimes] = mglGetKeyEvent(0.1,1);
    end

    % run adjustment until esc
    mglClearScreen;
    tic;
    while ~myscreen.userHitEsc  
        
        % display blue trajectory ring and fixation
        mglMetalArcs([0;0;0], [0; 0; 1; 1], [ecc_r-0.1; ecc_r+0.1], [0;2*pi], 1);
        mglMetalArcs([0;0;0], [1;1;1; 1], [pointer_r+0.1;pointer_r+0.3],[0;2*pi], 1);
        mglMetalDots([0;0;0], [0.5+0.5*rand(3,1);1], [pointer_r;pointer_r], 1, 1);

        % update parameters
        keystate = mglGetKeys;
        
        if use_numeric_keyboard
            codes_0123 = 83:86;
        else
            codes_0123 = [30,19:21];
        end

        if any(keystate(codes_0123)) %keys: 0,1,2,3
            orderarray  = [0,1,2,3];
            orders_pressed = orderarray(keystate(codes_0123));
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
        [ux, uy, mouse0] = cursor_update(myscreen,mouse0);
        state = ou_update_state(state, -1*ux/ ecc_r , wheel_params, dt);

        % display cursor
        theta = state(1);
        mglMetalDots([ecc_r * cos(theta); ecc_r * sin(theta); 0], ...
                     [1;0;0;1], [pointer_r;pointer_r], 1, 1);
        
        fprintf('r = %0.3e; ux = %0.5e; ',ecc_r, ux);
        fprintf('theta = %0.5e; order= %i; decay= %0.3e; int= %0.3e \n',...
            theta, wheel_params.maxorder,...
            wheel_params.(['invtaus_decay', num2str(wheel_params.maxorder)]),...
            wheel_params.(['taus_int', num2str(wheel_params.maxorder)]));
        
        % clear screen
        if toc > screencleartime
            mglClearScreen;
            tic;
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