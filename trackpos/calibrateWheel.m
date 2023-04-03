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
    
    if isfield(stimulus, 'pointerR')
        pointer_r = stimulus.pointerR;
    else
        pointer_r           = 0.4;
    end
    
    wordy = false;
    
    maxorder                = 4;
    wheel_params            = ornstein_uhlenbeck(maxorder, 0);
    wheel_params.maxorder   = order_init;
    
    state                   = zeros(maxorder,1);

    step    = 0.01;
    dt      = 1/myscreen.framesPerSecond;
    wheel_params.taus_int0 = dt;
    mouse0  = [0,0];

    % reset if wheel_params is specified
    if isfield(stimulus,'wheel_params')
        wheel_params = stimulus.wheel_params;
    end

    % set mouse position
    [x_screen,y_screen] = deg2screen(0, 0, myscreen);
    mglSetMousePosition(ceil(x_screen),floor(y_screen), myscreen.screenNumber);

    % show instructions
    mglClearScreen;
    mglTextDraw('Calibrating powerwheel...',[0 3]);
    % mglTextDraw('Press 0/1/2/3/4 to adjust the order of dynamics (0:pos, 1:vel, etc..) ...',[0 2]);
    % mglTextDraw('Press up/down to adjust decay constant...',[0 1]);
    mglTextDraw('Press left/right to adjust integration constant...',[0 0]);
    mglTextDraw('Press up/down to adjust sensitivity of the adjustments...',[0 -1]);
    mglTextDraw('Press <space> to continue...',[0 -3]);
    mglFlush;
%    [keyCodes, ~] = mglGetKeyEvent(0.1,1);
    keystate = mglGetKeys;
    while ~keystate(50) && ~keystate(13)
        if keystate(13)
          mglClearScreen; myscreen = tickScreen(myscreen,[]);
          return
        end
        keystate = mglGetKeys;
        % [keyCodes, ~] = mglGetKeyEvent(0.1,1);
    end

    % run adjustment until esc
    mglClearScreen;
    keystate = mglGetKeys;
    tic; str2print = '';
    while ~keystate(13)
        % display blue trajectory ring and fixation
        mglMetalRing_wlines(ecc_r, pointer_r/2, [0,0,1], 600)
        
        % mglMetalArcs([0;0;0], [0; 0; 1; 1], [ecc_r-0.1; ecc_r+0.1], [0;2*pi], 1);
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
        
%         if keystate(127) % up key => less decay
%             wheel_params.(['invtaus_decay', num2str(wheel_params.maxorder)]) = ...
%                 update_param(wheel_params.(['invtaus_decay', num2str(wheel_params.maxorder)]), 0, 1/dt, -1* step);
%         elseif keystate(126) % down key => more decay
%             wheel_params.(['invtaus_decay', num2str(wheel_params.maxorder)]) = ...
%                 update_param(wheel_params.(['invtaus_decay', num2str(wheel_params.maxorder)]), 0, 1/dt, step);
%         end

        if keystate(125) % right key => instantaneous integration
            wheel_params.(['taus_int', num2str(wheel_params.maxorder)]) = ...
                update_param_scale(wheel_params.(['taus_int', num2str(wheel_params.maxorder)]), 0.01*dt, 10, 1-step);
        elseif keystate(124) % left key => longer integration
            wheel_params.(['taus_int', num2str(wheel_params.maxorder)]) = ...
                update_param_scale(wheel_params.(['taus_int', num2str(wheel_params.maxorder)]), 0.01*dt, 10, 1+step);
        end
        
        if keystate(126) % down: decrease step size
            step = update_param_scale(step, 0.01, 10, 0.95);
        elseif keystate(127) %up: increase step size
            step = update_param_scale(step, 0.01, 10, 1.05);
        end

        if keystate(44) % < decrease step size
            ecc_r = ecc_r * 0.95;
        elseif keystate(48) %> increase step size
            ecc_r = ecc_r * 1.05;
        end

        % move cursor        
        [ux, uy, mouse0] = cursor_update(myscreen,mouse0);
        if isfield(stimulus, 'exp') && isfield(stimulus.exp, 'controlMethod') && strcmp(stimulus.exp.controlMethod, 'mouse_circ')
            tangent_vec_angle = atan2(uy,ux);
            theta = state(1); % current angle
            dtheta = (tangent_vec_angle-pi/2) - theta; % assume polar angle from tangent vector angle
            state = ou_update_state(state, dtheta , wheel_params, dt);
        else
            % use wheel
            state = ou_update_state(state, -1*ux/ ecc_r , wheel_params, dt);
        end

        % display cursor
        theta = state(1);
        mglMetalDots([ecc_r * cos(theta); ecc_r * sin(theta); 0], ...
                     [1;0;0;1], [pointer_r;pointer_r], 1, 1);
        
        if wordy
            fprintf(repmat('\b',1,numel(str2print)))
            str2print = sprintf('r = %0.3e; ux = %0.3e; theta = %0.3e; order= %i; decay= %0.3e; int= %0.3e',...
                ecc_r, ux,...
                theta, wheel_params.maxorder,...
                wheel_params.(['invtaus_decay', num2str(wheel_params.maxorder)]),...
                wheel_params.(['taus_int', num2str(wheel_params.maxorder)]));
            disp(str2print);
        end
        
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