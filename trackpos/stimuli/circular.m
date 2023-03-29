classdef circular < trackingTask
    % The target moves in a circle, with the radius doing a random walk
    % main stimulus and background

    % stim_dyngroup
        % 0,1,2: indicate which order has
        % cv: constant velocity stim_noiseStd gets used as velocity instead
    
properties
    % task parameter
    name            = 'circular'
    numTrials       = 0;    % total number of trials
    
    % cells to be initalized
    stimulus        = cell(1,1);  % mx1 cell of images to blt e.g. (mainstim*, background*, stim2,...)
    theta           = cell(1,1);   % mx1 cell {(numTrials x T), .... } of starting positions of stimuli

    % trial parameters
    t0              = 0; % start time

    movepointer     = false;
    doTrack         = true;    % indicates whether we should start recording tracking variables
    displayFix      = true; % display fixation at current segmention
    
    mousestate      = [0,0];

    % fixed parameters
    nonvarparams    = {'waitsecs', 'maxtrialtime', 'iti', 'trialpause'};
    waitsecs;               % allow movement after this many seconds
    maxtrialtime;           % maximum trial time in seconds
    iti;
    trialpause; % need to press backtick after each trial

    % variable parameters
    varparams   = {'backLum', ...
        'stimLum', 'stimStd', 'stimColor', 'stim_dyngroup', 'stim_noiseStd', 'stim_vel',...
        'pointLum', 'pointStd', 'pointColor', 'point_dyngroup', 'point_noiseStd', ...
        'ecc_r'};
    backLum;
    ecc_r;

    stimLum;
    stimStd;
    stimColor; 
    stim_dyngroup; 
    stim_noiseStd;
    stim_vel;
    
    pointLum;
    pointStd;
    pointColor;
    point_dyngroup; 
    point_noiseStd;

    % define segments
    segments = {'init', 'cue', 'fix', 'track', 'iti'};
end
    
    
methods
    % Constructor
    function obj = circular(myscreen, varargin)
        %% parse inputs
        % set to vector if you want to randomize parameter

        % initialize other parameters
        p = inputParser;
        p.KeepUnmatched = true;      
        p.addParameter('numTrials', 20, @(x)(isnumeric(x)));
        p.addParameter('waitsecs', 0, @(x)(isnumeric(x)));
        p.addParameter('maxtrialtime', 25, @(x)(isnumeric(x)));

        p.addParameter('iti', 3, @(x)(isnumeric(x)));
        p.addParameter('trialpause', false);

        p.addParameter('backLum', 0.7, @(x)(isnumeric(x)));
        p.addParameter('ecc_r', 5, @(x)(isnumeric(x) && all(x < 50))) ;

        p.addParameter('stimLum', 0.8, @(x)(isnumeric(x)));
        p.addParameter('stimStd', 1, @(x)(isnumeric(x))) ;
        p.addParameter('stimColor', 'k', @(x)(ischar(x) || iscell(x))) ;
        p.addParameter('stim_dyngroup', 0, @(x)(isnumeric(x) || iscell(x))) ;
        p.addParameter('stim_noiseStd', 1, @(x)(isnumeric(x))) ;
        p.addParameter('stim_vel', 0, @(x)(isnumeric(x))) ;

        p.addParameter('pointLum', 1, @(x)(isnumeric(x))) ;
        p.addParameter('pointStd', 0.2, @(x)(isnumeric(x))) ;
        p.addParameter('pointColor', 'r', @(x)(ischar(x) || iscell(x))) ;
        p.addParameter('point_dyngroup', '0', @(x)(ischar(x) || iscell(x))) ;
        p.addParameter('point_noiseStd', 0, @(x)(isnumeric(x)))
        
        p.parse(varargin{:})
        obj.initialize_params(p)
        
    end

    % return task object that can be run on trackpos.m
    function thistask = configureExperiment(obj, task, myscreen, stimulus) 
        thistask        = struct();
        inittime = 1;
        cue = 1;
        fix = 0.5;
        iti = max(0, obj.iti-inittime-fix-cue);
        
        thistask.segmin = [inittime, cue, fix, obj.maxtrialtime, iti];
        thistask.segmax = [inittime, cue, fix, obj.maxtrialtime, iti];
        
        if stimulus.exp.debug
            thistask.numTrials          = 5;
        else
            thistask.numTrials          = obj.numTrials;
        end
        
        thistask.getResponse        = [0, 0, 0];
        if obj.trialpause
            thistask.synchToVol     = [0, 0, 1];
        else
            thistask.synchToVol     = [0, 0, 0];
        end
        thistask.waitForBacktick    = 1;
        

        thistask.random             = 1; % randomized parameters by default
        
        for param = obj.varparams
            eval(['thistask.parameter.' param{1} ' = obj.' param{1}  ';'])
        end
    end

    % trial update?
    function [task, stimulus]  = initTrial(obj, task, myscreen, stimulus)
        dt              = 1/myscreen.framesPerSecond;
        
        
        obj.ecc_r       = task.thistrial.ecc_r;
        stimulus.ecc_r  = task.thistrial.ecc_r;

        target = struct();
        for param = {'stimStd', 'stimLum', 'stimColor'}
           eval(['target.' lower(param{1}(5:end)) ' = task.thistrial.' param{1} ';'])
        end
        stimulus.target     = trackposInitStimulus(target,myscreen);
        noiseStd            = task.thistrial.stim_noiseStd / sqrt(dt) / obj.ecc_r; % angular noise
        stim_dynparams      = obj.parameter_group(task.thistrial.stim_dyngroup, noiseStd);

        stimulus.target.dynparams = stim_dynparams;
        T = (obj.maxtrialtime + 1)*myscreen.framesPerSecond;
        dt = 1/myscreen.framesPerSecond;
        if task.thistrial.stim_dyngroup == 10
            state           = zeros(T, 1);
            state(1,1)      = task.thistrial.stim_vel / task.thistrial.ecc_r / dt; % linear velocity to angular velocity input (deg/s^2)
            stimulus.target.positions_trial = ou_simulate_full(stim_dynparams, T, dt,'state', state);
        else
            stimulus.target.positions_trial = ou_simulate_full(stim_dynparams, T, dt);
        end
        stimulus.target.position = obj.polar2cart(obj.ecc_r, stimulus.target.positions_trial(1));
        
        % initialize the pointer
        pointer             = struct();
        for param = {'pointStd', 'pointLum', 'pointColor'}
           eval(['pointer.' lower(param{1}(6:end)) ' = task.thistrial.' param{1} ';'])
        end
        stimulus.pointer    = pointer; %trackposInitStimulus(pointer,myscreen); 

        if isfield(task.thistrial, 'point_noiseStd')
            pnoiseStd = task.thistrial.point_noiseStd / sqrt(dt)  / obj.ecc_r ;
        else
            pnoiseStd = 0;
        end

        if isfield(stimulus,'wheel_params') % calibrated.
            pointer_dynparams   = stimulus.wheel_params;
            tau_noise = pointer_dynparams.(['taus_int', num2str(pointer_dynparams.maxorder)]);
            pointer_dynparams.(['thetaStd', num2str(pointer_dynparams.maxorder)]) = pnoiseStd * tau_noise;
        elseif isfield(task.thistrial,'point_dyngroup')
            pointer_dynparams   = obj.parameter_group(task.thistrial.point_dyngroup, pnoiseStd);
        else
            pointer_dynparams   = obj.parameter_group('0', pnoiseStd);
        end

        stimulus.pointer.dynparams  = pointer_dynparams;
        stimulus.pointer.position   = obj.polar2cart(obj.ecc_r, stimulus.target.positions_trial(1));
        stimulus.pointer.state      = [stimulus.target.positions_trial(1); zeros(pointer_dynparams.maxorder,1)];

        obj.t0 = tic; % start trial
    end

    
    function stimulus = startSegment(obj, task, myscreen, stimulus)       
        if strcmp(obj.segments{task.thistrial.thisseg}, 'tracking')
            disp('starting tracking segment')
            % task.thistrial.thisseg == 2 % tracking
            % set mouse position to the center of the screen
            [x_screen,y_screen] = deg2screen(0, 0, myscreen);
            mglSetMousePosition(ceil(x_screen),floor(y_screen), myscreen.screenNumber);
            obj.mousestate = [0,0];
            
            if ~stimulus.exp.showMouse, mglDisplayCursor(0);, end %hide cursor        

            obj.movepointer  = true;
            obj.doTrack     = true;
            obj.displayFix 	= true;
        elseif task.thistrial.thisseg == task.numsegs % ITI or strcmp(obj.segments(task.thistrial.thisseg), 'tracking')
            obj.movepointer  = false;
            obj.doTrack     = false;
            obj.displayFix 	= true;
            
            stimulus.target = [];
            stimulus.pointer = [];
        else 
            obj.movepointer  = false;
            obj.doTrack     = false;
            obj.displayFix 	= true;
        end   

    end

    % frame update
    % need to define an update function for the stimulus
    function [task, stimulus]  = update(obj, task, myscreen, stimulus) 
        framenum = task.thistrial.framecount;

        if obj.doTrack
            if toc(obj.t0) < obj.waitsecs
                obj.movepointer = false;

                % set mouse position to the middle. 
                [x_screen,y_screen] = deg2screen(0, 0, myscreen);
                mglSetMousePosition(ceil(x_screen),floor(y_screen), myscreen.screenNumber);
            else
                obj.movepointer = true;
            end

            if obj.movepointer
                if strcmp(stimulus.exp.controlMethod, 'wheel')
                    % see how far the mouse as moved
                    [dx, dy, obj.mousestate] = cursor_update(myscreen, obj.mousestate);

                    stimulus.pointer.state = ou_update_state(stimulus.pointer.state, ...
                        -1*dx/task.thistrial.ecc_r, stimulus.pointer.dynparams, 1/myscreen.framesPerSecond);
                    stimulus.pointer.position = obj.polar2cart(obj.ecc_r, stimulus.pointer.state(1));
                elseif strcmp(stimulus.exp.controlMethod, 'mouse')
                    mInfo = mglGetMouse(myscreen.screenNumber);
                    [x,y] = screen2deg(mInfo.x, mInfo.y, myscreen);
                    stimulus.pointer.position = atan2(y,x);
                end
            end

            if ~stimulus.exp.showMouse, mglDisplayCursor(0);, end %hide cursor

            stimulus.target.position = obj.polar2cart(obj.ecc_r, stimulus.target.positions_trial(framenum));

            task.thistrial.trackStim(task.thistrial.framecount) = ...
                task.thistrial.ecc_r * stimulus.target.positions_trial(framenum); 
            task.thistrial.trackResp(task.thistrial.framecount) = ...
                task.thistrial.ecc_r * stimulus.pointer.state(1);
            task.thistrial.trackTime(task.thistrial.framecount) = mglGetSecs(stimulus.t0); 
        end

        if strcmp(obj.segments{task.thistrial.thisseg}, 'cue')
            r0 = stimulus.pointerR;
            mglMetalArcs([stimulus.pointer.position, 0]', [1;0;0; 0.3], [2*r0; 3.5*r0],[0;2*pi], 1);
        end
    end

    function param = parameter_group(obj, groupname, noisestd)
        if groupname == 0 % white velocity
            param = ornstein_uhlenbeck(0, noisestd);
        elseif groupname == 1 % white acceleration
            param = ornstein_uhlenbeck(1, noisestd);
        elseif groupname == 2 % white jerk
            param = ornstein_uhlenbeck(2, noisestd);
        elseif groupname == 10 % constant velocity
            param = struct();
            param.maxorder              = 1;
            param.('thetaStd0')         = noisestd;
            param.('invtaus_decay0')    = 0;
            param.('taus_int0')         = 1;
            
            param.('thetaStd1')         = 0;
            param.('invtaus_decay1')    = 0;
            param.('taus_int1')         = 1;
        end
    end

    function pos_cart = polar2cart(obj, r, theta)
        pos_cart = [r*cos(theta), r*sin(theta)];
    end
    
end

end