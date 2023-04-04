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
    varparams0  = {}; % parameters to move to randvars
    varparams   = {'backLum', ...
        'stimType', 'stimLum', 'stimStd', 'stimColor', 'stim_dyngroup', 'stim_noiseStd', 'stim_vel',...
        'pointType', 'pointLum', 'pointStd', 'pointColor', 'point_dyngroup', 'point_noiseStd', 'point_vel', ...
        'ecc_r'};
    
    backLum;
    ecc_r;

    stimType;
    stimLum;
    stimStd;
    stimColor; 
    stim_dyngroup; 
    stim_noiseStd;
    stim_vel;
    
    pointType;
    pointLum;
    pointStd;
    pointColor;
    point_dyngroup; 
    point_noiseStd;
    point_vel;

    experiment; % for grouping parameters. defaults to ''  if not used
    experiment_paramset; % for grouping parameters. defaults to -1 if not used

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

        p.addParameter('stimType', 'gaussian', @(x)(ischar(x)));
        p.addParameter('stimLum', 0.8, @(x)(isnumeric(x)));
        p.addParameter('stimStd', 1, @(x)(isnumeric(x))) ;
        p.addParameter('stimColor', 'k', @(x)(ischar(x) || iscell(x))) ;
        p.addParameter('stim_dyngroup', 0, @(x)(isnumeric(x) || iscell(x))) ;
        p.addParameter('stim_noiseStd', 1, @(x)(isnumeric(x))) ;
        p.addParameter('stim_vel', 0, @(x)(isnumeric(x))) ;

        p.addParameter('pointType', 'dot', @(x)(ischar(x)));
        p.addParameter('pointLum', 1, @(x)(isnumeric(x))) ;
        p.addParameter('pointStd', 0.2, @(x)(isnumeric(x))) ;
        p.addParameter('pointColor', 'r', @(x)(ischar(x) || iscell(x))) ;
        p.addParameter('point_dyngroup', 0, @(x)(ischar(x) || iscell(x))) ;
        p.addParameter('point_noiseStd', 0, @(x)(isnumeric(x)))
        p.addParameter('point_vel', 0, @(x)(isnumeric(x)))

        p.addParameter('experiment', '', @(x)(ischar(x) || iscell(x)))
        p.addParameter('experiment_paramset', -1, @(x)(isnumeric(x)))
        
        p.parse(varargin{:})
        obj.initialize_params(p)

        if isempty(p.Results.experiment) || any(p.Results.experiment_paramset < 0)
            for param = p.Parameters
                eval(['obj.' param{1} ' = p.Results.' param{1} ';'])
            end
        else
            obj.varparams0 = obj.varparams;
            obj.varparams = {'experiment', 'experiment_paramset'};
            for param = p.Parameters
                if ~any(strcmp(param{1}, obj.varparams0)) % do not add variable parameters to the object -- this is now defined by experiment/experiment_paramset
                    eval(['obj.' param{1} ' = p.Results.' param{1} ';'])
                end
            end
        end
    end

    % return task object that can be run on trackpos.m
    function thistask = configureExperiment(obj, task, myscreen, stimulus) 
        thistask        = struct();
        
        inittime        = 1;
        cue             = 1;
        fix             = 0.5;
        iti             = max(0, obj.iti-inittime-fix-cue);
        
        thistask.segmin = [inittime, cue, fix, obj.maxtrialtime, iti];
        thistask.segmax = [inittime, cue, fix, obj.maxtrialtime, iti];
        
        if stimulus.exp.debug
            thistask.numTrials          = 5;
        else
            thistask.numTrials          = obj.numTrials;
        end
        
        thistask.getResponse        = [0, 0, 0, 0,0];
        if obj.trialpause
            thistask.synchToVol     = [0, 0, 0, 0, 1];
        else
            thistask.synchToVol     = [0, 0, 0, 0, 0];
        end
        thistask.waitForBacktick    = 1;
        

        thistask.random             = 1; % randomized parameters by default
        for param = obj.varparams
            eval(['thistask.parameter.' param{1} ' = obj.' param{1}  ';'])
        end

        for param = obj.varparams0
            if contains(lower(param{1}),'type')
                eval(['thistask.randVars.calculated.' param{1} '= ''somestring'';'])
            else
                eval(['thistask.randVars.calculated.' param{1} '= nan;'])
            end
        end
    end


    % trial update?
    function [task, stimulus]  = initTrial(obj, task, myscreen, stimulus)
        dt              = 1/myscreen.framesPerSecond;
        
        task.thistrial.ecc_r       = task.thistrial.ecc_r;
        stimulus.ecc_r  = task.thistrial.ecc_r;

        if ~isempty(task.thistrial.experiment) && ~isempty(obj.varparams0)
            paramset = obj.parameter_set(task.thistrial.experiment, task.thistrial.experiment_paramset);
            for param = obj.varparams0
                eval(['task.thistrial.' param{1} '= paramset.' param{1} ';']);
            end
        end

        target = stimulus.target;
        reinit = 0;
        for param = {'stimType', 'stimStd', 'stimLum', 'stimColor'}
            if ~isfield(target, lower(param{1}(5:end))) || (target.(lower(param{1}(5:end))) ~= task.thistrial.(param{1}))
                eval(['target.' lower(param{1}(5:end)) ' = task.thistrial.' param{1} ';'])
                reinit = 1;
            end
        end
        stimulus.target     = trackposInitStimulus(target,myscreen, 'reinit_img', reinit);
        noiseStd            = task.thistrial.stim_noiseStd / sqrt(dt) / task.thistrial.ecc_r; % angular noise
        stim_dynparams      = obj.dynparam_group(task.thistrial.stim_dyngroup, noiseStd);

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

        stimulus.target.position = obj.polar2cart(task.thistrial.ecc_r, stimulus.target.positions_trial(1));
        
        % initialize the pointer
        pointer             = stimulus.pointer;
        reinit = 0;
        for param = {'pointType', 'pointStd', 'pointLum', 'pointColor'}
            if ~isfield(pointer, lower(param{1}(6:end))) || (pointer.(lower(param{1}(6:end))) ~= task.thistrial.(param{1}))
                eval(['pointer.' lower(param{1}(6:end)) ' = task.thistrial.' param{1} ';'])
                reinit = 1;
            end
        end
        stimulus.pointer    = trackposInitStimulus(pointer, myscreen, 'reinit_img', reinit); 

        if isfield(task.thistrial, 'point_noiseStd')
            pnoiseStd = task.thistrial.point_noiseStd / sqrt(dt)  / task.thistrial.ecc_r ;
        else
            pnoiseStd = 0;
        end

        if isfield(stimulus,'wheel_params') % calibrated.
            pointer_dynparams   = stimulus.wheel_params;
            tau_noise = pointer_dynparams.(['taus_int', num2str(pointer_dynparams.maxorder)]);
            pointer_dynparams.(['thetaStd', num2str(pointer_dynparams.maxorder)]) = pnoiseStd * tau_noise;
        elseif isfield(task.thistrial,'point_dyngroup')
            pointer_dynparams   = obj.dynparam_group(task.thistrial.point_dyngroup, pnoiseStd);
        else
            pointer_dynparams   = obj.dynparam_group(0, pnoiseStd);
        end

        stimulus.pointer.dynparams  = pointer_dynparams;
        stimulus.pointer.position   = obj.polar2cart(task.thistrial.ecc_r, stimulus.target.positions_trial(1));
        stimulus.pointer.state      = [stimulus.target.positions_trial(1); zeros(pointer_dynparams.maxorder,1)];

        obj.t0 = tic; % start trial
    end

    
    function stimulus = startSegment(obj, task, myscreen, stimulus)       
        if strcmp(obj.segments{task.thistrial.thisseg}, 'track')
            % task.thistrial.thisseg == 2 % tracking
            % set mouse position to the center of the screen
            [x_screen,y_screen] = deg2screen(0, 0, myscreen);
            mglSetMousePosition(ceil(x_screen),floor(y_screen), myscreen.screenNumber);
            obj.mousestate = [0,0];
            
            if ~stimulus.exp.showMouse, mglDisplayCursor(0);, end %hide cursor        

            obj.movepointer  = true;
            obj.doTrack     = true;
            obj.displayFix 	= true;
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
                    
                    stimulus.pointer.position = obj.polar2cart(task.thistrial.ecc_r, stimulus.pointer.state(1));
                elseif strcmp(stimulus.exp.controlMethod, 'mouse_circ')
                    [ux, uy, obj.mousestate] = cursor_update(myscreen,obj.mousestate);
                    dtheta  = atan2(uy,ux);
                    stimulus.pointer.state = ou_update_state(stimulus.pointer.state, ...
                        dtheta, stimulus.pointer.dynparams, 1/myscreen.framesPerSecond);

                    stimulus.pointer.position = obj.polar2cart(task.thistrial.ecc_r, stimulus.pointer.state(1));
                elseif strcmp(stimulus.exp.controlMethod, 'mouse')
                    mInfo = mglGetMouse(myscreen.screenNumber);
                    [x,y] = screen2deg(mInfo.x, mInfo.y, myscreen);
                    stimulus.pointer.position = atan2(y,x);
                end
            end

            if ~stimulus.exp.showMouse, mglDisplayCursor(0);, end %hide cursor

            stimulus.target.position = obj.polar2cart(task.thistrial.ecc_r, stimulus.target.positions_trial(framenum));

            task.thistrial.trackStim(task.thistrial.framecount) = ...
                task.thistrial.ecc_r * stimulus.target.positions_trial(framenum); 
            task.thistrial.trackResp(task.thistrial.framecount) = ...
                task.thistrial.ecc_r * stimulus.pointer.state(1);
            task.thistrial.trackTime(task.thistrial.framecount) = mglGetSecs(stimulus.t0); 
        end

        if strcmp(obj.segments{task.thistrial.thisseg}, 'cue')
            r0 = stimulus.pointer.std;
            cuecolor = stimulus.pointer.color;
            mglMetalArcs([stimulus.pointer.position, 0]', [cuecolor; 0.3], [2*r0; 3.5*r0],[0;2*pi], 1);
        end
    end

    function param = dynparam_group(obj, groupname, noisestd)
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


    function params = parameter_set(obj,setname, setnum)
        params = struct();
        params.backLum = 0.7;
        params.ecc_r = 10;

        params.stimType = {'gaussian'};
        params.stimLum = 0.4;
        params.stimStd = 1;
        params.stimColor = 'k'; 
        params.stim_dyngroup = 0; 
        params.stim_noiseStd = 1;
        params.stim_vel = 0;
    
        params.pointType = {'dot'};
        params.pointLum = 1;   
        params.pointStd = 0.4;
        params.pointColor ='r';
        params.point_dyngroup = 0; 
        params.point_noiseStd = 0;
        params.point_vel = 0;

        if strcmp(setname, 'ecc') % effect of eccentricity
            ecc_r_list               = [5, 10, 15, 20]; % eccentricity

            if setnum < 10 % brownian motion
                params.stim_noiseStd       = 1; % in dva per second (linear velocity) 
                params.stim_dyngroup       = 0; % noise order, same size as stimStdList % 10: constant velocity
                params.stim_vel            = 0; 
            elseif setnum >= 10 % constant velocity, no brownain
                params.stim_noiseStd       = 0; % in dva per second (linear velocity) 
                params.stim_dyngroup       = 10; % noise order, same size as stimStdList % 10: constant velocity
                params.stim_vel            = 1; 
            end

            setorder = mod(setnum,length(ecc_r_list));
            if setorder == 0
                setorder = length(ecc_r_list);
            end
            params.ecc_r = ecc_r_list(setorder);

        elseif strcmp(setname,'mn') % effect of motor noise            
            stim_highlum    = {'dot', 1, 0.4,'r'}; % lum, size, color
            stim_lowlum     = {'gaussian', 0.4, 1,'k'}; % lum, size, color

            if setnum <= 12
                if setnum <= 6
                    % stimulus: low luminance
                    % pointer: high luminance
                    [params.stimType, params.stimLum, params.stimStd, params.stimColor]      = deal(stim_lowlum{:});
                    [params.pointType, params.pointLum, params.pointStd, params.pointColor]   = deal(stim_highlum{:});
                else
                    % stimulus: low luminance
                    % pointer: high luminance
                    [params.stimType, params.stimLum, params.stimStd, params.stimColor]      = deal(stim_highlum{:});
                    [params.pointType, params.pointLum, params.pointStd, params.pointColor]   = deal(stim_lowlum{:});
                end
                
                if mod(setnum,3) == 1
                    params.point_noiseStd = 0;
                elseif mod(setnum,3) == 2
                    params.point_noiseStd = 1;
                elseif mod(setnum,3) == 0
                    params.point_noiseStd = 2;
                end

                if mod(setnum,2) == 1
                    params.stim_noiseStd = 1;
                elseif mod(setnum,2) == 0
                    params.stim_noiseStd = 2;
                end
            elseif (setnum > 12) && (setnum <= 18) 
                % stabilization task - no stimulus dynamics noise
                params.stim_noiseStd = 0;
                
                if setnum <= 15
                    % stimulus: low luminance
                    % pointer: high luminance
                    [params.stimType, params.stimLum, params.stimStd, params.stimColor]      = deal(stim_lowlum{:});
                    [params.pointType, params.pointLum, params.pointStd, params.pointColor]   = deal(stim_highlum{:});
                else
                    % stimulus: low luminance
                    % pointer: high luminance
                    [params.stimType, params.stimLum, params.stimStd, params.stimColor]      = deal(stim_highlum{:});
                    [params.pointType, params.pointLum, params.pointStd, params.pointColor]   = deal(stim_lowlum{:});
                end
                
                % change only pointer noise -- 3 pointer noise conditions
                if mod(setnum,3) == 1
                    params.point_noiseStd = 1;
                elseif mod(setnum,3) == 2
                    params.point_noiseStd = 2;
                elseif mod(setnum,3) == 0
                    params.point_noiseStd = 10;
                end
            end
        elseif strcmp(setname, 'cv')
             % moving stimulus, change stimulus velocity
             params.point_noiseStd = 0;
             
             vel_list               = [1,2,3];
             ecc_r_list             = [5, 10, 15, 20]; % eccentricity
             stimLum_list           = [0.2, 0.5, 0.8];
             if setnum <=9
                 start = 0;
                 params.ecc_r = ecc_r_list(1);
             elseif setnum <= 18
                 start = 10;
                 params.ecc_r = ecc_r_list(2);
             elseif setnum <= 27
                 start = 19;
                 params.ecc_r = ecc_r_list(3);
             elseif setnum <= 36
                 start = 28;
                 params.ecc_r = ecc_r_list(4);
             end
                 
             % change velocity 
            if (setnum > start+0) && (setunm <= start+3)
                params.stim_vel = 1;
            elseif (setnum > start+3) && (setunm <= start+6)
                params.stim_vel = 2;
            elseif (setnum > start+6) && (setunm <= start+9)
                params.stim_vel = 3;
            end

            % change luminance
            if mod(setnum,3) == 1
                params.stimLum = stimLum_list(1);
            elseif mod(setnum,3) == 2
                params.stimLum = stimLum_list(2);
            elseif mod(setnum,3) == 3
                params.stimLum = stimLum_list(3);
            end
        end
    end
    
end        
    
end
