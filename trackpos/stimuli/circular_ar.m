classdef circular_ar < trackingTask
    % The target moves in a circle, with the radius doing a random walk
    % main stimulus and background

    % stim_dyngroup
        % 0,1,2: indicate which order has
        % cv: constant velocity stim_noiseStd gets used as velocity instead
    
properties
    % task parameter
    name            = 'circular_ar'
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
        'stimType', 'stimLum', 'stimStd', 'stimColor', 'stim_noiseStd', 'stim_noiseAR', 'stim_vel',...
        'pointType', 'pointLum', 'pointStd', 'pointColor', 'point_noiseStd', 'point_noiseAR', 'point_vel', ...
        'ecc_r'};
    
    backLum;
    ecc_r;

    stimType;
    stimLum;
    stimStd;
    stimColor; 
    stim_noiseStd;
    stim_noiseAR;
    stim_vel;
    
    pointType;
    pointLum;
    pointStd;
    pointColor;
    point_noiseStd;
    point_noiseAR; 
    point_vel;

    experiment; % for grouping parameters. defaults to ''  if not used
    experiment_paramset; % for grouping parameters. defaults to -1 if not used

    % define segments
    segments = {'init', 'cue', 'fix', 'track', 'iti'};
end
    
    
methods
    % Constructor
    function obj = circular_ar(myscreen, varargin)
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
        p.addParameter('stim_noiseStd', 1, @(x)(isnumeric(x)));
        p.addParameter('stim_noiseAR', 0.2, @(x)(isnumeric(x)))
        p.addParameter('stim_vel', 0, @(x)(isnumeric(x))) ;

        p.addParameter('pointType', 'dot', @(x)(ischar(x)));
        p.addParameter('pointLum', 1, @(x)(isnumeric(x))) ;
        p.addParameter('pointStd', 0.2, @(x)(isnumeric(x))) ;
        p.addParameter('pointColor', 'r', @(x)(ischar(x) || iscell(x))) ;
        p.addParameter('point_noiseStd', 0, @(x)(isnumeric(x)))
        p.addParameter('point_noiseAR', 0.2, @(x)(isnumeric(x)))
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
                % do not add variable parameters to the object
                % this is now defined by experiment/experiment_paramset
                if ~any(strcmp(param{1}, obj.varparams0)) 
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

        max_seg_frames = ceil(thistask.segmax(1)*myscreen.framesPerSecond) + 20;
        thistask.randVars.calculated.randomSeed   = nan;
        thistask.randVars.calculated.initStim     = [nan nan];
        
        thistask.randVars.calculated.trackTNoiseAR  = nan(1, max_seg_frames);
        thistask.randVars.calculated.trackTNoiseW   = nan(1, max_seg_frames);
        thistask.randVars.calculated.trackPNoiseAR  = nan(1, max_seg_frames);
        thistask.randVars.calculated.trackPNoiseW   = nan(1, max_seg_frames);
        
        thistask.randVars.calculated.trackStim      = nan(1,max_seg_frames);
        thistask.randVars.calculated.trackResp      = nan(1,max_seg_frames);
        thistask.randVars.calculated.trackEye     = nan(max_seg_frames,2);
        thistask.randVars.calculated.trackTime    = nan(1,max_seg_frames);
        thistask.randVars.calculated.trackEyeTime = nan(1,max_seg_frames); % for referencing edf file

    end


    % trial update?
    function [task, stimulus]  = initTrial(obj, task, myscreen, stimulus)
        dt                          = 1/myscreen.framesPerSecond;
        
        task.thistrial.ecc_r        = task.thistrial.ecc_r;
        stimulus.ecc_r              = task.thistrial.ecc_r;

        T = (obj.maxtrialtime + 1)*myscreen.framesPerSecond;
        dt = 1/myscreen.framesPerSecond;

        % load parameters from set
        if ~isempty(task.thistrial.experiment) && ~isempty(obj.varparams0)
            paramset = obj.parameter_set(task.thistrial.experiment, task.thistrial.experiment_paramset);
            for param = obj.varparams0
                eval(['task.thistrial.' param{1} '= paramset.' param{1} ';']);
            end
        end

        % set up target
        target = stimulus.target;
        reinit = 0;
        for param = {'stimType', 'stimStd', 'stimLum', 'stimColor'}
            if ~isfield(target, lower(param{1}(5:end))) || (target.(lower(param{1}(5:end))) ~= task.thistrial.(param{1}))
                eval(['target.' lower(param{1}(5:end)) ' = task.thistrial.' param{1} ';'])
                reinit = 1;
            end
        end
        stimulus.target     = trackposInitStimulus(target,myscreen, 'reinit_img', reinit);
        
        stimulu.target.ar   = task.thistrial.stim_noiseAR;
        noiseStd            = task.thistrial.stim_noiseStd * sqrt(dt) / task.thistrial.ecc_r * sqrt(prod(1-stimulu.target.ar.^2)); % angular noise
        [noise, wnoise]     = ar(T, noiseStd, stimulu.target.ar, 'plotfigs', false);
        task.thistrial.trackTNoiseAR    = task.thistrial.ecc_r * noise;
        task.thistrial.trackTNoiseW     = task.thistrial.ecc_r * wnoise;

        vel = noise + task.thistrial.stim_vel;
        stimulus.target.positions_trial = cumsum(vel) + 2*pi*rand();
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

        stimulu.pointer.ar   = task.thistrial.point_noiseAR;
        pnoiseStd            = task.thistrial.point_noiseStd * sqrt(dt) / task.thistrial.ecc_r * sqrt(prod(1-stimulu.pointer.ar .^2)); % angular noise

        [noise, wnoise]      = ar(T, pnoiseStd, stimulu.pointer.ar, 'plotfigs', false);
        task.thistrial.trackPNoiseAR    = task.thistrial.ecc_r * noise;
        task.thistrial.trackPNoiseW     = task.thistrial.ecc_r * wnoise;

        vel = noise + task.thistrial.point_vel;
        stimulus.pointer.inputs_trial = vel;
        stimulus.pointer.position = stimulus.target.position;

        if isfield(stimulus,'wheel_params') % calibrated.
            pointer_dynparams   = stimulus.wheel_params;
            pointer_dynparams.(['thetaStd', num2str(pointer_dynparams.maxorder)]) = 0;
        else
            pointer_dynparams   = ornstein_uhlenbeck(0, 0);
        end
        stimulus.pointer.dynparams  = pointer_dynparams;
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

                    [stimulus.pointer.state, noise] = ou_update_state(stimulus.pointer.state, ...
                        -1*dx/task.thistrial.ecc_r, stimulus.pointer.dynparams, 1/myscreen.framesPerSecond);
                    
                    stimulus.pointer.state(1) = stimulus.pointer.state(1) + stimulus.pointer.inputs_trial(framenum);
                    stimulus.pointer.position = obj.polar2cart(task.thistrial.ecc_r, stimulus.pointer.state(1));

                elseif strcmp(stimulus.exp.controlMethod, 'mouse_circ')
                    [ux, uy, obj.mousestate] = cursor_update(myscreen,obj.mousestate);
                    dtheta  = atan2(uy,ux);
                    [stimulus.pointer.state, noise] = ou_update_state(stimulus.pointer.state, ...
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
            
            % eye tracking
            if stimulus.exp.trackEye
                % mouse version for testing with no eyetracker
                [pos,postime] = mglEyelinkGetCurrentEyePos; % is this in image coordinates?
                task.thistrial.trackEye(task.thistrial.framecount,:)    = pos;
                task.thistrial.trackEyeTime(task.thistrial.framecount)  = postime;
            end

        end

        if strcmp(obj.segments{task.thistrial.thisseg}, 'cue')
            r0 = stimulus.pointer.std;
            cuecolor = stimulus.pointer.color;
            mglMetalArcs([stimulus.pointer.position, 0]', [cuecolor; 0.3], [2*r0; 3.5*r0],[0;2*pi], 1);
        end
    end

    
    function pos_cart = polar2cart(obj, r, theta)
        pos_cart = [r*cos(theta), r*sin(theta)];
    end


    function params = parameter_set(obj, setname, setnum)
        params = struct();
        params.backLum = 0.7;
        params.ecc_r = 10;

        stim_highlum    = {'dot', 1, 0.4,'r'}; % lum, size, color
        stim_lowlum     = {'gaussian', 0.4, 1,'k'}; % lum, size, color
        
        [params.stimType, params.stimLum, params.stimStd, params.stimColor]      ...
            = deal(stim_lowlum{:});

        params.stim_noiseStd    = 1;
        params.stim_noiseAR     = 0;
        params.stim_vel         = 0;

        [params.pointType, params.pointLum, params.pointStd, params.pointColor]  ...
            = deal(stim_highlum{:});
        params.point_noiseStd = 1;
        params.point_noiseAR = 0;
        params.point_vel = 0;

        if strcmp(setname,'mn') % effect of motor noise 
            if mod(setnum,3) == 1
                ar1 = 0.9;
                ar2 = 0;
            elseif mod(setnum,3) == 2
                ar1 = 0.5;
                ar2 = 0;
            elseif mod(setnum,6) == 3
                ar1 = 0;
                ar2 = 0;
            end

            if setnum <= 6
                noisestd1 = 1;
                noisestd2 = 0;
            elseif setnum <=12
                noisestd1 = 1;
                noisestd2 = 1;
            end

            if setnum <= 3 || (setnum > 6 &  setnum <=9)
                % stimulus: high luminance
                % pointer: low luminance,
                [params.stimType, params.stimLum, params.stimStd, params.stimColor]         = deal(stim_highlum{:});
                [params.pointType, params.pointLum, params.pointStd, params.pointColor]     = deal(stim_lowlum{:});

                params.stim_noiseStd    = noisestd2;
                params.point_noiseStd   = noisestd1;

                params.stim_noiseAR     = ar2;
                params.point_noiseAR    = ar1;
            elseif setnum <= 6 || (setnum > 9 &  setnum <=12)
                % stimulus: low luminance, moving
                % pointer: high luminance
                [params.stimType, params.stimLum, params.stimStd, params.stimColor]         = deal(stim_lowlum{:});
                [params.pointType, params.pointLum, params.pointStd, params.pointColor]     = deal(stim_highlum{:});
                
                params.stim_noiseStd    = noisestd1;
                params.point_noiseStd   = noisestd2;

                params.stim_noiseAR     = ar1;
                params.point_noiseAR    = ar2;
            end
        end
    end
    
end        
    
end
