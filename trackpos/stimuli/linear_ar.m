classdef linear_ar < trackingTask
    % The target moves in a circle, with the radius doing a random walk
    % main stimulus and background

    % stim_dyngroup
        % 0,1,2: indicate which order has
        % cv: constant velocity stim_noiseStd gets used as velocity instead
    
properties
    % task parameter
    name            = 'linear'
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
    waitsecs;       % allow movement after this many seconds
    maxtrialtime;   % maximum trial time in seconds
    iti;
    trialpause; % need to press backtick after each trial

    dyn_noise_phase;    % phase number to pull dynamics from 
    switch_tpnoise;
    
    % variable parameters
    varparams0  = {}; % parameters to move to randvars
    varparams   = {'backLum', 'ecc_r', 'pa', 'ori',...
        'stimType', 'stimLum', 'stimStd', 'stimColor', 'stim_noiseStd', 'stim_noiseTau', 'stim_vel',...
        'pointType', 'pointLum', 'pointStd', 'pointColor', 'point_noiseStd', 'point_noiseTau', 'point_vel'};
    
    backLum;
    ecc_r;
    pa; 
    ori;

    stimType;
    stimLum;
    stimStd;
    stimColor; 
    stim_noiseStd;
    stim_noiseTau;
    stim_vel;
    
    pointType;
    pointLum;
    pointStd;
    pointColor;
    point_noiseStd;
    point_noiseTau; 
    point_vel;

    experiment; % for grouping parameters. defaults to ''  if not used
    experiment_paramset; % for grouping parameters. defaults to -1 if not used

    % define segments
    segments = {'init', 'cue', 'fix', 'track', 'iti'};
end
    
    
methods
    % Constructor
    function obj = linear_ar(myscreen, varargin)
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
        
        p.addParameter('switch_tpnoise', false);
        p.addParameter('dyn_noise_phase', -1, @(x)(isnumeric(x))) ;

        p.addParameter('backLum', 0.7, @(x)(isnumeric(x)));
        p.addParameter('ecc_r', 5, @(x)(isnumeric(x) && all(x < 50))) ;
        p.addParameter('pa', 0, @(x)(isnumeric(x)));
        p.addParameter('ori', 0, @(x)(isnumeric(x) && all(x < 50))) ;

        p.addParameter('stimType', 'gaussian', @(x)(ischar(x)));
        p.addParameter('stimLum', 0.8, @(x)(isnumeric(x)));
        p.addParameter('stimStd', 1, @(x)(isnumeric(x))) ;
        p.addParameter('stimColor', 'k', @(x)(ischar(x) || iscell(x))) ;
        p.addParameter('stim_noiseStd', 1, @(x)(isnumeric(x)));
        p.addParameter('stim_noiseTau', 10/60, @(x)(isnumeric(x)))
        p.addParameter('stim_vel', 0, @(x)(isnumeric(x))) ;

        p.addParameter('pointType', 'dot', @(x)(ischar(x)));
        p.addParameter('pointLum', 1, @(x)(isnumeric(x))) ;
        p.addParameter('pointStd', 0.2, @(x)(isnumeric(x))) ;
        p.addParameter('pointColor', 'r', @(x)(ischar(x) || iscell(x))) ;
        p.addParameter('point_noiseStd', 0, @(x)(isnumeric(x)))
        p.addParameter('point_noiseTau', 10/60, @(x)(isnumeric(x)))
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

        thistask.numTrials          = obj.numTrials;
        
        thistask.getResponse        = [0, 0, 0, 0,0];
        if obj.trialpause
            thistask.synchToVol     = [0, 0, 0, 0, 1];
        else
            thistask.synchToVol     = [0, 0, 0, 0, 0];
        end
        thistask.waitForBacktick    = 0;
        

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
        
        % all the tracked variables are in linear space (not angular)
        thistask.randVars.calculated.trackTNoiseAR  = nan(1, max_seg_frames);
        thistask.randVars.calculated.trackTNoiseW   = nan(1, max_seg_frames);
        thistask.randVars.calculated.trackTNoise    = nan(1, max_seg_frames);
        
        thistask.randVars.calculated.trackPNoiseAR  = nan(1, max_seg_frames);
        thistask.randVars.calculated.trackPNoiseW   = nan(1, max_seg_frames);
        thistask.randVars.calculated.trackPNoise    = nan(1, max_seg_frames);
        
        thistask.randVars.calculated.trackControl   = nan(1, max_seg_frames);
        
        thistask.randVars.calculated.trackStim      = nan(1,max_seg_frames);
        thistask.randVars.calculated.trackResp      = nan(1,max_seg_frames);
        thistask.randVars.calculated.trackEye       = nan(max_seg_frames,2);
        thistask.randVars.calculated.trackTime      = nan(1,max_seg_frames);
        thistask.randVars.calculated.trackEyeTime   = nan(1,max_seg_frames); % for referencing edf file

    end


    % trial update?
    function [task, stimulus]  = initTrial(obj, task, myscreen, stimulus)
        dt                          = 1/myscreen.framesPerSecond;
        phaseNum                    = task.thistrial.thisphase;

        T       = (obj.maxtrialtime + 1)*myscreen.framesPerSecond;
        dt      = 1/myscreen.framesPerSecond;
        
        % load parameters from set
        if ~isempty(task.thistrial.experiment) && ~isempty(obj.varparams0)
            paramset = obj.parameter_set(task.thistrial.experiment, task.thistrial.experiment_paramset);
            for param = obj.varparams0
                eval(['task.thistrial.' param{1} '= paramset.' param{1} ';']);
            end
        end
        
        ecc_r                       = task.thistrial.ecc_r;
        pa                          = task.thistrial.pa;
        ori                         = task.thistrial.ori;
        stimulus.ecc_r              = ecc_r;
        
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

      
        if isfield(stimulus,'dyn_noise') && obj.dyn_noise_phase > 0 
            dnpn    = obj.dyn_noise_phase; % dyn_noise_phaseNum
            tn      = task.trialnum;
            % load noise, in linear space
            if obj.switch_tpnoise
                point_noiseW     = stimulus.dyn_noise(dnpn).stim_noiseW{tn};
                point_noiseAR    = stimulus.dyn_noise(dnpn).stim_noiseAR{tn};
                point_noise      = stimulus.dyn_noise(dnpn).stim_noise{tn};
            
                stim_noiseW     = stimulus.dyn_noise(dnpn).point_noiseW{tn};
                stim_noiseAR    = stimulus.dyn_noise(dnpn).point_noiseAR{tn};
                stim_noise      = stimulus.dyn_noise(dnpn).point_noise{tn};
            else
                stim_noiseW     = stimulus.dyn_noise(dnpn).stim_noiseW{tn};
                stim_noiseAR    = stimulus.dyn_noise(dnpn).stim_noiseAR{tn};
                stim_noise      = stimulus.dyn_noise(dnpn).stim_noise{tn};
            
                point_noiseW     = stimulus.dyn_noise(dnpn).point_noiseW{tn};
                point_noiseAR    = stimulus.dyn_noise(dnpn).point_noiseAR{tn};
                point_noise      = stimulus.dyn_noise(dnpn).point_noise{tn};
            end

            task.thistrial.trackTNoiseAR    = stim_noiseAR;
            task.thistrial.trackTNoiseW     = stim_noiseW;
            task.thistrial.trackTNoise      = stim_noise;
            
            task.thistrial.trackPNoiseAR    = point_noiseAR;
            task.thistrial.trackPNoiseW     = point_noiseW;
            task.thistrial.trackPNoise      = point_noise;

            % compute angular space variables
            stimulus.target.positions_trial = cumsum(stim_noise);    % angular position

            r = stimulus.target.positions_trial(1);
            stimulus.target.position = obj.r2position(ecc_r, pa, r, ori);

            stimulus.pointer.inputs_trial = point_noise;    % angular velocity inputs
            stimulus.pointer.position = stimulus.target.position;
            
        else
            % generate angular noise for target
            phi_t               = 1 - dt/task.thistrial.stim_noiseTau;
            noiseStd            = task.thistrial.stim_noiseStd * sqrt(dt) * sqrt(prod(1-phi_t.^2)); % angular noise
            [noise, wnoise]     = ar(T, noiseStd, phi_t, 'plotfigs', false);
            task.thistrial.trackTNoiseAR    = noise;
            task.thistrial.trackTNoiseW     = wnoise;

            vel = noise + task.thistrial.stim_vel;
            stimulus.target.positions_trial = cumsum(vel);

            r = stimulus.target.positions_trial(1);
            stimulus.target.position = obj.r2position(ecc_r, pa, r, ori);

            % generate angular noise for pointer
            phi_p                   = 1 - dt/task.thistrial.point_noiseTau;
            pnoiseStd               = task.thistrial.point_noiseStd * sqrt(dt) * sqrt(prod(1-phi_p.^2)); % angular noise
            [noise, wnoise]                 = ar(T, pnoiseStd, phi_p, 'plotfigs', false);
            task.thistrial.trackPNoiseAR    = noise;
            task.thistrial.trackPNoiseW     = wnoise;

            vel = noise + task.thistrial.point_vel;
            stimulus.pointer.inputs_trial = vel;
            stimulus.pointer.position = stimulus.target.position;
        end

        % wheel parameters
        if isfield(stimulus,'wheel_params') % calibrated.
            pointer_dynparams   = stimulus.wheel_params;
            pointer_dynparams.(['thetaStd', num2str(pointer_dynparams.maxorder)]) = 0;
        else
            pointer_dynparams   = ornstein_uhlenbeck(0, 0);
        end
        stimulus.pointer.dynparams  = pointer_dynparams;
        stimulus.pointer.state      = [stimulus.target.positions_trial(1); zeros(pointer_dynparams.maxorder,1)];
        
        if strcmp(task.thistrial.stimType, 'dot')
            stimulus.targetfirst = true;
        elseif strcmp(task.thistrial.pointType, 'dot')
            stimulus.targetfirst = false;
        end

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
        
        if ~stimulus.exp.showMouse, mglDisplayCursor(0);, end %hide cursor

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
                    % see how far the mouse has moved
                    [dx, dy, obj.mousestate] = cursor_update(myscreen, obj.mousestate);

                    prev_state = stimulus.pointer.state(1);
                    [stimulus.pointer.state, noise] = ou_update_state(stimulus.pointer.state, ...
                        dx, stimulus.pointer.dynparams, 1/myscreen.framesPerSecond);
                    
                    control_u_theta = stimulus.pointer.state(1) - prev_state; % angular perturbations
                    
                    stimulus.pointer.state(1) = stimulus.pointer.state(1) + stimulus.pointer.inputs_trial(framenum);
                    stimulus.pointer.position = obj.r2position(task.thistrial.ecc_r, task.thistrial.pa, stimulus.pointer.state(1), task.thistrial.ori);

                elseif strcmp(stimulus.exp.controlMethod, 'mouse_circ')
                    [ux, uy, obj.mousestate] = cursor_update(myscreen,obj.mousestate);
                    dtheta  = atan2(uy,ux);
                    [stimulus.pointer.state, noise] = ou_update_state(stimulus.pointer.state, ...
                        dtheta, stimulus.pointer.dynparams, 1/myscreen.framesPerSecond);

                    stimulus.pointer.position = obj.r2position(task.thistrial.ecc_r, task.thistrial.pa, stimulus.pointer.state(1), task.thistrial.ori);

                    control_u_theta = dtheta; % angular perturbations
                    
                elseif strcmp(stimulus.exp.controlMethod, 'mouse')
                    mInfo = mglGetMouse(myscreen.screenNumber);
                    [x,y] = screen2deg(mInfo.x, mInfo.y, myscreen);
                    stimulus.pointer.position = atan2(y,x);
                    
                    control_u = nan;
                end
            end

            if ~stimulus.exp.showMouse, mglDisplayCursor(0);, end %hide cursor

            if isfield(stimulus,'dyn_noise') % update pointer noise
                stimulus.dyn_noise(task.thistrial.thisphase).point_noise{task.trialnum}(framenum) = ...
                    (control_u_theta + stimulus.pointer.inputs_trial(framenum)) * task.thistrial.ecc_r;
            end
            
            stimulus.target.position = obj.r2position(task.thistrial.ecc_r, task.thistrial.pa, stimulus.target.positions_trial(framenum), task.thistrial.ori);

            task.thistrial.trackControl(task.thistrial.framecount) = ...
                control_u_theta; 
            task.thistrial.trackStim(task.thistrial.framecount) = ...
                stimulus.target.positions_trial(framenum); 
            task.thistrial.trackResp(task.thistrial.framecount) = ...
                stimulus.pointer.state(1);
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
            if strcmp(stimulus.pointer.color, '*')
                cuecolor = [1;0;0];
            else
                cuecolor = stimulus.pointer.color;
            end
            mglMetalArcs([stimulus.pointer.position, 0]', [cuecolor; 0.3], [2*r0; 3*r0],[0;2*pi], 1);

            da = task.thistrial.ori;
            arrow_length = 1.5; % length of arrow in deg (visual angle)
            arm_ratio = 1/3;
            arm_angle = pi/6;
            arrowidth = r0/5;

            mglMetalArrow(0.5*cos(da),0.5*sin(da),...
                da,arrow_length,arm_ratio, arm_angle, arrowidth, [1;0;0]);
            mglMetalArrow(0.5*cos(da+pi),0.5*sin(da+pi),...
                da+pi,arrow_length, arm_ratio, arm_angle, arrowidth, cuecolor);
        end
    end
    

    function pos = r2position(obj, ecc, pa, r, ori)
        center  = [ecc*cos(pa), ecc*sin(pa)];
        relpos  = [r*cos(ori), r*sin(ori)];
        pos     = center + relpos;
    end


    function params = parameter_set(obj, setname, setnum)
        params              = struct();
        params.backLum      = 0.7;
        params.ecc_r        = 15;

        stim_randcol        = {'dot', 1, 0.4, '*'}; % name, lum, size, color
        stim_lowlum         = {'gaussian', 0.4, 1, 'k'}; % name, lum, size, color
        
        [params.stimType, params.stimLum, params.stimStd, params.stimColor]      ...
            = deal(stim_randcol{:});
        params.stim_noiseStd    = 0;
        params.stim_noiseTau    = 0;
        params.stim_vel         = 0;

        [params.pointType, params.pointLum, params.pointStd, params.pointColor]  ...
            = deal(stim_lowlum{:});
        params.point_noiseStd = 0.7;
        params.point_noiseTau = 10/60;
        params.point_vel = 0;
        
        if strcmp(setname,'pa') 
            %% effect of polar angle            
            if setnum <= 3
                params.pa = - pi/2;
            elseif setnum <=6
                params.pa = 0;
            elseif setnum <= 9
                params.pa = pi/4;
            elseif setnum <= 12
                params.pa = pi/2;
            end

            if mod(setnum,3) == 1
                params.ori = 0;
            elseif mod(setnum,3) == 2
                params.ori = pi/4;
            elseif mod(setnum,3) == 0
                params.ori = pi/2;
            end
            
            if setnum == 13
                params.pa = pi/4;
                params.ori = -pi/4;
            end

        end
    end
    
end        
    
end


