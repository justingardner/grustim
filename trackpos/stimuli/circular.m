classdef circular < trackingTask
    % The target moves in a circle, with the radius doing a random walk
    % main stimulus and background
    
properties
    % task parameter
    name            = 'circular'
    numTrials       = 0;    % total number of trials
    pos_start       = cell(2,1);   % mx1 cell {(numTrials x 2), .... } of starting positions of stimuli
    
    % trial parameters
    stimulus        = cell(2,1);  % mx1 cell of images to Blt e.g. (mainstim*, background*, stim2,...)
    positions       = cell(2,1);  % mx1 cell of 4x1 position of stimulus [xpos ypos width height]. 
    
    % [r_target, theta_target, dtheta_target, r_pointer, theta_pointer]
    state;
    A;
    W;
    pidx;
    cidx = [4,5];
    
    movecursor      = 0;    
    cursor_steady   = 0; % frames for which the cursor is steady
    t0              = 0; % frames for which the cursor is steady
    
    bgfile          = '/Users/gru/data/trackpos/trackpos.mat';

    % fixed parameters
    nonvarparams   = {'steady_thresh_frame' 'steady_thresh_deg' 'waitsecs' 'maxtrialtime', 'maxtrials'};
    steady_thresh_frame;    % if subject is near the target for this long, go to next trial
    steady_thresh_deg;      % if the pointer is within this threshold of target, count frame as "steady"
    waitsecs;               % allow movement after this many seconds
    maxtrialtime;           % maximum trial time in seconds
    maxtrials;
      
    % variable parameters
    varparams   = {'backLum' 'noiseLum' 'stimLum' 'stimStd', 'thetaStep', ...
                    'thetaStd0', 'thetaStd1','thetaStd2','rStd', 'r_logSpace'};
    backLum;
    noiseLum;
    stimLum;
    stimStd;
    thetaStep;  % constant angular velocity input
    thetaStd0;  % angular velocity noise
    thetaStd1;  % angular acceleration noise 
    thetaStd2;  % angular jerk noise
    rStd;       % noise terms for r
    r_logSpace;
end
    
    
methods
    % Constructor
    function obj = circular(myscreen, varargin)
        %% parse inputs
        p = inputParser;
        p.KeepUnmatched = true;
        
        if any(cellfun(@(x)(strcmp(x,'pos_start')), varargin))
            % if start_pos exist as a parameter
            p.addParameter('pos_start', cell(1,1), @()(isnumeric(x)))
        else
            % otherwise generate test poitns around circles
            % radius can't be too big
            p.addParameter('r', [4], @(x)(isnumeric(x) && all(x < 20))) 
            p.addParameter('angle', 0:pi/2:(2*pi-pi/2), @(x)(isnumeric(x)))
        end
        
        % time threshold for going into next trial in frames
        % user needs to be steady for this amount of frames.
        p.addParameter('steady_thresh_frame', floor(1*myscreen.framesPerSecond), @(x)(isnumeric(x)))
        p.addParameter('steady_thresh_deg', 0.001, @(x)(isnumeric(x)))
        p.addParameter('waitsecs', 2, @(x)(isnumeric(x))) 
        p.addParameter('maxtrialtime', 20, @(x)(isnumeric(x)))
        p.addParameter('maxtrials', 10, @(x)(isnumeric(x)))
        
        p.addParameter('backLum', 90, @(x)(isnumeric(x))) 
        p.addParameter('noiseLum', 0, @(x)(isnumeric(x))) 
        p.addParameter('stimLum', 255, @(x)(isnumeric(x))) 
        p.addParameter('stimStd', 1, @(x)(isnumeric(x))) 
        
        p.addParameter('thetaStep', [0], @(x)(isnumeric(x))) 
        p.addParameter('thetaStd0', 0, @(x)(isnumeric(x))) 
        p.addParameter('thetaStd1', [0], @(x)(isnumeric(x)))  % (pi/30)^2
        p.addParameter('thetaStd2', [0], @(x)(isnumeric(x))) 

        p.addParameter('rStd', [0] , @(x)(isnumeric(x))) 
        p.addParameter('r_logSpace', false, @(x)(islogical(x))) 
        
        p.addParameter('randomize_order', true, @(x) (islogical(x)))
        p.parse(varargin{:})
                
        %% initialize
        % initial positions
        if isfield(p.Results, 'pos_start')
            obj.pos_start = p.Results.pos_start;
            obj.numTrials = min(p.Results.maxtrials, size(p.Results.pos_start{1},1));
        else
            % define starting position for main stimulus
            for r = p.Results.r
                for a = p.Results.angle
                    obj.pos_start{1} = [obj.pos_start{1}; r* cos(a), r * sin(a)];
                    obj.numTrials = obj.numTrials + 1;
                end
            end
            obj.numTrials = min(p.Results.maxtrials, obj.numTrials);
        end
        
        if p.Results.randomize_order
            permidx = randperm(obj.numTrials);
            obj.pos_start{1} = obj.pos_start{1}(permidx,:);
        end
        
        obj.initialize_params(p)
        
        % [r_target, theta_target, dtheta_target, r_pointer, theta_pointer]
        obj.state  = nan(5,1);       % current state vector 
        obj.A  = [1,0,0,0,0; 0,1,1,0,0; 0,0,1,0,0; 0,0,0,1,0; 0,0,0,0,1]; % dynamics update matrix
        obj.W  = [1;0;0;0;0];  % dynamics noise

    end
    
    % return task object that can be run on trackpos.m
    function thistask = configureExperiment(obj, task, myscreen, stimulus) 
        thistask        = struct();
        thistask.segmin = [obj.maxtrialtime];
        thistask.segmax = [obj.maxtrialtime];
        
        if stimulus.exp.debug
            thistask.numTrials          = 5;
        else
            thistask.numTrials          = obj.numTrials;
        end
        
        thistask.getResponse        = [0];
        thistask.synchToVol         = [0];
        thistask.waitForBacktick    = 1;
        
        for param = obj.varparams
            eval(['thistask.parameter.' param{1} ' = obj.' param{1}  ';'])
        end
    end

    % trial update?
    function task  = initTrial(obj, task, myscreen, stimulus)
        trackpos_stim = trackposInitStimulus(obj,myscreen);
        
        % initialize stimulus position
        obj.stimulus{1}    = trackpos_stim.gaussian;
        obj.positions{1}   = [obj.pos_start{1}(task.trialnum,:), [], []];
        
        % initialize background position
        if ~isempty(obj.bgfile) && task.thistrial.noiseLum > 0
            nframes = myscreen.framesPerSecond*task.segmax(1) + 20;%/downsample_timeRes; 
            task.thistrial.bgpermute(1:nframes) = randi(length(stimulus.backnoise),nframes,1);
            obj.stimulus{2}  = stimulus.backnoise{1}{task.thistrial.bgpermute(1)};
            obj.positions{2} = [0,0,myscreen.imageWidth, myscreen.imageHeight]; 
        end
        
        % set mouse position to the stimulus position. 
        x_img = obj.pos_start{1}(task.trialnum,1);  y_img = obj.pos_start{1}(task.trialnum,2);
        x_screen = x_img*myscreen.screenWidth/myscreen.imageWidth + myscreen.screenWidth/2;
        y_screen = y_img*myscreen.screenHeight/myscreen.imageHeight + myscreen.screenHeight/2;
        mglSetMousePosition(ceil(x_screen),floor(y_screen), myscreen.screenNumber); % correct for screen resolution???
        if ~stimulus.exp.showMouse, mglDisplayCursor(0);, end %hide cursor

        % trial terminal conditions
        obj.cursor_steady = 0; 
        obj.t0 = tic; % start trial
        
        % initialize state
        dtheta      = task.thistrial.thetaStep / myscreen.framesPerSecond;
        pos         = obj.cart2polar(obj.positions{1}(1:2));
                
        p0      = task.thistrial.thetaStep / myscreen.framesPerSecond;
        p_std   = task.thistrial.thetaStd0 / sqrt(myscreen.framesPerSecond);
        v_std   = task.thistrial.thetaStd1 / sqrt(myscreen.framesPerSecond);
        a_std   = task.thistrial.thetaStd2 / sqrt(myscreen.framesPerSecond);
        r_std   = task.thistrial.rStd / sqrt(myscreen.framesPerSecond);
        
        [A_theta, W_theta, ns0, ss0]  = NewtonianStateMatrix('p0', p0, 'p_std', p_std, 'v_std', v_std, 'a_std', a_std);
        [A_r, W_r, ns1, ss1]          = NewtonianStateMatrix('p_std', r_std);
        
        obj.A  = blkdiag(A_theta, A_r, eye(2,2)); 
        obj.W  = [blkdiag(W_theta, W_r); zeros(2, size(W_theta,2)+ size(W_r,2))];
        obj.W(:, all(obj.W==0,1)) = [];
        
        obj.pidx = [1, ss0+1];
        obj.cidx = [ss0+ss1+1, ss0+ss1+2];
        
        if p0==0
            obj.state   = [pos(1); zeros(ss0-1,1); ...
                   pos(2); zeros(ss1-1,1); pos'];
        else
            obj.state   = [pos(1); zeros(ss0-2,1); p0;...
                           pos(2); zeros(ss1-1,1); pos'];
        end
    end
    
    function task = startSegment(obj, task, myscreen, stimulus)
        % set mouse position to the stimulus position. 
        x_img = obj.pos_start{1}(task.trialnum,1);  y_img = obj.pos_start{1}(task.trialnum,2);
        x_screen = x_img*myscreen.screenWidth/myscreen.imageWidth + myscreen.screenWidth/2;
        y_screen = y_img*myscreen.screenHeight/myscreen.imageHeight + myscreen.screenHeight/2;
        mglSetMousePosition(ceil(x_screen),floor(y_screen), myscreen.screenNumber); % correct for screen resolution???
        if ~stimulus.exp.showMouse, mglDisplayCursor(0);, end %hide cursor        
        
    end

    % frame update
    % need to define an update function for the stimulus
    function task  = update(obj, task, myscreen, stimulus)       
        % background luminance
        mglClearScreen(task.thistrial.backLum/255);
        
        % blt all stimuli (including background)
        for stimidx = 1:length(obj.stimulus)
            mglBltTexture(obj.stimulus{stimidx}, obj.positions{stimidx})
        end
        
        % pointer updates
        if toc(obj.t0) < obj.waitsecs
            obj.movecursor = false;
            
            % set mouse position to the initial position. 
            x_img = obj.pos_start{1}(task.trialnum,1);  y_img = obj.pos_start{1}(task.trialnum,2);
            x_screen = x_img*myscreen.screenWidth/myscreen.imageWidth + myscreen.screenWidth/2;
            y_screen = y_img*myscreen.screenHeight/myscreen.imageHeight + myscreen.screenHeight/2;
            mglSetMousePosition(ceil(x_screen),floor(y_screen), myscreen.screenNumber); % correct for screen resolution???
            if ~stimulus.exp.showMouse, mglDisplayCursor(0);, end %hide cursor

        else
            obj.movecursor = true;
        end
        
        % cursorsteady, if the cursor is moving
        if obj.movecursor
            r = sqrt(sum((obj.positions{1} - stimulus.pointer).^2));
            if r < obj.steady_thresh_deg
                obj.cursor_steady = obj.cursor_steady + 1;
            else
                obj.cursor_steady =0;
            end

            if obj.cursor_steady > obj.steady_thresh_frame
                task = jumpSegment(task); 
            end
        end
        
        % update stimuli position
        obj.updateStimulus(myscreen,stimulus)
        
        % update background
        if ~isempty(obj.bgfile) && task.thistrial.noiseLum > 0
            obj.stimulus{2}  = stimulus.backnoise{1}{task.thistrial.bgpermute(task.thistrial.framecount+1)};        
        end
        
        % update fixation
        if stimulus.exp.fixateCenter == 1 % fixation below others.
            if obj.movecursor
                mglGluAnnulus(0,0,0.2,0.3,[0 0.7 0],60,1);
            else
                mglGluAnnulus(0,0,0.2,0.3,[1 1 1],60,1);
            end
            mglGluDisk(0,0,0.1,rand(1,3),60,1);
        end
    end
    
    function updateStimulus(obj, myscreen, stimulus)
        % update state
        noise       = obj.W * normrnd(0,1,size(obj.W,2),1);
        newstate    = obj.A * obj.state + noise;
        newpos      = obj.polar2cart(newstate(obj.pidx));
        
        % subtract back if out of bounds
        [horz_out, vert_out]    = check_oob(newpos, myscreen, stimulus);    
        newstate(obj.pidx(1))   = newstate(obj.pidx(1)) - horz_out * noise(obj.pidx(1));
        newstate(obj.pidx(2))   = newstate(obj.pidx(2)) - vert_out * noise(obj.pidx(2));
        
        % update position
        obj.positions{1}(1:2)   = obj.polar2cart(newstate(obj.pidx)); 
        obj.state               = newstate;
        obj.state(obj.cidx)     = obj.cart2polar(stimulus.pointer);
    end

    function pos_polar = cart2polar(obj, pos_cart)
        r = sqrt(sum(pos_cart.^2));
        if obj.r_logSpace
            r = log(r);
        end
        theta = atan2(pos_cart(2), pos_cart(1));
        pos_polar = [theta, r];
    end
    
    function pos_cart = polar2cart(obj, pos_polar)
        r = pos_polar(2);
        if obj.r_logSpace
            r = exp(r);
        end
        theta = pos_polar(1);
        pos_cart = [r*cos(theta), r*sin(theta)];
    end
    
end

end