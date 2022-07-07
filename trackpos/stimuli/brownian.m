classdef brownian < trackingTask
    
properties
    % task parameter
    name            = 'brownian'
    numTrials       = 0;    % total number of trials
    pos_start       = cell(2,1);   % mx1 cell {(numTrials x 2), .... } of starting positions of stimuli
    
    % trial parameters
    stimulus        = cell(2,1);  % mx1 cell of images to Blt e.g. (mainstim*, background, stim2,...)
                       % dimensions: 
                       % background: rgb
    positions       = cell(2,1);  % mx1 cell of 4x1 position of stimulus [xpos ypos width height]. 
    
    % [target_x, pointer_x, ...(higher order target representation)..., target_y, pointer_y]
    state;       % current state vector 
    A;           % dynamics update matrix
    W;           % dynamics noise
    cidx;        % controllable states   
    pidx;        % position index of target

    movecursor      = 0;    
    bgfile          = '/Users/gru/data/trackpos/trackpos.mat';

    % stillblob fixed parameters
    nonvarparams   = {'iti', 'maxtrialtime', 'trialpause'};
    iti;
    maxtrialtime;
    trialpause; % need to press backtick after each trial
      
    % stillblob variable parameters
    varparams   = {'backLum' 'noiseLum' 'stimLum' 'stimStd' 'stepStd' 'dynamics_order'};
    backLum;
    noiseLum;
    stimLum;
    stimStd;
    stepStd;
    dynamics_order;
end
    
    
methods
    % Constructor
    function obj = stillblob(myscreen, varargin)
        % {pos_start} or {r and angle}
        % randomize_order
        
        %% parse inputs
        p = inputParser;
        p.KeepUnmatched = true;
        
        if any(cellfun(@(x)(strcmp(x,'pos_start')), varargin))
            % if start_pos exist as a parameter
            p.addParameter('pos_start', cell(1,1), @()(isnumeric(x)))
        end
        p.addParameter('iti', 2, @(x)(isnumeric(x))) 
        p.addParameter('maxtrialtime', 15, @(x)(isnumeric(x))) 
        p.addParameter('trialpause', False, @(x)(islogical(x))) 
       
        p.addParameter('backLum', 90, @(x)(isnumeric(x))) ;
        p.addParameter('noiseLum', 0, @(x)(isnumeric(x))) ;
        p.addParameter('stimLum', 255, @(x)(isnumeric(x))) ;
        p.addParameter('stimStd', 0.4, @(x)(isnumeric(x))) ;
        p.addParameter('stepStd', 1, @(x)(isnumeric(x))) ;
        p.addParameter('dynamics_order', 0, @(x)(isinteger(x))) ;
                
        p.parse(varargin{:})
                
        %% initialize
        % initial positions
        if isfield(p.Results, 'pos_start')
            obj.pos_start = p.Results.pos_start;
            obj.numTrials = min(p.Results.maxtrials, size(p.Results.pos_start{1},1));
        else
            obj.numTrials = p.Results.maxtrials;
            x_img = 1/3*myscreen.imageWidth*(2*rand(obj.numTrials,1)-1); 
            y_img = 1/3*myscreen.imageWidth*(2*rand(obj.numTrials,1)-1);
            obj.pos_start{1} = [x_img, y_img];
        end
        
        obj.initialize_params(obj,p)
        
        obj.state = nan(4,1);       % current state vector 
        obj.A     = nan(4,4);           % dynamics update matrix
        obj.W     = nan(4,2);           % dynamics noise
        obj.cidx  = nan(4,1);        % controllable states   
    end
    
    % return task object that can be run on trackpos.m
    function thistask = configureExperiment(obj, task, myscreen, stimulus) 
        thistask        = struct();
        thistask.segmin = [obj.maxtrialtime, obj.iti];
        thistask.segmax = [obj.maxtrialtime, obj.iti];
        
        thistask.numTrials          = obj.numTrials;
        thistask.getResponse        = [0,0];
        if obj.trialpause
            thistask.synchToVol     = [0, 1];
        else
            thistask.synchToVol     = [0, 0];
        end
        thistask.waitForBacktick    = 1;
        
        for param = obj.varparams
            eval(['thistask.parameter.' param{1} ' = obj.' param{1}])
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
            obj.stimulus{2}  = stimulus.backnoise{task.thistrial.thisphase}{task.thistrial.bgpermute(1)};
            obj.positions{2} = [0,0,myscreen.imageWidth, myscreen.imageHeight]; 
        end
        
        % set mouse position to the stimulus position. 
        x_img = obj.pos_start{1}(task.trialnum,1);  y_img = obj.pos_start{1}(task.trialnum,2);
        x_screen = x_img*myscreen.screenWidth/myscreen.imageWidth + myscreen.screenWidth/2;
        y_screen = y_img*myscreen.screenHeight/myscreen.imageHeight + myscreen.screenHeight/2;
        mglSetMousePosition(ceil(x_screen),floor(y_screen), myscreen.screenNumber); % correct for screen resolution???
        if ~stimulus.exp.showmouse, mglDisplayCursor(0);, end %hide cursor
        
        % trial terminal conditions
        obj.cursor_steady = 0; 
        
        obj.t0 = tic; % start trial        
        
        % define dynamics
        ss = 2 + (task.thistrial.dynamics_order + 1)*2;
        do = task.thistrial.dynamics_order;
        
        obj.state = [obj.positions{1}(1:2)'; zeros(1,do*2); obj.positions{1}(1:2)']; % current state vector 
        obj.A     = blkdiag(triu(ones(do+1)), triu(ones(do+1)), zeros(2,2));         % dynamics update matrix
        obj.pidx  = [1, 2+do];
        obj.cidx  = [ss-2,ss-1];
        obj.W     = zeros(ss,2); % dynamics noise
        obj.W(do+1,1)   = task.thistrial.stepStd/sqrt(myscreen.framesPerSecond);
        obj.W(2*do+2,1) = task.thistrial.stepStd/sqrt(myscreen.framesPerSecond);
    end
    
    function task = startSegment(obj, task, myscreen, stimulus)
        obj.movecursor = true;
    end

    % frame update
    % need to define an update function for the stimulus
    function task  = update(obj, task, myscreen, stimulus)       
        % background luminance
        mglClearScreen(task.thistrial.backLum/255);
        
        if task.thistrial.thisseg == 1
            % blt all stimuli (including background)
            for stimidx = 1:length(obj.stimulus)
                mglBltTexture(obj.stimulus{stimidx}, obj.positions{stimidx})
            end

            % pointer updates
            if toc(obj.t0) < obj.waitsecs
                obj.movecursor = false;
            else

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
            obj.updateStimulus(stimulus)

            % update background
            if ~isempty(obj.bgfile) && task.thistrial.noiseLum > 0
                obj.stimulus{2}  = stimulus.backnoise{task.thistrial.bgpermute(task.thistrial.framecount+1)};        
            end
        end
        
        % update fixation
        if stimulus.exp.fixateCenter == 1 % fixation below others.
            if obj.movecursor
                mglGluAnnulus(0,0,0.2,0.3,[1 1 1],60,1);
            else
                mglGluAnnulus(0,0,0.2,0.3,[0 1 0],60,1);
            end
            mglGluDisk(0,0,0.1,rand(1,3),60,1);
        end
    end
    
        
    function updateStimulus(obj, stimulus)
        % update state
        obj.state               = obj.A * obj.state + obj.W * rand(size(obj.W,2),1);
        obj.state(obj.cidx)     = stimulus.position;
        
        % update position
        obj.positions{1}(1:2)   = obj.state(obj.pidx);
    end


end
    
end