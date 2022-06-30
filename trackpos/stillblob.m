classdef stillblob < trackingTask
    
properties
    % task parameter
    name        = 'stillblob'
    numTrials  = 0;    % total number of trials
    pos_start   = {};   % {(n x 2), .... } starting positions of stimulus
    
    % trial paremters
    stimulus    = {};  % mx1 cell of images to Blt e.g. (mainstim*, background, stim2,...)
                       % dimensions: 
                       % background: rgb
    positions   = {};  % mx1 cell of 4x1 position of stimulus [xpos ypos width height]. 
    
    state       = zeros(2,1);      % current state vector 
    A           = eye(2,2);      % dynamics update matrix
    W           = zeros(2,1);      % dynamics noise
    
    cursorsteady = 0; % frames for which the cursor is steady

    % stillblob fixed parameters
    nonvarparams   = {'steady_thresh_frame' 'steady_thresh_deg' 'waitsecs' 'maxtrialtime', 'maxtrials'};
    steady_thresh_frame;    % if subject is near the target for this long, go to next trial
    steady_thresh_deg;      % if the pointer is within this threshold of target, count frame as "steady"
    waitsecs;               % allow movement after this many seconds
    maxtrialtime;
    maxtrials;
      
    % stillblob variable parameters
    varparams   = {'backLum' 'noiseLum' 'stimLum' 'stimStd'};
    backLum;
    noiseLum;
    stimLum;
    stimStd;
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
            p.addParameter('pos_start', [], @()(isnumeric(x)))
        else
            % otherwise generate test poitns around circles
            % radius can't be too big
            p.addParameter('r', [1, 2, 3, 4, 5, 8], @(x)(isnumeric(x) && all(x < 20))) 
            p.addParameter('angle', 0:pi/4:(2*pi-pi/4), @(x)(isnumeric(x)))
        end
        
        % time threshold for going into next trial in frames
        % user needs to be steady for this amount of frames.
        p.addParameter('steady_thresh_frame', floor(2 * myscreen.framesPerSecond), @(x)(isnumeric(x)))
        p.addParameter('steady_thresh_deg', 0.2, @(x)(isnumeric(x)))
        p.addParameter('waitsecs', 2, @(x)(isnumeric(x))) 
        p.addParameter('maxtrialtime', 10, @(x)(isnumeric(x))) % in seconds
        p.addParameter('maxtrials', 300, @(x)(isnumeric(x)))
        
        p.addParameter('backLum', 90, @(x)(isnumeric(x))) 
        p.addParameter('noiseLum', 0, @(x)(isnumeric(x))) 
        p.addParameter('stimLum', 255, @(x)(isnumeric(x))) 
        p.addParameter('stimStd', 0.4, @(x)(isnumeric(x))) 
                
        p.addParameter('randomize_order', true, @(x) (islogical(x)))
        p.parse(varargin{:})
                
        %% initialize
        % initial positions
        if isfield(p.Results, 'pos_start')
            assert(size(p.Results.pos_start,2)==2, "check initial positions") 
            obj.pos_start{1} = p.Results.pos_start;
            obj.numTrials = min(p.Results.maxtrials, size(p.Results.pos_start,1));
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
            obj.pos_start{1} = obj.pos_start{1}(randperm(stillblob.n),:);
        end
        
        for param = [varparams, nonvarparams]
            eval(['obj.' param{1} ' = p.Results.' param{1}])
        end
    end
    
    % return task object that can be run on trackpos.m
    function thistask    = configureExperiment(obj, stimulus, task, myscreen) 
        thistask        = struct();
        thistask.segmin = [obj.maxtrialtime];
        thistask.segmax = [obj.maxtrialtime];
        
        thistask.numTrials          = obj.numTrials;
        thistask.getResponse        = [0];
        thistask.synchToVol         = [0];
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
        obj.positions{1}   = [obj.pos_start(task.trialnum,2), [], []];
        
        % initialize background position
        if task.thistrial.noiseLum > 0
            nframes = myscreen.framesPerSecond*task.segmax(1) + 20;%/downsample_timeRes; 
            task.thistrial.bgpermute(1:nframes) = randi(length(stimulus.backnoise),nframes,1);
            obj.stimulus{2}  = stimulus.backnoise{task.thistrial.bgpermute(1)};
            obj.positions{2} = [0,0,myscreen.imageWidth, myscreen.imageHeight]]; 
        end
        
        % trial terminal conditions
        obj.cursor_steady = 0; 
        
    end
    
    function task = startSegment(obj, task, myscreen, stimulus)
        tic
    end

    % frame update
    % need to define an update function for the stimulus
    function task  = update(obj, task, myscreen, stimulus)       
        % background luminance
        mglClearScreen(task.thistrial.backLum/255);
        
        % blt stimuli
        for stimidx = 1:length(obj.stimulus)
            mglBltTexture(obj.stimulus{stimidx},obj.positions{stimidx})
        end
        
        % pointer updates
        if tic < obj.waitsecs
            task.thistrial.movecursor = false;
        else
            task.thistrial.movecursor = true;
        end
        
        % cursorsteady
        if sum(abs(obj.positions{1}(1:2) - stimulus.pointer)) < stimulus.stimStd/2
            obj.cursor_steady = obj.cursor_steady + 1;
        end

        if obj.cursor_steady > obj.steady_thresh
            task = jumpSegment(task); 
        end
        
        % update framecount
        task.thistrial.framecount = task.thistrial.framecount + 1;
        
        % update stimuli position
            
        % update background
        obj.stimulus{2}  = stimulus.backnoise{task.thistrial.bgpermute(task.thistrial.framecount+1)};        
    end
    
end

end