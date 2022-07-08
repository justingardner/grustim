classdef joy_calib < trackingTask
    % todo: getresponse, play with buttons and joy2vel functional forms,
        % stimulus types, export settings, change stimulus.exp if needed
    % todo: automated calibration during tracking; x
    
properties
    % task parameter
    name            = 'joy_calib'
    numTrials       = 1;    % total number of trials
    pos_start       = cell(1,1);   % mx1 cell {(numTrials x 2), .... } of starting positions of stimuli
    
    % trial parameters
    stimulus        = cell(2,1);  % mx1 cell of images to Blt e.g. (mainstim*, background, stim2,...)
    positions       = cell(2,1);  % mx1 cell of 4x1 position of stimulus [xpos ypos width height]. 
    state           = nan(1,1);       % current state vector [target_x, pointer_x, target_y, pointer_y]
    A               = nan(1,1);         % dynamics update matrix
    W               = nan(1,1);       % dynamics noise
    cidx            = nan;
    
    movecursor      = 1;
    
    bgfile          = [];

    % fixed parameters
    nonvarparams   = {};    
    pointertype;  % dot, blob (feed in init stimulus struct)
      
    % stillblob variable parameters
    varparams   = {};
end
    
    
methods
    % Constructor
    function obj = joy_calib(pointertype)
        obj.pointertype = pointertype;
        obj.initialize_params(obj,p)
    end
    
    % return task object that can be run on trackpos.m
    function thistask = configureExperiment(obj, task, myscreen, stimulus) 
        thistask        = struct();
        task.seglen     = [inf]; 
        
        thistask.numTrials          = obj.numTrials;
        thistask.getResponse        = [1]; % get response to end trial
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
        if obj.pointertype == 'blob'
            obj.stimulus{1}    = trackpos_stim.gaussian;
            obj.positions{1}   = [obj.pos_start{1}(task.trialnum,:), [], []];
        elseif obj.pointertype == 'dot'
            obj.stimulus{1}    = 0;
            obj.positions{1}   = [obj.pos_start{1}(task.trialnum,:), [], []];
        end
        
        
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
    end
    
    function task = startSegment(obj, task, myscreen, stimulus)
        
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
                    
        % update background
        if ~isempty(obj.bgfile) && task.thistrial.noiseLum > 0
            obj.stimulus{2}  = stimulus.backnoise{task.thistrial.bgpermute(task.thistrial.framecount+1)};        
        end
        
        % update fixation
        if stimulus.exp.fixateCenter == 1 % fixation below others.
            mglGluAnnulus(0,0,0.2,0.3,[0 1 0],60,1);
            mglGluDisk(0,0,0.1,rand(1,3),60,1);
        end
    end

    end
    
end