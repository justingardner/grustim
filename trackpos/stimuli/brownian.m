classdef brownian < trackingTask
    %todo: checkoob
    
properties
    % task parameter
    name            = 'brownian'
    numTrials       = 0;    % total number of trials
    pos_start       = cell(3,1);   % mx1 cell {(numTrials x 2), .... } of starting positions of stimuli
    
    % trial parameters
    stimulus        = cell(3,1);  % mx1 cell of images to Blt e.g. (mainstim*, background, stim2,...)
                       % dimensions: 
                       % background: rgb
    positions       = cell(3,1);  % mx1 cell of 4x1 position of stimulus [xpos ypos width height]. 
    
    % [target_x, pointer_x, ...(higher order target representation)..., target_y, pointer_y]
    state;       % current state vector 
    A;           % dynamics update matrix
    W;           % dynamics noise
    pidx;        % position index of target
    cidx;        % controllable states   

    movecursor      = 0;    
    t0              = 0;
    bgfile          = '/Users/gru/data/trackpos/trackpos.mat';

    % fixed parameters
    nonvarparams   = {'iti', 'maxtrialtime', 'trialpause' 'maxtrials'};
    iti;
    maxtrialtime;
    trialpause; % need to press backtick after each trial
    maxtrials;
      
    % variable parameters
    varparams   = {'backLum' 'noiseLum' 'stimLum' 'stimColor' 'stimStd' 'stepStd' 'dynamics_order' ...
                   'pointLum' 'pointColor' 'pointStd' 'pointStepStd'};
    backLum;
    noiseLum;
    stimLum;
    stimStd;
    stepStd;
    stimColor; 
    dynamics_order;
    pointLum;
    pointStd;
    pointStepStd;
    pointColor;
end
    
    
methods
    % Constructor
    function obj = brownian(myscreen, varargin)
        % {pos_start} or {r and angle}
        % randomize_order
        
        %% parse inputs
        p = inputParser;
        p.KeepUnmatched = true;
        
        if any(cellfun(@(x)(ischar(x) && strcmp(x,'pos_start')), varargin))
            % if start_pos exist as a parameter
            p.addParameter('pos_start', cell(1,1), @(x)(isnumeric(x)))
        end
        p.addParameter('iti', 2, @(x)(isnumeric(x))) 
        p.addParameter('maxtrials', 10, @(x)(isnumeric(x))) 
        p.addParameter('maxtrialtime', 15, @(x)(isnumeric(x))) 
        p.addParameter('trialpause', false, @(x)(islogical(x))) 
       
        p.addParameter('backLum', 90, @(x)(isnumeric(x))) ;
        p.addParameter('noiseLum', 0, @(x)(isnumeric(x))) ;
        p.addParameter('stimLum', 255, @(x)(isnumeric(x))) ;
        p.addParameter('stimColor', 'k');
        p.addParameter('stimStd', 0.4, @(x)(isnumeric(x))) ;
        p.addParameter('stepStd', 1, @(x)(isnumeric(x))) ;
        p.addParameter('dynamics_order', 0, @(x)(isinteger(x))) ;
        p.addParameter('pointLum', 255, @(x)(isnumeric(x)));
        p.addParameter('pointColor', 'r'); %todo: add checks
        p.addParameter('pointStd', 0, @(x)(isnumeric(x)));
        p.addParameter('pointStepStd', 0, @(x)(isnumeric(x)));
                
        p.parse(varargin{:})
                
        %% initialize
        % initial positions
        if isfield(p.Results, 'pos_start')
            obj.pos_start = p.Results.pos_start;
            obj.numTrials = min(p.Results.maxtrials, size(p.Results.pos_start{1},1));
        else
            obj.numTrials = p.Results.maxtrials;
            x_img = 1/3*myscreen.imageWidth*(2*rand(obj.numTrials,1)-1); 
            y_img = 1/3*myscreen.imageHeight*(2*rand(obj.numTrials,1)-1);
            obj.pos_start{1} = [x_img, y_img];
        end
        
        obj.initialize_params(p)
        
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
        
        if stimulus.exp.debug
            thistask.numTrials          = 5;
        else
            thistask.numTrials          = obj.numTrials;
        end
        
        thistask.getResponse        = [0,0];
        if obj.trialpause
            thistask.synchToVol     = [0, 1];
        else
            thistask.synchToVol     = [0, 0];
        end
        thistask.waitForBacktick    = 1;
        
        for param = obj.varparams
            eval(['thistask.parameter.' param{1} ' = obj.' param{1} ';'])
        end
        
        obj.turnOffDefaultPointer();
%         if obj.pointLum > 0 && obj.pointStd > 0
%             obj.turnOffDefaultPointer();
%         end
    end
    
    % trial update?
    function task  = initTrial(obj, task, myscreen, stimulus)
        
        % initialize stimulus position
        target = struct();
        for param = {'stimStd', 'stepStd', 'stimLum', 'stimColor'}
           eval(['target.' param{1} ' = task.thistrial.' param{1} ';'])
        end
        trackpos_stim      = trackposInitStimulus(target,myscreen);
        obj.stimulus{1}    = trackpos_stim.gaussian;
        obj.positions{1}   = [obj.pos_start{1}(task.trialnum,:), [], []];
        
        % initialize background position
        if ~isempty(obj.bgfile) && task.thistrial.noiseLum > 0
            nframes = myscreen.framesPerSecond*task.segmax(1) + 20;%/downsample_timeRes; 
            task.thistrial.bgpermute(1:nframes) = randi(length(stimulus.backnoise),nframes,1);
            obj.stimulus{2}  = stimulus.backnoise{1}{task.thistrial.bgpermute(1)};
            obj.positions{2} = [0,0,myscreen.imageWidth, myscreen.imageHeight]; 
        end
        
        % initialize the pointer
        pointer             = struct();
%         for param = {'pointStd', 'pointStepStd', 'pointLum', 'pointColor'}
%            eval(['pointer.stim' param{1}[5:] ' = task.thistrial.' param{1} ';'])
%         end
        pointer.stimLum     = task.thistrial.pointLum;
        pointer.stimStd     = task.thistrial.pointStd;
        pointer.stepStd     = task.thistrial.pointStepStd;
        pointer.stimColor   = task.thistrial.pointColor;
        pointer_stim        = trackposInitStimulus(pointer,myscreen);
        obj.stimulus{3}     = pointer_stim.gaussian;
        obj.positions{3}    = [obj.pos_start{1}(task.trialnum,:), [], []];
        
        % set mouse position to the stimulus position. 
        x_img = obj.pos_start{1}(task.trialnum,1);  y_img = obj.pos_start{1}(task.trialnum,2);
        x_screen = x_img*myscreen.screenWidth/myscreen.imageWidth + myscreen.screenWidth/2;
        y_screen = y_img*myscreen.screenHeight/myscreen.imageHeight + myscreen.screenHeight/2;
        mglSetMousePosition(ceil(x_screen),floor(y_screen), myscreen.screenNumber); % correct for screen resolution???
        if ~stimulus.exp.showMouse, mglDisplayCursor(0);, end %hide cursor
        
        % keep track of start trial        
        obj.t0 = tic; 
        
        % define dynamics
        ss = 2 + (task.thistrial.dynamics_order + 1)*2;
        do = task.thistrial.dynamics_order;
        
        obj.state = [obj.positions{1}(1:2)'; zeros(1,do*2); obj.positions{1}(1:2)']; % current state vector 
        obj.A     = blkdiag(triu(ones(do+1)), triu(ones(do+1)), zeros(2,2));         % dynamics update matrix
        obj.pidx  = [1, 2+do];
        obj.cidx  = [ss-1,ss];
        obj.W     = zeros(ss,4); % dynamics noise
        obj.W((do+1),1)   = task.thistrial.stepStd/sqrt(myscreen.framesPerSecond);
        obj.W(2*(do+1),2) = task.thistrial.stepStd/sqrt(myscreen.framesPerSecond);
        obj.W(2*(do+1)+1,3) = task.thistrial.pointStepStd/sqrt(myscreen.framesPerSecond);
        obj.W(2*(do+1)+2,4) = task.thistrial.pointStepStd/sqrt(myscreen.framesPerSecond);
    end
  
    
    function task = startSegment(obj, task, myscreen, stimulus)
        if task.thistrial.thisseg == 1
            obj.movecursor = true;
        else
            obj.movecursor = false;
            task.thistrial.framecount = []; % don't write into calculated variables
        end   
    end

    % frame update
    % need to define an update function for the stimulus
    function task  = update(obj, task, myscreen, stimulus)       
        % background luminance
        mglClearScreen(task.thistrial.backLum/255);
        
        % tracking during segment 1
        if task.thistrial.thisseg == 1
            % blt all stimuli (including background)
            for stimidx = 1:length(obj.stimulus)
                if isstruct(obj.stimulus{stimidx}) % if we get null image, return a red annulus
                    mglBltTexture(obj.stimulus{stimidx}, obj.positions{stimidx})
                end
            end
            
            % overlay red disk
%             for stimidx = 1:length(obj.stimulus)
%                 if ~isstruct(obj.stimulus{stimidx}) % if we get null image, return a red annulus
%                     mglGluDisk(obj.positions{stimidx}(1), obj.positions{stimidx}(2), 0.2, [1 0 0])
%                 end
%             end
            
            % update stimuli position
            obj.updateStimulus(myscreen,stimulus)

            % update background
            if ~isempty(obj.bgfile) && task.thistrial.noiseLum > 0
                obj.stimulus{2}  = stimulus.backnoise{1}{task.thistrial.bgpermute(task.thistrial.framecount+1)};        
            end
        else
            task.thistrial.framecount = [];
        end
        
        % update fixation
        if stimulus.exp.fixateCenter == 1 % fixation below others.
            mglGluAnnulus(0,0,0.2,0.3,[1 1 1],60,1);
            mglGluDisk(0,0,0.1,rand(1,3),60,1);
        end
    end
     
    function updateStimulus(obj, myscreen, stimulus)
        % update state
        noise       = obj.W * normrnd(0,1,size(obj.W,2),1);
        newstate    = obj.A * obj.state + noise;
        
        % update mouse
        if strcmp(stimulus.exp.controlMethod, 'mouse')
            mInfo = mglGetMouse(myscreen.screenNumber);
            stimulus.pointer(1) = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
            stimulus.pointer(2) = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;
            
            % perturb current position
            newstate(obj.cidx) = stimulus.pointer' + noise(obj.cidx);
            
            % pointer: subtract back if out of bounds
            [horz_out, vert_out] = check_oob(newstate(obj.cidx), myscreen, stimulus);    
            newstate(obj.cidx(1))  = newstate(obj.cidx(1)) - horz_out * noise(obj.cidx(1));
            newstate(obj.cidx(2))  = newstate(obj.cidx(2)) - vert_out * noise(obj.cidx(2));
            stimulus.pointer       = newstate(obj.cidx)';
            
            x_img = stimulus.pointer(1);  y_img = stimulus.pointer(2);
            x_screen = x_img*myscreen.screenWidth/myscreen.imageWidth + myscreen.screenWidth/2;
            y_screen = y_img*myscreen.screenHeight/myscreen.imageHeight + myscreen.screenHeight/2;
            mglSetMousePosition(ceil(x_screen),floor(y_screen), myscreen.screenNumber); % correct for screen resolution???
        else
            newstate(obj.cidx) = stimulus.pointer' + noise(obj.cidx);
            % pointer: subtract back if out of bounds
            [horz_out, vert_out] = check_oob(newstate(obj.cidx), myscreen, stimulus);    
            newstate(obj.cidx(1))  = newstate(obj.cidx(1)) - horz_out * noise(obj.cidx(1));
            newstate(obj.cidx(2))  = newstate(obj.cidx(2)) - vert_out * noise(obj.cidx(2));
            stimulus.pointer       = newstate(obj.cidx)';
        end
        
        % target: subtract back if out of bounds
        [horz_out, vert_out] = check_oob(newstate(obj.pidx), myscreen, stimulus);    
        newstate(obj.pidx(1))  = newstate(obj.pidx(1)) - horz_out * noise(obj.pidx(1));
        newstate(obj.pidx(2))  = newstate(obj.pidx(2)) - vert_out * noise(obj.pidx(2));
        
        % update position
        obj.positions{1}(1:2)   = newstate(obj.pidx); 
        obj.positions{3}(1:2)   = stimulus.pointer;
        obj.state               = newstate;
    end

end

methods(Static)
    function turnOffDefaultPointer
        global stimulus;
        if stimulus.exp.dispPointer
            stimulus.exp.dispPointer = 0;
        end
    end
end

end