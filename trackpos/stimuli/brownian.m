classdef brownian < trackingTask
    %todo: checkoob
    
properties
    % task parameter
    name            = 'brownian'
    numTrials       = 0;    % total number of trials
    pos_start       = cell(3,1);   % mx1 cell {(numTrials x 2), .... } of starting positions of stimuli
    bgfile          = '/Users/gru/data/trackpos/trackpos.mat';
    
    % trial parameters    
    state;       % current state vector 
    A;           % dynamics update matrix
    W;           % dynamics noise matrix
    pidx;        % position index of target
    cidx;        % position index of pointer

    movecursor      = false;
    doTrack         = true;    % indicates whether we should start recording tracking variables

    % fixed parameters
    nonvarparams   = {'iti', 'maxtrialtime', 'trialpause' 'maxtrials'};
    iti;
    maxtrialtime;
    trialpause; % need to press backtick after each trial
    maxtrials;
      
    % variable parameters
    varparams   = {'backLum' 'noiseLum' 'stimLum' 'stimColor' 'stimStd' 'stimStepStd' 'dynamics_order' ...
                   'pointLum' 'pointColor' 'pointStd' 'pointStepStd'};
    backLum;
    noiseLum;
    stimLum;
    stimStd;
    stimStepStd;
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
        p.addParameter('bgfile', '/Users/gru/data/trackpos/trackpos.mat') 
       
        p.addParameter('backLum', 90, @(x)(isnumeric(x))) ;
        p.addParameter('noiseLum', 0, @(x)(isnumeric(x))) ;
        p.addParameter('stimLum', 255, @(x)(isnumeric(x))) ;
        p.addParameter('stimColor', 'k');
        p.addParameter('stimStd', 0.4, @(x)(isnumeric(x))) ;
        p.addParameter('stimStepStd', 1, @(x)(isnumeric(x))) ;
        p.addParameter('dynamics_order', 0, @(x)(isinteger(x))) ;
        p.addParameter('pointLum', 255, @(x)(isnumeric(x)));
        p.addParameter('pointColor', 'r'); %todo: add checks
        p.addParameter('pointStd', 0, @(x)(isnumeric(x)));
        p.addParameter('pointStepStd', 0, @(x)(isnumeric(x)));
                
        p.parse(varargin{:})
        
        obj.initialize_params(p)
                
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

        %initialize states
        obj.state = nan(4,1);       % current state vector 
        obj.A     = nan(4,4);       % dynamics update matrix
        obj.W     = nan(4,2);       % dynamics noise
        obj.pidx  = nan(2,1);       % controllable states   
        obj.cidx  = nan(2,1);       % controllable states   
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
    function stimulus = initTrial(obj, task, myscreen, stimulus)
        % initialize stimulus position
        target = struct();
        for param = {'stimStd', 'stimStepStd', 'stimLum', 'stimColor'}
           eval(['target.' lower(param{1}(5:end)) ' = task.thistrial.' param{1} ';'])
        end
        stimulus.target             = trackposInitStimulus(target,myscreen);
        initpos = obj.pos_start{1}(task.trialnum,:);
        stimulus.target.position    = initpos;
        
        % initialize the pointer
        pointer             = struct();
        for param = {'pointStd', 'pointStepStd', 'pointLum', 'pointColor'}
           eval(['pointer.stim' lower(param{1}(6:end)) ' = task.thistrial.' param{1} ';'])
        end
        stimulus.pointer            = trackposInitStimulus(pointer,myscreen);
        stimulus.pointer.position   = initpos; % initialize pointer position to the stimulus position
        
        % define dynamics
        % state: target position (x,y), target velocity + (higher order)
        % (x,y), pointer position (x,y) and velocity
        do = task.thistrial.dynamics_order;
        ss = (do + 1)*2 + 4; % state size
        
        % todo: consider dt in dynamics matrix
        A11         = triu(ones(do+1)) .* triu(ones(do+1),-1)'; % first order approximation to position/vel/etc...
        A11(end,:)  = zeros(1,do+1);
        Aee         = triu(ones(2));
        Aee(end,:)  = zeros(1,2);
        obj.A     = blkdiag(A11, A11, Aee, Aee); % dynamics update matrix
        obj.pidx  = [1, 2+do]; % index of target position
        obj.cidx  = [ss-3,ss-1]; % index of pointer position
        obj.W     = zeros(ss,4); % dynamics noise
        obj.W((do+1),1)     = task.thistrial.stimStepStd/sqrt(myscreen.framesPerSecond); % last element affected by noise
        obj.W(2*(do+1),2)   = task.thistrial.stimStepStd/sqrt(myscreen.framesPerSecond);
        obj.W(2*(do+1)+2,3) = task.thistrial.pointStepStd/sqrt(myscreen.framesPerSecond);
        obj.W(2*(do+1)+4,4) = task.thistrial.pointStepStd/sqrt(myscreen.framesPerSecond);
        
        % initialize state
        obj.state = zeros(ss,1);
        obj.state(obj.pidx) = initpos;
        obj.state(obj.cidx) = initpos;
    end
  
    
    function task = startSegment(obj, task, myscreen, stimulus)
        if task.thistrial.thisseg == 1 % tracking
            obj.movecursor  = true;
            obj.doTrack     = true;
        else % ITI
            obj.movecursor = false;
            obj.doTrack    = false;
        end   
    end

    % frame update
    % need to define an update function for the stimulus
    function stimulus  = update(obj, task, myscreen, stimulus)       
        % background luminance
        mglClearScreen(task.thistrial.backLum/255);
        
        % tracking during segment 1
        if task.thistrial.thisseg == 1
            % update stimuli position
            stimulus = obj.updateStimulus(myscreen,stimulus);
        end
    end
     
    function stimulus = updateStimulus(obj, myscreen, stimulus)
        % update state
        noise                       = obj.W * normrnd(0,1,size(obj.W,2),1);
        newstate                    = obj.A * obj.state + noise;
        
        % update pointer
        % perturb current position
        pointer_newpos              = stimulus.pointer.position' + noise(obj.cidx);

        % pointer: subtract back if out of bounds
        [horz_out, vert_out]        = check_oob(pointer_newpos, myscreen, stimulus.pointer.std);    
        newstate(obj.cidx(1))       = pointer_newpos(1) - horz_out * noise(obj.cidx(1));
        newstate(obj.cidx(2))       = pointer_newpos(2) - vert_out * noise(obj.cidx(2));
        stimulus.pointer.position   = newstate(obj.cidx)';

        % target: subtract back if out of bounds
        [horz_out, vert_out]        = check_oob(newstate(obj.pidx), myscreen, stimulus.target.std);    
        newstate(obj.pidx(1))       = newstate(obj.pidx(1)) - horz_out * noise(obj.pidx(1));
        newstate(obj.pidx(2))       = newstate(obj.pidx(2)) - vert_out * noise(obj.pidx(2));
        stimulus.target.position    = newstate(obj.pidx)';
        
        % update position
        obj.state                   = newstate;
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