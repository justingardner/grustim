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
    tpidx;       % position index of target
    ppidx;       % position index of pointer
    tnidx;       % noise index of target
    pnidx;       % noise index of pointer

    movecursor      = false;
    doTrack         = true;    % indicates whether we should start recording tracking variables
    displayFix      = true; % display fixation at current segmention

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
        p.addParameter('dynamics_order', 1, @(x)(isinteger(x))) ;
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
        obj.tpidx  = nan(2,1);       % controllable states   
        obj.ppidx  = nan(2,1);       % controllable states   
        obj.tnidx  = nan(2,1);       % controllable states   
        obj.pnidx  = nan(2,1);       % controllable states   
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
           eval(['pointer.' lower(param{1}(6:end)) ' = task.thistrial.' param{1} ';'])
        end
        stimulus.pointer            = trackposInitStimulus(pointer,myscreen);
        stimulus.pointer.position   = initpos; % initialize pointer position to the stimulus position
        
        obj.displayFix = false;
        
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
        obj.tpidx  = [1, 2+do];      % index of target position
        obj.ppidx  = [ss-3,ss-1];    % index of pointer position
        obj.tnidx  = [do+1, 2*(do+1)]; % target noise index
        obj.pnidx  = [2*(do+1)+2,2*(do+1)+4]; % pointer noise index
        obj.W     = zeros(ss,4); % dynamics noise
        obj.W(obj.tnidx(1),1)     = task.thistrial.stimStepStd/sqrt(myscreen.framesPerSecond); % last element affected by noise
        obj.W(obj.tnidx(2),2)   = task.thistrial.stimStepStd/sqrt(myscreen.framesPerSecond);
        obj.W(obj.pnidx(1),3) = task.thistrial.pointStepStd/sqrt(myscreen.framesPerSecond);
        obj.W(obj.pnidx(2),4) = task.thistrial.pointStepStd/sqrt(myscreen.framesPerSecond);
        
        % initialize state
        obj.state = zeros(ss,1);
        obj.state(obj.tpidx) = initpos;
        obj.state(obj.ppidx) = initpos;
    end
  
    
    function stimulus = startSegment(obj, task, myscreen, stimulus)
        if task.thistrial.thisseg == 1 % tracking
            obj.movecursor  = true;
            obj.doTrack     = true;
            obj.displayFix 	= true;
        else % ITI
            obj.movecursor  = false;
            obj.doTrack     = false;
            obj.displayFix 	= true;
            
            stimulus.target = [];
            stimulus.pointer = [];
        end   
    end

    % frame update
    % need to define an update function for the stimulus
    function stimulus  = update(obj, task, myscreen, stimulus)       
        % tracking during segment 1
        if task.thistrial.thisseg == 1
            % update stimuli position
            stimulus = obj.updateStimulus(myscreen,stimulus);
        end
    end
     
    function stimulus = updateStimulus(obj, myscreen, stimulus)
        % update state
        noise                       = obj.W * normrnd(0,1,size(obj.W,2),1);
        obj.state(obj.ppidx(1))     = stimulus.pointer.position(1);
        obj.state(obj.ppidx(2))     = stimulus.pointer.position(2);
        
        nanindex = isnan(obj.state);
        if any(nanindex)
            obj.state(nanindex) = 0; % replace NaN with 0s for proper matrix multiplication
        end
        
        newstate = obj.A * obj.state + noise;
        
        % update pointer
        % pointer: subtract back if out of bounds
        if any(noise(obj.pnidx)~=0) && all(~nanindex(obj.ppidx)) % if we added noise.
            [horz_out, vert_out]        = check_oob(newstate(obj.ppidx), myscreen, stimulus.pointer.std);
            newstate(obj.ppidx(1))      = (1-horz_out) * newstate(obj.ppidx(1)) + horz_out * obj.state(obj.ppidx(1));
            newstate(obj.ppidx(2))      = (1-vert_out) * newstate(obj.ppidx(2)) + vert_out * obj.state(obj.ppidx(2));
        end
%       disp(['pointer pos (brownian): ' num2str(stimulus.pointer.position)]) % todo: delete this line after checking

        % target: subtract back if out of bounds
        if any(noise(obj.tnidx)~=0) && all(~nanindex(obj.tpidx))
            [horz_out, vert_out]        = check_oob(newstate(obj.tpidx), myscreen, stimulus.target.std);    
            newstate(obj.tpidx(1))      = (1-horz_out) * newstate(obj.tpidx(1)) + horz_out * obj.state(obj.tpidx(1));
            newstate(obj.tpidx(2))      = (1-vert_out) * newstate(obj.tpidx(2)) + vert_out * obj.state(obj.tpidx(2));
        end
        
        % replace nans back
        newstate(nanindex) = NaN;
        
        % update position
        stimulus.pointer.position   = newstate(obj.ppidx)';
        stimulus.target.position    = newstate(obj.tpidx)';
        obj.state                   = newstate;
    end

end

methods(Static)

end

end