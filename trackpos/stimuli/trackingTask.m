% abstract class for defining a tracking task
% defines how the stimulus should update during the task
% 

classdef (Abstract) trackingTask < handle
    properties (Abstract)
        % task parameters
        name;       % name of task
        numTrials;  % number of trials
        pos_start;  % starting position of stimulus         
                        
        % trial parameters        
        state;      % current state vector 
        A;          % dynamics update matrix
        W;          % dynamics noise
        tpidx;       % position index of target
        ppidx;       % position index of pointer
        tnidx;       % noise index of target
        pnidx;       % noise index of pointer
        
        movecursor; % indicates whether we can move cursor during this trial
        doTrack;    % indicates whether we should start recording tracking variables
        displayFix; % display fixation at current segmention
        
        bgfile;     % background file to preload at task initialization
        
        % task parameters
        nonvarparams;   % cell of parameter names to be included as parameters in mgl task.
        varparams;      % cell of parameter names to be included as parameters in mgl task.
    end

    methods (Abstract)
        % return task object that can be run on trackpos.m
        thistask    = configureExperiment(obj,task, myscreen, stimulus) 
        
        % trial update
        % if there are parameters being randomized by mgl on trial by trial
        % basis, this function should interact with it 
        task        = initTrial(obj,task, myscreen, stimulus)
        
        task        = startSegment(obj, task, myscreen, stimulus)

        % update function for the stimulus
        % set background luminance
        % blt all images (stimuli and background)
        % allow pointer to move
        % update framecount
        % update stimuli position
        % update background image
        task = update(obj, task, myscreen, stimulus)        
    end
    
    methods
        function initialize_params(obj, parserOut)
            for param = parserOut.Parameters
                eval(['obj.' param{1} ' = parserOut.Results.' param{1} ';'])
            end
        end
    end
end