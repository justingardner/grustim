% abstract class for defining a tracking task
% defines how the stimulus should update during the task
% 

classdef (Abstract) trackingTask < handle
    properties (Abstract)
        % task parameters
        name;       % name of task
        numTrials;          % number of trials
        pos_start;  % starting position of stimulus         
                        
        % trial parameters
        stimulus;  % mx1 cell of images to Blt e.g. (mainstim*, background, stim2,...)
                           % dimensions: 
                           % background: rgb
        positions;  % mx1 cell of 4x1 position of stimulus [xpos ypos width height]. 
        
        state;      % current state vector 
        A;          % dynamics update matrix
        W;          % dynamics noise
        movecursor; % indicates whether we can move cursor
        
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
        function pos = returnstimpos(obj)
            pos = obj.positions(1,1:2); % first stimulus should be the main stimulus to track
        end
        
    end
    
end