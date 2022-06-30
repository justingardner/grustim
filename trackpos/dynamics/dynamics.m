classdef dynamics
    properties
        % trackpos properties
        backLum         = 90;
        stimLum         = 255; %full luminance
        stimStd         = [0.4]; % size of stimulus
        time            = [30]; % trial maximum time.
        maxtrials       = [5]; % maxmimum number of trials
        randomize_order = true; % randomize different conditions
        
        x0 % initial state

        % motorcalib properties
        steady_thresh
        waitsecs
    end
    
    methods
        function new_state = update(state, n, curr_task, myscreen,stimulus)
            error('NotImplemented')
            new_state = [];
        end
    end
end