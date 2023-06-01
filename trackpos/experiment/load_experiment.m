function cps = load_experiment(myscreen, stimulusname, experiment, exp_num, varargin)
    getArgs(varargin, {'debugmode=false', 'shuffle_set=false'});
    if strcmp(stimulusname, 'circular_ar')
        cps = load_experiment_circular_ar(myscreen, experiment, exp_num, shuffle_set,debugmode);
    elseif strcmp(stimulusname, 'linear_ar')
        cps = load_experiment_linear(myscreen, experiment, exp_num, shuffle_set,debugmode);
    end
end

function cps = load_experiment_linear(myscreen, experiment_name, exp_num, shuffle_set,debugmode);
    blocked             = true; % make things blocked
    
    cps = {};
    ntrial_learn        = 4;  % learning phase at full luminance, not analyzed
    nblocks_learn       = 1;
    nblocks_per_cond    = 3;  % number of same blocks for each condition
    trials_per_block    = 5;  % number of trials per block
    maxtrialtime        = 20; % seconds

    if debugmode, ntrial_learn= 2; trials_per_block = 1; maxtrialtime=2; end

    experiment          = {experiment_name};

    if strcmp(experiment_name, 'pa')        
        if exp_num == 1 % radial
            experiment_paramset = [3,4,8,12]; 
        elseif exp_num == 2 % tangential
            experiment_paramset = [1,6,13,10]; 
        elseif exp_num == 10 % horizontal tracking
            experiment_paramset = 1:3:12; 
        elseif exp_num == 20 % diagonal tracking
            experiment_paramset = 2:3:12;
        elseif exp_num == 30 % vertical tracking
            experiment_paramset = 3:3:12;
        end
        
        if shuffle_set
            experiment_paramset = experiment_paramset(randperm(length(experiment_paramset)));
        end

        Nconds              = length(experiment_paramset);

        if blocked
            for epset = experiment_paramset
                for n = 1:nblocks_learn
                    % learning phase -- not analyzed
                    cps{end+1} = linear_ar(myscreen, 'numTrials', ntrial_learn, 'maxtrialtime', maxtrialtime, ...
                        'dyn_noise_phase', length(cps)+1,  'switch_tpnoise', false, ...
                        'experiment', experiment, 'experiment_paramset', epset);
                end

                % tracking
                for b = 1:nblocks_per_cond
                    cps{end+1} = linear_ar(myscreen, 'numTrials', trials_per_block, 'maxtrialtime', maxtrialtime, ...
                        'dyn_noise_phase', length(cps)+1,  'switch_tpnoise', false, ...
                        'experiment', experiment, 'experiment_paramset', epset);
                end
            end
        else
            nblocks             = nblocks_per_cond * Nconds;

            for n = 1:nblocks_learn
                % learning phase -- not analyzed
                cps{end+1} = linear_ar(myscreen, 'numTrials', ntrial_learn, 'maxtrialtime', maxtrialtime, ...
                    'dyn_noise_phase', length(cps)+1,  'switch_tpnoise', false, ...
                    'experiment', experiment, 'experiment_paramset', experiment_paramset);
            end

            % tracking
            for b = 1:nblocks
                cps{end+1} = linear_ar(myscreen, 'numTrials', trials_per_block, 'maxtrialtime', maxtrialtime, ...
                    'dyn_noise_phase', length(cps)+1,  'switch_tpnoise', false, ...
                    'experiment', experiment, 'experiment_paramset', experiment_paramset);
            end
        end
    end

end

function cps = load_experiment_circular_ar(myscreen, experiment_name, exp_num, shuffle_set, debugmode)
    cps = {};
    maxtrialtime        = 20; % seconds

    if debugmode, ntrial_learn= 1; nblocks=1; trials_per_block = 1; maxtrialtime=2; end

    if strcmp(experiment_name, 'mn')
        %% effect of ar dynamics
        experiment          = {'mn'};
        ntrial_learn        = 4;  % learning phase at full luminance, not analyzed
        nblocks             = 3;  % number of same blocks for each condition
        trials_per_block    = 5;  % number of trials per block

        
        if exp_num == 1 
            experiment_paramset = [6, 5, 4, 3, 2, 1]; %., 9, 8, 7, 12,11,10];
            
            if shuffle_set
                experiment_paramset = experiment_paramset(randperm(length(experiment_paramset)));
            end
            Nconds              = length(experiment_paramset);

            for epset = experiment_paramset
                % learning phase -- not analyzed
                cps{end+1} = circular_ar(myscreen, 'numTrials', ntrial_learn, 'maxtrialtime', maxtrialtime, ...
                    'dyn_noise_phase', length(cps)+1,  'switch_tpnoise', false, ...
                    'experiment', experiment, 'experiment_paramset', epset);

                % tracking
                for b = 1:nblocks
                    cps{end+1} = circular_ar(myscreen, 'numTrials', trials_per_block, 'maxtrialtime', maxtrialtime, ...
                        'dyn_noise_phase', length(cps)+1,  'switch_tpnoise', false, ...
                        'experiment', experiment, 'experiment_paramset', epset);
                end
            end
        end
    elseif strcmp(experiment_name, 'ecc') 
        %% effect of eccentricity
        experiment          = {'ecc'};
        ntrial_learn        = 4;  % learning phase at full luminance, not analyzed
        nblocks             = 3;  % number of same blocks for each condition
        trials_per_block    = 5;  % number of trials per block
    	
        if exp_num == 1 
            experiment_paramset = [10, 7, 11, 8, 12, 9]; 
            
            if shuffle_set
                experiment_paramset = experiment_paramset(randperm(length(experiment_paramset)));
            end
            Nconds = length(experiment_paramset);

            for epset = experiment_paramset
                % learning phase 
                cps{end+1} = circular_ar(myscreen, 'numTrials', ntrial_learn, 'maxtrialtime', maxtrialtime, ...
                    'dyn_noise_phase', length(cps)+1,  'switch_tpnoise', false, ...
                    'experiment', experiment, 'experiment_paramset', epset);

                % tracking
                for b = 1:nblocks
                    cps{end+1} = circular_ar(myscreen, 'numTrials', trials_per_block, 'maxtrialtime', maxtrialtime, ...
                        'dyn_noise_phase', length(cps)+1,  'switch_tpnoise', false, ...
                        'experiment', experiment, 'experiment_paramset', epset);
                end
            end
            
        elseif exp_num == 2
            experiment_paramset = [4,1,5,2,6,3]; %., 9, 8, 7, 12,11,10];
            
            if shuffle_set
                experiment_paramset = experiment_paramset(randperm(length(experiment_paramset)));
            end
            Nconds = length(experiment_paramset);

            for epset = experiment_paramset
                % learning phase 
                cps{end+1} = circular_ar(myscreen, 'numTrials', ntrial_learn, 'maxtrialtime', maxtrialtime, ...
                    'dyn_noise_phase', length(cps)+1,  'switch_tpnoise', false, ...
                    'experiment', experiment, 'experiment_paramset', epset);

                % tracking
                for b = 1:nblocks
                    cps{end+1} = circular_ar(myscreen, 'numTrials', trials_per_block, 'maxtrialtime', maxtrialtime, ...
                        'dyn_noise_phase', length(cps)+1,  'switch_tpnoise', false, ...
                        'experiment', experiment, 'experiment_paramset', epset);
                end
            end
        end
    elseif strcmp(experiment_name, 'pert')         
        %% perturbation test
        experiment          = {'pert'};
        
        ntrial_learn        = 4;  % learning phase at full luminance, not analyzed
        nblocks             = 3;  % number of same blocks for each condition
        trials_per_block    = 5;  % number of trials per block
        
        if exp_num == 1
            experiment_paramset = [1,2,3,4,5,6]; %., 9, 8, 7, 12,11,10];
            
            if shuffle_set
                experiment_paramset = experiment_paramset(randperm(length(experiment_paramset)));
            end
            Nconds = length(experiment_paramset);
            
            setnums         = [];
            phasenums       = []; currphasenum = 0;
            blocknums       = [];

            for epset = [1,2,3,4,5,6]
                currphasenum = currphasenum + 1; 
                setnums(end+1)   = epset;
                phasenums(end+1) = currphasenum;
                blocknums(end+1) = 0;
                % learning phase -- not analyzed
                cps{end+1} = circular_ar(myscreen, 'numTrials', ntrial_learn, 'maxtrialtime', maxtrialtime, ...
                    'dyn_noise_phase', currphasenum, 'switch_tpnoise', false, ...
                    'experiment', experiment, 'experiment_paramset', epset);

                % tracking
                for b = 1:nblocks
                    currphasenum        = currphasenum + 1; 
                    setnums(end+1)      = epset;
                    phasenums(end+1)    = currphasenum;
                    blocknums(end+1)    = b;
                    
                    cps{end+1} = circular_ar(myscreen, 'numTrials', trials_per_block, 'maxtrialtime', maxtrialtime, ...
                        'dyn_noise_phase', currphasenum,  'switch_tpnoise', false, ...
                        'experiment', experiment, 'experiment_paramset', epset);
                end
            end
            
            % switch 
            for phase = 1:length(phasenums)
                if blocknums(phase) == 0
                    continue
                end
                
                if setnums(phase) == 5
                    idx = (setnums == 2) & (blocknums == blocknums(phase));
                    cps{phase}.dyn_noise_phase = find(idx);
                    cps{phase}.switch_tpnoise = true;
                elseif setnums(phase) == 6
                    idx = (setnums == 4) & (blocknums == blocknums(phase));
                    cps{phase}.dyn_noise_phase = find(idx);
                    cps{phase}.switch_tpnoise = true;
                end
            end
        end
        
    elseif strcmp(experiment_name, 'ind')
        %% test target pointer independence
        experiment          = {'ind'};
        
        nblocks_learn       = 3;
        ntrial_learn        = 4;  % learning phase at full luminance, not analyzed
        nblocks             = 5;
        trials_per_block    = 5;  % number of trials per block
        maxtrialtime        = 20; % seconds
    
        if debugmode, nblocks_learn=1; ntrial_learn= 1; nblocks=1; trials_per_block = 1; maxtrialtime=2; end

        if exp_num == 1 
            experiment_paramset = [1,2,3,4];
            if shuffle_set
                experiment_paramset = experiment_paramset(randperm(length(experiment_paramset)));
            end
            Nconds              = length(experiment_paramset);

            for epset = experiment_paramset
                % learning phase -- not analyzed
                for b = 1:nblocks_learn
                    cps{end+1} = circular_ar(myscreen, 'numTrials', ntrial_learn, 'maxtrialtime', maxtrialtime, ...
                        'dyn_noise_phase', length(cps)+1,  'switch_tpnoise', false, ...
                        'experiment', experiment, 'experiment_paramset', epset);
                end

                % tracking
                for b = 1:nblocks
                    cps{end+1} = circular_ar(myscreen, 'numTrials', trials_per_block, 'maxtrialtime', maxtrialtime, ...
                        'dyn_noise_phase', length(cps)+1,  'switch_tpnoise', false, ...
                        'experiment', experiment, 'experiment_paramset', epset);
                end
            end
        end
    end
end