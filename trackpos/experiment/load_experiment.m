function cps = load_experiment(myscreen, stimulusname, experiment, exp_num, varargin)
    getArgs(varargin, {'debugmode=false'});
    if strcmp(stimulusname, 'circular_ar')
        cps = load_experiment_circular_ar(myscreen, experiment, exp_num, debugmode);
    end
end


function cps = load_experiment_circular_ar(myscreen, experiment_name, exp_num, debugmode)
    cps = {};
    ntrial_learn        = 4;  % learning phase at full luminance, not analyzed
    nblocks             = 3;  % number of same blocks for each condition
    trials_per_block    = 5;  % number of trials per block
    maxtrialtime        = 20; % seconds

    if debugmode, ntrial_learn= 1; trials_per_block = 1; maxtrialtime=2; end

    if strcmp(experiment_name, 'mn')
        %% effect of ar dynamics
        experiment          = {'mn'};
        
        if exp_num == 1 
            experiment_paramset = [6, 5, 4, 3, 2, 1]; %., 9, 8, 7, 12,11,10];
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
    	
        if exp_num == 1 
            experiment_paramset = [10, 7, 11, 8, 12, 9]; 
            Nconds = length(experiment_paramset);

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
            
        elseif exp_num == 2
            experiment_paramset = [4,1,5,2,6,3]; %., 9, 8, 7, 12,11,10];
            Nconds = length(experiment_paramset);

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
    elseif strcmp(experiment_name, 'pert')         
        experiment          = {'pert'};

        if exp_num == 1
            experiment_paramset = [1,2,3,4,5,6]; %., 9, 8, 7, 12,11,10];
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

    end
end