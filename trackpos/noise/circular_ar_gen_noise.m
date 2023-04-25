function dyn_noise = circular_ar_gen_noise(varargin)
    getArgs(varargin, ...
        {'taskcfg', 'not_specified',...
        'myscreen', struct('framesPerSecond',60), ...
        'savenoise',false, ...
        'experiment', 'pert', 'exp_num',1});
    
    if strcmp(taskcfg,'not_specified')
        taskcfg = load_experiment(myscreen, 'circular_ar', experiment, exp_num, 'debugmode', false);
    end
    
    dyn_noise = generate_noise_struct(myscreen, taskcfg);
    if savenoise
        savefile = sprintf('circular_ar_%s_%s.mat',experiment, num2str(exp_num));
        save(fullfile(find_root_dir, 'proj/grustim/trackpos/noise', savefile), 'dyn_noise');
    end
end

function noiseStruct = generate_noise_struct(myscreen, taskcfg)
    for n = 1:length(taskcfg)
        thistask    = taskcfg{n};
        T           = ceil(thistask.maxtrialtime * myscreen.framesPerSecond) + 100; % with some additional frames
        trialN      = thistask.numTrials;
        paramset    = thistask.parameter_set(thistask.experiment, thistask.experiment_paramset);

        noiseStruct(n) = generate_ar_sequence_cond(myscreen, trialN, T, ...
            paramset.stim_noiseStd, paramset.stim_noiseTau, paramset.stim_vel, ...
            paramset.point_noiseStd, paramset.point_noiseTau, paramset.point_vel);
    end
end


function noiseStruct = generate_ar_sequence_cond(myscreen, trialN, T, ...
    stim_noiseStd, stim_noiseTau, stim_vel, ...
    point_noiseStd, point_noiseTau, point_vel)

    % std:      in linear velocity in deg
    % noiseTau: time constant with respect to dt
    % vel:      constant linear velocity to be added
    
    % check parameters including trialN,T

    dt = 1/myscreen.framesPerSecond;

    noiseStruct = struct();
    
    noiseStruct.framesPerSecond = myscreen.framesPerSecond;
    noiseStruct.trialN = trialN;
    noiseStruct.T = T;
    
    noiseStruct.stim_noiseStd   = stim_noiseStd;
    noiseStruct.stim_noiseTau   = stim_noiseTau; % in tau
    noiseStruct.stim_vel        = stim_vel;
    
    noiseStruct.point_noiseStd  = point_noiseStd;
    noiseStruct.point_noiseTau  = point_noiseTau;
    noiseStruct.point_vel       = point_vel;
    
    noiseStruct.stim_noiseW     = cell(trialN,1);
    noiseStruct.stim_noiseAR    = cell(trialN,1);
    noiseStruct.stim_noise      = cell(trialN,1);
    
    noiseStruct.point_noiseW     = cell(trialN,1);
    noiseStruct.point_noiseAR    = cell(trialN,1);
    noiseStruct.point_noise      = cell(trialN,1);
    
    t_phi = 1 - dt/stim_noiseTau; 
    p_phi = 1 - dt/point_noiseTau; 
    
    tnoiseStd = stim_noiseStd  * sqrt(dt) * sqrt(prod(1-t_phi.^2));
    pnoiseStd = point_noiseStd * sqrt(dt) * sqrt(prod(1-p_phi.^2));
    
    for tr = 1:trialN
        [t_noise, t_wnoise]             = ar(T, tnoiseStd, t_phi, 'plotfigs', false);
        noiseStruct.stim_noiseW{tr}     = t_wnoise;
        noiseStruct.stim_noiseAR{tr}    = t_noise;
        noiseStruct.stim_noise{tr}      = t_noise + stim_vel;

        [p_noise, p_wnoise]             = ar(T, pnoiseStd, p_phi, 'plotfigs', false);
        noiseStruct.point_noiseW{tr}     = p_wnoise;
        noiseStruct.point_noiseAR{tr}    = p_noise;
        noiseStruct.point_noise{tr}      = p_noise + point_vel;
    end
end