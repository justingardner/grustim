function test()
    T = 200;
    noisestd = 1;
    maxorder = 0;
    dt = 1/30;
    dynparams = generate_param(maxorder, noisestd);

    tic
    rng(0)
    state = initialize_full_trajectory(dynparams, T, dt);
    toc

    tic
    rng(0)
    state_t = zeros(maxorder+1,1);
    state1 = zeros(T,1);
    for t=1:T 
        state_t = update_state(state_t, ones(maxorder+1,1), dynparams, dt);
        state1(t) = state_t(1);
    end
    toc

    figure; hold on;
    plot(state); plot(state1);



end

function param = generate_param(maxorder, noisestd)
    param = struct();
    param.('thetaStd0')      = 0;
    param.('invtaus_decay0') = 0;
    param.('taus_int0')      = 1;

    param.maxorder = maxorder;

    for ord = 1:param.maxorder % velocity and acceleration has decays
        param.(['thetaStd', num2str(ord)])      = 0;
        param.(['invtaus_decay', num2str(ord)]) = 0.1;
        param.(['taus_int', num2str(ord)])      = 1;
    end
        
    param.(['thetaStd', num2str(param.maxorder)])   = noisestd;

end


function state = initialize_full_trajectory(dynparams, T, dt)
    state = zeros(T,1);
    for revord = 0:dynparams.maxorder
        ord = dynparams.maxorder - revord;
        noisestd        = dynparams.(['thetaStd', num2str(ord)]);
        invtau_decay    = dynparams.(['invtaus_decay', num2str(ord)]);
        tau_int         = dynparams.(['taus_int', num2str(ord)]);

        input               = state + noisestd * randn(T,1);
        integrator          = cumprod([dt/tau_int; (1 - invtau_decay*dt) * ones(T-1,1)]);
        integrated_noise    = conv(input,integrator, "full");
        state               = integrated_noise(1:T);
    end
end

function state_t = update_state(state_t, input_t, dynparams, dt)
    % state_t, input_t: array of order x 1
    state_order = 0;
    for revord = 0:dynparams.maxorder
        ord = dynparams.maxorder - revord;
        noisestd        = dynparams.(['thetaStd', num2str(ord)]);
        invtau_decay    = dynparams.(['invtaus_decay', num2str(ord)]);
        tau_int         = dynparams.(['taus_int', num2str(ord)]);

        input_all_t         = state_order + input_t(ord+1) + noisestd * randn();
        state_order         = (1-invtau_decay*dt)* state_t(ord+1) + dt/tau_int * input_all_t;
        state_t(ord+1)      = state_order;
    end
end
