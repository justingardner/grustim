function state = ou_simulate_full(dynparams, T, dt, varargin)
    % state: is added to the last order
    
    getArgs(varargin, {'state', zeros(T,1), 'b_input', zeros(T,dynparams.maxorder+1)});
    for revord = 0:dynparams.maxorder
        ord = dynparams.maxorder - revord;
        noisestd        = dynparams.(['thetaStd', num2str(ord)]);
        invtau_decay    = dynparams.(['invtaus_decay', num2str(ord)]);
        tau_int         = dynparams.(['taus_int', num2str(ord)]);

        input               = state + noisestd * randn(T,1);
        integrator          = cumprod([dt/tau_int; (1 - invtau_decay*dt) * ones(T-1,1)]);
        integrated_input    = conv(input,integrator, "full");
        state               = integrated_input(1:T) + b_input(:,ord+1);
    end
    
    % add random offset
    % state = state + 2*pi*rand();
end