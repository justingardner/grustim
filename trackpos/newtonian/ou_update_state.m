function [state_t, noise_t] = ou_update_state(state_t, input_t, dynparams, dt, varargin)
% N = (dynparams.maxorder + 1) defines the dimensions of the state
% dyn.params.maxorder = 0 =>  position representation only. 
%                             velocity (order 1) is the input to the position
% dyn.params.maxorder = 1 =>  position, velocity state representation. 
%                             velocity (order 1) and acceleration (order 2)
%                             are inputs

% state dimension >= N. All dimensions greater than N are ignored, not updated.
% input dimension >= N
    % if input dimension  = 1, the input is applied to the last dimension

    getArgs(varargin, {'noise_t', randn(dynparams.maxorder+1,1)});
    
    if numel(input_t) == 1 && numel(input_t) < (dynparams.maxorder+1)
        input_t = [zeros(dynparams.maxorder,1); input_t];
    end

    % state_t, input_t: array of order x 1
    state_order = 0;
    for revord = 0:dynparams.maxorder
        ord = dynparams.maxorder - revord;
        noisestd        = dynparams.(['thetaStd', num2str(ord)]);
        invtau_decay    = dynparams.(['invtaus_decay', num2str(ord)]);
        tau_int         = dynparams.(['taus_int', num2str(ord)]);

        input_all_t         = state_order + input_t(ord+1) + noisestd * noise_t(ord+1);
        state_order         = (1-invtau_decay*dt)* state_t(ord+1) + dt/tau_int * input_all_t;
        state_t(ord+1)      = state_order;
    end
end
    

