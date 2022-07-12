function [A, W, ns, ss] = NewtonianStateMatrix(varargin)
    p = inputParser;
    p.KeepUnmatched = true;
    
    p.addParameter('p0', 0, @(x)(isnumeric(x))) ; % constant input to position
    p.addParameter('v0', 0, @(x)(isnumeric(x))) ; % constant input to velocity 
    p.addParameter('a0', 0, @(x)(isnumeric(x))) ; % constant input to acceleration
    p.addParameter('p_std', 0, @(x)(isnumeric(x))); % position perturbations
    p.addParameter('v_std', 0, @(x)(isnumeric(x))); % velocity perturbations
    p.addParameter('a_std', 0, @(x)(isnumeric(x))); % acceleration perturbations

    p.parse(varargin{:})
    
    ss = 1; % state size
    if p.Results.a0 >0 || p.Results.a_std > 0
        ss = ss+2; 
    elseif p.Results.v0 >0 || p.Results.v_std > 0
        ss = ss+1; 
    end
    ns = ss;  % state size without constant inputs
    if p.Results.p0 >0, ss = ss+1;end
    if p.Results.v0 >0, ss = ss+1;end
    if p.Results.a0 >0, ss = ss+1;end
    
    Ac0 = eye(ns);
    Ac = Ac0(:,1:ss-ns); % constant inputs
    A = [triu(ones(ns)), Ac; zeros(ss-ns,ns), eye(ss-ns)];
    stds = [p.Results.p_std,p.Results.v_std,p.Results.a_std];
    W = [diag(stds(1:ns)); zeros(ss-ns,ns)];
end