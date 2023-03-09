function param = ornstein_uhlenbeck(maxorder, noisestd)
    
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

