function bestfit = fitNormCDF(testVals,pcorr, varargin)
    fittype = 'std';
    if strcmp(fittype, 'std')
        minParams = [0];
        maxParams = [inf];
        initParams = [1];
    end

    bestfit = fitCumulativeGaussian2(testVals,pcorr,...
        'minParams', minParams, 'maxParams', maxParams, 'initParams', initParams,...
        'fitType',fittype,varargin{:});
end