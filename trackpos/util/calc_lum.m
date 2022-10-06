function [maxStimLum, maxBackLum, maxAlpha] = calc_lum(stimLum, backLum, calib)
    alpha       = backLum + stimLum .* (1-backLum);
    maxLum      = interp1(calib.uncorrected.outputValues,calib.uncorrected.luminance,[stimLum, backLum, alpha],'linear');
    maxStimLum  = maxLum(1);
    maxBackLum  = maxLum(2);
    maxAlpha    = maxLum(3);
end

