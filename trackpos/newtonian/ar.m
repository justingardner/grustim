function [noise, white_noise] = ar(T, sigma_white, ar_coeff, varargin)
    getArgs(varargin, {'plotfigs', false});

    white_noise = normrnd(0, sigma_white, T, 1);
    noise = white_noise;
    for t = 1:T
        for ar_idx = 1:length(ar_coeff)
            if t-ar_idx > 0
                noise(t) = noise(t) + ar_coeff(ar_idx) * noise(t-ar_idx);
            end
        end
    end

    if plotfigs
        plot_ar(noise)
        % figure out variance
        disp(num2str(std(noise)));
        disp(num2str(sqrt(sigma_white^2./ prod(1-ar_coeff.^2))))
    end
end

function plot_ar(noise)
    shz = 60;
    figure;
    subplot(2,1,1)
    plot((1:length(noise))/shz, noise)
    title('Noise')
    
    subplot(2,1,2)
    [kernel, lags, n]= recoverFirstOrderKernel((1:length(noise))/shz,noise,noise,shz);
    plot(lags, kernel)
    title('Covariance')

end