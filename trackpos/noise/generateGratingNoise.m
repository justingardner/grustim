function generateGratingNoise()
    myscreen = initScreen();
    stimulus = struct();
    stimulus.lum        = 1;
    stimulus.std        = 0.1;    
    stimulus.color      = 'k';    
    stimulus.backLum    = 0.4;
    stimulus = trackposInitStimulus(stimulus,myscreen);
    
    w = 20;% w: width in degs
    h = w; % h: height in degs
    T = 1; % T: time in seconds
    savefile = '/Users/gru/proj/grustim/trackpos/noise/grating.mat';

    gratingNoise(myscreen, stimulus, w, h,T, savefile);
    mglClose;
end


function gratingNoise(myscreen, stimulus, w, h,T, savefile)
    % w: width in degs
    % h: height in degs
    % T: time in seconds
    T_frame = ceil(T * myscreen.framesPerSecond);
    [h_pixel, w_pixel] = size(mglMakeGrating(w,h,0,0,0));
    
    backgroundnoise_rgb  = ones(4,h_pixel,w_pixel,T_frame); %0.1165 s
    for idx2 = 1:T_frame
        angle  = rand(1)*180;
        sf     = exp(normrnd(1,1)); % cycles per degree
        grating = mglMakeGrating(w,h,sf,angle,0);
        backgroundnoise_rgb(4,:,:,idx2) = grating'/max(grating(:));  % normalize contrast %0.025s
    end
    
    save(savefile, 'backgroundnoise_rgb','-v7.3')
end


function phaseScrambleStim_old(myscreen,stimulus)
    % generate a lot of noise
    % takes about 12 minutes to generate on the stimulus test computer. 
    savefile = '/Users/gru/data/trackpos/trackpos.mat';
    
    %% noise
    % background noise
    downsample_timeRes  = 1;
    downsample_spatRes  = 10;
    
    %ntrials             = 1; %7*9; %7 trials * 9 conditions
    nframes             = 10000; %myscreen.framesPerSecond*30; % 30s; %/downsample_timeRes; 
    
    %for idx = 1:ntrials
    tic
    if ~isfield(stimulus,'gaussianFFT')
        xsize_deg          = round(myscreen.imageWidth/downsample_spatRes);
        ysize_deg          = round(myscreen.imageHeight/downsample_spatRes);
    
        backgaussian = mglMakeGaussian(xsize_deg,ysize_deg,...
            stimulus.stimStd/downsample_spatRes,stimulus.stimStd/downsample_spatRes)*255;
        stimulus.gaussianFFT = getHalfFourier(backgaussian);
    end 
    
    for idx2 = 1:nframes
        back                        = stimulus.gaussianFFT; %0.02s
        back.phase                  = rand(size(back.mag))*2*pi; % scramble phase % 0.02s
        backgroundnoise             = round(reconstructFromHalfFourier(back));   %0.04s
        if idx2 == 1 % to save time; only allocate memory once.
            backgroundnoise_rgb         = 255*ones(4,size(backgroundnoise,2),size(backgroundnoise,1),nframes,'uint8'); %0.1165 s
        end % otherwise just overwrite.
        backgroundnoise_rgb(4,:,:,idx2)  = backgroundnoise'/max(max(backgroundnoise))*stimulus.noiseLum;  % normalize contrast %0.025s
    end
    backgroundnoise_rgb = uint8(backgroundnoise_rgb);
    %eval(['backgroundnoise_rgb =uint8(backgroundnoise_rgb);']); %0.02s 
    toc
    %end
    save(savefile, 'backgroundnoise_rgb*','-v7.3')

end