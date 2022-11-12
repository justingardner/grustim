
function a = uniformrandnoise(myscreen, w,h,T)
    % x: width in degs
    % h: height in degs
    % T: time in seconds
    w_pixel = ceil(w*myscreen.screenWidth/myscreen.imageWidth);
    h_pixel = ceil(h*myscreen.screenHeight/myscreen.imageHeight);

    T_frame = ceil(T * myscreen.framesPerSecond);

    backgroundnoise_rgb = ones(4,w_pixel,h_pixel,T_frame);    
    backgroundnoise_rgb(4,:,:,:)  = rand(w_pixel,h_pixel,T_frame);

    savefile = '/Users/gru/proj/grustim/trackpos/noise/white0.mat';
    save(savefile, 'backgroundnoise_rgb','-v7.3')
end