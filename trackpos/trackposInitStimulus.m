%% Initialize stimulus, initialize blob
function stimulus = trackposInitStimulus(obj,myscreen)      
    
    stimulus = obj;
    % stimulus size
    if ~isfield(stimulus,'stimStd'), stimulus.stimStd = 1;,end %unit: imageX, in deg. 
    %GardnerLab: stimstd = 1; CSNL stimStd = 0.4.. (why...?)
    stimulus.patchsize = min(6*stimulus.stimStd,min(myscreen.imageWidth,myscreen.imageHeight));
    
    if ~isfield(stimulus,'position'), 
        %stimulus initial position. uniform distribution across the screen
        x_img = min(3*stimulus.stimStd,1/3*myscreen.imageWidth)*(2*rand(1)-1); 
        y_img = min(3*stimulus.stimStd,1/3*myscreen.imageWidth)*(2*rand(1)-1);
        stimulus.position = [x_img, y_img];
        stimulus.velocity = [0,0];
    end
    
    % stimulus speed
    % this might change based on effective sampling rate.
    if ~isfield(stimulus,'stepStd'), stimulus.stepStd = 3/myscreen.framesPerSecond;,end %unit: cm/s to deg/frame
        
    % stimulus luminance
    if ~isfield(stimulus,'stimLum'), stimulus.stimLum = 122;,end %unit: luminance
            
    % background noise
    if ~isfield(stimulus,'noiseLum'), stimulus.noiseLum = 122;,end; % unit: luminance
    
    % background luminance
    if ~isfield(stimulus,'backLum'), stimulus.backLum = 32;,end; % unit: luminance

    % generate stimulus image
    if isfield(stimulus,'gaussian'), mglDeleteTexture(stimulus.gaussian);, end 
    gaussian    =  mglMakeGaussian(stimulus.patchsize,stimulus.patchsize,...
        stimulus.stimStd,stimulus.stimStd)*(stimulus.stimLum);
    gaussian_rgb           = 255*ones(4,size(gaussian,2),size(gaussian,1),'uint8');
    gaussian_rgb(4,:,:)    = round(gaussian');
    gaussian_rgb           = uint8(gaussian_rgb);
    stimulus.gaussian = mglCreateTexture(gaussian_rgb);
    
    %% todo: move to main task?
    % pointer position
    stimulus.pointer            = [0, 0]; % pointer position 
    
    % fixation cross colors
    stimulus.fixColors.response = [1 1 1];
    stimulus.fixColors.stim     = [0 1 0]; % green
    stimulus.fixColors.est      = [1 0 0]; % red
    stimulus.fixColors.afc      = [0 0 1]; % blue
end