%% Initialize blob, initialize blob
function blob = trackposInitStimulus(obj,myscreen, varargin)    

getArgs(varargin, {'reinit_img', 1});

% todo: make this only about generating blob image
    if isempty(obj)
        obj = struct();
    end
    
    blob = struct();

    if ~isfield(obj,'type') && ~isprop(obj,'type') 
        blob.type = 'gaussian';
    else
        blob.type = obj.type;
    end

    % blob size
    if ~isfield(obj,'std') && ~isprop(obj,'std')
        blob.std = 1;
    else
        blob.std = obj.std;
    end %unit: imageX, in deg.
    %GardnerLab: stimstd = 1; CSNL std = 0.4.. 
    blob.patchsize = min(6*blob.std,min(myscreen.imageWidth,myscreen.imageHeight));
    
    if ~isfield(obj,'position') && ~isprop(obj,'position') 
        %blob initial position. uniform distribution across the screen
        x_img = min(3*blob.std,1/3*myscreen.imageWidth)*(2*rand(1)-1); 
        y_img = min(3*blob.std,1/3*myscreen.imageWidth)*(2*rand(1)-1);
        blob.position = [x_img, y_img];
        blob.velocity = [0,0];
    else
       blob.position = obj.position; 
    end
    
    % blob speed
    % this might change based on effective sampling rate.
    if ~isfield(obj,'stepstd') && ~isprop(obj,'stepstd') 
        blob.stepstd = 1; % in deg/sec
    else
       blob.stepstd = obj.stepstd; 
    end %unit: cm/s to deg/frame
        
    % blob luminance
    if ~isfield(obj,'lum') && ~isprop(obj,'lum') 
        blob.lum = 122;
    else
       blob.lum = obj.lum; 
    end %unit: luminance
    
    if ~isfield(obj,'color') && ~isprop(obj,'color') 
       blob.color = [1;1;1];
    else
        if iscell(obj.color)
            obj.color = obj.color{1};
        end
        if isscalar(obj.color) && obj.color > double('a') && obj.color > double('z') 
            blob.color = char(obj.color);
        end

        if ischar(obj.color)
            if obj.color == 'k'
                blob.color = [1;1;1];
            elseif obj.color == 'r'
                blob.color = [1;0;0];
            elseif obj.color == 'b'
                blob.color = [0;0;1];
            elseif obj.color == 'g'
                blob.color = [0;1;0];
            elseif obj.color == '*'
                blob.color = '*';
            else
                print('check blob color property')
                blob.color = [1;1;1];
            end
        elseif isvector(obj.color)  && (length(obj.color) == 3 || length(obj.color) == 4)
            blob.color = obj.color;
        else
            print('check input color property')
            blob.color = [1;1;1];
        end
    end 
        
    % generate blob image
    if ~isfield(blob,'img') && reinit_img && isfield(blob, 'type') && strcmp(blob.type,'gaussian')
        if isfield(blob,'img')
            mglDeleteTexture(blob.img)
        end
        if blob.lum == 0 || blob.std == 0
            blob.img = [];
        else
            gaussian                = mglMakeGaussian(blob.patchsize,blob.patchsize,blob.std,blob.std)*(blob.lum);
            tableSize = mglPrivateSetGammaTable;
            if mglResolution().bitDepth ~= 32 % todo check bits for older systems.
                gaussian_rgb            = 255*repmat([blob.color; 1], 1, size(gaussian,2),size(gaussian,1));
                gaussian_rgb(4,:,:)     = round(gaussian'); 
                gaussian_rgb            = uint8(gaussian_rgb);
                blob.img                = mglCreateTexture(gaussian_rgb);
            elseif tableSize == 1024
                % e.g. argb2101010
                gaussian_rgb        = repmat([blob.color], 1, size(gaussian,2),size(gaussian,1)); % 3,x,y
                gaussian_rgb(4,:,:) = gaussian'; % 4,x,y
                gaussian_rgb        = permute(gaussian_rgb,[3,2,1]); % should be: y,x,4
                blob.img            = mglMetalCreateTexture(gaussian_rgb);
            else
                gaussian_rgb        = repmat([blob.color], 1, size(gaussian,2),size(gaussian,1)); % 3,x,y
                gaussian_rgb(4,:,:) = gaussian'; % 4,x,y
                gaussian_rgb        = permute(gaussian_rgb,[3,2,1]); % should be: y,x,4
                blob.img            = mglMetalCreateTexture(gaussian_rgb);
            end
        end
    end
end