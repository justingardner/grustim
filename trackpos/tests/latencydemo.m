function latencydemo(varargin)

%     stim_highlum    = {'gaussian', 1, 1, 'k'}; 
%     stim_lowlum     = {'gaussian', 0.4, 1, 'k'}; 
%     stim_lowlum     = {'gaussian', 0.4, 1, 'r'}; 
%     stim_lowlum     = {'dot', 0.4, 1, '*'}; 
    objparams       = {...
        {'gaussian', 1, 1, 'k'}, ...
        {'gaussian', 0.2, 1, 'k'}, ...
        {'gaussian', 0.1, 1, 'k'}, ...
        {'gaussian', 0.05, 1, 'k'}};
% 
%     objparams       = {...
%         {'gaussian', 0.1, 1, 'g'}, ...
%         {'gaussian', 0.1, 1, 'r'}, ...
%         {'gaussian', 0.1, 1, 'k'}, ...
%         {'gaussian', 0.1, 1, 'b'}};

%     objparams       = {...
%         {'gaussian', 0.1, 1, 'g'}, ...
%         {'gaussian', 0.1, 1, 'r'}, ...
%         {'gaussian', 0.1, 1, 'b'}};

    movement_type = 'linear';

    vel     = 3; % angular velocity deg/s

    ecc_r   = 20;

    x_lim   = 5;

    ar          = 0.8;
    noiselevel  = 1;

    fixate              = false;
    fixation_size       = 0.4;
    randcolorsFile      = '/Users/jryu/proj/grustim/trackpos/util/labcolors.mat';
    randcolors          = load(randcolorsFile);

    N = length(objparams);
    objects         = cell(N,1);

    % init screen
    myscreen = setup_screen_jryu(); 
    myscreen = initScreen(myscreen);
    mglMetalSetViewColorPixelFormat(4);     % set to argb2101010 pixel format

    % 
    dt = 1/myscreen.framesPerSecond;

    for stimnum = 1:N
        stim = objparams{stimnum};
        [params.type, params.lum, params.std, params.color] = deal(stim{:});
        objects{stimnum,1}          = trackposInitStimulus(params, myscreen);

        if strcmp(movement_type,'linear')
            objects{stimnum,1}.y        = -ecc_r + stimnum/(N+1) * 2 * ecc_r ; 
            objects{stimnum,1}.position = [0, objects{stimnum,1}.y];
        elseif strcmp(movement_type, 'circular')
            objects{stimnum,1}.phase    = 2*pi*stimnum/N;
            objects{stimnum,1}.position = [ecc_r*cos(objects{stimnum,1}.phase), ecc_r*sin(objects{stimnum,1}.phase)];
        end
    end

    keystate = mglGetKeys;

    curr_phase      = 0;
    curr_cumnoise   = 0;

    while ~keystate(54) % esc
        keystate = mglGetKeys;
        mglClearScreen(0.7);

        % fixation
        if fixate
            mglMetalArcs([0;0;0], [1;1;1; 1], [fixation_size+0.1;fixation_size+0.3],[0;2*pi], 1);
            fixcolors = randcolors.colors(randi(size(randcolors.colors,1)),:)';
            mglMetalDots([0;0;0], [fixcolors; 1], [fixation_size;fixation_size], 1, 1);
        end

        % objects
        curr_cumnoise   = ar*curr_cumnoise + noiselevel*randn() * sqrt(dt);        
        curr_phase      = curr_phase + curr_cumnoise + vel*dt;
        if curr_phase > x_lim
            vel = -1*abs(vel); 
        elseif curr_phase < -1* x_lim
            vel = abs(vel);
        end

        for stimnum = 1:N
            if strcmp(movement_type,'linear')
                objects{stimnum,1}.position = [curr_phase objects{stimnum,1}.y];
            elseif strcmp(movement_type, 'circular')
                objects{stimnum,1}.phase = objects{stimnum,1}.phase + vel*dt;
                objects{stimnum,1}.position = [ecc_r*cos(objects{stimnum,1}.phase), ecc_r*sin(objects{stimnum,1}.phase)];
            end

            obj = objects{stimnum,1};
            
            if isfield(obj, 'img') && ~isempty(obj.img)
                mglBltTexture(obj.img, obj.position);
            else
                if strcmp(obj.color,'*') 
                    colors = randcolors.colors(randi(size(randcolors.colors,1)),:)';
                else
                    colors = obj.color;
                end
                mglMetalDots([obj.position(1);obj.position(2);0],[colors;1], [obj.std; obj.std], 1, 1);
            end
        end

        mglFlush();
    end
end

