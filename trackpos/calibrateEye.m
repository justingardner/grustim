function [positions_target, positions_eye] = calibrateEye(myscreen, stimulus, recalibrateEyelink)
    % todo: adaptive calibration
    % check which part of screen has big error and finely sample the area.
    % todo: option to skip calibration point


    % calibration parameters
    r = [0.5, 3, 6, 9, 12];
    a = 0:pi/6:2*pi;
    steady_time_thresh = 30; % in frames. have to fixate for this long.
    steady_deg_thresh = 0.1; % in degrees. threshold for eye being fixated
    
    disp(' (calibrateEye) calibrating eye tracker ...');
    
    % run eyelink calibration first
    if recalibrateEyelink
        myscreen = eyeCalibDisp(myscreen); 
    end

    [rlist, alist] = ndgrid(r,a);
    positions_target = [rlist(:), alist(:)];
    idx = randperm(size(positions_target,1));
    positions_target = positions_target(idx,:);
    positions_target = [0,0; positions_target]; 
    positions_target = [positions_target(:,1) .* cos(positions_target(:,2)),...
                        positions_target(:,1) .* sin(positions_target(:,2))];
    positions_eye = nan(size(positions_target));
    n = 1; 
    
    [pos, postime] = mglEyelinkGetCurrentEyePos; % is this in image coordinates?
    currpos = pos;
    steady= 0;
    
    mglClearScreen;
    
    while ~myscreen.userHitEsc
        mglClearScreen;
        
        if n == 1
            mglTextDraw('(calibrateEye) Calibrating eye tracker...',[0 -5]);
            mglTextDraw('Please fixate on the white dot (inside the red dot).',[0 -7]);
            mglTextDraw('Press esc to skip calibration',[0 -9]);
        end
        
        if steady > steady_time_thresh
            positions_eye(n,:) = pos; % record final position
            steady = 0;
            error = hypot(pos(1)-positions_target(n,1), pos(2)-positions_target(n,2));
            sprintf('target position (%2.2f, %2.2f); eye position: (%2.2f, %2.2f); error: %2.2f',...
                positions_target(n,1), positions_target(n,2), pos(1), pos(2), error);
            n = n+1; % go to next stimulus.
            if n > size(positions_target,1)
                break % done.
            end
        end
        
        [pos, postime] = mglEyelinkGetCurrentEyePos; % is this in image coordinates?
        mglGluDisk(pos(1), pos(2), 0.2, [0 0 1, 0.2]) % faint blue dot to track eye
        mglGluDisk(positions_target(n,1), positions_target(n,2), 0.2, [1 0 0, 1])
        mglGluDisk(positions_target(n,1), positions_target(n,2), 0.1, [1 1 1, 1])
        % disp(pos);
        
        if hypot(pos(1)-currpos(1),pos(2)-currpos(2)) < steady_deg_thresh
            steady = steady+1;
        else
            steady = 0;
            currpos = pos;
            % disp('not steady -- reseting center') 
        end
        
        myscreen = tickScreen(myscreen,[]);     % flip screen
    end

end