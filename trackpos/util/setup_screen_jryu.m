function myscreen = setup_screen_jryu()
    myscreen = struct();
    if isempty(mglGetSID)
        mglSetSID(-1);
    end

    myscreen.subjectID  = mglGetSID;

    rmpath(genpath('/Users/gru/proj/mgl'))
    addpath(genpath('/Users/gru/proj/mgl_jryu'))

    % myscreen.screenWidth = 860; myscreen.screenHeight = 600;
    myscreen.hideCursor         = 1;
    myscreen.displayName        = 'vpixx';
    myscreen.calibType          = 'Specify particular calibration';
    myscreen.calibFilename      = '0001_dn0a221834_221005.mat';
    myscreen.calibFullFilename  = '/Users/gru/proj/mgl/task/displays/0001_dn0a221834_221005';
    myscreen.saveData           = 1; % save stimfile to data directory
    myscreen.datadir            = '/Users/gru/data/';
    myscreen                    = initScreen(myscreen);

    % set to argb2101010 pixel format
    mglMetalSetViewColorPixelFormat(4);
end