function myscreen = setup_screen_jryu()
    myscreen = struct();
    if isempty(mglGetSID)
        mglSetSID(-1);
    end

    myscreen.subjectID  = mglGetSID;

    [ret, name] = system('hostname');
    
    if contains(name, 'Joshuas-MacBook-Air.local') ||contains(name, 'DN2t9suv.SUNet')  || contains(name, 'DNa80fac0.SUNet') 
        myscreen.hideCursor         = 1;
        myscreen.displayName        = 'fulldisp';
        myscreen.saveData           = 1; % save stimfile to data directory
        myscreen.datadir            = '/Users/jryu/Dropbox/GardnerLab/data/';
        
    else
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
    end

end