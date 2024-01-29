function myscreen = setup_screen_jryu()
% set hostname: sudo scutil --set HostName NAME

    myscreen = struct();
    if isempty(mglGetSID)
        mglSetSID(-1);
    end

    myscreen.subjectID  = mglGetSID;

    [ret, name] = system('hostname');
    
    if contains(name, 'Joshuas-MacBook-Air.local') ||contains(name, 'Joshuas-MacBook-Air')
        myscreen.hideCursor         = 1;
        myscreen.displayName        = 'fulldisp';
        myscreen.saveData           = 1; % save stimfile to data directory
        myscreen.datadir            = '/Users/jryu/Dropbox/GardnerLab/data/';
        
    else
        w = warning ('off','all');
        % rmpath(genpath('/Users/gru/proj/mgl'))
        % addpath(genpath('/Users/gru/proj/mgl_jryu'))
        % addpath(genpath('/Users/gru/proj/mgl'))
        w = warning ('on','all');

        % myscreen.screenWidth = 860; myscreen.screenHeight = 600;
        myscreen.hideCursor         = 1;
        myscreen.displayName        = 'vpixx_close';
        myscreen.calibType          = 'Specify particular calibration';
        myscreen.calibFilename      = '/Users/gru/proj/mgl/task/displays/0002_psy-g7mhty2hjd_240124.mat';
        myscreen.calibFullFilename  = '/Users/gru/proj/mgl/task/displays/0002_psy-g7mhty2hjd_240124.mat';
        myscreen.saveData           = 1; % save stimfile to data directory
        myscreen.datadir            = '/Users/gru/data/';
    end

    warning('off','MATLAB:RandStream:ActivatingLegacyGenerators')
    warning('off','MATLAB:RandStream:ReadingInactiveLegacyGeneratorState') 

end