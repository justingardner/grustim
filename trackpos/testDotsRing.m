
function testDotsRing(pointer_r, ring_r, ecc_r, reinit_screen)

if reinit_screen == 1 
    setup_screen();
elseif reinit_screen == 2
    mglOpen;
    mglScreenCoordinates;
end

%% plot stuff

mglClearScreen;
% draw stuff
mglMetalArcs([0;0;0], [0; 0; 1; 1], ...
    [ecc_r-0.1; ecc_r+0.1], [0;2*pi], 1);

% mglMetalArcs([0;0;0], [1;1;1; 1], [pointer_r+0.1;pointer_r+0.3],[0;2*pi], 1);
% mglMetalDots([0;0;0], [0.5+0.5*rand(3,1);1], [pointer_r;pointer_r], 1, 1);

for theta = 0:1:(2*pi)
%     mglGluDisk(ecc_r * cos(theta), ecc_r * sin(theta), ...
%         pointer_r, [1,1,0],60,1);
    mglMetalDots([ecc_r * cos(theta); ecc_r * sin(theta); 0], ...
        [1;0;0;1], [pointer_r;pointer_r], 1, 1);
end

mglFlush;

end

function setup_screen()
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