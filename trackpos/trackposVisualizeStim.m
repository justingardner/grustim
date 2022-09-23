function rmsContrast = trackposVisualizeStim()
    
% myscreen.calibType = 'None';
% myscreen.calibFilename = [];
% myscreen.calibFullFilename = [];

    stimLums    = linspace(0.05,1,10); %[0.05,0.1,0.2,0.3];
    backLums    = 0 * ones(length(stimLums),1);
    stimStds    = 1 * ones(length(stimLums),1);
    
    % myscreen = initScreen('gardnerlab420');
    myscreen = initScreen('gardnerlab420_uncalib');
    
    if isfield(myscreen, 'calibFullFilename') && ~isempty(myscreen.calibFullFilename)
        a = load(myscreen.calibFullFilename);
    else
        a = load('/Users/JRyu/github/mgl/task/displays/0001_dn0a22167c_220912.mat');
    end
    calib = a.calib;
    mglMetalSetViewColorPixelFormat(4);
    
    xpos = -myscreen.imageWidth/3;
    ypos = myscreen.imageHeight/3 - max(stimStds)*2;
    
    mglClearScreen(backLums(1));
    
    for cidx = 1:4
        if cidx == 1
            currcolor = 'k';
        elseif cidx == 2
            currcolor = 'r';
        elseif cidx == 3
            currcolor = 'g';
        elseif cidx == 4
            currcolor = 'b';
        end
    for i = 1:length(stimLums)
        obj         = struct();
        obj.lum     = stimLums(i);
        obj.std     = stimStds(i);
        obj.color   = currcolor;
        blob = trackposInitStimulus(obj,myscreen);  
        
        gaussian        = mglMakeGaussian(blob.patchsize,blob.patchsize,blob.std,blob.std);
        alpha           = backLums(i) + stimLums(i).*gaussian.*(1-backLums(i));
        gaussian_lums   = interp1(calib.uncorrected.outputValues,calib.uncorrected.luminance,alpha(:),'linear');
        backLum_lum     = interp1(calib.uncorrected.outputValues,calib.uncorrected.luminance,backLums(i),'linear');
        rmsContrast(i)  = rms(gaussian_lums - backLum_lum);
        maximumLuminance(i) = max(gaussian_lums(:));
        
        xpos = xpos + obj.std*2; 
        if xpos > myscreen.imageWidth/2
            ypos = ypos - max(stimStds)*6; 
            xpos = -myscreen.imageWidth/3 + obj.std*2; 
        end
        mglBltTexture(blob.img,[xpos,ypos]);
        xpos = xpos + obj.std*2; 
    end
    ypos = ypos - max(stimStds)*6; 
    xpos = -myscreen.imageWidth/3; 
    end

    disp('rms contrasts: ')
    disp(rmsContrast);
    disp('maximum luminance (cd/m^2): ')
    disp(maximumLuminance);
    mglFlush; 
    
    dbstep
    
    if false
        % check alpha compositing
        ypos = ypos - max(stimStds)*6; 
        mglClearScreen(0);            
        for i = 1:length(stimLums)
            alpha       = backLums(i) + stimLums(i).*gaussian.*(1-backLums(i));
            obj         = struct();
            obj.lum     = max(alpha(:));
            obj.std     = stimStds(i);
            blob = trackposInitStimulus(obj,myscreen);  

            xpos = xpos + blob.std*2; 
            if xpos > myscreen.imageWidth/2
                ypos = ypos - max(stimStds)*6; 
                xpos = -myscreen.imageWidth/3 + blob.std*2; 
            end
            mglBltTexture(blob.img,[xpos,ypos]);
            xpos = xpos + blob.std*2; 
        end
        mglFlush;  
    end
    
    mglClose;
    
end

function calibration

calib = moncalib('numRepeats=3','stepsize=1/127', 'initWaitTime=60','screenNumber=2',...
                 'spectrum=1', 'gamma=1', 'gammaEachChannel=1',...
                 'tableTest=1',...        
                 'bitTest=1','bitTestType=2','bitTestScreenMode=4','bitTestNumRepeats=3');

end

function colorTest()

myscreen = initScreen();
mglMetalSetViewColorPixelFormat(4);

pointer             = struct();
pointer.std = 0.1; 
pointer.lum = 1;
pointer.color = 'r';
stimulus.pointer            = trackposInitStimulus(pointer,myscreen);

mglClearScreen(0.4);

mglGluDisk(-3,0,0.2,[0,0,1],60,1);      
mglGluDisk(-2,0,0.2,[0,1,0],60,1);      
mglGluDisk(-1,0,0.2,[1,0,0],60,1);      
mglBltTexture(stimulus.pointer.img, [1,0]);

mglFlush;

end

function checkPrivateGammaTable

myscreen = initScreen('gardnerlab420');
% myscreen = initScreen('gardnerlab420_uncalib');

currcalib= load(calibFilename);
myscreen.gammaTable = currcalib.calib.table;
myscreen.initScreenGammaTable;

% compare gammaTable with initScreenGammabTable
myscreen.gammaTable
% myscreen.initScreenGammaTable % from mglGetGammaTable
tt = mglPrivateGetGammaTable;

figure;
subplot(2,1,1);hold on;
plot(tt.redTable);
plot(tt.greenTable);
plot(tt.blueTable);
plot(myscreen.gammaTable);
title('mglPrviateGetGammatable')
legend('r','g','b','inputTable')

subplot(2,1,2);hold on;
plot(myscreen.initScreenGammaTable.redTable);
plot(myscreen.initScreenGammaTable.greenTable);
plot(myscreen.initScreenGammaTable.blueTable);
plot(myscreen.gammaTable);
title('myscreen.initScreenGammaTable')
legend('r','g','b','inputTable')

myscreen = initScreen('gardnerlab420_uncalib');
tt1 = mglPrivateGetGammaTable;
mglClose;
myscreen = initScreen('gardnerlab420');
tt = mglPrivateGetGammaTable;
mglClose;

figure;hold on;
plot(tt.redTable,'r','LineWidth',2);
plot(tt1.redTable,'r:','LineWidth',3);
plot(tt.greenTable,'g','LineWidth',2);
plot(tt1.greenTable,'g:','LineWidth',3);
plot(tt.blueTable,'b','LineWidth',2);
plot(tt1.blueTable,'b:','LineWidth',3);
plot(myscreen.gammaTable,'k:','LineWidth',2);
title('mglPrviateGetGammatable')
legend('r-calib','r-nocalib','g-calib','g-nocalib','b-calib','b-nocalib','calib inputTable')

end

function linearizeColorCalib

figure; hold on; 
plot(a.calib.uncorrectedEachChannel.R.luminance/max(a.calib.uncorrectedEachChannel.R.luminance),'r'); 
plot(a.calib.uncorrectedEachChannel.G.luminance/max(a.calib.uncorrectedEachChannel.G.luminance),'g'); 
plot(a.calib.uncorrectedEachChannel.B.luminance/max(a.calib.uncorrectedEachChannel.B.luminance),'b');
plot(a.calib.table,'k');

a = load('/Users/JRyu/github/mgl/task/displays/0001_dn0a22167c_220912.mat');
tableSize = 1024; %mglPrivateSetGammaTable;
desiredOutput_r = min(a.calib.uncorrectedEachChannel.R.luminance):(max(a.calib.uncorrectedEachChannel.R.luminance)-min(a.calib.uncorrectedEachChannel.R.luminance))/(tableSize-1):max(a.calib.uncorrectedEachChannel.R.luminance);
desiredOutput_g = min(a.calib.uncorrectedEachChannel.G.luminance):(max(a.calib.uncorrectedEachChannel.G.luminance)-min(a.calib.uncorrectedEachChannel.G.luminance))/(tableSize-1):max(a.calib.uncorrectedEachChannel.G.luminance);
desiredOutput_b = min(a.calib.uncorrectedEachChannel.B.luminance):(max(a.calib.uncorrectedEachChannel.B.luminance)-min(a.calib.uncorrectedEachChannel.B.luminance))/(tableSize-1):max(a.calib.uncorrectedEachChannel.B.luminance);
t = struct();
outputValues_r = [];outputValues_g = []; outputValues_b=[];
for idx = 1:length(a.calib.uncorrectedEachChannel.R.outputValues) 
outputValues_r(idx) = a.calib.uncorrectedEachChannel.R.outputValues{idx}(1);
outputValues_g(idx) = a.calib.uncorrectedEachChannel.G.outputValues{idx}(2);
outputValues_b(idx) = a.calib.uncorrectedEachChannel.B.outputValues{idx}(3);
end
% fit exponent here. the measured values are not unique. 
redfit          = fitExponent(outputValues_r,a.calib.uncorrectedEachChannel.R.luminance,1);
greenfit        = fitExponent(outputValues_g(50:end),a.calib.uncorrectedEachChannel.G.luminance(50:end),1); % change to log loss to fit the small parts
bluefit         = fitExponent(outputValues_b,a.calib.uncorrectedEachChannel.B.luminance,1);

t.redTable      = abs(interp1(redfit.fit,outputValues_r,desiredOutput_r,'linear'));
t.greenTable    = abs(interp1(real(greenfit.fit),outputValues_g(50:end),desiredOutput_g,'linear'));
t.blueTable     = abs(interp1(bluefit.fit,outputValues_b,desiredOutput_b,'linear'));

t.redMin = 0;
t.redMax = 1;
t.redGamma = -1/redfit.tau;
t.greenMin = 0;
t.greenMax = 1;
t.greenGamma = 1/abs(greenfit.tau);
t.blueMin= 0;
t.blueMax= 1;
t.blueGamma= -1/bluefit.tau;

calib = a.calib;
calib.table = t;
save('/Users/JRyu/github/mgl/task/displays/0001_dn0a22167c_220912_rgb.mat','calib');


b = load('/Users/JRyu/github/mgl/task/displays/0003_dn0a2213b5_220901.mat') ;
c = load('/Users/JRyu/github/mgl/task/displays/0002_dn0a2213b5_220901.mat') ;
d = load('/Users/JRyu/github/mgl/task/displays/0001_dn0a2213b5_220901.mat') ;

end
