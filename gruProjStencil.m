function gruProjStencil(myscreen)
%MGLPROJSTENCIL Generate a mask stencil for the GE projecter @ Stanford
%   Builds a stencil based on the input list, which corresponds to the
%   visual eccentricity at angles [90, -90, 0, 15, 30, 45, 60, 75, -15, -30, -45, -60, -75]
% use by adding mglStencilSelect(1); and then mglStencilSelect(0); to turn
% off

% This approximates the full screen (needs to be adjusted):

% botangs1 = [-15 -30];
% botangs2 = [-45 -60 -75];
% transEccs = [10./cos(deg2rad(abs(botangs1))) 10./cos(deg2rad(90+botangs2))];
% eccArray = [15,15,12.5,27,24,21,18,15, transEccs];
% eccArray = eccArray;

% This is the actual recordings we made for what you can see in the 32
% channel coil:

% eccArray = [5.5, 7.5, 17, 18, 19, 11, 8, 6,14 , 15, 11, 8, 7.5];

eccArray = [11, 14, 21, 17, 15, 13, 12, 11, 25, 23, 16,15, 14]*32/21;


angles = [90, -90, 0, 15, 30, 45, 60, 75, -15, -30, -45, -60, -75];
angles2 = [180-angles(3:8) angles(9:end)-90];
angles = [angles angles2];
angles = mod(angles-90,360);
eccArray = [eccArray eccArray(3:8) fliplr(eccArray(9:end))];

if length(eccArray) ~= length(angles)
    error('Eccentricity array has incorrect length.');
end

% we will draw partial disks at every 3 deg angle, and using the linear
% interpolated position based on the eccArray

xs = [];
ys = [];
eccs = [];
sAs = [];

for curBlock = 0:15:345
    for withinB = 0:0.5:14.5
        xs = [xs 0];
        ys = [ys 0];
        cSA = curBlock+withinB;
        sAs = [sAs cSA];
        orig_ang = eccArray(angles==curBlock);
        final_ang = eccArray(angles==mod(curBlock+15,360));
        interpEcc = (final_ang-orig_ang) * withinB / 15 + orig_ang;
        %         disp(sprintf('%i: %.03f',curBlock,interpEcc));
        eccs = [eccs interpEcc];
    end
end

mglGluPartialDisk(xs,ys,zeros(size(xs)),eccs,sAs,ones(size(xs))*2,ones(3,size(xs,2)));
mglGluPartialDisk(xs,ys,zeros(size(xs)),eccs,sAs,-ones(size(xs))*2,ones(3,size(xs,2)));



end

