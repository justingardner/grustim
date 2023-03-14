
function testDotsRing(pointer_r, ring_w, ecc_r, reinit_screen)

if reinit_screen == 1 
    setup_screen_jryu();
elseif reinit_screen == 2
    mglOpen;
    mglScreenCoordinates;
elseif reinit_screen == 3
    mglOpen(2);
    mglVisualAngleCoordinates(57,[52.1 29.5]);
end

mglClearScreen;

plotradius = false;
if plotradius
    mglMetalArcs([0;0;0], [0; 0; 1; 1], ...
    [ecc_r-ring_w; ecc_r+ring_w], [0;2*pi], 1);
else
    mglMetalArcs([0;0;0], [0; 0; 1; 1], ...
    [2*(ecc_r-ring_w); 2*(ecc_r+ring_w)], [0;2*pi], 1);
end

thetas = 0:0.1:(2*pi);        
zlt = zeros(size(thetas)); % zeros like thetas
mglMetalDots([ecc_r * cos(thetas); ecc_r * sin(thetas); zlt], ...
        [zlt;zlt;zlt+1;zlt+1], [zlt + pointer_r/2; zlt + pointer_r/2], [zlt+1], [zlt+1]);
        
% 
% for theta = 0:0.3:(2*pi)
%     mglMetalDots([ecc_r * cos(theta); ecc_r * sin(theta); 0], ...
%         [1;0;0;1], [pointer_r;pointer_r], 1, 1);
% end

mglFlush;

end