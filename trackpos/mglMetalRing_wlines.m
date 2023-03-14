function mglMetalRing_wlines(r, lineWidth, color, res)
% with resolution 3; 0.0030 s

if isrow(color)
    color   = color(1:3)';
end

if size(color,2) == 1
    color = repmat(color(1:3), 1, res-1);
end
    
    thetas  = linspace(0,(2*pi),res);        
    x = r * cos(thetas);
    y = r * sin(thetas);
    mglMetalLines(x(1:end-1), y(1:end-1), x(2:end), y(2:end), ...
        repmat(lineWidth,1,res-1), color);
end

