function mglMetalArrow(x0,y0,da,dr,arm_ratio, arm_angle, lineWidth, color)
    if isrow(color)
        color = color';
    end

    x = x0 + dr*cos(da);
    y = y0 + dr*sin(da);

    x_from  = [x0, x, x];
    y_from  = [y0, y, y];

    x_to    = [x, ...
        x + dr*arm_ratio * cos(da-pi+arm_angle), ...
        x + dr*arm_ratio * cos(da-pi-arm_angle)];
    
    y_to    = [y, ...
        y + dr*arm_ratio * sin(da-pi+arm_angle), ...
        y + dr*arm_ratio * sin(da-pi-arm_angle)];

    mglMetalLines(x_from,y_from, x_to, y_to, repmat(lineWidth,1,3), repmat(color(1:3)', 1, 3))
end

