function mglMetalCircArrow(r0, pa, arrow_length, arm_ratio, arm_angle, lineWidth, color, varargin)
    getArgs(varargin,'res=10');
    
    % draws circular arrow
    if isrow(color)
        color = color';
    end

    da_vec = linspace(0,arrow_length/r0,res);
    x = r0 * cos(pa + da_vec);
    y = r0 * sin(pa + da_vec);

    x_from  = [x, x(end)];
    y_from  = [y, y(end)];

    x_to    = [x(2:end), ...
        x(end) + arrow_length*arm_ratio * cos(pa+da_vec(end)-pi/2+arm_angle), ...
        x(end) + arrow_length*arm_ratio * cos(pa+da_vec(end)-pi/2-arm_angle)];
    
    y_to    = [y(2:end), ...
        y(end) + arrow_length*arm_ratio * sin(pa+da_vec(end)-pi/2+arm_angle), ...
        y(end) + arrow_length*arm_ratio * sin(pa+da_vec(end)-pi/2-arm_angle)];

    mglMetalLines(x_from, y_from, x_to, y_to, ...
        repmat(lineWidth,1,numel(x_from)), ...
        repmat(color(1:3), 1, numel(x_from)));
end