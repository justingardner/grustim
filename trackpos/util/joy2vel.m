function [vx, vy] = joy2vel(joy, prms, myscreen)
    % convert the axes of the joystick to velocity
    % prms: struct with transducer parameters
        % prms.maxv: maximum velocity as a fraction of screen radius
        % prms.deadzone: width of deadzone in fraction. 
        % prms.sensivity: order of polynomial (positive real) s = s^{2*sens + 1}
        % prms.poly_order: order of polynomial for transducer (positive real) 
    % todo: try rotation in z;
    % todo: turn throttle off
    % todo: use buttons to control sensitivity and transducer exponents
        
    dt      = 1/myscreen.framesPerSecond;
    maxr    = prms.maxv * min(myscreen.imageWidth/2,myscreen.imageHeight/2);
    maxv    = maxr; % maxmimum velocity takes our pointer off the screen within one time frame
    
    a = axis(joy); 
    a = a - [prms.bias 0];
    % make forward go up; make sensitivity increase as throttle is more
    % forward
    theta = -1*a(3)*prms.theta_range;  % -1 to 1 =>  -pi/6 to pi/6 % 
    x0 = a(1)* prms.xyratio; y0=-1* a(2); 
    r0 = sqrt(x0^2 + y0^2);
    theta0 = atan2(y0,x0);
    x = r0*cos(theta0+theta);
    y = r0*sin(theta0+theta);
    s = -1* a(4); 

    sigma = polytrans(convert_range(s),prms.sensitivity);
    vx = lintrans(...
            polytrans(...
                deadzone(x,prms.deadzone),...
                prms.poly_order),...
            sigma*maxv);
    vy = lintrans(...
            polytrans(...
                deadzone(y,prms.deadzone),...
                prms.poly_order),...
            sigma*maxv);
end

function y = convert_range(x)
    % convert range from -1 to 1
    % to 0 to 1;
    y = (x+1)/2;

end

function y = deadzone(x, deadzone_frac)
% truncate 
    idx = (x < deadzone_frac) & (x > -1 * deadzone_frac);
    y = x;
    y(idx) = 0;
    % rescale everything else to new range
    y(y>0) = 1/(1-deadzone_frac)*(y(y>0) - deadzone_frac); 
    y(y<0) = 1/(1-deadzone_frac)*(deadzone_frac + y(y<0)); 
end

function y = lintrans(x,ymax)
% truncated linear 
% x is from -1 to 1. 
% y is from -max to max
    y = ymax * x;
end

function y = polytrans(x,order)    
    y = x;
    y(x>=0) = 1 * abs(x(x>=0)).^(order);
    y(x<0) = -1 * abs(x(x<0)).^(order);
end
