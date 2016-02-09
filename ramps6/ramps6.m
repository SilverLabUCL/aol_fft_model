function ramps = ramps6(z, v)

vel_angle = angle(v(1) + 1i * v(2));
if is_0_to_60(vel_angle)
    theta = vel_angle;
    a = @ramp0;
    b = @ramp60;
elseif is_60_to_120(vel_angle)
    theta = vel_angle - pi/3;
    a = @ramp60;
    b = @ramp120;
else 
    theta = vel_angle - 2*pi/3;
    a = @ramp120;
    b = @(z,v) ramp0(z, -v);
end

const = norm(v) / sin(2*pi/3); % use sine rule
v_a = const * sin(theta); 
v_b = const * sin(pi/3 - theta);
tot = sin(theta) + sin(pi/3 - theta);
ramps = a(z * tot/sin(theta), v_a) + b(z * tot/sin(pi/3 - theta), v_b);
end

function bool = is_0_to_60(vel_angle)
    pos = vel_angle > 0 && vel_angle <= pi/3;
    neg = vel_angle > -pi && vel_angle <= - 2*pi/3;
    bool = pos || neg;
end
function bool = is_60_to_120(vel_angle)
    pos = vel_angle > pi/3 && vel_angle <= 2*pi/3;
    neg = vel_angle > -2*pi/3 && vel_angle <= -pi/3;
    bool = pos || neg;
end
