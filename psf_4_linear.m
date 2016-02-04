function res = psf_4_linear(aol, time, xy, z, v, w3, w4, w5, ws, wf, plot)
    z_x = z + aol.spacing;
    z_y = z;
    L_x = 2 * aol.spacing;
    L_y = 2 * aol.spacing;   
    wavelen = (2*pi) ./ aol.k;
    pdr = 0;
    
    df = - aol.V ./ wavelen .* xy(:) ./ (pdr .* [(L_x + z_x); (L_y + z_y)] + [z_x; z_y]);
    a = 0e6 + [pdr .* df(1,:); pdr .* df(2,:); - df(1,:); - df(2,:)];
    
    b = - aol.V.^2 ./ wavelen .* [...
        (1 + v(1)) ./ (L_x .* (1 + v(1)) + 2 * z_x);...
        (1 + v(2)) ./ (L_y .* (1 + v(2)) + 2 * z_y);... 
        (1 - v(1)) ./ (2 * z_x);...
        (1 - v(2)) ./ (2 * z_y)];
    
    w1 = a * (aol.half_width ./ aol.V)^1 / 1;
    w2 = b * (aol.half_width ./ aol.V)^2 / 2;

    x1 = [w1(1), w2(1), w3, w4(1), w5]; % 36 waves of w2f gives 1 m focal length or 64 um after obj. 
    y1 = [w1(2), w2(2), w3, w4(2), w5];    
    x2 = [w1(3), w2(3), w3, w4(3), w5];
    y2 = [w1(4), w2(4), w3, w4(4), w5];   
    aod_drives_in_waves = {x1, y1, x2, y2};
    
    res = get_psf(aol, time, ws, wf, aod_drives_in_waves, plot);
end
