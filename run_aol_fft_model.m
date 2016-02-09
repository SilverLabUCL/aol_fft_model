function res = run_aol_fft_model(num_aods, use_drive_eqs, aol, times, xy_def, z_depth, xy_vel, w3, w4, w5, ws, wf, plot)
    if use_drive_eqs
        waves_aods = indirect_waves_to_aods(aol, num_aods, xy_def, z_depth, xy_vel, w3, w4, w5);
    else 
        waves_aods = direct_waves_to_aods(num_aods, xy_def, z_depth, xy_vel, w3, w4, w5);
    end 

    waves = build_waves(waves_aods, wf, ws);
    ma = make_model_analyser(aol, waves, times);
    
    res = ma.analyse_and_plot_focus(plot);
end

function waves_aods = indirect_waves_to_aods(aol, num_aods, xy_def, z_depth, xy_vel, w3, w4, w5)   
    if num_aods == 4
        waves_aods = drives4(aol, xy_def, z_depth, xy_vel, w3, w4, w5);
    else
        waves_aods = drives6(aol, xy_def, z_depth, xy_vel, w3, w4, w5);
    end
end

function waves_aods = direct_waves_to_aods(num_aods, w1, w2f, w2s, w3, w4, w5)
    if num_aods == 4
        x1 = [w1, w2f-w2s, w3, w4, w5]; % 36 waves of w2f gives 1 m focal length or 130 um after obj. 
        y1 = [w1, w2f    , w3, w4, w5];    
        x2 = [w1, w2f+w2s, w3, w4, w5];
        y2 = [w1, w2f    , w3, w4, w5];   
        waves_aods = {x1, y1, x2, y2};
    else
        a1 = [w1, w2f-w2s, w3, w4, w5];
        a2 = [w1, w2f    , w3, w4, w5];    
        a3 = [w1, w2f    , w3, w4, w5];
        a4 = [w1, w2f+w2s, w3, w4, w5];    
        a5 = [w1, w2f    , w3, w4, w5];    
        a6 = [w1, w2f    , w3, w4, w5];    
        waves_aods = {a1, a2, a3, a4, a5, a6};
    end
end

function waves = build_waves(waves_aods, wf, ws)
    waves = struct();
    waves.aods = waves_aods;
    waves.focus = wf;
    waves.spherical = ws;
end

function ma = make_model_analyser(aol, waves, times)
    ma = AolFftModelAnalyser(aol, waves);    
    ma.time_range = times;
    
    wavelength = 2*pi / aol.wavevector;
    if wavelength == 800e-9
        ma.wavelengths = [800, 800-1.75, 800+1.75, 800-2.5, 800+2.5] * 1e-9; 
    elseif wavelength == 920e-9
        ma.wavelengths = [920, 920-2.5, 920+2.5, 920-3.5, 920+3.5] * 1e-9;
    else
        ma.wavelengths = wavelength; 
    end
end

function waves4 = drives4(aol, xy_def, z_depth, xy_vel, w3, w4, w5)
    wavelength = 2*pi / aol.wavevector;
    z_x = z_depth + aol.aod_spacing;
    z_y = z_depth;
    L_x = 2 * aol.aod_spacing;
    L_y = 2 * aol.aod_spacing;   
    pdr = 0;
    
    df = - aol.acoustic_velocity ./ wavelength .* xy_def(:) ./ (pdr .* [(L_x + z_x); (L_y + z_y)] + [z_x; z_y]);
    a = [pdr .* df(1,:); pdr .* df(2,:); - df(1,:); - df(2,:)];
    
    b = - aol.acoustic_velocity.^2 ./ wavelength .* [...
        (1 + xy_vel(1)) ./ (L_x .* (1 + xy_vel(1)) + 2 * z_x);...
        (1 + xy_vel(2)) ./ (L_y .* (1 + xy_vel(2)) + 2 * z_y);... 
        (1 - xy_vel(1)) ./ (2 * z_x);...
        (1 - xy_vel(2)) ./ (2 * z_y)];
    
    w1 = a * (aol.aod_half_aperture_width ./ aol.acoustic_velocity)^1 / 1;
    w2 = b * (aol.aod_half_aperture_width ./ aol.acoustic_velocity)^2 / 2;

    x1 = [w1(1), w2(1), w3, w4(1), w5]; % 36 waves of w2f gives 1 m focal length or 64 um after obj. 
    y1 = [w1(2), w2(2), w3, w4(2), w5];    
    x2 = [w1(3), w2(3), w3, w4(3), w5];
    y2 = [w1(4), w2(4), w3, w4(4), w5];   
    waves4 = {x1, y1, x2, y2};
end

function waves6 = drives6(aol, xy, z, v, w3, w4, w5)
    wavelength = (2*pi) ./ aol.wavevector;
    shift = - aol.acoustic_velocity ./ wavelength .* xy; % (-) for -1 mode 
    a = [0; 0; 0;
        shift(2) ./ sqrt(3) ./ (z + aol.aod_spacing*2); 
       -shift(1) ./ (z + aol.aod_spacing);
       -shift(2) ./ sqrt(3) ./ z];
   
    if aol.aod_spacing == 0
        vx = v(1);
        vy = v(2);
        r = [    
                1 - 3/4*vx + sqrt(3)/2*vy 
                1 - 3/4*vx
                1 - 3/4*vx - sqrt(3)/2*vy
                1 + 3/4*vx - sqrt(3)/2*vy
                1 + 3/4*vx
                1 + 3/4*vx + sqrt(3)/2*vy
            ] / (3*z);
    elseif aol.aod_spacing == 4.554e-2
        r = ramps6(z, v);
    else
        error('6aod drive equations unknown for given spacing')
    end
    b = - aol.acoustic_velocity.^2 ./ wavelength .* r;
    
    w1 = a * (aol.aod_half_aperture_width ./ aol.acoustic_velocity)^1 / 1;
    w2 = b * (aol.aod_half_aperture_width ./ aol.acoustic_velocity)^2 / 2;

    a1 = [w1(1), w2(1), w3, w4(1), w5];
    a2 = [w1(2), w2(2), w3, w4(2), w5];    
    a3 = [w1(3), w2(3), w3, w4(3), w5];
    a4 = [w1(4), w2(4), w3, w4(4), w5];    
    a5 = [w1(5), w2(5), w3, w4(5), w5];    
    a6 = [w1(6), w2(6), w3, w4(6), w5];    
    waves6 = {a1, a2, a3, a4, a5, a6};
end
