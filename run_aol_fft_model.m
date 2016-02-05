function res = run_aol_fft_model(num_aods, use_drive_eqs, aol, times, xy_def, z_depth, xy_vel, w3, w4, w5, ws, wf, plot)
    if use_drive_eqs
        waves_aods = indirect_waves_to_aods(aol, num_aods, xy_def, z_depth, xy_vel, w3, w4, w5);
    else 
        waves_aods = direct_waves_to_aods(num_aods, xy_def, z_depth, xy_vel, w3, w4, w5);
    end 

    waves = build_waves(waves_aods, wf, ws);
    ma = make_model_analyser(aol, waves, times, z_depth);
    
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

function ma = make_model_analyser(aol, waves, times, z_depth)
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
    a = 0e6 + [pdr .* df(1,:); pdr .* df(2,:); - df(1,:); - df(2,:)];
    
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
    a = 0e6 + [0; 0; 0; -shift(2,:) ./ sqrt(3) ./ (z + aol.spacing*2); shift(1,:) ./ (z + aol.spacing); shift(2,:) ./ sqrt(3) ./ z];

    fac = 100000.*(160000000000000000000000000000000000000000.*z.^8 - 14572800000000000000000000000000000000000.*z.^7.*v + 99580800000000000000000000000000000000000.*z.^7 + 331822656000000000000000000000000000000.*z.^6.*v.^2 - 7576617312000000000000000000000000000000.*z.^6.*v + 26536595184000000000000000000000000000000.*z.^6 + 135686017043280000000000000000000000000.*z.^5.*v.^2 - 1643553286103520000000000000000000000000.*z.^5.*v + 3951999537393600000000000000000000000000.*z.^5 + 200714563865692800000000000000000000.*z.^4.*v.^3 + 22052317975196176800000000000000000000.*z.^4.*v.^2 - 192563322410924942400000000000000000000.*z.^4.*v + 359487217756191571200000000000000000000.*z.^4 + 54788839447099736088000000000000000.*z.^3.*v.^3 + 1810806508915961684688000000000000000.*z.^3.*v.^2 - 13129758774547585039584000000000000000.*z.^3.*v + 20428755008473152779328000000000000000.*z.^3 + 18737905508871124612410000000000.*z.^2.*v.^4 + 5450820578828526927195720000000000.*z.^2.*v.^3 + 78301044153351213892295400000000000.*z.^2.*v.^2 - 519401755515190682849172480000000000.*z.^2.*v + 707123870471722746706202880000000000.*z.^2 + 2818555746644394564198712200000.*z.*v.^4 + 231196795398228695070011130000000.*z.*v.^3 + 1660564224525949743596091328800000.*z.*v.^2 - 11001839737339016029887428174400000.*z.*v + 13605832459795647768150598214400000.*z + 105991788852562457586692572281.*v.^4 + 3443384562216748937737147375848.*v.^3 + 13218597203822813286369339129720.*v.^2 - 96087210765642630292884913412064.*v + 111255268037351854413413765243664).^(1./2);
    b = - aol.acoustic_velocity.^2 ./ wavelength .* [... 
    (549973834719289200000.*v - 80802964519200000000000.*z + fac + 28727027970300000000000.*z.*v + 1141210700190000000000.*z.*v.^2 + 421691292000000000000000.*z.^2.*v + 1821600000000000000000000.*z.^3.*v - 1631461392000000000000000.*z.^2 - 13662000000000000000000000.*z.^3 - 40000000000000000000000000.*z.^4 + 84527115437486700000.*v.^2 - 1387081003857555600000)./(18216.*(11385000000000000000.*z.^3.*v - 200000000000000000000.*z.^4 - 68689500000000000000.*z.^3 + 2644211790000000000.*z.^2.*v - 8347413690000000000.*z.^2 + 5902813966500000.*z.*v.^2 + 179511131403450000.*z.*v - 419187240717300000.*z + 439063108456203.*v.^2 + 3436838663017716.*v - 7249353394646484))
    (836708925955993200000.*v - 79281350252280000000000.*z + fac + 44126813740680000000000.*z.*v + 511577210430000000000.*z.*v.^2 + 712036116000000000000000.*z.^2.*v + 3643200000000000000000000.*z.^3.*v - 1615331124000000000000000.*z.^2 - 13662000000000000000000000.*z.^3 - 40000000000000000000000000.*z.^4 + 39127392658341900000.*v.^2 - 1354159567456304400000)./(36432.*(6831000000000000000.*z.^3.*v - 100000000000000000000.*z.^4 - 31878000000000000000.*z.^3 + 1313464680000000000.*z.^2.*v - 3593305237500000000.*z.^2 + 590281396650000.*z.*v.^2 + 82393444949062500.*z.*v - 172263787589025000.*z + 44802358005735.*v.^2 + 1654700422345146.*v - 2962929276112608))
    (348462784489050000000.*v - 58450975632720000000000.*z + fac + 20961548263260000000000.*z.*v + 432873024210000000000.*z.*v.^2 + 352561572000000000000000.*z.^2.*v + 1821600000000000000000000.*z.^3.*v - 1304247384000000000000000.*z.^2 - 12144000000000000000000000.*z.^3 - 40000000000000000000000000.*z.^4 + 32556380150834100000.*v.^2 - 895117926022729200000)./(36432.*(2500000000.*z.^2 + 303600000.*z + 9793377).*(1821600000.*z - 57032019.*v - 683100000.*z.*v + 10000000000.*z.^2 + 86412150))
    (77549858155440000000000.*z - 674225707588527600000.*v - fac - 32557298366340000000000.*z.*v + 432873024210000000000.*z.*v.^2 - 449343180000000000000000.*z.^2.*v - 1821600000000000000000000.*z.^3.*v + 1613026800000000000000000.*z.^2 + 13662000000000000000000000.*z.^3 + 40000000000000000000000000.*z.^4 + 32556380150834100000.*v.^2 + 1252209312794365200000)./(291456.*(625000000.*z.^2 + 71156250.*z + 1728243).*(3529350000.*z - 84683907.*v - 1138500000.*z.*v + 20000000000.*z.^2 + 146324574))
    (85262868405000000000000.*z + 169253352466110000000.*v - fac + 5876579237760000000000.*z.*v + 275464651770000000000.*z.*v.^2 + 48390804000000000000000.*z.^2.*v + 1679852196000000000000000.*z.^2 + 13662000000000000000000000.*z.^3 + 40000000000000000000000000.*z.^4 + 21206449456047900000.*v.^2 + 1462747208489463600000)./(4554.*(800000000000000000000.*z.^4 + 224664000000000000000.*z.^3 + 483908040000000000.*z.^2.*v + 20819567340000000000.*z.^2 + 55224103997700000.*z.*v + 770688881285400000.*z + 26881414803441.*v.^2 + 1358009251551612.*v + 9898998034037508))
    1./(3.*z)];
    
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