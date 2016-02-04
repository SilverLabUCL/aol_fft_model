function res = psf_6_waves(aol, time, w1, w2f, w2s, w3, w4, w5, ws, wf)
    a1 = [w1, w2f-w2s, w3, w4, w5];
    a2 = [w1, w2f    , w3, w4, w5];    
    a3 = [w1, w2f    , w3, w4, w5];
    a4 = [w1, w2f+w2s, w3, w4, w5];    
    a5 = [w1, w2f    , w3, w4, w5];    
    a6 = [w1, w2f    , w3, w4, w5];    
    aod_drives_in_waves = {a1, a2, a3, a4, a5, a6};
    res = get_psf(aol, time, ws, wf, aod_drives_in_waves, true);
end
