function res = psf_4_waves(aol, time, w1, w2f, w2s, w3, w4, w5, ws, wf)   
    x1 = [w1, w2f-w2s, w3, w4, w5]; % 36 waves of w2f gives 1 m focal length or 130 um after obj. 
    y1 = [w1, w2f    , w3, w4, w5];    
    x2 = [w1, w2f+w2s, w3, w4, w5];
    y2 = [w1, w2f    , w3, w4, w5];   
    aod_drives_in_waves = {x1, y1, x2, y2};
    res = get_psf(aol, time, ws, wf, aod_drives_in_waves, true);
end
    

