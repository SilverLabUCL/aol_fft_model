% Note that the colormap and some parameters are slightly different to those used in figure 6 of 
% Dynamic wavefront shaping with an acousto-optic lens for laser scanning microscopy.

NA = 0.7;
z_res = 199; 
z_um = 0;
x_um = 0;

aol = AolFftModel();
aol.fft_adjustment = 400;
aol.fft_number_of_samples = 2^8 - 1;
aol.wavevector = 2*pi/920e-9;
aol.aod_spacing = 0e-2;
aol.aod_half_aperture_width = NA * 1.1e-2;
aol.beam_width = NA * 0.7e-2; 

aol.z_range = linspace(z_um-60, z_um+20, z_res) * 1e-6;
z_aol = -130 / (z_um + 1e-12);
x_aol = z_aol * x_um * 1e-4; % assumes 0.8 mag relay

for do_plot = [0.5, 2]
    res_not_aberrated =     run_aol_fft_model(4, 1, aol, 0, [x_aol, 0], z_aol, [0,0], 0, zeros(1,4),      0,   0, 0, do_plot); 
    res_aberrated =         run_aol_fft_model(4, 1, aol, 0, [x_aol, 0], z_aol, [0,0], 0, zeros(1,4),      0, -10, 0, do_plot); 
    res_4_aod_corrected =   run_aol_fft_model(4, 1, aol, 0, [x_aol, 0], z_aol, [0,0], 0, ones(1,4)*3.9,   0, -10, 0, do_plot);
    res_6_aod_corrected =   run_aol_fft_model(6, 1, aol, 0, [x_aol, 0], z_aol, [0,0], 0, ones(1,6)*4.45,  0, -10, 0, do_plot);
end