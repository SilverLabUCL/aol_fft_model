NA = 0.6;
z_res = 199; 
z_um = 0;
x_um = 0;

aol = AolFftModel();
aol.fft_adjustment = 400;
aol.fft_number_of_samples = 2^8 - 1;
aol.wavevector = 2*pi/920e-9;
aol.aod_spacing = 0e-2;
aol.aod_half_aperture_width = NA * 1e-2;
aol.beam_width = NA/1.5 * 1e-2; 
aol.z_range = linspace(z_um-50, z_um+50, z_res) * 1e-6;

fprintf('%f %f\n', [x_um z_um])

z_aol = -130 / (z_um + 1e-12);
x_aol = z_aol * x_um * 1e-4; % assumes 0.8 mag relay
%res = run_aol_fft_model(4, 1, aol, [-10e-6 0 10e-6], [x_aol, 0], z_aol, [0,0], 3, [0,0,0,0,0,0], 0, 0, 0, true);
%res2 = run_aol_fft_model(6, 1, aol, [-10e-6 0 10e-6], [x_aol, 0], z_aol, [0,0], 2, [0,0,0,0,0,0], 0.1, 0, 0, true);
%res2 = run_aol_fft_model(6, 1, aol, [-10e-6 0 10e-6], [x_aol, 0], z_aol, [0,0], 2.5, [0,0,0,0,0,0], 0, 0, 0, true);
res2 = run_aol_fft_model(6, 1, aol, [-10e-6 0 10e-6], [x_aol, 0], z_aol, [0,0], 0, [0,0,0,0,0,0], 0, 0, 0, true);