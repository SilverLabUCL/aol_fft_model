NA = 0.7;
z_res = 199; 
z_um = 0;
x_um = 0;
do_plot = true;

aol = AolFftModel();
aol.fft_adjustment = 400;
aol.fft_number_of_samples = 2^8 - 1;
aol.wavevector = 2*pi/920e-9;
aol.aod_spacing = 0e-2;
aol.aod_half_aperture_width = NA * 1.1e-2;
aol.beam_width = NA * 0.7e-2; 

res_4_aod_corrected = [];
res_6_aod_corrected = [];
res_no_correction = [];

z_step = 100; % 10 slow but good resolution, 100 fast but poor resolution
z_list = -400:z_step:400;
for z_um = z_list
    do_plot = false;

    aol.z_range = linspace(z_um-20, z_um+30, z_res) * 1e-6;
    z_aol = -130 / (z_um + 1e-12);
    x_aol = z_aol * x_um * 1e-4; % assumes 0.8 mag relay
    fprintf('%f %f\n', [x_um z_um])

    res_4_aod_corrected = [res_4_aod_corrected; run_aol_fft_model(4, 1, aol, 0, [x_aol, 0], z_aol, [0,0], 0, ones(1,4)*4.2/400 * z_um, 0, 0, 0, do_plot)]; 
    res_6_aod_corrected = [res_6_aod_corrected; run_aol_fft_model(6, 1, aol, 0, [x_aol, 0], z_aol, [0,0], 0, ones(1,6)*4.8/400 * z_um, 0, 0, 0, do_plot)];
    res_no_correction =   [res_no_correction;   run_aol_fft_model(4, 1, aol, 0, [x_aol, 0], z_aol, [0,0], 0, ones(1,4)*0,              0, 0, 0, do_plot)];
end

figure; hold on; plot(z_list, res_4_aod_corrected(:, 4)); plot(z_list, res_6_aod_corrected(:, 4),'k'); plot(z_list, res_no_correction(:, 4), 'r');