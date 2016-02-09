aol = AolFftModel();
aol.fft_adjustment = 400;
aol.fft_number_of_samples = 2^9 - 1;
aol.wavevector = 2*pi/920e-9;
aol.aod_spacing = 0;

z_list = [1e-9 50:50:250];
x_list = 0:50:250;
[x_mesh, z_mesh] = meshgrid(x_list, z_list);

[x_fwhm6, z_fwhm6, max_fl6, total_fl6] = plot_na(aol, 0.6, x_mesh, z_mesh);
[x_fwhm7, z_fwhm7, max_fl7, total_fl7] = plot_na(aol, 0.7, x_mesh, z_mesh);
[x_fwhm8, z_fwhm8, max_fl8, total_fl8] = plot_na(aol, 0.8, x_mesh, z_mesh);

x_mesh = [fliplr(-x_mesh), x_mesh; fliplr(-x_mesh), x_mesh];
z_mesh = [flipud(-z_mesh), flipud(-z_mesh); z_mesh, z_mesh];

