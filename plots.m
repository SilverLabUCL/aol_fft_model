aol = aol_fft();
aol.adjustment = 400;
aol.number_of_samples = 2^9 - 1;
aol.k = 2*pi/920e-9;
aol.spacing = 0;

[x_fwhm6, z_fwhm6] = plot_na(aol, 0.6);
[x_fwhm8, z_fwhm8] = plot_na(aol, 0.8);

% aol.k = 2*pi/800.1e-9;
% 
% [x_fwhm6c, z_fwhm6c] = plot_na(aol, 0.6);
% [x_fwhm8c, z_fwhm8c] = plot_na(aol, 0.8);