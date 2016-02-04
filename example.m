aol = aol_fft();
aol.adjustment = 400;
aol.number_of_samples = 2^9 - 1;
aol.k = 2*pi/920e-9;
aol.spacing = 0;

NA = 0.6;
aol.half_width = NA * 1e-2;
aol.beam_sigma = NA/2 * 1e-2; 

z_res = 199; 

z = 1e-8;
x = 250;

fprintf('%f %f\n', [x z])
aol.z_list = linspace(z-25,z+25,z_res) * 1e-6;
z_aol = -130 / z;
x_aol = z_aol * x * 1e-4; % assumes 0.8 mag relay
res = psf_4_linear(aol, 0, [x_aol, 0], z_aol, [0,0], 0, [0,0,0,0], 0, 0, 0, true)
