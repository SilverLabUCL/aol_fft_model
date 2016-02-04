function [x_fwhm, z_fwhm] = plot_na(aol, NA)

aol.half_width = NA * 1e-2;
aol.beam_sigma = NA/2 * 1e-2; 

z_res = 199; 

z_list = [1e-9 50:50:250];
x_list = 0:50:250;
[x_mesh, z_mesh] = meshgrid(x_list, z_list);
points = size(x_mesh);
z_fwhm = zeros(points);
x_fwhm = zeros(points);

for m = 1:points(1)
    for n = 1:points(2)
        z = z_mesh(m,n);
        x = x_mesh(m,n);
        fprintf('%f %f\n', [x z])
        aol.z_list = linspace(z-25,z+25,z_res) * 1e-6;
        z_aol = -130 / z;
        x_aol = z_aol * x * 1e-4; % assumes 0.8 mag relay
        res = psf_4_linear(aol, 0, [x_aol, 0], z_aol, [0,0], 0, [0,0,0,0], 0, 0, 0, false);
        x_fwhm(m,n) = res(3);
        z_fwhm(m,n) = res(4);
    end
end
x_mesh = [fliplr(-x_mesh), x_mesh; fliplr(-x_mesh), x_mesh];
z_mesh = [flipud(-z_mesh), flipud(-z_mesh); z_mesh, z_mesh];
x_fwhm = [flipud(fliplr(x_fwhm)), flipud(x_fwhm); fliplr(x_fwhm), x_fwhm];
z_fwhm = [flipud(fliplr(z_fwhm)), flipud(z_fwhm); fliplr(z_fwhm), z_fwhm];

ellipse_grid(x_mesh, z_mesh, x_fwhm, z_fwhm);
end