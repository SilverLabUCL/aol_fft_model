function [x_fwhm, z_fwhm, max_fl, total_fl] = plot_na(aol, numerical_aperture, x_mesh, z_mesh)

z_res = 199; 

aol.aod_half_aperture_width = numerical_aperture * 1e-2;
aol.beam_width = numerical_aperture/1.5 * 1e-2; 

points = size(x_mesh);
z_fwhm = zeros(points);
x_fwhm = zeros(points);
max_fl = zeros(points);
total_fl = zeros(points);

for m = 1:points(1)
    for n = 1:points(2)
        z = z_mesh(m,n);
        x = x_mesh(m,n);
        fprintf('%f %f\n', [x z])
        aol.z_range = linspace(z-25,z+25,z_res) * 1e-6;
        z_aol = -130 / z;
        x_aol = z_aol * x * 1e-4; % assumes 0.8 mag relay
        res = run_aol_fft_model(4, 1, aol, 0, [x_aol, 0], z_aol, [0,0], 0, [0,0,0,0], 0, 0, 0, false);
        x_fwhm(m,n) = res(3);
        z_fwhm(m,n) = res(4);
        max_fl(m,n) = res(5);
        total_fl(m,n) = res(6);
    end
end

x_fwhm = expand(x_fwhm);
z_fwhm = expand(z_fwhm);
max_fl = expand(10.^(max_fl - max(max_fl(:))));
total_fl = expand(10.^(total_fl - max(total_fl(:))));

ellipse_grid(x_mesh, z_mesh, x_fwhm, z_fwhm);
ellipse_grid(x_mesh, z_mesh, max_fl, total_fl);
end

function expanded = expand(mat)
    expanded = [rot90(mat,2), flipud(mat); fliplr(mat), mat];
end