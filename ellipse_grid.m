function ellipse_grid(x_mesh, y_mesh, x_fwhm, y_fwhm)

x_fwhm = min(x_fwhm, 3.132);

% t = linspace(0, 2*pi);
% figure; hold on;
% for n = 1:numel(x_mesh)
%     plot(x_mesh(n) + x_fwhm(n)/2*5 * cos(t), y_mesh(n) + y_fwhm(n)/2*5 * sin(t));
% end
% axis equal
% excess = 40;
% xlim([min(x_mesh(:))-excess, max(x_mesh(:))+excess])
% ylim([min(y_mesh(:))-excess, max(y_mesh(:))+excess])
% 
% figure
hold on
plot(mean(x_mesh), mean(x_fwhm))
%plot(x_mesh, x_fwhm, 'rx')
plot(mean(y_mesh'), mean(y_fwhm'))
%plot(y_mesh, y_fwhm, 'bx')

x_u = unique(x_mesh);
y_u = unique(y_mesh);
wavelength = 920e-3;
wavelength_fwhm = 7e-3 / sqrt(2);
focal_fwhm = abs(y_u) ./ wavelength .* wavelength_fwhm;
angular_fwhm = abs(x_u) ./ wavelength .* wavelength_fwhm;
plot(x_u, angular_fwhm, 'k--')
plot(y_u, focal_fwhm, 'k--')
ylim([0, 10])

end

