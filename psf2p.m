% theoretical psf dimensions from Zipfel 2003, fig4

na = 0.5:0.01:1;
lambda = 920e-9;

wx = 0.32  * lambda / sqrt(2) ./ na;
wxx = 0.325 * lambda / sqrt(2) ./ na.^0.91;
wz = 0.532 * lambda / sqrt(2) ./ (1.33 - sqrt(1.33^2 - na.^2));

% w = sqrt(2) sigma
% fwhm = 2 * sqrt(ln2)
figure
hold on
plot(na, wx * 2 * sqrt(log(2)),'b-')
plot(na, wxx * 2 * sqrt(log(2)),'b.')
plot(na, wz * 2 * sqrt(log(2)),'r')
xlabel('numerical aperture')
ylabel('fwhm (um)')