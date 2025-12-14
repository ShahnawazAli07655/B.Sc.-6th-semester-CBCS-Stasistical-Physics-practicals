h = 6.626e-34
c = 3e8
kB = 1.381e-23

lambda = linspace(100e-9, 3000e-9, 500)
lambda_nm = lambda * 1e9

Tvals = [3000, 6000]

for i = 1:length(Tvals)
    T = Tvals(i);

    exponent = h * c ./ (lambda * kB * T);
    

    exponent(exponent > 700) = 700

    B = ((2 * h * c^2) ./ (lambda.^5)).*(1 ./ (exp(exponent) - 1));
    B_RJ = (2 * c * kB * T) ./ (lambda.^4);

    B = B ./ max(B);
    B_RJ = B_RJ ./ max(B_RJ);

    subplot(2,1,i);
    plot(lambda_nm, B, 'b-');     // Planck
    plot(lambda_nm, B_RJ, 'r-'); // Rayleigh-Jeans
    title("Spectral Radiance at T = " + string(T) + " K");
    xlabel("Wavelength (nm)");
    ylabel("Normalized Radiance");
    legend("Plancks Law", "Rayleigh-Jeans Law");
end
