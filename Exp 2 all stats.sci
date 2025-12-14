kB = 1.380649e-23
epsilon = 1e-21
E = [0, 1] * epsilon
N = 2
T = linspace(10, 1000, 100)
U_MB = zeros(1, length(T)); U_FD = zeros(1, length(T)); U_BE = zeros(1, length(T));
Cv_MB = zeros(1, length(T)); Cv_FD = zeros(1, length(T)); Cv_BE = zeros(1, length(T));
deltaE_MB = zeros(1, length(T)); deltaE_FD = zeros(1, length(T)); deltaE_BE = zeros(1, length(T));
ratio_MB = zeros(1, length(T)); ratio_FD = zeros(1, length(T)); ratio_BE = zeros(1, length(T));

function states = gen_states_FD_BE(N, type)
    states = [];
    for n1 = 0:N
        n2 = N - n1;
        if type == "FD" & (n1 > 1 | n2 > 1) then
            continue;
        end
        states = [states; n1, n2];
    end
endfunction

for i = 1:length(T)
    temp = T(i);
    beta = 1 / (kB * temp);
    
    // ===== MAXWELL-BOLTZMANN =====
    z1 = sum(exp(-beta * E));
    Z_MB = z1^N;
    probs = exp(-beta * E) / z1;
    avg_e = sum(E .* probs);
    avg_e2 = sum(E.^2 .* probs);
    U_MB(i) = N * avg_e;
    deltaE_MB(i) = sqrt(N * (avg_e2 - avg_e^2));
    Cv_MB(i) = kB * beta^2 * N * (avg_e2 - avg_e^2)
    ratio_MB(i) = probs(2) / probs(1);
    
    // ===== FERMI-DIRAC =====
    states_FD = gen_states_FD_BE(N, "FD");
    Z_FD = 0; EU_FD = 0; EU2_FD = 0; n1_FD = 0; n2_FD = 0;
    for s = 1:size(states_FD, 1)
        occ = states_FD(s, :);
        Etot = sum(occ .* E);
        w = exp(-beta * Etot);
        Z_FD = Z_FD + w;
        EU_FD = EU_FD + Etot * w;
        EU2_FD = EU2_FD + Etot^2 * w;
        n1_FD = n1_FD + occ(1) * w;
        n2_FD = n2_FD + occ(2) * w;
    end
    U_FD(i) = EU_FD / Z_FD;
    deltaE_FD(i) = sqrt(EU2_FD/Z_FD - (U_FD(i))^2);
    Cv_FD(i) = kB * beta^2 * (EU2_FD/Z_FD - (U_FD(i))^2)
    ratio_FD(i) = (n2_FD/Z_FD) / (n1_FD/Z_FD);
    
    // ===== BOSE-EINSTEIN =====
    states_BE = gen_states_FD_BE(N, "BE");
    Z_BE = 0; EU_BE = 0; EU2_BE = 0; n1_BE = 0; n2_BE = 0;
    for s = 1:size(states_BE, 1)
        occ = states_BE(s, :);
        Etot = sum(occ .* E);
        w = exp(-beta * Etot);
        Z_BE = Z_BE + w;
        EU_BE = EU_BE + Etot * w;
        EU2_BE = EU2_BE + Etot^2 * w;
        n1_BE = n1_BE + occ(1) * w;
        n2_BE = n2_BE + occ(2) * w;
    end
    U_BE(i) = EU_BE / Z_BE;
    deltaE_BE(i) = sqrt(EU2_BE/Z_BE - (U_BE(i))^2);
    Cv_BE(i) = kB * beta^2 * (EU2_BE/Z_BE - (U_BE(i))^2)
    ratio_BE(i) = (n2_BE/Z_BE) / (n1_BE/Z_BE);
end

subplot(2,2,1)
plot(T, [U_MB; U_FD; U_BE]')
legend("MB", "FD", "BE", pos="lower right")
xlabel("Temperature (K)"); 
ylabel("⟨E⟩ (J)"); 
title("Average Energy");

subplot(2,2,2)
plot(T, [deltaE_MB; deltaE_FD; deltaE_BE]')
legend("MB", "FD", "BE", pos="upper left")
xlabel("Temperature (K)"); 
ylabel("ΔE (J)"); 
title("Energy Fluctuations");

subplot(2,2,3)
plot(T, [Cv_MB; Cv_FD; Cv_BE]')
legend("MB", "FD", "BE", pos="upper left")
xlabel("Temperature (K)"); 
ylabel("C_v (J/K)"); 
title("Specific Heat");

subplot(2,2,4)
plot(T, [ratio_MB; ratio_FD; ratio_BE]')
legend("MB", "FD", "BE", pos="lower right")
xlabel("Temperature (K)"); 
ylabel("P₂/P₁"); 
title("Occupation Ratio");
