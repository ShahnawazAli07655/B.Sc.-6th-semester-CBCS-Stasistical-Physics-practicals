E = linspace(0, 2, 500)
kB = 8.617e-5
Ef = 1
Tval = [300, 600, 900]
subplot(3,1,1)
title("Maxwell-Boltzmann Distribution")
xlabel("Energy (eV) -->")
ylabel("f(E) -->")
for i = 1:length(Tval)
    T = Tval(i)
    f_MB = exp(-E ./ (kB * T))
    plot(E, f_MB, i, leg="T = " + string(T) + "K", style=i)
end
legend("T = 300K", "T = 600K", "T = 900K")
subplot(3,1,2)
title("Fermi-Dirac Distribution")
xlabel("Energy (eV) -->")
ylabel("f(E) -->")
for i = 1:length(Tval)
    T = Tval(i)
    f_FD = 1 ./ (exp((E - Ef) ./ (kB * T)) + 1)
    plot(E, f_FD, i, leg="T = " + string(T) + "K", style=i)
end
legend("T = 300K", "T = 600K", "T = 900K")
subplot(3,1,3)
title("Bose-Einstein Distribution")
xlabel("Energy (eV) -->")
ylabel("f(E) -->")
for i = 1:length(Tval)
    T = Tval(i)
    denom = exp((E - Ef) ./ (kB * T)) - 1
    f_BE = 1 ./ denom
    f_BE(find(denom <= 0)) = %nan
    plot(E, f_BE, i, leg="T = " + string(T) + "K", style=i)
end
legend("T = 300K", "T = 600K", "T = 900K")
