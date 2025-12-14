R=8.314
h=6.626*10.^-34
k=1.38*10.^-23
v=15*10^12
for T=10:10:1000
DP=3*R
a=(h*v)/(k*T)
E=(((3*R)*(a^2)*(exp(a)))/((exp(a)-1)^2))
function c=f(x)
    c=(x^4)*exp(x)/((exp(x)-1)^2)
endfunction
integral=intg(0,a,f)
D=(9*R*((1/a)^3)*integral)
plot(T,DP,'r+')
plot(T,E,'bx')
plot(T,D,'g*')
end
xlabel('T  -->','fontsize',4)
ylabel('Cv  -->','fontsize',4)
legend('Dulong-Petit','Einstein','Debye')
title('Specific heat Vs. Temperature Plot','fontsize',4)
