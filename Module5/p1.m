%% Finite strip PO Scattered field
phi = linspace(0,2*pi,721);
phi_0 = 45*pi/180;
phi_RSB = pi - phi_0;
phi_ISB = phi_0 + pi;

E0 = 1;
freq = 1e9;
Z = 120*pi;

L = 10*3e8/freq;
k = 2*pi*freq;

E = -(k*L*E0*sin(phi_0)/sqrt(2*pi*k))*exp(1i*pi/4)*exp(-1i*k)*(exp(1i*k*(cos(phi_0) + cos(phi))*L) - 1)./(1i*k*(cos(phi_0 + cos(phi))));

figure;
plot(phi,20*log10(abs(E)));