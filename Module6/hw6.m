x0 = linspace(-5,5,1001);
y0 = 5;
k = 2*pi;
rho = sqrt(x0.^2 + y0.^2);

a = 90 + atand(x0/y0);
phi = 360 - a;


E_incident = exp(1i*k*rho);

E_GO = E_incident;
E_GO(x0 > 0) = 0;

figure;
plot(x0,abs(E_GO),'k','DisplayName','|E_{Go}|')

xlabel('x(\lambda)')
ylabel('|E|')
