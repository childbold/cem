x0 = linspace(-5,5,1001);
y0 = 5;
k = 2*pi;
rho = sqrt(x0.^2 + y0.^2);

a = 90 + atan2d(x0,y0);
phi = 360 - a;

E_incident = exp(1i*k*rho.*sind(phi));

% Geometric Optics Field
E_GO = E_incident;
E_GO(x0 > 0) = 0;

% Diffracted Field
E_diffacted = zeros(length(rho),2);

for i = 1:length(E_diffacted)
    [DH,DV] = wdca(rho(i),phi(i),90,90,2,1);
    E_diffacted(i,1) = DH*exp(-1i*k*rho(i))/sqrt(rho(i));
    E_diffacted(i,2) = DV*exp(-1i*k*rho(i))/sqrt(rho(i));
end

E_total = E_diffacted + E_GO.';

%Plotting
figure; hold on; grid on
plot(x0,abs(E_GO),'k','LineWidth',2,'DisplayName','|E_{Go}|')
plot(x0,abs(E_diffacted(:,1)),'b--','LineWidth',2,'DisplayName','|E_d|(hard)')
plot(x0,abs(E_diffacted(:,2)),'r--','LineWidth',2,'DisplayName','|E_d|(soft)')
plot(x0,abs(E_total(:,1)),'b','LineWidth',2,'DisplayName','|E_t|(hard)')
plot(x0,abs(E_total(:,2)),'r','LineWidth',2,'DisplayName','|E_t|(soft)')

xlabel('x(\lambda)')
ylabel('|E|')
yticks(0:.2:1.4)
legend()
title('GO and Corrected Fields')
saveas(gcf,'p6_2','png')
