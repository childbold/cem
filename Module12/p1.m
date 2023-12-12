% close all;
clear;
clc;
f = 3e9;
lambda = 3e8/f;
a = 5; % in wavelengths
dphi = 0.1;
N = floor(2*pi*a/dphi);
thi = 315*(pi/180);
k0 = 2*pi;
Z0 = 120*pi;
phi = linspace(0,2*pi,N);
V = zeros(N,1);
Z = zeros(N,N);
for n=1:N
    for m=1:N
        if n==m
            Z(n,m) = 0.5 - (dphi/(4*pi)) + 1j*(k0*a).^2*(dphi).^3./192;
        else
            arg = 2*a*k0*abs(sin((phi(n) - phi(m))/2));
            Z(n,m) = 1j*k0*a*dphi*abs(sin((phi(n) - phi(m))/2))*besselh(1,2,arg)./4;
        end
    end
    V(n) = 2*cos(k0*a*sin(phi(m))*sin(thi))*exp(-1j*k0*a*cos(phi(m)*cos(thi)));
end
J = V.'/Z;

% Jpo = zeros(N,1);
% 
% for i=N/4+1:N/2
%     Jpo(i) = (2/Z0)*cos(thi)*exp(-1i*k0*(x(i)*sin(thi)-y(i)*cos(thi)));
% end
% 
% for i=N/2+1:3*N/4
%     Jpo(i) = (2/Z0)*sin(thi)*exp(-1i*k0*(x(i)*sin(thi)-y(i)*cos(thi)));
% end

figure; hold on;
plot(phi,abs(J),'linewidth',3);
% plot(l*lambda,abs(Jpo),'linewidth',3);
legend({'MoM','PO'})
% xlim([0 max(l*lambda)])
% ylim([0 max(abs(J))])
grid on
box on
xlabel('l-m')
ylabel('|J| (A/m)')

% dths = 0.1;
% Nths = floor(360/dths)+1;
% ths = linspace(0,360,Nths);
% Es = zeros(Nths,1);
% Espo = zeros(Nths,1);
% 
% for i=1:Nths
%     sum = 0;
%     for j=1:N
%         sum = sum + J(j)*exp(1i*k0*(x(j)*cosd(ths(i))+y(j)*sind(ths(i))));
%     end
%     Es(i) = -(k0*Z0/4)*sqrt(2/(pi*k0))*exp(1i*pi/4)*dphi*sum;
% end
% 
% for i=1:Nths
%     sum = 0;
%     for j=1:N
%         sum = sum + Jpo(j)*exp(1i*k0*(x(j)*cosd(ths(i))+y(j)*sind(ths(i))));
%     end
%     Espo(i) = -(k0*Z0/4)*sqrt(2/(pi*k0))*exp(1i*pi/4)*dphi*sum;
% end
% 
% figure(2); hold on;
% plot(ths,20*log10(abs(Es)),'linewidth',3);
% plot(ths,20*log10(abs(Espo)),'linewidth',3);
% legend({'MoM','PO'})
% xlim([0 360])
% ylim([-70 20])
% grid on
% box on
% xlabel('\theta_s^o')
% ylabel('|E_s|-dB')
% set(gca,'fontsize',24)


