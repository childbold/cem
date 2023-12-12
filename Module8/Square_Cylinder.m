% close all;
clear;
clc;
f = 3e9;
lambda = 3e8/f;
w = 1/lambda;
L = 4*w;
dl = 0.1;
N = floor(L/dl);
thi = 45*(pi/180);
a = 1.781;
e = 2.7183;
k0 = 2*pi;
Z0 = 120*pi;
x = zeros(N,1);
y = zeros(N,1);
for i = 1:N
    if i<=N/4
        x(i) = w/2;
        y(i) = -w/2 + (i-0.5)*dl;
    elseif i>N/4&&i<=N/2
        j = i - N/4;
        x(i) = w/2 - (j-0.5)*dl;
        y(i) = w/2;
    elseif i>N/2&&i<=3*N/4
        j= i - N/2;
        x(i) = -w/2;
        y(i) = w/2 - (j-0.5)*dl;
    else
        j = i - 3*N/4;
        x(i) = -w/2 + (j-0.5)*dl;
        y(i) = -w/2;
    end
    
end

l = linspace(dl/2,L-dl/2,N);
V = zeros(N,1);
for i=1:N
    V(i) = exp(-1i*k0*(x(i)*sin(thi)-y(i)*cos(thi)));
end
Z = zeros(N,N);
for i=1:N
    for j=1:N
        if i==j
            Z(i,j) = (k0*Z0*dl/4)*(1-1i*(2/pi)*log(a*k0*dl/(4*e)));
        else
            arg = k0*sqrt((x(j)-x(i))^2+(y(j)-y(i))^2);
            Z(i,j) = (k0*Z0*dl/4)*besselh(0,2,arg);
        end
    end
end
J = V.'/Z;

Jpo = zeros(N,1);

for i=N/4+1:N/2
    Jpo(i) = (2/Z0)*cos(thi)*exp(-1i*k0*(x(i)*sin(thi)-y(i)*cos(thi)));
end

for i=N/2+1:3*N/4
    Jpo(i) = (2/Z0)*sin(thi)*exp(-1i*k0*(x(i)*sin(thi)-y(i)*cos(thi)));
end

figure(1); hold on;
plot(l*lambda,abs(J),'linewidth',3);
plot(l*lambda,abs(Jpo),'linewidth',3);
legend({'MoM','PO'})
xlim([0 max(l*lambda)])
ylim([0 max(abs(J))])
grid on
box on
xlabel('l-m')
ylabel('|J| (A/m)')
set(gca,'fontsize',24)

dths = 0.1;
Nths = floor(360/dths)+1;
ths = linspace(0,360,Nths);
Es = zeros(Nths,1);
Espo = zeros(Nths,1);

for i=1:Nths
    sum = 0;
    for j=1:N
        sum = sum + J(j)*exp(1i*k0*(x(j)*cosd(ths(i))+y(j)*sind(ths(i))));
    end
    Es(i) = -(k0*Z0/4)*sqrt(2/(pi*k0))*exp(1i*pi/4)*dl*sum;
end

for i=1:Nths
    sum = 0;
    for j=1:N
        sum = sum + Jpo(j)*exp(1i*k0*(x(j)*cosd(ths(i))+y(j)*sind(ths(i))));
    end
    Espo(i) = -(k0*Z0/4)*sqrt(2/(pi*k0))*exp(1i*pi/4)*dl*sum;
end

figure(2); hold on;
plot(ths,20*log10(abs(Es)),'linewidth',3);
plot(ths,20*log10(abs(Espo)),'linewidth',3);
legend({'MoM','PO'})
xlim([0 360])
ylim([-70 20])
grid on
box on
xlabel('\theta_s^o')
ylabel('|E_s|-dB')
set(gca,'fontsize',24)


