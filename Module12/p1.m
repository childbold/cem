close all;
clear;
clc;

% Parameters
f = 3e9;
a = .5;            % radius in meters;
delta = 0.1;        % spacing in wavelengths
thi = 315*(pi/180); % Incident angle

c = 3e8;
a = a*f/c;          % radius in wavelengths
N = floor(2*pi*a/delta);
lambda = 3e8/f;
k0 = 2*pi;
Z0 = 120*pi;
phi = linspace(0,2*pi,N);

% Initialization
V = zeros(N,1);
Z = zeros(N,N);

% Solving for Z matrix
for n=1:N
    for m=1:N
        if n==m
            Z(n,m) = 0.5 - (delta/(4*pi*a)) - 1j*(k0*a).^2*(delta/4/a - 0.5*sin(delta/2/a))/2;
        else
            arg = sin((phi(m) - phi(n))/2);
            Z(n,m) = -1j*k0*delta*arg*besselh(1,2,2*k0*a*arg)./4;
        end
    end
    V(n) = exp(-1j*k0*a*cos(phi(n) - thi)) + exp(-1j*k0*a*cos(phi(n) + thi));
end
% Hz = V.'/Z;
% Hz = inv(Z)*V;
Hz = Z/V.';

figure; hold on; grid on;
plot(phi*180/pi,abs(V))

figure; hold on; grid on;
plot(phi*180/pi,abs(sum(Z,2)))

figure; hold on; grid on;
plot(phi*180/pi,abs(Hz))

figure; hold on; grid on;
plot(phi*180/pi,abs(V)./abs(sum(Z,2)))



dths = 0.1;
Nths = floor(360/dths)+1;
ths = linspace(0,2*pi,Nths);
Es = zeros(Nths,1);

for i=1:Nths
    sum = 0;
    for j=1:N
        sum = sum + Hz(j)*exp(1i*k0*(a*cos(ths(i) - phi(j))))*cos(ths(i) - phi(j));
    end
    Es(i) = -(k0/4)*sqrt(2/(pi*k0))*exp(1i*pi/4)*delta*sum;
end


figure; hold on;
plot(ths,20*log10(abs(Es)),'linewidth',3);



