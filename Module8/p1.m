% Parameters
f = 3e9;
c = 3e8;
lambda = c/f;
th0 = 45;
dl = 1/10;
w = 1/lambda;
N = floor(4*w/dl);
k0 = 2*pi;
Z0 = 120*pi;
alpha = 1.781;

% Defining surface
x = zeros(N,1);
y = zeros(N,1);
for n = 1:N
    l_n = (n - 1/2)*dl;
    if n <= N/4
        x(n) = w/2;
        y(n) = l_n - w/2;
    elseif n > N/4 && n <= N/2
        x(n) = 3*w/2 - l_n;
        y(n) = w/2;
    elseif n > N/2 && n <= 3*N/4
        x(n) = -w/2;
        y(n) = 5*w/2 - l_n;
    elseif n > 3*N/4
        x(n) = l_n - 7*w/2;
        y(n) = -w/2;
    end
end

% Incident field
l = linspace(dl/2,4*w-dl/2,N);
V = zeros(N,1);
for n=1:N
    V(n) = exp(-1i*k0*(x(n)*sind(th0)-y(n)*cosd(th0)));
end

Z = zeros(N);
for n = 1:N

    for m = 1:N
        R = sqrt((x(m) - x(n))^2 + (y(m) - y(n))^2);
        if n == m
            
            Z(n,m) = (k0*Z0/4)*dl*(1 - 2j*log(alpha*k0*dl/(4*exp(1)))/pi);

        else

            Z(n,m) = (k0*Z0/4)*dl*besselh(0,2,k0*R);

        end
    end
end

J = V.'/Z;

Jpo = zeros(N,1);
for i=N/4+1:N/2
    Jpo(i) = (2/Z0)*cosd(th0)*exp(-1i*k0*(x(i)*sind(th0)-y(i)*cosd(th0)));
end
for i=N/2+1:3*N/4
    Jpo(i) = (2/Z0)*sind(th0)*exp(-1i*k0*(x(i)*sind(th0)-y(i)*cosd(th0)));
end

figure; hold on; grid on;
plot(l*lambda,abs(J),'linewidth',2);
plot(l*lambda,abs(Jpo),'linewidth',2);
legend({'MoM','PO'})
xlabel('l [m]')
ylabel('|J| [A/m]')
saveas(gcf,'InducedCurrents.png')

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

figure; hold on; grid on;
plot(ths,20*log10(abs(Es)),'linewidth',2);
plot(ths,20*log10(abs(Espo)),'linewidth',2);
legend({'MoM','PO'})
xlim([0 360])
xlabel('\phi_s')
ylabel('|E^s_Z|')
saveas(gcf,'ScatteredField.png')