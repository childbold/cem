clear all;
close all;
clc;

% Parameters
maxTime = 50e-6;
line_length = 12000;      %Length of transmission line (m)
R = 0.001;
a = 4.5e-3;
b = 7e-3;
epsr = 2.36;

c = 299792458;
eps0 = 8.85419e-12;
mu0 = 1.25663706143592e-06;
L = (mu0/2/pi)*log(b/a);
C = (2*pi*eps0*epsr)/log(b/a);
dz = 10;                %Spatial sampling
vp = 1/sqrt(L*C);       % Phase velocity
dt = dz/vp;             %Temporal sampling
dt = dz/c;

Nz = line_length/dz;
km = floor(Nz/2);
Nt = round(maxTime/dt) + 1;
t = linspace(0,maxTime,Nt);

% Initializatations
V = zeros(Nz+1,Nt);
I = zeros(Nz,Nt);

% Surge current parameters
Ip = 300;
eta = 0.85;
T1 = 1e-6;
T2 = 10e-6;
m = 5;
gamma = t./T1;
I0 = (Ip/eta).*(gamma.^m./(1 + gamma.^m)).*exp(-t/T2);




Zv = linspace(0,line_length,Nz+1);

Zi = linspace(dz/2,line_length-dz/2,Nz);

for n=2:Nt
    
    V(1,n) = V(2,n-1);
    V(Nz+1,n) = V(Nz,n-1);


    for k = 2:Nz
        V(k,n) = V(k,n-1) - (dt/(C*dz))*(I(k,n-1) - I(k-1,n-1));
        if k == km
            V(k,n) = V(k,n) + I0(n);
        end
    end
    for k = 1:Nz
        I(k,n) = I(k,n-1) - (V(k+1,n) - V(k,n))./(dz*((L/dt) - (R/2)));
    end
end
figure; grid on; hold on;
plot(t*1e6,I0)
xlabel('Time (\mus)')
ylabel('Surge Current Source (A)')

figure; hold on; grid on;
plot(Zv/1e3,V(:,round(10e-6/dt)),'.','DisplayName','10 \mus','LineWidth',2)
plot(Zv/1e3,V(:,round(20e-6/dt)),'.','DisplayName','20 \mus','LineWidth',2)
plot(Zv/1e3,V(:,round(30e-6/dt)),'.','DisplayName','30 \mus','LineWidth',2)
plot(Zv/1e3,V(:,round(40e-6/dt)),'.','DisplayName','40 \mus','LineWidth',2)
legend()
xlabel('Distance (km)')
ylabel('Voltage (V)')

figure; hold on; grid on;
plot(Zi/1e3,I(:,round(10e-6/dt)),'.','DisplayName','10 \mus','LineWidth',2)
plot(Zi/1e3,I(:,round(20e-6/dt)),'.','DisplayName','20 \mus','LineWidth',2)
plot(Zi/1e3,I(:,round(30e-6/dt)),'.','DisplayName','30 \mus','LineWidth',2)
plot(Zi/1e3,I(:,round(40e-6/dt)),'.','DisplayName','40 \mus','LineWidth',2)
legend()
xlabel('Distance (km)')
ylabel('Current (A)')

% Voltage animation
figure; hold on; grid on;
subplot(2,1,1)
vPlot = plot(Zv/1e3,V(:,1),'LineWidth',2);
axis([0, 12, 0 250])
hText = text(6.5, 200, '');

subplot(2,1,2)
iPlot = plot(Zi/1e3,I(:,1),'LineWidth',2);
axis([0, 12, -15 15])

for ix = 1:10:Nt
    set(vPlot,'YData',V(:,ix))
    set(iPlot,'YData',I(:,ix))
    title(['Time: ' num2str(ix*dt*1e6) ' (\mus)'])
%     set(hText, 'String', ['Time: ' num2str(ix*dt*1e6) ' (\mus)']);
    pause(0.1)
end