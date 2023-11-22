%This script solves Problem 3 in HW7
close all
clear 
clc

k = 2*pi;         %Lengths are normalized by lamda
E0 = 1;            

a = 5;            %Disc radius
th0 = 4;

th = th0:0.01:180-th0;   
EGO = zeros(length(th),1);
for i=1:length(th)
    if th(i)<=90
        EGO(i) = E0*cos((pi/2)*cos(th(i)*pi/180))/sind(th(i));
    end
end
b0d = 90;
n = 2;
nu = 1;
Ed1 = zeros(length(th),1);
Ed2 = zeros(length(th),1);
psip1 = 0;
psip2 = 0;
for m=1:length(th)
    psi1 = th(m)+90;
    if th(m)>=th0&&th(m)<=90
        psi2 = 90-th(m);
    elseif th(m)>90&&th(m)<=180-th0
        psi2 = 450-th(m);
    end
    [DH1,DV1] = wdca(a,psi1,psip1,b0d,n,nu);
    [DH2,DV2] = wdca(a,psi2,psip2,b0d,n,nu);
    Ed1(m) = E0*(exp(-1j*k*a)/sqrt(a))*(exp(1j*k*a*sind(th(m)))/sqrt(sind(th(m))))*DV1;
    Ed2(m) = E0*(exp(-1j*k*a)/sqrt(a))*(exp(-1j*k*a*sind(th(m)))/sqrt(-sind(th(m))))*DV2;
end
Et = Ed1+Ed2+EGO;
figure(1)
polarplot(th*pi/180, abs(Et),'b','linewidth',2)
hold on
polarplot(th*pi/180, abs(EGO),'r','linewidth',2)
rlim([0 1.5]);
set(gca,'FontSize',14,'FontWeight','bold','ThetaDir','clockwise','ThetaZeroLocation','top')

EGO = EGO+min(abs(Et));
Etmin = min(20*log10(abs(Et)));
Egomin = min(20*log10(abs(EGO)));

figure(2)
polarplot(th*pi/180, 20*log10(abs(Et))-Etmin,'b','linewidth',2)
hold on;
polarplot(th*pi/180, 20*log10(abs(EGO))-Egomin,'r','linewidth',2)
set(gca,'FontSize',14,'FontWeight','bold','ThetaDir','clockwise','ThetaZeroLocation','top')

EGO2 = flip(EGO);
th2 = 180+th0:0.01:360-th0;
th3 = cat(2,th,th2);

figure(1)
polarplot(th2*pi/180, abs(EGO2),'r','linewidth',2)
Et2 = flip(Et);
polarplot(th2*pi/180, abs(Et2),'b','linewidth',2)
grid on;

figure(2)
polarplot(th2*pi/180, 20*log10(abs(Et2))-Etmin,'b','linewidth',2)
hold on;
polarplot(th2*pi/180, 20*log10(abs(EGO2))-Egomin,'r','linewidth',2)
set(gca,'FontSize',14,'FontWeight','bold','ThetaDir','clockwise','ThetaZeroLocation','top')
grid on;