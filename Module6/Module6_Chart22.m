%This script reproduces the result on Chart 22 from the Lecture Notes of
%Module 6
% close all
% clear 
% clc

k = 2*pi;         %Lengths are normalized by lamda
E = 1;            

h = 0.5;          %h is normalized by lambda
w = 2;            %w is normalized by lambda
s = sqrt(h^2 + w^2/4);
p = s;
a = atand(2*h/w); %alpha
ph = 0:0.1:360;   

Ego = zeros(length(ph),1);
for i=1:length(ph)
    if (ph(i)>=a)&&(ph(i)<=180-a) 
        Ego(i) = 1j*2*E*sin(k*h*sind(ph(i)));
    elseif ((ph(i)>=0)&&(ph(i)<=a))||((ph(i)>=360-a)&&(ph(i)<=360))||((ph(i)>=180-a)&&(ph(i)<=180+a))
        Ego(i) = E*exp(1j*k*h*sind(ph(i)));
    end
end

R = s; 
b0d = 90;
n = 2;
nu = 1;
Ed1 = zeros(length(ph),1);
Ed2 = zeros(length(ph),1);
for m=1:length(ph)
    if ph(m)>=0&&ph(m)<=180
        psi1 = 180-ph(m);
    else
        psi1 = 540-ph(m);
    end
    psi2 = ph(m);
    [DH1,DV1] = wdca(R,psi1,a,b0d,n,nu);
    [DH2,DV2] = wdca(R,psi2,a,b0d,n,nu);
    Ed1(m) = E*(exp(-1j*k*s)/sqrt(s))*DH1*exp(1j*k*(w/2)*cosd(ph(m)));
    Ed2(m) = E*(exp(-1j*k*s)/sqrt(s))*DH2*exp(-1j*k*(w/2)*cosd(ph(m)));
end
Et = Ed1+Ed2+Ego;
figure(1)
polarplot(ph*pi/180, abs(Et), 'b')
hold on
polarplot(ph*pi/180, abs(Ego), 'r')
rlim([0 5]);
hold off
Ego = Ego+min(abs(Et));
Etmin = min(20*log10(abs(Et)));
Egomin = min(20*log10(abs(Ego)));
figure(2)
polarplot(ph*pi/180, 20*log10(abs(Et))-Etmin, 'b')
hold on;
polarplot(ph*pi/180, 20*log10(abs(Ego))-Egomin, 'r')