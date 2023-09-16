% close all;
% clear;
% clc;
pol = 'TEz';                %Polarization
f = 1.0e9;                 %Frequency(Hz)
a = 0.25;                    %Radius(m)
k = 2*pi*f/3.0e8;           %Wavenumber
ka = k*a;                   %Wavenumber times radius
nm = 200;                   %Number of modes
nph = 361;                  %Number of phi samples
ph = linspace(0,2*pi,nph);
current = zeros(nph,1);
current_po = zeros(nph,1);
sum = zeros(nm,1);
if strcmp(pol,'TEz')
    sum(1) = -1/besselh(1,2,ka);
elseif strcmp(pol,'TMz')
    sum(1) = 1/besselh(0,2,ka);
end
if strcmp(pol,'TEz')
    H0=1/(120*pi);
    for i=1:nph
        for j=2:nm   
            n = j-1;
            CosineT=2*cos(n*ph(i));
            CT=1j^(-n);
            HT=(n/ka)*besselh(n,2,ka)-besselh(n+1,2,ka);
            sum(j)=sum(j-1)+CT*CosineT/HT;
        end
        current(i)=1j*(2*H0/(pi*ka))*sum(nm);
        if (ph(i)<pi/2)||(ph(i)>3*pi/2)
            current_po(i)=0;
        else
            current_po(i)=-2*H0*exp(-1j*ka*cos(ph(i)));
        end
    end
elseif strcmp(pol,'TMz')
    E0=1;
    for i=1:nph
        for j=2:nm
            n = j-1;
            CosineT=2*cos(n*ph(i));
            CT=1j^(-n);
            HT=besselh(n,2,ka);
            sum(j)=sum(j-1)+CT*CosineT/HT;
        end
        current(i)=(2*E0/(pi*ka*120*pi))*sum(nm);
        if (ph(i)<pi/2)||(ph(i)>3*pi/2)
            current_po(i)=0;
        else
            current_po(i)=-(2.0*E0/(120*pi))*cos(ph(i))*exp(-1j*ka*cos(ph(i)));
        end
    end
end
figure(1);
plot(ph*180/pi,abs(current),'linewidth',3);
hold on;
plot(ph*180/pi,abs(current_po),'g--','linewidth',3);
xlim([0 360]);
xlabel('\phi^o')
if strcmp(pol,'TMz')
    ylabel('|J_z| (A/m)');
else
    ylabel('|J_\phi| (A/m)');
end
grid on
% box on
% set(gca,'fontsize',36)
% keyboard
        
        