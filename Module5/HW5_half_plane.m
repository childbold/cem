%Edu: Solution via uniform theory of diffraction (UTD)
%Edk: Solution via Keller's geometrical theory of diffraction (GTD)
%Epo: PO solution for the half-plane obtained by adding the GO reflected 
%field due to the stationary-point contributio, the end-point (diffraction)
%contribution,and the incident field
%Epo: PO solution for the finite strip of length L

close all;
clear;
clc;
k=2*pi;
php = 45;
dph = 0.01*pi/180;
ph = 0:dph:2*pi;
nph = length(ph);
btd = 90;
n = 2;
r = 10;
ph0 = php*pi/180;
sph0 = sin(ph0);
cph0 = cos(ph0);
L = 10;
eps = dph/8;
Edu = zeros(nph,1);        
Edk = zeros(nph,1);
Edkd = zeros(nph,1);    
Edp = zeros(nph,1);        
Epo = zeros(nph,1);
Epod = zeros(nph,1);
Epoi = zeros(nph,1);
Epof = zeros(nph,1);       
Ea = zeros(nph,1);         
Eex = zeros(nph,1);
% UTD or GTD Solution
for i=1:nph
    if ph(i)<=ph0+pi
        Ui = 1;
    else
        Ui = 0;
    end
    if ph(i)<=pi-ph0
        Ur = 1;
    else
        Ur = 0;     
    end
    u0i = exp(1j*k*r*cos(ph(i) - ph0));
    u0r = exp(1j*k*r*cos(ph(i) + ph0));
    [DHU,DVU] = wdca(r,ph(i)*180/pi,php,btd,n,1);
    [DHK,DVK] = wdca(r,ph(i)*180/pi,php,btd,n,2);
    Edu(i)=(exp(-1j*k*r)/sqrt(r))*DHU + Ui*u0i - Ur*u0r;
    Edk(i)=(exp(-1j*k*r)/sqrt(r))*DHK + Ui*u0i - Ur*u0r;
    Edkd(i) = (exp(-1j*k*r)/sqrt(r))*DHK;
end
% Physical optics solution
for i=1:nph
    if ph(i)<=ph0+pi
        Ui = 1;
    else
        Ui = 0;
    end
    if ph(i)<=pi-ph0
        Ur = 1;
    else
        Ur = 0;     
    end
    u0i = exp(1j*k*r*cos(ph(i) - ph0));
    u0r = -exp(1j*k*r*cos(ph(i) + ph0));
    Epod(i) = (exp(-1j*pi/4)/sqrt(2*pi*k))*(exp(-1j*k*r)/sqrt(r))...
            *(sph0/(cph0+cos(ph(i))));
    Epoi(i) = u0i*Ui+u0r*Ur;
    Epo(i)=Epod(i)+Epoi(i);
    term = k*(L/2)*(cph0+cos(ph(i)));
    Epof(i) = k*L*sqrt(1/(2*pi*k))*sph0*exp(1j*pi/4)...
        *(exp(-1j*k*r)/sqrt(r))*exp(1j*term)*(sin(term)/term);
    Ea(i) = sqrt(2/(k*pi))*(sin(ph0/2)*sin(ph(i)/2)/(cph0+cos(ph(i))))...
        *exp(-1j*pi/4)*(exp(-1j*k*r)/sqrt(r));
end
% Exact Solution
for i=1:nph
    ai = sqrt(2)*cos(0.5*(ph(i) - ph0));
    ar = sqrt(2)*cos(0.5*(ph(i) + ph0));
    if ai>=0
        ei = 1;
        Ui = 1;
    else
        ei = -1;
        Ui = 0;
    end
    if ar>=0
        er = 1;
        Ur = 1;
    else
        er = -1;
        Ur = 0;
    end
    u0i = exp(1j*k*r*cos(ph(i) - ph0));
    u0r = exp(1j*k*r*cos(ph(i) + ph0));
    xi = abs(ai)*sqrt(k*r);
    xr = abs(ar)*sqrt(k*r);
    yi = sqrt(2/pi)*xi;
    yr = sqrt(2/pi)*xr;
    [FCi,FSi]= fresnel(yi);
    [FCr,FSr]= fresnel(yr);
    Kni = (sqrt(1j)/sqrt(2))*exp(1j*xi^2)...
        *((exp(-1j*pi/4)/sqrt(2))-(FCi-1j*FSi));
    Knr = (sqrt(1j)/sqrt(2))*exp(1j*xr^2)...
        *((exp(-1j*pi/4)/sqrt(2))-(FCr-1j*FSr));
    udi = -ei*Kni*exp(-1j*k*r);
    udr = -er*Knr*exp(-1j*k*r);
    Eex(i) = Ui*u0i-Ur*u0r+udi-udr-u0i+u0i;
end
figure(3);hold on;
plot(ph*180/pi,20*log10(abs(Eex)),'r','linewidth',2)
plot(ph*180/pi,20*log10(abs(Edu)),'g--','linewidth',2)
plot(ph*180/pi,20*log10(abs(Edk)),'b--','linewidth',2)
plot(ph*180/pi,20*log10(abs(Epo)),'k','linewidth',2)
xlim([0 360]);
ylim([-50 20]);
xlabel('\phi^o')
ylabel('|E_z| (dB)')
grid on;
box on
legend({'Exact','UTD','GTD','PO-Equation (46) in Solution'});
set(gca,'fontsize',24);

figure(2);hold on;
plot(ph*180/pi,20*log10(abs(Epod)),'linewidth',2)
plot(ph*180/pi,20*log10(abs(Edkd)),'linewidth',2)
xlim([0 360]);
ylim([-50 20]);
xlabel('\phi^o')
ylabel('|E_z| (dB)')
grid on;
box on
legend({'Diffracted Field - PO: Equation (44) Solution','Diffracted Field - GTD'})
set(gca,'fontsize',24);

figure(1);hold on;
plot(ph*180/pi,20*log10(abs(Epof)),'linewidth',2)
xlim([0 360]);
ylim([-50 20]);
xlabel('\phi^o')
ylabel('|E_z| (dB)')
grid on;
box on
set(gca,'fontsize',24);

keyboard


    
      

     

    