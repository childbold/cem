% Paramters
clear
close all

E0 = 1;
a = 10;
B0 = 90;
n = 2;
v = 1;

th0 = 4;
theta = th0:0.01:180-th0;

% Geometric Optics Field
E_GO = E0*cosd(90*cosd(theta))./sind(theta);
E_GO(theta >= 90 & theta <= 270) = 0;

% Diffracted Field
E_D1 = zeros(size(theta));
E_D2 = zeros(size(theta));

for i = 1:length(theta)

    % Edge 1
    phd = mod(90 + theta(i),360);

    [DH1, DV1] = wdca(a,phd,0,B0,n,v);
    E_D1(i) = E0*DV1*exp(-1i*2*pi*a)*exp(1i*2*pi*a*sind(theta(i)))./sqrt(a*sind(theta(i)));

    % Edge 2
    phd = mod(90 - theta(i),360);
 
    [DH2, DV2] = wdca(a,phd,0,B0,n,v);
    E_D2(i) = E0*DV2*exp(-1i*2*pi*a)*exp(-1i*2*pi*a*sind(theta(i)))./sqrt(-a*sind(theta(i)));

end

E_T = E_GO + E_D1 + E_D2;

% Plotting
f1=figure;
set(gcf,'Position',[387 565 560 420])
ax1 = polaraxes;
polarplot(ax1,theta*pi/180,abs(E_GO),'r','LineWidth',2)
ax1.ThetaZeroLocation = 'top';
hold on
polarplot(ax1,theta*pi/180,abs(E_T),'b','LineWidth',2)

E_GO = E_GO+min(abs(E_T));
E_T_min = min(20*log10(abs(E_T)));
E_GO_min = min(20*log10(abs(E_GO)));

f2=figure;
set(gcf,'Position',[968 565 560 420])
ax2 = polaraxes;
polarplot(ax2,theta*pi/180, 20*log10(abs(E_T))-E_T_min,'b','linewidth',2)
ax2.ThetaZeroLocation = 'top';
hold on;
polarplot(ax2,theta*pi/180, 20*log10(abs(E_GO))-E_GO_min,'r','linewidth',2)

EGO2 = flip(E_GO);
th2 = 180+th0:0.01:360-th0;

polarplot(ax1,th2*pi/180, abs(EGO2),'r','linewidth',2)
Et2 = flip(E_T);
polarplot(ax1,th2*pi/180, abs(Et2),'b','linewidth',2)
title(ax1,sprintf('a = %d\\lambda',a))
saveas(f1,sprintf('fields_1_a_%d',a),'png')

polarplot(ax2,th2*pi/180, 20*log10(abs(Et2))-E_T_min,'b','linewidth',2)
polarplot(ax2,th2*pi/180, 20*log10(abs(EGO2))-E_GO_min,'r','linewidth',2)
title(sprintf('a = %d\\lambda',a))
saveas(f2,sprintf('fields_2_a_%d',a),'png')
