%% Script to plot the radiation pattern of a current dipole in the prescence of a PEC sphere
close all; clear all
load('../constants.mat')
lambda = c*2*pi/(omega);

I = 1;
L = 1;
r = 1;

a_vec = [0.05, 0.25, 2, 10, 20];
% a_vec = [2, 10, 75, 200, 500];

theta_vec = linspace(0,2*pi,361);
power = zeros(numel(theta_vec),numel(a_vec));

N = 300;
a_idx = 1;

for a = a_vec

    coeff = (I*L)^2./(32*pi^2*r^2*k0^2*(a*lambda)^4);
    coeff = 1;
   
    p_idx = 1;
    for theta = theta_vec

        sum = 0;

        for n = 1:N
    
            temp = (1i^n)*(2*n+1);
            
            Pnm = legendre(n, cos(theta));
            
            P1n = Pnm(2, :);
            
            hank = sqrt(pi*k0*a*lambda/8)*(-besselh(n+1.5,2,k0*a*lambda) + ((n+1)/(k0*a*lambda))*besselh(n+.5,2,k0*a*lambda));

            if isnan(hank)
                invH = 0;
            else
                invH = 1/hank;
            end
            sum = sum+temp*P1n*invH;

        end

        power(p_idx,a_idx) = abs(sum).^2;
        p_idx = p_idx+1;

    end

    power(:,a_idx) = power(:,a_idx)*coeff;
    power(:,a_idx) = power(:,a_idx)./max(power(:,a_idx));
    a_idx = a_idx+1;
        
end

ax = polaraxes;
ax.ThetaZeroLocation = 'top';
hold on; set(gcf,'Position',[-1208 280 752 559]);
polarplot(ax,theta_vec, power(:,1),'b','LineWidth',2,'DisplayName','a = 0.05\lambda')
polarplot(ax,theta_vec, power(:,2),'r','LineWidth',2,'DisplayName','a = 0.25\lambda')
polarplot(ax,theta_vec, power(:,3),'g','LineWidth',2,'DisplayName','a = 2\lambda')
polarplot(ax,theta_vec, power(:,4),'m','LineWidth',2,'DisplayName','a = 10\lambda')
polarplot(ax,theta_vec, power(:,5),'k','LineWidth',2,'DisplayName','a = 20\lambda')
legend()

