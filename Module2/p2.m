% Initial parameters
theta_0 = 90;
phi_vec = linspace(0,2*pi,361);
E_0 = 1;
freqs = [1,5,10,20];
a = 0.25;
eta = 120*pi;    %377;
n_vec = 0:200;

currents = zeros(numel(phi_vec),numel(freqs));
currents_po = zeros(size(phi_vec));

figure; grid on; hold on;
title('Surface Currents from a PEC Cylinder')
for freq_idx = 1:length(freqs)

    k = 2*pi*freqs(freq_idx)*1e9./3e8; %Wavenumber
    coeff = 2*E_0/(pi*k*a*eta);

    for phi_idx = 1:length(phi_vec)
        
           summation = 1/besselh(0,2,k*a);
        for n_idx = 2:length(n_vec)
    
            cosine = 2*cos(n_vec(n_idx)*phi_vec(phi_idx));
            hank = besselh(n_vec(n_idx),2,k*a);
    
            summation = summation + 1i^(-n_vec(n_idx))*cosine/hank;
    
        end

       currents(phi_idx,freq_idx) = coeff*summation;
       currents_po(phi_idx) = (-2/eta)*cos(phi_vec(phi_idx))*E_0*exp(-1i*k*a*cos(phi_vec(phi_idx)));

       currents_po(phi_vec < pi/2 | phi_vec > 3*pi/2) = 0;

    end

    plot(phi_vec*180/pi,abs(currents(:,freq_idx)),'DisplayName',sprintf('%d GHz',freqs(freq_idx)),'LineWidth',2)
    xlim([0 360])
    xlabel('\phi')
    ylabel('|J_z| A/m')
end
saveas(gcf,'p3_exact','png')

plot(phi_vec*180/pi,abs(currents_po),'--','DisplayName','PO Currents','LineWidth',2)
legend()

saveas(gcf,'p3_po','png')
        
