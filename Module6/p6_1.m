% Parameters
E0 = 1;
k = 2*pi;
h = 0.5;
w = 2;
a = atan2d(h,w/2);

RSB1 = a;
ISB1 = 360 - a;
RSB2 = 180 - a;
ISB2 = 180 + a;

phi = linspace(0,360,721);
rho = 1; % Set to 1 because we only care about the radiation pattern

%% Geometric Optics Solution
E_GO = zeros(size(phi));

for i = 1:length(phi)

    % Region II
    if (phi(i) < RSB1 || phi(i) > ISB1)
        E_GO(i) = E0*exp(1i*k*h*sind(phi(i)));
    % Region IV
    elseif (phi(i) > RSB2 && phi(i) < ISB2)
        E_GO(i) = E0*exp(1i*k*h*sind(phi(i)));
    % Region III
    elseif (phi(i) > RSB1 && phi(i) < RSB2)
        E_GO(i) = 2i*E0*sin(k*h*sind(phi(i)));
    end
end


%% GTD Solution
E_D1 = zeros(size(phi));
E_D2 = zeros(size(phi));

L1 = sqrt(h.^2 + w.^2/4);
L2 = L1;
s = L1;
phpd = a;
b0 = 90;
n = 2;
nu = 1; % GTD

for i = 1:length(phi)
    %Edge 1
    if phi(i) < 180
        phd = 180 - phi(i);
    else
        phd = 540 - phi(i);
    end
    [DH1, ~] = wdca(L1,phd,phpd,b0,n,nu);
    E_D1(i) = E0*exp(-1i*k*s)*DH1*exp(1i*k*(w/2)*cosd(phi(i)))/sqrt(s);

    %Edge 2
    phd = phi(i);
    [DH2, ~] = wdca(L2,phd,phpd,b0,n,nu);
    E_D2(i) = E0*exp(-1i*k*s)*DH2*exp(-1i*k*(w/2)*cosd(phi(i)))/sqrt(s);

end

E_UTD = E_D1 + E_D2 + E_GO;

% Plotting
figure;
polarplot(phi*pi/180,abs(E_GO),'r','LineWidth',2,'DisplayName','|E_{GO}|')
hold on
polarplot(phi*pi/180,abs(E_UTD),'b','LineWidth',2,'DisplayName','|E_{UTD}|')
title('Plots of GO and Diffracted Fields')
legend()
saveas(gcf,'p6_1','png')

E_GO = E_GO+min(abs(E_UTD));
E_UTD_min = min(20*log10(abs(E_UTD)));
E_GO_min = min(20*log10(abs(E_GO)));

figure;
polarplot(phi*pi/180,20*log10(abs(E_GO)) - E_GO_min,'r','LineWidth',2,'DisplayName','|E_{GO}|')
hold on
polarplot(phi*pi/180,20*log10(abs(E_UTD)) - E_UTD_min,'b','LineWidth',2,'DisplayName','|E_{UTD}|')
title('Plots of GO and Diffracted Fields - Normalized to UTD')
legend()
saveas(gcf,'p6_1_norm','png')