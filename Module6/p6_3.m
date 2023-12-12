% Parameters

E0 = 1;
n = 1.5;
k = 2*pi;
rho = 50;

boundary = (2 - n)*180;

phi = linspace(0,270,1001);

E_incident = E0*exp(-1i*k*rho*cosd(phi - n*180/2));

RSB1 = 180 - n*180/2;
RSB2 = n*180 - RSB1;

%% GO Field
E_GO = zeros(size(phi));

for i = 1:length(phi)

    if phi(i) < RSB1
        E_GO(i) = -E0*exp(1i*k*rho*cosd(phi(i) + n*180/2));
    elseif phi(i) > RSB2
        E_GO(i) = -E0*exp(-1i*k*rho*cosd(phi(i) + n*180/2));
    end
end

%% Diffracted Field
L = rho;
phpd = n*180/2;
b0 = 90;
nu = 1;

E_D_h = zeros(size(phi));
E_D_v = zeros(size(phi));
for i = 1:length(phi)

    phd = phi(i);
    [DH,DV] = wdca(L,phd,phpd,b0,n,nu);
    E_D_h(i) = DH*exp(-1i*k*L)/sqrt(L);
    E_D_v(i) = -DV*exp(-1i*k*L)/sqrt(L);
end

%% Plotting
figure; hold on
plot(phi,20*log10(abs(E_GO+E_D_h)),'b','LineWidth',2)
plot(phi,20*log10(abs(E_GO+E_D_v)),'r','LineWidth',2)