%Input
%nm: number of Mie (modal) series modes to be summed

nm = 50;

close all

abwl = linspace(0,4,4001);

RCS_exact = zeros(length(abwl),1);

RCS_po = zeros(length(abwl),1);

for i =1:length(abwl)
    ka = 2*pi*abwl(i);
    sum = 0;
    for n=1:nm
        H = ka*sqrt(pi/(2*ka))*besselh(n+0.5,2,ka);
        Hd = (n+1).*sqrt(pi/(2*ka))*besselh(n+0.5,2,ka)...
            -ka*sqrt(pi/(2*ka))*besselh(n+1.5,2,ka);
        sum = sum + (-1)^n*(2*n+1)/(H*Hd);
    end
    RCS_exact(i) = (1/(4*(pi^2)))*(1/abwl(i)^2)*abs(sum)^2;
    RCS_po(i) = 1-2*cos(ka)*(sin(ka)/ka)+(sin(ka)/ka)^2;
    
end
figure(1);
hold on;
plot(abwl,10*log10(RCS_exact),'linewidth',3)
plot(abwl,10*log10(RCS_po),'linewidth',3)
grid on
box on
xlabel('a/\lambda')
ylabel('RCS/\pia^2 (dB)')
% set(gca,'fontsize',36) 

