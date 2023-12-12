f = 3e9;
c = 3e8;
lambda = c/f;
k = 2*pi/lambda;
N = 100;

a_over_lambda = linspace(0,4,4001);
rcs_bs = zeros(length(a_over_lambda),1);

for i = 1:length(a_over_lambda)
ka = 2*pi*a_over_lambda(i);

summation = zeros(N,1);
bess1 = zeros(N,1);
bess2 = zeros(N,1);

for n = 1:N
    bess1(n) = ka*sqrt(pi/(2*ka))*besselh(n+0.5,2,ka);
    bess2(n) = (n+1).*sqrt(pi/(2*ka))*besselh(n+0.5,2,ka)-ka*sqrt(pi/(2*ka))*besselh(n+1.5,2,ka);
    summation(n) = (-1).^n*(2*n+1)/(bess1(n)*bess2(n));
end

rcs_bs(i) = lambda.^2/(4*pi)*abs(sum(summation)).^2;

end

figure;
plot(a_over_lambda,10*log10(rcs_bs))

