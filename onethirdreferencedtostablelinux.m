clear

t = 0:0.001:10;
nt = length(t);

% thetaonethird = [1/3, 1, 3 * sqrt(3) / 8, 0];
% stablepdfonethird = stablepdf(t, thetaonethird, 1);

onethirddist = makedist('Stable', 'alpha', 1/3, 'beta', 1, 'gam', 3 * sqrt(3) / 8, 'delta', 3/8);
pdfonethird = pdf(onethirddist, t);

twentysixb = heaviside(t) / 3^(1/3) ./ t.^(4/3) .* airy(1 ./ (3 * t).^(1/3));

figure(1)
plot(t, pdfonethird, '-',  t, twentysixb, '--',   'linewidth', 2)
legend('matlab stable y = 1/3', 'airy y = 1/3')
title('y = 1/3, referenced to PLWE paper')

scale = 1;
% thetaonethirdscaleone = [1/3, 1, scale, 0];
% stablepdfonethirdscaleone = stablepdf(t, thetaonethirdscaleone, 1);
a = 3 * sqrt(3) / 8 / scale;
twentysixbscaled = a * heaviside(t) / 3^(1/3) ./ (a * t).^(4/3) .* airy(1 ./ (a * 3 * t).^(1/3));

approxairy = a * heaviside(t) .* (3 * a * t).^(-1/4) / (2 * sqrt(pi)) .* exp(- 2 * (a * t) .* sqrt(a * t) / (3 * sqrt(3)));
besselKexact =  a * heaviside(t) .* (a * t).^(-3/2) / (3 * pi) .* besselk(1/3, 2 ./ (3 * sqrt(3 * (a * t))));

figure(2)
plot( t, twentysixbscaled, '--', ...
    t, real(besselKexact), ':', t, real(approxairy), ':',  'linewidth', 2)
legend('airy y = 1/3', 'bessel k', 'approx airy')
title('y = 1/3, referenced to stable')
v = axis;
v(3) = 0;
v(4) = 1.2;
axis(v)
grid

figure(3)
airyfunction = a * heaviside(t) ./ (3 * (a * t).^4).^(1/3) .* airy(1 ./ (3 * (a * t)).^(1/3));
airyJFKapprox = a * heaviside(t) ./ (3 * (a * t).^4).^(1/3) .* ...
    (3 ./ (3 * (a * t)).^(1/3)) .^ (-1/4) / (2 * sqrt(pi)) .* exp(2 * a * t .* sqrt(a * t) / 3 / sqrt(3));
% airy(1 ./ (3 * (a * t)).^(1/3));
besselKexact =  a * heaviside(t) .* (a * t).^(-3/2) / (3 * pi) .* besselk(1/3, 2 ./ (3 * sqrt(3 * a * t)));
scale = 1;
% thetaonethirdscaleone = [1/3, 1, 1, 0];
% stablepdfonethirdscaleone = stablepdf(t, thetaonethirdscaleone, 1);
plot(t, besselKexact, '--', t, airyfunction, ':', t, airyJFKapprox, ':', 'linewidth', 2)
legend('bessel k', 'airy function', 'airy approx')
title('y = 1/3, \beta = 1')
v = axis;
v(1) = -0.1;
v(2) = 1;
v(3) = 0;
v(4) = 2;
axis(v)
grid
set(gca, 'fontsize', 16)

figure(4)
airysmall = t.^(-3/2) * gamma(1/3) .* t.^(1/6)/ sqrt(3 * pi) / 2;
airyapprox = heaviside(t) .* (3 * t).^(-1/4) / (2 * sqrt(pi)) .*...
    exp(- 2 * t .* sqrt(t) / (3 * sqrt(3)));
plot(t, airyapprox, '-', t, airy(t), '--', t, airysmall, '--', 'linewidth', 2)
legend('approx', 'airy', 'airysmall')
title('y = 1/3, \beta = 1')
v = axis;
v(1) = -0.1;
v(2) = 1;
v(3) = 0;
v(4) = 2;
axis(v)
grid

scale = 1;
a = 3 * sqrt(3) / 8 / scale;
figure(5)
airyfunction = a * heaviside(t) ./ (3 * (a * t).^4).^(1/3) .* airy(1 ./ (3 * (a * t)).^(1/3));
% airyjmj = a * heaviside(t) ./ (3 * (a * t).^4).^(1/3) .* airyjianmingjin(1 ./ (3 * (a * t)).^(1/3));
airyjmj = a * heaviside(t) ./ (3 * (a * t).^4).^(1/3) .* airyjmjcode(1 ./ (3 * (a * t)).^(1/3));

% airyJFKapprox = a * heaviside(t) ./ (3 * (a * t).^4).^(1/3) .* ...
%     (3 ./ (3 * (a * t)).^(1/3)) .^ (-1/4) / (2 * sqrt(pi)) .* exp(2 * a * t .* sqrt(a * t) / 3 / sqrt(3));
% airy(1 ./ (3 * (a * t)).^(1/3));
besselKexact =  a * heaviside(t) .* (a * t).^(-3/2) / (3 * pi) .* ...
    besselk(1/3, 2 ./ (3 * sqrt(3 * a * t)));
% thetaonethirdscaleone = [1/3, 1, scale, 0];
% stablepdfonethirdscaleone = stablepdf(t, thetaonethirdscaleone, 1);
plot( t, besselKexact, '-', t, airyfunction, '--',  t, airyjmj, ':', 'linewidth', 2)
legend('bessel k', 'airy function', 'airy jmj')
title('y = 1/3, \beta = 1')
v = axis;
v(1) = -0.1;
v(2) = 1;
v(3) = 0;
v(4) = 1.4;
axis(v)
grid
set(gca, 'fontsize', 16)

figure(6)
airyfunction = a * heaviside(t) ./ (3 * (a * t).^4).^(1/3) .* airy(1 ./ (3 * (a * t)).^(1/3));
airyjmj = a * heaviside(t) ./ (3 * (a * t).^4).^(1/3) .* airyjianmingjin(1 ./ (3 * (a * t)).^(1/3));
besselKexact =  a * heaviside(t) .* (a * t).^(-3/2) / (3 * pi) .* ...
    besselk(1/3, 2 ./ (3 * sqrt(3 * a * t)));
semilogy( t, abs(airyfunction - besselKexact), '-', t, ... 
    abs(airyfunction -  airyjmj), '--', 'linewidth', 2)
legend('bessel k',  'airy jmj')
title('y = 1/3, \beta = 1')
v = axis;
v(1) = 0;
v(2) = 2;
% v(3) = 0;
% v(4) = 1.5;
axis(v)
grid
set(gca, 'fontsize', 16)

