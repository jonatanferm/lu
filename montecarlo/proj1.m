clear
clc
load powercurve_V90
%%
n = 1000;
kl=[10.6, 9.7,9.2, 8.0,7.8,8.1,7.8,8.1,9.1,9.9,10.6,10.6;
    2.0,2.0,2.0,1.9,1.9,1.9,1.9,1.9,2.0,1.9,2.0,2.0];
m = 1;
colors = [
    241,163,64
153,142,195
253,141,60
230,85,13
166,54,3]./255;

%% Values CRUDE MC
N = 1000000;
means = zeros(1, size(kl, 2));
cis = zeros(1, size(kl, 2));
for cm = 1:size(kl, 2)
    wind_speeds = wblrnd(kl(1, cm), kl(2, cm), N, 1);
    %wind_speeds = wbl_trunc(kl(1, cm), kl(2, cm), N);
    powers = P(wind_speeds);
    clear wind_speeds
    E_P = mean(powers); %.*(wblcdf(25,kl(1,cm),kl(2,cm)) - wblcdf(4,kl(1,cm),kl(2,cm)));
    STD_P = std(powers);
    CI = norminv(0.975)*STD_P./sqrt(N);
    means(cm) = E_P;
    cis(cm) = CI;
end
%% IS
nn = 10000;
xrange = linspace(4, 25, nn);
options = optimset('MaxFunEvals', 100000, 'MaxIter', 100000);
fphi = P(xrange)'.*wblpdf(xrange, mean(kl(1,:), 2), mean(kl(2,:), 2));
fphi = fphi .* (nn/(sum(fphi) * 25-4));
[sm, v2] = fminsearch(@(x) var(fphi ./ gampdf_trunc(xrange, x(1), x(2))), [8.3007 1.4869], options);
%[sm2, v22] = fminsearch(@(x) var(fphi ./ normpdf_trunc(xrange, x(1), x(2))), [11 4], options);es
N = 1000000;
means = zeros(1, size(kl, 2));
cis = zeros(1, size(kl, 2));

for cm = 1:size(kl, 2)
    cm
    X = gamrnd_trunc(sm(1), sm(2), N);
    cm
    omega = (wblpdf(X, kl(1,cm), kl(2,cm))./gampdf_trunc(X, sm(1), sm(2)))';
    tau = P(X).*omega;
    E_P = mean(tau);
    STD_P = std(tau);
    CI = norminv(0.975)*STD_P./sqrt(N);
    means(cm) = E_P;
    cis(cm) = CI;
end
%%
X = gamrnd_trunc(sm(1), sm(2), n);
omega = (wblpdf(X, kl(1,m), kl(2,m))./gampdf_trunc(X, sm(1), sm(2)))';

tau = cumsum(P(X).*omega)./(1:n)';

phiomega = P(X).*omega;
phiomega_mean = cumsum(phiomega)./(1:numel(phiomega))';
tau_std = sqrt((sum(triu(phiomega(2:end)-phiomega_mean(2:end)').^2)./(2:n)));

CMC = cumsum(P(wblrnd(kl(1, m), kl(2, m), n, 1)))./(1:n)';


CI = norminv(0.975)*tau_std./sqrt(2:n);
mean_ci = [tau(2:end) + CI', tau(2:end) - CI'];
figure;
hold on;
plot(mean_ci, 'r--');
plot(tau(2:end), 'b-');
plot(CMC, 'k--')
hold off;
tau_std(end)
%% Antithetic sampling values
N = 500000;
means = zeros(1, size(kl, 2));
cis = zeros(1, size(kl, 2));
for cm = 1:size(kl, 2)
    U = rand(1, N);
    ws1 = wblinv(U, kl(1, cm), kl(2, cm));
    ws2 = wblinv(1-U, kl(1, cm), kl(2, cm));
    powers = (P(ws1) + P(ws2))/2;
    E_P = mean(powers);
    STD_P = std(powers);
    CI = norminv(0.975)*STD_P./sqrt(N);
    means(cm) = E_P;
    cis(cm) = CI;
end
%% sampling values TRUNC
N = 1e7;
EP = zeros(1, size(kl, 2));
EP_cis = zeros(1, size(kl, 2));
for cm = 1:size(kl, 2)
    [ws1, ws2] = wbl_trunc_antithetic(kl(1,cm), kl(2,cm), N);
    powers = (P(ws1) + P(ws2))/2;
    powers = powers.*(wblcdf(25,kl(1,cm),kl(2,cm)) - wblcdf(4,kl(1,cm),kl(2,cm)));
    E_P = mean(powers);
    STD_P = std(powers);
    CI = norminv(0.975)*STD_P./sqrt(N);
    EP(cm) = E_P;
    EP_cis(cm) = CI;
end
%% plots
%% Antithetic sampling truncated
[ws1, ws2] = wbl_trunc_antithetic(kl(1,m), kl(2,m), n/2);
powers = (P(ws1) + P(ws2))/2;
powers_mean = cumsum(powers)./(1:numel(powers))' ...
    .*(wblcdf(25,kl(1,m),kl(2,m)) - wblcdf(4,kl(1,m),kl(2,m)));
powers_std = sqrt(sum(triu(powers(2:end)-powers_mean(2:end)').^2)./(1:n/2-1));
CI = norminv(0.975)*powers_std./sqrt(2:n/2);
mean_ci = [powers_mean(2:end) + CI', powers_mean(2:end) - CI'];
figure;
hold on;
plot(mean_ci, 'r--');
plot(powers_mean, 'b-');
%plot(mean_ci(:,1) - mean_ci(:,2), 'g-')
hold off;
powers_std(end)
EP = powers_mean(end);
%% IS anti trunc
N = 500000;
means = zeros(1, size(kl, 2));
cis = zeros(1, size(kl, 2));
for cm = 1:size(kl, 2)
    [x1, x2] = gamrnd_trunc_antitethic(sm(1), sm(2), N);
    omega1 = (wblpdf(x1, kl(1,cm), kl(2,cm))./gampdf_trunc(x1, sm(1), sm(2)))';
    omega2 = (wblpdf(x2, kl(1,cm), kl(2,cm))./gampdf_trunc(x2, sm(1), sm(2)))';
    powers = (P(x1).*omega1 + P(x2).*omega2)/2;
    %powers = powers.*(wblcdf(25,kl(1,cm),kl(2,cm)) - wblcdf(4,kl(1,cm),kl(2,cm)));
    E_P = mean(powers);
    STD_P = std(powers);
    CI = norminv(0.975)*STD_P./sqrt(N);
    means(cm) = E_P;
    cis(cm) = CI;
end
%% Probabilities that power will be generated
wblcdf(25,kl(1,:), kl(2,:)) - wblcdf(4,kl(1,:), kl(2,:))
%% Explicit PTOT values
explicit_P_tot = [
         6.16928436328976800555221540683e6   % jan
         4.72750934302877747648340332285e6   % feb
         4.03348712551416846516540794645e6   % mar
         2.80770364105932821998696428878e6   % apr
         2.60234644193747329189698022885e6   % may
         2.91431412247697353234002400118e6   % jun
         2.60234644193747329189698022885e6   % jul
         2.91431412247697353234002400118e6   % aug
         3.90338483020264527572424595168e6   % sep
         5.32092194378168967290064738755e6   % oct
         6.16928436328976800555221540683e6    % nov
         6.16928436328976800555221540683e6    % dec
    ];

%% P TOT ESTIMATE
N = 1e7;
P_tot = @(v) .5 * 1.225 * pi * 90^2 * .25 * v.^3;

eptots = nan(1, 12);
cis = nan(1, 12);
for cm = 1:12
    [ws1, ws2] = wbl_antithetic(kl(1,cm), kl(2,cm), N);
    tau = (P_tot(ws1) + P_tot(ws2))/2;
    E_P_tot = mean(tau);
    STD_P_tot = std(tau);
    CI = norminv(0.975)*STD_P_tot./sqrt(N);
    eptots(cm) = E_P_tot;
    eptots_cis(cm) = CI;
end
%%
plot(xrange, fphi ./ gampdf_trunc(xrange, sm(1), sm(2)))
%%
N = 1e6;
for cm = 1:size(kl, 2)
    cm
    X = gamrnd_trunc(sm(1), sm(2), N);
    omega = (wblpdf(X, kl(1,cm), kl(2,cm))./gampdf_trunc(X, sm(1), sm(2)))';
    tau = phi(X)'.*omega;
    E_P = mean(tau);
    STD_P = std(tau);
    CI = norminv(0.975)*STD_P./sqrt(N);
    means(cm) = E_P;
    cis(cm) = CI;
end
%% antithetic 3
k = 1.96;
lamb = 9.13;
U = rand(1, 10000000);
ws1 = wblinv(U, lamb, k);
ws2 = wblinv(1-U, lamb, k);
tau = (P(ws1) + P(ws2))/2;
EP = mean(tau);
STDP = std(tau);
CI = norminv(0.975)*STDP./sqrt(10000000);

%%
[v1, v2] = meshgrid(linspace(4, 25), linspace(4, 25));
options = optimset('MaxFunEvals', 30000, 'MaxIter', 30000);
fphi = (reshape(P(v1), size(v1, 1), size(v1, 2)).* ...
    reshape(P(v2), size(v1, 1), size(v1, 2))).*biwblpdf(v1, v2);
fphi = fphi .* (size(v2, 1)/(sum(fphi, 'all') * 30-0));
[sm, vv2] = fminsearch(@(x) var( fphi ./ reshape(mvnpdf([v1(:) v2(:)], x(1), get_sigma(x(2), x(3))), size(v1)), 0, 'all'),...
    [12 10 -10], options);
%[sm, vv2] = fminsearch(@(x) var( fphi ./ reshape(mvtpdf([v1(:)+x(4) v2(:)+x(4)], get_sigma(x(2), x(3)), abs(x(1))), size(v1)), 0, 'all'),...
%    [12 10 -10 -12], options);
%[sm, vv2] = fminsearch(@(x) norm( reshape(fphi, [], 1) - mvnpdf([v1(:) v2(:)], x(1), get_sigma(x(2),x(2)))),...
%    [12 1 0], options);
%[sm, vv2] = fminsearch(@(x) var( reshape(fphi, [], 1) ./ gmpdf([v1(:) v2(:)], x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8))),...
%    [11.1765, -16.3261, 35.7846, 3.7291, 6.3852, 0.5215, 1.7827, 0.9084], options);

%%
Z = mvnpdf([v1(:) v2(:)], sm(1), get_sigma(sm(2), sm(3)));
%%
Z = reshape(Z, size(v1));
mesh(v1, v2, fphi./Z);
%% bivar IS values
N = 10000000;
X = mvnrnd([sm(1) sm(1)], [sm(2) sm(3); sm(3) sm(2)], N);
omega = biwblpdf(X(:,1), X(:,2))./mvnpdf(X, sm(1), get_sigma(sm(2), sm(3)));
tau = (P(X(:,1)).*P(X(:,2))).*omega;
E_P = mean(tau);
STD_P = std(tau);
CI = norminv(0.975)*STD_P./sqrt(N);
EPP = [E_P+CI, E_P-CI];
%%

X = mvnrnd([sm(1) sm(1)], [sm(2) sm(3); sm(3) sm(2)], n);

omega = biwblpdf(X(:,1), X(:,2))./mvnpdf(X, sm(1), [sm(2) sm(3); sm(3) sm(2)]);
PX = [P(X(:,1)), P(X(:,2))];

covar = mean((prod(PX, 2).*omega)) - EP.^2
tau_std = (1/sqrt(n))*std(prod(PX, 2).*omega);


% integrate (.5*1.225*pi*(90^2)/4)*(((1.96)/(9.13))*(x/9.13)^(1.96-1))*e^(-((x/9.13)^1.96))^3 from 0 to inf

phiomega = (P(X(:,1))+P(X(:,2))).*omega;
phiomega_mean = cumsum(phiomega)./(1:numel(phiomega))';
tau_std = sqrt((sum(triu(phiomega(2:end)-phiomega_mean(2:end)').^2)./(2:n)));

PG3 = (sum(PX, 2) > 3*10^6).*omega;
PL3 = (sum(PX, 2) < 3*10^6).*omega;

PG3_CI = [mean(PG3) - std(PG3)/sqrt(n), mean(PG3) + std(PG3)/sqrt(n)]
PL3_CI = [mean(PL3) - std(PL3)/sqrt(n), mean(PL3) + std(PL3)/sqrt(n)]
%%

%CMC = cumsum(P(wblrnd(kl(1, m), kl(2, m), n, 1)))./(1:n)';

CI = norminv(0.975)*tau_std./sqrt(2:n);
mean_ci = [tau(2:end) + CI', tau(2:end) - CI'];
figure;
hold on;
plot(mean_ci, 'r--');
plot(tau(2:end), 'b-');
%plot(CMC, 'k--')
hold off;
mean_ci(end,1) - mean_ci(end,2)
%% antithetic VAR
k = 1.96;
lamb = 9.13;
U = rand(1, 10000000);
ws1 = wblinv(U, lamb, k);
ws2 = wblinv(1-U, lamb, k);
tau = (P(ws1).*P(ws1) + P(ws2).*P(ws2))/2;
EP = mean(tau);
STDP = std(tau);
CI = norminv(0.975)*STDP./sqrt(10000000);

%% d leq
[v1, v2] = meshgrid(linspace(0, 17), linspace(0, 17));
options = optimset('MaxFunEvals', 30000, 'MaxIter', 30000);
fphi = (reshape(P(v1), size(v1, 1), size(v1, 2)) +  ...
    reshape(P(v2), size(v1, 1), size(v1, 2)) < 3e6).*biwblpdf(v1, v2);

fphi = fphi .* (size(v2, 1)/(sum(fphi, 'all')));
%[sm, vv2] = fminsearch(@(x) norm( fphi - reshape(mvnpdf([v1(:) v2(:)], x(1), get_sigma(x(2), x(3))), size(v1))),...
%    [12 10 -10], options);
[sm, vv2] = fminsearch(@(x) var( fphi ./ reshape(mvnpdf([v1(:) v2(:)], x(1), get_sigma(x(2), x(3))), size(v1)), 0, 'all'),...
    [12 10 -10], options);
N = 10000000;
X = mvnrnd([sm(1) sm(1)], get_sigma(sm(2), sm(3)), N);
omega = biwblpdf(X(:,1), X(:,2))./mvnpdf(X, sm(1), get_sigma(sm(2), sm(3)));
tau = (P(X(:,1))+P(X(:,2)) < 3e6).*omega;
PL3 = mean(tau);
STDPL3 = std(tau);
CIPL3 = norminv(0.975)*STDPL3./sqrt(N);
PL3_CI = [PL3+CIPL3, PL3-CIPL3];
%% d geq
[v1, v2] = meshgrid(linspace(10, 27), linspace(10, 27));
options = optimset('MaxFunEvals', 30000, 'MaxIter', 30000);
fphi = (reshape(P(v1), size(v1, 1), size(v1, 2)) +  ...
    reshape(P(v2), size(v1, 1), size(v1, 2)) > 3e6).*biwblpdf(v1, v2);
fphi = fphi .* (size(v2, 1)/(sum(fphi, 'all')));
[sm, vv2] = fminsearch(@(x) var( fphi ./ reshape(mvnpdf([v1(:) v2(:)], x(1), get_sigma(x(2), x(3))), size(v1)), 0, 'all'),...
    [12 10 -10], options);
N = 10000000;
X = mvnrnd([sm(1) sm(1)], get_sigma(sm(2), sm(3)), N);
omega = biwblpdf(X(:,1), X(:,2))./mvnpdf(X, sm(1), get_sigma(sm(2), sm(3)));
tau = (P(X(:,1))+P(X(:,2)) > 3e6).*omega;
PG3 = mean(tau);
STDPG3 = std(tau);
CIPG3 = norminv(0.975)*STDPG3./sqrt(N);
PG3_CI = [PG3+CIPG3, PG3-CIPG3];
%% d eq
[v1, v2] = meshgrid(linspace(0, 25), linspace(0, 25));
options = optimset('MaxFunEvals', 30000, 'MaxIter', 30000);
fphi = (reshape(P(v1), size(v1, 1), size(v1, 2)) +  ...
    reshape(P(v2), size(v1, 1), size(v1, 2)) == 3e6).*biwblpdf(v1, v2);
fphi = fphi .* (size(v2, 1)/(sum(fphi, 'all')));
[sm, vv2] = fminsearch(@(x) var( fphi ./ reshape(mvnpdf([v1(:) v2(:)], x(1), get_sigma(x(2), x(3))), size(v1)), 0, 'all'),...
    [12 10 -10], options);
N = 10000000;
X = mvnrnd([sm(1) sm(1)], get_sigma(sm(2), sm(3)), N);
omega = biwblpdf(X(:,1), X(:,2))./mvnpdf(X, sm(1), get_sigma(sm(2), sm(3)));
tau = (P(X(:,1))+P(X(:,2)) == 3e6).*omega;
Peq3 = mean(tau);
STDPeq3 = std(tau);
CIPeq3 = norminv(0.975)*STDPeq3./sqrt(N);
Peq3_CI = [Peq3+CIPeq3, Peq3-CIPeq3];
%%
PG3_CI + PL3_CI + Peq3_CI + (1-biwblcdf(25, 25))
%%
function sigma = get_sigma(s, c0)
    c = s*(atan(c0)/(.5001*pi));
    sigma = [s+1e-8 c; c s+1e-8];
end
function gm = get_gm(mu1, mu21, mu22, s1, s2, c1, c2, p0)
    mu = [mu1 mu1; mu21 mu22; mu22 mu21];
    s12 = s1.^2;
    s22 = s2.^2;
    cp1 = s12.*atan(c1)/(.5001*pi);
    cp2 = s22.*atan(c2)/(.5001*pi);
    sigma = cat(3,[s12 cp1; cp1 s12],...
                  [s22 cp2; cp2 s22],...
                  [s22 cp2; cp2 s22]);
    p = [p0.^2 1 1];
    p = p/sum(p);
    gm = gmdistribution(mu,sigma,p);
end
function Y = gmpdf(x, mu1, mu21, mu22, s1, s2, c1, c2, p0)
    gm = get_gm(mu1, mu21, mu22, s1, s2, c1, c2, p0); 
    Y = pdf(gm,x);
end
function P = biwblcdf(v1, v2)
    bk = 1.96;
    kl = 9.13;
    alpha = 0.638;
    p = 3;
    q = 1.5;
    F1 = wblcdf(v1, kl, bk);
    F2 = wblcdf(v2, kl, bk);
    P = F1.*F2.*(1 + alpha .* ((1-F1.^p).^q) .* ((1-F2.^p).^q));
end
function Y = biwblpdf(v1, v2)
    bk = 1.96;
    kl = 9.13;
    alpha = 0.638;
    p = 3;
    q = 1.5;
    f1 = wblpdf(v1, kl, bk);
    f2 = wblpdf(v2, kl, bk);
    F1 = wblcdf(v1, kl, bk);
    F2 = wblcdf(v2, kl, bk);
    Y = f1.*f2.*(1 + alpha .* ((1-F1.^p).^(q-1)) .* ((1-F2.^p).^(q-1)) .* ...
        ((F1.^p).*(1+p.*q)-1) .* ((F2.^p).*(1+p.*q)-1));
end
function x = normrnd_trunc(mu, sigma, n)
    u = rand(1, n);
    x = norminv(...
            u.*(normcdf(25, mu, sigma)-normcdf(4, mu, sigma)) + normcdf(4, mu, sigma),...
            mu, sigma);
end
function D = normpdf_trunc(x, mu, sigma)
    D = normpdf(x, mu, sigma) ./ ...
        (normcdf(25, mu, sigma)-normcdf(4, mu, sigma));
    D(x>25 | x<4) = 0;
end
function D = gampdf_trunc(x, A, B)
    D = abs((gampdf(x, A, B)) ./ ...
        (gamcdf(25, A, B)-gamcdf(4, A, B)));
    D(x>25 | x<4) = 0;
end
function x = gamrnd_trunc(A, B, n)
    u = rand(1, n);
    gcdf4 = gamcdf(4, A, B);
    c = (gamcdf(25, A, B)-gcdf4);
    x = gaminv(u.*c + gcdf4, A, B);
end
function [x1 x2] = gamrnd_trunc_antitethic(A, B, n)
    u = rand(1, n);
    gcdf4 = gamcdf(4, A, B);
    c = (gamcdf(25, A, B)-gcdf4);
    x1 = gaminv(u.*c + gcdf4, A, B);
    x2 = gaminv((1-u).*c + gcdf4, A, B);
end
function x = wbl_trunc(A, B, n)
    u = rand(1, n);
    x = wblinv(...
            u.*(wblcdf(25, A, B)-wblcdf(4, A, B)) + wblcdf(4, A, B),...
            A, B);
end
function [x1 x2] = wbl_trunc_antithetic(A, B, n)
    u = rand(1, n);
    x1 = wblinv(...
            u.*(wblcdf(25, A, B)-wblcdf(4, A, B)) + wblcdf(4, A, B),...
            A, B);
    x2 = wblinv(...
            (1-u).*(wblcdf(25, A, B)-wblcdf(4, A, B)) + wblcdf(4, A, B),...
            A, B);
end
function [x1 x2] = wbl_antithetic(A, B, n)
    u = rand(1, n);
    x1 = wblinv(u,A, B);
    x2 = wblinv((1-u), A, B);
end
