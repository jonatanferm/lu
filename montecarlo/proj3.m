clear;
clc;
load('coal_mine_disasters.mat');
%% colors for plots;
cl = [166,206,227
    31,120,180
    178,223,138
    51,160,44
    251,154,153
    227,26,28
    253,191,111
    255,127,0
    202,178,214
    106,61,154
    255,255,153
    177,89,40]/255;

ct = [141,211,199
    255,255,179
    190,186,218
    251,128,114
    128,177,211
    253,180,98
    179,222,105
    252,205,229
    217,217,217
    188,128,189
    204,235,197
    255,237,111]/255;

cline = [102,194,165
    252,141,98
    141,160,203
    231,138,195
    166,216,84
    255,217,47]/255;
%%
%for psi = logspace(-5, 5, 11)
%for d = 2:12
for rho = logspace(-5, 3, 7)
    d = 5;
    %rho = 1;
    psi = 10;
    N = 1e5;
    burn = 0;
    target_rates = [.25]';
    tic
    bigt = T(1:end);
    [thetas, lambdas, breakpoints, accs] = hybrid_sampler(bigt, psi, rho, d, N, burn, target_rates);
    toc
    %%
    accfig = figure('Name', sprintf('Acceptance rate d=%i', d),'Visible','Off');
    plot(movmean(accs, 500, 2)', 'LineWidth', 2)
    set(gca,'FontSize',18)
    ylabel('Moving acceptance rate, 500 samples.','FontSize', 18, 'Interpreter', 'latex')
    xlabel('Sample number','FontSize', 18, 'Interpreter', 'latex')
    bins = floor(N/250);
    %saveas(accfig,sprintf('d%i_acc.png', d))
    %saveas(accfig,sprintf('psi%f_acc.png', psi))
    saveas(accfig,sprintf('rho%f_acc.png', rho))
    %%
    lambfig = figure('Name', sprintf('Lambdas d=%i', d),'Visible','Off');
    xl = [0 10];
    p2s = .98-.9*(1/d);
    for i=1:d
        subplot(d,1, i)
        histogram(lambdas(:,i), bins,'facecolor', [115,115,115]/255,'facealpha',1,'edgecolor','none', 'Normalization','probability');
        xlim(xl)
        %pos = get(gca, 'Position');
        %pos(1) = 0.055;
        %pos(2) = p2s;
        %pos(3) = 0.9;
        %pos(4) = .9*(1/d);
        %p2s = pos(2)-1.02*pos(4);
        %set(gca, 'Position', pos)
        yticks([])
        xticks([])
        ylabel(sprintf('$\\lambda_%i$', i),'FontSize', 18, 'Interpreter', 'latex')
    end
    xticks(xl(1):xl(2))
    set(gca,'TickLength',[0 0]);
    set(gca,'FontSize',16)
    %saveas(lambfig,sprintf('d%i_lambda.png', d))
    %saveas(lambfig,sprintf('psi%f_lambda.png', psi))
    saveas(lambfig,sprintf('rho%f_lambda.png', rho))
    tfig = figure('Name', sprintf('Break points d=%i', d),'Visible','Off');
    hold on
    xl = [1658 1980];
    for i=1:d-1
        %for i=2
        subplot(d-1,1, i)
        histogram(breakpoints(:,i+1), bins,'facecolor', [115,115,115]/255,'facealpha',1,'edgecolor','none', 'Normalization','probability')
        xlim(xl)
        %pos = get(gca, 'Position');
        %pos(1) = 0.055;
        %pos(2) = pos(2)-0.05;
        %pos(3) = 0.9;
        %pos(4) = .9*(1/(d-1));
        %set(gca, 'Position', pos)
        yticks([])
        xticks([])
        %xtck = xticks;
        ylabel(sprintf('$t_%i$', i+1),'FontSize', 18, 'Interpreter', 'latex')
    end
    xticks(xl)
    xlabel('year','FontSize', 18, 'Interpreter', 'latex')
    set(gca,'FontSize',16)
    %saveas(tfig,sprintf('d%i_t.png', d))
    %saveas(tfig,sprintf('psi%f_t.png', psi))
    saveas(tfig,sprintf('rho%f_t.png', rho))
end
%%

%% 2

clear;
clc;
load('atlantic.txt')
T=3*14*100;
%Bootstrap
N=1000;
A = atlantic;
[beta,mu]=est_gumbel(A);
tau_hat = mean(A);
boot = zeros(1,N);
Ai_star=zeros(582,1);
Aistar=zeros(582,1);
betabootvector=zeros(1,N);
mubootvector=zeros(1,N);
for i=1:N
    for j=1:582
        u=rand(1);
        Aistar(j)=inversgumb(u,mu,beta);
    end
    [beta_boot, mu_boot] =est_gumbel(Aistar);
    betabootvector(i)=beta_boot;
    mubootvector(i)=mu_boot;
    
end
betaboot=mean(betabootvector);
muboot=mean(mubootvector);

delta = sort(mubootvector -mu); % sorting to obtain quantiles
alpha = 0.05; % CB level
Lmu = mu - delta(ceil((1 - alpha/2)*N)); % constructing CB
Umu = mu - delta(ceil(alpha*N/2));

delta = sort(betabootvector -beta); % sorting to obtain quantiles
alpha = 0.05; % CB level
Lbe = beta - delta(ceil((1 - alpha/2)*N)); % constructing CB
Ube = beta - delta(ceil(alpha*N/2));

%%
%2c

T=3*14*100;
U=1-1/T;
N=1000;

Upperboundary=0;
maxwave = zeros(length(betabootvector),1);
for p=1:length(betabootvector)
    maxwave(p)=inversgumb(U, mubootvector(p), betabootvector(p));
end

[beta,mu]=est_gumbel(A);
maxwavestart=inversgumb(U, mu, beta);

delta = sort(maxwave -maxwavestart); % sorting to obtain quantiles
alpha = 0.05; % CB level
hmm1= maxwavestart - delta(ceil((1 - alpha)*N)); % constructing CB
hmm2= maxwavestart - delta(ceil(alpha*N));
%%
function [thetas, lambdas, t, accs] = hybrid_sampler(T, psi, rhos, d, N, burn, target_rates)
    if isscalar(rhos)
        rhos = ones(d-1, 1)*rhos;
    end
    if isscalar(target_rates)
        target_rates = ones(d-1, 1)*target_rates;
    end
    bs = 50;
    num_samples = N+burn;
    tenprcnt = ceil(num_samples/10);
    thetas = nan(num_samples, 1);
    thetas(1) = gamrnd(2, 1/psi);
    lambdas = nan(num_samples, d);
    lambdas(1, :) = gamrnd(2, 1/thetas(1), 1, d);
    % Initiate breakpoints at points, more likely at places where intensity
    % changes
    t1 = 1658;
    td1 = 1980;
    at = abs(diff(diff(T)));
    [~, ind] = histc((rand(1,d-1)+(0:(d-2)))/d, [-inf; cumsum(at)/sum(at); inf]);
    brpts = [t1 T(ind)' td1];
    brpts = [brpts; nan(size(brpts) + [num_samples-2 0])];
    %
    accs = nan(d-1, num_samples-1);
    for i = 2:num_samples
        thetas(i) = sample_theta(lambdas(i-1,:), psi);
        [brpts(i,:), accs(:,i-1)] = MH(brpts(i-1,:), lambdas(i-1,:), T, rhos);
        lambdas(i,:) = sample_lambdas(thetas(i), brpts(i,:), T);
        if mod(i, bs) == 0 %% Varying rho from book. Unsure how to make it
            %work for beta proposal
            rhos = rhos .* exp( .5*sign(  mean(accs(:,i-bs+1:i-1), 2)-target_rates  )*min(1e-2, 1/sqrt(i)) );
        end
        if mod(i, tenprcnt) == 0
            %ACCEPT_RATE = mean(accs(:,i-bs+1:i-1), 2)'
            %CURRENT_RHOS = rhos'
            PRCNT_DONE = i/num_samples
        end
    end
    thetas = thetas(burn+1:end);
    lambdas = lambdas(burn+1:end,:);
    t = brpts(burn+1:end,:);
    accs = accs(:,burn+1:end);
end

function [brpts, accs] = MH(old_brpts, lambdas, T, rhos)
    brpts = old_brpts;
    accs = nan(1, numel(brpts)-2);
    for i = 2:numel(brpts)-1
        brpts_star = brpts;
        %brpts_star(i) = brpts(i-1) + betarnd(rhos(i-1), rhos(i-1))*(brpts(i+1)-brpts(i-1));
        brpts_star(i) = brpts(i) + rhos(i-1)*(brpts(i+1)-brpts(i-1))*(2*rand()-1);
        if sum(sign(diff(brpts_star))-1) == 0
            fq = log_ft_posterior(brpts_star, lambdas, T, i) - ...
                log_ft_posterior(brpts, lambdas, T, i);
            %rq  = log_r(brpts, rhos(i-1), i) - log_r(brpts_star, rhos(i-1), i);
            if log(rand()) < fq
                brpts = brpts_star;
                accs(i-1) = 1;
            else
                accs(i-1) = 0;
            end
        else
            accs(i-1) = 0;
        end
    end
end

function lr = log_r(brpts, rho, i)
    lr = gammaln(2*rho) - 2*gammaln(rho) ...
        + (rho-1)*( log(brpts(i)-brpts(i-1))+log(brpts(i+1)-brpts(i)) ) ...
        - (2*rho-1)*log(brpts(i+1)-brpts(i-1));
end
function ld = log_ft_posterior(brpts, lambdas, T, ii)
    t1mt = (brpts(2:end) - brpts(1:end-1));
    s1 = -sum(t1mt.*lambdas);
    s2 = sum((nt(brpts, T)+1).*log(lambdas));
    s3 = sum(log(t1mt));
    ld = s1+s2+s3;
end
function ld = log_ft_posterior_ii(brpts, lambdas, T, ii)
    nts = nti(brpts, T, [ii-1, ii]);
    s1 = (lambdas(ii)-lambdas(ii-1))*brpts(ii);
    s2 = sum((nts+1)'.*log(lambdas(ii-1:ii)));
    s3 = sum(log(brpts(ii:ii+1)-brpts(ii-1:ii)));
    ld = s1+s2+s3;
end
function lambdas = sample_lambdas(theta, brpts, T)
    alpha = nt(brpts, T) + 2;
    beta = brpts(2:end) - brpts(1:end-1) + theta;
    lambdas = gamrnd(alpha, 1./beta);
end
function theta = sample_theta(lambdas, psi)
    alpha = 2*numel(lambdas)+2;
    beta = sum(lambdas)+psi;
    theta = gamrnd(alpha, 1/beta);
end
function x = nt(breakpoints, T)
    N = histc(T, breakpoints);
    x=N(1:end-1)';
end
function x = nti(breakpoints, T, i)
    N = histc(T, [-inf breakpoints(2:end-1) inf]);
    x=N(i);
end
function [ out ] = inversgumb(u,mu, beta )
%INVERSGUMB Summary of this function goes here
%   Detailed explanation goes here
out=mu-beta.*(log(-log(u)));

end
