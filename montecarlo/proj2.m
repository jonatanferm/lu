clear;
clc;
%% True values from oeis.org
true_d2 = [1, 4, 12, 36, 100, 284, 780, 2172, 5916, 16268, 44100, 120292, 324932, 881500, 2374444, 6416596, 17245332, 46466676, 124658732, 335116620, 897697164, 2408806028, 6444560484, 17266613812, 46146397316, 123481354908, 329712786220, 881317491628];
true_d3 = [1, 6, 30, 150, 726, 3534, 16926, 81390, 387966, 1853886, 8809878, 41934150, 198842742, 943974510, 4468911678, 21175146054, 100121875974, 473730252102, 2237723684094, 10576033219614, 49917327838734, 235710090502158, 1111781983442406, 5245988215191414, 24730180885580790, 116618841700433358, 549493796867100942, 2589874864863200574, 12198184788179866902, 57466913094951837030, 270569905525454674614];
c_2_71 = 4190893020903935054619120005916;
cn2 = nan(1, 71);
cn2(1)=4; 
cn2(2)=12;
cn2(3)=36;
cn2(4)=100;
cn2(5)=284;
cn2(6)=780;
cn2(7)=2172;
cn2(8)=5916;
cn2(9)=16268;
cn2(10)=44100;
cn2(11)=120292;
cn2(12)=324932;
cn2(13)=881500;
cn2(14)=2374444;
cn2(15)=6416596;
cn2(16)=17245332;
cn2(17)=46466676;
cn2(18)=124658732;
cn2(19)=335116620;
cn2(20)=897697164;
cn2(21)=2408806028;
cn2(22)=6444560484;
cn2(23)=17266613812;
cn2(24)=46146397316;
cn2(25)=123481354908;
cn2(26)=329712786220;
cn2(27)=881317491628;
cn2(28)=2351378582244;
cn2(29)=6279396229332;
cn2(30)=16741957935348;
cn2(31)=44673816630956;
cn2(32)=119034997913020;
cn2(33)=317406598267076;
cn2(34)=845279074648708;
cn2(35)=2252534077759844;
cn2(36)=5995740499124412;
cn2(37)=15968852281708724;
cn2(38)=42486750758210044;
cn2(39)=113101676587853932;
cn2(40)=300798249248474268;
cn2(41)=800381032599158340;
cn2(42)=2127870238872271828;
cn2(43)=5659667057165209612;
cn2(44)=15041631638016155884;
cn2(45)=39992704986620915140;
cn2(46)=106255762193816523332;
cn2(47)=282417882500511560972;
cn2(48)=750139547395987948108;
cn2(49)=1993185460468062845836;
cn2(50)=5292794668724837206644;
cn2(51)=14059415980606050644844;
cn2(52)=37325046962536847970116;
cn2(53)=99121668912462180162908;
cn2(54)=263090298246050489804708;
cn2(55)=698501700277581954674604;
cn2(56)=1853589151789474253830500;
cn2(57)=4920146075313000860596140;
cn2(58)=13053884641516572778155044;
cn2(59)=34642792634590824499672196;
cn2(60)=91895836025056214634047716;
cn2(61)=243828023293849420839513468;
cn2(62)=646684752476890688940276172;
cn2(63)=1715538780705298093042635884;
cn2(64)=4549252727304405545665901684;
cn2(65)=12066271136346725726547810652;
cn2(66)=31992427160420423715150496804;
cn2(67)=84841788997462209800131419244;
cn2(68)=224916973773967421352838735684;
cn2(69)=596373847126147985434982575724;
cn2(70)=1580784678250571882017480243636;
cn2(71)=4190893020903935054619120005916;
%%
d = 2;
N = 1e7;
n_list = [1,2,3,5,10,20,30,40,50,60,70];
%% 3
taus = zeros(size(n_list));
for n = n_list
    c = random_walks_sis_naive(N, n, d);
    taus(n_list == n) = c;
end
taus3 = taus;
%%
make_table(n_list, taus, cn2(n_list))
%% 4
taus = zeros(size(n_list));
for n = n_list
    [weights, ~] = random_walks_sis(N, n, d);
    taus(n == n_list) = mean(weights);
end
taus4 = taus;
%%
[w_70, ~] = random_walks_sis(1e7, 70, 2);
%%
w_70_f = w_70;
w_70_f(w_70 == 0) = 1;
[n, xout] = hist(log(w_70_f));
bar(xout, n, 'barwidth', 1, 'basevalue', 1);
set(gca,'YScale','log')
xlabel('log($w$)', 'FontSize', 18,'Interpreter', 'latex')
sum(w_70==0)
%%
w_70_sorted = sort(w_70, 'descend');
y = cumsum(w_70_sorted)./sum(w_70_sorted);
figure;
hold on;
plot(y, 'k','LineWidth', 1.5)
set(gca,'XScale','log')
refline(0, .95)
hold off;
%%
make_table(n_list, taus, cn2(n_list))
%% 5
taus = zeros(size(n_list));
for n = n_list
    n
    [weights, ~] = random_walks_sisr(N, n, d);  
    taus(n == n_list) = prod(mean(weights));
end
taus5 = taus;
%%
[w_70_sisr, ~] = random_walks_sisr(1e7, 70, 2);
%%
w_70_f(w_70 == 0) = 1;
[n, xout] = hist(reshape(w_70_sisr, [], 1));
bar(xout, n, 'barwidth', 1, 'basevalue', 1);
%set(gca,'YScale','log')
xlabel('$w$', 'FontSize', 18,'Interpreter', 'latex')
sum(w_70==0)
%%
make_table(n_list, taus, cn2(n_list))
%%
clc
make_table(n_list, taus5, cn2(n_list))
%%
mmms = nan(1000, 1);
for i = 1:numel(mmms)
    i
    [weights, ~]= random_walks_sisr(1e4, 10, 2);
    mmms(i) = prod(mean(weights));
end
%%
mean(mmms)
std(mmms)*norminv(.975)*sqrt(1000)^(-.5)
%%
%%
%% 6 
%% Estimates from start
maxn = 30;
lncn = zeros(maxn, 1);
lnn = log(1:maxn)';
for n = 1:maxn
    [weights, ~] = random_walks_sisr(1e4, n, 2);
    lncn(n) = log(prod(mean(weights)));
end
regr = [ones(n, 1) (1:n)' lnn];
estimates = zeros(3, maxn);
for n = 1:maxn
    estimates(:,n) = regr(1:n,:)\lncn(1:n,:);
    estimates(1:2,n) = exp(estimates(1:2,n));
    estimates(3,n) = estimates(3,n)+1;
end
%% PLOT ESTIMATES
colors = [27,158,119
217,95,2
117,112,179]/255;

figure;
plot(estimates', 'LineWidth', 2)
legend('A_2','\mu_2', '\gamma_2', 'FontSize',18)
legend boxoff
xlabel('Number of steps, n.', 'FontSize', 18)
ylabel('Parameter estimate.', 'FontSize', 18)
%% Several estimates not from zero
step = 1;
number_of_steps = 20;
nn = step*number_of_steps;
ests_per_step = 10;
skip = 30;
estimates = zeros(3, 100);
tic
for estno = 1:100
    toc
    estno
    tic
    lncn = zeros(number_of_steps*ests_per_step, 1);
    lnn = repmat(log((skip+1):step:(skip+nn))', ests_per_step, 1);
    nn_arr = repmat(((skip+1):step:(skip+nn))', ests_per_step, 1);
    c = 0;
    for n = nn_arr'
        c = c + 1;
        [weights, ~] = random_walks_sisr(1e3, n, 2);
        lncn(c) = log(prod(mean(weights)));
    end
    regr = [ones(number_of_steps*ests_per_step, 1) nn_arr lnn];
    estimates2 = regr\lncn;
    estimates2(1:2) = exp(estimates2(1:2));
    estimates2(3) = estimates2(3)+1;
    estimates(:,estno) = estimates2
end
save('paramestimatesd2.mat', 'estimates')
%%
gammas = zeros(50, 1);
num_steps = 10;
for start = 1:50
    cend = start+num_steps;
    logcn_true = log(cn2(start:cend))';
    regr_true = [ones(size(logcn_true)) (start:cend)' log((start:cend))'];
    true_ests = regr_true\logcn_true;
    gammas(start) = true_ests(3)+1;
end
figure;
hold on;
plot(gammas, '-r', 'LineWidth', 1.5)
refline(0, 43/32)
ylim([1.2 1.4])
legend('$\hat{\gamma}_2$', '$\gamma_2 = 43/32$', 'FontSize', 18, 'Interpreter', 'latex')
xlabel('Number of skipped steps.', 'FontSize', 18)
hold off;
%%
CI = norminv(0.975)*std(estimates, 0, 2)*sqrt(100)^-1;
[mean(estimates, 2) + CI, mean(estimates, 2) - CI]
mean(estimates, 2);
%%
means = mean(estimates, 2);
stds = std(estimates, 0, 2)
%% 9
step = 1;
number_of_steps = 20;
nn = step*number_of_steps;
ests_per_step = 10;
skip = 30;
estimates = zeros(2, 100);
tic
for estno = 1:100
    toc
    estno
    tic
    lncn = zeros(number_of_steps*ests_per_step, 1);
    nn_arr = repmat(((skip+1):step:(skip+nn))', ests_per_step, 1);
    c = 0;
    for n = nn_arr'
        c = c + 1;
        [weights, ~] = random_walks_sisr(1e3, n, 5);
        lncn(c) = log(prod(mean(weights)));
    end
    regr = [ones(number_of_steps*ests_per_step, 1) nn_arr];
    estimates2 = regr\lncn;
    estimates2(1:2) = exp(estimates2(1:2));
    estimates(:,estno) = estimates2
end
save('paramestimatesd2.mat', 'estimates')

%%
CI = norminv(0.975)*std(estimates, 0, 2)*sqrt(100)^-1;
[mean(estimates, 2) + CI, mean(estimates, 2) - CI]
means = mean(estimates, 2)
stds = std(estimates, 0, 2)
%%
(9-means(2))/stds(2)
%%
means(1) - norminv(0.9994)*stds(1)*sqrt(100)^-1
%%
D = 5;
2*D - 1 - 1/(2*D) - 3/((2*D)^2) - 16/((2*D)^3) 
%% Estimates from tau5
regr = [ones(size(n_list))' n_list' log(n_list)'];
skp = 9;
tau5params = regr(skp:end,:)\log(taus5(skp:end)');
A = exp(tau5params(1))
mu = exp(tau5params(2))
gamma = tau5params(3)+1
%%
A_save = A;
mu_save = mu;
gamma_save = gamma;
%% 7
%%
print_c(n_list, cn2(n_list))
function [weights, walks] = random_walks_sisr(N, n, d)
    weights = ones(N, n);
    base = n+1;
    visited = zeros(N, n+1);
    steps = [(base * (0:d-1) + 1)'; -(base * (0:d-1) + 1)'];
    Nl = 1:N;
    for cn = 2:n+1        
        possteps = reshape(visited(:, cn-1)+steps', N, 1, []); % all possible steps
        possteps(sum(possteps - visited(:,1:cn-1) == 0, 2) ~= 0) = NaN; % remove visited
        possteps = reshape(possteps, N, d*2);
        possteps = sort(possteps, 2, 'MissingPlacement', 'last'); % Nan last
        num_pos_steps = sum(~isnan(possteps), 2); % number of possible steps
        weights(:,cn-1) = num_pos_steps;
        ind = ceil(rand(N, 1) .* num_pos_steps); % pick random steps
        visited(ind ~= 0,cn) = possteps(sub2ind(size(possteps), Nl(ind ~= 0), (ind(ind ~= 0))'))';

        CW=cumsum([0 weights(:,cn-1)']);
        [~,ind] = histc((rand(1,N)+(0:(N-1)))/N,CW/CW(end));
 %       [~,ind] = histc(rand(1,N),CW/CW(end));
        %weights(:, cn) = weights(ind, cn);
        visited = visited(ind,:);
    end
    walks = visited;
end
function [weights, walks] = random_walks_sis(N, n, d)
    weights = ones(N, 1);
    base = n+1;
    visited = zeros(N, n+1);
    steps = [(base * (0:d-1) + 1)'; -(base * (0:d-1) + 1)'];
    Nl = 1:N;
    for cn = 2:n+1
        possteps = reshape(visited(:, cn-1)+steps', N, 1, []);
        possteps(sum(possteps - visited(:,1:cn-1) == 0, 2) ~= 0) = NaN;
        possteps = reshape(possteps, N, d*2);
        possteps = sort(possteps, 2, 'MissingPlacement', 'last');
        num_pos_steps = sum(~isnan(possteps), 2);
        weights = weights .* num_pos_steps;
        ind = ceil(rand(N, 1) .* num_pos_steps);
        visited(ind ~= 0,cn) = possteps(sub2ind(size(possteps), Nl(ind ~= 0), (ind(ind ~= 0))'))';
    end
    walks = visited;
end
function c = random_walks_sis_naive(N, n, d)
    base = n+1;
    visited = zeros(N, n+1);
    steps = [(base * (0:d-1) + 1)'; -(base * (0:d-1) + 1)'];
    for cn = 2:n+1
        visited(:, cn) = visited(:, cn-1)+steps(randi(2*d, N, 1));
    end
    sorted_v = sort(visited, 2);
    c = ((2*d)^n)*mean(sum(sorted_v(:,2:end) - sorted_v(:,1:end-1) == 0, 2) == 0);
end
function make_plots(steps, chat, c)
    figure;
    hold on;
    plot(steps, chat, 'xr');
    plot(steps, c, 'ok');
    legend('$\hat{c}_n(2)$','$c_{n}(2)$', 'FontSize',18, 'Location', 'northwest', 'Interpreter', 'latex')
    legend boxoff
    xlabel('Number of steps, n.', 'FontSize', 18)
    xticks(steps)
    hold off
    set(gca, 'YScale', 'log')
    figure;
    stem(steps, abs((chat-c)))
    xlabel('Number of steps, n.', 'FontSize', 18)
    ylabel('$|c_{n}(2) - \hat{c}_n(2)|$', 'FontSize',18, 'Interpreter', 'latex')
    xticks(steps)
    set(gca, 'YScale', 'log')
    figure;
    stem(steps, abs((chat-c))./c )
    xlabel('Number of steps, n.', 'FontSize', 18)
    ylabel('$|c_{n}(2) - \hat{c}_n(2)|/c_{n}(2)$', 'FontSize',18, 'Interpreter', 'latex')
    xticks(steps)
    set(gca, 'YScale', 'log')
end
function print_c(steps, c)
    formatSpec = 'c_{%0.0f}(2) = & \\num[group-separator={,}]{%0.0f} \\\\ \n';
    fprintf(formatSpec,[steps; c]) 
end
function make_table(steps, chat, c)
    formatSpec = '$%0.6g$ & $%0.6g$ & $%0.6g$ & $%0.6g$ \\\\ \n'; 
    fprintf(formatSpec,[steps; chat; c; abs(chat-c)]) 
end