s = 24;
C = [1 zeros(1, 23) -.77];
A = conv([1 -.97], [1 zeros(1, s-1) -1]);
[A, C] = equalLength(A, C);
k = 8;
[F, G] = deconv( conv([1 zeros(1, k-1)], C), A)
%%
s = 24;
C = [1 zeros(1, 23) -.77];
A = conv([1 -.97], [1 zeros(1, s-1) -1]);
[A, C] = equalLength(A, C);
k = 26;
[F, G] = deconv( conv([1 zeros(1, k-1)], C), A);
%%
s = 24;
C = [1 zeros(1, 23) -.77];
A = conv([1 -.97], [1 zeros(1, s-1) -1]);
[A, C] = equalLength(A, C);
k = 1;
[~, G1] = deconv( conv([1 zeros(1, k-1)], C), A);
k = 8;
[~, G8] = deconv( conv([1 zeros(1, k-1)], C), A);
k = 26;
[~, G26] = deconv( conv([1 zeros(1, k-1)], C), A);

N = 1e5;
Nsim = 100;

burn = 300;
et = sqrt(.32).*randn(Nsim, N+burn);
Tt = filter(C, A, et, [], 2);
Tt = Tt(:,burn+1:end);

pred1 = filter(G1, C, Tt, [], 2);
pred8 = filter(G8, C, Tt, [], 2);
pred26 = filter(G26, C, Tt, [], 2);
%%
% Probability that the absolute value of the one step prediction
% is bigger than 2C
k = 1;
EG2p1 = abs((Tt(:,k+1:end)-pred1(:,k+1:end))) > 2;
PEG2p1 = mean(EG2p1, 2);
mean_PEG2p1 = mean(PEG2p1) %Estimated probability
mean_PEG2p1_CI = norminv(.975)*std(PEG2p1)/sqrt(Nsim) %95% confidence
%%
% Probability that the absolute value of the 8 step prediction
% is bigger than 2C
k = 8;
EG2p8 = abs((Tt(:,k+1:end)-pred8(:,k+1:end))) > 2;
PEG2p8 = mean(EG2p8, 2);
mean_PEG2p8 = mean(PEG2p8)
mean_PEG2p8_CI = norminv(.975)*std(PEG2p8)/sqrt(Nsim)
%%
% Variance estimate of the 26 step prediction error
k = 26;
vars = var((Tt(:,k+1:end)-pred26(:,k+1:end)), 0, 2);
pred26vars = mean(vars) % Variance estimate
pred26vars_CI = norminv(.975)*std(vars)/sqrt(Nsim) %with 95% confidence
%%
colors = [27,158,119
          217,95,2
          117,112,179
          231,41,138
          102,166,30
          230,171,2
          166,118,29]/255;
figure(1)
set(gca, 'Color', [222,235,247]/255);
hold on;
plot(Tt(1,:), 'LineWidth', 3, 'Color', colors(1,:))
plot(pred1(1,:), '-', 'LineWidth', 1.2, 'Color', colors(4,:))
plot(pred8(1,:), '-', 'LineWidth', 1.2, 'Color', colors(2,:))
plot(pred26(1,:), '-', 'LineWidth', 1.2, 'Color', colors(3,:))
legend({'Simulated temperate', 'k=1', 'k=8', 'k=26'})
hold off;