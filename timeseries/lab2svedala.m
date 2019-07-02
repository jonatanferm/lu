clear;
load("./tsadl/data/svedala");
y = svedala;
A = [1, -1.79, 0.84];
C = [1, -0.18, -0.11];

k = 3;
[CS, AS] = equalLength(C, A);
[Fk, Gk] = deconv(conv([1,zeros(1, k-1)], CS), AS);

yhat_k = filter(Gk, C, y);
ve = var(filter(A, C, y));
%Q6 estimated variance is 0.3895;
e = y - yhat_k; %Do I need to filter this?
me = mean(abs(e));
%Q7 est mean k3 = 1.2383 k26 = 2.6427 Seems low?
%Q8 theoretical Ve = (1+f_1 ^2 +...+ f_k ^2)sigma^2
%sum(Fk .^ 2) * 0.3895 = 13.1014; for k3 = 2.8528
%var(e) = 10.7563, for k3 = 2.6481
figure;
hold on;
plot(y)
plot(yhat_k)
%plot(y(k+1:end))
%plot(yhat_k(1:end-k))
hold off
%analyzets(e) Behaves like MA
var(e)
mean(e)