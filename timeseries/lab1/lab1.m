A1 = [1 -1.79 0.84];
C1 = [1 -0.18 -0.11];

A2 = [1 -1.79];
C2 = [1 -0.18 -0.11];

armap1 = idpoly(A1, [], C1);
armap2 = idpoly(A2, [], C2);

N = 200;
sigma2 = 1.5;
e = sqrt(sigma2) * randn(N, 1);
y1 = filter(armap1.c, armap1.a, e);
y2 = filter(armap2.c, armap2.a, e);
%%
plot(y1)
%%
plot(y2)
%%
m = 100;
rtheo = kovarians(armap1.c, armap1.a, m);
stem(0:m, rtheo*sigma2)
hold on
rest = covf(y1, m+1);
stem(0:m, rest, 'r')
%%
normplot(y1)
%%
figure
hold on
plot(acf(y1, 5))
plot(pacf(y1, 5))
refline(0, 2/sqrt(N))
refline(0, -2/sqrt(N))
hold off
%%
data = iddata(y1, [], 1);
%%
ar_model = arx(y1, 3);
%%
arma_model = armax(y1, [2 0]);
%%
e_hat = filter(arma_model.A, arma_model.C, y1);
figure;
hold on;
plot(acf(e_hat, 20));
plot(pacf(e_hat, 20));
refline(0, 2/sqrt(N))
refline(0, -2/sqrt(N))
hold off
%%
clear;
%%
n = 5000;
A = [1 -1.35 0.43];
sigma2 = 4;
noise = sqrt(sigma2) * randn(n, 1);
y = filter(1, A, noise);
y = y(101:end);
subplot(211); plot(y);
subplot(212); plot(noise);
%%
n_est = floor(2/3 * n);
y_est = iddata(y(1:n_est));
y_val = iddata(y(n_est + 1:end));
NN = (1:10)';
%%
n_orders = zeros(100, 1);
n_aics = zeros(100, 1);
%%
for i = 1:100
    noise = sqrt(sigma2) * randn(n, 1);
    y = filter(1, A, noise);
    y = y(101:end);
    y_est = iddata(y(1:n_est));
    y_val = iddata(y(n_est + 1:end));
    V = arxstruc(y_est, y_val, NN);
    n_orders(i) = selstruc(V, 0);
    n_aics(i) = selstruc(V, 'aic');
end
%%
subplot(211); hist(n_orders);
subplot(212); hist(n_aics);
%%
ar_model = arx(y, n_orders(end));
ar_model.NoiseVariance
ar_model.CovarianceMatrix
present(ar_model)
%%
load('/Users/ferm/lu/timeseries/tsadl/data/data.dat');
idata = iddata(data);
%%
arma11 = armax(idata, [1 1]);
arma22 = armax(idata, [2 2]);
rar11 = resid(arma11, idata);
rar22 = resid(arma22, idata);
%%
A = [1 -1.5 0.7];
C = [1 zeros(1, 11) - 0.5];
A12 = [1 zeros(1, 11) -1];
A_star = conv(A, A12);
N = 600;
e = randn(N + 100, 1);
y = filter(C, A_star, e);
y = y(100:end);
%%
plot(y)
%%
normplot(y)
%%
figure
hold on
plot(acf(y, 50), 'b')
plot(pacf(y, 50), 'g')
refline(0, 2/sqrt(N))
refline(0, -2/sqrt(N))
hold off
%%
y_s = filter(A12, 1, y);
data = iddata(y_s);
%%
normplot(y_s)
%%
figure
hold on
plot(acf(y_s, 50), 'b')
plot(pacf(y_s, 50), 'g')
refline(0, 2/sqrt(N))
refline(0, -2/sqrt(N))
hold off
%%
model_init = idpoly([1 0 0], [], []);
model_armax = pem(data, model_init);
%%
res=resid(model_armax, data)
%%
normplot(res.OutputData)
%%
figure
hold on
plot(acf(res.OutputData, 50), 'b')
plot(pacf(res.OutputData, 50), 'g')
refline(0, 2/sqrt(N))
refline(0, -2/sqrt(N))
hold off
%%
model_init = idpoly([1 0 0], [], [1 zeros(1, 12)]);
model_init.Structure.c.Free = [zeros(1, 12) 1];
model_armax = pem(data, model_init);
%%
res=resid(model_armax, data);
%%
normplot(res.OutputData)
%%
figure
hold on
plot(acf(res.OutputData, 50), 'b')
plot(pacf(res.OutputData, 50), 'g')
refline(0, 2/sqrt(N))
refline(0, -2/sqrt(N))
hold off
%%
clear;
%%
load('/Users/ferm/lu/timeseries/tsadl/data/svedala.mat');
N = size(svedala, 1);
%%
plot(svedala);
%%
figure
hold on
plot(acf(svedala, 50), 'b')
plot(pacf(svedala, 50), 'g')
refline(0, 2/sqrt(N))
refline(0, -2/sqrt(N))
hold off
%% Def sesonality, seems to be 24
s_s = conv(svedala, [1 zeros(1, 23) -1]);
N_s = size(s_s, 1);
%%
figure
hold on
plot(acf(s_s, 10), 'b')
plot(pacf(s_s, 10), 'g')
refline(0, 2/sqrt(N_s))
refline(0, -2/sqrt(N_s))
hold off
%% Try modelling AR(2)
s_s = iddata(s_s);
ar2 = armax(s_s, [2 0]);
res = resid(ar2, s_s);
%%
figure
hold on
plot(acf(res.OutputData, 50), 'b')
plot(pacf(res.OutputData, 50), 'g')
refline(0, 2/sqrt(N_s))
refline(0, -2/sqrt(N_s))
hold off
%% AR(24) looks strong
model_init = idpoly([1 zeros(1, 2)], [], [1 zeros(1, 23) 1]);
model_init.Structure.a.Free = [1 1 1];
model_init.Structure.c.Free = [zeros(1, 24) 1];
model_armax = pem(s_s, model_init);
res = resid(model_armax, s_s);
%%
model_init = idpoly([1 zeros(1, 25)], [], [1 zeros(1, 23) 0 0]);
model_init.Structure.a.Free = [1 1 1 zeros(1, 20) 0 1 0];
model_init.Structure.c.Free = [zeros(1, 23) 0 0 0];
model_armax = pem(svedala, model_init);
res = resid(model_armax, svedala);
%%
analyzets(res)