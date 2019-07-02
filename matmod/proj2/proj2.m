clear;
clc;
%%
data = xlsread('pk',1);
data = data(:,2:end);
data = permute(reshape(data', [3, 10, 10]), [2 1 3]);
%%
set(groot,'defaultAxesColorOrder',jet(10));
plot(squeeze(data(:,1,:)), squeeze(data(:,2,:)), 'x-','LineWidth', 1.5);
%%
C = @(t) @(x) sum(x(end/2+1:end).*exp( t*-abs(x(1:end/2)) ), 2)';
t = data(:,1,2);
Ct = C(t);
Cthat = data(:,2,1)';

fun = @(x) sum((Ct(x) - Cthat).^2, 'all');
%%
options = optimset('MaxIter', 5000,'MaxFunEvals',10000,'PlotFcns',@optimplotfval, 'Display', 'iter');
x = fminsearch(fun, [0 0 0 0 0 0 0 0 0 0], options);
%%
tt = linspace(0, 100)';
Ct2 = C(tt);
figure;
hold on;
plot(tt, Ct2(x));
plot(t, data(:,2,1));
hold off;
