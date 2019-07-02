clear;
clc;
x = linspace(-10, 10);
k = abs(x) < 1;
y = 1 ./ (sqrt(1 + x.^2));
figure;
hold on;
plot(x, y);
plot(x, k);
hold off;
figure;
plot(conv(k, y));
figure;
plot(atan(x + .5) - atan(x-.5))