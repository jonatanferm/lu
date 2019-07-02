%% 1.1
clear;
clc;
n = 10000;

s_b = [1 zeros(1, n-1)];

f_b = 1:n;
f_b = 1./f_b;
f_b(2:2:n) = 0;
 
t_b = 1:n;
t_b = 1./(t_b.^2);
t_b(2:2:n) = 0;
t_b(3:4:n) = -t_b(3:4:n);

saw_b = 1:n;
saw_b = 1./saw_b;

S = @(b, t, w) sum(b.*sin( (1:n).*w.*t' ), 2);

t = linspace(-10, 10, 100000);

plot(t, S(saw_b(1:n), t, 1));

%% 2.2
clear;
clc;
n = 2; %spalter
b = 6; %bredd spalt
d = 20; %avstånd mellan spalter
pad = 100;
N = 10000;

f = [zeros(1, d) repmat([ones(1, b) zeros(1, d)], 1, n)];
f = [zeros(1, pad), f, zeros(1, pad)];
fn = length(f);

ft = abs(fftshift(fft(f, N))).^2;

f = f.*max(ft)/2;
f = interp1(1:fn, f, (1:N)*(fn/N), 'next', 'extrap');

figure;
hold on;
plot(f)
plot(ft)
hold off;

