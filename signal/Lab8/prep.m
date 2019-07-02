
Fs = 1000;
N = 75;
xt = @(t) sin(2*pi*200*t) + sin(2*pi*220*t);
n = 0:(1/Fs):((N-1)/Fs);
xn = xt(n);
plot(abs(fft(xn)))