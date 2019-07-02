clear;
clc;
load('projekt.mat')
%%
fs = 256;
f = (1:(length(acc)))./fs;
%% fft acc
N = 5000; %Resolution of ft
ff = (1:N/2).*(fs/(N)); %frequencies
f_acc = fft(acc, N, 2); %calculate ft
f_acc = f_acc(:,1:end/2)'; %discard mirrored half and flip it.
%% plot of accelerometer frequencies
figure(1)
plot(ff, abs(f_acc), 'Linewidth', 1.5)
legend(['$A_x$'; '$A_y$'; '$A_z$'], 'fontsize', 18, 'interpreter', 'latex')
axis([0 3 0 1e7])
xlabel('Hz', 'interpreter', 'latex', 'fontsize', 18)
%% fft ppg
f_ppg = fft(ppg, N, 2); %calculate ft
f_ppg = f_ppg(:,1:end/2)'; %discard mirrored half and flip it.
%% plot of ppg frequencies
figure(2)
plot(ff, abs(f_ppg), 'Linewidth', 1.5)
%legend(['$F], 'fontsize', 18, 'interpreter', 'latex')
axis([0 3 0 1e7])
xlabel('Hz', 'interpreter', 'latex', 'fontsize', 18)
ylabel('F', 'interpreter', 'latex', 'fontsize', 18)
%% filter accelerometer data
r = .99;
w0 = [.5630 .9728];
z = exp (-1j*2*pi*[-w0 w0]/fs);
p = exp (-1j*2*pi*[-w0 w0]/fs)*r;
B = poly(z);
A = poly(p);
%% plot filter frequency response
[H,W] = freqz(B,A,N);
figure(3)
plot(W, abs(H), 'LineWidth', 1.5)
axis([0 .2 0 1.1])
%% plot impulse response
[H,T] = impz(B,A,N);
figure(4)
plot(T(2:1000), H(2:1000), 'LineWidth', 1.5)
%set(gca, 'YScale', 'log')
%set(gca, 'XScale', 'log')
%axis([0 .2 0 1.1])
%% filter acc data
y = filter(B, A, acc, [], 2);
%% fft fitlered accdata
f_y = fft(y, N, 2); %calculate ft
f_y = f_y(:,1:end/2)'; %discard mirrored half and flip it.
%% plot of filtered accelerometer frequencies
figure(5)
plot(ff, abs(f_y), 'Linewidth', 1.5)
legend(['$A_x$'; '$A_y$'; '$A_z$'], 'fontsize', 18, 'interpreter', 'latex')
axis([0 3 0 1e7])
xlabel('Hz', 'interpreter', 'latex', 'fontsize', 18)
%% filter ppg data
py = filter(B, A, ppg, [], 2);
%% fft fitlered ppg data
f_py = fft(py, N); %calculate ft
f_py = f_py(:,1:end/2)'; %discard mirrored half and flip it.
%% plot of filtered ppg frequencies
figure(6)
hold on
plot(ff, abs(f_py), 'Linewidth', 1.5)
plot(ff, abs(f_ppg), 'Linewidth', 1.5)
hold off
legend({'Filtered signal', 'Unfiltered signal'})
axis([0 3 0 1e7])
xlabel('Hz', 'interpreter', 'latex', 'fontsize', 18)
%% pzplot
o = linspace(0, 2*pi, 1000);
figure(7)
hold on;
pzmap(tf(A, B))
plot(sin(o), cos(o))
axis([0.99 1.0001 -0.03 0.03])
hold off
%%
%% plot fft of ecg
f_ecg = fft(ecg, N, 2); %calculate ft
f_ecg = f_ecg(:,1:end/2)'; %discard mirrored half and flip it.
figure(8)
plot(ff, abs(f_ecg), 'Linewidth', 1.5)
legend(['$A_x$'; '$A_y$'; '$A_z$'], 'fontsize', 18, 'interpreter', 'latex')
axis([0 3 0 1e7])
xlabel('Hz', 'interpreter', 'latex', 'fontsize', 18)

