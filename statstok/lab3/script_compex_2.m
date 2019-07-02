%% Stationary stochastic processes, computer exercise 2
%% 2.1 Study of AR(1)-processes
close all
clear all
clc

% Setup
a1 = -0.5;
C = 1;
A = [1 a1];

% Spectrum
[H,w] = freqz(C,A);
R = abs(H).^2;
figure
plot(w/2/pi,R)
title('Spectral density')
xlabel('Frequency')

% Covariance function
H = freqz(C,A,512,'whole');
Rd = abs(H).^2;
r = ifft(Rd);
figure
stem([0:49],r(1:50))
title('Covariance function')
xlabel('Time lag')

% Simulate AR(1)-process
n = 400;
m = 0;
sigma = 1;
e = normrnd(m, sigma, 1, n);
x = filter(C, A, e);
figure
plot(x)
title('Realisation')
xlabel('Sample')

%% 2.2 Study of an ARMA(2)-process
close all
clear all
clc

% Setup
A = [1 -1 0.5];
C = [1 1 0.5];

% Look at pole-zero plot
P = roots(A);
Z = roots(C);
figure
zplane(Z,P)

% Spectrum (for verification)
[H,w] = freqz(C,A);
R = abs(H).^2;
figure
plot(w/2/pi,R)
title('Spectral density')
xlabel('Frequency')

%% 3
close all force
clear all
clc

%armagui

%% 4.1
close all force
clear all
clc

% Audio file
[x,Fs] = audioread('fa.wav');
%sound(x,Fs)

n = length(x);                 
t = (0:n-1)/Fs;
figure
plot(t,x)
title('Speech signal')
xlabel('Time (s)')

dt = 20*10^(-3);                    % Set length of time-windows (20 ms)
dn = Fs*dt;                         % Sample length of each time-window
N_sec = floor(n/dn);                % Number of sections with length 20ms

M = 20;                             % Set AR model order
ar_param = zeros(N_sec,M+1);

for i = 1:N_sec                           % For each time frame:
    x_temp = x((i-1)*dn+1:i*dn);          % Pick out the i:th  20 ms section
    
    % Student code
    U = zeros(size(x_temp,1)-20,20);
    for k = 1:20
        U(:,k) = -x_temp(21-k:end-k);
    end
    ar_temp = (U'*U) \ U'*x_temp(21:end);
    Q = (x_temp(21:end) - U*ar_temp)' * (x_temp(21:end) - U*ar_temp);
    e_temp = Q / size(U,1);
    
    ar_param(i,:) = [ar_temp ; e_temp];     % and save in ar_param. Why do we not need to save the first one?
    
end

whos x ar_param

%% 4.2

% Calculate spectrum for some 20 ms
i = floor(N_sec/2);                % Pick the middle section.
x_temp = x((i-1)*dn+1:i*dn);
x_rec_temp = filter(1,[1 ar_param(i,1:end-1)],sqrt(ar_param(i,end))*randn(dn,1)); 
N_fft = 1024;

f = (0:N_fft-1)/N_fft;
Px = abs(fft(x_temp,N_fft)).^2/N_sec ;

w  = exp(2i*pi*f);
Pa = (ar_param(i,end)) ./abs( polyval([1 ar_param(i,1:end-1)],w).' ).^2;

figure
semilogy(f,Px);
hold on
semilogy(f,Pa,'r');
legend('Speech','AR Reconstruction');

%% 4.3

% Recreate the audio in each time frame using the AR parameters
x_rec = zeros(n,1);
for jj = 1:N_sec
   x_rec((jj-1)*dn+1:jj*dn) = filter(1,[1 ar_param(jj,1:end-1)],sqrt(ar_param(jj,end))*randn(dn,1)); 
end

figure
subplot(211);
plot(t,x);
title('Original sound')
xlabel('Time (s)')
subplot(212);
plot(t,x_rec);
title('Reconstructed sound')
xlabel('Time (s)')

% If stable, play sound
%sound(x_rec,Fs)