%% Stationary stochastic processes, computer exercise 3
% 1.1 Question answers 
%
% 1) To receive spectral density values for larger amount of frequencies,
% i.e increase the resolution of the visualization.
%
% 2) The covariance function is r(tau)=E[A^2]/2 cos(2 pi f_0 tau),
% tau=0,+-1, +-2.... Using the periodogram will change the
% expected value of the estimated covariance function (and with the same 
% effect of each realization), to 
% r(tau)*(1-|tau|/n), tau=0,+-1,...+-(n-1), see p. 240.
%
% 3) A data window has the advantage of better sidelobe suppression, which
% implies less leakage and better estimation of the expected value for
% spectra with large dynamics. The drawback is the increased mainlobe
% width.
% 
% 4) The variance is V[Rhat_x(f)] approx R^2_x(f). The Welch method cuts
% the data into shorter sequences (non-overlapping or partly overlapping)
% which results in approximately un-correlated  spectrum estimates, 
% which are averaged for the final Welch spectrum estimate. The variance is
% reduced a factor K where K is the number of averaged spectrum estimates.
%
% 5) The window length of each window of the Welch method is the largest integer L<=2n/(K+1)
% 
% 6) R_XY(f)=H(f)R_XX(f), kappa^2(f)=1
%


%% 2 The periodogram and zero-padding
close all
%clear all
clc

load unknowndata.mat

% FFT without zero-padding and without subtracting the mean!
X = fft(data);
n = length(data);
Rhat = (X .* conj(X)) / n;
f = (0:n-1)/n;
figure
plot(f,Rhat)
xlabel('Normalised frequency')
title('Periodogram without zero-padding of the non-zero mean signal')

% FFT without zero-padding of the zero-mean signal

x=data-mean(data);

X= fft(x);
n = length(x);
Rhat = (X .* conj(X)) / n;
f = (0:n-1)/n;
figure
plot(f,Rhat)
xlabel('Normalised frequency')
title('Periodogram of the zero-mean signal without zero-padding')

% The use of matlab fftshift

figure
plot(f-0.5,fftshift(Rhat))
xlabel('Normalised frequency')
title('Periodogram using fftshift')

% FFT with zero-padding
NFFT = 512;
X = fft(x,NFFT);
Rhat = (X .* conj(X)) / n;
f = (0:NFFT-1)/NFFT;
figure
plot(f,Rhat)
xlabel('Normalised frequency')
title('Periodogram with zero-padding')

% The true frequency is f_0=0.2 



% Covariance function
rhat = ifft(Rhat);
figure
plot([0:length(rhat)-1],rhat)

%For explanation, see preparation exercise 2

%% 3 Investigation of different spectrum estimation techniques
close all force
%clear all
clc

% Define model and look at poles, zeros and spectral density
modell.A = [1 -2.39 3.35 -2.34 0.96];
modell.C = [1 0 1];
%armagui

% Visualize periodogram
NFFT = 4096;
e = randn(500,1);
x = filter(modell.C,modell.A,e);

% Periodogram with logarithmic scale
figure
periodogram(x,[],NFFT);

% Windowed periodogram using a Hanning window
figure
periodogram(x,hanning(500),NFFT);

% The windowed periodogram captures the zero at f=0.25, where the leakage of the
% periodogram causes severe bias of the zero. 


% Plotting the output in linear or logarithmic scale

[Rhat,f]=periodogram(x,[],NFFT,1);
figure
plot(f,Rhat)
figure
semilogy(f,Rhat)

% Welch method (window length L=90)
figure
pwelch(x,hanning(90),[],NFFT);
title('Welch method, window length 90, K=10')
xlabel('Normalised frequency')

Rhate=periodogram(e,[],NFFT);
Rhatew=pwelch(e,hanning(90),[],NFFT);

%The ratio of the variance of the periodogram and the Welch method (should
%be around K=10)

var(Rhate)./var(Rhatew)

%The Thomson estimator (you can tell them that his name is not Thompson as
%said in the title of the Matlab-produced plot!)

figure
pmtm(x,(10+1)/2,NFFT);





%% 4.1 Spectral estimates
close all force
clear all
clc

load threeprocessdata.mat

% Plot outputs
figure
subplot(311)
plot(y1)
title('y1')
subplot(312)
plot(y2)
title('y2')
subplot(313)
plot(y3)
title('y3')

% Plot covariance functions
L = 50;
r1 = covf(y1,L+1);
r2 = covf(y2,L+1);
r3 = covf(y3,L+1);
figure
subplot(131)
plot(0:L,r1)
title('Covariance funtion for y1')
xlabel('Lags')
subplot(132)
plot(0:L,r2)
title('Covariance funtion for y2')
xlabel('Lags')
subplot(133)
plot(0:L,r3)
title('Covariance funtion for y3')
xlabel('Lags')

% Periodograms
NFFT = 4096;
[Rhat1,f] = periodogram(y1,[],NFFT,1);
Rhat2 = periodogram(y2,[],NFFT);
Rhat3 = periodogram(y3,[],NFFT);
figure
subplot(131)
semilogy(f,Rhat1)
title('Periodogram for y1')
xlabel('Normalised frequency')
subplot(132)
semilogy(f,Rhat2)
title('Periodogram for y2')
xlabel('Normalised frequency')
subplot(133)
semilogy(f,Rhat3)
title('Periodogram for y3')
xlabel('Normalised frequency')

% Welch method

% K=10, L=90

[Rhat1w,f] = pwelch(y1,hanning(90),[],NFFT,1);
    Rhat2w = pwelch(y2,hanning(90),[],NFFT);
    Rhat3w = pwelch(y3,hanning(90),[],NFFT);

figure
subplot(131)
semilogy(f,Rhat1w)
title('Welch method for y1')
xlabel('Normalised frequency')
subplot(132)
semilogy(f,Rhat2w)
title('Welch method for y2')
xlabel('Normalised frequency')
subplot(133)
semilogy(f,Rhat3w)
title('Welch method for y3')
xlabel('Normalised frequency')

%% 4.2 Cross-spectrum and coherence spectrum
close all
clc

% Coherence spectra using Welch estimator
[Cxy1,f] = mscohere(x1,y1,hanning(90),[],NFFT,1);
Cxy3 = mscohere(x3,y3,hanning(90),[],NFFT);
figure
subplot(121)
plot(f,Cxy1)
title('Coh. sp. x1 and y1')
xlabel('Normalised frequency')

subplot(122)
plot(f,Cxy3)
title('Coh. sp. x3 and y3')
xlabel('Normalised frequency')

[Cxy1,f] = mscohere(x3,y1,hanning(90),[],NFFT,1);
Cxy3 = mscohere(x1,y3,hanning(90),[],NFFT);
figure
subplot(121)
plot(f,Cxy1)
title('Coh. sp. x3 and y1')
xlabel('Normalised frequency')

subplot(122)
plot(f,Cxy3)
title('Coh. sp. x1 and y3')
xlabel('Normalised frequency')


%% Extra!  Transfer/amplitude function

%Correcting labels of x1 and x3
x1ny=x3;
x3ny=x1;

Rxy1 = cpsd(x1ny,y1,hanning(90),[],NFFT);
Rxy2 = cpsd(x2,y2,hanning(90),[],NFFT);
Rxy3 = cpsd(x3ny,y3,hanning(90),[],NFFT);
Rxx1 = pwelch(x1ny,hanning(90),[],NFFT);
Rxx2 = pwelch(x2,hanning(90),[],NFFT);
Rxx3 = pwelch(x3ny,hanning(90),[],NFFT);
A1 = abs(Rxy1./Rxx1);
A2 = abs(Rxy2./Rxx2);
A3 = abs(Rxy3./Rxx3);
figure
subplot(131)
plot(f,20*log10(A1))
title('|H(f)|, system 1')
xlabel('Normalised frequency')
subplot(132)
plot(f,20*log10(A2))
title('|H(f)|, system 2')
xlabel('Normalised frequency')
subplot(133)
plot(f,20*log10(A3))
title('|H(f)|, system 3')
xlabel('Normalised frequency')

% Impulse response
H1 = tfestimate(x1new,y1,hanning(500),[],NFFT,'twosided');
H2 = tfestimate(x2,y2,hanning(500),[],NFFT,'twosided');
H3 = tfestimate(x3new,y3,hanning(500),[],NFFT,'twosided');
h1 = ifft(H1);
h2 = ifft(H2);
h3 = ifft(H3);
figure
subplot(131)
plot([0:length(h1)-1],h1)
axis([0 300 min(h1) max(h1)])
title('h(t), system 1')
xlabel('Samples')
subplot(132)
plot([0:length(h2)-1],h2)
axis([0 300 min(h2) max(h2)])
title('h(t), system 2')
xlabel('Samples')
subplot(133)
plot([0:length(h3)-1],h3)
axis([0 300 min(h3) max(h3)])
title('h(t), system 3')
xlabel('Samples')

% Redo calculation of transfer function for system 2
H2 = tfestimate(x2(1:end-64),y2(65:end),hanning(90),[],NFFT);
figure
plot(f,20*log10(abs(H2)))
title('Estimated amplitude function, system 2')
xlabel('Normalised frequency')