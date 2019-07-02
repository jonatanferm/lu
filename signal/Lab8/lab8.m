clear;
clc;
load EEGOnsetSeizure
%%
figure;
plot(X(:,3));
%%
eeg = X(1001:4100,3);
eeg = eeg - mean(eeg);
%%
Fs=256;
STFT=[];
Overlap=0.5;
N=128;
M = 1 + floor((length(eeg)-N)/(N-floor(N*Overlap)));
w=hamming(N); % Definiera ett Hammingf¨onster av l¨angd N
for i=1:M
    STFT = [STFT abs(fft(w.*eeg(1+(i-1)*N*Overlap:N+(i-1)*N*Overlap),N))];
end;
%%
figure
subplot(211)
t = (0:length(eeg)-1)/Fs;
plot(t,eeg)
xlim([0 max(t)]);
subplot(212)
t = (0:M-1)/M*length(eeg)/Fs;
f = (0:N-1)/N*Fs;
imagesc(t,f,STFT)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(gca,'ylim',[0 25])
% Den verkar inte varierande frekvensinnehåll
%%
Overlap=0.1;
N=256;
WV=WVD(eeg',N,1,ceil(Overlap*N),1,0);
[N,M]=size(WV);
figure
subplot(211)
t = (0:length(eeg)-1)/Fs;
plot(t,eeg)
xlim([0 max(t)]);
subplot(212)
t = (0:M-1)/M*length(eeg)/Fs;
f = (0:N-1)/N*Fs/2;
imagesc(t,f,abs(WV))
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(gca,'ylim',[0 25])
% betydligt högre upplösning
%%
delfin=audioread('delfin1.wav');
%%
Fs=44100;
STFT=[];
Overlap=0.5;
%N=1024;
N = 10;
M = floor(length(delfin)/N/Overlap-1); % Hur m˚anga f¨onster f˚ar plats i signalen
w=hamming(N); % Definiera ett Hammingf¨onster av l¨angd N
for i=1:M
    STFT = [STFT abs(fft(w.* delfin(1+(i-1)*N*Overlap:N+(i-1)*N*Overlap),N))];
end;
figure
subplot(211)
t = (0:length(delfin)-1)/Fs;
plot(t, delfin)
xlim([0 max(t)]);
subplot(212)
t = (0:M-1)/M*length(delfin)/Fs;
f = (0:N-1)/N*Fs;
imagesc(t,f,STFT)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(gca,'ylim',[0 15000])
%%

Overlap=0.5;
N= 1024;
WV=WVD(delfin',N,1,ceil(Overlap*N),1,0);
[N,M]=size(WV);
figure
subplot(211)
t = (0:length(delfin)-1)/Fs;
plot(t, delfin)
xlim([0 max(t)]);
subplot(212)
t = (0:M-1)/M*length(delfin)/Fs;
f = (0:N-1)/N*Fs/2;
imagesc(t,f,abs(WV))
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(gca,'ylim',[0 15000])
%%
Overlap=0.5;
N= 1024;
sigma = 1/1e9;
[WV,WV_smooth,A,A_smooth,k]=WVD(delfin',N,1,ceil(Overlap*N),1,sigma);
[N,M]=size(WV);
figure;
subplot(231);
imagesc(t,f,abs(WV));
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(gca,'ylim',[0 15000])
set(gca,'clim',[0 10])
title('Wigner-Ville');
subplot(233);
imagesc(t,f,abs(WV_smooth));
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(gca,'ylim',[0 15000])
set(gca,'clim',[0 10])
title('Wigner-Ville smoothed');
subplot(234);
imagesc(abs(A));
title('Ambiguity');
subplot(235);
imagesc(k);
title('Kernel');
subplot(236);
imagesc(abs(A_smooth));
title('Ambiguity x Kernel');