clear;
clc;
load ekg1
%%
size(ekg1)
figure
plot(ekg1)
%%
FT=1000;
N=10000;
n=0:N-1;
t=n/FT;
plot(t,ekg1)
xlabel('tid (s)')
ylabel('amplitud (mikroV)')
%%
M=10000;
ekgspek=fft(ekg1,M);
%%
figure
plot(abs(ekgspek))
%%
f=(0:M-1)/M*FT;
plot(f,abs(ekgspek))
xlabel('frequency (Hz)')
axis([0 20 0 150000])
%%
figure
subplot(211)
plot(t,ekg1)
h=1/15*ones(1,15);
y=conv(ekg1,h);
y=y(15:length(y));
subplot(212)
plot(t,y)
%%
ekgspek_y=fft(y,M);
figure
f=(0:M-1)/M*FT;
plot(f,abs(ekgspek))
xlabel('frequency (Hz)')
axis([0 20 0 150000])
%%
load man1
load kvinna1
%%
FT=8000;
sound(man1)
sound(kvinna1)
figure
subplot(211)
plot(man1)
subplot(212)
plot(kvinna1)
%%
e_man1=man1(10401:11200);
i_man1=man1(8201:9000);
e_kvi1=kvinna1(6201:7000);
o_kvi1=kvinna1(3801:4600);
%%
sound(e_man1,FT)
%%
sound(i_man1,FT)
%%
sound(e_kvi1,FT)
%%
sound(o_kvi1,FT)
%%
figure
N=800;
subplot(221)
plot((0:N-1)/N*FT,abs(fft(e_man1,N)))
axis([0 1500 0 300])
subplot(222)
plot((0:N-1)/N*FT,abs(fft(i_man1,N)))
axis([0 1500 0 300])
subplot(223)
plot((0:N-1)/N*FT,abs(fft(e_kvi1,N)))
axis([0 1500 0 300])
subplot(224)
plot((0:N-1)/N*FT,abs(fft(o_kvi1,N)))
axis([0 1500 0 300])
%%
lowpass=[1 1 1 1];
highpass=[1 -1 1 -1];
denominator=1; %Ej IIR
man1_l=filter(lowpass,denominator,man1);
man1_h=filter(highpass,denominator,man1);
kvinna1_l=filter(lowpass,denominator,kvinna1);
kvinna1_h=filter(highpass,denominator,kvinna1);
%%
sound([man1;man1_l;man1_h])
sound([kvinna1;kvinna1_l;kvinna1_h])
%%
ekofilter=[1 zeros(1,250) 0.8 zeros(1,250) 0.6];
%%
klang_kvinna1=filter(ekofilter,1,kvinna1);
sound([kvinna1;klang_kvinna1])
%%
clear;
clc;
load drugsignals.mat
%%
figure;
subplot(511); plot(ECG1)
subplot(512); plot(ECG2)
subplot(513); plot(ECG3)
subplot(514); plot(ECG4)
subplot(515); plot(ECG5)
%%
fs = 1000;
M = 1000;
f=(0:M-1)/M*fs;
SP_ECG1=fft(ECG1,M);
SP_ECG2=fft(ECG2,M);
SP_ECG3=fft(ECG3,M);
SP_ECG4=fft(ECG4,M);
SP_ECG5=fft(ECG5,M);
figure;
hold on;
plot(f, abs(SP_ECG1), 'LineWidth', 1.5, 'DisplayName','1')
plot(f, abs(SP_ECG2), 'LineWidth', 1.5, 'DisplayName','2')
plot(f, abs(SP_ECG3), 'LineWidth', 1.5, 'DisplayName','3')
plot(f, abs(SP_ECG4), 'LineWidth', 1.5, 'DisplayName','4')
plot(f, abs(SP_ECG5), 'LineWidth', 1.5, 'DisplayName','5')
legend
axis([0 100 0 15000])
hold off;
%%
clear;
clc;
load drugresponse.mat
%%
plot(drugresponse)
%%
h=0.1*ones(1,10);
hzeros=roots(h);
figure;
zplane(hzeros);
%%
freqz(h)
%%
[H,W] = freqz(h);
figure;
subplot(211);
plot(W/(2*pi),abs(H));
subplot(212);
plot(W/(2*pi),angle(H));
%%
figure
stem(abs(fft(h)))
%%
figure
subplot(211)
plot(abs(fft([h zeros(1,10)])))
subplot(212)
plot(abs(fft([h zeros(1,100)])))
%%
y = filter(h,1, drugresponse);
figure()
hold on
plot(drugresponse)
hold on
plot(y)
%%
h2 = 0.1 * [1, 1, 1, 1, 1, -1, -1, -1, -1, -1];
y2 = filter(h2,1, drugresponse);
plot(y2)