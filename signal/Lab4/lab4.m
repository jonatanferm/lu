clear;
clc;
load senddata1
load kanalfilter1
load ikanalfilter1
%%
figure
subplot(311)
stem(kanalfilter1)
subplot(312)
stem(ikanalfilter1)
subplot(313)
stem(conv(kanalfilter1,ikanalfilter1))
%%
figure
subplot(311)
stem(senddata1)
mottagen=conv(senddata1,kanalfilter1);
subplot(312)
stem(mottagen)
subplot(313)
avkodad=conv(mottagen,ikanalfilter1);
stem(avkodad)
%%
clear;
clc;
load kanal1
load data2
%%
N=64;
n=0:N-1;
f=[3.2;4.3];
s=sqrt(2/N)*sin(2*pi*(f/N)*n)';
%%
figure
subplot(211)
plot(s)
title('Two waveforms')
subplot(212)
plot((0:N-1)/N,abs(fft(s,N)));
title('Spectrum')
xlabel('Normalized frequency')
axis([0 0.5 0 6])
%%
skalprod=s(:,1)'*s(:,2)
%%
f=[1;2;3;4;5;6;7;8;9;10];
s=sqrt(2/N)*sin(2*pi*(f/N)*n)';
%%
figure
subplot(211)
plot(s)
title('Ten waveforms')
subplot(212)
plot((0:N-1)/N,abs(fft(s,N)));
title('Spectrum')
xlabel('Normalized frequency')
axis([0 0.5 0 6])
%%
data1 = [5 4 3 2 1 1 2 3 4 5];
x=data1*s';
%%
figure
plot((0:N-1)/N,abs(fft(x,N)));
title('Spectrum')
xlabel('Normalized frequency')
axis([0 0.5 0 40])
%%
x=[x x];
%%
mottaget=conv(x,kanal1);
mottaget=mottaget(N+1:2*N);
%%
avkodaddata=mottaget*s;
%%
kanalskattning = avkodaddata./data1;
%%
figure
plot(kanalskattning)
xlabel('Waveform #')
axis([1 10 0 2])
%%
x=data2*s';
x=[x x];
mottaget=conv(x,kanal1);
mottaget=mottaget(N+1:2*N);
avkodaddata=mottaget*s;
%%
skattaddata2=avkodaddata./kanalskattning
data2