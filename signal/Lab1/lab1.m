%% 1.2
[Y,FS] = audioread('delfin1.wav');
%[Y,FS] = audioread('jamesbond.wav');
%%
soundsc(Y,FS);
%%
Y2=decimate(Y,2);
%%
soundsc(Y2,FS/2)
%%
Y3=Y(1:2:end);
soundsc(Y3,FS/2)
%%
Y2=decimate(Y,3);
Y3=Y(1:3:end);
%%
soundsc([Y2;Y3],FS/3)
%% 1.3
load piano1
load shout1
load drum1
load hi1
%%
size(piano1)
size(shout1)
size(drum1)
size(hi1)
%%
figure(1)
plot(drum1)
figure(2)
plot(shout1)
figure(3)
plot(drum1)
figure(4)
plot(hi1)
%%
FT = 44100;
C4=261.63;
A3=220;
F3=174.61;
G3=196;
%%
Ntakt=44100;
n=(0:Ntakt-1)';
%%
ton1=sin(2*pi*C4/FT*n);
ton2=sin(2*pi*A3/FT*n);
ton3=sin(2*pi*F3/FT*n);
ton4=sin(2*pi*G3/FT*n);
%%
bakgrundsslinga=[ton1;ton2;ton3;ton4];
sound(bakgrundsslinga,44100);
%%
figure(5);
plot(bakgrundsslinga);
%%
dft = fft(bakgrundsslinga);
dft = dft(1:length(dft)/2);
plot(abs(dft))
xlim([400 1200])
%%
komptakt=[drum1;hi1;hi1;hi1];
komp=[komptakt;komptakt;komptakt;komptakt];
sound(komp,44100);
%%
sound(0.5*bakgrundsslinga+komp,44100);
%%
piano2=resample(piano1,6,5);
piano3=resample(piano1,6,4);
piano4=resample(piano1,4,3);
%%
pianoC4=[piano1(1:5512);piano1(1:5513);piano1(1:5512);piano1(1:5513);
piano1(1:5512);piano1(1:5513);piano1(1:5512);piano1(1:5513)];
pianoA3=[piano2(1:5512);piano2(1:5513);piano2(1:5512);piano2(1:5513);
piano2(1:5512);piano2(1:5513);piano2(1:5512);piano2(1:5513)];
pianoF3=[piano3(1:5512);piano3(1:5513);piano3(1:5512);piano3(1:5513);
piano3(1:5512);piano3(1:5513);piano3(1:5512);piano3(1:5513)];
pianoG3=[piano4(1:5512);piano4(1:5513);piano4(1:5512);piano4(1:5513);
piano4(1:5512);piano4(1:5513);piano4(1:5512);piano4(1:5513)];
pianoslinga=[pianoC4;pianoA3;pianoF3;pianoG3];
%%
sound(pianoslinga,44100);
%%
sound(0.5*bakgrundsslinga+komp+pianoslinga,44100);
%%
shoutslinga=[zeros(size(shout1));shout1;zeros(size(shout1));shout1];
sound(shoutslinga,44100);
%%
sound(0.5*bakgrundsslinga+komp+pianoslinga+shoutslinga,44100);
%%
load bild1
figure
imagesc(bild1)
colormap('gray')
%%
figure
subplot(211)
plot(bild1(4:144,110))
h=[1 1 1 1 1 1 1 1];
kolonn110filt=filter(h,1,bild1(4:144,110));
subplot(212)
plot(kolonn110filt)
%%
figure
subplot(221)
imagesc(bild1)
colormap('gray')
bild2=filter(h,1,bild1);
subplot(222)
imagesc(bild2)
bild3=filter(h,1,bild1)';
subplot(223)
imagesc(bild3)
h2 = h'*h;
bild4=filter2(h2,bild1);
subplot(224)
imagesc(bild4)
%%
h=[1 -1];
h2 = h'*h;
bild5=filter2(h2,bild1);
figure
imagesc(abs(bild5))
colormap('gray')