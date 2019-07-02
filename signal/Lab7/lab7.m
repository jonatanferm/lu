%% 7.1
x1 = [1 1 1 1;
      0 0 0 0;
      1 1 1 1;
      0 0 0 0];
x2 = [1 0 1 0;
      1 0 1 0;
      1 0 1 0;
      1 0 1 0];

fft2(x1)
fft2(x2)
%%
x = sin(2*pi*0.25*[1:40])';
X = repmat(x,[1 40]);
figure;
subplot(211);
surf(X);
colormap('gray');
subplot(212);
imagesc(X);
colormap('gray');
%%
[a,b] = size(X);
fx = (0:a-1)/a;
fy = (0:b-1)/b;
figure;
imagesc(fx,fy,abs(fft2(X)));
colormap('gray');
xlabel('f_x');
ylabel('f_y');
%%
T = toeplitz(x);
figure;
subplot(211);
imagesc(T);
colormap('gray');
subplot(212);
imagesc(fx,fy,abs(fft2(T)));
colormap('gray');
%%
figure;
subplot(211);
imagesc(X+X');
colormap('gray');
subplot(212);
imagesc(fx,fy,abs(fft2(X+X')));
colormap('gray');
%%
bild=imread('mars','bmp');
bild=double(bild);
set(0,'DefaultFigureColormap',feval('gray'));
figure;
imagesc(bild);
figure;
surf(abs(fft2(bild)))
%%
bild_mv=filter2(1/9*ones(3),bild);
figure;
imagesc(bild_mv)
figure;
surf(abs(fft2(bild_mv)))
%% 
brusigbild=bild+mean(mean(bild))*0.1*randn(size(bild));
bild_mv=filter2(1/9*ones(3),brusigbild);
%bild_mv2=wiener2(brusigbild, [5 5]);
figure;
imagesc(bild_mv);
figure;
imagesc(brusigbild);
figure;
imagesc(bild);
%%
% Frekvensdom¨an l˚agpassfiltrering med cirkul¨ar mask;
N=length(brusigbild);
[f1,f2]=freqspace(N);
[f1,f2]=meshgrid(f1,f2);
cutoff=.5;
Hlp=zeros(N,N);d=find(sqrt(f1.^2+f2.^2) < cutoff);
Hlp(d)=ones(size(d));
Filtered=fft2(brusigbild).*fftshift(Hlp);
filtered=abs(ifft2(Filtered));
imagesc(filtered)
%Kanske lite tråkigt att höga frekvenser från den riktiga bilden försvinner
%med
%%
h_x=[-1 -2 -1;0 0 0;1 2 1];
h_y=[-1 0 1;-2 0 2;-1 0 1];
filtered_x=filter2(h_x,bild);
filtered_y=filter2(h_y,bild);
filtered=sqrt(filtered_x.^2+filtered_y.^2);
imagesc(filtered)
% Tycker den gör det ganska bra
%%
imagesc(edge(bild))
%%
bild = 10*ones(13);
bild(5:9,5:9) = 200;
imagesc(bild);
%%
% Atm. turb. 5x5
bx=20;
by=20;
K=10;
h_AT=K*exp(-(-2:2).^2./(2*bx^2))'*exp(-(-2:2).^2./(2*by^2));
%%
% Fel. kam. 5x5
ax=0.2;
ay=0.5;
K=10;
h_FK=K*sinc(ax*(-2:2))'*sinc(ay*(-2:2));
%%
% Spatiell filtrering av testbild
filtered_AT=filter2(h_AT,bild);
filtered_FK=filter2(h_FK,bild);
figure;
imagesc(filtered_AT);
figure;
imagesc(filtered_FK);
%%
% Filtrering (i frekvensplanet) av testbild och inversfiltrering
h_AT13=zeros(13);
h_AT13(5:9,5:9)=h_AT; % Lika stort filter
filtered=abs(ifft2(fft2(bild).*fft2(h_AT13)));
invfiltered=abs(ifft2(fft2(filtered)./fft2(h_AT13)));
imagesc(invfiltered);
%%
h_AT13=zeros(13);
h_AT13(5:9,5:9)=h_FK; % Lika stort filter
filtered=abs(ifft2(fft2(bild).*fft2(h_AT13)));
invfiltered=abs(ifft2(fft2(filtered)./fft2(h_AT13)));
imagesc(invfiltered);
%%
clear all
bild=imread('fran','bmp');
bild=double(bild);
%%
bx=20;
by=20;
K=10;
h_AT=K*exp(-(-12.5:1:12.5).^2./(2*bx^2))'*exp(-(-12.5:1:12.5).^2./(2*by^2));
h_AT256=zeros(256);
h_AT256(116:141,116:141)=h_AT;
ax=0.2;
ay=0.2;
K=10;
h_FK=K*sinc(ax*(-12.5:1:12.5))'*sinc(ay*(-12.5:1:12.5));
h_FK256=zeros(256);
h_FK256(116:141,116:141)=h_FK;
filtered=abs(ifft2(fft2(bild).*fft2(h_AT256).*fft2(h_FK256)));
brus=mean(mean(bild))*250000;
brusfiltered=filtered+brus*randn(size(bild));
figure
subplot(131)
imagesc(bild)
colormap('gray')
subplot(132)
imagesc(filtered)
subplot(133)
imagesc(brusfiltered)
% A den blir ju väldigt distorted
%%
Hinv=fft2(h_AT256).*fft2(h_FK256);
invfiltered=abs(ifft2(fft2(brusfiltered).*conj(Hinv)./ ...
                (Hinv.*conj(Hinv)+brus^2./(bild.^2))));
imagesc(invfiltered)
% Ganska långt ifrån fortfarance
%%
bild1=double(imread('m100blurred','bmp'));
bild2=double(imread('m100sharp','bmp'));
bx=2;
by=2;
K=10;
h_AT=K*exp(-(-22.5:1:22.5).^2./(2*bx^2))'*exp(-(-22.5:1:22.5).^2./(2*by^2));
h_AT256=zeros(256);
h_AT256(106:151,106:151)=h_AT;
invfiltered=abs(ifft2(fft2(bild1)./(fft2(h_AT256)+50)));
figure;
subplot(131)
imagesc(bild1);
colormap('gray')
title('Original');
subplot(132)
imagesc(abs(invfiltered))
title('Inverse filtered');
subplot(133)
imagesc(bild2)
% Sfäriska abberationen är inte samma problem som innan?
title('Sharp');