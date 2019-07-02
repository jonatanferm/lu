xn = 0:5;
n = 0:5;
k = 0:5;
N = length(xn);

Xk = @(k) sum(xn .* exp(-pi*2i/N).*k'.*n, 2);
Vk = @(k) 2*sum(xn .* cos(pi/N * (n + .5) .* k' ), 2);
%%
Xk(k)
Vk(k)
%%
clear all
f=imread('man256','bmp');
f=double(f);
figure
imagesc(f)
colormap('gray')
%%
DCTMat=dct_basis(8);
coeffsC=pic_coeff(f,DCTMat,8);
coeffsC=coeff_mask(coeffsC,5,8);
bild_recC=pic_reconstr(coeffsC,DCTMat,8);
lseC=sum(sum((f-bild_recC).^2))
%%
DFTMat=dft_basis(8);
coeffsF=pic_coeff(f,DFTMat,8);
coeffsF=coeff_mask(coeffsF,5,8);
bild_recF=pic_reconstr(coeffsF,DFTMat,8);
lseF=sum(sum((f-bild_recF).*conj(f-bild_recF)))
%%
figure
subplot(321)
imagesc(f)
set(gca,'clim',[0 255])
subplot(322)
stem(abs([lseC lseF]))
set(gca,'xtick',[1 2],'xticklabel',['DCT';'DFT'])
xlim([0 3]);
ylabel('Error');
colormap('gray')
subplot(323)
imagesc(coeffsC)
subplot(324)
imagesc(abs(coeffsF))
subplot(325)
imagesc(bild_recC)
title('C')
set(gca,'clim',[0 255])
subplot(326)
imagesc(abs(bild_recF))
title('F')
set(gca,'clim',[0 255])
%%
clear;
clc;
%%
block=[2 4 8 16];
cpu_start = cputime;
for p=1:length(block)
DCTMat=dct_basis(block(p));
coeffsC=pic_coeff(f,DCTMat,block(p));
coeffsC=coeff_mask(coeffsC,10,block(p));
bild_recC{p}=pic_reconstr(coeffsC,DCTMat,block(p));
cpu(p)= cputime;
lseC(p)=sum(sum((f-bild_recC{p}).^2));
end;
cpu = cpu - cpu_start;
%%
figure('Name', 'sq diff and cpu')
subplot(211)
plot(lseC)
set(gca,'xtick',[1 2 3 4],'xticklabel',[' 2';' 4';' 8';'16'], 'YScale', 'log')
subplot(212)
plot(cpu)
set(gca,'xtick',[1 2 3 4],'xticklabel',[' 2';' 4';' 8';'16'])
%%
figure('name', 'recons')
subplot(221)
imagesc(abs(bild_recC{1}))
set(gca,'clim',[0 255])
subplot(222)
imagesc(abs(bild_recC{2}))
set(gca,'clim',[0 255])
subplot(223)
imagesc(abs(bild_recC{3}))
set(gca,'clim',[0 255])
subplot(224)
imagesc(abs(bild_recC{4}))
set(gca,'clim',[0 255])
%%
FPlab
%%
fpextractdemo