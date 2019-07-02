filename = '412.mat'
X = imresize(imread('./data/coll_2/HE/12HE.jpg'), 1);
%Y = imresize(imread('./data/coll_1/p63AMACR/6.0.bmp'), 1);

t = Tiff('./data/coll_2/TRF/12TRF.tif','r');
Y = read(t);

Xg = im2single(rgb2gray(X));
%Yg = im2single(rgb2gray(Y));
Yg = histeq(im2single(Y)); %single(rgb2gray(Y));
%Yg = 1 - Yg;

figure;
ax1 = subplot(121);
imshow(Xg);
hold on;
ax2 = subplot(122);
imshow(Yg);
hold on;

npoints = 5;

[X1 X2] = getpts(ax1);
[Y1 Y2] = getpts(ax2);

xps = [X1 X2];
yps = [Y1 Y2];

save(filename, 'xps', 'yps');


[d, Z, t] = procrustes(xps,yps);
M = [[t.b * t.T [0; 0]]; [t.c(1,:) 1]]; % get the transformation matrix
mt = affine2d(M);
%Ty = imwarp(Y, mt); % transform image
Ty = imwarp(histeq(im2single(Y)), mt, 'OutputView', imref2d(size(X)));
imsize = [2*1024, 2*1024];
fused = imfuse(imresize(Ty, imsize, 'bicubic'), imresize(X, imsize), 'blend');
figure; imshow(fused)