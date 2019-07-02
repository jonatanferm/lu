run('./vlfeat/toolbox/vl_setup');


%for file_name = dir('./data/coll_1/HE/*.bmp')'

%file_name.name = '8.4.bmp';
%X = imresize(imread(strcat('./data/coll_1/HE/', file_name.name)), 1);
%Y = imresize(imread(strcat('./data/coll_1/p63AMACR/', file_name.name)), 1);

%file_name.name = '8.3.bmp';
file_names = [];
for i = [12]
X = imread(strcat(strcat('./data/coll_2/HE/', num2str(i)), 'HE.jpg'));
t = Tiff(strcat(strcat('./data/coll_2/TRF/', num2str(i)), 'TRF.tif'),'r');
Y = read(t);

%Y = imresize(imread('./data/coll_2/TRF/1TRF.tif'), 1);

%X = imcrop(X, [max(size(Y, 2)/2 - size(X, 2)/2, 0), max(size(Y, 1)/2 - size(X, 1)/2, 0), size(X, 2)/2 + size(Y, 2)/2, size(X, 1)/2 + size(Y, 1)/2]);
%Y = imcrop(Y, [max(size(Y, 2)/2 - size(X, 2)/2, 0), max(size(Y, 1)/2 - size(X, 1)/2, 0), size(X, 2)/2 + size(Y, 2)/2, size(X, 1)/2 + size(Y, 1)/2]);

Xg = im2single(rgb2gray(X));
%Yg = single(rgb2gray(Y));
Yg = histeq(im2single(Y)); %single(rgb2gray(Y));
Yg = 1 - Yg;

Xg = imgaussfilt(Xg, 3);
Xg = imsharpen(Xg,'Radius',2,'Amount',2);

Yg = imgaussfilt(Yg, 3);
Yg = imsharpen(Yg,'Radius',2,'Amount',2);

[Fx, Dx] = vl_sift(Xg);
[Fy, Dy] = vl_sift(Yg);

Fxc = Fx(1:2,:);
Fyc = Fy(1:2,:);

[m, s] = vl_ubcmatch(Dx, Dy, 1.5);

mx = Fxc(:,m(1,:));
my = Fyc(:,m(2,:));
clear Fyc Fxc

best_indices = find_best_indices(mx, my, true);
[d, Z, t] = procrustes(mx(:,best_indices)',my(:,best_indices)');
size(best_indices, 1)
if t.b < 0.01
  t.b = 1;
end
M = [[t.b * t.T [0; 0]]; [t.c(1,:) 1]]; % get the transformation matrix
mt = affine2d(M);
ax = mx(:,best_indices);
ay = my(:,best_indices);
save(strcat(strcat('pro2_', num2str(i)), '.mat'), 'ax', 'ay');
%Ty = imwarp(Y, mt); % transform image
Ty = imwarp(histeq(im2single(Y)), mt, 'OutputView', imref2d(size(X)));

%%transform points
imsize = [2*1024, 2*1024];
text_str = strcat(strcat(strcat('rotation:', num2str((acos(t.T(1,1)) / (2*pi)) * 360, 4)) ...
    , strcat('translation:',num2str(t.c(1,:)))), ...
    strcat('scale:',num2str(t.b)));

final_image = imfuse(Ty, X, 'blend');
final_image = insertText(final_image,[0 0],text_str,'FontSize',32,'BoxColor',...
    'black','BoxOpacity',0.4,'TextColor','white');
%imwrite(final_image, strcat('c2_',strcat(num2str(i),'.png')));
figure; imshow(final_image)
end

%end
