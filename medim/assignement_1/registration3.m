run('./vlfeat/toolbox/vl_setup');

pars = [];

%for file_name = dir('./data/coll_1/HE/*.bmp')'
%for file_name2 = ['4.7.bmp' '6.3.bmp' '8.0.bmp']'
file_name.name = '1.1.bmp'; 

X = imresize(imread(strcat('./data/coll_1/HE/', file_name.name)), 1);
Y = imresize(imread(strcat('./data/coll_1/p63AMACR/', file_name.name)), 1);

Xg = im2single(rgb2gray(X));
Yg = im2single(rgb2gray(Y));

%Xg = imgaussfilt(Xg, 2);
%Xg = imsharpen(Xg,'Radius',2,'Amount',5);
%Xg = im2single(histeq(Xg, 25));

%Yg = imgaussfilt(Yg, 2);
%Yg = imsharpen(Yg,'Radius',2,'Amount',5);
%Yg = im2single(histeq(Yg, 25));


[Fx, Dx] = vl_sift(Xg);
[Fy, Dy] = vl_sift(Yg);

Fxc = Fx(1:2,:);
Fyc = Fy(1:2,:);

[m, s] = vl_ubcmatch(Dx, Dy, 1.3);

mx = Fxc(:,m(1,:));
my = Fyc(:,m(2,:));
clear Fyc Fxc
%imshow(insertMarker(Xg, mx'));
%figure; imshow(insertMarker(Yg, my'));

best_indices = find_best_indices(mx, my, false);
ax = mx(:,best_indices);
ay = my(:,best_indices);
save(strcat('pro_', file_name.name), 'ax', 'ay');
[d, Z, t] = procrustes(mx(:,best_indices)',my(:,best_indices)');
size(best_indices, 1)
if t.b < 0.01
  t.b = 1;
end
M = [[t.b * t.T [0; 0]]; [t.c(1,:) 1]]; % get the transformation matrix
mt = affine2d(M);
%Ty = imwarp(Y, mt); % transform image
Ty = imwarp(histeq(im2single(Y)), mt, 'OutputView', imref2d(size(X)));

%%transform points
imsize = [800, 800];
text_str = strcat(strcat('rotation:', num2str((acos(t.T(1,1)) / (2*pi)) * 360, 4)) ...
    , strcat('translation:',num2str(t.c(1,:))));
final_image = ...
    [imresize(X, imsize), ...
    imresize(Y, imsize), ... 
    imfuse(imresize(Ty, imsize, 'bicubic'), imresize(X, imsize), 'blend')];
final_image = insertText(final_image,[0 0],text_str,'FontSize',32,'BoxColor',...
    'black','BoxOpacity',0.4,'TextColor','white');
%imwrite(final_image, strcat(extractBefore(file_name.name, '.bmp'), '.png'));
imshow(final_image)
%figure;imshow(final_image)

%end
