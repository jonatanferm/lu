% VL_DEMO_SIFT_MATCH  Demo: SIFT: basic matching

pfx = fullfile(vl_root,'figures','demo') ;
randn('state',0) ;
rand('state',0) ;
figure(1) ; clf ;

% --------------------------------------------------------------------
%                                                    Create image pair
% --------------------------------------------------------------------
sz = [1000 1000];
Ia = imread('./data/coll_2/HE/1HE.jpg') ;
Ia = im2single(rgb2gray(Ia));
t = Tiff('./data/coll_2/TRF/1TRF.tif','r');
Ib = read(t);
Ib = histeq(im2single(Ib)); %single(rgb2gray(Y));
Ib = 1 - Ib;
imwrite(Ib, './1TRF.png');
Ia = imresize(Ia, sz);
Ib = imresize(Ib, sz);

Ia = imgaussfilt(Ia, 3);
Ia = imsharpen(Ia,'Radius',2,'Amount',2);

Ib = imgaussfilt(Ib, 3);
Ib = imsharpen(Ib,'Radius',2,'Amount',2);
%Ib = imresize(imread('./data/coll_2/TRF/1.1.bmp'), sz) ;

% --------------------------------------------------------------------
%                                           Extract features and match
% --------------------------------------------------------------------

[fa,da] = vl_sift(im2single(Ia)) ;
[fb,db] = vl_sift(im2single(Ib)) ;

[matches, scores] = vl_ubcmatch(da,db, 1.5) ;

[drop, perm] = sort(scores, 'descend') ;
matches = matches(:, perm) ;
scores  = scores(perm) ;

matches = matches(:,:);
scores = scores(:,:);


imshow(cat(2, Ia, Ib)) ;

xa = fa(1,matches(1,:)) ;
xb = fb(1,matches(2,:)) + size(Ia,2) ;
ya = fa(2,matches(1,:)) ;
yb = fb(2,matches(2,:)) ;

hold on ;
h = line([xa ; xb], [ya ; yb]) ;
set(h,'linewidth', 1, 'color', 'b') ;

vl_plotframe(fa(:,matches(1,:))) ;
fb(1,:) = fb(1,:) + size(Ia,2) ;
vl_plotframe(fb(:,matches(2,:))) ;
axis image off ;