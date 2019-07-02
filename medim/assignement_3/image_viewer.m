load(fullfile('databases','hep_proper_mask'));
X1_masks = Y1;
X2_masks = Y2;
load(fullfile('databases','hep_proper'));
%%
i = 15;
im = X1(:,:,1,i);
%im = im/255;
%im(find(~X1_masks(:,:,1,i))) = 0;
figure;
imagesc(im);
%imhist(im, 16);