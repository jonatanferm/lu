addpath("./code");
%%
%[im, info] = mydicomreadfolder("images/MR-thorax-transversal");
[im, info] = mydicomreadfolder("images/MR-carotid-coronal");
%[info, im] = mydicomread("images/MR-heart-single.dcm");
%[info, im] = mydicomread("images/CT-thorax-single.dcm");
%[info, im] = mydicomread("images/MR-knee-single.dcm");

im = im.*info.RescaleSlope + info.RescaleIntercept;
%%
imagesc(do_speed_image(im, 125, 125));
colorbar
%%
disp_image(fliplr(permute(im, [2 1])))
%%
sz = size(im);
disp_image(fliplr(squeeze(im(round(sz(1)/2),:,:))), [info.PixelSpacing(1) info.SliceThickness])
disp_image(squeeze(im(:,round(sz(2)/2),:))', [info.SliceThickness info.PixelSpacing(2)])
image = squeeze(im(:,:,round(sz(3)/2)))';
%imgrad = imgradient(image);
%imgrad = max(imgrad, [], [1 2]) - imgrad;
%disp_image(imgrad)
disp_image(image, info.PixelSpacing)
%%
im_m1 = fliplr(squeeze(max(im,[],1)));
im_m2 = squeeze(max(im,[],2));
im_m3 = squeeze(max(im,[],3));
disp_image(im_m1, [info.PixelSpacing(1) info.SliceThickness]);
disp_image(im_m2', [info.SliceThickness info.PixelSpacing(2)]);
disp_image(im_m3', info.PixelSpacing);
%%

%%
interactive
%%
[i1 i2] = mydicomread("images/MR-knee-single.dcm");

%i2 = imgaussfilt(i2, 2);
%i2 = wiener2(i2, [3 3]);
%i2(1,:) = 0;

%i2 = abs((i2(120, 120) - i2));
%[~, thresh] = edge(i2, 'log');
%i2 = double(edge(i2, 'log', thresh*0.8));

%i2 = imgradientxy(i2);
%i2 = imgradientxy(i2);
%i2 = imgradientxy(i2);
%i2 = imgradientxy(i2);
%i2 = stdfilt(i2);
i2 = stdfilt(i2, ones(9));
i2 = rangefilt(i2);
i2 = imgaussfilt(i2, 3);
disp_image(i2);
%[info, im] = mydicomread("./images/MR-heart-single.dcm");
%[info, im] = mydicomread("./images/CT-thorax-single.dcm");
%%
imagesc(do_speed_image(im, 125, 125));
colorbar
%%
function sim = do_speed_image(image, x, y)
    fi = image;
    [Gx,Gy] = imgradientxy(fi);
    [Gmag,Gdir] = imgradient(Gx,Gy);
    fi = Gmag ./ stdfilt(Gdir, true(3)).^2;
    fi = min(fi, max(fi(:))./20);
    fi = max(fi, max(fi(:))./100);
    fi = fi-min(fi(:));

    xcs = zeros(size(fi));
    xcs(1:x,:) = cumsum(fi(1:x,:), 1, 'reverse');
    xcs(x:end,:) = cumsum(fi(x:end,:), 1);
    ycs = zeros(size(fi));
    ycs(:,1:y) = cumsum(fi(:,1:y), 2, 'reverse');
    ycs(:,y:end) = cumsum(fi(:,y:end), 2);
    
    sim = max(xcs, ycs);
    sim = imbinarize(sim);
    sim = imerode(sim, strel('disk', 2));
    sim = imclose(sim, strel('disk', 10));
    sim = imclose(sim, strel('disk', 8));
    sim = 1-sim;
end

%%
function [] = disp_image(image, aspect)
    figure;
    imagesc(image)
    sz = size(image);
    set(gca,'xcolor',[1 1 1],'ycolor',[1 1 1]);
    set(gca, 'Xgrid', true);
    set(gca, 'Ygrid', true);
    set(gca,'ticklength',[0.01 0.01]);
    set(gca,'XTick', 0:(10/aspect(2)):sz(2));
    set(gca,'YTick', 0:(10/aspect(1)):sz(1));
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    colorbar;
    pbaspect(fliplr([1 size(image).*aspect]))
    colormap(gray)
    sz.*aspect
end