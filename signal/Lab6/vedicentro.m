% image acquisition
a=imread('37_7.bmp');
% rotate image
a=imrotate(a,0);
% write to a BitMap image the rotated mage
imwrite(a,'pippo.bmp');

% if image sizes are not 8*K where K is an integer
% image is resized
img=double(a);
imgN=size(img,1);
imgM=size(img,2);
modN=mod(imgN,8);
modM=mod(imgM,8);
img=img(modN+1:imgN,modM+1:imgM);

% call "centralizing" function
[out,xc,yc]=centralizing(img,0);

% show image and the center point determinated by
% the previous function
imshow(a);
hold on;
plot(xc,yc,'O');
hold off;