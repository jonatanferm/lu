clc;
clear;
addpath("./code");
%%
[im, info] = mydicomreadfolder("images/MR-carotid-coronal");
%%
%[~,i] = max(im(:));
%[x, y, z] = ind2sub(size(im),i);
x = 51;
y = 129;
z = 31;
fi = im;
%{
fi = imgaussfilt3(fi);
[Gmag,Gazimuth,Gelevation] = imgradient3(im);
gs = Gmag ./ (stdfilt(Gazimuth) + stdfilt(Gelevation));
gs(isnan(gs)) = max(gs(:));
gs = min(gs, max(gs(:))./3);
gs = max(gs, max(gs(:))./100);
gs = gs-min(gs(:));
gs = (gs).^2;

gs = gs./max(gs(:));
gs = 1-gs;
%}

speed = abs(fi-fi(x, y, z));
speed = speed.^2;
speed = speed./max(speed(:));
speed = 1-speed;



ts = speed; %+ 0.5*gs; 
ts(ts<1e-2) = 1e-2;
%%
T=msfm3d(ts, [x;y;z], true, true);

%msfm2d(GUI.SPEED,[GUI.YSeed;GUI.XSeed],true,true)


%%
r = -1:0.1:1;
[a, b, c] = meshgrid(r, r, r);
foo = a.^2 + b.^2 + c.^2;

%%
im3d = T;
[x,y,z] = size(im3d);
[X,Y,Z] = meshgrid(1:x,1:y,1:z);
is = isosurface(X,Y,Z,im3d,50);
p = patch(is);
isonormals(X,Y,Z,im3d,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1 1 1])
view(3); 
axis tight
camlight 
lighting gouraud
%%
