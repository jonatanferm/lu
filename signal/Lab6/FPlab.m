close all

%--------------------------------------
% Läs in bild
%--------------------------------------

namefile = input('Ange namnet på bilden som ska läsas in: ')

[img,map]=imread(namefile);

if isa(img,'uint8')
    graylevmax=2^8-1;
end
if isa(img,'uint16')
    graylevmax=2^16-1;
end
if isa(img,'uint32')
    graylevmax=2^32-1;
end

%--------------------------------------
% Skala om bilden
%--------------------------------------

% Om bilden är N x M stor med mod(N,8)~=0 eller mod(M,8)~=0
% skalas den inlästa bilden om
imgN=size(img,1);
imgM=size(img,2);
modN=mod(imgN,8);
modM=mod(imgM,8);

img=img(modN+1:imgN,modM+1:imgM);

%--------------------------------------
% Plotta originalbilden
%--------------------------------------
figure
img = double(img)/graylevmax;

imshow(img)
title('Originalbild')

fingerprint = img*graylevmax;

%--------------------------------------
% Beräkna mittpunkten
%--------------------------------------

disp('Beräknar mittpunkt')
disp(' ')

[BinarizedPrint,XofCenter,YofCenter] = centralizing(fingerprint,0);

%--------------------------------------
% Plotta mittpunkten i originalbilden
%--------------------------------------

disp('Tryck enter för att visa mittpunkt')
disp(' ')
pause

hold on;
plot(XofCenter,YofCenter,'o');
hold off;

%--------------------------------------
% Fortsätta eller inte?
%--------------------------------------

ctrl = input('Vill du fortsätta? (J/N) ');
disp(' ')
if ctrl == 'N'
    return;
end

%--------------------------------------
% Klipp ut relevant region av original
% bilden
%--------------------------------------

disp('Klipper ut relevant region av bilden')
disp(' ')

[CroppedPrint]=cropping(XofCenter,YofCenter,fingerprint);
CroppedPrintCpy = CroppedPrint;

%--------------------------------------
% Plotta den klippta bilden
%--------------------------------------

figure
imshow(CroppedPrint/graylevmax)
title('Klippt bild')

%--------------------------------------
% Fortsätta eller inte?
%--------------------------------------

ctrl = input('Vill du fortsätta? (J/N) ');
disp(' ')
if ctrl == 'N'
    return;
end

%--------------------------------------
% Beräkna sektorerna
%--------------------------------------

disp('Beräknar sektorer')
disp(' ')

for (i=1:1:175*175)
    tmp=CroppedPrint(i);
    CroppedPrint(i)=whichsector(i);
    if (CroppedPrint(i)==36 | CroppedPrint(i)==37)
        CroppedPrint(i)=tmp/graylevmax;
    else 
        CroppedPrint(i)=CroppedPrint(i)/64;
    end
    
end

%--------------------------------------
% Plotta sektorerna
%--------------------------------------

figure
imshow(CroppedPrint)
title('Sektorer')

%--------------------------------------
% Fortsätta eller inte?
%--------------------------------------

ctrl = input('Vill du fortsätta? (J/N) ');
disp(' ')
if ctrl == 'N'
    return;
end

%--------------------------------------
% Normalisera bilden
%--------------------------------------

disp('Normaliserar bilden')
disp(' ')

[NormalizedPrint,vector] = sector_norm( CroppedPrintCpy , 0 , 0);
NormalizedPrint = NormalizedPrint/100;


%--------------------------------------
% Plotta den normaliserade bilden
%--------------------------------------

figure
imshow(NormalizedPrint)
title('Normaliserad bild')

%--------------------------------------
% Fortsätta eller inte?
%--------------------------------------

ctrl = input('Vill du fortsätta? (J/N) ');
disp(' ')
if ctrl == 'N'
    return;
end

%--------------------------------------
% Beräkna och plotta Gaborfilterna
%--------------------------------------

disp('Beräknar och plottar Gaborfilter')
disp(' ')

figure
num_disk=8;

filter = struct([]);
for (angle=0:1:num_disk-1)
    filter{angle+1}.gabor=gabor2d_sub(angle,num_disk);
    
    %gabor=gabor*128;
    switch angle<num_disk
        case (angle==0),
            tmp = filter{angle+1}.gabor;
            subplot(241)
            mesh(tmp)
            title('Gabor filter (0^{o})')
        case (angle==1),
            tmp = filter{angle+1}.gabor;
            subplot(242)
            mesh(tmp)
            title('Gabor filter (22.5^{o})')
        case (angle==2),
            tmp = filter{angle+1}.gabor;
            subplot(243)
            mesh(tmp)
            title('Gabor filter (45^{o})')
        case (angle==3),
            tmp = filter{angle+1}.gabor;
            subplot(244)
            mesh(tmp)
            title('Gabor filter (67.5^{o})')
        case (angle==4),
            tmp = filter{angle+1}.gabor;
            subplot(245)
            mesh(tmp)
            title('Gabor filter (90^{o})')
        case (angle==5),
            tmp = filter{angle+1}.gabor;
            subplot(246)
            mesh(tmp)
            title('Gabor filter (112.5^{o})')
        case (angle==6),
            tmp = filter{angle+1}.gabor;
            subplot(247)
            mesh(tmp)
            title('Gabor filter (135^{o})')
        case (angle==7),
            tmp = filter{angle+1}.gabor;
            subplot(248)
            mesh(tmp)
            title('Gabor filter (157.5^{o})')
        otherwise 
            error('Nothing !');
    end
end

%--------------------------------------
% Fortsätta eller inte?
%--------------------------------------

ctrl = input('Vill du fortsätta? (J/N) ');
disp(' ')
if ctrl == 'N'
    return;
end

%--------------------------------------
% Falta bilden med Gaborfilter
%--------------------------------------

disp('Faltar bilden med Gaborfilter')
disp(' ')

figure

feature = struct([]);

for (angle=0:1:num_disk-1)
    
    z2 = filter{angle+1}.gabor;
    z1=NormalizedPrint*100;
    z1x=size(z1,1);
    z1y=size(z1,2);
    z2x=size(z2,1);
    z2y=size(z2,2);    
    ComponentPrint=real(ifft2(fft2(z1,z1x+z2x-1,z1y+z2y-1).*fft2(z2,z1x+z2x-1,z1y+z2y-1)));    
    px=((z2x-1)+mod((z2x-1),2))/2;
    py=((z2y-1)+mod((z2y-1),2))/2;
    ComponentPrint=ComponentPrint(px+1:px+z1x,py+1:py+z1y);
    
    
    [disk,vector]=sector_norm(ComponentPrint,1,0);
    feature{angle+1}.disk = disk;
    feature{angle+1}.vector = vector;
    %img = double(ComponentPrint)/graylevmax;
    img = ComponentPrint/graylevmax;
    
    % FIXA TILARNA PÅ BILDEN
    switch angle<num_disk
        case (angle==0),
            subplot(2,4,1)
            imagesc(img)
			colormap gray
            title('0^{o} komponent')
        case (angle==1),
            subplot(2,4,2)
            imagesc(img)
			colormap gray
            title('22.5^{o} komponent')
        case (angle==2),            
            subplot(2,4,3)
            imagesc(img)
			colormap gray
            title('45^{o} komponent')
        case (angle==3),            
            subplot(2,4,4)
            imagesc(img)
 			colormap gray
           title('67.5^{o} komponent')
        case (angle==4),            
            subplot(245)
            imagesc(img)
			colormap gray
            title('90^{o} komponent')
        case (angle==5),            
            subplot(246)
            imagesc(img)
 			colormap gray
           title('112.5^{o} komponent')
        case (angle==6),            
            subplot(2,4,7)
            imagesc(img)
			colormap gray
            title('135^{o} komponent')
        case (angle==7),            
            subplot(2,4,8)
            imagesc(img)
			colormap gray
            title('157.5^{o} komponent')
        otherwise 
            error('Nothing !');
    end
end

%--------------------------------------
% Fortsätta eller inte?
%--------------------------------------

ctrl = input('Vill du fortsätta? (J/N) ');
disp(' ')
if ctrl == 'N'
    return;
end

%--------------------------------------
% Extrahera egenskaper
%--------------------------------------

figure

disp('Extraherar egenskaper')
disp(' ')

for (angle=0:1:num_disk-1)
    switch angle<num_disk
        case (angle==0),
            img1 = feature{angle+1}.disk/51200;
			subplot(241)
            imagesc(img1)
			colormap bone
            title('Egenskaper (0^{o})')
        case (angle==1),
            img1 = feature{angle+1}.disk/51200;
          subplot(242)
            imagesc(img1)
  			colormap bone
          title('Egenskaper (22.5^{o})')
        case (angle==2),           
            img1 = feature{angle+1}.disk/51200;
            subplot(243)
            imagesc(img1)
   			colormap bone
         title('Egenskaper (45^{o})')
        case (angle==3),            
            img1 = feature{angle+1}.disk/51200;
            subplot(244)
            imagesc(img1)
  			colormap bone
          title('Egenskaper (67.5^{o})')
        case (angle==4),            
            img1 = feature{angle+1}.disk/51200;
            subplot(245)
            imagesc(img1)
   			colormap bone
         title('Egenskaper (90^{o})')
        case (angle==5),            
            img1 = feature{angle+1}.disk/51200;
            subplot(246)
            imagesc(img1)
			colormap bone
            title('Egenskaper (112.5^{o})')
        case (angle==6),            
            img1 = feature{angle+1}.disk/51200;
            subplot(247)
            imagesc(img1)
			colormap bone
            title('Egenskaper (135^{o})')
        case (angle==7),            
            img1 = feature{angle+1}.disk/51200;
            subplot(248)
            imagesc(img1)
			colormap bone
            title('Egenskaper (157.5^{o})')
        otherwise 
            error('Nothing !');
    end
end