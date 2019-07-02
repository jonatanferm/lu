function [disk,vector] = sector_norm( image , mode , mix)
% Modified by Luigi Rosa

% N=175 size of cropped image (175 x 175)
N=175;
% Number of sectors
M=38;

size_m=N*N;


mean_s=zeros(M,1);
varn_s=zeros(M,1);
num_s=zeros(M,1);
image1=zeros(175,175);
Mo=50;
Vo=50;

for ( i=1:1:size_m)
    tmp=whichsector(i);
    tmp=tmp+1;
    if (tmp>=1)
        mean_s(tmp)= mean_s(tmp)+image(i);
        num_s(tmp)=num_s(tmp)+1;
    end   
end   

for (i=1:1:M)
    mean_s(i)=mean_s(i)/num_s(i);
end   

for ( i=1:1:size_m)
    tmp=whichsector(i);
    tmp=tmp+1;
    if (tmp>=1)
        varn_s(tmp)= varn_s(tmp) + (image(i)- mean_s(tmp))^2;
    end   
end   

for (i=1:1:M)
    varn_s(i)= varn_s(i) / num_s(i);
end   

if (mix==0 | mix==1)
    for (i=1:1:size_m)
        tmp=whichsector(i);
        tmp=tmp+1;			
        image1(i)=varn_s(tmp);
    end
end

if ( mode == 0 )
    for ( i=1:1:size_m)
        tmp=whichsector(i);
        tmp=tmp+1;
        if (tmp>=1 & abs(varn_s(tmp))>1)
            if ((image(i) - mean_s(tmp))<0)
                if (tmp==37 | tmp==38) & mix==0
                    image1(i)=50;
                else
                    image1(i)=Mo - (Vo/varn_s(tmp)*((image(i) - mean_s(tmp))^2))^0.5;
                end
            else 
                if (tmp==37 | tmp==38) & mix==0
                    image1(i)=50;
                else
                    image1(i)=Mo + (Vo/varn_s(tmp)*((image(i) - mean_s(tmp))^2))^0.5;
                end
            end
        else
            image1(i)=Mo;
        end
    end
    
    disk=image1;
    vector=varn_s;   
else   
    disk=image1;
    vector=varn_s;    
end