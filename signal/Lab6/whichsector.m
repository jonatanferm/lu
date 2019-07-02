function [sector_num] = whichsector(index)
% Modiofied by Luigi Rosa
% index is the index of current pixel of cropped image ( cropped image is
% 175 x 175 ); sector_num is the output and represents what is the
% corresponding sector.

length = 175;
x = rem( index , length );
y = floor(index / length);

x = x - floor(length / 2);
y = y - floor(length / 2);

rad = (x*x) + (y*y);
if rad < 144        % innerest radius = 12 (144=12*12)
    sector_num = 36;
    sector_num;
    return
end

if rad >= 5184       % outtest radius = 72 (5184=72*72)
    sector_num = 37;
    sector_num;
    return
end   

if x ~= 0
    theta = atan( y / x );
else 
    if y > 0
        theta = pi/2;
    else
        theta = -pi/2;
    end
end   

if x < 0
    theta = theta + pi;
else
    if theta < 0
        theta = theta + 2*pi;
    end
end

if theta < 0
    theta = theta + 2*pi;
end

r = floor(rad ^ 0.5);
ring = floor(( r-12 )/20);
arc = floor(theta /(pi/6));

sector_num = ring * 12 + arc;