function [F, STR] = get_features(image, mask)
%keyboard;
F(1) = mean(image(find(mask)));
STR{1} = 'mean intensity';

F(2) = std(image(find(mask)));
STR{2} = 'std dev';

F(3) = mean(mean(mask)) + mean(mean(mask, 2));
STR{3} = 'mask size';

F(4) = max(image(find(mask)));
STR{4} = 'max val';

%pow = sum(real(fft(xcorr2(image))), 2);
%range = abs(pow - pow(1)/10);
%[~, range] = min(range);

cov = [mean(xcorr2(image)) 0];
lc = numel(cov)/2;
cov = sum([fliplr(cov(1:lc)); cov(lc+1:end)]);
range = abs(cov - cov(1)/10);
[~, range] = min(range);

F(4) = range;
STR{4} = 'range';

%F(5) = cov(2);
%STR{5} = 'cov2';
%F(6) = cov(2);
%STR{6} = 'cov2';
%F(7) = cov(3);
%STR{7} = 'cov3';
%F(8) =cov(4);
%STR{8} = 'cov4';
%{
F(9) = cov(5);
STR{9} = 'cov5';
F(10) = cov(6);
STR{10} = 'cov6';
F(11) =cov(7);
STR{11} = 'cov7';
F(12) = cov(8);
STR{12} = 'cov8';
%}
% Need to name all features.
assert(numel(F) == numel(STR));