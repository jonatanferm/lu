clear;
clc;
%%
A = -2;
AW = 5;
WB = 5;
B = 3;
n = [1 1];

N = 300;
wp = linspace(A, B, N)';
paths = [repmat(A, N, 1) wp repmat(B, N, 1)];
path_lengths = sum((diff(paths, 1, 2).^2 + [AW^2 WB^2]).*n, 2);

%%Coloring
colors = flipud([255,245,240
    254,224,210
    252,187,161
    252,146,114
    251,106,74
    239,59,44
    203,24,29
    165,15,21
    103,0,13]/255);
best_line_color = [66,146,198]/255;
num_colors = size(colors, 1);

[min_length, min_index] = min(path_lengths);
path_diff = path_lengths - min_length;
path_colors = round((num_colors-1)*(path_diff/max(path_diff)))+1;

figure;
hold on;
set(gca,'ColorOrder',colors(path_colors,:));
X = [0, AW, AW+WB];
plot(X, paths', 'LineWidth', 2);
plot(X, paths(min_index,:)', 'LineWidth', 3, 'Color', best_line_color);
hold off;
%%
D = 1;
d = 5;
n = 1.1;
h = linspace(-D, D);
t = (sqrt(D^2 + d^2)*n - ...
    sqrt(d^2 * n^2 + h.^2 * n^2 - 2*sqrt(D^2 + d^2)*d*n + D^2 + d^2 - h.^2) - d)/(n^2 - 1);
plot([-t' t'], h', 'k')
grid on
daspect([1 1 1])