clear;
diffs = [];

load('pro_1.1.mat');
load('311.mat');
[ad, aZ, at] = procrustes(ax', ay', 'scaling', false);
[md, mZ, mt] = procrustes(xps, yps, 'scaling', false);
aM = affine2d([[at.T [0; 0]]; [at.c(1,:) 1]]);
ayp = transformPointsForward(aM, yps);
mM = affine2d([[mt.T [0; 0]]; [mt.c(1,:) 1]]);
myp = transformPointsForward(mM, yps);
diffs = horzcat(diffs, sqrt(mean(sum((ayp - myp) .^ 2, 2))));

load('pro_4.5.mat');
load('345.mat');
[ad, aZ, at] = procrustes(ax', ay', 'scaling', false);
[md, mZ, mt] = procrustes(xps, yps, 'scaling', false);
aM = affine2d([[at.T [0; 0]]; [at.c(1,:) 1]]);
ayp = transformPointsForward(aM, yps);
mM = affine2d([[mt.T [0; 0]]; [mt.c(1,:) 1]]);
myp = transformPointsForward(mM, yps);
diffs = horzcat(diffs, sqrt(mean(sum((ayp - myp) .^ 2, 2))));

load('pro_6.0.mat');
load('360.mat');
[ad, aZ, at] = procrustes(ax', ay', 'scaling', false);
[md, mZ, mt] = procrustes(xps, yps, 'scaling', false);
aM = affine2d([[at.T [0; 0]]; [at.c(1,:) 1]]);
ayp = transformPointsForward(aM, yps);
mM = affine2d([[mt.T [0; 0]]; [mt.c(1,:) 1]]);
myp = transformPointsForward(mM, yps);
diffs = horzcat(diffs, sqrt(mean(sum((ayp - myp) .^ 2, 2))));

load('pro_6.3.mat');
load('363.mat');
[ad, aZ, at] = procrustes(ax', ay', 'scaling', false);
[md, mZ, mt] = procrustes(xps, yps, 'scaling', false);
aM = affine2d([[at.T [0; 0]]; [at.c(1,:) 1]]);
ayp = transformPointsForward(aM, yps);
mM = affine2d([[mt.T [0; 0]]; [mt.c(1,:) 1]]);
myp = transformPointsForward(mM, yps);
diffs = horzcat(diffs, sqrt(mean(sum((ayp - myp) .^ 2, 2))));

load('pro_8.3.mat');
load('383.mat');
[ad, aZ, at] = procrustes(ax', ay', 'scaling', false);
[md, mZ, mt] = procrustes(xps, yps, 'scaling', false);
aM = affine2d([[at.T [0; 0]]; [at.c(1,:) 1]]);
ayp = transformPointsForward(aM, yps);
mM = affine2d([[mt.T [0; 0]]; [mt.c(1,:) 1]]);
myp = transformPointsForward(mM, yps);
diffs = horzcat(diffs, sqrt(mean(sum((ayp - myp) .^ 2, 2))));

load('pro2_3.mat');
load('43.mat');
[ad, aZ, at] = procrustes(ax', ay', 'scaling', true);
[md, mZ, mt] = procrustes(xps, yps, 'scaling', true);
aM = affine2d([[at.b * at.T [0; 0]]; [at.c(1,:) 1]]);
ayp = transformPointsForward(aM, yps);
mM = affine2d([[at.b * mt.T [0; 0]]; [mt.c(1,:) 1]]);
myp = transformPointsForward(mM, yps);
diffs = horzcat(diffs, sqrt(mean(sum((ayp - myp) .^ 2, 2))));

load('pro2_7.mat');
load('47.mat');
[ad, aZ, at] = procrustes(ax', ay', 'scaling', true);
[md, mZ, mt] = procrustes(xps, yps, 'scaling', true);
aM = affine2d([[at.b * at.T [0; 0]]; [at.c(1,:) 1]]);
ayp = transformPointsForward(aM, yps);
mM = affine2d([[at.b * mt.T [0; 0]]; [mt.c(1,:) 1]]);
myp = transformPointsForward(mM, yps);
diffs = horzcat(diffs, sqrt(mean(sum((ayp - myp) .^ 2, 2))));

load('pro2_9.mat');
load('49.mat');
[ad, aZ, at] = procrustes(ax', ay', 'scaling', true);
[md, mZ, mt] = procrustes(xps, yps, 'scaling', true);
aM = affine2d([[at.b * at.T [0; 0]]; [at.c(1,:) 1]]);
ayp = transformPointsForward(aM, yps);
mM = affine2d([[at.b * mt.T [0; 0]]; [mt.c(1,:) 1]]);
myp = transformPointsForward(mM, yps);
diffs = horzcat(diffs, sqrt(mean(sum((ayp - myp) .^ 2, 2))));

load('pro2_10.mat');
load('410.mat');
[ad, aZ, at] = procrustes(ax', ay', 'scaling', true);
[md, mZ, mt] = procrustes(xps, yps, 'scaling', true);
aM = affine2d([[at.b * at.T [0; 0]]; [at.c(1,:) 1]]);
ayp = transformPointsForward(aM, yps);
mM = affine2d([[at.b * mt.T [0; 0]]; [mt.c(1,:) 1]]);
myp = transformPointsForward(mM, yps);
diffs = horzcat(diffs, sqrt(mean(sum((ayp - myp) .^ 2, 2))));

load('pro2_12.mat');
load('412.mat');
[ad, aZ, at] = procrustes(ax', ay', 'scaling', true);
[md, mZ, mt] = procrustes(xps, yps, 'scaling', true);
aM = affine2d([[at.b * at.T [0; 0]]; [at.c(1,:) 1]]);
ayp = transformPointsForward(aM, yps);
mM = affine2d([[at.b * mt.T [0; 0]]; [mt.c(1,:) 1]]);
myp = transformPointsForward(mM, yps);
diffs = horzcat(diffs, sqrt(mean(sum((ayp - myp) .^ 2, 2))));