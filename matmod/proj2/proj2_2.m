clear;
clc;
%%
data = xlsread('pk',1);
data = data(:,2:end);
data = permute(reshape(data', [3, 10, 10]), [2 1 3]);
%%
set(groot,'defaultAxesColorOrder',winter(10));
plot(squeeze(data(:,1,:)), squeeze(data(:,2,:)), 'x-','LineWidth', 1.5);
%%

function C = Ct(t, CL, alpha, f)
    A = [1-alpha, 0, 0; alpha, 1-CL-f, f; 0, f, 1-f];
    
end