clear;
load('proj18.mat');
addpath('/Users/ferm/lu/timeseries/tsadl/matlab')
%%
y = ElGeneina.rain_org;
t_y = ElGeneina.rain_org_t;
t_x = ElGeneina.rain_t;
[lambda, ~, offset] = bcNormPlot(y, false);
y = (lambda^-1)*((y+offset).^lambda-1);
model = arx(y, 1);
res = resid(model, y);
%%
a1 = model.A(2);
%%
ap = [-(a1+1) -(a1+1) -a1];
xmodel = idpoly([1 ap]);
covX = covM(sim(xmodel, randn(10000, 1)), 4, true);
covX = (covX./covX(1, 1)) * var(y);

covX = [covX(1,2:end); covX(3,1) covX(2,3:end);covX(4,1:end-1);];
%%
A = ([ap' [1 0; 0 1; 0 0]]')^3;
V_wt=  10.002395;
V_et = var(res.OutputData) * 10;
r_1 = V_et+V_wt*(1 + a1^2);
r_2 = V_wt*(1+a1);
r_3 = V_wt*(a1);
R_e = diag(repmat(r_1, 1, 3), 0) + ...
    diag(repmat(r_2, 1, 2), 1) + diag(repmat(r_2, 1, 2), -1) + ...
    diag(repmat(r_3, 1, 1), 2) + diag(repmat(r_3, 1, 1), -2);
R_e = R_e^3;
R_w = V_wt; 
B = 0;
C = [1 1 1];
%x_hat_1 = repmat(mean(y), [3 1]);
x_hat_1 = zeros(3, 1);
R_xx_1 = eye(3) * var(y);
R_yy_1 = var(y);

%R_yy_1 = var(y);
x = nan(numel(t_x), 1);
x(1) = 0;
x(end) = 0;

ypm = interp1(t_y, y, t_x);
ypm(1) = 0;
ypm(end) = 0;

R_xx_3 = covX;
R_yy_1_inv = 1/var(y);
K_t = (R_xx_3*C')/R_yy_1_inv;
%%
for i = 1:numel(y)
    x_hat = x_hat_1 + K_t*(y(i) - C * x_hat_1);
    x(i*3-2:i*3) = x_hat;
    x_hat_1 = A*x_hat;
end
ypn = interp1(t_y, y, t_x, 'next');
ypn(1) = 0;
ypn(end) = 0;
xy = x.*(ypn>0);
xy(xy<0) = 0;
xy = (xy./(sum(xy))).*sum(y);
figure('Name', 'y/3 and normrecon')
hold on;
plot(t_x, ypn./3, 'r');
plot(t_x, xy, 'k')
hold off;

figure('Name', 'ypm and recon')
hold on;
plot(t_x, ypm./3, 'r');
plot(t_x, x, 'k')
hold off;

xys = xy(1:end-3);
xys = permute(sum(reshape(xys, 3, [])), [2 1]);
figure('Name', 'Summed normrecond and org')
hold on;
plot(t_y, y, 'r');
plot(t_y(1:end-1), xys, 'k')
hold off;
var(y(1:end-1)-xys)

Elgeneina_rain_recon = xy;
save('Elgeneina_rain_recon.mat', 'Elgeneina_rain_recon')