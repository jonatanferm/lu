clear;
load('proj18.mat');
addpath('/Users/ferm/lu/timeseries/tsadl/matlab')
%%
y = ElGeneina.rain_org;
t_y = ElGeneina.rain_org_t;
t_x = ElGeneina.rain_t;
[lambda, ~, offset] = bcNormPlot(y, false);
y = (lambda^-1)*((y+offset).^lambda-1);
%y = y-mean(y);
%%
Ar = [1 zeros(1, 10)];
Cr = [1 zeros(1, 10)];
Ar_free = Ar;
Cr_free = Cr;
Ar_free([2]) = 1;
Cr_free([]) = 1;
model_poly = idpoly(Ar, [], Cr);
model_poly.Structure.A.Free = Ar_free;
model_poly.Structure.C.Free = Cr_free;
model_rain = pem(y, model_poly);
res = resid(model_rain, y);
%analyzets(res);
present(model_rain)

[AS, CS] = equalLength(model_rain.a, model_rain.c);
xa = conv(model_rain.a, [1 1 1]);
[Fyk, Gyk] = deconv(xa, CS);

A = [-xa(2:6)' [eye(4); zeros(1, 4)]];
V_wt =  1.5037e-32;
V_et = var(res.OutputData);
r_0 = mar(model_rain.c, 0, V_et) + mar(model_rain.a, 0, V_wt);
r_1 = mar(model_rain.c, 1, V_et) + mar(model_rain.a, 1, V_wt);
r_2 = mar(model_rain.c, 2, V_et) + mar(model_rain.a, 2, V_wt);
r_3 = mar(model_rain.c, 3, V_et) + mar(model_rain.a, 3, V_wt);
r_4 = mar(model_rain.c, 4, V_et) + mar(model_rain.a, 4, V_wt);
R_e = diag(repmat(r_0, 1, 5), 0) + ...
    diag(repmat(r_1, 1, 4), 1) + diag(repmat(r_1, 1, 4), -1) + ...
    diag(repmat(r_2, 1, 3), 2) + diag(repmat(r_2, 1, 3), -2) + ...
    diag(repmat(r_3, 1, 2), 3) + diag(repmat(r_3, 1, 2), -3) + ...
    diag(repmat(r_4, 1, 1), 4) + diag(repmat(r_4, 1, 1), -4);
R_w = V_wt; 

C = [1 1 1 0 0];
x_hat_1 = zeros(5, 1);
R_xx_1 = eye(5) * var(y) * 0.001;
R_yy_1 = var(y) * 0.001;
padding = 5;
x = zeros(numel(t_x) + padding, 1);

ypm = interp1(t_y, y, t_x);
ypm(1) = 0;
ypm(end) = 0;
residual = zeros(size(y));
for i = 1:numel(t_x)
    K_t = (R_xx_1*C')/R_yy_1;
    %ypm(i) = mean(yp(i:i+2));
    x_hat = x_hat_1;
    
    x_hat = x_hat + K_t*(ypm(i) - C * x_hat_1);
    %x(i) = x_hat(1);
    
    yi = find(abs(t_y - t_x(i)) < 0.001, 1, 'first');
    if(yi)
        residual(yi) = ypm(i) - sum(x_hat(1:3));
    %    x_hat = x_hat + K_t*(y(yi) - C * x_hat_1);
        x(i-2+padding:i+padding) = x_hat(1:3);
    end
    
    x_hat_1 = A*x_hat;
    R_xx = (eye(5) - K_t*C)*R_xx_1;
    R_xx_1 = A*R_xx*A'+R_e;
    R_yy_1 = C*R_xx_1*C'+R_w;
end
x = x(padding+1:end);
ypn = interp1(t_y, y, t_x, 'next');
xy = x.*(ypn>0);
xy(xy<0) = 0;
nf = sum(y)/sum(xy)
%xy = xy.*nf;
figure('Name', 'y/3 and normrecon')
hold on;
plot(t_x, ypn./3, 'r');
plot(t_x, xy, 'k')
hold off;

figure('Name', 'ypm/3 and recon')
hold on;
plot(t_x, ypm./3, 'r');
plot(t_x, x, 'k')
hold off;

xys = xy(3:end-1);
xys = permute(sum(reshape(xys(3:end-1), 3, [])), [2 1]);
figure('Name', 'Summed normrecond and org')
hold on;
plot(t_y, y, 'r');
plot(t_y(2:end-1), xys, 'k')
hold off;
var(residual)
%%
mar(C, 7, 1)
%%
function rk = mar(C, k, sigma2)
    ak = abs(k) + 1;
    mc = [1 C];
    if (ak > length(mc))
        rk = 0;
    elseif (ak == length(mc))
        rk = sigma2*mc(ak);
    else
        rk = sigma2*sum(mc(ak:end).*mc(1:end-ak+1));
    end
end