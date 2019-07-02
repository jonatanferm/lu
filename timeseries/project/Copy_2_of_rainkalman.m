clear;
load('proj18.mat');
addpath('/Users/ferm/lu/timeseries/tsadl/matlab')
%%
%y = Kassala.rain_org;
%t_y = Kassala.rain_org_t;
%t_x = Kassala.rain_t;

y = ElGeneina.rain_org;
t_y = ElGeneina.rain_org_t;
t_x = ElGeneina.rain_t;
my = mean(y);
y_org = y;
[lambda, ~, offset] = bcNormPlot(y, false);
y = (lambda^-1)*((y+offset).^lambda-1);
ynz = y(y~=0);
model = arx(y, 1);
res = resid(model, y);
%%

%%
xoffs = 5;

a1 = model.A(2);
ap = [-1 -1 -a1 -a1 -a1];
A = [ap; eye(4), zeros(4, 1)];

V_wt =  0.1;
V_et = var(res.OutputData);
R_w = V_wt;

C = [1 1 1 0 0];
x_hat_1 = zeros(5, 1);
R_xx_1 = eye(5) * var(y);
%R_e = eye(5) * V_et;
R_e = (diag(ones(2, 1), 3)+diag(ones(2, 1), -3))*a1*V_wt + eye(5)*(V_wt+V_et);

R_yy_1 = var(y);

x = nan(numel(t_x) + xoffs, 1);
x(1) = 0;
x(end) = 0;
w = zeros(size(t_y));
for i = 1:(numel(x)-xoffs)
    R_yy_1 = C*R_xx_1*C'+R_w;
    K_t = (R_xx_1*C')/R_yy_1;
    x_hat = x_hat_1;
    %w_1(i+xoffs) = ypm(i) - C*x_hat;
    yi = find(abs(t_y - t_x(i)) < 0.001, 1, 'first');
    if(yi)
         x_hat = x_hat + K_t*(y(yi) - C * x_hat);
         w(yi) = y(yi) - C * x_hat;
         x(i-4+xoffs:i+xoffs) = x_hat;
         %x_hat_1 = x_hat_1 + B*w(i-1+xoffs:i+1+xoffs);
    end
    R_xx = (eye(5) - K_t*C)*R_xx_1;
    R_xx_1 = A*R_xx*A'+R_e;
    x_hat_1 = A*x_hat;
end
y_kl_recon = [zeros(2, 1); sum([x(:), [x(2:end); 0], [x(3:end); 0; 0]], 2)];
y_kl_recon = [y_kl_recon(6:end-3); 0];
x = x(xoffs+1:end);
y_linear_recon = interp1(t_y, y, t_x);
y_linear_recon(1) = 0;
y_linear_recon(end) = 0;

colors = [197,27,125   
77,146,33]/255;
figure('Name', 'Linear interpolation and Kalman reconsctruction')
hold on;
set(gca,'Color',[247,247,247]/255)
set(gca, 'ColorOrder', colors);
plot(t_x, y_kl_recon, '-','LineWidth', 1.5);
plot(t_x, y_linear_recon, '--', 'LineWidth', 1)
legend({'Kalman reconstruction', 'Linear Interpolation'},'Location','northwest')
hold off;

rang = 4:3:numel(t_x);
var(w)
sum(y)./sum(x)
%%
Elgeneina_rain_recon = y_kl_recon;

%Kassala_rain_recon = max(x, 0);
save('Elgeneina_rain_recon_nt.mat', 'Elgeneina_rain_recon')
%save('Kassala_rain_recon.mat', 'Kassala_rain_recon')