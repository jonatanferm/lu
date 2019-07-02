clear;
load('proj18.mat');
addpath('/Users/ferm/lu/timeseries/tsadl/matlab')
%%
y = ElGeneina.rain_org;
t_y = ElGeneina.rain_org_t;
t_x = ElGeneina.rain_t;
my = mean(y);
y_org = y;
[lambda, ~, offset] = bcNormPlot(y, false);
y = (lambda^-1)*((y+offset).^lambda-1);
%%

%a1 = model.A(2);
%ap = [-(a1+1) -(a1+1) -a1];

V_wt = 1;
V_et = 10;


sum_quote = sumdiff(ind)
x_Variance = xvar(ind)
a1 = a1space(ind)
y_kl_recon = [zeros(2, 1); sum([x(:), [x(2:end); 0], [x(3:end); 0; 0]], 2)];
y_kl_recon = [y_kl_recon(xoffs+1:end-3); 0];
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
%%
Elgeneina_rain_recon = y_kl_recon;

%Kassala_rain_recon = max(x, 0);
save('Elgeneina_rain_recon_nt.mat', 'Elgeneina_rain_recon')
%save('Kassala_rain_recon.mat', 'Kassala_rain_recon')
%%
f = get_norm_fun(y, t_x, t_y);
options = optimset('Display','iter','PlotFcns',@optimplotfval);
fminsearch(f,[0.1, 10], options)
%%
function f = get_norm_fun(y, t_x, t_y)
function n = residual_norm(arg)
    V_wt = abs(arg(1));
    V_et = abs(arg(2));
    xoffs = 2;
    R_w = V_wt;
    C = [1 1 1];
    R_xx_1 = eye(3) * var(y);
    R_e = eye(3) * V_et;
    a1space = linspace(-.99, -.01);
    wnorm = zeros(size(a1space));
    c = 1;
    for a1=a1space
        ap = [0 0 -a1];
        A = [ap; eye(2), zeros(2, 1)];
        x_hat_1 = zeros(3, 1);
        x = nan(numel(t_x) + xoffs, 1);
        x(1) = 0;
        x(end) = 0;
        w = zeros(size(t_y));
        for i = 1:(numel(x)-xoffs)
            R_yy_1 = C*R_xx_1*C'+R_w;
            K_t = (R_xx_1*C')/R_yy_1;
            x_hat = x_hat_1;
            yi = find(abs(t_y - t_x(i)) < 0.001, 1, 'first');
            if(yi)
                w(yi) = y(yi) - C * x_hat;
                x_hat = x_hat + K_t*(y(yi) - C * x_hat);
                x(i-2+xoffs:i+xoffs) = x_hat;
            end
            R_xx = (eye(3) - K_t*C)*R_xx_1;
            R_xx_1 = A*R_xx*A'+R_e;
            x_hat_1 = A*x_hat;
        end
        
        wnorm(c) = norm(w);
        c = c+1;
    end
    [n, ~] = min(wnorm);
end
    f = @residual_norm;
end