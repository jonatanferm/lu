clear;
load('proj18.mat');
load('Kassala_rain_recon.mat')
addpath('/Users/ferm/lu/timeseries/tsadl/matlab')
load('Elgeneina_model.mat');
load('data_transform.mat');
%%
yraw = Kassala.nvdi;
y_t = Kassala.nvdi_t;
uraw = Kassala_rain_recon;
u_t = Kassala.rain_t;
start_index = find(u_t == y_t(1));
uraw = uraw(start_index:end);
u_t = u_t(start_index:end);
assert(numel(find(abs(u_t - y_t) > 0.00001)) == 0)
%%
y = yraw;
u = uraw;
%%
y_offset = - mean(y);
u_offset = - mean(u);
y = y + y_offset;
u = u + u_offset;
[lambda, ~, offset] = bcNormPlot(y, false);
%u = (lambda^-1)*((uraw+offset).^lambda-1);
y = (lambda^-1)*((y+offset).^lambda-1);
y = y - mean(y);
figure;
hold on;
plot(u)
plot(y)
hold off
%%
y_validation = y;
u_validation = u;
nbr_of_validation_points = 0;
y = y(1:end-nbr_of_validation_points);
u = u(1:end-nbr_of_validation_points);
%%
z = iddata(y, u);
ehat = resid(ElGeneina_model, z);
analyzets(ehat)
present(ElGeneina_model)
MboxJ = ElGeneina_model;
%%
retrained = pem(z, ElGeneina_model);
ehatr = resid(ElGeneina_model, z);
analyzets(ehatr)
present(retrained)
MboxJ = retrained;
%MboxJ = ElGeneina_model;
%%
A =  MboxJ.D;
C =  MboxJ.C;

A2 =  MboxJ.F;
B =  MboxJ.B;

k_range = [1 4 6 8 10];
yhat_k = zeros(size(y_validation, 1), numel(k_range));
for ik = 1:numel(k_range)
    k = k_range(ik);
    [AS, CS] = equalLength(A, C);
    [Fyk, Gyk] = deconv(conv([1,zeros(1, k-1)], CS), AS);

    [BFS, CS] = equalLength(conv(B, Fyk), C);
    [Fuk, Guk] = deconv(conv([1,zeros(1, k-1)], BFS), CS);

    yhat_k(:,ik) = filter(Gyk, C, y_validation) + filter(Guk, C, u_validation);
end
colors = flipud([
    254,237,222
253,190,133
253,141,60
230,85,13
166,54,3]./255);
figure;
hold on;
set(gca,'Color',[107,174,214]/255)
%ylim([-3 5])
xlim([y_t(end-72) y_t(end)])
set(gca, 'ColorOrder', colors);
plot(y_t, yhat_k, 'LineWidth', 1.5)
plot(y_t, y_validation, 'k', 'LineWidth', 2)
legend({'k=1','k=4','k=6','k=8','k=10','observed'},'Location','northwest')
legend('boxoff')
hold off;
%%
alim = 20;
yhat_acf = zeros(alim+1, size(yhat_k, 2));
for i=1:size(yhat_k, 2)
    yhat_acf(:,i) = acf(y_validation-yhat_k(:,i), alim);
end
figure;
hold on;
set(gca,'Color',[107,174,214]/255)
%ylim([-3 5])
xlim([0 alim])
set(gca, 'ColorOrder', colors);
plot(0:alim, yhat_acf, 'LineWidth', 1.5)
legend({'k=1','k=4','k=6','k=8','k=10'},'Location','northeast')
legend('boxoff')
hold off;