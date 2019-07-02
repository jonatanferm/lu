clear;
load('proj18.mat');
addpath('/Users/ferm/lu/timeseries/tsadl/matlab')
%%
y = ElGeneina.nvdi;
y_t = ElGeneina.nvdi_t;
y_validation = y;
yv_t = y_t;
nbr_of_validation_points = 36*2;
y = y(1:end-nbr_of_validation_points);
y_t = y_t(1:end-nbr_of_validation_points);
%%
figure; analyzets(y);
%%
A = [1 zeros(1, 80)];
C = [1 zeros(1, 80)];
Af = [0 zeros(1, 80)];
Cf = [0 zeros(1, 80)];
Af([2 73]) = 1;
Cf([3]) = 1;

model_poly = idpoly(A, [], C);
model_poly.Structure.A.Free = Af;
model_poly.Structure.C.Free = Cf;
model = pem(y, model_poly);
present(model)
res = resid(model, y);

figure; analyzets(res, 40);

r_acf = acf(res.Outputdata(1:end), floor(numel(res.OutputData)/4));
dagostinoK2test(r_acf(2:end), 0.05)
%%
figure;
hold on;
acf(res.OutputData, 40, 0.05, true);
ylim([-.12 .12])
hold off;
%%
[CS, AS] = equalLength(model.c, model.a);

%k_range = 1:1;
k_range = [1 4 6 8 10];
yhat_k = zeros(size(y_validation, 1), numel(k_range));
for ik = 1:numel(k_range)
    k = k_range(ik);
    [Fk, Gk] = deconv(conv([1,zeros(1, k-1)], CS), AS);
    yhat = filter(Gk, model.c, y_validation);
    yhat_k(:,ik) = yhat;
    
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
xlim([y_t(end-nbr_of_validation_points) y_t(end)])
set(gca, 'ColorOrder', colors);
plot(yv_t, yhat_k, 'LineWidth', 1.5)
plot(yv_t, y_validation, 'k', 'LineWidth', 2)
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
%%

