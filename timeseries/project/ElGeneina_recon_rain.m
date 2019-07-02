clear;
load('proj18.mat');
load('Elgeneina_rain_recon.mat')
addpath('/Users/ferm/lu/timeseries/tsadl/matlab')
%%
u = Elgeneina_rain_recon;
u_t = ElGeneina.rain_t;
y = ElGeneina.nvdi;
y_t = ElGeneina.nvdi_t;
yv_t = y_t;
%%
start_index = find(u_t == y_t(1));
u = u(start_index:end);
u_t = u_t(start_index:end);
assert(numel(find(abs(u_t - y_t) > 0.00001)) == 0)
%%
y_validation = y;
u_validation = u;
nbr_of_validation_points = 36*2;
y = y(1:end-nbr_of_validation_points);
y_t = y_t(1:end-nbr_of_validation_points);
u = u(1:end-nbr_of_validation_points);
u_t = u_t(1:end-nbr_of_validation_points);
%%
y_offset = - mean(y);
u_offset = - mean(u);
y = y + y_offset;%log(yraw);
u = u + u_offset;
[lambda, ~, offset] = bcNormPlot(y, false);
%u = (lambda^-1)*((uraw+offset).^lambda-1);
y = (lambda^-1)*(((y+offset).^lambda)-1);
%y = log(y-min(y)+1);
%y = y-mean(y);
%y = sqrt(y);
%scale_factor = std(y)/std(u);
%u = u*scale_factor;
%save('data_transform.mat', 'scale_factor', 'y_offset', 'u_offset')
figure;
hold on;
plot(u)
plot(y)
hold off
%%
u_s = filter([1 zeros(1, 35) -1], 1, u);
%%
analyzets(u_s);
%%
Ar = [1 zeros(1, 40)];
Cr = [1 zeros(1, 40)];
Ar_free = Ar;
Cr_free = Cr;
Ar_free([2 13 37 38]) = 1;
Cr_free([2 3 4 5]) = 1;
model_poly = idpoly(Ar, [], Cr);
model_poly.Structure.A.Free = Ar_free;
model_poly.Structure.C.Free = Cr_free;
model_rain = pem(u_s, model_poly);
res = resid(model_rain, u_s);
analyzets(res);
present(model_rain)
%%
r_acf = acf(res.Outputdata(1:end), floor(numel(res.OutputData)/4));
dagostinoK2test(r_acf(2:end), 0.05)
%% predicting rain
A = conv(model_rain.a, [1 zeros(1, 35) -1]);
[CS, AS] = equalLength(model_rain.c, A);
k = 1;
[Fk, Gk] = deconv(conv([1,zeros(1, k-1)], CS), AS);
u_k = filter(Gk, model_rain.c, u);
%u_k = u_k(length(Gk)+1:end);
figure;
hold on;
plot(u)
plot([300, 300+k], [1.2*max(u) 1.2*max(u)], 'k', 'LineWidth', 1.5)
plot(u_k, '--')
hold off
%%
ypw = filter(A, model_rain.c, y);
upw = filter(A, model_rain.c, u);
M = 80;
stem(-M:M, crosscorr(upw, ypw, M));
hold on;
refline([0 2/sqrt(size(u, 1))])
refline([0 -2/sqrt(size(u, 1))])
hold off
%%
%%
A2 = [1 zeros(1, 10)];
A2_free = A2;
A2_free([1 2 3]) = 1;
d = 7;
dr = 0;
B = [0 zeros(1, 10)];
B_free = B;
B_free(d:d+dr) = 1;
Mi = idpoly(1, B, [], [], A2);
Mi.Structure.F.Free = A2_free;
Mi.Structure.B.Free = B_free;
zpw = iddata(ypw, upw);
Mba2 = pem(zpw, Mi);
vhat = resid(Mba2, zpw);
analyzets(vhat)
present(Mba2)
%%
df = 4;
dfr = 0;
B_f = zeros(1, 10);
B_ff = B_f;
B_ff(df:df+dfr) = 1;

rc_f = [1 zeros(1, 40)];
rc_ff = zeros(size(rc_f));
rc_ff(find(model_rain.c~=0)) = 1;
rc_ff([2 5]) = 0;
%rc_ff([2 3]) = 1;
%rc_ff([2]) = 0;

ra_f = [1 zeros(1, 80)];
ra_ff = zeros(size(ra_f));
%ra_ff([2 3 35 36 37 38]) = 1;
ra_ff(find(A~=0)) = 1;
ra_ff([38 49 74]) = 0;

A2_f = [1 zeros(1, 10)];
A2_ff = zeros(size(A2_f));
A2_ff([1 3]) = 1;

Mi = idpoly(1, B_f, rc_f, ra_f, A2_f);
Mi.Structure.F.Free = A2_ff;
Mi.Structure.C.Free = rc_ff;
Mi.Structure.D.Free = ra_ff;
Mi.Structure.B.Free = B_ff;
z = iddata(y, u);
MboxJ = pem(z, Mi);
present(MboxJ);
ehat = resid(MboxJ, z);
analyzets(ehat)
%%
whitenessTest(ehat.Outputdata)
ehacf = acf(ehat.OutputData, 100);
dagostinoK2test(ehacf(2:end), 0.05)
dagostinoK2test(ehat.OutputData, 0.05)
%%
ElGeneina_model = MboxJ;
save('ElGeneina_model.mat', 'ElGeneina_model')
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