clear;
load('proj18.mat');
load('Elgeneina_rain_recon_nt.mat')
addpath('/Users/ferm/lu/timeseries/tsadl/matlab')
%%
yraw = ElGeneina.nvdi;
y_t = ElGeneina.nvdi_t;
uraw = Elgeneina_rain_recon;
u_t = ElGeneina.rain_t;
%%
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
%scale_factor = std(y)/std(u);
%u = u*scale_factor;
%save('data_transform.mat', 'scale_factor', 'y_offset', 'u_offset')
y_validation = y;
u_validation = u;
nbr_of_validation_points = 36*2;
y = y(1:end-nbr_of_validation_points);
u = u(1:end-nbr_of_validation_points);
figure;
hold on;
plot(u)
plot(y)
hold off
%%
u_s = filter([1 zeros(1, 35) -1],1, u);
%%
analyzets(u_s);
%%
Ar = [1 zeros(1, 40)];
Cr = [1 zeros(1, 20)];
Ar_free = Ar;
Cr_free = Cr;
Ar_free([2 3 4 13 37 38 39]) = 1;
Cr_free([2 3 4 5]) = 1;
model_poly = idpoly(Ar, [], Cr);
model_poly.Structure.A.Free = Ar_free;
model_poly.Structure.C.Free = Cr_free;
model_rain = pem(u_s, model_poly);
res = resid(model_rain, u_s);
analyzets(res);
present(model_rain)
%% predicting rain
A = conv(model_rain.a, [1 zeros(1, 35) -1]);
[CS, AS] = equalLength(model_rain.c, A);
k = 10;
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
%ypw = ypw(find(~A, 1,'last'):end);
%upw = upw(find(~A, 1,'last'):end);
%%
%ypw = ypw(40:end);
%upw = upw(40:end);
%%
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
d = 4;
dr = 0;
B = zeros(1, 10);
B_free = B;
B_free(d:d+dr) = 1;
Mi = idpoly(1, B, [], [], A2);
Mi.Structure.F.Free = A2_free;
Mi.Structure.B.Free = B_free;
zpw = iddata(ypw, upw);
Mba2 = pem(zpw, Mi);
present(Mba2)
vhat = resid(Mba2, zpw);
analyzets(vhat)% Looks like ARMA, which is good I think?
%%
df = 4;
dfr = 0;
B_f = zeros(1, 10);
B_ff = B_f;
B_ff(df:df+dfr) = 1;

rc_f = Cr;
rc_ff = Cr_free;
rc_ff([]) = 1;
%rc_ff([3 5]) = 0;
rc_ff([4]) = 0;

ra_f = [1 zeros(1, 80)];
ra_ff = zeros(size(ra_f));
%ra_ff([2 3 35 36 37 38]) = 1;
ra_ff(find(A~=0)) = 1;
ra_ff([]) = 1;
%ra_ff([3 13 37 39 40 49 73 74 75]) = 0;
ra_ff([49 39 3]) = 0;

A2_f = [1 zeros(1, 10)];
A2_ff = A2_free;
A2_ff([2]) = 0;



Mi = idpoly(1, B_f, rc_f, ra_f, A2_f);
Mi.Structure.F.Free = A2_ff;
Mi.Structure.C.Free = rc_ff;
Mi.Structure.D.Free = ra_ff;
Mi.Structure.B.Free = B_ff;
z = iddata(y, u);
MboxJ = pem(z, Mi);
present(MboxJ);
ehat = resid(MboxJ, z);
%analyzets(ehat)
%%
ElGeneina_model = MboxJ;
save('ElGeneina_model.mat', 'ElGeneina_model')
%%
Af =  MboxJ.D;
Cf =  MboxJ.C;

A2f =  MboxJ.F;
Bf =  MboxJ.B;

%k_range = 1:1;
k_range = [1 4 6 8 10];
yhat_k = zeros(size(y_validation, 1), numel(k_range));
%bu = u;
%bu(numel(u)/2:end) = 0;
for ik = 1:numel(k_range)
    k = k_range(ik);
    [AS, CS] = equalLength(Af, Cf);
    [Fyk, Gyk] = deconv(conv([1,zeros(1, k-1)], CS), AS);

    [BFS, CS] = equalLength(conv(Bf, Fyk), Cf);
    [Fuk, Guk] = deconv(conv([1,zeros(1, k-1)], BFS), CS);

    yhat_k(:,ik) = filter(Gyk, Cf, y_validation) + filter(Guk, Cf, u_validation);
end
colors = flipud([254,237,222
253,190,133
253,141,60
230,85,13
166,54,3]./255);
%colors = colors(end-2,:);
figure;
hold on;
%caxis([0 1])
%ticks = (1/numel(k_range))*(k_range-k_range(1)) + (1/numel(k_range))/2;
%colorbar('Ticks',ticks,...
%     'TickLabels',k_range,...
%     'FontSize', 15);
%colormap(colors);
set(gca,'Color',[107,174,214]/255)
ylim([-3 5])
xlim([y_t(end-nbr_of_validation_points) y_t(end)])
set(gca, 'ColorOrder', colors);
plot(y_t, yhat_k, 'LineWidth', 1.5)
plot(y_t, y_validation, 'k', 'LineWidth', 2)
legend({'k=1','k=4','k=6','k=8','k=10','observed'},'Location','northwest')
legend('boxoff')
hold off;