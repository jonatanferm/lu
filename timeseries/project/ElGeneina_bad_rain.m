clear;
load('proj18.mat');
addpath('/Users/ferm/lu/timeseries/tsadl/matlab')
%%
y_raw = ElGeneina.nvdi;
y_t = ElGeneina.nvdi_t;
u_raw = ElGeneina.rain;
u_t = ElGeneina.rain_t;
%%
y_raw = y_raw - mean(y_raw);
u_raw = u_raw - mean(u_raw);
offset = min([u_raw; y_raw]) - 1;
y_raw = y_raw - offset;
u_raw = u_raw - offset;
%%
start_index = find(u_t == y_t(1));
u_raw = u_raw(start_index:end);
u_t = u_t(start_index:end);
%%
lambda_u = -0.2230;
u_tform = (1/lambda_u)*(u_raw.^lambda_u-1);
lambda_y = -0.6571;
y_tform = (1/lambda_y)*(y_raw.^lambda_y-1);
%y_tform = y_raw;
figure; 
hold on;
plot(u_tform);
plot(y_tform);
hold off;
%%
myt = mean(y_tform);
y_tform = (y_tform - myt) * (std(u_tform) / std(y_tform)) + myt;
%% remove period of 3*12
period = [1 zeros(1, 3*12-2) 0 -1];
y = y_tform;%filter(period, 1, y_tform);
u = u_tform;%filter(period, 1, u_tform);
%y = y(37:end);
%u = u(37:end);
%%
figure;
hold on;
plot(u);
plot(y);
hold off;
%%
analyzets(u, 100);
%%
ra = [1 zeros(1, 40)];
ra_free = zeros(size(ra));
ra_free([2 15 16 36 37 38]) = 1;

rc = [1 zeros(1, 40)];
rc_free = zeros(size(rc));
rc_free([1 2 3]) = 1;

model_poly = idpoly(ra, [], rc);
model_poly.Structure.a.Free = ra_free;
model_poly.Structure.c.Free = rc_free;
model_rain = pem(u, model_poly);
res = resid(model_rain, u);
analyzets(res);
figure; pzmap(model_rain.a, model_rain.c);
%%
[CS, AS] = equalLength(model_rain.c, model_rain.a);
k = 1;
[Fk, Gk] = deconv(conv([1,zeros(1, k-1)], CS), AS);
u_k = filter(Gk, model_rain.c, u);
%u_k = u_k(length(Gk)+1:end);
figure;
hold on;
plot(u)
plot(u_k, '--')
hold off
%%
ypw = filter(model_rain.a, model_rain.c, y);
upw = filter(model_rain.a, model_rain.c, u);
ypw = ypw(25:end);
upw = upw(25:end);
%%
figure;
hold on;
plot(upw)
plot(ypw)
hold off;
%%
analyzets(upw)
%%
M = 80;
stem(-M:M, crosscorr(upw, ypw, M));
hold on;
refline([0 2/sqrt(size(u, 1))])
refline([0 -2/sqrt(size(u, 1))])
hold off
%%
A2 = [1 0];
A2_free = [1 1];
d = 4;
dr = 1;
B =      [0 0 0 0 0 0 0];
B_free = [0 0 0 0 0 1 1];
Mi = idpoly(1, B, [], [], A2);
Mi.Structure.F.Free = A2_free;
Mi.Structure.B.Free = B_free;
zpw = iddata(ypw, upw);
Mba2 = pem(zpw, Mi);
vhat = resid(Mba2, zpw);
analyzets(vhat)% Looks like ARMA, which is good I think?
%%
%% Final model
B_f = zeros(1, 7);
B_ff = zeros(1, 7);
B_ff([6]) = 1;

rc_f = [1 zeros(1, 40)];
rc_ff = zeros(size(rc));
rc_ff([1 2 3]) = 1;

ra_f = [1 zeros(1, 40)];
ra_ff = zeros(size(ra));
ra_ff([2]) = 1;

A2_f = [1 zeros(1, 10)];
A2_ff = zeros(size(A2_f));
A2_ff([1 2]) = 1;

Mi = idpoly(1, B_f, rc_f, ra_f, A2_f);
Mi.Structure.F.Free = A2_ff;
Mi.Structure.C.Free = rc_ff;
Mi.Structure.D.Free = ra_ff;
Mi.Structure.B.Free = B_ff;
z = iddata(y, u);
MboxJ = pem(z, Mi);
%present(MboxJ);
ehat = resid(MboxJ, z);
analyzets(ehat)
%%
A = 3;
k = 3;
[AS, CS] = equalLength(A, C);
[Fyk, Gyk] = deconv(conv([1,zeros(1, k-1)], CS), AS);

[BFS, CS] = equalLength(conv(B, Fyk), C);
[Fuk, Guk] = deconv(conv([1,zeros(1, k-1)], BFS), CS);

yhat_k = filter(Gyk, C, y) + filter(Guk, C, u);

