addpath("tsadl/matlab");
%%
%Generate data following BJ
n = 500;
A1 = [1 -.65];
A2 = [1 .9 .78];
C = 1;
B = [0 0 0 .4];
e = sqrt(1.5) * randn(n + 100, 1);
w = sqrt(2) * randn(n + 200, 1);
A3 = [1 .5];
C3 = [1 -.3 .2];
u = filter(C3, A3, w);
u = u(101:end);
y = filter(C, A1, e) + filter(B, A2, u);
u = u(101:end);
y = y(101:end);
clear A1 A2 C B e w A3 C3;
%% Task 1, ARMA model for u
%analyzets(u) %looks like it can be AR(2)
%ar2_model = arx(u, 2);
ar2_model = armax(u, [1 2]);
ar2_res = resid(ar2_model, u);
%upw = ar2_res.OutputData;
%analyzets(ar2_res); %Does look pretty good, everything is sort of on the edge
%whitenessTest(ar2_res.OutputData) %Super white
%normplot(ar2_res.OutputData); %Super normal
%ar2_model % prints A(z) = 1 + 0.7699 z^-1 + 0.1194 z^-2
%Q1: AR(2) was most suitable. Think it is reasonably close. 
%True is ARMA(1, 2), Higher A estimates the missing C. 
%% Pre white y
ypw = filter(ar2_model.a,ar2_model.c, y); %??? Is this right?
upw = filter(ar2_model.a,ar2_model.c, u);
%ypw = resid(ar2_model, y);
%ypw = ypw.OutputData;
%% plot crosscorr
M = 40;
stem(-M:M, crosscorr(upw, ypw, M));
hold on;
refline(0, 2/sqrt(size(ypw, 1)))
refline(0, -2/sqrt(size(ypw, 1)))
hold off %d = 2, s - r + 1 = 9, ordA2 = 2 ??
%% First none zero at 3. Looks like 8 values of random stuff and maybe then MA(1)
A2 = [1 0];
B = [0 0 0 1];
Mi = idpoly(1, B, [], [], A2);
zpw = iddata(ypw, upw);
Mba2 = pem(zpw, Mi);
%present(Mba2);
vhat = resid(Mba2, zpw);
%whitenessTest(vhat.OutputData)
%dagostinoK2test(vhat.OutputData)
%Q2 I found the above orders most suitable. Should it be white? I think it
%should be ARMA?
%% form x plot crosscorr. Crosscorr like unrelated for k>= 0. It should be wierd before 0 right?
x = y - filter(Mba2.b, Mba2.f, u);
M = 40;
stem(-M:M, crosscorr(x, u, M));
hold on;
refline(0, 2/sqrt(size(ypw, 1)))
refline(0, -2/sqrt(size(ypw, 1)))
hold off %d = 2, s - r + 1 = 9, ordA2 = 2 ??
%% Analyze x
%analyzets(x); % Maybe AR(4)
ar_4_x = arx(x, 4);
xres = resid(ar_4_x, x);
analyzets(xres);% looks good!
%% Final model
A2 = [1 0 0];
A1 = [1 0];
B = [0 0 0 1];
C = [1];
Mi = idpoly(1, B, C, A1, A2);
z = iddata(y, u);
MboxJ = pem(z, Mi);
%present(MboxJ);
ehat = resid(MboxJ, z);
%analyzets(ehat)
%%