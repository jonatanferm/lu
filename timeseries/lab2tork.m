clear;
load("/Users/ferm/lu/timeseries/tsadl/data/tork.dat");
tork = tork - mean(tork);
y = tork(:,1);
u = tork(:,2);
z = iddata(y, u);
%plot(z)
%%
%analyzets(u) %looks like AR(1)
arumodel = arx(u, 1);
ures = resid(arumodel, u);
%analyzets(ures)% looks good. Really not gaussian though
%prewhite
ypw = filter(arumodel.a, 1, y);
upw = filter(arumodel.a, 1, u);
%% plot crosscorr. Looks like d = 2 maybe r = 1 s =2
M = 40;
stem(-M:M, crosscorr(upw, ypw, M));
hold on;
refline(0, 2/sqrt(size(ypw, 1)))
refline(0, -2/sqrt(size(ypw, 1)))
hold off
%%
A2 = [1 0];
B = [0 0 0 0 0 0];
Mi = idpoly(1, B, [], [], A2);
Mi.Structure.f.Free = [0 1];
Mi.Structure.b.Free = [0 0 0 1 1 1];
zpw = iddata(ypw, upw);
Mba2 = pem(zpw, Mi);
vhat = resid(Mba2, zpw);
analyzets(vhat)% Looks like ARMA, which is good I think?
%% form x and plot xcorr
x = y - filter(Mba2.b, Mba2.f, u);
M = 40;
stem(-M:M, crosscorr(x, u, M));
hold on;
refline(0, 2/sqrt(size(ypw, 1)))
refline(0, -2/sqrt(size(ypw, 1)))
hold off %This is not white :( Maybe it should be since data not normal?
%%
%analyzets(x); %Look like AR(1)
ar_x = arx(x, 1);
xres = resid(ar_x, x);
analyzets(xres);% looks good!
%% Final model
A2 = [1 0];
A1 = [1 0];
B = [0 0 0 1 1 1];
C = [1];
Mi = idpoly(1, B, C, A1, A2);
z = iddata(y, u);
MboxJ = pem(z, Mi);
%present(MboxJ);
ehat = resid(MboxJ, z);
analyzets(ehat) %crosscorr does not look good
%%
