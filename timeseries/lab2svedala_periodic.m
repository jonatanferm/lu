load("./tsadl/data/svedala");
%y = svedala;
per = [1 zeros(1, 23) -1];
y = filter(per, 1, svedala);

model_init = idpoly([1 zeros(1, 2)], [], [1 zeros(1, 23) 1]);
model_init.Structure.a.Free = [1 1 1];
model_init.Structure.c.Free = [zeros(1, 24) 1];
model_armax = pem(y, model_init);
res = resid(model_armax, y);

%analyzets(res);

A = model_armax.a;
C = model_armax.c;

A = conv(A, per);

k = 3;
[CS, AS] = equalLength(C, A);
[Fk, Gk] = deconv(conv([1,zeros(1, k-1)], CS), AS);

svedala_k = filter(Gk, C, svedala);
figure;
hold on;
plot(svedala)
plot(svedala_k)
%plot(y(k+1:end))
%plot(yhat_k(1:end-k))
hold off

e = svedala - svedala_k;
mean(e)
var(e)
