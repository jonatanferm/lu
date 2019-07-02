A = [1 -1.5 0.7];
C = [1 zeros(1, 11) - 0.5];
A12 = [1 zeros(1, 11) -1];
A_star = conv(A, A12) ;
e = randn(600, 1);
y = filter(C, A_star, e);
y = y(100:end);
%%
plot(y)
%%
S5 = [1 zeros(1, 4) -1];
y_star = conv(y, S5);
model_init = idpoly([1 zeros(1, 3)], [], [1 zeros(1, 12)]);
model_init.Structure.a.Free = [0 1 0 1];
model_init.Structure.c.Free = [0 1 zeros(1, 10) 1];
model_armax = pem(y_star, model_init);
a1 = model_armax.a(2)
a3 = model_armax.a(4)
c1 = model_armax.c(2)
c12 = model_armax.c(13)
