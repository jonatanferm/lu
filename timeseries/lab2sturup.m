load("./tsadl/data/sturup")
load("./tsadl/data/svedala");
y = svedala;
u = sturup;
A = [1 -1.49 0.57];
B = [0 0 0 .28 -.26];
C = 1;
%Q12 3

k = 3;
[AS, CS] = equalLength(A, C);
[Fyk, Gyk] = deconv(conv([1,zeros(1, k-1)], CS), AS);

[BFS, CS] = equalLength(conv(B, Fyk), C);
[Fuk, Guk] = deconv(conv([1,zeros(1, k-1)], BFS), CS);

yhat_k = filter(Gyk, C, y) + filter(Guk, C, u);

var(y - yhat_k)
%Q13 variance is  2.0016