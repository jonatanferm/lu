load('tsadl/data/tar2.dat');
load('tsadl/data/thx.dat');
%%
figure;plot(tar2)
figure; plot(thx)
%%
model = [2];
lambda = 0.95;
%%
[Aest, yhat, covAest, yprev] = rarx(tar2, model, 'ff', lambda);
%%
figure;
hold on;
plot(thx)
plot(Aest)
hold off
%%
n = 100;
lambda_line = linspace(0.85, 1, n);
ls2 = zeros(n, 1);
for i = 1:length(lambda_line)
    [Aest, yhat, covAest, ~] = rarx(tar2, model, 'ff', lambda_line(i));
    ls2(i) = sum((tar2-yhat).^2);
end
%%
[y, i] = min(ls2);
%%

y = ?;
N = length(y);
A = [? ?; ? ?];
Re = [? ?; ? ?];
Rw = ?;
Rxx_1 = ? * eye(2);
xtt_1 = [? ?]';
xsave = zeros(2,N);
for k=3:N
C = [? ?];
Ryy = ?;
Kt = ?;
xtt = ?;
Rxx = ?;
xsave(:,k) = ?;
Rxx_1 = ?;
xtt_1 = ?;
end



