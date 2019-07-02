clear
sz = [60 60];

[X Y] = meshgrid(0:sz(1)-1, 0:sz(2)-1);
coords = cat(3, X, Y);
coords = reshape(coords, [], 2);
dists = distance_matrix(coords);
mu = 2;
N = size(coords, 1);

Sigma = matern_covariance(dists, 0.3, 0.2, 2);
%R = chol(Sigma); % Calculate the Cholesky factorisation
R = chol(Sigma + eye(size(Sigma))*1e-5);
eta = mu+R'*randn(N,1); % Simulate a sample
eta_image = reshape(eta,sz); %reshape the column to an image
figure; imagesc(eta_image);
se = 0.4;
y = eta + randn(N, 1)*se;
z=y-mu*1.1;
%figure; plot(Sigma, dists, '.k');
%%

%I_obs = (rand(sz)<=0.01);
I_obs = false(sz);
sd = 5;
I_obs(1:sd:end,1:sd:end) = true;


[xc_obc, yc_obc] = meshgrid(0:sd:sz(1)-1, 0:sd:sz(2)-1);
obs_coords = cat(3, xc_obc, yc_obc);
obs_coords = reshape(obs_coords, [], 2);
obs_dists = distance_matrix(obs_coords);

zobs = z(I_obs)
figure; plot(obs_dists, zobs*zobs','.k');

bins = 10;
md = 100;
[rhatr,s2hatr,m,n,d]=covest_nonparametric(obs_dists,zobs,bins,md);
rhatm = zeros(100, size(rhatr, 2));
for i = 1:100
   [rhat,s2hat,m,n,d]=covest_nonparametric(obs_dists, ...
       zobs(randperm(size(zobs, 1)),1),bins,md);
   rhatm(i,:) = rhat;
end

figure; plot(d,rhatr,'-', d, rhatm, '.',0,s2hat,'o');


par = covest_ml(obs_coords ,z(I_obs), 'matern', []);

Sigma_hat = matern_covariance(dists, par(1), par(2), par(3));

Sigma_yy = Sigma_hat;
Sigma_uu = Sigma_yy(~I_obs, ~I_obs);
Sigma_uo = Sigma_yy(~I_obs, I_obs);
Sigma_oo = Sigma_yy(I_obs, I_obs);
y_o = y(I_obs);
y_u = y(~I_obs);
X = ones(prod(sz),1);
X_u = X(~I_obs);
X_o = X(I_obs);
y_rec = nan(sz);
y_rec(I_obs) = y_o;
beta_hat = (X_o' * (Sigma_oo \ X_o)) \ (X_o' * (Sigma_oo \ y_o));
y_rec(~I_obs) = X_u * beta_hat + Sigma_uo * (Sigma_oo \ ...
    (y_o - X_o * beta_hat));
y_rec_image = reshape(y_rec,sz); %reshape the column to an image
figure; imagesc(y_rec_image);

