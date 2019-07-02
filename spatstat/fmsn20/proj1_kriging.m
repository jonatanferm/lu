y_hat = X * beta;
y_hat(y_hat<0) = 0;
e = [Y(:,1) - y_hat, Y(:,3) Y(:,4)];
%{
figure; plot(distance_matrix(e(:,2:3)), e(:,1)*e(:,1)','.k');
bins = 10;
md = 70;
[rhatr,s2hatr,m,n,d]=covest_nonparametric(e(:,2:3),e(:,1),bins,md);
nbr_of_perms = 1000;
rhatm = zeros(nbr_of_perms, size(rhatr, 2));
for i = 1:nbr_of_perms
   [rhat,s2hat,m,n,d]=covest_nonparametric(e(:,2:3),e(randperm(size(e, 1)),1),bins,md);
   rhatm(i,:) = rhat;
end

figure;
plot(d, rhatr, '-', ...
    d, prctile(rhatm, 95), '-r', ...
    d, prctile(rhatm, 5), '-r')

%}
dists_oo_est = distance_matrix(e(:,2:3));
dists_oo_est(dists_oo_est > 70) = 10^6;
[par, ~] = covest_ml(dists_oo_est,e(:,1), 'matern', [0 0 0.8 0]);
rhat_p = matern_covariance(d, par(1), par(2), par(3));
figure;
plot(d, rhatr, '-b', ...
    0, s2hatr, 'ob', ... 
    d, rhat_p, '-r')

%[par, ~] = covest_ml(e(:,2:3),e(:,1), 'cauchy', []);

swiss_coords = cat(3, swissX, swissY);
swiss_coords = reshape(swiss_coords(~isnan(swiss_coords)), [], 2);
dists_uu = distance_matrix(swiss_coords);

dists_oo = distance_matrix(e(:,2:3));
dists_uo = distance_matrix(swiss_coords, e(:,2:3));
dists_ou = distance_matrix(e(:,2:3), swiss_coords);

dists_ov = distance_matrix(e(:,2:3), Yvalid(:,3:4));
dists_vo = distance_matrix(Yvalid(:,3:4), e(:,2:3));
dists_vv = distance_matrix(Yvalid(:,3:4));

Sigma_uu = matern_covariance(dists_uu, par(1), par(2), par(3)) + eye(size(dists_uu, 1)) * par(4);
Sigma_uo = matern_covariance(dists_uo, par(1), par(2), par(3));
Sigma_ou = matern_covariance(dists_ou, par(1), par(2), par(3));
Sigma_oo = matern_covariance(dists_oo, par(1), par(2), par(3)) + eye(size(dists_oo, 1)) * par(4);
Sigma_vo = matern_covariance(dists_vo, par(1), par(2), par(3));
Sigma_ov = matern_covariance(dists_ov, par(1), par(2), par(3));
Sigma_vv = matern_covariance(dists_vv, par(1), par(2), par(3)) + eye(size(dists_vv, 1)) * par(4);

y_o = Y(:,1);
X_u = Xgrid;
X_o = X;

xst = X_u' - X_o' * (Sigma_oo \ Sigma_ou);
V_y = Sigma_uu - Sigma_uo * (Sigma_oo \ Sigma_ou) + ...
    xst' * ((X_o' * (Sigma_oo \ X_o)) \ xst);
V_y = diag(V_y);
V_pic = nan(size(swissElevation));
V_pic(~isnan(swissElevation)) = V_y;

xst_v = X_v' - X_o' * (Sigma_oo \ Sigma_ov);
V_yv = Sigma_vv - Sigma_vo * (Sigma_oo \ Sigma_ov) + ...
    xst_v' * ((X_o' * (Sigma_oo \ X_o)) \ xst_v);
V_yv = diag(V_yv);

beta_hat = (X_o' * (Sigma_oo \ X_o)) \ (X_o' * (Sigma_oo \ y_o));

y_rec = X_u * beta_hat + Sigma_uo * (Sigma_oo \ ...
    (y_o - X_o * beta_hat));
y_rec_valid = X_v * beta_hat + Sigma_vo * (Sigma_oo \ ...
    (y_o - X_o * beta_hat));
y_rec(y_rec < 0) = 0;

y_rec_image = nan(size(swissElevation));
y_rec_image(~isnan(swissElevation)) = y_rec;
figure;
imagesc([0 max(swissX(:))], [0 max(swissY(:))], y_rec_image, ...
        'alphadata', ~isnan(y_rec_image));
hold on
plot(swissBorder(:,1), swissBorder(:,2),'k')
scatter(Y(:,3), Y(:,4), 200, Y(:,1), 'filled','markeredgecolor','r', 'LineWidth',.5)
scatter(Yvalid(:,3), Yvalid(:,4), 200, Yvalid(:,1), 'filled','markeredgecolor','w')
axis xy tight; hold off; colorbar
title('sqrt of rainfall predictions')

figure;
hold on
imagesc([0 max(swissX(:))], [0 max(swissY(:))], V_pic, ...
        'alphadata', ~isnan(V_pic));
plot(swissBorder(:,1), swissBorder(:,2),'k');
scatter(Y(:,3), Y(:,4), 20, 'markeredgecolor','r', 'LineWidth',.5)
axis xy tight; hold off; colorbar
%title('Variance of estimates')

validation_error = abs((Yvalid(:,1) - y_rec_valid)) ./ Yvalid(:,1);
mean_ve = mean(validation_error)
mean_vv = mean(V_yv)
par



