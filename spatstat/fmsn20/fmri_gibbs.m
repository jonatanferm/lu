clear
clc
load fmri.mat

sz = size(img);
beta = X\colstack(img)';
beta = reshape(beta', sz(1), sz(2), []);
[y_beta, ~, P_beta] = pca(colstack(beta));

modes = 2:4;
mdim = numel(modes);
y_beta_flat = y_beta(:,modes);

y_beta_sq = reshape(y_beta_flat, sz(1), sz(2), mdim);

K = 3;
[theta, prior, p_flat]=normmix_gibbs(y_beta_flat,K);

p_sq = reshape(p_flat, sz(1), sz(2), K);

N1 = [0 1 0;
      1 0 1;
      0 1 0];
N2 = [1 1 1;
      1 0 1;
      1 1 1]; 
N = N2;
beta_prior = 10;
beta = zeros(1,length(theta));
alpha = log(prior/prior(1));
alpha_post = mrf_gaussian_post(alpha, theta, y_beta_sq);

[z,Mz,Mf,Mzf] = mrf_sim(p_sq, N, alpha_post, beta, 100);
%figure('Name', 'mrf')
%imagesc(z2img(z));
zf = 0;
figure('Name', 'dynamic')
for iter = 1:1000
    theta0 = cell(1,K);
    zf_new = reshape(z, [], K);
    if norm(zf_new - zf) < 1
        break
    end
    zf = zf_new;
    ap_flat = reshape(alpha_post, [], K);
    for k=1:K
        [theta0{k}.mu, theta0{k}.Sigma] = gibbs_mu_sigma(y_beta_flat(zf(:,k) == 1, :));
    end
    alpha_post = mrf_gaussian_post(alpha, theta0, y_beta_sq);
    z = mrf_sim(z, N, alpha_post, beta, 1);
    [alpha0, beta, acc] =  gibbs_alpha_beta(alpha(2:end), beta, z, N, beta_prior);
    alpha = [0; alpha0];
    %alpha_post = mrf_gaussian_post(alpha, theta0, y_beta_sq);
    %z = mrf_sim(z, N, alpha_post, beta, 1);
    image(rgbimage(z))
    drawnow
end
figure('Name', 'normmix')
imagesc(z2img(p_sq))
%figure('Name', 'mrf2')
%imagesc(z2img(z));
function img = z2img(z)
    [~, img] = max(z, [], 3);
end

