% Load data
clear
clc
load fmri.mat

%size of data
sz = size(img);

%Option 1: regress onto indicator functions
beta = X\colstack(img)';
%reshape back to an image where the beta-coefs are in each "color"-layer
beta = reshape(beta', sz(1), sz(2), []);
[y_beta, ~, P_beta] = pca(colstack(beta));
y_beta = reshape(y_beta, sz(1), sz(2), []);
modes = 2:4;
y_beta = y_beta(:,:,modes);
y_beta_flat = reshape(y_beta, [], size(y_beta, 3));

rgb = [ ...
    103,0,31
178,24,43
214,96,77
244,165,130
253,219,199
247,247,247
209,229,240
146,197,222
67,147,195
33,102,172
5,48,97  ] / 255;
K = 6;
bsz = size(beta);
km = kmeans(y_beta_flat, K);
figure;
imagesc(reshape(km, sz(1:2)));
colormap(rgb)

%%
K = 6;
[theta,prior, p]=normmix_gibbs(rs,K);
[pm, pi] = max(p, [], 2);
p = reshape(p, [sz(1:2) K]);
pm = reshape(pm, [sz(1:2) 1]);
pi = reshape(pi, [sz(1:2) 1]);
figure;
imagesc(pi);
colormap(rgb)
colorbar
figure;
imagesc(pm, [0 1]);
colormap(gray())
colorbar
%%
%N = [0 -1 0; -1 4+0.1^2 -1; 0 -1 0];
N = [0 1 0; 1 0 1; 0 1 1];
%Q = gmrfprec([size(beta, 1) size(beta, 2)],N);
%s = Q \ sparse(size(beta, 2)*size(beta, 1)/2+size(beta, 2)/2,1,1,size(Q,1),1);
%s = reshape(s, [size(beta, 1), size(beta, 2)]);
%%
y_beta = reshape(y_beta, sz(1), sz(2), []);
beta = zeros(1,length(theta));
%beta = 1;
alpha = log(prior/prior(1));
alpha_post = mrf_gaussian_post(alpha, theta, y_beta);
[z,Mz,Mf,Mzf] = mrf_sim(p, 10, alpha_post, beta, 1);
[~, zi] = max(z, [], 3);

gibbs_mu_sigma(z);

%{
figure()
axis tight
colormap(rgb)
colorbar
subplot(1, 3, 1)
imagesc(zi)
subplot(1, 3, 2)
imagesc(pi)
subplot(1, 3, 3)
imagesc(zi-pi)
%}
%%
for i=1:10
    alpha_post = mrf_gaussian_post(alpha_post, theta, y_beta);
    [z,Mz,Mf,Mzf]= mrf_sim(z, 3, alpha_post, beta, 1000);
    [alpha, beta, ~] = gibbs_alpha_beta([], beta, z, 3, 10);
    image(rgbimage(z))
    %imagesc(mat2gray(alpha_post))
    drawnow
end
%%
alpha_post = mrf_gaussian_post(beta, theta, y_beta);
[alpha, beta, acc] = gibbs_alpha_beta([], beta, alpha_post, N, 10);
z = mrf_sim(alpha_post2, N, [0; alpha], beta, 100);
imagesc(z(:,:,3));
%%
%help gibbs_mu_sigma
%%
    N = [0 1 0;1 4 1;0 1 0];
    z = mrf_sim(zeros(100,120,3),N,log([0.2 0.3 0.5]),0.5,100);
    for iter=1:200
      z = mrf_sim(z,N,log([0.2 0.2 0.6]),[0.9 0.9 0],1);
      image(rgbimage(z))
      drawnow
    end

%%
%Option 2: Compute SVD directly on the data (this is essentially the SVD
%example from lecture 9)
[y,V,P] = pca(colstack(img));
y = reshape(y, sz(1), sz(2), []);

%study the temporal components to find those with 20s periodicity
figure('Name','PCA directly')
subplot(3,4,1)
psp = P/sum(P);
semilogy(1:20,psp(1:20))
axis tight
for i=1:11
  subplot(3,4,i+1)
  plot(V(:,i))
  axis tight
  title(i)
end
%%
figure;
imagesc(y(:,:,11))
colormap(gray())
%%
%alternatively using a periodogram (stationary/timeseries course)
figure('Name','PCA periodogram')
for i=1:8
  subplot(4,4,2*i-1)
  plot(V(:,i))
  axis tight
  title(i)
  subplot(4,4,2*i)
	periodogram(V(:,i))
	line([0.1 0.1], get(gca,'YLim'));
	xlabel(''); ylabel(''); title(i);
end

%and plot a few leading components
figure('Name','leading components')
subplot(3,4,1)
semilogy(P/sum(P))
axis tight
for i=1:11
  subplot(3,4,i+1)
  imagesc(y(:,:,i))
  title(i)
end