%load data
load HA2_forest
load observations.mat

%figure; imagesc( bei_counts )
%figure; imagesc( reshape(I(:,1), size(bei_elev)) )
%figure; imagesc( bei_elev )
%figure; imagesc( bei_grad )

%size of the grid
sz = size(bei_counts);
%observations
Y = bei_counts(:);
%missing data
I_val = I(:,2);
I_test = I(:,3);
I = I(:,1);

%create Q-matrix
[u1, u2] = ndgrid(1:sz(1),1:sz(2));
[C,G,G2] = matern_prec_matrices([u1(:) u2(:)]);
%mean value-vector (might not need all)
%Bgrid = [ones(prod(sz),1) bei_elev(:) bei_grad(:)];
%Bgrid = [ones(prod(sz),1) bei_elev(:)            ];
%Bgrid = [                 bei_elev(:) bei_grad(:)];
%Bgrid = [ones(prod(sz),1)             bei_grad(:)];
%Bgrid = [ones(prod(sz),1)                        ];
%Bgrid = [                 bei_elev(:)            ];
%Bgrid = [                             bei_grad(:)];
%Bgrid = [];
Bgrid = [ones(prod(sz),1) bei_elev(:) bei_elev(:).^2];

%and observation matrix for the grid
Agrid = speye(prod(sz));

%G2 is the most dense of the matrices, lets reorder
p = amd(G2);
%{
figure
subplot(121)
spy(G2)
subplot(122)
spy(G2(p,p))
%}
%reorder precision matrices
C = C(p,p);
G = G(p,p);
G2 = G2(p,p);
%and observation matrix
Agrid = Agrid(:,p);

%create A tilde matrix
Atilde = [Agrid Bgrid]; 
%Atilde = [Agrid]; 
%Atilde = [Bgrid]; 

%we need a global variable for x_mode to reuse
%between optimisation calls
global x_mode;
x_mode = [];

Qmodel = 'CAR';
%Qmodel = 'SAR';
%Qmodel = 'OSC';
qbeta = 1e-6;
%subset Y and Atilde to observed points
par = fminsearch( @(x) gmrf_negloglike(x, Y(I), Atilde(I,:), C, G, G2, qbeta, Qmodel), [0 0]);
%conditional mean is given by the mode  
E_xy = x_mode;
%and reconstruction (field+covariates)
E_zy = Atilde*x_mode;
%imagesc( reshape(E_zy,sz) )

tau = exp(par(1));
kappa2 = exp(par(2));
%xi = exp(par(3));

%comput Q for a CAR(1) or SAR(1) process
switch Qmodel
    case 'CAR'
        Q_x = tau * (kappa2 * speye(size(G, 1)) + G);
    case 'SAR'
        Q_x = tau * (kappa2^2*speye(size(G, 1)) + 2*kappa2*G + G2);
    case 'OSC'
        Q_x = tau * (kappa2^2*speye(size(G, 1)) + 2*xi*kappa2*G + G2);
end

%combine Q_x and qbeta
Nbeta = size(Atilde,2) - size(Q_x,1);
Qtilde = blkdiag(Q_x, qbeta*speye(Nbeta));

%reuse taylor expansion to compute posterior precision
[~, ~, Q_xy] = gmrf_taylor_po(E_xy, Y(I), Atilde(I,:), Qtilde);

%1000 samples from the approximate posterior
Rxy = chol(Q_xy);
x_samp = bsxfun(@plus, E_xy, Rxy\randn(size(Rxy,1),1000));
stand_err = sqrt(var(Atilde*x_samp, 0, 2));
%%
[y, x] = ind2sub(sz, find(I));
figure; 
hold on
imagesc(reshape(stand_err, sz));
plot(x + 1i*y, '.k');
hold off
%%
figure; plot(reshape(bei_elev, [], 1), stand_err, '.');

AO = [Agrid, zeros(size(Bgrid))];
OB = [zeros(size(Agrid)), Bgrid];

E_spat = reshape(AO*x_mode,sz);
E_mc = reshape(OB*x_mode,sz);

figure; imagesc(reshape(E_spat, sz));
figure; imagesc(reshape(E_mc, sz));
%{
figure;
subplot(311)
imagesc( E_spat )
caxis([-5 5]);
subplot(312)
imagesc( E_mc )
caxis([-5 5]);
subplot(313)
imagesc(reshape(E_zy,sz))
caxis([-5 5]);

e = [zeros(size(Q_xy,1)-size(Bgrid,2), size(Bgrid,2)); eye(size(Bgrid,2))];
V_beta0 = diag(e'*(Q_xy\e));

mean((Y(I_val) - exp(E_zy(I_val))).^2)
%mean(stand_err(I_val))
%}
%%
figure;
colormap(hot(25))
hist2([reshape(bei_counts(I_test), [], 1), round(exp(E_zy(I_test)))], 5)


