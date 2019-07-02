function negloglike = gmrf_negloglike(theta, y, Atilde, C, G, G2, qbeta, Qmodel)
% GMRF_NEGLOGLIKE_SKELETON  Calculate the GMRF data likelihood, non-Gaussian observations
%
% negloglike = GMRF_negloglike_skeleton(theta, y, A, B, C, G, G2, qbeta, is_CAR)
%
% theta = [log(tau2) log(kappa2)]
% y = the data vector, as a column with n elements
% Atilde = the observation matrix, sparse n-by-(N+Nbeta)
% C,G,G2 = matrices used to build a Matern-like precision,
%          see matern_prec_matrices, sparse N-by-N
% qbeta = prior precision the regression parameters scalar
% isCAR = use CAR model (false=SAR)
%
% This is only a skeleton for Home Assignment 2.

% $Id: gmrf_negloglike_skeleton.m 4454 2011-10-02 18:29:12Z johanl $


tau = exp(theta(1));
kappa2 = exp(theta(2));
%xi = exp(theta(3));

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

%declare x_mode as global so that we start subsequent optimisations from
%the previous mode (speeds up nested optimisation).
global x_mode;
if isempty(x_mode)
  %no existing mode, compute a rough initial guess assuming Gaussian errors of log(y+0.1)
  x_mode = (Qtilde + Atilde'*Atilde)\(Atilde'*log(y+.1));
end

%find mode using Newton-Raphson
x_mode = fminNR(@(x) gmrf_taylor_po(x, y, Atilde, Qtilde), x_mode);

%find the Laplace approximation of the denominator
[g, ~, Q_xy] = gmrf_taylor_po(x_mode, y, Atilde, Qtilde);
%HINT: What does g contain?

%Compute choleskey factors
[R_x,p_x] = chol(Q_x);
[R_xy,p_xy] = chol(Q_xy);
if p_x~=0 || p_xy~=0
  %choleskey factor fail -> (almost) semidefinite matrix -> 
  %-> det(Q) ~ 0 -> log(det(Q)) ~ -inf -> negloglike ~ inf
  %Set negloglike to a REALLY big value
  %printf('Chol fail');
  negloglike = realmax; 
  return;
end

negloglike = g + sum(log(diag(R_xy))) - sum(log(diag(R_x)));

%print diagnostic/debug information (optimization progress)
fprintf(1, 'Theta: %11.4e %11.4e %11.4e; fval: %11.4e;\n', theta(1), theta(2),0, negloglike);