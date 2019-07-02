function [ev,mult] = heltalsev(A,tol)
% HELTALSEV Find integer eigenvalues for A, with the tolerance tol
% defualt tolerance is 2
% By Jonatan Ferm jo1640fe-s 880302-4911
if nargin < 2, tol = 2; end
[~, D] = eig(A); %Find eigenvalues for A
D = round(diag(D) * 10^tol) / 10^tol;
%Checks if the eigenvalues are integers with our tolerance
%If they are not just return zeros.
if sum(floor(D) - D) ~= 0
    ev = 0;
    mult = 0;
else
    %Sort the eigenvalues and find the last index of each eigenvalue
    [ev, i] = unique(round(sort(D)), 'last');
    %The index of an eigenvalue minus the index of the last eigenvalue
    %will be the multiplicty
    mult = i - [0; i(1:end-1)];
end
end