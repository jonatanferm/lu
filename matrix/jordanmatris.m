function J = jordanmatris(A,tol)
% JORDANMATRIS finds the jordan normal form of a matrix A
% Must have integer (with tolerance tol) eigenvalues or this will return 0
% By Jonatan Ferm jo1640fe-s 880302-4911
if nargin < 2, tol = 2; end
[ev,mult] = heltalsev(A,tol);
function block = createblock(r)
% CREATEBLOCK creates a jordan block
% r should be an array with the first value as the eigenvalue
% and the second value as the multiplicity
    h = r(1);
    m = r(2);
    block = eye(m)*h + [zeros(m-1, 1), eye(m-1); zeros(1, m)];
end
if ev ~= 0
    %create a cell for each eigenvalue/multiplicty pair
    c = num2cell([ev, mult].', 1);
    %call the create block function for each pair
    blocks = cellfun(@createblock, c, 'un', 0);
    %create a matrix of all the blocks.
    J = blkdiag(blocks{:});
else
    J = 0;
end
end
