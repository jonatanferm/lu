C = diag([Cii(0),Cii(1),Cii(2),Cii(3),Cii(4), Cii(5)])
G = [G]

function v = Cii(i)
    v = integral( @(x) triangularPulse(i-1,i+1,x), 0, 5);
end
function v = Gij(i, j)
    v = integral( @(x) triangularPulse(i-1,i+1,x)*triangularPulse(j-1,j+1,x), 0, 5);
end