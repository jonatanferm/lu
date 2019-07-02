% pred.m	% make a separate function of the equation
function ndot=pred(t,n)
global k r a b d;
if n(1)<0
n(1)=0;
end
if n(2)<0
n(2)=0;
end
ndot=[r.*n(1).*(1- n(1)./k)-a.*n(1).*n(2);b.*a.*n(1).*n(2)-d*n(2)];