	% icomp.m
	% functions for n(1) species 1 and n(2) species 2, as a function of time

	function ndot=icomp(t,n)
	global r k alfa beta;
	ndot=[r(1).*n(1).*((k(1)-n(1)-alfa*n(2))/k(1));r(2).*n(2).*((k(2)-n(2)-beta*n(1))/k(2))];