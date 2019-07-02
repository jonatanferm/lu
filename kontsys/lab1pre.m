n = 8; % 256 punkter
N=2^n;
K=1:(N/2-1);

L=1;
delta=L/(N-1);
X=0:delta:(N-1)*delta;

funktion1 = 'max(sin(X*2*pi/L),0)';
funktion2 = '(X < L/2) - (X > L/2)' ;
funktion3 = 'X.*(X < L/4) + (L-X).*(X >= L/4)/3' ;
funktion4 ='(X<L/4).*X - (X>2*L/3).*(X<3*L/4)';
funktion5='sin(10*pi*X).*(X>0.4*L).*(X<0.5*L)'; 


anoll = 0;
akoeff = zeros(size(K));
bkoeff = ones(size(K))./K; 