r = 1;
mu = .2;
eps = .1;
a0 = 10;
k = 10;

a = @(x) a0 - .5.*k.*x.^2;
f = @(x1, x2) mu.*(a(x1)./a(x2)-1);

[x1, x2] = meshgrid(linspace(-1, 1));

pcolor(sign(f(x1, x2)))
colorbar
xlabel('x1')
ylabel('x2')