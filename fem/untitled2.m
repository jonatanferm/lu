N = 1000;
n = 1000; %hits

b = .05;
inc = 1;
f = b*(1+inc) * ones(1, n);

o = zeros(N, n);

for i = 1:N
    o(i,:) = rand(1, n) < f;
end

[r, c, v] = find(cumsum(o) == 20);