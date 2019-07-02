

N = 10000;
gens = nan(N, 1);

for i = 1:N
    gens(i) = Perceptron;
end


std(gens)
CI = norminv(.975)*std(gens)/sqrt(N);
[mean(gens) - CI, mean(gens) + CI]

    