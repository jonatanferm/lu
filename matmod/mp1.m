
stds = [44.6 20.4 42.5;
        52.0 45.4 21.3;
        33.9 42.7 38.5;
        29.1 23.7 37.5];
means = [54.1 18.3 34.9;
        47.9  96.1 40.6;
         43.3 49.7 57.8;
         101.7 67.6  80.9];
    
    
counts = [9 8 8;
          8 9 8;
          8 8 9;
          8 8 8];
    
    
tot_stds = sqrt(sum(stds.^2 .* counts, 2)./sum(counts, 2));
tot_means = sum(means .* counts, 2)./sum(counts, 2);
CI = tot_stds .* 2./sqrt(sum(counts, 2));
tot_ci = [tot_means + CI, tot_means - CI]


plot([20, 40, 80, 160], tot_ci, 'LineWidth', 1.5, 'color', 'k')
set(gca, 'XScale', 'log')
xticks([20, 40, 80, 160])
xlabel('Dos mg')
ylabel('TSS')

