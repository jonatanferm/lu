load HA2_forest
I = ~isnan(bei_counts(:));
plot(reshape(bei_elev(I), [], 1), reshape(bei_counts(I), [], 1), '.')
figure; plot(reshape(bei_grad(I), [], 1), reshape(bei_counts(I), [], 1), '.')
