load HA2_forest
Y = bei_counts(:);
I = ~isnan(Y);

obs_indices = find(I);
vt_indices = randperm(size(obs_indices, 1), size(obs_indices, 1)*.2);

test_observations = sort(obs_indices(vt_indices(1:size(vt_indices, 2)/2)));
val_obervations = sort(obs_indices(vt_indices(size(vt_indices, 2)/2+1:end)));

test_I = zeros(size(I));
val_I = zeros(size(I));

test_I(test_observations) = 1;
val_I(val_obervations) = 1;

I = boolean([I - test_I - val_I, val_I, test_I]);
%save('observations.mat', 'I')
