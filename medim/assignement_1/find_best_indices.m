function best_indices = find_best_indices(mx, my, use_scaling)

best_indices = [];
best_error = 100000;
for n = 1:2500
    i = randperm(size(mx, 2), 2); % get random indices
    %i = [6 7 8];
    [~, ~, t] = procrustes(mx(:,i)',my(:, i)', 'scaling', use_scaling); % do prcrustes with random indices
    if t.b < 0.01
        t.b = 1;
    end
    scale = 1;
    if use_scaling
        scale = t.b;
    end
    mt = affine2d([[scale * t.T [0; 0]]; [t.c(1,:) 1]]);
    et = mx' - transformPointsForward(mt, my');
    good_indices = find(sum(et .^ 2, 2) < 15); % number of points under treshhold
    good_error_tot = abs(sum(et(good_indices)));
    if (size(good_indices, 1) > size(best_indices, 1)) || (size(good_indices, 1) > 2 && size(good_indices, 1) == size(best_indices, 1) && best_error > good_error_tot && good_error_tot > 0)
        best_indices = good_indices;
        best_error = good_error_tot;
    end
end