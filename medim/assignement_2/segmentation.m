
kk = 16;
%for i = 1:5
%for i = 6:10
%for i = 11:15
for i = 16:20
    max_modes = 15;
    image = dmsa_images(:,:,i);
    tresh_hold = max(max(image)) / 4;
    
    [left_kidney, right_kidney, middle] = split_image(image, tresh_hold);
    left_kidney = fliplr(left_kidney);
    kidneys = {left_kidney, right_kidney};
    for cur_kid_index = 1:2
        cur_kid = kidneys{cur_kid_index};
        fig = figure('Visible', 'Off');
        hold on;
        set(gca, 'ColorOrder', interpolate_colors(get_colomap(2), max_modes));
        imagesc(cur_kid);
        colormap(gray(255))
        %colormap(interpolate_colors(get_colomap(4), 255));
        for num_modes = 1:max_modes
            cur_basis = full_basis(:,end-num_modes+1:end);
            [fit_shape, initial_b] = shape_fitter(cur_kid, tresh_hold, mean_shape, cur_basis);
            fit_shape_closed = [fit_shape;fit_shape(1)];
            plot(fit_shape_closed, 'LineWidth', 1.5)
         
            
            %%Things to validate models
            %val_model = val_models(:,i-20*(cur_kid_index-2));
            %fit_shape + middle
            
            %val_poly = polyshape(real(val_model - middle), imag(val_model));
            %fit_poly = polyshape(real(fit_shape), imag(fit_shape));
            
            %eval_m(i, 1, num_modes) = area(intersect(val_poly, fit_poly)) / area(val_poly);
            %eval_m(i, 2, num_modes) = area(xor(val_poly, fit_poly)) / area(val_poly);
            %%plot fit thing
            %figure; hold on;
            %imagesc(cur_kid);
            %plot(val_model - middle, 'r');
            %plot(fit_shape, 'k');
            %colormap(interpolate_colors(get_colomap(2), 255));
            %set(gca, 'YLim', [min(imag(fit_shape))*0.8, max(imag(fit_shape)) * 1.2]);
            %set(gca, 'XLim', [min(real(fit_shape))*0.8, max(real(fit_shape)) * 1.2]);
            %hold off;
            
        end
        plot(initial_b, '-', 'Color', [49,130,189]/255, 'LineWidth', 2);
        set(gca, 'YLim', [min(imag(initial_b))*0.9, max(imag(initial_b)) * 1.1]);
        set(gca, 'XLim', [min(real(initial_b))*0.9, max(real(initial_b)) * 1.1]);
        axis off;
        hold off;
        saveas(fig,"segmentation_i" + i + "_k"+cur_kid_index +".png");
        close(fig);
    end
end

function [f_shape, initial_b] = shape_fitter(image2segment, tresh_hold, base_shape, basis)
    model_size = size(base_shape) / 2;
    [x, y, center] = find_initial_boundry(image2segment, tresh_hold);
    area = polyarea(x, y);
    [p_axis, center] = find_parameters(x, y, center);
    x = x - center(1);
    y = y - center(2);
    xy_c = x + 1i*y;
    initial_b = xy_c + center(1) + 1i*center(2);
    xy_steps = abs(xy_c(2:end)' - xy_c(1:end-1)')' * ...
        sum(abs(xy_c(2:end)' - xy_c(1:end-1)'), 2).^-1;
    xy_resampled = interp1(cumsum([0;xy_steps]), xy_c, ...
        linspace(0, 1, 1000), 'spline');

    ms_rot_c = rotate_shape(p_axis, area, base_shape);
    indices = zeros(size(ms_rot_c, 1), 1);
    deltas = zeros(size(ms_rot_c, 1), 1);
    %temp = zeros(size(ms_rot_c, 1), 1);
    for ind = 1:size(ms_rot_c, 1)
        [m, ~] = min(abs(xy_resampled - ms_rot_c(ind)));
        lm = ms_rot_c(ind);
        lm_m1 = ms_rot_c(mod(ind-2, size(ms_rot_c, 1)) + 1);
        lm_p1 = ms_rot_c(mod(ind, size(ms_rot_c, 1)) + 1);
        %ort = 2*lm - lm_m1 - lm_p1;
        dp = (lm - lm_p1);
        dm = (lm - lm_m1);
        ort =  dp / abs(dp) + dm / abs(dm);
        [m1, li1] = min(abs(xy_resampled - ms_rot_c(ind) - (ort/abs(ort)) * m));
        [m2, li2] = min(abs(xy_resampled - ms_rot_c(ind) + (ort/abs(ort)) * m));
        if m1 < m2
            li = li1;
            %temp(ind) = ms_rot_c(ind) + (ort/abs(ort)) * m;
        else
            li = li2;       
            %temp(ind) = ms_rot_c(ind) - (ort/abs(ort)) * m;
        end
        %+1 k mod(k, 5)+1]
        indices(ind) = li;
        deltas(ind) = xy_resampled(li) - lm;
    end
    %figure; hold on;
    %plot(xy_resampled)
    xy_resampled = permute(xy_resampled(indices), [2 1]);
    %plot(temp, 'ok');
    %plot(xy_resampled, 'x');
    %plot(ms_rot_c, 'or');
    %hold off;

    b = zeros(size(basis, 2), 1);
    f_shape = base_shape;
    for iteration_number = 1:1000
        [~,rotated_f_shape] = procrustes([real(xy_resampled), imag(xy_resampled)], ...
            [f_shape(1:model_size), f_shape(model_size+1:end)]);
        f_shape = rotated_f_shape(:,1) + 1i*rotated_f_shape(:,2);

        deltas_r = xy_resampled - f_shape;
        db = basis'*[real(deltas_r); imag(deltas_r)];
        if (norm(db)< 0.0001)
            break
        end
        b = b + db;
        f_shape = base_shape + basis*b;
    end

    f_shape = f_shape + center(1) + 1i*center(2);
end

function rotated_shape = rotate_shape(axis, area, shape)
    model_size = size(shape, 1) / 2;
    ms_x = shape(1:model_size,:);
    ms_y = shape(model_size+1:end,:);
    [ms_axis, ms_center] = find_parameters(ms_x, ms_y, nan);

    ms_x = ms_x - ms_center(1);
    ms_y = ms_y - ms_center(2);
    ms_area = polyarea([ms_x; ms_x(end)], [ms_y; ms_y(end)]);
    
    scale = sqrt((area / ms_area));
    cos_theta = dot(axis, ms_axis);
    sin_theta = sqrt(1-cos_theta^2);
    R = [cos_theta, -sin_theta; sin_theta, cos_theta];
    mr_rot = (scale*R)*[ms_x, ms_y]';

    rotated_shape = mr_rot(1,:)' + 1i*mr_rot(2,:)';
    
end

function [li, ri, middle] = split_image(image, tresh_hold)
    ts_image = zeros(size(image));
    ts_image(image > tresh_hold) = 1;
    c_rows = sum(ts_image);
    first_c_row = find(c_rows, 1);
    last_c_row = find(~c_rows(first_c_row:end), 1);
    second_c_row = find(c_rows(first_c_row+last_c_row:end), 1);
    
    middle = last_c_row + first_c_row + floor(second_c_row/2);

    ri = image(:,middle:end);
    li = image(:,1:middle);
end

function [x, y, center] = find_initial_boundry(image_to_split, tresh_hold)
    ts_image = zeros(size(image_to_split));
    ts_image(image_to_split > tresh_hold) = 1;
    [y, x] = ind2sub(size(ts_image), find(ts_image > 0));
    center = [mean(x), mean(y)];
    k = boundary(x, y);
    x = x(k);
    y = y(k);
end

function [axis, center] = find_parameters(x, y, center)
    if isnan(center)
        center = [mean(x), mean(y)];
    end
    y_hat = y - center(2);
    x_hat = x - center(1);
    I = [sum(x_hat .^ 2, 'all'), sum(x_hat .* y_hat, 'all');...
        sum(x_hat .* y_hat, 'all'), sum(y_hat .^ 2, 'all')];
    [v, ~] = eig(I);
    axis = v(1,:);
end

function colormap = get_colomap(n)
    switch n
        case 1
            colormap = [229,245,224;
                161,217,155;
                49,163,84];
        case 2
            colormap = [254,230,206;
                253,174,107;
                230,85,13];
        case 3
            colormap = [189,215,231;
                107,174,214;
                49,130,189;
                8,81,156];
        otherwise
            colormap = [166,97,26;
                223,194,125;
                245,245,245;
                128,205,193;
                1,133,113];
    end
end

function colors = interpolate_colors(colormap, n)

    noc = size(colormap, 1);
    colors = [interp1(1:noc, colormap(:,1), 1:((noc-1)/(n-1)):noc)' ...
          interp1(1:noc, colormap(:,2), 1:((noc-1)/(n-1)):noc)' ...
          interp1(1:noc, colormap(:,3), 1:((noc-1)/(n-1)):noc)'];
    colors = (1/255)*colors;
end