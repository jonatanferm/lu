clear;
load('./data/dmsa_images');
load('./data/man_seg');
load('./data/models');

%% Resample reference image
%val_models = models;
%models = models(:, 6:40);
%models = models(:,[1:5 11:40]);
%models = models(:,[1:10 16:40]);
%models = models(:,[1:15 21:40]);
closed_models = [models; models(1,:)];
closed_manual = [man_seg; man_seg(1,:)];
avg_model_step = abs(closed_models(2:end,:)' - closed_models(1:end-1,:)')' * sum(abs(closed_models(2:end,:)' - closed_models(1:end-1,:)'), 2).^-1 / size(models, 2);
man_steps = abs(closed_manual(2:end)' - closed_manual(1:end-1)')' * sum(abs(closed_manual(2:end)' - closed_manual(1:end-1)'), 2).^-1;
resampled = interp1(cumsum([0;man_steps]), [man_seg; man_seg(1,:)], cumsum([0;avg_model_step]), 'spline');
resampled = resampled(1:end-1,:);
%% Plot resampling
%{
fig = figure('Visible','Off');
hold on;
plot(resampled, 'xr');
pbaspect([1 1 1])
plot(man_seg, '.k');
legend({'Resampled points','Original points'},'Location','southwest')
hold off;
saveas(fig,"resampledpoints.png")
%}
%% Align images
model_size = size(models, 1);
y = [real(resampled(:,1)), imag(resampled(:,1))];
x = num2cell(reshape([real(models); imag(models)], model_size, 2, []), [1 2]);
aligner = get_aligner(y);
tform = cellfun(aligner, x, 'UniformOutput', false);
x_aligned = cell2mat(...
    cellfun(@transformPointsForward, tform, x, 'UniformOutput', false));
x_aligned = [squeeze(x_aligned(:,1,:)); squeeze(x_aligned(:,2,:))]; %Make row
%x_aligned = squeeze(x_aligned(:,1,:) + 1i*x_aligned(:,2,:));
x_aligned_c = x_aligned(1:model_size,:) + 1i*x_aligned(model_size+1:end,:);
mean_shape = squeeze(sum(x_aligned,2)) / size(x_aligned, 2);
complex_mean_shape = mean_shape(1:model_size) + 1i*mean_shape(model_size+1:end); %make mean complex form
%% Plot aligned images
%{
fig = figure;%('Visible', 'Off');
hold on;
pbaspect([1 1 1]);
colors = (1/255)*[216,179,101;90,180,172];
set(gca, 'ColorOrder', colors);
set(gca,'Color',(1/255)*[245,245,245])
plot(permute(x_aligned_c, [2 1]), '.'); %plot locations of landmarks
mean = plot(complex_mean_shape, 'kx');
legend(mean, {'Mean image landmarks'},'Location','southwest');
saveas(fig,"alignedlandmarks.png")
hold off;
%}
%% Find a basis
x_d = x_aligned - mean_shape;
%[U, v] = eig(x_d*x_d');
[full_basis, v] = eig(cov(x_d'));
modes = 10;%sum(cumsum(flip(diag(v))) / sum(diag(v)) < 0.95) + 1;
basis = full_basis(:,end-modes+1:end);
%% plot mode powers
%{
fig = figure('Visible','Off');
hold on
title("Cumlative power in modes.", 'FontSize', 18);
xlabel("Number of modes.", 'FontSize', 18) 
ylabel("Part of total power.", 'FontSize', 18)
axis([0 15 0 1])
stem((1:size(v, 1))',cumsum(flip(diag(v)) / sum(diag(v))), '-')
refline([0, .95]);
saveas(fig,"modepowers.png")
hold off
close(fig)
%}
%% Vizualize a mode.
%{
max_val = 3;
step = 0.01;
row_steps = -max_val:step:max_val;
nbr_of_rows = size(row_steps, 2);
colors = [166,97,26;223,194,125;245,245,245;128,205,193;1,133,113];
colors = [interp1(1:size(colors, 1), colors(:,1), 1:(4/(nbr_of_rows-1)):size(colors, 1))' ...
          interp1(1:size(colors, 1), colors(:,2), 1:(4/(nbr_of_rows-1)):size(colors, 1))' ...
          interp1(1:size(colors, 1), colors(:,3), 1:(4/(nbr_of_rows-1)):size(colors, 1))'];
colors = (1/255)*colors;
closed_complex_mean = [complex_mean_shape; complex_mean_shape(1,:)]; %Make mean closed figure
sq_lambdas = sqrt(v(end-modes+1:end,end-modes+1:end));
for mode = 1:modes
    b = [zeros(nbr_of_rows,modes-mode), ...
        row_steps', ...
        zeros(nbr_of_rows, mode-1)];
    b = b*sq_lambdas;
    rec_image = mean_shape + basis * b';
    rec_image_temp = rec_image(1:model_size, :) + 1i*rec_image(model_size+1:end, :);%Make complex form
    rec_image = [rec_image_temp; rec_image_temp(1,:)]; %Make closed figure
    fig = figure('Visible','Off');
    hold on;
    axis off;
    colorbar('Ticks',(-max_val:1:max_val) * sq_lambdas(modes-mode+1, modes-mode+1), ...
        'TickLabels',{'-3\lambda^{0.5}','-2\lambda^{0.5}','-\lambda^{0.5}','0','\lambda^{0.5}', '2\lambda^{0.5}', '3\lambda^{0.5}'},...
        'FontSize', 15);
    pbaspect([1 1 1])
    set(gca, 'ColorOrder', colors);
    colormap(colors);
    caxis([-max_val, max_val] * sq_lambdas(modes-mode+1, modes-mode+1));
    plot(rec_image);
    plot([rec_image(:,1) rec_image(:,end)], 'k');
    plot(closed_complex_mean, 'k', 'LineWidth', 3)
    plot(permute([rec_image(:,1) rec_image(:,end)], [2 1]), 'k','LineWidth', 1.5);
    title("Variantions of mode " + num2str(mode) + ". \lambda^{0.5} = " + num2str(sq_lambdas(modes-mode+1, modes-mode+1), 2), 'FontSize', 15)
    hold off
    saveas(fig,"mode"+mode+".png")
    close(fig)
end
%}
%%
function colors = get_colormap(n)
    colors_b = [166,97,26;223,194,125;245,245,245;128,205,193;1,133,113];
    colors = [interp1(1:size(colors_b, 1), colors_b(:,1), 1:(4/(n-1)):size(colors_b, 1))' ...
          interp1(1:size(colors_b, 1), colors_b(:,2), 1:(4/(n-1)):size(colors_b, 1))' ...
          interp1(1:size(colors_b, 1), colors_b(:,3), 1:(4/(n-1)):size(colors_b, 1))'];
    colors = (1/255)*colors;
end

function f = get_aligner(y)
    function tform = get_trans(x)
        tform = fitgeotrans(x,y,'NonreflectiveSimilarity');
    end
    f = @get_trans;
end