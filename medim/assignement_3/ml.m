
m1 = 0;
m2 = 1;
samples_per_class = 10;
test_samples_per_class = 1000;
c = 0;


sd = 5;
x = linspace(min(m1, m2) - sd*3, max(m1, m2) + sd*3, 1000);
y1s = zeros(numel(x), numel(0.01:0.01:3));
%fxs = zeros(numel(x), numel(0.01:0.01:3));
fxs = zeros(numel(x), numel(0.01:0.01:3));
%for sd = [0.5 5]
    training_data = [sd*randn(samples_per_class, 1) + m1; ...
                     sd*randn(samples_per_class, 1)+m2];
    training_classes = [-ones(samples_per_class, 1); ones(samples_per_class, 1)];
    test_data_1 = sd*randn(test_samples_per_class, 1) + m1;
    test_data_2 = sd*randn(test_samples_per_class, 1) + m2;
    i = 0;
    windows = 0.01:0.01:3;
  %  for window_size = windows
        i = i+1;
window_size = 0.01;
        c = c+1;
        
        y = parzen(training_data, x, window_size);
        y1 = parzen(training_data(1:samples_per_class), x, window_size);
        y1s(:,i) = y1;
        y2 = parzen(training_data(samples_per_class+1:end), x, window_size);
        fx = regr(training_classes, training_data, x, window_size);
        fxs(:,i) = fx;
        %[~,i] = min(abs(fx-(m1+m2)/2));
        %{
        figure;
        hold on;
        scatter(training_data, zeros(size(training_data)), 'x', 'LineWidth', 2);
        plot(x, min(y, 1),'m', 'LineWidth', 1);
        plot(x, min(y1, 1),'r', 'LineWidth', 1);
        plot(x, min(y2, 1),'b', 'LineWidth', 1);
        plot(x, normdens(sd, x, m1), 'r:')
        plot(x, normdens(sd, x, m2), 'b:')
        plot(x, 0.5*(normdens(sd, x, m1)+normdens(sd, x, m2)), 'm:')
        %plot(x(i), fx(i), 'ko', 'LineWidth', 2)
        %plot(x, y-normdens(sd, x, m2))
        hold off;
        figure;
        hold on;
        plot(x, fx,'-', 'LineWidth', 2)
        hold off;
        %}
                plot(x, fx,'--', 'LineWidth', 2)

        [~, xindices1] = min(abs(test_data_1' - x'));
        [~, xindices2] = min(abs(test_data_2' - x'));
        %window_size
        e1 = mean(fx(xindices1) > 0)
        e2 = mean(fx(xindices2) < 0)
        
   % end
%{
window_size = 1;
sd = 0.5;
figure;
hold on;
set(gca, 'ColorOrder', summer(numel(windows)));
colormap(summer(numel(windows)));
colorbar('Ticks',[0 1], ...
        'TickLabels',{'0.01','3'},...
        'FontSize', 15)
plot(x, y1s)
scatter(training_data(1:samples_per_class), ...
    ones(size(training_data(1:samples_per_class))),...
    'kx', 'LineWidth', 2)
hold off;

figure;
hold on;
set(gca, 'ColorOrder', summer(numel(windows)));
colormap(summer(numel(windows)));
colorbar('Ticks',[0 1], ...
        'TickLabels',{'0.01','3'},...
        'FontSize', 15)
plot(x, fxs)
hold off;
%}

function fx = regr(training_classes, d, x, r)
    nd = normdens(r, x, d);
    fx = sum(training_classes.*nd)./sum(nd);
end

function y = parzen(d, x, r)
    y = mean(normdens(r, x, d));
end

function p = normdens(sigma, x, mu)
    p = (2*pi*sigma^2)^(-.5) * exp((-1/(2*sigma^2))*(x-mu).^2);
end
