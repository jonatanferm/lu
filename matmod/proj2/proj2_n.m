clear all
data = xlsread('pk',1);
con = data(:,3); % blood concentration
side = data(:,4); % side effect
%%
% clf
% hold on
% p = zeros(10,10);
p = zeros(10,10);
measure_times = [0.5 1 2 4 8 12 18 24 48 96];
interpolated_time = [0.5:0.5:96];
interpolated_data = zeros(size(p,1),size(interpolated_time,2));
for i = 1:10
    p1 = data((i-1)*10+1:(i-1)*10+10,3);
    p(i,:) = p1;
    interpolated_data(i,:) = interp1(measure_times, p1, interpolated_time);
%     plot(measure_times,p1,'x-')
end
p_mean = mean(p,1);
p_median = median(p,1);
%%
med_const = [0.005 0.004 0.007 0.004];
interpolated_median = interp1(measure_times, p_median, interpolated_time);
interpolated_data2 = reshape(interpolated_data', [1 size(interpolated_data,1)*size(interpolated_data,2)]);
med_const = fminsearch('Qfunc2', med_const, [], interpolated_time, interpolated_median)
%% Simulate medication
% med_const = [0.005 0.004 0.007 0.004]; 
% med_const = test;
m = 20;
human_const = [m 0.07*m];
doses = [15 0];
[x, t] = simulate_med(med_const, human_const, doses, 300, 0);
doses = [2 2];
[x_5, t_5] = simulate_med(med_const, human_const, doses, 300, 12);
% Plot simulation
clf
hold on
plot([0 96], [1 1], 'r-.', 'linewidth',1) % 
plot([0 96], [2 2], 'b-.', 'linewidth',1) % 1
plot([0 96], [3 3], 'b-.', 'linewidth',1) % 2
plot([0 96], [5 5], 'b-.', 'linewidth',1) % 3
plot(t,x(2,:),'k-', 'LineWidth',2)
plot(t_5,x_5(2,:),'b-')
%plot(measure_times, p_mean,'g.-', 'linewidth',2)
%plot(measure_times, p_median,'m.-', 'linewidth',2)

for i = 1:size(p,1)
   plot(measure_times,p(i,:))
end
%%
clf
p_median