function err = Qfunc2(constants, time, data)
m = 70;
med_const = constants;
human_const = [m 0.07*m];
doses = [15 0];
indx = time*60+1;
[x, ~] = simulate_med(med_const, human_const, doses, 96, 0);
pred = x(2,:);
% err = sum((pred(indx)-data).^2);
err = 0;
N = length(time);
for i = 1:length(data)/N
    err = err + sum((pred(indx) - data((i-1)*N + 1:(i-1)*N + N)).^2);
%     err = err + sum((pred(indx) - data((i-1)*N+1:(i-1)*N+N)).^2);
end
%if(sum(constants <= 10^-3) > 0)
%    err = 100000000000;
%end
end