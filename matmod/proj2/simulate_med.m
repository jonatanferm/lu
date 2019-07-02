function [x_hist, time] = simulate_med(med_const, human_const, doses, end_time, dose_time)
% dose_time in hours
% end_time in hours

CL = med_const(1); % clearing
a = med_const(2); % absorption
k_bo = med_const(3); % blood to other (fat)
k_ob = med_const(4); % other (fat) to blood
m = human_const(1); % human mass
V_b = human_const(2); % human blood volume 
init_dose = [doses(1) 0 0 0]'; % inital dose
cont_dose = [doses(2) 0 0 0]'; % continuous dose.

x = init_dose;
x_hist = x;

K = [1-a 0 0 0; a/V_b 1-CL-k_bo k_ob/V_b 0; 0 k_bo*V_b 1-k_ob 0; 0 V_b*CL 0 1];

time = 0:1/60:end_time;
for i = 1/60:1/60:end_time
    if(rem(i,dose_time) == 0)
        x = K*x + cont_dose;
    else
        x = K*x;
    end
    x_hist = [x_hist x];
end