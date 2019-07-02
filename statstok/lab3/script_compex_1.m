%% Stationary stochastic processes, computer exercise 1
%% 2.1
close all
clear all
clc

% Load data
load data1.mat

% Plot
figure
plot(data1.x)
xlabel('Measurements')

% Estimated expected value
m = mean(data1.x)

% Confidence interval
lower = m - 1.96*std(data1.x)/sqrt(100);
upper = m + 1.96*std(data1.x)/sqrt(100);
conf_int = [lower upper]




%% 2.2
close all
clear all
clc

% Load and decide max time shift
load covProc
n = 3;

% Plot y(t) against y(t-k)
for k = 1:n
    figure
    plot(covProc(1:end-k), covProc(k+1:end),'.')
    title(['Plot for k = ' num2str(k)])
    xlabel('y(t)')
    ylabel(['y(t-' num2str(k) ')'])
end

% Look at covariance function
r = covf(covProc,n+1)

%% 2.3
close all force
clear all
clc

% Create data
simuleraenkelsumma

% Open interface
%
spekgui

%% 3.1
close all force
clear all
clc

% Load data
load cello.mat
load trombone.mat

% Open interface
spekgui

%% 3.2
close all force
clear all
clc

% Load data
load cello.mat
load trombone.mat

% Down-sample cello
n = 2;
cello2.x = cello.x(1:n:end);
cello2.dt = cello.dt*n;

% Down-sample trombone
k = 4;
tromboneD.x = trombone.x(1:k:end);
tromboneD.dt = trombone.dt*k;

% Low-pass filter then down-sample
cello2_dec.x = decimate(cello.x,2);
cello2_dec.dt = cello.dt*2;

% Open interface
spekgui


