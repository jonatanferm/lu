% pop_growth.m
% Parameters and start values used in this excercise:
n0=100; % size of population at time 0
r=2; % population growth rate
t=0:0.2:10; % time period (in years)
nt=n0*exp(r.*t); % the equation for pop. growth in continuos time
figure % opens a new figure-window
plot(t,nt) %plots the figure
%% VARYING K
% intra_comp
% parameters and start values used in this excercise:
n0=1; % size of population at time 0
r=1; % population growth rate
t=0:0.2:20; % timeperiod (in years) 

k=(1:1000)'; %carrying capacity
nt=k./(1+(k./n0-1).*exp(-r.*t)); % the solution (equation 16) to the logistic equation
figure
set(groot,'defaultAxesColorOrder',summer(1000))
plot(t,nt)
colormap(summer(1000))
%axes('Clim', [0 10])
caxis([0 1])
colorbar('XtickLabel', {'K=0','K=500', 'K=1000'}, 'Xtick', [0, .5, 1])
%% VARYING r
n0=10000; % size of population at time 0
r=logspace(-.5, .5, 5000)'; % population growth rate
t=0:0.2:20; % timeperiod (in years) 

k=1000; %carrying capacity
nt=k./(1+(k./n0-1).*exp(-r.*t)); % the solution (equation 16) to the logistic equation
figure
set(groot,'defaultAxesColorOrder',summer(numel(r)))
plot(t,nt)
colormap(summer(numel(r)))
ylim([0 1e4])
%axes('Clim', [0 10])
caxis([0 1])
colorbar('XtickLabel', {'r~.3','r~1', 'r~3'}, 'Xtick', [0, .5, 1])
%% m;g
% intra_comp
% parameters and start values used in this excercise:
n0=1001; % size of population at time 0
r=-1; % population growth rate
t=0:0.2:20; % timeperiod (in years)
k=1000; %carrying capacity
nt=k./(1+(k./n0-1)*exp(-r.*t)); % the solution (equation 16) to the logistic equation
figure
plot(t,nt)
%%
% parameters
n=0:10:1000; % initial and end population size
r=1; % pop. growth rate
k=1000; %carrying capacity
dndt=r.*n.*(1-n/k); % the logistic equation
figure
plot(n,dndt)
%%
% discrete.m
% parameters
clear all;
n0=1; %start value
r=1; %rate of increase
k=1000; %carrying capacity
%loop:
n(1) = n0; %assign a value
for t = 1:20 %loop from 1 to end
%equation 16, calculates values into matrix n:
n(t+1)=n(t)+(r.*n(t).*(1-n(t)./k));
end
t=1:21; %final step 21 because of t+1
figure % figure window
hold on;
plot(t,n); % plot n

n0=1; % size of population at time 0
r=1; % population growth rate
t=0:0.2:20; % timeperiod (in years)
k=1000; %carrying capacity
nt=k./(1+(k./n0-1)*exp(-r.*t)); % the solution (equation 16) to the logistic equation
plot(t,nt)
hold off;
%%
clear all;
n0=1; %start value
r=2.1; %rate of increase
k=1000; %carrying capacity
%loop:
n(1) = n0; %assign a value
N = 1000;
for t = 1:N %loop from 1 to end
%equation 16, calculates values into matrix n:
n(t+1)=n(t)+(r.*n(t).*(1-n(t)./k));
end
t=1:N+1; %final step 21 because of t+1
figure % figure window
plot(t,n); % plot n
%% PRED
% Predation.m
% Parameters
% n1 = density of prey population
% n2 = density of predator population

clear all;
global r a b k d;
r=0.5;         	 % Growth rate of prey population
a=0.02;        	% predator kill rate
b=0.1;        	% predator conversion coefficient  
d=0.1;         	% predator mortality in abscence of prey
k=2000;        	% carrying capacity of prey population
n0=[10;10];    	% density of prey and predator at t=0
t0=0; tend=1000; 	% simulation start and stop

dat=[ ];       	% Vector for storing the result
t=[ ];		% time vector

% Calculate population densities with ODE45, a MATLAB application which 
% uses the differential equation below (pred.m), solves the equation numerically 

[t,dat]=ode45('pred',[t0 tend],n0); 	%’pred’ is the equation below
figure	% opens a new figure window
plot(t,dat);	% draws the figure
legend('prey','pred');	% makes a legend
title(sprintf('r=%.3f a=%.3f  b=%.3f  d=%.3f  k=%.0f', r, a, b, d, k))
%%
%Add to the file predation.m
figure
plot(dat(:,1),dat(:,2));
ylabel('predator');
xlabel('prey');
title(sprintf('r=%.3f a=%.3f  b=%.3f  d=%.3f  k=%.0f', r, a, b, d, k))
%%
% inter_comp.m 
% Interspecific competition 
clear all;
close all;
global r k alfa beta;
% Parameters, initial values
n0=[100;1]; 	% initial densities of the population n(1)=1 & n(2)=1
r=[1;.2];  	% growth rate of the populations r(1)=1 & r(2)=1;
k=[300;300]; 	% carrying capacity of the system k(1)=1000 och K(2)=1000 
alfa=1.1;      	% per capita effect of n(2) on n(1)
beta=1.1;   	% per capita effect of n(1) on n(2)

t0=0;tend=300;	% time step 0 - 300
dat=[];	% vector for storing the result
t=[];		% vector for storing the timesteps

% Function, equation, see separate file ’icomp.m’ 
% ode-routine
[t,dat]=ode45('icomp',[t0 tend],n0);	
% solves the equation numerically
% making the plot
figure;
subplot(2,1,1);
plot(t,dat);ylabel('Population density');
xlabel('Time');
title(sprintf('$\\alpha=%.3f, \\beta=%.3f, r_1=%.3f, r_2=%.3f, k_1=%.0f, k_2=%.0f$', ...
    alfa, beta, r(1), r(2), k(1), k(2)), 'Interpreter','latex')
	% phase diagram
[t,dat]=ode45('icomp',[t0 tend], n0);
subplot(2,1,2);
plot(dat(:,1),dat(:,2))
xlabel ('species 1')
ylabel ('species 2')

%add starting point as a star
hold on;
plot(dat(1,1),dat(1,2),'*')
%make dn1dt=0 isocline
y1=linspace(0,k(1)/alfa,10);
x1=k(1)-alfa*y1;

%make dn2dt=0 isocline
x2=linspace(0,k(2)/beta,10);
y2=k(2)-beta*x2;

%plot the isoclines
plot(x1,y1,'r');
plot(x2,y2,'g');

%make legend and set axis
legend('Population sizes','Starting point','dn1dt=0','dn2dt=0')
axis([0 max([x1 x2]) 0 max([y1 y2])]);