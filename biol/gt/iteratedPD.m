function iteratedPD
% An iterated Prisoner's dilemma game
% Pure strategies:
% C = Cooperate (coded as 1)
% D = Defect (coded as 2)
% n = number of rounds
% Payoff matrix:
% U = [R, S; T, P]
% Strategies are defined as a sequence of pure strategies, 
% e.g. [C, C, D, C, D], which is coded as the row vector [1, 1, 2, 1, 2]

T = 10;
R = 5;
P = 2;
S = 0;
U = [R S; T P];

n = 3; % number of rounds
N = 500; % population size
% Start with a population of allC, i.e. always cooperate:
pop = ones(N,n); % each row is an individual

tmax = 100; % number of generations
Pmutation = 1e-2; % Mutation probability
meanStrategies = zeros(tmax,n); % storage of the population mean
meanStrategies(1,:) = mean(pop); % First generation mean
for t=2:tmax
	% First calculate the fitness of each individual by letting it play
	% against all other individuals in the population:
	fitness = zeros(1,N);
	for i = 1:N-1
		for j = i+1:N
			[wi,wj] = payoff(pop(i,:), pop(j,:), n, U);
			fitness(i) = fitness(i) + wi;
			fitness(j) = fitness(j) + wj;
		end
	end
	% Next, create the next generation
	nextpop = zeros(N,n);
	cumsumW = cumsum(fitness);
	for i=1:N
		% Select a parent randomly, with probabilities weighted by the individuals' fitness:
		x = rand*cumsumW(end);
		parent = find(x<cumsumW,1,'first'); 
		nextpop(i,:) = pop(parent,:); % copy the selected parents strategy
		if rand<Pmutation % mutate the strategy with a small probability
			mutatedPosition = ceil(rand*n);
			nextpop(i,mutatedPosition) = 3 - nextpop(i,mutatedPosition); % This turns a 1 into a 2 and vice versa
		end
	end
	pop = nextpop;
	meanStrategies(t,:) = mean(pop); % Save the generation mean
	if mod(t,10)==0 % every 10th generation, do some plotting
		figure(1), clf
        for i=1:n
            subplot(n,1,i)
    		plot(meanStrategies(1:t,i))
            axis([0 tmax 1 2])
    		xlabel('time (generations)')
    		title(['Mean strategy, round ' num2str(i)])
        end
		drawnow
	end
end