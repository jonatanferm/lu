%% Divergens
a = [150 200 250]; % avstånd laser -> lins (mm)
b = [55 76 97]; % stålbredd (mm)

divs = atan(b./(2*a)) * (360/(2*pi)) %divergens vinklar
mean(divs);

%% Vinklar för prismat

x = 350; %x avstånd från prisma till mätpunkt
y = 255; %y avstånd
thetas2 = [40 42 44 46 48 50 52 54 56 58];
thetas4 = [40 44 48 52 56 60 64 68 72 74]; %Infallsvinklar
thetasm = [40 42 44 46 48 50 52 54 56 58 60 62 64 66 68 70];
deltasm = [59.4594 56.4424 54.1205 52.7652 51.8091 51.1080 50.8774 50.7626 50.8774 51.1080 51.4568 51.9272 52.6443 53.3753 54.1205 54.8800];
neg_y_deltas = [17 82 115 125 126 120 107 86 57 16]; %Avsånt på pappret
deltas = atan((y+ 250 - neg_y_deltas)./x)*360./(pi*2);

n2 = fminsearch(@(x) delta_squared_diff(thetas2, deltas, 60, 1, x(1)), [1]);
n4 = fminsearch(@(x) delta_squared_diff(thetas4, deltas, 60, 1, x(1)), [1]);
nm = fminsearch(@(x) delta_squared_diff(thetasm, deltasm, 60, 1, x(1)), [1]);
t_thetas = 40:.1:74;
figure(1);
hold on;
plot(thetas2, atan((y+ 250 - neg_y_deltas)./x)*360./(pi*2), 'ob')
plot(thetasm, deltasm, 'or')
plot(t_thetas, find_delta(t_thetas, 60, 1, n2) ,'-b', 'LineWidth', 1.2);
plot(t_thetas, find_delta(t_thetas, 60, 1, nm) ,'-r', 'LineWidth', 1.2);
legend('Mätvärden 641nm', 'Mätvärden 404nm', 'Anpassad kurva 641nm', 'Anpassad kurva 404nm')
xlabel('Infallsvinklar', 'FontSize', 14)
ylabel('Avlänkningsvinklar', 'FontSize', 14)

figure(2);
hold on;
plot(thetas4, atan((y+ 250 - neg_y_deltas)./x)*360./(pi*2), 'ob')
plot(thetasm, deltasm, 'or')
plot(t_thetas, find_delta(t_thetas, 60, 1, n4) ,'-b', 'LineWidth', 1.2);
plot(t_thetas, find_delta(t_thetas, 60, 1, nm) ,'-r', 'LineWidth', 1.2);
legend('Mätvärden 641nm', 'Mätvärden 404nm', 'Anpassad kurva 641nm', 'Anpassad kurva 404nm')
xlabel('Infallsvinklar', 'FontSize', 14)
ylabel('Avlänkningsvinklar', 'FontSize', 14)
%%

fminsearch(@(x) delta_squared_diff(thetas, deltas, 60, 1, x(1)), [1])

function t_delta = find_delta(theta1,alpha, n1, n2)
    theta2 = asind( (n1./n2) .* sind(theta1) );
    theta3 = alpha - theta2;
    theta4 = asind( (n2./n1) .* sind(theta3) );
    t_delta = theta1 - theta2 - theta3 + theta4;
end
function sd = delta_squared_diff(theta1, delta,alpha, n1, n2)
    sd = sum((delta - find_delta(theta1,alpha, n1, n2)).^2);
end