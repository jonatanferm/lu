alpha = [pi/3];
theta_1 = linspace(0, pi/2, 1000);
n1 = 1;
n2 = 1.5;

theta_2 = asin( (n1./n2)*sin(theta_1) );
theta_3 = alpha-theta_2;
theta_4 = asin( (n2./n1)*sin(theta_3) );
theta_4(imag(theta_4) ~= 0) = nan;

delta = theta_1 - theta_2 - theta_3 + theta_4;
delta(imag(delta) ~= 0) = nan;

figure(1)
plot(theta_1, delta', 'k');
xlim([theta_1(1) theta_1(end)]);
xlabel('Infallsvinkel')
ylabel('Avlänkningsvinkel')
%%
lambdas = linspace(380, 780, 1000);
colors = squeeze(spectrumRGB(lambdas));

mA = [1.517, 10.72;
    1.653, 10.27];
A = mA(1,:);
A(2) = 100;
n2 = A(1) + A(2)./lambdas;
theta_2 = asin( (n1./n2)'.*sin(theta_1) );
theta_3 = alpha-theta_2;
theta_4 = asin( (n2./n1)'.*sin(theta_3) );
theta_4(imag(theta_4) ~= 0) = nan;
delta = theta_1 - theta_2 - theta_3 + theta_4;
delta = delta';
delta(imag(delta) ~= 0) = nan;

figure(2)
set(groot,'defaultAxesColorOrder',colors)
plot(theta_1, delta);
xlim([theta_1(1) theta_1(end)]);
xlabel('Infallsvinkel')
ylabel('Avlänkningsvinkel')
%%
delta = delta - delta(:,end);
plot(theta_1, delta);