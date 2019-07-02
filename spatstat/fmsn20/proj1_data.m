clear
%%
%load data
load swissRainfall.mat
%%
% extract covariates and reshape to images to columns
swissGrid = [swissElevation(:) swissX(:) swissY(:)];
%remove points outside of Switzerland
swissGrid = swissGrid( ~isnan(swissGrid(:,1)),:);
%%
%plot observations
figure(1)
subplot(221)
scatter(swissRain(:,3), swissRain(:,4), 20, swissRain(:,1), 'filled')
hold on
plot(swissBorder(:,1), swissBorder(:,2),'k')
axis xy tight; hold off; colorbar
title('Precip. (mm)')
%plot elevation (prediction surface + at observations sites)
subplot(222)
imagesc([0 max(swissX(:))], [0 max(swissY(:))], swissElevation, ...
        'alphadata', ~isnan(swissElevation))
hold on
plot(swissBorder(:,1), swissBorder(:,2),'k')
scatter(swissRain(:,3), swissRain(:,4), 20, swissRain(:,2), ...
        'filled','markeredgecolor','r')
%replace markeredgecolor with markeredge on older version of matlab.
axis xy tight; hold off; colorbar
title('Elevation (km)')
%plot X and Y coordinates
subplot(223)
imagesc([0 max(swissX(:))], [0 max(swissY(:))], swissX, ...
        'alphadata', ~isnan(swissX))
hold on
plot(swissBorder(:,1), swissBorder(:,2),'k')
scatter(swissRain(:,3), swissRain(:,4), 20, swissRain(:,3), ...
        'filled','markeredgecolor','r')
axis xy tight; hold off; colorbar
title('X-dist (km)')
subplot(224)
imagesc([0 max(swissX(:))], [0 max(swissY(:))], swissY, ...
        'alphadata', ~isnan(swissX))
hold on
plot(swissBorder(:,1), swissBorder(:,2),'k')
scatter(swissRain(:,3), swissRain(:,4), 20, swissRain(:,4), ...
        'filled','markeredgecolor','r')
axis xy tight; hold off; colorbar
title('Y-dist (km)')
%%
%Extract observations and split data into obsrvation and validation set
Y = swissRain(swissRain(:,5)==0,:);
Yvalid = swissRain(swissRain(:,5)==1,:);

%transform data
Y(:,1) = sqrt(Y(:,1));
Yvalid(:,1) = sqrt(Yvalid(:,1));

%perform a linear regression on elevation
%X = [ones(size(Y,1),1)];
%X = [ones(size(Y,1),1), Y(:,2)];
%X = [ones(size(Y, 1), 1), Y(:,2), Y(:,2).^2 ];
%X = [ones(size(Y, 1), 1), Y(:,2), Y(:,2).^2, Y(:,2).^3];
%X = [ones(size(Y,1),1), Y(:,3) * 0.01, Y(:,4) * 0.01];
X = [ones(size(Y,1),1), Y(:,2), Y(:,3) * 0.01];
%X = [ones(size(Y,1),1), Y(:,2), Y(:,2).^2, Y(:,3) * 0.01];
%X = [ones(size(Y,1),1), Y(:,2), Y(:,2).^2, Y(:,2).^3, Y(:,3) * 0.01];
%X = [ones(size(Y,1),1), Y(:,2), Y(:,3) * 0.01, Y(:,4) * 0.01];
%X = [ones(size(Y, 1), 1), Y(:,2), Y(:,2).^2, Y(:,3) * 0.01, Y(:,4) * 0.01];
%X_v = [ones(size(Yvalid, 1), 1)];
%X_v = [ones(size(Yvalid, 1), 1), Yvalid(:,2)];
%X_v = [ones(size(Yvalid, 1), 1), Yvalid(:,2),Yvalid(:,2).^2];
%X_v = [ones(size(Yvalid, 1), 1), Yvalid(:,2),Yvalid(:,2).^2,Yvalid(:,2).^3 ];
%X_v = [ones(size(Yvalid, 1), 1), Yvalid(:,3) * 0.01, Yvalid(:,4) * 0.01];
X_v = [ones(size(Yvalid, 1), 1), Yvalid(:,2), Yvalid(:,3) * 0.01];
%X_v = [ones(size(Yvalid, 1), 1), Yvalid(:,2),Yvalid(:,2).^2, Yvalid(:,3) * 0.01];
%X_v = [ones(size(Yvalid, 1), 1), Yvalid(:,2), Yvalid(:,2).^2, Yvalid(:,2).^3, Yvalid(:,3) * 0.01];
%X_v = [ones(size(Yvalid, 1), 1), Yvalid(:,2), Yvalid(:,3) * 0.01, Yvalid(:,4) * 0.01];
%X_v = [ones(size(Yvalid, 1), 1), Yvalid(:,2), Yvalid(:,2).^2, Yvalid(:,3) * 0.01, Yvalid(:,4) * 0.01];
beta = X\Y(:,1);

%extract covariates (elevation) for the prediction locations
%Xgrid = [ones(size(swissGrid,1),1)];
%Xgrid = [ones(size(swissGrid,1),1), swissGrid(:,1)];
%Xgrid = [ones(size(swissGrid,1),1), swissGrid(:,1), swissGrid(:,1).^2];
%Xgrid = [ones(size(swissGrid,1),1), swissGrid(:,1), swissGrid(:,1).^2, swissGrid(:,1).^3];
%Xgrid = [ones(size(swissGrid,1),1), swissGrid(:,2) * 0.01, swissGrid(:,3) * 0.01];
Xgrid = [ones(size(swissGrid,1),1), swissGrid(:,1), swissGrid(:,2) * 0.01];
%Xgrid = [ones(size(swissGrid,1),1), swissGrid(:,1),swissGrid(:,1).^2, swissGrid(:,2) * 0.01];
%Xgrid = [ones(size(swissGrid,1),1), swissGrid(:,1), swissGrid(:,1).^2, swissGrid(:,1).^3 swissGrid(:,2) * 0.01];
%Xgrid = [ones(size(swissGrid,1),1), swissGrid(:,1), swissGrid(:,2) * 0.01, swissGrid(:,3) * 0.01];
%Xgrid = [ones(size(swissGrid,1),1) swissGrid(:,1), swissGrid(:,1).^2, swissGrid(:,2) * 0.01, swissGrid(:,3) * 0.01];

%recall that swissGrid contains all relevant locations (i.e. not nan values)
%we need to recreate the full matrix
mu = nan( size(swissElevation) );
%and place the predictions at the correct locations
mu_pred = Xgrid*beta;
validation_pred = X_v*beta;
mu_pred(mu_pred<0) = 0;
%mu_pred = mu_pred .^ 2;
mu( ~isnan(swissElevation) ) = mu_pred;
%%
%and plot the result
figure(2)
imagesc([0 max(swissX(:))], [0 max(swissY(:))], mu, ...
        'alphadata', ~isnan(mu))
hold on
plot(swissBorder(:,1), swissBorder(:,2),'k')
scatter(Y(:,3), Y(:,4), 200, Y(:,1), 'filled','markeredgecolor','r')
axis xy tight; hold off; colorbar
title('sqrt of rainfall and predictions')

validation_error = abs((Yvalid(:,1) - validation_pred)) ./ Yvalid(:,1);
mean_ve = mean(validation_error)
max_ve = max(validation_error);
beta;


%sigma_e_valid = sum((Yvalid(:,1) - X_v*beta) .^ 2) / (size(X, 1) - size(X, 2));
%V_beta = sigma_e*inv(X'*X);
sigma_e = sum((Y(:,1) - X*beta) .^ 2) / (size(X, 1) - size(X, 2));
V_yhat_rec = sigma_e * (1 + sum((Xgrid / (X'*X)) .* Xgrid,2));
V_yhat_valid = sigma_e * (1 + sum((X_v / (X'*X)) .* X_v,2));

V_yhat_rec_2d = nan( size(swissElevation) );
V_yhat_rec_2d( ~isnan(swissElevation) ) = V_yhat_rec;

figure(3)
imagesc([0 max(swissX(:))], [0 max(swissY(:))], V_yhat_rec_2d, ...
        'alphadata', ~isnan(V_yhat_rec_2d))
hold on
plot(swissBorder(:,1), swissBorder(:,2),'k')
%scatter(Y(:,3), Y(:,4), 20, Y(:,1), 'filled','markeredgecolor','r')
axis xy tight; hold off; colorbar
title('Variance')
mean_vv = mean(V_yhat_valid)

%%
