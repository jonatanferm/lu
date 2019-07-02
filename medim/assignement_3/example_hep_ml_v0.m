%% Load data (training data X1 - labels Y1)
% test data X2 - labels Y2

clear all
load(fullfile('databases','hep_proper_mask'));
X1_masks = Y1;
X2_masks = Y2;
load(fullfile('databases','hep_proper'));

%% Use hand made features

nr_of_training_images = size(X1,4);
for i = 1:nr_of_training_images,
   [fv,str]=get_features(X1(:,:,1,i),X1_masks(:,:,1,i));
   X1f(i,:)=fv;
end

nr_of_test_images = size(X2,4);
for i = 1:nr_of_test_images,
   [fv,str]=get_features(X2(:,:,1,i),X2_masks(:,:,1,i));
   X2f(i,:)=fv;
end

%% Visualize the data
figure;
dims = size(str, 2);
blorp = dims*(dims-1)/2;
f = floor(sqrt(blorp));
g = ceil(blorp/f);
c = 0;
for i = 1:dims-1
    for j = i+1:dims
        %[i j]
        c = c+1;
        subplot(f,  g, c)
        gscatter(X1f(:,i), X1f(:,j), Y1 ,'rgbrgb','ooo+++');
        xlabel(str{i});
        ylabel(str{j});
    end
end

%% Train machine learning models to data

%% Fit a decision tree model
%{
model1 = fitctree(X1f(:,1:2),Y1,'PredictorNames',str);
% Test the classifier on the training set
[Y_result1,node_1] = predict(model1,X1f);
accuracy1 = sum(Y_result1 == Y1)/numel(Y_result1);
disp(['The accuracy on the training set: ' num2str(accuracy1)]);
% Test the classifier on the test set
[Y_result2,node_2] = predict(model1,X2f);
accuracy2 = sum(Y_result2 == Y2)/numel(Y_result2);
disp(['The accuracy on the test set: ' num2str(accuracy2)]);
%}

%% Fit a random forest model
model2 = TreeBagger(50,X1f,Y1,'OOBPrediction','On',...
    'Method','classification')
% Test the classifier on the training set
[Y_result1,node_1] = predict(model2,X1f);
accuracy1 = sum(Y_result1 == Y1)/numel(Y_result1);
disp(['The accuracy on the training set: ' num2str(accuracy1)]);
% Test the classifier on the test set
[Y_result2,node_2] = predict(model2,X2f);
accuracy2 = sum(Y_result2 == Y2)/numel(Y_result2);
disp(['The accuracy on the test set: ' num2str(accuracy2)]);

%% Fit a random forest model

model3 = fitcecoc(X1f,Y1);
% Test the classifier on the training set
[Y_result1,node_1] = predict(model3,X1f);
accuracy1 = sum(Y_result1 == Y1)/numel(Y_result1);
disp(['The accuracy on the training set: ' num2str(accuracy1)]);
% Test the classifier on the test set
[Y_result2,node_2] = predict(model3,X2f);
accuracy2 = sum(Y_result2 == Y2)/numel(Y_result2);
disp(['The accuracy on the test set: ' num2str(accuracy2)]);
%}
%% Fit a k-nearest neighbour model
model4 = fitcknn(X1f,Y1, 'NumNeighbors', 4);
% Test the classifier on the training set
[Y_result1,node_1] = predict(model4,X1f);
accuracy1 = sum(Y_result1 == Y1)/numel(Y_result1);
disp(['The accuracy on the training set: ' num2str(accuracy1)]);
% Test the classifier on the test set
[Y_result2,node_2] = predict(model4,X2f);
accuracy2 = sum(Y_result2 == Y2)/numel(Y_result2);
disp(['The accuracy on the test set: ' num2str(accuracy2)]);

%model5 = fitcecoc(X1f,Y1);
%[Y_result1,node_1] = predict(model5,X1f);
%accuracy1 = sum(Y_result1 == Y1)/numel(Y_result1);
%disp(['The accuracy on the training set: ' num2str(accuracy1)]);
% Test the classifier on the test set
%[Y_result2,node_2] = predict(model5,X2f);
%accuracy2 = sum(Y_result2 == Y2)/numel(Y_result2);
%disp(['The accuracy on the test set: ' num2str(accuracy2)]);



