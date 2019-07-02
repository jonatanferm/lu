%% 1. Load mnist data
[train_im,train_classes,train_angles]=digitTrain4DArrayData;
[test_im,test_classes,test_angles]=digitTest4DArrayData;
%% 2. Select deep learning architecture
%layers = [
%    imageInputLayer([28 28 1]) % Specify input sizes
%    fullyConnectedLayer(10)    % Fully connected is a affine map from 28^2 pixels to 10 numbers
%    softmaxLayer               % Convert to 'probabilities'
%    classificationLayer];      % Specify output layer
drop_prob = 0.5;
layers = [
    imageInputLayer([28 28 1])
    convolution2dLayer(7,36, 'Stride', 2, 'Padding', 'same')
    batchNormalizationLayer
    dropoutLayer(drop_prob)
    reluLayer
    maxPooling2dLayer(2,'Stride',2)
    convolution2dLayer(4,2*36, 'Stride', 1, 'Padding', 'same')
    batchNormalizationLayer
    dropoutLayer(drop_prob)
    reluLayer
    maxPooling2dLayer(2,'Stride',2)
    fullyConnectedLayer(100)
    dropoutLayer(drop_prob)
    fullyConnectedLayer(50)
    dropoutLayer(drop_prob)
    fullyConnectedLayer(50)
    fullyConnectedLayer(10)
    softmaxLayer
    classificationLayer];
%% 3. Train deep learning network
miniBatchSize = 512;       
max_epochs = 300;           % Specify how long we should optimize
learning_rate = 0.0003;     % Try different learning rates 
options = trainingOptions( 'adam',...
    'MaxEpochs',max_epochs,...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.8, ...
    'LearnRateDropPeriod',25, ...
    'InitialLearnRate',learning_rate, ...
    'ValidationData',{test_im,test_classes}, ...
    'ValidationFrequency',30, ...
    'Plots', 'training-progress');
net = trainNetwork(train_im, train_classes, layers, options);
%% 4. Test the classifier on the test set
[Y_result2,scores2] = classify(net,test_im);
accuracy2 = sum(Y_result2 == test_classes)/numel(Y_result2);
disp(['The accuracy on the test set: ' num2str(accuracy2)]);
%% Test the classifier on the training set
[Y_result1,scores1] = classify(net,train_im);
accuracy1 = sum(Y_result1 == train_classes)/numel(Y_result1);
disp(['The accuracy on the training set: ' num2str(accuracy1)]);

