load(fullfile('databases','hep_proper'));

shuffle = randperm(size(X1(:,:,1,:), 4));
train_im = X1(:,:,1,:);
train_classes = Y1;
train_im = train_im(:,:,1,shuffle);
train_classes = train_classes(shuffle);

shuffle = randperm(size(X2(:,:,1,:), 4));
test_im = X2(:,:,1,:);
test_classes = Y2;
test_im = test_im(:,:,1,shuffle);
test_classes = test_classes(shuffle);

imageSize = [64 64 1];
pixelRange = [-4 4];
%    'RandRotation', [0 180], ...

imageAugmenter = imageDataAugmenter( ...
    'RandXReflection',true, ...
    'RandYReflection',true, ...
    'RandXTranslation',pixelRange, ...
    'RandYTranslation',pixelRange);
augimdsTrain = augmentedImageDatastore(imageSize,train_im,train_classes, ...
    'DataAugmentation',imageAugmenter, ...
    'OutputSizeMode','randcrop');

drop_prob = 0.5;
netWidth = 64;
layers = [
    imageInputLayer([64 64 1])
    dropoutLayer(0.2)   
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
    fullyConnectedLayer(60)
    dropoutLayer(drop_prob)
    fullyConnectedLayer(30)
    dropoutLayer(drop_prob)
    fullyConnectedLayer(30)
    fullyConnectedLayer(6)
    softmaxLayer
    classificationLayer
];

% slightly deeper 

%% Train deep learning network
%lgraph = layerGraph(layers);
%lgraph = connectLayers(lgraph,'relu1','coll');
%lgraph = connectLayers(lgraph,'relu2','coll');
%lgraph = connectLayers(lgraph,'relu3','coll');

%plot(lgraph);

miniBatchSize = 128;
max_epochs = 300;
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

net = trainNetwork(augimdsTrain, layers, options);


%% Test the classifier on the test set

[Y_result2,scores2] = classify(net,test_im);
accuracy2 = sum(Y_result2 == test_classes)/numel(Y_result2);
disp(['The accuracy on the test set: ' num2str(accuracy2)]);

%% Test the classifier on the training set

[Y_result1,scores1] = classify(net,train_im);
accuracy1 = sum(Y_result1 == train_classes)/numel(Y_result1);
disp(['The accuracy on the training set: ' num2str(accuracy1)]);


function c = convUnit(numF, stride, tag)
    c = [ ...
    convolution2dLayer(3,numF,'Padding','same','Stride',stride,'Name',[tag,'conv1']) ...
    dropoutLayer(.5,'Name',[tag,'DO1'])  ...
    maxPooling2dLayer(2,'Stride',2, 'Padding','same','Name',[tag,'MP1']) ...
    batchNormalizationLayer('Name',[tag,'BN1']) ...
    reluLayer('Name',[tag,'relu1'])...
    ];
end

function layers = convolutionalUnit(numF,stride,tag)
layers = [
    convolution2dLayer(3,numF,'Padding','same','Stride',stride,'Name',[tag,'conv1'])
    dropoutLayer(.5,'Name',[tag,'DO1']) 
    maxPooling2dLayer(2,'Stride',2, 'Padding','same','Name',[tag,'MP1']) 
    %batchNormalizationLayer('Name',[tag,'BN1']) 
    reluLayer('Name',[tag,'relu1']) 
    convolution2dLayer(3,numF,'Padding','same','Name',[tag,'conv2']) 
    dropoutLayer(.5,'Name',[tag,'DO2']) 
    maxPooling2dLayer(2,'Stride',2, 'Padding','same','Name',[tag,'MP2']) 
    %batchNormalizationLayer('Name',[tag,'BN2'])
    ];
end

