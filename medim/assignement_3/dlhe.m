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
pixelRange = [-8 8];
imageAugmenter = imageDataAugmenter( ...
    'RandXReflection',true, ...
    'RandYReflection',true, ...
    'RandRotation', [0 180], ...
    'RandXTranslation',pixelRange, ...
    'RandYTranslation',pixelRange);
augimdsTrain = augmentedImageDatastore(imageSize,train_im,train_classes, ...
    'DataAugmentation',imageAugmenter, ...
    'OutputSizeMode','randcrop');

drop_prob = 0.5;
netWidth = 8;
layers = [
    imageInputLayer([64 64 1],'Name','inp')
    convolution2dLayer(5,netWidth, 'Stride',1, 'Padding','same','Name','convInp')
    batchNormalizationLayer('Name','BNInp')
    reluLayer('Name','reluInp')
    
    convolutionalUnit(netWidth,1,'S1U1')
    additionLayer(2,'Name','add11')
    reluLayer('Name','relu11')
    convolutionalUnit(netWidth,1,'S1U2')
    additionLayer(2,'Name','add12')
    reluLayer('Name','relu12')
    
    convolutionalUnit(2*netWidth,2,'S2U1')
    additionLayer(2,'Name','add21')
    reluLayer('Name','relu21')
    convolutionalUnit(2*netWidth,1,'S2U2')
    additionLayer(2,'Name','add22')
    reluLayer('Name','relu22')
    
    convolutionalUnit(4*netWidth,2,'S3U1')
    additionLayer(2,'Name','add31')
    reluLayer('Name','relu31')
    convolutionalUnit(4*netWidth,1,'S3U2')
    additionLayer(2,'Name','add32')
    reluLayer('Name','relu32')
    
    averagePooling2dLayer(4,'Name','globalPool')
    %dropoutLayer(.5,'Name','DOF') 
    fullyConnectedLayer(6,'Name','fcFinal')
    softmaxLayer('Name','softmax')
    classificationLayer('Name','class')];

% slightly deeper 

%% Train deep learning network
lgraph = layerGraph(layers);

lgraph = connectLayers(lgraph,'reluInp','add11/in2');
lgraph = connectLayers(lgraph,'relu11','add12/in2');
skip1 = [
    convolution2dLayer(1,2*netWidth,'Stride',2,'Name','skipConv1')
    batchNormalizationLayer('Name','skipBN1')];
lgraph = addLayers(lgraph,skip1);
lgraph = connectLayers(lgraph,'relu12','skipConv1');
lgraph = connectLayers(lgraph,'skipBN1','add21/in2');
lgraph = connectLayers(lgraph,'relu21','add22/in2');

skip2 = [
    convolution2dLayer(1,4*netWidth,'Stride',2,'Name','skipConv2')
    batchNormalizationLayer('Name','skipBN2')];
lgraph = addLayers(lgraph,skip2);
lgraph = connectLayers(lgraph,'relu22','skipConv2');
lgraph = connectLayers(lgraph,'skipBN2','add31/in2');
lgraph = connectLayers(lgraph,'relu31','add32/in2');


plot(lgraph);

miniBatchSize = 128;
max_epochs = 300;
learning_rate = 0.001;
options = trainingOptions( 'sgdm',...
    'Momentum', 0.9, ...
    'MaxEpochs',max_epochs,...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.8, ...
    'LearnRateDropPeriod',10, ...
    'Shuffle','every-epoch', ...
    'InitialLearnRate',learning_rate, ...
    'ValidationData',{test_im,test_classes}, ...
    'ValidationFrequency',10, ...
    'Plots', 'training-progress');

net = trainNetwork(augimdsTrain, lgraph, options);


%% Test the classifier on the test set

[Y_result2,scores2] = classify(net,test_im);
accuracy2 = sum(Y_result2 == test_classes)/numel(Y_result2);
disp(['The accuracy on the test set: ' num2str(accuracy2)]);

%% Test the classifier on the training set

[Y_result1,scores1] = classify(net,train_im);
accuracy1 = sum(Y_result1 == train_classes)/numel(Y_result1);
disp(['The accuracy on the training set: ' num2str(accuracy1)]);

function layers = convolutionalUnit(numF,stride,tag)
layers = [
    convolution2dLayer(3,numF,'Padding','same','Stride',stride,'Name',[tag,'conv1'])
    %dropoutLayer(.5,'Name',[tag,'DO1']) 
    %maxPooling2dLayer(2,'Stride',2, 'Padding','same','Name',[tag,'MP1']) 
    batchNormalizationLayer('Name',[tag,'BN1']) 
    reluLayer('Name',[tag,'relu1']) 
    convolution2dLayer(3,numF,'Padding','same','Name',[tag,'conv2']) 
    %dropoutLayer(.5,'Name',[tag,'DO2']) 
    %maxPooling2dLayer(2,'Stride',2, 'Padding','same','Name',[tag,'MP2']) 
    batchNormalizationLayer('Name',[tag,'BN2'])
    ];
end

