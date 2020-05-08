%% Load maps as ForceMap object into 1xNmaps cell.



prompt = 'Enter how many force maps you want to load for training:';
dlgtitle = 'Number of training maps';
dims = [1 35]; 
definput = {'2'};
Nmaps = str2double(inputdlg(prompt,dlgtitle,dims,definput));

% Try loading force maps from previous session. If that fails, let user
% choose to load from eighter an existing .mat ForceMap-object or from a
% folder with .csv files

try
    load('ForceMapFolders.mat');
    for i=1:length(objectfolders)
        FM{i} = ForceMap(objectfolders{i});
    end
catch ME
    warning('Could not find existing ForceMap objects. Proceeding with manual loading')
    rethrow(ME);
    for i=1:Nmaps
        FM{i} = ForceMap();
        objectfolders{i} = FM{i}.folder;
    end
    save('ForceMapFolders','objectfolders');
end

%% Prepare and manually select and label every single force curve.
% If you loaded already existing .mat files for the ForceMap objects on
% which the following operations have already been done, you can skip this
% section.

for i=20:Nmaps
    FM{i}.level_height_map();
end
for i=20:(Nmaps)
    FM{i}.choose_fibril(0.8);
end
for i=20:Nmaps
    FM{i}.base_and_tilt();
end
for i=20:Nmaps
    FM{i}.manual_CP();
end
for i=1:Nmaps
    FM{i}.CP_RoV();
end
for i=1:Nmaps
    FM{i}.old_CP;
end

%% Split the data into a training and a validation set with ratio 3:1
%  and set the NN layers and training options

ImgSize = 128; %bigger sizes improve results marginally but significantly
               %increase traning time. The trainer might even run out
               %of GPU memory... not really worth it

%[X,Y] = CP_CNN_batchprep(FM,ImgSize);
[X,Y] = CP_CNN_batchprep_alt(FM,ImgSize);

Val_idx = randperm(size(X,4),floor(size(X,4)/4));
XTrain = X;
YTrain = Y;
XValidation = X(:,:,:,Val_idx);
XTrain(:,:,:,Val_idx) = [];
YValidation = Y(Val_idx,:);
YTrain(Val_idx,:) = [];

%Defining shape of the CNN.
layers = [
    imageInputLayer([ImgSize ImgSize 1],"Name","imageinput",...
    'Normalization','zscore')
    convolution2dLayer([5 5],32,"Name","conv_1","Padding","same","Stride",[1 1])
    reluLayer("Name","relu_1")
    maxPooling2dLayer([5 5],"Name","maxpool_1","Padding","same","Stride",[3 3])
    convolution2dLayer([3 3],32,"Name","conv_2","Padding","same")
    reluLayer("Name","relu_2")
    maxPooling2dLayer([5 5],"Name","maxpool_2","Padding","same","Stride",[3 3])
    dropoutLayer(0.01,'Name','dropout_1')
    fullyConnectedLayer(2048,"Name","fc_1")
    reluLayer('Name','relu_3')
    dropoutLayer(0.01,'Name','dropout_2')
    fullyConnectedLayer(1024,'Name','fc_2')
    reluLayer('Name','relu_4')
    fullyConnectedLayer(2,"Name","fc_3")
    regressionLayer("Name","regressionoutput")];

resnet14 = layerGraph();

tempLayers = [
    imageInputLayer([ImgSize ImgSize 1],"Name","imageinput","Normalization","zscore")
    convolution2dLayer([7 7],64,"Name","conv_11","Padding","same","Stride",[2 2])
    batchNormalizationLayer("Name","batchnorm_9")
    reluLayer("Name","relu_9")
    maxPooling2dLayer([3 3],"Name","maxpool","Padding","same","Stride",[2 2])];
resnet14 = addLayers(resnet14,tempLayers);

tempLayers = [
    batchNormalizationLayer("Name","batchnorm_1")
    reluLayer("Name","relu_1")
    convolution2dLayer([3 3],64,"Name","conv_1","Padding","same")
    batchNormalizationLayer("Name","batchnorm_2")
    reluLayer("Name","relu_2")
    convolution2dLayer([3 3],64,"Name","conv_2","Padding","same")];
resnet14 = addLayers(resnet14,tempLayers);

tempLayers = additionLayer(2,"Name","addition_1");
resnet14 = addLayers(resnet14,tempLayers);

tempLayers = [
    batchNormalizationLayer("Name","batchnorm_3")
    reluLayer("Name","relu_3")
    convolution2dLayer([3 3],64,"Name","conv_3","Padding","same")
    batchNormalizationLayer("Name","batchnorm_4")
    reluLayer("Name","relu_4")
    convolution2dLayer([3 3],64,"Name","conv_4","Padding","same")];
resnet14 = addLayers(resnet14,tempLayers);

tempLayers = additionLayer(2,"Name","addition_2");
resnet14 = addLayers(resnet14,tempLayers);

tempLayers = [
    batchNormalizationLayer("Name","batchnorm_5")
    reluLayer("Name","relu_5")
    convolution2dLayer([3 3],128,"Name","conv_5","Padding","same","Stride",[2 2])
    batchNormalizationLayer("Name","batchnorm_6")
    reluLayer("Name","relu_6")
    convolution2dLayer([3 3],128,"Name","conv_6","Padding","same")];
resnet14 = addLayers(resnet14,tempLayers);

tempLayers = convolution2dLayer([1 1],128,"Name","conv_9","Padding","same","Stride",[2 2]);
resnet14 = addLayers(resnet14,tempLayers);

tempLayers = additionLayer(2,"Name","addition_3");
resnet14 = addLayers(resnet14,tempLayers);

tempLayers = convolution2dLayer([1 1],128,"Name","conv_10","Padding","same");
resnet14 = addLayers(resnet14,tempLayers);

tempLayers = [
    batchNormalizationLayer("Name","batchnorm_7")
    reluLayer("Name","relu_7")
    convolution2dLayer([3 3],128,"Name","conv_7","Padding","same")
    batchNormalizationLayer("Name","batchnorm_8")
    reluLayer("Name","relu_8")
    convolution2dLayer([3 3],128,"Name","conv_8","Padding","same")];
resnet14 = addLayers(resnet14,tempLayers);

tempLayers = additionLayer(2,"Name","addition_4");
resnet14 = addLayers(resnet14,tempLayers);

tempLayers = [
    batchNormalizationLayer("Name","batchnorm_10")
    reluLayer("Name","relu_10")
    convolution2dLayer([3 3],256,"Name","conv_12","Padding","same","Stride",[2 2])
    batchNormalizationLayer("Name","batchnorm_11")
    reluLayer("Name","relu_11")
    convolution2dLayer([3 3],256,"Name","conv_13","Padding","same")];
resnet14 = addLayers(resnet14,tempLayers);

tempLayers = convolution2dLayer([1 1],256,"Name","conv_16","Padding","same","Stride",[2 2]);
resnet14 = addLayers(resnet14,tempLayers);

tempLayers = additionLayer(2,"Name","addition_5");
resnet14 = addLayers(resnet14,tempLayers);

tempLayers = [
    batchNormalizationLayer("Name","batchnorm_12")
    reluLayer("Name","relu_12")
    convolution2dLayer([3 3],256,"Name","conv_14","Padding","same")
    batchNormalizationLayer("Name","batchnorm_13")
    reluLayer("Name","relu_13")
    convolution2dLayer([3 3],256,"Name","conv_15","Padding","same")];
resnet14 = addLayers(resnet14,tempLayers);

tempLayers = convolution2dLayer([1 1],256,"Name","conv_17","Padding","same");
resnet14 = addLayers(resnet14,tempLayers);

tempLayers = [
    additionLayer(2,"Name","addition_6")
    globalAveragePooling2dLayer("Name","gapool")
    fullyConnectedLayer(2,"Name","fc")
    regressionLayer("Name","regressionoutput")];
resnet14 = addLayers(resnet14,tempLayers);

resnet14 = connectLayers(resnet14,"maxpool","batchnorm_1");
resnet14 = connectLayers(resnet14,"maxpool","addition_1/in2");
resnet14 = connectLayers(resnet14,"conv_2","addition_1/in1");
resnet14 = connectLayers(resnet14,"addition_1","batchnorm_3");
resnet14 = connectLayers(resnet14,"addition_1","addition_2/in2");
resnet14 = connectLayers(resnet14,"conv_4","addition_2/in1");
resnet14 = connectLayers(resnet14,"addition_2","batchnorm_5");
resnet14 = connectLayers(resnet14,"addition_2","conv_9");
resnet14 = connectLayers(resnet14,"conv_6","addition_3/in1");
resnet14 = connectLayers(resnet14,"conv_9","addition_3/in2");
resnet14 = connectLayers(resnet14,"addition_3","conv_10");
resnet14 = connectLayers(resnet14,"addition_3","batchnorm_7");
resnet14 = connectLayers(resnet14,"conv_10","addition_4/in2");
resnet14 = connectLayers(resnet14,"conv_8","addition_4/in1");
resnet14 = connectLayers(resnet14,"addition_4","batchnorm_10");
resnet14 = connectLayers(resnet14,"addition_4","conv_16");
resnet14 = connectLayers(resnet14,"conv_16","addition_5/in2");
resnet14 = connectLayers(resnet14,"conv_13","addition_5/in1");
resnet14 = connectLayers(resnet14,"addition_5","batchnorm_12");
resnet14 = connectLayers(resnet14,"addition_5","conv_17");
resnet14 = connectLayers(resnet14,"conv_15","addition_6/in1");
resnet14 = connectLayers(resnet14,"conv_17","addition_6/in2");

options = trainingOptions('adam','Plots','training-progress',...
    'ValidationData',{XValidation,YValidation},...
    'ValidationFrequency',20,...
    'MiniBatchSize',20,...
    'MaxEpochs',500,...
    'InitialLearnRate',0.002,...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropPeriod',10,...
    'LearnRateDropFactor',0.5,...
    'Shuffle','every-epoch');

% analyzeNetwork(resnet14);
%% Start actual training of the NN. GPU usage for training is recommended.
%  For this option to be available, the Parallel-Computing-Toolbox has to be
%  installed in MATLAB

CP_CNN = trainNetwork(XTrain,YTrain,lgraph_2,options);

%% Evaluate your model looking at the models predictions

% predict outcome of XValidation
YPredicted = predict(CP_CNN,XValidation); 

% show random curve out of force map i and draw the manual CP in green and
% the CNNs CP in red

for i=1:Nmaps
FM{i}.CP_CNN_predict;
end

i = randi(Nmaps);
idxs = find(FM{i}.selected_curves);
rand = randperm(length(idxs),1);
randidx = idxs(rand);
plot(FM{i}.thapp{randidx},(FM{i}.basedapp{randidx}));
manual = drawpoint('Position',FM{i}.Man_CP(randidx,:),'Color','green');
net = drawpoint('Position',FM{i}.CP(randidx,:),'Color','red');
[CP_RoV(1,1),RoVidx] = max(FM{i}.CP_RoV{randidx});
CP_RoV(1,2) = FM{i}.basedapp{randidx}(RoVidx);
RoV = drawpoint('Position',CP_RoV,'Color','blue');

fig = figure('Name','Evaluate the model','Position',[949 79 971 915]);
idx = randi(length(XValidation));
imshow(XValidation(:,:,1,idx),'InitialMagnification','fit');
manpoint =drawpoint('Position',[YValidation(idx,1)*ImgSize,...
    (1-YValidation(idx,2))*ImgSize],'Color','green');
predpoint = drawpoint('Position',[YPredicted(idx,1)*ImgSize,...
    (1-YPredicted(idx,2))*ImgSize],'Color','red');

%% Create dropout model for uncertainty estimation during inference. Extract and apply weights
%  from previously trained model

DropParam = 0.2;

DropoutNet = layerGraph();

tempLayers = [
    CP_CNN.Layers(1)
    dropoutLayer(DropParam,"Name","dropout_1")
    CP_CNN.Layers(2)
    CP_CNN.Layers(3)
    reluLayer("Name","relu_9")
    CP_CNN.Layers(5)];
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = [
    CP_CNN.Layers(6)
    reluLayer("Name","relu_1")
    dropoutLayer(DropParam,"Name","dropout_2")
    CP_CNN.Layers(8)
    CP_CNN.Layers(9)
    reluLayer("Name","relu_2")
    dropoutLayer(DropParam,"Name","dropout_3")
    CP_CNN.Layers(11)];
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = additionLayer(2,"Name","addition_1");
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = [
    CP_CNN.Layers(13)
    reluLayer("Name","relu_3")
    dropoutLayer(DropParam,"Name","dropout_4")
    CP_CNN.Layers(15)
    CP_CNN.Layers(16)
    reluLayer("Name","relu_4")
    dropoutLayer(DropParam,"Name","dropout_5")
    CP_CNN.Layers(18)];
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = additionLayer(2,"Name","addition_2");
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = CP_CNN.Layers(26);
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = [
    CP_CNN.Layers(20)
    reluLayer("Name","relu_5")
    dropoutLayer(DropParam,"Name","dropout_6")
    CP_CNN.Layers(22)
    CP_CNN.Layers(23)
    reluLayer("Name","relu_6")
    dropoutLayer(DropParam,"Name","dropout_7")
    CP_CNN.Layers(25)];
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = additionLayer(2,"Name","addition_3");
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = CP_CNN.Layers(28);
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = [
    CP_CNN.Layers(29)
    reluLayer("Name","relu_7")
    dropoutLayer(DropParam,"Name","dropout_8")
    CP_CNN.Layers(31)
    CP_CNN.Layers(32)
    reluLayer("Name","relu_8")
    dropoutLayer(DropParam,"Name","dropout_9")
    CP_CNN.Layers(34)];
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = additionLayer(2,"Name","addition_4");
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = CP_CNN.Layers(42);
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = [
    CP_CNN.Layers(36)
    reluLayer("Name","relu_10")
    dropoutLayer(DropParam,"Name","dropout_10")
    CP_CNN.Layers(38)
    CP_CNN.Layers(39)
    reluLayer("Name","relu_11")
    dropoutLayer(DropParam,"Name","dropout_11")
    CP_CNN.Layers(41)];
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = additionLayer(2,"Name","addition_5");
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = [
    CP_CNN.Layers(44)
    reluLayer("Name","relu_12")
    dropoutLayer(DropParam,"Name","dropout_12")
    CP_CNN.Layers(46)
    CP_CNN.Layers(47)
    reluLayer("Name","relu_13")
    dropoutLayer(DropParam,"Name","dropout_13")
    CP_CNN.Layers(49)];
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = CP_CNN.Layers(50);
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = [
    additionLayer(2,"Name","addition_6")
    globalAveragePooling2dLayer("Name","gapool")
    CP_CNN.Layers(53)];
DropoutNet = addLayers(DropoutNet,tempLayers);

% clean up helper variable
clear tempLayers;

DropoutNet = connectLayers(DropoutNet,"maxpool","batchnorm_1");
DropoutNet = connectLayers(DropoutNet,"maxpool","addition_1/in2");
DropoutNet = connectLayers(DropoutNet,"conv_2","addition_1/in1");
DropoutNet = connectLayers(DropoutNet,"addition_1","batchnorm_3");
DropoutNet = connectLayers(DropoutNet,"addition_1","addition_2/in2");
DropoutNet = connectLayers(DropoutNet,"conv_4","addition_2/in1");
DropoutNet = connectLayers(DropoutNet,"addition_2","conv_9");
DropoutNet = connectLayers(DropoutNet,"addition_2","batchnorm_5");
DropoutNet = connectLayers(DropoutNet,"conv_9","addition_3/in2");
DropoutNet = connectLayers(DropoutNet,"conv_6","addition_3/in1");
DropoutNet = connectLayers(DropoutNet,"addition_3","conv_10");
DropoutNet = connectLayers(DropoutNet,"addition_3","batchnorm_7");
DropoutNet = connectLayers(DropoutNet,"conv_10","addition_4/in2");
DropoutNet = connectLayers(DropoutNet,"conv_8","addition_4/in1");
DropoutNet = connectLayers(DropoutNet,"addition_4","conv_16");
DropoutNet = connectLayers(DropoutNet,"addition_4","batchnorm_10");
DropoutNet = connectLayers(DropoutNet,"conv_16","addition_5/in2");
DropoutNet = connectLayers(DropoutNet,"conv_13","addition_5/in1");
DropoutNet = connectLayers(DropoutNet,"addition_5","batchnorm_12");
DropoutNet = connectLayers(DropoutNet,"addition_5","conv_17");
DropoutNet = connectLayers(DropoutNet,"conv_17","addition_6/in2");
DropoutNet = connectLayers(DropoutNet,"conv_15","addition_6/in1");

dlXValidation = dlarray(XValidation,'SSCB');
DropoutNet = dlnetwork(DropoutNet);
forward(DropoutNet,dlXValidation(:,:,:,1));
predict(CP_CNN,dlXValidation(:,:,:,1));

%% Define and train CNN, that tries to predict the regression error 
%  from the previous Network in order to get some kind of confidence metric
%  for our Network

% predict outcome of XValidation
YTrainPredicted = predict(CP_CNN,XTrain);
YValidationPredicted = predict(CP_CNN,XValidation);

% Define new target-values ConfY instead of Y. X stayes the same
ConfXTrain = XTrain;
ConfXValidation = XValidation;
ConfYTrain = zeros(length(YTrain),1);
ConfYValidation = zeros(length(YValidation),1);
for i=1:length(YTrain)
    ConfYTrain(i) = norm(YTrain(i,:)-YTrainPredicted(i,:));
end
for i=1:length(YValidation)
    ConfYValidation(i) = norm(YValidation(i,:)-YValidationPredicted(i,:));
end


% Defining shape of the CNN. This is the same architecture as the first one
% used above
ConfLayers = [
    imageInputLayer([ImgSize ImgSize 1],"Name","imageinput",...
    'Normalization','zscore')
    convolution2dLayer([5 5],32,"Name","conv_1","Padding","same","Stride",[1 1])
    reluLayer("Name","relu_1")
    maxPooling2dLayer([5 5],"Name","maxpool_1","Padding","same","Stride",[3 3])
    convolution2dLayer([3 3],32,"Name","conv_2","Padding","same")
    reluLayer("Name","relu_2")
    maxPooling2dLayer([5 5],"Name","maxpool_2","Padding","same","Stride",[3 3])
    dropoutLayer(DropParam,'Name','dropout_1')
    fullyConnectedLayer(2048,"Name","fc_1")
    reluLayer('Name','relu_3')
    dropoutLayer(DropParam,'Name','dropout_2')
    fullyConnectedLayer(1024,'Name','fc_2')
    reluLayer('Name','relu_4')
    fullyConnectedLayer(1,"Name","fc_3")
    regressionLayer("Name","regressionoutput")];

ConfResnet14 = layerGraph();

tempLayers = [
    imageInputLayer([ImgSize ImgSize 1],"Name","imageinput","Normalization","zscore")
    convolution2dLayer([7 7],64,"Name","conv_11","Padding","same","Stride",[2 2])
    batchNormalizationLayer("Name","batchnorm_9")
    reluLayer("Name","relu_9")
    maxPooling2dLayer([3 3],"Name","maxpool","Padding","same","Stride",[2 2])];
ConfResnet14 = addLayers(ConfResnet14,tempLayers);

tempLayers = [
    batchNormalizationLayer("Name","batchnorm_1")
    reluLayer("Name","relu_1")
    convolution2dLayer([3 3],64,"Name","conv_1","Padding","same")
    batchNormalizationLayer("Name","batchnorm_2")
    reluLayer("Name","relu_2")
    convolution2dLayer([3 3],64,"Name","conv_2","Padding","same")];
ConfResnet14 = addLayers(ConfResnet14,tempLayers);

tempLayers = additionLayer(2,"Name","addition_1");
ConfResnet14 = addLayers(ConfResnet14,tempLayers);

tempLayers = [
    batchNormalizationLayer("Name","batchnorm_3")
    reluLayer("Name","relu_3")
    convolution2dLayer([3 3],64,"Name","conv_3","Padding","same")
    batchNormalizationLayer("Name","batchnorm_4")
    reluLayer("Name","relu_4")
    convolution2dLayer([3 3],64,"Name","conv_4","Padding","same")];
ConfResnet14 = addLayers(ConfResnet14,tempLayers);

tempLayers = additionLayer(2,"Name","addition_2");
ConfResnet14 = addLayers(ConfResnet14,tempLayers);

tempLayers = [
    batchNormalizationLayer("Name","batchnorm_5")
    reluLayer("Name","relu_5")
    convolution2dLayer([3 3],128,"Name","conv_5","Padding","same","Stride",[2 2])
    batchNormalizationLayer("Name","batchnorm_6")
    reluLayer("Name","relu_6")
    convolution2dLayer([3 3],128,"Name","conv_6","Padding","same")];
ConfResnet14 = addLayers(ConfResnet14,tempLayers);

tempLayers = convolution2dLayer([1 1],128,"Name","conv_9","Padding","same","Stride",[2 2]);
ConfResnet14 = addLayers(ConfResnet14,tempLayers);

tempLayers = additionLayer(2,"Name","addition_3");
ConfResnet14 = addLayers(ConfResnet14,tempLayers);

tempLayers = convolution2dLayer([1 1],128,"Name","conv_10","Padding","same");
ConfResnet14 = addLayers(ConfResnet14,tempLayers);

tempLayers = [
    batchNormalizationLayer("Name","batchnorm_7")
    reluLayer("Name","relu_7")
    convolution2dLayer([3 3],128,"Name","conv_7","Padding","same")
    batchNormalizationLayer("Name","batchnorm_8")
    reluLayer("Name","relu_8")
    convolution2dLayer([3 3],128,"Name","conv_8","Padding","same")];
ConfResnet14 = addLayers(ConfResnet14,tempLayers);

tempLayers = additionLayer(2,"Name","addition_4");
ConfResnet14 = addLayers(ConfResnet14,tempLayers);

tempLayers = [
    batchNormalizationLayer("Name","batchnorm_10")
    reluLayer("Name","relu_10")
    convolution2dLayer([3 3],256,"Name","conv_12","Padding","same","Stride",[2 2])
    batchNormalizationLayer("Name","batchnorm_11")
    reluLayer("Name","relu_11")
    convolution2dLayer([3 3],256,"Name","conv_13","Padding","same")];
ConfResnet14 = addLayers(ConfResnet14,tempLayers);

tempLayers = convolution2dLayer([1 1],256,"Name","conv_16","Padding","same","Stride",[2 2]);
ConfResnet14 = addLayers(ConfResnet14,tempLayers);

tempLayers = additionLayer(2,"Name","addition_5");
ConfResnet14 = addLayers(ConfResnet14,tempLayers);

tempLayers = [
    batchNormalizationLayer("Name","batchnorm_12")
    reluLayer("Name","relu_12")
    convolution2dLayer([3 3],256,"Name","conv_14","Padding","same")
    batchNormalizationLayer("Name","batchnorm_13")
    reluLayer("Name","relu_13")
    convolution2dLayer([3 3],256,"Name","conv_15","Padding","same")];
ConfResnet14 = addLayers(ConfResnet14,tempLayers);

tempLayers = convolution2dLayer([1 1],256,"Name","conv_17","Padding","same");
ConfResnet14 = addLayers(ConfResnet14,tempLayers);

tempLayers = [
    additionLayer(2,"Name","addition_6")
    globalAveragePooling2dLayer("Name","gapool")
    fullyConnectedLayer(1,"Name","fc")
    regressionLayer("Name","regressionoutput")];
ConfResnet14 = addLayers(ConfResnet14,tempLayers);

ConfResnet14 = connectLayers(ConfResnet14,"maxpool","batchnorm_1");
ConfResnet14 = connectLayers(ConfResnet14,"maxpool","addition_1/in2");
ConfResnet14 = connectLayers(ConfResnet14,"conv_2","addition_1/in1");
ConfResnet14 = connectLayers(ConfResnet14,"addition_1","batchnorm_3");
ConfResnet14 = connectLayers(ConfResnet14,"addition_1","addition_2/in2");
ConfResnet14 = connectLayers(ConfResnet14,"conv_4","addition_2/in1");
ConfResnet14 = connectLayers(ConfResnet14,"addition_2","batchnorm_5");
ConfResnet14 = connectLayers(ConfResnet14,"addition_2","conv_9");
ConfResnet14 = connectLayers(ConfResnet14,"conv_6","addition_3/in1");
ConfResnet14 = connectLayers(ConfResnet14,"conv_9","addition_3/in2");
ConfResnet14 = connectLayers(ConfResnet14,"addition_3","conv_10");
ConfResnet14 = connectLayers(ConfResnet14,"addition_3","batchnorm_7");
ConfResnet14 = connectLayers(ConfResnet14,"conv_10","addition_4/in2");
ConfResnet14 = connectLayers(ConfResnet14,"conv_8","addition_4/in1");
ConfResnet14 = connectLayers(ConfResnet14,"addition_4","batchnorm_10");
ConfResnet14 = connectLayers(ConfResnet14,"addition_4","conv_16");
ConfResnet14 = connectLayers(ConfResnet14,"conv_16","addition_5/in2");
ConfResnet14 = connectLayers(ConfResnet14,"conv_13","addition_5/in1");
ConfResnet14 = connectLayers(ConfResnet14,"addition_5","batchnorm_12");
ConfResnet14 = connectLayers(ConfResnet14,"addition_5","conv_17");
ConfResnet14 = connectLayers(ConfResnet14,"conv_15","addition_6/in1");
ConfResnet14 = connectLayers(ConfResnet14,"conv_17","addition_6/in2");

options = trainingOptions('adam','Plots','training-progress',...
    'ValidationData',{ConfXValidation,ConfYValidation},...
    'ValidationFrequency',50,...
    'MiniBatchSize',50,...
    'MaxEpochs',1000,...
    'InitialLearnRate',0.002,...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropPeriod',50,...
    'LearnRateDropFactor',0.5,...
    'Shuffle','every-epoch');


Confidence_CP_CNN = trainNetwork(ConfXTrain,ConfYTrain,ConfResnet14,options);