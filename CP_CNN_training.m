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
        objectfolders{i} = FM{i}.Folder;
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
    FM{i}.manual_cp();
end
for i=1:Nmaps
    FM{i}.cp_rov();
end
for i=1:Nmaps
    FM{i}.old_cp;
end

%% Split the data into a training and a validation set with ratio 3:1
%  and set the NN layers and training options

ImgSize = 128; %bigger sizes improve results marginally but significantly
               %increase traning time. The trainer might even run out
               %of GPU memory... not really worth it

k = 1;
for i = A
objcell{k} = E.FM{i};
k = k + 1;
end               

[X,Y] = CP_CNN_batchprep(objcell,ImgSize);

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
DropParam = 0.2;

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

tempLayers = convolution2dLayer([1 1],128,"Name","conv_9","Padding","same","Stride",[2 2],...
    'BiasInitializer','zeros','BiasLearnRateFactor',0);
resnet14 = addLayers(resnet14,tempLayers);

tempLayers = additionLayer(2,"Name","addition_3");
resnet14 = addLayers(resnet14,tempLayers);

tempLayers = convolution2dLayer([1 1],128,"Name","conv_10","Padding","same",...
    'BiasInitializer','zeros','BiasLearnRateFactor',0);
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

tempLayers = convolution2dLayer([1 1],256,"Name","conv_16","Padding","same","Stride",[2 2],...
    'BiasInitializer','zeros','BiasLearnRateFactor',0);
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

tempLayers = convolution2dLayer([1 1],256,"Name","conv_17","Padding","same",...
    'BiasInitializer','zeros','BiasLearnRateFactor',0);
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

MonteCarlo14 =  layerGraph();
DropParam = 0.2;

tempLayers = [
    imageInputLayer([128 128 1],"Name","imageinput","Normalization","zscore")
    dropoutLayer(DropParam,"Name","dropout_1")
    convolution2dLayer([7 7],64,"Name","conv_11","Padding","same","Stride",[2 2])
    batchNormalizationLayer("Name","batchnorm_9")
    reluLayer("Name","relu_9")
    maxPooling2dLayer([3 3],"Name","maxpool","Padding","same","Stride",[2 2])];
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = [
    reluLayer("Name","relu_1")
    dropoutLayer(DropParam,"Name","dropout_2")
    convolution2dLayer([3 3],64,"Name","conv_1","Padding","same")
    batchNormalizationLayer("Name","batchnorm_2")
    reluLayer("Name","relu_2")
    dropoutLayer(DropParam,"Name","dropout_3")
    convolution2dLayer([3 3],64,"Name","conv_2","Padding","same")
    ];
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = additionLayer(2,"Name","addition_1");
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = [
    batchNormalizationLayer("Name","batchnorm_3")
    reluLayer("Name","relu_3")
    dropoutLayer(DropParam,"Name","dropout_4")
    convolution2dLayer([3 3],64,"Name","conv_3","Padding","same")
    batchNormalizationLayer("Name","batchnorm_4")
    reluLayer("Name","relu_4")
    dropoutLayer(DropParam,"Name","dropout_5")
    convolution2dLayer([3 3],64,"Name","conv_4","Padding","same")
    ];
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = additionLayer(2,"Name","addition_2");
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = [
    batchNormalizationLayer("Name","batchnorm_5")
    reluLayer("Name","relu_5")
    dropoutLayer(DropParam,"Name","dropout_6")
    convolution2dLayer([3 3],128,"Name","conv_5","Padding","same","Stride",[2 2])
    batchNormalizationLayer("Name","batchnorm_6")
    reluLayer("Name","relu_6")
    dropoutLayer(DropParam,"Name","dropout_7")
    convolution2dLayer([3 3],128,"Name","conv_6","Padding","same")
    ];
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = [
    dropoutLayer(DropParam,"Name","dropout_8")
    convolution2dLayer([1 1],128,"Name","conv_9","Padding","same","Stride",[2 2],...
    'BiasInitializer','zeros','BiasLearnRateFactor',0)
    ];
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = additionLayer(2,"Name","addition_3");
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = [
    batchNormalizationLayer("Name","batchnorm_7")
    reluLayer("Name","relu_7")
    dropoutLayer(DropParam,"Name","dropout_9")
    convolution2dLayer([3 3],128,"Name","conv_7","Padding","same")
    batchNormalizationLayer("Name","batchnorm_8")
    reluLayer("Name","relu_8")
    dropoutLayer(DropParam,"Name","dropout_10")
    convolution2dLayer([3 3],128,"Name","conv_8","Padding","same")
    ];
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = [
    dropoutLayer(DropParam,"Name","dropout_11")
    convolution2dLayer([1 1],128,"Name","conv_10","Padding","same",...
    'BiasInitializer','zeros','BiasLearnRateFactor',0)
    ];
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = additionLayer(2,"Name","addition_4");
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = [
    batchNormalizationLayer("Name","batchnorm_10")
    reluLayer("Name","relu_10")
    dropoutLayer(DropParam,"Name","dropout_12")
    convolution2dLayer([3 3],256,"Name","conv_12","Padding","same","Stride",[2 2])
    batchNormalizationLayer("Name","batchnorm_11")
    reluLayer("Name","relu_11")
    dropoutLayer(DropParam,"Name","dropout_14")
    convolution2dLayer([3 3],256,"Name","conv_13","Padding","same")
    ];
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = [
    dropoutLayer(DropParam,"Name","dropout_13")
    convolution2dLayer([1 1],256,"Name","conv_16","Padding","same","Stride",[2 2],...
    'BiasInitializer','zeros','BiasLearnRateFactor',0)
    ];
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = additionLayer(2,"Name","addition_5");
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = [
    dropoutLayer(DropParam,"Name","dropout_16")
    convolution2dLayer([1 1],256,"Name","conv_17","Padding","same",...
    'BiasInitializer','zeros','BiasLearnRateFactor',0)];
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = [
    batchNormalizationLayer("Name","batchnorm_12")
    reluLayer("Name","relu_12")
    dropoutLayer(DropParam,"Name","dropout_15")
    convolution2dLayer([3 3],256,"Name","conv_14","Padding","same")
    batchNormalizationLayer("Name","batchnorm_13")
    reluLayer("Name","relu_13")
    dropoutLayer(DropParam,"Name","dropout_17")
    convolution2dLayer([3 3],256,"Name","conv_15","Padding","same")
    ];
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = [
    additionLayer(2,"Name","addition_6")
    globalAveragePooling2dLayer("Name","gapool")
    dropoutLayer(DropParam,"Name","dropout_18")
    fullyConnectedLayer(2,"Name","fc")
    regressionLayer("Name","regressionoutput")];
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

% clean up helper variable
clear tempLayers;

MonteCarlo14 = connectLayers(MonteCarlo14,"maxpool","relu_1");
MonteCarlo14 = connectLayers(MonteCarlo14,"maxpool","addition_1/in1");
MonteCarlo14 = connectLayers(MonteCarlo14,"conv_2","addition_1/in2");
MonteCarlo14 = connectLayers(MonteCarlo14,"addition_1","batchnorm_3");
MonteCarlo14 = connectLayers(MonteCarlo14,"addition_1","addition_2/in1");
MonteCarlo14 = connectLayers(MonteCarlo14,"conv_4","addition_2/in2");
MonteCarlo14 = connectLayers(MonteCarlo14,"addition_2","dropout_8");
MonteCarlo14 = connectLayers(MonteCarlo14,"addition_2","batchnorm_5");
MonteCarlo14 = connectLayers(MonteCarlo14,"conv_9","addition_3/in2");
MonteCarlo14 = connectLayers(MonteCarlo14,"conv_6","addition_3/in1");
MonteCarlo14 = connectLayers(MonteCarlo14,"addition_3","batchnorm_7");
MonteCarlo14 = connectLayers(MonteCarlo14,"addition_3","dropout_11");
MonteCarlo14 = connectLayers(MonteCarlo14,"conv_8","addition_4/in1");
MonteCarlo14 = connectLayers(MonteCarlo14,"conv_10","addition_4/in2");
MonteCarlo14 = connectLayers(MonteCarlo14,"addition_4","batchnorm_10");
MonteCarlo14 = connectLayers(MonteCarlo14,"addition_4","dropout_13");
MonteCarlo14 = connectLayers(MonteCarlo14,"conv_16","addition_5/in2");
MonteCarlo14 = connectLayers(MonteCarlo14,"conv_13","addition_5/in1");
MonteCarlo14 = connectLayers(MonteCarlo14,"addition_5","batchnorm_12");
MonteCarlo14 = connectLayers(MonteCarlo14,"addition_5","dropout_16");
MonteCarlo14 = connectLayers(MonteCarlo14,"conv_17","addition_6/in2");
MonteCarlo14 = connectLayers(MonteCarlo14,"conv_15","addition_6/in1");


options = trainingOptions('adam','Plots','training-progress',...
    'ValidationData',{XValidation,YValidation},...
    'ValidationFrequency',50,...
    'MiniBatchSize',100,...
    'MaxEpochs',1000,...
    'InitialLearnRate',0.002,...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropPeriod',100,...
    'LearnRateDropFactor',0.5,...
    'Shuffle','every-epoch');

% analyzeNetwork(resnet14);
%% Start actual training of the NN. GPU usage for training is recommended.
%  For this option to be available, the Parallel-Computing-Toolbox has to be
%  installed in MATLAB

CP_CNN_29_10 = trainNetwork(XTrain,YTrain,MonteCarlo14,options);

%% Evaluate your model looking at the models predictions

% predict outcome of XValidation
YPredicted = predict(CP_CNN,XValidation); 
YPredictedMC = predict(CP_CNN_29_10,XValidation);

% show random curve out of force map i and draw the manual CP in green and
% the CNNs CP in red

for i=1:Nmaps
FM{i}.estimate_cp_cnn('Fast');
end

i = randi(Nmaps);
idxs = find(FM{i}.SelectedCurves);
rand = randperm(length(idxs),1);
randidx = idxs(rand);
plot(FM{i}.THApp{randidx},(FM{i}.BasedApp{randidx}));
manual = drawpoint('Position',FM{i}.Man_CP(randidx,:),'Color','green');
net = drawpoint('Position',FM{i}.CP(randidx,:),'Color','red');
[CP_RoV(1,1),RoVidx] = max(FM{i}.cp_rov{randidx});
CP_RoV(1,2) = FM{i}.BasedApp{randidx}(RoVidx);
RoV = drawpoint('Position',CP_RoV,'Color','blue');

fig = figure('Name','Evaluate the model','Position',[949 79 971 915]);
for i=1:1115
    Diff(i) = norm(YPredictedMC(i,:)-YValidation(i,:));
end
k=1;
% [~,IdxVec] = sort(Diff,'descend');
idx = randi(length(XValidation));
% idx = IdxVec(k);
imshow(XValidation(:,:,1,idx),'InitialMagnification','fit');
manpoint =drawpoint('Position',[YValidation(idx,1)*ImgSize,...
    (1-YValidation(idx,2))*ImgSize],'Color','green');
% predpoint = drawpoint('Position',[YPredicted(idx,1)*ImgSize,...
%     (1-YPredicted(idx,2))*ImgSize],'Color','red');
predpoint2 = drawpoint('Position',[YPredictedMC(idx,1)*ImgSize,...
    (1-YPredictedMC(idx,2))*ImgSize],'Color','yellow');
k = k + 1;


%% Create dropout model for uncertainty estimation during inference. Extract and apply weights
%  from previously trained model

DropoutNet = assemble_DON(CP_CNN,0.1);
MonteCarloDrop = assemble_MC(CP_CNN_MC);
% predict(DropoutNet,X(:,:,:,1))
% predict(CP_CNN,X(:,:,:,1))

