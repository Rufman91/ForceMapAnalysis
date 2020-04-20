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

for i=1:Nmaps
    FM{i}.level_height_map();
end
for i=1:(Nmaps)
    FM{i}.choose_fibril(0.8);
end
for i=1:Nmaps
    FM{i}.base_and_tilt();
end
for i=12:Nmaps
    FM{i}.manual_CP();
end


%% Split the data into a training and a validation set with ratio 3:1
%  and set the NN layers and training options

ImgSize = 128; %bigger sizes improve results marginally but significantly
               %increase traning time. The trainer might even run out
               %of GPU memory

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
    fullyConnectedLayer(2048,"Name","fc_1")
    reluLayer('Name','relu_3')
    fullyConnectedLayer(1024,'Name','fc_2')
    reluLayer('Name','relu_4')
    fullyConnectedLayer(2,"Name","fc_3")
    regressionLayer("Name","regressionoutput")];

lgraph = layerGraph();

tempLayers = [
    imageInputLayer([128 128 1],"Name","imageinput","Normalization","zscore")
    convolution2dLayer([7 7],64,"Name","conv_11","Padding","same","Stride",[2 2])
    batchNormalizationLayer("Name","batchnorm_9")
    reluLayer("Name","relu_9")
    maxPooling2dLayer([3 3],"Name","maxpool","Padding","same","Stride",[2 2])];
lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    batchNormalizationLayer("Name","batchnorm_1")
    reluLayer("Name","relu_1")
    convolution2dLayer([3 3],64,"Name","conv_1","Padding","same")
    batchNormalizationLayer("Name","batchnorm_2")
    reluLayer("Name","relu_2")
    convolution2dLayer([3 3],64,"Name","conv_2","Padding","same")];
lgraph = addLayers(lgraph,tempLayers);

tempLayers = additionLayer(2,"Name","addition_1");
lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    batchNormalizationLayer("Name","batchnorm_3")
    reluLayer("Name","relu_3")
    convolution2dLayer([3 3],64,"Name","conv_3","Padding","same")
    batchNormalizationLayer("Name","batchnorm_4")
    reluLayer("Name","relu_4")
    convolution2dLayer([3 3],64,"Name","conv_4","Padding","same")];
lgraph = addLayers(lgraph,tempLayers);

tempLayers = additionLayer(2,"Name","addition_2");
lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    batchNormalizationLayer("Name","batchnorm_5")
    reluLayer("Name","relu_5")
    convolution2dLayer([3 3],128,"Name","conv_5","Padding","same","Stride",[2 2])
    batchNormalizationLayer("Name","batchnorm_6")
    reluLayer("Name","relu_6")
    convolution2dLayer([3 3],128,"Name","conv_6","Padding","same")];
lgraph = addLayers(lgraph,tempLayers);

tempLayers = convolution2dLayer([1 1],128,"Name","conv_9","Padding","same","Stride",[2 2]);
lgraph = addLayers(lgraph,tempLayers);

tempLayers = additionLayer(2,"Name","addition_3");
lgraph = addLayers(lgraph,tempLayers);

tempLayers = convolution2dLayer([1 1],128,"Name","conv_10","Padding","same");
lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    batchNormalizationLayer("Name","batchnorm_7")
    reluLayer("Name","relu_7")
    convolution2dLayer([3 3],128,"Name","conv_7","Padding","same")
    batchNormalizationLayer("Name","batchnorm_8")
    reluLayer("Name","relu_8")
    convolution2dLayer([3 3],128,"Name","conv_8","Padding","same")];
lgraph = addLayers(lgraph,tempLayers);

tempLayers = additionLayer(2,"Name","addition_4");
lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    batchNormalizationLayer("Name","batchnorm_10")
    reluLayer("Name","relu_10")
    convolution2dLayer([3 3],256,"Name","conv_12","Padding","same","Stride",[2 2])
    batchNormalizationLayer("Name","batchnorm_11")
    reluLayer("Name","relu_11")
    convolution2dLayer([3 3],256,"Name","conv_13","Padding","same")];
lgraph = addLayers(lgraph,tempLayers);

tempLayers = convolution2dLayer([1 1],256,"Name","conv_16","Padding","same","Stride",[2 2]);
lgraph = addLayers(lgraph,tempLayers);

tempLayers = additionLayer(2,"Name","addition_5");
lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    batchNormalizationLayer("Name","batchnorm_12")
    reluLayer("Name","relu_12")
    convolution2dLayer([3 3],256,"Name","conv_14","Padding","same")
    batchNormalizationLayer("Name","batchnorm_13")
    reluLayer("Name","relu_13")
    convolution2dLayer([3 3],256,"Name","conv_15","Padding","same")];
lgraph = addLayers(lgraph,tempLayers);

tempLayers = convolution2dLayer([1 1],256,"Name","conv_17","Padding","same");
lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    additionLayer(2,"Name","addition_6")
    globalAveragePooling2dLayer("Name","gapool")
    fullyConnectedLayer(2,"Name","fc")
    regressionLayer("Name","regressionoutput")];
lgraph = addLayers(lgraph,tempLayers);

lgraph = connectLayers(lgraph,"maxpool","batchnorm_1");
lgraph = connectLayers(lgraph,"maxpool","addition_1/in2");
lgraph = connectLayers(lgraph,"conv_2","addition_1/in1");
lgraph = connectLayers(lgraph,"addition_1","batchnorm_3");
lgraph = connectLayers(lgraph,"addition_1","addition_2/in2");
lgraph = connectLayers(lgraph,"conv_4","addition_2/in1");
lgraph = connectLayers(lgraph,"addition_2","batchnorm_5");
lgraph = connectLayers(lgraph,"addition_2","conv_9");
lgraph = connectLayers(lgraph,"conv_6","addition_3/in1");
lgraph = connectLayers(lgraph,"conv_9","addition_3/in2");
lgraph = connectLayers(lgraph,"addition_3","conv_10");
lgraph = connectLayers(lgraph,"addition_3","batchnorm_7");
lgraph = connectLayers(lgraph,"conv_10","addition_4/in2");
lgraph = connectLayers(lgraph,"conv_8","addition_4/in1");
lgraph = connectLayers(lgraph,"addition_4","batchnorm_10");
lgraph = connectLayers(lgraph,"addition_4","conv_16");
lgraph = connectLayers(lgraph,"conv_16","addition_5/in2");
lgraph = connectLayers(lgraph,"conv_13","addition_5/in1");
lgraph = connectLayers(lgraph,"addition_5","batchnorm_12");
lgraph = connectLayers(lgraph,"addition_5","conv_17");
lgraph = connectLayers(lgraph,"conv_15","addition_6/in1");
lgraph = connectLayers(lgraph,"conv_17","addition_6/in2");

options = trainingOptions('adam','Plots','training-progress',...
    'ValidationData',{XValidation,YValidation},...
    'ValidationFrequency',50,...
    'MiniBatchSize',90,...
    'MaxEpochs',1000,...
    'InitialLearnRate',0.002,...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropPeriod',100,...
    'LearnRateDropFactor',0.5,...
    'Shuffle','every-epoch');

%% Start actual training of the NN. GPU usage for training is recommended.
%  For this option to be available, the Parallel-Computing-Toolbox has to be
%  installed in MATLAB

CP_CNN = trainNetwork(XTrain,YTrain,lgraph,options);

%% Evaluate your model looking at the models predictions

% show random curve out of force map i and draw the manual CP in green and
% the CNNs CP in red
i = randi(19);
CP_CNN_predict(FM{i},CP_CNN);
idxs = find(FM{i}.selected_curves);
rand = randperm(length(idxs),1);
randidx = idxs(rand);
plot(FM{i}.thapp{randidx},(FM{i}.basedapp{randidx}));
manual = drawpoint('Position',FM{i}.Man_CP(randidx,:),'Color','green');
net = drawpoint('Position',FM{i}.CP(randidx,:),'Color','red');

YPredicted = predict(CP_CNN,XValidation); 
fig = figure('Name','Evaluate the model','Position',[949 79 971 915])
idx = randi(length(XValidation));
imshow(XValidation(:,:,1,idx),'InitialMagnification','fit');
manpoint =drawpoint('Position',[YValidation(idx,1)*ImgSize , (1-YValidation(idx,2))*ImgSize],'Color','green');
predpoint = drawpoint('Position',[YPredicted(idx,1)*ImgSize , (1-YPredicted(idx,2))*ImgSize],'Color','red');