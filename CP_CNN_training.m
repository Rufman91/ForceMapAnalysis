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



options = trainingOptions('adam','Plots','training-progress',...
    'ValidationData',{XValidation,YValidation},...
    'ValidationFrequency',50,...
    'MiniBatchSize',60,...
    'MaxEpochs',1000,...
    'InitialLearnRate',0.002,...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropPeriod',100,...
    'LearnRateDropFactor',0.5,...
    'Shuffle','every-epoch');

%% Start actual training of the NN. GPU usage for training is recommended.
%  For this option to be available, the Parallel-Computing-Toolbox has to be
%  installed in MATLAB

CP_CNN = trainNetwork(XTrain,YTrain,layers,options);

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