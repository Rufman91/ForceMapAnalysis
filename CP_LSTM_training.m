%% Load maps as ForceMap object into 1xNmaps cell.
% Prepare and manually select and label every single force curve


prompt = 'Enter how many force maps you want to load for training:';
dlgtitle = 'Number of training maps';
dims = [1 35]; 
definput = {'2'};
Nmaps = str2double(inputdlg(prompt,dlgtitle,dims,definput));

% Load force maps, choose, basefit, and manually CP-label all the curves.
% Save after every ForceMap to keep progress in case of
% crash/error/boredom-related-scriptstop.
for i=1:Nmaps
    FM{i} = ForceMap();
    %FM{i}.choose_curves();
    FM{i}.base_and_tilt();
    FM{i}.manual_CP();
    quicksave = sprintf('%i_out_of_%i_ForceMaps',i,Nmaps);
    save(quicksave);
end

%% Prepare training data
k = 1;
for i =1:length(FM)
    jRange = find(FM{i}.selected_curves);
    for j=jRange'
        LSTMX(k,1) = {[(FM{i}.basedapp{j}-min(FM{i}.basedapp{j}))/range(FM{i}.basedapp{j}),...
            (FM{i}.thapp{j}-min(FM{i}.thapp{j}))/range(FM{i}.thapp{j})]'};
        LSTMY(k,:) = [(FM{i}.Man_CP(j,1)-min(FM{i}.thapp{j}))/range(FM{i}.thapp{j}),(FM{i}.Man_CP(j,2)-min(FM{i}.basedapp{j}))/range(FM{i}.basedapp{j})];
        k = k + 1;
    end
end
       
Val_idx = randperm(size(LSTMX,1),floor(size(LSTMX,1)/4));
LSTMXTrain = LSTMX;
LSTMYTrain = LSTMY;
LSTMXValidation = LSTMX(Val_idx,1);
LSTMXTrain(Val_idx,:) = [];
LSTMYValidation = LSTMY(Val_idx,:);
LSTMYTrain(Val_idx,:) = [];

%% Define Layers and traning options

layers = [
    sequenceInputLayer(2,"Name","sequence")
    bilstmLayer(256,"Name","bilstm","OutputMode","last")
    fullyConnectedLayer(2,"Name","fc")
    regressionLayer("Name","regressionoutput")];

options = trainingOptions('adam','Plots','training-progress',...
    'ValidationData',{LSTMXValidation,LSTMYValidation},...
    'ValidationFrequency',20,...
    'MiniBatchSize',50,...
    'MaxEpochs',100,...
    'InitialLearnRate',0.002,...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropPeriod',10,...
    'LearnRateDropFactor',0.5,...
    'Shuffle','every-epoch');

%% Start actual training of the NN. GPU usage for training is recommended.
%  For this option to be available, the Parallel-Computing-Toolbox has to be
%  installed in MATLAB

CP_LSTM = trainNetwork(LSTMXTrain,LSTMYTrain,layers,options);
