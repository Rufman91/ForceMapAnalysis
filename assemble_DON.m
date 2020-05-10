function DropoutNet = assemble_DON(CP_CNN, DropParam)

DropoutNet = layerGraph();

tempLayers = [
    CP_CNN.Layers(1)
    predictiondropoutLayer(DropParam,'preddrop_1')
    CP_CNN.Layers(2)
    CP_CNN.Layers(3)
    reluLayer("Name","relu_9")
    CP_CNN.Layers(5)];
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = [
    CP_CNN.Layers(6)
    reluLayer("Name","relu_1")
    predictiondropoutLayer(DropParam,'preddrop_2')
    CP_CNN.Layers(8)
    CP_CNN.Layers(9)
    reluLayer("Name","relu_2")
    predictiondropoutLayer(DropParam,'preddrop_3')
    CP_CNN.Layers(11)];
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = additionLayer(2,"Name","addition_1");
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = [
    CP_CNN.Layers(13)
    reluLayer("Name","relu_3")
    predictiondropoutLayer(DropParam,'preddrop_4')
    CP_CNN.Layers(15)
    CP_CNN.Layers(16)
    reluLayer("Name","relu_4")
    predictiondropoutLayer(DropParam,'preddrop_5')
    CP_CNN.Layers(18)];
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = additionLayer(2,"Name","addition_2");
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = [
    predictiondropoutLayer(DropParam,'preddrop_?1')
    CP_CNN.Layers(26)];
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = [
    CP_CNN.Layers(20)
    reluLayer("Name","relu_5")
    predictiondropoutLayer(DropParam,'preddrop_6')
    CP_CNN.Layers(22)
    CP_CNN.Layers(23)
    reluLayer("Name","relu_6")
    predictiondropoutLayer(DropParam,'preddrop_7')
    CP_CNN.Layers(25)];
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = additionLayer(2,"Name","addition_3");
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = [
    predictiondropoutLayer(DropParam,'preddrop_?2')
    CP_CNN.Layers(28)];
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = [
    CP_CNN.Layers(29)
    reluLayer("Name","relu_7")
    predictiondropoutLayer(DropParam,'preddrop_8')
    CP_CNN.Layers(31)
    CP_CNN.Layers(32)
    reluLayer("Name","relu_8")
    predictiondropoutLayer(DropParam,'preddrop_9')
    CP_CNN.Layers(34)];
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = additionLayer(2,"Name","addition_4");
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = [
    predictiondropoutLayer(DropParam,'preddrop_?3')
    CP_CNN.Layers(42)];
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = [
    CP_CNN.Layers(36)
    reluLayer("Name","relu_10")
    predictiondropoutLayer(DropParam,'preddrop_10')
    CP_CNN.Layers(38)
    CP_CNN.Layers(39)
    reluLayer("Name","relu_11")
    predictiondropoutLayer(DropParam,'preddrop_11')
    CP_CNN.Layers(41)];
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = additionLayer(2,"Name","addition_5");
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = [
    CP_CNN.Layers(44)
    reluLayer("Name","relu_12")
    predictiondropoutLayer(DropParam,'preddrop_12')
    CP_CNN.Layers(46)
    CP_CNN.Layers(47)
    reluLayer("Name","relu_13")
    predictiondropoutLayer(DropParam,'preddrop_13')
    CP_CNN.Layers(49)];
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = [
    predictiondropoutLayer(DropParam,'preddrop_?4')
    CP_CNN.Layers(50)];
DropoutNet = addLayers(DropoutNet,tempLayers);

tempLayers = [
    additionLayer(2,"Name","addition_6")
    predictiondropoutLayer(DropParam,'preddrop_?5')
    globalAveragePooling2dLayer("Name","gapool")
    CP_CNN.Layers(53)
    CP_CNN.Layers(54)];
DropoutNet = addLayers(DropoutNet,tempLayers);

% clean up helper variable
clear tempLayers;

DropoutNet = connectLayers(DropoutNet,"maxpool","batchnorm_1");
DropoutNet = connectLayers(DropoutNet,"maxpool","addition_1/in2");
DropoutNet = connectLayers(DropoutNet,"conv_2","addition_1/in1");
DropoutNet = connectLayers(DropoutNet,"addition_1","batchnorm_3");
DropoutNet = connectLayers(DropoutNet,"addition_1","addition_2/in2");
DropoutNet = connectLayers(DropoutNet,"conv_4","addition_2/in1");
DropoutNet = connectLayers(DropoutNet,"addition_2","preddrop_?1");
DropoutNet = connectLayers(DropoutNet,"addition_2","batchnorm_5");
DropoutNet = connectLayers(DropoutNet,"conv_9","addition_3/in2");
DropoutNet = connectLayers(DropoutNet,"conv_6","addition_3/in1");
DropoutNet = connectLayers(DropoutNet,"addition_3","preddrop_?2");
DropoutNet = connectLayers(DropoutNet,"addition_3","batchnorm_7");
DropoutNet = connectLayers(DropoutNet,"conv_10","addition_4/in2");
DropoutNet = connectLayers(DropoutNet,"conv_8","addition_4/in1");
DropoutNet = connectLayers(DropoutNet,"addition_4","preddrop_?3");
DropoutNet = connectLayers(DropoutNet,"addition_4","batchnorm_10");
DropoutNet = connectLayers(DropoutNet,"conv_16","addition_5/in2");
DropoutNet = connectLayers(DropoutNet,"conv_13","addition_5/in1");
DropoutNet = connectLayers(DropoutNet,"addition_5","batchnorm_12");
DropoutNet = connectLayers(DropoutNet,"addition_5","preddrop_?4");
DropoutNet = connectLayers(DropoutNet,"conv_17","addition_6/in2");
DropoutNet = connectLayers(DropoutNet,"conv_15","addition_6/in1");

DropoutNet = assembleNetwork(DropoutNet);
save('DropoutNet.mat','DropoutNet')
end