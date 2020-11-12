function MonteCarlo14 = assemble_MC(CP_CNN)

MonteCarlo14 =  layerGraph();
DropParam = CP_CNN.Layers(2).Probability;

tempLayers = [
    CP_CNN.Layers(1)
    predictiondropoutLayer(DropParam,"dropout_1")
    CP_CNN.Layers(3)
    CP_CNN.Layers(4)
    CP_CNN.Layers(5)
    CP_CNN.Layers(6)];
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = [
    CP_CNN.Layers(7)
    predictiondropoutLayer(DropParam,"dropout_2")
    CP_CNN.Layers(9)
    CP_CNN.Layers(10)
    CP_CNN.Layers(11)
    predictiondropoutLayer(DropParam,"dropout_3")
    CP_CNN.Layers(13)
    ];
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = additionLayer(2,"Name","addition_1");
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = [
    CP_CNN.Layers(15)
    CP_CNN.Layers(16)
    predictiondropoutLayer(DropParam,"dropout_4")
    CP_CNN.Layers(18)
    CP_CNN.Layers(19)
    CP_CNN.Layers(20)
    predictiondropoutLayer(DropParam,"dropout_5")
    CP_CNN.Layers(22)
    ];
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = additionLayer(2,"Name","addition_2");
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = [
    CP_CNN.Layers(24)
    CP_CNN.Layers(25)
    predictiondropoutLayer(DropParam,"dropout_6")
    CP_CNN.Layers(27)
    CP_CNN.Layers(28)
    CP_CNN.Layers(29)
    predictiondropoutLayer(DropParam,"dropout_7")
    CP_CNN.Layers(31)
    ];
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = [
    predictiondropoutLayer(DropParam,"dropout_8")
    CP_CNN.Layers(33)
    ];
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = additionLayer(2,"Name","addition_3");
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = [
    CP_CNN.Layers(35)
    CP_CNN.Layers(36)
    predictiondropoutLayer(DropParam,"dropout_9")
    CP_CNN.Layers(38)
    CP_CNN.Layers(39)
    CP_CNN.Layers(40)
    predictiondropoutLayer(DropParam,"dropout_10")
    CP_CNN.Layers(42)
    ];
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = [
    predictiondropoutLayer(DropParam,"dropout_11")
    CP_CNN.Layers(44)
    ];
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = additionLayer(2,"Name","addition_4");
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = [
    CP_CNN.Layers(46)
    CP_CNN.Layers(47)
    predictiondropoutLayer(DropParam,"dropout_12")
    CP_CNN.Layers(49)
    CP_CNN.Layers(50)
    CP_CNN.Layers(51)
    predictiondropoutLayer(DropParam,"dropout_14")
    CP_CNN.Layers(53)
    ];
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = [
    predictiondropoutLayer(DropParam,"dropout_13")
    CP_CNN.Layers(55)
    ];
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = additionLayer(2,"Name","addition_5");
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = [
    predictiondropoutLayer(DropParam,"dropout_16")
    CP_CNN.Layers(58)];
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = [
    CP_CNN.Layers(59)
    CP_CNN.Layers(60)
    predictiondropoutLayer(DropParam,"dropout_15")
    CP_CNN.Layers(62)
    CP_CNN.Layers(63)
    CP_CNN.Layers(64)
    predictiondropoutLayer(DropParam,"dropout_17")
    CP_CNN.Layers(66)
    ];
MonteCarlo14 = addLayers(MonteCarlo14,tempLayers);

tempLayers = [
    additionLayer(2,"Name","addition_6")
    CP_CNN.Layers(68)
    predictiondropoutLayer(DropParam,"dropout_18")
    CP_CNN.Layers(70)
    CP_CNN.Layers(71)];
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

MonteCarlo14 = assembleNetwork(MonteCarlo14);
end