function [MaxPeakMap,MinPeakMap,MaxPeakValueMap,DensityMap,minCoords, maxCoords,...
    subsetsGrid, boundariesGrid, subsetsCloud, boundariesCloud] =...
    PC2HM_point_cloud_to_height_map(PC, Resolution, zResolution,...
    MaxPointsPerGP, PartitionShrinkingFactor, InterpolationExpansionFactor,...
    Lambda, Sigma, Noise)

MinPeakDistance = .2/zResolution;

[NormPC,Scale,Offset] = PC2HM_normalize_pointcloud(PC,'KeepAspects');

PCWeights = ones(size(NormPC,1),1);
[GridPoints,minCoords, maxCoords] = PC2HM_createGrid(NormPC,Resolution);
GridWeights = ones(size(GridPoints,1),1).*zResolution;
MaxPeakMap = zeros(Resolution,Resolution);
MinPeakMap = MaxPeakMap;
MaxValueMap = MaxPeakMap;
NPeakMap = MaxPeakMap;

[subsetsGrid, boundariesGrid, subsetsCloud, boundariesCloud] = ...
    PC2HM_partitionPointCloud(GridPoints, GridWeights, NormPC, PCWeights,...
    MaxPointsPerGP, PartitionShrinkingFactor,InterpolationExpansionFactor);

MaxMean = [];
MaxSigma = [];
MaxPoints = [];
MinPoints = [];
MaxValuePoints = [];

h = waitbar(0,'setting up');
% FF = figure;
% Interpolate
for i=1:length(subsetsGrid)
    if isempty(subsetsCloud{i})
        continue
    end
    MinZ = min(subsetsCloud{i}(:,3));
    MaxZ = max(subsetsCloud{i}(:,3));
    if MinZ==MaxZ
        continue
    end
    Z = linspace(MinZ,MaxZ,zResolution);
    QueuePoints = [];
    for j=1:size(subsetsGrid{i},1)
        Temp = repmat(subsetsGrid{i}(j,:),[zResolution 1]);
        Temp(:,3) = Z;
        QueuePoints = [QueuePoints ; Temp];
    end
    if isempty(QueuePoints)
        continue
    end
    F = ones(size(subsetsCloud{i},1),1);
    
%     PC2HM_visualizePartition(subsetsCloud(i),boundariesGrid(i),'Figure',FF)
%     drawnow
    
    [InterpolatedMean,InterpolatedSigma] = PC2HM_predictGP_mean_model_zero(...
        subsetsCloud{i},QueuePoints,Sigma,Lambda,F,Noise,'CalculateSigma',false);
    
    DeflattenedIM = reshape(InterpolatedMean,zResolution,[]);
%     DeflattenedIS = reshape(InterpolatedSigma,size(subsetsGrid{i},1),[]);
    
    [TempMaxMean,MaxIdx] = max(DeflattenedIM',[],2);
    MaxMean = [MaxMean ; TempMaxMean];
    [MaxIdx, MinIdx, MaxValue, NPeaks] = PC2HM_findPeakVectors(DeflattenedIM,...
        'Threshold',0,'MinPeakDistance',MinPeakDistance,'MinPeakHeight',0);
%     TempMaxSigma = DeflattenedIS;
%     MaxSigma = [MaxSigma ; TempMaxSigma];
    TempMaxPoints = subsetsGrid{i};
    TempMaxPoints(:,3) = Z(MaxIdx);
    MaxPoints = [MaxPoints ; TempMaxPoints];
    TempMinPoints = subsetsGrid{i};
    TempMinPoints(:,3) = Z(MinIdx);
    MinPoints = [MinPoints ; TempMinPoints];
    TempMaxValuePoints = subsetsGrid{i};
    TempMaxValuePoints(:,3) = Z(MaxValue);
    MaxValuePoints = [MaxValuePoints ; TempMaxValuePoints];
    waitbar(i/length(subsetsGrid),h,sprintf('%i out of %i subsets processed',i,length(subsetsGrid)))
end

waitbar(length(subsetsGrid)/length(subsetsGrid),h,...
    sprintf('%i out of %i subsets processed',length(subsetsGrid),length(subsetsGrid)))

MaxPoints = PC2HM_backtransform_pointcloud(MaxPoints,Scale,Offset);
MinPoints = PC2HM_backtransform_pointcloud(MinPoints,Scale,Offset);
MaxValuePoints = PC2HM_backtransform_pointcloud(MaxValuePoints,Scale,Offset);

minCoords = PC2HM_backtransform_pointcloud(minCoords,Scale,Offset);
maxCoords = PC2HM_backtransform_pointcloud(maxCoords,Scale,Offset);


MaxPeakMap = PC2HM_pointcloud2grid(MaxPoints,Resolution);
MinPeakMap = PC2HM_pointcloud2grid(MinPoints,Resolution);
MaxPeakValueMap = PC2HM_pointcloud2grid(MaxValuePoints,Resolution);
MaxDensity = MaxPoints;
MaxDensity(:,3) = MaxMean;
DensityMap = PC2HM_pointcloud2grid(MaxDensity,Resolution);

% PC2HM_visualizePartition(subsetsCloud,boundariesGrid)
% figure
% PC2HM_visualizePartition(subsetsCloud,boundariesCloud)
% figure
% PC2HM_visualizePartition(subsetsGrid,boundariesGrid)

close(h)

end