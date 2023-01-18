l = logspace(-5,1,5);
s = logspace(-3,3,5);
n = logspace(-2,3,5);

[L,S,N] = meshgrid(l,s,n);
L = L(:);
S = S(:);
N = N(:);

for i=1:length(L)
Lambda = L(i);
Sigma = S(i);
Noise = N(i);
Resolution = 1024;
zResolution = 100;
MaxPointsPerGP = 10000;
PartitionShrinkingFactor = 2;
InterpolationExpansionFactor = 10;
[HeightMap,DensityMap,UncertaintyMap] = PC2HM_point_cloud_to_height_map(PC, Resolution, zResolution, MaxPointsPerGP, PartitionShrinkingFactor, InterpolationExpansionFactor,...
Lambda, Sigma, Noise);
FileName = ['ParamSweep_L' strrep(num2str(L(i)),'.','c') '_S' strrep(num2str(S(i)),'.','c') '_N' strrep(num2str(N(i)),'.','c')];
save([FileName '.mat'],'DensityMap','HeightMap','Lambda','Sigma','Noise','Resolution','zResolution','MaxPointsPerGP','InterpolationExpansionFactor','PartitionShrinkingFactor')
end


