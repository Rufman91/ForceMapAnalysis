function OutCell = data_aug_zoom(InCell,DataAugMult)
% OutCell = data_aug_zoom(InCell,DataAugMult)

DataAugMult = ceil(abs(DataAugMult));
NIn = length(InCell);
OutCell = cell(1,NIn*DataAugMult);
for i=1:DataAugMult
    for j = 1:NIn
        OutCell{(i-1)*NIn+j} = InCell{j}.copy;
    end
end
clear InCell

for i = (NIn+1):(NIn*DataAugMult)
    % max CutPercent restricted to 95% due to stability reasons
    CutPercent = rand(1,OutCell{i}.NCurves)*0.95;
    for j = 1:OutCell{i}.NCurves
        NNCPts = length(OutCell{i}.HHApp{j}(OutCell{i}.HHApp{j}<OutCell{i}.Man_CP(j,1)));
        OutCell{i}.HHApp{j}(1:floor(CutPercent(j)*NNCPts)) = [];
        OutCell{i}.BasedApp{j}(1:floor(CutPercent(j)*NNCPts)) = [];
    end
    clear CutPercent
end


end