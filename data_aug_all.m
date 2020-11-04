function OutCell = data_aug_all(InCell,DataAugMult,NInsertedPoints)
% OutCell = data_aug_all(InCell,DataAugMult,NInsertedPoints)

DataAugMult = ceil(abs(DataAugMult));
NInsertedPoints = floor(abs(NInsertedPoints));
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

for i = (NIn+1):(NIn*DataAugMult)
    for j = 1:OutCell{i}.NCurves
        % insert NInsertedPoints by linear interpolation between actual
        % curve points.
        N = randi([0 NInsertedPoints]);
        L = length(OutCell{i}.HHApp{j})*(N+1) - N;
        NoiseScale = range(OutCell{i}.BasedApp{j})/(100+randi([-20 20]));
        if N ~= 0
            NewHH = zeros(L,1);
            NewHH(1:N+1:end) = OutCell{i}.HHApp{j};
            NewBased = zeros(L,1);
            NewBased(1:N+1:end) = OutCell{i}.BasedApp{j};
            for k = 1:N
                NewHH(k+1:N+1:end+k-(N+1)) = NewHH(1:N+1:end-(N+1))+(NewHH(1+N+1:N+1:end) - NewHH(1:N+1:end-(N+1)))*(k/(N+2));
                NewBased(k+1:N+1:end+k-(N+1)) = NewBased(1:N+1:end-(N+1))+(NewBased(1+N+1:N+1:end) - NewBased(1:N+1:end-(N+1)))*(k/(N+2));
            end
            OutCell{i}.HHApp{j} = NewHH;
            OutCell{i}.BasedApp{j} = NewBased;
        end
        clear NewHH NewBased
        
        Noise = randn(L,1)*NoiseScale;
        OutCell{i}.BasedApp{j} = OutCell{i}.BasedApp{j} + Noise;
    end
end

end