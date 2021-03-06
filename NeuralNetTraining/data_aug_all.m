function OutCell = data_aug_all(InCell,DataAugMult,NInsertedPoints, RunMode)
% OutCell = data_aug_all(InCell,DataAugMult,NInsertedPoints)

if nargin < 4
    RunMode = 'all';
end

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


% Random Zoom
if isequal(lower(RunMode),'all') || isequal(lower(RunMode(1:4)),'zoom')
    for i = (NIn+1):(NIn*DataAugMult)
        % max CutPercent restricted to 90% for stability reasons
        CutPercent = rand(1,OutCell{i}.NCurves)*0.9;
        jRange = find(OutCell{i}.SelectedCurves)';
        for j = jRange
            NNCPts = length(OutCell{i}.HHApp{j}(OutCell{i}.HHApp{j}<OutCell{i}.Man_CP(j,1)));
            OutCell{i}.HHApp{j}(1:floor(CutPercent(j)*NNCPts)) = [];
            OutCell{i}.BasedApp{j}(1:floor(CutPercent(j)*NNCPts)) = [];
        end
        clear CutPercent
    end
end

% Random gaussian noise added + inserting additional points
if isequal(lower(RunMode),'all') || isequal(lower(RunMode(1:5)),'gauss')
    for i = (NIn+1):(NIn*DataAugMult)
        jRange = find(OutCell{i}.SelectedCurves)';
        for j = jRange
            % insert NInsertedPoints by linear interpolation between actual
            % curve points.
            N = randi([0 NInsertedPoints]);
            L = length(OutCell{i}.HHApp{j})*(N+1) - N;
            NoiseScale = range(OutCell{i}.BasedApp{j})/(100+randi([-40 40]));
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

% add colored noise
if isequal(lower(RunMode),'all') || isequal(lower(RunMode),'colored')
    FNLength = 1000000;
    cn = dsp.ColoredNoise('Color','brown','SamplesPerFrame',FNLength);
    FullColoredNoise = cn();
    for i = (NIn+1):(NIn*DataAugMult)
        jRange = find(OutCell{i}.SelectedCurves)';
        for j = jRange
            HeavyNoiseFlag = randi([0 1]);
            if HeavyNoiseFlag == 0
                continue
            end
            L = length(OutCell{i}.HHApp{j}(OutCell{i}.HHApp{j}<OutCell{i}.Man_CP(j,1)));
            NoiseScale = range(OutCell{i}.BasedApp{j})/(700+randi([-100 1000]));
            TanhDist = tanh(5*((L-1):-1:0)/L);
            Gauss = normpdf(((L-1):-1:0),randi(L),L/randi([1 10]));
            NoiseDist = normalize(TanhDist.*Gauss,'range');
            SampleRate = randi(20);
            NoiseWindow = randi(FNLength-SampleRate*L);
            ColoredNoise = FullColoredNoise(NoiseWindow:SampleRate:NoiseWindow+(SampleRate*L-1)).*NoiseScale.*NoiseDist';
            OutCell{i}.BasedApp{j}(1:L) = OutCell{i}.BasedApp{j}(1:L) + ColoredNoise;
        end
    end
end

end