function [X,Y,DummyForceMap] = generate_force_curves(NCurves)
% [X,Y,DummyForceMap] = generate_force_curves()
%
% create synthetic force curves for CNN training for deep CP estimates

% create ground truth force curve

% idealised material indentation response function with baseline
syms IndentResp(x,p1,p2)
IndentResp(x,p1,p2) = piecewise(x<0,0,x>=0,p2*x^(p1));

DummyForceMap = ForceMap('I dont care','leave britney alone','ID not needed','Dummy',NCurves);


for i=1:NCurves
    NPoints = floor(1800*rand+200);
    Base2ResponseRatio = 0.7*rand+0.3;
    xVal = [0:1/(NPoints-1):1] - Base2ResponseRatio;
    pr1 = 2*rand+1;
    pr2 = 10000*rand;
    DummyForceMap.HHApp{i} = xVal;
    DummyForceMap.BasedApp{i} = double(IndentResp(xVal,pr1,pr2));
    DummyForceMap.BasedRet{i} = zeros(1,length(xVal));
    DummyForceMap.HHRet{i} = xVal;
end

DummyForceMap.Man_CP = zeros(NCurves,2);
DummyForceMap.CPFlag.Manual = 1;
X = 1;
Y = 2;

% add noise to onset of the the force response curve using colored noise
% with a randomly varied gaussian as an envelope
FNLength = 1000000;
cn = dsp.ColoredNoise('Color','brown','SamplesPerFrame',FNLength);
FullColoredNoise = cn();

jRange = find(DummyForceMap.SelectedCurves)';
for j = jRange
    HeavyNoiseFlag = randi([0 1]);
    if HeavyNoiseFlag == 0
        continue
    end
    L = length(DummyForceMap.HHApp{j});
    NoiseScale = range(DummyForceMap.BasedApp{j})/(700+randi([-100 1000]));
    %         TanhDist = tanh(5*((L-1):-1:0)/L);
    Gauss = normpdf(1:L,length(DummyForceMap.HHApp{j}(DummyForceMap.HHApp{j} < 0)),L/randi([3 10]));
    NoiseDist = normalize(Gauss,'range');
    SampleRate = randi(20);
    NoiseWindow = randi(FNLength-SampleRate*L);
    ColoredNoise = FullColoredNoise(NoiseWindow:SampleRate:NoiseWindow+(SampleRate*L-1));
    ColoredNoise = abs(ColoredNoise);
    ColoredNoise = smoothdata(ColoredNoise,'movmedian');
    ColoredNoise = smoothdata(ColoredNoise,'movmean');
    ColoredNoise = ColoredNoise.*NoiseScale.*NoiseDist';
    DummyForceMap.BasedApp{j} = DummyForceMap.BasedApp{j} + ColoredNoise';
    plot(DummyForceMap.HHApp{j},DummyForceMap.BasedApp{j})
    hold on
    plot(DummyForceMap.HHApp{j},NoiseScale.*NoiseDist)
    plot(DummyForceMap.HHApp{j},ColoredNoise)
    hold off
    pause(1)
    clear L NoiseScale NoiseDist Gauss SampleRate NoiseWindow ColoredNoise
end

if NCurves <= 36
    tiledlayout(floor(sqrt(NCurves)+1),ceil(sqrt(NCurves)),'Padding','none')
else
    tiledlayout(6,6,'Padding','none')
end

for i=1:NCurves
    if i <= 36
        nexttile
        plot(DummyForceMap.HHApp{i},DummyForceMap.BasedApp{i})
    end
end

end