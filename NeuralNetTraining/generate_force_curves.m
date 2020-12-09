function [X,Y,DummyForceMap] = generate_force_curves(NCurves)
% [X,Y,DummyForceMap] = generate_force_curves()
%
% create synthetic force curves for CNN training for deep CP estimates

% create ground truth force curve

% idealised material indentation response function with baseline
syms IndentResp(x,p1,p2)
IndentResp(x,p1,p2) = piecewise(x<0,0,x>=0,p2*x^(p1));

DummyForceMap = ForceMap('I dont care','leave britney alone','Dummy',NCurves);

if NCurves <= 36
    tiledlayout(floor(sqrt(NCurves)+1),ceil(sqrt(NCurves)),'Padding','none')
else
    tiledlayout(6,6,'Padding','none')
end

for i=1:NCurves
    NPoints = floor(1800*rand+200);
    Base2ResponseRatio = 0.7*rand+0.3;
    xVal = [0:1/(NPoints-1):1] - Base2ResponseRatio;
    pr1 = 2*rand+1;
    pr2 = 10000*rand;
    if i <= 36
        nexttile
    end
    DummyForceMap.BasedApp{i} = xVal;
    DummyForceMap.App{i} = IndentResp(xVal,pr1,pr2);
    plot(xVal,DummyForceMap.App{i})
end

DummyForceMap.Man_CP = zeros()
X = 1;
Y = 2;

end