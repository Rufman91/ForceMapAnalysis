function [X,Y,DummyForceMap] = generate_force_curves(NCurves)
% [X,Y,DummyForceMap] = generate_force_curves()
%
% create synthetic force curves for CNN training for deep CP estimates

% create ground truth force curve

% idealised material indentation response function with baseline
syms IndentResp(x,p1,p2)
IndentResp(x,p1,p2) = piecewise(x<0,0,x>=0,p2*x^(p1));

for i=1:NCurves
    NPoints = floor(1800*rand+200);
    Base2ResponseRatio = rand;
    xVal = [0:1/(NPoints-1):1] - Base2ResponseRatio;
    pr1 = 2*rand+1;
    pr2 = 10000*rand;
end
plot(xVal,IndentResp(xVal,pr1,pr2))

end