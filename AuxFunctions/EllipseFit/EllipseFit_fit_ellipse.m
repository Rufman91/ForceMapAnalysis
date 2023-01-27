function EllFit =  EllipseFit_fit_ellipse(X,Y,Axis)
RangeX = range(X);
RangeY = range(Y);
CompoundScale = sqrt(RangeX^2 + RangeY^2);
XScaled = X./CompoundScale;
YScaled = Y./CompoundScale;

YScaled = YScaled-max(YScaled)/2;
YMirrored = -YScaled;

YCombo = [YScaled ; YMirrored];
XCombo = [XScaled ; XScaled];

if nargin > 2
    cla(Axis)
    EllFit = fit_ellipse(XCombo,YCombo,Axis);
    hold on
    plot(XCombo,YCombo,'bO');
    axis equal
    hold off
else
    EllFit = fit_ellipse(XCombo,YCombo);
end



end

% b = max(YScaled)/2;
% 
% % Define the elliptic equation
% FitFunction = ['sqrt(1 - (x-c_x).^2./a^2).*' num2str(b) ' + ' num2str(b)];
% 
% s = fitoptions('Method','NonlinearLeastSquares',...
%     'Lower',[range(XScaled)/2 -range(XScaled)/2],...
%     'Upper',[inf range(XScaled)/2],...
%     'MaxIter',400,...
%     'StartPoint',[range(XScaled) 0]);
% f = fittype(FitFunction,'options',s);
% try
%     [EllFit,GoF] = fit(XScaled,...
%         YScaled,f);
%     a = EllFit.a*CompoundScale;
%     b = EllFit.b*CompoundScale;
%     c_x = EllFit.c_x*CompoundScale;
% catch ME
%     a = nan;
%     b = nan;
%     c_x = nan;
%     GoF = nan;
% end
