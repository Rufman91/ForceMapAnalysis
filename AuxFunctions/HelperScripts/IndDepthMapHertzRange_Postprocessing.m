%% Fill these in correctly
UpperForceCutoff = 1;
LowerForceCutoff = 1;
UpperCurveFraction = 1;
LowerCurveFraction = 1;

%%
for k=1:E.NumForceCurves
    Hertzfit = obj.HertzFitValues;
    iRange = find(E.FM{k}.SelectedCurves);
    E.FM{k}.IndentationDepthHertzFitRange = zeros(E.FM{k}.NCurves,1);
    MaxFitRange = zeros(E.FM{k}.NCurves,1);
    for i=iRange'
        [App,HHApp] = E.FM{k}.get_force_curve_data(iRange(i),'AppRetSwitch',0,...
            'BaselineCorrection',1,'TipHeightCorrection',0,...
            'Sensitivity','corrected',...
            'Unit','N');
        CP = E.FM{k}.CP_Old(i,:);
        if SortHeightDataForFit
            HHApp = sort(HHApp,'ascend');
        end
        
        force = App - CP;
        tip_h = (HHApp - CP(1,1)) - force/E.FM{k}.SpringConstant;
        force(1:(length(force)-length(tip_h))) = [];
        if length(tip_h) < 2
            continue
        end
        Max = max(tip_h);
        % Apply absolute force cut offs
        if ~isempty(UpperForceCutoff)
            tip_h(force>UpperForceCutoff) = [];
            force(force>UpperForceCutoff) = [];
        end
        if ~isempty(LowerForceCutoff)
            tip_h(force<LowerForceCutoff) = [];
            force(force<LowerForceCutoff) = [];
        end
        % delete everything below curve_percent of the maximum
        % force
        force(force<(LowerCurveFraction)*max(force)) = [];
        force(force>(UpperCurveFraction)*max(force)) = [];
        tip_h(1:(length(tip_h)-length(force))) = [];
        % Allocate the maximum indentation for later
        % calculation of max indent. depth of fitrange
        MaxFR = max(tip_h);
        MaxFitRange(i) = MaxFR(1);
    end
    
    for i=1:E.FM{k}.NCurves
        E.FM{k}.IndentationDepthHertzFitRange(iRange(i)) = MaxFitRange(i)+Hertzfit{i}.b;
    end
    for i=1:E.FM{k}.NumPixelsX
        for j=1:E.FM{k}.NumPixelsY
            IndDepMapHertzFitRange(i,j) = E.FM{k}.IndentationDepthHertzFitRange(E.FM{k}.Map2List(i,j));
        end
    end
    Channel = E.FM{k}.create_standard_channel(IndDepMapHertzFitRange,'Indentation Depth Hertz Fit Range','m');
    E.FM{k}.add_channel(Channel,~KeepOldResults)
end