
figure('name', 'Height dependence'); hold on
box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
CsEModHertz_mean_l = nan(1,total_cs);
CsHeight_mean_l = nan(1,total_cs);
for i = 1:length(sample_idx)
    for j = 1:counter(i)
        clr_counter = j + sum(counter(1:i-1));
        c = color_map(clr_counter,:);
        CsEModHertz_mean_l(clr_counter) = CsEModHertz_mean{i,j}.*1e-3;
        CsHeight_mean_l(clr_counter) = CsHeight_mean{i,j}.*1e9;

%         scatter(CsHeight_mean_l(clr_counter), CsEModHertz_mean_l(clr_counter), 50, [0 0.4470 0.7410], "filled");
    end
end


hold on
[xData, yData] = prepareCurveData( CsHeight_mean_l, CsEModHertz_mean_l );

% Set up fittype and options.
ft = fittype( 'power2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
% opts.StartPoint = [6133207.15536504 -2.05322611301745 -35.6361239539309]; % w/o modifications
opts.StartPoint = [2058429.67363441 -2.04285486042691 -8.67399822133641]; % thin film 

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
h = plot( fitresult, xData, yData ); legend boxoff
scatter(xData, yData, 50, [0 0.4470 0.7410], "filled"); 
ylabel('Indentation modulus [kPa]'); 
xlabel('Centrosome height [nm]'); 
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off')
xlim([0 1000])
ylim([0 300])
