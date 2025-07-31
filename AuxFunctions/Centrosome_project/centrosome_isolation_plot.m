total_prot_mg_ml = [2.673692946
0.034688797
0.077842324
0.047966805
0.02473029
0.112697095
0.052946058
0.069543568
0.081161826
0.099419087
0.122655602
0.24879668
0.79153527
1.438838174
1.853775934
1.910207469
2.046307054
2.026390041
1.97659751
2.006473029
2.218921162
2.029709544
2.092780083
1.923485477
2.023070539
2.182406639]; 

centrosome_ml = [6.72E+07
5.97E+06
2.99E+06
3.29E+07
2.33E+08
2.63E+08
1.87E+08
3.73E+07
3.14E+07
2.84E+07
1.19E+07
7.47E+06
1.49E+06
0.00E+00
0.00E+00
0.00E+00
0.00E+00
0.00E+00
0.00E+00
0.00E+00
0.00E+00
0.00E+00
0.00E+00
0.00E+00
0.00E+00
0.00E+00]; 

fractions = 0:1:25; 

figure(); hold on

% Set axis properties
set(gca, 'FontSize', 19, 'LineWidth', 1.5);
xlim([0 16]);
set(gca, 'xtick', 0:1:16);

% Left Y-axis (Centrosomes/ml)
yyaxis left 
p1 = plot(fractions, centrosome_ml, 'LineWidth', 3, ...
          'MarkerSize', 8, 'MarkerFaceColor', [0 0.447 0.741]);
ylabel('Centrosomes/ml', 'FontSize', 19, 'FontWeight', 'bold');
set(gca, 'YColor', [0 0.447 0.741]); % Match y-axis color to plot

% Right Y-axis (Protein concentration)
yyaxis right
p2 = plot(fractions, total_prot_mg_ml, 'LineWidth', 3, ...
          'MarkerSize', 8, 'MarkerFaceColor', [0.85 0.325 0.098]);
ylabel('Total protein [mg/ml]', 'FontSize', 19, 'FontWeight', 'bold');
set(gca, 'YColor', [0.85 0.325 0.098]); % Match y-axis color to plot

% X-axis label
xlabel('Fraction Number', 'FontSize', 19, 'FontWeight', 'bold');

% Enrichment factor line
factor = 374.7955257;
xline(4, '--', 'LineWidth', 2, 'Color', [0.5 0 0], ...
      'Alpha', 0.7, 'Label', sprintf('Enrichment factor = %.1f', factor), ...
      'FontSize', 16, 'LabelOrientation', 'horizontal');
legend([p1 p2], {'Centrosomes', 'Total protein'}, ...
       'FontSize', 16, 'Box', 'off');
xlim([0 16]);
set(gca, 'XTick', 0:1:16); % Keep x-ticks at the bottom

yyaxis left 
hold on
% tick_positions = [0, 0.5e8, 1e8, 1.5e8, 2e8, 2.5e8, 3e8];  % Adjust to your data
% tick_labels = {'0', '0.5', '1', '1.5', '2', '2.5', '3'};  % Unicode superscripts

% Apply custom ticks
% yticks(tick_positions);
% yticklabels(tick_labels);
