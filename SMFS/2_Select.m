close all

%% Allocate data

xxA0=h_approach;
xxA1=h_retraction;

yyA0=f_approach;
yyA1=f_retraction;


%% Define variables

ii=1; % ii corresponds to the number of subplots

jj=1; % jj corresponds to the number of force curves 

kk=0;


%% Define the margins of the subtightplot 
    %%% "subtightplot.m" is needed
    subplot = @(m,n,p) subtightplot (m, n, p,...
        [0.01 0.005],...% gap [vertical horizontal]
        [0.01 0.01],...% margin height [lower upper]
        [0.005 0.005]); % margin width [left right]

%% Figure
for ii=1:4  
    fig=figure(ii);
    h_fig1=gcf; % Defines the handle for the current figure 

    h_fig1.Color='white'; % changes the background color of the figure
    h_fig1.Units='normalized'; % Defines the units 
    h_fig1.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
    h_fig1.PaperOrientation='landscape';
    % h_fig1.WindowState='fullscreen'; % Shows the figure in fullscreen mode

    %%% Axes 
    %h_axes1=gca; % Defines the handle for the current axes
    box on
    h_axes1=axes;
    h_axes1.FontSize = 20;
    h_axes1.XLabel.String = 'Tip-sample seperation  (nm)';
    h_axes1.XLabel.FontSize = 20;
    h_axes1.YLabel.String = 'Force (nN)';
    h_axes1.YLabel.FontSize = 20;

    %%% Plots
    % Plot the subplot
    for jj= (1:25)+25*kk             
        subplot(5,5,jj-25*kk)
        hold on
        grid on
        plot(xxA0(:,jj),yyA0(:,jj)); % Plots all row entries in column jj
        plot(xxA1(:,jj),yyA1(:,jj)); % Plots all row entries in column jj
  
        %%% Axes
        box on
        h_axes1=axes;
        h_axes1.FontSize = 20;
        h_axes1.XLabel.String = 'Tip-sample seperation  (nm)';
        h_axes1.XLabel.FontSize = 20;
        h_axes1.YLabel.String = 'Force (nN)';
        h_axes1.YLabel.FontSize = 20;
    end
    
    %% Dialog boxes
        %%% "bttnChoiseDialog.m" is needed
        inputOptions={'Select all', 'Select none', 'Select all - except of', 'Select none - except of'}; % Define the input arguments
        % 'Select all' = 1
        % 'Select none' = 2
        % 'Select all - except of' = 3
        % 'Select none - except of' = 4

        defSelection=inputOptions{1}; % Default selection; Defined selection if the window is closed without choosing a selection possibility

        SelectionData(k).subplotselect{j,1}=bttnChoiseDialog(inputOptions, 'Force curve selection', defSelection,...
        'Please choose the appropriate button ...'); % Selected button number
        
        % Case 3:
        if isequal(SelectionData(k).subplotselect{j,1},3)
            prompt = {'Enter the force curve number you do not want for further analysis'};
            definput = {''};
            opts.Interpreter = 'tex';
            SelectionData(k).numbnonselectfc(j,1) = inputdlg(prompt,'Select all except of ...',[1 100],definput,opts);        
            %%% Convert the 'except of' force curves in the SelectionData structure   
            SelectionData(k).numbnonselectfc{j,1} = str2num(SelectionData(k).numbnonselectfc{j,1}); 
                
        % Case 4:      
        elseif isequal(SelectionData(k).subplotselect{j,1},4)
            prompt = {'Enter the force curve number you do not want for further analysis'};
            definput = {''};
            opts.Interpreter = 'tex';
            SelectionData(k).numbselectfc(j,1) = inputdlg(prompt,'Select none except of ...',[1 100],definput,opts);        
            %%% Convert the 'except of' force curves in the SelectionData structure   
           SelectionData(k).numbselectfc{j,1} = str2num(SelectionData(k).numbselectfc{j,1}); 
        
        end
    
    kk=ii;
end