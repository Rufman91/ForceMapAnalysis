%% 
clear
close all
clc

%% Define RGB colors

RGB1=[0 26 255]./255;  % Blue 
RGB2=[255 119 0]./255; % Orange
RGB7=[255 230 0]./255; % Yellow

RGB3=[80 220 100]./255; % Emerald
RGB4=[200 81 160]./255; % Compl to emerald
RGB5=[81 172 200]./255; 
RGB6=[200 108 81]./255;

%% Allocate data

x0=mapsData(1).load_data{1,1}(:,1);
y0=mapsData(1).load_data{1,1}(:,2);

x1=mapsData(1).unload_data{1,1}(:,1);
y1=mapsData(1).unload_data{1,1}(:,2);

%% Define variables

ii=1; % ii corresponds to the number of subplots

jj=1; % jj corresponds to the number of force curves 

kk=0;

%% Figure
%%
for ii=1:4  
    fig=figure(1);
    h_fig1=gcf; % Defines the handle for the current figure 

    h_fig1.Color='white'; % changes the background color of the figure
    h_fig1.Units='normalized'; % Defines the units 
    h_fig1.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
    h_fig1.PaperOrientation='landscape';

    %% Plotting the tiles
    t = tiledlayout(5,5);
    %t.TileSpacing = 'compact';
    %t.Padding = 'compact';
    t.TileSpacing = 'none'; % To reduce the spacing between the tiles
    t.Padding = 'none'; % To reduce the padding of perimeter of a tile

        for jj=1:25  
            % Tile jj
            ax=nexttile
            hold on
            grid on
            plot(x0,y0);
            plot(x1,y1);
            % Legend, x- and y-labels and title
            legend('Approach','Retraction','Location','best')
            xlabel('Tip-sample seperation  (nm)','FontSize',11,'Interpreter','latex');
            ylabel('Force (nN)','FontSize',11,'Interpreter','latex');
            title('Sample 1')
        end
        
        %% Dialog boxes
        %%% "bttnChoiseDialog.m" file is needed
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
        
        %% Colour highlighting of the force curves regarding the choosen answer and storage in a structure
           
        if isequal(SelectionData(k).subplotselect{j,1},1)
            for i = (1:25)+25*l % If all force curve from a SubMap are selected to be analyzed
            subplot(5,5,i-25*l)
            title(num2str(indent(i,1)),'Color','b') 
            SelectionData(k).selectfc_load(i,1)=mapsData(k).load_data(i,1); % Store the selected Dataset of mapsData.load in a cell
            SelectionData(k).selectfc_unload(i,1)=mapsData(k).unload_data(i,1);
            end
             
        elseif isequal(SelectionData(k).subplotselect{j,1},2)
            for i = (1:25)+25*l % If none force curve from a SubMap are selected to be analyzed
            subplot(5,5,i-25*l)
            title(num2str(indent(i,1)),'Color','r') 
            SelectionData(k).nonselectfc_load(i,1)=mapsData(k).load_data(i,1); % Store the selected Dataset of mapsData.Approach in a cell
            SelectionData(k).nonselectfc_unload(i,1)=mapsData(k).unload_data(i,1);
            end

        elseif isequal(SelectionData(k).subplotselect{j,1},3)        
            for i = (1:25)+25*l 
        %%% Step 1: Select all force curves of the SubMap
            subplot(5,5,i-25*l)
            title(num2str(indent(i,1)),'Color','b')
            SelectionData(k).selectfc_load(i,1)=mapsData(k).load_data(i,1); % Store the selected Dataset of mapsData.Approach in a cell
            SelectionData(k).selectfc_unload(i,1)=mapsData(k).unload_data(i,1);
            end   
         
        %%% Step 2: Correct for the non-selected force curves of selection 3            
            for i=SelectionData(k).numbnonselectfc{j,1} % Deselect the choosen force curves
            subplot(5,5,i-25*l)
            title(num2str(indent(i,1)),'Color','r')
            SelectionData(k).selectfc_load(i,1)={[]}; % Delete the non-selected force curves in the selectfc cell
            SelectionData(k).selectfc_unload(i,1)={[]};
            SelectionData(k).nonselectfc_load(i,1)=mapsData(k).load_data(i,1); % Store the selected Dataset of mapsData.Approach in a cell
            SelectionData(k).nonselectfc_unload(i,1)=mapsData(k).unload_data(i,1);
            end
         
        elseif isequal(SelectionData(k).subplotselect{j,1},4)        
            for i = (1:25)+25*l 
        %%% Step 1: Select all force curves of the SubMap
            subplot(5,5,i-25*l)
            title(num2str(indent(i,1)),'Color','r')
            SelectionData(k).nonselectfc_load(i,1)=mapsData(k).load_data(i,1); % Store the selected Dataset of mapsData.Approach in a cell
            SelectionData(k).nonselectfc_unload(i,1)=mapsData(k).unload_data(i,1);
            end   
         
        %%% Step 2: Correct for the non-selected force curves of selection 3            
            for i=SelectionData(k).numbselectfc{j,1} % Deselect the choosen force curves
            subplot(5,5,i-25*l)
            title(num2str(indent(i,1)),'Color','b')
            SelectionData(k).nonselectfc_load(i,1)={[]}; % Delete the selected force curves in the selectfc cell
            SelectionData(k).nonselectfc_unload(i,1)={[]};
            SelectionData(k).selectfc_load(i,1)=mapsData(k).load_data(i,1); % Store the selected Dataset of mapsData.Approach in a cell
            SelectionData(k).selectfc_unload(i,1)=mapsData(k).unload_data(i,1);
            end    
        end
    kk=ii;
end
