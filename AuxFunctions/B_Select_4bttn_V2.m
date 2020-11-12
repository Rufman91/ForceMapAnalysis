% clc
% clear all
% close all

%% Define RGB colors

% RGB1=[0 26 255]./255;  % Blue 
% RGB2=[255 119 0]./255; % Orange

%% Pre-define needed variables
    l=0;
    i=1;
 
%% preallocate a structure
    if i == 1
    SelectionData = struct('MapNo',cell(i,1),...
                      'MapName',cell(i,1),...               
                      'subplotselect',{cell(i(i,1),1)},...
                      'numbselectfc',{cell(i(i,1),1)},...
                      'numbnonselectfc',{cell(i(i,1),1)},...
                      'selectfc_load',{cell(num_of_force_curves(i,1),1)},...
                      'selectfc_unload',{cell(num_of_force_curves(i,1),1)},...
                      'nonselectfc_load',{cell(num_of_force_curves(i,1),1)},...  
                      'nonselectfc_unload',{cell(num_of_force_curves(i,1),1)}...  
                      );
     % it is enough to preallocate once per force map, that is why we do not need to do it for the following 99 force curves
    end
    
for k=1:num_of_files
    %% Debugging: Give current Force Map Position
    sprintf('Force Map No. %d of %d',k,num_of_files)
    
    %% Force map number in the current analysis
        SelectionData(k).MapNo=sprintf('Force Map No. %d',k);
        
    %% Name the force maps in the 'SelectionData' structure
        SelectionData(k).MapName=originname{k,1};
        
    %% Figure and Subplots
    for j=1:4 % j corresponds to the number of subplots -> how many subplots are plotted
   
    %%% Define the names for the figure title    
        partname=sprintf('-Submap-%d',j);    
        fullname=sprintf('%s%s',originname{k,1}{1,1},partname);    
       
    %%% Define the margins of the subtightplot 
    %%% "subtightplot.m" is needed
    subplot = @(m,n,p) subtightplot (m, n, p,...
        [0.01 0.005],...% gap [vertical horizontal]
        [0.01 0.01],...% margin height [lower upper]
        [0.005 0.005]); % margin width [left right]
 
    %%% Plot the figure and subplots   
        figure(j)
        fig=gcf;     
        fig.Name=fullname; % Title for the whole figure
        fig.NumberTitle='off'; % switches off the figure number before the titlefor i=SelectionData(k).numbselectfc{j,1} % Deselect the choosen force curves
        fig.Units='normalized'; % changes to normalized unit settings, necessary to receive the full screen size in the next line
        fig.Color='white'; % changes the background color of the figure
     %  fig.WindowState='fullscreen';% to change the size to full screen % DID
     %  NOT WORK!?
        fig.OuterPosition=[0 0 1 1];% changes the size of the figure to full screen
        fig.PaperOrientation='landscape';

            for i = (1:25)+25*l % i corresponds to the number of force curves -> how many force curves are plotted              
                %%%% Define some variables
                dist=100e-9; % define a distance of 100nm
                x100=mapsData(k).load_data{i,1}(end,1)+dist;
                %%%% Plot the subplot
                subplot(5,5,i-25*l)
                hold on
                grid on
                plot(mapsData(k).load_data{i,1}(:,1),mapsData(k).load_data{i,1}(:,2),...
                    'Color',RGB1)
                plot(mapsData(k).unload_data{i,1}(:,1),mapsData(k).unload_data{i,1}(:,2),...
                    'Color',RGB2)
                line([x100 x100], ylim,... % Draws a line into the plot 
                    'Color',RGB8);
                % x100 ... at the defined x-position 
                % ylim ... Within the min and max y-value plotted       
                t=title(num2str(indent(i,1))); % Title for each Subplot
                t.Units='normalized'; % Set units to 'normalized'  
                t.Position=[0.5,0.9]; % Position the subplot title within the subplot
                clear dist x100
            end
            
        %%%% Add axes labels and legend only to the last Force curve of the submap
        legend('Load','Unload','Location','best')
        xlabel('Height ($\mu$m)','FontSize',11,'Interpreter','latex');
        ylabel('Deflection ($\mu$m)','FontSize',11,'Interpreter','latex');

    %% Dialog boxes

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
        
    %% Save all figures 

    %%% Create a folder "Figures" in datadir for saving the results
     mkdir(datadir,'B_Figures'); 
     Figuresdir=fullfile(datadir,'B_Figures');
          
     %%% Save the current figure
     %%% "export_fig folder" is needed          
     figfilename=fullfile(Figuresdir,fullname);  
     export_fig(figure(j),'-pdf', figfilename)
     
     %%% Close the current figure
     close(figure(j));
        
    %%% Adopt the variable 'l'    
    l=l+1;
    end
l=0;   % Reset to starting position for the next force map
end    

%% Save the created variables as .mat-file

    readdataname=strcat('B_Select_',date);
    readdatafilename=fullfile(Analysisdir,readdataname);
    save(readdatafilename);