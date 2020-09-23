classdef Experiment < matlab.mixin.Copyable
    
    properties
        ExperimentName
        ExperimentFolder
        ForceMapNames
        ForceMapFolders
        SurfacePotentialMapFolders
        SurfacePotentialMapNames
        NumFiles
        NumSpec
        NumMeas
        FM
        SPM
        FMFlag
        SPMFlag
        GroupFM
        GroupSPM
        DropoutNet
        CP_CNN
        CantileverTip
        
    end
    
    methods
        % constructor method and methods related with data handling
        
        function obj = Experiment()
            
            % Force Maps + KPFM or only one of them?
            answer = questdlg('What kind of measurements were done?', ...
                'Experiment Type',...
                'Surface Potential Maps','Indentation Force Maps','Both','Indentation Force Maps');
            % Handle response
            switch answer
                case 'Surface Potential Maps'
                    WhichFiles = 1;
                case 'Indentation Force Maps'
                    WhichFiles = 2;
                case 'Both'
                    WhichFiles = 0;
            end
            
            % How many Specimens were tested? Multiple measurements per
            % specimen?
            prompt = {'Enter Number of tested specimen','Enter number of measurements per specimen','How would you like to name the Experiment?'};
            dlgtitle = 'Experiment Layout';
            dims = [1 35];
            definput = {'20','2','MGO-Glycation'};
            answer = inputdlg(prompt,dlgtitle,dims,definput);
            
            
            
            obj.NumSpec = str2double(answer{1});
            obj.NumMeas = str2double(answer{2});
            obj.NumFiles = obj.NumMeas*obj.NumSpec;
            obj.ExperimentName = answer{3};
            obj.ExperimentFolder = uigetdir('Choose a Folder where the Experiment is to be saved');
            N = obj.NumFiles;
            obj.NumFiles = N;
            MapFullFile = {};
            k = 1;
            while length(MapFullFile) < N
                Title = sprintf('Choose one or more .jpk-force-map files. %i/%i',length(MapFullFile),N);
                [TempFile,TempPath] = uigetfile('*.jpk-force-map',Title,'MultiSelect','on');
                if  ~iscell(TempFile)
                    MapFullFile{k} = fullfile(TempPath,TempFile);
                    k = k + 1;
                else
                    for i=1:length(TempFile)
                        MapFullFile{k} = fullfile(TempPath,TempFile{i});
                        k = k + 1;
                    end
                end
                clear TempFile
            end
            obj.FM = cell(N,1);
            obj.SPM = cell(N,1);
            for i=1:N
                if WhichFiles == 2 || WhichFiles == 0
                    obj.FM{i} = ForceMap(MapFullFile{i},obj.ExperimentFolder);
                    obj.ForceMapFolders{i} = obj.FM{i}.Folder;
                    obj.ForceMapNames{i} = obj.FM{i}.Name;
                elseif WhichFiles == 1 || WhichFiles == 0
                    obj.SPM{i} = SurfacePotentialMap();
                    obj.SurfacePotentialMapFolders{i} = obj.SPM{i}.Folder;
                    obj.SurfacePotentialMapNames{i} = obj.SPM{i}.Name;
                end
            end
            obj.FMFlag.Analysis = zeros(N,1);
            obj.FMFlag.Grouping = 0;
            obj.SPMFlag.Analysis = zeros(N,1);
            obj.SPMFlag.Grouping = 0;
            
            if WhichFiles == 2 || WhichFiles == 0
                obj.grouping_force_map();
            elseif WhichFiles == 1 || WhichFiles == 0
                obj.grouping_surface_potential_map();
            end
            
            Temp = load('DropoutNet.mat');
            obj.DropoutNet = Temp.DropoutNet;
            Temp2 = load('CP_CNN.mat');
            obj.CP_CNN = Temp2.CP_CNN;
            
            obj.save_experiment();
        end
        
        function load_data(obj)
            for i=1:obj.NumFiles
                obj.FM{i} = ForceMap(obj.ForceMapFolders{i},obj.ForceMapNames{i});
                obj.SPM{i} = SurfacePotentialMap(obj.SurfacePotentialMapFolders{i},obj.SurfacePotentialMapNames{i});
            end
        end
        
        function save_data(obj)
            for i=1:obj.NumFiles
                disp('')
                obj.FM{i}.save();
                obj.ForceMapFolders{i} = obj.FM{i}.Folder;
                obj.ForceMapNames{i} = obj.FM{i}.Name;
                obj.SPM{i}.save();
                obj.SurfacePotentialMapFolders{i} = obj.SPM{i}.Folder;
                obj.SurfacePotentialMapNames{i} = obj.SPM{i}.Name;
            end
        end
        
        function save_experiment(obj)
            current = what();
            cd(obj.ExperimentFolder)
            savename = sprintf('%s.mat',obj.ExperimentName);
            save(savename,'obj')
            cd(current.path)
            savemsg = sprintf('Changes to Experiment %s saved to %s',obj.ExperimentName,obj.ExperimentFolder);
            disp(savemsg);
        end
    end
    
    methods
       % methods for sequential data analysis mostly looping over child-classes methods
       
       function force_map_analysis_fibril(obj,CPOption,EModOption)
           % CPOption = 'Fast' ...(Default) contact point estimation through single
           % pass through CNN
           % CPOption = 'Dropout' ... contact point estimation through
           % averaging over multiple passes through monte carlo dropout
           % net. NPasses (Default=100) times slower than 'Fast'
           % CPOption = 'Old' ... old method for contact point estimation
           %
           % EModOption = 'Hertz' ... E-Modulus calculation through Hertz-Sneddon
           % method
           % EModOption = 'Oliver' ... E-Modulus calculation through
           % Oliver-Pharr-like method (O. Andriotis 2014)
           
           h = waitbar(0,'setting up','Units','normalized','Position',[0.4 0.3 0.2 0.1]);
           NLoop = length(obj.ForceMapNames);
           if sum(obj.FMFlag.Analysis) >= 1
               KeepFlagged = questdlg(sprintf('Some maps have been processed already.\nDo you want to skip them and keep old results?'),...
                   'Processing Options',...
                   'Yes',...
                   'No',...
                   'No');
           else
               KeepFlagged = 'No';
           end
           % Preprocessing everything, that needs user input
           answer = questdlg('Do you want to skip manual exclusion of problematic areas?',...
               'Manual Exclusion',...
               'Yes',...
               'No','No');
           for i=1:NLoop
               if isequal(KeepFlagged,'Yes') && obj.FMFlag.Analysis(i) == 1
                   continue
               end
               obj.FM{i}.create_and_level_height_map();
               obj.FM{i}.create_fibril_mask();
               if isequal(answer,'Yes')
                   continue
               end
               obj.FM{i}.manual_exclusion();
           end
           % Main loop for contact point estimation, Fibril Diameter and
           % Fibril E-Modulus calculation
           for i=1:NLoop
               if isequal(KeepFlagged,'Yes') && obj.FMFlag.Analysis(i) == 1
                   continue
               end
               waitbar(i/NLoop,h,sprintf('Processing Fibril %i/%i\nFitting Base Line',i,NLoop));
               obj.FM{i}.base_and_tilt('linear');
               obj.FM{i}.calculate_fib_diam();
               waitbar(i/NLoop,h,sprintf('Processing Fibril %i/%i\nFinding Contact Point',i,NLoop));
               if isequal(lower(CPOption),'dropout')
                   obj.FM{i}.estimate_cp_cnn(obj.DropoutNet,'Dropout',100);
               elseif isequal(lower(CPOption),'old')
                   obj.FM{i}.estimate_cp_old_oliverpharr;
               else
                   obj.FM{i}.estimate_cp_cnn(obj.CP_CNN,'Fast')
               end
               waitbar(i/NLoop,h,sprintf('Processing Fibril %i/%i\nCalculating E-Modulus',i,NLoop));
               if isequal(lower(EModOption),'hertz')
                   obj.FM{i}.calculate_e_mod_hertz('cnn','parabolic',1);
               else
                   obj.FM{i}.calculate_e_mod_oliverpharr(obj.CantileverTip.ProjArea,0.75);
               end
               waitbar(i/NLoop,h,sprintf('Processing Fibril %i/%i\nWrapping Up And Saving',i,NLoop));
               obj.FM{i}.show_height_map();
               obj.FM{i}.save();
               obj.FMFlag.Analysis(i) = 1;
           end
           close(h);
       end
       
       function surface_potential_analysis_fibril(obj)
           
       end
       
       function force_map_analysis_general(obj)
           % WORK IN PROGRESS
       end
       
       function surface_potential_analysis_general(obj)
           % WORK IN PROGRESS
       end
       
       function statistical_analysis_force_maps(obj)
           % Basic statistical analysis of the force map experiment. First get the grouping
           % of data from the user and then perform several tests, plots
           % etc. 
           N = obj.NumFiles;
           obj.grouping_force_map();
           
           % Ask user which groups are to be compared
           prompt = 'Specialize which groups are to be tested against which';
           definput = {'1 2 ; 3 4'};
           answer = inputdlg(prompt,'Pairings',[1 50],definput);
           
           % Define Test Matrix
           TestMat = [str2num(answer{1})];
           
           % Write Data into local variables and
           % replace extreme outlier values (median +- 2.5*IQR) with NaN
           % also replace data points that lie on the exclusion mask with
           % NaN
           for i=1:N
            DataOP(i,:) = obj.FM{i}.EModOliverPharr(obj.FM{i}.RectApexIndex);
            DataHS(i,:) = obj.FM{i}.EModHertz(obj.FM{i}.RectApexIndex);
            for j=1:length(obj.FM{i}.RectApexIndex)
                if DataOP(i,j) > (nanmedian(DataOP(i,:))+2.5*iqr(DataOP(i,:))) || ...
                        DataOP(i,j) < (nanmedian(DataOP(i,:))-2.5*iqr(DataOP(i,:))) || ...
                        obj.FM{i}.ExclMask(obj.FM{i}.List2Map(obj.FM{i}.RectApexIndex(j),1),obj.FM{i}.List2Map(obj.FM{i}.RectApexIndex(j),2)) == 0
                    DataOP(i,j) = NaN;
                elseif DataOP(i,j) < 0
                    DataOP(i,j) = NaN;
                end
                if DataHS(i,j) > (nanmedian(DataHS(i,:))+2.5*iqr(DataHS(i,:))) || ...
                        DataHS(i,j) < (nanmedian(DataHS(i,:))-2.5*iqr(DataHS(i,:))) || ...
                        obj.FM{i}.ExclMask(obj.FM{i}.List2Map(obj.FM{i}.RectApexIndex(j),1),obj.FM{i}.List2Map(obj.FM{i}.RectApexIndex(j),2)) == 0
                    DataHS(i,j) = NaN;
                elseif DataHS(i,j) < 0
                    DataHS(i,j) = NaN;
                end
            end
           end
           
           DataMeansOP = nanmean(DataOP,2);
           DataMeansHS = nanmean(DataHS,2);
           
           figure('Name','OliverPharr vs HertzSneddon','Color','w');
           plot(1:N,DataMeansOP,'bO',1:N,DataMeansHS,'rO')
           legend('E-Mod Oliver-Pharr','E-Mod Hertz-Sneddon')
           xlim([0,N+1])
           xlabel('Force Maps')
           ylabel('E-Mod [Pa]')
           
           % loop over all rows of the Test Matrix, doing paired ttests 
           for i=1:size(TestMat,1)
               % Statistics for Oliver-Pharr Method
               [hOP(i),pOP(i)] = ...
                   ttest(DataMeansOP(obj.GroupFM(TestMat(i,2)).Indices),...
                   DataMeansOP(obj.GroupFM(TestMat(i,1)).Indices),'Tail','right');
               
               figure('Name','Paired Right Tailed T-Test','Units','normalized','Position',[0.2 0.2 0.5 0.5],'Color','w')
               boxplot([DataMeansOP(obj.GroupFM(TestMat(i,1)).Indices) DataMeansOP(obj.GroupFM(TestMat(i,2)).Indices)],'Notch','on')
               title('Paired Right Tailed T-Test for Oliver-Pharr Method')
               xticklabels({obj.GroupFM(TestMat(i,1)).Name,obj.GroupFM(TestMat(i,2)).Name})
               xlabel('Test Group')
               ylabel('E-Mod Oliver-Pharr[Pa]')
               DeltaMean = mean(DataMeansOP(obj.GroupFM(TestMat(i,2)).Indices)) - mean(DataMeansOP(obj.GroupFM(TestMat(i,1)).Indices));
               Sigma = std(DataMeansOP(obj.GroupFM(TestMat(i,2)).Indices) - DataMeansOP(obj.GroupFM(TestMat(i,1)).Indices));
               Beta = sampsizepwr('t',[0 Sigma],DeltaMean,[],length(obj.GroupFM(TestMat(i,2)).Indices));
               Stats = {sprintf('\\DeltaMean = %.2fMPa',DeltaMean*1e-6),...
                   sprintf('P-Value = %.4f%',pOP(i)),...
                   sprintf('Power \\beta = %.2f%%',Beta*100),...
                   sprintf('Number of Specimen = %i',length(obj.GroupFM(TestMat(i,2)).Indices))};
               text(0.5,0.8,Stats,...
                   'Units','normalized',...
                   'FontSize',12,...
                   'HorizontalAlignment','center')
               
               %Statistics for Hertz-Sneddon Method
               [hHS(i),pHS(i)] = ...
                   ttest(DataMeansHS(obj.GroupFM(TestMat(i,2)).Indices),...
                   DataMeansHS(obj.GroupFM(TestMat(i,1)).Indices),'Tail','right');
               
               figure('Name','Paired Right Tailed T-Test','Units','normalized','Position',[0.2 0.2 0.5 0.5],'Color','w')
               boxplot([DataMeansHS(obj.GroupFM(TestMat(i,1)).Indices) DataMeansHS(obj.GroupFM(TestMat(i,2)).Indices)],'Notch','on')
               title('Paired Right Tailed T-Test for Hertz-Sneddon Method')
               xticklabels({obj.GroupFM(TestMat(i,1)).Name,obj.GroupFM(TestMat(i,2)).Name})
               xlabel('Test Group')
               ylabel('E-Mod Hertz-Sneddon[Pa]')
               DeltaMean = mean(DataMeansHS(obj.GroupFM(TestMat(i,2)).Indices)) - mean(DataMeansHS(obj.GroupFM(TestMat(i,1)).Indices));
               Sigma = std(DataMeansHS(obj.GroupFM(TestMat(i,2)).Indices) - DataMeansHS(obj.GroupFM(TestMat(i,1)).Indices));
               Beta = sampsizepwr('t',[0 Sigma],DeltaMean,[],length(obj.GroupFM(TestMat(i,2)).Indices));
               Stats = {sprintf('\\DeltaMean = %.2f MPa',DeltaMean*1e-6),...
                   sprintf('P-Value = %.4f%',pHS(i)),...
                   sprintf('Power \\beta = %.2f%%',Beta*100),...
                   sprintf('Number of Specimen = %i',length(obj.GroupFM(TestMat(i,2)).Indices))};
               text(0.5,0.8,Stats,...
                   'Units','normalized',...
                   'FontSize',12,...
                   'HorizontalAlignment','center')
           end
           
           % Now test, if the difference in one pair of groups is
           % statistically different from the other in a two sample t test
           DiffControlOP = DataMeansOP(11:20) - DataMeansOP(1:10);
           DiffMGOOP = DataMeansOP(31:40) - DataMeansOP(21:30);
           DiffControlHS = DataMeansHS(11:20) - DataMeansHS(1:10);
           DiffMGOHS = DataMeansHS(31:40) - DataMeansHS(21:30);
           
           % For Oliver-Pharr
           [hOP,pOP,ciOP,statsOP] = ttest2(DiffMGOOP,DiffControlOP);
           figure('Name','Two Sample T-Test','Units','normalized','Position',[0.2 0.2 0.5 0.5],'Color','w')
               yyaxis left
               boxplot([DiffControlOP DiffMGOOP],'Notch','on')
               ax = gca;
               YLim = ax.YLim;
               ylabel('Difference Before-After E-Mod [Pa]')
               DeltaMean = mean(DiffMGOOP) - mean(DiffControlOP);
               PooledSTD = statsOP.sd;
               yyaxis right
               errorbar(1.5,DeltaMean,ciOP(2)-DeltaMean,'O');
               ylim(YLim)
               xticks([1 1.5 2])
               title('Two Sample T-Test for E-Mod Oliver-Pharr Method')
               ax = gca;
               ax.TickLabelInterpreter = 'tex';
               xticklabels({sprintf('%s - %s',obj.GroupFM(2).Name,obj.GroupFM(1).Name),...
                   '\DeltaMean with CI',...
                   sprintf('%s - %s',obj.GroupFM(4).Name,obj.GroupFM(3).Name)})
               ylabel('Difference of Differences [Pa]')
               Beta = sampsizepwr('t2',[mean(DiffControlOP) PooledSTD],mean(DiffMGOOP),[],length(DiffControlOP),'Ratio',length(DiffMGOOP)/length(DiffControlOP));
               Stats = {sprintf('\\DeltaMean = %.2f MPa',DeltaMean*1e-6),...
                   sprintf('P-Value = %.4f%',pOP),...
                   sprintf('Power \\beta = %.2f%%',Beta*100),...
                   sprintf('Degrees of freedom df = %i',statsOP.df)};
               text(0.5,0.8,Stats,...
                   'Units','normalized',...
                   'FontSize',12,...
                   'HorizontalAlignment','center')
           
           % For Hertz-Sneddon
           [hHS,pHS,ciHS,statsHS] = ttest2(DiffMGOHS,DiffControlHS);
           figure('Name','Two Sample T-Test','Units','normalized','Position',[0.2 0.2 0.5 0.5],'Color','w')
               yyaxis left
               boxplot([DiffControlHS DiffMGOHS],'Notch','on')
               ax = gca;
               YLim = ax.YLim;
               ylabel('Difference Before-After E-Mod [Pa]')
               DeltaMean = mean(DiffMGOHS) - mean(DiffControlHS);
               PooledSTD = statsHS.sd;
               yyaxis right
               errorbar(1.5,DeltaMean,ciHS(2)-DeltaMean,'O');
               ylim(YLim)
               xticks([1 1.5 2])
               title('Two Sample T-Test for E-Mod Hertz-Sneddon Method')
               ax = gca;
               ax.TickLabelInterpreter = 'tex';
               xticklabels({sprintf('%s - %s',obj.GroupFM(2).Name,obj.GroupFM(1).Name),...
                   '\DeltaMean with CI',...
                   sprintf('%s - %s',obj.GroupFM(4).Name,obj.GroupFM(3).Name)})
               ylabel('Difference of Differences [Pa]')
               Beta = sampsizepwr('t2',[mean(DiffControlHS) PooledSTD],mean(DiffMGOHS),[],length(DiffControlHS),'Ratio',length(DiffMGOHS)/length(DiffControlHS));
               Stats = {sprintf('\\DeltaMean = %.2f MPa',DeltaMean*1e-6),...
                   sprintf('P-Value = %.4f%',pHS),...
                   sprintf('Power \\beta = %.2f%%',Beta*100),...
                   sprintf('Degrees of freedom df = %i',statsHS.df)};
               text(0.5,0.8,Stats,...
                   'Units','normalized',...
                   'FontSize',12,...
                   'HorizontalAlignment','center')
           
       end
       
       function statistical_analysis_surface_potential(obj)
           % Statistical analysis of the surface potential map experiment.
           % First get the grouping of data from the user and then perform
           % several tests, plots etc.
           
           obj.grouping_surface_potential_map();
           
           %%%%%%% DISCLAIMER: just works for specific cases at the moment %%%%%%%
           
           N = obj.NumFiles;
           
           % Ask user which groups are to be compared
           prompt = 'Specialize which groups are to be tested against which';
           definput = {'1 2 ; 3 4'};
           answer = inputdlg(prompt,'Pairings',[1 50],definput);
           
           % Define Test Matrix
           TestMat = [str2num(answer{1})];
           % Get Data
           for i=1:N
               FibPot(i) = obj.SPM{i}.FibPot;
           end
           % Do Paired T-Test of Surface Potential between paired groups
           for i=1:size(TestMat,1)
               [h(i),p(i)] = ...
                   ttest(FibPot(obj.GroupFM(TestMat(i,2)).Indices),...
                   FibPot(obj.GroupFM(TestMat(i,1)).Indices));
               
               figure('Name','Paired T-Test for Surface Potential Changes','Units','normalized','Position',[0.2 0.2 0.5 0.5],'Color','w')
               boxplot([FibPot(obj.GroupFM(TestMat(i,1)).Indices)' FibPot(obj.GroupFM(TestMat(i,2)).Indices)'])
               title('Paired T-Test for Surface Potential Changes')
               xticklabels({obj.GroupFM(TestMat(i,1)).Name,obj.GroupFM(TestMat(i,2)).Name})
               xlabel('Test Group')
               ylabel('Surface Potential [V]')
               DeltaMean = mean(FibPot(obj.GroupFM(TestMat(i,2)).Indices)) - mean(FibPot(obj.GroupFM(TestMat(i,1)).Indices));
               Sigma = std(FibPot(obj.GroupFM(TestMat(i,2)).Indices) - FibPot(obj.GroupFM(TestMat(i,1)).Indices));
               Beta = sampsizepwr('t',[0 Sigma],DeltaMean,[],length(obj.GroupFM(TestMat(i,2)).Indices));
               Stats = {sprintf('\\DeltaMean = %.2fmV',DeltaMean*1e3),...
                   sprintf('P-Value = %.2e%',p(i)),...
                   sprintf('Power \\beta = %.2f%%',Beta*100),...
                   sprintf('Number of Specimen = %i',length(obj.GroupFM(TestMat(i,2)).Indices))};
               text(0.5,0.8,Stats,...
                   'Units','normalized',...
                   'FontSize',12,...
                   'HorizontalAlignment','center')
           end
           
           % Now test, if the difference in one pair of groups is
           % statistically different from the other in a two sample t test
           DiffControl = (FibPot(11:20) - FibPot(1:10))';
           DiffMGO = (FibPot(31:40) - FibPot(21:30))';
           
           
           [h,p,ci,stats] = ttest2(DiffMGO,DiffControl);
           figure('Name','Two Sample T-Test','Units','normalized','Position',[0.2 0.2 0.5 0.5],'Color','w')
               yyaxis left
               boxplot([DiffControl DiffMGO],'Notch','on')
               ax = gca;
               YLim = ax.YLim;
               ylabel('Difference Before-After Potential [V]')
               DeltaMean = mean(DiffMGO) - mean(DiffControl);
               PooledSTD = stats.sd;
               yyaxis right
               errorbar(1.5,DeltaMean,ci(2)-DeltaMean,'O');
               ylim(YLim)
               xticks([1 1.5 2])
               title('Two Sample T-Test for Surface Potential Changes')
               ax = gca;
               ax.TickLabelInterpreter = 'tex';
               xticklabels({sprintf('%s - %s',obj.GroupSPM(2).Name,obj.GroupSPM(1).Name),...
                   '\DeltaMean with CI',...
                   sprintf('%s - %s',obj.GroupSPM(4).Name,obj.GroupSPM(3).Name)})
               ylabel('Difference of Differences [V]')
               Beta = sampsizepwr('t2',[mean(DiffControl) PooledSTD],mean(DiffMGO),[],length(DiffControl),'Ratio',length(DiffMGO)/length(DiffControl));
               Stats = {sprintf('\\DeltaMean = %.2fmV',DeltaMean*1e3),...
                   sprintf('P-Value = %.2e%',p),...
                   sprintf('Power \\beta = %.2f%%',Beta*100),...
                   sprintf('Degrees of freedom df = %i',stats.df)};
               text(0.5,0.8,Stats,...
                   'Units','normalized',...
                   'FontSize',12,...
                   'HorizontalAlignment','center')
       end
       
       function statistical_analysis_swelling(obj)
           %%%%%%% DISCLAIMER: just works for specific cases at the moment %%%%%%%
           
           N = obj.NumFiles;
           
           % Ask user which groups are to be compared
           prompt = 'Specialize which groups are to be tested against which';
           definput = {'1 2 ; 3 4'};
           answer = inputdlg(prompt,'Pairings',[1 50],definput);
           
           % Define Test Matrix
           TestMat = [str2num(answer{1})];
           % Get Data
           for i=1:N
               Dry(i) = obj.SPM{i}.FibDiam;
               Wet(i) = obj.FM{i}.FibDiam;
               RelChange(i) = Wet(i)/Dry(i);
           end
           % Do Paired T-Test of relative swelling between paired groups
           for i=1:size(TestMat,1)
               [h(i),p(i)] = ...
                   ttest(RelChange(obj.GroupFM(TestMat(i,2)).Indices),...
                   RelChange(obj.GroupFM(TestMat(i,1)).Indices));
               
               figure('Name','Paired T-Test for relative Swelling Changes','Units','normalized','Position',[0.2 0.2 0.5 0.5],'Color','w')
               boxplot([RelChange(obj.GroupFM(TestMat(i,1)).Indices)' RelChange(obj.GroupFM(TestMat(i,2)).Indices)'])
               title('Paired T-Test for relative Swelling Changes')
               xticklabels({obj.GroupFM(TestMat(i,1)).Name,obj.GroupFM(TestMat(i,2)).Name})
               xlabel('Test Group')
               ylabel('Relative Swelling')
               DeltaMean = mean(RelChange(obj.GroupFM(TestMat(i,2)).Indices)) - mean(RelChange(obj.GroupFM(TestMat(i,1)).Indices));
               Sigma = std(RelChange(obj.GroupFM(TestMat(i,2)).Indices) - RelChange(obj.GroupFM(TestMat(i,1)).Indices));
               Beta = sampsizepwr('t',[0 Sigma],DeltaMean,[],length(obj.GroupFM(TestMat(i,2)).Indices));
               Stats = {sprintf('\\DeltaMean = %.2f',DeltaMean),...
                   sprintf('P-Value = %.3f%',p(i)),...
                   sprintf('Power \\beta = %.2f%%',Beta*100),...
                   sprintf('Number of Specimen = %i',length(obj.GroupFM(TestMat(i,2)).Indices))};
               text(0.5,0.8,Stats,...
                   'Units','normalized',...
                   'FontSize',12,...
                   'HorizontalAlignment','center')
           end
           
           % Now test, if the difference in one pair of groups is
           % statistically different from the other in a two sample t test
           DiffControl = (RelChange(11:20) - RelChange(1:10))';
           DiffMGO = (RelChange(31:40) - RelChange(21:30))';
           
           
           [h,p,ci,stats] = ttest2(DiffMGO,DiffControl);
           figure('Name','Two Sample T-Test','Units','normalized','Position',[0.2 0.2 0.5 0.5],'Color','w')
               yyaxis left
               boxplot([DiffControl DiffMGO],'Notch','on')
               ax = gca;
               YLim = ax.YLim;
               ylabel('Difference Before-After relative Swelling')
               DeltaMean = mean(DiffMGO) - mean(DiffControl);
               PooledSTD = stats.sd;
               yyaxis right
               errorbar(1.5,DeltaMean,ci(2)-DeltaMean,'O');
               ylim(YLim)
               xticks([1 1.5 2])
               title('Two Sample T-Test for relative Swelling')
               ax = gca;
               ax.TickLabelInterpreter = 'tex';
               xticklabels({sprintf('%s - %s',obj.GroupSPM(2).Name,obj.GroupSPM(1).Name),...
                   '\DeltaMean with CI',...
                   sprintf('%s - %s',obj.GroupSPM(4).Name,obj.GroupSPM(3).Name)})
               ylabel('Difference of Differences')
               Beta = sampsizepwr('t2',[mean(DiffControl) PooledSTD],mean(DiffMGO),[],length(DiffControl),'Ratio',length(DiffMGO)/length(DiffControl));
               Stats = {sprintf('\\DeltaMean = %.2f',DeltaMean),...
                   sprintf('P-Value = %.2e%',p),...
                   sprintf('Power \\beta = %.2f%%',Beta*100),...
                   sprintf('Degrees of freedom df = %i',stats.df)};
               text(0.5,0.8,Stats,...
                   'Units','normalized',...
                   'FontSize',12,...
                   'HorizontalAlignment','center')
       end
       
       function statistical_analysis_d_banding(obj)
                      %%%%%%% DISCLAIMER: just works for specific cases at the moment %%%%%%%
           
           N = obj.NumFiles;
           
           % Ask user which groups are to be compared
           prompt = 'Specialize which groups are to be tested against which';
           definput = {'1 2 ; 3 4'};
           answer = inputdlg(prompt,'Pairings',[1 50],definput);
           
           % Define Test Matrix
           TestMat = [str2num(answer{1})];
           % Get Data
           for i=1:N
               DBanding(i) = obj.SPM{i}.DBanding;
           end
           % Do Paired T-Test of relative swelling between paired groups
           for i=1:size(TestMat,1)
               [h(i),p(i)] = ...
                   ttest(DBanding(obj.GroupFM(TestMat(i,2)).Indices),...
                   DBanding(obj.GroupFM(TestMat(i,1)).Indices));
               
               figure('Name','Paired T-Test for D-Banding Changes','Units','normalized','Position',[0.2 0.2 0.5 0.5],'Color','w')
               boxplot([DBanding(obj.GroupFM(TestMat(i,1)).Indices)' DBanding(obj.GroupFM(TestMat(i,2)).Indices)'])
               title('Paired T-Test for D-Banding Changes')
               xticklabels({obj.GroupFM(TestMat(i,1)).Name,obj.GroupFM(TestMat(i,2)).Name})
               xlabel('Test Group')
               ylabel('D-Banding [m]')
               DeltaMean = mean(DBanding(obj.GroupFM(TestMat(i,2)).Indices)) - mean(DBanding(obj.GroupFM(TestMat(i,1)).Indices));
               Sigma = std(DBanding(obj.GroupFM(TestMat(i,2)).Indices) - DBanding(obj.GroupFM(TestMat(i,1)).Indices));
               Beta = sampsizepwr('t',[0 Sigma],DeltaMean,[],length(obj.GroupFM(TestMat(i,2)).Indices));
               Stats = {sprintf('\\DeltaMean = %.2fnm',DeltaMean*1e9),...
                   sprintf('P-Value = %.3f%',p(i)),...
                   sprintf('Power \\beta = %.2f%%',Beta*100),...
                   sprintf('Number of Specimen = %i',length(obj.GroupFM(TestMat(i,2)).Indices))};
               text(0.5,0.8,Stats,...
                   'Units','normalized',...
                   'FontSize',12,...
                   'HorizontalAlignment','center')
           end
           
           % Now test, if the difference in one pair of groups is
           % statistically different from the other in a two sample t test
           DiffControl = (DBanding(11:20) - DBanding(1:10))';
           DiffMGO = (DBanding(31:40) - DBanding(21:30))';
           
           
           [h,p,ci,stats] = ttest2(DiffMGO,DiffControl);
           figure('Name','Two Sample T-Test','Units','normalized','Position',[0.2 0.2 0.5 0.5],'Color','w')
               yyaxis left
               boxplot([DiffControl DiffMGO],'Notch','on')
               ax = gca;
               YLim = ax.YLim;
               ylabel('Difference Before-After in D-Banding [m]')
               DeltaMean = mean(DiffMGO) - mean(DiffControl);
               PooledSTD = stats.sd;
               yyaxis right
               errorbar(1.5,DeltaMean,ci(2)-DeltaMean,'O');
               ylim(YLim)
               xticks([1 1.5 2])
               title('Two Sample T-Test for D-Banding Changes')
               ax = gca;
               ax.TickLabelInterpreter = 'tex';
               xticklabels({sprintf('%s - %s',obj.GroupSPM(2).Name,obj.GroupSPM(1).Name),...
                   '\DeltaMean with CI',...
                   sprintf('%s - %s',obj.GroupSPM(4).Name,obj.GroupSPM(3).Name)})
               ylabel('Difference of Differences [m]')
               Beta = sampsizepwr('t2',[mean(DiffControl) PooledSTD],mean(DiffMGO),[],length(DiffControl),'Ratio',length(DiffMGO)/length(DiffControl));
               Stats = {sprintf('\\DeltaMean = %.2fnm',DeltaMean*1e9),...
                   sprintf('P-Value = %.2e%',p),...
                   sprintf('Power \\beta = %.2f%%',Beta*100),...
                   sprintf('Degrees of freedom df = %i',stats.df)};
               text(0.5,0.8,Stats,...
                   'Units','normalized',...
                   'FontSize',12,...
                   'HorizontalAlignment','center')
       end
       
       function grouping_force_map(obj)
           % A series of input dialogues that determine the structure and
           % relations of the force map files in the experiment
           
           if obj.FMFlag.Grouping == 1
               disp('Force Maps already have a Grouping assigned')
               return
           elseif obj.SPMFlag.Grouping == 1
               answer = questdlg('Is the Grouping and File order the same as for Surface Potential Maps?');
               if isequal(lower(answer),'yes')
                   obj.GroupFM = obj.GroupSPM;
                   obj.FMFlag.Grouping = 1;
                   disp('Same Grouping as Surface Potential Maps assigned to Force Maps')
                   return
               else
               end
           end
           
           % How many statistical groups are there?
           answer1 = inputdlg('How many different groupings of data are there? (E.g. Control-Before and Control-After count as two separate groups)',...
               'Statistical Groupings',[1 50]);
           NGroups = str2num(answer1{1});
           
           % create the appropriate inputdlg for assigning the groups show
           % a table with numbered map-names in background
           Names = obj.ForceMapNames;
           Fig = figure('Units', 'Normalized', 'Position',[0, 0, 0.4, 1],'Color','w');
           T = table(Names');
           uitable('Data',T{:,:},'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
           
           for i=1:2:2*NGroups
               prompts{i} = sprintf('Whats the Name of Group %i?',(i+1)/2);
               prompts{i+1} = sprintf('Which Indices belong to Group %i?',(i+1)/2);
               definput{i} = sprintf('Group %i Name',(i+1)/2);
               definput{i+1} = 'e.g. 1 2 3 4 8 9 10 or 1:4 8:10';
           end
           
           dims = [1 50]; 
           opts.WindowStyle = 'normal';
           dlgtitle = 'Choose and name Groups';
           answer2 = inputdlg(prompts,dlgtitle,dims,definput,opts);
           
           for i=1:NGroups
               obj.GroupFM(i).Name = answer2{2*i-1};
               obj.GroupFM(i).Indices = str2num(answer2{2*i});
           end
           obj.FMFlag.Grouping = 1;
           close(Fig);
           
       end
       
       function grouping_surface_potential_map(obj)
           % A series of input dialogues that determine the structure and
           % relations of the force map files in the experiment
           
           if obj.SPMFlag.Grouping == 1
               disp('Surface Potential Maps already have a Grouping assigned')
               return
           elseif obj.FMFlag.Grouping == 1
               answer = questdlg('Is the Grouping and File order the same as for Force Maps?');
               if isequal(lower(answer),'yes')
                   obj.GroupSPM = obj.GroupFM;
                   disp('Same Grouping as Force Maps assigned to Surface Potential Maps')
                   obj.SPMFlag.Grouping = 1;
                   return
               else
               end
           end
           
           % How many statistical groups are there?
           answer1 = inputdlg('How many different groupings of data are there? (E.g. Control-Before and Control-After count as two separate groups)',...
               'Statistical Groupings',[1 50]);
           NGroups = str2num(answer1{1});
           
           % create the appropriate inputdlg for assigning the groups show
           % a table with numbered map-names in background
           Names = obj.SurfacePotentialMapNames;
           Fig = figure('Units', 'Normalized', 'Position',[0, 0, 0.4, 1],'Color','w');
           T = table(Names');
           uitable('Data',T{:,:},'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
           
           for i=1:2:2*NGroups
               prompts{i} = sprintf('Whats the Name of Group %i?',(i+1)/2);
               prompts{i+1} = sprintf('Which Indices belong to Group %i?',(i+1)/2);
               definput{i} = sprintf('Group %i Name',(i+1)/2);
               definput{i+1} = 'e.g. 1 2 3 4 8 9 10 or 1:4 8:10';
           end
           
           dims = [1 50]; 
           opts.WindowStyle = 'normal';
           dlgtitle = 'Choose and name Groups';
           answer2 = inputdlg(prompts,dlgtitle,dims,definput,opts);
           
           for i=1:NGroups
               obj.GroupSPM(i).Name = answer2{2*i-1};
               obj.GroupSPM(i).Indices = str2num(answer2{2*i});
           end
           
           obj.SPMFlag.Grouping = 1;
           
           close(Fig);
           
       end
       
    end
    
    methods
       % auxiliary methods
       
       function RadiusNM = calculate_tip_radius(obj,TipDepthNM)
            if nargin < 2
                TipDepthNM = 20;
            end
            m = 1;
            Niter = 200;
            Cumulative = zeros(Niter,1);
            for n=1:Niter
                k = 1;
                SizeArray = size(obj.CantileverTip.data);
                for i=1:SizeArray(1)
                    for j=1:SizeArray(2)
                        if obj.CantileverTip.data(i,j) > -TipDepthNM*1e-9
                            X(k) = (i-1)*obj.CantileverTip.XaxisSizeUM*1e-6/SizeArray(1);
                            Y(k) = (j-1)*obj.CantileverTip.YaxisSizeUM*1e-6/SizeArray(2);
                            Z(k) = obj.CantileverTip.data(i,j);
                            k = k + 1;
                        end
                    end
                end
                
                %             SurfFit = fit([X', Y'],Z','poly33');
                %             syms f(x,y)
                %             % use for 'poly33'
                %             f(x,y) = SurfFit.p00 + SurfFit.p10*x...
                %                 + SurfFit.p20*x^2 + SurfFit.p11*x*y...
                %                 + SurfFit.p02*y^2 + SurfFit.p30*x^3 ...
                %                 + SurfFit.p21*x^2*y + SurfFit.p12*x*y^2 ...
                %                 + SurfFit.p03*y^3 ;
                %             % Use for 'poly22'
                % %             f(x,y) = SurfFit.p00 + SurfFit.p10*x...
                % %                 + SurfFit.p20*x^2 + SurfFit.p11*x*y...
                % %                 + SurfFit.p02*y^2;
                %         % The mean curvature of a surface (weighted integral over all
                %         % possible 1D curves on the surface) for a surface of the form
                %         % z=f(x,y) is given as half the negative divergence of the surface
                %         % unit normal at any given point, which expands to the formula
                %         % given here:
                %             MeanCurvature(x,y) = 1/2*abs(((1+diff(f,x)^2)*diff(f,y,y)...
                %                 -2*diff(f,x)*diff(f,y)*diff(f,x,y)...
                %                 +(1+diff(f,y)^2)*diff(f,x,x))...
                %                 /(1+diff(f,x)^2+diff(f,y)^2)^(3/2));
                %             RadiusNM = double(1/MeanCurvature(8.966*1e-7,1.124.*1e-6)*1e9)
                PC(:,1) = X*1e6;
                PC(:,2) = Y*1e6;
                PC(:,3) = Z*1e6;
                PC=pointCloud(PC);
                Sph = pcfitsphere(PC,1e-1,'Confidence',99.9,'MaxNumTrials',10000);
                plot(Sph)
                axis equal
                hold on
                scatter3(PC.Location(:,1),PC.Location(:,2),PC.Location(:,3))
                hold off
                pause(0.01)
                if n == 1
                    BestSph = Sph;
                elseif Sph.Radius < BestSph.Radius
                    BestSph = Sph;
                end
                RadiusNM = Sph.Radius*1e3;
                Cumulative(m) = RadiusNM;
                m = m + 1;
                clear PC
            end
            Sorted = sort(Cumulative);
            RadiusNM = mean(Sorted(1:round(Niter/2)));
            obj.CantileverTip.RadiusNM = RadiusNM;
            for i=1:obj.NumFiles
                obj.FM{i}.TipRadius = RadiusNM;
            end
       end
       
       function show_tip_data(obj)
           figure('Name','Cantilever Tip',...
               'Units','normalized','Position',[0.7 0.1 0.3 0.8],'Color','w');
           subplot(3,1,1)
           surf(obj.CantileverTip.data,'LineStyle','none','FaceLighting','gouraud','FaceColor','interp')
           light('Style','local')
           subplot(3,1,2)
           plot(1:length(obj.CantileverTip.ProjArea),obj.CantileverTip.ProjArea)
           xlabel('Indentation Depth [nm]')
           ylabel('Projected Area [m^2]')
           subplot(3,1,3)
           plot(1:30,obj.CantileverTip.ProjArea(1:30))
           xlabel('Typical Indentation Depth [nm]')
           ylabel('Projected Area [m^2]')
       end
        
    end
end