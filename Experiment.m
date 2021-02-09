classdef Experiment < matlab.mixin.Copyable
    
    properties
        ExperimentName      % Shows name of the experiment
        ExperimentFolder    % Shows Folder where Experiment is saved
        HostOS              % Shows the current operating system
        HostName            % Shows the current Host (User of running machine)
        ForceMapNames       % Shows names of the force maps
        ForceMapFolders
        SurfacePotentialMapFolders
        SurfacePotentialMapNames
        NumFiles
        FM
        EMod
        RefFM
        SPM
        SurfPot
        FMFlag
        SPMFlag
        GroupFM
        GroupSPM
        DropoutNet
        CP_CNN
        CantileverTip
        CantileverTipFlag
        idxSubstrate
        idxEnvCond
        
    end
    
    methods
        % constructor method and methods related with Experiment-file handling
        
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
            prompt = {'Enter Number of force maps','How would you like to name the Experiment?'};
            dlgtitle = 'Experiment Layout';
            dims = [1 35];
            definput = {'10','YourExperimentName'};
            answer = inputdlg(prompt,dlgtitle,dims,definput);
            
            % ReferenceSlopeMaps
            answerRef = questdlg('Are there any force maps from reference material (glass, mica)?','Reference Slope Data','Yes','No','No');
            
            if isequal(answerRef,'Yes')
                prompt = 'Enter Number of reference force maps';
                dlgtitle = 'Experiment Layout';
                dims = [1 35];
                definput = {'1'};
                answerRefN = inputdlg(prompt,dlgtitle,dims,definput);
                NRef = str2double(answerRefN{1});
            else
                NRef = 0;
            end
            
            % Set HostOS and HostName properties
            FullOS = computer;
            OS = FullOS(1:3);
            if isequal(OS,'PCW')
                Host = getenv('COMPUTERNAME');
            elseif isequal(OS,'GLN')
                Host = getenv('HOSTNAME');
            elseif isequal(OS,'MAC')
                Host = getenv('HOSTNAME');
            end
            obj.HostOS = OS;
            obj.HostName = Host;
            
            obj.NumFiles = str2double(answer{1});
            obj.ExperimentName = answer{2};
            current = what();
            ParentFolder = uigetdir(current.path,'Choose a Folder where the Experiment is to be saved');
            mkdir(ParentFolder,obj.ExperimentName);
            obj.ExperimentFolder = fullfile(ParentFolder,obj.ExperimentName,filesep);
            N = obj.NumFiles;
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
            
            while length(MapFullFile) < N + NRef
                Title = sprintf('Choose one or more REFERENCE .jpk-force-map files. %i/%i',length(MapFullFile)-N,NRef);
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
            FM = cell(N,1);
            SPM = cell(N,1);
            ExperimentName = obj.ExperimentName;
            ExperimentFolder = obj.ExperimentFolder;
            if contains(struct2array(ver), 'Parallel Computing Toolbox')
                parfor i=1:N
                    % for i=1:N Debugging
                    if WhichFiles == 2 || WhichFiles == 0
                        TempID = sprintf('%s-%i',ExperimentName,i);
                        FM{i} = ForceMap(MapFullFile{i},ExperimentFolder,TempID);
                    elseif WhichFiles == 1 || WhichFiles == 0
                        SPM{i} = SurfacePotentialMap();
                    end
                end
            else
                for i=1:N
                    % for i=1:N Debugging
                    if WhichFiles == 2 || WhichFiles == 0
                        TempID = sprintf('%s-%i',ExperimentName,i);
                        FM{i} = ForceMap(MapFullFile{i},ExperimentFolder,TempID);
                    elseif WhichFiles == 1 || WhichFiles == 0
                        SPM{i} = SurfacePotentialMap();
                    end
                end
            end
            
            % Assign the objects created in the parfor loop to the
            % Experiment object
            obj.FM = FM;
            obj.SPM = SPM;
            for i=1:N
                if WhichFiles == 2 || WhichFiles == 0
                    obj.ForceMapFolders{i} = obj.FM{i}.Folder;
                    obj.ForceMapNames{i} = obj.FM{i}.Name;
                elseif WhichFiles == 1 || WhichFiles == 0
                    obj.SurfacePotentialMapFolders{i} = obj.SPM{i}.Folder;
                    obj.SurfacePotentialMapNames{i} = obj.SPM{i}.Name;
                end
            end
            
            for i=1:NRef
                if WhichFiles == 2 || WhichFiles == 0
                    obj.RefFM{i} = ForceMap(MapFullFile{N+i},obj.ExperimentFolder);
                end
            end
            obj.FMFlag.FibrilAnalysis = zeros(N,1);
            obj.FMFlag.ForceMapAnalysis = zeros(N,1);
            obj.FMFlag.Preprocessed = zeros(N,1);
            obj.FMFlag.Grouping = 0;
            obj.SPMFlag.FibrilAnalysis = zeros(N,1);
            obj.SPMFlag.Grouping = 0;
            obj.CantileverTipFlag = 0;
            
%             if WhichFiles == 2 || WhichFiles == 0
%                 obj.grouping_force_map();
%             elseif WhichFiles == 1 || WhichFiles == 0
%                 obj.grouping_surface_potential_map();
%             end
            
            Temp = load('DropoutNetFinal.mat');
            obj.DropoutNet = Temp.MC14_Drop;
            Temp2 = load('CP_CNN_Final.mat');
            obj.CP_CNN = Temp2.CNN;
            
            obj.check_for_new_host();
            
            obj.save_experiment();
        end
        
        function Out = add_data(obj)
            
            warning('Note that for this function to work properly you need to assign your Experiment in workspace to itself e.g. ">> E = E.add_data" otherwise the updated Experiment will be stored in the temporary variable "ans" and inevitably overwritten at some point!')
            
            % create save copy to restore if function produces errors
            SaveCopy = obj.copy_experiment;
            
            try
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
            prompt = {'Enter Number of additional force maps'};
            dlgtitle = 'How many files to add?';
            dims = [1 35];
            definput = {'5'};
            answer = inputdlg(prompt,dlgtitle,dims,definput);
            
            NOld = SaveCopy.NumFiles;
            NumNewFiles = str2double(answer{1});
            SaveCopy.NumFiles = NOld + NumNewFiles;
            N = NumNewFiles;
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
            
            % Load maps into Experiment with parfor-loop
            FM = cell(SaveCopy.NumFiles,1);
            SPM = cell(SaveCopy.NumFiles,1);
            FM(1:NOld) = SaveCopy.FM;
            SPM(1:NOld) = SaveCopy.SPM;
            ExperimentName = obj.ExperimentName;
            ExperimentFolder = obj.ExperimentFolder;
            if contains(struct2array(ver), 'Parallel Computing Toolbox')
                parfor i=1:N
                    if WhichFiles == 2 || WhichFiles == 0
                        TempID = sprintf('%s-%i',ExperimentName,NOld+i);
                        FM{NOld+i} = ForceMap(MapFullFile{i},ExperimentFolder,TempID);
                    elseif WhichFiles == 1 || WhichFiles == 0
                        SPM{NOld+i} = SurfacePotentialMap();
                    end
                end
            else
                for i=1:N
                    if WhichFiles == 2 || WhichFiles == 0
                        TempID = sprintf('%s-%i',ExperimentName,NOld+i);
                        FM{NOld+i} = ForceMap(MapFullFile{i},ExperimentFolder,TempID);
                    elseif WhichFiles == 1 || WhichFiles == 0
                        SPM{NOld+i} = SurfacePotentialMap();
                    end
                end
                
            end
            % Assign the objects created in the parfor loop to the
            % Experiment object
            SaveCopy.FM = FM;
            SaveCopy.SPM = SPM;
            for i=1:N
                if WhichFiles == 2 || WhichFiles == 0
                    SaveCopy.ForceMapFolders{NOld+i} = SaveCopy.FM{NOld+i}.Folder;
                    SaveCopy.ForceMapNames{NOld+i} = SaveCopy.FM{NOld+i}.Name;
                elseif WhichFiles == 1 || WhichFiles == 0
                    SaveCopy.SurfacePotentialMapFolders{NOld+i} = SaveCopy.SPM{NOld+i}.Folder;
                    SaveCopy.SurfacePotentialMapNames{NOld+i} = SaveCopy.SPM{NOld+i}.Name;
                end
            end
            
            SaveCopy.FMFlag.FibrilAnalysis(NOld+1:NOld+N) = zeros(N,1);
            SaveCopy.FMFlag.ForceMapAnalysis(NOld+1:NOld+N) = zeros(N,1);
            SaveCopy.FMFlag.Preprocessed(NOld+1:NOld+N) = zeros(N,1);
            SaveCopy.FMFlag.Grouping = 0;
            SaveCopy.SPMFlag.FibrilAnalysis(NOld+1:NOld+N) = zeros(N,1);
            SaveCopy.SPMFlag.Grouping = 0;
            
%             if WhichFiles == 2 || WhichFiles == 0
%                 obj.grouping_force_map();
%             elseif WhichFiles == 1 || WhichFiles == 0
%                 obj.grouping_surface_potential_map();
%             end

            Out = SaveCopy;
            Out.save_experiment
            warning('Did you read the warning above?');
            catch ME
                disp('data adding failed. restored original experiment object')
                fclose('all')
                cd(obj.ExperimentFolder)
                
                % Check for and remove any temporary folders that were
                % created during failed file-add attempt
                CheckDir = dir('Temp*');
                if 0<length(CheckDir)
                    for i=1:length(CheckDir)
                        if CheckDir(i).isdir
                            rmdir(CheckDir(i).name)
                        end
                    end
                end
                
                Out = obj;
            end
            
        end
        
        function load_data(obj)
            for i=1:obj.NumFiles
                obj.FM{i} = ForceMap(obj.ForceMapFolders{i},obj.ForceMapNames{i});
                obj.SPM{i} = SurfacePotentialMap(obj.SurfacePotentialMapFolders{i},obj.SurfacePotentialMapNames{i});
            end
        end
        
        function save_data(obj)
            disp('saving');
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
            disp('saving...');
            save(savename,'obj','-v7.3')
            cd(current.path)
            savemsg = sprintf('Changes to Experiment %s saved to %s',obj.ExperimentName,obj.ExperimentFolder);
            disp(savemsg);
        end
        
        function ExperimentCopy = copy_experiment(obj)
            % ExperimentCopy = copy_experiment(obj)
            %
            % makes a proper copy of the object, so that also contained
            % handle objects, such as ForceMaps, SurfacePotentialMaps etc.
            % are copied and not only referenced to
            
            ExperimentCopy = obj.copy;
            for i=1:obj.NumFiles
                if i<=length(obj.FM)
                    MCFM = metaclass(obj.FM{i});
                end
                if i<=length(obj.SPM)
                    MCSPM = metaclass(obj.SPM{i});
                end
                if i<=length(obj.FM) && ~isempty(MCFM.SuperclassList) && isequal(MCFM.SuperclassList.Name,'matlab.mixin.Copyable')
                    ExperimentCopy.FM{i} = obj.FM{i}.copy;
                end
                if i<=length(obj.SPM) && ~isempty(MCSPM.SuperclassList) && isequal(MCSPM.SuperclassList.Name,'matlab.mixin.Copyable')
                    ExperimentCopy.SPM{i} = obj.SPM{i}.copy;
                end
            end
        end
    end
    methods(Static)
        % Static methods related with Experiment-file handling
        
        function E = load()
            % E = load()
            %
            % recommended way of loading an existing Experiment() from its
            % folder. Checks, if path has changed and adjusts object
            % properties. Also checks if running on different system,
            % and, if so, updating object properties and setting
            % CPFlag.CNNOpt = 0
            
            [File,Path] = uigetfile('*.mat','Choose Experiment .mat from folder');
            Fullfile = fullfile(Path,File);
            disp('Loading Experiment... this can take a while for larger Experiments')
            load(Fullfile);
            
            E = obj;
            clear obj
            
            E.check_for_new_host();
            FMFolder = fullfile(Path,filesep,'ForceData');
            for i=1:E.NumFiles
                if ~isempty(E.FM{i})
                    E.FM{i}.check_for_new_host();
                    E.FM{i}.Folder = FMFolder;
                    E.ForceMapFolders{i} = FMFolder;
                end
            end
            
            E.ExperimentFolder = Path;
            
        end
        
    end
    
    methods
        % methods for sequential data analysis mostly looping over child-classes methods
        
        function preprocessing(obj)
            % preprocessing(obj)
            % 
            % bare minimum of preprocessing steps to prepare data for
            % custom further data processing. (.base_and_tilt())
            
            h = waitbar(0,'setting up','Units','normalized','Position',[0.4 0.3 0.2 0.1]);
            NLoop = length(obj.ForceMapNames);
            if sum(obj.FMFlag.Preprocessed) >= 1
                KeepFlagged = questdlg(sprintf('Some maps have been processed already.\nDo you want to skip them and keep old results?'),...
                    'Processing Options',...
                    'Yes',...
                    'No',...
                    'No');
            else
                KeepFlagged = 'No';
            end
            
            for i=1:NLoop
                if isequal(KeepFlagged,'Yes') && obj.FMFlag.Preprocessed(i) == 1
                    continue
                end
                obj.FM{i}.create_and_level_height_map();
            end
            
            
            % Main loop for base and tilt
            for i=1:NLoop
                if isequal(KeepFlagged,'Yes') && obj.FMFlag.Preprocessed(i) == 1
                    continue
                end
                waitbar(i/NLoop,h,sprintf('Preprocessing ForceMap %i/%i\nFitting Base Lines',i,NLoop));
                obj.FM{i}.base_and_tilt('linear');
                waitbar(i/NLoop,h,sprintf('Preprocessing ForceMap %i/%i\nWrapping Up And Saving',i,NLoop));
                
                obj.FM{i}.save();
                obj.FMFlag.Preprocessed(i) = 1;
            end
            
            obj.save_experiment;
            
            close(h);
        end
        
        function force_map_analysis_fibril(obj,CPOption,EModOption)
            % force_map_analysis_fibril(obj,CPOption,EModOption)
            %
            % CPOption = 'Fast' ...(Default) contact point estimation through single
            % pass through CNN
            % CPOption = 'Dropout' ... contact point estimation through
            % averaging over multiple passes through monte carlo dropout
            % net. NPasses (Default=100) times slower than 'Fast'
            % CPOption = 'Zoom' ... CNN-based method 
            % CPOption = 'ZoomDropout' ... CNN-based method (still in development)
            %%%%%% RECOMMENDED %%%%%%
            % CPOption = 'Zoomsweep' ... CNN-based method
            %%%%%% RECOMMENDED %%%%%%
            % CPOption = 'Old' ... old method for contact point estimation
            % CPOption = 'RoV' ... RoV method for contact point estimation
            % CPOption = 'GoF' ... GoF method for contact point estimation
            % CPOption = 'Combo' ... RoV and GoF combined method for contact point estimation
            % CPOption = 'Manual' ... go through manual CP determination for contact point estimation
            %
            % EModOption = 'Hertz' ... E-Modulus calculation through Hertz-Sneddon
            % method
            % EModOption = 'Oliver' ... E-Modulus calculation through
            % Oliver-Pharr-like method (O. Andriotis 2014)
            
            h = waitbar(0,'setting up','Units','normalized','Position',[0.4 0.3 0.2 0.1]);
            NLoop = length(obj.ForceMapNames);
            if sum(obj.FMFlag.FibrilAnalysis) >= 1
                KeepFlagged = questdlg(sprintf('Some maps have been processed already.\nDo you want to skip them and keep old results?'),...
                    'Processing Options',...
                    'Yes',...
                    'No',...
                    'No');
            else
                KeepFlagged = 'No';
            end
            
            if obj.CantileverTipFlag == 1
                KeepTip = questdlg(sprintf('There already exists data from a deconvoluted tip\nDo you want to skip tip deconvolution and keep old tip data?'),...
                    'Processing Options',...
                    'Yes',...
                    'No',...
                    'No');
            else
                KeepTip = 'No';
            end
            
            if isequal(KeepTip,'No')
                waitbar(0,h,'deconvoluting cantilever tip...')
                obj.deconvolute_cantilever_tip;
            elseif isequal(KeepTip,'Yes')
            end
            
            % Preprocessing everything, that needs user input
            answer = questdlg('Do you want to skip manual exclusion of problematic areas?',...
                'Manual Exclusion',...
                'Yes',...
                'No','No');
            for i=1:NLoop
                if isequal(KeepFlagged,'Yes') && obj.FMFlag.FibrilAnalysis(i) == 1
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
                if isequal(KeepFlagged,'Yes') && obj.FMFlag.FibrilAnalysis(i) == 1
                    continue
                end
                waitbar(i/NLoop,h,sprintf('Processing Fibril %i/%i\nFitting Base Line',i,NLoop));
                obj.FM{i}.base_and_tilt('linear');
                obj.FM{i}.calculate_fib_diam();
                waitbar(i/NLoop,h,sprintf('Processing Fibril %i/%i\nFinding Contact Point',i,NLoop));
                
                % contact point estimation happens here
                obj.cp_option_converter(CPOption,i);
                
                waitbar(i/NLoop,h,sprintf('Processing Fibril %i/%i\nCalculating E-Modulus',i,NLoop));
                if isequal(lower(EModOption),'hertz')
                    obj.FM{i}.calculate_e_mod_hertz(CPOption,'parabolic',1);
                else
                    obj.FM{i}.calculate_e_mod_oliverpharr(obj.CantileverTip.ProjArea,0.75);
                end
                waitbar(i/NLoop,h,sprintf('Processing Fibril %i/%i\nWrapping Up And Saving',i,NLoop));
                
                if i > 1
                    close(Fig{i-1})
                end
                Fig{i} = obj.FM{i}.show_analyzed_fibril();
                obj.FM{i}.save();
                obj.FMFlag.FibrilAnalysis(i) = 1;
            end
            
            % Assign the Apex curves EMod and exclude +-2.5*IQR and curves
            % from ExclMask
            for i=1:NLoop
                if isequal(lower(EModOption(1:5)),'hertz')
                    EMods = obj.FM{i}.EModHertz;
                elseif isequal(lower(EModOption(1:6)),'oliver')
                    EMods = obj.FM{i}.EModOliverPharr;
                else
                    EMods = obj.FM{i}.EModOliverPharr;
                end
                obj.EMod.Apex(i,1:length(obj.FM{i}.RectApexIndex)) = EMods(obj.FM{i}.RectApexIndex);
                for j=1:length(obj.FM{i}.RectApexIndex)
                    if obj.EMod.Apex(i,j) > (nanmedian(obj.EMod.Apex(i,:))+2.5*iqr(obj.EMod.Apex(i,:))) || ...
                            obj.EMod.Apex(i,j) < (nanmedian(obj.EMod.Apex(i,:))-2.5*iqr(obj.EMod.Apex(i,:))) || ...
                            obj.FM{i}.ExclMask(obj.FM{i}.List2Map(obj.FM{i}.RectApexIndex(j),1),obj.FM{i}.List2Map(obj.FM{i}.RectApexIndex(j),2)) == 0
                        obj.EMod.Apex(i,j) = NaN;
                    elseif obj.EMod.Apex(i,j) < 0
                        obj.EMod.Apex(i,j) = NaN;
                    end
                end
            end
            
            obj.EMod.Mean = nanmean(obj.EMod.Apex,2);
            obj.EMod.STD = nanstd(obj.EMod.Apex,[],2);
            obj.save_experiment;
            
            close(h);
        end
        
        function force_map_analysis_general(obj,CPOption,EModOption)
            % force_map_analysis_general(obj,CPOption,EModOption)
            %
            % CPOption = 'Fast' ...(Default) contact point estimation through single
            % pass through CNN
            % CPOption = 'Dropout' ... contact point estimation through
            % averaging over multiple passes through monte carlo dropout
            % net. NPasses (Default=100) times slower than 'Fast'
            % CPOption = 'Zoom' ... CNN-based method (still in development)
            % CPOption = 'ZoomDrop' ... CNN-based method (still in development)
            %%%%%% RECOMMENDED %%%%%%
            % CPOption = 'Zoomsweep' ... CNN-based method
            %%%%%% RECOMMENDED %%%%%%
            % CPOption = 'Old' ... old method for contact point estimation
            % CPOption = 'RoV' ... RoV method for contact point estimation
            % CPOption = 'GoF' ... GoF method for contact point estimation
            % CPOption = 'Combo' ... RoV and GoF combined method for contact point estimation
            % CPOption = 'Manual' ... go through manual CP determination for contact point estimation
            %
            % EModOption = 'Hertz' ... E-Modulus calculation through Hertz-Sneddon
            % method
            % EModOption = 'Oliver' ... E-Modulus calculation through
            % Oliver-Pharr-like method (O. Andriotis 2014)
            
            h = waitbar(0,'setting up','Units','normalized','Position',[0.4 0.3 0.2 0.1]);
            NLoop = length(obj.ForceMapNames);
            if sum(obj.FMFlag.ForceMapAnalysis) >= 1
                KeepFlagged = questdlg(sprintf('Some maps have been processed already.\nDo you want to skip them and keep old results?'),...
                    'Processing Options',...
                    'Yes',...
                    'No',...
                    'No');
            else
                KeepFlagged = 'No';
            end
            if isequal(lower(EModOption),'oliver')
                if obj.CantileverTipFlag == 1
                    KeepTip = questdlg(sprintf('There already exists data from a deconvoluted tip\nDo you want to skip tip deconvolution and keep old tip data?'),...
                        'Processing Options',...
                        'Yes',...
                        'No',...
                        'No');
                else
                    KeepTip = 'No';
                end
                
                if isequal(KeepTip,'No')
                    waitbar(0,h,'deconvoluting cantilever tip...')
                    obj.deconvolute_cantilever_tip;
                elseif isequal(KeepTip,'Yes')
                end
            end
            % Preprocessing everything, that needs user input
            answer = questdlg('Do you want to skip manual exclusion of problematic areas?',...
                'Manual Exclusion',...
                'Yes',...
                'No','No');
            for i=1:NLoop
                if isequal(KeepFlagged,'Yes') && obj.FMFlag.ForceMapAnalysis(i) == 1
                    continue
                end
                obj.FM{i}.create_and_level_height_map();
                if isequal(answer,'Yes')
                    continue
                end
                obj.FM{i}.manual_exclusion();
            end
            
            % preprocess reference ForceMaps
            if isequal(lower(EModOption),'oliver')
                ExternalRefSlope = true;
                waitbar(0,h,'Processing reference ForceMaps')
                for i=1:length(obj.RefFM)
                    obj.RefFM{i}.base_and_tilt('linear');
                    obj.cp_option_converter(CPOption,i,ExternalRefSlope);
                    obj.RefFM{i}.calculate_reference_slope;
                end
            end
            
            % Main loop for contact point estimation and E-Modulus calculation
            for i=1:NLoop
                if isequal(KeepFlagged,'Yes') && obj.FMFlag.ForceMapAnalysis(i) == 1
                    continue
                end
                waitbar(i/NLoop,h,sprintf('Processing ForceMap %i/%i\nFitting Base Line',i,NLoop));
                obj.FM{i}.base_and_tilt('linear');
                waitbar(i/NLoop,h,sprintf('Processing ForceMap %i/%i\nFinding Contact Point',i,NLoop));
                
                % contact point estimation happens here
                obj.cp_option_converter(CPOption,i);
                
                waitbar(i/NLoop,h,sprintf('Processing ForceMap %i/%i\nCalculating E-Modulus',i,NLoop));
                if isequal(lower(EModOption),'hertz')
                    obj.FM{i}.calculate_e_mod_hertz(CPOption,'parabolic',1);
                else
                    obj.FM{i}.RefSlope = obj.RefFM{1}.RefSlope;
                    obj.FM{i}.calculate_e_mod_oliverpharr(obj.CantileverTip.ProjArea,0.75,ExternalRefSlope);
                end
                waitbar(i/NLoop,h,sprintf('Processing ForceMap %i/%i\nWrapping Up And Saving',i,NLoop));
                
                obj.FM{i}.save();
                obj.FMFlag.ForceMapAnalysis(i) = 1;
            end
            
            obj.save_experiment;
            
            close(h);
        end
        
        function surface_potential_analysis_fibril(obj)
            
        end
        
        function surface_potential_analysis_general(obj)
            % WORK IN PROGRESS
        end
        
        function cp_option_converter(obj,CPOption,i,RefFM)
            % cp_option_converter(CPOption,i,RefFM)
            %
            % aux-function for CPOption choice
            
            NumPasses = 20;
            
            if nargin < 4
                RefFM = false;
            end
            if RefFM == false
                if isequal(lower(CPOption),'rov')
                    obj.FM{i}.estimate_cp_rov();
                end
                if isequal(lower(CPOption),'gof')
                    obj.FM{i}.estimate_cp_gof();
                end
                if isequal(lower(CPOption),'old')
                    obj.FM{i}.estimate_cp_old();
                end
                if isequal(lower(CPOption),'combo')
                    if obj.FM{i}.CPFlag.RoV == 0
                        obj.FM{i}.estimate_cp_rov();
                    end
                    if obj.FM{i}.CPFlag.GoF == 0
                        obj.FM{i}.estimate_cp_gof();
                    end
                    obj.FM{i}.estimate_cp_combined();
                end
                if isequal(lower(CPOption),'manual')
                    obj.FM{i}.estimate_cp_manually;
                end
                if isequal(lower(CPOption),'fast')
                    obj.FM{i}.estimate_cp_cnn(obj.CP_CNN,'Fast');
                end
                if isequal(lower(CPOption),'dropout')
                    obj.FM{i}.estimate_cp_cnn(obj.DropoutNet,'Dropout',NumPasses);
                end
                if isequal(lower(CPOption),'zoom')
                    obj.FM{i}.estimate_cp_cnn(obj.CP_CNN,'Zoom');
                end
                if isequal(lower(CPOption),'zoomdropout')
                    obj.FM{i}.estimate_cp_cnn(obj.CP_CNN,'zoomdropout',NumPasses);
                end
                if isequal(lower(CPOption),'zoomsweep')
                    obj.FM{i}.estimate_cp_cnn(obj.CP_CNN,'Zoomsweep',NumPasses);
                end
            elseif RefFM == true
                if isequal(lower(CPOption),'rov')
                    obj.RefFM{i}.estimate_cp_rov();
                end
                if isequal(lower(CPOption),'gof')
                    obj.RefFM{i}.estimate_cp_gof();
                end
                if isequal(lower(CPOption),'old')
                    obj.RefFM{i}.estimate_cp_old();
                end
                if isequal(lower(CPOption),'combo')
                    obj.RefFM{i}.estimate_cp_rov();
                    obj.RefFM{i}.estimate_cp_gof();
                    obj.RefFM{i}.estimate_cp_combined();
                end
                if isequal(lower(CPOption),'manual')
                    obj.RefFM{i}.estimate_cp_manually;
                end
                if isequal(lower(CPOption),'fast')
                    obj.RefFM{i}.estimate_cp_cnn(obj.CP_CNN,'Fast');
                end
                if isequal(lower(CPOption),'dropout')
                    obj.RefFM{i}.estimate_cp_cnn(obj.DropoutNet,'Dropout',NumPasses);
                end
                if isequal(lower(CPOption),'zoom')
                    obj.RefFM{i}.estimate_cp_cnn(obj.CP_CNN,'Zoom');
                end
                if isequal(lower(CPOption),'zoomdropout')
                    obj.RefFM{i}.estimate_cp_cnn(obj.CP_CNN,'zoomdropout',NumPasses);
                end
                if isequal(lower(CPOption),'zoomsweep')
                    obj.RefFM{i}.estimate_cp_cnn(obj.CP_CNN,'Zoomsweep',NumPasses);
                end
            end
                
        end
            
        function Results = results_readout(obj)
            % result_readout(obj)
            %
            % Interactive function to specify how and which results to
            % return
            %
            % Under Development!!!!!
            
            if obj.FMFlag.Grouping == 0
                obj.grouping_force_map;
            end
            
            N = length(obj.GroupFM);
            
            for i=1:N
                Index = regexp(obj.GroupFM(i).Name,'\w');
                FieldName = obj.GroupFM(i).Name(Index);
                Results.(FieldName) = [];
                
                clear Index
            end
            
        end
        
        function create_moving_dot_gif_emod_vs_surfpot(obj,NFrames,NStartFrames,NEndFrames)
            % create_moving_dot_gif_emod_vs_surfpot(obj,NFrames,FreezeFrames)
            %
            % This gif-creator has been written for a specific use case and
            % has to be manually recoded to fit other Experiments. Videos
            % generated run at 30 FPS so e.g. NFrames=60, NEndFrames=30
            % and NStartFrames=30 creates a 4 second video
            
            % Get Fibril Data
            X = obj.statistical_analysis_force_maps*1e-6;
            Y = obj.statistical_analysis_surface_potential*1e3;
            Y = Y';
            close all
            
            %
            DX1 = X(11:20) - X(1:10); 
            DX2 = X(31:40) - X(21:30);
            DY1 = Y(11:20) - Y(1:10);
            DY2 = Y(31:40) - Y(21:30);
            
            Steps = 0:1/(NFrames-1):1;
            
            % Create and and open .avi for write-in of frames
            FullFile = fullfile(obj.ExperimentFolder,filesep,'EModvsSurfPot_gif');
            
            Vid = VideoWriter(FullFile,'MPEG-4');
            Vid.open
            
            
            % Create figure and loop over NFrames to write into Vid
            Fig = figure('Color','white',...
                'Units','normalized','Position',[0 0 1 1]);
            for i=1:NFrames
                hold off
                XTemp1 = X(1:10) + Steps(i)*DX1;
                XTemp2 = X(21:30) + Steps(i)*DX2;
                YTemp1 = Y(1:10) + Steps(i)*DY1;
                YTemp2 = Y(21:30) + Steps(i)*DY2;
                MeanX1 = mean(XTemp1);
                ErrX1 = std(XTemp1)/sqrt(10);
                MeanX2 = mean(XTemp2);
                ErrX2 = std(XTemp2)/sqrt(10);
                MeanY1 = mean(YTemp1);
                ErrY1 = std(YTemp1)/sqrt(10);
                MeanY2 = mean(YTemp2);
                ErrY2 = std(YTemp2)/sqrt(10);
                plot(XTemp1,YTemp1,'square','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',10)
                ax = Fig.CurrentAxes;
                xlim([min(X)-0.1*range(X) max(X)+0.1*range(X)])
                ylim([min(Y)-0.1*range(Y) max(Y)+0.1*range(Y)])
                hold on
                errorbar(MeanX1,MeanY1,-ErrY1,ErrY1,-ErrX1,ErrX1,'bsquare','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',20)
                plot(XTemp2,YTemp2,'diamond','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',10)
                errorbar(MeanX2,MeanY2,-ErrY1,ErrY1,-ErrX1,ErrX1,'rdiamond','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',20)
                legend({'Control Group','Mean Control +- SE','MGO Group','Mean MGO +- SE'},'FontSize',22)
                xlabel('Indentation Modulus [MPa]','FontSize',22)
                ylabel('Rel. Surface Potential [mV]','FontSize',22)
                ax.FontSize = 22;
                title('...','FontSize',42)
                
                if i==1
                    title('Before PBS/MGO','FontSize',42)
                    Frame = getframe(Fig);
                    for j=1:NStartFrames
                        Vid.writeVideo(Frame)
                    end
                    title('...','FontSize',42)
                end
                
                Frame = getframe(Fig);
                
                Vid.writeVideo(Frame);
            end
            
            for i=1:NEndFrames
                % Lingers on the last frame for NEndFrames/30 seconds
                title('After PBS/MGO','FontSize',42)
                Frame = getframe(Fig);
                Vid.writeVideo(Frame)
            end
            
            Vid.close
            close(Fig)
        end
        
        function min_batch(obj)
            
            for ii=1:obj.NumFiles
               obj.FM{ii}.base_and_tilt('linear');
               obj.FM{ii}.min_force; 
            end
            obj.save_experiment
        end
        
        function SMFS_print(obj)
         
            % Change into the Folder of Interest
            cd(obj.ExperimentFolder) % Move into the folder 
            % Create folders for saving the produced figures
            foldername='FM_raw_Fig';    % Defines the folder name
            mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
            currpath=fullfile(obj.ExperimentFolder,foldername);
            cd(currpath); 
            
            % Loop over the imported force maps
            for ii=1:obj.NumFiles
            %for ii=3:5 % Debugging
               % Command window output
               sprintf('Force Map No. %d of %d',ii,obj.NumFiles) % Gives current Force Map Position
               % Run the chosen functions
               obj.FM{ii}.estimate_cp_hardsurface
               obj.FM{ii}.fc_print;     
               obj.save_experiment;        % Save immediately after each force curve
            end    
        end
        
        function SMFS_print_sort(obj,StartDate,EndDate)
         % Comment: Date format is: 'YYYY.MM.DD'
            
         % Change into the Folder of Interest
            cd(obj.ExperimentFolder) % Move into the folder 
            if nargin<2
                StartDate='0000.00.00';
                EndDate='2999.00.00';
            elseif nargin<3
                EndDate='2999.00.00';
            end
            % Loop over the imported force maps
             for ii=1:obj.NumFiles
                % Needed function
                obj.FM{ii}.fc_chipprop
         
                if ~obj.FM{ii}.FlagPrintSort
                % if obj.FM{ii}.FlagPrintSort    
                    continue
                end
                % Remove the dots in the dates               
                FMDate=split(obj.FM{ii}.Date,'.');
                StartDateSplit=split(StartDate,'.');
                EndDateSplit=split(EndDate,'.');
                % if conditions for the time selection (start date and end date)              
                if ~(str2num(FMDate{1}) > str2num(StartDateSplit{1}) || ... % Verifies if the year of the FM object Date is greater than the start date
                        (str2num(FMDate{1}) == str2num(StartDateSplit{1})... % Verifies if the year of the FM object Date is equal with the start date
                        && str2num(FMDate{2}) > str2num(StartDateSplit{2})) || ...  % Verifies if the month of the FM object Date is greater than the start date
                    (str2num(FMDate{1}) == str2num(StartDateSplit{1})...    
                    &&(str2num(FMDate{2}) == str2num(StartDateSplit{2}))... % Verifies if the month of the FM object Date is equal than the start date
                        && (str2num(FMDate{3}) >= str2num(StartDateSplit{3})))) % Verifies if the day of the FM object Date is greater than or equal with the start date
                    continue
                elseif ~(str2num(FMDate{1}) < str2num(EndDateSplit{1}) || ... 
                        (str2num(FMDate{1}) == str2num(EndDateSplit{1})...
                        && str2num(FMDate{2}) < str2num(EndDateSplit{2})) || ...
                    (str2num(FMDate{1}) == str2num(EndDateSplit{1})...    
                    &&(str2num(FMDate{2}) == str2num(EndDateSplit{2}))...
                        && (str2num(FMDate{3}) <= str2num(EndDateSplit{3})))) 
                    continue
                end  
                % Define variables for the folder name
                VelocityConvert=num2str(obj.FM{ii}.Velocity*1e+9); % Convert into nm
                StartDateMod=strrep(StartDate,'.','');
                EndDateMod=strrep(EndDate,'.','');
                foldername=append('FM_',VelocityConvert,obj.FM{ii}.Substrate,'_',obj.FM{ii}.EnvCond,'_','_',StartDateMod,'-',EndDateMod); % Defines the folder name
                warning('off','all');
                mkdir(foldername);
                warning('on','all');
                cd(foldername)         
               % Run the chosen functions
               obj.FM{ii}.estimate_cp_hardsurface      
               obj.FM{ii}.fc_print
               cd(obj.ExperimentFolder) % Move into the folder 
               obj.FM{ii}.FlagPrintSort=1;
            end 
            %obj.save_experiment        % Save immediately after each force curve
        end
        
        function SMFS_selection(obj)
         
            % Change into the Folder of Interest
            cd(obj.ExperimentFolder) % Move into the folder 
            % Create folders for saving the produced figures
            foldername='FM_Fig';    % Defines the folder name
            mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
            currpath=fullfile(obj.ExperimentFolder,foldername);
            cd(currpath); 
            
            % Loop over the imported force maps
            for ii=1:obj.NumFiles
            %for ii=3:5 % Debugging
               % Command window output
               sprintf('Force Map No. %d of %d',ii,obj.NumFiles) % Gives current Force Map Position
               % Run the chosen functions
               obj.FM{ii}.estimate_cp_hardsurface
               obj.FM{ii}.fc_selection;     
               %obj.save_experiment;        % Save immediately after each force curve
            end    
        end
        
        
    end
    methods
        %%%   WARNING! %%%
        % The following methods were programmed for specific use cases
        % and are yet to be generalized! However, you can of course adjust
        % them for your own use, though its probably easier to just take
        % the processed raw data from your Experiment and do your own
        % statistics
        
        function [DataMeansOP,DataMeansHS,DataOP,DataHS] = statistical_analysis_force_maps(obj)
            % Basic statistical analysis of the force map experiment. First get the grouping
            % of data from the user and then perform several tests, plots
            % etc.
            
            warning('The following methods were programmed for specific use cases and are yet to be generalized! However, you can of course adjust them for your own use, though its probably easier to just take the processed raw data from your Experiment and do your own statistics')
            
            N = obj.NumFiles;
            obj.grouping_force_map();
            
            % Ask user which groups are to be compared
            prompt = 'Specialize which groups are to be tested against which';
            definput = {'1 2 ; 3 4'};
            answer = inputdlg(prompt,'Pairings',[1 50],definput);
            
            % Define Test Matrix
            TestMat = [str2num(answer{1})];
            
            % Write Data into local variables and
            % replace extreme outlier values with NaN
            % also replace data points that lie on the exclusion mask with
            % NaN
            for i=1:N
                DataOP(i,:) = obj.FM{i}.EModOliverPharr(obj.FM{i}.RectApexIndex);
                DataHS(i,:) = obj.FM{i}.EModHertz(obj.FM{i}.RectApexIndex);
                OutliersOP = isoutlier(DataOP(i,:));
                OutliersHS = isoutlier(DataHS(i,:));
                for j=1:length(obj.FM{i}.RectApexIndex)
                    if DataOP(i,j) > (nanmedian(DataOP(i,:))+2.5*iqr(DataOP(i,:))) || ...
                            DataOP(i,j) < (nanmedian(DataOP(i,:))-2.5*iqr(DataOP(i,:))) || ...
                            obj.FM{i}.ExclMask(obj.FM{i}.List2Map(obj.FM{i}.RectApexIndex(j),1),obj.FM{i}.List2Map(obj.FM{i}.RectApexIndex(j),2)) == 0 ||...
                            OutliersOP(j) == 1
                        DataOP(i,j) = NaN;
                    elseif DataOP(i,j) < 0
                        DataOP(i,j) = NaN;
                    end
                    if DataHS(i,j) > (nanmedian(DataHS(i,:))+2.5*iqr(DataHS(i,:))) || ...
                            DataHS(i,j) < (nanmedian(DataHS(i,:))-2.5*iqr(DataHS(i,:))) || ...
                            obj.FM{i}.ExclMask(obj.FM{i}.List2Map(obj.FM{i}.RectApexIndex(j),1),obj.FM{i}.List2Map(obj.FM{i}.RectApexIndex(j),2)) == 0 ||...
                            OutliersHS(j) == 1
                        DataHS(i,j) = NaN;
                    elseif DataHS(i,j) < 0
                        DataHS(i,j) = NaN;
                    end
                end
            end
            
            DataMeansOP = nanmean(DataOP,2);
            DataMeansHS = nanmean(DataHS,2);
            
            figure('Name','OliverPharr vs HertzSneddon','Color','w');
            plot(DataMeansHS,DataMeansOP,'bO')
            legend(sprintf('E-Mod Hertz vs. Oliver-Pharr (R=%.3f)',corr(DataMeansHS,DataMeansOP)))
%             xlim([0,N+1])
            xlabel('E-Mod Hertz [Pa]')
            ylabel('E-Mod Oliver-Pharr [Pa]')
            
            % loop over all rows of the Test Matrix, doing paired ttests
            for i=1:size(TestMat,1)
                % Statistics for Oliver-Pharr Method
                [hOP(i),pOP(i)] = ...
                    ttest(DataMeansOP(obj.GroupFM(TestMat(i,2)).Indices),...
                    DataMeansOP(obj.GroupFM(TestMat(i,1)).Indices),'Tail','right');
                
                figure('Name','Paired Right Tailed T-Test','Units','normalized','Position',[0.2 0.2 0.5 0.5],'Color','w')
                boxplot([DataMeansOP(obj.GroupFM(TestMat(i,1)).Indices) DataMeansOP(obj.GroupFM(TestMat(i,2)).Indices)]*1e-6)
                title('Paired Right Tailed T-Test for Oliver-Pharr Method')
                xticklabels({obj.GroupFM(TestMat(i,1)).Name,obj.GroupFM(TestMat(i,2)).Name})
                xlabel('Test Group')
                ylabel('Indentation Modulus [MPa]')
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
                
%                 %Statistics for Hertz-Sneddon Method
%                 [hHS(i),pHS(i)] = ...
%                     ttest(DataMeansHS(obj.GroupFM(TestMat(i,2)).Indices),...
%                     DataMeansHS(obj.GroupFM(TestMat(i,1)).Indices),'Tail','right');
%                 
%                 figure('Name','Paired Right Tailed T-Test','Units','normalized','Position',[0.2 0.2 0.5 0.5],'Color','w')
%                 boxplot([DataMeansHS(obj.GroupFM(TestMat(i,1)).Indices) DataMeansHS(obj.GroupFM(TestMat(i,2)).Indices)])
%                 title('Paired Right Tailed T-Test for Hertz-Sneddon Method')
%                 xticklabels({obj.GroupFM(TestMat(i,1)).Name,obj.GroupFM(TestMat(i,2)).Name})
%                 xlabel('Test Group')
%                 ylabel('E-Mod Hertz-Sneddon[Pa]')
%                 DeltaMean = mean(DataMeansHS(obj.GroupFM(TestMat(i,2)).Indices)) - mean(DataMeansHS(obj.GroupFM(TestMat(i,1)).Indices));
%                 Sigma = std(DataMeansHS(obj.GroupFM(TestMat(i,2)).Indices) - DataMeansHS(obj.GroupFM(TestMat(i,1)).Indices));
%                 Beta = sampsizepwr('t',[0 Sigma],DeltaMean,[],length(obj.GroupFM(TestMat(i,2)).Indices));
%                 Stats = {sprintf('\\DeltaMean = %.2f MPa',DeltaMean*1e-6),...
%                     sprintf('P-Value = %.4f%',pHS(i)),...
%                     sprintf('Power \\beta = %.2f%%',Beta*100),...
%                     sprintf('Number of Specimen = %i',length(obj.GroupFM(TestMat(i,2)).Indices))};
%                 text(0.5,0.8,Stats,...
%                     'Units','normalized',...
%                     'FontSize',12,...
%                     'HorizontalAlignment','center')
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
            boxplot([DiffControlOP DiffMGOOP])
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
            
%             % For Hertz-Sneddon
%             [hHS,pHS,ciHS,statsHS] = ttest2(DiffMGOHS,DiffControlHS);
%             figure('Name','Two Sample T-Test','Units','normalized','Position',[0.2 0.2 0.5 0.5],'Color','w')
%             yyaxis left
%             boxplot([DiffControlHS DiffMGOHS])
%             ax = gca;
%             YLim = ax.YLim;
%             ylabel('Difference Before-After E-Mod [Pa]')
%             DeltaMean = mean(DiffMGOHS) - mean(DiffControlHS);
%             PooledSTD = statsHS.sd;
%             yyaxis right
%             errorbar(1.5,DeltaMean,ciHS(2)-DeltaMean,'O');
%             ylim(YLim)
%             xticks([1 1.5 2])
%             title('Two Sample T-Test for E-Mod Hertz-Sneddon Method')
%             ax = gca;
%             ax.TickLabelInterpreter = 'tex';
%             xticklabels({sprintf('%s - %s',obj.GroupFM(2).Name,obj.GroupFM(1).Name),...
%                 '\DeltaMean with CI',...
%                 sprintf('%s - %s',obj.GroupFM(4).Name,obj.GroupFM(3).Name)})
%             ylabel('Difference of Differences [Pa]')
%             Beta = sampsizepwr('t2',[mean(DiffControlHS) PooledSTD],mean(DiffMGOHS),[],length(DiffControlHS),'Ratio',length(DiffMGOHS)/length(DiffControlHS));
%             Stats = {sprintf('\\DeltaMean = %.2f MPa',DeltaMean*1e-6),...
%                 sprintf('P-Value = %.4f%',pHS),...
%                 sprintf('Power \\beta = %.2f%%',Beta*100),...
%                 sprintf('Degrees of freedom df = %i',statsHS.df)};
%             text(0.5,0.8,Stats,...
%                 'Units','normalized',...
%                 'FontSize',12,...
%                 'HorizontalAlignment','center')
            
        end
        
        function FibPot = statistical_analysis_surface_potential(obj)
            % Statistical analysis of the surface potential map experiment.
            % First get the grouping of data from the user and then perform
            % several tests, plots etc.
            
            warning('The following methods were programmed for specific use cases and are yet to be generalized! However, you can of course adjust them for your own use, though its probably easier to just take the processed raw data from your Experiment and do your own statistics')
            
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
                boxplot([FibPot(obj.GroupFM(TestMat(i,1)).Indices)' FibPot(obj.GroupFM(TestMat(i,2)).Indices)']*1e3)
                title('Paired T-Test for Surface Potential Changes')
                xticklabels({obj.GroupFM(TestMat(i,1)).Name,obj.GroupFM(TestMat(i,2)).Name})
                xlabel('Test Group')
                ylabel('Surface Potential [mV]')
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
            boxplot([DiffControl DiffMGO])
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
            
            warning('The following methods were programmed for specific use cases and are yet to be generalized! However, you can of course adjust them for your own use, though its probably easier to just take the processed raw data from your Experiment and do your own statistics')
            
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
            boxplot([DiffControl DiffMGO])
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
            
            warning('The following methods were programmed for specific use cases and are yet to be generalized! However, you can of course adjust them for your own use, though its probably easier to just take the processed raw data from your Experiment and do your own statistics')
            
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
            boxplot([DiffControl DiffMGO])
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
            Fig = figure('Units', 'Normalized', 'Position',[0.2 0.2 0.4 8],'Color','w');
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
        
        function deconvolute_cantilever_tip(obj)
            %deconvolute_tip.m
            %Deconvolutes an image of TGT1 grating to obtain the AFM tip using an
            %envelope algorithm. Goes through all the ibw files in the folder and
            %calculates the tip shape. The tip shape is then saved in a new ibw file
            %called "probe_tip_NAME OF THE SCAN.ibw"
            %The ibw files used to find the tip is required to have a ZSensor record as
            %channel 4
            %! ! !ONLY THE ZSENSOR IMAGE IN THE NEW FILE REPRESENTS THE TIP ! ! !
            
            %requires readibw.m, getinfo.m, scalescan.m, cone.m, getpeak.m,
            %minimiseW.m, writedata.m
            
            system_chosen='W';
            
            %Step 1: AFM tip reconstruction
            
            % ---------------------------- FOR SINGLE FILE ----------------------------
            [filename,pathname]=uigetfile('*.jpk','Filedirectory of AFM tip image');
            [~,name,~] = fileparts(filename);
            varname=cell(1,1);
            ImDir=fullfile(pathname,filename);
            % for i = 1:length(flist); ********  FOR MULTIPLE FILES ********
            i=1;
            info=imfinfo(ImDir);
            k=1;
            for j = 2:size(info,1)
                istring(k,1)=find([info(j).UnknownTags.ID]==32851);
                s_name=info(j).UnknownTags(istring(k,1)).Value;
                % ONLY THE MEASURED HEIGHT RETRACE IMAGES ARE STORED FOR ANALYSIS
                % OLD version
                % s1: is a character of the tagged image which specifies the image channel
                % and the image direction i.e. trace or retrace.
                s1={['channel.name : capacitiveSensorHeight' char(10) 'channel.type : channel' char(10) 'retrace : true' char(10) 'type : channel-retrace' char(10) '']};
                % NEW version
                % By NEW version we mean the JPK Software Version 5.0.72 which was updated
                % after the hard drive replacement. With this version the tags are changed.
                % The measurd height and image is located in a different position of the
                % UnkwownTag.
                % s2: is a character as is in the newer version that speficies the image
                % channel and the imaging direction, i.e. trace or retrace.
                s2={['channel.name : measuredHeight' char(10) 'channel.type : channel' char(10) 'retrace : true' char(10) 'type : channel-retrace' char(10) '']};
                % You will probably need this line for a further updated version!
                
                if isequal(s_name,s1{1,1})==1
                    pos(i,1)=j;
                    
                elseif isequal(s_name,s2{1,1})==1
                    pos(i,1)=j;
                    
                end
                k=k+1;
            end
            % end ********  FOR MULTIPLE FILES ********
            
            %
            % i=3
            % for i = 1:length(flist);********  FOR MULTIPLE FILES ********
            
            data=double(imread(ImDir,pos(i,1)));  %read height data from each file
            
            
            
            
            %finding the scansize and the image pixelsize
            %scansize
            iscansizex=find([info(1).UnknownTags.ID]==32834); % location of the scan size in um
            iscansizey=find([info(1).UnknownTags.ID]==32835); % location of the scan size in um
            Scansizex=info(1).UnknownTags(iscansizex).Value*10^6; %in um
            Scansizey=info(1).UnknownTags(iscansizey).Value*10^6; %in um
            %pixelsize
            %        pixels_x = info(2).Width;
            %        pixels_yt = info(2).Height;
            
            %finding the multiplier and the offset value
            imult=find([info(pos(i,1)).UnknownTags.ID]==33028);
            ioffset=find([info(pos(i,1)).UnknownTags.ID]==33029);
            
            
            mult=info(pos(i,1)).UnknownTags(imult).Value; %calls multiplier of the jpk file
            offset=info(pos(i,1)).UnknownTags(ioffset).Value;   %calls offset of the jpk file
            himage = (offset + data.*mult); %Calculate real height data by applying offset and mult to the var data
            
            
            
            s = strrep(ImDir,'.jpk','');%define variable 'filename'
            name_list(i,1)=cellstr(s);
            varname{i,1}=himage;
            % end;********  FOR MULTIPLE FILES ********
            clear flist i himage s  imult ioffset iscansizex iscansizey ...
                data offset mult info s1 s2 pos ans istring j k s_name ext
            %clears temporary variables
            
            
            varname{1,1} = obj.planefit_tgt1(varname{1,1});
            
            
            for i = 1:length(name_list) %loop through all the ibw files
                
                %     if size(varname,3)<4
                %         disp(['Error! Channel 4 is missing in ' name_list{i} '.ibw']);
                %         disp('Please ensure that the Z Sensor data is recorded on channel 4');
                %         disp('Proceeding to next file');
                %     else
                [z, height]=obj.scalescan(varname{i,1}); %put the ZSensor scan image in the positive range
                pixelz_x=size(z,1);
                pixelz_y=size(z,2);
                [s,size_pixel_x,size_pixel_y] = obj.cone(pixelz_x, pixelz_y, height,filename); %generate the perfect cone of TGT1
                [max_x, max_y] = obj.getpeak(s); %find location of the tip peak
                pixels_x=size(s,1);
                pixels_y=size(s,2);
                
                v = genvarname(['probe_tip']);
                g = strcat(name_list{i},'.jpk');
                eval([v ' = obj.minimiseW(z,s,pixelz_x,pixelz_y,pixels_x,pixels_y,max_x, max_y);']); %calculate the tip geometry
            end
            
            
            [depth_num(:,1), projected_area(:,1)] = obj.proj_area(...
                probe_tip,Scansizex*10^-6);
            % depth_num in NANOMETERS
            % projected_area in METERS
            % Store the tip data in a structure variable
            Tip.data = probe_tip;
            Tip.ProjArea = projected_area; % in METERS
            Tip.HeightData = depth_num; % in NANOMETERS
            Tip.XaxisSizeUM = Scansizex;
            Tip.YaxisSizeUM = Scansizey;
            
            % write relevant tip deconv data into experiment property
            % obj.CantileverTip
            
            obj.CantileverTip = Tip;
            obj.CantileverTipFlag = 1;
            
            tipdataname = sprintf('%s\%s.mat',obj.ExperimentFolder,name);
            save(tipdataname)
        end
        
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
        
        function check_for_new_host(obj)
            % check_for_new_host(obj)
            %
            % Checks, if the system environment has changed and, if so,
            % adjusts object properties
            
            FullOS = computer;
            OS = FullOS(1:3);
            if isequal(OS,'PCW')
                Host = getenv('COMPUTERNAME');
            elseif isequal(OS,'GLN')
                Host = getenv('HOSTNAME');
            elseif isequal(OS,'MAC')
                Host = getenv('HOSTNAME');
            end
            
            obj.HostOS = OS;
            obj.HostName = Host;
        end
        
    end
    methods(Static)
        % Static auxilary methods mainly for tip deconvolution (code by Orestis Andriotis)
        % and Experiment() loading
        
        function [depth_num, projected_area] = proj_area(shape,image_x_axis)
            % shape=probe_tip_tgt1_run12_0002;
            % image_x_axis = 1.5e-6
            % close all
            % figure;
            % meshc(shape)
            
            % % convert to gray scale
            I_gray=mat2gray(shape);
            %
            % % mask
            % % At what threshold of the maximum height do you want to mask the data
            
            y=0.1;
            I_mask = (I_gray > y);
            boundary=bwboundaries((I_gray > y));
            
            
            
            figure;
            subplot(1,2,2)
            imshow(I_mask)
            title('Binary mask')
            subplot(1,2,1)
            imshow(I_gray)
            title('Image & Masked area')
            
            hold on
            p=patch(boundary{1}(:,2),boundary{1}(:,1),'g','EdgeColor', [0.8 1 0.4]);
            set(p,'FaceAlpha',0.2)
            hold off
            
            
            % steps
            
            hmax = max(max(shape)).*10^9; % NANOMETERS
            hmin = min(min(shape)).*10^9; % NANOMETERS
            
            h_abs = hmax-hmin; % absolute height NANOMETERS
            
            prompt = {...
                'Step size (in nanometer):'};
            dlg_title = 'Input data';
            num_lines = 1;
            def = {'1'};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            step = str2double(cell2mat(answer(1,1)));
            % Stores the reference slope, the radius
            clear prompt dlg_title num_lines def answer;
            % step = 1; % 1 nanometer step
            % I_tip_bin=(I_gray > y);
            % I_tip = I_gray.*I_tip_bin;
            
            
            
            
            
            pixel_step = step./h_abs;
            clear depth_points
            
            depth_points(:,1) = 1:-pixel_step:y;
            
            % estimate area with BWAREA matlab function.
            % clear area_bin
            k=1;
            for i=1:-pixel_step:y;
                
                area_bin(k,1) = bwarea((I_gray > i));
                
                k=k+1;
            end
            
            
            size_pixel_x=image_x_axis/length(shape); % in METERS
            area_per_pixel=size_pixel_x.^2; % in SQUARE METERS
            
            depth_num(:,1)= 0:step:(1-y)*h_abs; % in NANOMETERS
            
            projected_area=area_per_pixel.*area_bin; % in SQUARE METERS
            
        end
        
        function wdata = getinfo(filename)
            %getinfo.m
            %reads header info and extracts wave data from ibw file
            fid = fopen(filename);
            fseek(fid,80,'bof');
            type = fread(fid,1,'int16');%read type
            fseek(fid,50,'cof');
            x = fread(fid,1,'int32');%read width
            y = fread(fid,1,'int32');%read height
            n = fread(fid,1,'int32');%read no of images
            if type == 2
                d = 'float32';%32 bit float
            elseif type == 3
                d = 'float64';%64 bit float
            elseif type == 4
                d = 'int8';%8 bit signed integer
            else disp('type unknown');
            end
            fseek(fid,384,'bof');%start of wave data
            wdata = zeros(x,y,n);%preallocate wdata array for speed
            for i = 1:n
                wdata(:,:,i) = fread(fid,[x,y],d);%read wave data into 3d array
            end
            fclose(fid);
        end
        
        function [leveled_scan,height]= scalescan(scan)
            %scalescan.m
            %moves the image in the positive scale so that its minimum height is 0
            
            leveled_scan=scan-min(scan(:)); %min(scan(:)) is faster than min(min(scan))
            height=max(leveled_scan(:)); %idem
        end
        
        function [shape,size_pixel_x,size_pixel_y] = cone(pixel_x,pixel_y,height,name_scan)
            %cone.m, version 1.1
            %The file creates a TGT1 grating surface with a cone peaking at the centre
            %of the sample. The angle of the is assumed to be 50 degrees and the height
            %is taken to be equal to the scan height. The radius of the cone tip can be
            %changed so that different levels of accuracy may be achieved.
            %Ask the user for the size of the scan.
            prompt = {['What is the scanning range in micrometres for  ' name_scan '.ibw ?']};
            dlg_title = 'Scan Size';
            num_lines= 1;
            def     = {'1.5'};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            size_scan= (str2num(cell2mat(answer(1,1))))*1e-6;
            
            %Asks the user for the radius of the grating peak.
            prompt = {'What is the radius of the grating tip in nanometres?'};
            dlg_title = 'Grating Tip Radius';
            num_lines= 1;
            def     = {'5'};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            peak_radius= (str2num(cell2mat(answer(1,1))))*1e-9;
            
            %Variables
            height_cone=height;
            angle_cone=50;
            height_loss=(peak_radius*cosd(angle_cone/2))/(tand(angle_cone/2))+(peak_radius*sind(angle_cone/2))-peak_radius;
            %Calculates the height that is lost from a perfect tip when a rounded tip
            %is used instead.
            height_tip=height_cone+height_loss;%The ideal tip is derived from the required
            %height of the cone (from the scan height) and the amount of height loss
            %experienced for the desired grating radius. This value is then used to
            %generate the ideal cone with a perfect tip such that when the curved tip
            %is added the height of the cone is equal to the scan height.
            
            radius_cone=tand(angle_cone/2)*height_tip; %Calculates the cone radius.
            shape=zeros(pixel_x,pixel_y); %initiates a flat surface of size equal to scan.
            size_pixel_x=size_scan/pixel_x;
            size_pixel_y=size_scan/pixel_y;
            shape(floor(pixel_x/2), floor(pixel_y/2))=height_tip; % Positions the cone.
            max_pixel_movement_x=floor(radius_cone/size_pixel_x); %Radius of cone divided
            %by size of a pixel to find the maximal number of pixels in line in the
            %cone radius.
            max_pixel_movement_y=floor(radius_cone/size_pixel_y);
            
            %Determine the limits of the cone whether it fits completely in the image
            %or not. Done for each dimension and limit using the centre of the image as
            %a reference.
            if pixel_x/2-max_pixel_movement_x>=1
                limit_x_1=pixel_x/2-max_pixel_movement_x;
            else limit_x_1=1;
            end
            
            if pixel_x/2+max_pixel_movement_x<=pixel_x
                limit_x_2=pixel_x/2+max_pixel_movement_x;
            else limit_x_2=pixel_x;
            end
            
            if pixel_y/2-max_pixel_movement_y>=1
                limit_y_1=pixel_y/2-max_pixel_movement_y;
            else limit_y_1=1;
            end
            
            if pixel_y/2+max_pixel_movement_y<=pixel_y
                limit_y_2=pixel_y/2+max_pixel_movement_y;
            else limit_y_2=pixel_y;
            end
            
            
            curve_start_height = height_tip - ((peak_radius*cosd(angle_cone/2))/(tand(angle_cone/2)));
            %Calculates the hieght at which the cone leaves from a constant gradient
            %into the curved profile.
            
            %Generates the cone
            for i=limit_x_1:limit_x_2
                for j=limit_y_1:limit_y_2
                    distance=sqrt(((pixel_x/2-i)*size_pixel_x)^2+((pixel_y/2-j)*size_pixel_y)^2);
                    %Distance of point i,j with reference from the centre of the image.
                    
                    curve_height=sqrt(((peak_radius)^2)-(distance^2))-(peak_radius*sind(angle_cone/2));
                    %The absolute hieght of each point of the curved tip.
                    
                    if distance<=radius_cone;
                        shape(i,j)=(radius_cone-distance)*height_tip/radius_cone;
                        %If the dstance is smaller than than the radius then the
                        %constant slope of the cone is generated.
                    end
                    
                    if distance<=peak_radius*cosd(angle_cone/2);
                        shape(i,j) = curve_height + curve_start_height;
                        %If the distance is smaller than the radius of the curved peak
                        %radius then the curved peak is assumed.
                    end
                    
                end
            end
        end
        
        function [max_x, max_y] = getpeak(scan)
            %getpeak.m
            %Returns the coordinates of the peak of the image
            
            [c,i]=max(scan);
            [c,max_y]=max(c);
            max_x=i(max_y);
        end
        
        function tip=minimiseW(z,s,pixelz_x,pixelz_y,pixels_x,pixels_y,max_x,max_y)
            %minimiseW.m
            %Minimises the w function w(x,y,x',y')=z(x',y')-s(x-x',y-y') to erode the
            %image. Returns the real image of the sample or the tip shape depending on
            %the input s. To obtain the tip shape, s must be the sample. To
            %obtain the eroded image, s must be the tip upside down.
            
            tip=ones(pixelz_x,pixelz_y); %creates the empty image array
            
            h = waitbar(0,'Please wait, Processing...');
            
            %loops over all the elements and find the minimum value of w and allocate it
            for j=1:pixelz_y %loops over points in image output
                waitbar(j/pixelz_y);
                s_ymin=max(-j,-max_y); %determines the allowed range for tip scanning
                s_ymax=min(pixels_y-max_y,pixelz_y-j)-1; %idem
                for i=1:pixelz_x %loops over points in other direction in image output
                    s_xmin=max(-max_x,-i); %determines allowed range for tip scanning
                    s_xmax=min(pixels_x-max_x,pixelz_x-i)-1; %idem
                    %need to add 1 for matrix starts on row 1 and not 0
                    minimum=z(i+s_xmin+1,j+s_ymin+1)-s(s_xmin+max_x+1,s_ymin+max_y+1);
                    for k=s_xmin:s_xmax %loop over points in tip
                        for l=s_ymin:s_ymax %idem
                            %need to add 1 for matrix starts on row 1 and not 0
                            temp=z(i+k+1,j+l+1)-s(k+max_x+1,l+max_y+1); %calculate w.
                            minimum=min(temp,minimum); %checks if w is minimum
                        end
                    end
                    tip(i,j)=minimum; %allocates the minimum value
                end
            end
            close(h);
        end
        
        function image_pfit = planefit_tgt1(image)
            
            
            %
            % image=varname{1,1};
            leng = length(image);
            width = length(image);
            
            Igray = mat2gray(image); %convert to gray scale
            
            % select ROI to mask
            
            figure;
            subplot(1,2,1)
            hold on
            title('median filter')
            % imagesc(himage_filt2)
            imagesc(Igray)
            axis([0 width 0 leng])
            axis square
            % colorbar
            caxis([0 1])
            set(gca,'FontSize',8)
            %create mask of interest
            maskROI=roipoly;
            % Invert mask
            maskROI = (maskROI ==0);
            % plot the mask region of interest
            
            subplot(1,2,2)
            hold on
            
            title('mask')
            imagesc(flip(maskROI,1))
            axis([0 width 0 leng])
            axis square
            caxis([0 1])
            set(gca,'FontSize',8)
            
            
            
            backgroundI=image.*maskROI;
            
            
            %
            %prealocate matrices and vectors for planefit
            %matrix for x-coordinate vector
            A=zeros(leng,width);
            %matrix for y-coordinate vector
            B=zeros(leng,width);
            %matrix including height data of background
            C=backgroundI;
            %matrices for creating vectors including x,y coordinates of original image
            %size
            %matrix including x positions
            D=zeros(leng,width);
            %matrix including y positions
            E=zeros(leng,width);
            %matrix including height values of the fitted plane
            F=zeros(leng,width);
            %filling matrix with the specific x and y values
            for i=1:leng
                A(i,:)=1:width;
                B(i,:)=i*ones(1,width);
                D(i,:)=1:width;
                E(i,:)=1:width;
                F(i,:)=i*ones(1,width);
            end
            
            %apply mask on all the matrices
            A=A.*maskROI;
            B=B.*maskROI;
            
            %creating vectors out of the matrices
            xx=reshape(A',[numel(A),1]);
            yy=reshape(B',[numel(B),1]);
            zz=reshape(C',[numel(C),1]);
            xo=reshape(E',[numel(E),1]);
            yo=reshape(F',[numel(F),1]);
            zo=zeros(leng*width,1);
            
            %deleting all positions where height was set zero from masking
            X=[xx,yy,zz];
            X(~any(X,2),:)=[];  %rows where x,y,z are empty should be taken out
            %rebuild xx,yy,zz vectors of Matrix X
            xx=X(:,1);
            yy=X(:,2);
            zz=X(:,3);
            
            %plnefitting by fit a linear function to each data line of background image
            linefit=D; %this matrix should constist of flattened lines
            for l=1:leng
                ins=find(yy==l,1,'first');
                inl=find(yy==l,1,'last');
                P=polyfit(xx(ins:inl,1),zz(ins:inl,1),1);
                linefit(l,:)=polyval(P,D(l,:));
            end
            image_pfit=image-linefit;
            
            %
            figure;
            
            subplot(1,2,1)
            hold on
            grid on
            % plot original image
            mesh(image,'facealpha',0.6);
            axis tight
            axis square
            view([-32,20])
            % plot original image
            mesh(linefit,'FaceColor',[0.5 0.5 0.3],'FaceAlpha',0.3,'EdgeColor','none');
            colormap jet
            set(gca,'FontSize',8)
            zlabel('Height (m)')
            
            subplot(1,2,2)
            
            hold on
            grid on
            % plot planefited image
            mesh(image_pfit,'facealpha',0.8);
            axis tight
            axis square
            colormap jet
            view([-32,20])
            set(gca,'FontSize',8)
            zlabel('Height (m)')
            
%             waitforbuttonpress
            close all
            
            
            clear maskROI Igray backgroundI leng width
        end
        
    end
end