classdef Experiment < matlab.mixin.Copyable & matlab.mixin.SetGet
    
    properties
        % Essential properties for File and subclass management
        ExperimentName      % Shows name of the experiment
        ExperimentFolder    % Shows Folder where Experiment is saved
        HostOS              % Shows the current operating system
        HostName            % Shows the current Host (User of running machine)
        BigDataFlag         % Determines how Experiment loads ForceMap objects. If true, 
                            % force voluime data are not loaded into RAM
                            % but read from unpacked jpk data container
        CurrentLogFile
        FM                  % Cellarray containing Force/QI Maps
        NumForceMaps
        ForceMapNames       
        ForceMapFolders
        RefFM               % Cellarray containing reference Force/QI Maps
        NumReferenceForceMaps
        ReferenceForceMapNames
        ReferenceForceMapFolders
        I                   % Cellarray containing AFM-Images
        NumAFMImages
        AFMImageNames
        AFMImageFolders
        SPM                 % Cellarray containging Surface Potential Maps (to be made obsolete by AFMImage in the future)
        NumSurfacePotentialMaps
        SurfacePotentialMapFolders
        SurfacePotentialMapNames
        CantileverTips      % Cellarray containing images of Cantilever Tips from TGT-1 grating for deconvolution
        NumCantileverTips
        CantileverTipNames
        CantileverTipFolders
    end
    properties
        % Non-essential properties
        EMod
        WhichRefMap
        WhichTip
        SurfPot
        GroupFM
        GroupSPM
        DropoutNet
        CP_CNN
        idxSubstrate
        idxEnvCond
        MinFM
        NumFcUncorrupt % SMFS property
    end
    properties
        % All the Flags
        FMFlag
        SPMFlag
        SMFSFlag
        SMFSFlagDown
        DebugFlag
        ReferenceSlopeFlag
        AssignedReferenceMaps
        CantileverTipFlag
        AssignedCantileverTips
    end
    
    % SM related
    properties
        SMFSFMParameters
        SMFSResults
        SMFSResultsParameters
        WLCFit
        SMFSLillieAdhMaxApp
        SMFSLillieAdhMaxRet
        SMFSLillieAdhUnbinding
        SMFSLillieAdhEneApp
        SMFSLillieAdhEneRet
        SMFSLilliePullingLength
        SMFSLillieSnapInLength
        SMFSWilcoxonAdhMaxApp
        SMFSWilcoxonAdhMaxRet
        SMFSWilcoxonAdhUnbinding
        SMFSWilcoxonAdhEneApp
        SMFSWilcoxonAdhEneRet
        SMFSWilcoxonPullingLength
        SMFSWilcoxonSnapInLength
    end
    
    methods
        % constructor method and methods related with Experiment-file handling
        
        function obj = Experiment()
            % Author: Manuel Rufin
            % Description:
            % This is the constructor method of the class. The function creates an Experiment-object by calling different functions that leads to the construction of the object. 
            % First, host settings are determined ('obj.check_for_new_host'). 
            % Second, file name, file type - .mat-file - and storage location of the newly created object are defined. 
            % Third, a GUI window allows to assign the raw data according to the jpk-data-type and the corresponding data paths are determined ('obj.constructor_user_input_parser'). 
            % Fourth, the raw data is loaded and stored in a defined object structure ('obj.take_paths_and_load_files'). 
            % Fifth, flags are initialized (obj.initialize_flags) which are used in the further course of the analysis. 
            % Sixth, the newly created object is saved in the defined storage location (obj.save_experiment).

            %% Function body
            % Set HostOS and HostName properties
            obj.check_for_new_host
            
            % set Name and choose Experiment folder
            [ExperimentName, ParentFolder] = uiputfile('*.mat','Select Experiment Name and Parent Folder');
            obj.ExperimentName = ExperimentName(1:end-4);
            mkdir(ParentFolder,obj.ExperimentName);
            obj.ExperimentFolder = fullfile(ParentFolder,obj.ExperimentName,filesep);
            
            % get Experiment filenames and paths from user input
            [FileTypes, FullFileStruct, IsValid, BigData] = obj.constructor_user_input_parser(obj.ExperimentName,obj.HostOS);
            
            if ~IsValid
                obj = [];
                return
            end
            
            % get paths of requested files and load them in
            obj.take_paths_and_load_files(FileTypes,FullFileStruct,true,BigData)
            
            obj.initialize_flags
            
            Temp = load('DropoutNetFinal.mat');
            obj.DropoutNet = Temp.MC14_Drop;
            Temp2 = load('CP_CNN_Final.mat');
            obj.CP_CNN = Temp2.CNN;
            
            obj.save_experiment();
        end
        
        function take_paths_and_load_files(obj,FileTypes,FullFileStruct,isNew,BigData)
            % Author: Manuel Rufin
            % Description: This function loads the raw data and stores it
            % in a defined object structure.

            %% Function body
            if nargin<4
                isNew = false;
            end
            if isNew
                % set everything to zero if creating new experiment.
                % otherwise skip and just add new entries
                obj.FM = cell(0,0);
                obj.NumForceMaps = 0;
                obj.ForceMapNames = cell(0,0);
                obj.ForceMapFolders = cell(0,0);
                obj.RefFM = cell(0,0);
                obj.NumReferenceForceMaps = 0;
                obj.ReferenceForceMapNames = cell(0,0);
                obj.ReferenceForceMapFolders = cell(0,0);
                obj.I = cell(0,0);
                obj.NumAFMImages = 0;
                obj.AFMImageNames = cell(0,0);
                obj.AFMImageFolders = cell(0,0);
                obj.SPM = cell(0,0);
                obj.NumSurfacePotentialMaps = 0;
                obj.SurfacePotentialMapFolders = cell(0,0);
                obj.SurfacePotentialMapNames = cell(0,0);
                obj.CantileverTips = cell(0,0);
                obj.NumCantileverTips = 0;
                obj.CantileverTipNames = cell(0,0);
                obj.CantileverTipFolders = cell(0,0);
                obj.BigDataFlag = BigData;
            end
            
            % Need to assign something to *FullFile. Otherwise parfor will
            % crash
            for i=1:5
                if isempty(FullFileStruct(i).FullFile{1})
                    NumFiles(i) = 0;
                else
                    NumFiles(i) = length(FullFileStruct(i).FullFile);
                end
            end
            
            FMFullFile = cell(max(NumFiles),1);
            RefFMFullFile = cell(max(NumFiles),1);
            IFullFile = cell(max(NumFiles),1);
            SPMFullFile = cell(max(NumFiles),1);
            CantTipFullFile = cell(max(NumFiles),1);
            
            FMFullFile(1:NumFiles(1)) = FullFileStruct(1).FullFile;
            RefFMFullFile(1:NumFiles(2)) = FullFileStruct(2).FullFile;
            IFullFile(1:NumFiles(3)) = FullFileStruct(3).FullFile;
            SPMFullFile(1:NumFiles(4)) = FullFileStruct(4).FullFile;
            CantTipFullFile(1:NumFiles(5)) = FullFileStruct(5).FullFile;
            
            
            StartID = obj.NumAFMImages + obj.NumCantileverTips + obj.NumForceMaps + obj.NumReferenceForceMaps + obj.NumSurfacePotentialMaps + 1;
            IDs = StartID:(StartID + sum(NumFiles)-1);
            TempCell = cell(max(NumFiles),5);
            TempID = cell(max(NumFiles),5);
            L = max(NumFiles);
            
            if contains(struct2array(ver), 'Parallel Computing Toolbox') && ((sum(NumFiles(1:2)) > 1) || (sum(NumFiles) > 20))
                for i=1:5
                    parfor j=1:L
                        if (j+sum(NumFiles(1:i))-NumFiles(i)) <= length(IDs)
                            TempID{j,i} = sprintf('%s-%i',obj.ExperimentName,IDs(j+sum(NumFiles(1:i))-NumFiles(i)));
                        end
                        if i == 1 && FileTypes(i) && (j<=NumFiles(i))
                            TempCell{j,i} = ForceMap(FMFullFile{j},obj.ExperimentFolder,TempID{j,i},obj.BigDataFlag);
                        end
                        if i == 2 && FileTypes(i) && (j<=NumFiles(i))
                            TempCell{j,i} = ForceMap(RefFMFullFile{j},obj.ExperimentFolder,TempID{j,i},obj.BigDataFlag);
                        end
                        if i == 3 && FileTypes(i) && (j<=NumFiles(i))
                            TempCell{j,i} = AFMImage(IFullFile{j},obj.ExperimentFolder,TempID{j,i});
                        end
                        if i == 4 && FileTypes(i) && (j<=NumFiles(i))
                            TempCell{j,i} = SurfacePotentialMap(SPMFullFile{j},TempID{j,i});
                        end
                        if i == 5 && FileTypes(i) && (j<=NumFiles(i))
                            TempCell{j,i} = AFMImage(CantTipFullFile{j},obj.ExperimentFolder,TempID{j,i});
                        end
                    end
                end
            else
                for i=1:5
                    for j=1:NumFiles(i)
                        TempID{j,i} = sprintf('%s-%i',obj.ExperimentName,IDs(j+sum(NumFiles(1:i))-NumFiles(i)));
                        if i == 1 && FileTypes(i)
                            TempCell{j,i} = ForceMap(FMFullFile{j},obj.ExperimentFolder,TempID{j,i},obj.BigDataFlag);
                        end
                        if i == 2 && FileTypes(i)
                            TempCell{j,i} = ForceMap(RefFMFullFile{j},obj.ExperimentFolder,TempID{j,i},obj.BigDataFlag);
                        end
                        if i == 3 && FileTypes(i)
                            TempCell{j,i} = AFMImage(IFullFile{j},obj.ExperimentFolder,TempID{j,i});
                        end
                        if i == 4 && FileTypes(i)
                            TempCell{j,i} = SurfacePotentialMap(SPMFullFile{j},TempID{j,i});
                        end
                        if i == 5 && FileTypes(i)
                            TempCell{j,i} = AFMImage(CantTipFullFile{j},obj.ExperimentFolder,TempID{j,i});
                        end
                    end
                end
            end
            
            % Now write everything into the obj
            for i=1:5
                for j=1:NumFiles(i)
                    if i == 1 && FileTypes(i)
                        Idx = obj.NumForceMaps+j;
                        obj.FM{Idx,1} = TempCell{j,i};
                        obj.ForceMapNames{Idx,1} = obj.FM{Idx}.Name;
                        obj.ForceMapFolders{Idx,1} = obj.FM{Idx}.Folder;
                    end
                    if i == 2 && FileTypes(i)
                        Idx = obj.NumReferenceForceMaps+j;
                        obj.RefFM{Idx,1} = TempCell{j,i};
                        obj.ReferenceForceMapNames{Idx,1} = obj.RefFM{Idx}.Name;
                        obj.ReferenceForceMapFolders{Idx,1} = obj.RefFM{Idx}.Folder;
                    end
                    if i == 3 && FileTypes(i)
                        Idx = obj.NumAFMImages+j;
                        obj.I{Idx,1} = TempCell{j,i};
                        obj.AFMImageNames{Idx,1} = obj.I{Idx}.Name;
                        obj.AFMImageFolders{Idx,1} = 'Placeholder';
                    end
                    if i == 4 && FileTypes(i)
                        Idx = obj.NumSurfacePotentialMaps+j;
                        obj.SPM{Idx,1} = TempCell{j,i};
                        obj.SurfacePotentialMapNames{Idx,1} = obj.SPM{Idx}.Name;
                        obj.SurfacePotentialMapFolders{Idx,1} = obj.SPM{Idx}.Folder;
                    end
                    if i == 5 && FileTypes(i)
                        Idx = obj.NumCantileverTips+j;
                        obj.CantileverTips{Idx,1} = TempCell{j,i};
                        obj.CantileverTipNames{Idx,1} = obj.CantileverTips{Idx}.Name;
                        obj.CantileverTipFolders{Idx,1} = 'Placeholder';
                    end
                end
            end
            
            obj.NumForceMaps = obj.NumForceMaps + NumFiles(1);
            obj.NumReferenceForceMaps = obj.NumReferenceForceMaps + NumFiles(2);
            obj.NumAFMImages = obj.NumAFMImages + NumFiles(3);
            obj.NumSurfacePotentialMaps = obj.NumSurfacePotentialMaps + NumFiles(4);
            obj.NumCantileverTips = obj.NumCantileverTips + NumFiles(5);
            
        end
        
        function Out = add_data(obj)
            
            warning('Note that for this function to work properly you need to assign your Experiment in workspace to itself e.g. ">> E = E.add_data" otherwise the updated Experiment will be stored in the temporary variable "ans" and inevitably overwritten at some point!')
            
            % create save copy to restore if function produces errors
            SaveCopy = obj.copy_experiment;
            
            try
                
                % Set HostOS and HostName properties
                SaveCopy.check_for_new_host
                
                % get Experiment name and layout from user
                isNew = false;
                [FileTypes, FullFileStruct, isValid] = SaveCopy.constructor_user_input_parser(obj.ExperimentName,obj.HostOS);
                
                if ~isValid
                    SaveCopy = [];
                    return
                end
                
                % get paths of requested files and load them in
                SaveCopy.take_paths_and_load_files(FileTypes,FullFileStruct,isNew)
                
                % SaveCopy.initialize_flags % What to do with this?
                
                Out = SaveCopy;
                Out.save_experiment
                warning('Did you read the warning above?');
            catch ME
                warning(ME.getReport)
                disp('data adding failed. restored original experiment object')
                fclose('all');
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
        
        function save_experiment(obj)
            % Author: Manuel Rufin
            % Description: 
            % This funciton saves the newly created object in the defined storage location 

            %% Function body
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
            for i=1:max([obj.NumAFMImages obj.NumCantileverTips obj.NumForceMaps obj.NumReferenceForceMaps obj.NumSurfacePotentialMaps])
                if i<=length(obj.FM)
                    MCFM = metaclass(obj.FM{i});
                end
                if i<=length(obj.RefFM)
                    MCRefFM = metaclass(obj.RefFM{i});
                end
                if i<=length(obj.I)
                    MCI = metaclass(obj.I{i});
                end
                if i<=length(obj.SPM)
                    MCSPM = metaclass(obj.SPM{i});
                end
                if i<=length(obj.CantileverTips)
                    MCCantileverTips = metaclass(obj.CantileverTips{i});
                end
                if i<=length(obj.FM) && ~isempty(MCFM.SuperclassList) && isequal(MCFM.SuperclassList.Name,'matlab.mixin.Copyable')
                    ExperimentCopy.FM{i} = obj.FM{i}.copy;
                end
                if i<=length(obj.RefFM) && ~isempty(MCRefFM.SuperclassList) && isequal(MCRefFM.SuperclassList.Name,'matlab.mixin.Copyable')
                    ExperimentCopy.RefFM{i} = obj.RefFM{i}.copy;
                end
                if i<=length(obj.I) && ~isempty(MCI.SuperclassList) && isequal(MCI.SuperclassList.Name,'matlab.mixin.Copyable')
                    ExperimentCopy.I{i} = obj.I{i}.copy;
                end
                if i<=length(obj.SPM) && ~isempty(MCSPM.SuperclassList) && isequal(MCSPM.SuperclassList.Name,'matlab.mixin.Copyable')
                    ExperimentCopy.SPM{i} = obj.SPM{i}.copy;
                end
                if i<=length(obj.CantileverTips) && ~isempty(MCCantileverTips.SuperclassList) && isequal(MCCantileverTips.SuperclassList.Name,'matlab.mixin.Copyable')
                    ExperimentCopy.CantileverTips{i} = obj.CantileverTips{i}.copy;
                end
            end
        end
        
        function delete_experiment(obj)
            
            if ~isequal(obj.HostOS,'PCW')
                disp('This method is only available on Windows systems at the moment');
                return
            end
            
            Folder = string(obj.ExperimentFolder);
            
            answer = questdlg('Caution! This PERMANENTLY deletes the experiment folder, experiment file and all subfolders', ...
                sprintf('Deletion of %s',Folder),'Yes, delete all', ...
                'Abort','Abort');
            % Handle response
            switch answer
                case 'Yes, delete all'
                    disp(sprintf('Deleting %s. This may take from a few to tens of minutes',obj.ExperimentName))
                    
                    current = what();
                    cd(Folder)
                    cd ..
                    
                    cmd1_1 = string('DEL /F/Q/S ');
                    cmd1_2 = string(' > nul');
                    cmd2_1 = string('RMDIR /Q/S ');
                    
                    CMD1 = strcat(cmd1_1,Folder,cmd1_2);
                    CMD2 = strcat(cmd2_1,Folder);
                    
                    system(CMD1);
                    system(CMD2);
                    
                    cd(current.path)
                case 'Abort'
                    return
            end
            
        end
    end
    methods(Static)
        % Static methods related with Experiment-file handling
        
        function E = load(Fullfile)
            % E = load()
            %
            % recommended way of loading an existing Experiment() from its
            % folder. Checks, if path has changed and adjusts object
            % properties. Also checks if running on different system,
            % and, if so, updating object properties and setting
            % CPFlag.CNNOpt = 0
            
            if nargin < 1
                [File,Path] = uigetfile('*.mat','Choose Experiment .mat from folder');
                Fullfile = fullfile(Path,File);
                disp('Loading Experiment... this can take a while for larger Experiments')
            end
            load(Fullfile);
            
            E = obj;
            clear obj
            
            E.check_for_new_host();
            FMFolder = fullfile(Path,filesep,'ForceData');
            for i=1:E.NumForceMaps
                if ~isempty(E.FM{i})
                    E.FM{i}.check_for_new_host();
                    OldDataStore = E.FM{i}.DataStoreFolder;
                    Split = strsplit(OldDataStore,filesep);
                    E.FM{i}.DataStoreFolder = fullfile(Path,Split{end-1});
                    E.FM{i}.Folder = FMFolder;
                    E.ForceMapFolders{i} = FMFolder;
                end
            end
            for i=1:E.NumAFMImages
                if ~isempty(E.I{i})
                    E.I{i}.check_for_new_host();
                end
            end
            for i=1:E.NumReferenceForceMaps
                if ~isempty(E.RefFM{i})
                    E.RefFM{i}.check_for_new_host();
                end
            end
            for i=1:E.NumCantileverTips
                if ~isempty(E.CantileverTips{i})
                    E.CantileverTips{i}.check_for_new_host();
                end
            end
            
            E.ExperimentFolder = Path;
            
        end
        
        function delete_folderstructure(FolderPath)
            
            
            FullOS = computer;
            OS = FullOS(1:3);
            
            if ~isequal(OS,'PCW')
                disp('This method is only available on Windows systems at the moment');
                return
            end
            
            Folder = FolderPath;
            
            answer = questdlg({'Caution! This PERMANENTLY deletes the given folder, its files and all subfolders:',Folder}, ...
                sprintf('Deletion of %s',Folder),'Yes, delete all', ...
                'Abort','Abort');
            % Handle response
            switch answer
                case 'Yes, delete all'
                    disp(sprintf('Deleting %s. This may take from a few to tens of minutes',Folder))
                    
                    current = what();
                    cd(Folder)
                    cd ..
                    
                    cmd1_1 = string('DEL /F/Q/S ');
                    cmd1_2 = string(' > nul');
                    cmd2_1 = string('RMDIR /Q/S ');
                    
                    CMD1 = strcat(cmd1_1,Folder,cmd1_2);
                    CMD2 = strcat(cmd2_1,Folder);
                    
                    system(CMD1);
                    system(CMD2);
                    
                    cd(current.path)
                case 'Abort'
                    return
            end
            
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
                
                obj.FMFlag.Preprocessed(i) = 1;
            end
            
            obj.save_experiment;
            
            close(h);
        end
        
        function force_map_analysis_fibril(obj,CPOption,EModOption,BaseLineCorrectBool,TemporaryLoadInBool,UseTipInHertzBool)
            % force_map_analysis_fibril(obj,CPOption,EModOption)
            %
            % CPOption = 'Snap-In' ... Preferred Option for data with
            % snap-in effect
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
            %                           method
            % EModOption = 'HertzCorrected' ... E-Modulus calculation through Hertz-Sneddon
            %                               method with RefSlope-corrected
            %                               senstitivity
            % EModOption = 'Oliver' ... E-Modulus calculation through
            % Oliver-Pharr-like method (O. Andriotis 2014)
            % EModOption = 'Both'   ... 'Oliver' and 'Hertz'
            
            if nargin < 4
                BaseLineCorrectBool = false;
                TemporaryLoadInBool = true;
                UseTipInHertzBool = true;
            elseif nargin < 5
                TemporaryLoadInBool = true;
                UseTipInHertzBool = true;
            end
            
            obj.write_to_log_file('Analysis Function','force_map_analysis_fibril()','start')
            obj.write_to_log_file('Contact Point Option',CPOption)
            obj.write_to_log_file('EMod Option',EModOption)
            obj.write_to_log_file('BaseLineCorrectBool',BaseLineCorrectBool)
            
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
            
            % Preprocessing everything, that needs user input
            answer = questdlg('Do you want to skip manual exclusion of problematic areas?',...
                'Manual Exclusion',...
                'Yes',...
                'No','No');
            for i=1:NLoop
                if isequal(KeepFlagged,'Yes') && obj.FMFlag.FibrilAnalysis(i) == 1
                    continue
                end
                obj.FM{i}.create_fibril_mask();
                if isequal(answer,'Yes')
                    continue
                end
                obj.FM{i}.manual_exclusion();
            end
            
            % Setting and calculating preferred method of reference slope
            obj.reference_slope_parser(5) % The input argument sets the default refslope method to AutomaticFibril
            
            if obj.ReferenceSlopeFlag.SetAllToValue
                RefSlopeOption = 'SetAllToValue';
            elseif obj.ReferenceSlopeFlag.UserInput
                RefSlopeOption = 'UserInput';
            elseif obj.ReferenceSlopeFlag.FromRefFM
                RefSlopeOption = 'FromRefFM';
            elseif obj.ReferenceSlopeFlag.FromArea
                RefSlopeOption = 'FromArea';
            elseif obj.ReferenceSlopeFlag.AutomaticFibril
                RefSlopeOption = 'AutomaticFibril';
            elseif obj.ReferenceSlopeFlag.Automatic
                RefSlopeOption = 'Automatic';
            end
            obj.write_to_log_file('Reference Slope Option',RefSlopeOption)
            
            % Deconvoluting cantilever tip(s)
            if isequal(lower(EModOption),'oliver') || isequal(lower(EModOption),'both') || UseTipInHertzBool
                if obj.NumCantileverTips == 0
                    if isequal(class(obj.CantileverTips{1}),'AFMImage')
                        obj.NumCantileverTips = length(obj.CantileverTips);
                        obj.CantileverTipFlag = 1;
                        obj.AssignedCantileverTips = 1;
                        obj.WhichTip = ones(obj.NumForceMaps,1);
                    else
                        Warn = warndlg('You need to load in TGT-1 images of your cantilever for this kind of analysis');
                        uiwait(Warn);
                        IsValid = false;
                        while ~IsValid
                            UsrInput = inputdlg('How many tips were used in this Experiment?');
                            NumTips = str2num(UsrInput{1});
                            IsValid = isnumeric(NumTips)&&~isempty(NumTips);
                        end
                        obj.get_paths_and_load_files([0 0 0 0 1],[0 0 0 0 ceil(NumTips)],false);
                    end
                end
                if ~obj.AssignedCantileverTips
                    obj.assign_cantilever_tips
                end
                % Check if all needed tips are deconvoluted, if not, do it
                for i=1:obj.NumForceMaps
                    if ~obj.CantileverTips{obj.WhichTip(i)}.hasDeconvolutedCantileverTip
                        obj.CantileverTips{obj.WhichTip(i)}.deconvolute_cantilever_tip;
                    end
                end
            end
            
            % Main loop for contact point estimation, Fibril Diameter and
            % Fibril E-Modulus calculation
            for i=1:NLoop
                if isequal(KeepFlagged,'Yes') && obj.FMFlag.FibrilAnalysis(i) == 1
                    continue
                end
                
                waitbar(i/NLoop,h,sprintf('Processing Fibril %i/%i\nReading out data...',i,NLoop));
                if TemporaryLoadInBool && obj.BigDataFlag
                    obj.FM{i}.temporary_data_load_in(true);
                end
                
                waitbar(i/NLoop,h,sprintf('Processing Fibril %i/%i\nCreating and levelling Height Map',i,NLoop));
                obj.FM{i}.create_and_level_height_map
                obj.FM{i}.create_automatic_background_mask(1)
                
                Thresh = 1/2;
                AppRetSwitch = 2;
                obj.FM{i}.unselect_curves_by_fraction_of_max_data_points(Thresh,AppRetSwitch)
                
                waitbar(i/NLoop,h,sprintf('Processing Fibril %i/%i\nFitting Base Line',i,NLoop));
                if ~obj.FM{i}.BaseAndTiltFlag
                    obj.FM{i}.base_and_tilt('linear');
                end
                if i == 1
                    obj.write_to_log_file('Baseline and Tilt option','linear')
                end
                
                obj.FM{i}.calculate_fib_diam();
                
                % contact point estimation happens here
                waitbar(i/NLoop,h,sprintf('Processing Fibril %i/%i\nFinding Contact Point',i,NLoop));
                if BaseLineCorrectBool
                    if ~obj.FM{i}.CPFlag.CNN
                        obj.cp_option_converter('Fast',i);
                    end
                    FractionBeforeCP = .7;
                    obj.FM{i}.base_and_tilt_using_cp(FractionBeforeCP)
                    obj.write_to_log_file('FractionBeforeCP',FractionBeforeCP);
                end
                
                if ~obj.FM{i}.CPFlag.CNNZoomSweep
                obj.cp_option_converter(CPOption,i);
                end
                
                % reference slope calculation happens here
                waitbar(i/NLoop,h,sprintf('Processing Fibril %i/%i\nProcessing and calculating Reference Slope',i,NLoop));
                obj.reference_slope_calculator(i);
                
                waitbar(i/NLoop,h,sprintf('Processing Fibril %i/%i\nCalculating E-Modulus',i,NLoop));
                if isequal(lower(EModOption),'hertz') || isequal(lower(EModOption),'both') || isequal(lower(EModOption),'hertzcorrected')
                    if isequal(lower(EModOption),'hertzcorrected') || isequal(lower(EModOption),'both')
                        CorrectSens = true;
                    else
                        CorrectSens = false;
                    end
                    AllowXShift = true;
                    obj.FM{i}.calculate_e_mod_hertz(CPOption,'parabolic',1,AllowXShift,CorrectSens,UseTipInHertzBool,0,obj.CantileverTips{obj.WhichTip(i)});
                    if i == 1
                        obj.write_to_log_file('Hertzian Tip-Shape','parabolic')
                        obj.write_to_log_file('Hertzian CurvePercent','1')
                        obj.write_to_log_file('Allow X-Shift',AllowXShift)
                    end
                end
                if isequal(lower(EModOption),'oliver') || isequal(lower(EModOption),'both')
                    obj.FM{i}.calculate_e_mod_oliverpharr(obj.CantileverTips{obj.WhichTip(i)}.ProjectedTipArea,0.75);
                    if i == 1
                        obj.write_to_log_file('OliverPharr CurvePercent','0.75')
                    end
                end
                
                obj.FM{i}.calculate_adhesion_energy_and_length(2);
                obj.write_to_log_file('Adhesion Energy STD-Threshold Multiplier','2')
                obj.FM{i}.calculate_adhesion_force;
                obj.FM{i}.calculate_dissipated_and_elastic_energy;
                obj.FM{i}.calculate_peak_indentation_angle(.5);
                obj.write_to_log_file('Upper Percent of Curve considered for Peak Indentation','50%')
                obj.FM{i}.create_and_level_height_map_by_current_cp;
                
                waitbar(i/NLoop,h,sprintf('Processing Fibril %i/%i\nWrapping Up And Saving',i,NLoop));
                
                if exist('Fig')
                    for j=1:(i-1)
                        if ishandle(Fig{j})
                            close(Fig{j})
                        end
                    end
                end
                
                if TemporaryLoadInBool && obj.BigDataFlag
                    obj.FM{i}.temporary_data_load_in(false);
                    if i < NLoop
%                         obj.save_experiment;
                    end
                end
                
                if i==1
                    for k=2:NLoop
                        obj.FM{k}.CPFlag.CNNOpt = 1;
                        obj.FM{k}.MiniBatchSize = obj.FM{1}.MiniBatchSize;
                    end
                end
                
                Fig{i} = obj.FM{i}.show_analyzed_fibril();
                obj.FMFlag.FibrilAnalysis(i) = 1;
            end
            
            % Assign the Apex curves EMod and exclude +-2.5*IQR and curves
            % from ExclMask
            for i=1:NLoop
                if isequal(lower(EModOption),'hertz') || isequal(lower(EModOption),'hertzcorrected')
                    EMods = obj.FM{i}.EModHertz;
                elseif isequal(lower(EModOption),'oliver')
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
            obj.write_to_log_file('','','end')
        end
        
        function force_map_analysis_general(obj,CPOption,EModOption,BaseLineCorrectBool,TemporaryLoadInBool,UseTipInHertzBool,TiltCorrectionBool)
            % force_map_analysis_general(obj,CPOption,EModOption)
            %
            % CPOption = 'Snap-In' ... Preferred Option for data with
            % snap-in effect
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
            % EModOption = 'HertzCorrected' ... E-Modulus calculation through Hertz-Sneddon
            %                               method with RefSlope-corrected
            %                               senstitivity
            % EModOption = 'Oliver' ... E-Modulus calculation through
            % Oliver-Pharr-like method (O. Andriotis 2014)
            
            if nargin < 4
                BaseLineCorrectBool = false;
                TemporaryLoadInBool = true;
                UseTipInHertzBool = true;
                TiltCorrectionBool = true;
            elseif nargin < 5
                TemporaryLoadInBool = true;
                UseTipInHertzBool = true;
                TiltCorrectionBool = true;
            elseif nargin < 6
                UseTipInHertzBool = true;
                TiltCorrectionBool = true;
            end
            
            obj.write_to_log_file('Analysis Function','force_map_analysis_general()','start')
            obj.write_to_log_file('Contact Point Option',CPOption)
            obj.write_to_log_file('EMod Option',EModOption)
            obj.write_to_log_file('BaseLineCorrectBool',BaseLineCorrectBool)
            
            h = waitbar(0,'setting up','Units','normalized','Position',[0.4 0.3 0.2 0.1]);
            NLoop = obj.NumForceMaps;
            if sum(obj.FMFlag.ForceMapAnalysis) >= 1
                KeepFlagged = questdlg(sprintf('Some maps have been processed already.\nDo you want to skip them and keep old results?'),...
                    'Processing Options',...
                    'Yes',...
                    'No',...
                    'No');
            else
                KeepFlagged = 'No';
            end
            
            % Setting and calculating preferred method of reference slope
            obj.reference_slope_parser(1)
            
            if obj.ReferenceSlopeFlag.SetAllToValue
                RefSlopeOption = 'SetAllToValue';
            elseif obj.ReferenceSlopeFlag.UserInput
                RefSlopeOption = 'UserInput';
            elseif obj.ReferenceSlopeFlag.FromRefFM
                RefSlopeOption = 'FromRefFM';
            elseif obj.ReferenceSlopeFlag.FromArea
                RefSlopeOption = 'FromArea';
            elseif obj.ReferenceSlopeFlag.AutomaticFibril
                RefSlopeOption = 'AutomaticFibril';
            elseif obj.ReferenceSlopeFlag.Automatic
                RefSlopeOption = 'Automatic';
            end
            obj.write_to_log_file('Reference Slope Option',RefSlopeOption)
            
            % Deconvoluting cantilever tip(s)
            if isequal(lower(EModOption),'oliver') || isequal(lower(EModOption),'both') || UseTipInHertzBool
                if obj.NumCantileverTips == 0
                    if isequal(class(obj.CantileverTips{1}),'AFMImage')
                        obj.NumCantileverTips = length(obj.CantileverTips);
                        obj.CantileverTipFlag = 1;
                        obj.AssignedCantileverTips = 1;
                        obj.WhichTip = ones(obj.NumForceMaps,1);
                    else
                        Warn = warndlg('You need to load in TGT-1 images of your cantilever for this kind of analysis');
                        uiwait(Warn);
                        IsValid = false;
                        while ~IsValid
                            UsrInput = inputdlg('How many tips were used in this Experiment?');
                            NumTips = str2num(UsrInput{1});
                            IsValid = isnumeric(NumTips)&&~isempty(NumTips);
                        end
                        obj.get_paths_and_load_files([0 0 0 0 1],[0 0 0 0 ceil(NumTips)],false);
                    end
                end
                if ~obj.AssignedCantileverTips
                    obj.assign_cantilever_tips
                end
                % Check if all needed tips are deconvoluted, if not, do it
                for i=1:obj.NumForceMaps
                    if ~obj.CantileverTips{obj.WhichTip(i)}.hasDeconvolutedCantileverTip
                        obj.CantileverTips{obj.WhichTip(i)}.deconvolute_cantilever_tip;
                    end
                end
            end
            
            % Main loop for contact point estimation and E-Modulus calculation
            for i=1:NLoop
                if isequal(KeepFlagged,'Yes') && obj.FMFlag.ForceMapAnalysis(i) == 1
                    continue
                end
                
                waitbar(i/NLoop,h,sprintf('Processing Fibril %i/%i\nReading out data...',i,NLoop));
                if TemporaryLoadInBool && obj.BigDataFlag
                    obj.FM{i}.temporary_data_load_in(true);
                end
                
                waitbar(i/NLoop,h,sprintf('Processing Fibril %i/%i\nCreating and levelling Height Map',i,NLoop));
                obj.FM{i}.create_and_level_height_map
                obj.FM{i}.create_automatic_background_mask(1)
                
                Thresh = 1/2;
                AppRetSwitch = 2;
                obj.FM{i}.unselect_curves_by_fraction_of_max_data_points(Thresh,AppRetSwitch)
                
                waitbar(i/NLoop,h,sprintf('Processing ForceMap %i/%i\nFitting Base Line',i,NLoop));
                if ~obj.FM{i}.BaseAndTiltFlag
                    obj.FM{i}.base_and_tilt('linear',TiltCorrectionBool);
                if i == 1
                    obj.write_to_log_file('Baseline and Tilt option','linear')
                end
                end
                
                % contact point estimation happens here
                waitbar(i/NLoop,h,sprintf('Processing ForceMap %i/%i\nFinding Contact Point',i,NLoop));
                if BaseLineCorrectBool
                    if ~obj.FM{i}.CPFlag.CNN
                        obj.cp_option_converter('Fast',i);
                    end
                    FractionBeforeCP = .7;
                    obj.FM{i}.base_and_tilt_using_cp(FractionBeforeCP)
                    obj.write_to_log_file('FractionBeforeCP',FractionBeforeCP);
                end
                
                if ~obj.FM{i}.CPFlag.CNNZoomSweep
                obj.cp_option_converter(CPOption,i);
                end
                
                % reference slope calculation happens here
                waitbar(i/NLoop,h,sprintf('Processing ForceMap %i/%i\nProcessing and calculating Reference Slope',i,NLoop));
                obj.reference_slope_calculator(i);
                
                waitbar(i/NLoop,h,sprintf('Processing Fibril %i/%i\nCalculating E-Modulus',i,NLoop));
                if isequal(lower(EModOption),'hertz') || isequal(lower(EModOption),'both') || isequal(lower(EModOption),'hertzcorrected')
                    if isequal(lower(EModOption),'hertzcorrected') || isequal(lower(EModOption),'both')
                        CorrectSens = true;
                    else
                        CorrectSens = false;
                    end
                    AllowXShift = true;
                    obj.FM{i}.calculate_e_mod_hertz(CPOption,'parabolic',.9,AllowXShift,CorrectSens,UseTipInHertzBool,0,obj.CantileverTips{obj.WhichTip(i)});
                    if i == 1
                        obj.write_to_log_file('Hertzian Tip-Shape','parabolic')
                        obj.write_to_log_file('Hertzian CurvePercent','1')
                        obj.write_to_log_file('Allow X-Shift',AllowXShift)
                    end
                end
                if isequal(lower(EModOption),'oliver') || isequal(lower(EModOption),'both')
                    obj.FM{i}.calculate_e_mod_oliverpharr(obj.CantileverTips{obj.WhichTip(i)}.ProjectedTipArea,0.75);
                    if i == 1
                        obj.write_to_log_file('OliverPharr CurvePercent','0.75')
                    end
                end
                
                obj.FM{i}.calculate_adhesion_energy_and_length(2);
                obj.write_to_log_file('Adhesion Energy STD-Threshold Multiplier','2')
                obj.FM{i}.calculate_adhesion_force;
                obj.FM{i}.calculate_dissipated_and_elastic_energy;
                obj.FM{i}.calculate_peak_indentation_angle(.5);
                obj.write_to_log_file('Upper Percent of Curve considered for Peak Indentation','50%')
                try
                    obj.FM{i}.create_and_level_height_map_by_current_cp;
                catch
                end
                
                waitbar(i/NLoop,h,sprintf('Processing ForceMap %i/%i\nWrapping Up And Saving',i,NLoop));
                
                if TemporaryLoadInBool && obj.BigDataFlag
                    obj.FM{i}.temporary_data_load_in(false);
                    if i < NLoop
%                         obj.save_experiment;
                    end
                end
                
                obj.FMFlag.ForceMapAnalysis(i) = 1;
                if i==1
                    for k=2:NLoop
                        obj.FM{k}.CPFlag.CNNOpt = 1;
                        obj.FM{k}.MiniBatchSize = obj.FM{1}.MiniBatchSize;
                    end
                end
            end
            
            obj.save_experiment;
            
            close(h);
            obj.write_to_log_file('','','end')
        end
        
        function image_analysis_base_on_even_background(obj,UpperLim,NIter)
            
            if nargin < 2
                UpperLim = 1;
                NIter = 1;
            end
            %main loop
            h = waitbar(0,'setting up...');
            for i=1:obj.NumAFMImages
                waitbar(i/obj.NumAFMImages,h,{sprintf('Processing %i/%i:',i,obj.NumAFMImages),sprintf('%s',obj.I{i}.Name)});
                [Processed,Index] = obj.I{i}.get_channel('ProcessedSimple');
                Height = obj.I{i}.get_channel('Height (Trace)');
                if isempty(Processed)
                    Processed = Height;
                    Processed.Name = 'ProcessedSimple';
                    Index = length(obj.I{i}.Channel)+1;
                    obj.I{i}.NumChannels = Index;
                end
                Processed.Image = obj.I{i}.subtract_line_fit_hist(Height.Image, UpperLim);
                for j=1:NIter
                    Processed.Image = obj.I{i}.subtract_line_fit_hist(Processed.Image, UpperLim);
                end
                obj.I{i}.Channel(Index) = Processed;
                obj.I{i}.hasProcessed = 1;
            end
            close(h)
        end
        
        function image_analysis_flatten_on_even_background_automatic(obj,WindowSize,NIter)
            
            if nargin < 2
                WindowSize = .2;
                NIter = 3;
            end
            %main loop
            h = waitbar(0,'setting up...');
            for i=1:obj.NumAFMImages
                waitbar(i/obj.NumAFMImages,h,{sprintf('Processing %i/%i:',i,obj.NumAFMImages),sprintf('%s',obj.I{i}.Name)});
                [Processed,Index] = obj.I{i}.get_channel('Processed');
                Height = obj.I{i}.get_channel('Height (Trace)');
                if isempty(Processed)
                    Processed = Height;
                    Processed.Name = 'Processed';
                    Index = length(obj.I{i}.Channel)+1;
                    obj.I{i}.NumChannels = Index;
                end
               % Processed.Image = AFMImage.subtract_line_fit_hist(Height.Image, .5);
                for j=1:NIter
                    Processed.Image = AFMImage.subtract_line_fit_vertical_rov(Processed.Image, WindowSize,false);
                end
                obj.I{i}.Channel(Index) = Processed;
                obj.I{i}.hasProcessed = 1;
            end
            close(h)
        end
        
        function image_analysis_flatten_and_combine_trace_retrace_automatic(obj,WindowSize,NIter)
            
            if nargin < 2
                WindowSize = .2;
                NIter = 3;
            end
            
            %main loop
            h = waitbar(0,'setting up...');
            for i=1:obj.NumAFMImages
                waitbar(i/obj.NumAFMImages,h,{sprintf('Processing %i/%i:',i,obj.NumAFMImages),sprintf('%s',obj.I{i}.Name)});
                [Processed,Index] = obj.I{i}.get_channel('R-T Combined');
                T = obj.I{i}.get_channel('Height (Trace)');
                RT = obj.I{i}.get_channel('Height (Retrace)');
                if isempty(Processed)
                    Processed = T;
                    Processed.Name = 'R-T Combined';
                    Index = length(obj.I{i}.Channel)+1;
                    obj.I{i}.NumChannels = Index;
                end
                T.Image = AFMImage.subtract_line_fit_hist(T.Image,.5);
                RT.Image = AFMImage.subtract_line_fit_hist(RT.Image,.5);
                for j=1:NIter
                    T.Image = AFMImage.subtract_line_fit_vertical_rov(T.Image,WindowSize,false);
                    RT.Image = AFMImage.subtract_line_fit_vertical_rov(RT.Image,WindowSize,false);
                end
                Processed.Image = min(T.Image,RT.Image);
                for j=1:NIter
                    Processed.Image = AFMImage.subtract_line_fit_vertical_rov(Processed.Image,WindowSize,false);
                end
                obj.I{i}.Channel(Index) = Processed;
            end
            close(h)
        end
        
        function image_analysis_mask_background(obj,Thresh)
            
            if nargin < 2
                Thresh = 0;
                Automatic = 'on';
            else
                Automatic = 'off';
            end
            
            h = waitbar(0,'setting up...');
            for i=1:obj.NumAFMImages
                waitbar(i/obj.NumAFMImages,h,{sprintf('Processing %i/%i:',i,obj.NumAFMImages),sprintf('%s',obj.I{i}.Name)});
                InChannel = obj.I{i}.get_channel('Processed');
                OutChannel = InChannel;
                try
                    OutChannel.Image = AFMImage.mask_background_by_threshold(InChannel.Image,Thresh,Automatic);
                catch
                    warning("Couldn't determine threshold automatically. Taking default of 10%, instead.")
                    OutChannel.Image = AFMImage.mask_background_by_threshold(InChannel.Image,10,'off');
                end
                OutChannel.Name = 'Background Mask';
                OutChannel.Unit = 'Categorical';
                obj.I{i}.Channel(end+1) = OutChannel;
                obj.I{i}.hasBackgroundMask = true;
            end
            close(h)
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
                if isequal(lower(CPOption),'snap-in')
                    obj.FM{i}.estimate_cp_snap_in();
                end
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
                if isequal(lower(CPOption),'snap-in')
                    obj.FM{i}.estimate_cp_snap_in();
                end
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
    end     
        
    methods  % method block concerning SMFS (single molecule force spectroscopy) related data analysis 
        
        function SMFS_initialize_arrays(obj)
            % Author: Andreas Rohatschek
            % Description: 
            % The function initializes an array showing parameters of analysis selections.       
            
            %% Function body
            TableSize=[1 9];
            VarTypes = {'double','string','string','string','string','string','double','double','double'};
            VarNames = {'SMFSResults Idx','Chipbox number','Chip cantilever number','Linker','Substrate','Medium','Extend velocity','Retraction velocity','Holding time'};
            obj.SMFSResultsParameters=table('Size',TableSize,'VariableTypes',VarTypes,'VariableNames',VarNames);

        end
        
        function SMFS_initialize_flags(obj)
            % Author: Andreas Rohatschek
            % Description: 
            % The function initializes flags relevant for SM analysis.

            %% Function body
            obj.initialize_flags
            NFM = obj.NumForceMaps;

             obj.SMFSFlag.PropertiesParameters = false(NFM,1);
                         for Fm=1:obj.NumForceMaps
                              obj.FM{Fm}.initialize_flags
                         end

            NFM = obj.NumForceMaps;

            obj.SMFSFlagDown.AnalysedPreSelected = false(NFM,1);
            obj.SMFSFlagDown.AnalysedPostSelected = false(NFM,1);

            if isempty(obj.FMFlag)
                obj.SMFSFlag.AnalysedPreSelected = false(NFM,1);
                obj.SMFSFlag.AnalysedPostSelected = false(NFM,1);
            else
                PrevNFM = length(obj.FMFlag.FibrilAnalysis);
                NFM = obj.NumForceMaps;
                DiffFM = NFM - PrevNFM;
                obj.SMFSFlag.AnalysedPreSelected = false(DiffFM,1);
                obj.SMFSFlag.AnalysedPostSelected = false(DiffFM,1);
            end

        end


        function SMFS_properties_parameters(obj,SetNumFcValue)
            % Author: Andreas Rohatschek
            % Description:       
            % The function verifies the number of force curves per force map and flags each force map based on the determined value compared to the input variable 'SetNumFcValue'. 
            % If the 'SMFSFlag.NumForceCurves' is set for a force map the following part of the function is executed. There, properties and parameters, i.e. probe box, substrate, medium, velocities, from the jpk file name of the force map file are read out. 
            % Force maps that have passed through are flagged at the end of the function with the 'SMFSFlag.PropertiesParameters'.
            % Required input variables:
            % SetNumFcValue ... Defines the number of force curves per
            % force map: double
            % Flags:
            % SMFSFlag.PropertiesParameters
            % SMFSFlag.NumForceCurves
                        
            %% Function body
            h = waitbar(0,'setting up','Units','normalized','Position',[0.4 0.3 0.2 0.1]);
            NLoop = length(obj.ForceMapNames);
            if sum(obj.SMFSFlag.PropertiesParameters) >= 1
                KeepFlagged = questdlg(sprintf('Some maps have been processed already.\nDo you want to skip them and keep old results?'),...
                    'Processing Options',...
                    'Yes',...
                    'No',...
                    'No');
            else
                KeepFlagged = 'No';
            end
            for Fm=1:obj.NumForceMaps
            %  for Fm=45:obj.NumForceMaps % debugging
                 if isequal(KeepFlagged,'Yes') && obj.SMFSFlag.PropertiesParameters(Fm) == 1
                    continue
                 end
                waitbar(Fm/NLoop,h,sprintf('Preprocessing ForceMap %i/%i\nProcessing force curves',Fm,NLoop));  
                obj.FM{Fm}.fc_measurement_prop
                if obj.FM{Fm}.NCurves==SetNumFcValue
                    obj.SMFSFlag.NumForceCurves(Fm)=1;
                else
                    obj.SMFSFlag.NumForceCurves(Fm)=0;
                end 
                obj.SMFSFMParameters(Fm,:)={Fm,obj.FM{Fm}.ID,obj.FM{Fm}.Name,obj.FM{Fm}.Date,obj.FM{Fm}.Time,obj.FM{Fm}.ExtendVelocity,obj.FM{Fm}.RetractVelocity,obj.FM{Fm}.HoldingTime,obj.FM{Fm}.Linker,obj.FM{Fm}.Substrate,obj.FM{Fm}.EnvCond,obj.FM{Fm}.ChipCant,obj.FM{Fm}.Chipbox};
                waitbar(Fm/NLoop,h,sprintf('Preprocessing ForceMap %i/%i\nWrapping Up And Saving',Fm,NLoop));
            % Set flag
                obj.SMFSFlag.PropertiesParameters(Fm) = 1;
            end
            close(h);
        end
        
        function SMFS_preprocessing(obj)
            % Author: Andreas Rohatschek
            % Description:       
            % The function applies various preprocessing steps to the data.
            % First, the general cantilever sensitivity measured before starting an AFM SM experiment with the cantilever is corrected with the individual cantilever sensitivity determined from each individual force curve. 
            % Typically, in AFM experiments it is not necessary to determine the individual sensitivity but in case of the force curves detected during AFM SM experiments it was necessary to ensure most accurate results possible. 
            % Second, to correct for artifacts which occurred in some force curves and the general tilt, the approach data is fitted. Two fits, a sinoidal fit and linear fit, are applied separately onto the same amount of data and the better fit is selected. 
            % Additionally, a second linear fit is applied to correct for the tilt in the data. 
            % Third, to optimize the tilt correction for the retraction part, a linear fit of selected retraction data is applied. 
            % Fourth, the z-piezo height data is corrected with the cantilever deflection. 
            % Fifth, the contact point is estimated. 
            % Force maps that have passed through are flagged at the end of the function with the 'SMFSFlag.Preprocessed'.
            % Flags:
            % SMFSFlag.Preprocessed
                             
            %% Function body
            % Dialog box
            h = waitbar(0,'setting up','Units','normalized','Position',[0.4 0.3 0.2 0.1]);
            NLoop = length(obj.ForceMapNames);
            if sum(obj.SMFSFlag.Preprocessed) >= 1
                KeepFlagged = questdlg(sprintf('Some maps have been processed already.\nDo you want to skip them and keep old results?'),...
                    'Processing Options',...
                    'Yes',...
                    'No',...
                    'No');
            else
                KeepFlagged = 'No';
            end
            % Figure visibility
            set(groot,'defaultFigureVisible','off')      
            % set(groot,'defaultFigureVisible','on')  
            % Change into the Folder of Interest
            cd(obj.ExperimentFolder) % Move into the folder
            % Create folders for saving the produced figures
            foldername='SMFS_preprocessing';               
            mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
            currpath=fullfile(obj.ExperimentFolder,foldername);
            cd(currpath);           
            % force map loop
            
            for Fm=1:obj.NumForceMaps   
            %for Fm=263:obj.NumForceMaps  % debugging
                if isequal(KeepFlagged,'Yes') && obj.SMFSFlag.Preprocessed(Fm) == 1
                    KeepFlagged = questdlg(sprintf('Some maps have been processed already.\nDo you want to skip them and keep old results?'),...
                    'Processing Options',...
                    'Yes',...
                    'No',...
                    'No');
                else
                    KeepFlagged = 'No';
                end
                waitbar(Fm/NLoop,h,sprintf('Preprocessing ForceMap %i/%i\nProcessing force curves',Fm,NLoop));            
                obj.FM{Fm}.fc_sensitivity_correction
                obj.FM{Fm}.fc_sinoidal_fit
                obj.FM{Fm}.fc_linear_fit
                obj.FM{Fm}.fc_TipHeight_calculation
                obj.FM{Fm}.fc_estimate_cp_hardsurface
                waitbar(Fm/NLoop,h,sprintf('Preprocessing ForceMap %i/%i\nWrapping Up And Saving',Fm,NLoop));
                % Set flag
                obj.SMFSFlag.Preprocessed(Fm) = 1;
            end
            close(h);
        end
                  
      
        function SMFS_preparation(obj,xLimit1,xLimit2,xLimit3,AppThreshValue,RetThreshValue)
            % Author: Andreas Rohatschek
            % Description: 
            % In this function force curves are presorted to allow further
            % analysis. Only force maps with set SMFSFlag.Preprocessed are
            % processed. 
            % First, indices of data points corresponding to
            % three lengths defined by the input variables
            % (xLimit1,xLimit2,xLimit3) are determined.
            % Second, thresholds defined by the input variables (AppThreshValue,RetThreshValue) are used to flag force curves.
            % The flags indicate if interactions on the approach part
            % and/or retraction part are present.
            % Force maps that have passed through are flagged at the end of the function with the 'SMFSFlag.Presorted'.
            % Required input variables:
            % xLimit1 ... 1st retraction x-value limit in meters (m): double
            % xLimit2 ... Approach x-value limit in meters (m): double
            % xLimit3 ... 2nd retraction x-value limit in meters (m): double
            % AppThreshValue ... Approach threshold force value in
            % Newton (N): double
            % RetThreshValue ... Retraction threshold force value in
            % Newton (N): double
            % Flags:
            % SMFSFlag.Presorted
            % SMFSFlag.SelectFM
            % Typical input variables:  
            % Limits used for the analysis of TC            
            % xLimit1=20e-9;    % 20 nm
            % xLimit2=200e-9;   % 200 nm
            % xLimit3=300e-9;   % 300 nm
            % AppThreshValue=15e-12;     % 15 pN 
            % RetThreshValue=25e-12;     % 25 pN 
            
            %% Function body
            h = waitbar(0,'setting up','Units','normalized','Position',[0.4 0.3 0.2 0.1]);
            NLoop = length(obj.ForceMapNames);
            if sum(obj.SMFSFlag.Presorted) >= 1
                KeepFlagged = questdlg(sprintf('Some maps have been processed already.\nDo you want to skip them and keep old results?'),...
                    'Processing Options',...
                    'Yes',...
                    'No',...
                    'No');
            else
                KeepFlagged = 'No';
            end                        
            % Loop over the imported force maps
            for Fm=1:obj.NumForceMaps
            %for Fm=1:76 % Debugging
                if isequal(KeepFlagged,'Yes') && ~obj.SMFSFlag.Preprocessed(Fm)
                    continue
                end   
                waitbar(Fm/NLoop,h,sprintf('Preprocessing ForceMap %i/%i\nProcessing force curves',Fm,NLoop));
                waitbar(Fm/NLoop,h,sprintf('Preprocessing ForceMap %i/%i\nWrapping Up And Saving',Fm,NLoop));
                sprintf('Force Map No. %d of %d',Fm,obj.NumForceMaps) % Gives current Force Map Position               
                obj.FM{Fm}.fc_xLimit_idx(xLimit1,xLimit2,xLimit3)
                obj.FM{Fm}.fc_selection_threshold(AppThreshValue,RetThreshValue)
                    if nnz(obj.FM{Fm}.SMFSFlag.RetMinCrit)<20 % Only if more than 20 force curves fulfil the citeria the whole force map is considered successfully functionalized
                        obj.SMFSFlag.SelectFM(Fm)=0;
                    else
                        obj.SMFSFlag.SelectFM(Fm)=1;
                    end
                    obj.SMFSFlag.Presorted(Fm) = 1;
            end
            close(h);
        end
      
        function SMFS_visual_selection_pre(obj,XMin,XMax,YMin,YMax,NumFcMax)
            % Author: Andreas Rohatschek
            % Description: 
            % In this function force curves are plotted for visual selection via GUI (graphical user interface). 
            % Only force maps with set 'SMFSFlag.Preprocessed' are processed. 
            % Required input variables:
            % XMin ... Minimum limit of the X-axis in nanometers (nm):
            % double
            % XMax ... Maximum limit of the X-axis in nanometers (nm):
            % double
            % YMin ... Minimum limit of the Y-axis in nanoNewtons (nN):
            % double
            % YMax ... Maximum limit of the Y-axis in nanoNewtons (nN):
            % double
            % NumFcMax ... Maximum number of force curves per figure:
            % double
            % Typical input variables
            % XMin= -inf;
            % XMax= 800;     
            % YMin= -1.5;     
            % YMax= 0.100;         
            % NumFcMax = 25;  

            %% Function body
            % Figure visibility
            % set(groot,'defaultFigureVisible','off')      
            set(groot,'defaultFigureVisible','on')  
            % Change into the Folder of Interest
            cd(obj.ExperimentFolder) % Move into the folder 
            % Create folders for saving the produced figures
      %      foldername='SMFS_visual_selection';    % OLD
            foldername='SMFS_visual_selection_pre';    % Defines the folder name
            mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
            currpath=fullfile(obj.ExperimentFolder,foldername);
            cd(currpath);             
            % Loop over the imported force maps
            for fm=1:obj.NumForceMaps
            %for fm=38:53 % Debugging
                if ~obj.SMFSFlag.Preprocessed(fm)
                    continue
                end
               % Command window output
               sprintf('Force Map No. %d of %d',fm,obj.NumForceMaps) % Gives current Force Map Position
               % Call the function
               obj.FM{fm}.fc_visual_selection_flag_pre(XMin,XMax,YMin,YMax,NumFcMax)
            end    
        end
  
        function SMFS_visual_selection_post(obj,XMin,XMax,YMin,YMax,NumFcMax)
            % Author: Andreas Rohatschek
            % Description: 
            % In this function force curves are plotted, including determined force curve features, for visual selection via GUI (graphical user interface). 
            % Required input variables:
            % XMin ... Minimum limit of the X-axis in nanometers (nm):
            % double
            % XMax ... Maximum limit of the X-axis in nanometers (nm):
            % double
            % YMin ... Minimum limit of the Y-axis in nanoNewtons (nN):
            % double
            % YMax ... Maximum limit of the Y-axis in nanoNewtons (nN):
            % double
            % NumFcMax ... Maximum number of force curves per figure: double (Only
            % natural numbers are allowed)
            % Typical input variables:
            % XMin= -inf;
            % XMax= 800;     
            % YMin= -1.5;     
            % YMax= 0.100;         
            % NumFcMax = 25;

            %% Function body
            % Figure visibility
            % set(groot,'defaultFigureVisible','off')      
            set(groot,'defaultFigureVisible','on')  
            % Change into the Folder of Interest
            cd(obj.ExperimentFolder) % Move into the folder 
            % Create folders for saving the produced figures
      %      foldername='SMFS_visual_selection_analysed';    % OLD
            foldername='SMFS_visual_selection_post';    % Defines the folder name
            mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
            currpath=fullfile(obj.ExperimentFolder,foldername);
            cd(currpath);           
            % Loop over the imported force maps
            for Fm=1:obj.NumForceMaps
            %for Fm=27 % Debugging
               % Command window output
               sprintf('Force Map No. %d of %d',Fm,obj.NumForceMaps) % Gives current Force Map Position
               % Run the chosen functions
                obj.FM{Fm}.fc_visual_selection_post(XMin,XMax,YMin,YMax,NumFcMax)
            end    
        end
 
        function SMFS_analysis(obj,SelectedFlagStatus)           
            % Author: Andreas Rohatschek
            % Description: 
            % In this function the determination of the force curves features takes place. 
            % First, a possible snap-in length on the approach part is
            % determined.
            % Second, the pull-off length on the retraction part is
            % determined.
            % Third, the maximum adhesion force of the approach and
            % retraction part are determined. Additionally, the maximum
            % adhesion force of the last unbinding event is determined.
            % Fourth, the adhesion energies of the approach and retraction
            % part are determined.
            % Force maps that have passed through are flagged at the end of the function. The flagging is based on the input variable 'SelectedFlagStatus'. Either 'SMFSFlag.AnalysedPreSelected' or 'SMFSFlag.AnalysedPostSelected' is set.
            % Required input variables:
            % SelectedFlagStatus ... Selected flag status indicates if only 'SMFS_visual_selection_pre' or also 'SMFS_visual_selection_post' have been executed: string, either 'pre' or 'post' 
            % Input variable options
            % pre ... 'SMFS_visual_selection_pre' has been executed 
            % post ... 'SMFS_visual_selection_post' has been executed 
            % Flags:
            % SMFSFlag.AnalysedPreSelected
            % SMFSFlag.AnalysedPostSelected
     
            %% Function body
            h = waitbar(0,'setting up','Units','normalized','Position',[0.4 0.3 0.2 0.1]);
            NLoop = length(obj.ForceMapNames);
            switch upper(SelectedFlagStatus)
                case 'PRE'
                if sum(obj.SMFSFlag.AnalysedPreSelected) >= 1
                    KeepFlagged = questdlg(sprintf('Some maps have been processed already.\nDo you want to skip them and keep old results?'),...
                    'Processing Options',...
                    'Yes',...
                    'No',...
                    'No');
                else
                KeepFlagged = 'No';
                end    
                case 'POST'
                if sum(obj.SMFSFlag.AnalysedPostSelected) >= 1
                KeepFlagged = questdlg(sprintf('Some maps have been processed already.\nDo you want to skip them and keep old results?'),...
                    'Processing Options',...
                    'Yes',...
                    'No',...
                    'No');
                else
                KeepFlagged = 'No';
                end   
            end
            % Figure visibility
            set(groot,'defaultFigureVisible','off')      
            % set(groot,'defaultFigureVisible','on') 
            %% Loop
            for Fm=1:obj.NumForceMaps
            %for Fm=27:36 % Debugging    
            if isequal(KeepFlagged,'Yes') 
                    continue
            end   
               waitbar(Fm/NLoop,h,sprintf('Preprocessing ForceMap %i/%i\nProcessing force curves',Fm,NLoop));
               waitbar(Fm/NLoop,h,sprintf('Preprocessing ForceMap %i/%i\nWrapping Up And Saving',Fm,NLoop));
               sprintf('Force Map No. %d of %d',Fm,obj.NumForceMaps) % Gives current Force Map Position   
               % Print force curves containing label for the pulling length
               % and colored area for the adhesion energy                              
               % Snap-in
               obj.FM{Fm}.fc_snap_in_length_MAD(SelectedFlagStatus)
               % Pulling length
               obj.FM{Fm}.fc_pulling_length_MAD(SelectedFlagStatus)
               % Maximum adhesion force
               obj.FM{Fm}.fc_adh_force_max(SelectedFlagStatus)
               % Adhesion energy
               obj.FM{Fm}.fc_adhesion_energy_idxlength(SelectedFlagStatus)    
               % Flag
               %% General variables 1
           if strcmpi(SelectedFlagStatus,'Pre')
               obj.SMFSFlag.AnalysedPreSelected(Fm) = 1;
           elseif strcmpi(SelectedFlagStatus,'Post')
               obj.SMFSFlag.AnalysedPostSelected(Fm) = 1;
           end
            end
            close(h);
        end
                 
     
                       
        function SMFS_print_raw(obj,XMin,XMax,YMin,YMax)
            % Author: Andreas Rohatschek  
            % Description:
            % In this function all force curves of all
            % force maps are plotted in their raw form.
            % Required input variables:
            % XMin ... Minimum limit of the X-axis in meters (nm): double
            % XMax ... Maximum limit of the X-axis in meters (nm): double
            % YMin ... Minimum limit of the Y-axis in Newtons (nN): double
            % YMax ... Maximum limit of the Y-axis in Newtons (nN): double
            % Typical input variables:
            % XMin= -10e-9;
            % XMax= 500e-9;
            % YMin= -1e-9;
            % YMax= 100e-12;

            %% Function body
            if nargin<2
                XMin= -inf;     % Limit of the X-axis in meters (m)
                XMax= inf;      % Limit of the X-axis in meters (m)
                YMin= -inf;     % Limit of the Y-axis in Newtons (N)
                YMax= inf;      % Limit of the Y-axis in Newtons (N)
            else

                % Figure visibility
                set(groot,'defaultFigureVisible','off')
                % set(groot,'defaultFigureVisible','on')
                % Change into the Folder of Interest
                % Change into the Folder of Interest
                cd(obj.ExperimentFolder) % Move into the folder
                % Create folders for saving the produced figures
                foldername='SMFS_print_raw';    % Defines the folder name
                mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
                currpath=fullfile(obj.ExperimentFolder,foldername);
                cd(currpath);

                % Loop over the imported force maps
                for Fm=1:obj.NumForceMaps
                    %for Fm=1:50 % Debugging
                    % Command window output
                    sprintf('Force Map No. %d of %d',Fm,obj.NumForceMaps) % Gives current Force Map Position
                    % Run the chosen functions
                    obj.FM{Fm}.fc_print_raw(XMin,XMax,YMin,YMax);
                end
            end
        end

        function SMFS_print_sort(obj,StartDate,EndDate,XMin,XMax,YMin,YMax,Flags)
            % Author: Andreas Rohatschek
            % Description: 
            % In this function force maps ca be selected by date and time and plotted 
            % Required input variables:
            % StartDate ... Starting date of the selection (Input format:
            % 'YYYY.MM.DD')
            % EndDate ... Ending date of the selection (Input format:
            % 'YYYY.MM.DD')
            % XMin ... Minimum limit of the X-axis in nanometers (nm):
            % double
            % XMax ... Maximum limit of the X-axis in nanometers (nm):
            % double
            % YMin ... Minimum limit of the Y-axis in nanoNewtons (nN):
            % double
            % YMax ... Maximum limit of the Y-axis in nanoNewtons (nN):
            % double
            % Flags ... Flags that have been set already: string either
            % 'presorted' or 'fit' 
            % Typical input variables:
            % StartDate= '0000.00.00'
            % EndDate= '2050.00.00'
            % XMin= -500;
            % XMax= 10;     
            % YMin= -1.5;     
            % YMax= 0.1;         
            % Flags = 'fit';

            %% Function body
            % Figure visibility
            set(groot,'defaultFigureVisible','off')      
            % set(groot,'defaultFigureVisible','on')  
            % Change into the Folder of Interest
            cd(obj.ExperimentFolder) % Move into the folder 
            % Create folders for saving the produced figures
            foldername='SMFS_print_sort';    % Defines the folder name
            mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
            currpath=fullfile(obj.ExperimentFolder,foldername);
            cd(currpath); 
            % Input variable adaptation
            if nargin<2
                StartDate='0000.00.00';
                EndDate='2999.00.00';
                XMin= -inf;     % Limit of the X-axis in nanometers (m)  
                XMax= inf;      % Limit of the X-axis in nanometers (m)
                YMin= -inf;     % Limit of the Y-axis in nanoNewtons (N)   
                YMax= inf;      % Limit of the Y-axis in nanoNewtons (N)
            elseif nargin<6    
                StartDate='0000.00.00';
                EndDate='2999.00.00';
            end
            
            % Output time and date for the dairy
            datetime('now')
            
            % Loop over the imported force maps
             for Fm=1:obj.NumForceMaps
             %for Fm=1
                 % Needed function               
                %if ~obj.SMFSFlag(ii)     % Selects all flagged 1 force maps
                %if obj.SMFSFlag(ii)     % Selects all flagged 0 force maps
                %    continue
                %end
                % Remove the dots in the dates               
                FMDate=split(obj.FM{Fm}.Date,'.');
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
                StartDateMod=strrep(StartDate,'.','');
                EndDateMod=strrep(EndDate,'.','');
                %foldername=append('FM_Flag',SMFSFlagConvert,'_',VelocityConvert,'_',obj.FM{ii}.Substrate,'_',obj.FM{ii}.EnvCond,'_',StartDateMod,'-',EndDateMod); % Defines the folder name
                %foldername=append(obj.FM{ii}.Substrate,'_',obj.FM{ii}.EnvCond,'_',StartDateMod,'-',EndDateMod); % Defines the folder name 
                foldername=append(StartDateMod,'-',EndDateMod,'_',obj.FM{Fm}.Chipbox,'_',obj.FM{Fm}.ChipCant); % Defines the folder name 
                warning('off','all'); % Set warning off to prevent from showing the same warning each loop 
                mkdir(foldername);
                warning('on','all');
                cd(foldername)         
               % Run the chosen function  
               obj.FM{Fm}.fc_print_fitted(XMin,XMax,YMin,YMax,Flags)
               cd(obj.ExperimentFolder) % Move into the folder                             
            end            
        end
        
        function SMFS_print_analysed_fc(obj,XMin,XMax,YMin,YMax,NumFcMax)
            % Author: Andreas Rohatschek
            % Description: 
            % In this function all analysed force curves and the corresponding determined
            % features are plotted
            % Required input variables:
            % XMin ... Minimum limit of the X-axis in nanometers (nm):
            % double
            % XMax ... Maximum limit of the X-axis in nanometers (nm):
            % double
            % YMin ... Minimum limit of the Y-axis in nanoNewtons (nN):
            % double
            % YMax ... Maximum limit of the Y-axis in nanoNewtons (nN):
            % double
            % NumFcMax ... Maximum number of force curves per figure: double (Only
            % natural numbers are allowed)
            % Typical input variables:
            % XMin= -10;
            % XMax= 500;     
            % YMin= -1.0;     
            % YMax= 0.1;
            % NumFcMax=25

            %% Function body
            if nargin<2
                XMin= -inf;     % Limit of the X-axis in meters (m)
                XMax= 50e-9;      % Limit of the X-axis in meters (m)
                YMin= -inf;     % Limit of the Y-axis in Newtons (N)
                YMax= 100e-12;      % Limit of the Y-axis in Newtons (N)
                NumFcMax = 25;   % Maximum number of force curves per figure
            elseif nargin<3
                XMin= -inf;     % Limit of the X-axis in meters (m)
                XMax= inf;      % Limit of the X-axis in meters (m)
                YMin= -inf;     % Limit of the Y-axis in Newtons (N)
                YMax= inf;      % Limit of the Y-axis in Newtons (N)
            end

            % Output time and date for the dairy
            datetime('now')
            
            % Figure visibility
            set(groot,'defaultFigureVisible','off')
        %    set(groot,'defaultFigureVisible','on')
            cd(obj.ExperimentFolder) % Move into the folder
            % Create folders for saving the produced figures
            %foldername='FM_Test';    % for debugging
            foldername='SMFS_print_analysed_fc';    % Defines the folder name
            mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
            currpath=fullfile(obj.ExperimentFolder,foldername);
            cd(currpath);
            %% loop
            for Fm=1:obj.NumForceMaps
            %for Fm=7:10 % Debugging
            %sprintf('Force Map No. %d of %d',hh,obj.NumForceMaps) % Gives current Force Map Position   
            if ~obj.SMFSFlag.Preprocessed(Fm)
                    continue
            end   
            % Determine needed input variable
               NumFcUncorrupt(Fm)=nnz(obj.FM{Fm}.SMFSFlag.Uncorrupt); % Determine the number of uncorrupted force curves     
               obj.FM{Fm}.fc_print_analysed(XMin,XMax,YMin,YMax,NumFcMax)         
            end
            obj.NumFcUncorrupt=NumFcUncorrupt;
        end
        

        function SMFS_results_structure(obj,ChipboxValue,ChipCantValue,LinkerValue,SubstrateValue,EnvCondValue,ExtVelocityValue,RetVelocityValue,HoldingTimeValue)
            % Author: Andreas Rohatschek
            % Description: 
            % In this function a results structure entry is created based on the
            % input variables (properties and parameters)
            % Required input variables:
            % ChipboxValue ... Chipbox number: string
            % ChipCantValue ... Chip and cantilever number and letter:
            % string
            % LinkerValue ... Linker type: string , either 'long' or 'short'
            % SubstrateValue ... Subtstrate value: string
            % EnvCondValue ... Media (environmental conditions) of the
            % experiment: string
            % ExtVelocityValue ... Velocity value of the approach part
            % (extend): double
            % RetVelocityValue ... Velocity value of the retraction part: double 
            % HoldingTimeValue ... Holding time (dwell time) value: double
            % Typical input variables:
            % ChipboxValue ... 'XXX'
            % ChipCantValue ... '5E'
            % LinkerValue ... 'short'
            % SubstrateValue ... 'mica'
            % EnvCondValue ... 'water'
            % ExtVelocityValue ... 0 (selects all velocities)
            % RetVelocityValue ...  0 (selects all velocities)  
            % HoldingTimeValue ... -1 (selects all holding times) 

            %% Function body
            % Define variables
            FMIdxArray=[];
            jj=1;
            DateFormat='yyyy-MM-dd HH-mm-ss-SSS';
            for Fm=1:obj.NumForceMaps
                %% Debugging
                %for Fm=65:117 % for debugging
                sprintf('Force map No. %d',Fm) % Gives current force map
                %% Force map selection criteria
                if ~obj.SMFSFlag.AnalysedPostSelected(Fm) || ~obj.SMFSFlag.NumForceCurves(Fm)    % Exclude force map if analysis has not been done
%                if ~obj.SMFSFlag.Analysed(Fm) || ~obj.SMFSFlag.NumForceCurves(Fm)    % Exclude force map if analysis has not been done
                    continue
                end
                % Parameters
                if ((round(obj.FM{Fm}.ExtendVelocity,8)==ExtVelocityValue || ExtVelocityValue==0) ...
                        && (round(obj.FM{Fm}.RetractVelocity,8)==RetVelocityValue || RetVelocityValue==0) ...
                        && (round(obj.FM{Fm}.HoldingTime,2)==HoldingTimeValue || HoldingTimeValue==-1) ...    % Round holding time value to correct for tiny deviations of decimal places after the comma origin from instrument (AFM)
                        && (strcmpi(obj.FM{Fm}.Substrate,SubstrateValue) || strcmpi(SubstrateValue,'All')) ...
                        && (strcmpi(obj.FM{Fm}.EnvCond,EnvCondValue) || strcmpi(EnvCondValue,'All')) ...
                        && (strcmpi(obj.FM{Fm}.ChipCant,ChipCantValue) || strcmpi(ChipCantValue,'All')) ...
                        && (strcmpi(obj.FM{Fm}.Chipbox,ChipboxValue) || strcmpi(ChipboxValue,'All')) ...
                        && (strcmpi(obj.FM{Fm}.Linker,LinkerValue) || strcmpi(LinkerValue,'All')))
                    % Define variables for the if condition
                    FMIdxArray(jj,1)=Fm;
                    % Adjust variable
                    jj=jj+1;
                end
            end
            % If condition to handle an empty index array
            if isempty(FMIdxArray)
                return
            else
            end
            % Preallocate
            ConcateArrayAdhMaxApp=zeros(length(FMIdxArray)*obj.FM{FMIdxArray(1)}.NCurves,1);
            ConcateArrayAdhMaxRet=zeros(length(FMIdxArray)*obj.FM{FMIdxArray(1)}.NCurves,1);
            ConcateArrayAdhMaxRetUnbinding=zeros(length(FMIdxArray)*obj.FM{FMIdxArray(1)}.NCurves,1);
            ConcateArrayAdhEneApp=zeros(length(FMIdxArray)*obj.FM{FMIdxArray(1)}.NCurves,1);
            ConcateArrayAdhEneRet=zeros(length(FMIdxArray)*obj.FM{FMIdxArray(1)}.NCurves,1);
            ConcateArrayPullLength=zeros(length(FMIdxArray)*obj.FM{FMIdxArray(1)}.NCurves,1);
            ConcateArraySnapIn=zeros(length(FMIdxArray)*obj.FM{FMIdxArray(1)}.NCurves,1);
            yAdhMaxAppAll=zeros(obj.FM{FMIdxArray(1)}.NCurves,length(FMIdxArray));
            FCAnalysedyAdhMaxApp=zeros(1,length(FMIdxArray));
            yAdhMaxRetAll=zeros(obj.FM{FMIdxArray(1)}.NCurves,length(FMIdxArray));
            FCAnalysedyAdhMaxRet=zeros(1,length(FMIdxArray));
            yAdhUnbindingAll=zeros(obj.FM{FMIdxArray(1)}.NCurves,length(FMIdxArray));
            FCAnalysedyAdhUnbind=zeros(1,length(FMIdxArray));
            yAdhEneAppAll=zeros(obj.FM{FMIdxArray(1)}.NCurves,length(FMIdxArray));
            FCAnalysedyAdhEneApp=zeros(1,length(FMIdxArray));
            yAdhEneRetAll=zeros(obj.FM{FMIdxArray(1)}.NCurves,length(FMIdxArray));
            FCAnalysedyAdhEneRet=zeros(1,length(FMIdxArray));
            yPullingLengthAll=zeros(obj.FM{FMIdxArray(1)}.NCurves,length(FMIdxArray));
            FCAnalysedyPullingLength=zeros(1,length(FMIdxArray));
            ySnapInLengthAll=zeros(obj.FM{FMIdxArray(1)}.NCurves,length(FMIdxArray));
            FCAnalysedySnapInLength=zeros(1,length(FMIdxArray));            
            FMID=cell(length(FMIdxArray),1);
            FMNum=zeros(length(FMIdxArray),1);
            FMExtVelocity=zeros(length(FMIdxArray),1);
            FMRetVelocity=zeros(length(FMIdxArray),1);
            FMHoldingTime=zeros(length(FMIdxArray),1);
            FMSubstrate=cell(length(FMIdxArray),1);
            FMEnvCond=cell(length(FMIdxArray),1);
            FMChipCant=cell(length(FMIdxArray),1);
            FMChipbox=cell(length(FMIdxArray),1);
            FMLinker=cell(length(FMIdxArray),1);
            FMDateTime=cell(length(FMIdxArray),1);
            FMDateTimeNumber=zeros(length(FMIdxArray),1);
            TotalNumFc=zeros(length(FMIdxArray),1);            
            SMFSFlagUncorrupt=zeros(length(FMIdxArray),1);
            SMFSFlagSelected=zeros(length(FMIdxArray),1);
            SMFSFlagAppMinCrit=zeros(length(FMIdxArray),1);
            SMFSFlagRetMinCrit=zeros(length(FMIdxArray),1);
            SMFSFlagLengthRequisite=zeros(length(FMIdxArray),1);
            SMFSFlagFit=zeros(length(FMIdxArray),1);
            SMFSFlagFitLinear=zeros(length(FMIdxArray),1);
            SMFSFlagFitSinoidal=zeros(length(FMIdxArray),1);
            SMFSFlagSnapIn=zeros(length(FMIdxArray),1);
            SMFSFlagPullingLength=zeros(length(FMIdxArray),1);
            AdhMaxAppMinArray=zeros(length(FMIdxArray),1);
            AdhMaxAppMinFcArray=zeros(length(FMIdxArray),1);
            AdhMaxAppMaxArray=zeros(length(FMIdxArray),1);
            AdhMaxAppMaxFcArray=zeros(length(FMIdxArray),1);
            AdhMaxRetMinArray=zeros(length(FMIdxArray),1);
            AdhMaxRetMinFcArray=zeros(length(FMIdxArray),1);
            AdhMaxRetMaxArray=zeros(length(FMIdxArray),1);
            AdhMaxRetMaxFcArray=zeros(length(FMIdxArray),1);
            AdhMaxRetUnbindingMinArray=zeros(length(FMIdxArray),1);
            AdhMaxRetUnbindingMinFcArray=zeros(length(FMIdxArray),1);
            AdhMaxRetUnbindingMaxArray=zeros(length(FMIdxArray),1);
            AdhMaxRetUnbindingMaxFcArray=zeros(length(FMIdxArray),1);
            AdhEneAppMinArray=zeros(length(FMIdxArray),1);
            AdhEneAppMinFcArray=zeros(length(FMIdxArray),1);
            AdhEneAppMaxArray=zeros(length(FMIdxArray),1);
            AdhEneAppMaxFcArray=zeros(length(FMIdxArray),1);
            AdhEneRetMinArray=zeros(length(FMIdxArray),1);
            AdhEneRetMinFcArray=zeros(length(FMIdxArray),1);
            AdhEneRetMaxArray=zeros(length(FMIdxArray),1);
            AdhEneRetMaxFcArray=zeros(length(FMIdxArray),1);
            PullLengthMinArray=zeros(length(FMIdxArray),1);
            PullLengthMinFcArray=zeros(length(FMIdxArray),1);
            PullLengthMaxArray=zeros(length(FMIdxArray),1);
            PullLengthMaxFcArray=zeros(length(FMIdxArray),1);
            SnapInMinArray=zeros(length(FMIdxArray),1);
            SnapInMinFcArray=zeros(length(FMIdxArray),1);
            SnapInMaxArray=zeros(length(FMIdxArray),1);
            SnapInMaxFcArray=zeros(length(FMIdxArray),1);
            % Loop
            for ff=1:length(FMIdxArray)
                %% Debugging
                %for ff=6 % for debugging
                %sprintf('Index array row No. %d',ff) % Gives current Force curve
                % Allocate data
                yAdhMaxApp=obj.FM{FMIdxArray(ff)}.AdhForceMaxApp;
                yAdhMaxRet=obj.FM{FMIdxArray(ff)}.AdhForceMaxRet;
                yAdhUnbinding=obj.FM{FMIdxArray(ff)}.AdhForceUnbinding;
                yAdhEneApp=obj.FM{FMIdxArray(ff)}.AppAdhEnergy_IdxMethod;
                yAdhEneRet=obj.FM{FMIdxArray(ff)}.RetAdhEnergy_IdxMethod;
                yPullingLength=obj.FM{FMIdxArray(ff)}.PullingLength;
                ySnapInLength=obj.FM{FMIdxArray(ff)}.SnapInLength;
                FMID{ff,1}=obj.FM{FMIdxArray(ff)}.ID;
                FMNum(ff,1)=ff;
                FMExtVelocity(ff,1)=round(obj.FM{FMIdxArray(ff)}.ExtendVelocity,8); % Round holding time value to correct for tiny deviations of decimal places after the comma origin from instrument (AFM)
                FMRetVelocity(ff,1)=round(obj.FM{FMIdxArray(ff)}.RetractVelocity,8); % Round holding time value to correct for tiny deviations of decimal places after the comma origin from instrument (AFM)
                FMHoldingTime(ff,1)=round(obj.FM{FMIdxArray(ff)}.HoldingTime,2); % Round holding time value to correct for tiny deviations of decimal places after the comma origin from instrument (AFM)
                FMSubstrate{ff,1}=obj.FM{FMIdxArray(ff)}.Substrate;
                FMEnvCond{ff,1}=obj.FM{FMIdxArray(ff)}.EnvCond;
                FMChipCant{ff,1}=obj.FM{FMIdxArray(ff)}.ChipCant;
                FMChipbox{ff,1}=obj.FM{FMIdxArray(ff)}.Chipbox;
                FMLinker{ff,1}=obj.FM{FMIdxArray(ff)}.Linker;
                FMDate=obj.FM{FMIdxArray(ff)}.Date;
                FMTime=obj.FM{FMIdxArray(ff)}.Time;
                DateTimeStr=[FMDate,' ',FMTime];
                FMDateTime{ff,1}=cellstr(datetime(DateTimeStr,'InputFormat',DateFormat,'Format',DateFormat));
                FMDateTimeNumber(ff,1)=datenum(datetime(DateTimeStr,'InputFormat',DateFormat,'Format',DateFormat));
                TotalNumFc(ff,1)=obj.FM{FMIdxArray(ff)}.NCurves;
                % SMFS Flags
                SMFSFlagUncorrupt(ff,1)=nnz(obj.FM{FMIdxArray(ff)}.SMFSFlag.Uncorrupt');
                SMFSFlagSelected(ff,1)=nnz(obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected');
                SMFSFlagAppMinCrit(ff,1)=nnz(obj.FM{FMIdxArray(ff)}.SMFSFlag.AppMinCrit');
                SMFSFlagRetMinCrit(ff,1)=nnz(obj.FM{FMIdxArray(ff)}.SMFSFlag.RetMinCrit');
                SMFSFlagLengthRequisite(ff,1)=nnz(obj.FM{FMIdxArray(ff)}.SMFSFlag.LengthRequisite');
                SMFSFlagFit(ff,1)=nnz(obj.FM{FMIdxArray(ff)}.SMFSFlag.Fit');
                SMFSFlagFitLinear(ff,1)=nnz(obj.FM{FMIdxArray(ff)}.SMFSFlag.FitLinear');
                SMFSFlagFitSinoidal(ff,1)=nnz(obj.FM{FMIdxArray(ff)}.SMFSFlag.FitSinoidal');
                SMFSFlagSnapIn(ff,1)=nnz(obj.FM{FMIdxArray(ff)}.SMFSFlag.SnapIn');
                SMFSFlagPullingLength(ff,1)=nnz(obj.FM{FMIdxArray(ff)}.SMFSFlag.PullingLength');
                %% Concatenate arrays
                % FCs of each FM in seperate column
                yAdhMaxAppAll(:,ff)=yAdhMaxApp'.*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected';
                yAdhMaxAppAll(yAdhMaxAppAll==0)=nan; % Replace zero entries by nan´s
                FCAnalysedyAdhMaxApp(1,ff)=nnz(~isnan(yAdhMaxAppAll(:,ff)));
                yAdhMaxRetAll(:,ff)=yAdhMaxRet'.*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected';
                yAdhMaxRetAll(yAdhMaxRetAll==0)=nan; % Replace zero entries by nan´s
                FCAnalysedyAdhMaxRet(1,ff)=nnz(~isnan(yAdhMaxRetAll(:,ff)));
                yAdhUnbindingAll(:,ff)=yAdhUnbinding'.*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected';
                yAdhUnbindingAll(yAdhUnbindingAll==0)=nan; % Replace zero entries by nan´s
                FCAnalysedyAdhUnbind(1,ff)=nnz(~isnan(yAdhUnbindingAll(:,ff)));
                yAdhEneAppAll(:,ff)=yAdhEneApp'.*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected';
                yAdhEneAppAll(yAdhEneAppAll==0)=nan; % Replace zero entries by nan´s
                FCAnalysedyAdhEneApp(1,ff)=nnz(~isnan(yAdhEneAppAll(:,ff)));
                yAdhEneRetAll(:,ff)=yAdhEneRet'.*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected';
                yAdhEneRetAll(yAdhEneRetAll==0)=nan; % Replace zero entries by nan´s
                FCAnalysedyAdhEneRet(1,ff)=nnz(~isnan(yAdhEneRetAll(:,ff)));
                yPullingLengthAll(:,ff)=yPullingLength'.*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected';
                yPullingLengthAll(yPullingLengthAll==0)=nan; % Replace zero entries by nan´s
                FCAnalysedyPullingLength(1,ff)=nnz(~isnan(yPullingLengthAll(:,ff)));
                ySnapInLengthAll(:,ff)=ySnapInLength'.*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected';
                ySnapInLengthAll(ySnapInLengthAll==0)=nan; % Replace zero entries by nan´s
                FCAnalysedySnapInLength(1,ff)=nnz(~isnan(ySnapInLengthAll(:,ff)));
                % All FCs of all FM in one column                   
                if ~isempty(yAdhMaxApp)
                    % Determine the number of rows per force map
                    ArrayLength=length(yAdhMaxApp); % Define the length of the array
                    row_start = ((ff-1) * ArrayLength) + 1; % Define the appropriate row start to append the new data
                    row_end   = ff * ArrayLength; % Define the appropriate row end to append the new data
                    % Concatenated data
                    ConcateArrayAdhMaxApp(row_start:row_end,:)=yAdhMaxApp'; % Append the new data into the concatenated vector
                    ConcateArrayAdhMaxApp(row_start:row_end,:)=ConcateArrayAdhMaxApp(row_start:row_end,:).*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected'; % Set non-selected force curves from the concatenated arrays to zero
                    ConcateArrayAdhMaxApp(ConcateArrayAdhMaxApp==0)=nan; % Replace zero entries by nan´s
                    % Allocate parameters
                    FMIDArray(row_start:row_end,:)={FMID}; % Allocate the FM ID to each row
                    FMIndexArray(row_start:row_end,:)=FMIdxArray(ff);
                    FMNumArray(row_start:row_end,:)=ff;
                    FcNumArray(row_start:row_end,:)=(row_start:row_end)';
                    FMExtVelocityArray(row_start:row_end,:)=FMExtVelocity(ff,1);
                    FMRetVelocityArray(row_start:row_end,:)=FMRetVelocity(ff,1);
                    FMHoldingTimeArray(row_start:row_end,:)=FMHoldingTime(ff,1);
                    FMSubstrateArray(row_start:row_end,:)=FMSubstrate(ff,1);
                    FMEnvCondArray(row_start:row_end,:)=FMEnvCond(ff,1);
                    FMChipCantArray(row_start:row_end,:)=FMChipCant(ff,1);
                    FMChipboxArray(row_start:row_end,:)=FMChipbox(ff,1);
                    FMLinkerArray(row_start:row_end,:)=FMLinker(ff,1);
                    FMDateTimeFcArray(row_start:row_end,:)=datetime(DateTimeStr,'InputFormat',DateFormat,'Format',DateFormat);
                    FMDateTimeNumberFcArray(row_start:row_end,:)=datenum(datetime(DateTimeStr,'InputFormat',DateFormat,'Format',DateFormat));
                else
                end
                if ~isempty(yAdhMaxRet)
                    % Determine the number of rows per force map
                    ArrayLength=length(yAdhMaxRet); % Define the length of the array
                    row_start = ((ff-1) * ArrayLength) + 1; % Define the appropriate row start to append the new data
                    row_end   = ff * ArrayLength; % Define the appropriate row end to append the new data
                    % Concatenated data
                    ConcateArrayAdhMaxRet(row_start:row_end,:)=yAdhMaxRet'; % Append the new data into the concatenated vector
                    ConcateArrayAdhMaxRet(row_start:row_end,:)=ConcateArrayAdhMaxRet(row_start:row_end,:).*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected'; % Set non-selected force curves from the concatenated arrays to zero
                    ConcateArrayAdhMaxRet(ConcateArrayAdhMaxRet==0)=nan; % Replace zero entries by nan´s
                else
                end
                if ~isempty(yAdhUnbinding)
                    % Determine the number of rows per force map
                    ArrayLength=length(yAdhUnbinding); % Define the length of the array
                    row_start = ((ff-1) * ArrayLength) + 1; % Define the appropriate row start to append the new data
                    row_end   = ff * ArrayLength; % Define the appropriate row end to append the new data
                    % Concatenated data
                    ConcateArrayAdhMaxRetUnbinding(row_start:row_end,:)=yAdhUnbinding'; % Append the new data into the concatenated vector
                    ConcateArrayAdhMaxRetUnbinding(row_start:row_end,:)=ConcateArrayAdhMaxRetUnbinding(row_start:row_end,:).*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected'; % Set non-selected force curves from the concatenated arrays to zero
                    ConcateArrayAdhMaxRetUnbinding(ConcateArrayAdhMaxRetUnbinding==0)=nan; % Replace zero entries by nan´s
                else
                end
                if ~isempty(yAdhEneApp)
                    % Determine the number of rows per force map
                    ArrayLength=length(yAdhEneApp); % Define the length of the array
                    row_start = ((ff-1) * ArrayLength) + 1; % Define the appropriate row start to append the new data
                    row_end   = ff * ArrayLength; % Define the appropriate row end to append the new data
                    % Concatenated data
                    ConcateArrayAdhEneApp(row_start:row_end,:)=yAdhEneApp'; % Append the new data into the concatenated vector
                    ConcateArrayAdhEneApp(row_start:row_end,:)=ConcateArrayAdhEneApp(row_start:row_end,:).*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected'; % Set non-selected force curves from the concatenated arrays to zero
                    ConcateArrayAdhEneApp(ConcateArrayAdhEneApp==0)=nan; % Replace zero entries by nan´s
                else
                end
                if ~isempty(yAdhEneRet)
                    % Determine the number of rows per force map
                    ArrayLength=length(yAdhEneRet); % Define the length of the array
                    row_start = ((ff-1) * ArrayLength) + 1; % Define the appropriate row start to append the new data
                    row_end   = ff * ArrayLength; % Define the appropriate row end to append the new data
                    % Concatenated data
                    ConcateArrayAdhEneRet(row_start:row_end,:)=yAdhEneRet'; % Append the new data into the concatenated vector
                    ConcateArrayAdhEneRet(row_start:row_end,:)=ConcateArrayAdhEneRet(row_start:row_end,:).*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected'; % Set non-selected force curves from the concatenated arrays to zero
                    ConcateArrayAdhEneRet(ConcateArrayAdhEneRet==0)=nan; % Replace zero entries by nan´s
                else
                end
                if ~isempty(yPullingLength)
                    % Determine the number of rows per force map
                    ArrayLength=length(yPullingLength); % Define the length of the array
                    row_start = ((ff-1) * ArrayLength) + 1; % Define the appropriate row start to append the new data
                    row_end   = ff * ArrayLength; % Define the appropriate row end to append the new data
                    % Concatenated data
                    ConcateArrayPullLength(row_start:row_end,:)=yPullingLength'; % Append the new data into the concatenated vector
                    ConcateArrayPullLength(row_start:row_end,:)=ConcateArrayPullLength(row_start:row_end,:).*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected'; % Set non-selected force curves from the concatenated arrays to zero
                    ConcateArrayPullLength(ConcateArrayPullLength==0)=nan; % Replace zero entries by nan´s
                else
                end
                if ~isempty(ySnapInLength)
                    % Determine the number of rows per force map
                    ArrayLength=length(ySnapInLength); % Define the length of the array
                    row_start = ((ff-1) * ArrayLength) + 1; % Define the appropriate row start to append the new data
                    row_end   = ff * ArrayLength; % Define the appropriate row end to append the new data
                    % Concatenated data
                    ConcateArraySnapIn(row_start:row_end,:)=ySnapInLength'; % Append the new data into the concatenated vector
                    ConcateArraySnapIn(row_start:row_end,:)=ConcateArraySnapIn(row_start:row_end,:).*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected'; % Set non-selected force curves from the concatenated arrays to zero
                    ConcateArraySnapIn(ConcateArraySnapIn==0)=nan; % Replace zero entries by nan´s
                else
                end
                % Min max values and location
                if ff==1
                    [AdhMaxAppMinArray(ff,1),AdhMaxAppMinFcArray(ff,1)]=min(ConcateArrayAdhMaxApp(1:100),[],'omitnan');
                    [AdhMaxAppMaxArray(ff,1),AdhMaxAppMaxFcArray(ff,1)]=max(ConcateArrayAdhMaxApp(1:100),[],'omitnan');
                    [AdhMaxRetMinArray(ff,1),AdhMaxRetMinFcArray(ff,1)]=min(ConcateArrayAdhMaxRet(1:100),[],'omitnan');
                    [AdhMaxRetMaxArray(ff,1),AdhMaxRetMaxFcArray(ff,1)]=max(ConcateArrayAdhMaxRet(1:100),[],'omitnan');
                    [AdhMaxRetUnbindingMinArray(ff,1),AdhMaxRetUnbindingMinFcArray(ff,1)]=min(ConcateArrayAdhMaxRetUnbinding(1:100),[],'omitnan');
                    [AdhMaxRetUnbindingMaxArray(ff,1),AdhMaxRetUnbindingMaxFcArray(ff,1)]=max(ConcateArrayAdhMaxRetUnbinding(1:100),[],'omitnan');
                    [AdhEneAppMinArray(ff,1),AdhEneAppMinFcArray(ff,1)]=min(ConcateArrayAdhEneApp(1:100),[],'omitnan');
                    [AdhEneAppMaxArray(ff,1),AdhEneAppMaxFcArray(ff,1)]=max(ConcateArrayAdhEneApp(1:100),[],'omitnan');
                    [AdhEneRetMinArray(ff,1),AdhEneRetMinFcArray(ff,1)]=min(ConcateArrayAdhEneRet(1:100),[],'omitnan');
                    [AdhEneRetMaxArray(ff,1),AdhEneRetMaxFcArray(ff,1)]=max(ConcateArrayAdhEneRet(1:100),[],'omitnan');
                    [PullLengthMinArray(ff,1),PullLengthMinFcArray(ff,1)]=min(ConcateArrayPullLength(1:100),[],'omitnan');
                    [PullLengthMaxArray(ff,1),PullLengthMaxFcArray(ff,1)]=max(ConcateArrayPullLength(1:100),[],'omitnan');
                    [SnapInMinArray(ff,1),SnapInMinFcArray(ff,1)]=min(ConcateArraySnapIn(1:100),[],'omitnan');
                    [SnapInMaxArray(ff,1),SnapInMaxFcArray(ff,1)]=max(ConcateArraySnapIn(1:100),[],'omitnan');
                else
                    [AdhMaxAppMinArray(ff,1),AdhMaxAppMinFcArray(ff,1)]=min(ConcateArrayAdhMaxApp(((ff-1)*100)+1:ff*100),[],'omitnan');
                    [AdhMaxAppMaxArray(ff,1),AdhMaxAppMaxFcArray(ff,1)]=max(ConcateArrayAdhMaxApp(((ff-1)*100)+1:ff*100),[],'omitnan');
                    [AdhMaxRetMinArray(ff,1),AdhMaxRetMinFcArray(ff,1)]=min(ConcateArrayAdhMaxRet(((ff-1)*100)+1:ff*100),[],'omitnan');
                    [AdhMaxRetMaxArray(ff,1),AdhMaxRetMaxFcArray(ff,1)]=max(ConcateArrayAdhMaxRet(((ff-1)*100)+1:ff*100),[],'omitnan');
                    [AdhMaxRetUnbindingMinArray(ff,1),AdhMaxRetUnbindingMinFcArray(ff,1)]=min(ConcateArrayAdhMaxRetUnbinding(((ff-1)*100)+1:ff*100),[],'omitnan');
                    [AdhMaxRetUnbindingMaxArray(ff,1),AdhMaxRetUnbindingMaxFcArray(ff,1)]=max(ConcateArrayAdhMaxRetUnbinding(((ff-1)*100)+1:ff*100),[],'omitnan');
                    [AdhEneAppMinArray(ff,1),AdhEneAppMinFcArray(ff,1)]=min(ConcateArrayAdhEneApp(((ff-1)*100)+1:ff*100),[],'omitnan');
                    [AdhEneAppMaxArray(ff,1),AdhEneAppMaxFcArray(ff,1)]=max(ConcateArrayAdhEneApp(((ff-1)*100)+1:ff*100),[],'omitnan');
                    [AdhEneRetMinArray(ff,1),AdhEneRetMinFcArray(ff,1)]=min(ConcateArrayAdhEneRet(((ff-1)*100)+1:ff*100),[],'omitnan');
                    [AdhEneRetMaxArray(ff,1),AdhEneRetMaxFcArray(ff,1)]=max(ConcateArrayAdhEneRet(((ff-1)*100)+1:ff*100),[],'omitnan');
                    [PullLengthMinArray(ff,1),PullLengthMinFcArray(ff,1)]=min(ConcateArrayPullLength(((ff-1)*100)+1:ff*100),[],'omitnan');
                    [PullLengthMaxArray(ff,1),PullLengthMaxFcArray(ff,1)]=max(ConcateArrayPullLength(((ff-1)*100)+1:ff*100),[],'omitnan');
                    [SnapInMinArray(ff,1),SnapInMinFcArray(ff,1)]=min(ConcateArraySnapIn(((ff-1)*100)+1:ff*100),[],'omitnan');
                    [SnapInMaxArray(ff,1),SnapInMaxFcArray(ff,1)]=max(ConcateArraySnapIn(((ff-1)*100)+1:ff*100),[],'omitnan');
                end
            end
            % Sorting the force maps and force curves based on measurement time
            [~, FMSortDateTimeIdxFc]=sort(FMDateTimeNumber); % Sort force maps based on the time of measurement
            FMIndexChrono=FMIdxArray(FMSortDateTimeIdxFc); % Order the force maps chronologically
            [FMSortDateTimeFcArray, FMSortDateTimeIdxFcArray]=sort(FMDateTimeNumberFcArray);
            FcNumArrayChrono=FcNumArray(FMSortDateTimeIdxFcArray);
            FMIndexChronoArray=FMIndexArray(FMSortDateTimeIdxFcArray);
            % Statistics
            AdhMaxAppMean=mean(ConcateArrayAdhMaxApp,'omitnan');
            AdhMaxAppMedian=median(ConcateArrayAdhMaxApp,'omitnan');
            AdhMaxAppStd=std(ConcateArrayAdhMaxApp,'omitnan');
            [AdhMaxAppMin,AdhMaxAppMinIdx]=min(AdhMaxAppMinArray);
            [AdhMaxAppMax,AdhMaxAppMaxIdx]=max(AdhMaxAppMaxArray);
            AdhMaxAppMinFM=FMIdxArray(AdhMaxAppMinIdx);
            AdhMaxAppMaxFM=FMIdxArray(AdhMaxAppMaxIdx);
            AdhMaxAppMinFc=AdhMaxAppMinFcArray(AdhMaxAppMinIdx,1);
            AdhMaxAppMaxFc=AdhMaxAppMaxFcArray(AdhMaxAppMaxIdx,1);
            AdhMaxRetMean=mean(ConcateArrayAdhMaxRet,'omitnan');
            AdhMaxRetMedian=median(ConcateArrayAdhMaxRet,'omitnan');
            AdhMaxRetStd=std(ConcateArrayAdhMaxRet,'omitnan');
            [AdhMaxRetMin,AdhMaxRetMinIdx]=min(AdhMaxRetMinArray);
            [AdhMaxRetMax,AdhMaxRetMaxIdx]=max(AdhMaxRetMaxArray);
            AdhMaxRetMinFM=FMIdxArray(AdhMaxRetMinIdx);
            AdhMaxRetMaxFM=FMIdxArray(AdhMaxRetMaxIdx);
            AdhMaxRetMinFc=AdhMaxRetMinFcArray(AdhMaxRetMinIdx,1);
            AdhMaxRetMaxFc=AdhMaxRetMaxFcArray(AdhMaxRetMaxIdx,1);
            AdhMaxRetUnbindingMean=mean(ConcateArrayAdhMaxRetUnbinding,'omitnan');
            AdhMaxRetUnbindingMedian=median(ConcateArrayAdhMaxRetUnbinding,'omitnan');
            AdhMaxRetUnbindingStd=std(ConcateArrayAdhMaxRetUnbinding,'omitnan');
            [AdhMaxRetUnbindingMin,AdhMaxRetUnbindingMinIdx]=min(AdhMaxRetUnbindingMinArray);
            [AdhMaxRetUnbindingMax,AdhMaxRetUnbindingMaxIdx]=max(AdhMaxRetUnbindingMaxArray);
            AdhMaxRetUnbindingMinFM=FMIdxArray(AdhMaxRetUnbindingMinIdx);
            AdhMaxRetUnbindingMaxFM=FMIdxArray(AdhMaxRetUnbindingMaxIdx);
            AdhMaxRetUnbindingMinFc=AdhMaxRetUnbindingMinFcArray(AdhMaxRetUnbindingMinIdx,1);
            AdhMaxRetUnbindingMaxFc=AdhMaxRetUnbindingMaxFcArray(AdhMaxRetUnbindingMaxIdx,1);
            AdhEneAppMean=mean(ConcateArrayAdhEneApp,'omitnan');
            AdhEneAppMedian=median(ConcateArrayAdhEneApp,'omitnan');
            AdhEneAppStd=std(ConcateArrayAdhEneApp,'omitnan');
            [AdhEneAppMin,AdhEneAppMinIdx]=min(AdhEneAppMinArray);
            [AdhEneAppMax,AdhEneAppMaxIdx]=max(AdhEneAppMaxArray);
            AdhEneAppMinFM=FMIdxArray(AdhEneAppMinIdx);
            AdhEneAppMaxFM=FMIdxArray(AdhEneAppMaxIdx);
            AdhEneAppMinFc=AdhEneAppMinFcArray(AdhEneAppMinIdx,1);
            AdhEneAppMaxFc=AdhEneAppMaxFcArray(AdhEneAppMaxIdx,1);
            AdhEneRetMean=mean(ConcateArrayAdhEneRet,'omitnan');
            AdhEneRetMedian=median(ConcateArrayAdhEneRet,'omitnan');
            AdhEneRetStd=std(ConcateArrayAdhEneRet,'omitnan');
            [AdhEneRetMin,AdhEneRetMinIdx]=min(AdhEneRetMinArray);
            [AdhEneRetMax,AdhEneRetMaxIdx]=max(AdhEneRetMaxArray);
            AdhEneRetMinFM=FMIdxArray(AdhEneRetMinIdx);
            AdhEneRetMaxFM=FMIdxArray(AdhEneRetMaxIdx);
            AdhEneRetMinFc=AdhEneRetMinFcArray(AdhEneRetMinIdx,1);
            AdhEneRetMaxFc=AdhEneRetMaxFcArray(AdhEneRetMaxIdx,1);
            PullLengthMean=mean(ConcateArrayPullLength,'omitnan');
            PullLengthMedian=median(ConcateArrayPullLength,'omitnan');
            PullLengthStd=std(ConcateArrayPullLength,'omitnan');
            [PullLengthMin,PullLengthMinIdx]=min(PullLengthMinArray);
            [PullLengthMax,PullLengthMaxIdx]=max(PullLengthMaxArray);
            PullLengthMinFM=FMIdxArray(PullLengthMinIdx);
            PullLengthMaxFM=FMIdxArray(PullLengthMaxIdx);
            PullLengthMinFc=PullLengthMinFcArray(PullLengthMinIdx,1);
            PullLengthMaxFc=PullLengthMaxFcArray(PullLengthMaxIdx,1);
            SnapInMean=mean(ConcateArraySnapIn,'omitnan');
            SnapInMedian=median(ConcateArraySnapIn,'omitnan');
            SnapInStd=std(ConcateArraySnapIn,'omitnan');
            [SnapInMin,SnapInMinIdx]=min(SnapInMinArray);
            [SnapInMax,SnapInMaxIdx]=max(SnapInMaxArray);
            SnapInMinFM=FMIdxArray(SnapInMinIdx);
            SnapInMaxFM=FMIdxArray(SnapInMaxIdx);
            SnapInMinFc=SnapInMinFcArray(PullLengthMinIdx,1);
            SnapInMaxFc=SnapInMaxFcArray(PullLengthMaxIdx,1);
            %% Data selection based on input parameters
            % Determine unique entries
            ExtVelocityValues=unique(FMExtVelocity)';
            RetVelocityValues=unique(FMRetVelocity)';
            HoldingTimeValues=unique(FMHoldingTime)';
            SubstrateValues=unique(FMSubstrate)';
            EnvCondValues=unique(FMEnvCond)';
            ChipCantValues=unique(FMChipCant)';
            ChipboxValues=unique(FMChipbox)';
            LinkerValues=unique(FMLinker)';
            % Preallocate
            ExtVelocityFMIdx=zeros(length(FMExtVelocity),length(ExtVelocityValues));
            ExtVelocityConcateIdx=zeros(length(FMIndexArray),length(ExtVelocityValues));
            RetVelocityFMIdx=zeros(length(FMRetVelocity),length(RetVelocityValues));
            RetVelocityConcateIdx=zeros(length(FMIndexArray),length(RetVelocityValues));
            HoldingTimeFMIdx=zeros(length(FMHoldingTime),length(HoldingTimeValues));
            HoldingTimeConcateIdx=zeros(length(FMIndexArray),length(HoldingTimeValues));
            SubstrateFMIdx=zeros(length(FMSubstrate),length(SubstrateValues));
            SubstrateConcateIdx=zeros(length(FMIndexArray),length(SubstrateValues));
            EnvCondFMIdx=zeros(length(FMEnvCond),length(EnvCondValues));
            EnvCondConcateIdx=zeros(length(FMIndexArray),length(EnvCondValues));
            ChipCantFMIdx=zeros(length(FMChipCant),length(ChipCantValues));
            ChipCantConcateIdx=zeros(length(FMIndexArray),length(ChipCantValues));
            ChipboxFMIdx=zeros(length(FMChipbox),length(ChipboxValues));
            ChipboxConcateIdx=zeros(length(FMIndexArray),length(ChipboxValues));
            LinkerFMIdx=zeros(length(FMLinker),length(LinkerValues));
            LinkerConcateIdx=zeros(length(FMIndexArray),length(LinkerValues));
            for Fm=1:length(ExtVelocityValues)
                [~, Loc] =ismember(FMExtVelocity,ExtVelocityValues(Fm)); % Identify the force maps that have the same value
                Loc=find(Loc); % Use find to get ride of nonzeros
                FMIdx=FMIdxArray(Loc); % Identify same entries in the concatenate data of the results table
                [~, ConcateLoc]=ismember(FMIndexArray,FMIdx); % Identify same entries in the concatenate data of the results table
                ConcateLoc=find(ConcateLoc); % Use find to get ride of nonzeros
                % Allocate the data
                ExtVelocityFMIdx(1:length(Loc),Fm)=FMIdx;
                ExtVelocityConcateIdx(1:length(ConcateLoc),Fm)=ConcateLoc;
            end
            for Fm=1:length(RetVelocityValues)
                [~, Loc] =ismember(FMRetVelocity,RetVelocityValues(Fm)); % Identify the force maps that have the same value
                Loc=find(Loc); % Use find to get ride of nonzeros
                FMIdx=FMIdxArray(Loc); % Identify same entries in the concatenate data of the results table
                [~, ConcateLoc]=ismember(FMIndexArray,FMIdx); % Identify same entries in the concatenate data of the results table
                ConcateLoc=find(ConcateLoc); % Use find to get ride of nonzeros
                % Allocate the data
                RetVelocityFMIdx(1:length(Loc),Fm)=FMIdx;
                RetVelocityConcateIdx(1:length(ConcateLoc),Fm)=ConcateLoc;
            end
            for Fm=1:length(HoldingTimeValues)
                [~, Loc] =ismember(FMHoldingTime,HoldingTimeValues(Fm)); % Identify the force maps that have the same value
                Loc=find(Loc); % Use find to get ride of nonzeros
                FMIdx=FMIdxArray(Loc); % Identify same entries in the concatenate data of the results table
                [~, ConcateLoc]=ismember(FMIndexArray,FMIdx); % Identify same entries in the concatenate data of the results table
                ConcateLoc=find(ConcateLoc); % Use find to get ride of nonzeros
                % Allocate the data
                HoldingTimeFMIdx(1:length(Loc),Fm)=FMIdx;
                HoldingTimeConcateIdx(1:length(ConcateLoc),Fm)=ConcateLoc;
            end
            for Fm=1:length(SubstrateValues)
                [~, Loc] =ismember(FMSubstrate,SubstrateValues(Fm)); % Identify the force maps that have the same value
                Loc=find(Loc); % Use find to get ride of nonzeros
                FMIdx=FMIdxArray(Loc); % Identify same entries in the concatenate data of the results table
                [~, ConcateLoc]=ismember(FMIndexArray,FMIdx); % Identify same entries in the concatenate data of the results table
                ConcateLoc=find(ConcateLoc); % Use find to get ride of nonzeros
                % Allocate the data
                SubstrateFMIdx(1:length(Loc),Fm)=FMIdx;
                SubstrateConcateIdx(1:length(ConcateLoc),Fm)=ConcateLoc;
            end
            for Fm=1:length(EnvCondValues)
                [~, Loc] =ismember(FMEnvCond,EnvCondValues(Fm)); % Identify the force maps that have the same value
                Loc=find(Loc); % Use find to get ride of nonzeros
                FMIdx=FMIdxArray(Loc); % Identify same entries in the concatenate data of the results table
                [~, ConcateLoc]=ismember(FMIndexArray,FMIdx); % Identify same entries in the concatenate data of the results table
                ConcateLoc=find(ConcateLoc); % Use find to get ride of nonzeros
                % Allocate the data
                EnvCondFMIdx(1:length(Loc),Fm)=FMIdx;
                EnvCondConcateIdx(1:length(ConcateLoc),Fm)=ConcateLoc;
            end
            for Fm=1:length(ChipCantValues)
                [~, Loc] =ismember(FMChipCant,ChipCantValues(Fm)); % Identify the force maps that have the same value
                Loc=find(Loc); % Use find to get ride of nonzeros
                FMIdx=FMIdxArray(Loc); % Identify same entries in the concatenate data of the results table
                [~, ConcateLoc]=ismember(FMIndexArray,FMIdx); % Identify same entries in the concatenate data of the results table
                ConcateLoc=find(ConcateLoc); % Use find to get ride of nonzeros
                % Allocate the data
                ChipCantFMIdx(1:length(Loc),Fm)=FMIdx;
                ChipCantConcateIdx(1:length(ConcateLoc),Fm)=ConcateLoc;
            end
            for Fm=1:length(ChipboxValues)
                [~, Loc] =ismember(FMChipbox,ChipboxValues(Fm)); % Identify the force maps that have the same value
                Loc=find(Loc); % Use find to get ride of nonzeros
                FMIdx=FMIdxArray(Loc); % Identify same entries in the concatenate data of the results table
                [~, ConcateLoc]=ismember(FMIndexArray,FMIdx); % Identify same entries in the concatenate data of the results table
                ConcateLoc=find(ConcateLoc); % Use find to get ride of nonzeros
                % Allocate the data
                ChipboxFMIdx(1:length(Loc),Fm)=FMIdx;
                ChipboxConcateIdx(1:length(ConcateLoc),Fm)=ConcateLoc;
            end
            for Fm=1:length(LinkerValues)
                [~, Loc] =ismember(FMLinker,LinkerValues(Fm)); % Identify the force maps that have the same value
                Loc=find(Loc); % Use find to get ride of nonzeros
                FMIdx=FMIdxArray(Loc); % Identify same entries in the concatenate data of the results table
                [~, ConcateLoc]=ismember(FMIndexArray,FMIdx); % Identify same entries in the concatenate data of the results table
                ConcateLoc=find(ConcateLoc); % Use find to get ride of nonzeros
                % Allocate the data
                LinkerFMIdx(1:length(Loc),Fm)=FMIdx;
                LinkerConcateIdx(1:length(ConcateLoc),Fm)=ConcateLoc;
            end
            %% SMFS Results structure
            % Check entry
            if ~isempty(obj.SMFSResults)
                jj=length(obj.SMFSResults)+1;
            else
                jj=1;
            end
            % Debugging
            % jj=11
            % Allocate data
            obj.SMFSResults{jj,1}.Concatenate(1).FMID=FMIDArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMIndex=FMIndexArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMNum=FMNumArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMIndexChrono=FMIndexChronoArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FcNum=FcNumArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FcNumChrono=FcNumArrayChrono;
            obj.SMFSResults{jj,1}.Concatenate(1).FMExtVelocity=FMExtVelocityArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMRetVelocity=FMRetVelocityArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMHoldingTime=FMHoldingTimeArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMSubstrate=FMSubstrateArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMEnvCond=FMEnvCondArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMChipCant=FMChipCantArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMChipbox=FMChipboxArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMLinker=FMLinkerArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMDateTime=FMDateTimeFcArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMDateTimeNumber=FMDateTimeNumberFcArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMDateTimeNumberSort=FMSortDateTimeFcArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMDateTimeNumberSortIdx=FMSortDateTimeIdxFcArray;
            obj.SMFSResults{jj,1}.Concatenate(1).AdhMaxApp=ConcateArrayAdhMaxApp;
            obj.SMFSResults{jj,1}.Concatenate(1).AdhMaxRet= ConcateArrayAdhMaxRet;
            obj.SMFSResults{jj,1}.Concatenate(1).AdhUnbinding=ConcateArrayAdhMaxRetUnbinding;
            obj.SMFSResults{jj,1}.Concatenate(1).AdhEneApp=ConcateArrayAdhEneApp;
            obj.SMFSResults{jj,1}.Concatenate(1).AdhEneRet=ConcateArrayAdhEneRet;
            obj.SMFSResults{jj,1}.Concatenate(1).PullingLength=ConcateArrayPullLength;
            obj.SMFSResults{jj,1}.Concatenate(1).SnapInLength=ConcateArraySnapIn;
            obj.SMFSResults{jj,1}.Flags(1).SMFSFlagUncorrupt=sum(SMFSFlagUncorrupt);
            obj.SMFSResults{jj,1}.Flags(1).UncorruptPct=(sum(SMFSFlagUncorrupt)/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Flags(1).SMFSFlagSelected=sum(SMFSFlagSelected);
            obj.SMFSResults{jj,1}.Flags(1).SelectedPct=(sum(SMFSFlagSelected)/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Flags(1).SMFSFlagAppMinCrit=sum(SMFSFlagAppMinCrit);
            obj.SMFSResults{jj,1}.Flags(1).AppMinCritPct=(sum(SMFSFlagAppMinCrit)/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Flags(1).SMFSFlagRetMinCrit=sum(SMFSFlagRetMinCrit);
            obj.SMFSResults{jj,1}.Flags(1).RetMinCritPct=(sum(SMFSFlagRetMinCrit)/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Flags(1).SMFSFlagLengthRequisite=sum(SMFSFlagLengthRequisite);
            obj.SMFSResults{jj,1}.Flags(1).LengthRequisitePct=(sum(SMFSFlagLengthRequisite)/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Flags(1).SMFSFlagFit=sum(SMFSFlagFit);
            obj.SMFSResults{jj,1}.Flags(1).FitPct=(sum(SMFSFlagFit)/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Flags(1).SMFSFlagFitLinear=sum(SMFSFlagFitLinear);
            obj.SMFSResults{jj,1}.Flags(1).FitLinearPct=(sum(SMFSFlagFitLinear)/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Flags(1).SMFSFlagFitSinoidal=sum(SMFSFlagFitSinoidal);
            obj.SMFSResults{jj,1}.Flags(1).FitSinoidalPct=(sum(SMFSFlagFitSinoidal)/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Flags(1).SMFSFlagSnapIn=sum(SMFSFlagSnapIn);
            obj.SMFSResults{jj,1}.Flags(1).SnapInPct=(sum(SMFSFlagSnapIn)/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Flags(1).SMFSFlagPullingLength=sum(SMFSFlagPullingLength);
            obj.SMFSResults{jj,1}.Flags(1).PullingLengthPct=(sum(SMFSFlagPullingLength)/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Data(1).FMIndex=FMIdxArray;
            obj.SMFSResults{jj,1}.Data(1).FMIndexChrono=FMIndexChrono;
            obj.SMFSResults{jj,1}.Data(1).FMID=FMID;
            obj.SMFSResults{jj,1}.Data(1).FMNum=FMNum;
            obj.SMFSResults{jj,1}.Data(1).FMExtVelocity=FMExtVelocity;
            obj.SMFSResults{jj,1}.Data(1).FMRetVelocity=FMRetVelocity;
            obj.SMFSResults{jj,1}.Data(1).FMHoldingTime=FMHoldingTime;
            obj.SMFSResults{jj,1}.Data(1).FMSubstrate=FMSubstrate;
            obj.SMFSResults{jj,1}.Data(1).FMEnvCond=FMEnvCond;
            obj.SMFSResults{jj,1}.Data(1).FMChipCant=FMChipCant;
            obj.SMFSResults{jj,1}.Data(1).FMChipbox=FMChipbox;
            obj.SMFSResults{jj,1}.Data(1).FMLinker=FMLinker;
            obj.SMFSResults{jj,1}.Data(1).FMDateTime=FMDateTime;
            obj.SMFSResults{jj,1}.Data(1).FMDateTimeSortIdx=FMSortDateTimeIdxFc;
            obj.SMFSResults{jj,1}.Data(1).TotalNumFc=sum(TotalNumFc);
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedAdhMaxApp=nnz(~isnan(ConcateArrayAdhMaxApp));
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedAdhMaxAppPct=(nnz(~isnan(ConcateArrayAdhMaxApp))/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedAdhMaxRet=nnz(~isnan(ConcateArrayAdhMaxRet));
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedAdhMaxRetPct=(nnz(~isnan(ConcateArrayAdhMaxRet))/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedAdhUnbinding=nnz(~isnan(ConcateArrayAdhMaxRetUnbinding));
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedAdhUnbindingPct=(nnz(~isnan(ConcateArrayAdhMaxRetUnbinding))/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedAdhEneApp=nnz(~isnan(ConcateArrayAdhEneApp));
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedAdhEneAppPct=(nnz(~isnan(ConcateArrayAdhEneApp))/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedAdhEneRet=nnz(~isnan(ConcateArrayAdhEneRet));
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedAdhEneRetPct=(nnz(~isnan(ConcateArrayAdhEneRet))/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedyPullingLength=nnz(~isnan(ConcateArrayPullLength));
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedPullingLengthPct=(nnz(~isnan(ConcateArrayPullLength))/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedySnapInLength=nnz(~isnan(ConcateArraySnapIn));
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedySnapInLengthPct=(nnz(~isnan(ConcateArraySnapIn))/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Data(1).NumFcAnalysedAdhMaxApp=FCAnalysedyAdhMaxApp;
            obj.SMFSResults{jj,1}.Data(1).NumFcAnalysedAdhMaxRet=FCAnalysedyAdhMaxRet;
            obj.SMFSResults{jj,1}.Data(1).NumFcAnalysedAdhUnbinding=FCAnalysedyAdhUnbind;
            obj.SMFSResults{jj,1}.Data(1).NumFcAnalysedAdhEneApp=FCAnalysedyAdhEneApp;
            obj.SMFSResults{jj,1}.Data(1).NumFcAnalysedAdhEneRet=FCAnalysedyAdhEneRet;
            obj.SMFSResults{jj,1}.Data(1).NumFcAnalysedyPullingLength=FCAnalysedyPullingLength;
            obj.SMFSResults{jj,1}.Data(1).NumFcAnalysedySnapInLength=FCAnalysedySnapInLength;
            obj.SMFSResults{jj,1}.Data(1).AdhMaxApp=yAdhMaxAppAll;
            obj.SMFSResults{jj,1}.Data(1).AdhMaxRet=yAdhMaxRetAll;
            obj.SMFSResults{jj,1}.Data(1).AdhUnbinding=yAdhUnbindingAll;
            obj.SMFSResults{jj,1}.Data(1).AdhEneApp=yAdhEneAppAll;
            obj.SMFSResults{jj,1}.Data(1).AdhEneRet=yAdhEneRetAll;
            obj.SMFSResults{jj,1}.Data(1).PullingLength=yPullingLengthAll;
            obj.SMFSResults{jj,1}.Data(1).SnapInLength=ySnapInLengthAll;
            obj.SMFSResults{jj,1}.Parameters(1).ExtendVelocity=ExtVelocityValue;
            obj.SMFSResults{jj,1}.Parameters(1).RetractVelocity=RetVelocityValue;
            obj.SMFSResults{jj,1}.Parameters(1).HoldingTime=HoldingTimeValue;
            obj.SMFSResults{jj,1}.Parameters(1).Substrate=SubstrateValue;
            obj.SMFSResults{jj,1}.Parameters(1).Medium=EnvCondValue;
            obj.SMFSResults{jj,1}.Parameters(1).ChipCantilever=ChipCantValue;
            obj.SMFSResults{jj,1}.Parameters(1).Chipbox=ChipboxValue;
            obj.SMFSResults{jj,1}.Parameters(1).Linker=LinkerValue;
            obj.SMFSResults{jj,1}.Selection(1).ExtVelocityParameters=ExtVelocityValues;
            obj.SMFSResults{jj,1}.Selection(1).ExtVelocityFMIdx=ExtVelocityFMIdx;
            obj.SMFSResults{jj,1}.Selection(1).ExtVelocityConcateIdx=ExtVelocityConcateIdx;
            obj.SMFSResults{jj,1}.Selection(1).RetVelocityParameters=RetVelocityValues;
            obj.SMFSResults{jj,1}.Selection(1).RetVelocityFMIdx=RetVelocityFMIdx;
            obj.SMFSResults{jj,1}.Selection(1).RetVelocityConcateIdx=RetVelocityConcateIdx;
            obj.SMFSResults{jj,1}.Selection(1).HoldingTimeParameters=HoldingTimeValues;
            obj.SMFSResults{jj,1}.Selection(1).HoldingTimeFMIdx=HoldingTimeFMIdx;
            obj.SMFSResults{jj,1}.Selection(1).HoldingTimeConcateIdx=HoldingTimeConcateIdx;
            obj.SMFSResults{jj,1}.Selection(1).SubstrateParameters=SubstrateValues;
            obj.SMFSResults{jj,1}.Selection(1).SubstrateFMIdx=SubstrateFMIdx;
            obj.SMFSResults{jj,1}.Selection(1).SubstrateConcateIdx=SubstrateConcateIdx;
            obj.SMFSResults{jj,1}.Selection(1).EnvCondParameters=EnvCondValues;
            obj.SMFSResults{jj,1}.Selection(1).EnvCondFMIdx=EnvCondFMIdx;
            obj.SMFSResults{jj,1}.Selection(1).EnvCondConcateIdx=EnvCondConcateIdx;
            obj.SMFSResults{jj,1}.Selection(1).ChipCantParameters=ChipCantValues;
            obj.SMFSResults{jj,1}.Selection(1).ChipCantFMIdx=ChipCantFMIdx;
            obj.SMFSResults{jj,1}.Selection(1).ChipCantConcateIdx=ChipCantConcateIdx;
            obj.SMFSResults{jj,1}.Selection(1).ChipboxParameters=ChipboxValues;
            obj.SMFSResults{jj,1}.Selection(1).ChipboxFMIdx=ChipboxFMIdx;
            obj.SMFSResults{jj,1}.Selection(1).ChipboxConcateIdx=ChipboxConcateIdx;
            obj.SMFSResults{jj,1}.Selection(1).LinkerParameters=LinkerValues;
            obj.SMFSResults{jj,1}.Selection(1).LinkerFMIdx=LinkerFMIdx;
            obj.SMFSResults{jj,1}.Selection(1).LinkerConcateIdx=LinkerConcateIdx;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppMean=AdhMaxAppMean;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppMedian=AdhMaxAppMedian;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppStd=AdhMaxAppStd;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppMinArray=AdhMaxAppMinArray;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppMin=AdhMaxAppMin;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppMaxArray=AdhMaxAppMaxArray;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppMax=AdhMaxAppMax;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppMinFM=AdhMaxAppMinFM;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppMaxFM=AdhMaxAppMaxFM;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppMinFcArray=AdhMaxAppMinFcArray;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppMinFc=AdhMaxAppMinFc;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppMaxFcArray=AdhMaxAppMaxFcArray;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppMaxFc=AdhMaxAppMaxFc;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetMean=AdhMaxRetMean;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetMedian=AdhMaxRetMedian;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetStd=AdhMaxRetStd;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetMinArray=AdhMaxRetMinArray;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetMin=AdhMaxRetMin;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetMaxArray=AdhMaxRetMaxArray;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetMax=AdhMaxRetMax;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetMinFM=AdhMaxRetMinFM;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetMaxFM=AdhMaxRetMaxFM;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetMinFcArray=AdhMaxRetMinFcArray;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetMinFc=AdhMaxRetMinFc;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetMaxFcArray=AdhMaxRetMaxFcArray;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetMaxFc=AdhMaxRetMaxFc;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingMean=AdhMaxRetUnbindingMean;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingMedian=AdhMaxRetUnbindingMedian;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingStd=AdhMaxRetUnbindingStd;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingMinArray=AdhMaxRetUnbindingMinArray;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingMin=AdhMaxRetUnbindingMin;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingMaxArray=AdhMaxRetUnbindingMaxArray;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingMax=AdhMaxRetUnbindingMax;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingMinFM=AdhMaxRetUnbindingMinFM;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingMaxFM=AdhMaxRetUnbindingMaxFM;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingMinFcArray=AdhMaxRetUnbindingMinFcArray;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingMinFc=AdhMaxRetUnbindingMinFc;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingMaxFcArray=AdhMaxRetUnbindingMaxFcArray;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingMaxFc=AdhMaxRetUnbindingMaxFc;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppMean=AdhEneAppMean;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppMedian=AdhEneAppMedian;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppStd=AdhEneAppStd;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppMinArray=AdhEneAppMinArray;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppMin=AdhEneAppMin;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppMaxArray=AdhEneAppMaxArray;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppMax=AdhEneAppMax;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppMinFM=AdhEneAppMinFM;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppMaxFM=AdhEneAppMaxFM;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppMinFcArray=AdhEneAppMinFcArray;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppMinFc=AdhEneAppMinFc;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppMaxFcArray=AdhEneAppMaxFcArray;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppMaxFc=AdhEneAppMaxFc;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetMean=AdhEneRetMean;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetMedian=AdhEneRetMedian;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetStd=AdhEneRetStd;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetMinArray=AdhEneRetMinArray;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetMin=AdhEneRetMin;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetMaxArray=AdhEneRetMaxArray;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetMax=AdhEneRetMax;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetMinFM=AdhEneRetMinFM;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetMaxFM=AdhEneRetMaxFM;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetMinFcArray=AdhEneRetMinFcArray;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetMinFc=AdhEneRetMinFc;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetMaxFcArray=AdhEneRetMaxFcArray;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetMaxFc=AdhEneRetMaxFc;
            obj.SMFSResults{jj,1}.Results(1).PullLengthMean=PullLengthMean;
            obj.SMFSResults{jj,1}.Results(1).PullLengthMedian=PullLengthMedian;
            obj.SMFSResults{jj,1}.Results(1).PullLengthStd=PullLengthStd;
            obj.SMFSResults{jj,1}.Results(1).PullLengthMinArray=PullLengthMinArray;
            obj.SMFSResults{jj,1}.Results(1).PullLengthMin=PullLengthMin;
            obj.SMFSResults{jj,1}.Results(1).PullLengthMinFM=PullLengthMinFM;
            obj.SMFSResults{jj,1}.Results(1).PullLengthMinFcArray=PullLengthMinFcArray;
            obj.SMFSResults{jj,1}.Results(1).PullLengthMinFc=PullLengthMinFc;
            obj.SMFSResults{jj,1}.Results(1).PullLengthMaxArray=PullLengthMaxArray;
            obj.SMFSResults{jj,1}.Results(1).PullLengthMax=PullLengthMax;
            obj.SMFSResults{jj,1}.Results(1).PullLengthMaxFM=PullLengthMaxFM;
            obj.SMFSResults{jj,1}.Results(1).PullLengthMaxFcArray=PullLengthMaxFcArray;
            obj.SMFSResults{jj,1}.Results(1).PullLengthMaxFc=PullLengthMaxFc;
            obj.SMFSResults{jj,1}.Results(1).SnapInMean=SnapInMean;
            obj.SMFSResults{jj,1}.Results(1).SnapInMedian=SnapInMedian;
            obj.SMFSResults{jj,1}.Results(1).SnapInStd=SnapInStd;
            obj.SMFSResults{jj,1}.Results(1).SnapInMinArray=SnapInMinArray;
            obj.SMFSResults{jj,1}.Results(1).SnapInMin=SnapInMin;
            obj.SMFSResults{jj,1}.Results(1).SnapInMinFM=SnapInMinFM;
            obj.SMFSResults{jj,1}.Results(1).SnapInMinFcArray=SnapInMinFcArray;
            obj.SMFSResults{jj,1}.Results(1).SnapInMinFc=SnapInMinFc;
            obj.SMFSResults{jj,1}.Results(1).SnapInMaxArray=SnapInMaxArray;
            obj.SMFSResults{jj,1}.Results(1).SnapInMax=SnapInMax;
            obj.SMFSResults{jj,1}.Results(1).SnapInMaxFM=SnapInMaxFM;
            obj.SMFSResults{jj,1}.Results(1).SnapInMaxFcArray=SnapInMaxFcArray;
            obj.SMFSResults{jj,1}.Results(1).SnapInMaxFc=SnapInMaxFc;
            % Paramater list
            obj.SMFSResultsParameters(jj,:)={jj,ChipboxValue,ChipCantValue,LinkerValue,SubstrateValue,EnvCondValue,ExtVelocityValue,RetVelocityValue,HoldingTimeValue};
 
        end

        function SMFS_results_structure_FM(obj,Var)
            % Author: Andreas Rohatschek
            % Description: 
            % In this function a results structure entry is created based on the
            % input variable (Variant).
            % Required input variables:
            % Var ... Variant: double, either 1 or 2
            % Typical input variables:
            % Var ... 1

            %% Function body
            % Define variables
            DateFormat='yyyy-MM-dd HH-mm-ss-SSS';
            ExtVelocityValue=0;
            RetVelocityValue=0;
            HoldingTimeValue=-1;
            SubstrateValue='Various';
            EnvCondValue='Various';
            ChipCantValue='Various';
            ChipboxValue='Various';
            LinkerValue='Various';
            
            if Var==1
            % Input dialog
            prompt = {'Enter the force map number you do not want to have included in the "SMFSResults"-structure (For multiple selections just use the space key to separeat entries)'};
            definput = {''};
            opts.Interpreter = 'tex';
            FMIdxArray=inputdlg(prompt,'Select all - except of ...',[1 150],definput,opts);
            FMIdxArray=str2num(FMIdxArray{1}); % Convert the cell array to numerals
            elseif Var==2
            warning('This variation requires that a none parameter based selection row in the "SMFSResults"-structure exists (Typically row 1)')
            dlgtitle='Enter the "SMFSResults"-structure row the chronologcial data will be taken from';
            prompt = {'ResultsRow'};
            definput={'1'};
            dims=[1 150];
            ResultsRow = inputdlg(prompt,dlgtitle,dims,definput); 
            ResultsRow=str2num(ResultsRow{1});
            dlgtitle='Enter the chronological force map index of the first and last force map you want to have included in the "SMFSResults"-structure';
            prompt = {'First','Last'};
            dims=[1 150; 1 150];
            FMIdxChrono = inputdlg(prompt,dlgtitle,dims); 
            FMIdxChrono1=str2num(FMIdxChrono{1});
            FMIdxChrono2=str2num(FMIdxChrono{2});           
            % Allocate the force maps in chronological order as defined by
            % the input dialog box
            FMIdxArray=obj.SMFSResults{ResultsRow}.Data.FMDateTimeSortIdx(FMIdxChrono1:FMIdxChrono2);
            end
            % If condition to handle an empty index array
            if isempty(FMIdxArray)
                return
            else
            end
            % Preallocate
            ConcateArrayAdhMaxApp=zeros(length(FMIdxArray)*obj.FM{FMIdxArray(1)}.NCurves,1);
            ConcateArrayAdhMaxRet=zeros(length(FMIdxArray)*obj.FM{FMIdxArray(1)}.NCurves,1);
            ConcateArrayAdhMaxRetUnbinding=zeros(length(FMIdxArray)*obj.FM{FMIdxArray(1)}.NCurves,1);
            ConcateArrayAdhEneApp=zeros(length(FMIdxArray)*obj.FM{FMIdxArray(1)}.NCurves,1);
            ConcateArrayAdhEneRet=zeros(length(FMIdxArray)*obj.FM{FMIdxArray(1)}.NCurves,1);
            ConcateArrayPullLength=zeros(length(FMIdxArray)*obj.FM{FMIdxArray(1)}.NCurves,1);
            ConcateArraySnapIn=zeros(length(FMIdxArray)*obj.FM{FMIdxArray(1)}.NCurves,1);
            yAdhMaxAppAll=zeros(obj.FM{FMIdxArray(1)}.NCurves,length(FMIdxArray));
            FCAnalysedyAdhMaxApp=zeros(1,length(FMIdxArray));
            yAdhMaxRetAll=zeros(obj.FM{FMIdxArray(1)}.NCurves,length(FMIdxArray));
            FCAnalysedyAdhMaxRet=zeros(1,length(FMIdxArray));
            yAdhUnbindingAll=zeros(obj.FM{FMIdxArray(1)}.NCurves,length(FMIdxArray));
            FCAnalysedyAdhUnbind=zeros(1,length(FMIdxArray));
            yAdhEneAppAll=zeros(obj.FM{FMIdxArray(1)}.NCurves,length(FMIdxArray));
            FCAnalysedyAdhEneApp=zeros(1,length(FMIdxArray));
            yAdhEneRetAll=zeros(obj.FM{FMIdxArray(1)}.NCurves,length(FMIdxArray));
            FCAnalysedyAdhEneRet=zeros(1,length(FMIdxArray));
            yPullingLengthAll=zeros(obj.FM{FMIdxArray(1)}.NCurves,length(FMIdxArray));
            FCAnalysedyPullingLength=zeros(1,length(FMIdxArray));
            ySnapInLengthAll=zeros(obj.FM{FMIdxArray(1)}.NCurves,length(FMIdxArray));
            FCAnalysedySnapInLength=zeros(1,length(FMIdxArray));            
            FMID=cell(length(FMIdxArray),1);
            FMNum=zeros(length(FMIdxArray),1);
            FMExtVelocity=zeros(length(FMIdxArray),1);
            FMRetVelocity=zeros(length(FMIdxArray),1);
            FMHoldingTime=zeros(length(FMIdxArray),1);
            FMSubstrate=cell(length(FMIdxArray),1);
            FMEnvCond=cell(length(FMIdxArray),1);
            FMChipCant=cell(length(FMIdxArray),1);
            FMChipbox=cell(length(FMIdxArray),1);
            FMLinker=cell(length(FMIdxArray),1);
            FMDateTime=cell(length(FMIdxArray),1);
            FMDateTimeNumber=zeros(length(FMIdxArray),1);
            TotalNumFc=zeros(length(FMIdxArray),1);            
            SMFSFlagUncorrupt=zeros(length(FMIdxArray),1);
            SMFSFlagSelected=zeros(length(FMIdxArray),1);
            SMFSFlagAppMinCrit=zeros(length(FMIdxArray),1);
            SMFSFlagRetMinCrit=zeros(length(FMIdxArray),1);
            SMFSFlagLengthRequisite=zeros(length(FMIdxArray),1);
            SMFSFlagFit=zeros(length(FMIdxArray),1);
            SMFSFlagFitLinear=zeros(length(FMIdxArray),1);
            SMFSFlagFitSinoidal=zeros(length(FMIdxArray),1);
            SMFSFlagSnapIn=zeros(length(FMIdxArray),1);
            SMFSFlagPullingLength=zeros(length(FMIdxArray),1);
            AdhMaxAppMinArray=zeros(length(FMIdxArray),1);
            AdhMaxAppMinFcArray=zeros(length(FMIdxArray),1);
            AdhMaxAppMaxArray=zeros(length(FMIdxArray),1);
            AdhMaxAppMaxFcArray=zeros(length(FMIdxArray),1);
            AdhMaxRetMinArray=zeros(length(FMIdxArray),1);
            AdhMaxRetMinFcArray=zeros(length(FMIdxArray),1);
            AdhMaxRetMaxArray=zeros(length(FMIdxArray),1);
            AdhMaxRetMaxFcArray=zeros(length(FMIdxArray),1);
            AdhMaxRetUnbindingMinArray=zeros(length(FMIdxArray),1);
            AdhMaxRetUnbindingMinFcArray=zeros(length(FMIdxArray),1);
            AdhMaxRetUnbindingMaxArray=zeros(length(FMIdxArray),1);
            AdhMaxRetUnbindingMaxFcArray=zeros(length(FMIdxArray),1);
            AdhEneAppMinArray=zeros(length(FMIdxArray),1);
            AdhEneAppMinFcArray=zeros(length(FMIdxArray),1);
            AdhEneAppMaxArray=zeros(length(FMIdxArray),1);
            AdhEneAppMaxFcArray=zeros(length(FMIdxArray),1);
            AdhEneRetMinArray=zeros(length(FMIdxArray),1);
            AdhEneRetMinFcArray=zeros(length(FMIdxArray),1);
            AdhEneRetMaxArray=zeros(length(FMIdxArray),1);
            AdhEneRetMaxFcArray=zeros(length(FMIdxArray),1);
            PullLengthMinArray=zeros(length(FMIdxArray),1);
            PullLengthMinFcArray=zeros(length(FMIdxArray),1);
            PullLengthMaxArray=zeros(length(FMIdxArray),1);
            PullLengthMaxFcArray=zeros(length(FMIdxArray),1);
            SnapInMinArray=zeros(length(FMIdxArray),1);
            SnapInMinFcArray=zeros(length(FMIdxArray),1);
            SnapInMaxArray=zeros(length(FMIdxArray),1);
            SnapInMaxFcArray=zeros(length(FMIdxArray),1);
            % Loop
            for ff=1:length(FMIdxArray)
                %% Debugging
                %for ff=6 % for debugging
                %sprintf('Index array row No. %d',ff) % Gives current Force curve
                % Allocate data
                yAdhMaxApp=obj.FM{FMIdxArray(ff)}.AdhForceMaxApp;
                yAdhMaxRet=obj.FM{FMIdxArray(ff)}.AdhForceMaxRet;
                yAdhUnbinding=obj.FM{FMIdxArray(ff)}.AdhForceUnbinding;
                yAdhEneApp=obj.FM{FMIdxArray(ff)}.AppAdhEnergy_IdxMethod;
                yAdhEneRet=obj.FM{FMIdxArray(ff)}.RetAdhEnergy_IdxMethod;
                yPullingLength=obj.FM{FMIdxArray(ff)}.PullingLength;
                ySnapInLength=obj.FM{FMIdxArray(ff)}.SnapInLength;
                FMID{ff,1}=obj.FM{FMIdxArray(ff)}.ID;
                FMNum(ff,1)=ff;
                FMExtVelocity(ff,1)=round(obj.FM{FMIdxArray(ff)}.ExtendVelocity,8); % Round holding time value to correct for tiny deviations of decimal places after the comma origin from instrument (AFM)
                FMRetVelocity(ff,1)=round(obj.FM{FMIdxArray(ff)}.RetractVelocity,8); % Round holding time value to correct for tiny deviations of decimal places after the comma origin from instrument (AFM)
                FMHoldingTime(ff,1)=round(obj.FM{FMIdxArray(ff)}.HoldingTime,2); % Round holding time value to correct for tiny deviations of decimal places after the comma origin from instrument (AFM)
                FMSubstrate{ff,1}=obj.FM{FMIdxArray(ff)}.Substrate;
                FMEnvCond{ff,1}=obj.FM{FMIdxArray(ff)}.EnvCond;
                FMChipCant{ff,1}=obj.FM{FMIdxArray(ff)}.ChipCant;
                FMChipbox{ff,1}=obj.FM{FMIdxArray(ff)}.Chipbox;
                FMLinker{ff,1}=obj.FM{FMIdxArray(ff)}.Linker;
                FMDate=obj.FM{FMIdxArray(ff)}.Date;
                FMTime=obj.FM{FMIdxArray(ff)}.Time;
                DateTimeStr=[FMDate,' ',FMTime];
                FMDateTime{ff,1}=cellstr(datetime(DateTimeStr,'InputFormat',DateFormat,'Format',DateFormat));
                FMDateTimeNumber(ff,1)=datenum(datetime(DateTimeStr,'InputFormat',DateFormat,'Format',DateFormat));
                TotalNumFc(ff,1)=obj.FM{FMIdxArray(ff)}.NCurves;
                % SMFS Flags
                SMFSFlagUncorrupt(ff,1)=nnz(obj.FM{FMIdxArray(ff)}.SMFSFlag.Uncorrupt');
                SMFSFlagSelected(ff,1)=nnz(obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected');
                SMFSFlagAppMinCrit(ff,1)=nnz(obj.FM{FMIdxArray(ff)}.SMFSFlag.AppMinCrit');
                SMFSFlagRetMinCrit(ff,1)=nnz(obj.FM{FMIdxArray(ff)}.SMFSFlag.RetMinCrit');
                SMFSFlagLengthRequisite(ff,1)=nnz(obj.FM{FMIdxArray(ff)}.SMFSFlag.LengthRequisite');
                SMFSFlagFit(ff,1)=nnz(obj.FM{FMIdxArray(ff)}.SMFSFlag.Fit');
                SMFSFlagFitLinear(ff,1)=nnz(obj.FM{FMIdxArray(ff)}.SMFSFlag.FitLinear');
                SMFSFlagFitSinoidal(ff,1)=nnz(obj.FM{FMIdxArray(ff)}.SMFSFlag.FitSinoidal');
                SMFSFlagSnapIn(ff,1)=nnz(obj.FM{FMIdxArray(ff)}.SMFSFlag.SnapIn');
                SMFSFlagPullingLength(ff,1)=nnz(obj.FM{FMIdxArray(ff)}.SMFSFlag.PullingLength');
                %% Concatenate arrays
                % FCs of each FM in seperate column
                yAdhMaxAppAll(:,ff)=yAdhMaxApp'.*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected';
                yAdhMaxAppAll(yAdhMaxAppAll==0)=nan; % Replace zero entries by nan´s
                FCAnalysedyAdhMaxApp(1,ff)=nnz(~isnan(yAdhMaxAppAll(:,ff)));
                yAdhMaxRetAll(:,ff)=yAdhMaxRet'.*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected';
                yAdhMaxRetAll(yAdhMaxRetAll==0)=nan; % Replace zero entries by nan´s
                FCAnalysedyAdhMaxRet(1,ff)=nnz(~isnan(yAdhMaxRetAll(:,ff)));
                yAdhUnbindingAll(:,ff)=yAdhUnbinding'.*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected';
                yAdhUnbindingAll(yAdhUnbindingAll==0)=nan; % Replace zero entries by nan´s
                FCAnalysedyAdhUnbind(1,ff)=nnz(~isnan(yAdhUnbindingAll(:,ff)));
                yAdhEneAppAll(:,ff)=yAdhEneApp'.*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected';
                yAdhEneAppAll(yAdhEneAppAll==0)=nan; % Replace zero entries by nan´s
                FCAnalysedyAdhEneApp(1,ff)=nnz(~isnan(yAdhEneAppAll(:,ff)));
                yAdhEneRetAll(:,ff)=yAdhEneRet'.*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected';
                yAdhEneRetAll(yAdhEneRetAll==0)=nan; % Replace zero entries by nan´s
                FCAnalysedyAdhEneRet(1,ff)=nnz(~isnan(yAdhEneRetAll(:,ff)));
                yPullingLengthAll(:,ff)=yPullingLength'.*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected';
                yPullingLengthAll(yPullingLengthAll==0)=nan; % Replace zero entries by nan´s
                FCAnalysedyPullingLength(1,ff)=nnz(~isnan(yPullingLengthAll(:,ff)));
                ySnapInLengthAll(:,ff)=ySnapInLength'.*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected';
                ySnapInLengthAll(ySnapInLengthAll==0)=nan; % Replace zero entries by nan´s
                FCAnalysedySnapInLength(1,ff)=nnz(~isnan(ySnapInLengthAll(:,ff)));
                % All FCs of all FM in one column                   
                if ~isempty(yAdhMaxApp)
                    % Determine the number of rows per force map
                    ArrayLength=length(yAdhMaxApp); % Define the length of the array
                    row_start = ((ff-1) * ArrayLength) + 1; % Define the appropriate row start to append the new data
                    row_end   = ff * ArrayLength; % Define the appropriate row end to append the new data
                    % Concatenated data
                    ConcateArrayAdhMaxApp(row_start:row_end,:)=yAdhMaxApp'; % Append the new data into the concatenated vector
                    ConcateArrayAdhMaxApp(row_start:row_end,:)=ConcateArrayAdhMaxApp(row_start:row_end,:).*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected'; % Set non-selected force curves from the concatenated arrays to zero
                    ConcateArrayAdhMaxApp(ConcateArrayAdhMaxApp==0)=nan; % Replace zero entries by nan´s
                    % Allocate parameters
                    FMIDArray(row_start:row_end,:)={FMID}; % Allocate the FM ID to each row
                    FMIndexArray(row_start:row_end,:)=FMIdxArray(ff);
                    FMNumArray(row_start:row_end,:)=ff;
                    FcNumArray(row_start:row_end,:)=(row_start:row_end)';
                    FMExtVelocityArray(row_start:row_end,:)=FMExtVelocity(ff,1);
                    FMRetVelocityArray(row_start:row_end,:)=FMRetVelocity(ff,1);
                    FMHoldingTimeArray(row_start:row_end,:)=FMHoldingTime(ff,1);
                    FMSubstrateArray(row_start:row_end,:)=FMSubstrate(ff,1);
                    FMEnvCondArray(row_start:row_end,:)=FMEnvCond(ff,1);
                    FMChipCantArray(row_start:row_end,:)=FMChipCant(ff,1);
                    FMChipboxArray(row_start:row_end,:)=FMChipbox(ff,1);
                    FMLinkerArray(row_start:row_end,:)=FMLinker(ff,1);
                    FMDateTimeFcArray(row_start:row_end,:)=datetime(DateTimeStr,'InputFormat',DateFormat,'Format',DateFormat);
                    FMDateTimeNumberFcArray(row_start:row_end,:)=datenum(datetime(DateTimeStr,'InputFormat',DateFormat,'Format',DateFormat));
                else
                end
                if ~isempty(yAdhMaxRet)
                    % Determine the number of rows per force map
                    ArrayLength=length(yAdhMaxRet); % Define the length of the array
                    row_start = ((ff-1) * ArrayLength) + 1; % Define the appropriate row start to append the new data
                    row_end   = ff * ArrayLength; % Define the appropriate row end to append the new data
                    % Concatenated data
                    ConcateArrayAdhMaxRet(row_start:row_end,:)=yAdhMaxRet'; % Append the new data into the concatenated vector
                    ConcateArrayAdhMaxRet(row_start:row_end,:)=ConcateArrayAdhMaxRet(row_start:row_end,:).*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected'; % Set non-selected force curves from the concatenated arrays to zero
                    ConcateArrayAdhMaxRet(ConcateArrayAdhMaxRet==0)=nan; % Replace zero entries by nan´s
                else
                end
                if ~isempty(yAdhUnbinding)
                    % Determine the number of rows per force map
                    ArrayLength=length(yAdhUnbinding); % Define the length of the array
                    row_start = ((ff-1) * ArrayLength) + 1; % Define the appropriate row start to append the new data
                    row_end   = ff * ArrayLength; % Define the appropriate row end to append the new data
                    % Concatenated data
                    ConcateArrayAdhMaxRetUnbinding(row_start:row_end,:)=yAdhUnbinding'; % Append the new data into the concatenated vector
                    ConcateArrayAdhMaxRetUnbinding(row_start:row_end,:)=ConcateArrayAdhMaxRetUnbinding(row_start:row_end,:).*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected'; % Set non-selected force curves from the concatenated arrays to zero
                    ConcateArrayAdhMaxRetUnbinding(ConcateArrayAdhMaxRetUnbinding==0)=nan; % Replace zero entries by nan´s
                else
                end
                if ~isempty(yAdhEneApp)
                    % Determine the number of rows per force map
                    ArrayLength=length(yAdhEneApp); % Define the length of the array
                    row_start = ((ff-1) * ArrayLength) + 1; % Define the appropriate row start to append the new data
                    row_end   = ff * ArrayLength; % Define the appropriate row end to append the new data
                    % Concatenated data
                    ConcateArrayAdhEneApp(row_start:row_end,:)=yAdhEneApp'; % Append the new data into the concatenated vector
                    ConcateArrayAdhEneApp(row_start:row_end,:)=ConcateArrayAdhEneApp(row_start:row_end,:).*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected'; % Set non-selected force curves from the concatenated arrays to zero
                    ConcateArrayAdhEneApp(ConcateArrayAdhEneApp==0)=nan; % Replace zero entries by nan´s
                else
                end
                if ~isempty(yAdhEneRet)
                    % Determine the number of rows per force map
                    ArrayLength=length(yAdhEneRet); % Define the length of the array
                    row_start = ((ff-1) * ArrayLength) + 1; % Define the appropriate row start to append the new data
                    row_end   = ff * ArrayLength; % Define the appropriate row end to append the new data
                    % Concatenated data
                    ConcateArrayAdhEneRet(row_start:row_end,:)=yAdhEneRet'; % Append the new data into the concatenated vector
                    ConcateArrayAdhEneRet(row_start:row_end,:)=ConcateArrayAdhEneRet(row_start:row_end,:).*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected'; % Set non-selected force curves from the concatenated arrays to zero
                    ConcateArrayAdhEneRet(ConcateArrayAdhEneRet==0)=nan; % Replace zero entries by nan´s
                else
                end
                if ~isempty(yPullingLength)
                    % Determine the number of rows per force map
                    ArrayLength=length(yPullingLength); % Define the length of the array
                    row_start = ((ff-1) * ArrayLength) + 1; % Define the appropriate row start to append the new data
                    row_end   = ff * ArrayLength; % Define the appropriate row end to append the new data
                    % Concatenated data
                    ConcateArrayPullLength(row_start:row_end,:)=yPullingLength'; % Append the new data into the concatenated vector
                    ConcateArrayPullLength(row_start:row_end,:)=ConcateArrayPullLength(row_start:row_end,:).*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected'; % Set non-selected force curves from the concatenated arrays to zero
                    ConcateArrayPullLength(ConcateArrayPullLength==0)=nan; % Replace zero entries by nan´s
                else
                end
                if ~isempty(ySnapInLength)
                    % Determine the number of rows per force map
                    ArrayLength=length(ySnapInLength); % Define the length of the array
                    row_start = ((ff-1) * ArrayLength) + 1; % Define the appropriate row start to append the new data
                    row_end   = ff * ArrayLength; % Define the appropriate row end to append the new data
                    % Concatenated data
                    ConcateArraySnapIn(row_start:row_end,:)=ySnapInLength'; % Append the new data into the concatenated vector
                    ConcateArraySnapIn(row_start:row_end,:)=ConcateArraySnapIn(row_start:row_end,:).*obj.FM{FMIdxArray(ff)}.SMFSFlag.Selected'; % Set non-selected force curves from the concatenated arrays to zero
                    ConcateArraySnapIn(ConcateArraySnapIn==0)=nan; % Replace zero entries by nan´s
                else
                end
                % Min max values and location
                if ff==1
                    [AdhMaxAppMinArray(ff,1),AdhMaxAppMinFcArray(ff,1)]=min(ConcateArrayAdhMaxApp(1:100),[],'omitnan');
                    [AdhMaxAppMaxArray(ff,1),AdhMaxAppMaxFcArray(ff,1)]=max(ConcateArrayAdhMaxApp(1:100),[],'omitnan');
                    [AdhMaxRetMinArray(ff,1),AdhMaxRetMinFcArray(ff,1)]=min(ConcateArrayAdhMaxRet(1:100),[],'omitnan');
                    [AdhMaxRetMaxArray(ff,1),AdhMaxRetMaxFcArray(ff,1)]=max(ConcateArrayAdhMaxRet(1:100),[],'omitnan');
                    [AdhMaxRetUnbindingMinArray(ff,1),AdhMaxRetUnbindingMinFcArray(ff,1)]=min(ConcateArrayAdhMaxRetUnbinding(1:100),[],'omitnan');
                    [AdhMaxRetUnbindingMaxArray(ff,1),AdhMaxRetUnbindingMaxFcArray(ff,1)]=max(ConcateArrayAdhMaxRetUnbinding(1:100),[],'omitnan');
                    [AdhEneAppMinArray(ff,1),AdhEneAppMinFcArray(ff,1)]=min(ConcateArrayAdhEneApp(1:100),[],'omitnan');
                    [AdhEneAppMaxArray(ff,1),AdhEneAppMaxFcArray(ff,1)]=max(ConcateArrayAdhEneApp(1:100),[],'omitnan');
                    [AdhEneRetMinArray(ff,1),AdhEneRetMinFcArray(ff,1)]=min(ConcateArrayAdhEneRet(1:100),[],'omitnan');
                    [AdhEneRetMaxArray(ff,1),AdhEneRetMaxFcArray(ff,1)]=max(ConcateArrayAdhEneRet(1:100),[],'omitnan');
                    [PullLengthMinArray(ff,1),PullLengthMinFcArray(ff,1)]=min(ConcateArrayPullLength(1:100),[],'omitnan');
                    [PullLengthMaxArray(ff,1),PullLengthMaxFcArray(ff,1)]=max(ConcateArrayPullLength(1:100),[],'omitnan');
                    [SnapInMinArray(ff,1),SnapInMinFcArray(ff,1)]=min(ConcateArraySnapIn(1:100),[],'omitnan');
                    [SnapInMaxArray(ff,1),SnapInMaxFcArray(ff,1)]=max(ConcateArraySnapIn(1:100),[],'omitnan');
                else
                    [AdhMaxAppMinArray(ff,1),AdhMaxAppMinFcArray(ff,1)]=min(ConcateArrayAdhMaxApp(((ff-1)*100)+1:ff*100),[],'omitnan');
                    [AdhMaxAppMaxArray(ff,1),AdhMaxAppMaxFcArray(ff,1)]=max(ConcateArrayAdhMaxApp(((ff-1)*100)+1:ff*100),[],'omitnan');
                    [AdhMaxRetMinArray(ff,1),AdhMaxRetMinFcArray(ff,1)]=min(ConcateArrayAdhMaxRet(((ff-1)*100)+1:ff*100),[],'omitnan');
                    [AdhMaxRetMaxArray(ff,1),AdhMaxRetMaxFcArray(ff,1)]=max(ConcateArrayAdhMaxRet(((ff-1)*100)+1:ff*100),[],'omitnan');
                    [AdhMaxRetUnbindingMinArray(ff,1),AdhMaxRetUnbindingMinFcArray(ff,1)]=min(ConcateArrayAdhMaxRetUnbinding(((ff-1)*100)+1:ff*100),[],'omitnan');
                    [AdhMaxRetUnbindingMaxArray(ff,1),AdhMaxRetUnbindingMaxFcArray(ff,1)]=max(ConcateArrayAdhMaxRetUnbinding(((ff-1)*100)+1:ff*100),[],'omitnan');
                    [AdhEneAppMinArray(ff,1),AdhEneAppMinFcArray(ff,1)]=min(ConcateArrayAdhEneApp(((ff-1)*100)+1:ff*100),[],'omitnan');
                    [AdhEneAppMaxArray(ff,1),AdhEneAppMaxFcArray(ff,1)]=max(ConcateArrayAdhEneApp(((ff-1)*100)+1:ff*100),[],'omitnan');
                    [AdhEneRetMinArray(ff,1),AdhEneRetMinFcArray(ff,1)]=min(ConcateArrayAdhEneRet(((ff-1)*100)+1:ff*100),[],'omitnan');
                    [AdhEneRetMaxArray(ff,1),AdhEneRetMaxFcArray(ff,1)]=max(ConcateArrayAdhEneRet(((ff-1)*100)+1:ff*100),[],'omitnan');
                    [PullLengthMinArray(ff,1),PullLengthMinFcArray(ff,1)]=min(ConcateArrayPullLength(((ff-1)*100)+1:ff*100),[],'omitnan');
                    [PullLengthMaxArray(ff,1),PullLengthMaxFcArray(ff,1)]=max(ConcateArrayPullLength(((ff-1)*100)+1:ff*100),[],'omitnan');
                    [SnapInMinArray(ff,1),SnapInMinFcArray(ff,1)]=min(ConcateArraySnapIn(((ff-1)*100)+1:ff*100),[],'omitnan');
                    [SnapInMaxArray(ff,1),SnapInMaxFcArray(ff,1)]=max(ConcateArraySnapIn(((ff-1)*100)+1:ff*100),[],'omitnan');
                end
            end
            % Sorting the force maps and force curves based on measurement time
            [~, FMSortDateTimeIdxFc]=sort(FMDateTimeNumber); % Sort force maps based on the time of measurement
            FMIndexChrono=FMIdxArray(FMSortDateTimeIdxFc); % Order the force maps chronologically
            [FMSortDateTimeFcArray, FMSortDateTimeIdxFcArray]=sort(FMDateTimeNumberFcArray);
            FcNumArrayChrono=FcNumArray(FMSortDateTimeIdxFcArray);
            FMIndexChronoArray=FMIndexArray(FMSortDateTimeIdxFcArray);
            % Statistics
            AdhMaxAppMean=mean(ConcateArrayAdhMaxApp,'omitnan');
            AdhMaxAppMedian=median(ConcateArrayAdhMaxApp,'omitnan');
            AdhMaxAppStd=std(ConcateArrayAdhMaxApp,'omitnan');
            [AdhMaxAppMin,AdhMaxAppMinIdx]=min(AdhMaxAppMinArray);
            [AdhMaxAppMax,AdhMaxAppMaxIdx]=max(AdhMaxAppMaxArray);
            AdhMaxAppMinFM=FMIdxArray(AdhMaxAppMinIdx);
            AdhMaxAppMaxFM=FMIdxArray(AdhMaxAppMaxIdx);
            AdhMaxAppMinFc=AdhMaxAppMinFcArray(AdhMaxAppMinIdx,1);
            AdhMaxAppMaxFc=AdhMaxAppMaxFcArray(AdhMaxAppMaxIdx,1);
            AdhMaxRetMean=mean(ConcateArrayAdhMaxRet,'omitnan');
            AdhMaxRetMedian=median(ConcateArrayAdhMaxRet,'omitnan');
            AdhMaxRetStd=std(ConcateArrayAdhMaxRet,'omitnan');
            [AdhMaxRetMin,AdhMaxRetMinIdx]=min(AdhMaxRetMinArray);
            [AdhMaxRetMax,AdhMaxRetMaxIdx]=max(AdhMaxRetMaxArray);
            AdhMaxRetMinFM=FMIdxArray(AdhMaxRetMinIdx);
            AdhMaxRetMaxFM=FMIdxArray(AdhMaxRetMaxIdx);
            AdhMaxRetMinFc=AdhMaxRetMinFcArray(AdhMaxRetMinIdx,1);
            AdhMaxRetMaxFc=AdhMaxRetMaxFcArray(AdhMaxRetMaxIdx,1);
            AdhMaxRetUnbindingMean=mean(ConcateArrayAdhMaxRetUnbinding,'omitnan');
            AdhMaxRetUnbindingMedian=median(ConcateArrayAdhMaxRetUnbinding,'omitnan');
            AdhMaxRetUnbindingStd=std(ConcateArrayAdhMaxRetUnbinding,'omitnan');
            [AdhMaxRetUnbindingMin,AdhMaxRetUnbindingMinIdx]=min(AdhMaxRetUnbindingMinArray);
            [AdhMaxRetUnbindingMax,AdhMaxRetUnbindingMaxIdx]=max(AdhMaxRetUnbindingMaxArray);
            AdhMaxRetUnbindingMinFM=FMIdxArray(AdhMaxRetUnbindingMinIdx);
            AdhMaxRetUnbindingMaxFM=FMIdxArray(AdhMaxRetUnbindingMaxIdx);
            AdhMaxRetUnbindingMinFc=AdhMaxRetUnbindingMinFcArray(AdhMaxRetUnbindingMinIdx,1);
            AdhMaxRetUnbindingMaxFc=AdhMaxRetUnbindingMaxFcArray(AdhMaxRetUnbindingMaxIdx,1);
            AdhEneAppMean=mean(ConcateArrayAdhEneApp,'omitnan');
            AdhEneAppMedian=median(ConcateArrayAdhEneApp,'omitnan');
            AdhEneAppStd=std(ConcateArrayAdhEneApp,'omitnan');
            [AdhEneAppMin,AdhEneAppMinIdx]=min(AdhEneAppMinArray);
            [AdhEneAppMax,AdhEneAppMaxIdx]=max(AdhEneAppMaxArray);
            AdhEneAppMinFM=FMIdxArray(AdhEneAppMinIdx);
            AdhEneAppMaxFM=FMIdxArray(AdhEneAppMaxIdx);
            AdhEneAppMinFc=AdhEneAppMinFcArray(AdhEneAppMinIdx,1);
            AdhEneAppMaxFc=AdhEneAppMaxFcArray(AdhEneAppMaxIdx,1);
            AdhEneRetMean=mean(ConcateArrayAdhEneRet,'omitnan');
            AdhEneRetMedian=median(ConcateArrayAdhEneRet,'omitnan');
            AdhEneRetStd=std(ConcateArrayAdhEneRet,'omitnan');
            [AdhEneRetMin,AdhEneRetMinIdx]=min(AdhEneRetMinArray);
            [AdhEneRetMax,AdhEneRetMaxIdx]=max(AdhEneRetMaxArray);
            AdhEneRetMinFM=FMIdxArray(AdhEneRetMinIdx);
            AdhEneRetMaxFM=FMIdxArray(AdhEneRetMaxIdx);
            AdhEneRetMinFc=AdhEneRetMinFcArray(AdhEneRetMinIdx,1);
            AdhEneRetMaxFc=AdhEneRetMaxFcArray(AdhEneRetMaxIdx,1);
            PullLengthMean=mean(ConcateArrayPullLength,'omitnan');
            PullLengthMedian=median(ConcateArrayPullLength,'omitnan');
            PullLengthStd=std(ConcateArrayPullLength,'omitnan');
            [PullLengthMin,PullLengthMinIdx]=min(PullLengthMinArray);
            [PullLengthMax,PullLengthMaxIdx]=max(PullLengthMaxArray);
            PullLengthMinFM=FMIdxArray(PullLengthMinIdx);
            PullLengthMaxFM=FMIdxArray(PullLengthMaxIdx);
            PullLengthMinFc=PullLengthMinFcArray(PullLengthMinIdx,1);
            PullLengthMaxFc=PullLengthMaxFcArray(PullLengthMaxIdx,1);
            SnapInMean=mean(ConcateArraySnapIn,'omitnan');
            SnapInMedian=median(ConcateArraySnapIn,'omitnan');
            SnapInStd=std(ConcateArraySnapIn,'omitnan');
            [SnapInMin,SnapInMinIdx]=min(SnapInMinArray);
            [SnapInMax,SnapInMaxIdx]=max(SnapInMaxArray);
            SnapInMinFM=FMIdxArray(SnapInMinIdx);
            SnapInMaxFM=FMIdxArray(SnapInMaxIdx);
            SnapInMinFc=SnapInMinFcArray(PullLengthMinIdx,1);
            SnapInMaxFc=SnapInMaxFcArray(PullLengthMaxIdx,1);
            %% Data selection based on input parameters
            % Determine unique entries
            ExtVelocityValues=unique(FMExtVelocity)';
            RetVelocityValues=unique(FMRetVelocity)';
            HoldingTimeValues=unique(FMHoldingTime)';
            SubstrateValues=unique(FMSubstrate)';
            EnvCondValues=unique(FMEnvCond)';
            ChipCantValues=unique(FMChipCant)';
            ChipboxValues=unique(FMChipbox)';
            LinkerValues=unique(FMLinker)';
            % Preallocate
            ExtVelocityFMIdx=zeros(length(FMExtVelocity),length(ExtVelocityValues));
            ExtVelocityConcateIdx=zeros(length(FMIndexArray),length(ExtVelocityValues));
            RetVelocityFMIdx=zeros(length(FMRetVelocity),length(RetVelocityValues));
            RetVelocityConcateIdx=zeros(length(FMIndexArray),length(RetVelocityValues));
            HoldingTimeFMIdx=zeros(length(FMHoldingTime),length(HoldingTimeValues));
            HoldingTimeConcateIdx=zeros(length(FMIndexArray),length(HoldingTimeValues));
            SubstrateFMIdx=zeros(length(FMSubstrate),length(SubstrateValues));
            SubstrateConcateIdx=zeros(length(FMIndexArray),length(SubstrateValues));
            EnvCondFMIdx=zeros(length(FMEnvCond),length(EnvCondValues));
            EnvCondConcateIdx=zeros(length(FMIndexArray),length(EnvCondValues));
            ChipCantFMIdx=zeros(length(FMChipCant),length(ChipCantValues));
            ChipCantConcateIdx=zeros(length(FMIndexArray),length(ChipCantValues));
            ChipboxFMIdx=zeros(length(FMChipbox),length(ChipboxValues));
            ChipboxConcateIdx=zeros(length(FMIndexArray),length(ChipboxValues));
            LinkerFMIdx=zeros(length(FMLinker),length(LinkerValues));
            LinkerConcateIdx=zeros(length(FMIndexArray),length(LinkerValues));
            for Fm=1:length(ExtVelocityValues)
                [~, Loc] =ismember(FMExtVelocity,ExtVelocityValues(Fm)); % Identify the force maps that have the same value
                Loc=find(Loc); % Use find to get ride of nonzeros
                FMIdx=FMIdxArray(Loc); % Identify same entries in the concatenate data of the results table
                [~, ConcateLoc]=ismember(FMIndexArray,FMIdx); % Identify same entries in the concatenate data of the results table
                ConcateLoc=find(ConcateLoc); % Use find to get ride of nonzeros
                % Allocate the data
                ExtVelocityFMIdx(1:length(Loc),Fm)=FMIdx;
                ExtVelocityConcateIdx(1:length(ConcateLoc),Fm)=ConcateLoc;
            end
            for Fm=1:length(RetVelocityValues)
                [~, Loc] =ismember(FMRetVelocity,RetVelocityValues(Fm)); % Identify the force maps that have the same value
                Loc=find(Loc); % Use find to get ride of nonzeros
                FMIdx=FMIdxArray(Loc); % Identify same entries in the concatenate data of the results table
                [~, ConcateLoc]=ismember(FMIndexArray,FMIdx); % Identify same entries in the concatenate data of the results table
                ConcateLoc=find(ConcateLoc); % Use find to get ride of nonzeros
                % Allocate the data
                RetVelocityFMIdx(1:length(Loc),Fm)=FMIdx;
                RetVelocityConcateIdx(1:length(ConcateLoc),Fm)=ConcateLoc;
            end
            for Fm=1:length(HoldingTimeValues)
                [~, Loc] =ismember(FMHoldingTime,HoldingTimeValues(Fm)); % Identify the force maps that have the same value
                Loc=find(Loc); % Use find to get ride of nonzeros
                FMIdx=FMIdxArray(Loc); % Identify same entries in the concatenate data of the results table
                [~, ConcateLoc]=ismember(FMIndexArray,FMIdx); % Identify same entries in the concatenate data of the results table
                ConcateLoc=find(ConcateLoc); % Use find to get ride of nonzeros
                % Allocate the data
                HoldingTimeFMIdx(1:length(Loc),Fm)=FMIdx;
                HoldingTimeConcateIdx(1:length(ConcateLoc),Fm)=ConcateLoc;
            end
            for Fm=1:length(SubstrateValues)
                [~, Loc] =ismember(FMSubstrate,SubstrateValues(Fm)); % Identify the force maps that have the same value
                Loc=find(Loc); % Use find to get ride of nonzeros
                FMIdx=FMIdxArray(Loc); % Identify same entries in the concatenate data of the results table
                [~, ConcateLoc]=ismember(FMIndexArray,FMIdx); % Identify same entries in the concatenate data of the results table
                ConcateLoc=find(ConcateLoc); % Use find to get ride of nonzeros
                % Allocate the data
                SubstrateFMIdx(1:length(Loc),Fm)=FMIdx;
                SubstrateConcateIdx(1:length(ConcateLoc),Fm)=ConcateLoc;
            end
            for Fm=1:length(EnvCondValues)
                [~, Loc] =ismember(FMEnvCond,EnvCondValues(Fm)); % Identify the force maps that have the same value
                Loc=find(Loc); % Use find to get ride of nonzeros
                FMIdx=FMIdxArray(Loc); % Identify same entries in the concatenate data of the results table
                [~, ConcateLoc]=ismember(FMIndexArray,FMIdx); % Identify same entries in the concatenate data of the results table
                ConcateLoc=find(ConcateLoc); % Use find to get ride of nonzeros
                % Allocate the data
                EnvCondFMIdx(1:length(Loc),Fm)=FMIdx;
                EnvCondConcateIdx(1:length(ConcateLoc),Fm)=ConcateLoc;
            end
            for Fm=1:length(ChipCantValues)
                [~, Loc] =ismember(FMChipCant,ChipCantValues(Fm)); % Identify the force maps that have the same value
                Loc=find(Loc); % Use find to get ride of nonzeros
                FMIdx=FMIdxArray(Loc); % Identify same entries in the concatenate data of the results table
                [~, ConcateLoc]=ismember(FMIndexArray,FMIdx); % Identify same entries in the concatenate data of the results table
                ConcateLoc=find(ConcateLoc); % Use find to get ride of nonzeros
                % Allocate the data
                ChipCantFMIdx(1:length(Loc),Fm)=FMIdx;
                ChipCantConcateIdx(1:length(ConcateLoc),Fm)=ConcateLoc;
            end
            for Fm=1:length(ChipboxValues)
                [~, Loc] =ismember(FMChipbox,ChipboxValues(Fm)); % Identify the force maps that have the same value
                Loc=find(Loc); % Use find to get ride of nonzeros
                FMIdx=FMIdxArray(Loc); % Identify same entries in the concatenate data of the results table
                [~, ConcateLoc]=ismember(FMIndexArray,FMIdx); % Identify same entries in the concatenate data of the results table
                ConcateLoc=find(ConcateLoc); % Use find to get ride of nonzeros
                % Allocate the data
                ChipboxFMIdx(1:length(Loc),Fm)=FMIdx;
                ChipboxConcateIdx(1:length(ConcateLoc),Fm)=ConcateLoc;
            end
            for Fm=1:length(LinkerValues)
                [~, Loc] =ismember(FMLinker,LinkerValues(Fm)); % Identify the force maps that have the same value
                Loc=find(Loc); % Use find to get ride of nonzeros
                FMIdx=FMIdxArray(Loc); % Identify same entries in the concatenate data of the results table
                [~, ConcateLoc]=ismember(FMIndexArray,FMIdx); % Identify same entries in the concatenate data of the results table
                ConcateLoc=find(ConcateLoc); % Use find to get ride of nonzeros
                % Allocate the data
                LinkerFMIdx(1:length(Loc),Fm)=FMIdx;
                LinkerConcateIdx(1:length(ConcateLoc),Fm)=ConcateLoc;
            end
            %% SMFS Results structure
            % Check entry
            if ~isempty(obj.SMFSResults)
                jj=length(obj.SMFSResults)+1;
            else
                jj=1;
            end
            % Debugging
            % jj=11
            % Allocate data
            obj.SMFSResults{jj,1}.Concatenate(1).FMID=FMIDArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMIndex=FMIndexArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMNum=FMNumArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMIndexChrono=FMIndexChronoArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FcNum=FcNumArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FcNumChrono=FcNumArrayChrono;
            obj.SMFSResults{jj,1}.Concatenate(1).FMExtVelocity=FMExtVelocityArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMRetVelocity=FMRetVelocityArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMHoldingTime=FMHoldingTimeArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMSubstrate=FMSubstrateArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMEnvCond=FMEnvCondArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMChipCant=FMChipCantArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMChipbox=FMChipboxArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMLinker=FMLinkerArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMDateTime=FMDateTimeFcArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMDateTimeNumber=FMDateTimeNumberFcArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMDateTimeNumberSort=FMSortDateTimeFcArray;
            obj.SMFSResults{jj,1}.Concatenate(1).FMDateTimeNumberSortIdx=FMSortDateTimeIdxFcArray;
            obj.SMFSResults{jj,1}.Concatenate(1).AdhMaxApp=ConcateArrayAdhMaxApp;
            obj.SMFSResults{jj,1}.Concatenate(1).AdhMaxRet= ConcateArrayAdhMaxRet;
            obj.SMFSResults{jj,1}.Concatenate(1).AdhUnbinding=ConcateArrayAdhMaxRetUnbinding;
            obj.SMFSResults{jj,1}.Concatenate(1).AdhEneApp=ConcateArrayAdhEneApp;
            obj.SMFSResults{jj,1}.Concatenate(1).AdhEneRet=ConcateArrayAdhEneRet;
            obj.SMFSResults{jj,1}.Concatenate(1).PullingLength=ConcateArrayPullLength;
            obj.SMFSResults{jj,1}.Concatenate(1).SnapInLength=ConcateArraySnapIn;
            obj.SMFSResults{jj,1}.Flags(1).SMFSFlagUncorrupt=sum(SMFSFlagUncorrupt);
            obj.SMFSResults{jj,1}.Flags(1).UncorruptPct=(sum(SMFSFlagUncorrupt)/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Flags(1).SMFSFlagSelected=sum(SMFSFlagSelected);
            obj.SMFSResults{jj,1}.Flags(1).SelectedPct=(sum(SMFSFlagSelected)/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Flags(1).SMFSFlagAppMinCrit=sum(SMFSFlagAppMinCrit);
            obj.SMFSResults{jj,1}.Flags(1).AppMinCritPct=(sum(SMFSFlagAppMinCrit)/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Flags(1).SMFSFlagRetMinCrit=sum(SMFSFlagRetMinCrit);
            obj.SMFSResults{jj,1}.Flags(1).RetMinCritPct=(sum(SMFSFlagRetMinCrit)/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Flags(1).SMFSFlagLengthRequisite=sum(SMFSFlagLengthRequisite);
            obj.SMFSResults{jj,1}.Flags(1).LengthRequisitePct=(sum(SMFSFlagLengthRequisite)/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Flags(1).SMFSFlagFit=sum(SMFSFlagFit);
            obj.SMFSResults{jj,1}.Flags(1).FitPct=(sum(SMFSFlagFit)/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Flags(1).SMFSFlagFitLinear=sum(SMFSFlagFitLinear);
            obj.SMFSResults{jj,1}.Flags(1).FitLinearPct=(sum(SMFSFlagFitLinear)/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Flags(1).SMFSFlagFitSinoidal=sum(SMFSFlagFitSinoidal);
            obj.SMFSResults{jj,1}.Flags(1).FitSinoidalPct=(sum(SMFSFlagFitSinoidal)/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Flags(1).SMFSFlagSnapIn=sum(SMFSFlagSnapIn);
            obj.SMFSResults{jj,1}.Flags(1).SnapInPct=(sum(SMFSFlagSnapIn)/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Flags(1).SMFSFlagPullingLength=sum(SMFSFlagPullingLength);
            obj.SMFSResults{jj,1}.Flags(1).PullingLengthPct=(sum(SMFSFlagPullingLength)/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Data(1).FMIndex=FMIdxArray;
            obj.SMFSResults{jj,1}.Data(1).FMIndexChrono=FMIndexChrono;
            obj.SMFSResults{jj,1}.Data(1).FMID=FMID;
            obj.SMFSResults{jj,1}.Data(1).FMNum=FMNum;
            obj.SMFSResults{jj,1}.Data(1).FMExtVelocity=FMExtVelocity;
            obj.SMFSResults{jj,1}.Data(1).FMRetVelocity=FMRetVelocity;
            obj.SMFSResults{jj,1}.Data(1).FMHoldingTime=FMHoldingTime;
            obj.SMFSResults{jj,1}.Data(1).FMSubstrate=FMSubstrate;
            obj.SMFSResults{jj,1}.Data(1).FMEnvCond=FMEnvCond;
            obj.SMFSResults{jj,1}.Data(1).FMChipCant=FMChipCant;
            obj.SMFSResults{jj,1}.Data(1).FMChipbox=FMChipbox;
            obj.SMFSResults{jj,1}.Data(1).FMLinker=FMLinker;
            obj.SMFSResults{jj,1}.Data(1).FMDateTime=FMDateTime;
            obj.SMFSResults{jj,1}.Data(1).FMDateTimeSortIdx=FMSortDateTimeIdxFc;
            obj.SMFSResults{jj,1}.Data(1).TotalNumFc=sum(TotalNumFc);
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedAdhMaxApp=nnz(~isnan(ConcateArrayAdhMaxApp));
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedAdhMaxAppPct=(nnz(~isnan(ConcateArrayAdhMaxApp))/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedAdhMaxRet=nnz(~isnan(ConcateArrayAdhMaxRet));
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedAdhMaxRetPct=(nnz(~isnan(ConcateArrayAdhMaxRet))/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedAdhUnbinding=nnz(~isnan(ConcateArrayAdhMaxRetUnbinding));
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedAdhUnbindingPct=(nnz(~isnan(ConcateArrayAdhMaxRetUnbinding))/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedAdhEneApp=nnz(~isnan(ConcateArrayAdhEneApp));
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedAdhEneAppPct=(nnz(~isnan(ConcateArrayAdhEneApp))/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedAdhEneRet=nnz(~isnan(ConcateArrayAdhEneRet));
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedAdhEneRetPct=(nnz(~isnan(ConcateArrayAdhEneRet))/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedyPullingLength=nnz(~isnan(ConcateArrayPullLength));
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedPullingLengthPct=(nnz(~isnan(ConcateArrayPullLength))/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedySnapInLength=nnz(~isnan(ConcateArraySnapIn));
            obj.SMFSResults{jj,1}.Data(1).SumNumFcAnalysedySnapInLengthPct=(nnz(~isnan(ConcateArraySnapIn))/sum(TotalNumFc))*100;
            obj.SMFSResults{jj,1}.Data(1).NumFcAnalysedAdhMaxApp=FCAnalysedyAdhMaxApp;
            obj.SMFSResults{jj,1}.Data(1).NumFcAnalysedAdhMaxRet=FCAnalysedyAdhMaxRet;
            obj.SMFSResults{jj,1}.Data(1).NumFcAnalysedAdhUnbinding=FCAnalysedyAdhUnbind;
            obj.SMFSResults{jj,1}.Data(1).NumFcAnalysedAdhEneApp=FCAnalysedyAdhEneApp;
            obj.SMFSResults{jj,1}.Data(1).NumFcAnalysedAdhEneRet=FCAnalysedyAdhEneRet;
            obj.SMFSResults{jj,1}.Data(1).NumFcAnalysedyPullingLength=FCAnalysedyPullingLength;
            obj.SMFSResults{jj,1}.Data(1).NumFcAnalysedySnapInLength=FCAnalysedySnapInLength;
            obj.SMFSResults{jj,1}.Data(1).AdhMaxApp=yAdhMaxAppAll;
            obj.SMFSResults{jj,1}.Data(1).AdhMaxRet=yAdhMaxRetAll;
            obj.SMFSResults{jj,1}.Data(1).AdhUnbinding=yAdhUnbindingAll;
            obj.SMFSResults{jj,1}.Data(1).AdhEneApp=yAdhEneAppAll;
            obj.SMFSResults{jj,1}.Data(1).AdhEneRet=yAdhEneRetAll;
            obj.SMFSResults{jj,1}.Data(1).PullingLength=yPullingLengthAll;
            obj.SMFSResults{jj,1}.Data(1).SnapInLength=ySnapInLengthAll;
            obj.SMFSResults{jj,1}.Parameters(1).ExtendVelocity=ExtVelocityValue;
            obj.SMFSResults{jj,1}.Parameters(1).RetractVelocity=RetVelocityValue;
            obj.SMFSResults{jj,1}.Parameters(1).HoldingTime=HoldingTimeValue;
            obj.SMFSResults{jj,1}.Parameters(1).Substrate=SubstrateValue;
            obj.SMFSResults{jj,1}.Parameters(1).Medium=EnvCondValue;
            obj.SMFSResults{jj,1}.Parameters(1).ChipCantilever=ChipCantValue;
            obj.SMFSResults{jj,1}.Parameters(1).Chipbox=ChipboxValue;
            obj.SMFSResults{jj,1}.Parameters(1).Linker=LinkerValue;
            obj.SMFSResults{jj,1}.Selection(1).ExtVelocityParameters=ExtVelocityValues;
            obj.SMFSResults{jj,1}.Selection(1).ExtVelocityFMIdx=ExtVelocityFMIdx;
            obj.SMFSResults{jj,1}.Selection(1).ExtVelocityConcateIdx=ExtVelocityConcateIdx;
            obj.SMFSResults{jj,1}.Selection(1).RetVelocityParameters=RetVelocityValues;
            obj.SMFSResults{jj,1}.Selection(1).RetVelocityFMIdx=RetVelocityFMIdx;
            obj.SMFSResults{jj,1}.Selection(1).RetVelocityConcateIdx=RetVelocityConcateIdx;
            obj.SMFSResults{jj,1}.Selection(1).HoldingTimeParameters=HoldingTimeValues;
            obj.SMFSResults{jj,1}.Selection(1).HoldingTimeFMIdx=HoldingTimeFMIdx;
            obj.SMFSResults{jj,1}.Selection(1).HoldingTimeConcateIdx=HoldingTimeConcateIdx;
            obj.SMFSResults{jj,1}.Selection(1).SubstrateParameters=SubstrateValues;
            obj.SMFSResults{jj,1}.Selection(1).SubstrateFMIdx=SubstrateFMIdx;
            obj.SMFSResults{jj,1}.Selection(1).SubstrateConcateIdx=SubstrateConcateIdx;
            obj.SMFSResults{jj,1}.Selection(1).EnvCondParameters=EnvCondValues;
            obj.SMFSResults{jj,1}.Selection(1).EnvCondFMIdx=EnvCondFMIdx;
            obj.SMFSResults{jj,1}.Selection(1).EnvCondConcateIdx=EnvCondConcateIdx;
            obj.SMFSResults{jj,1}.Selection(1).ChipCantParameters=ChipCantValues;
            obj.SMFSResults{jj,1}.Selection(1).ChipCantFMIdx=ChipCantFMIdx;
            obj.SMFSResults{jj,1}.Selection(1).ChipCantConcateIdx=ChipCantConcateIdx;
            obj.SMFSResults{jj,1}.Selection(1).ChipboxParameters=ChipboxValues;
            obj.SMFSResults{jj,1}.Selection(1).ChipboxFMIdx=ChipboxFMIdx;
            obj.SMFSResults{jj,1}.Selection(1).ChipboxConcateIdx=ChipboxConcateIdx;
            obj.SMFSResults{jj,1}.Selection(1).LinkerParameters=LinkerValues;
            obj.SMFSResults{jj,1}.Selection(1).LinkerFMIdx=LinkerFMIdx;
            obj.SMFSResults{jj,1}.Selection(1).LinkerConcateIdx=LinkerConcateIdx;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppMean=AdhMaxAppMean;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppMedian=AdhMaxAppMedian;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppStd=AdhMaxAppStd;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppMinArray=AdhMaxAppMinArray;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppMin=AdhMaxAppMin;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppMaxArray=AdhMaxAppMaxArray;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppMax=AdhMaxAppMax;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppMinFM=AdhMaxAppMinFM;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppMaxFM=AdhMaxAppMaxFM;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppMinFcArray=AdhMaxAppMinFcArray;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppMinFc=AdhMaxAppMinFc;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppMaxFcArray=AdhMaxAppMaxFcArray;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppMaxFc=AdhMaxAppMaxFc;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetMean=AdhMaxRetMean;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetMedian=AdhMaxRetMedian;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetStd=AdhMaxRetStd;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetMinArray=AdhMaxRetMinArray;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetMin=AdhMaxRetMin;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetMaxArray=AdhMaxRetMaxArray;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetMax=AdhMaxRetMax;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetMinFM=AdhMaxRetMinFM;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetMaxFM=AdhMaxRetMaxFM;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetMinFcArray=AdhMaxRetMinFcArray;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetMinFc=AdhMaxRetMinFc;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetMaxFcArray=AdhMaxRetMaxFcArray;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetMaxFc=AdhMaxRetMaxFc;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingMean=AdhMaxRetUnbindingMean;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingMedian=AdhMaxRetUnbindingMedian;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingStd=AdhMaxRetUnbindingStd;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingMinArray=AdhMaxRetUnbindingMinArray;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingMin=AdhMaxRetUnbindingMin;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingMaxArray=AdhMaxRetUnbindingMaxArray;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingMax=AdhMaxRetUnbindingMax;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingMinFM=AdhMaxRetUnbindingMinFM;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingMaxFM=AdhMaxRetUnbindingMaxFM;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingMinFcArray=AdhMaxRetUnbindingMinFcArray;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingMinFc=AdhMaxRetUnbindingMinFc;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingMaxFcArray=AdhMaxRetUnbindingMaxFcArray;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingMaxFc=AdhMaxRetUnbindingMaxFc;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppMean=AdhEneAppMean;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppMedian=AdhEneAppMedian;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppStd=AdhEneAppStd;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppMinArray=AdhEneAppMinArray;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppMin=AdhEneAppMin;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppMaxArray=AdhEneAppMaxArray;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppMax=AdhEneAppMax;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppMinFM=AdhEneAppMinFM;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppMaxFM=AdhEneAppMaxFM;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppMinFcArray=AdhEneAppMinFcArray;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppMinFc=AdhEneAppMinFc;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppMaxFcArray=AdhEneAppMaxFcArray;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppMaxFc=AdhEneAppMaxFc;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetMean=AdhEneRetMean;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetMedian=AdhEneRetMedian;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetStd=AdhEneRetStd;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetMinArray=AdhEneRetMinArray;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetMin=AdhEneRetMin;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetMaxArray=AdhEneRetMaxArray;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetMax=AdhEneRetMax;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetMinFM=AdhEneRetMinFM;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetMaxFM=AdhEneRetMaxFM;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetMinFcArray=AdhEneRetMinFcArray;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetMinFc=AdhEneRetMinFc;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetMaxFcArray=AdhEneRetMaxFcArray;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetMaxFc=AdhEneRetMaxFc;
            obj.SMFSResults{jj,1}.Results(1).PullLengthMean=PullLengthMean;
            obj.SMFSResults{jj,1}.Results(1).PullLengthMedian=PullLengthMedian;
            obj.SMFSResults{jj,1}.Results(1).PullLengthStd=PullLengthStd;
            obj.SMFSResults{jj,1}.Results(1).PullLengthMinArray=PullLengthMinArray;
            obj.SMFSResults{jj,1}.Results(1).PullLengthMin=PullLengthMin;
            obj.SMFSResults{jj,1}.Results(1).PullLengthMinFM=PullLengthMinFM;
            obj.SMFSResults{jj,1}.Results(1).PullLengthMinFcArray=PullLengthMinFcArray;
            obj.SMFSResults{jj,1}.Results(1).PullLengthMinFc=PullLengthMinFc;
            obj.SMFSResults{jj,1}.Results(1).PullLengthMaxArray=PullLengthMaxArray;
            obj.SMFSResults{jj,1}.Results(1).PullLengthMax=PullLengthMax;
            obj.SMFSResults{jj,1}.Results(1).PullLengthMaxFM=PullLengthMaxFM;
            obj.SMFSResults{jj,1}.Results(1).PullLengthMaxFcArray=PullLengthMaxFcArray;
            obj.SMFSResults{jj,1}.Results(1).PullLengthMaxFc=PullLengthMaxFc;
            obj.SMFSResults{jj,1}.Results(1).SnapInMean=SnapInMean;
            obj.SMFSResults{jj,1}.Results(1).SnapInMedian=SnapInMedian;
            obj.SMFSResults{jj,1}.Results(1).SnapInStd=SnapInStd;
            obj.SMFSResults{jj,1}.Results(1).SnapInMinArray=SnapInMinArray;
            obj.SMFSResults{jj,1}.Results(1).SnapInMin=SnapInMin;
            obj.SMFSResults{jj,1}.Results(1).SnapInMinFM=SnapInMinFM;
            obj.SMFSResults{jj,1}.Results(1).SnapInMinFcArray=SnapInMinFcArray;
            obj.SMFSResults{jj,1}.Results(1).SnapInMinFc=SnapInMinFc;
            obj.SMFSResults{jj,1}.Results(1).SnapInMaxArray=SnapInMaxArray;
            obj.SMFSResults{jj,1}.Results(1).SnapInMax=SnapInMax;
            obj.SMFSResults{jj,1}.Results(1).SnapInMaxFM=SnapInMaxFM;
            obj.SMFSResults{jj,1}.Results(1).SnapInMaxFcArray=SnapInMaxFcArray;
            obj.SMFSResults{jj,1}.Results(1).SnapInMaxFc=SnapInMaxFc;
            % Paramater list
            obj.SMFSResultsParameters(jj,:)={jj,ChipboxValue,ChipCantValue,LinkerValue,SubstrateValue,EnvCondValue,ExtVelocityValue,RetVelocityValue,HoldingTimeValue};        
        end


        function SMFS_fine_figure_publication(obj,XMin,XMax,YMin,YMax,FmoI,FcoI,Linker,CArea,Axis)
            % Author: Andreas Rohatschek
            % Description: 
            % In this function force curves can be plotted in publication quality
            % Required input variables:
            % XMin ... Minimum x-axis value in nanometers (nm): double
            % XMax ... Maximum x-axis value in nanometers (nm): double
            % YMin ... Minimum y-axis value in nanoNewton (nN): double
            % YMax ... Maximum y-axis value in nanoNewton (nN): double
            % FmoI ... Force map of interest: double
            % FcoI ... Force curve of interest: double 
            % Linker ... Linker type: string 'Short' or 'Long' or 'No' 
            % CArea ... Color range (shows theoretical Linker-TC complex range): string 'yes' or
            % 'no'
            % Axis ... Show axis in graph: string 'Yes' or 'No'

            %% Function body
            if nargin < 3
                XMin= -inf;
                XMax= inf;
                YMin= -inf;
                YMax= inf;
            end
            % Figure visibility
            %set(groot,'defaultFigureVisible','off')
            set(groot,'defaultFigureVisible','on')
            % Set figure position
            set(groot,'defaultFigurePaperPositionMode','auto')
            % Change into the Folder of Interest
            cd(obj.ExperimentFolder) % Move into the folder            
            % Create folders for saving the produced figures
            foldername='SMFS_fine_figure_publication';    % Defines the folder name
            %foldername='SMFS_fine_figure_publication_sens';    % Defines the folder name
            mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
            currpath=fullfile(obj.ExperimentFolder,foldername);
            cd(currpath);
            % Run the chosen functions  
            obj.FM{FmoI}.fc_fine_figure_publication(XMin,XMax,YMin,YMax,FmoI,FcoI,Linker,CArea,Axis)
            
        end
              
        
       function SMFS_results_gramm_boxplot2_publication(obj,ResultsRow,Linker,xArg,MarkerArg,LegendArg,CBar,Var)
            % Author: Andreas Rohatschek
            % Description:
            % This function allows to plot the analysis results using the
            % gramm toolbox
            % (https://joss.theoj.org/papers/10.21105/joss.00568). The
            % different features (snap-in length, pull-off length) determined are plotted as boxplots to
            % the corresponding force map.
            % Required input variables:
            % ResultsRow ... row in the results structure: double
            % Linker ... Linker type: string , either 'long' or 'short'
            % xArg ... x-axis argument: string, either 'Index' or 'DateTime'
            % MarkerArg ... Marker argument: string, either 'Y' or 'N'
            % LegendArg ... Legend argument: string, either 'Y' or 'N'
            % CBar ... Color Bar: string, either 'Y' or 'N'
            % Var ... Variant: double, Either 1, 2, 3, 4 or 5

            %% Function body
            % Define variables
            yLimMaxFactor=0.22;
            yLimMinFactor=0.25;
            yAxisMinFactor=0.02;
            yAxisMaxFactor=0.25;
            NumColor=6;
            NumLightness=6;
            Res=[1 1 2560 1250]; % Define the figure resolution
            MarkerStyle={'d' 's' 'v' 'o'};
            MarkerSize=10;
            BaseFontSize=32;
            LabelScaling=1.2;
            LegendScaling=1;
            LineWidth=1.5;
            % Define color bar positions
            xCBar1=[0 2 2 0];
            xCBar2=[2 21 21 2];
            xCBar3=[21 121 121 21];
            xCBar4=[121 222 222 121];
            % Color and color maps
            Ochreish=[253 174 97]./255;
            SteelBlue=[116 173 209]./255;
            ColorBrewerMap1=[[253 174 97]./255; % Ochreish
                [116 173 209]./255]; % Steel blueish
            ColorBarMap=[[54 163 0]./255; % Dark green HEX 8DB600
                [206 22 32]./255; % Fire Engine Red HEX CE162
                [0 24 204]./255; % Blue HEX 8DB600
                [135 0 224]./255]; % Violet HEX 8F00FF
            % Change into the Folder of Interest
            cd(obj.ExperimentFolder) % Move into the folder
            % Create folders for saving the produced figures
            foldername='SMFS_results_gramm_boxplot2_publication';    % Defines the folder name
            mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
            currpath=fullfile(obj.ExperimentFolder,foldername);
            cd(currpath);
            % Input variables
            % Linker
            if strcmpi(Linker,'Long')
                LimitLengthRet1=[0 378];
                LimitLengthRet2=[378 522];
                LimitLengthApp=[50 120];
            elseif strcmpi(Linker,'Short')
                LimitLengthRet1=[0 308];
                LimitLengthRet2=[308 463];
                LimitLengthApp=[50 120];
            end
            % xArg
            if strcmpi(xArg,'DateTime')
                LegendxAxis='Date and Time';
                xData=obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSort;
                xDataMin=min(xData);
                xDataMax=max(xData);
                xAxisCorr=(xDataMax-xDataMin)*0.005;
                BoxplotWidth=0.8;
                BoxplotDodge=2;
            elseif strcmpi(xArg,'Index')
               % LegendxAxis='Chronological force set index';
                LegendxAxis='Number of cycles (x100)';
                xData=obj.SMFSResults{ResultsRow}.Concatenate.FMNum;
                xDataMin=min(xData);
                xDataMax=max(xData);
                xAxisCorr=(xDataMax-xDataMin)*0.005; % For BoxplotWidth=10
                %xAxisCorr=(xDataMax-xDataMin)*0.02; % For BoxplotWidth=2 and BoxplotWidth=0.02
                % BoxplotWidth=18;
                % BoxplotWidth=10;
                BoxplotWidth=2.5;
                BoxplotDodge=1;
            end
            % Variant
            if Var==1
                ColorName='Medium';
                LightnessName='Substrate';
                MarkerName='ChipCantilever';
                LightnessData=obj.SMFSResults{ResultsRow}.Concatenate.FMSubstrate(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx);
                ColorData=obj.SMFSResults{ResultsRow}.Concatenate.FMEnvCond(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx);
                MarkerData=obj.SMFSResults{ResultsRow}.Concatenate.FMChipCant(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx);
            elseif Var==2
                ColorName='Approach speed ($\mu$m/s)';
                LightnessName='Retraction speed ($\mu$m/s)';
                MarkerName='Dwell Time (s)';
                FMExtVeloData=obj.SMFSResults{ResultsRow}.Concatenate.FMExtVelocity(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*1e6;
                FMRetVeloData=obj.SMFSResults{ResultsRow}.Concatenate.FMRetVelocity(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*1e6;
                ColorData=FMExtVeloData;
                LightnessData=FMRetVeloData;
                MarkerData=obj.SMFSResults{ResultsRow}.Concatenate.FMHoldingTime(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx);
            elseif Var==3
                ColorName='Dwell Time (s)';
                LightnessName='Retraction speed ($\mu$m/s)';
                MarkerName='Approach speed ($\mu$m/s)';
                FMExtVeloData=obj.SMFSResults{ResultsRow}.Concatenate.FMExtVelocity(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*1e6;
                FMRetVeloData=obj.SMFSResults{ResultsRow}.Concatenate.FMRetVelocity(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*1e6;
                ColorData=obj.SMFSResults{ResultsRow}.Concatenate.FMHoldingTime(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx);
                LightnessData=FMRetVeloData;
                MarkerData=FMExtVeloData;
            elseif Var==4
                ColorName='Retraction speed ($\mu$m/s)';
                LightnessName='Medium';
                MarkerName='Substrate';
                FMRetVeloData=obj.SMFSResults{ResultsRow}.Concatenate.FMRetVelocity(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*1e6;
                ColorData=FMRetVeloData;
                LightnessData=obj.SMFSResults{ResultsRow}.Concatenate.FMEnvCond(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx);
                MarkerData=obj.SMFSResults{ResultsRow}.Concatenate.FMSubstrate(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx);
             elseif Var==5
                ColorName='Dwell Time (s)';
                LightnessName='Medium';
                MarkerName='Substrate';
                ColorData=obj.SMFSResults{ResultsRow}.Concatenate.FMHoldingTime(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx);
                LightnessData=obj.SMFSResults{ResultsRow}.Concatenate.FMEnvCond(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx);
                MarkerData=obj.SMFSResults{ResultsRow}.Concatenate.FMSubstrate(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx);

            end
            %%
            if obj.SMFSResults{ResultsRow}.Parameters.ExtendVelocity==0
                ExtVelocityValueStr='All';
            else
                ExtVelocityValueStr=num2str(round(obj.SMFSResults{ResultsRow}.Parameters.ExtendVelocity*1e9));
            end
            if obj.SMFSResults{ResultsRow}.Parameters.RetractVelocity==0
                RetVelocityValueStr='All';
            else
                RetVelocityValueStr=num2str(round(obj.SMFSResults{ResultsRow}.Parameters.RetractVelocity*1e9));
            end
            if obj.SMFSResults{ResultsRow}.Parameters.HoldingTime==-1
                HoldingTimeValueStr='All';
            else
                HoldingTimeValueStr=num2str(obj.SMFSResults{ResultsRow}.Parameters.HoldingTime);
            end
            % General names
            FigNamePt1=sprintf('SMFSResultRow%d_',ResultsRow);
            FigNamePt2=strcat(ExtVelocityValueStr,{'_'},RetVelocityValueStr,{'_'},HoldingTimeValueStr,{'_'},obj.SMFSResults{ResultsRow}.Parameters.Substrate,{'_'},obj.SMFSResults{ResultsRow}.Parameters.Medium,{'_'},obj.SMFSResults{ResultsRow}.Parameters.ChipCantilever,{'_'},obj.SMFSResults{ResultsRow}.Parameters.Chipbox,{'_'},obj.SMFSResults{ResultsRow}.Parameters.Linker);
            FigNamePt2=char(FigNamePt2);
            FigNamePt3='_Boxplot2';
            %% Gramm object 1
            % Define variables
            Plottitle1=sprintf('%d Force Maps containing %d Force Curves selected',length(obj.SMFSResults{ResultsRow,1}.Data(1).FMIndex),obj.SMFSResults{ResultsRow,1}.Data(1).SumNumFcAnalysedAdhMaxApp);
            LegendyAxis1='Adhesion force (nN)';
            NameSuffix1='_MaxAdhesionForceApproach';
            % Allocate data
            yData1=obj.SMFSResults{ResultsRow}.Concatenate.AdhMaxApp(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*-1e9;
            yData1Min=min(yData1);
            yData1Max=max(yData1);
            yData1diff=yData1Max-yData1Min;
            % Create a gramm object
            if strcmpi(MarkerArg,'Y')
                g1=gramm('x',xData,'y',yData1,...
                    'color',ColorData,...
                    'lightness',LightnessData,...
                    'marker',MarkerData);
            elseif strcmpi(MarkerArg,'N')
                g1=gramm('x',xData,'y',yData1,...
                    'color',ColorData,...
                    'lightness',LightnessData);
            end
            % Plot data
            g1.stat_boxplot('notch',true,...
                'width',BoxplotWidth,...
                'dodge',BoxplotDodge); % Plot data in boxplot
            % Set options
            g1.axe_property('LineWidth',LineWidth,'xlim',[xDataMin-xAxisCorr xDataMax+xAxisCorr],'ylim',[yData1Min-yData1diff*yAxisMinFactor yData1Max+yData1diff*yAxisMaxFactor],'XTick',[0:10:220],'TickDir','out','TickLength',[0.005 0.005]);
            if strcmpi(xArg,'DateTime')
                g1.set_datetick('x',0,'keeplimits') % Format x-axis
            end
            % Color Bar
            if strcmpi(CBar,'Y')
            y1CBar=[yData1Max+yData1diff*yLimMinFactor yData1Max+yData1diff*yLimMinFactor yData1Max+yData1diff*yLimMaxFactor yData1Max+yData1diff*yLimMaxFactor];
            g1.geom_polygon('x',{xCBar1;xCBar2;xCBar3;xCBar4},'y',{y1CBar;y1CBar;y1CBar;y1CBar},'color',ColorBarMap,'alpha',1);
            else
            end
            g1.set_point_options('base_size',MarkerSize)
            %g1.set_title(Plottitle1) %Set figure title
            if strcmpi(MarkerArg,'Y')
                g1.set_names('x',LegendxAxis,'y',LegendyAxis1,'color',ColorName,'lightness',LightnessName,'marker',MarkerName)
            elseif strcmpi(MarkerArg,'N')
                g1.set_names('x',LegendxAxis,'y',LegendyAxis1,'color',ColorName,'lightness',LightnessName)
            end
            g1.set_color_options('n_color',NumColor,...
                'n_lightness',NumLightness,...
                'legend','separate_gray')
            g1.set_text_options('font','Arial','base_size',BaseFontSize,'label_scaling',LabelScaling,'legend_scaling',LegendScaling)
            if strcmpi(LegendArg,'Y')
                g1.set_layout_options("legend",1) % Don't show legend
            elseif strcmpi(MarkerArg,'N')
                g1.set_layout_options("legend",0) % Show legend
            end
            % Figure
            h_fig1=figure(1);
            h_fig1.Color='white'; % changes the background color of the figure
            h_fig1.Units='pixel'; % Defines the units
            h_fig1.OuterPosition=Res;
            h_fig1.PaperOrientation='landscape';
            h_fig1.Name=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix1);
            % The actual plotting
            g1.draw()
            % Save the figure
            FullName1=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix1);
            print(h_fig1,FullName1,'-dpng');
            exportgraphics(h_fig1,[FullName1,'.pdf'],'ContentType','vector');

            %% Gramm object 2
            % Define variables
            Plottitle2=sprintf('%d Force Maps containing %d Force Curves selected',length(obj.SMFSResults{ResultsRow,1}.Data(1).FMIndex),obj.SMFSResults{ResultsRow,1}.Data(1).SumNumFcAnalysedAdhMaxRet);
            LegendyAxis2='Adhesion force (nN)';
            NameSuffix2='_MaxAdhesionForceRetract';
            % Allocate data
            yData2=obj.SMFSResults{ResultsRow}.Concatenate.AdhMaxRet(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*-1e9;
            yData2Min=min(yData2);
            yData2Max=max(yData2);
            yData2diff=yData2Max-yData2Min;
            % Create a gramm object
            if strcmpi(MarkerArg,'Y')
                g2=gramm('x',xData,'y',yData2,...
                    'color',ColorData,...
                    'lightness',LightnessData,...
                    'marker',MarkerData);
            elseif strcmpi(MarkerArg,'N')
                g2=gramm('x',xData,'y',yData2,...
                    'color',ColorData,...
                    'lightness',LightnessData);
            end           
            % Plot data           
            g2.stat_boxplot('notch',true,...
                'width',BoxplotWidth,...
                'dodge',BoxplotDodge); % Plot data in boxplot
            % Set options
            g2.axe_property('LineWidth',LineWidth,'xlim',[xDataMin-xAxisCorr xDataMax+xAxisCorr],'ylim',[yData2Min-yData2diff*yAxisMinFactor yData2Max+yData2diff*yAxisMaxFactor],'XTick',[0:10:220],'TickDir','out','TickLength',[0.005 0.005]);
            if strcmpi(xArg,'DateTime')
                g2.set_datetick('x',0,'keeplimits') % Format x-axis
            end
            % Color Bar
            if strcmpi(CBar,'Y')
            y2CBar=[yData2Max+yData2diff*yLimMinFactor yData2Max+yData2diff*yLimMinFactor yData2Max+yData2diff*yLimMaxFactor yData2Max+yData2diff*yLimMaxFactor];
            g2.geom_polygon('x',{xCBar1;xCBar2;xCBar3;xCBar4},'y',{y2CBar;y2CBar;y2CBar;y2CBar},'color',ColorBarMap,'alpha',1);
            else
            end
            g2.set_point_options('base_size',MarkerSize)
            %g2.set_title(Plottitle2) %Set figure title
            if strcmpi(MarkerArg,'Y')
                g2.set_names('x',LegendxAxis,'y',LegendyAxis2,'color',ColorName,'lightness',LightnessName,'marker',MarkerName)
            elseif strcmpi(MarkerArg,'N')
                g2.set_names('x',LegendxAxis,'y',LegendyAxis2,'color',ColorName,'lightness',LightnessName)
            end
            g2.set_color_options('n_color',NumColor,...
                'n_lightness',NumLightness,...
                'legend','separate_gray')
            g2.set_text_options('font','Arial','base_size',BaseFontSize,'label_scaling',LabelScaling,'legend_scaling',LegendScaling)
            if strcmpi(LegendArg,'Y')
                g2.set_layout_options("legend",1) % Don't show legend
            elseif strcmpi(MarkerArg,'N')
                g2.set_layout_options("legend",0) % Show legend
            end
            % Figure
            h_fig2=figure(2);
            h_fig2.Color='white'; % changes the background color of the figure
            h_fig2.Units='pixel'; % Defines the units
            h_fig2.OuterPosition=Res;
            h_fig2.PaperOrientation='landscape';
            h_fig2.Name=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix2);
            % The actual plotting
            g2.draw()
            % Save the figure
            FullName2=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix2);
            print(h_fig2,FullName2,'-dpng');
            exportgraphics(h_fig2,[FullName2,'.pdf'],'ContentType','vector');

            %% Gramm object 3
            % Define variables
            Plottitle3=sprintf('%d Force Maps containing %d Force Curves selected',length(obj.SMFSResults{ResultsRow,1}.Data(1).FMIndex),obj.SMFSResults{ResultsRow,1}.Data(1).SumNumFcAnalysedAdhUnbinding);
            LegendyAxis3='Adhesion force (nN)';
            NameSuffix3='_AdhForceUnbinding';
            % Allocate data
            yData3=obj.SMFSResults{ResultsRow}.Concatenate.AdhUnbinding(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*-1e9;
            yData3Min=min(yData3);
            yData3Max=max(yData3);
            yData3diff=yData3Max-yData3Min;
            % Create a gramm object
            if strcmpi(MarkerArg,'Y')
                g3=gramm('x',xData,'y',yData3,...
                    'color',ColorData,...
                    'lightness',LightnessData,...
                    'marker',MarkerData);
            elseif strcmpi(MarkerArg,'N')
                g3=gramm('x',xData,'y',yData3,...
                    'color',ColorData,...
                    'lightness',LightnessData);
            end
            % Plot data
            g3.stat_boxplot('notch',true,...
                'width',BoxplotWidth,...
                'dodge',BoxplotDodge); % Plot data in boxplot
            % Set options
            g3.axe_property('LineWidth',LineWidth,'xlim',[xDataMin-xAxisCorr xDataMax+xAxisCorr],'ylim',[yData3Min-yData3diff*yAxisMinFactor yData3Max+yData3diff*yAxisMaxFactor],'XTick',[0:10:220],'TickDir','out','TickLength',[0.005 0.005]);
            if strcmpi(xArg,'DateTime')
                g3.set_datetick('x',0,'keeplimits') % Format x-axis
            end
            % Color Bar
            if strcmpi(CBar,'Y')
            y3CBar=[yData3Max+yData3diff*yLimMinFactor yData3Max+yData3diff*yLimMinFactor yData3Max+yData3diff*yLimMaxFactor yData3Max+yData3diff*yLimMaxFactor];
            g3.geom_polygon('x',{xCBar1;xCBar2;xCBar3;xCBar4},'y',{y3CBar;y3CBar;y3CBar;y3CBar},'color',ColorBarMap,'alpha',1);
            else
            end
            g3.set_point_options('base_size',MarkerSize)
            %g3.set_title(Plottitle3) %Set figure title
            if strcmpi(MarkerArg,'Y')
                g3.set_names('x',LegendxAxis,'y',LegendyAxis3,'color',ColorName,'lightness',LightnessName,'marker',MarkerName)
            elseif strcmpi(MarkerArg,'N')
                g3.set_names('x',LegendxAxis,'y',LegendyAxis3,'color',ColorName,'lightness',LightnessName)
            end
            g3.set_color_options('n_color',NumColor,...
                'n_lightness',NumLightness,...
                'legend','separate_gray')
            g3.set_text_options('font','Arial','base_size',BaseFontSize,'label_scaling',LabelScaling,'legend_scaling',LegendScaling)
            if strcmpi(LegendArg,'Y')
                g3.set_layout_options("legend",1) % Don't show legend
            elseif strcmpi(MarkerArg,'N')
                g3.set_layout_options("legend",0) % Show legend
            end
            % Figure
            h_fig3=figure(3);
            h_fig3.Color='white'; % changes the background color of the figure
            h_fig3.Units='pixel'; % Defines the units
            h_fig3.OuterPosition=Res;
            h_fig3.PaperOrientation='landscape';
            h_fig3.Name=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix3);
            % The actual plotting
            g3.draw()
            % Save the figure
            FullName3=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix3);
            print(h_fig3,FullName3,'-dpng');
            exportgraphics(h_fig3,[FullName3,'.pdf'],'ContentType','vector');

            %% Gramm object 4
            % Define variables
            Plottitle4=sprintf('%d Force Maps containing %d Force Curves selected',length(obj.SMFSResults{ResultsRow,1}.Data(1).FMIndex),obj.SMFSResults{ResultsRow,1}.Data(1).SumNumFcAnalysedAdhEneApp);
            LegendyAxis4='Work of adhesion (aJ)';
            NameSuffix4='_AdhEnergyApproach';
            % Allocate data
            yData4=obj.SMFSResults{ResultsRow}.Concatenate.AdhEneApp(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*-1e18;
            yData4Min=min(yData4);
            yData4Max=max(yData4);
            yData4diff=yData4Max-yData4Min;
            % Create a gramm object
            if strcmpi(MarkerArg,'Y')
                g4=gramm('x',xData,'y',yData4,...
                    'color',ColorData,...
                    'lightness',LightnessData,...
                    'marker',MarkerData);
            elseif strcmpi(MarkerArg,'N')
                g4=gramm('x',xData,'y',yData4,...
                    'color',ColorData,...
                    'lightness',LightnessData);
            end
            % Plot data            
            g4.stat_boxplot('notch',true,...
                'width',BoxplotWidth,...
                'dodge',BoxplotDodge); % Plot data in boxplot
            % Set options
            g4.axe_property('LineWidth',LineWidth,'xlim',[xDataMin-xAxisCorr xDataMax+xAxisCorr],'ylim',[yData4Min-yData4diff*yAxisMinFactor yData4Max+yData4diff*yAxisMaxFactor],'XTick',[0:10:220],'TickDir','out','TickLength',[0.005 0.005]);
            if strcmpi(xArg,'DateTime')
                g4.set_datetick('x',0,'keeplimits') % Format x-axis
            end
            % Color Bar
            if strcmpi(CBar,'Y')
            y4CBar=[yData4Max+yData4diff*yLimMinFactor yData4Max+yData4diff*yLimMinFactor yData4Max+yData4diff*yLimMaxFactor yData4Max+yData4diff*yLimMaxFactor];
            g4.geom_polygon('x',{xCBar1;xCBar2;xCBar3;xCBar4},'y',{y4CBar;y4CBar;y4CBar;y4CBar},'color',ColorBarMap,'alpha',1);
            else
            end
            g4.set_point_options('base_size',MarkerSize)
            %g4.set_title(Plottitle4) %Set figure title
            if strcmpi(MarkerArg,'Y')
                g4.set_names('x',LegendxAxis,'y',LegendyAxis4,'color',ColorName,'lightness',LightnessName,'marker',MarkerName)
            elseif strcmpi(MarkerArg,'N')
                g4.set_names('x',LegendxAxis,'y',LegendyAxis4,'color',ColorName,'lightness',LightnessName)
            end
            g4.set_color_options('n_color',NumColor,...
                'n_lightness',NumLightness,...
                'legend','separate_gray')
            g4.set_text_options('font','Arial','base_size',BaseFontSize,'label_scaling',LabelScaling,'legend_scaling',LegendScaling)
            if strcmpi(LegendArg,'Y')
                g4.set_layout_options("legend",1) % Don't show legend
            elseif strcmpi(MarkerArg,'N')
                g4.set_layout_options("legend",0) % Show legend
            end
            % Figure
            h_fig4=figure(4);
            h_fig4.Color='white'; % changes the background color of the figure
            h_fig4.Units='pixel'; % Defines the units
            h_fig4.OuterPosition=Res;
            h_fig4.PaperOrientation='landscape';
            h_fig4.Name=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix4);
            % The actual plotting
            g4.draw()
            % Save the figure
            FullName4=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix4);
            print(h_fig4,FullName4,'-dpng');
            exportgraphics(h_fig4,[FullName4,'.pdf'],'ContentType','vector');

            %% Gramm object 5
            % Define variables
            Plottitle5=sprintf('%d Force Maps containing %d Force Curves selected',length(obj.SMFSResults{ResultsRow,1}.Data(1).FMIndex),obj.SMFSResults{ResultsRow,1}.Data(1).SumNumFcAnalysedAdhEneRet);
            LegendyAxis5='Work of adhesion (aJ)';
            NameSuffix5='_AdhEnergyRetract';
            % Allocate data
            yData5=obj.SMFSResults{ResultsRow}.Concatenate.AdhEneRet(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*-1e18;
            yData5Min=min(yData5);
            yData5Max=max(yData5);
            yData5diff=yData5Max-yData5Min;
            % Create a gramm object
            if strcmpi(MarkerArg,'Y')
                g5=gramm('x',xData,'y',yData5,...
                    'color',ColorData,...
                    'lightness',LightnessData,...
                    'marker',MarkerData);
            elseif strcmpi(MarkerArg,'N')
                g5=gramm('x',xData,'y',yData5,...
                    'color',ColorData,...
                    'lightness',LightnessData);
            end
            % Plot data
            g5.stat_boxplot('notch',true,...
                'width',BoxplotWidth,...
                'dodge',BoxplotDodge); % Plot data in boxplot
            % Set options
            g5.axe_property('LineWidth',LineWidth,'xlim',[xDataMin-xAxisCorr xDataMax+xAxisCorr],'ylim',[yData5Min-yData5diff*yAxisMinFactor yData5Max+yData5diff*yAxisMaxFactor],'XTick',[0:10:220],'TickDir','out','TickLength',[0.005 0.005]);
            if strcmpi(xArg,'DateTime')
                g5.set_datetick('x',0,'keeplimits') % Format x-axis
            end
            % Color Bar
            if strcmpi(CBar,'Y')
            y5CBar=[yData5Max+yData5diff*yLimMinFactor yData5Max+yData5diff*yLimMinFactor yData5Max+yData5diff*yLimMaxFactor yData5Max+yData5diff*yLimMaxFactor];
            g5.geom_polygon('x',{xCBar1;xCBar2;xCBar3;xCBar4},'y',{y5CBar;y5CBar;y5CBar;y5CBar},'color',ColorBarMap,'alpha',1);
            else
            end
            g5.set_point_options('base_size',MarkerSize)
            %g5.set_title(Plottitle5) %Set figure title
            if strcmpi(MarkerArg,'Y')
                g5.set_names('x',LegendxAxis,'y',LegendyAxis5,'color',ColorName,'lightness',LightnessName,'marker',MarkerName)
            elseif strcmpi(MarkerArg,'N')
                g5.set_names('x',LegendxAxis,'y',LegendyAxis5,'color',ColorName,'lightness',LightnessName)
            end
            g5.set_color_options('n_color',NumColor,...
                'n_lightness',NumLightness,...
                'legend','separate_gray')
            g5.set_text_options('font','Arial','base_size',BaseFontSize,'label_scaling',LabelScaling,'legend_scaling',LegendScaling)
            if strcmpi(LegendArg,'Y')
                g5.set_layout_options("legend",1) % Don't show legend
            elseif strcmpi(MarkerArg,'N')
                g5.set_layout_options("legend",0) % Show legend
            end
            %g5.set_layout_options('legend_position',[0.75,0.4,0.35,0.6]) %[left bottom width height]
            % Figure
            h_fig5=figure(5);
            h_fig5.Color='white'; % changes the background color of the figure
            h_fig5.Units='pixel'; % Defines the units
            h_fig5.OuterPosition=Res;
            h_fig5.PaperOrientation='landscape';
            h_fig5.Name=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix5);
            % The actual plotting
            g5.draw()
            % Save the figure
            FullName5=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix5);
            print(h_fig5,FullName5,'-dpng');
            exportgraphics(h_fig5,[FullName5,'.pdf'],'ContentType','vector')

            %% Gramm object 6
            % Define variables
            Plottitle6=sprintf('%d Force Maps containing %d Force Curves selected',length(obj.SMFSResults{ResultsRow,1}.Data(1).FMIndex),obj.SMFSResults{ResultsRow,1}.Data(1).SumNumFcAnalysedyPullingLength);
            LegendyAxis6='Pull-off length (nm)';
           % LegendyAxis6='L_{Pull-off} (nm)';
            NameSuffix6='_Pullinglength';
            % Allocate data
            yData6=obj.SMFSResults{ResultsRow}.Concatenate.PullingLength(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*1e9;
            yData6Min=min(yData6);
            yData6Max=max(yData6);
            yData6diff=yData6Max-yData6Min;
            % Create a gramm object
            if strcmpi(MarkerArg,'Y')
                g6=gramm('x',xData,'y',yData6,...
                    'color',ColorData,...
                    'lightness',LightnessData,...
                    'marker',MarkerData);
            elseif strcmpi(MarkerArg,'N')
                g6=gramm('x',xData,'y',yData6,...
                    'color',ColorData,...
                    'lightness',LightnessData);
            end
            % Plot data
            g6.stat_boxplot('notch',true,...
                'width',BoxplotWidth,...
                'dodge',BoxplotDodge); % Plot data in boxplot
            % Set options
            g6.axe_property('LineWidth',LineWidth,'xlim',[xDataMin-xAxisCorr xDataMax+xAxisCorr],'ylim',[yData6Min-yData6diff*yAxisMinFactor yData6Max+yData6diff*yAxisMaxFactor],'XTick',[0:10:220],'TickDir','out','TickLength',[0.005 0.005]);
            % Reference length color box
            g6.geom_polygon('y',{LimitLengthRet2},'color',SteelBlue);    
            if strcmpi(xArg,'DateTime')
                g6.set_datetick('x',0,'keeplimits') % Format x-axis
            end
            % Color bar
            if strcmpi(CBar,'Y')
            y6CBar=[yData6Max+yData6diff*yLimMinFactor yData6Max+yData6diff*yLimMinFactor yData6Max+yData6diff*yLimMaxFactor yData6Max+yData6diff*yLimMaxFactor];
            g6.geom_polygon('x',{xCBar1;xCBar2;xCBar3;xCBar4},'y',{y6CBar;y6CBar;y6CBar;y6CBar},'color',ColorBarMap,'alpha',1);
            else
            end
            g6.set_point_options('base_size',MarkerSize)
            %g6.set_title(Plottitle6) %Set figure title
            if strcmpi(MarkerArg,'Y')
                g6.set_names('x',LegendxAxis,'y',LegendyAxis6,'color',ColorName,'lightness',LightnessName,'marker',MarkerName)
            elseif strcmpi(MarkerArg,'N')
                g6.set_names('x',LegendxAxis,'y',LegendyAxis6,'color',ColorName,'lightness',LightnessName)
            end
            g6.set_color_options('n_color',NumColor,...
                'n_lightness',NumLightness,...
                'legend','separate_gray')
            g6.set_text_options('font','Arial','base_size',BaseFontSize,'label_scaling',LabelScaling,'legend_scaling',LegendScaling)
            if strcmpi(LegendArg,'Y')
                g6.set_layout_options("legend",1) % Don't show legend
            elseif strcmpi(MarkerArg,'N')
                g6.set_layout_options("legend",0) % Show legend
            end
            % Figure
            h_fig6=figure(6);
            h_fig6.Color='white'; % changes the background color of the figure
            h_fig6.Units='pixel'; % Defines the units
            h_fig6.OuterPosition=Res;
            h_fig6.PaperOrientation='landscape';
            h_fig6.Name=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix6);
            % The actual plotting
            g6.draw()
            % Save the figure
            FullName6=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix6);
            print(h_fig6,FullName6,'-dpng');
            exportgraphics(h_fig6,[FullName6,'.pdf'],'ContentType','vector')
            
            %% Gramm object 7
            % Define variables
            Plottitle7=sprintf('%d Force Maps containing %d Force Curves selected',length(obj.SMFSResults{ResultsRow,1}.Data(1).FMIndex),obj.SMFSResults{ResultsRow,1}.Data(1).SumNumFcAnalysedySnapInLength);
            LegendyAxis7='Snap-In length (nm)';
            NameSuffix7='_SnapInLength';
            % Allocate data
            yData7=obj.SMFSResults{ResultsRow}.Concatenate.SnapInLength(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*1e9;
            yData7Min=min(yData7);
            yData7Max=max(yData7);
            yData7diff=yData7Max-yData7Min;
            % Create a gramm object
            if strcmpi(MarkerArg,'Y')
                g7=gramm('x',xData,'y',yData7,...
                    'color',ColorData,...
                    'lightness',LightnessData,...
                    'marker',MarkerData);
            elseif strcmpi(MarkerArg,'N')
                g7=gramm('x',xData,'y',yData7,...
                    'color',ColorData,...
                    'lightness',LightnessData);
            end
            % Plot data
            g7.stat_boxplot('notch',true,...
                'width',BoxplotWidth,...
                'dodge',BoxplotDodge); % Plot data in boxplot
            % Set options
            g7.axe_property('LineWidth',LineWidth,'xlim',[xDataMin-xAxisCorr xDataMax+xAxisCorr],'ylim',[yData7Min-yData7diff*yAxisMinFactor yData7Max+yData7diff*yAxisMaxFactor],'XTick',[0:10:220],'TickDir','out','TickLength',[0.005 0.005]);
            g7.geom_polygon('y',{LimitLengthApp},'color',Ochreish);
            if strcmpi(xArg,'DateTime')
                g7.set_datetick('x',0,'keeplimits') % Format x-axis
            end
            if strcmpi(CBar,'Y')
            y7CBar=[yData7Max+yData7diff*yLimMinFactor yData7Max+yData7diff*yLimMinFactor yData7Max+yData7diff*yLimMaxFactor yData7Max+yData7diff*yLimMaxFactor];
            g7.geom_polygon('x',{xCBar1;xCBar2;xCBar3;xCBar4},'y',{y7CBar;y7CBar;y7CBar;y7CBar},'color',ColorBarMap,'alpha',1);
            else
            end
            g7.set_point_options('base_size',MarkerSize)
            %g7.set_title(Plottitle6) %Set figure title
            if strcmpi(MarkerArg,'Y')
                g7.set_names('x',LegendxAxis,'y',LegendyAxis7,'color',ColorName,'lightness',LightnessName,'marker',MarkerName)
            elseif strcmpi(MarkerArg,'N')
                g7.set_names('x',LegendxAxis,'y',LegendyAxis7,'color',ColorName,'lightness',LightnessName)
            end
            g7.set_color_options('n_color',NumColor,...
                'n_lightness',NumLightness,...
                'legend','separate_gray')
            g7.set_text_options('font','Arial','base_size',BaseFontSize,'label_scaling',LabelScaling,'legend_scaling',LegendScaling)
            if strcmpi(LegendArg,'Y')
                g7.set_layout_options("legend",1) % Don't show legend
            elseif strcmpi(MarkerArg,'N')
                g7.set_layout_options("legend",0) % Show legend
            end
            % Figure
            h_fig7=figure(7);
            h_fig7.Color='white'; % changes the background color of the figure
            h_fig7.Units='pixel'; % Defines the units
            h_fig7.OuterPosition=Res;
            h_fig7.PaperOrientation='landscape';
            h_fig7.Name=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix7);
            % The actual plotting
            g7.draw()
            % Save the figure
            FullName7=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix7);
            print(h_fig7,FullName7,'-dpng');
            exportgraphics(h_fig7,[FullName7,'.pdf'],'ContentType','vector')

            % House keeping
            close all
        end

       function SMFS_results_gramm_plot2_publication(obj,ResultsRow,Linker,xArg,MarkerArg,LegendArg,CBar,Var,yDataArg)
           % Author: Andreas Rohatschek
           % Description:
           % This function allows to plot the analysis results using the
           % gramm toolbox (https://joss.theoj.org/papers/10.21105/joss.00568). The
           % different features (snap-in length, pull-off length) determined are plotted as data points to
           % the corresponding force map.
           % Required input variables:
           % ResultsRow ... row in the results structure: double
           % Linker ... Linker type: string , either 'long' or 'short'
           % xArg ... x-axis argument: string, either 'Index' or 'DateTime'
           % MarkerArg ... Marker argument: string, either 'Y' or 'N'
           % LegendArg ... Legend argument: string, either 'Y' or 'N'
           % CBar ... Color Bar: string, either 'Y' or 'N'
           % Var ... Variant: double, Either 1, 2 or 3
           % yDataArg ... y-data argument (To be able to use the fuction also for trials which are not completely analysable): string, either 'old' or 'new' 

           %% Function body
           % Define variables
           yLimMaxFactor=0.25;
           yLimMinFactor=0.2;
           %yLimMaxFactor=0.3;
           %yLimMinFactor=0.25;
           yAxisMinFactor=0.09;
           yAxisMaxFactor=0.25;
           Res=[1 1 2560 1250]; % Define the figure resolution
           %MarkerSize=10;
           MarkerSize=20;
           BaseFontSize=46;
           LabelScaling=1.2;
           LegendScaling=1;
           LineWidth=1.5;
           FontName='Arial';
            % Define color bar
            xCBar1=[0 134 134 0]; % Native
            xCBar2=[134 2100 2100 134]; % Sliding
            xCBar3=[2200 12100 12100 2200]; % Unraveled
            xCBar4=[12100 22200 22200 12100]; % Dissociated
            % Color and Color maps
            Ochreish=[253 174 97]./255;
            SteelBlue=[116 173 209]./255;
           ColorBrewerMap1=[[253 174 97]./255; % Ochreish
               [116 173 209]./255]; % Steel blueish
           ColorBarMap=[[54 163 0]./255; % Dark green HEX 8DB600
               [206 22 32]./255; % Fire Engine Red HEX CE162
               [0 24 204]./255; % Blue HEX 8DB600
               [135 0 224]./255]; % Violet HEX 8F00FF
%             ColorMarkerMap=[[255 194 191]./255; % Light rose
%                [255 194 191]./255;  % Light rose
%                [209 217 161]./255;  % Light green
%                [209 217 161]./255; % Light green
%                [0 136 255]./255; % Light blue
%                [0 136 255]./255; % Light blue
%                [200 0 200]./255; % Pink
%                [200 0 200]./255]; % Pink
            ColorMarkerMap='parula'; 
           % Change into the Folder of Interest
           cd(obj.ExperimentFolder) % Move into the folder
           % Create folders for saving the produced figures
           foldername='SMFS_results_gramm_plot2_publication';    % Defines the folder name
           mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
           currpath=fullfile(obj.ExperimentFolder,foldername);
           cd(currpath);
           %% Input variables
           % Linker
           if strcmpi(Linker,'Long')
               LimitLengthRet1=[0 378];
               LimitLengthRet2=[378 522];
               LimitLengthApp=[50 120];
           elseif strcmpi(Linker,'Short')
               LimitLengthRet1=[0 308];
               LimitLengthRet2=[308 463];
               LimitLengthApp=[50 120];
           end
           % xArg
           if strcmpi(xArg,'DateTime')
               LegendxAxis='Date and Time';
               xData=obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSort;
               xDataMin=min(xData);
               xDataMax=max(xData);
               xAxisCorr=(xDataMax-xDataMin)*0.05;
           elseif strcmpi(xArg,'Index')
               LegendxAxis='Number of cycles';
               xData=obj.SMFSResults{ResultsRow}.Concatenate.FcNum;
               xDataMin=min(xData);
               xDataMax=max(xData);
               xAxisCorr=(xDataMax-xDataMin)*0.005;
           end
           % Var
           if Var==1
               ColorName='Medium';
               LightnessName='Substrate';
               MarkerName='ChipCantilever';
               LightnessData=obj.SMFSResults{ResultsRow}.Concatenate.FMSubstrate(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx);
               ColorData=obj.SMFSResults{ResultsRow}.Concatenate.FMEnvCond(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx);
               MarkerData=obj.SMFSResults{ResultsRow}.Concatenate.FMChipCant(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx);
               NumColor=2;
               NumLightness=2;
           elseif Var==2
               ColorName='Approach speed ($\mu$m/s)';
               LightnessName='Retraction speed ($\mu$m/s)';
               MarkerName='Dwell Time (s)';
               FMExtVeloData=obj.SMFSResults{ResultsRow}.Concatenate.FMExtVelocity(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*1e6;
               FMRetVeloData=obj.SMFSResults{ResultsRow}.Concatenate.FMRetVelocity(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*1e6;
               ColorData=FMExtVeloData;
               LightnessData=FMRetVeloData;
               MarkerData=obj.SMFSResults{ResultsRow}.Concatenate.FMHoldingTime(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx);
               NumColor=6;
               NumLightness=6;
           elseif Var==3
                ColorName='Dwell Time (s)';
                LightnessName='Retraction velocity ($\mu$m/s)';
                MarkerName='Approach velocity ($\mu$m/s)';
                FMExtVeloData=obj.SMFSResults{ResultsRow}.Concatenate.FMExtVelocity(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*1e6;
                FMRetVeloData=obj.SMFSResults{ResultsRow}.Concatenate.FMRetVelocity(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*1e6;
                ColorData=obj.SMFSResults{ResultsRow}.Concatenate.FMHoldingTime(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx);
                LightnessData=FMRetVeloData;
                MarkerData=FMExtVeloData;
                NumColor=2;
                NumLightness=2;
           end
           % Define parameter term for the figure name 
           if obj.SMFSResults{ResultsRow}.Parameters.ExtendVelocity==0
               ExtVelocityValueStr='All';
           else
               ExtVelocityValueStr=num2str(round(obj.SMFSResults{ResultsRow}.Parameters.ExtendVelocity*1e9));
           end
           if obj.SMFSResults{ResultsRow}.Parameters.RetractVelocity==0
               RetVelocityValueStr='All';
           else
               RetVelocityValueStr=num2str(round(obj.SMFSResults{ResultsRow}.Parameters.RetractVelocity*1e9));
           end
           if obj.SMFSResults{ResultsRow}.Parameters.HoldingTime==-1
               HoldingTimeValueStr='All';
           else
               HoldingTimeValueStr=num2str(obj.SMFSResults{ResultsRow}.Parameters.HoldingTime);
           end
           % General names
           FigNamePt1=sprintf('SMFSResultRow%d_',ResultsRow);
           FigNamePt2=strcat(ExtVelocityValueStr,{'_'},RetVelocityValueStr,{'_'},HoldingTimeValueStr,{'_'},obj.SMFSResults{ResultsRow}.Parameters.Substrate,{'_'},obj.SMFSResults{ResultsRow}.Parameters.Medium,{'_'},obj.SMFSResults{ResultsRow}.Parameters.ChipCantilever,{'_'},obj.SMFSResults{ResultsRow}.Parameters.Chipbox,{'_'},obj.SMFSResults{ResultsRow}.Parameters.Linker);
           FigNamePt2=char(FigNamePt2);
           FigNamePt3='_Plot2';
           %% Gramm object 1
           % Define variables
           Plottitle1=sprintf('%d Force Maps containing %d Force Curves selected',length(obj.SMFSResults{ResultsRow,1}.Data(1).FMIndex),obj.SMFSResults{ResultsRow,1}.Data(1).SumNumFcAnalysedAdhMaxApp);
           LegendyAxis1='Adhesion force (nN)';
           NameSuffix1='_MaxAdhesionForceApproach';
           % Allocate data
           if strcmpi(yDataArg,'new')
           yData1=obj.SMFSResults{ResultsRow}.Concatenate.AdhMaxApp(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*-1e9;
           elseif strcmpi(yDataArg,'old')
           yData1=obj.SMFSResults{ResultsRow}.Data.AdhMaxAppConcat(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*-1;   
           end
           yData1Min=min(yData1);
           yData1Max=max(yData1);
           yData1diff=yData1Max-yData1Min;
           % Create a gramm object
           if strcmpi(MarkerArg,'Y')
               g1=gramm('x',xData,'y',yData1,...
                   'color',ColorData,...
                   'lightness',LightnessData,...
                   'marker',MarkerData);
           elseif strcmpi(MarkerArg,'N')
               g1=gramm('x',xData,'y',yData1,...
                   'color',ColorData,...
                   'lightness',LightnessData);
           end
           % Plot data
           g1.geom_point();
           % Set options
           g1.axe_property('LineWidth',LineWidth,'xlim',[xDataMin-xAxisCorr xDataMax+xAxisCorr],'ylim',[yData1Min-yData1diff*yAxisMinFactor yData1Max+yData1diff*yAxisMaxFactor],'TickDir','out','TickLength',[0.005 0.005]);
           if strcmpi(xArg,'DateTime')
               g1.set_datetick('x',0,'keeplimits') % Format x-axis
           end
           % Color Bar
           if strcmpi(CBar,'Y')
           y1CBar=[yData1Max+yData1diff*yLimMinFactor yData1Max+yData1diff*yLimMinFactor yData1Max+yData1diff*yLimMaxFactor yData1Max+yData1diff*yLimMaxFactor];
           g1.geom_polygon('x',{xCBar1;xCBar2;xCBar3;xCBar4},'y',{y1CBar;y1CBar;y1CBar;y1CBar},'color',ColorBarMap,'alpha',1);
           else
           end
           g1.set_point_options('base_size',MarkerSize)
           %g1.set_title(Plottitle1) %Set figure title
           if strcmpi(MarkerArg,'Y')
               g1.set_names('x',LegendxAxis,'y',LegendyAxis1,'color',ColorName,'lightness',LightnessName,'marker',MarkerName)
           elseif strcmpi(MarkerArg,'N')
               g1.set_names('x',LegendxAxis,'y',LegendyAxis1,'color',ColorName,'lightness',LightnessName)
           end
           g1.set_color_options('Map',ColorMarkerMap,...
               'n_color',NumColor,...
               'n_lightness',NumLightness,...
               'legend','separate_gray')
           g1.set_text_options('font',FontName,'base_size',BaseFontSize,'label_scaling',LabelScaling,'legend_scaling',LegendScaling)
           if strcmpi(LegendArg,'Y')
               g1.set_layout_options("legend",1) % Show legend
           elseif strcmpi(LegendArg,'N')
               g1.set_layout_options("legend",0) % Don't show legend
           end
           % Figure
           h_fig1=figure(1);
           h_fig1.Color='white'; % changes the background color of the figure
           h_fig1.Units='pixel'; % Defines the units
           h_fig1.OuterPosition=Res;
           h_fig1.PaperOrientation='landscape';
            h_fig1.Name=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix1);
            % The actual plotting
            g1.draw()
            % Save the figure
            FullName1=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix1);
            print(h_fig1,FullName1,'-dpng');
            exportgraphics(h_fig1,[FullName1,'.pdf'],'ContentType','vector')

           %% Gramm object 2
           % Define variables
           Plottitle2=sprintf('%d Force Maps containing %d Force Curves selected',length(obj.SMFSResults{ResultsRow,1}.Data(1).FMIndex),obj.SMFSResults{ResultsRow,1}.Data(1).SumNumFcAnalysedAdhMaxRet);
           LegendyAxis2='Adhesion force (nN)';
           NameSuffix2='_MaxAdhesionForceRetract';
           % Allocate data
           if strcmpi(yDataArg,'new')
           yData2=obj.SMFSResults{ResultsRow}.Concatenate.AdhMaxRet(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*-1e9;
           elseif strcmpi(yDataArg,'old')
            yData2=obj.SMFSResults{ResultsRow}.Data.AdhMaxRetConcat(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*-1e9;   
           end     
           yData2Min=min(yData2);
           yData2Max=max(yData2);
           yData2diff=yData2Max-yData2Min;
           % Create a gramm object
           if strcmpi(MarkerArg,'Y')
               g2=gramm('x',xData,'y',yData2,...
                   'color',ColorData,...
                   'lightness',LightnessData,...
                   'marker',MarkerData);
           elseif strcmpi(MarkerArg,'N')
               g2=gramm('x',xData,'y',yData2,...
                   'color',ColorData,...
                   'lightness',LightnessData);
           end
           % Plot data
           g2.geom_point();
           % Set options
            g2.axe_property('LineWidth',LineWidth,'xlim',[xDataMin-xAxisCorr xDataMax+xAxisCorr],'ylim',[yData2Min-yData2diff*yAxisMinFactor yData2Max+yData2diff*yAxisMaxFactor],'TickDir','out','TickLength',[0.005 0.005]);
            if strcmpi(xArg,'DateTime')
                g2.set_datetick('x',0,'keeplimits') % Format x-axis
            end
            % Color Bar
            if strcmpi(CBar,'Y')
            y2CBar=[yData2Max+yData2diff*yLimMinFactor yData2Max+yData2diff*yLimMinFactor yData2Max+yData2diff*yLimMaxFactor yData2Max+yData2diff*yLimMaxFactor];
            g2.geom_polygon('x',{xCBar1;xCBar2;xCBar3;xCBar4},'y',{y2CBar;y2CBar;y2CBar;y2CBar},'color',ColorBarMap,'alpha',1);
            else
            end
           g2.set_point_options('base_size',MarkerSize)
           %g2.set_title(Plottitle2) %Set figure title
           if strcmpi(MarkerArg,'Y')
               g2.set_names('x',LegendxAxis,'y',LegendyAxis2,'color',ColorName,'lightness',LightnessName,'marker',MarkerName)
           elseif strcmpi(MarkerArg,'N')
               g2.set_names('x',LegendxAxis,'y',LegendyAxis2,'color',ColorName,'lightness',LightnessName)
           end
           g2.set_color_options('Map',ColorMarkerMap,...
               'n_color',NumColor,...
               'n_lightness',NumLightness,...
               'legend','separate_gray')
           g2.set_text_options('font',FontName,'base_size',BaseFontSize,'label_scaling',LabelScaling,'legend_scaling',LegendScaling)
           if strcmpi(LegendArg,'Y')
               g2.set_layout_options("legend",1) % Show legend
           elseif strcmpi(LegendArg,'N')
               g2.set_layout_options("legend",0) % Don't show legend
           end
           % Figure
           h_fig2=figure(2);
           h_fig2.Color='white'; % changes the background color of the figure
           h_fig2.Units='pixel'; % Defines the units
           h_fig2.OuterPosition=Res;
           h_fig2.PaperOrientation='landscape';          
            h_fig2.Name=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix2);
            % The actual plotting
            g2.draw()
            % Save the figure
            FullName2=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix2);
            print(h_fig2,FullName2,'-dpng');
            exportgraphics(h_fig2,[FullName2,'.pdf'],'ContentType','vector')

           %% Gramm object 3
           % Define variables
           Plottitle3=sprintf('%d Force Maps containing %d Force Curves selected',length(obj.SMFSResults{ResultsRow,1}.Data(1).FMIndex),obj.SMFSResults{ResultsRow,1}.Data(1).SumNumFcAnalysedAdhUnbinding);
           LegendyAxis3='Adhesion force (nN)';
           NameSuffix3='_AdhForceUnbinding';
           % Allocate data
           if strcmpi(yDataArg,'new')
           yData3=obj.SMFSResults{ResultsRow}.Concatenate.AdhUnbinding(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*-1e9;
           elseif strcmpi(yDataArg,'old')
           yData3=obj.SMFSResults{ResultsRow}.Data.AdhUnbindingConcat(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*-1e9;    
           end     
           yData3Min=min(yData3);
           yData3Max=max(yData3);
           yData3diff=yData3Max-yData3Min;
           % Create a gramm object
           if strcmpi(MarkerArg,'Y')
               g3=gramm('x',xData,'y',yData3,...
                   'color',ColorData,...
                   'lightness',LightnessData,...
                   'marker',MarkerData);
           elseif strcmpi(MarkerArg,'N')
               g3=gramm('x',xData,'y',yData3,...
                   'color',ColorData,...
                   'lightness',LightnessData);
           end
           % Plot data
           g3.geom_point()
            % Set options
            g3.axe_property('LineWidth',LineWidth,'xlim',[xDataMin-xAxisCorr xDataMax+xAxisCorr],'ylim',[yData3Min-yData3diff*yAxisMinFactor yData3Max+yData3diff*yAxisMaxFactor],'TickDir','out','TickLength',[0.005 0.005]);
            if strcmpi(xArg,'DateTime')
                g3.set_datetick('x',0,'keeplimits') % Format x-axis
            end
            % Color Bar
            if strcmpi(CBar,'Y')
            y3CBar=[yData3Max+yData3diff*yLimMinFactor yData3Max+yData3diff*yLimMinFactor yData3Max+yData3diff*yLimMaxFactor yData3Max+yData3diff*yLimMaxFactor];
            g3.geom_polygon('x',{xCBar1;xCBar2;xCBar3;xCBar4},'y',{y3CBar;y3CBar;y3CBar;y3CBar},'color',ColorBarMap,'alpha',1);
            else
            end
           g3.set_point_options('base_size',MarkerSize)
           %g3.set_title(Plottitle3) %Set figure title
           if strcmpi(MarkerArg,'Y')
               g3.set_names('x',LegendxAxis,'y',LegendyAxis3,'color',ColorName,'lightness',LightnessName,'marker',MarkerName)
           elseif strcmpi(MarkerArg,'N')
               g3.set_names('x',LegendxAxis,'y',LegendyAxis3,'color',ColorName,'lightness',LightnessName)
           end
           g3.set_color_options('Map',ColorMarkerMap,...
               'n_color',NumColor,...
               'n_lightness',NumLightness,...
               'legend','separate_gray')
           g3.set_text_options('font',FontName,'base_size',BaseFontSize,'label_scaling',LabelScaling,'legend_scaling',LegendScaling)
           if strcmpi(LegendArg,'Y')
               g3.set_layout_options("legend",1) % Show legend
           elseif strcmpi(LegendArg,'N')
               g3.set_layout_options("legend",0) % Don't show legend
           end
           % Figure
           h_fig3=figure(3);
           h_fig3.Color='white'; % changes the background color of the figure
           h_fig3.Units='pixel'; % Defines the units
           h_fig3.OuterPosition=Res;
           h_fig3.PaperOrientation='landscape';          
           h_fig3.Name=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix3);
           % The actual plotting
           g3.draw()
           % Save the figure
           FullName3=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix3);
           print(h_fig3,FullName3,'-dpng');
           exportgraphics(h_fig3,[FullName3,'.pdf'],'ContentType','vector')

           %% Gramm object 4
           % Define variables
           Plottitle4=sprintf('%d Force Maps containing %d Force Curves selected',length(obj.SMFSResults{ResultsRow,1}.Data(1).FMIndex),obj.SMFSResults{ResultsRow,1}.Data(1).SumNumFcAnalysedAdhEneApp);
           LegendyAxis4='Work of adhesion (aJ)';
           NameSuffix4='_AdhEnergyApproach';
           % Allocate data
           if strcmpi(yDataArg,'new')
           yData4=obj.SMFSResults{ResultsRow}.Concatenate.AdhEneApp(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*-1e18;
           elseif strcmpi(yDataArg,'old')
           yData4=obj.SMFSResults{ResultsRow}.Data.AdhEneAppConcat(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*-1e18;     
           end     
           yData4Min=min(yData4);
           yData4Max=max(yData4);
           yData4diff=yData4Max-yData4Min;
           % Create a gramm object
           if strcmpi(MarkerArg,'Y')
               g4=gramm('x',xData,'y',yData4,...
                   'color',ColorData,...
                   'lightness',LightnessData,...
                   'marker',MarkerData);
           elseif strcmpi(MarkerArg,'N')
               g4=gramm('x',xData,'y',yData4,...
                   'color',ColorData,...
                   'lightness',LightnessData);
           end
           % Plot data
           g4.geom_point()
            % Set options
            g4.axe_property('LineWidth',LineWidth,'xlim',[xDataMin-xAxisCorr xDataMax+xAxisCorr],'ylim',[yData4Min-yData4diff*yAxisMinFactor yData4Max+yData4diff*yAxisMaxFactor],'TickDir','out','TickLength',[0.005 0.005]);
            if strcmpi(xArg,'DateTime')
                g4.set_datetick('x',0,'keeplimits') % Format x-axis
            end
            % Color Bar
            if strcmpi(CBar,'Y')
            y4CBar=[yData4Max+yData4diff*yLimMinFactor yData4Max+yData4diff*yLimMinFactor yData4Max+yData4diff*yLimMaxFactor yData4Max+yData4diff*yLimMaxFactor];
            g4.geom_polygon('x',{xCBar1;xCBar2;xCBar3;xCBar4},'y',{y4CBar;y4CBar;y4CBar;y4CBar},'color',ColorBarMap,'alpha',1);
            else
            end
           g4.set_point_options('base_size',MarkerSize)
           %g4.set_title(Plottitle4) %Set figure title
           if strcmpi(MarkerArg,'Y')
               g4.set_names('x',LegendxAxis,'y',LegendyAxis4,'color',ColorName,'lightness',LightnessName,'marker',MarkerName)
           elseif strcmpi(MarkerArg,'N')
               g4.set_names('x',LegendxAxis,'y',LegendyAxis4,'color',ColorName,'lightness',LightnessName)
           end
           g4.set_color_options('Map',ColorMarkerMap,...
               'n_color',NumColor,...
               'n_lightness',NumLightness,...
               'legend','separate_gray')
           g4.set_text_options('font',FontName,'base_size',BaseFontSize,'label_scaling',LabelScaling,'legend_scaling',LegendScaling)
           if strcmpi(LegendArg,'Y')
               g4.set_layout_options("legend",1) % Show legend
           elseif strcmpi(LegendArg,'N')
               g4.set_layout_options("legend",0) % Don't show legend
           end
           % Figure
           h_fig4=figure(4);
           h_fig4.Color='white'; % changes the background color of the figure
           h_fig4.Units='pixel'; % Defines the units
           h_fig4.OuterPosition=Res;
           h_fig4.PaperOrientation='landscape';
            h_fig4.Name=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix4);
            % The actual plotting
            g4.draw()
            % Save the figure
            FullName4=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix4);
            print(h_fig4,FullName4,'-dpng');
            exportgraphics(h_fig4,[FullName4,'.pdf'],'ContentType','vector')

           %% Gramm object 5
           % Define variables
           Plottitle5=sprintf('%d Force Maps containing %d Force Curves selected',length(obj.SMFSResults{ResultsRow,1}.Data(1).FMIndex),obj.SMFSResults{ResultsRow,1}.Data(1).SumNumFcAnalysedAdhEneRet);
           LegendyAxis5='Work of adhesion (aJ)';
           NameSuffix5='_AdhEnergyRetract';
           % Allocate data
           if strcmpi(yDataArg,'new')
           yData5=obj.SMFSResults{ResultsRow}.Concatenate.AdhEneRet(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*-1e18;
           elseif strcmpi(yDataArg,'old')
           yData5=obj.SMFSResults{ResultsRow}.Data.AdhEneRetConcat(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*-1e18;     
           end     
           yData5Min=min(yData5);
           yData5Max=max(yData5);
           yData5diff=yData5Max-yData5Min;
           % Create a gramm object
           if strcmpi(MarkerArg,'Y')
               g5=gramm('x',xData,'y',yData5,...
                   'color',ColorData,...
                   'lightness',LightnessData,...
                   'marker',MarkerData);
           elseif strcmpi(MarkerArg,'N')
               g5=gramm('x',xData,'y',yData5,...
                   'color',ColorData,...
                   'lightness',LightnessData);
           end
           % Plot data
           g5.geom_point()
            % Set options
            g5.axe_property('LineWidth',LineWidth,'xlim',[xDataMin-xAxisCorr xDataMax+xAxisCorr],'ylim',[yData5Min-yData5diff*yAxisMinFactor yData5Max+yData5diff*yAxisMaxFactor],'TickDir','out','TickLength',[0.005 0.005]);
            if strcmpi(xArg,'DateTime')
                g5.set_datetick('x',0,'keeplimits') % Format x-axis
            end
            % Color Bar
            if strcmpi(CBar,'Y')
            y5CBar=[yData5Max+yData5diff*yLimMinFactor yData5Max+yData5diff*yLimMinFactor yData5Max+yData5diff*yLimMaxFactor yData5Max+yData5diff*yLimMaxFactor];
            g5.geom_polygon('x',{xCBar1;xCBar2;xCBar3;xCBar4},'y',{y5CBar;y5CBar;y5CBar;y5CBar},'color',ColorBarMap,'alpha',1);
            else
            end
           g5.set_point_options('base_size',MarkerSize)
           %g5.set_title(Plottitle5) %Set figure title
           if strcmpi(MarkerArg,'Y')
               g5.set_names('x',LegendxAxis,'y',LegendyAxis5,'color',ColorName,'lightness',LightnessName,'marker',MarkerName)
           elseif strcmpi(MarkerArg,'N')
               g5.set_names('x',LegendxAxis,'y',LegendyAxis5,'color',ColorName,'lightness',LightnessName)
           end
           g5.set_color_options('Map',ColorMarkerMap,...
               'n_color',NumColor,...
               'n_lightness',NumLightness,...
               'legend','separate_gray')
           g5.set_text_options('font',FontName,'base_size',BaseFontSize,'label_scaling',LabelScaling,'legend_scaling',LegendScaling)
           if strcmpi(LegendArg,'Y')
               g5.set_layout_options("legend",1) % Show legend
           elseif strcmpi(LegendArg,'N')
               g5.set_layout_options("legend",0) % Don't show legend           
           end
           % Figure
           h_fig5=figure(5);
           h_fig5.Color='white'; % changes the background color of the figure
           h_fig5.Units='pixel'; % Defines the units
           h_fig5.OuterPosition=Res;
           h_fig5.PaperOrientation='landscape';
            h_fig5.Name=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix5);
            % The actual plotting
            g5.draw()
            % Save the figure
            FullName5=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix5);
            print(h_fig5,FullName5,'-dpng');
            exportgraphics(h_fig5,[FullName5,'.pdf'],'ContentType','vector')

           %% Gramm object 6
           % Define variables
           Plottitle6=sprintf('%d Force Maps containing %d Force Curves selected',length(obj.SMFSResults{ResultsRow,1}.Data(1).FMIndex),obj.SMFSResults{ResultsRow,1}.Data(1).SumNumFcAnalysedyPullingLength);
           LegendyAxis6='Pull-off length (nm)';
           NameSuffix6='_Pullinglength';
           % Allocate data
           if strcmpi(yDataArg,'new')
           yData6=obj.SMFSResults{ResultsRow}.Concatenate.PullingLength(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*1e9;
           elseif strcmpi(yDataArg,'old')
           yData6=obj.SMFSResults{ResultsRow}.Data.yPullingLengthConcat(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*1e9;     
           end           
           yData6Min=min(yData6);
           yData6Max=max(yData6);
           yData6diff=yData6Max-yData6Min;
           % Create a gramm object
           if strcmpi(MarkerArg,'Y')
               g6=gramm('x',xData,'y',yData6,...
                   'color',ColorData,...
                   'lightness',LightnessData,...
                   'marker',MarkerData);
           elseif strcmpi(MarkerArg,'N')
               g6=gramm('x',xData,'y',yData6,...
                   'color',ColorData,...
                   'lightness',LightnessData);
           end
           % Plot data
           g6.geom_point();
            % Set options
            g6.axe_property('LineWidth',LineWidth,'xlim',[xDataMin-xAxisCorr xDataMax+xAxisCorr],'ylim',[yData6Min-yData6diff*yAxisMinFactor yData6Max+yData6diff*yAxisMaxFactor],'TickDir','out','TickLength',[0.005 0.005]);
            % Reference length color box
            g6.geom_polygon('y',{LimitLengthRet2},'color',SteelBlue);    
            if strcmpi(xArg,'DateTime')
                g6.set_datetick('x',0,'keeplimits') % Format x-axis
            end
            % Color bar
            if strcmpi(CBar,'Y')
            y6CBar=[yData6Max+yData6diff*yLimMinFactor yData6Max+yData6diff*yLimMinFactor yData6Max+yData6diff*yLimMaxFactor yData6Max+yData6diff*yLimMaxFactor];
            g6.geom_polygon('x',{xCBar1;xCBar2;xCBar3;xCBar4},'y',{y6CBar;y6CBar;y6CBar;y6CBar},'color',ColorBarMap,'alpha',1);
            else
            end
           g6.set_point_options('base_size',MarkerSize)
           %g6.set_title(Plottitle6) %Set figure title
           if strcmpi(MarkerArg,'Y')
               g6.set_names('x',LegendxAxis,'y',LegendyAxis6,'color',ColorName,'lightness',LightnessName,'marker',MarkerName)
           elseif strcmpi(MarkerArg,'N')
               g6.set_names('x',LegendxAxis,'y',LegendyAxis6,'color',ColorName,'lightness',LightnessName)
           end
           g6.set_color_options('Map',ColorMarkerMap,...
               'n_color',NumColor,...
               'n_lightness',NumLightness,...
               'legend','separate_gray')
           g6.set_text_options('font',FontName,'base_size',BaseFontSize,'label_scaling',LabelScaling,'legend_scaling',LegendScaling)
           if strcmpi(LegendArg,'Y')
               g6.set_layout_options("legend",1) % Show legend
           elseif strcmpi(LegendArg,'N')
               g6.set_layout_options("legend",0) % Don't show legend
           end
           %g6.set_layout_options('legend_position',[0.75,0.4,0.35,0.6]) %[left bottom width height]
           % Figure
           h_fig6=figure(6);
           h_fig6.Color='white'; % changes the background color of the figure
           h_fig6.Units='pixel'; % Defines the units
           h_fig6.OuterPosition=Res;
           h_fig6.PaperOrientation='landscape';
            h_fig6.Name=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix6);
            % The actual plotting
            g6.draw()
            % Save the figure
            FullName6=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix6);
            print(h_fig6,FullName6,'-dpng');
            exportgraphics(h_fig6,[FullName6,'.pdf'],'ContentType','vector')

           %% Gramm object 7
           % Define variables
           if strcmpi(yDataArg,'new')     
           Plottitle7=sprintf('%d Force Maps containing %d Force Curves selected',length(obj.SMFSResults{ResultsRow,1}.Data(1).FMIndex),obj.SMFSResults{ResultsRow,1}.Data(1).SumNumFcAnalysedySnapInLength);
           LegendyAxis7='Snap-In length (nm)';
           NameSuffix7='_SnapInLength';
           % Allocate data
           yData7=obj.SMFSResults{ResultsRow}.Concatenate.SnapInLength(obj.SMFSResults{ResultsRow}.Concatenate.FMDateTimeNumberSortIdx)*1e9;
           yData7Min=min(yData7);
           yData7Max=max(yData7);
           yData7diff=yData7Max-yData7Min;
           % Create a gramm object
           if strcmpi(MarkerArg,'Y')
               g7=gramm('x',xData,'y',yData7,...
                   'color',ColorData,...
                   'lightness',LightnessData,...
                   'marker',MarkerData);
           elseif strcmpi(MarkerArg,'N')
               g7=gramm('x',xData,'y',yData7,...
                   'color',ColorData,...
                   'lightness',LightnessData);
           end
           % Plot data
           g7.geom_polygon('y',{LimitLengthApp},'color',Ochreish);
           g7.geom_point();
            % Set options
            g7.axe_property('LineWidth',LineWidth,'xlim',[xDataMin-xAxisCorr xDataMax+xAxisCorr],'ylim',[yData7Min-yData7diff*yAxisMinFactor yData7Max+yData7diff*yAxisMaxFactor],'TickDir','out','TickLength',[0.005 0.005]);
            g7.geom_polygon('y',{LimitLengthApp},'color',Ochreish);
            if strcmpi(xArg,'DateTime')
                g7.set_datetick('x',0,'keeplimits') % Format x-axis
            end
            if strcmpi(CBar,'Y')
            y7CBar=[yData7Max+yData7diff*yLimMinFactor yData7Max+yData7diff*yLimMinFactor yData7Max+yData7diff*yLimMaxFactor yData7Max+yData7diff*yLimMaxFactor];
            g7.geom_polygon('x',{xCBar1;xCBar2;xCBar3;xCBar4},'y',{y7CBar;y7CBar;y7CBar;y7CBar},'color',ColorBarMap,'alpha',1);
            else
            end
           g7.set_point_options('base_size',MarkerSize)
           %g7.set_title(Plottitle7) %Set figure title
           if strcmpi(MarkerArg,'Y')
               g7.set_names('x',LegendxAxis,'y',LegendyAxis7,'color',ColorName,'lightness',LightnessName,'marker',MarkerName)
           elseif strcmpi(MarkerArg,'N')
               g7.set_names('x',LegendxAxis,'y',LegendyAxis7,'color',ColorName,'lightness',LightnessName)
           end
           g7.set_color_options('Map',ColorMarkerMap,...
               'n_color',NumColor,...
               'n_lightness',NumLightness,...
               'legend','separate_gray')
           g7.set_text_options('font',FontName,'base_size',BaseFontSize,'label_scaling',LabelScaling,'legend_scaling',LegendScaling)
           if strcmpi(LegendArg,'Y')
               g7.set_layout_options("legend",1) % Show legend
           elseif strcmpi(LegendArg,'N')
               g7.set_layout_options("legend",0) % Don't show legend
           end
           % Figure
           h_fig7=figure(7);
           h_fig7.Color='white'; % changes the background color of the figure
           h_fig7.Units='pixel'; % Defines the units
           h_fig7.OuterPosition=Res;
           h_fig7.PaperOrientation='landscape';
            h_fig7.Name=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix7);
            % The actual plotting
            g7.draw()
            % Save the figure
            FullName7=strcat(FigNamePt1,FigNamePt2,FigNamePt3,'_',xArg,'_','Var',num2str(Var),NameSuffix7);
            print(h_fig7,FullName7,'-dpng');
            exportgraphics(h_fig7,[FullName7,'.pdf'],'ContentType','vector')
          
           elseif strcmpi(yDataArg,'old')
          
           end      
           % House keeping
           close all
       end


    
       function SMFS_force_landscape_app(obj,Var)
            % Author: Manuel Rufin, Andreas Rohatschek
            % Description: 
            % In this function a force landscape of the approach part is plotted. This function is based on "2024_11_14_AveragingCorrection.m" file in the "FCLandscape"-folder and slightly modified.
            % Description: 
            % In this function  
            % Required input variables:
            % Var ... Variant: double, either 1 or 2
            % Typical input variables:
            % Var ... 1

            %% Function body
            if Var==1
                % Input dialog
                prompt = {'Enter the force map number you do not want to have included in the "SMFSResults"-structure (For multiple selections just use the space key to separeat entries)'};
                definput = {''};
                opts.Interpreter = 'tex';
                FMIdxArray=inputdlg(prompt,'Select all - except of ...',[1 150],definput,opts);
                FMIdxArray=str2num(FMIdxArray{1}); % Convert the cell array to numerals
            elseif Var==2
                warning('This variation requires that a none parameter based selection row in the "SMFSResults"-structure exists (Typically row 1)')
                dlgtitle='Enter the "SMFSResults"-structure row the chronologcial data will be taken from';
                prompt = {'ResultsRow'};
                definput={'1'};
                dims=[1 150];
                ResultsRow = inputdlg(prompt,dlgtitle,dims,definput);
                ResultsRow=str2num(ResultsRow{1});
                dlgtitle='Enter the chronological force map index of the first and last force map you want to have included in the "SMFSResults"-structure';
                prompt = {'First','Last'};
                dims=[1 150; 1 150];
                FMIdxChrono = inputdlg(prompt,dlgtitle,dims);
                FMIdxChrono1=str2num(FMIdxChrono{1});
                FMIdxChrono2=str2num(FMIdxChrono{2});
                % Allocate the force maps in chronological order as defined by
                % the input dialog box
                FMIdxArray=obj.SMFSResults{ResultsRow}.Data.FMDateTimeSortIdx(FMIdxChrono1:FMIdxChrono2);
            end
            % If condition to handle an empty index array
            if isempty(FMIdxArray)
                return
            else
            end

            % Allocate the force maps in chronological order as defined by
            % the input dialog box
            FMIdxArray=obj.SMFSResults{ResultsRow}.Data.FMDateTimeSortIdx(FMIdxChrono1:FMIdxChrono2);

            % Change into the Folder of Interest
            cd(obj.ExperimentFolder) % Move into the folder
            % Create folders for saving the produced figures
            foldername='SMFS_force_landscape';    % Defines the folder name
            mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
            currpath=fullfile(obj.ExperimentFolder,foldername);
            cd(currpath);

           %% Figures
           FigNamePt1=sprintf('SMFSResultRow%d_',ResultsRow);
           FigNamePt2=sprintf('FM%d_to_FM%d_',FMIdxChrono1,FMIdxChrono2);
           FigNamePt3='App';
           FigNamePt4='ForceCurveLandscape';

            % Define variables
            FcCount = 0; % Force curve count
            for ii=FMIdxArray'
                FcCount = FcCount + obj.FM{ii}.NCurves;
            end
            ResY = FcCount; % y-axis resolution
            ResX = 512; % x-axis resolution
            CritLength=200*1e-9; % Max. length on x-axis
            % Allocate data
            MaxRangeApp = 0;
            SkippedCurves = 0;
            k=1;
            DwellTime = [];
            Phase = [];
            ExtSpeed = [];
            RetSpeed = [];
            AbsCycleIndex = [];
            m = 0;
            % Read out data
            for j=obj.SMFSResults{ResultsRow}.Data.FMIndexChrono'
                for i=1:obj.FM{j}.NCurves
                    m=m+1;
                     if ~obj.FM{j}.SMFSFlag.Selected(i)
                         SkippedCurves = SkippedCurves + 1;
                         continue
                     end
                    DwellTime(k,1) = obj.SMFSResults{ResultsRow}.Data.FMHoldingTime(j);
                    AbsCycleIndex(k,1) = m;
                    % Assigning the phase/state
                    if any(j == [1:2])
                        Phase(k,1) = 1;
                    elseif any(j == [3:22])
                        Phase(k,1) = 2;
                    elseif any(j == [23:121])
                        Phase(k,1) = 3;
                    elseif any(j == [121:obj.NumForceMaps])
                        Phase(k,1) = 4;
                    end
                    % Allocating data
                    vDefApp{k} = obj.FM{j}.BasedApp{i};
                    THApp{k} = obj.FM{j}.THApp{i};
                    THApp{k} = THApp{k} - max(THApp{k});
                    THApp{k}=THApp{k}*-1;
                    AppIdx=find(THApp{k}>CritLength,1,'last'); % Find all entries longer than the defined length
                    THApp{k}(1:AppIdx)=[];
                    vDefApp{k}(1:AppIdx)=[];

                    MaxRangeApp = max(range(THApp{k}),MaxRangeApp);

                    k = k + 1;
                end
                    if j == FMIdxChrono2 % Leave loop when the entry FM is reached
                        break
                    end   
            end

            % Define axes properties
            xmax=200*1e-9;
            XQApp = linspace(0,xmax,ResX);
            ResY = ResY - SkippedCurves; % y-axis resolution
            % Preallocate
            FcMapApp = zeros(ResY, ResX);
            % Interp Data and fill FCMap
            for i=1:ResY
                FCMapApp(i,:) = interp1(THApp{i},vDefApp{i},XQApp);
            end
            ApproachForceCurveMap = FCMapApp;

            % Font Sizes
            FS = 18;

            % Processed Image Preparation
            ProcAppMap = ApproachForceCurveMap * 1e9;
            YPixels = 512;
            XPixels = 512;
            %img = imresize(-ProcRetMap, [YPixels XPixels], 'bilinear');
            img=-ProcAppMap;
            img = imfilter(img,ones(round(size(-ProcAppMap,1)/YPixels),1)./round(size(ProcAppMap,1)/YPixels)); % Applying a rollling average filter for smearing out data
            img = imresize(img, [YPixels XPixels], 'nearest'); % Downsampling
            ProcDwellTime = imresize(DwellTime, [YPixels 1], 'nearest');
            ProcPhase = imresize(Phase, [YPixels 1], 'nearest');
            xmin = 0;
            xmax = 200;
            ymax = size(ApproachForceCurveMap, 1) / 1000;
            MaxForce = 0.2;
            MinForce = 0;
            LineWidth = 1.5;

            % Figure and Image Display
            Fig = figure('Color', 'w');
            imshow(img, [MinForce MaxForce]);
            Fig.Units='normalized';
            % Axis
            ax = gca;  % Get the current axes handle
            axis on;   % Turn on the axis visibility
            xlabel('Tip-sample distance (nm)');
            ylabel('Number of cycles (x100)');
            xlim([0 XPixels]);  % Set X-axis limits
            ylim([0 YPixels]);  % Set Y-axis limits
            set(ax, 'XTick', [0:50:200] ./ xmax * XPixels);  % Custom X-axis ticks
            YTickSpacing = [0:0.5:floor(ymax)] ./ ymax;
            set(ax, 'YTick', YTickSpacing * YPixels);  % Custom Y-axis ticks
            set(ax, 'XTickLabel', {string(xmin), string(50), string(100), string(150), string(200)});  % Custom X-axis tick labels
            % Custom Y-axis tick labels
            for i=1:length(YTickSpacing)
                if i == 1
                    YTickLabels{i} = '0';
                    continue
                end
               YTickLabels{i} =string(round(AbsCycleIndex(YTickSpacing(i)*ymax*1000)/100,1));
             %   YTickLabels{i} = string(round(AbsCycleIndex(YTickSpacing(i)*ymax*1000),1));
            end
            set(ax, 'YTickLabel',YTickLabels)              
            set(ax, 'FontSize', FS);  % Set font size for axes
            set(ax, 'FontName', 'Arial');  % Set font type for axes
            % Display axes only on the left and bottom
            ax.XAxisLocation = 'bottom';  % X-axis at the bottom
            ax.YAxisLocation = 'left';    % Y-axis on the left
            ax.Box = 'off';  % Turn off the box around the axes
            ax.TickDir = 'out';  % Ticks pointing outwards
            % Ensure the top and right axes are not shown
            ax.XColor = 'k';  % X-axis color (bottom)
            ax.YColor = 'k';  % Y-axis color (left)
            ax.XRuler.Axle.Visible = 'off';  % Hide the top axis
            ax.YRuler.Axle.Visible = 'off';  % Hide the right axis
            ax.LineWidth = LineWidth;
            % Adjust Axes Position to Make Space for Rectangles
            ax.Position = [0.15 0.11 0.7 0.8];  % Adjust the axis position
            LeftAxExpand = 0.07;

            %% Adding Colored Rectangles
            % Left-hand side rectangles for dwell time colored gray 
            ScalingParam = 0.18*LeftAxExpand;
            UpperPos = 0;
            for i=1:YPixels
                DwellSize = ProcDwellTime(i) + 0.4;
                Thickness = ScalingParam * DwellSize * XPixels;
                Height = 1;
                rectangle('Position', [-Thickness, UpperPos,...
                    Thickness, Height],...
                    'FaceColor', 0.5.*[1 1 1], 'EdgeColor', 'none');
                UpperPos = UpperPos + Height;
            end

            % Add Bar for cycle axis
            rectangle('Position', [-LeftAxExpand * XPixels, 0,...
                0.01 * XPixels, YPixels],...
                'FaceColor', [0 0 0], 'EdgeColor', 'none');

            ColorBarMap=[[54 163 0]./255; % Dark green HEX 8DB600
                [206 22 32]./255; % Fire Engine Red HEX CE162
                [0 24 204]./255; % Blue HEX 8DB600
                [135 0 224]./255]; % Violet HEX 8F00FF

            UpperPos = 0;
            for i=1:YPixels
                Color = ColorBarMap(ProcPhase(i),:);
                Thickness = 0.05 * XPixels;
                Height = 1;
                rectangle('Position', [XPixels, UpperPos,...
                    Thickness, Height],...
                    'FaceColor', Color, 'EdgeColor', 'none');
                UpperPos = UpperPos + Height;
            end

            % Adjust the axis limits to avoid cutting off the rectangles
            xlim([-LeftAxExpand * XPixels XPixels * 1.05]);

            % Adding and Customizing Colorbar
            c = colorbar;
            c.Location = 'northoutside';
            c.Label.String = 'Attractive force (nN)';
            c.Ticks = [MinForce, MaxForce/3, MaxForce*2/3, MaxForce];
            c.FontSize = FS;
            c.FontName = 'Arial';
            c.Orientation = 'horizontal';
            c.LineWidth = LineWidth;
            colormap(turbo);

            % Final Adjustments
            c.Position = [0.1939    0.8650    0.6239    0.030];

            % Save figure
            fullname=sprintf('%s%s%s%s',FigNamePt1,FigNamePt2,FigNamePt3,FigNamePt4);
            print(gcf,fullname,'-dpng');
            exportgraphics(gcf,[fullname,'.pdf'],'ContentType','vector')

            %% House keeping
            close all

        end

       function SMFS_force_landscape_ret(obj,Var)
            % Author: Manuel Rufin, Andreas Rohatschek
            % Description: 
            % In this function a force landscape of the retraction part is plotted. This function is based on "2024_11_14_AveragingCorrection.m" file in the "FCLandscape"-folder and slightly modified. 
            % Required input variables:
            % Var ... Variant: double, either 1 or 2
            % Typical input variables:
            % Var ... 1

            %% Function body
            if Var==1
                % Input dialog
                prompt = {'Enter the force map number you do not want to have included in the "SMFSResults"-structure (For multiple selections just use the space key to separeat entries)'};
                definput = {''};
                opts.Interpreter = 'tex';
                FMIdxArray=inputdlg(prompt,'Select all - except of ...',[1 150],definput,opts);
                FMIdxArray=str2num(FMIdxArray{1}); % Convert the cell array to numerals
            elseif Var==2
                warning('This variation requires that a none parameter based selection row in the "SMFSResults"-structure exists (Typically row 1)')
                dlgtitle='Enter the "SMFSResults"-structure row the chronologcial data will be taken from';
                prompt = {'ResultsRow'};
                definput={'1'};
                dims=[1 150];
                ResultsRow = inputdlg(prompt,dlgtitle,dims,definput);
                ResultsRow=str2num(ResultsRow{1});
                dlgtitle='Enter the chronological force map index of the first and last force map you want to have included in the "SMFSResults"-structure';
                prompt = {'First','Last'};
                dims=[1 150; 1 150];
                FMIdxChrono = inputdlg(prompt,dlgtitle,dims);
                FMIdxChrono1=str2num(FMIdxChrono{1});
                FMIdxChrono2=str2num(FMIdxChrono{2});
                % Allocate the force maps in chronological order as defined by
                % the input dialog box
                FMIdxArray=obj.SMFSResults{ResultsRow}.Data.FMDateTimeSortIdx(FMIdxChrono1:FMIdxChrono2);
            end
            % If condition to handle an empty index array
            if isempty(FMIdxArray)
                return
            else
            end

            % Allocate the force maps in chronological order as defined by
            % the input dialog box
            FMIdxArray=obj.SMFSResults{ResultsRow}.Data.FMDateTimeSortIdx(FMIdxChrono1:FMIdxChrono2);

            % Change into the Folder of Interest
            cd(obj.ExperimentFolder) % Move into the folder
            % Create folders for saving the produced figures
            foldername='SMFS_force_landscape';    % Defines the folder name
            mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
            currpath=fullfile(obj.ExperimentFolder,foldername);
            cd(currpath);
 
           %% Figures
           FigNamePt1=sprintf('SMFSResultRow%d_',ResultsRow);
           FigNamePt2=sprintf('FM%d_to_FM%d_',FMIdxChrono1,FMIdxChrono2);
           FigNamePt3='Ret';
           FigNamePt4='ForceCurveLandscape';

            % Force Landscape variables
            FcCount = 0; % Force curve count
            for ii=FMIdxArray'
                FcCount = FcCount + obj.FM{ii}.NCurves;
            end
            ResY = FcCount; % y-axis resolution
            ResX = 1024; % x-axis resolution
            % Allocate data
            MaxRangeRet = 0;
            SkippedCurves = 0;
            k=1;
            DwellTime = [];
            Phase = [];
            ExtSpeed = [];
            RetSpeed = [];
            AbsCycleIndex = [];
            m = 0;
            % Read ou data
            for j=obj.SMFSResults{ResultsRow}.Data.FMIndexChrono'
                for i=1:obj.FM{j}.NCurves
                    m=m+1;
                    if ~obj.FM{j}.SMFSFlag.Selected(i)
                        SkippedCurves = SkippedCurves + 1;
                        continue
                    end

                    DwellTime(k,1) = obj.SMFSResults{ResultsRow}.Data.FMHoldingTime(j);
                    ExtSpeed(k,1) = obj.SMFSResults{ResultsRow}.Data.FMExtVelocity(j);
                    RetSpeed(k,1) = obj.SMFSResults{ResultsRow}.Data.FMRetVelocity(j);
                    AbsCycleIndex(k,1) = m;

                    if any(j == [1:2])
                        Phase(k,1) = 1;
                    elseif any(j == [3:22])
                        Phase(k,1) = 2;
                    elseif any(j == [23:121])
                        Phase(k,1) = 3;
                    elseif any(j == [121:obj.NumForceMaps])
                        Phase(k,1) = 4;
                    end

                    vDefRet{k} = obj.FM{j}.BasedRet{i};
                    THRet{k} = obj.FM{j}.THRet{i};
                    THRet{k} = THRet{k} - max(THRet{k});
                    MaxRangeRet = max(range(THRet{k}),MaxRangeRet);

                    k = k + 1;
                end
            end

            % Define axes properties
            XQRet = linspace(0,-MaxRangeRet,ResX); % Query points on x-axis
            ResY = ResY - SkippedCurves; % y-axis resolution
            % Preallocate
            FcMapRet = zeros(ResY, ResX);
            % Interp Data and fill FCMap
            for i=1:ResY
                FCMapRet(i,:) = interp1(THRet{i},vDefRet{i},XQRet);
            end
            RetractForceCurveMap = FCMapRet;
            RetractMaxRange = MaxRangeRet;

            %% Plot: Dwell Time only            
            % Font Sizes
            %FS = 18;
            FS = 7;
           
            % Processed Image Preparation
            ProcRetMap = RetractForceCurveMap * 1e9;
            YPixels = 2048;
            XPixels = 1024;
            img=-ProcRetMap;
            img = imfilter(img,ones(round(size(-ProcRetMap,1)/YPixels),1)./round(size(ProcRetMap,1)/YPixels)); % Applying a rollling average filter for smearing out data
            img = imresize(img, [YPixels XPixels], 'nearest'); % Downsampling
            ProcDwellTime = imresize(DwellTime, [YPixels 1], 'nearest');
            ProcPhase = imresize(Phase, [YPixels 1], 'nearest');
            xmin = 0;
            xmax = RetractMaxRange * 1e9;
            ymax = size(RetractForceCurveMap, 1) / 1000;
            MaxForce = 0.3;
            MinForce = 0;         
            LineWidth = 1.5;
            
            % Figure and Image Display
            Fig = figure('Color', 'w');
            imshow(img, [MinForce MaxForce]);
            
            ax = gca;  % Get the current axes handle
            axis on;   % Turn on the axis visibility
            xlabel('Tip-sample distance (nm)');
            ylabel('Number of cycles (x100)');
            xlim([0 XPixels]);  % Set X-axis limits
            ylim([0 YPixels]);  % Set Y-axis limits         
            set(ax, 'XTick', [0:200:1000] ./ xmax * XPixels);  % Custom X-axis ticks
            YTickSpacing = [0:2:floor(ymax)] ./ ymax;
            set(ax, 'YTick', YTickSpacing * YPixels);  % Custom Y-axis ticks
            set(ax, 'XTickLabel', {string(xmin), string(200), string(400), string(600), string(800), string(1000)});  % Custom X-axis tick labels
            % Custom Y-axis tick labels
            for i=1:length(YTickSpacing)
                if i == 1
                    YTickLabels{i} = '0';
                    continue
                end
               YTickLabels{i} =string(round(AbsCycleIndex(YTickSpacing(i)*ymax*1000)/100,1));
             %   YTickLabels{i} = string(round(AbsCycleIndex(YTickSpacing(i)*ymax*1000),1));
            end
            set(ax, 'YTickLabel',YTickLabels)           
            set(ax, 'FontSize', FS);  % Set font size for axes
            set(ax, 'FontName', 'Arial');  % Set font type for axes
            
            % Display axes only on the left and bottom
            ax.XAxisLocation = 'bottom';  % X-axis at the bottom
            ax.YAxisLocation = 'left';    % Y-axis on the left
            ax.Box = 'off';  % Turn off the box around the axes
            ax.TickDir = 'out';  % Ticks pointing outwards
            
            % Ensure the top and right axes are not shown
            ax.XColor = 'k';  % X-axis color (bottom)
            ax.YColor = 'k';  % Y-axis color (left)
            ax.XRuler.Axle.Visible = 'off';  % Hide the top axis
            ax.YRuler.Axle.Visible = 'off';  % Hide the right axis
            ax.LineWidth = LineWidth;
            
            % Adjust Axes Position to Make Space for Rectangles
            ax.Position = [0.15 0.11 0.7 0.8];  % Adjust the axis position
            
            % Lines and boxes demarking critical lengths
            LinesFlag = true;
            if LinesFlag
                XPos1 = 250;
                LineX1 = ones(1,2)*XPos1 ./ xmax * XPixels;
                LineY1 = [.735 1] * YPixels;
                hold on
                plot(LineX1,LineY1,'w--','LineWidth',LineWidth)
                % text(LineX1(1) + 10,mean(LineY1),string(XPos1),'Color','w','FontName','Arial','FontSize',FS)
            
                % XPos2 = 430;
                % LineX2 = ones(1,2)*XPos2 ./ xmax * XPixels;
                % LineY2 = [.55 .68] * YPixels;
                % plot(LineX2,LineY2,'w--','LineWidth',LineWidth)
                % text(LineX2(1) + 10,mean(LineY2),string(XPos2),'Color','w','FontName','Arial','FontSize',FS)
            
                XPos3 = 670;
                LineX3 = ones(1,2)*XPos3 ./ xmax * XPixels;
                % LineY3 = [.29 .543] * YPixels;
                % plot(LineX3,LineY3,'w--','LineWidth',LineWidth)
                % text(LineX3(1) + 10,mean(LineY3),string(XPos3),'Color','w','FontName','Arial','FontSize',FS)
                
                XPos4 = 820;
                LineX4 = ones(1,2)*XPos4 ./ xmax * XPixels;
                % LineY4 = [.29 .543] * YPixels;
                % plot(LineX4,LineY4,'w--','LineWidth',LineWidth)
                % text(LineX4(1) + 10,mean(LineY4),string(XPos4),'Color','w','FontName','Arial','FontSize',FS)          
            
                Band1X1 = 670 ./ xmax * XPixels;
                Band1X2 = 820 ./ xmax * XPixels;
                Band1Pos = [Band1X1 .29*YPixels Band1X2-Band1X1 (.543 - .29)*YPixels];
                Band1Rect =rectangle('Position', Band1Pos,...
                    'FaceColor', 'none', 'EdgeColor', 'w','LineStyle','--','LineWidth',LineWidth);
                       
                Band2X1 = 250 ./ xmax * XPixels;
                Band2X2 = 580 ./ xmax * XPixels;
                Band2Pos = [Band2X1 .56*YPixels Band2X2-Band2X1 (.68 - .56)*YPixels];
                Band2Rect =rectangle('Position', Band2Pos,...
                    'FaceColor', 'none', 'EdgeColor', 'w','LineStyle','--','LineWidth',LineWidth);
            
                %     text(Band2X1(1) + 10,mean([.56 .68].*YPixels),string(250),'Color','w','FontName','Arial','FontSize',FS)
                %     text(Band2X2(1) + 10,mean([.56 .68].*YPixels),string(580),'Color','w','FontName','Arial','FontSize',FS)
            end
            
            % White band (demarking theoretical linker-molecule complex
            % range)
            BandFlag = true;
            if BandFlag
                BandX1 = 308 ./ xmax * XPixels;
                BandX2 = 463 ./ xmax * XPixels;
                BandPos = [BandX1 -2 BandX2-BandX1 YPixels];
                BandRect =rectangle('Position', BandPos,...
                    'FaceColor', [1 1 1 .3], 'EdgeColor', 'none');           
            
            end
            LeftAxExpand = 0.07;
            
            % Adding Colored Rectangles
            % Dwell time (Gray rectangles on the left-hand side with thicknesses according to the dwell times)            
            DwellSizes = [0.2 1.2 2.2 3.2];
            ScalingParam = 0.18*LeftAxExpand;
            UpperPos = 0;
            for i=1:YPixels
                DwellSize = ProcDwellTime(i) + 0.2;
                Thickness = ScalingParam * DwellSize * XPixels;
                Height = 1;
                rectangle('Position', [-Thickness, UpperPos,...
                    Thickness, Height],...
                    'FaceColor', 0.5.*[1 1 1], 'EdgeColor', 'none');
                UpperPos = UpperPos + Height;
            end
            
            % Add Bar for cycle axis
            rectangle('Position', [-LeftAxExpand * XPixels, 0,...
                0.01 * XPixels, YPixels],...
                'FaceColor', [0 0 0], 'EdgeColor', 'none');
            
            % State/phase rectangles on the right-hand side         
            ColorBarMap=[[54 163 0]./255; % Dark green HEX 8DB600
                [206 22 32]./255; % Fire Engine Red HEX CE162
                [0 24 204]./255; % Blue HEX 8DB600
                [135 0 224]./255]; % Violet HEX 8F00FF
            
            UpperPos = 0;
            for i=1:YPixels
                Color = ColorBarMap(ProcPhase(i),:);
                Thickness = 0.05 * XPixels;
                Height = 1;
                rectangle('Position', [XPixels, UpperPos,...
                    Thickness, Height],...
                    'FaceColor', Color, 'EdgeColor', 'none');
                UpperPos = UpperPos + Height;
            end
            
            % Adjust the axis limits to avoid cutting off the rectangles
            xlim([-LeftAxExpand * XPixels XPixels * 1.05]);
            
            % Adding and Customizing Colorbar
            c = colorbar;
            c.Location = 'northoutside';
            c.Label.String = 'Adhesive force (nN)';
            c.Ticks = [MinForce, MaxForce/3, MaxForce*2/3, MaxForce];
            c.FontSize = FS;
            c.FontName = 'Arial';
            c.Orientation = 'horizontal';
            c.LineWidth = LineWidth;
            colormap(turbo);
            
            % Final Adjustments            
            Fig.Position = [1858 68 686 1277];
            c.Position = [0.1939    0.8446    0.6239    0.0167];            
            set(Fig, 'Renderer', 'Painters');
            
            % Save figure
            fullname=sprintf('%s%s%s%s',FigNamePt1,FigNamePt2,FigNamePt3,FigNamePt4);
            print(gcf,fullname,'-dpng');
            exportgraphics(gcf,[fullname,'.pdf'],'ContentType','vector')

            %% House keeping
            close all

        end


       function SMFS_analysis_flag_status(obj)
            % Author: Andreas Rohatschek
            % Description: 
            % The function analysis flags relevant for SMFS analysis and gives their status. In detail the number of flags not set (down flags) which correspong to not processed force maps, is determined and given          
            
            %% Function body
            obj.SMFSFlagDown.SelectFM=find(~obj.SMFSFlag.SelectFM);
            obj.SMFSFlagDown.PropertiesParameters=find(~obj.SMFSFlag.PropertiesParameters);
            obj.SMFSFlagDown.Preprocessed=find(~obj.SMFSFlag.Preprocessed);
            obj.SMFSFlagDown.Presorted=find(~obj.SMFSFlag.Presorted);
            obj.SMFSFlagDown.NumForceCurves=find(~obj.SMFSFlag.NumForceCurves);            
            obj.SMFSFlagDown.AnalysedPreSelected=find(~obj.SMFSFlag.AnalysedPreSelected);
            obj.SMFSFlagDown.AnalysedPostSelected=find(~obj.SMFSFlag.AnalysedPostSelected);
     %       obj.SMFSFlagDown.Analysed=find(~obj.SMFSFlag.Analysed);
            for Fm=1:obj.NumForceMaps
            %for Fm=122:133
            obj.FM{Fm}.fc_flag_status          
            end
        end
  
      
        end
   
    
    methods
        % Methods for data visualization spanning all the data
       
        function show_image(obj)
            % TODO: implement ui elements for customization
            
            h.ColorMode(1).Background = 'k';
            h.ColorMode(1).Profile1 = [219 21 223]./255; %[189 0 96]./255; % 'b';
            h.ColorMode(1).Profile2 = 'c';
            h.ColorMode(1).Text = 'w';
            
            
            h.ColorMode(2).Background = 'w';
            h.ColorMode(2).Profile1 = [219 21 223]./255; %[94 170 170]./255; % [189 0 96]./255; %'b';
            h.ColorMode(2).Profile2 = [0 181 26]./255; % alternatives %[80 200 204]./255;%[0,0.870588235294118,0.407843137254902];
            h.ColorMode(2).Text = 'k';
            
            h.ColorIndex = 1;
            h.ReferenceFontSize = 24;
            h.ProfileLineWidth = 3;
            
            h.Fig = figure('Name',sprintf('%s',obj.ExperimentName),...
                'Units','pixels',...
                'Position',[200 200 1024 512],...
                'Color',h.ColorMode(h.ColorIndex).Background);
            
            h.B(1) = uicontrol('style','togglebutton',...
                'String','Cross Section',...
                'units','normalized',...
                'position',[.85 .45 .1 .05],...
                'Callback',@cross_section_toggle);
            
            [ClassPopUp,ClassIndex] = obj.string_of_existing_class_instances();
            h.NumClasses = length(ClassPopUp);
            Class{1} = obj.get_class_instance(ClassIndex(1,:));
            Class{2} = obj.get_class_instance(ClassIndex(1,:));
            PopUp = Class{1}.string_of_existing();
            
            h.B(4) = uicontrol('style','text',...
                'String','Channel 1',...
                'units','normalized',...
                'position',[.85 .95 .15 .04]);
            
            
            h.B(16) = uicontrol('style','popupmenu',...
                'String',ClassPopUp,...
                'units','normalized',...
                'position',[.85 .9 .15 .05],...
                'Callback',@draw_channel_1);
            
            h.B(2) = uicontrol('style','popupmenu',...
                'String',PopUp,...
                'units','normalized',...
                'position',[.85 .85 .15 .05],...
                'Callback',@draw_channel_1);
            
            h.B(5) = uicontrol('style','text',...
                'String','Channel 2',...
                'units','normalized',...
                'position',[.85 .7 .15 .04]);
            
            h.B(17) = uicontrol('style','popupmenu',...
                'String',ClassPopUp,...
                'units','normalized',...
                'position',[.85 .65 .15 .05],...
                'Callback',@draw_channel_2);
            
            h.B(3) = uicontrol('style','popupmenu',...
                'String',PopUp,...
                'units','normalized',...
                'position',[.85 .6 .15 .05],...
                'Callback',@draw_channel_2);
            
            h.B(6) = uicontrol('style','pushbutton',...
                'String','Save Figure',...
                'units','normalized',...
                'position',[.85 .1 .1 .05],...
                'Callback',@save_figure_to_file);
            
            h.B(7) = uicontrol('style','checkbox',...
                'String','...with white background',...
                'units','normalized',...
                'position',[.85 .05 .1 .04],...
                'Callback',@changed_color);
            
            h.B(8) = uicontrol('style','slider',...
                'Value',1,...
                'Units','normalized',...
                'Position',[.85 .83 .1 .02],...
                'Callback',@changed_slider);
            
            h.B(9) = uicontrol('style','slider',...
                'Value',0,...
                'Units','normalized',...
                'Position',[.85 .81 .1 .02],...
                'Callback',@changed_slider);
            
            h.B(10) = uicontrol('style','slider',...
                'Value',1,...
                'Units','normalized',...
                'Position',[.85 .58 .1 .02],...
                'Callback',@changed_slider);
            
            h.B(11) = uicontrol('style','slider',...
                'Value',0,...
                'Units','normalized',...
                'Position',[.85 .56 .1 .02],...
                'Callback',@changed_slider);
            
            h.B(12) = uicontrol('style','text',...
                'String','Max',...
                'Units','normalized',...
                'Position',[.95 .83 .03 .02]);
            
            h.B(13) = uicontrol('style','text',...
                'String','Min',...
                'Units','normalized',...
                'Position',[.95 .81 .03 .02]);
            
            h.B(14) = uicontrol('style','text',...
                'String','Max',...
                'Units','normalized',...
                'Position',[.95 .58 .03 .02]);
            
            h.B(15) = uicontrol('style','text',...
                'String','Min',...
                'Units','normalized',...
                'Position',[.95 .56 .03 .02]);
            
            h.B(18) = uicontrol('Style','checkbox',...
                'String','Both Channels',...
                'Value',0,...
                'Tooltip','Green, if both Channels have the same size scaling',...
                'Units','normalized',...
                'Position',[.85 .42 .1 .03],...
                'Callback',@checked_both_cross_sections);
            
            h.B(19) = uicontrol('style','checkbox',...
                'String','Upscale Images',...
                'units','normalized',...
                'position',[.85 .15 .1 .04],...
                'Callback',@upscale_images);
            
            h.B(20) = uicontrol('style','checkbox',...
                'String','Lock Channels',...
                'units','normalized',...
                'position',[.85 .3   .1 .04],...
                'Callback',@lock_channels);
            
            h.B(21) = uicontrol('style','checkbox',...
                'String','Lock Scalebars',...
                'units','normalized',...
                'position',[.85 .20 .1 .04],...
                'Callback',@lock_scalebars);
            
            h.B(22) = uicontrol('style','checkbox',...
                'String','Statistical CMapping',...
                'units','normalized',...
                'position',[.85 .25 .1 .04],...
                'Callback',@statistical_cmapping);
            
            h.B(23) = uicontrol('style','edit',...
                'String','',...
                'units','normalized',...
                'position',[.85 .75 .035 .05],...
                'Callback',@set_scale);
            
            h.B(24) = uicontrol('style','edit',...
                'String','',...
                'units','normalized',...
                'position',[.925 .75 .035 .05],...
                'Callback',@set_scale);
            
            h.B(25) = uicontrol('style','edit',...
                'String','',...
                'units','normalized',...
                'position',[.85 .5 .035 .05],...
                'Callback',@set_scale);
            
            h.B(26) = uicontrol('style','edit',...
                'String','',...
                'units','normalized',...
                'position',[.925 .5 .035 .05],...
                'Callback',@set_scale);
            
            h.B(27) = uicontrol('style','text',...
                'String',{'Min','[]'},...
                'Units','normalized',...
                'Position',[.885 .75 .04 .05]);
            
            h.B(28) = uicontrol('style','text',...
                'String',{'Max','[]'},...
                'Units','normalized',...
                'Position',[.96 .75 .04 .05]);
            
            h.B(29) = uicontrol('style','text',...
                'String',{'Min','[]'},...
                'Units','normalized',...
                'Position',[.885 .5 .04 .05]);
            
            h.B(30) = uicontrol('style','text',...
                'String',{'Max','[]'},...
                'Units','normalized',...
                'Position',[.96 .5 .04 .05]);
            
            h.B(31) = uicontrol('Style','checkbox',...
                'String','Use Overlay',...
                'Value',0,...
                'Tooltip','Green, if both Channels share an overlay',...
                'Units','normalized',...
                'Position',[.85 .39 .1 .03]);
            
            h.Channel1Max = 1;
            h.Channel1Min = 0;
            h.Channel2Max = 1;
            h.Channel2Min = 0;
            
            h.CurChannel1Idx = h.B(2).Value;
            h.CurChannel2Idx = h.B(3).Value;
            h.RelativeChannelIndex = 0;
            
            h.MainLine = [];
            h.ChildLine = [];
            h.hasCrossSection = 0;
            h.hasChannel2 = 0;
            h.isUpscaled = false;
            h.lockedChannels = h.B(20).Value;
            [~,DefIndex] = Class{1}.get_channel('Processed');
            if isempty(DefIndex)
                DefIndex = 2;
            else
                DefIndex = DefIndex + 1;
            end
            h.B(2).Value = DefIndex;
            draw_channel_1
            
            function cross_section_toggle(varargin)
                h.hasCrossSection = ~h.hasCrossSection;
                try
                    delete(h.ImAx(3));
                catch
                end
                draw_channel_1
                draw_channel_2
            end
            
            function draw_channel_1(varargin)
                LeftRight = 'Left';
                if h.hasChannel2 && h.hasCrossSection
                    FullPart = 'PartTwo';
                elseif h.hasChannel2 && ~h.hasCrossSection
                    FullPart = 'FullTwo';
                elseif ~h.hasChannel2 && h.hasCrossSection
                    FullPart = 'PartOne';
                elseif ~h.hasChannel2 && ~h.hasCrossSection
                    FullPart = 'FullOne';
                end
                h.hasChannel1 = true;
                draw_image(LeftRight,FullPart)
                if isequal(h.Channel{1},'none')
                    h.hasChannel1 = false;
                end
                if h.hasChannel2 && ~h.OnePass
                    h.OnePass = true;
                    draw_channel_2
                end
                h.OnePass = false;
            end
            
            function draw_channel_2(varargin)
                LeftRight = 'Right';
                if h.hasChannel1 && h.hasCrossSection
                    FullPart = 'PartTwo';
                elseif h.hasChannel1 && ~h.hasCrossSection
                    FullPart = 'FullTwo';
                elseif ~h.hasChannel1 && h.hasCrossSection
                    FullPart = 'PartOne';
                elseif ~h.hasChannel1 && ~h.hasCrossSection
                    FullPart = 'FullOne';
                end
                h.hasChannel2 = true;
                draw_image(LeftRight,FullPart)
                if isequal(h.Channel{2},'none')
                    h.hasChannel2 = false;
                end
                if h.hasChannel1 && ~h.OnePass
                    h.OnePass = true;
                    draw_channel_1
                end
                h.OnePass = false;
            end
            
            function moving_cross_section_channel_1(src,evt)
                evname = evt.EventName;
                if ~get(h.B(1),'Value')
                    return
                end
                delete(h.ImAx(3))
                h.ImAx(3) = [];
                Pos1 = [h.MainLine.Position(1,1) h.MainLine.Position(1,2)];
                Pos2 = [h.MainLine.Position(2,1) h.MainLine.Position(2,2)];
                MainProfile = improfile(h.Image{h.MainIndex},[Pos1(1) Pos2(1)],[Pos1(2) Pos2(2)],'bicubic');
                Len = norm(Pos1-Pos2)/h.NumPixelsX(h.MainIndex)*h.ScanSizeX(h.MainIndex);
                Points = [0:1/(length(MainProfile)-1):1].*Len;
                [MainMultiplierY,UnitY,~] = AFMImage.parse_unit_scale(range(MainProfile),h.BaseUnit{h.MainIndex},1);
                [MultiplierX,UnitX,~] = AFMImage.parse_unit_scale(range(Points),'m',1);
                h.ImAx(3) = subplot(10,10,[71:78 81:88 91:98]);
                P = plot(Points.*MultiplierX,Profile.*MultiplierY);
                grid on
                CurrentAxHeight = round(h.Fig.Position(4)*h.ImAx(h.MainIndex).Position(4));
                h.ImAx(3).Color = h.ColorMode(h.ColorIndex).Background;
                h.ImAx(3).LineWidth = 1;
                h.ImAx(3).FontSize = round(h.ReferenceFontSize*(CurrentAxHeight/756));
                h.ImAx(3).XColor = h.ColorMode(h.ColorIndex).Text;
                h.ImAx(3).YColor = h.ColorMode(h.ColorIndex).Profile1;
                h.ImAx(3).GridColor = h.ColorMode(h.ColorIndex).Text;
                xlabel(sprintf('[%s]',UnitX))
                ylabel(sprintf('%s [%s]',h.Channel{1},UnitY))
                xlim([0 Points(end).*MultiplierX])
                h.P.LineWidth = 2;
                h.P.Color = h.ColorMode(h.ColorIndex).Profile1;
                if h.hasBothCrossSections && (h.hasChannel2 && h.hasChannel1)
                    if ~isempty(h.ChildLine)
                        if ~isvalid(h.ChildLine)
                            h.ChildLine = [];
                        end
                    end
                    h.ChildLine.Visible = 'off';
                    
                    if h.B(31).Value && ...
                            Class{h.MainIndex}.OverlayGroup.hasOverlayGroup &&...
                            Class{h.ChildIndex}.OverlayGroup.hasOverlayGroup &&...
                            isequal(Class{h.MainIndex}.OverlayGroup.Names,Class{h.ChildIndex}.OverlayGroup.Names)
                        XDiff = Class{h.MainIndex}.Channel(1).OriginX - Class{h.ChildIndex}.Channel(1).OriginX;
                        SizePerPixelX = Class{h.ChildIndex}.Channel(1).ScanSizeX./Class{h.ChildIndex}.Channel(1).NumPixelsX;
                        XDiff = XDiff/SizePerPixelX;
                        YDiff = Class{h.MainIndex}.Channel(1).OriginY - Class{h.ChildIndex}.Channel(1).OriginY;
                        SizePerPixelY = Class{h.ChildIndex}.Channel(1).ScanSizeY./Class{h.ChildIndex}.Channel(1).NumPixelsY;
                        YDiff = YDiff/SizePerPixelY;
                        AngleDiff = Class{h.MainIndex}.Channel(1).ScanAngle - Class{h.ChildIndex}.Channel(1).ScanAngle;
                        AngleDiff = deg2rad(-AngleDiff);
                        
                        InitPos = h.MainLine.Position;
                        
                        CPos1 = [InitPos(1,1) InitPos(1,2)];
                        CPos2 = [InitPos(2,1) InitPos(2,2)];
                        
                        ImCenter = [Class{h.ChildIndex}.Channel(1).NumPixelsX/2 Class{h.ChildIndex}.Channel(1).NumPixelsY/2];
                        
                        TempCP1 = CPos1 - ImCenter;
                        TempCP2 = CPos2 - ImCenter;
                        
                        RotationMatrix = [cos(AngleDiff) -sin(AngleDiff);sin(AngleDiff) cos(AngleDiff)];
                        
                        CPos1 = [RotationMatrix*TempCP1']' + ImCenter + [XDiff -YDiff];
                        CPos2 = [RotationMatrix*TempCP2']' + ImCenter + [XDiff -YDiff];
                        
                        InitPos = [CPos1(1) CPos1(2); CPos2(1) CPos2(2)];
                    else
                        InitPos = h.MainLine.Position;
                    end
                    h.ChildLine = drawline('Position',InitPos,...
                        'Parent',h.ImAx(h.ChildIndex),'Color',h.ColorMode(h.ColorIndex).Profile2,...
                        'LineWidth',h.ProfileLineWidth);
                    addlistener(h.ChildLine,'MovingROI',@moving_cross_section);
                    addlistener(h.ChildLine,'ROIMoved',@moving_cross_section);
                    CPos1 = [h.ChildLine.Position(1,1) h.ChildLine.Position(1,2)];
                    CPos2 = [h.ChildLine.Position(2,1) h.ChildLine.Position(2,2)];
                    ChildProfile = improfile(h.Image{h.ChildIndex},[CPos1(1) CPos2(1)],[CPos1(2) CPos2(2)],'bicubic');
                    ChildPoints = [0:1/(length(ChildProfile)-1):1].*Len;
                    [ChildMultiplierY,UnitY,~] = AFMImage.parse_unit_scale(range(ChildProfile),h.BaseUnit{h.ChildIndex},1);
                    yyaxis right
                    h.CP = plot(ChildPoints.*MultiplierX,ChildProfile.*ChildMultiplierY);
                    grid on
                    CurrentAxHeight = round(h.Fig.Position(4)*h.ImAx(h.ChildIndex).Position(4));
                    h.ImAx(3).Color = h.ColorMode(h.ColorIndex).Background;
                    h.ImAx(3).LineWidth = 1;
                    h.ImAx(3).FontSize = round(h.ReferenceFontSize*(CurrentAxHeight/756));
                    h.ImAx(3).XColor = h.ColorMode(h.ColorIndex).Text;
                    h.ImAx(3).YColor = h.ColorMode(h.ColorIndex).Profile2;
                    h.ImAx(3).GridColor = h.ColorMode(h.ColorIndex).Text;
                    xlabel(sprintf('[%s]',UnitX))
                    ylabel(sprintf('%s [%s]',h.Channel{h.ChildIndex},UnitY))
                    xlim([0 ChildPoints(end).*MultiplierX])
                    h.CP.LineWidth = 2;
                    h.CP.Color = h.ColorMode(h.ColorIndex).Profile2;
                    if isequal(h.BaseUnit{1},h.BaseUnit{2})
                        yyaxis left
                        ylim([min([min(ChildProfile)*ChildMultiplierY min(MainProfile)*MainMultiplierY])...
                            max([max(ChildProfile)*ChildMultiplierY max(MainProfile)*MainMultiplierY])]);
                        yyaxis right
                        ylim([min([min(ChildProfile)*ChildMultiplierY min(MainProfile)*MainMultiplierY])...
                            max([max(ChildProfile)*ChildMultiplierY max(MainProfile)*MainMultiplierY])]);
                    else
                        ylim([min(ChildProfile)*ChildMultiplierY max(ChildProfile)*ChildMultiplierY]);
                    end
                    
%                     % Temporary
%                     legend({'Before','After'},'Location','northwest',...
%                         'FontSize',h.ReferenceFontSize);
                    
                end
                hold off
            end
            
            function moving_cross_section_channel_2(src,evt)
                evname = evt.EventName;
                if ~get(h.B(1),'Value')
                    return
                end
                Pos1 = [h.Line.Position(1,1) h.Line.Position(1,2)];
                Pos2 = [h.Line.Position(2,1) h.Line.Position(2,2)];
                if norm(Pos1-Pos2)==0
                    get_and_draw_profile;
                    return
                end
                Profile = improfile(h.Image{2},[Pos1(1) Pos2(1)],[Pos1(2) Pos2(2)]);
                Len = norm(Pos1-Pos2)/Class{2}.NumPixelsX*Class{2}.ScanSizeX;
                Points = [0:1/(length(Profile)-1):1].*Len;
                [MultiplierY,UnitY,~] = AFMImage.parse_unit_scale(range(Profile),h.BaseUnit{2},1);
                [MultiplierX,UnitX,~] = AFMImage.parse_unit_scale(range(Points),'m',1);
                h.ImAx(3) = subplot(10,10,[71:78 81:88 91:98]);
                P = plot(Points.*MultiplierX,Profile.*MultiplierY);
                grid on
                CurrentAxHeight = round(h.Fig.Position(4)*h.ImAx(2).Position(4));
                h.ImAx(3).Color = 'k';
                h.ImAx(3).LineWidth = 1;
                h.ImAx(3).FontSize = round(22*(CurrentAxHeight/756));
                h.ImAx(3).XColor = 'w';
                h.ImAx(3).YColor = 'w';
                h.ImAx(3).GridColor = 'w';
                xlabel(sprintf('[%s]',UnitX))
                ylabel(sprintf('%s [%s]',h.Channel{2},UnitY))
                xlim([0 Points(end).*MultiplierX])
                P.LineWidth = 2;
                P.Color = 'b';
            end
            
            function get_and_draw_profile_channel_1(varargin)
                if ~get(h.B(1),'Value')
                    return
                end
                if ~isempty(h.Line)
                    if ~isvalid(h.Line)
                        h.Line = [];
                    end
                end
                h.Line.Visible = 'off';
                h.Line = drawline('Color','b','Parent',h.ImAx(1));
                addlistener(h.Line,'MovingROI',@moving_cross_section_channel_1);
                addlistener(h.Line,'ROIMoved',@moving_cross_section_channel_1);
                Pos1 = [h.Line.Position(1,1) h.Line.Position(1,2)];
                Pos2 = [h.Line.Position(2,1) h.Line.Position(2,2)];
                if norm(Pos1-Pos2)==0
                    get_and_draw_profile_channel_1;
                    return
                end
                Profile = improfile(h.Image{1},[Pos1(1) Pos2(1)],[Pos1(2) Pos2(2)]);
                Len = norm(Pos1-Pos2)/Class{1}.NumPixelsX*Class{1}.ScanSizeX;
                Points = [0:1/(length(Profile)-1):1].*Len;
                [MultiplierY,UnitY,~] = AFMImage.parse_unit_scale(range(Profile),h.BaseUnit{1},1);
                [MultiplierX,UnitX,~] = AFMImage.parse_unit_scale(range(Points),'m',1);
                h.ImAx(3) = subplot(10,10,[71:78 81:88 91:98]);
                P = plot(Points.*MultiplierX,Profile.*MultiplierY);
                grid on
                CurrentAxHeight = round(h.Fig.Position(4)*h.ImAx(1).Position(4));
                h.ImAx(3).Color = 'k';
                h.ImAx(3).LineWidth = 1;
                h.ImAx(3).FontSize = round(22*(CurrentAxHeight/756));
                h.ImAx(3).XColor = 'w';
                h.ImAx(3).YColor = 'w';
                h.ImAx(3).GridColor = 'w';
                xlabel(sprintf('[%s]',UnitX))
                ylabel(sprintf('%s [%s]',h.Channel{1},UnitY))
                xlim([0 Points(end).*MultiplierX])
                P.LineWidth = 2;
                P.Color = 'b';
            end
            
            function get_and_draw_profile_channel_2(varargin)
                if ~get(h.B(1),'Value')
                    return
                end
                if ~isempty(h.Line)
                    if ~isvalid(h.Line)
                        h.Line = [];
                    end
                end
                h.MainLine = drawline('Color',h.ColorMode(h.ColorIndex).Profile1,...
                    'Parent',varargin{1}.Parent,'LineWidth',h.ProfileLineWidth);
                addlistener(h.MainLine,'MovingROI',@moving_cross_section);
                addlistener(h.MainLine,'ROIMoved',@moving_cross_section);
                Pos1 = [h.MainLine.Position(1,1) h.MainLine.Position(1,2)];
                Pos2 = [h.MainLine.Position(2,1) h.MainLine.Position(2,2)];
                if norm(Pos1-Pos2)==0
                    get_and_draw_profile_channel_2;
                    return
                end
                MainProfile = improfile(h.Image{h.MainIndex},[Pos1(1) Pos2(1)],[Pos1(2) Pos2(2)],'bicubic');
                Len = norm(Pos1-Pos2)/h.NumPixelsX(h.MainIndex)*h.ScanSizeX(h.MainIndex);
                Points = [0:1/(length(MainProfile)-1):1].*Len;
                [MainMultiplierY,UnitY,~] = AFMImage.parse_unit_scale(range(MainProfile),h.BaseUnit{h.MainIndex},1);
                [MultiplierX,UnitX,~] = AFMImage.parse_unit_scale(range(Points),'m',1);
                h.ImAx(3) = subplot(10,10,[71:78 81:88 91:98]);
                P = plot(Points.*MultiplierX,Profile.*MultiplierY);
                grid on
                CurrentAxHeight = round(h.Fig.Position(4)*h.ImAx(h.MainIndex).Position(4));
                h.ImAx(3).Color = h.ColorMode(h.ColorIndex).Background;
                h.ImAx(3).LineWidth = 1;
                h.ImAx(3).FontSize = round(h.ReferenceFontSize*(CurrentAxHeight/756));
                h.ImAx(3).XColor = h.ColorMode(h.ColorIndex).Text;
                h.ImAx(3).YColor = h.ColorMode(h.ColorIndex).Profile1;
                h.ImAx(3).GridColor = h.ColorMode(h.ColorIndex).Text;
                xlabel(sprintf('[%s]',UnitX))
                ylabel(sprintf('%s [%s]',h.Channel{2},UnitY))
                xlim([0 Points(end).*MultiplierX])
                h.P.LineWidth = 2;
                h.P.Color = h.ColorMode(h.ColorIndex).Profile1;
                if h.hasBothCrossSections && (h.hasChannel2 && h.hasChannel1)
                    
                    if h.B(31).Value && ...
                            Class{h.MainIndex}.OverlayGroup.hasOverlayGroup &&...
                            Class{h.ChildIndex}.OverlayGroup.hasOverlayGroup &&...
                            isequal(Class{h.MainIndex}.OverlayGroup.Names,Class{h.ChildIndex}.OverlayGroup.Names)
                        XDiff = Class{h.MainIndex}.Channel(1).OriginX - Class{h.ChildIndex}.Channel(1).OriginX;
                        SizePerPixelX = Class{h.ChildIndex}.Channel(1).ScanSizeX./Class{h.ChildIndex}.Channel(1).NumPixelsX;
                        XDiff = XDiff/SizePerPixelX;
                        YDiff = Class{h.MainIndex}.Channel(1).OriginY - Class{h.ChildIndex}.Channel(1).OriginY;
                        SizePerPixelY = Class{h.ChildIndex}.Channel(1).ScanSizeY./Class{h.ChildIndex}.Channel(1).NumPixelsY;
                        YDiff = YDiff/SizePerPixelY;
                        AngleDiff = Class{h.MainIndex}.Channel(1).ScanAngle - Class{h.ChildIndex}.Channel(1).ScanAngle;
                        AngleDiff = deg2rad(-AngleDiff);
                        
                        InitPos = h.MainLine.Position;
                        
                        CPos1 = [InitPos(1,1) InitPos(1,2)];
                        CPos2 = [InitPos(2,1) InitPos(2,2)];
                        
                        ImCenter = [Class{h.ChildIndex}.Channel(1).NumPixelsX/2 Class{h.ChildIndex}.Channel(1).NumPixelsY/2];
                        
                        TempCP1 = CPos1 - ImCenter;
                        TempCP2 = CPos2 - ImCenter;
                        
                        RotationMatrix = [cos(AngleDiff) -sin(AngleDiff);sin(AngleDiff) cos(AngleDiff)];
                        
                        CPos1 = [RotationMatrix*TempCP1']' + ImCenter + [XDiff -YDiff];
                        CPos2 = [RotationMatrix*TempCP2']' + ImCenter + [XDiff -YDiff];
                        
                        InitPos = [CPos1(1) CPos1(2); CPos2(1) CPos2(2)];
                    else
                        InitPos = h.MainLine.Position;
                    end
                    h.ChildLine = drawline('Position',InitPos,...
                        'Parent',h.ImAx(h.ChildIndex),'Color',h.ColorMode(h.ColorIndex).Profile2,...
                        'LineWidth',h.ProfileLineWidth);
                    addlistener(h.ChildLine,'MovingROI',@moving_cross_section);
                    addlistener(h.ChildLine,'ROIMoved',@moving_cross_section);
                    CPos1 = [h.ChildLine.Position(1,1) h.ChildLine.Position(1,2)];
                    CPos2 = [h.ChildLine.Position(2,1) h.ChildLine.Position(2,2)];
                    ChildProfile = improfile(h.Image{h.ChildIndex},[CPos1(1) CPos2(1)],[CPos1(2) CPos2(2)],'bicubic');
                    ChildPoints = [0:1/(length(ChildProfile)-1):1].*Len;
                    [ChildMultiplierY,UnitY,~] = AFMImage.parse_unit_scale(range(ChildProfile),h.BaseUnit{h.ChildIndex},1);
                    yyaxis right
                    h.CP = plot(ChildPoints.*MultiplierX,ChildProfile.*ChildMultiplierY);
                    grid on
                    CurrentAxHeight = round(h.Fig.Position(4)*h.ImAx(h.MainIndex).Position(4));
                    h.ImAx(3).Color = h.ColorMode(h.ColorIndex).Background;
                    h.ImAx(3).LineWidth = 1;
                    h.ImAx(3).FontSize = round(h.ReferenceFontSize*(CurrentAxHeight/756));
                    h.ImAx(3).XColor = h.ColorMode(h.ColorIndex).Text;
                    h.ImAx(3).YColor = h.ColorMode(h.ColorIndex).Profile2;
                    h.ImAx(3).GridColor = h.ColorMode(h.ColorIndex).Text;
                    xlabel(sprintf('[%s]',UnitX))
                    ylabel(sprintf('%s [%s]',h.Channel{h.ChildIndex},UnitY))
                    xlim([0 ChildPoints(end).*MultiplierX])
                    if isequal(h.BaseUnit{1},h.BaseUnit{2})
                        yyaxis left
                        ylim([min([min(ChildProfile)*ChildMultiplierY min(MainProfile)*MainMultiplierY])...
                            max([max(ChildProfile)*ChildMultiplierY max(MainProfile)*MainMultiplierY])]);
                        yyaxis right
                        ylim([min([min(ChildProfile)*ChildMultiplierY min(MainProfile)*MainMultiplierY])...
                            max([max(ChildProfile)*ChildMultiplierY max(MainProfile)*MainMultiplierY])]);
                    else
                        ylim([min(ChildProfile)*ChildMultiplierY max(ChildProfile)*ChildMultiplierY]);
                    end
                    h.CP.LineWidth = 2;
                    h.CP.Color = h.ColorMode(h.ColorIndex).Profile2;
                end
                hold off
            end
            
            function draw_image(LeftRight,FullPart)
                if isequal(LeftRight,'Left')
                    Index = 1;
                elseif isequal(LeftRight,'Right')
                    Index = 2;
                end
                BarToImageRatio = 1/5;
                try
                    delete(h.ImAx(Index));
                    delete(h.I(Index));
                catch
                end
                if isequal(FullPart,'FullOne')
                    h.ImAx(Index) = axes(h.Fig,'Position',[0.1 0.1 .6 .8]);
                elseif isequal(FullPart,'FullTwo')
                    if isequal(LeftRight,'Left')
                    h.ImAx(Index) = axes(h.Fig,'Position',[0.12 0.1 .3 .8]);
                    else
                    h.ImAx(Index) = axes(h.Fig,'Position',[.47 0.1 .3 .8]);
                    end
                elseif isequal(FullPart,'PartOne')
                    h.ImAx(Index) = axes(h.Fig,'Position',[0.1 .35 .6 .6]);
                elseif isequal(FullPart,'PartTwo')
                    if isequal(LeftRight,'Left')
                    h.ImAx(Index) = axes(h.Fig,'Position',[0.12 .35 .3 .6]);
                    else
                    h.ImAx(Index) = axes(h.Fig,'Position',[0.47 .35 .3 .6]);
                    end
                end
                
                
                if h.lockedChannels && Index==1
                    h.CurChannel1Idx = h.B(16).Value;
                    h.B(17).Value = mod(h.CurChannel1Idx + h.RelativeChannelIndex,h.NumClasses+1);
                    h.CurChannel2Idx = h.B(17).Value;
                    CurIndex = h.CurChannel1Idx;
                elseif h.lockedChannels && Index==2
                    h.CurChannel2Idx = h.B(17).Value;
                    h.B(16).Value = mod(h.CurChannel2Idx - h.RelativeChannelIndex,h.NumClasses+1);
                    h.CurChannel1Idx = h.B(16).Value;
                    CurIndex = h.CurChannel2Idx;
                else
                    CurIndex = h.B(15+Index).Value;
                    h.CurChannel2Idx = h.B(17).Value;
                    h.CurChannel1Idx = h.B(16).Value;
                end
                
                h.RelativeChannelIndex = h.CurChannel2Idx - h.CurChannel1Idx;
                
                Class{Index} = obj.get_class_instance(ClassIndex(CurIndex,:));
                CurrentChannelName = h.B(1+Index).String{h.B(1+Index).Value};
                PopUp = Class{Index}.string_of_existing();
                set(h.B(1+Index),'String',PopUp)
                
                if h.B(1+Index).Value > (length(Class{Index}.Channel) + 1)
                    set(h.B(1+Index),'Value',1)
                end
                
                h.Channel{Index} = h.B(1+Index).String{h.B(1+Index).Value};
                if isequal(h.Channel{Index},'none')
                    try
                        delete(h.ImAx(Index));
                        delete(h.I(Index));
                    catch
                    end
                    if Index == 1
                        h.hasChannel1 = 0;
                    elseif Index == 2
                        h.hasChannel2 = 0;
                    end
                    return
                else
                    [Channel,ChannelIndex] = Class{Index}.get_channel(h.Channel{Index});
                    if h.isUpscaled
                        Channel.Image = fillmissing(Channel.Image,'linear','EndValues','nearest');
                        Channel = AFMImage.resize_channel(Channel,1,1920,true);
                    end
                    h.Image{Index} = fillmissing(Channel.Image,'linear','EndValues','nearest');
                    h.BaseUnit{Index} = Channel.Unit;
                    h.ScanSizeX(Index) = Channel.ScanSizeX;
                    h.ScanSizeY(Index) = Channel.ScanSizeY;
                    h.NumPixelsX(Index) = Channel.NumPixelsX;
                    h.NumPixelsY(Index) = Channel.NumPixelsY;
                    ColorPattern = Class{Index}.CMap;
                end
                
                if h.B(21).Value && h.hasChannel2 && h.hasChannel1 && isequal(h.BaseUnit{1},h.BaseUnit{2})
                    CurImage = h.Image{Index};
                    Range = range(CurImage,'all');
                    OtherImage = h.Image{mod(Index,2)+1};
                    OtherRange = range(OtherImage,'all');
                    
                    [FinalRange,FinalIndex] = max([Range OtherRange]);
                    
                    if FinalIndex==1
                        Min = min(CurImage,[],'all');
                    else
                        Min = min(OtherImage,[],'all');
                    end
                    
                    if ~h.B(22).Value
                        CutMax = FinalRange*h.Channel1Max + Min;
                        CutMin = FinalRange*h.Channel1Min + Min;
                    else
                        CutMax = quantile(h.Image{FinalIndex},h.Channel1Max,'all') + Min;
                        CutMin = quantile(h.Image{FinalIndex},h.Channel1Min,'all') + Min;
                    end
                else
                    CurImage = h.Image{Index};
                    FinalRange = range(CurImage,'all');
                    if Index==1
                        if ~h.B(22).Value
                            CutMax = FinalRange*h.Channel1Max + min(CurImage,[],'all');
                            CutMin = FinalRange*h.Channel1Min + min(CurImage,[],'all');
                        else
                            CutMax = quantile(CurImage,h.Channel1Max,'all') + min(CurImage,[],'all');
                            CutMin = quantile(CurImage,h.Channel1Min,'all') + min(CurImage,[],'all');
                        end
                    elseif Index==2
                        if ~h.B(22).Value
                            CutMax = FinalRange*h.Channel2Max + min(CurImage,[],'all');
                            CutMin = FinalRange*h.Channel2Min + min(CurImage,[],'all');
                        else
                            CutMax = quantile(CurImage,h.Channel2Max,'all') + min(CurImage,[],'all');
                            CutMin = quantile(CurImage,h.Channel2Min,'all') + min(CurImage,[],'all');
                        end
                    end
                end
                
                MinIndex = 21+2*Index;
                MaxIndex = 22+2*Index;
                
                [h.Multiplier{Index},h.Unit{Index},~] = AFMImage.parse_unit_scale(FinalRange,h.BaseUnit{Index},1);
                
                if ~isempty(get(h.B(MinIndex),'String')) && ~isempty(get(h.B(MaxIndex),'String'))
                    MinValue = str2double(get(h.B(MinIndex),'String'))/h.Multiplier{Index};
                    MaxValue = str2double(get(h.B(MaxIndex),'String'))/h.Multiplier{Index};
                    if isnan(MinValue) || isnan(MaxValue)
                    else
                        CutMin = MinValue;
                        CutMax = MaxValue;
                    end
                end
                
                CurImage(CurImage>CutMax) = CutMax;
                CurImage(CurImage<CutMin) = CutMin;
                
                if isempty(CurImage(CurImage>=CutMax))
                    CurImage(1,1) = CutMax;
                end
                if isempty(CurImage(CurImage<=CutMin))
                    CurImage(end,end) = CutMin;
                end
                
                FinalRange = abs(CutMax - CutMin);
                
                [h.Multiplier{Index},h.Unit{Index},~] = AFMImage.parse_unit_scale(FinalRange,h.BaseUnit{Index},1);
                h.I(Index) = imshow(CurImage*h.Multiplier{Index},[],'Colormap',ColorPattern);
                h.I(Index).ButtonDownFcn = @get_and_draw_profile;
                hold on
                CurrentAxHeight = round(h.Fig.Position(4)*h.ImAx(Index).Position(4));
                CurrentAxWidth = round(h.Fig.Position(3)*h.ImAx(Index).Position(3));
                AFMImage.draw_scalebar_into_current_image(Channel.NumPixelsX,Channel.NumPixelsY,Channel.ScanSizeX,BarToImageRatio,CurrentAxHeight,CurrentAxWidth);
                c = colorbar('northoutside');
                c.FontSize = round(h.ReferenceFontSize*(CurrentAxHeight/756));
                c.Color = h.ColorMode(h.ColorIndex).Text;
                c.Label.String = sprintf('%s [%s]',h.Channel{Index},h.Unit{Index});
                c.Label.FontSize = round(h.ReferenceFontSize*(CurrentAxHeight/756));
                c.Label.Color = h.ColorMode(h.ColorIndex).Text;
                
                set(h.B(MinIndex+4),'String',{'Min',sprintf('[%s]',h.Unit{Index})});
                set(h.B(MaxIndex+4),'String',{'Max',sprintf('[%s]',h.Unit{Index})});
                
                try
                    if (h.ScanSizeX(1) == h.ScanSizeX(2)) &&...
                            (h.ScanSizeY(1) == h.ScanSizeY(2))
                        h.B(18).BackgroundColor = 'g';
                    else
                        h.B(18).BackgroundColor = 'r';
                    end
                catch
                end
            end
            
            function changed_slider(varargin)
                C1Max = get(h.B(8),'value');
                C1Min = get(h.B(9),'value');
                C2Max = get(h.B(10),'value');
                C2Min = get(h.B(11),'value');
                
                if C1Max <= C1Min
                    C1Min = C1Max*0.9;
                    set(h.B(8),'value',C1Max);
                    set(h.B(9),'value',C1Min);
                end
                if C2Max <= C2Min
                    C2Min = C2Max*0.9;
                    set(h.B(10),'value',C2Max);
                    set(h.B(11),'value',C2Min);
                end
                
                h.Channel1Max = C1Max;
                h.Channel1Min = C1Min;
                h.Channel2Max = C2Max;
                h.Channel2Min = C2Min;
                
                draw_channel_1
                draw_channel_2
            end
            
            function save_figure_to_file(varargin)
                
                filter = {'*.png';'*.tif'};
                [file, path] = uiputfile(filter);
                FullFile = fullfile(path,file);
                exportgraphics(h.Fig,FullFile,'Resolution',300,'BackgroundColor','current')
            end
            
            function changed_color(varargin)
                
                if ~h.B(7).Value
                    h.ColorIndex = 1;
                else
                    h.ColorIndex = 2;
                end
                
                h.Fig.Color = h.ColorMode(h.ColorIndex).Background;
                
                draw_channel_1
                draw_channel_2
                
            end
            
            function upscale_images(varargin)
                
                h.isUpscaled = h.B(19).Value;
                
                draw_channel_1
                draw_channel_2
                
            end
            
            function lock_channels(varargin)
                
                h.lockedChannels = h.B(20).Value;
                
                h.RelativeChannelIndex = h.CurChannel2Idx - h.CurChannel1Idx;
                draw_channel_1
                draw_channel_2
                
            end
            
            function lock_scalebars(varargin)
                
                if h.B(21).Value
                    h.B(10).Visible = 'off';
                    h.B(11).Visible = 'off';
                    h.B(14).Visible = 'off';
                    h.B(15).Visible = 'off';
                else
                    h.B(10).Visible = 'on';
                    h.B(11).Visible = 'on';
                    h.B(14).Visible = 'on';
                    h.B(15).Visible = 'on';
                end
                
                draw_channel_1
                draw_channel_2
            end
            
            function set_scale(varargin)
                
                draw_channel_1
                draw_channel_2
            end
            
            function statistical_cmapping(varargin)
                draw_channel_1
                draw_channel_2
            end
            
            uiwait(h.Fig)
        end
        
        function [Fig,DataMat] = visualize_listed_data(obj,Property,YAxisName,BaseUnit,Method,ErrorBars,UseGrouping,ListOfIndizes)
            % [Fig,DataMat] = visualize_listed_data(obj,Property,YAxisName,BaseUnit,Method,ErrorBars,UseGrouping,ListOfIndizes)
            %
            % visualize_listed_data(varargin)
            % Input arbitrary ForceMap listed Property
            % e.g. 'EModOliverPharr' or 'IndentationDepth'
            % Methods are 'Boxplot', 'Barplot', 'BarMedian'
            % ErrorBars, 'std'...standard deviation
            %            'ste'...standard error
            %            'ci'...confidence interval
            
            if nargin < 8
                ListOfIndizes = 1:obj.NumForceMaps;
            end
            
            k = 1;
            if UseGrouping
                for i=1:length(obj.GroupFM)
                    TempArray = [];
                    for j=1:length(obj.GroupFM(i).Indices)
                        TempArray = cat(1,TempArray,reshape(obj.FM{obj.GroupFM(i).Indices(j)}.get(Property),[],1));
                    end
                    Data(k).Values = TempArray;
                    [Data(k).Multiplier,Data(k).Unit] = AFMImage.parse_unit_scale(range(Data(k).Values),BaseUnit,5);
                    Data(k).Name = obj.GroupFM(i).Name;
                    Data(k).Length = length(Data(k).Values);
                    if isequal(ErrorBars,'std')
                        Data(k).Bars = nanstd(Data(k).Values);
                    elseif isequal(ErrorBars,'ste')
                        Data(k).Bars = nanstd(Data(k).Values)/sqrt(length(Data(k).Values(~isnan(Data(k).Values))));
                    elseif isequal(ErrorBars,'ci')
                        [~,~,Data(k).Bars,~] = ttest(Data(k).Values(~isnan(Data(k).Values)),nanmean(Data(k).Values));
                    end
                    k = k + 1;
                end
            else
                for i=ListOfIndizes
                    Data(k).Values = obj.FM{i}.get(Property);
                    [Data(k).Multiplier,Data(k).Unit] = AFMImage.parse_unit_scale(range(Data(k).Values),BaseUnit,5);
                    Data(k).Name = obj.FM{i}.Name;
                    Data(k).Length = length(Data(k).Values);
                    if isequal(ErrorBars,'std')
                        Data(k).Bars = nanstd(Data(k).Values);
                    elseif isequal(ErrorBars,'ste')
                        Data(k).Bars = nanstd(Data(k).Values)/sqrt(length(Data(k).Values(~isnan(Data(k).Values))));
                    elseif isequal(ErrorBars,'ci')
                        [~,~,Data(k).Bars,~] = ttest(Data(k).Values(~isnan(Data(k).Values)),nanmean(Data(k).Values));
                    end
                    k = k + 1;
                end
            end
            
            [MaxMult,MaxIdx] = max([Data(:).Multiplier],[],'all','linear');
            MaxLen = max([Data(:).Length]);
            
            Scale = MaxMult;
            Unit = Data(MaxIdx).Unit;
            YAxisLabel = sprintf('%s [%s]',YAxisName,Unit);
            XTickLabels = {Data(:).Name};
            
            Fig = figure('Color','w');
            DataMat = nan(MaxLen,length(Data));
            for i=1:length(Data)
                DataMat(1:Data(i).Length,i) = Data(i).Values;
                Data(i).Bars = Data(i).Bars.*Scale;
            end
            DataMat = DataMat.*Scale;
            if isequal(Method,'Boxplot')
                boxplot(DataMat)
                title(sprintf('Boxplot of %s',YAxisName))
            elseif isequal(Method,'Barplot')
                bar(nanmean(DataMat,1));
                hold on
                for i=1:length(Data)
                    errorbar(i,nanmean(DataMat(:,i)),-Data(i).Bars(1),Data(i).Bars(end),'Color','k','LineWidth',1.5);
                end
                title(sprintf('Barplot of %s',YAxisName))
            elseif isequal(Method,'BarMedian')
                bar(nanmedian(DataMat,1));
                hold on
                for i=1:length(Data)
                    errorbar(i,nanmedian(DataMat(:,i)),-Data(i).Bars(1),Data(i).Bars(end),'Color','k','LineWidth',1.5);
                end
                title(sprintf('Median Barplot of %s',YAxisName))
            end
            ax = gca;
            ax.LineWidth = 1;
            ax.FontSize = 16;
            ax.TickLabelInterpreter = 'latex';
            ax.XTickLabel = XTickLabels;
            ylabel(YAxisLabel)
            xlim([.5 length(Data)+.5])
        end
        
        function visualize_depth_dependend_tip_radius(obj,CantileverTipIndex)
            
            TipClass = obj.CantileverTips{CantileverTipIndex};
            DepthRadii = TipClass.DepthDependendTipRadius;
            
            
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
            
            N = obj.NumForceMaps;
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
            
            for i=1:obj.NumForceMaps
                obj.FM{i}.FibrilEModOliverPharr = DataMeansOP(i);
                obj.FM{i}.FibrilEModHertz = DataMeansHS(i);
            end
            
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
                
                %Statistics for Hertz-Sneddon Method
                [hHS(i),pHS(i)] = ...
                    ttest(DataMeansHS(obj.GroupFM(TestMat(i,2)).Indices),...
                    DataMeansHS(obj.GroupFM(TestMat(i,1)).Indices),'Tail','right');
                
                figure('Name','Paired Right Tailed T-Test','Units','normalized','Position',[0.2 0.2 0.5 0.5],'Color','w')
                boxplot([DataMeansHS(obj.GroupFM(TestMat(i,1)).Indices) DataMeansHS(obj.GroupFM(TestMat(i,2)).Indices)])
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
            
            % For Hertz-Sneddon
            [hHS,pHS,ciHS,statsHS] = ttest2(DiffMGOHS,DiffControlHS);
            figure('Name','Two Sample T-Test','Units','normalized','Position',[0.2 0.2 0.5 0.5],'Color','w')
            yyaxis left
            boxplot([DiffControlHS DiffMGOHS])
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
        
        function FibPot = statistical_analysis_surface_potential(obj)
            % Statistical analysis of the surface potential map experiment.
            % First get the grouping of data from the user and then perform
            % several tests, plots etc.
            
            warning('The following methods were programmed for specific use cases and are yet to be generalized! However, you can of course adjust them for your own use, though its probably easier to just take the processed raw data from your Experiment and do your own statistics')
            
            obj.grouping_surface_potential_map();
            
            %%%%%%% DISCLAIMER: just works for specific cases at the moment %%%%%%%
            
            N = obj.NumSurfacePotentialMaps;
            
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
            
            N = obj.NumForceMaps;
            
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
            
            N = obj.NumSurfacePotentialMaps;
            
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
            Fig = figure('Units', 'Normalized', 'Position',[0.2 0.2 0.2 .7],'Color','w');
            Names = reshape(Names,[],1);
            T = table(Names);
            Widths = cell(1,length(Names));
            for i=1:length(Widths)
                Widths{i} = 10000;
            end
            uitable('Data',T{:,:},'Units', 'Normalized', 'Position',[0, 0, 1, 1],'ColumnWidth', Widths);
            
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
        
        function RadiusNM = calculate_tip_radius(obj,TipDepthNM,TipIndex)
            if nargin < 2
                TipDepthNM = 20;
            end
            Chan = obj.CantileverTips{TipIndex}.get_channel('Eroded Tip');
            TipImage = -Chan.Image;
            m = 1;
            Niter = 100;
            Cumulative = zeros(Niter,1);
            for n=1:Niter
                k = 1;
                SizeArray = size(TipImage);
                for i=1:SizeArray(1)
                    for j=1:SizeArray(2)
                        if TipImage(i,j) > -TipDepthNM*1e-9
                            X(k) = (i-1)*obj.CantileverTips{TipIndex}.ScanSizeX/SizeArray(1);
                            Y(k) = (j-1)*obj.CantileverTips{TipIndex}.ScanSizeY/SizeArray(2);
                            Z(k) = TipImage(i,j);
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
%             obj.CantileverTip.RadiusNM = RadiusNM;
%             for i=1:obj.NumForceMaps
%                 obj.FM{i}.TipRadius = RadiusNM;
%             end
        end
        
        function check_for_new_host(obj)
            % Author: Manuel Rufin
            % Description:
            % This function verifies if the system environment has changed and, if so,
            % adjusts object properties.
            
            %% Main function
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
        
        function reference_slope_parser(obj,DefaultOption)
            
            if nargin < 2
                DefaultOption = 1;
            end
            
            Methods = false(6,1);
            Methods(DefaultOption) = true;
            [ChosenMethod, AppRetSwitch] = obj.reference_slope_parser_gui(Methods);
            Methods = false(6,1);
            Methods(ChosenMethod) = true;
            
            obj.ReferenceSlopeFlag.SetAllToValue = Methods(1);
            obj.ReferenceSlopeFlag.UserInput = Methods(2);
            obj.ReferenceSlopeFlag.FromRefFM = Methods(3);
            obj.ReferenceSlopeFlag.FromArea = Methods(4);
            obj.ReferenceSlopeFlag.AutomaticFibril = Methods(5);
            obj.ReferenceSlopeFlag.Automatic = Methods(6);
            obj.ReferenceSlopeFlag.AppRetSwitch = AppRetSwitch;
            
            % Process all the RefSlope-Methods that can be done frontloaded
            if obj.ReferenceSlopeFlag.SetAllToValue
                GetValue = inputdlg('To which value should the reference slopes be set?','Choose a value',1);
                Value = str2double(GetValue{1});
                for i=1:obj.NumForceMaps
                    obj.FM{i}.set_reference_slope_to_value(Value)
                end
            elseif obj.ReferenceSlopeFlag.UserInput
                for i=1:obj.NumForceMaps
                    obj.FM{i}.set_reference_slope_to_user_input
                end
            elseif obj.ReferenceSlopeFlag.FromArea
                for i=1:obj.NumForceMaps
                    Mask = obj.FM{i}.create_mask_general;
                    obj.FM{i}.calculate_reference_slope_from_area(Mask)
                end
            end
        end
        
        function reference_slope_calculator(obj,Index)
            
            i = Index;
            if obj.ReferenceSlopeFlag.FromRefFM
                if isempty(obj.RefFM)
                    Warn = warndlg({'There are no reference Force Maps in your Experiment','Please load in one or more reference files'});
                    uiwait(Warn)
                    obj.load_in_reference_maps
                end
                
                if length(obj.RefFM) > 1
                    if ~obj.AssignedReferenceMaps
                        obj.assign_reference_force_map
                    end
                    RefIdx = obj.WhichRefMap(i);
                    if ~obj.RefFM{RefIdx}.HasRefSlope
                        if ~obj.RefFM{RefIdx}.BaseAndTiltFlag
                            obj.RefFM{RefIdx}.base_and_tilt
                        end
                        obj.RefFM{1}.estimate_cp_old
                        Mask = ones(obj.RefFM{RefIdx}.NumProfiles,obj.RefFM{RefIdx}.NumPoints);
                        obj.RefFM{RefIdx}.calculate_reference_slope_from_area(Mask);
                    end
                    obj.FM{i}.RefSlope = obj.RefFM{RefIdx}.RefSlope;
                    obj.FM{i}.HasRefSlope = true;
                elseif length(obj.RefFM) == 1
                    if ~obj.RefFM{1}.HasRefSlope
                        if ~obj.RefFM{1}.BaseAndTiltFlag
                            obj.RefFM{1}.base_and_tilt
                        end
                        obj.RefFM{1}.estimate_cp_old
                        Mask = ones(obj.RefFM{1}.NumProfiles,obj.RefFM{1}.NumPoints);
                        obj.RefFM{1}.calculate_reference_slope_from_area(Mask);
                    end
                    obj.FM{i}.RefSlope = obj.RefFM{1}.RefSlope;
                    obj.FM{i}.HasRefSlope = true;
                end
            elseif obj.ReferenceSlopeFlag.AutomaticFibril
                Mask = obj.FM{i}.BackgroundMask;
                obj.FM{i}.calculate_reference_slope_from_area(Mask,obj.ReferenceSlopeFlag.AppRetSwitch)
            elseif obj.ReferenceSlopeFlag.Automatic
                obj.FM{i}.create_automatic_background_mask(.8)
                Mask = obj.FM{i}.BackgroundMask;
                obj.FM{i}.calculate_reference_slope_from_area(Mask,obj.ReferenceSlopeFlag.AppRetSwitch)
            elseif obj.ReferenceSlopeFlag.FromArea
                obj.FM{i}.calculate_reference_slope_from_area(obj.FM{i}.RefSlopeMask,obj.ReferenceSlopeFlag.AppRetSwitch)
            end
        end
        
        function assign_reference_force_map(obj,DefaultValues)
            
            
            obj.WhichRefMap = zeros(obj.NumForceMaps,1);
            
            NGroups = length(obj.RefFM);
            
            if nargin < 2
                for i=1:NGroups
                    DefaultValues{i} = 'e.g. 1 2 3 4 8 9 10 or 1:4 8:10';
                end
            end
            
            % create the appropriate inputdlg for assigning the groups show
            % a table with numbered map-names in background
            Names = obj.ForceMapNames;
            Fig = figure('Name','Names and coresponding numbers of your Force Maps','Units', 'Normalized', 'Position',[0.1 0.1 0.3 0.4],'Color','w');
            T = table(Names');
            uitable('Data',T{:,:},'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
            
            for i=1:NGroups
                prompts{i} = sprintf('Which Force Maps belong to Reference Force Map %i?',i);
                definput{i} = DefaultValues{i};
            end
            
            dims = [1 50];
            opts.WindowStyle = 'normal';
            dlgtitle = 'Assign a the Force Maps to their respective Reference Maps';
            answer = inputdlg(prompts,dlgtitle,dims,definput,opts);
            
            for i=1:NGroups
                obj.WhichRefMap(str2num(answer{i})) = i;
            end
            
            if  ~isempty(obj.WhichRefMap(obj.WhichRefMap == 0))
                Warn = warndlg('You need to assign exactly one Reference Map to every Force Map','Parsing Error');
                close(Fig);
                uiwait(Warn);
                obj.assign_reference_force_map(answer)
                return
            end
            
            obj.AssignedReferenceMaps = true;
            close(Fig);
            
        end
        
        function assign_cantilever_tips(obj,DefaultValues)
            
            obj.WhichTip = zeros(obj.NumForceMaps,1);
            
            NGroups = obj.NumCantileverTips;
            
            if nargin < 2
                for i=1:NGroups
                    DefaultValues{i} = 'e.g. 1 2 3 4 8 9 10 or 1:4 8:10';
                end
            end
            
            % create the appropriate inputdlg for assigning the groups show
            % a table with numbered map-names in background
            Names = obj.ForceMapNames;
            Fig = figure('Name','Names and corresponding numbers of your Force Maps','Units', 'pixels', 'Position',[100 200 400 800],'Color','w');
            T = table(Names');
            uitable('Data',T{:,:},'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
            
            for i=1:NGroups
                prompts{i} = sprintf('Which Force Maps belong to Cantilever Tip %i, %s?',i,obj.CantileverTipNames{i});
                definput{i} = DefaultValues{i};
            end
            
            dims = [1 50];
            opts.WindowStyle = 'normal';
            dlgtitle = 'Assign the Force Maps to their respective Cantilever Tips';
            answer = inputdlg(prompts,dlgtitle,dims,definput,opts);
            
            for i=1:NGroups
                obj.WhichTip(str2num(answer{i})) = i;
            end
            
            if  ~isempty(obj.WhichTip(obj.WhichTip == 0))
                Warn = warndlg('You need to assign exactly one Tip to every Force Map','Parsing Error');
                close(Fig);
                uiwait(Warn);
                obj.assign_cantilever_tips(answer)
                return
            end
            
            obj.AssignedCantileverTips = true;
            close(Fig);
            
        end
        
        function load_in_reference_maps(obj)
            
            prompt = 'Enter Number of reference force maps';
            dlgtitle = 'Experiment Layout';
            dims = [1 35];
            definput = {'1'};
            answerRefN = inputdlg(prompt,dlgtitle,dims,definput);
            NRef = str2double(answerRefN{1});
            
            MapFullFile = {};
            k = 1;
            while length(MapFullFile) < NRef
                Title = sprintf('Choose one or more REFERENCEs .jpk-force-map or .jpk-qi-data files. %i/%i',length(MapFullFile),NRef);
                [TempFile,TempPath] = uigetfile({'*.jpk-force-map;*.jpk-qi-data',...
                    'Valid Types (*.jpk-force-map,*.jpk-qi-data)'},...
                    Title,'MultiSelect','on');
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
            
            for i=1:NRef
                TempID = sprintf('ReferenceMap-%i',i);
                TempRefFM{i,1} = ForceMap(MapFullFile{i},obj.ExperimentFolder,TempID);
            end
            
            obj.RefFM = TempRefFM;
            
        end
        
        function initialize_flags(obj)
            % Author: Manuel Rufin
            % Description: This function initializes flags which are used in the further course of the analysis. 

            %% Function body
            NFM = obj.NumForceMaps;
                obj.SMFSFlagDown.SelectFM = false(NFM,1);
                obj.SMFSFlagDown.PropertiesParameters = false(NFM,1);
                obj.SMFSFlagDown.Preprocessed = false(NFM,1);
                obj.SMFSFlagDown.Presorted = false(NFM,1);
                obj.SMFSFlagDown.AnalysedPreSelected = false(NFM,1);
                obj.SMFSFlagDown.AnalysedPostSelected = false(NFM,1);
                obj.SMFSFlagDown.NumForceCurves = false(NFM,1);

            if isempty(obj.FMFlag)
                NFM = obj.NumForceMaps;
                obj.FMFlag.FibrilAnalysis = false(NFM,1);
                obj.FMFlag.ForceMapAnalysis = false(NFM,1);
                obj.FMFlag.Preprocessed = false(NFM,1);
                obj.FMFlag.Grouping = false;
                NSPM = obj.NumSurfacePotentialMaps;
                obj.SPMFlag.FibrilAnalysis = false(NSPM,1);
                obj.SPMFlag.Grouping = false;
                % SMFSFlag
                obj.SMFSFlag.SelectFM = false(NFM,1);
                obj.SMFSFlag.PropertiesParameters = false(NFM,1);               
                obj.SMFSFlag.Preprocessed = false(NFM,1);
                obj.SMFSFlag.Presorted = false(NFM,1);
                obj.SMFSFlag.AnalysedPreSelected = false(NFM,1);
                obj.SMFSFlag.AnalysedPostSelected = false(NFM,1);
                obj.SMFSFlag.NumForceCurves = false(NFM,1);
                
                obj.DebugFlag.Plot= false(NFM,1);
                obj.CantileverTipFlag = false;
                if obj.NumCantileverTips > 0
                    obj.CantileverTipFlag = true;
                end
                obj.AssignedCantileverTips = false;
                if obj.NumCantileverTips == 1
                    obj.WhichTip = true(obj.NumForceMaps,1);
                    obj.AssignedCantileverTips = true;
                end
                
                obj.AssignedReferenceMaps = false;
                obj.ReferenceSlopeFlag.SetAllToValue = false;
                obj.ReferenceSlopeFlag.UserInput = false;
                obj.ReferenceSlopeFlag.FromRefFM = false;
                obj.ReferenceSlopeFlag.FromArea = false;
                obj.ReferenceSlopeFlag.AutomaticFibril = false;
                obj.ReferenceSlopeFlag.Automatic = false;
                if obj.NumReferenceForceMaps > 0
                    obj.ReferenceSlopeFlag.FromRefFM = true;
                end
            else
                PrevNFM = length(obj.FMFlag.FibrilAnalysis);
                NFM = obj.NumForceMaps;
                DiffFM = NFM - PrevNFM;
                obj.FMFlag.FibrilAnalysis(end+1:NFM) = false(DiffFM,1);
                obj.FMFlag.ForceMapAnalysis(end+1:NFM) = false(DiffFM,1);
                obj.FMFlag.Preprocessed(end+1:NFM) = false(DiffFM,1);
                obj.FMFlag.Grouping = false;
                PrevNSPM = length(obj.SPMFlag.FibrilAnalysis);
                NSPM = obj.NumSurfacePotentialMaps;
                DiffSPM = PrevNSPM - NSPM;
                obj.SPMFlag.FibrilAnalysis(end+1:NSPM) = false(DiffSPM,1);
                obj.SPMFlag.Grouping = false;
                obj.SMFSFlag.SelectFM(end+1:NFM) = false(DiffFM,1);
                obj.SMFSFlag.PropertiesParameters = false(DiffFM,1);  
                obj.SMFSFlag.Preprocessed(end+1:NFM) = false(DiffFM,1);
                obj.SMFSFlag.Presorted(end+1:NFM) = false(DiffFM,1);
                obj.SMFSFlag.AnalysedPreSelected = false(DiffFM,1);
                obj.SMFSFlag.AnalysedPostSelected = false(DiffFM,1);

                
                obj.CantileverTipFlag = false;
                if obj.NumCantileverTips > 0
                    obj.CantileverTipFlag = true;
                end
                obj.AssignedCantileverTips = false;
                if obj.NumCantileverTips == 1
                    obj.WhichTip = true(NFM,1);
                    obj.AssignedCantileverTips = true;
                end
                
                obj.AssignedReferenceMaps = false;
                obj.ReferenceSlopeFlag.SetAllToValue = false;
                obj.ReferenceSlopeFlag.UserInput = false;
                obj.ReferenceSlopeFlag.FromRefFM = false;
                obj.ReferenceSlopeFlag.FromArea = false;
                obj.ReferenceSlopeFlag.AutomaticFibril = false;
                obj.ReferenceSlopeFlag.Automatic = false;
                if obj.NumReferenceForceMaps > 0
                    obj.ReferenceSlopeFlag.FromRefFM = true;
                end
            end
        end
        
        function write_to_log_file(obj,Name, Value, StartEnd)
            % creates a Log-file that tracks all important analysis
            % parameters aswell as the current git branch and hash of the
            % current commit. All data analysis done this way should then
            % be reproducible down the line by downloading the
            % corresponding program version from GitHub and running with
            % the same parameters. 
            % This function is meant to be called multiple times in a
            % Experiment wrapper method (e.g. force_map_analysis_general())
            % Arguments:
            % Name: Name of the stored parameter in the log-file
            % Value: Value of the stored parameter (duh)
            % StartEnd = none, 'start', 'end'
            %       'start'... will create the logfile and write the
            %       GitInfo and the current date and time into the file
            %       header. The first function call should log the called
            %       wrapper-method
            %       none ... Leaving out the StartEnd argument will just
            %       append the Name Value pair to the existing file
            %       'end'... will append a last comment to the file,
            %       confirming that the analysis has terminated
            %       successfully and should therefore only be called as the
            %       very last Line in a wrapper-method.
            %
            %   NOTE: calling with 'end' will IGNORE whatever
            %   was written into Name and Value. However, you still need to
            %   Pass something in for those arguments, e.g.:
            %       obj.write_to_log_file('foo','bar','end')
            
            Current = what();
            cd(obj.ExperimentFolder)
            
            if nargin == 4
                if isequal(lower(StartEnd),'start')
                    DateTime = datestr(now,30);
                    obj.CurrentLogFile = strcat('AnalysisLog',DateTime,'.txt');
                    fid = fopen( obj.CurrentLogFile, 'wt' );
                    if fid == -1
                        error('Cannot open log file.');
                    end
                    fclose(fid);
                    
                    % Change to source code folder, get GitInfo and then
                    % change back to ExperimentFolder
                    Source = which('ForceMap');
                    Source = Source(1:end-11);
                    cd(Source)
                    GitInfo = getGitInfo();
                    cd(obj.ExperimentFolder);
                    
                    obj.write_to_log_file('Date and Time',datestr(now,0))
                    if isempty(GitInfo)
                        obj.write_to_log_file('GitHub Error','Could not determine GitInfo')
                    else
                        obj.write_to_log_file('Branch',GitInfo.branch)
                        obj.write_to_log_file('Hash of current commit',GitInfo.hash)
                        obj.write_to_log_file('Remote',GitInfo.remote)
                        obj.write_to_log_file('URL',GitInfo.url)
                    end
                elseif isequal(lower(StartEnd),'end')
                    obj.write_to_log_file('Analysis ended successfully',':D')
                    obj.CurrentLogFile = [];
                    cd(Current.path);
                    return
                end
            end
            
            if islogical(Value)
                if Value
                    Value = 'true';
                else
                    Value = 'false';
                end
            end
            Value = char(Value);
            
            fid = fopen(obj.CurrentLogFile, 'a');
            if fid == -1
                error('Cannot open log file.');
            end
            
            fprintf(fid, '%s: %s\n', Name, Value);
            
            fclose(fid);
            cd(Current.path);
            
        end
        
        function [PopUp,ClassIndex] = string_of_existing_class_instances(obj)
            
            k = 1;
            for i=1:obj.NumAFMImages
                PopUp{k} = sprintf('Image: %s%',obj.I{i}.Name);
                ClassIndex(k,:) = [1 i];
                k = k + 1;
            end
            for i=1:obj.NumForceMaps
                PopUp{k} = sprintf('Force Map: %s%',obj.FM{i}.Name);
                ClassIndex(k,:) = [2 i];
                k = k + 1;
            end
            for i=1:obj.NumReferenceForceMaps
                PopUp{k} = sprintf('Ref. Force Map: %s%',obj.RefFM{i}.Name);
                ClassIndex(k,:) = [3 i];
                k = k + 1;
            end
            for i=1:obj.NumCantileverTips
                PopUp{k} = sprintf('Cant. Tip Image: %s%',obj.CantileverTips{i}.Name);
                ClassIndex(k,:) = [4 i];
                k = k + 1;
            end
        end
        
        function Class = get_class_instance(obj,ClassIndex)
            
            if ClassIndex(1) == 1
                Class = obj.I{ClassIndex(2)};
            end
            if ClassIndex(1) == 2
                Class = obj.FM{ClassIndex(2)};
            end
            if ClassIndex(1) == 3
                Class = obj.RefFM{ClassIndex(2)};
            end
            if ClassIndex(1) == 4
                Class = obj.CantileverTips{ClassIndex(2)};
            end
            
        end
        
        function cast_volume_data_to_single_precision(obj)
            for i=1:obj.NumForceMaps
                obj.FM{i}.THApp = [];
                obj.FM{i}.THRet = [];
                for j=1:obj.FM{i}.NCurves
                    obj.FM{i}.App{j} = single(obj.FM{i}.App{j});
                    obj.FM{i}.Ret{j} = single(obj.FM{i}.Ret{j});
                    obj.FM{i}.HHApp{j} = single(obj.FM{i}.HHApp{j});
                    obj.FM{i}.HHRet{j} = single(obj.FM{i}.HHRet{j});
                    obj.FM{i}.BasedApp{j} = single(obj.FM{i}.BasedApp{j});
                    obj.FM{i}.BasedRet{j} = single(obj.FM{i}.BasedRet{j});
                end
            end
        end
        
    end
    methods(Static)
        % Static auxilary methods mainly for tip deconvolution (code by Orestis Andriotis)
        % and Experiment() loading
        
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
        
        function [Out, AppRetSwitch] = reference_slope_parser_gui(Methods)
            % Create figure
            left = 0.3;
            bottom = 0.4;
            width = 0.4;
            height = 0.15;
            h.f = figure('units','normalized','position',[left bottom width height],...
                'toolbar','none','menu','none');
            
            h.bg = uibuttongroup('Visible','off',...
                  'Position',[0.1 0.15 0.5 0.8]);
            
            % Create checkboxes
            c(1) = uicontrol(h.bg,'style','radiobutton','units','normalized',...
                'position',[0.1 0.8 0.8 1/7],'string','Set all RefSlopes to a chosen value');
            c(2) = uicontrol(h.bg,'style','radiobutton','units','normalized',...
                'position',[0.1 0.65 0.8 1/7],'string','Get RefSlopes from user input for each Force Map individually');
            c(3) = uicontrol(h.bg,'style','radiobutton','units','normalized',...
                'position',[0.1 0.5 0.8 1/7],'string','Get RefSlopes from Reference Force Maps');
            c(4) = uicontrol(h.bg,'style','radiobutton','units','normalized',...
                'position',[0.1 0.35 0.8 1/7],'string','Get RefSlopes from chosen area on Force Map');
            c(5) = uicontrol(h.bg,'style','radiobutton','units','normalized',...
                'position',[0.1 0.2 0.8 1/7],'string','Get RefSlope from automatically identified glass background around fibril');
            c(6) = uicontrol(h.bg,'style','radiobutton','units','normalized',...
                'position',[0.1 0.05 0.8 1/7],'string','Get RefSlope from automatically identified glass background around object');
            
            h.bg.SelectedObject = c(find(Methods));
            
            h.bg.Visible = 'on';
            
            
            h.bgAR = uibuttongroup('Visible','off',...
                  'Position',[0.65 0.15 0.3 0.8]);
            
            % Create checkboxes
            cAR(1) = uicontrol(h.bgAR,'style','radiobutton','units','normalized',...
                'position',[0.05 0.6 0.9 1/7],'string','From Approach');
            cAR(2) = uicontrol(h.bgAR,'style','radiobutton','units','normalized',...
                'position',[0.05 0.3 0.9 1/7],'string','From Retract');
            
            h.bgAR.Visible = 'on';
            
            % Create OK pushbutton
            h.p = uicontrol('style','pushbutton','units','normalized',...
                'position',[0.4 0.05 0.2 0.1],'string','OK',...
                'callback',@p_close);
            % Pushbutton callback
%             function p_call(varargin)
%                 vals = get(h.c,'Value');
%                 set(h.c(3),'Value',1)
%                 checked = find([vals{:}]);
%                 disp(checked)
%             end

            function p_close(varargin)
                vals = get(c,'Value');
                checked = find([vals{:}]);
                Out = checked;
                valsAR = get(cAR,'Value');
                checkedAR = find([valsAR{:}]);
                AppRetSwitch = checkedAR - 1;
                close(h.f)
            end
            uiwait(h.f)
        end
        
        function [FileTypes,OutStruct,IsValid,BigData] = constructor_user_input_parser(ExperimentName,OS)
            % Author: Manuel Rufin
            % Description:
            % This function opens a GUI window which allows to assign the raw data according to the jpk-data-type and the corresponding data paths are determined. 

            %% Function body
            IsValid = false;
            OutStruct = struct('FullFile',{cell(1,1),cell(1,1),cell(1,1),cell(1,1),cell(1,1)});
            FileTypes = zeros(1,5);
            
            Initial = what();
            h.LastFolder = Initial.path;
            
            % Create figure
            left = 0.2;
            bottom = 0.1;
            width = .6;
            height = 0.8;
            h.f = figure('Name',sprintf('%s: Data loader',ExperimentName),'units','normalized','position',[left bottom width height],...
                'toolbar','none','menu','none');
            
            % Create texttitles for filetypes
            c(1) = uicontrol(h.f,'style','text','units','normalized',...
                'position',[.05 .93 .2 .05],'string','Force/QI Maps',...
                'FontSize',14);
            c(2) = uicontrol(h.f,'style','text','units','normalized',...
                'position',[.05 .75 .2 .05],'string','Reference Force Maps',...
                'FontSize',14);
            c(3) = uicontrol(h.f,'style','text','units','normalized',...
                'position',[.05 .58 .2 .05],'string','AFM Image files',...
                'FontSize',14);
            c(4) = uicontrol(h.f,'style','text','units','normalized',...
                'position',[.05 .39 .2 .05],'string',"Surface Potential Maps",...
                'FontSize',14);
            c(5) = uicontrol(h.f,'style','text','units','normalized',...
                'position',[.05 .21 .2 .05],'string','Cantilever Tip data',...
                'FontSize',14);
            
            c(6) = uicontrol(h.f,'style','pushbutton','units','normalized',...
                'position',[.05 .88 .1 .05],'string','Open Browser',...
                'Callback',@open_browser);
            c(7) = uicontrol(h.f,'style','pushbutton','units','normalized',...
                'position',[.05 .7 .1 .05],'string','Open Browser',...
                'Callback',@open_browser);
            c(8) = uicontrol(h.f,'style','pushbutton','units','normalized',...
                'position',[.05 .53 .1 .05],'string','Open Browser',...
                'Callback',@open_browser);
            c(9) = uicontrol(h.f,'style','pushbutton','units','normalized',...
                'position',[.05 .34 .1 .05],'string',"Open Browser",...
                'Callback',@open_browser);
            c(10) = uicontrol(h.f,'style','pushbutton','units','normalized',...
                'position',[.05 .16 .1 .05],'string','Open Browser',...
                'Callback',@open_browser);
            
            c(11) = uicontrol(h.f,'style','pushbutton','units','normalized',...
                'position',[.15 .88 .1 .05],'string','Delete Selected',...
                'Callback',@delete_selected);
            c(12) = uicontrol(h.f,'style','pushbutton','units','normalized',...
                'position',[.15 .7 .1 .05],'string','Delete Selected',...
                'Callback',@delete_selected);
            c(13) = uicontrol(h.f,'style','pushbutton','units','normalized',...
                'position',[.15 .53 .1 .05],'string','Delete Selected',...
                'Callback',@delete_selected);
            c(14) = uicontrol(h.f,'style','pushbutton','units','normalized',...
                'position',[.15 .34 .1 .05],'string',"Delete Selected",...
                'Callback',@delete_selected);
            c(15) = uicontrol(h.f,'style','pushbutton','units','normalized',...
                'position',[.15 .16 .1 .05],'string','Delete Selected',...
                'Callback',@delete_selected);
            c(16) = uicontrol(h.f,'style','checkbox','units','normalized',...
                'position',[.3 .025 .2 .05],'string','Big Data mode',...
                'tooltip','Recomme2nded when processing >~100000 force curves',...
                'FontSize',18);
            
            HeightPos = [.82 .64 .46 .28 .1];
            
            for i=1:5
                h.ListBox(i) = uicontrol(h.f,...
                    'Style','listbox',...
                    'Max',1000000,'Min',1,...
                    'Units','normalized',...
                    'Position',[.3 HeightPos(i) .65 .16]);
            end
            
%             if isequal(upper(OS),'PCW')
                % Get back the java component associated to the axis
                % NB1: See §3.7.2 of Undocumented Secrets of Matlab Java Programming
                % NB2: or use findjobj, or javaObjectEDT for drop support onto other component types
                warning('off')
                jFrame = get(handle(h.f), 'JavaFrame');
                warning ('on')
                jAxis = jFrame.getAxisComponent();
                % Add listener for drop operations
                DropListener(jAxis, ... % The component to be observed
                    'DropFcn', @(s, e)onDrop(h.f, s, e)); % Function to call on drop operation
                h.DragNDrop = uicontrol(h.f,'style','text','units','normalized',...
                    'position',[.5 .025 .5 .05],...
                    'string',sprintf('Drag and Drop files into boxes\n (excl. to Windows and Linux) or load from browser'),...
                    'FontSize',16);
%             end
            
            % Create OK pushbutton
            h.p = uicontrol('style','pushbutton','units','normalized',...
                'position',[.05 .025 .1 .05],'string','Confirm',...
                'callback',@p_close);
            % Create cancel button
            h.Cancel = uicontrol('style','pushbutton','units','normalized',...
                'position',[.15 .025 .1 .05],'string','Cancel',...
                'callback',@pushed_cancel);
            
            function pushed_cancel(varargin)
                OutStruct = [];
                FileTypes = 'Cancel';
                IsValid = false;
                close(h.f)
            end
            
            function p_close(varargin)
                for i=1:5
                    if isempty(OutStruct(i).FullFile{1})
                        FileTypes(i) = 0;
                    else
                        FileTypes(i) = 1;
                    end
                end
                IsValid = true;
                BigData = c(16).Value;
                close(h.f)
            end
            
            function onDrop(fig, listener, evtArg) %#ok<INUSL>
                %[
                % Get back the dropped data
                Screen = get(groot);
                ScreenSize = Screen.ScreenSize;
                PointerPosition = get(0, 'PointerLocation');
                RelPointerPosition = PointerPosition./[ScreenSize(3) ScreenSize(4)];
                Index = which_listbox(RelPointerPosition);
                if isequal(Index,'invalid')
                    evtArg.DropComplete(false);
                    return
                end
                data = evtArg.GetTransferableData();
                
                % Is it transferable as a list of files
                if (data.IsTransferableAsFileList)
                    
                    TempFileCell = data.TransferAsFileList;
                    
                    TempFileCell = sort(TempFileCell);

                    k = 1;
                    DelIdx = [];
                    for i=1:length(TempFileCell)
                        SplitString = split(TempFileCell{i},filesep);
                        TempTempFile{i} = SplitString{end};
                        SplitName = split(TempTempFile{i},'.');
                        FileExtension = SplitName{end};
                        if ((Index == 1) || (Index == 2)) &&...
                                (isequal(FileExtension,'jpk-force-map') || isequal(FileExtension,'jpk-qi-data'))
                            % All Good
                        elseif ((Index == 3) || (Index == 5)) &&...
                                (isequal(FileExtension,'jpk') || isequal(FileExtension,'jpk-qi-image'))
                            % All Good
                        elseif (Index == 4)&&...
                                (isequal(FileExtension,'sdf'))
                            % All Good
                        else
                            DelIdx(k) = i;
                            k = k + 1;
                        end
                    end
                    
                    if ~isempty(DelIdx)
                        TempTempFile(DelIdx) = [];
                        TempFileCell(DelIdx) = [];
                    end
                    
                    k = 1;
                    if  ~iscell(TempTempFile)
                        TempFile{1} = TempTempFile;
                        TempFileCell{k} = fullfile(TempPath,TempFile{1});
                        k = k + 1;
                    else
                        TempFile = TempTempFile';
                    end
                    OldFiles = h.ListBox(Index).String;
                    if ~iscell(OldFiles)
                        NewFiles = TempFile;
                    else
                        NewFiles = cell(length(OldFiles)+length(TempFile),1);
                        for i=1:length(OldFiles)
                            NewFiles{i} = OldFiles{i};
                        end
                        for i=1:length(TempFile)
                            NewFiles{i+length(OldFiles)} = TempFile{i};
                        end
                    end
                    
                    if (length(OutStruct(Index).FullFile) == 1) && isempty(OutStruct(Index).FullFile{1})
                        OldLength = 0;
                    else
                        OldLength = length(OutStruct(Index).FullFile);
                    end
                    
                    for i=1:length(TempFile)
                        OutStruct(Index).FullFile{i+OldLength} = TempFileCell{i};
                    end
                    
                    
                    set(h.ListBox(Index),'String',NewFiles)
                    
                    % Indicate to the source that drop has completed
                    evtArg.DropComplete(true);
                    
                elseif (data.IsTransferableAsString)
                    
                    % Not interested
                    evtArg.DropComplete(false);
                    
                else
                    
                    % Not interested
                    evtArg.DropComplete(false);
                    
                end
                %]
            end
            
            function open_browser(varargin)
                
                ButtonPos = get(varargin{2}.Source,'InnerPosition');
                switch ButtonPos(2)
                    case .16
                        Index = 5;
                        AllowedFiles = {'*.jpk;*.jpk-qi-image',...
                            'Valid Types (*.jpk,*.jpk-qi-image)'};
                    case .34
                        Index = 4;
                        AllowedFiles = {'*.sdf',...
                            'Valid Types (*.sdf)'};
                    case .53
                        Index = 3;
                        AllowedFiles = {'*.jpk;*.jpk-qi-image',...
                            'Valid Types (*.jpk,*.jpk-qi-image)'};
                    case .7
                        Index = 2;
                        AllowedFiles = {'*.jpk-force-map;*.jpk-qi-data',...
                            'Valid Types (*.jpk-force-map,*.jpk-qi-data)'};
                    case .88
                        Index = 1;
                        AllowedFiles = {'*.jpk-force-map;*.jpk-qi-data',...
                            'Valid Types (*.jpk-force-map,*.jpk-qi-data)'};
                end
                
                current = what();
                cd(h.LastFolder)
                [TempTempFile,TempPath] = uigetfile(AllowedFiles,...
                    'MultiSelect','on');
                if isempty(TempTempFile)
                    return
                end
                k = 1;
                if  ~iscell(TempTempFile)
                    TempFile{1} = TempTempFile;
                    TempFileCell{k} = fullfile(TempPath,TempFile{1});
                    k = k + 1;
                else
                    TempFile = TempTempFile';
                    for i=1:length(TempFile)
                        TempFileCell{k} = fullfile(TempPath,TempFile{i});
                        k = k + 1;
                    end
                end
                h.LastFolder = TempPath;
                cd(current.path)
                OldFiles = h.ListBox(Index).String;
                if ~iscell(OldFiles)
                    NewFiles = TempFile;
                else
                    NewFiles = cell(length(OldFiles)+length(TempFile),1);
                    for i=1:length(OldFiles)
                        NewFiles{i} = OldFiles{i};
                    end
                    for i=1:length(TempFile)
                        NewFiles{i+length(OldFiles)} = TempFile{i};
                    end
                end
                
                if (length(OutStruct(Index).FullFile) == 1) && isempty(OutStruct(Index).FullFile{1})
                    OldLength = 0;
                else
                    OldLength = length(OutStruct(Index).FullFile);
                end
                
                for i=1:length(TempFile)
                    OutStruct(Index).FullFile{i+OldLength} = TempFileCell{i};
                end
                
                
                set(h.ListBox(Index),'String',NewFiles)
                OutStruct(Index).FullFile
            end
            
            function delete_selected(varargin)
                
                ButtonPos = get(varargin{2}.Source,'InnerPosition');
                switch ButtonPos(2)
                    case .16
                        Index = 5;
                    case .34
                        Index = 4;
                    case .53
                        Index = 3;
                    case .7
                        Index = 2;
                    case .88
                        Index = 1;
                end
                
                OldString = get(h.ListBox(Index),'String');
                DeleteIdx = get(h.ListBox(Index),'Value');
                OldString(DeleteIdx) = [];
                OutStruct(Index).FullFile(DeleteIdx) = [];
                
                set(h.ListBox(Index),'Value',1); 
                set(h.ListBox(Index),'String',OldString);    
                OutStruct(Index).FullFile
            end
            
            function FileTypeIndex = which_listbox(PPos)
                
                FigPos = get(h.f,'Position');
                
                for i=1:5
                    Pos = get(h.ListBox(i),'InnerPosition');
                    Pos(1) = Pos(1)*FigPos(3) + FigPos(1);
                    Pos(2) = Pos(2)*FigPos(4) + FigPos(2);
                    Pos(3) = Pos(1)*FigPos(3) + FigPos(1);
                    Pos(4) = Pos(2)*FigPos(4) + FigPos(2);
                    if (PPos(1) <= (Pos(1)+Pos(3)*FigPos(3))) &&...
                            (PPos(2) <= (Pos(2)+Pos(4)*FigPos(4))) &&...
                            (PPos(1) >= Pos(1)) &&...
                            (PPos(2) >= Pos(2))
                        FileTypeIndex = i;
                        return
                    else
                        FileTypeIndex = 'invalid';
                    end
                end
                
            end
            
            uiwait(h.f)
        end
        
    end
end