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
    
    % SMFS related
    properties
        SMFSResults
    end
    
    methods
        % constructor method and methods related with Experiment-file handling
        
        function obj = Experiment()
            
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
        
        %% SMFS section                     
        
        function SMFS_preprocessing(obj)
            % SMFS_preprocessing: A function to run a bundle of other 
            % typically required functions for further analysis
            
            
            % Output time and date for the dairy
            datetime('now')
            
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
            foldername='SMFS_preprocessing';    % for debugging                
            mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
            currpath=fullfile(obj.ExperimentFolder,foldername);
            cd(currpath);
            
            % force map loop
            for ii=1:obj.NumForceMaps   
            %for ii=139:obj.NumForceMaps % debugging
                if isequal(KeepFlagged,'Yes') && obj.SMFSFlag.Preprocessed(ii) == 1
                    continue
                end
                waitbar(ii/NLoop,h,sprintf('Preprocessing ForceMap %i/%i\nProcessing force curves',ii,NLoop));            
                obj.FM{ii}.fc_measurement_prop
                obj.FM{ii}.fc_sinoidal_fit
                obj.FM{ii}.fc_linear_fit
                obj.FM{ii}.fc_TipHeight_calculation
                obj.FM{ii}.fc_estimate_cp_hardsurface
                waitbar(ii/NLoop,h,sprintf('Preprocessing ForceMap %i/%i\nWrapping Up And Saving',ii,NLoop));
                                
                
                obj.SMFSFlag.Preprocessed(ii) = 1;
            end
            close(h);
        end
                  
      
        function SMFS_preparation(obj)
            % SMFS_presorting: This function allows to conduct an automated presorting of the force curves 
            % The function flags force curves and whole force maps that are
            % non-functionalize
            % Needed function: obj.preprocessing
            
            % Output time and date for the dairy
            datetime('now')
            
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
            for ii=1:obj.NumForceMaps
             %for ii=1:43 % Debugging
                if isequal(KeepFlagged,'Yes') && obj.SMFSFlag.Preprocessed(ii) == 1
                    continue
                end   
                waitbar(ii/NLoop,h,sprintf('Preprocessing ForceMap %i/%i\nProcessing force curves',ii,NLoop));
                waitbar(ii/NLoop,h,sprintf('Preprocessing ForceMap %i/%i\nWrapping Up And Saving',ii,NLoop));
                sprintf('Force Map No. %d of %d',ii,obj.NumForceMaps) % Gives current Force Map Position               
                obj.FM{ii}.fc_xLimit_idx
                obj.FM{ii}.fc_selection_threshold
                    if nnz(obj.FM{ii}.SMFSFlag.RetMinCrit)<20 % Only if more than 20 force curves fulfil the citeria the whole force map is considered successfully functionalized
                        obj.SMFSFlag.SelectFM(ii)=0;
                    else
                        obj.SMFSFlag.SelectFM(ii)=1;
                    end
                    obj.SMFSFlag.Presorted(ii) = 1;
            end
            close(h);
        end
      
        function SMFS_visual_selection(obj,XMin,XMax,YMin,YMax)
            % 
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
            % set(groot,'defaultFigureVisible','off')      
            set(groot,'defaultFigureVisible','on')  
            % Change into the Folder of Interest
            cd(obj.ExperimentFolder) % Move into the folder 
            % Create folders for saving the produced figures
            foldername='SMFS_visual_selection';    % Defines the folder name
            mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
            currpath=fullfile(obj.ExperimentFolder,foldername);
            cd(currpath); 
            
            % Loop over the imported force maps
            %for ii=1:obj.NumForceMaps
            for ii=1 % Debugging
               % Command window output
               sprintf('Force Map No. %d of %d',ii,obj.NumForceMaps) % Gives current Force Map Position
               % Run the chosen functions
            %   obj.FM{ii}.fc_visual_selection_all(XMin,XMax,YMin,YMax);  
               obj.FM{ii}.fc_visual_selection_flag_Selected(XMin,XMax,YMin,YMax,NumFcMax);
            %    obj.FM{ii}.fc_visual_selection_flag_Uncorrupt(obj,XMin,XMax,YMin,YMax,NumFcMax)
               %obj.save_experiment;        % Save immediately after each force curve
            end    
        end
      
        function SMFS_analysis(obj)
            % This function allows to analyse different force curve
            % criteria, i.e. pulling length, adhesion energy. 
           
            % Output time and date for the dairy
            datetime('now')
            
            % Figure visibility
            set(groot,'defaultFigureVisible','off')      
            % set(groot,'defaultFigureVisible','on') 
            %% Loop
            for hh=1:obj.NumForceMaps
            %for hh=4 % Debugging     
               sprintf('Force Map No. %d of %d',hh,obj.NumForceMaps) % Gives current Force Map Position   
               % Print force curves containing label for the pulling length
               % and colored area for the adhesion energy                              
               % Snap-in
               obj.FM{hh}.fc_snap_in_length_MAD
               % Pulling length
               obj.FM{hh}.fc_pulling_length_MAD
               % Maximum adhesion force
               obj.FM{hh}.fc_adh_force_max
               % Adhesion energy
               obj.FM{hh}.fc_adhesion_energy_idxlength                    
            end
            
        end
                            
        function SMFS_print_analysed_fc(obj,XMin,XMax,YMin,YMax,NumFcMax,NumFcUncorrupt,hh)
            %Furthermore, all analysed force curves are plotted and the determined
            % criteria are plotted for visual inspection
            % Input variable adaptation
            % IMPORTANT: Input variable NumFcMax - Only natural numbers are allowed
            % that result in natural numbers after square root extraction.
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
            foldername='FM_analysed';    % Defines the folder name
            mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
            currpath=fullfile(obj.ExperimentFolder,foldername);
            cd(currpath);
            %% loop
            for hh=1:obj.NumForceMaps
            %for hh=544:obj.NumForceMaps % Debugging
            %sprintf('Force Map No. %d of %d',hh,obj.NumForceMaps) % Gives current Force Map Position   
            % Determine needed input variable
               NumFcUncorrupt(hh)=nnz(obj.FM{hh}.SMFSFlag.Uncorrupt); % Determine the number of uncorrupted force curves     
               obj.FM{hh}.fc_print_properties(XMin,XMax,YMin,YMax,NumFcMax,NumFcUncorrupt,hh)         
            end
            obj.NumFcUncorrupt=NumFcUncorrupt;
        end
                
        function SMFS_print_raw(obj,XMin,XMax,YMin,YMax)
            % SMFS_print: A function to simply plot all force curves of all
            % force maps loaded and calssified based on the SMFS Flag
            % Needed function: obj.presorting
            
            % Show time and date for the dairy
            datetime('now')
            
            if nargin<2
                XMin= -inf;     % Limit of the X-axis in meters (m)  
                XMax= inf;      % Limit of the X-axis in meters (m)
                YMin= -inf;     % Limit of the Y-axis in Newtons (N)   
                YMax= inf;      % Limit of the Y-axis in Newtons (N)
            else
                
            end
            % Figure visibility
            set(groot,'defaultFigureVisible','off')      
            % set(groot,'defaultFigureVisible','on')  
            % Change into the Folder of Interest
            % Change into the Folder of Interest
            cd(obj.ExperimentFolder) % Move into the folder 
            % Create folders for saving the produced figures
            %foldername='FM_test';    % for debugging
            foldername='SMFS_print_raw';    % Defines the folder name
            mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
            currpath=fullfile(obj.ExperimentFolder,foldername);
            cd(currpath); 
            
            % Loop over the imported force maps
            for ii=1:obj.NumForceMaps
            %for ii=463 % Debugging
               % Command window output
               sprintf('Force Map No. %d of %d',ii,obj.NumForceMaps) % Gives current Force Map Position
               % Run the chosen functions
               obj.FM{ii}.fc_print_raw(XMin,XMax,YMin,YMax);     
            end    
        end
                         
        function SMFS_print_sort(obj,StartDate,EndDate,XMin,XMax,YMin,YMax)
            % SMFS_print_sort: A function to plot all force curves of all
            % force maps sorted by different properties 
            % Comment: Date format is: 'YYYY.MM.DD'
            
            % Figure visibility
            set(groot,'defaultFigureVisible','off')      
            % set(groot,'defaultFigureVisible','on')  
            % Change into the Folder of Interest
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
                XMin= -inf;     % Limit of the X-axis in meters (m)  
                XMax= inf;      % Limit of the X-axis in meters (m)
                YMin= -inf;     % Limit of the Y-axis in Newtons (N)   
                YMax= inf;      % Limit of the Y-axis in Newtons (N)
            elseif nargin<6    
                StartDate='0000.00.00';
                EndDate='2999.00.00';
            end
            
            % Output time and date for the dairy
            datetime('now')
            
            % Loop over the imported force maps
             for ii=1:obj.NumForceMaps
             % for ii=1:43
                 % Needed function               
                %if ~obj.SMFSFlag(ii)     % Selects all flagged 1 force maps
                %if obj.SMFSFlag(ii)     % Selects all flagged 0 force maps
                %    continue
                %end
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
                StartDateMod=strrep(StartDate,'.','');
                EndDateMod=strrep(EndDate,'.','');
                %foldername=append('FM_Flag',SMFSFlagConvert,'_',VelocityConvert,'_',obj.FM{ii}.Substrate,'_',obj.FM{ii}.EnvCond,'_',StartDateMod,'-',EndDateMod); % Defines the folder name
                foldername=append(obj.FM{ii}.Substrate,'_',obj.FM{ii}.EnvCond,'_',StartDateMod,'-',EndDateMod); % Defines the folder name 
                warning('off','all'); % To not showing the warning that the same folder is created each loop
                mkdir(foldername);
                warning('on','all');
                cd(foldername)         
               % Run the chosen functions    
               obj.FM{ii}.fc_print_raw(XMin,XMax,YMin,YMax)
               cd(obj.ExperimentFolder) % Move into the folder                             
            end 
            %obj.save_experiment        % Save immediately after each force curve
        end
        
        
        function SMFS_analysis_flag_status(obj)
            % SMFS_sinoidal_fit: This function allows to conduct an automated presorting of the force curves 
            % The function flags force curves and whole force maps that are
            % non-functionalize
            % Needed function: obj.preprocessing
         
            % Loop over the imported force maps
            
            %for ii=464:obj.NumForceMaps % Debugging
            
            % Find not processed force maps
            obj.SMFSFlagDown.Fit=find(~obj.SMFSFlag.Fit);
            obj.SMFSFlagDown.Preprocessed=find(~obj.SMFSFlag.Preprocessed);
            obj.SMFSFlagDown.Presorted=find(~obj.SMFSFlag.Presorted);
            obj.SMFSFlagDown.SelectFM=find(~obj.SMFSFlag.SelectFM);
          
            for ii=1:obj.NumForceMaps
            obj.FM{ii}.fc_flag_status          
            end
        end
                       
        
        function SMFS_analysis_dashboard(obj,ExtVelocityValue,RetVelocityValue,HoldingTimeValue,SubstrateValue,EnvCondValue,ChipCantValue,ChipboxValue,LinkerValue)
            % I all velocities should be selected use input variable: 0
            
            % Output time and date for the dairy
            datetime('now')
            % Write to log file
            obj.write_to_log_file('Function: SMFS_analysis_dashboard','Trial15','start')
            obj.write_to_log_file('Extend Velocity',num2str(ExtVelocityValue))
            obj.write_to_log_file('Retention Velocity',num2str(RetVelocityValue))
            obj.write_to_log_file('Holding Time',num2str(HoldingTimeValue))
            obj.write_to_log_file('Substrate',SubstrateValue)
            obj.write_to_log_file('Linker',EnvCondValue)
            obj.write_to_log_file('Chip',ChipCantValue)
            obj.write_to_log_file('Chipbox',ChipboxValue)
            obj.write_to_log_file('Linker',LinkerValue)
            obj.write_to_log_file('','','end')
            % Define variables
            jj=1;
            IdxArray=[];
            for ii=1:obj.NumForceMaps
            %% Debugging
            %for ii=1:10 % for debugging
                % sprintf('Force curve No. %d',ii) % Gives current Force curve
                % for debugging
                if ((obj.FM{ii}.ExtendVelocity==ExtVelocityValue || ExtVelocityValue==0) ...
                        && (obj.FM{ii}.RetractVelocity==RetVelocityValue || RetVelocityValue==0) ...
                        && (obj.FM{ii}.HoldingTime==HoldingTimeValue || HoldingTimeValue==-1) ...
                        && (strcmpi(obj.FM{ii}.Substrate,SubstrateValue) || strcmpi(SubstrateValue,'All')) ...
                        && (strcmpi(obj.FM{ii}.EnvCond,EnvCondValue) || strcmpi(EnvCondValue,'All')) ...
                        && (strcmpi(obj.FM{ii}.ChipCant,ChipCantValue) || strcmpi(ChipCantValue,'All')) ...
                        && (strcmpi(obj.FM{ii}.Chipbox,ChipboxValue) || strcmpi(ChipboxValue,'All')) ...
                        && (strcmpi(obj.FM{ii}.Linker,LinkerValue) || strcmpi(LinkerValue,'All')))
                    % Define variables for the if condition
                    IdxArray(jj,1)=ii;
                    % Adjust variable
                    jj=jj+1;
                end     
            end
            % If condition to handle an empty index array
            if isempty(IdxArray)
                return
            else    
            end
            % Define variables
            ExtVelocityValueStr=num2str(ExtVelocityValue);
            RetVelocityValueStr=num2str(RetVelocityValue);
            HoldingTimeValueStr=num2str(HoldingTimeValue);
            RGB11=[200 255 150]./255;
            % Change into the Folder of Interest
            cd(obj.ExperimentFolder) % Move into the folder
            % Create folders for saving the produced figures
            %foldername='FM_test';    % for debugging
            foldername='FM_analysis_dashboard';    % Defines the folder name
            mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
            currpath=fullfile(obj.ExperimentFolder,foldername);
            cd(currpath);
            % Define names
            figname=strcat(obj.ExperimentName,{'_'},ExtVelocityValueStr,{'_'},RetVelocityValueStr,{'_'},HoldingTimeValueStr,{'_'},SubstrateValue,{'_'},EnvCondValue,{'_'},ChipCantValue,{'_'},ChipboxValue,{'_'},LinkerValue);
            figname=char(figname);
            parttitle1='Parameters';
            %% Figure
            % Figure 1
            h_fig1=figure(1);
            h_fig1.Color='white'; % changes the background color of the figure
            h_fig1.Units='normalized'; % Defines the units
            h_fig1.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
            h_fig1.PaperOrientation='landscape';
            h_fig1.Name=figname;
            %% Plotting the tiles
            t = tiledlayout(3,3);
            t.TileSpacing = 'compact';
            t.Padding = 'compact';
            %t.TileSpacing = 'none'; % To reduce the spacing between the tiles
            %t.Padding = 'none'; % To reduce the padding of perimeter of a tile            
            % Define variables
            ConcateArray1=zeros(1,1);
            ConcateArray2=zeros(1,1);
            ConcateArray3=zeros(1,1);
            ConcateArray4=zeros(1,1);
            ConcateArray5=zeros(1,1);
            ConcateArray6=zeros(1,1);
            ConcateArray7=zeros(1,1);
                       
            % for loop
            for ff=1:length(IdxArray)
            %% Debugging
            %for ff=6 % for debugging
             sprintf('Force curve No. %d',ff) % Gives current Force curve
            % for debugging
                % Allocate data
                yAdhMaxApp=obj.FM{IdxArray(ff)}.AdhForceMaxRet;
                yAdhMaxRet=obj.FM{IdxArray(ff)}.AdhForceMaxApp;
                yAdhUnbinding=obj.FM{IdxArray(ff)}.AdhForceUnbinding;
                yAdhEneApp=obj.FM{IdxArray(ff)}.AppAdhEnergy_IdxMethod;
                yAdhEneRet=obj.FM{IdxArray(ff)}.RetAdhEnergy_IdxMethod;
                yPullingLength=obj.FM{IdxArray(ff)}.PullingLength;
                ySnapInLength=obj.FM{IdxArray(ff)}.SnapInLength;
                % Concatenate condition
                if ~isempty(yAdhMaxApp)
                    % Determine the number of rows per force map
                    ArrayLength=length(yAdhMaxApp); % Define the length of the array
                    row_start = ((ff-1) * ArrayLength) + 1; % Define the appropriate row start to append the new data
                    row_end   = ff * ArrayLength; % Define the appropriate row end to append the new data
                    % Concatenated data
                    ConcateArray1(row_start:row_end,:)=yAdhMaxApp'; % Append the new data into the concatenated vector
                else                   
                end                
                if ~isempty(yAdhMaxRet)
                    % Determine the number of rows per force map
                    ArrayLength=length(yAdhMaxRet); % Define the length of the array
                    row_start = ((ff-1) * ArrayLength) + 1; % Define the appropriate row start to append the new data
                    row_end   = ff * ArrayLength; % Define the appropriate row end to append the new data
                    % Concatenated data
                    ConcateArray2(row_start:row_end,:)=yAdhMaxRet'; % Append the new data into the concatenated vector
                else                  
                end
                if ~isempty(yAdhUnbinding)
                    % Determine the number of rows per force map
                    ArrayLength=length(yAdhUnbinding); % Define the length of the array
                    row_start = ((ff-1) * ArrayLength) + 1; % Define the appropriate row start to append the new data
                    row_end   = ff * ArrayLength; % Define the appropriate row end to append the new data
                    % Concatenated data
                    ConcateArray3(row_start:row_end,:)=yAdhUnbinding'; % Append the new data into the concatenated vector
                else                  
                end
                if ~isempty(yAdhEneApp)
                    % Determine the number of rows per force map
                    ArrayLength=length(yAdhEneApp); % Define the length of the array
                    row_start = ((ff-1) * ArrayLength) + 1; % Define the appropriate row start to append the new data
                    row_end   = ff * ArrayLength; % Define the appropriate row end to append the new data
                    % Concatenated data
                    ConcateArray4(row_start:row_end,:)=yAdhEneApp'; % Append the new data into the concatenated vector
                else                  
                end
                if ~isempty(yAdhEneRet)
                    % Determine the number of rows per force map
                    ArrayLength=length(yAdhEneRet); % Define the length of the array
                    row_start = ((ff-1) * ArrayLength) + 1; % Define the appropriate row start to append the new data
                    row_end   = ff * ArrayLength; % Define the appropriate row end to append the new data
                    % Concatenated data
                    ConcateArray5(row_start:row_end,:)=yAdhEneRet'; % Append the new data into the concatenated vector
                else                  
                end
                if ~isempty(yPullingLength)
                    % Determine the number of rows per force map
                    ArrayLength=length(yPullingLength); % Define the length of the array
                    row_start = ((ff-1) * ArrayLength) + 1; % Define the appropriate row start to append the new data
                    row_end   = ff * ArrayLength; % Define the appropriate row end to append the new data
                    % Concatenated data
                     ConcateArray6(row_start:row_end,:)=yPullingLength'; % Append the new data into the concatenated vector
                else                  
                end
                if ~isempty(ySnapInLength)
                    % Determine the number of rows per force map
                    ArrayLength=length(ySnapInLength); % Define the length of the array
                    row_start = ((ff-1) * ArrayLength) + 1; % Define the appropriate row start to append the new data
                    row_end   = ff * ArrayLength; % Define the appropriate row end to append the new data
                    % Concatenated data
                     ConcateArray7(row_start:row_end,:)=ySnapInLength'; % Append the new data into the concatenated vector
                else                  
                end  
                
                      
                % Statistics
                AdhMaxAppSelMean=mean(ConcateArray1);
                AdhMaxAppSelStd=std(ConcateArray1);
                AdhMaxRetSelMean=mean(ConcateArray2);
                AdhMaxRetSelStd=std(ConcateArray2);
                AdhMaxRetSelUnbindingMean=mean(ConcateArray3);
                AdhMaxRetSelUnbindingStd=std(ConcateArray3);
                AdhEneAppSelMean=mean(ConcateArray4);
                AdhEneAppSelStd=std(ConcateArray4);
                AdhEneRetSelMean=mean(ConcateArray5);
                AdhEneRetSelStd=std(ConcateArray5);
                PullLengthMedian=median(ConcateArray6);
                PullLengthMin=min(ConcateArray6);
                PullLengthMax=max(ConcateArray6);
                SnapInMedian=median(ConcateArray7);
                SnapInMin=min(ConcateArray7);
                SnapInMax=max(ConcateArray7);
                         
                % Plot
                % Tile 1 - Max. adhesion force approach
                ax1=nexttile(1);
                hold on
                plot(yAdhMaxApp,'o')
                % Title for each Subplot
                ti=title('Max. Adhesion Force Approach');                                     
                ti.Units='normalized'; % Set units to 'normalized'  
                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot    
                % Tile 2 - Max. adhesion force retract
                ax2=nexttile(2);
                hold on
                plot(yAdhMaxRet,'o')
                ti=title('Max. Adhesion Force Retraction');                                     
                ti.Units='normalized'; % Set units to 'normalized'  
                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot    
                % Tile 3 - Unbinding adhesion force (retract)
                ax3=nexttile(3);
                hold on
                plot(yAdhUnbinding,'o')
                ti=title('Unbinding Adhesion Force (Retraction)');                                     
                ti.Units='normalized'; % Set units to 'normalized'  
                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot    
                % Tile 4 - Adhesion energy approach
                ax4=nexttile(4);
                hold on
                plot(yAdhEneApp,'o')
                ti=title('Adhesion Energy Approach');                                     
                ti.Units='normalized'; % Set units to 'normalized'  
                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot    
                % Tile 5 - Adhesion energy retraction
                ax5=nexttile(5);
                hold on
                plot(yAdhEneRet,'o')
                ti=title('Adhesion Energy Retraction');                                     
                ti.Units='normalized'; % Set units to 'normalized'  
                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot    
                % Tile 6 - Pulling length
                ax6=nexttile(6);
                hold on
                plot(yPullingLength,'o')
                ti=title('Pulling Length');                                     
                ti.Units='normalized'; % Set units to 'normalized'  
                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot    
                % Tile 7 - Snap-in length
                ax7=nexttile(7);
                hold on
                plot(ySnapInLength,'o')
                ti=title('Snap-In Length');                                     
                ti.Units='normalized'; % Set units to 'normalized'  
                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot    
                % Tile 8
                ax8=nexttile(8);
                ax8.Color=RGB11;
                ax8.Box='on';
                ax8.LineWidth = 10;
                ax8.XTick=[];
                ax8.XTickLabel=[];
                ax8.YTick=[];
                ax8.YTickLabel=[];
                % Tile 9
                ax9=nexttile(9);
                ax9.Color=RGB11;
                ax9.Box='on';
                ax9.LineWidth = 10;
                ax9.XTick=[];
                ax9.XTickLabel=[];
                ax9.YTick=[];
                ax9.YTickLabel=[];
            end
            % Axes
            ax1.FontSize = 12;
            ax1.XLabel.String = 'Index (1)';
            ax1.XLabel.FontSize = 10;
            ax1.YLabel.String = 'Adhesion Force (N)';
            ax1.YLabel.FontSize = 10;
            ax2.FontSize = 12;
            ax2.XLabel.String = 'Index (1)';
            ax2.XLabel.FontSize = 10;
            ax2.YLabel.String = 'Adhesion Force (N)';
            ax2.YLabel.FontSize = 10;
            ax3.FontSize = 12;
            ax3.XLabel.String = 'Index (1)';
            ax3.XLabel.FontSize = 10;
            ax3.YLabel.String = 'Adhesion Force (N)';
            ax3.YLabel.FontSize = 10;
            ax4.FontSize = 12;
            ax4.XLabel.String = 'Index (1)';
            ax4.XLabel.FontSize = 10;
            ax4.YLabel.String = 'Adhesion Energy (J)';
            ax4.YLabel.FontSize = 10;
            ax5.FontSize = 12;
            ax5.XLabel.String = 'Index (1)';
            ax5.XLabel.FontSize = 10;
            ax5.YLabel.String = 'Adhesion Energy (J)';
            ax5.YLabel.FontSize = 10;
            ax6.FontSize = 12;
            ax6.XLabel.String = 'Index (1)';
            ax6.XLabel.FontSize = 10;
            ax6.YLabel.String = 'Pulling length (m)';
            ax6.YLabel.FontSize = 10;
            ax7.FontSize = 12;
            ax7.XLabel.String = 'Index (1)';
            ax7.XLabel.FontSize = 10;
            ax7.YLabel.String = 'Snap-in length (m)';
            ax7.YLabel.FontSize = 10;
            % Add text to plot
            % Tile 1 - Max. adhesion force approach
            ax1TextPos = [max(ax1.XLim) max(ax1.YLim)]; % Define the position in the plot
            partstr11a='Mean=';
            partstr12a=num2str(AdhMaxAppSelMean);
            partstr11b='Std=';
            partstr12b=num2str(AdhMaxAppSelStd);
            partstr13=' N';
            fullstr1a=strcat(partstr11a,partstr12a,partstr13); % Define the string that shall be shown in the plot
            fullstr1b=strcat(partstr11b,partstr12b,partstr13); % Define the string that shall be shown in the plot
            te1=text(ax1,ax1TextPos(1), ax1TextPos(2),{fullstr1a, fullstr1b}, 'VerticalAlignment','top', 'HorizontalAlignment','right');
            te1.FontSize = 12;
            te1.BackgroundColor=RGB11;
            % Tile 2 - Max. adhesion force retract
            ax2TextPos = [max(ax2.XLim) max(ax2.YLim)]-[diff(ax2.XLim) diff(ax2.YLim)]*0.01; % Define the position in the plot
            partstr21a='Mean=';
            partstr22a=num2str(AdhMaxRetSelMean);
            partstr21b='Std=';
            partstr22b=num2str(AdhMaxRetSelStd);
            fullstr2a=strcat(partstr21a,partstr22a,partstr13); % Define the string that shall be shown in the plot
            fullstr2b=strcat(partstr21b,partstr22b,partstr13); % Define the string that shall be shown in the plot
            te2=text(ax2,ax2TextPos(1), ax2TextPos(2),{fullstr2a, fullstr2b}, 'VerticalAlignment','top', 'HorizontalAlignment','right');
            te2.FontSize = 12;
            te2.BackgroundColor=RGB11;
            % Tile 3 - Unbinding adhesion force (retract)
            ax3TextPos = [max(ax3.XLim) max(ax3.YLim)]-[diff(ax3.XLim) diff(ax3.YLim)]*0.01; % Define the position in the plot
            partstr31a='Mean=';
            partstr32a=num2str(AdhMaxRetSelUnbindingMean);
            partstr31b='Std=';
            partstr32b=num2str(AdhMaxRetSelUnbindingStd);
            fullstr3a=strcat(partstr31a,partstr32a,partstr13); % Define the string that shall be shown in the plot
            fullstr3b=strcat(partstr31b,partstr32b,partstr13); % Define the string that shall be shown in the plot
            te3=text(ax3,ax3TextPos(1), ax3TextPos(2),{fullstr3a, fullstr3b}, 'VerticalAlignment','top', 'HorizontalAlignment','right');
            te3.FontSize = 12;
            te3.BackgroundColor=RGB11;
            % Tile 4 - Adhesion energy approach
            ax4TextPos = [max(ax4.XLim) max(ax4.YLim)]; % Define the position in the plot
            partstr41a='Mean=';
            partstr42a=num2str(AdhEneAppSelMean);        
            partstr41b='Std=';
            partstr42b=num2str(AdhEneAppSelStd);
            partstr43=' J';
            fullstr4a=strcat(partstr41a,partstr42a,partstr43); % Define the string that shall be shown in the plot
            fullstr4b=strcat(partstr41b,partstr42b,partstr43); % Define the string that shall be shown in the plot
            te4=text(ax4,ax4TextPos(1), ax4TextPos(2),{fullstr4a, fullstr4b}, 'VerticalAlignment','top', 'HorizontalAlignment','right');
            te4.FontSize = 12;
            te4.BackgroundColor=RGB11;
            % Tile 5 - Adhesion energy retraction
            ax5TextPos = [max(ax5.XLim) max(ax5.YLim)]; % Define the position in the plot
            partstr51a='Mean=';
            partstr52a=num2str(AdhEneRetSelMean);
            partstr51b='Std=';
            partstr52b=num2str(AdhEneRetSelStd);
            fullstr5a=strcat(partstr51a,partstr52a,partstr43); % Define the string that shall be shown in the plot
            fullstr5b=strcat(partstr51b,partstr52b,partstr43); % Define the string that shall be shown in the plot
            te5=text(ax5,ax5TextPos(1), ax5TextPos(2),{fullstr5a, fullstr5b}, 'VerticalAlignment','top', 'HorizontalAlignment','right');         
            te5.FontSize = 12;
            te5.BackgroundColor=RGB11;
            % Tile 6 - Pulling length
            ax6TextPos = [max(ax6.XLim) max(ax6.YLim)]; % Define the position in the plot
            partstr61a='Median=';
            partstr62a=num2str(PullLengthMedian);
            partstr61b='Min=';
            partstr72b=num2str(PullLengthMin);
            partstr71c='Max=';
            partstr72c=num2str(PullLengthMax);
            partstr63=' m';
            fullstr6a=strcat(partstr61a,partstr62a,partstr63); % Define the string that shall be shown in the plot
            fullstr6b=strcat(partstr61b,partstr72b,partstr63); % Define the string that shall be shown in the plot
            fullstr6c=strcat(partstr71c,partstr72c,partstr63); % Define the string that shall be shown in the plot
            te6=text(ax6,ax6TextPos(1), ax6TextPos(2),{fullstr6a, fullstr6b, fullstr6c}, 'VerticalAlignment','top', 'HorizontalAlignment','right');
            te6.FontSize = 12;
            te6.BackgroundColor=RGB11;
            % Tile 7 - Snap-in length
            ax7TextPos = [max(ax7.XLim) max(ax7.YLim)]-[diff(ax7.XLim) diff(ax7.YLim)]*0.01; % Define the position in the plot
            partstr71a='Median=';
            partstr72a=num2str(SnapInMedian);
            partstr71b='Min=';
            partstr72b=num2str(SnapInMin);
            partstr71c='Max=';
            partstr72c=num2str(SnapInMax);
            fullstr7a=strcat(partstr71a,partstr72a,partstr63); % Define the string that shall be shown in the plot
            fullstr7b=strcat(partstr71b,partstr72b,partstr63); % Define the string that shall be shown in the plot
            fullstr7c=strcat(partstr71c,partstr72c,partstr63); % Define the string that shall be shown in the plot
            te7=text(ax7,ax7TextPos(1), ax7TextPos(2),{fullstr7a, fullstr7b, fullstr7c}, 'VerticalAlignment','top', 'HorizontalAlignment','right');
            te7.FontSize = 12;
            te7.BackgroundColor=RGB11;
            % Tile 8
            ax8TextPos = [max(ax8.XLim)*0.5 max(ax8.YLim)*0.5];  % Define the position in the plot
            partstr81a='Extend Velocity: ';
            partstr82a=num2str(ExtVelocityValueStr);
            partstr83a='m*s^{-1}';
            fullstr8a=strcat(partstr81a,partstr82a,partstr83a); % Define the string that shall be shown in the plot
            partstr81b='Retract Velocity: ';
            partstr82b=num2str(RetVelocityValueStr);
            partstr83b='m*s^{-1}';
            fullstr8b=strcat(partstr81b,partstr82b,partstr83b); % Define the string that shall be shown in the plot
            partstr81c='Holding time: ';
            partstr82c=num2str(HoldingTimeValueStr);
            partstr83c='s';
            fullstr8c=strcat(partstr81c,partstr82c,partstr83c); % Define the string that shall be shown in the plot
            partstr81d='Substrate: ';
            partstr82d=num2str(SubstrateValue);
            fullstr8d=strcat(partstr81d,partstr82d); % Define the string that shall be shown in the plot          
            partstr81e='Medium: ';
            partstr82e=num2str(EnvCondValue);
            fullstr8e=strcat(partstr81e,partstr82e); % Define the string that shall be shown in the plot     
            partstr81f='Chip & Cantilever: ';
            partstr82f=num2str(ChipCantValue);
            fullstr8f=strcat(partstr81f,partstr82f); % Define the string that shall be shown in the plot      
            partstr81g='Chipbox: ';
            partstr82g=num2str(ChipboxValue);
            fullstr8g=strcat(partstr81g,partstr82g); % Define the string that shall be shown in the plot           
            partstr81h='Linker: ';
            partstr82h=num2str(ChipboxValue);
            fullstr8h=strcat(partstr81h,partstr82h); % Define the string that shall be shown in the plot   
            te8=text(ax8,ax8TextPos(1), ax8TextPos(2),{fullstr8a, fullstr8b, fullstr8c, fullstr8d, fullstr8e, fullstr8f, fullstr8g, fullstr8h}, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
            te8.FontSize = 18;
            % Tile 9     
            ax9TextPos = [max(ax9.XLim)*0.5 max(ax9.YLim)*0.5];  % Define the position in the plot
            fullstr9a='Number of force curves analysed: ';
            if ~isempty(yAdhMaxApp)
                 partstr91b='Max. Adhesion Force Approach: ';
                 partstr92b=num2str(length(yAdhMaxApp));
                 fullstr9b=strcat(partstr91b,partstr92b); % Define the string that shall be shown in the plot
                else                   
            end             
            if ~isempty(yAdhMaxRet)
                 partstr91c='Max. Adhesion Force Retraction:';
                 partstr92c=num2str(length(yAdhMaxRet));
                 fullstr9c=strcat(partstr91c,partstr92c); % Define the string that shall be shown in the plot
                else                   
            end 
            if ~isempty(yAdhUnbinding)
                 partstr91d='Unbinding Adhesion Force (Retraction): ';
                 partstr92d=num2str(length(yAdhUnbinding));
                 fullstr9d=strcat(partstr91d,partstr92d); % Define the string that shall be shown in the plot
                else                   
            end 
            if ~isempty(yAdhEneApp)
                 partstr91e='Adhesion Energy Approach: ';
                 partstr92e=num2str(length(yAdhEneApp));
                 fullstr9e=strcat(partstr91e,partstr92e); % Define the string that shall be shown in the plot
                else                   
            end 
            if ~isempty(yAdhEneRet)
                 partstr91d='Adhesion Energy Retraction: ';
                 partstr92d=num2str(length(yAdhEneRet));
                 fullstr9d=strcat(partstr91d,partstr92d); % Define the string that shall be shown in the plot
                else                   
            end 
            if ~isempty(yPullingLength)
                 partstr91e='Pulling Length: ';
                 partstr92e=num2str(length(yPullingLength));
                 fullstr9e=strcat(partstr91e,partstr92e); % Define the string that shall be shown in the plot
                else                   
            end 
             if ~isempty(ySnapInLength)
                 partstr91f='Snap-In Length: ';
                 partstr92f=num2str(length(ySnapInLength));
                 fullstr9f=strcat(partstr91f,partstr92f); % Define the string that shall be shown in the plot
                else                   
            end 
            te9=text(ax9,ax9TextPos(1), ax9TextPos(2),{fullstr9a, fullstr9b, fullstr9c, fullstr9d, fullstr9e, fullstr9f},'VerticalAlignment','middle', 'HorizontalAlignment','center');
            te9.FontSize = 18;            
            % Save figure
            %%% Define the name for the figure title
            partname='Dashboard';
            % fullname=sprintf('%s%s',figname,partname);
            fullname1=strcat(figname,{'_'},parttitle1,{'_'},partname);
            fullname1=char(fullname1);
            %%% Save the current figure in the current folder
            print(h_fig1,fullname1,'-dpng');
            
            %% SMFS Results structure 
            % Check entry
            if ~isempty(obj.SMFSResults)
                jj=length(obj.SMFSResults)+1;
            else
                jj=1;
            end
            
            % Append the new data into the concatenated array
            obj.SMFSResults{jj,1}.Data(1).AdhMaxApp=ConcateArray1;
            obj.SMFSResults{jj,1}.Data(1).AdhMaxRet= ConcateArray2;
            obj.SMFSResults{jj,1}.Data(1).AdhUnbinding=ConcateArray3;
            obj.SMFSResults{jj,1}.Data(1).AdhEneApp=ConcateArray4;
            obj.SMFSResults{jj,1}.Data(1).AdhEneRet=ConcateArray5;
            obj.SMFSResults{jj,1}.Data(1).AdhPullingLength=ConcateArray6;
            obj.SMFSResults{jj,1}.Properties(1).ExtendVelocity=ExtVelocityValue;
            obj.SMFSResults{jj,1}.Properties(1).RetractVelocity=RetVelocityValue;
            obj.SMFSResults{jj,1}.Properties(1).Substrate=SubstrateValue;
            obj.SMFSResults{jj,1}.Properties(1).Medium=EnvCondValue;
            obj.SMFSResults{jj,1}.Properties(1).ChipCantilever=ChipCantValue;
            obj.SMFSResults{jj,1}.Properties(1).Chipbox=ChipboxValue;
            obj.SMFSResults{jj,1}.Properties(1).Linker=LinkerValue;
            obj.SMFSResults{jj,1}.FMIndex=IdxArray;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppMean=AdhMaxAppSelMean;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppMean=AdhMaxAppSelMean;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxAppStd=AdhMaxAppSelStd;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetMean=AdhMaxRetSelMean;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetStd=AdhMaxRetSelStd;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingMean=AdhMaxRetSelUnbindingMean;
            obj.SMFSResults{jj,1}.Results(1).AdhMaxRetUnbindingStd=AdhMaxRetSelUnbindingStd;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppMean=AdhEneAppSelMean;
            obj.SMFSResults{jj,1}.Results(1).AdhEneAppStd=AdhEneAppSelStd;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetMean=AdhEneRetSelMean;
            obj.SMFSResults{jj,1}.Results(1).AdhEneRetStd=AdhEneRetSelStd;
            obj.SMFSResults{jj,1}.Results(1).PullLengthMedian=PullLengthMedian;
            obj.SMFSResults{jj,1}.Results(1).PullLengthMin=PullLengthMin;
            obj.SMFSResults{jj,1}.Results(1).PullLengthMax=PullLengthMax;
            obj.SMFSResults{jj,1}.Results(1).SnapInMedian=SnapInMedian;
            obj.SMFSResults{jj,1}.Results(1).SnapInMin=SnapInMin;
            obj.SMFSResults{jj,1}.Results(1).SnapInMax=SnapInMax;
                       
            %% House keeping
            close all
        end
        
      
        function SMFS_analysis_dashboard2(obj,ExtVelocityValue,RetVelocityValue,HoldingTimeValue,SubstrateValue,EnvCondValue,ChipCantValue,ChipboxValue,LinkerValue)
            % I all velocities should be selected use input variable: 0
            
            % Output time and date for the dairy
            datetime('now')
            
            % Define variables
            j=1;
            IdxArray=[];
            %for ii=1:obj.NumForceMaps
            %% Debugging
            for ii=1:10 % for debugging
                % sprintf('Force curve No. %d',ii) % Gives current Force curve
                % for debugging
                if ((obj.FM{ii}.ExtendVelocity==ExtVelocityValue || ExtVelocityValue==0) ...
                        && (obj.FM{ii}.RetractVelocity==RetVelocityValue || RetVelocityValue==0) ...
                        && (obj.FM{ii}.HoldingTime==HoldingTimeValue || HoldingTimeValue==0) ...
                        && (strcmpi(obj.FM{ii}.Substrate,SubstrateValue) || strcmpi(SubstrateValue,'All')) ...
                        && (strcmpi(obj.FM{ii}.EnvCond,EnvCondValue) || strcmpi(EnvCondValue,'All')) ...
                        && (strcmpi(obj.FM{ii}.ChipCant,ChipCantValue) || strcmpi(ChipCantValue,'All')) ...
                        && (strcmpi(obj.FM{ii}.Chipbox,ChipboxValue) || strcmpi(ChipboxValue,'All')) ...
                        && (strcmpi(obj.FM{ii}.Linker,LinkerValue) || strcmpi(LinkerValue,'All')))
                    % Define variables for the if condition
                    IdxArray(j,1)=ii;
                    % Adjust variable
                    j=j+1;
                end     
            end
            % If condition to handle an empty index array
            if isempty(IdxArray)
                return
            else    
            end
            % Define variables
            ExtVelocityValueStr=num2str(ExtVelocityValue);
            RetVelocityValueStr=num2str(RetVelocityValue);
            HoldingTimeValueStr=num2str(HoldingTimeValue);
            RGB11=[200 255 150]./255;
            % Change into the Folder of Interest
            cd(obj.ExperimentFolder) % Move into the folder
            % Create folders for saving the produced figures
            %foldername='FM_test';    % for debugging
            foldername='FM_analysis';    % Defines the folder name
            mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
            currpath=fullfile(obj.ExperimentFolder,foldername);
            cd(currpath);
            % Define names
            figname=strcat(obj.ExperimentName,{'_'},ExtVelocityValueStr,{'_'},RetVelocityValueStr,{'_'},HoldingTimeValueStr,{'_'},SubstrateValue,{'_'},EnvCondValue,{'_'},ChipCantValue,{'_'},ChipboxValue,{'_'},LinkerValue);
            figname=char(figname);
            parttitle1='Parameters';
            %% Figure
            % Figure 1
            h_fig1=figure(1);
            h_fig1.Color='white'; % changes the background color of the figure
            h_fig1.Units='normalized'; % Defines the units
            h_fig1.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
            h_fig1.PaperOrientation='landscape';
            h_fig1.Name=figname;
            %% Plotting the tiles
            t = tiledlayout(3,3);
            t.TileSpacing = 'compact';
            t.Padding = 'compact';
            %t.TileSpacing = 'none'; % To reduce the spacing between the tiles
            %t.Padding = 'none'; % To reduce the padding of perimeter of a tile            
            % Define variables
            ConcateArray1=zeros(1,1);
            ConcateArray2=zeros(1,1);
            ConcateArray3=zeros(1,1);
            ConcateArray4=zeros(1,1);
            ConcateArray5=zeros(1,1);
            ConcateArray6=zeros(1,1);
            ConcateArray7=zeros(1,1);
                       
            % for loop
            for ff=1:length(IdxArray)
                % Allocate data
                yAdhMaxRet=obj.FM{IdxArray(ff)}.AdhForceMaxApp;
                yAdhMaxApp=obj.FM{IdxArray(ff)}.AdhForceMaxRet;
                yAdhUnbinding=obj.FM{IdxArray(ff)}.AdhForceUnbinding;
                yAdhEneApp=obj.FM{IdxArray(ff)}.AppAdhEnergy_IdxMethod;
                yAdhEneRet=obj.FM{IdxArray(ff)}.RetAdhEnergy_IdxMethod;
                yPullingLength=obj.FM{IdxArray(ff)}.PullingLength;
                ySnapInLength=obj.FM{IdxArray(ff)}.SnapInLength;
                % Determine the number of rows per force map
                ArrayLength=length(yAdhMaxRet); % Define the length of the array
                row_start = ((ff-1) * ArrayLength) + 1; % Define the appropriate row start to append the new data
                row_end   = ff * ArrayLength; % Define the appropriate row end to append the new data
                % Concatenated data
                ConcateArray1(row_start:row_end,:)=yAdhMaxApp; % Append the new data into the concatenated vector
                ConcateArray2(row_start:row_end,:)=yAdhMaxRet; % Append the new data into the concatenated vector
                ConcateArray3(row_start:row_end,:)=yAdhUnbinding; % Append the new data into the concatenated vector
                ConcateArray4(row_start:row_end,:)=yAdhEneApp; % Append the new data into the concatenated vector
                ConcateArray5(row_start:row_end,:)=yAdhEneRet; % Append the new data into the concatenated vector
                ConcateArray6(row_start:row_end,:)=yPullingLength; % Append the new data into the concatenated vector
                ConcateArray7(row_start:row_end,:)=ySnapInLength; % Append the new data into the concatenated vector     
                % Statistics
                AdhMaxAppSelMean=mean(ConcateArray1);
                AdhMaxAppSelStd=std(ConcateArray1);
                AdhMaxRetSelMean=mean(ConcateArray2);
                AdhMaxRetSelStd=std(ConcateArray2);
                AdhMaxRetSelUnbindingMean=mean(ConcateArray3);
                AdhMaxRetSelUnbindingStd=std(ConcateArray3);
                AdhEneAppSelMean=mean(ConcateArray4);
                AdhEneAppSelStd=std(ConcateArray4);
                AdhEneRetSelMean=mean(ConcateArray5);
                AdhEneRetSelStd=std(ConcateArray5);
                PullLengthMedian=median(ConcateArray6);
                PullLengthMin=min(ConcateArray6);
                PullLengthMax=max(ConcateArray6);
                SnapInMedian=median(ConcateArray7);
                SnapInMin=min(ConcateArray7);
                SnapInMax=max(ConcateArray7);                        
                % Plot
                % Tile 1 - Max. adhesion force approach
                ax1=nexttile(1);
                hold on
                plot(yAdhMaxApp,'o')
                % Title for each Subplot
                ti=title('Max. Adhesion Force Approach');                                     
                ti.Units='normalized'; % Set units to 'normalized'  
                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot    
                % Tile 2 - Max. adhesion force retract
                ax2=nexttile(2);
                hold on
                plot(yAdhMaxRet,'o')
                ti=title('Max. Adhesion Force Retraction');                                     
                ti.Units='normalized'; % Set units to 'normalized'  
                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot    
                % Tile 3 - Unbinding adhesion force (retract)
                ax3=nexttile(3);
                hold on
                plot(yAdhUnbinding,'o')
                ti=title('Unbinding Adhesion Force (Retraction)');                                     
                ti.Units='normalized'; % Set units to 'normalized'  
                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot    
                % Tile 4 - Adhesion energy approach
                ax4=nexttile(4);
                hold on
                plot(yAdhEneApp,'o')
                ti=title('Adhesion Energy Approach');                                     
                ti.Units='normalized'; % Set units to 'normalized'  
                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot    
                % Tile 5 - Adhesion energy retraction
                ax5=nexttile(5);
                hold on
                plot(yAdhEneRet,'o')
                ti=title('Adhesion Energy Retraction');                                     
                ti.Units='normalized'; % Set units to 'normalized'  
                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot    
                % Tile 6 - Pulling length
                ax6=nexttile(6);
                hold on
                plot(yPullingLength,'o')
                ti=title('Pulling Length');                                     
                ti.Units='normalized'; % Set units to 'normalized'  
                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot    
                % Tile 7 - Snap-in length
                ax7=nexttile(7);
                hold on
                plot(ySnapInLength,'o')
                ti=title('Snap-In Length');                                     
                ti.Units='normalized'; % Set units to 'normalized'  
                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot    
                % Tile 8
                ax8=nexttile(8);
                ax8.Color=RGB11;
                ax8.Box='on';
                ax8.LineWidth = 10;
                ax8.XTick=[];
                ax8.XTickLabel=[];
                ax8.YTick=[];
                ax8.YTickLabel=[];
            end
            % Axes
            ax1.FontSize = 12;
            ax1.XLabel.String = 'Index (1)';
            ax1.XLabel.FontSize = 10;
            ax1.YLabel.String = 'Adhesion Force (N)';
            ax1.YLabel.FontSize = 10;
            ax2.FontSize = 12;
            ax2.XLabel.String = 'Index (1)';
            ax2.XLabel.FontSize = 10;
            ax2.YLabel.String = 'Adhesion Force (N)';
            ax2.YLabel.FontSize = 10;
            ax3.FontSize = 12;
            ax3.XLabel.String = 'Index (1)';
            ax3.XLabel.FontSize = 10;
            ax3.YLabel.String = 'Adhesion Force (N)';
            ax3.YLabel.FontSize = 10;
            ax4.FontSize = 12;
            ax4.XLabel.String = 'Index (1)';
            ax4.XLabel.FontSize = 10;
            ax4.YLabel.String = 'Adhesion Energy (J)';
            ax4.YLabel.FontSize = 10;
            ax5.FontSize = 12;
            ax5.XLabel.String = 'Index (1)';
            ax5.XLabel.FontSize = 10;
            ax5.YLabel.String = 'Adhesion Energy (J)';
            ax5.YLabel.FontSize = 10;
            ax6.FontSize = 12;
            ax6.XLabel.String = 'Index (1)';
            ax6.XLabel.FontSize = 10;
            ax6.YLabel.String = 'Pulling length (m)';
            ax6.YLabel.FontSize = 10;
            ax7.FontSize = 12;
            ax7.XLabel.String = 'Index (1)';
            ax7.XLabel.FontSize = 10;
            ax7.YLabel.String = 'Snap-in length (m)';
            ax7.YLabel.FontSize = 10;
            % Add text to plot
            % Tile 1 - Max. adhesion force approach
            ax1TextPos = [max(ax1.XLim) max(ax1.YLim)]; % Define the position in the plot
            partstr11a='Mean=';
            partstr12a=num2str(AdhMaxAppSelMean);
            partstr11b='Std=';
            partstr12b=num2str(AdhMaxAppSelStd);
            partstr13=' N';
            fullstr1a=strcat(partstr11a,partstr12a,partstr13); % Define the string that shall be shown in the plot
            fullstr1b=strcat(partstr11b,partstr12b,partstr13); % Define the string that shall be shown in the plot
            te1=text(ax1,ax1TextPos(1), ax1TextPos(2),{fullstr1a, fullstr1b}, 'VerticalAlignment','top', 'HorizontalAlignment','right');
            te1.FontSize = 12;
            te1.BackgroundColor=RGB11;
            % Tile 2 - Max. adhesion force retract
            ax2TextPos = [max(ax2.XLim) max(ax2.YLim)]-[diff(ax2.XLim) diff(ax2.YLim)]*0.01; % Define the position in the plot
            partstr21a='Mean=';
            partstr22a=num2str(AdhMaxRetSelMean);
            partstr21b='Std=';
            partstr22b=num2str(AdhMaxRetSelStd);
            fullstr2a=strcat(partstr21a,partstr22a,partstr13); % Define the string that shall be shown in the plot
            fullstr2b=strcat(partstr21b,partstr22b,partstr13); % Define the string that shall be shown in the plot
            te2=text(ax2,ax2TextPos(1), ax2TextPos(2),{fullstr2a, fullstr2b}, 'VerticalAlignment','top', 'HorizontalAlignment','right');
            te2.FontSize = 12;
            te2.BackgroundColor=RGB11;
            % Tile 3 - Unbinding adhesion force (retract)
            ax3TextPos = [max(ax3.XLim) max(ax3.YLim)]-[diff(ax3.XLim) diff(ax3.YLim)]*0.01; % Define the position in the plot
            partstr31a='Mean=';
            partstr32a=num2str(AdhMaxRetSelUnbindingMean);
            partstr31b='Std=';
            partstr32b=num2str(AdhMaxRetSelUnbindingStd);
            fullstr3a=strcat(partstr31a,partstr32a,partstr13); % Define the string that shall be shown in the plot
            fullstr3b=strcat(partstr31b,partstr32b,partstr13); % Define the string that shall be shown in the plot
            te3=text(ax3,ax3TextPos(1), ax3TextPos(2),{fullstr3a, fullstr3b}, 'VerticalAlignment','top', 'HorizontalAlignment','right');
            te3.FontSize = 12;
            te3.BackgroundColor=RGB11;
            % Tile 4 - Adhesion energy approach
            ax4TextPos = [max(ax4.XLim) max(ax4.YLim)]; % Define the position in the plot
            partstr41a='Mean=';
            partstr42a=num2str(AdhEneAppSelMean);        
            partstr41b='Std=';
            partstr42b=num2str(AdhEneAppSelStd);
            partstr43=' J';
            fullstr4a=strcat(partstr41a,partstr42a,partstr43); % Define the string that shall be shown in the plot
            fullstr4b=strcat(partstr41b,partstr42b,partstr43); % Define the string that shall be shown in the plot
            te4=text(ax4,ax4TextPos(1), ax4TextPos(2),{fullstr4a, fullstr4b}, 'VerticalAlignment','top', 'HorizontalAlignment','right');
            te4.FontSize = 12;
            te4.BackgroundColor=RGB11;
            % Tile 5 - Adhesion energy retraction
            ax5TextPos = [max(ax5.XLim) max(ax5.YLim)]; % Define the position in the plot
            partstr51a='Mean=';
            partstr52a=num2str(AdhEneRetSelMean);
            partstr51b='Std=';
            partstr52b=num2str(AdhEneRetSelStd);
            fullstr5a=strcat(partstr51a,partstr52a,partstr43); % Define the string that shall be shown in the plot
            fullstr5b=strcat(partstr51b,partstr52b,partstr43); % Define the string that shall be shown in the plot
            te5=text(ax5,ax5TextPos(1), ax5TextPos(2),{fullstr5a, fullstr5b}, 'VerticalAlignment','top', 'HorizontalAlignment','right');         
            te5.FontSize = 12;
            te5.BackgroundColor=RGB11;
            % Tile 6 - Pulling length
            ax6TextPos = [max(ax6.XLim) max(ax6.YLim)]; % Define the position in the plot
            partstr61a='Median=';
            partstr62a=num2str(PullLengthMedian);
            partstr61b='Min=';
            partstr72b=num2str(PullLengthMin);
            partstr71c='Max=';
            partstr72c=num2str(PullLengthMax);
            partstr63=' m';
            fullstr6a=strcat(partstr61a,partstr62a,partstr63); % Define the string that shall be shown in the plot
            fullstr6b=strcat(partstr61b,partstr72b,partstr63); % Define the string that shall be shown in the plot
            fullstr6c=strcat(partstr71c,partstr72c,partstr63); % Define the string that shall be shown in the plot
            te6=text(ax6,ax6TextPos(1), ax6TextPos(2),{fullstr6a, fullstr6b, fullstr6c}, 'VerticalAlignment','top', 'HorizontalAlignment','right');
            te6.FontSize = 12;
            te6.BackgroundColor=RGB11;
            % Tile 7 - Snap-in length
            ax7TextPos = [max(ax7.XLim) max(ax7.YLim)]-[diff(ax7.XLim) diff(ax7.YLim)]*0.01; % Define the position in the plot
            partstr71a='Median=';
            partstr72a=num2str(SnapInMedian);
            partstr71b='Min=';
            partstr72b=num2str(SnapInMin);
            partstr71c='Max=';
            partstr72c=num2str(SnapInMax);
            fullstr7a=strcat(partstr71a,partstr72a,partstr63); % Define the string that shall be shown in the plot
            fullstr7b=strcat(partstr71b,partstr72b,partstr63); % Define the string that shall be shown in the plot
            fullstr7c=strcat(partstr71c,partstr72c,partstr63); % Define the string that shall be shown in the plot
            te7=text(ax7,ax7TextPos(1), ax7TextPos(2),{fullstr7a, fullstr7b, fullstr7c}, 'VerticalAlignment','top', 'HorizontalAlignment','right');
            te7.FontSize = 12;
            te7.BackgroundColor=RGB11;
            % Tile 8
            ax8TextPos = [max(ax8.XLim)*0.5 max(ax8.YLim)*0.5];  % Define the position in the plot
            partstr81a='Extend Velocity: ';
            partstr82a=num2str(ExtVelocityValueStr);
            partstr83a='m*s^{-1}';
            fullstr8a=strcat(partstr81a,partstr82a,partstr83a); % Define the string that shall be shown in the plot
            partstr81b='Retract Velocity: ';
            partstr82b=num2str(RetVelocityValueStr);
            partstr83b='m*s^{-1}';
            fullstr8b=strcat(partstr81b,partstr82b,partstr83b); % Define the string that shall be shown in the plot
            partstr81c='Holding time: ';
            partstr82c=num2str(HoldingTimeValueStr);
            partstr83c='s';
            fullstr8c=strcat(partstr81c,partstr82c,partstr83c); % Define the string that shall be shown in the plot
            partstr81d='Substrate: ';
            partstr82d=num2str(SubstrateValue);
            fullstr8d=strcat(partstr81d,partstr82d); % Define the string that shall be shown in the plot          
            partstr81e='Medium: ';
            partstr82e=num2str(EnvCondValue);
            fullstr8e=strcat(partstr81e,partstr82e); % Define the string that shall be shown in the plot     
            partstr81f='Chip & Cantilever: ';
            partstr82f=num2str(ChipCantValue);
            fullstr8f=strcat(partstr81f,partstr82f); % Define the string that shall be shown in the plot      
            partstr81g='Chipbox: ';
            partstr82g=num2str(ChipboxValue);
            fullstr8g=strcat(partstr81g,partstr82g); % Define the string that shall be shown in the plot           
            partstr81h='Linker: ';
            partstr82h=num2str(ChipboxValue);
            fullstr8h=strcat(partstr81h,partstr82h); % Define the string that shall be shown in the plot   
            partstr81i='Number of force curves analysed: ';
            partstr82i=num2str(length(ConcateArray1));
            fullstr8i=strcat(partstr81i,partstr82i); % Define the string that shall be shown in the plot
            partstr81j='Number of force maps analysed: ';
            partstr82j=num2str(length(IdxArray));
            fullstr8j=strcat(partstr81j,partstr82j); % Define the string that shall be shown in the plot
            te8=text(ax8,ax8TextPos(1), ax8TextPos(2),{fullstr8a, fullstr8b, fullstr8c, fullstr8d, fullstr8e, fullstr8f, fullstr8g, fullstr8h, fullstr8i, fullstr8j}, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
            te8.FontSize = 18;
    
            % Save figure
            %%% Define the name for the figure title
            partname='Dashboard';
            % fullname=sprintf('%s%s',figname,partname);
            fullname1=strcat(figname,{'_'},parttitle1,{'_'},partname);
            fullname1=char(fullname1);
            %%% Save the current figure in the current folder
            print(h_fig1,fullname1,'-dpng');
            
            %% SMFSResults Structure
            % Check entry 
            if length(obj.SMFSResults)>0
                jj=length(obj.SMFSResults)+1;
            else
                jj=1;
            end
            % Append the new data into the concatenated array
            obj.SMFSResults{jj,1}.Data(1).AdhMaxApp=ConcateArray1;
            obj.SMFSResults{jj,1}.Data(1).AdhMaxRet=ConcateArray2;
            obj.SMFSResults{jj,1}.Data(1).AdhUnbinding=ConcateArray3;
            obj.SMFSResults{jj,1}.Data(1).AdhEneApp=ConcateArray4;
            obj.SMFSResults{jj,1}.Data(1).AdhEneRet=ConcateArray5;
            obj.SMFSResults{jj,1}.Data(1).PullingLength=ConcateArray6;
            obj.SMFSResults{jj,1}.Data(1).SnapInLength=ConcateArray7;
            obj.SMFSResults{jj,1}.Properties(1).ExtendVelocity=ExtVelocityValue;
            obj.SMFSResults{jj,1}.Properties(1).RetractVelocity=RetVelocityValue;
            obj.SMFSResults{jj,1}.Properties(1).HoldingTime=HoldingTimeValue;
            obj.SMFSResults{jj,1}.Properties(1).Substrate=SubstrateValue;
            obj.SMFSResults{jj,1}.Properties(1).Medium=EnvCondValue;
            obj.SMFSResults{jj,1}.Properties(1).ChipCantilever=ChipCantValue;
            obj.SMFSResults{jj,1}.Properties(1).Chipbox=ChipboxValue;
            obj.SMFSResults{jj,1}.Properties(1).Linker=LinkerValue;
            obj.SMFSResults{jj,1}.FMIndex=IdxArray;
            
            
            %% House keeping
            close all
        end
        
       
        
        function SMFS_analysis_selction_fit(obj,VelocityValue,SubstrateValue,EnvCondValue,ChipCantValue,ChipboxValue,LinkerValue)
                
            % Output time and date for the dairy
            datetime('now')
            
                % Define colors
                RGB8=[80 200 204]./255; % Turquoise
                RGB9=[0 181 26]./255; % Green
                RGB10=[69 22 113]./255; % Violet
                RGB11=[124 223 124]./255; % Light Green
                % Loop
                for ii=1:obj.NumForceMaps
                    if (strcmpi(obj.FM{ii}.VelocityConvert,VelocityValue) || strcmpi(VelocityValue,'All')) ...
                            && (strcmpi(obj.FM{ii}.Substrate,SubstrateValue) || strcmpi(SubstrateValue,'All')) ...
                            && (strcmpi(obj.FM{ii}.EnvCond,EnvCondValue) || strcmpi(EnvCondValue,'All')) ...
                            && (strcmpi(obj.FM{ii}.ChipCant,ChipCantValue) || strcmpi(ChipCantValue,'All')) ...
                            && (strcmpi(obj.FM{ii}.Chipbox,ChipboxValue) || strcmpi(ChipboxValue,'All')) ...
                            && (strcmpi(obj.FM{ii}.Linker,LinkerValue) || strcmpi(LinkerValue,'All'))
                        % Define variables for the if condition 
                        IdxArray(ii,1)=ii;
                        IdxNonzero=find(IdxArray,1,'first');
                        
                        if ~isempty(IdxArray) && IdxNonzero==ii
                            obj.FM{ii}.PullingLength(obj.FM{ii}.PullingLength==0)=nan;
                            ConcatArrayPullingLength1=obj.FM{ii}.PullingLength(:);
                            ConcatArrayPullingLength2=obj.FM{ii}.PullingLength(:);
                            obj.FM{ii}.RetAdhEnergy_IdxMethod(obj.FM{ii}.RetAdhEnergy_IdxMethod==0)=nan;
                            ConcatArrayRetAdhEnergy1=obj.FM{ii}.RetAdhEnergy_IdxMethod(:);
                            ConcatArrayRetAdhEnergy2=obj.FM{ii}.RetAdhEnergy_IdxMethod(:);
                            ConcatArrayRetForce=obj.FM{ii}.MinRet;
                        else
                            obj.FM{ii}.PullingLength(obj.FM{ii}.PullingLength==0)=nan;
                            ConcatArrayPullingLength1=horzcat(ConcatArrayPullingLength1,obj.FM{ii}.PullingLength(:));
                            ConcatArrayPullingLength2=vertcat(ConcatArrayPullingLength2,obj.FM{ii}.PullingLength(:));
                            obj.FM{ii}.RetAdhEnergy_IdxMethod(obj.FM{ii}.RetAdhEnergy_IdxMethod==0)=nan;
                            ConcatArrayRetAdhEnergy1=horzcat(ConcatArrayRetAdhEnergy1,obj.FM{ii}.RetAdhEnergy_IdxMethod(:));
                            ConcatArrayRetAdhEnergy2=vertcat(ConcatArrayRetAdhEnergy2,obj.FM{ii}.RetAdhEnergy_IdxMethod(:));
                            ConcatArrayRetForce=obj.FM{ii}.MinRet;
                        end
                    end
                end
                ConcatArrayRetAdhEnergy1=ConcatArrayRetAdhEnergy1*-1; % Multiply with -1 to bring the ahdesion energy values into quadrant I
                ConcatArrayRetAdhEnergy2=ConcatArrayRetAdhEnergy2*-1; % Multiply with -1 to bring the ahdesion energy values into quadrant I
                
                % Change into the Folder of Interest
                cd(obj.ExperimentFolder) % Move into the folder
                % Create folders for saving the produced figures
                foldername='FM_test';    % for debugging
                %foldername='FM_analysis_fit';    % Defines the folder name
                mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
                currpath=fullfile(obj.ExperimentFolder,foldername);
                cd(currpath);
                
                % Define names
                figname=strcat(obj.ExperimentName,{'_'},VelocityValue,{'_'},SubstrateValue,{'_'},EnvCondValue,{'_'},ChipCantValue,{'_'},ChipboxValue,{'_'},LinkerValue,{'_'});
                figname=char(figname);
                parttitle1='PullingLength';
                parttitle2='AdhesionEnergy';
                fulltitle1=strcat(parttitle1,{'_'},figname);
                fulltitle2=strcat(parttitle2,{'_'},figname);
                
                % Define variables
                pdname='Lognormal'; % Chooses the probability distribution 
                
                % Calculate statistic values
                MaxPullingLength=max(ConcatArrayPullingLength2);
                MinPullingLength=min(ConcatArrayPullingLength2);
                MeanArrayRetForce=mean(ConcatArrayRetForce);
                
                NumFC=sum(~isnan(ConcatArrayRetAdhEnergy2))
                             
                % Fitting the data
                distLognormalRetAdh=fitdist(ConcatArrayRetAdhEnergy2,pdname); % Fits probability distribution object to data
                distLognormalPullLen=fitdist(ConcatArrayPullingLength2,pdname); % Fits probability distribution object to data
                % Calculate the mode (peak) of the log-normal distribution 
                LogNormalModeRetAdh=exp(distLognormalRetAdh.mu-distLognormalRetAdh.sigma^2);
                LogNormalModePullLen=exp(distLognormalPullLen.mu-distLognormalPullLen.sigma^2);
                % Calculate the mode (peak) of the log-normal distribution                
                LogNormalMedianPullLen=exp(distLognormalPullLen.mu);
                % Calculate the variance of the log-normal distribution 
                LogNormalVarRetAdh=exp((distLognormalRetAdh.sigma^2)-1)*exp(2*distLognormalRetAdh.mu+distLognormalRetAdh.sigma^2);
                LogNormalVarPullLen=exp((distLognormalPullLen.sigma^2)-1)*exp(2*distLognormalPullLen.mu+distLognormalPullLen.sigma^2);
                % Determine the probability density function for the plot
                pdfLognormalRetAdh=pdf(distLognormalRetAdh,linspace(0,max(ConcatArrayRetAdhEnergy2),1000));    
                pdfLognormalPullLen=pdf(distLognormalPullLen,linspace(0,max(ConcatArrayPullingLength2),1000));    
                % pdfLognormal=distLognormal.pdf(linspace(0,max(ConcatArrayRetAdhEnergy2),1000))
         
%                 % Figure 1
%                 h_fig1=figure(1);
%                 h_fig1.Color='white'; % changes the background color of the figure
%                 h_fig1.Units='normalized'; % Defines the units
%                 h_fig1.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
%                 h_fig1.PaperOrientation='landscape';
%                 h_fig1.Name=figname;
%                 % Histogram
%                  % Histogram
%                 hold on
%                 yyaxis left               
%                 h1=histogram(ConcatArrayPullingLength2(~isnan(ConcatArrayPullingLength2)));
%                 yyaxis right
%                 plot(linspace(0,max(ConcatArrayPullingLength2),1000),pdfLognormalPullLen)
%                 xlim([0 max(ConcatArrayPullingLength2)])
%               %  h3.BinWidth=20e-9;
%                 h1.FaceAlpha=1;
%                 h1.FaceColor='b';
%                 h1.EdgeColor='k';
%                 % title 
%                 t1=title(fulltitle1);
%                 % axes
%                 ax1=gca;
%                 ax1.FontSize = 16;
%                 ax1.XLabel.String = 'PullingLength (m)';
%                 ax1.XLabel.FontSize = 20;
%                 ax1.YLabel.String = 'Frequency count (1)';
%                 ax1.YLabel.FontSize = 20;               
%                 hold on;
%                 % text 
%                 NE = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.2;   % Define the position in the plot    
%                 SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.1;
%                 partstrA1='median=';
%                 partstrA2=num2str( LogNormalMedianPullLen);
%                 partstrB1='min=';
%                 partstrB2=num2str(MinPullingLength);
%                 partstrC1='max=';
%                 partstrC2=num2str(MaxPullingLength);
%                 fullstrA=strcat(partstrA1,partstrA2); % Define the string that shall be shown in the plot
%                 fullstrB=strcat(partstrB1,partstrB2); % Define the string that shall be shown in the plot
%                 fullstrC=strcat(partstrC1,partstrC2); % Define the string that shall be shown in the plot
%                 fullstrBC=strcat(fullstrB,{'  '},fullstrC);
%                 te11=text(NE(1), NE(2),fullstrA, 'VerticalAlignment','top', 'HorizontalAlignment','right');
%                 te11.FontSize = 22;
%                 te2=text(SE(1), SE(2),fullstrBC, 'VerticalAlignment','top', 'HorizontalAlignment','right');
%                 te2.FontSize = 22;
%                 % Save figure
%                 %%% Define the name for the figure title
%                 partname='histo';
%                 % fullname=sprintf('%s%s',figname,partname);
%                 fullname3=strcat(figname,{'_'},parttitle1,{'_'},partname);
%                 fullname3=char(fullname3);
%                 %%% Save the current figure in the current folder
%                 print(gcf,fullname3,'-dpng');
%                 
%                 % Figure 2
%                 h_fig2=figure(2);
%                 h_fig2.Color='white'; % changes the background color of the figure
%                 h_fig2.Units='normalized'; % Defines the units
%                 h_fig2.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
%                 h_fig2.PaperOrientation='landscape';
%                 h_fig2.Name=figname;
%                 % Histogram
%                 hold on
%                 yyaxis left               
%                 h2=histogram(ConcatArrayRetAdhEnergy2(~isnan(ConcatArrayRetAdhEnergy2))); 
%                 h2.BinWidth=20e-19;
%                 yyaxis right
%                 plot(linspace(0,max(ConcatArrayRetAdhEnergy2),1000),pdfLognormalRetAdh)
%                 xlim([0 max(ConcatArrayRetAdhEnergy2)])               
%                 h2.FaceAlpha=1;
%                 h2.FaceColor='g';
%                 h2.EdgeColor='k';
%                 % title
%                 t2=title(fulltitle2);
%                 % axes
%                 ax2=gca;
%                 ax2.FontSize = 16;
%                 %ax4.YScale='log';
%                 ax2.XLabel.String = 'Adhesion Energy (J)';
%                 ax2.XLabel.FontSize = 20;
%                 ax2.YLabel.String = 'Frequency count (1)';
%                 ax2.YLabel.FontSize = 20;               
%                 hold on;
%                 % text 
%                 NE = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.2;   % Define the position in the plot
%                 SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.1;   % Define the position in the plot       
%                 partstrA1='Mode=';
%                 partstrA2=num2str(LogNormalModeRetAdh);                                 
%                 partstrB1='Variance=';
%                 partstrB2=num2str(LogNormalVarRetAdh);
%                 partstrC1='Num fc analysed=';
%                 partstrC2=num2str(NumFC);           
%                 fullstrA=strcat(partstrA1,partstrA2); % Define the string that shall be shown in the plot
%                 fullstrB=strcat(partstrB1,partstrB2); % Define the string that shall be shown in the plot               
%                 fullstrAB=strcat(fullstrA,{'  '},fullstrB);
%                 fullstrC=strcat(partstrC1,partstrC2); % Define the string that shall be shown in the plot               
%                 te21=text(NE(1), NE(2),fullstrAB, 'VerticalAlignment','top', 'HorizontalAlignment','right');
%                 te21.FontSize = 22;
%                 te22=text(SE(1), SE(2),fullstrC, 'VerticalAlignment','top', 'HorizontalAlignment','right');
%                 te22.FontSize = 22;
%                 
                % Fine Figure
                h_fig10=figure(10);
                h_fig10.Color='white'; % changes the background color of the figure
                h_fig10.Units='normalized'; % Defines the units
                h_fig10.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
                h_fig10.PaperOrientation='landscape';
                h_fig10.Name=figname;
                % Histogram
                hold on
                yyaxis left
                h10=histogram(ConcatArrayRetAdhEnergy2(~isnan(ConcatArrayRetAdhEnergy2)),'LineWidth',6); 
                %h10.BinWidth=5e-19;
                h10.BinWidth=15e-19;
                yyaxis right
                ax10=gca
                plot(linspace(0,max(ConcatArrayRetAdhEnergy2),1000),pdfLognormalRetAdh,'Color',RGB10,'LineWidth',6)
                xlim([0 max(ConcatArrayRetAdhEnergy2)])               
                h10.FaceAlpha=1;
                h10.FaceColor=RGB11;
                h10.EdgeColor='k';
                % legend
                le10=legend('Frequency count','Log normal fit','Location','best');
                le10.FontSize = 48;      
                le10.EdgeColor='w';
                % axes
                yyaxis left
                ax10.FontSize = 48;
                ax10.LineWidth = 5;
                ax10.XLim =[1.5e-17 7e-17];
                ax10.YAxis(1).Color='k';
                ax10.YAxis(2).Color='w';
                ax10.YAxis(2).TickLabel=[];
                ax10.TickDir='out';
                ax10.XLabel.String = 'Adhesion Energy (J)';
                ax10.XLabel.FontSize = 52;
                ax10.YLabel.String = 'Frequency count (1)';
                ax10.YLabel.FontSize = 52;               
                hold on;
                        
                % Save figure
                %%% Define the name for the figure title
                partname='histo';
                % fullname=sprintf('%s%s',figname,partname);
                fullname4=strcat(figname,{'_'},parttitle2,{'_'},partname);
                fullname4=char(fullname4);
                %%% Save the current figure in the current folder
                print(gcf,fullname4,'-dpng');
                %% House keeping
                close all
        end
       
        
        function SMFS_boxplot_pulllength(obj,XMin,XMax,YMin,YMax) % fc ... force curve
            %
            % Output time and date for the dairy
            datetime('now')
            
            if nargin < 2
                XMin= -inf;
                XMax= inf;
                YMin= -inf;
                YMax= inf;
            end
            % Define variables for the figure name
%             VelocityConvert=num2str(obj.Velocity*1e+9); % Convert into nm
%             % Classification criteria
%             figname=strcat(obj.DateAdapt,{'_'},obj.TimeAdapt,{'_'},obj.ID,{'_'},obj.Substrate,{'_'},obj.EnvCond,{'_'},VelocityConvert,{'_'},obj.Chipbox,{'_'},obj.ChipCant);
%             figname=char(figname);
%             % Define variables for the plot loop
             
            %  figname=strcat(obj.DateAdapt,{'_'},obj.TimeAdapt,{'_'})
             %  figname=char(figname);
                % Figure
                h_fig=figure;
                h_fig.Color='white'; % changes the background color of the figure
                h_fig.Units='normalized'; % Defines the units
                h_fig.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
                h_fig.PaperOrientation='landscape';
          %      h_fig.Name=figname;
                    %% Plot loop
                   
                    for jj=1:obj.NumForceMaps
                        % Define variables
                        PlotTitle=obj.FM{jj}.ID;                   
                            ax=nexttile;
                            ax.XLim = [XMin XMax];
                            ax.YLim = [YMin YMax];                             
                            %ax.YLim = [obj.MinPullingLength obj.MaxPullingLength];
                            hold on
                            grid on                           
                            boxplot(nonzeros(obj.FM{jj}.PullingLength));
                           
                            % Title for each Subplot
                            %ti=title(sprintf('%i',jj),'Color','k');
                            ti=title(num2str(PlotTitle));
                            ti.Units='normalized'; % Set units to 'normalized'
                            ti.Position=[0.5,1]; % Position the subplot title within the subplot                     
                    end                  
              
%                 %% Save figures
%                 %%% Define the name for the figure title
%                 partname=sprintf('-p%d',ii);
%                 % fullname=sprintf('%s%s',figname,partname);
%                 fullname=sprintf('%s%s',figname,partname);
%                 %%% Save the current figure in the current folder
%                 print(gcf,fullname,'-dpng');
           
        %    close all
        
        
           
        end
             
       
       

        function SMFS_fine_figure(obj,XMin,XMax,YMin,YMax,ii)
            % Function to plot individual fine figures for publication
            if nargin < 2
                XMin= -inf;
                XMax= inf;
                YMin= -inf;
                YMax= inf;
            end
            ii=65
            % Change into the Folder of Interest
            cd(obj.ExperimentFolder) % Move into the folder
            % Create folders for saving the produced figures
            foldername='FineFigures';    % Defines the folder name
            mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
            currpath=fullfile(obj.ExperimentFolder,foldername);
            cd(currpath);
            % Load FM function           
            obj.FM{ii}.fc_fine_figure(XMin,XMax,YMin,YMax,ii)
        end
                
        function SMFS_statistics(obj)
           
            % Uncorrupt force curves
            SumFcUncorrupt=sum(obj.NumFcUncorrupt);
            SumFcCorrupt=obj.NumForceMaps*obj.FM{1}.NCurves;
            PercentUncorrupt=SumFcUncorrupt/SumFcCorrupt*100;
            PercentCorrupt=100-PercentUncorrupt;
        end
        
        % Individual ForceMap function related
        
        function SMFS_print_pulllength(obj)
            % SMFS_print: A function to simply plot all force curves of all
            % force maps loaded and calssified based on the SMFS Flag
            % Needed function: obj.presorting

            
            % Figure visibility
            set(groot,'defaultFigureVisible','off')      
            %set(groot,'defaultFigureVisible','on')           
            % Change into the Folder of Interest
            cd(obj.ExperimentFolder) % Move into the folder 
            % Create folders for saving the produced figures
            %foldername='FM_test';    % for debugging
            foldername='FM_Pulllength_MAD';    % Defines the folder name
            mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
            currpath=fullfile(obj.ExperimentFolder,foldername);
            cd(currpath); 
            
            % Loop over the imported force maps
            %for ii=1:obj.NumForceMaps
            for ii=2 % Debugging
               % Command window output
               sprintf('Force Map No. %d of %d',ii,obj.NumForceMaps) % Gives current Force Map Position
               % Run the chosen functions
               obj.FM{ii}.fc_pulling_length_MAD;     
            end    
           end
               
        function SMFS_snap_in_length_MAD(obj)            
            
            % Figure visibility
            %set(groot,'defaultFigureVisible','off')      
             set(groot,'defaultFigureVisible','on')           
            % Change into the Folder of Interest
            cd(obj.ExperimentFolder) % Move into the folder 
            % Create folders for saving the produced figures
            %foldername='FM_test';    % for debugging
            foldername='FM_SnapIn_MAD';    % Defines the folder name
            mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
            currpath=fullfile(obj.ExperimentFolder,foldername);
            cd(currpath); 
            
            % Loop over the imported force maps
            %for ii=1:obj.NumForceMaps
            for ii=2 % Debugging
               % Command window output
               sprintf('Force Map No. %d of %d',ii,obj.NumForceMaps) % Gives current Force Map Position
               % Run the chosen functions
               obj.FM{ii}.fc_snap_in_length_MAD;
            end 
        end
    
        function SMFS_adh_force_max(obj)            
            
            % Figure visibility
            set(groot,'defaultFigureVisible','off')      
            %set(groot,'defaultFigureVisible','on')           
            % Change into the Folder of Interest
            cd(obj.ExperimentFolder) % Move into the folder 
            % Create folders for saving the produced figures
            %foldername='FM_test';    % for debugging
            foldername='FM_Adh_Force';    % Defines the folder name
            mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
            currpath=fullfile(obj.ExperimentFolder,foldername);
            cd(currpath); 
            
            % Loop over the imported force maps
            %for ii=1:obj.NumForceMaps
            for ii=11:13 % Debugging
               % Command window output
               sprintf('Force Map No. %d of %d',ii,obj.NumForceMaps) % Gives current Force Map Position
               % Run the chosen functions
               obj.FM{ii}.fc_adh_force_max;
            end 
        end
   

        function SMFS_min_max(obj)
            

            for ii=1:obj.NumForceMaps
                %    obj.FM{ii}.base_and_tilt('linear');
                %obj.FM{ii}.fc_min_max_values;
                
                if ii==1
                    ConcatArrayMax=obj.FM{ii}.FMPullingLengthMax;
                 %   ConcatArrayAdhEnergy=obj.FM{ii}.obj.MinRet;
                else               
                    ConcatArrayMax=horzcat(ConcatArrayMax,obj.FM{ii}.FMPullingLengthMax);
                  %  ConcatArrayAdhEnergy=horzcat(ConcatArrayMax,obj.FM{ii}.FMPullingLengthMax);
                end                
            end
             ExpPullingLengthMax=max(ConcatArrayMax)
        end
        
        function [m,n,NumFigures] = adjust_tiled_layout(obj,NumFcMax)
            
            if nargin < 2
                NumFcMax=25; % The maximum of allowed plots per figure
            end
            
            for ii=1:obj.NumForceMaps
                %for ii=1:8 %for debugging
                NumFcUncorrupt(ii)=nnz(obj.FM{ii}.SMFSFlag.Uncorrupt); % Determine the number of uncorrupted force curves
                if ~any(NumFcUncorrupt(ii))    
                    continue
                end
                NumFigures=ceil(NumFcUncorrupt(ii)./NumFcMax); % Determine the number of figures
                Remainder=mod(NumFcUncorrupt(ii),NumFcMax); % Check for remainder
                if Remainder ~= 0
                    m(ii)=floor(sqrt(Remainder)); % Determine the number of rows in the figure
                    n(ii)=ceil(sqrt(Remainder)); % Determine the number of columns in the figure
                else
                    m(ii)=sqrt(NumFcMax);
                    n(ii)=m(ii);
                end
            end
        end
        
        
         function SMFS_initialize_flags(obj)
             
             obj.initialize_flags
             for ii=1:obj.NumForceMaps
                  obj.FM{ii}.initialize_flags
             end
         end
         
        
        function SMFS_testing_function(obj)
            % Function to quickly loop over all force maps for testing and
            % debugging
            
            % Figure visibility
            %set(groot,'defaultFigureVisible','off')      
             set(groot,'defaultFigureVisible','on')  
            
             
             for ii=1:obj.NumForceMaps
          %  for ii=135

          
               obj.FM{ii}.fc_TipHeight_calculation
          %    obj.FM{ii}.initialize_flags
          %
          
            
                 % obj.write_to_log_file('Test','variable ii','start')
                 % aa=33;
                  
                 % obj.write_to_log_file('','','end')
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
                
            NFM = obj.NumForceMaps;
                obj.SMFSFlagDown.SelectFM = false(NFM,1);
                obj.SMFSFlagDown.Preprocessed = false(NFM,1);
                obj.SMFSFlagDown.Presorted = false(NFM,1);
                obj.SMFSFlagDown.Fit= false(NFM,1);
                
            if isempty(obj.FMFlag)
                NFM = obj.NumForceMaps;
                obj.FMFlag.FibrilAnalysis = false(NFM,1);
                obj.FMFlag.ForceMapAnalysis = false(NFM,1);
                obj.FMFlag.Preprocessed = false(NFM,1);
                obj.FMFlag.Grouping = false;
                NSPM = obj.NumSurfacePotentialMaps;
                obj.SPMFlag.FibrilAnalysis = false(NSPM,1);
                obj.SPMFlag.Grouping = false;
                obj.SMFSFlag.SelectFM = false(NFM,1);
                obj.SMFSFlag.Preprocessed = false(NFM,1);
                obj.SMFSFlag.Presorted = false(NFM,1);
                obj.SMFSFlag.Uncorrupt= false(NFM,1);
                obj.SMFSFlag.AppMinCrit= false(NFM,1);
                obj.SMFSFlag.RetMinCrit= false(NFM,1);
                obj.SMFSFlag.Length= false(NFM,1);
                obj.SMFSFlag.Fit= false(NFM,1);
                obj.SMFSFlag.SnapIn= false(NFM,1);
                obj.SMFSFlag.PullingLength= false(NFM,1);
                
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
                obj.SMFSFlag.Preprocessed(end+1:NFM) = false(DiffFM,1);
                obj.SMFSFlag.Presorted(end+1:NFM) = false(DiffFM,1);
                
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

% funciton graveyard
%      %   function SMFS_analysis_classified(obj)
%             % Trial 10: Mica-HOH
%             % Cantilever 3
%             ii=235;
%             jj=236;
%                 C3=cat(2,obj.FM{ii}.RetAdhEnergy,obj.FM{jj}.RetAdhEnergy);
%                 C3Mean=mean(C3);
%                 C3Std=std(C3);
%                 C3Pull=cat(2,obj.FM{ii}.PullingLength ,obj.FM{jj}.PullingLength);
%                 C3PullMean=mean(C3Pull);
%                 C3PullStd=std(C3Pull);
%             % Cantilever 4
%             kk=241;
%             ll=242;
%                 C4=cat(2,obj.FM{kk}.RetAdhEnergy,obj.FM{ll}.RetAdhEnergy);
%                 C4Mean=mean(C4);
%                 C4Std=std(C4);
%                 C4Pull=cat(2,obj.FM{kk}.PullingLength ,obj.FM{ll}.PullingLength);
%                 C4PullMean=mean(C4Pull);
%                 C4PullStd=std(C4Pull);
% %             % Cantilever 5
% %             mm=245;
% %             nn=246;
% %                 C5=cat(2,obj.FM{mm}.RetAdhEnergy,obj.FM{nn}.RetAdhEnergy);
% %                 C5Mean=mean(C5);
% %                 C5Std=std(C5);
% %                 C5Pull=cat(2,obj.FM{mm}.PullingLength ,obj.FM{nn}.PullingLength);
% %                 C5PullMean=mean(C5Pull);
% %                 C5PullStd=std(C5Pull);
%             % Cantilever 7
%             oo=261;
%                 C7=obj.FM{oo}.RetAdhEnergy;
%                 C7Mean=mean(C7);
%                 C7Std=std(C7);
%                 C7Pull=obj.FM{oo}.PullingLength;
%                 C7PullMean=mean(C7Pull);
%                 C7PullStd=std(C7Pull);
%              % Concatenate arrays
%              CantRetAdhEnergyMean=cat(1,C3Mean,C4Mean,C7Mean);
%              CantRetAdhEnergyStd=cat(1,C3Std,C4Std,C7Std);
%              CantPullingLengthMean=cat(1,C3PullMean,C4PullMean,C7PullMean);
%              CantPullingLengthStd=cat(1,C3PullStd,C4PullStd,C7PullStd);
%              CantRetAdhCorrPullinglength=CantRetAdhEnergyMean./CantPullingLengthMean;
%              
%             % figure 1
%             h_fig=figure(1);
%             h_fig.Color='white'; % changes the background color of the figure
%             h_fig.Units='normalized'; % Defines the units
%             h_fig.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
%             h_fig.PaperOrientation='landscape';
%             [columns,~]=size(CantRetAdhEnergyMean);
%             errlow=CantRetAdhEnergyStd*-1;
%             errhigh=CantRetAdhEnergyStd;
%             xbar=1:1:columns;
%             % bar chart
%             bar(xbar,CantRetAdhEnergyMean)
%             hold on
%             er = errorbar(xbar,CantRetAdhEnergyMean,errhigh,errlow);    
%             er.Color = [0 0 0];                            
%             er.LineStyle = 'none';
%             hold off
%             %%% Axes 
%             ax = gca; % current axes
%             ax.FontSize = 20;
%             ax.XLabel.String = 'Cantilever';
%             ax.XLabel.FontSize = 20;
%             ax.YLabel.String = 'Adhesion Energy (N*m)';
%             ax.YLabel.FontSize = 20;
%             % figure 2
%             h_fig=figure(2);
%             h_fig.Color='white'; % changes the background color of the figure
%             h_fig.Units='normalized'; % Defines the units
%             h_fig.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
%             h_fig.PaperOrientation='landscape';
%             [columns,~]=size(CantPullingLengthMean);
%             errlow=CantPullingLengthStd*-1;
%             errhigh=CantPullingLengthStd;
%             xbar=1:1:columns;
%             % bar chart
%             bar(xbar,CantPullingLengthMean)
%             hold on
%             er = errorbar(xbar,CantPullingLengthMean,errhigh,errlow);    
%             er.Color = [0 0 0];                            
%             er.LineStyle = 'none';
%             hold off
%             %%% Axes 
%             ax = gca; % current axes
%             ax.FontSize = 20;
%             ax.XLabel.String = 'Cantilever';
%             ax.XLabel.FontSize = 20;
%             ax.YLabel.String = 'Pulling length (m)';
%             ax.YLabel.FontSize = 20;
%             % figure 3
%             h_fig=figure(3);
%             h_fig.Color='white'; % changes the background color of the figure
%             h_fig.Units='normalized'; % Defines the units
%             h_fig.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
%             h_fig.PaperOrientation='landscape';
%             [columns,~]=size(CantRetAdhCorrPullinglength);
%             xbar=1:1:columns;
%             % bar chart
%             bar(xbar,CantRetAdhCorrPullinglength)
%             %%% Axes 
%             ax = gca; % current axes
%             ax.FontSize = 20;
%             ax.XLabel.String = 'Cantilever';
%             ax.XLabel.FontSize = 20;
%             ax.YLabel.String = 'Ahesion Energy / Pulling length (N)';
%             ax.YLabel.FontSize = 20;
%         
%         end