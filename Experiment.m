classdef Experiment < matlab.mixin.Copyable & matlab.mixin.SetGet
    
    properties
        % Essential properties for File and subclass management
        ExperimentName = ''      % Shows name of the experiment
        ExperimentFolder = ''    % Shows Folder where Experiment is saved
        HostOS = ''              % Shows the current operating system
        HostName = ''            % Shows the current Host (User of running machine)
        BigDataFlag = true         % Determines how Experiment loads ForceMap objects. If true, 
                            % force voluime data are not loaded into RAM
                            % but read from unpacked jpk data container
        PythonLoaderFlag = true    % All ForceMap objects raw data are loaded by the 
                            % Python.zipfile module. This avoids the giant
                            % folder structures that get built in
                            % BigData-mode, while keeping load on memory
                            % AND disk relatively low and the programm
                            % almost as quick as in original
                            % 'all-data-to-memory'-mode
        KeepPythonFilesOpen = false % Decides whether to preload all PythonLoader Files into memory
                            % all the time
        FractionedSaveFiles = true
        CurrentLogFile = ''
        DummyAFMIMageFile = "TGT-1_MSCT-F-Box33-Nr10-fit-2022.05.13-09.50.09.883.jpk-qi-image"
        ShowImageSettings = Experiment.set_default_show_image_settings()
        ForceMapAnalysisOptions = Experiment.set_default_fma_options()
        CurrentFMA_ID = 'none'
        FMAExitErrorStack
        FiberSegmentProcessingOptions = Experiment.set_default_fsp_options()
        GrammOptions = Experiment.set_default_gramm_options()
        CustomCantileverTipOptions = Experiment.set_custom_cantilever_tip_options()
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
        ReferenceSlopeFlag
        AssignedReferenceMaps
        CantileverTipFlag
        AssignedCantileverTips
    end
    
    methods
        % constructor method and methods related with Experiment-file handling
        
        function obj = Experiment()
            % function obj = Experiment()
            %
            % The Experiment function initializes an experiment object by
            % gathering user input, creating an experiment folder, and
            % loading various types of files (Force Maps, Reference Force
            % Maps, AFM Images, Surface Potential Maps, and Cantilever
            % Tips) into the object. It also initializes flags and
            % properties, loads pretrained neural networks, and saves the
            % experiment object.
            %
            % Input: - There are no input parameters for this function, as
            % the necessary information is gathered through user input
            % dialogs.
            %
            % Output: - obj: The initialized experiment object with loaded
            % data, neural networks, and updated properties.
            
            % Set HostOS and HostName properties
            obj.check_for_new_host
            
            % set Name and choose Experiment folder
            [ExperimentName, ParentFolder] = uiputfile('*.mat','Select Experiment Name and Parent Folder');
            obj.ExperimentName = ExperimentName(1:end-4);
            mkdir(ParentFolder,obj.ExperimentName);
            obj.ExperimentFolder = fullfile(ParentFolder,obj.ExperimentName,filesep);
            
            % get Experiment filenames and paths from user input
            [FileTypes, FullFileStruct, IsValid, BigData,PythonLoaderFlag,KeepPythonFilesOpen] = obj.constructor_user_input_parser(obj.ExperimentName,obj.HostOS);
            
            if ~IsValid
                obj = [];
                return
            end
            
            % get paths of requested files and load them in
            obj.take_paths_and_load_files(FileTypes,FullFileStruct,true,BigData,PythonLoaderFlag,KeepPythonFilesOpen)
            
            obj.initialize_flags
            
            Temp = load('DropoutNetFinal.mat');
            obj.DropoutNet = Temp.MC14_Drop;
            Temp2 = load('CP_CNN_Final.mat');
            obj.CP_CNN = Temp2.CNN;
            
            obj.FractionedSaveFiles = true;
            obj.save_experiment();
        end
        
        function take_paths_and_load_files(obj,FileTypes,FullFileStruct,isNew,BigData,PythonLoaderFlag,KeepPythonFilesOpen)
            % function
            % take_paths_and_load_files(obj,FileTypes,FullFileStruct,isNew,BigData,PythonLoaderFlag,KeepPythonFilesOpen)
            %
            % This function processes and loads various types of files
            % (Force Maps, Reference Force Maps, AFM Images, Surface
            % Potential Maps, and Cantilever Tips) into an experiment
            % object. It sets up the experiment object properties according
            % to the input parameters and loads the data from the specified
            % file paths. It also supports parallel processing when loading
            % multiple files.
            %
            % Input: - obj: The experiment object to which the loaded data
            % will be added. - FileTypes (boolean array, size: 1x5):
            % Indicates which file types are to be loaded. The order is
            % [Force Maps, Reference Force Maps, AFM Images, Surface
            % Potential Maps, Cantilever Tips]. - FullFileStruct (struct):
            % A struct containing the full file paths for each of the five
            % file types. Each field in the struct should be named FullFile
            % and contain a cell array of file paths. - isNew (boolean,
            % optional): A flag indicating if this is a new experiment. If
            % true, existing data in the object will be reset. Defaults to
            % false. - BigData (boolean, optional): A flag indicating if
            % the data being loaded is large, affecting the memory
            % management. Defaults to false. - PythonLoaderFlag (boolean,
            % optional): A flag indicating if Python should be used for
            % data loading. Defaults to false. - KeepPythonFilesOpen
            % (boolean, optional): A flag indicating if the Python files
            % should be kept open after loading. Defaults to false.
            %
            % Output: - The function updates the obj properties with the
            % loaded data.
            
            obj.BigDataFlag = BigData;
            obj.PythonLoaderFlag = PythonLoaderFlag;
            obj.KeepPythonFilesOpen = KeepPythonFilesOpen;
            
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
            
            if contains(struct2array(ver), 'Parallel Computing Toolbox') && ((sum(NumFiles(1:2)) > 1) || (sum(NumFiles) > 20)) && ~obj.PythonLoaderFlag
                for i=1:5
                    parfor j=1:L
                        if (j+sum(NumFiles(1:i))-NumFiles(i)) <= length(IDs)
                            TempID{j,i} = sprintf('%s-%i',obj.ExperimentName,IDs(j+sum(NumFiles(1:i))-NumFiles(i)));
                        end
                        if i == 1 && FileTypes(i) && (j<=NumFiles(i))
                            TempCell{j,i} = ForceMap(FMFullFile{j},obj.ExperimentFolder,TempID{j,i},obj.BigDataFlag,obj.PythonLoaderFlag,obj.KeepPythonFilesOpen);
                        end
                        if i == 2 && FileTypes(i) && (j<=NumFiles(i))
                            TempCell{j,i} = ForceMap(RefFMFullFile{j},obj.ExperimentFolder,TempID{j,i},obj.BigDataFlag,obj.PythonLoaderFlag,obj.KeepPythonFilesOpen);
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
                            TempCell{j,i} = ForceMap(FMFullFile{j},obj.ExperimentFolder,TempID{j,i},obj.BigDataFlag,obj.PythonLoaderFlag,obj.KeepPythonFilesOpen);
                        end
                        if i == 2 && FileTypes(i)
                            TempCell{j,i} = ForceMap(RefFMFullFile{j},obj.ExperimentFolder,TempID{j,i},obj.BigDataFlag,obj.PythonLoaderFlag,obj.KeepPythonFilesOpen);
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
            % function Out = add_data(obj)
            %
            % The add_data function allows users to add new data to an
            % existing experiment object. The function prompts the user for
            % file selection, loads the new files, updates the experiment
            % object, and saves the updated object. If an error occurs
            % during this process, the function restores the original
            % experiment object.
            %
            % Input: - obj: The existing experiment object that new data
            % should be added to.
            %
            % Output: - Out: The updated experiment object after adding the
            % new data.
            
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
                SaveCopy.take_paths_and_load_files(FileTypes,FullFileStruct,isNew,SaveCopy.BigDataFlag,SaveCopy.PythonLoaderFlag,obj.KeepPythonFilesOpen)
                
                SaveCopy.initialize_flags % What to do with this?
                
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
        
        function create_custom_cantilever_tip(obj,varargin)
            % function create_custom_cantilever_tip(obj,varargin)
            %
            %
            % The create_custom_cantilever_tip function creates a custom
            % cantilever tip for an experiment object, optionally allowing
            % users to modify the tip shape options through a user
            % interface. The function then updates the properties of the
            % experiment object and saves the object if desired.
            %
            % Input: - obj: The existing experiment object to which a
            % custom cantilever tip should be added.
            %
            % - varargin (Optional Name-Value pairs):
            %
            % - "OpenUIOptions": A logical value indicating whether to open
            % the user interface for modifying tip shape options. Default
            % is true.
            %
            % - "SaveExperiment": A logical value indicating whether to
            % save the experiment object after adding the custom cantilever
            % tip. Default is true.
            %
            %
            % The function performs the following steps: 1. Parses the
            % input arguments, ensuring the validity of the experiment
            % object and optional name-value pairs.
            %
            % 2. If "OpenUIOptions" is true, opens the user interface for
            % modifying the tip shape options and updates the
            % CustomCantileverTipOptions property of the experiment object.
            %
            % 3. Updates the experiment object properties related to the
            % custom cantilever tip, such as CantileverTipFolders,
            % CantileverTips, CantileverTipFlag, CantileverTipNames, and
            % NumCantileverTips.
            %
            % 4. If "SaveExperiment" is true, saves the updated experiment
            % object.
            
            p = inputParser;
            p.FunctionName = "create_custom_cantilever_tip";
            p.CaseSensitive = false;
            p.PartialMatching = true;
            
            % Required inputs
            validobj = @(x)true;
            addRequired(p,"obj",validobj);
            
            % NameValue inputs
            defaultOpenUIOptions = true;
            defaultSaveExperiment = true;
            validOpenUIOptions = @(x)islogical(x);
            validSaveExperiment = @(x)islogical(x);
            addParameter(p,"OpenUIOptions",defaultOpenUIOptions,validOpenUIOptions);
            addParameter(p,"SaveExperiment",defaultSaveExperiment,validSaveExperiment);
            
            parse(p,obj,varargin{:});
            
            % Assign parsing results to named variables
            obj = p.Results.obj;
            OpenUIOptions = p.Results.OpenUIOptions;
            SaveExperiment = p.Results.SaveExperiment;
            
            if OpenUIOptions
                obj.CustomCantileverTipOptions = ui_set_struct_fields(obj.CustomCantileverTipOptions);
            end
            
            
            
            % Adjust Experiment properties and save if wanted
            obj.CantileverTipFolders{end+1} = TempPath;
            obj.CantileverTips{end+1} = tgt1_dummy;
            obj.CantileverTipFlag = 1;
            obj.CantileverTipNames{1} = TempFile;
            obj.NumCantileverTips = obj.NumCantileverTips + 1;
            
            if SaveExperiment
                obj.save_experiment
            end
            
        end
        
        function save_experiment(obj)
            % function save_experiment(obj)
            %
            % The save_experiment function saves the current state of an
            % experiment object to a .mat file in the specified
            % ExperimentFolder. If the FractionedSaveFiles property is set,
            % the function saves the individual AFM classes in separate
            % files within the ExperimentFolder before saving the
            % experiment object.
            %
            % Input: - obj: The experiment object to be saved.
            %
            % The function performs the following steps:
            %
            % 1. Changes the current directory to the ExperimentFolder.
            %
            % 2. Initializes temporary variables to store the OpenZipFile
            % property of force maps and reference force maps.
            %
            % 3. If FractionedSaveFiles is empty or false, saves the
            % experiment object as a single .mat file.
            %
            % 4. If FractionedSaveFiles is true, saves each AFM class (AFM
            % images, cantilever tips, force maps, and reference force
            % maps) in separate files within the ExperimentFolder, and then
            % saves the experiment object as a .mat file.
            %
            % 5. Restores the original current directory and OpenZipFile
            % properties of force maps and reference force maps.
            %
            % 6. Displays a message indicating that the experiment object
            % has been saved.
            
            
            current = what();
            cd(obj.ExperimentFolder)
            savename = sprintf('%s.mat',obj.ExperimentName);
            disp('saving...');
            TempZipFiles = cell(obj.NumForceMaps,1);
            TempRefZipFiles = cell(obj.NumReferenceForceMaps,1);
            for i=1:obj.NumForceMaps
                TempZipFiles{i} = obj.FM{i}.OpenZipFile;
                obj.FM{i}.OpenZipFile = [];
            end
            for i=1:obj.NumReferenceForceMaps
                TempRefZipFiles{i} = obj.RefFM{i}.OpenZipFile;
                obj.RefFM{i}.OpenZipFile = [];
            end
            if isempty(obj.FractionedSaveFiles) || ~obj.FractionedSaveFiles
                save(savename,'obj','-v7.3')
            elseif obj.FractionedSaveFiles
                for i=1:obj.NumAFMImages
                    obj.I{i}.save_afm_class(obj.ExperimentFolder);
                    obj.AFMImageFolders{i} = obj.I{i}.Folder;
                    obj.I{i}.clear_all_properties;
                end
                for i=1:obj.NumCantileverTips
                    obj.CantileverTips{i}.save_afm_class(obj.ExperimentFolder);
                    obj.CantileverTipFolders{i} = obj.CantileverTips{i}.Folder;
                    obj.CantileverTips{i}.clear_all_properties;
                end
                for i=1:obj.NumForceMaps
                    obj.FM{i}.save_afm_class(obj.ExperimentFolder);
                    obj.ForceMapFolders{i} = obj.FM{i}.Folder;
                    obj.FM{i}.clear_all_properties;
                end
                for i=1:obj.NumReferenceForceMaps
                    obj.RefFM{i}.save_afm_class(obj.ExperimentFolder);
                    obj.ReferenceForceMapFolders{i} = obj.RefFM{i}.Folder;
                    obj.RefFM{i}.clear_all_properties;
                end
                save(savename,'obj','-v7')
                obj.load_fractioned_afm_classes
            end
            cd(current.path)
            savemsg = sprintf('Changes to Experiment %s saved to %s',obj.ExperimentName,obj.ExperimentFolder);
            for i=1:obj.NumForceMaps
                obj.FM{i}.OpenZipFile = TempZipFiles{i};
            end
            for i=1:obj.NumReferenceForceMaps
                obj.RefFM{i}.OpenZipFile = TempRefZipFiles{i};
            end
            disp(savemsg);
        end
        
        function load_fractioned_afm_classes(obj)
            % function load_fractioned_afm_classes(obj)
            %
            % The load_fractioned_afm_classes function loads AFM class
            % properties for AFM images, cantilever tips, force maps, and
            % reference force maps from separate files within the
            % ExperimentFolder, assuming that the FractionedSaveFiles
            % property is true.
            %
            % Input: - obj: The experiment object containing the AFM
            % classes to be loaded.
            %
            % The function performs the following steps:
            %
            % 1. Iterates through the AFMImageFolders array and calls the
            % load_afm_class_properties method for each AFM image in the
            % array, passing the folder path as an argument.
            %
            % 2. Iterates through the CantileverTipFolders array and calls
            % the load_afm_class_properties method for each cantilever tip
            % in the array, passing the folder path as an argument.
            %
            % 3. Iterates through the ForceMapFolders array and calls the
            % load_afm_class_properties method for each force map in the
            % array, passing the folder path as an argument.
            %
            % 4. Iterates through the ReferenceForceMapFolders array and
            % calls the load_afm_class_properties method for each reference
            % force map in the array, passing the folder path as an
            % argument.
            %
            %
            % No output is returned. The function updates the AFM class
            % properties of the AFM images, cantilever tips, force maps,
            % and reference force maps within the experiment object.
            
            for i=1:obj.NumAFMImages
                obj.I{i}.load_afm_class_properties(obj.AFMImageFolders{i});
            end
            for i=1:obj.NumCantileverTips
                obj.CantileverTips{i}.load_afm_class_properties(obj.CantileverTipFolders{i});
            end
            for i=1:obj.NumForceMaps
                obj.FM{i}.load_afm_class_properties(obj.ForceMapFolders{i});
            end
            for i=1:obj.NumReferenceForceMaps
                obj.RefFM{i}.load_afm_class_properties(obj.ReferenceForceMapFolders{i});
            end
            
        end
        
        function ExperimentCopy = copy_experiment(obj)
            % ExperimentCopy = copy_experiment(obj)
            %
            %
            % The copy_experiment function creates a deep copy of the input
            % experiment object, including all its associated AFM images,
            % cantilever tips, force maps, reference force maps, and
            % surface potential maps.
            %
            % Input: - obj: The experiment object to be copied.
            %
            % Output: - ExperimentCopy: The deep copy of the input
            % experiment object.
            %
            % The function performs the following steps:
            %
            % 1. Calls the copy method of the input experiment object to
            % create a shallow copy named ExperimentCopy.
            %
            % 2. Iterates through the highest number of AFM images,
            % cantilever tips, force maps, reference force maps, and
            % surface potential maps within the input experiment object.
            %
            % 3. For each object (AFM image, cantilever tip, force map,
            % reference force map, and surface potential map), checks
            % whether it belongs to the 'matlab.mixin.Copyable' superclass.
            %
            % 4. If the object belongs to the 'matlab.mixin.Copyable'
            % superclass, the copy method of that object is called to
            % create a deep copy and it is added to the corresponding array
            % in the ExperimentCopy object.
            %
            %
            % Note: The input experiment object must belong to the
            % 'matlab.mixin.Copyable' superclass in order to be copied by
            % this function.

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
            %function delete_experiment(obj)
            %
            % The delete_experiment function permanently deletes the input
            % experiment object, its folder, experiment file, and all
            % subfolders. This method is only available on Windows systems.
            %
            % Input: - obj: The experiment object to be deleted.
            %
            % The function performs the following steps:
            %
            % 1. Checks whether the host operating system is Windows
            % ('PCW'). If not, displays a message that the method is only
            % available on Windows systems and returns.
            %
            % 2. Sets the Folder variable as the experiment folder.
            %
            % 3. Displays a warning dialog to confirm deletion, providing
            % the user with two choices: 'Yes, delete all' or 'Abort'.
            %
            % 4. If the user selects 'Yes, delete all', a second warning
            % dialog is displayed, asking the user to confirm once more:
            % 'YES, DO IT!' or 'Oh no, I changed my mind'.
            %
            % 5. If the user confirms deletion ('YES, DO IT!'), the
            % function changes the current directory to the parent of the
            % experiment folder, constructs the necessary command-line
            % commands, and calls the system function to execute the
            % commands, deleting the experiment folder and its contents.
            %
            % 6. If the user selects 'Abort' or 'Oh no, I changed my mind'
            % at any point, the function returns without deleting the
            % experiment.
            %
            %
            % Note: This function is designed to work only on Windows
            % systems. Use with caution, as the deletion is permanent and
            % cannot be undone.
            
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
                    answer = questdlg('Are you REALLY sure?', ...
                        sprintf('Deletion of %s',Folder),'YES, DO IT!', ...
                        'Oh no, i changed my mind','Oh no, i changed my mind');
                    switch answer
                        case 'YES, DO IT!'
                            
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
                        case 'Oh no, i changed my mind'
                            return
                    end
                case 'Abort'
                    return
            end
            
        end
        
        function update_absolute_paths(E,Path,JustExperiment)
            % function update_absolute_paths(E,Path,JustExperiment)
            %
            % The update_absolute_paths function updates the absolute file
            % paths of the experiment object and its associated objects,
            % such as force maps, AFM images, reference force maps, and
            % cantilever tips.
            %
            % Inputs: - E: The experiment object whose absolute paths need
            % to be updated.
            %
            % - Path: The new top-level path (string) to replace the old
            % top-level path in the experiment object and its associated
            % objects.
            %
            % - JustExperiment (optional, default = 0): A flag (0 or 1) to
            % indicate whether to update only the experiment object (1) or
            % also its associated objects (0).
            %
            % The function performs the following steps: 1. Checks the
            % number of input arguments. If JustExperiment is not provided,
            % it is set to 0 by default.
            %
            % 2. Replaces the file separators in the experiment folder with
            % the appropriate separators for the current system.
            %
            % 3. Iterates through the force maps, AFM images, reference
            % force maps, and cantilever tips associated with the
            % experiment object.
            %
            % 4. If JustExperiment is 0, the function updates the file
            % paths of the associated objects by replacing the old
            % top-level path with the new top-level path.
            %
            % 5. Updates the file paths of the experiment object itself.
            %
            % Note: This function should be used when the file paths of the
            % experiment and its associated objects need to be updated,
            % such as when moving the experiment to a new location.
            
            if nargin < 3
                JustExperiment = 0;
            end
            %             if isequal(E.HostOS,'GLN')
            %                 Path = [Path filesep];
            %             end
            
            E.ExperimentFolder = Experiment.replace_fileseps(E.ExperimentFolder);
            
            for i=1:E.NumForceMaps
                if ~isempty(E.FM{i})
                    if ~JustExperiment
                        E.FM{i}.check_for_new_host();
                        if E.BigDataFlag && ~E.PythonLoaderFlag
                            OldDataStore = ...
                                Experiment.replace_fileseps(E.FM{i}.DataStoreFolder);
                            Split = strsplit(OldDataStore,filesep);
                            E.FM{i}.DataStoreFolder = fullfile(Path,Split{end-1});
                        elseif E.BigDataFlag && E.PythonLoaderFlag
                            OldDataStore = ...
                                Experiment.replace_fileseps(E.FM{i}.DataStoreFolder);
                            Split = strsplit(OldDataStore,filesep);
                            E.FM{i}.DataStoreFolder = fullfile(Path,Split{end});
                            OldFilePath = split(...
                                Experiment.replace_fileseps(E.FM{i}.RawDataFilePath),...
                                filesep);
                            if iscell(E.FM{i}.RawDataFilePath)
                                E.FM{i}.RawDataFilePath = E.FM{i}.RawDataFilePath{1};
                            end
                            E.FM{i}.RawDataFilePath = fullfile(E.FM{i}.DataStoreFolder,OldFilePath{end});
                        end
                        E.FM{i}.Folder = E.switch_old_with_new_toplvl_path(...
                            Experiment.replace_fileseps(E.FM{i}.Folder),E.ExperimentFolder,Path);
                    end
                    E.ForceMapFolders{i} = E.switch_old_with_new_toplvl_path(...
                        Experiment.replace_fileseps(E.ForceMapFolders{i}),E.ExperimentFolder,Path);
                end
            end
            for i=1:E.NumAFMImages
                if ~isempty(E.I{i})
                    if ~JustExperiment
                        E.I{i}.check_for_new_host();
                        E.I{i}.Folder = E.switch_old_with_new_toplvl_path(...
                            Experiment.replace_fileseps(E.I{i}.Folder),E.ExperimentFolder,Path);
                    end
                    E.AFMImageFolders{i} = E.switch_old_with_new_toplvl_path(...
                        Experiment.replace_fileseps(E.AFMImageFolders{i}),E.ExperimentFolder,Path);
                end
            end
            for i=1:E.NumReferenceForceMaps
                if ~isempty(E.RefFM{i})
                    if ~JustExperiment
                        E.RefFM{i}.check_for_new_host();
                        if E.BigDataFlag && ~E.PythonLoaderFlag
                            E.RefFM{i}.DataStoreFolder = fullfile(Path,Split{end-1});
                            OldDataStore = Experiment.replace_fileseps(E.RefFM{i}.DataStoreFolder);
                            Split = strsplit(OldDataStore,filesep);
                        elseif E.BigDataFlag && E.PythonLoaderFlag
                            E.RefFM{i}.DataStoreFolder = fullfile(Path,Split{end});
                            OldFilePath = split(...
                                Experiment.replace_fileseps(E.RefFM{i}.RawDataFilePath),filesep);
                            E.RefFM{i}.RawDataFilePath = fullfile(E.RefFM{i}.DataStoreFolder,OldFilePath{end});
                        end
                        E.RefFM{i}.Folder = E.switch_old_with_new_toplvl_path(...
                            Experiment.replace_fileseps(E.RefFM{i}.Folder),E.ExperimentFolder,Path);
                    end
                    E.ReferenceForceMapFolders{i} = E.switch_old_with_new_toplvl_path(...
                        Experiment.replace_fileseps(E.ReferenceForceMapFolders{i}),E.ExperimentFolder,Path);
                end
            end
            for i=1:E.NumCantileverTips
                if ~isempty(E.CantileverTips{i})
                    if ~JustExperiment
                        E.CantileverTips{i}.check_for_new_host();
                        E.CantileverTips{i}.Folder = E.switch_old_with_new_toplvl_path(...
                            Experiment.replace_fileseps(E.CantileverTips{i}.Folder),E.ExperimentFolder,Path);
                    end
                    E.CantileverTipFolders{i} = E.switch_old_with_new_toplvl_path(...
                        Experiment.replace_fileseps(E.CantileverTipFolders{i}),E.ExperimentFolder,Path);
                end
            end
            
        end
        
        function add_dummy_afmimage_data(obj,Preprocessing)
            % function add_dummy_afmimage_data(obj,Preprocessing)
            %
            % The add_dummy_afmimage_data function adds an artificial AFM
            % image object to the experiment object based on a dummy AFM
            % image file.
            %
            % Inputs: - obj: The experiment object to which the artificial
            % AFM image object should be added.
            %
            % - Preprocessing (optional, default = false): A boolean flag
            % indicating whether preprocessing should be applied to the
            % artificial AFM image object.
            %
            %
            % The function performs the following steps:
            %
            % 1. Checks the number of input arguments. If Preprocessing is
            % not provided, it is set to false by default.
            %
            % 2. Retrieves the source file for the dummy AFM image.
            %
            % 3. Creates an AFMImage object from the source file and adds
            % it to the experiment object.
            %
            % 4. Updates the NumAFMImages property of the experiment object
            % to reflect the addition of the new AFMImage object.
            %
            % 5. Removes all file-specific information from the AFMImage
            % object, leaving only the structure and properties of an
            % artificial AFM image.
            %
            %
            % Note: This function is useful when you want to create a
            % placeholder AFM image for testing purposes or to simulate
            % data in an experiment object.
            
            if nargin <2
                Preprocessing = false;
            end
            
            SourceFile = which(obj.DummyAFMIMageFile);
            
            ImageClass = AFMImage(SourceFile,obj.ExperimentFolder,'ArtificialAFMImage-1',Preprocessing);
            
            if isempty(obj.NumAFMImages) || obj.NumAFMImages == 0
                obj.I{1} = ImageClass;
                obj.NumAFMImages = 1;
            else
                obj.I{end+1} = ImageClass;
                obj.NumAFMImages = obj.NumAFMImages + 1;
            end
            
            Index = obj.NumAFMImages;
            
            % Clean up the Instance of all file specific Information
            obj.I{Index}.Channel = [];
            obj.I{Index}.NumChannels = 0;
            obj.I{Index}.NumPixelsX = [];
            obj.I{Index}.NumPixelsY = [];
            obj.I{Index}.ScanAngle = [];
            obj.I{Index}.ScanSizeX = [];
            obj.I{Index}.ScanSizeY = [];
            obj.I{Index}.OriginX = [];
            obj.I{Index}.OriginY = [];
            obj.I{Index}.List2Map = [];
            obj.I{Index}.Map2List = [];
            obj.I{Index}.Name = '';
            obj.AFMImageNames{Index} = '';
        end
        
        function PC2HM_convert_pc_to_hm_and_add_afmimage(obj,PCName, PC, Resolution, zResolution,...
                MaxPointsPerGP, PartitionShrinkingFactor, InterpolationExpansionFactor,...
                Lambda, Sigma, Noise)
            %function PC2HM_convert_pc_to_hm_and_add_afmimage(obj,PCName,
            %PC, Resolution, zResolution,...
            %                 MaxPointsPerGP, PartitionShrinkingFactor,
            %                 InterpolationExpansionFactor,... Lambda,
            %                 Sigma, Noise)
            %
            % The PC2HM_convert_pc_to_hm_and_add_afmimage function converts
            % a point cloud to a height map and adds the height map as an
            % AFM image object to the experiment object.
            %
            % Inputs:
            %
            % - obj: The experiment object to which the height map should
            % be added as an AFM image.
            %
            % - PCName: A string representing the name of the point cloud.
            %
            % - PC: The point cloud matrix (Nx3) with columns representing
            % X, Y, and Z coordinates.
            %
            % - Resolution: The resolution of the height map (integer).
            %
            % - zResolution: The resolution in the Z direction (float).
            %
            % - MaxPointsPerGP: The maximum number of points per Gaussian
            % process (integer).
            %
            % - PartitionShrinkingFactor: The shrinking factor for
            % partitioning (float).
            %
            % - InterpolationExpansionFactor: The expansion factor for
            % interpolation (float).
            %
            % - Lambda: The regularization parameter for Gaussian processes
            % (float).
            %
            % - Sigma: The standard deviation for Gaussian processes
            % (float).
            %
            % - Noise: The noise level for Gaussian processes (float).
            %
            % The function performs the following steps:
            %
            % 1. Converts the point cloud to a height map using the
            % PC2HM_point_cloud_to_height_map function.
            %
            % 2. Adds a dummy AFM image object to the experiment object.
            %
            % 3. Sets the properties of the AFM image object based on the
            % input parameters and the output of the height map conversion.
            %
            % 4. Creates and adds channels for MaxPeakMap, MinPeakMap,
            % MaxPeakValueMap, and DensityMap to the AFM image object.
            %
            % 5. Constructs list-to-map relations for the AFM image object.
            %
            % 6. Sets the name of the AFM image object based on the input
            % parameters.
            %
            %
            % This function is useful when you want to convert a point
            % cloud representation of a surface into a height map and
            % incorporate it into an experiment as an AFM image object.

            [MaxPeakMap,MinPeakMap,MaxPeakValueMap,DensityMap,minCoords, maxCoords,...
                subsetsGrid, boundariesGrid, subsetsCloud, boundariesCloud] =...
                PC2HM_point_cloud_to_height_map(PC, Resolution, zResolution,...
                MaxPointsPerGP, PartitionShrinkingFactor, InterpolationExpansionFactor,...
                Lambda, Sigma, Noise);
            
            obj.add_dummy_afmimage_data(false)
            
            NewName = PCName;
            Counter = 1;
            while contains([obj.AFMImageNames{:}],NewName)
                NewName = sprintf('%s-%03i',PCName,Counter);
                Counter = Counter + 1;
            end
            
            CompoundName = NewName;
            
            obj.I{end}.Metadata.PC2HM.Resolution = Resolution;
            obj.I{end}.Metadata.PC2HM.zResolution = zResolution;
            obj.I{end}.Metadata.PC2HM.MaxPointsPerGP = MaxPointsPerGP;
            obj.I{end}.Metadata.PC2HM.PartitionShrinkingFactor = PartitionShrinkingFactor;
            obj.I{end}.Metadata.PC2HM.InterpolationExpansionFactor = InterpolationExpansionFactor;
            obj.I{end}.Metadata.PC2HM.Lambda = Lambda;
            obj.I{end}.Metadata.PC2HM.Sigma = Sigma;
            obj.I{end}.Metadata.PC2HM.Noise = Noise;
            
            obj.I{end}.NumPixelsX = Resolution;
            obj.I{end}.NumPixelsY = Resolution;
            obj.I{end}.ScanAngle = 0;
            obj.I{end}.ScanSizeX = maxCoords(1) - minCoords(1);
            obj.I{end}.ScanSizeY = maxCoords(2) - minCoords(2);
            obj.I{end}.OriginX = 0;
            obj.I{end}.OriginY = 0;
            MaxPeakChannel = obj.I{end}.create_standard_channel(MaxPeakMap,'MaxPeakMap','m');
            MinPeakChannel = obj.I{end}.create_standard_channel(MinPeakMap,'MinPeakMap','m');
            MaxPeakValueChannel = obj.I{end}.create_standard_channel(MaxPeakValueMap,'MaxPeakValueMap','m');
            DensityChannel = obj.I{end}.create_standard_channel(DensityMap,'DensityMap','m');
            obj.I{end}.Channel = MaxPeakChannel;
            obj.I{end}.Channel(2) = MinPeakChannel;
            obj.I{end}.Channel(3) = MaxPeakValueChannel;
            obj.I{end}.Channel(4) = DensityChannel;
            obj.I{end}.NumChannels = numel(obj.I{end}.Channel);
            obj.I{end}.construct_list_to_map_relations;
            obj.I{end}.Name = CompoundName;
            obj.AFMImageNames{end} = CompoundName;
            
        end
        
        function add_dummy_cantilever_tip_data(obj,Preprocessing,Name)
            % function add_dummy_cantilever_tip_data(obj,Preprocessing)
            %
            % The add_dummy_cantilever_tip_data function adds a dummy
            % cantilever tip object to the experiment object.
            %
            %
            % Inputs:
            %
            % - obj: The experiment object to which the dummy cantilever
            % tip should be added.
            %
            % - Preprocessing (optional): A boolean value indicating
            % whether preprocessing should be applied to the cantilever tip
            % data. Default is true.
            %
            % The function performs the following steps:
            %
            % 1. Locates the source file for the dummy cantilever tip data.
            %
            % 2. Creates an AFMImage object with the source file and the
            % experiment folder, using the specified preprocessing setting.
            %
            % 3. Adds the AFMImage object as a cantilever tip to the
            % experiment object's CantileverTips property, updating the
            % NumCantileverTips property accordingly.
            %
            % This function is useful when you want to add a dummy
            % cantilever tip object to an experiment for testing purposes,
            % without needing to load real data. The dummy cantilever tip
            % data can be used for simulating different aspects of the
            % experiment, such as analyzing the effect of the cantilever
            % tip on the acquired AFM images.
            
            if nargin <2
                Preprocessing = true;
                Name = 'ArtificialAFMImage-1';
            else
                Name = [Name '-1'];
            end
            
            SourceFile = which("TGT-1_MSCT-F-Box33-Nr10-fit-2022.05.13-09.50.09.883.jpk-qi-image");
            
            Tip = AFMImage(SourceFile,obj.ExperimentFolder,Name,Preprocessing);
            
            if isempty(obj.NumCantileverTips) || obj.NumCantileverTips == 0
                obj.CantileverTips{1} = Tip;
                obj.NumCantileverTips = 1;
            else
                obj.CantileverTips{end+1} = Tip;
                obj.NumCantileverTips = obj.NumCantileverTips + 1;
            end
        end
        
        function create_artificial_cantilever_tip(obj,Name,Shape,varargin)
            % function [OutChannel1,OutChannel2] = create_custom_height_topography(Shape,varargin)
            %
            % <FUNCTION DESCRIPTION HERE>
            %
            %
            % Required inputs
            % Shape ... <VARIABLE DESCRIPTION>
            %
            % Name-Value pairs
            % "TipRadius" ... <NAMEVALUE DESCRIPTION>
            % "TipHeight" ... <NAMEVALUE DESCRIPTION>
            % "TipTilt" ... <NAMEVALUE DESCRIPTION>
            % "Radius" ... <NAMEVALUE DESCRIPTION>
            % "HalfAngle" ... <NAMEVALUE DESCRIPTION>
            % "ImageResolution" ... <NAMEVALUE DESCRIPTION>
            % "StartFraction" ... <NAMEVALUE DESCRIPTION>
            % "EndFraction" ... <NAMEVALUE DESCRIPTION>
            % "FuseMethod" ... <NAMEVALUE DESCRIPTION>
            % "HeightDifference" ... <NAMEVALUE DESCRIPTION>
            % "TipApexChannel" ... <NAMEVALUE DESCRIPTION>
            
            p = inputParser;
            p.FunctionName = "create_custom_height_topography";
            p.CaseSensitive = false;
            p.PartialMatching = true;
            
            % Required inputs
            validShape = @(x)true;
            addRequired(p,"Shape",validShape);
            
            % NameValue inputs
            defaultTipRadius = [];
            defaultTipHeight = 1000e-9;
            defaultTipTilt = 0;
            defaultRadius = 20e-9;
            defaultHalfAngle = 15;
            defaultImageResolution = 256;
            defaultStartFraction = .8;
            defaultEndFraction = .99;
            defaultFuseMethod = 'linear';
            defaultHeightDifference = 0;
            defaultTipApexChannel = [];
            defaultProjectedAreaStepSize = 1e-9;
            validTipRadius = @(x)true;
            validTipHeight = @(x)true;
            validTipTilt = @(x)true;
            validRadius = @(x)true;
            validHalfAngle = @(x)true;
            validImageResolution = @(x)true;
            validStartFraction = @(x)true;
            validEndFraction = @(x)true;
            validFuseMethod = @(x)true;
            validHeightDifference = @(x)true;
            validTipApexChannel = @(x)true;
            validProjectedAreaStepSize = @(x)true;
            addParameter(p,"TipRadius",defaultTipRadius,validTipRadius);
            addParameter(p,"TipHeight",defaultTipHeight,validTipHeight);
            addParameter(p,"TipTilt",defaultTipTilt,validTipTilt);
            addParameter(p,"Radius",defaultRadius,validRadius);
            addParameter(p,"HalfAngle",defaultHalfAngle,validHalfAngle);
            addParameter(p,"ImageResolution",defaultImageResolution,validImageResolution);
            addParameter(p,"StartFraction",defaultStartFraction,validStartFraction);
            addParameter(p,"EndFraction",defaultEndFraction,validEndFraction);
            addParameter(p,"FuseMethod",defaultFuseMethod,validFuseMethod);
            addParameter(p,"HeightDifference",defaultHeightDifference,validHeightDifference);
            addParameter(p,"TipApexChannel",defaultTipApexChannel,validTipApexChannel);
            addParameter(p,"ProjectedAreaStepSize",defaultProjectedAreaStepSize,validProjectedAreaStepSize);
            
            parse(p,Shape,varargin{:});
            
            % Assign parsing results to named variables
            Shape = p.Results.Shape;
            TipRadius = p.Results.TipRadius;
            TipHeight = p.Results.TipHeight;
            TipTilt = p.Results.TipTilt;
            Radius = p.Results.Radius;
            HalfAngle = p.Results.HalfAngle;
            ImageResolution = p.Results.ImageResolution;
            StartFraction = p.Results.StartFraction;
            EndFraction = p.Results.EndFraction;
            FuseMethod = p.Results.FuseMethod;
            HeightDifference = p.Results.HeightDifference;
            TipApexChannel = p.Results.TipApexChannel;
            ProjectedAreaStepSize = p.Results.ProjectedAreaStepSize;
            
            obj.add_dummy_cantilever_tip_data(false,Name)
            
            CT = obj.CantileverTips{end};
            CT.Name = Name;
            obj.CantileverTipNames{end+1} = CT.Name;
            obj.CantileverTipFolders{end+1} = CT.Folder;
            
            Channel = AFMImage.create_custom_height_topography(Shape,...
                "TipRadius",TipRadius,...
                "TipHeight",TipHeight,...
                "TipTilt",TipTilt,...
                "Radius",Radius,...
                "HalfAngle",HalfAngle,...
                "ImageResolution",ImageResolution,...
                "StartFraction",StartFraction,...
                "EndFraction",EndFraction,...
                "FuseMethod",FuseMethod,...
                "HeightDifference",HeightDifference,...
                "TipApexChannel",TipApexChannel);
            
            Channel.Name = 'Processed';
            CT.Channel(end+1) = Channel;
            CT.hasProcessed = true;
            
            CT.Channel(1:CT.NumChannels) = [];
            
            Channel.Name = 'Eroded Tip';
            Channel.Image = Channel.Image - max(Channel.Image,[],'all');
            CT.Channel(end+1) = Channel;
            CT.hasDeconvolutedCantileverTip = true;
            
            CT.OriginX = Channel.OriginX;
            CT.OriginY = Channel.OriginY;
            CT.ScanAngle = Channel.ScanAngle;
            CT.NumPixelsX = Channel.NumPixelsX;
            CT.NumPixelsY = Channel.NumPixelsY;
            CT.ScanSizeX = Channel.ScanSizeX;
            CT.ScanSizeY = Channel.ScanSizeY;
            
            CT.NumChannels = 2;
            
            CT.create_pixel_difference_channel(Channel.Image);
            
            CT.ProjectedTipArea = AFMImage.calculate_projected_area_histcumsum(CT.get_channel('Eroded Tip'),ProjectedAreaStepSize);
            CT.DepthDependendTipRadius = CT.calculate_depth_dependend_tip_data(CT.ProjectedTipArea,75);
            
        end
        
        function load_force_map_analysis_options(obj)
            % function load_force_map_analysis_options(obj)
            %
            % The load_force_map_analysis_options function loads force map
            % analysis options from a .mat file and assigns them to the
            % experiment object.
            %
            % Inputs:
            %
            % - obj: The experiment object for which the force map analysis
            % options should be loaded and assigned.
            %
            % The function performs the following steps:
            
            % 1. Sets the default force map analysis options using
            % Experiment.set_default_fma_options.
            %
            % 2. Opens a file dialog to let the user choose a .mat file
            % containing the force map analysis options.
            %
            % 3. Loads the .mat file and checks if the file contains a
            % valid ForceMapAnalysis Options structure (FMA).
            %
            % 4. Checks if the loaded FMA structure has the same fields as
            % the default FMA structure.
            %
            % 5. If the loaded FMA structure is valid, assigns it to the
            % obj.ForceMapAnalysisOptions property.
            %
            % 6. Displays a message indicating the successful assignment of
            % the force map analysis options.
            %
            % This function is useful when you want to load a set of force
            % map analysis options from a .mat file and apply them to an
            % experiment object. The loaded options will then be used to
            % analyze force maps associated with the experiment. The
            % function checks the validity of the loaded options to ensure
            % compatibility with the current default options. If the loaded
            % options are not valid, an error is thrown.
            
            Default = Experiment.set_default_fma_options;
            
            
            [File,Path] = uigetfile('*.mat','Choose ForceMapAnalysis Options .mat from folder');
            Fullfile = fullfile(Path,File);
            
            disp('Loading ForceMapAnalysis Options...')
            load(Fullfile);
            
            if ~isvarname('FMA')
                error('The provided file is not a ForceMapAnalysis Options file')
            end
            
            if ~check_structs_for_same_fields_recursively(FMA,Default)
                error(sprintf('\n%s\n%s',...
                    'Error: The provided file is either not a ForceMapAnalysis',...
                    'Options file or the default struct has since been updated.'))
            end
            
            obj.ForceMapAnalysisOptions = FMA;
            disp('ForceMapAnalysis Options have been assigned successfully!')
            
        end
        
    end
    methods(Static)
        % Static methods related with Experiment-file handling
        
        function E = load(Fullfile)
            % E = load(Fullfile)
            %
            %
            % The load function reads an Experiment object from a .mat file
            % and loads it into memory. It also handles necessary updates
            % such as checking for a new host, updating absolute paths, and
            % loading/clearing zipped files using Python.
            %
            % Inputs:
            %
            % - Fullfile (optional): A string containing the full file path
            % of the .mat file to load. If not provided, a file dialog will
            % open for the user to choose the .mat file.
            %
            % Outputs:
            %
            % - E: The loaded Experiment object.
            %
            % The function performs the following steps:
            %
            % 1. Opens a file dialog to let the user choose a .mat file
            % containing the Experiment object if Fullfile is not provided.
            %
            % 2. Loads the .mat file and extracts the Experiment object
            % (obj).
            %
            % 3. Checks if the host has changed and updates the absolute
            % paths accordingly.
            %
            % 4. If FractionedSaveFiles is set, loads the fractioned AFM
            % classes.
            %
            % 5. Sets default values for PythonLoaderFlag,
            % KeepPythonFilesOpen, and BigDataFlag if they are empty.
            %
            % 6. Re-checks the host and updates absolute paths.
            %
            % 7. If PythonLoaderFlag is true, loads or clears zipped files
            % using Python, depending on the KeepPythonFilesOpen flag.
            %
            % 8. Sets the ExperimentFolder property to the loaded file's
            % path.
            %
            % This function is useful when you want to load an Experiment
            % object from a .mat file to work with it in the current
            % session. The loaded object is updated to ensure compatibility
            % with the current system and any necessary files are loaded or
            % cleared according to the set flags.
            
            if nargin < 1
                [File,Path] = uigetfile('*.mat','Choose Experiment .mat from folder');
                Fullfile = fullfile(Path,File);
            else
                [Path,FileName,Extension] = fileparts(Fullfile);
                Path = [Path filesep];
                File = [FileName Extension];
            end
            disp('Loading Experiment... this can take a while for larger Experiments')
            load(Fullfile);
            
            E = obj;
            clear obj
            
            E.check_for_new_host();
            E.update_absolute_paths(Path,true);
            
            if isempty(E.FractionedSaveFiles) || ~E.FractionedSaveFiles
            elseif E.FractionedSaveFiles
                E.load_fractioned_afm_classes
            end
            
            if isempty(E.PythonLoaderFlag)
                E.PythonLoaderFlag = false;
            end
            if isempty(E.KeepPythonFilesOpen)
                E.KeepPythonFilesOpen = false;
            end
            if isempty(E.BigDataFlag)
                E.BigDataFlag = false;
            end
            
            E.check_for_new_host();
            E.update_absolute_paths(Path);
            
            if E.PythonLoaderFlag
                if E.KeepPythonFilesOpen
                    for i=1:E.NumForceMaps
                        E.FM{i}.load_zipped_files_with_python
                    end
                    for i=1:E.NumReferenceForceMaps
                        E.RefFM{i}.load_zipped_files_with_python
                    end
                else
                    for i=1:E.NumForceMaps
                        E.FM{i}.clear_zipped_files_from_memory
                    end
                    for i=1:E.NumReferenceForceMaps
                        E.RefFM{i}.clear_zipped_files_from_memory
                    end
                end
            end
            
            E.ExperimentFolder = Path;
            
        end
        
        function delete_folderstructure(FolderPath)
            % function delete_folderstructure(FolderPath)
            %
            % The delete_folderstructure function permanently deletes the
            % specified folder, its files, and all subfolders on a Windows
            % system. This function should be used with caution as the
            % deletion is irreversible.
            %
            % Inputs:
            %
            % - FolderPath: A string containing the full path of the folder
            % to be deleted.
            %
            % Outputs: None.
            %
            % The function performs the following steps:
            %
            % 1. Checks if the current operating system is a Windows
            % system; if not, the function is not executed.
            %
            % 2. Displays a warning dialog asking the user to confirm the
            % folder deletion.
            %
            % 3. If the user confirms the deletion, the function:
            %
            % a. Changes the current working directory to the parent folder
            % of the target folder.
            %
            % b. Constructs and executes two command line commands (CMD1
            % and CMD2) to delete the folder and its contents.
            %
            % c. Restores the current working directory to its previous
            % state.
            %
            % Use this function when you need to permanently delete a
            % folder and its contents on a Windows system. Be aware that
            % this operation is irreversible, so use caution and make sure
            % you have backups of any important data before proceeding.
            
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
            %
            % The preprocessing function performs preprocessing steps on a
            % set of ForceMap objects. It creates and levels the height
            % maps, then fits base lines for each ForceMap. Optionally, the
            % user can skip ForceMaps that have already been preprocessed.
            %
            % Inputs:
            %
            % - obj: An instance of a class containing ForceMap objects
            % that need preprocessing.
            %
            % Outputs: None.
            %
            % The function performs the following steps:
            %
            % 1. Initializes a progress waitbar to indicate the
            % preprocessing progress.
            %
            % 2. Checks if any ForceMaps have already been preprocessed. If
            % so, the user is prompted to decide whether to skip them and
            % keep the old results or reprocess them.
            %
            % 3. Loops through each ForceMap in the obj, and for each:
            %
            % a. If the user chose to skip preprocessed ForceMaps and the
            % current ForceMap is flagged as preprocessed, the loop moves
            % to the next ForceMap.
            %
            % b. Creates and levels the height map for the current ForceMap
            % using the create_and_level_height_map method.
            %
            % 4. Loops through each ForceMap in the obj again, and for
            % each:
            %
            % a. If the user chose to skip preprocessed ForceMaps and the
            % current ForceMap is flagged as preprocessed, the loop moves
            % to the next ForceMap.
            %
            % b. Fits base lines for the current ForceMap using the
            % base_and_tilt method.
            %
            % c. Flags the current ForceMap as preprocessed.
            %
            % 5. Closes the progress waitbar upon completion.
            %
            % Use this function when you need to preprocess a set of
            % ForceMap objects by creating and leveling height maps and
            % fitting base lines. You can choose to skip ForceMaps that
            % have already been preprocessed to save time and keep previous
            % results.
            
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
                obj.FM{i}.base_and_tilt();
                waitbar(i/NLoop,h,sprintf('Preprocessing ForceMap %i/%i\nWrapping Up And Saving',i,NLoop));
                
                obj.FMFlag.Preprocessed(i) = 1;
            end
            
%             obj.save_experiment;
            
            close(h);
        end
        
        function force_map_analysis_fibril(obj,CPOption,EModOption,TemporaryLoadIn,UseTipInHertz)
            % This method has been deprecated! Should work as it is but
            % wont be updated anymore.
            %
            % force_map_analysis_fibril(obj, CPOption, EModOption,
            % TemporaryLoadIn, UseTipInHertz)
            %
            % This function is a deprecated method for analyzing force maps
            % of fibrils. Although it should still work, it will not be
            % updated in the future. The function handles preprocessing
            % tasks, contact point estimation, fibril diameter calculation,
            % and elastic modulus calculation for each force map in the
            % set. It also provides options to use Oliver or Hertz models
            % and the deconvolution of the cantilever tip if needed.
            %
            % Inputs:
            %
            % obj: An object containing the force maps to be analyzed.
            %
            % CPOption: String, representing the contact point estimation
            % method. It can be any valid option supported by the
            % cp_option_converter function (e.g., 'cnnzoomsweep').
            %
            % EModOption: String, specifying the elastic modulus
            % calculation method. Possible values are 'hertz', 'oliver',
            % 'both', or 'hertzcorrected'.
            %
            % TemporaryLoadIn: (Optional, default: true) Boolean, if true,
            % temporarily loads data in memory during analysis for cases
            % when the data is too big to be stored in memory all at once.
            %
            % UseTipInHertz: (Optional, default: true) Boolean, if true,
            % uses the deconvoluted cantilever tip in Hertz model
            % calculations.
            %
            % Outputs: No explicit outputs. The function updates the force
            % maps within the obj with analysis results. It also saves the
            % experiment and displays analyzed fibril data.
            
            
            if nargin < 5
                TemporaryLoadIn = true;
                UseTipInHertz = true;
            end
            
            obj.write_to_log_file('Analysis Function','force_map_analysis_fibril()','start')
            obj.write_to_log_file('Contact Point Option',CPOption)
            obj.write_to_log_file('EMod Option',EModOption)
            
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
            if isequal(lower(EModOption),'oliver') || isequal(lower(EModOption),'both') || UseTipInHertz
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
                if ~obj.KeepPythonFilesOpen && obj.PythonLoaderFlag && obj.BigDataFlag
                    obj.load_python_files_to_memory(i,[])
                end
                if TemporaryLoadIn && obj.BigDataFlag
                    obj.FM{i}.temporary_data_load_in(true);
                end
                
                waitbar(i/NLoop,h,sprintf('Processing Fibril %i/%i\nCreating and levelling Height Map',i,NLoop));
                obj.FM{i}.create_and_level_height_map
                obj.FM{i}.create_automatic_background_mask(1)
                
                Thresh = 1/2;
                AppRetSwitch = 2;
                obj.FM{i}.unselect_curves_by_fraction_of_max_data_points(Thresh,AppRetSwitch)
                
                waitbar(i/NLoop,h,sprintf('Processing Fibril %i/%i\nFitting Base Line',i,NLoop));
                if ~obj.FM{i}.BaseAndTiltFlag || ~obj.ForceMapAnalysisOptions.SkipAlreadyProcessedData
                    obj.FM{i}.base_and_tilt(obj.ForceMapAnalysisOptions.BaselineCorrection);
                    if i == 1
                        obj.write_to_log_file('Baseline and Tilt option','linear')
                    end
                end
                
                obj.FM{i}.calculate_fib_diam();
                
                % contact point estimation happens here
                waitbar(i/NLoop,h,sprintf('Processing Fibril %i/%i\nFinding Contact Point',i,NLoop));
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
                    if UseTipInHertz
                        obj.FM{i}.calculate_e_mod_hertz(CPOption,'parabolic',0,1,AllowXShift,CorrectSens,UseTipInHertz,0,obj.CantileverTips{obj.WhichTip(i)});
                    else
                        obj.FM{i}.calculate_e_mod_hertz(CPOption,'parabolic',0,1,AllowXShift,CorrectSens,UseTipInHertz,0);
                    end
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
                
                if TemporaryLoadIn && obj.BigDataFlag
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
                
                if ~obj.KeepPythonFilesOpen && obj.PythonLoaderFlag && obj.BigDataFlag
                    obj.clear_python_files_from_memory(i,[])
                end
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
                    if obj.EMod.Apex(i,j) > (median(obj.EMod.Apex(i,:),'omitnan')+2.5*iqr(obj.EMod.Apex(i,:))) || ...
                            obj.EMod.Apex(i,j) < (median(obj.EMod.Apex(i,:),'omitnan')-2.5*iqr(obj.EMod.Apex(i,:))) || ...
                            obj.FM{i}.ExclMask(obj.FM{i}.List2Map(obj.FM{i}.RectApexIndex(j),1),obj.FM{i}.List2Map(obj.FM{i}.RectApexIndex(j),2)) == 0
                        obj.EMod.Apex(i,j) = NaN;
                    elseif obj.EMod.Apex(i,j) < 0
                        obj.EMod.Apex(i,j) = NaN;
                    end
                end
            end
            
            obj.EMod.Mean = mean(obj.EMod.Apex,2,'omitnan');
            obj.EMod.STD = std(obj.EMod.Apex,[],2,'omitnan');
            obj.save_experiment;
            
            close(h);
            obj.write_to_log_file('','','end')
        end
        
        function force_map_analysis_general(obj,CPOption,EModOption,TemporaryLoadIn,UseTipInHertz,TiltCorrection)
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
            
            if nargin < 2
                current = what();
                cd(obj.ExperimentFolder)
                Time = char(datetime);
                Time = strrep(Time,' ','_');
                Time = strrep(Time,':','-');
                Time = strrep(Time,'.','_');
                FMA = obj.ForceMapAnalysisOptions;
                [~,FMA_ID] = fileparts(tempname);
                obj.broadcast_current_fma_id(FMA_ID);
                OptionsSnapshot = sprintf('ForceMapAnalysis-Options_%s.mat',Time);
                save(OptionsSnapshot,'FMA','FMA_ID');
                cd(current.path)
                CPOption = obj.ForceMapAnalysisOptions.ContactPointOption;
                EModOption = obj.ForceMapAnalysisOptions.EModOption.Type;
                TemporaryLoadIn = obj.ForceMapAnalysisOptions.TemporaryLoadIn;
                UseTipInHertz = obj.ForceMapAnalysisOptions.EModOption.Hertz.UseTipInHertz;
                TiltCorrection = obj.ForceMapAnalysisOptions.BaselineCorrection.TiltCorrection;
                obj.ReferenceSlopeFlag.Options = obj.ForceMapAnalysisOptions.SensitivityCorrection;
            elseif nargin < 4
                TemporaryLoadIn = true;
                UseTipInHertz = true;
                TiltCorrection = true;
            elseif nargin < 5
                UseTipInHertz = true;
                TiltCorrection = true;
            elseif nargin < 6
                TiltCorrection = true;
            end
            
            % Check if all FMAOptions are given. Else let user know how to
            % reset them
            [ContainsSame,IsSame] = check_structs_for_same_fields_recursively(obj.ForceMapAnalysisOptions,Experiment.set_default_fma_options);
            if ~ContainsSame && nargin < 2
                disp(sprintf('\n'))
                warning(sprintf(['Cannot run force map analysis!\n'...
                    'Your ForceMapAnalysisOptions do not contain all the necessary information!\n'...
                    'Most likely the settings have been updated and you need to reset and reapply your settings; try:\n\n'...
                    '>> E.ForceMapAnalysisOptions = Experiment.set_default_fma_options;\n\n'...
                    'and then choose your desired settings again \n\n'...
                    '>> E.set_force_map_analysis_options\n']));
                return
            end
            
            if obj.ForceMapAnalysisOptions.ResetFlagsSelectedAndCorrupted
                obj.reset_selected_and_corrupted_flags();
            end
            
            obj.write_to_log_file('Analysis Function','force_map_analysis_general()','start')
            obj.write_to_log_file('Contact Point Option',CPOption)
            obj.write_to_log_file('EMod Option',EModOption)
            
            h = waitbar(0,'setting up','Units','normalized','Position',[0.4 0.3 0.2 0.1]);
            NLoop = obj.NumForceMaps;
            if nargin < 2
                if obj.ForceMapAnalysisOptions.SkipAlreadyProcessedData
                    KeepFlagged = 'Yes';
                else
                    KeepFlagged = 'No';
                end
            else
                if sum(obj.FMFlag.ForceMapAnalysis) >= 1
                    KeepFlagged = questdlg(sprintf('Some maps have been processed already.\nDo you want to skip them and keep old results?'),...
                        'Processing Options',...
                        'Yes',...
                        'No',...
                        'No');
                else
                    KeepFlagged = 'No';
                end
            end
            
            % Setting and calculating preferred method of reference slope
            if nargin < 2
                obj.reference_slope_parser(1,true)
            else
                obj.reference_slope_parser(1,false)
            end
            
            if obj.ReferenceSlopeFlag.SetAllToValue
                RefSlopeOption = 'SetAllToValue';
            elseif obj.ReferenceSlopeFlag.UserInput
                RefSlopeOption = 'UserInput';
            elseif obj.ReferenceSlopeFlag.FromRefFM
                RefSlopeOption = 'FromRefFM';
            elseif obj.ReferenceSlopeFlag.FromArea
                RefSlopeOption = 'FromArea';
            elseif obj.ReferenceSlopeFlag.Adaptive
                RefSlopeOption = 'Adaptive';
            elseif obj.ReferenceSlopeFlag.KeepOriginal
                RefSlopeOption = 'KeepOriginal';
            elseif obj.ReferenceSlopeFlag.AutomaticFibril
                RefSlopeOption = 'AutomaticFibril';
            elseif obj.ReferenceSlopeFlag.Automatic
                RefSlopeOption = 'Automatic';
            end
            obj.write_to_log_file('Reference Slope Option',RefSlopeOption)
            
            % Deconvoluting cantilever tip(s)
            if isequal(lower(EModOption),'oliverpharr') || isequal(lower(EModOption),'both') || UseTipInHertz
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
            
            % Preallocate cell array to store error information
            errors = cell(1, NLoop);
            % Main loop for contact point estimation and E-Modulus calculation
            for i=1:NLoop
                try
                    if isequal(KeepFlagged,'Yes') && obj.FMFlag.ForceMapAnalysis(i) == 1
                        continue
                    end
                    
                    waitbar(i/NLoop,h,sprintf('Processing Force Map %i/%i\nReading out data...',i,NLoop));
                    
                    if ~obj.KeepPythonFilesOpen && obj.PythonLoaderFlag && obj.BigDataFlag
                        obj.load_python_files_to_memory(i,[])
                    end
                    if TemporaryLoadIn && obj.BigDataFlag
                        obj.FM{i}.temporary_data_load_in(true);
                    end
                    
                    waitbar(i/NLoop,h,sprintf('Processing Force Map %i/%i\nCreating and levelling Height Map',i,NLoop));
                    obj.FM{i}.create_and_level_height_map
                    obj.FM{i}.create_automatic_background_mask(1)
                    
                    Thresh = obj.ForceMapAnalysisOptions.UnselectCurveFragmentsThreshold;
                    AppRetSwitch = obj.ForceMapAnalysisOptions.UnselectCurveFragmentsAppRetSwitch;
                    obj.FM{i}.unselect_curves_by_fraction_of_max_data_points(Thresh,AppRetSwitch)
                    
                    waitbar(i/NLoop,h,sprintf('Processing Force Map %i/%i\nFitting Base Line',i,NLoop));
                    if ~obj.FM{i}.BaseAndTiltFlag || ~obj.ForceMapAnalysisOptions.SkipAlreadyProcessedData
                        obj.FM{i}.base_and_tilt(obj.ForceMapAnalysisOptions.BaselineCorrection);
                        if i == 1
                            obj.write_to_log_file('Baseline and Tilt option','linear')
                        end
                    end
                    
                    % contact point estimation happens here
                    waitbar(i/NLoop,h,sprintf('Processing Force Map %i/%i\nFinding Contact Point',i,NLoop));
                    obj.cp_option_converter(CPOption,i);
                    
                    % reference slope calculation happens here
                    waitbar(i/NLoop,h,sprintf('Processing Force Map %i/%i\nProcessing and calculating Reference Slope',i,NLoop));
                    obj.reference_slope_calculator(i);
                    
                    if obj.ForceMapAnalysisOptions.EModOption.Hertz.UseTopography
                        % Calculate local curvature maps
                        waitbar(i/NLoop,h,sprintf('Processing Force Map %i/%i\nProcessing and calculating Local Curvature Maps',i,NLoop));
                        if obj.ForceMapAnalysisOptions.EModOption.Hertz.UseTipInHertz
                            % At the moment Local Radius is determined by
                            % the mean cantilever depth dependend tip
                            % radius this is kind of a placeholder until
                            % proper contact area based radii can be fed to
                            % the routine.
                            if isempty(obj.FM{i}.search_channel('Indentation Depth'))
                                obj.FM{i}.create_and_level_height_map_by_current_cp(...
                                    obj.ForceMapAnalysisOptions.KeepOldResults);
                                [~,AllNames] = obj.FM{i}.search_channel('Contact Height');
                                CHChan = obj.FM{i}.get_channel(AllNames{end});
                                TmpHChan = obj.FM{i}.get_channel('Processed');
                                SIDChan = obj.FM{i}.create_standard_channel(...
                                    CHChan.Image - TmpHChan.Image,'Indentation Depth','m');
                            else
                                SIDChan = obj.FM{i}.get_channel(...
                                    obj.FM{i}.search_channel('Indentation Depth'));
                            end
                            THChan = obj.CantileverTips{obj.WhichTip(i)}.get_channel('Eroded Tip');
                            SHChan = obj.FM{i}.get_channel(...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.TopographyHeightChannel);
                            RadiiMap = obj.FM{i}.local_contact_area_ClassWrapper(...
                                THChan,SHChan,SIDChan,...
                                obj.ForceMapAnalysisOptions.KeepOldResults);
                            Radii = obj.FM{i}.convert_map_to_data_list(RadiiMap);
                            obj.FM{i}.localSurfaceFit_ClassWrapper(obj.ForceMapAnalysisOptions.EModOption.Hertz.TopographyHeightChannel,...
                                Radii,...
                                obj.ForceMapAnalysisOptions.KeepOldResults);
                        else
                            Radii = 0.25*obj.FM{i}.TipRadius;
                            obj.FM{i}.localSurfaceFit_ClassWrapper(obj.ForceMapAnalysisOptions.EModOption.Hertz.TopographyHeightChannel,...
                                Radii,...
                                obj.ForceMapAnalysisOptions.KeepOldResults);
                        end
                    end
                    
                    waitbar(i/NLoop,h,sprintf('Processing Force Map %i/%i\nCalculating E-Modulus',i,NLoop));
                    if isequal(lower(EModOption),'hertz') || isequal(lower(EModOption),'both') || isequal(lower(EModOption),'hertzcorrected')
                        if isequal(lower(EModOption),'hertzcorrected') || isequal(lower(EModOption),'both')
                            CorrectSens = true;
                        else
                            CorrectSens = false;
                        end
                        AllowXShift = obj.ForceMapAnalysisOptions.EModOption.Hertz.AllowXShift;
                        if UseTipInHertz
                            obj.FM{i}.calculate_e_mod_hertz(CPOption,...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.TipShape,...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.LowerCurveFraction,...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.UpperCurveFraction,...
                                AllowXShift,...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.CorrectSensitivity,...
                                UseTipInHertz,...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.UseTopography,...
                                obj.CantileverTips{obj.WhichTip(i)},...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.WeighPointsByInverseDistance,...
                                obj.ForceMapAnalysisOptions.SortHeightDataForFit,...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.FitDegreeForSneddonPolySurf,...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.LowerForceCutOff,...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.UpperForceCutOff,...
                                obj.ForceMapAnalysisOptions.SensitivityCorrection.SensitivityCorrectionMethod,...
                                obj.ForceMapAnalysisOptions.KeepOldResults,...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.TopographyHeightChannel,...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.ThinFilmMode,...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.ThinFilmThickness);
                        else
                            obj.FM{i}.calculate_e_mod_hertz(CPOption,...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.TipShape,...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.LowerCurveFraction,...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.UpperCurveFraction,...
                                AllowXShift,...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.CorrectSensitivity,...
                                UseTipInHertz,...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.UseTopography,...
                                [],...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.WeighPointsByInverseDistance,...
                                obj.ForceMapAnalysisOptions.SortHeightDataForFit,...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.FitDegreeForSneddonPolySurf,...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.LowerForceCutOff,...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.UpperForceCutOff,...
                                obj.ForceMapAnalysisOptions.SensitivityCorrection.SensitivityCorrectionMethod,...
                                obj.ForceMapAnalysisOptions.KeepOldResults,...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.TopographyHeightChannel,...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.ThinFilmMode,...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.ThinFilmThickness);
                        end
                        if i == 1
                            obj.write_to_log_file('Hertzian Tip-Shape',...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.TipShape)
                            obj.write_to_log_file('Hertzian UpperCurveFraction',...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.UpperCurveFraction)
                            obj.write_to_log_file('Hertzian LowerCurveFraction',...
                                obj.ForceMapAnalysisOptions.EModOption.Hertz.LowerCurveFraction)
                            obj.write_to_log_file('Allow X-Shift',AllowXShift)
                        end
                    end
                    if isequal(lower(EModOption),'oliverpharr') || isequal(lower(EModOption),'both')
                        obj.FM{i}.CP = obj.FM{i}.CP_HertzFitted; % DELETE ME AFTER THIS USE!
                        obj.FM{i}.calculate_e_mod_oliverpharr(obj.CantileverTips{obj.WhichTip(i)}.ProjectedTipArea,...
                            obj.ForceMapAnalysisOptions.EModOption.OliverPharr.CurvePercent,...
                            obj.ForceMapAnalysisOptions.KeepOldResults);
                        if i == 1
                            obj.write_to_log_file('OliverPharr UpperCurvePercent for linear fit',...
                                obj.ForceMapAnalysisOptions.EModOption.OliverPharr.CurvePercent)
                        end
                    end
                    
                    if obj.ForceMapAnalysisOptions.OptionalMetrics.AdhesionEnergy
                        obj.FM{i}.calculate_adhesion_energy_and_length(2,...
                            obj.ForceMapAnalysisOptions.KeepOldResults);
                        obj.write_to_log_file('Adhesion Energy STD-Threshold Multiplier','2')
                    end
                    if obj.ForceMapAnalysisOptions.OptionalMetrics.MaximumAdhesionForce
                        obj.FM{i}.calculate_adhesion_force(...
                            obj.ForceMapAnalysisOptions.KeepOldResults);
                    end
                    if obj.ForceMapAnalysisOptions.OptionalMetrics.DissipatedAndElasticEnergy
                        obj.FM{i}.calculate_dissipated_and_elastic_energy(...
                            obj.ForceMapAnalysisOptions.KeepOldResults);
                    end
                    if obj.ForceMapAnalysisOptions.OptionalMetrics.PeakIndentationAngle
                        obj.FM{i}.calculate_peak_indentation_angle(.5,...
                            obj.ForceMapAnalysisOptions.KeepOldResults);
                        obj.write_to_log_file('Upper Percent of Curve considered for Peak Indentation','50%')
                    end
                    if obj.ForceMapAnalysisOptions.OptionalMetrics.HeightMapByCurrentCP
                        try
                            obj.FM{i}.create_and_level_height_map_by_current_cp(...
                                obj.ForceMapAnalysisOptions.KeepOldResults);
                        catch
                        end
                    end
                    if obj.ForceMapAnalysisOptions.OptionalMetrics.CorrectedDZSlopes
                        if isequal(obj.ForceMapAnalysisOptions.OptionalMetrics.CorrectedDZSlopesOptions.FitRangeMode,'SameAsSensitivityCorrection')
                            obj.FM{i}.calculate_corrected_dzslopes(...
                                'SensitivityCorrectionMethod',obj.ForceMapAnalysisOptions.OptionalMetrics.CorrectedDZSlopesOptions.SensitivityCorrectionMethod,...
                                'FitRangeMode',obj.ForceMapAnalysisOptions.SensitivityCorrection.FitRangeMode,...
                                'FitRangeLowerValue',obj.ForceMapAnalysisOptions.SensitivityCorrection.FitRangeLowerValue,...
                                'FitRangeUpperValue',obj.ForceMapAnalysisOptions.SensitivityCorrection.FitRangeUpperValue,...
                                'FitRangeLowerFraction',obj.ForceMapAnalysisOptions.SensitivityCorrection.FitRangeLowerFraction,...
                                'FitRangeUpperFraction',obj.ForceMapAnalysisOptions.SensitivityCorrection.FitRangeUpperFraction,...
                                'AppRetSwitch',obj.ForceMapAnalysisOptions.SensitivityCorrection.AppRetSwitch,...
                                'KeepOldResults',obj.ForceMapAnalysisOptions.KeepOldResults...
                                );
                        else
                            obj.FM{i}.calculate_corrected_dzslopes(...
                                'SensitivityCorrectionMethod',obj.ForceMapAnalysisOptions.OptionalMetrics.CorrectedDZSlopesOptions.SensitivityCorrectionMethod,...
                                'FitRangeMode',obj.ForceMapAnalysisOptions.OptionalMetrics.CorrectedDZSlopesOptions.FitRangeMode,...
                                'FitRangeLowerValue',obj.ForceMapAnalysisOptions.OptionalMetrics.CorrectedDZSlopesOptions.FitRangeLowerValue,...
                                'FitRangeUpperValue',obj.ForceMapAnalysisOptions.OptionalMetrics.CorrectedDZSlopesOptions.FitRangeUpperValue,...
                                'FitRangeLowerFraction',obj.ForceMapAnalysisOptions.OptionalMetrics.CorrectedDZSlopesOptions.FitRangeLowerFraction,...
                                'FitRangeUpperFraction',obj.ForceMapAnalysisOptions.OptionalMetrics.CorrectedDZSlopesOptions.FitRangeUpperFraction,...
                                'AppRetSwitch',obj.ForceMapAnalysisOptions.OptionalMetrics.CorrectedDZSlopesOptions.AppRetSwitch,...
                                'KeepOldResults',obj.ForceMapAnalysisOptions.KeepOldResults...
                                );
                        end
                    end
                    if obj.ForceMapAnalysisOptions.EModOption.Hertz.PredictiveRSquare &&...
                            isequal(lower(EModOption),'hertz') || isequal(lower(EModOption),'both')
                        obj.FM{i}.calculate_predictive_rsquare_hertz(...
                            obj.ForceMapAnalysisOptions.EModOption.Hertz.CorrectSensitivity,...
                            obj.ForceMapAnalysisOptions.SensitivityCorrection.SensitivityCorrectionMethod,...
                            obj.ForceMapAnalysisOptions.SortHeightDataForFit,...
                            obj.ForceMapAnalysisOptions.EModOption.Hertz.AllowXShift,...
                            obj.ForceMapAnalysisOptions.KeepOldResults...
                            );
                    end
                    
                    waitbar(i/NLoop,h,sprintf('Processing Force Map %i/%i\nWrapping Up And Saving',i,NLoop));
                    
                    if TemporaryLoadIn && obj.BigDataFlag
                        obj.FM{i}.temporary_data_load_in(false);
                        if (i < NLoop &&...
                                obj.ForceMapAnalysisOptions.SaveAfterEachMap) ||...
                                (obj.ForceMapAnalysisOptions.SaveAfterEachMap &&...
                                ~obj.ForceMapAnalysisOptions.SaveWhenFinished)
                            obj.save_experiment;
                        end
                    end
                    
                    obj.FMFlag.ForceMapAnalysis(i) = 1;
                    if i==1
                        for k=2:NLoop
                            obj.FM{k}.CPFlag.CNNOpt = 1;
                            obj.FM{k}.MiniBatchSize = obj.FM{1}.MiniBatchSize;
                        end
                    end
                    if ~obj.KeepPythonFilesOpen && obj.PythonLoaderFlag && obj.BigDataFlag
                        obj.clear_python_files_from_memory(i,[])
                    end
                    
                catch ME
                    % If an error occurs, store it in the errors cell array
                    errors{i} = sprintf('Error in iteration %d ForceMap name %s: %s\nStack Trace:\n%s', i, obj.FM{i}.Name, ME.message, getReport(ME, 'extended', 'hyperlinks', 'off'));
                end
            end
            if obj.ForceMapAnalysisOptions.SaveWhenFinished
                obj.save_experiment;
            end
            
            close(h);
            obj.write_to_log_file('','','end')
            obj.FMAExitErrorStack = errors;
            for i = 1:NLoop
                if ~isempty(errors{i})
                    warning(errors{i});
                end
            end
            
            obj.broadcast_current_fma_id('none')
        end
        
        function fiber_segment_analysis(obj)
            
            
        end
        
        function OutStruct = characterize_fiber_like_polyline_segments(obj,varargin)
            % function OutStruct = characterize_fiber_like_polyline_segments(obj,varargin)
            %
            % BATCH PROCESS
            % Takes in SNAPPED polyline-segmented fiber-like structures and
            % determines several characterstics such as Height, FWHM, Area.
            %
            %
            % Required inputs
            % obj ... Experiment class object containing instances of 
            %         AFMBaseClass object that already have polyline
            %         segments and a 'Processed' channel. E.I and E.FM are
            %         processed
            %
            % Name-Value pairs
            % "WidthLocalWindowMeters" ... Width of the local profile
            %                           that is taken from the
            %                           orthogonal to the local
            %                           direction vector. Points within
            %                           this profile are further
            %                           processed. Unit meters e.g. for
            %                           200 nm type 200e-9
            % "SmoothingWindowSize" ... Integer>3. Number of data points in
            %                       the running window smoothing step.
            % "MinPeakDistanceMeters" ... Minimum Distance between two
            %                         peaks determined by findpeaks() in
            %                         order for them to actually count as
            %                         valid peaks. If distance is snmaller,
            %                         only the higher peak persists.
            %                         Unit meters e.g. for
            %                         50 nm type 50e-9
            % "LowerEndThreshold" ... Determines the lower end of the
            %                     fiber-like object. Can be given as a
            %                     fraction of fiber height or in meters. 
            % "ThresholdType" ... 'Fraction'(def), 'Meters'
            % "Verbose" ... logical; if true, the function will draw a
            %               visualization of what is happening.
            % "RecordMovieBool" ... logical. if true, the verbose figure
            %                       will be recorded to a videofile
            % "KeyFrames" ... Integer>0. Video will only record every
            %                 KeyFrames step.
            %                 FramesInVideo=TotalFrames/KeyFrames
            
            p = inputParser;
            p.FunctionName = "characterize_fiber_like_polyline_segments";
            p.CaseSensitive = false;
            p.PartialMatching = true;
            
            % Required inputs
            validobj = @(x)true;
            addRequired(p,"obj",validobj);
            
            % NameValue inputs
            defaultHeightChannelName = 'Processed';
            defaultWidthLocalWindowMeters = 800e-9;
            defaultSmoothingWindowSize = 41;
            defaultMinPeakDistanceMeters = 50e-9;
            defaultLowerEndThreshold = .1;
            defaultEllipseFitThreshold = .75;
            defaultThresholdType = 'Fraction';
            defaultVerbose = false;
            defaultRecordMovieBool = false;
            defaultKeyFrames = 3;
            validHeightChannelName = @(x)ischar(x);
            validWidthLocalWindowMeters = @(x)isnumeric(x)&&isscalar(x);
            validSmoothingWindowSize = @(x)isnumeric(x)&&mod(x,1)==0;
            validMinPeakDistanceMeters = @(x)isnumeric(x)&&isscalar(x);
            validLowerEndThreshold = @(x)isnumeric(x)&&isscalar(x);
            validEllipseFitThreshold = @(x)isnumeric(x)&&isscalar(x); 
            validThresholdType = @(x)any(validatestring(x,{'Fraction','Meters'}));
            validVerbose = @(x)islogical(x);
            validRecordMovieBool = @(x)islogical(x);
            validKeyFrames = @(x)isnumeric(x)&&mod(x,1)==0;
            addParameter(p,"HeightChannelName",defaultHeightChannelName,validHeightChannelName);
            addParameter(p,"WidthLocalWindowMeters",defaultWidthLocalWindowMeters,validWidthLocalWindowMeters);
            addParameter(p,"SmoothingWindowSize",defaultSmoothingWindowSize,validSmoothingWindowSize);
            addParameter(p,"MinPeakDistanceMeters",defaultMinPeakDistanceMeters,validMinPeakDistanceMeters);
            addParameter(p,"LowerEndThreshold",defaultLowerEndThreshold,validLowerEndThreshold);
            addParameter(p,"EllipseFitThreshold",defaultEllipseFitThreshold,validEllipseFitThreshold);
            addParameter(p,"ThresholdType",defaultThresholdType,validThresholdType);
            addParameter(p,"Verbose",defaultVerbose,validVerbose);
            addParameter(p,"RecordMovieBool",defaultRecordMovieBool,validRecordMovieBool);
            addParameter(p,"KeyFrames",defaultKeyFrames,validKeyFrames);
            
            parse(p,obj,varargin{:});
            
            % Assign parsing results to named variables
            obj = p.Results.obj;
            HeightChannelName = p.Results.HeightChannelName;
            WidthLocalWindowMeters = p.Results.WidthLocalWindowMeters;
            SmoothingWindowSize = p.Results.SmoothingWindowSize;
            MinPeakDistanceMeters = p.Results.MinPeakDistanceMeters;
            LowerEndThreshold = p.Results.LowerEndThreshold;
            ThresholdType = p.Results.ThresholdType;
            EllipseFitThreshold = p.Results.EllipseFitThreshold;
            Verbose = p.Results.Verbose;
            RecordMovieBool = p.Results.RecordMovieBool;
            KeyFrames = p.Results.KeyFrames;
            
            
            k = 1;
            for i=1:obj.NumForceMaps
                [OutStruct(k).Array,OutStruct(k).Struct,OutStruct(k).StructAll] = ...
                    obj.FM{i}.characterize_fiber_like_polyline_segments(...
                    'HeightChannelName',HeightChannelName,...
                    'WidthLocalWindowMeters',WidthLocalWindowMeters,...
                    'SmoothingWindowSize',SmoothingWindowSize,...
                    'MinPeakDistanceMeters',MinPeakDistanceMeters,...
                    'LowerEndThreshold',LowerEndThreshold,...
                    'EllipseFitThreshold',EllipseFitThreshold,...
                    'ThresholdType',ThresholdType,...
                    'Verbose',Verbose,...
                    'RecordMovieBool',RecordMovieBool,...
                    'KeyFrames',KeyFrames);
                if ~isempty(OutStruct(k).Array)
                    k = k + 1;
                end
            end
            for i=1:obj.NumAFMImages
                try
                [OutStruct(k).Array,OutStruct(k).Struct,...
                    OutStruct(k).StructAll] = ...
                    obj.I{i}.characterize_fiber_like_polyline_segments(...
                    'HeightChannelName',HeightChannelName,...
                    'WidthLocalWindowMeters',WidthLocalWindowMeters,...
                    'SmoothingWindowSize',SmoothingWindowSize,...
                    'MinPeakDistanceMeters',MinPeakDistanceMeters,...
                    'LowerEndThreshold',LowerEndThreshold,...
                    'EllipseFitThreshold',EllipseFitThreshold,...
                    'ThresholdType',ThresholdType,...
                    'Verbose',Verbose,...
                    'RecordMovieBool',RecordMovieBool,...
                    'KeyFrames',KeyFrames);
                if ~isempty(OutStruct(k).Array)
                    k = k + 1;
                end
                catch ME
                    warning(ME)
                end
            end
        end
        
        function Table = compile_experiment_content_to_table_sqldatabase(obj)
            % function Table =
            % compile_experiment_content_to_table_sqldatabase(obj)
            %
            % Cycle through ALL the contents of FM,I,CantileverTips and
            % create an entry of every single Pixel with columns being
            % every possible property with special rules for some
            % properties that are derived from
            % 'characterize_fiber_like_polyline_segments' or
            % 'create_overlay_group'
            %
            % This function compiles the content of an AFM experiment,
            % including force maps, AFM images, and cantilever tips, into a
            % single table for use in an SQL database. The output table has
            % columns for each property of the experiment content, with
            % specific rules applied to properties derived from
            % 'characterize_fiber_like_polyline_segments' and
            % 'create_overlay_group' functions.
            %
            % Input:
            %
            % obj: An object representing the AFM experiment content. It
            % should have properties for force maps (FM), AFM images (I),
            % and cantilever tips (CantileverTips), as well as other
            % properties required for the output table columns.
            %
            % Output:
            %
            % Table: A table with rows for each pixel in the AFM experiment
            % content and columns for each property. The table is
            % preallocated with the correct size and variable types based
            % on the input object properties. The table's properties
            % include variable units and descriptions derived from the
            % 'assign_units_and_fancy_names_to_properties' function.
            %
            % Properties and Data Types:
            %
            % The input object properties can be of the following data
            % types: float, integer, char, logical, string, or cell array
            % of strings. The properties are stored in the output table
            % with specific rules applied to derived properties, as
            % mentioned above.
            %
            % Note that 'single' and 'double' types are replaced with
            % 'singlenan' and 'doublenan' respectively, and 'char' is
            % replaced with 'string'.
            %
            % Blacklisted Properties:
            %
            % The following properties are blacklisted and will not be
            % included in the output table: 'FileVersion',
            % 'DataStoreFolder', 'RawDataFilePath', 'OpenZipFile',
            % 'Folder', 'HostOS', and 'HostName'.
            %
            % Usage: Table =
            % compile_experiment_content_to_table_sqldatabase(obj)
            
            h = waitbar(0,'setting up...','Name','Compiling SQL Table');
            
            obj.write_unrolled_channels_to_dynamic_properties(true);
            
            % Collect all the Column Categories:
            ColumnNames = [];
            TempColumnTypes = [];
            BlackList = {'FileVersion',...
                'DataStoreFolder',...
                'RawDataFilePath',...
                'OpenZipFile',...
                'Folder',...
                'HostOS',...
                'HostName',...
                };
            ExProps = {'FM','I','CantileverTips'};
            LoopNumber = [obj.NumForceMaps obj.NumAFMImages obj.NumCantileverTips];
            k = 1;
            for j=1:3
                if LoopNumber(j) == 0
                    continue
                end
                NPixels = obj.(ExProps{j}){1}.NumPixelsX.*obj.(ExProps{j}){1}.NumPixelsY;
                if NPixels > 0
                    PropertyStruct = obj.(ExProps{j}){1};
                    PropertyNames = fieldnames(PropertyStruct);
                    for i=1:length(PropertyNames)
                        if (isfloat(PropertyStruct.(PropertyNames{i})) ||...
                                isinteger(PropertyStruct.(PropertyNames{i})) ||...
                                ischar(PropertyStruct.(PropertyNames{i})) ||...
                                islogical(PropertyStruct.(PropertyNames{i})) ||...
                                isstring(PropertyStruct.(PropertyNames{i})) ||...
                                (iscellstr(PropertyStruct.(PropertyNames{i})) &&...
                                length(PropertyStruct.(PropertyNames{i}))==NPixels))
                            if ((isfloat(PropertyStruct.(PropertyNames{i})) ||...
                                    isinteger(PropertyStruct.(PropertyNames{i})) ||...
                                    islogical(PropertyStruct.(PropertyNames{i}))) &&...
                                    ~(isequal(size(PropertyStruct.(PropertyNames{i})),[1 NPixels]) ||...
                                    isequal(size(PropertyStruct.(PropertyNames{i})),[NPixels 1]) ||...
                                    isequal(size(PropertyStruct.(PropertyNames{i})),[1 1]))...
                                    )
                                continue
                            end
                            if any(strcmp(PropertyNames{i},BlackList))
                                continue
                            end
                            ColumnNames{k} = PropertyNames{i};
                            TempColumnTypes{k} = class(PropertyStruct.(PropertyNames{i}));
                            k = k + 1;
                        end
                    end
                end
            end
            
            % Delete non unique entries. This will require shared property
            % names between subclasses to be of the same type across
            % AFMImage and ForceMap
            [ColumnNames,iPre] = unique(ColumnNames);
            ColumnTypes = TempColumnTypes(iPre);
            
            % replace some of the types
            ColumnTypes = strrep(ColumnTypes,'single','singlenan');
            ColumnTypes = strrep(ColumnTypes,'double','doublenan');
            ColumnTypes = strrep(ColumnTypes,'char','string');
            
            % Preallocate the Table
            LoopNumber = [obj.NumForceMaps obj.NumAFMImages obj.NumCantileverTips];
            NumRowsTable = 0;
            for i=1:3
                if LoopNumber(i) == 0
                    continue
                end
                for j=1:LoopNumber(i)
                    NumRowsTable = NumRowsTable + ...
                        obj.(ExProps{i}){j}.NumPixelsX.*obj.(ExProps{i}){j}.NumPixelsY;
                end
            end
            [Units,FancyNames] = Experiment.assign_units_and_fancy_names_to_properties(ColumnNames);
            Table = table('Size',[NumRowsTable length(ColumnNames)],...
                'VariableTypes',ColumnTypes,'VariableNames',ColumnNames);
            
            Table.Properties.VariableUnits = Units;
            Table.Properties.VariableDescriptions = FancyNames;
            
            % Now loop over all the Classes and all their points and fill
            % up the Table
            m = 1;
            Barcounter = 1;
            for i=1:3
                if LoopNumber(i) == 0
                    continue
                end
                for j=1:LoopNumber(i)
                    waitbar(Barcounter/sum(LoopNumber),h,sprintf('writing Class %i/%i...',Barcounter,sum(LoopNumber)))
                    NumPixels = (obj.(ExProps{i}){j}.NumPixelsX*obj.(ExProps{i}){j}.NumPixelsY);
                    for k=1:length(ColumnNames)
                        if isprop(obj.(ExProps{i}){j},ColumnNames{k})
                            try
                                if isempty(obj.(ExProps{i}){j}.(ColumnNames{k})) && ~isequal(ColumnTypes{k},'cell')
                                    continue
                                elseif isempty(obj.(ExProps{i}){j}.(ColumnNames{k})) && isequal(ColumnTypes{k},'cell')
                                    DefaultCell = cell(NumPixels,1);
                                    [DefaultCell] = {''};
                                    Table(m:m+NumPixels-1,k) = DefaultCell;
                                elseif length(obj.(ExProps{i}){j}.(ColumnNames{k})) <= 1 &&...
                                        ~isequal(ColumnTypes{k},'string') &&...
                                        ~isempty(obj.(ExProps{i}){j}.(ColumnNames{k}))
                                    Table{m:m+NumPixels-1,k} = obj.(ExProps{i}){j}.(ColumnNames{k});
                                elseif isequal(ColumnTypes{k},'string') &&...
                                        ~isempty(obj.(ExProps{i}){j}.(ColumnNames{k}))
                                    Table{m:m+NumPixels-1,k} = string(obj.(ExProps{i}){j}.(ColumnNames{k}));
                                elseif isequal(ColumnTypes{k},'cell') &&...
                                        length(obj.(ExProps{i}){j}.(ColumnNames{k})) == NumPixels
                                    Table(m:m+NumPixels-1,k) = [obj.(ExProps{i}){j}.(ColumnNames{k})];
                                else
                                    Table{m:m+NumPixels-1,k} = reshape(obj.(ExProps{i}){j}.(ColumnNames{k}),[],1);
                                end
                            catch ME
                                i
                                j
                                k
                                m
                                rethrow(ME)
                            end
                        end
                    end
                    m = m + NumPixels;
                    Barcounter = Barcounter + 1;
                end
            end
            
            % Save Tbale as SQLDatabase
            
            close(h)
        end
        
        function image_analysis_base_on_even_background(obj,UpperLim,NIter)
            % function image_analysis_base_on_even_background(obj,
            % UpperLim, NIter)
            %
            % This function performs image analysis on AFM images by
            % subtracting the line fit histogram to obtain an even
            % background. The analysis is performed iteratively for a
            % specified number of iterations. If no 'ProcessedSimple'
            % channel exists, it will be created based on the 'Height
            % (Trace)' channel.
            %
            % Input:
            %
            % obj: An object representing the AFM experiment content. It
            % should have properties for AFM images (I) and their channels
            % (Channel), including 'Height (Trace)' and 'ProcessedSimple'
            % (if it exists).
            %
            % UpperLim (optional): A numeric value indicating the upper
            % limit for the line fit histogram subtraction. Default is 1.
            %
            % NIter (optional): An integer value representing the number of
            % iterations for the line fit histogram subtraction. Default is
            % 1.
            %
            % Output:
            %
            % The function does not return any output variables. It
            % modifies the input object by updating the 'ProcessedSimple'
            % channel (or creating it if it doesn't exist) with the
            % modified image data. Additionally, the 'hasProcessed'
            % property of the AFM image object is set to 1.
            %
            % Usage: image_analysis_base_on_even_background(obj)
            % image_analysis_base_on_even_background(obj, UpperLim)
            % image_analysis_base_on_even_background(obj, UpperLim, NIter)
            %
            % Example: To process AFM images with an upper limit of 1.5 and
            % 2 iterations, call the function as:
            % image_analysis_base_on_even_background(obj, 1.5, 2)
            
            
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
            % function
            % image_analysis_flatten_on_even_background_automatic(obj,
            % WindowSize, NIter)
            %
            % This function performs image analysis on AFM images by
            % iteratively flattening the image on an even background using
            % the 'subtract_line_fit_vertical_rov' function. If no
            % 'Processed' channel exists, it will be created based on the
            % 'Height (Trace)' channel.
            %
            % Input:
            %
            % obj: An object representing the AFM experiment content. It
            % should have properties for AFM images (I) and their channels
            % (Channel), including 'Height (Trace)' and 'Processed' (if it
            % exists).
            %
            % WindowSize (optional): A numeric value representing the
            % window size for the 'subtract_line_fit_vertical_rov'
            % function. Default is 0.2.
            %
            % NIter (optional): An integer value representing the number of
            % iterations for the image flattening process. Default is 3.
            %
            % Output:
            %
            % The function does not return any output variables. It
            % modifies the input object by updating the 'Processed' channel
            % (or creating it if it doesn't exist) with the modified image
            % data. Additionally, the 'hasProcessed' property of the AFM
            % image object is set to 1.
            %
            % Usage:
            %
            % image_analysis_flatten_on_even_background_automatic(obj)
            % image_analysis_flatten_on_even_background_automatic(obj,
            % WindowSize)
            % image_analysis_flatten_on_even_background_automatic(obj,
            % WindowSize, NIter)
            %
            % Example:
            %
            % To process AFM images with a window size of 0.3 and 4
            % iterations, call the function as:
            % image_analysis_flatten_on_even_background_automatic(obj, 0.3,
            % 4)
            
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
            % function
            % image_analysis_flatten_and_combine_trace_retrace_automatic(obj,
            % WindowSize, NIter)
            %
            % This function performs image analysis on AFM images by
            % flattening and combining the trace and retrace height
            % channels. It uses the 'subtract_line_fit_hist' and
            % 'subtract_line_fit_vertical_rov' functions for flattening and
            % creates a new channel 'R-T Combined' with the minimum of the
            % processed trace and retrace images.
            %
            % Input:
            %
            % obj: An object representing the AFM experiment content. It
            % should have properties for AFM images (I) and their channels
            % (Channel), including 'Height (Trace)' and 'Height (Retrace)'.
            %
            % WindowSize (optional): A numeric value representing the
            % window size for the 'subtract_line_fit_vertical_rov'
            % function. Default is 0.2.
            %
            % NIter (optional): An integer value representing the number of
            % iterations for the image flattening process. Default is 3.
            %
            % Output:
            %
            % The function does not return any output variables. It
            % modifies the input object by updating or creating the 'R-T
            % Combined' channel with the combined trace and retrace image
            % data.
            %
            % Usage:
            % image_analysis_flatten_and_combine_trace_retrace_automatic(obj)
            % image_analysis_flatten_and_combine_trace_retrace_automatic(obj,
            % WindowSize)
            % image_analysis_flatten_and_combine_trace_retrace_automatic(obj,
            % WindowSize, NIter)
            %
            % Example:
            %
            % To process AFM images with a window size of 0.3 and 4
            % iterations, call the function as:
            % image_analysis_flatten_and_combine_trace_retrace_automatic(obj,
            % 0.3, 4)
            
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
            % function image_analysis_mask_background(obj, Thresh)
            %
            % This function applies a background mask to the processed AFM
            % images by thresholding. The mask is generated using the
            % 'mask_background_by_threshold' function and creates a new
            % channel 'Background Mask' in the object's Channel property.
            %
            % Input:
            %
            % obj: An object representing the AFM experiment content. It
            % should have properties for AFM images (I) and their channels
            % (Channel), including a 'Processed' channel.
            %
            % Thresh (optional): A numeric value representing the threshold
            % for background masking. If not provided, the function will
            % attempt to determine the threshold automatically.
            %
            % Output:
            %
            % The function does not return any output variables. It
            % modifies the input object by creating the 'Background Mask'
            % channel with the masked image data.
            %
            % Usage:
            %
            % image_analysis_mask_background(obj)
            % image_analysis_mask_background(obj, Thresh)
            %
            % Example:
            %
            % To apply a background mask with a threshold of 5, call the
            % function as: image_analysis_mask_background(obj, 5)
            
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
            % aux-function for CPOption choice function
            % cp_option_converter(obj, CPOption, i, RefFM)
            %
            % This function is an auxiliary function for handling various
            % contact point (CP) estimation methods for AFM force maps.
            % Based on the input option (CPOption), it applies the selected
            % method to the force map and updates the contact point
            % information.
            %
            % Input:
            %
            % obj: An object representing the AFM experiment content. It
            % should have properties for force maps (FM), contact points
            % (CP), and various contact point estimation methods.
            %
            % CPOption: A string specifying the contact point estimation
            % method. Possible values are: 'snap-in', 'none', 'rov',
            % 'hardsurface', 'gof', 'old', 'combo', 'manual', 'fast',
            % 'dropout', 'zoom', 'zoomdropout', 'zoomsweep', and
            % 'curveorigin'.
            %
            % i: An integer representing the index of the force map to be
            % processed.
            %
            % RefFM (optional): A boolean flag indicating if the reference
            % force map (RefFM) should be used. If not provided, the
            % default value is false.
            %
            % Output:
            %
            % The function does not return any output variables. It
            % modifies the input object by updating the contact point
            % information based on the selected method.
            %
            % Usage: cp_option_converter(obj, CPOption, i)
            % cp_option_converter(obj, CPOption, i, RefFM)
            %
            % Example: To apply the 'snap-in' method to the first force
            % map, call the function as: cp_option_converter(obj,
            % 'snap-in', 1)
            
            NumPasses = 20;
            
            if nargin < 4
                RefFM = false;
            end
            if RefFM == false
                if isequal(lower(CPOption),'snap-in')
                    obj.FM{i}.estimate_cp_snap_in();
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                        obj.FM{i}.CP_SnapIn(:,2) = 0;
                    end
                end
                if isequal(lower(CPOption),'none')
                    obj.FM{i}.CP_None = zeros(obj.FM{i}.NCurves,2);
                    obj.FM{i}.CP = obj.FM{i}.CP_None;
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                        obj.FM{i}.CP_None(:,2) = 0;
                    end
                end
                if isequal(lower(CPOption),'rov')
                    obj.FM{i}.estimate_cp_rov();
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                        obj.FM{i}.CP_RoV(:,2) = 0;
                    end
                end
                if isequal(lower(CPOption),'hardsurface')
                    obj.FM{i}.estimate_cp_hardsurface();
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                        obj.FM{i}.CP_HardSurface(:,2) = 0;
                    end
                end
                if isequal(lower(CPOption),'gof')
                    obj.FM{i}.estimate_cp_gof();
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                        obj.FM{i}.CP_GoF(:,2) = 0;
                    end
                end
                if isequal(lower(CPOption),'old')
                    obj.FM{i}.estimate_cp_old();
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                        obj.FM{i}.CP_Old(:,2) = 0;
                    end
                end
                if isequal(lower(CPOption),'combo')
                    if obj.FM{i}.CPFlag.RoV == 0
                        obj.FM{i}.estimate_cp_rov();
                    end
                    if obj.FM{i}.CPFlag.GoF == 0
                        obj.FM{i}.estimate_cp_gof();
                    end
                    obj.FM{i}.estimate_cp_combined();
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                        obj.FM{i}.CP_Combo(:,2) = 0;
                    end
                end
                if isequal(lower(CPOption),'manual')
                    obj.FM{i}.estimate_cp_manually;
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                    end
                end
                if isequal(lower(CPOption),'fast')
                    obj.FM{i}.estimate_cp_cnn(obj.CP_CNN,'Fast');
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                    end
                end
                if isequal(lower(CPOption),'dropout')
                    obj.FM{i}.estimate_cp_cnn(obj.DropoutNet,'Dropout',NumPasses);
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                        obj.FM{i}.CP_Dropout(:,2) = 0;
                    end
                end
                if isequal(lower(CPOption),'zoom')
                    obj.FM{i}.estimate_cp_cnn(obj.CP_CNN,'Zoom');
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                    end
                end
                if isequal(lower(CPOption),'zoomdropout')
                    obj.FM{i}.estimate_cp_cnn(obj.CP_CNN,'zoomdropout',NumPasses);
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                    end
                end
                if isequal(lower(CPOption),'zoomsweep')
                    obj.FM{i}.estimate_cp_cnn(obj.CP_CNN,'Zoomsweep',NumPasses);
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                    end
                end
                if isequal(lower(CPOption),'curveorigin')
                    obj.FM{i}.estimate_cp_curve_origin();
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                        obj.FM{i}.CP_CurveOrigin(:,2) = 0;
                    end
                end
            elseif RefFM == true
                if isequal(lower(CPOption),'snap-in')
                    obj.FM{i}.estimate_cp_snap_in();
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                        obj.FM{i}.CP_SnapIn(:,2) = 0;
                    end
                end
                if isequal(lower(CPOption),'none')
                    obj.FM{i}.CP_None = zeros(obj.FM{i}.NCurves,2);
                    obj.FM{i}.CP = obj.FM{i}.CP_None;
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                        obj.FM{i}.CP_None(:,2) = 0;
                    end
                end
                if isequal(lower(CPOption),'hardsurface')
                    obj.FM{i}.estimate_cp_hardsurface();
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                        obj.FM{i}.CP_HardSurface(:,2) = 0;
                    end
                end
                if isequal(lower(CPOption),'rov')
                    obj.RefFM{i}.estimate_cp_rov();
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                        obj.FM{i}.CP_RoV(:,2) = 0;
                    end
                end
                if isequal(lower(CPOption),'gof')
                    obj.RefFM{i}.estimate_cp_gof();
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                        obj.FM{i}.CP_GoF(:,2) = 0;
                    end
                end
                if isequal(lower(CPOption),'old')
                    obj.RefFM{i}.estimate_cp_old();
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                        obj.FM{i}.CP_Old(:,2) = 0;
                    end
                end
                if isequal(lower(CPOption),'combo')
                    obj.RefFM{i}.estimate_cp_rov();
                    obj.RefFM{i}.estimate_cp_gof();
                    obj.RefFM{i}.estimate_cp_combined();
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                        obj.FM{i}.CP_Combo(:,2) = 0;
                    end
                end
                if isequal(lower(CPOption),'manual')
                    obj.RefFM{i}.estimate_cp_manually;
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                    end
                end
                if isequal(lower(CPOption),'fast')
                    obj.RefFM{i}.estimate_cp_cnn(obj.CP_CNN,'Fast');
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                    end
                end
                if isequal(lower(CPOption),'dropout')
                    obj.RefFM{i}.estimate_cp_cnn(obj.DropoutNet,'Dropout',NumPasses);
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                        obj.FM{i}.CP_Dropout(:,2) = 0;
                    end
                end
                if isequal(lower(CPOption),'zoom')
                    obj.RefFM{i}.estimate_cp_cnn(obj.CP_CNN,'Zoom');
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                    end
                end
                if isequal(lower(CPOption),'zoomdropout')
                    obj.RefFM{i}.estimate_cp_cnn(obj.CP_CNN,'zoomdropout',NumPasses);
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                    end
                end
                if isequal(lower(CPOption),'zoomsweep')
                    obj.RefFM{i}.estimate_cp_cnn(obj.CP_CNN,'Zoomsweep',NumPasses);
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                    end
                end
                if isequal(lower(CPOption),'curveorigin')
                    obj.FM{i}.estimate_cp_curve_origin();
                    if obj.ForceMapAnalysisOptions.ContactPointJustXAxis
                        obj.FM{i}.CP(:,2) = 0;
                        obj.FM{i}.CP_CurveOrigin(:,2) = 0;
                    end
                end
            end
                
            if obj.FM{i}.CPFlag.CNNopt
                obj.broadcast_CNN_optimization(i);
            end
            
        end
        
        function apply_segmentation_to_overlay_group(obj,GroupMemberInstance)
            % function apply_segmentation_to_overlay_group(obj,
            % GroupMemberInstance)
            %
            % This function applies a segmentation from a given group
            % member instance to all other members of the same overlay
            % group. It identifies the corresponding members within the AFM
            % experiment object (obj) and applies the segmentation to force
            % maps (FM) and AFM images (I).
            %
            % Input:
            %
            % obj: An object representing the AFM experiment content. It
            % should have properties for force maps (FM) and AFM images
            % (I), as well as their corresponding names (ForceMapNames,
            % AFMImageNames).
            %
            % GroupMemberInstance: An object representing a group member,
            % which includes the segmentation to be applied. It should have
            % properties for the overlay group (OverlayGroup) and its name
            % (Name).
            %
            % Output:
            %
            % The function does not return any output variables. It
            % modifies the input objects by applying the segmentation from
            % the given group member instance to all other members of the
            % same overlay group.
            %
            % Usage:
            %
            % apply_segmentation_to_overlay_group(obj, GroupMemberInstance)
            %
            % Example:
            %
            % To apply the segmentation from a specific group member
            % instance to all other members of the same overlay group, call
            % the function as: apply_segmentation_to_overlay_group(obj,
            % GroupMemberInstance)
            
            if ~GroupMemberInstance.OverlayGroup.hasOverlayGroup
                warning('The class instance you input has no overlay group')
                return
            end
            
            for i=1:GroupMemberInstance.OverlayGroup.Size
                if strcmp(GroupMemberInstance.OverlayGroup.Names{i},GroupMemberInstance.Name)
                    continue
                end
                Index = find(strcmp({GroupMemberInstance.OverlayGroup.Names{i}},obj.ForceMapNames));
                if ~isempty(Index)
                    GroupMemberInstance.apply_segmentation_to_other_baseclass(obj.FM{Index})
                end
                Index = find(strcmp({GroupMemberInstance.OverlayGroup.Names{i}},obj.AFMImageNames));
                if ~isempty(Index)
                    GroupMemberInstance.apply_segmentation_to_other_baseclass(obj.I{Index})
                end
            end
            
        end
        
        function create_proximity_mapping_for_overlay_group_segments(obj,GroupMemberInstance,AllToOneBool)
            % create_proximity_mapping_for_overlay_group_segments(obj,GroupMemberInstance,AllToOneBool)
            %
            % Takes an Overlay Group an tries to map Segments with same
            % Names and SubSectionNames onto each other. Mapping is done by
            % least distance of polyline vertices to the vertices of the
            % transformed segments. The mapping is injective but not
            % bijective.
            %
            % This function creates a proximity mapping for segments in an
            % overlay group. Segments with the same names and subsection
            % names are mapped to each other based on the least distance of
            % polyline vertices to the vertices of the transformed
            % segments. The mapping is injective but not bijective. The
            % function provides two modes: All-to-One mode or sequential
            % mode.
            %
            % Input:
            %
            % obj: An object representing the AFM experiment content. It
            % should have properties for force maps (FM) and AFM images
            % (I), as well as their corresponding names (ForceMapNames,
            % AFMImageNames). GroupMemberInstance: An object representing a
            % group member in the overlay group. It should have properties
            % for the overlay group (OverlayGroup) and its name (Name).
            %
            % AllToOneBool: A boolean value (0 or 1). When set to 1, the
            % function will map all group members to the first group member
            % (All-to-One mode). When set to 0, the function will map each
            % group member to the next one sequentially (sequential mode).
            % Default value is 0.
            %
            % Output:
            %
            % The function does not return any output variables. It
            % modifies the input objects by creating proximity mappings for
            % segments in the overlay group.
            %
            % Usage:
            %
            % create_proximity_mapping_for_overlay_group_segments(obj,
            % GroupMemberInstance, AllToOneBool)
            %
            % Example:
            %
            % To create a proximity mapping for segments in an overlay
            % group, call the function as:
            % create_proximity_mapping_for_overlay_group_segments(obj,
            % GroupMemberInstance, AllToOneBool)
            
            if ~GroupMemberInstance.OverlayGroup.hasOverlayGroup
                warning('The class instance you input has no overlay group')
                return
            end
            
            if nargin < 3
                AllToOneBool = 0;
            end
            
            h = waitbar(0,'Proximity Mapping');
            
            if AllToOneBool
                for i=1:GroupMemberInstance.OverlayGroup.Size
                    waitbar(i/(GroupMemberInstance.OverlayGroup.Size - 1),h,sprintf('Mapping Group Member %i to %i',i,1));
                    if strcmp(GroupMemberInstance.OverlayGroup.Names{i},GroupMemberInstance.Name)
                        continue
                    end
                    Index = find(strcmp({GroupMemberInstance.OverlayGroup.Names{i}},obj.ForceMapNames));
                    if ~isempty(Index)
                        GroupMemberInstance.create_proximity_mapping_for_segments(obj.FM{Index})
                    end
                    Index = find(strcmp({GroupMemberInstance.OverlayGroup.Names{i}},obj.AFMImageNames));
                    if ~isempty(Index)
                        GroupMemberInstance.create_proximity_mapping_for_segments(obj.I{Index})
                    end
                end
            else
                for i=1:GroupMemberInstance.OverlayGroup.Size - 1
                    waitbar(i/(GroupMemberInstance.OverlayGroup.Size - 1),h,sprintf('Mapping Group Member %i to %i',i+1,i));
                    Index1 = find(strcmp({GroupMemberInstance.OverlayGroup.Names{i}},obj.ForceMapNames));
                    Index2 = find(strcmp({GroupMemberInstance.OverlayGroup.Names{i+1}},obj.ForceMapNames));
                    if ~isempty(Index1)
                        ObjOne = obj.FM{Index1};
                    end
                    if ~isempty(Index2)
                        ObjTwo = obj.FM{Index2};
                    end
                    Index1 = find(strcmp({GroupMemberInstance.OverlayGroup.Names{i}},obj.AFMImageNames));
                    Index2 = find(strcmp({GroupMemberInstance.OverlayGroup.Names{i+1}},obj.AFMImageNames));
                    if ~isempty(Index1)
                        ObjOne = obj.I{Index1};
                    end
                    if ~isempty(Index2)
                        ObjTwo = obj.I{Index2};
                    end
                    ObjOne.create_proximity_mapping_for_segments(ObjTwo)
                end
            end
            close(h)
        end
        
        function reset_proximity_mapping_for_overlay_group_segments(obj,GroupMemberInstance)
            % function
            % reset_proximity_mapping_for_overlay_group_segments(obj,
            % GroupMemberInstance)
            %
            % This function resets the proximity mapping for segments in an
            % overlay group. It undoes the previously applied proximity
            % mappings for all group members in the overlay group. The
            % function iterates through all the members of the overlay
            % group and resets their proximity mappings.
            %
            % Input:
            %
            % obj: An object representing the AFM experiment content. It
            % should have properties for force maps (FM) and AFM images
            % (I), as well as their corresponding names (ForceMapNames,
            % AFMImageNames). GroupMemberInstance: An object representing a
            % group member in the overlay group. It should have
            %
            % properties for the overlay group (OverlayGroup) and its name
            % (Name).
            %
            % Output:
            %
            % The function does not return any output variables. It
            % modifies the input objects by resetting the proximity
            % mappings for segments in the overlay group.
            %
            % Usage:
            %
            % reset_proximity_mapping_for_overlay_group_segments(obj,
            % GroupMemberInstance)
            %
            % Example:
            %
            % To reset the proximity mapping for segments in an overlay
            % group, call the function as:
            % reset_proximity_mapping_for_overlay_group_segments(obj,
            % GroupMemberInstance)
            
            
            if ~GroupMemberInstance.OverlayGroup.hasOverlayGroup
                warning('The class instance you input has no overlay group')
                return
            end
            
            for i=1:GroupMemberInstance.OverlayGroup.Size
                Index = find(strcmp({GroupMemberInstance.OverlayGroup.Names{i}},obj.ForceMapNames));
                if ~isempty(Index)
                    obj.FM{Index}.reset_proximity_mapping;
                end
                Index = find(strcmp({GroupMemberInstance.OverlayGroup.Names{i}},obj.AFMImageNames));
                if ~isempty(Index)
                    obj.I{Index}.reset_proximity_mapping;
                end
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
        
        % SMFS section
                   
        function SMFS_testing_function(obj,XMin,XMax,YMin,YMax,ii)
            % Function to quickly loop over all force maps for testing and
            % debugging
             for ii=1:obj.NumForceMaps
            %for ii=2
             %   obj.FM{ii}.initialize_flags
                  % obj.FM{ii}.fc_fc_measurement_prop;
                %   obj.FM{ii}.fc_pulling_length_MAD
%                   obj.FM{ii}.fc_adhesion_energy_idxpulllength
%                   obj.FM{ii}.fc_adhesion_energy_threshold
              %      obj.FM{ii}.fc_find_idx
          %       obj.FM{ii}.fc_adh_force_max
          %  obj.FM{ii}.fc_sinoidal_fit    
        %  obj.FM{ii}.fc_fit_based_yData
            
              obj.FM{ii}.fc_print_raw
            
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
        
        function SMFS_preprocessing(obj)
            % SMFS_preprocessing: A function to run a bundle of other 
            % typically required functions for further analysis
            % obj.preprocessing
                       
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
            
            % force map loop
            for ii=1:obj.NumForceMaps
            % for ii=1:20 % debugging
                if isequal(KeepFlagged,'Yes') && obj.SMFSFlag.Preprocessed(ii) == 1
                    continue
                end
                waitbar(ii/NLoop,h,sprintf('Preprocessing ForceMap %i/%i\nProcessing force curves',ii,NLoop));            
                obj.FM{ii}.fc_measurement_prop             
                waitbar(ii/NLoop,h,sprintf('Preprocessing ForceMap %i/%i\nWrapping Up And Saving',ii,NLoop));
                                
                
                obj.SMFSFlag.Preprocessed(ii) = 1;
            end
            close(h);
        end
                  
        function SMFS_sinoidal_fit(obj)
            % SMFS_sinoidal_fit: This function allows to conduct an automated presorting of the force curves 
            % The function flags force curves and whole force maps that are
            % non-functionalize
            % Needed function: obj.preprocessing
            
             h = waitbar(0,'setting up','Units','normalized','Position',[0.4 0.3 0.2 0.1]);
            NLoop = length(obj.ForceMapNames);
            if sum(obj.SMFSFlag.Fit) >= 1
                KeepFlagged = questdlg(sprintf('Some maps have been processed already.\nDo you want to skip them and keep old results?'),...
                    'Processing Options',...
                    'Yes',...
                    'No',...
                    'No');
            else
                KeepFlagged = 'No';
            end
            
            % Change into the Folder of Interest
            cd(obj.ExperimentFolder) % Move into the folder
            % Create folders for saving the produced figures
            foldername='FM_sinoidal_fits';    % for debugging                
            mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
            currpath=fullfile(obj.ExperimentFolder,foldername);
            cd(currpath);
            
            % Loop over the imported force maps
            for ii=1:obj.NumForceMaps
            %for ii=1:2 % Debugging
                waitbar(ii/NLoop,h,sprintf('Preprocessing ForceMap %i/%i\nProcessing force curves',ii,NLoop));
                obj.FM{ii}.fc_sinoidal_fit
                obj.FM{ii}.fc_fit_based_yData
                waitbar(ii/NLoop,h,sprintf('Preprocessing ForceMap %i/%i\nWrapping Up And Saving',ii,NLoop));
                sprintf('Force Map No. %d of %d',ii,obj.NumForceMaps) % Gives current Force Map Position             
            end
            close(h);
        end
                       
        function SMFS_presorting(obj)
            % SMFS_presorting: This function allows to conduct an automated presorting of the force curves 
            % The function flags force curves and whole force maps that are
            % non-functionalize
            % Needed function: obj.preprocessing
            
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
            % for ii=3 % Debugging
                if isequal(KeepFlagged,'Yes') && obj.SMFSFlag.Preprocessed(ii) == 1
                    continue
                end   
                waitbar(ii/NLoop,h,sprintf('Preprocessing ForceMap %i/%i\nProcessing force curves',ii,NLoop));
                waitbar(ii/NLoop,h,sprintf('Preprocessing ForceMap %i/%i\nWrapping Up And Saving',ii,NLoop));
                sprintf('Force Map No. %d of %d',ii,obj.NumForceMaps) % Gives current Force Map Position               
                obj.FM{ii}.fc_estimate_cp_hardsurface
                obj.FM{ii}.fc_selection_threshold
                    if nnz(obj.FM{ii}.SMFSFlag.Min)<20 % Only if more than 20 force curves fulfil the citeria the whole force map is considered successfully functionalized
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
            elseif nargin<3
                XMin= -inf;     % Limit of the X-axis in meters (m)
                XMax= inf;      % Limit of the X-axis in meters (m)
                YMin= -inf;     % Limit of the Y-axis in Newtons (N)
                YMax= inf;      % Limit of the Y-axis in Newtons (N)
            end
            % Change into the Folder of Interest
            cd(obj.ExperimentFolder) % Move into the folder 
            % Create folders for saving the produced figures
            foldername='FM_Fig';    % Defines the folder name
            mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
            currpath=fullfile(obj.ExperimentFolder,foldername);
            cd(currpath); 
            
            % Loop over the imported force maps
            %for ii=1:obj.NumForceMaps
            for ii=1 % Debugging
               % Command window output
               sprintf('Force Map No. %d of %d',ii,obj.NumForceMaps) % Gives current Force Map Position
               % Run the chosen functions
               obj.FM{ii}.fc_visual_selection(XMin,XMax,YMin,YMax);     
               %obj.save_experiment;        % Save immediately after each force curve
            end    
        end
      
        function SMFS_analysis(obj,XMin,XMax,YMin,YMax,NumFcMax,NumFcUncorrupt,hh)
            % This function allows to analyse different force curve
            % criteria, i.e. pulling length, adhesion energy. Furthermore,
            % all analysed force curves are plotted and the determined
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
            %% Folder
            % Change into the Folder of Interest
            cd(obj.ExperimentFolder) % Move into the folder
            % Create folders for saving the produced figures
            %foldername='FM_fcAnalysis';    % for debugging
            foldername='FM_analysed';    % Defines the folder name
            mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
            currpath=fullfile(obj.ExperimentFolder,foldername);
            cd(currpath);
            %% loop
            for hh=1:obj.NumForceMaps
            %for hh=3:6 % Debugging
               sprintf('Force map No. %d',hh);
               % Print force curves containing label for the pulling length
               % and colored area for the adhesion energy                              
               % Baseline correction
               obj.FM{hh}.fc_based_ret_correction
               % Snap-in
               obj.FM{hh}.fc_snap_in_length_MAD
               % Pulling length
               obj.FM{hh}.fc_pulling_length_MAD
               % 50 nm limit index
               obj.FM{hh}.fc_find_idx
               % Maximum adhesion force
               obj.FM{hh}.fc_adh_force_max
               % Adhesion energy
               obj.FM{hh}.fc_adhesion_energy_idxpulllength
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
                     
            % Input variable adaptation
            if nargin < 2
                XMin= -inf;     % Limit of the X-axis in meters (m)  
                XMax= inf;      % Limit of the X-axis in meters (m)
                YMin= -inf;     % Limit of the Y-axis in Newtons (N)   
                YMax= inf;      % Limit of the Y-axis in Newtons (N)
            end
            if nargin < 4         
                XMax= 50e-9;      % Limit of the X-axis in meters (m)
                YMax= 100e-12;      % Limit of the Y-axis in Newtons (N)
            end
            
            % Change into the Folder of Interest
            cd(obj.ExperimentFolder) % Move into the folder 
            % Create folders for saving the produced figures
            %foldername='FM_test';    % for debugging
            foldername='FM_unsorted';    % Defines the folder name
            mkdir(obj.ExperimentFolder,foldername);  % Creates for each force map a folder where the corresponding figures are stored in
            currpath=fullfile(obj.ExperimentFolder,foldername);
            cd(currpath); 
            
            % Loop over the imported force maps
            for ii=1:obj.NumForceMaps
            %for ii=1:10 % Debugging
            % Presort condition 
              %  if ~obj.SMFSFlag(ii)   % Selects all flagged 1 force maps
                %if obj.SMFSFlag(ii)     % Selects all flagged 0 force maps
               %     continue
               %end
               % Command window output
               sprintf('Force Map No. %d of %d',ii,obj.NumForceMaps) % Gives current Force Map Position
               % Run the chosen functions
               obj.FM{ii}.fc_print_raw(XMin,XMax,YMin,YMax);     
             %  obj.save_experiment;        % Save immediately after each force curve
            end    
        end
        
        function SMFS_print_sort(obj,StartDate,EndDate,XMin,XMax,YMin,YMax)
            % SMFS_print_sort: A function to plot all force curves of all
            % force maps sorted by different properties 
            % Comment: Date format is: 'YYYY.MM.DD'
            % Required function: obj.FM{ii}.fc_measurement_prop, obj.FM{ii}.estimate_cp_hardsurface
         
            % Change into the Folder of Interest
            cd(obj.ExperimentFolder) % Move into the folder 
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
            % Loop over the imported force maps
             for ii=1:obj.NumForceMaps
             %for ii=3
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
           
        function SMFS_analysis_selction_fit(obj,VelocityValue,SubstrateValue,EnvCondValue,ChipCantValue,ChipboxValue,LinkerValue)
                
            
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
             
        function SMFS_analysis_selection(obj,VelocityValue,SubstrateValue,EnvCondValue,ChipCantValue,ChipboxValue,LinkerValue)
                
                
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
                        else
                            obj.FM{ii}.PullingLength(obj.FM{ii}.PullingLength==0)=nan;
                            ConcatArrayPullingLength1=horzcat(ConcatArrayPullingLength1,obj.FM{ii}.PullingLength(:));
                            ConcatArrayPullingLength2=vertcat(ConcatArrayPullingLength2,obj.FM{ii}.PullingLength(:));
                            obj.FM{ii}.RetAdhEnergy_IdxMethod(obj.FM{ii}.RetAdhEnergy_IdxMethod==0)=nan;
                            ConcatArrayRetAdhEnergy1=horzcat(ConcatArrayRetAdhEnergy1,obj.FM{ii}.RetAdhEnergy_IdxMethod(:));
                            ConcatArrayRetAdhEnergy2=vertcat(ConcatArrayRetAdhEnergy2,obj.FM{ii}.RetAdhEnergy_IdxMethod(:));
                        end
                    end
                end
                
                % Change into the Folder of Interest
                cd(obj.ExperimentFolder) % Move into the folder
                % Create folders for saving the produced figures
                %foldername='FM_test';    % for debugging
                foldername='FM_analysis';    % Defines the folder name
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
                
                %% Figure
                % Figure 1
                h_fig1=figure(1);
                h_fig1.Color='white'; % changes the background color of the figure
                h_fig1.Units='normalized'; % Defines the units
                h_fig1.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
                h_fig1.PaperOrientation='landscape';
                h_fig1.Name=figname;
                % boxplot
                bo1=boxplot(ConcatArrayPullingLength1);
                % title
                t1=title(fulltitle1);
                % Axes
                ax1 = gca; % current axes
                ax1.FontSize = 16;
                ax1.XLabel.String = 'Force map number (1)';
                ax1.XLabel.FontSize = 20;
                ax1.YLabel.String = 'PullingLength (m)';
                ax1.YLabel.FontSize = 20;
                % Save figure
                %%% Define the name for the figure title
                partname='boxplot';
                % fullname=sprintf('%s%s',figname,partname);
                fullname1=strcat(figname,{'_'},parttitle1,{'_'},partname);
                fullname1=char(fullname1);
                %%% Save the current figure in the current folder
                print(h_fig1,fullname1,'-dpng');
            
                % Figure 2
                h_fig2=figure(2);
                h_fig2.Color='white'; % changes the background color of the figure
                h_fig2.Units='normalized'; % Defines the units
                h_fig2.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
                h_fig2.PaperOrientation='landscape';
                h_fig2.Name=figname;
                % boxplot
                bo2= boxplot(ConcatArrayRetAdhEnergy1);
                % title
                t2=title(fulltitle2);
                % Axes
                ax2 = gca; % current axes
                ax2.FontSize = 16;
                ax2.XLabel.String = 'Force map number (1)';
                ax2.XLabel.FontSize = 20;
                ax2.YLabel.String = 'Adhesion Energy (J)';
                ax2.YLabel.FontSize = 20;
                % Save figure
                %%% Define the name for the figure title
                partname='boxplot';
                % fullname=sprintf('%s%s',figname,partname);
                fullname2=strcat(figname,{'_'},parttitle2,{'_'},partname);
                fullname2=char(fullname2);
                %%% Save the current figure in the current folder
                print(gcf,fullname2,'-dpng');
                
                % Figure 3
                h_fig3=figure(3);
                h_fig3.Color='white'; % changes the background color of the figure
                h_fig3.Units='normalized'; % Defines the units
                h_fig3.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
                h_fig3.PaperOrientation='landscape';
                h_fig3.Name=figname;
                % Histogram
                h3=histogram(ConcatArrayPullingLength2);
                h3.BinWidth=1e-9;
                h3.FaceAlpha=1;
                h3.FaceColor='b';
                h3.EdgeColor='k';
                % title
                t1=title(fulltitle1);
                % axes
                ax3=gca;
                ax3.FontSize = 16;
                ax3.XLabel.String = 'PullingLength (m)';
                ax3.XLabel.FontSize = 20;
                ax3.YLabel.String = 'Frequency count (1)';
                ax3.YLabel.FontSize = 20;               
                hold on;
                % Save figure
                %%% Define the name for the figure title
                partname='histo';
                % fullname=sprintf('%s%s',figname,partname);
                fullname3=strcat(figname,{'_'},parttitle1,{'_'},partname);
                fullname3=char(fullname3);
                %%% Save the current figure in the current folder
                print(gcf,fullname3,'-dpng');
                
                % Figure 4
                h_fig4=figure(4);
                h_fig4.Color='white'; % changes the background color of the figure
                h_fig4.Units='normalized'; % Defines the units
                h_fig4.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
                h_fig4.PaperOrientation='landscape';
                h_fig4.Name=figname;
                % Histogram
                h4=histogram(ConcatArrayRetAdhEnergy2);
                h4.BinWidth=1e-19;
                h4.FaceAlpha=1;
                h4.FaceColor='g';
                h4.EdgeColor='k';
                % title
                t2=title(fulltitle2);
                % axes
                ax4=gca;
                ax4.FontSize = 16;
                ax4.YScale='log';
                ax4.XLabel.String = 'Adhesion Energy (J)';
                ax4.XLabel.FontSize = 20;
                ax4.YLabel.String = 'Frequency count (1)';
                ax4.YLabel.FontSize = 20;               
                hold on;
                % Save figure
                %%% Define the name for the figure title
                partname='histo';
                % fullname=sprintf('%s%s',figname,partname);
                fullname4=strcat(figname,{'_'},parttitle2,{'_'},partname);
                fullname4=char(fullname4);
                %%% Save the current figure in the current folder
                print(gcf,fullname3,'-dpng');
                %% House keeping
                close all
        end
  
       

        function SMFS_fine_figure(obj,XMin,XMax,YMin,YMax,ii)
            % Function to plot individual fine figures for publication
            if nargin < 2
                XMin= -inf;
                XMax= inf;
                YMin= -inf;
                YMax= inf;
            end
            
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
        
     
            
    end
    methods
        % Methods for data visualization spanning all the data
        
        function show_image(obj)
            % TODO: implement ui elements for customization
            
            % Check if all FMAOptions are given. Else let user know how to
            % reset them
            [ContainsSame,IsSame] = check_structs_for_same_fields_recursively(obj.ShowImageSettings,Experiment.set_default_show_image_settings);
            if ~ContainsSame
                disp(sprintf('\n'))
                warning(sprintf(['Cannot run "show image"!\n'...
                    'Your ShowImageSettings do not contain all the necessary information!\n'...
                    'Most likely the settings have been updated.\n'...
                    'Setting everything to current default...']));
                obj.ShowImageSettings = Experiment.set_default_show_image_settings;
            end
            
            h.ColorMode = set_default_color_options();
            
            h.Fig = figure('Name',sprintf('%s',obj.ExperimentName),...
                'Units','pixels',...
                'Position',[200 100 1024 512],...
                'Color',h.ColorMode(obj.ShowImageSettings.ColorIndex).Background);
            
            set(gcf,'SizeChangedFcn',@(a,b)update_interface)
            
            h.Backdrop = uicontrol('style','text',...
                'String','',...
                'Units','normalized',...
                'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Backdrop,...
                'Position',[.845 .0 .2 1]);
            
            initialize_starting_class_and_popups
            
            set_main_tab_ui
            set_tab_ui
            set_volume_data_tab_ui
            set_results_tab_ui
            
            initialize_starting_state
            draw_channel_1
            
            function initialize_starting_class_and_popups()
                
                [h.ClassPopUp,h.ClassIndex] = obj.string_of_existing_class_instances();
                h.NumClasses = length(h.ClassPopUp);
                if h.NumClasses < obj.ShowImageSettings.DefaultChannel1Index
                    obj.ShowImageSettings.DefaultChannel1Index = 1;
                end
                if h.NumClasses < obj.ShowImageSettings.DefaultChannel2Index
                    obj.ShowImageSettings.DefaultChannel2Index = 1;
                end
                
                h.Class{1} = obj.get_class_instance(...
                    h.ClassIndex(obj.ShowImageSettings.DefaultChannel1Index,:));
                h.Class{2} = obj.get_class_instance(...
                    h.ClassIndex(obj.ShowImageSettings.DefaultChannel2Index,:));
                h.PopUp1 = h.Class{1}.string_of_existing();
                h.PopUp2 = h.Class{2}.string_of_existing();
                if length(h.PopUp1) < obj.ShowImageSettings.DefaultChannel1SubIndex
                    obj.ShowImageSettings.DefaultChannel1SubIndex = 2;
                end
                if length(h.PopUp2) < obj.ShowImageSettings.DefaultChannel2SubIndex
                    obj.ShowImageSettings.DefaultChannel2SubIndex = 2;
                end
                
            end
            
            function set_tab_ui()
                
                h.MenuTabs(1) = uicontrol('style','togglebutton',...
                    'String','Main',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'position',[.85 .42 .05 .06],...
                    'Tooltip','Main Menu and Line Profiles',...
                    'Value',obj.ShowImageSettings.MainMenuValue,...
                    'Callback',@switch_to_main_menu);
                h.MenuTabs(2) = uicontrol('style','togglebutton',...
                    'String','Volume Data',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'position',[.9 .42 .05 .06],...
                    'Tooltip','Volume Data such as Force-Distance curves',...
                    'Value',obj.ShowImageSettings.VolumeMenuValue,...
                    'Callback',@switch_to_volume_menu);
                h.MenuTabs(3) = uicontrol('style','togglebutton',...
                    'String','Results',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'position',[.95 .42 .05 .06],...
                    'Tooltip','Results and Statistics based on Segmentations',...
                    'Value',obj.ShowImageSettings.ResultsMenuValue,...
                    'Callback',@switch_to_results_menu);
                
            end
            
            function set_main_tab_ui()
                
                h.B(1) = uicontrol('style','togglebutton',...
                    'String','Cross Section',...
                    'units','normalized',...
                    'Visible',obj.ShowImageSettings.MainMenuValue,...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'position',[.875 .36 .1 .04],...
                    'Callback',@cross_section_toggle);
                
                h.B(4) = uicontrol('style','text',...
                    'String','Channel 1',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'position',[.85 .95 .15 .04]);
                
                
                h.B(16) = uicontrol('style','popupmenu',...
                    'String',h.ClassPopUp,...
                    'Value',obj.ShowImageSettings.DefaultChannel1Index,...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'position',[.85 .9 .15 .05],...
                    'Callback',@draw_channel_1);
                
                h.B(2) = uicontrol('style','popupmenu',...
                    'String',h.PopUp1,...
                    'Value',obj.ShowImageSettings.DefaultChannel1SubIndex,...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'position',[.85 .85 .15 .05],...
                    'Callback',@draw_channel_1);
                
                h.B(5) = uicontrol('style','text',...
                    'String','Channel 2',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'position',[.85 .7 .15 .04]);
                
                h.B(17) = uicontrol('style','popupmenu',...
                    'String',h.ClassPopUp,...
                    'Value',obj.ShowImageSettings.DefaultChannel2Index,...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'position',[.85 .65 .15 .05],...
                    'Callback',@draw_channel_2);
                
                h.B(3) = uicontrol('style','popupmenu',...
                    'String',h.PopUp2,...
                    'Value',obj.ShowImageSettings.DefaultChannel2SubIndex,...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'position',[.85 .6 .15 .05],...
                    'Callback',@draw_channel_2);
                
                h.B(6) = uicontrol('style','pushbutton',...
                    'String','Save Figure',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'Visible',1,...
                    'position',[.875 .08 .1 .04],...
                    'Callback',@save_figure_to_file);
                
                h.B(7) = uicontrol('style','checkbox',...
                    'String','...with white background',...
                    'Value',mod(obj.ShowImageSettings.ColorIndex-1,2),...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'Visible',1,...
                    'position',[.875 .05 .1 .03],...
                    'Callback',@changed_color);
                
                h.B(8) = uicontrol('style','slider',...
                    'Value',1,...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Position',[.85 .83 .12 .02],...
                    'Callback',@changed_slider);
                
                h.B(9) = uicontrol('style','slider',...
                    'Value',0,...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Position',[.85 .81 .12 .02],...
                    'Callback',@changed_slider);
                
                h.B(10) = uicontrol('style','slider',...
                    'Value',1,...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Position',[.85 .58 .12 .02],...
                    'Callback',@changed_slider);
                
                h.B(11) = uicontrol('style','slider',...
                    'Value',0,...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Position',[.85 .56 .12 .02],...
                    'Callback',@changed_slider);
                
                h.B(12) = uicontrol('style','text',...
                    'String','Max',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Position',[.97 .83 .03 .02]);
                
                h.B(13) = uicontrol('style','text',...
                    'String','Min',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Position',[.97 .81 .03 .02]);
                
                h.B(14) = uicontrol('style','text',...
                    'String','Max',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Position',[.97 .58 .03 .02]);
                
                h.B(15) = uicontrol('style','text',...
                    'String','Min',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Position',[.97 .56 .03 .02]);
                
                h.B(18) = uicontrol('Style','checkbox',...
                    'String','Both Channels',...
                    'Value',obj.ShowImageSettings.BothCrossSections,...
                    'Tooltip','Green, if both Channels have the same size scaling',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Visible',obj.ShowImageSettings.MainMenuValue,...
                    'Position',[.875 .33 .1 .03],...
                    'Callback',@checked_both_cross_sections);
                
                h.B(19) = uicontrol('style','checkbox',...
                    'String','Upscale Images',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'Value',obj.ShowImageSettings.IsUpscaled,...
                    'Visible',obj.ShowImageSettings.MainMenuValue,...
                    'position',[.875 .15 .1 .03],...
                    'Callback',@upscale_images);
                
                h.B(20) = uicontrol('style','checkbox',...
                    'String','Lock Channels',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'Value',obj.ShowImageSettings.LockChannels,...
                    'Visible',obj.ShowImageSettings.MainMenuValue,...
                    'position',[.875 .24 .1 .03],...
                    'Callback',@lock_channels);
                
                h.B(21) = uicontrol('style','checkbox',...
                    'String','Lock Scalebars',...
                    'Value',obj.ShowImageSettings.LockScalebars,...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'Visible',obj.ShowImageSettings.MainMenuValue,...
                    'position',[.875 .18 .1 .03],...
                    'Callback',@lock_scalebars);
                
                h.B(22) = uicontrol('style','checkbox',...
                    'String','Statistical CMapping',...
                    'Value',obj.ShowImageSettings.StatisticalCMap,...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'Visible',obj.ShowImageSettings.MainMenuValue,...
                    'position',[.875 .21 .1 .03],...
                    'Callback',@statistical_cmapping);
                
                h.B(23) = uicontrol('style','edit',...
                    'String','',...
                    'units','normalized',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'position',[.85 .75 .035 .05],...
                    'Callback',@set_scale);
                
                h.B(24) = uicontrol('style','edit',...
                    'String','',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'position',[.925 .75 .035 .05],...
                    'Callback',@set_scale);
                
                h.B(25) = uicontrol('style','edit',...
                    'String','',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'position',[.85 .5 .035 .05],...
                    'Callback',@set_scale);
                
                h.B(26) = uicontrol('style','edit',...
                    'String','',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'position',[.925 .5 .035 .05],...
                    'Callback',@set_scale);
                
                h.B(27) = uicontrol('style','text',...
                    'String',{'Min','[]'},...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Position',[.885 .75 .04 .05]);
                
                h.B(28) = uicontrol('style','text',...
                    'String',{'Max','[]'},...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Position',[.96 .75 .04 .05]);
                
                h.B(29) = uicontrol('style','text',...
                    'String',{'Min','[]'},...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Position',[.885 .5 .04 .05]);
                
                h.B(30) = uicontrol('style','text',...
                    'String',{'Max','[]'},...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Position',[.96 .5 .04 .05]);
                
                h.B(31) = uicontrol('Style','checkbox',...
                    'String','Use Overlay',...
                    'Value',obj.ShowImageSettings.UseOverlay,...
                    'Tooltip','Green, if both Channels share an overlay',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Visible',obj.ShowImageSettings.MainMenuValue,...
                    'Position',[.875 .3 .1 .03],...
                    'Callback',@changed_use_overlay);
                
                    h.ToggledMainMenuIndizes = [18 19 20 21 22 31 1];
                    
            end
            
            function set_volume_data_tab_ui()
                
                % Volume Data Buttons
                h.FVM(1) = uicontrol('style','radiobutton',...
                    'String','Baseline corrected',...
                    'Tooltip','Subtracts a baseline from the data',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Value',obj.ShowImageSettings.BaselineCorrected,...
                    'Visible',obj.ShowImageSettings.VolumeMenuValue,...
                    'Position',[.875 .33 .1 .03],...
                    'Callback',@changed_baseline_corrected);
                h.FVM(2) = uicontrol('style','radiobutton',...
                    'String','Tip Height',...
                    'Tooltip','Subtracts the baselined vDeflection from the height data',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Value',obj.ShowImageSettings.TipHeight,...
                    'Visible',obj.ShowImageSettings.VolumeMenuValue,...
                    'Position',[.875 .30 .1 .03],...
                    'Callback',@changed_tip_height);
                h.FVM(3) = uicontrol('style','radiobutton',...
                    'String','Contact Point Shifted',...
                    'Tooltip','Defines the currently chosen Contact Point as the origin',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Value',obj.ShowImageSettings.ContactPointShifted,...
                    'Visible',obj.ShowImageSettings.VolumeMenuValue,...
                    'Position',[.875 .27 .1 .03],...
                    'Callback',@changed_contact_point_shifted);
                h.FVM(4) = uicontrol('style','radiobutton',...
                    'String','Show Hertz Fit',...
                    'Tooltip','Additionally plots the fitted Hertz model',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Value',obj.ShowImageSettings.ShowHertzFit,...
                    'Visible',obj.ShowImageSettings.VolumeMenuValue,...
                    'Position',[.875 .24 .1 .03],...
                    'Callback',@changed_show_hertz_fit);
                h.FVM(5) = uicontrol('style','radiobutton',...
                    'String','Use Corrected Sens.',...
                    'Tooltip','Uses the corrected sensitivity from chosen reference slope method',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Value',obj.ShowImageSettings.UseCorrectedSensitivity,...
                    'Visible',obj.ShowImageSettings.VolumeMenuValue,...
                    'Position',[.875 .21 .1 .03],...
                    'Callback',@changed_use_corrected_sensitivity);
                h.FVM(6) = uicontrol('style','radiobutton',...
                    'String','Plot Time',...
                    'Tooltip','Plots the Deflection over Time',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Value',obj.ShowImageSettings.PlotTime,...
                    'Visible',obj.ShowImageSettings.VolumeMenuValue,...
                    'Position',[.875 .18 .1 .03],...
                    'Callback',@changed_plot_time);
                h.FVM(8) = uicontrol('style','popupmenu',...
                    'String',{'N','m','V'},...
                    'Tooltip','Choose in terms of which unit the deflection is to be plotted',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Value',obj.ShowImageSettings.vDeflectionUnitIndex,...
                    'Visible',obj.ShowImageSettings.VolumeMenuValue,...
                    'Position',[.925 .15 .05 .03],...
                    'Callback',@changed_vdeflection_unit);
                h.FVM(9) = uicontrol('style','text',...
                    'String','VDef. Unit',...
                    'Tooltip','Choose in terms of which unit the deflection is to be plotted',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Visible',obj.ShowImageSettings.VolumeMenuValue,...
                    'Position',[.875 .15 .05 .03]);
                h.FVM(7) = uicontrol('style','radiobutton',...
                    'String','Extended Information',...
                    'Tooltip','Writes additional info into the plot',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Value',obj.ShowImageSettings.ExtendedInformation,...
                    'Visible',obj.ShowImageSettings.VolumeMenuValue,...
                    'Position',[.875 .12 .1 .03],...
                    'Callback',@changed_extended_information);
                
                    h.ToggledVolumeMenuIndizes = [1 2 3 4 5 6 7 8 9];
                
            end
            
            function set_results_tab_ui()
                
                % Results Data Buttons
                h.Res(6) = uicontrol('style','pushbutton',...
                    'String','Export Data',...
                    'Tooltip','Export Data to workspace aswell as .mat-file into user chosen folder',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'Visible',obj.ShowImageSettings.ResultsMenuValue,...
                    'position',[.875 .12 .1 .04],...
                    'Callback',@export_data);
                % THIS NEEDS TO STAY NR. 8!
                h.Res(8) = uicontrol('style','text',...
                    'String','',...
                    'Units','normalized',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Backdrop,...
                    'Position',[0 0 .85 .34]);
                % THIS NEEDS TO STAY NR. 8!
                h.Res(1) = uicontrol('style','radiobutton',...
                    'String','Use Snapped',...
                    'Tooltip','Use algorithmically snapped segments instead of handdrawn',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Value',obj.ShowImageSettings.UseSnapped,...
                    'Visible',obj.ShowImageSettings.ResultsMenuValue,...
                    'Position',[.03 .26 .1 .03],...
                    'Callback',@changed_use_snapped);
                h.Res(2) = uicontrol('style','radiobutton',...
                    'String','Plot whole Group',...
                    'Tooltip',"Compare all members of Channel 1's Overlay Group",...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Value',obj.ShowImageSettings.PlotGroup,...
                    'Visible',obj.ShowImageSettings.ResultsMenuValue,...
                    'Position',[.03 .23 .1 .03],...
                    'Callback',@changed_plot_group);
                h.Res(3) = uicontrol('style','radiobutton',...
                    'String','Ignore Zeros',...
                    'Tooltip','Replace zeros in data set with NaNs',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Value',obj.ShowImageSettings.IgnoreZeros,...
                    'Visible',obj.ShowImageSettings.ResultsMenuValue,...
                    'Position',[.03 .20 .1 .03],...
                    'Callback',@changed_ignore_zeros);
                h.Res(4) = uicontrol('style','radiobutton',...
                    'String','Just Channel 1',...
                    'Tooltip','Plot just Channel 1',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Value',obj.ShowImageSettings.JustChannel1,...
                    'Visible',obj.ShowImageSettings.ResultsMenuValue,...
                    'Position',[.03 .17 .1 .03],...
                    'Callback',@changed_just_channel_1);
                h.Res(5) = uicontrol('style','radiobutton',...
                    'String','Just Channel 2',...
                    'Tooltip','Plot just Channel 2',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Units','normalized',...
                    'Value',obj.ShowImageSettings.JustChannel2,...
                    'Visible',obj.ShowImageSettings.ResultsMenuValue,...
                    'Position',[.03 .14 .1 .03],...
                    'Callback',@changed_just_channel_2);
                h.Res(9) = uicontrol('style','text',...
                    'String','',...
                    'Units','normalized',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Position',[.445 0.03 .002 .29]);
                h.Res(10) = uicontrol('style','text',...
                    'String','Data',...
                    'Units','normalized',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Position',[0.13 .3 .2 .03]);
                h.Res(11) = uicontrol('style','text',...
                    'String','Plot Settings',...
                    'Units','normalized',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Position',[0.55 .3 .2 .03]);
                h.Res(12) = uicontrol('style','text',...
                    'String','',...
                    'Units','normalized',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'Position',[0.843 0 .002 .34]);
                
                % Plotting settings
                h.Res(13) = uicontrol('style','pushbutton',...
                    'String','Gramm Plot Options',...
                    'Tooltip','Set Gramm Plot settings. Credit to https://github.com/piermorel/gramm',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'Visible',obj.ShowImageSettings.ResultsMenuValue,...
                    'position',[.55 .25 .2 .03],...
                    'Callback',@set_gramm_options);
                h.Res(14) = uicontrol('style','pushbutton',...
                    'String','Redraw Plot',...
                    'Tooltip','Redraw plot with last used data but current plot settings',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'Visible',obj.ShowImageSettings.ResultsMenuValue,...
                    'position',[.47 .05 .09 .03],...
                    'Callback',@redraw_gramm_plot);
                h.Res(15) = uicontrol('style','pushbutton',...
                    'String','Overlay Current Plot',...
                    'Tooltip','Overlay current plot with current data and settings',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'Visible',obj.ShowImageSettings.ResultsMenuValue,...
                    'position',[.6 .05 .09 .03],...
                    'Callback',@draw_over_current_gramm_plot);
                h.Res(7) = uicontrol('style','pushbutton',...
                    'String','Get data and plot',...
                    'Tooltip','Get data, assign according to settings and plot with current plotting settings',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'Visible',obj.ShowImageSettings.ResultsMenuValue,...
                    'position',[.73 .05 .09 .03],...
                    'Callback',@draw_results);
                
                % Data settings
                
                PossibleOptions = {'none','Channel 1','Channel 2','Channel Name',...
                    'Segment Name','AFM File Name','Height','CharFib: WidthHalfHeight',...
                    'CharFib: Prominence','CharFib: WidthHalfProminence',...
                    'CharFib: Area','CharFib: WidthBase','CharFib: AspectRatioHalfHeight',...
                    'CharFib: AspectRatioBaseHeight','CharFib: AreaDerivedDiameter'};
                
                h.Res(16) = uicontrol('style','text',...
                    'String','X',...
                    'Tooltip','Which data should be assigned as X data',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'Visible',obj.ShowImageSettings.ResultsMenuValue,...
                    'position',[.32 .26 .1 .03]);
                h.Res(17) = uicontrol('style','popupmenu',...
                    'Value',obj.ShowImageSettings.GrammX,...
                    'String',PossibleOptions,...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'Visible',obj.ShowImageSettings.ResultsMenuValue,...
                    'position',[.32 .23 .1 .03],...
                    'Callback',@changed_gramm_x);
                h.Res(18) = uicontrol('style','text',...
                    'String','Y',...
                    'Tooltip','Which data should be assigned as Y data',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'Visible',obj.ShowImageSettings.ResultsMenuValue,...
                    'position',[.32 .19 .1 .03]);
                h.Res(19) = uicontrol('style','popupmenu',...
                    'Value',obj.ShowImageSettings.GrammY,...
                    'String',PossibleOptions,...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'Visible',obj.ShowImageSettings.ResultsMenuValue,...
                    'position',[.32 .16 .1 .03],...
                    'Callback',@changed_gramm_y);
                h.Res(20) = uicontrol('style','text',...
                    'String','Group',...
                    'Tooltip','Which data should be assigned as Grouping data',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'Visible',obj.ShowImageSettings.ResultsMenuValue,...
                    'position',[.32 .12 .1 .03]);
                h.Res(21) = uicontrol('style','popupmenu',...
                    'Value',obj.ShowImageSettings.GrammGroup,...
                    'String',PossibleOptions,...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'Visible',obj.ShowImageSettings.ResultsMenuValue,...
                    'position',[.32 .09 .1 .03],...
                    'Callback',@changed_gramm_group);
                h.Res(22) = uicontrol('style','text',...
                    'String','Facets',...
                    'Tooltip','Which data should be assigned as Facet/Subplot data',...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'Visible',obj.ShowImageSettings.ResultsMenuValue,...
                    'position',[.32 .05 .1 .03]);
                h.Res(23) = uicontrol('style','popupmenu',...
                    'Value',obj.ShowImageSettings.GrammFacets,...
                    'String',PossibleOptions,...
                    'ForegroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText,...
                    'BackgroundColor',h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons,...
                    'units','normalized',...
                    'Visible',obj.ShowImageSettings.ResultsMenuValue,...
                    'position',[.32 .02 .1 .03],...
                    'Callback',@changed_gramm_facets);
                
                
                
                h.ToggledResultMenuIndizes = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23];
                
            end
            
            function initialize_starting_state()
                
                h.Channel1Max = 1;
                h.Channel1Min = 0;
                h.Channel2Max = 1;
                h.Channel2Min = 0;
                
                h.Channel{1} = 'none';
                h.Channel{2} = 'none';
                
                h.CurChannel1Idx = h.B(2).Value;
                h.CurChannel2Idx = h.B(3).Value;
                
                h.ResetOnNextZoom = false;
                h.MainLine = [];
                h.ChildLine = [];
                h.hasCrossSection = 0;
                h.OldMultiplier{1} = [];
                h.OldMultiplier{2} = [];
                h.SubChannelName{1} = [];
                h.SubChannelName{2} = [];
                h.VolumeStruct = struct('Point',[]);
                
                
                if obj.ShowImageSettings.DefaultChannel1SubIndex == 0 ||...
                        obj.ShowImageSettings.DefaultChannel2SubIndex == 0
                    obj.ShowImageSettings.HasChannel2 = 0;
                    [~,DefIndex] = h.Class{1}.get_channel('Processed');
                    if isempty(DefIndex)
                        DefIndex = 2;
                    else
                        DefIndex = DefIndex + 1;
                    end
                    h.B(2).Value = DefIndex;
                    obj.ShowImageSettings.DefaultChannel2SubIndex = DefIndex;
                    h.B(3).Value = 1;
                    obj.ShowImageSettings.DefaultChannel2SubIndex = 1;
                end
                
                h.CurrentClassName{1} = h.Class{1}.Name;
                h.CurrentClassName{2} = h.Class{1}.Name;
                h.OnePass = false;
                
                % Related to Results Struct
                h.ResultStruct.SelectedSegmentIndex = cell(1,2);
                h.ResultStruct.SelectedSegmentIndex{1} = 0;
                h.ResultStruct.SelectedSegmentIndex{2} = 0;
                h.ResultStruct.DisabledSegmentIndizes{1} = [];
                h.ResultStruct.DisabledSegmentIndizes{2} = [];
                h.ResultStruct.Results(1).Name = h.Class{1}.Name;
                h.ResultStruct.Results(2).Name = h.Class{1}.Name;
                h.ResultStruct.Results(1).ChannelType = 'none';
                h.ResultStruct.Results(2).ChannelType = 'none';
                h.ResultStruct.Results(1).Data = [];
                h.ResultStruct.Results(2).Data = [];
                h.ResultStruct.Results(1).SegmentNames = [];
                h.ResultStruct.Results(2).SegmentNames = [];
                h.ResultStruct.Results(1).Unit = [];
                h.ResultStruct.Results(2).Unit = [];
                h.ResultStruct.Results(1).SubSegmentCell = [];
                h.ResultStruct.Results(2).SubSegmentCell = [];
                h.ResultStruct.Results(1).SegmentCell = [];
                h.ResultStruct.Results(2).SegmentCell = [];
                h.ResultStruct.Results(1).Mask = [];
                h.ResultStruct.Results(2).Mask = [];
                
                if obj.ShowImageSettings.MainMenuValue
                        switch_to_main_menu
                elseif obj.ShowImageSettings.VolumeMenuValue
                        switch_to_volume_menu
                elseif obj.ShowImageSettings.ResultsMenuValue
                        switch_to_results_menu
                end
                
            end
            
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
                if obj.ShowImageSettings.HasChannel2 && ...
                        (h.hasCrossSection ||...
                        obj.ShowImageSettings.VolumeMenuValue ||...
                        obj.ShowImageSettings.ResultsMenuValue)
                    FullPart = 'PartTwo';
                elseif obj.ShowImageSettings.HasChannel2 && ...
                        ~(h.hasCrossSection ||...
                        obj.ShowImageSettings.VolumeMenuValue ||...
                        obj.ShowImageSettings.ResultsMenuValue)
                    FullPart = 'FullTwo';
                elseif ~obj.ShowImageSettings.HasChannel2 && ...
                        (h.hasCrossSection ||...
                        obj.ShowImageSettings.VolumeMenuValue ||...
                        obj.ShowImageSettings.ResultsMenuValue)
                    FullPart = 'PartOne';
                elseif ~obj.ShowImageSettings.HasChannel2 && ...
                        ~(h.hasCrossSection ||...
                        obj.ShowImageSettings.VolumeMenuValue ||...
                        obj.ShowImageSettings.ResultsMenuValue)
                    FullPart = 'FullOne';
                end
                obj.ShowImageSettings.DefaultChannel1Index = h.B(16).Value;
                obj.ShowImageSettings.DefaultChannel2Index = h.B(17).Value;
                obj.ShowImageSettings.DefaultChannel1SubIndex = h.B(2).Value;
                obj.ShowImageSettings.DefaultChannel2SubIndex = h.B(3).Value;
                h.hasChannel1 = true;
                draw_image(LeftRight,FullPart)
                if isequal(h.Channel{1},'none')
                    h.hasChannel1 = false;
                end
                if obj.ShowImageSettings.HasChannel2 && ~h.OnePass
                    h.OnePass = true;
                    draw_channel_2
                end
                h.OnePass = false;
            end
            
            function draw_channel_2(varargin)
                LeftRight = 'Right';
                if h.hasChannel1 && ...
                        (h.hasCrossSection ||...
                        obj.ShowImageSettings.VolumeMenuValue ||...
                        obj.ShowImageSettings.ResultsMenuValue)
                    FullPart = 'PartTwo';
                elseif h.hasChannel1 && ...
                        ~(h.hasCrossSection ||...
                        obj.ShowImageSettings.VolumeMenuValue ||...
                        obj.ShowImageSettings.ResultsMenuValue)
                    FullPart = 'FullTwo';
                elseif ~h.hasChannel1 && ...
                        (h.hasCrossSection ||...
                        obj.ShowImageSettings.VolumeMenuValue ||...
                        obj.ShowImageSettings.ResultsMenuValue)
                    FullPart = 'PartOne';
                elseif ~h.hasChannel1 && ...
                        ~(h.hasCrossSection ||...
                        obj.ShowImageSettings.VolumeMenuValue ||...
                        obj.ShowImageSettings.ResultsMenuValue)
                    FullPart = 'FullOne';
                end
                obj.ShowImageSettings.DefaultChannel1Index = h.B(16).Value;
                obj.ShowImageSettings.DefaultChannel2Index = h.B(17).Value;
                obj.ShowImageSettings.DefaultChannel1SubIndex = h.B(2).Value;
                obj.ShowImageSettings.DefaultChannel2SubIndex = h.B(3).Value;
                obj.ShowImageSettings.HasChannel2 = true;
                draw_image(LeftRight,FullPart)
                if isequal(h.Channel{2},'none')
                    obj.ShowImageSettings.HasChannel2 = false;
                end
                if h.hasChannel1 && ~h.OnePass
                    h.OnePass = true;
                    draw_channel_1
                end
                h.OnePass = false;
            end
            
            function moving_cross_section(src,evt)
                if ~get(h.B(1),'Value')
                    return
                end
                if isempty(h.MainLine)
                    return
                end
                if ~isvalid(h.MainLine)
                    h.MainLine = [];
                    h.MainLine = drawline('Position',h.MainProfilePosition,...
                        'Color',h.ColorMode(obj.ShowImageSettings.ColorIndex).Profile1,...
                        'Parent',h.ImAx(h.MainIndex),'LineWidth',obj.ShowImageSettings.ProfileLineWidth);
                    addlistener(h.MainLine,'MovingROI',@moving_cross_section);
                    addlistener(h.MainLine,'ROIMoved',@moving_cross_section);
                end
                delete(h.ImAx(3))
                h.ImAx(3) = [];
                h.MainProfilePosition = h.MainLine.Position;
                Pos1 = [h.MainLine.Position(1,1) h.MainLine.Position(1,2)];
                Pos2 = [h.MainLine.Position(2,1) h.MainLine.Position(2,2)];
                MainProfile = improfile(h.Image{h.MainIndex},[Pos1(1) Pos2(1)],[Pos1(2) Pos2(2)],'bicubic');
                Len = sqrt(((Pos1(1)-Pos2(1))*h.ScanSizeX(h.MainIndex)/h.NumPixelsY(h.MainIndex))^2 + ...
                ((Pos1(2)-Pos2(2))*h.ScanSizeY(h.MainIndex)/h.NumPixelsX(h.MainIndex))^2);
                Points = [0:1/(length(MainProfile)-1):1].*Len;
                [MainMultiplierY,UnitY,~] = AFMImage.parse_unit_scale(range(MainProfile),h.BaseUnit{h.MainIndex},1);
                [MultiplierX,UnitX,~] = AFMImage.parse_unit_scale(range(Points),'m',1);
                h.ImAx(3) = subplot(10,10,[71:78 81:88 91:98]);
                h.P = plot(Points.*MultiplierX,MainProfile.*MainMultiplierY);
                if obj.ShowImageSettings.BothCrossSections && (obj.ShowImageSettings.HasChannel2 && h.hasChannel1)
                    yyaxis left
                end
                hold on
                grid on
                CurrentAxHeight = round(h.Fig.Position(4)*h.ImAx(h.MainIndex).Position(4));
                h.ImAx(3).Color = h.ColorMode(obj.ShowImageSettings.ColorIndex).Background;
                h.ImAx(3).LineWidth = 1;
                h.ImAx(3).FontSize = round(obj.ShowImageSettings.ReferenceFontSize*(CurrentAxHeight/756));
                h.ImAx(3).XColor = h.ColorMode(obj.ShowImageSettings.ColorIndex).Text;
                h.ImAx(3).YColor = h.ColorMode(obj.ShowImageSettings.ColorIndex).Profile1;
                h.ImAx(3).GridColor = h.ColorMode(obj.ShowImageSettings.ColorIndex).Text;
                xlabel(sprintf('[%s]',UnitX))
                ylabel(sprintf('%s [%s]',h.Channel{h.MainIndex},UnitY))
                xlim([0 Points(end).*MultiplierX])
                h.P.LineWidth = 2;
                h.P.Color = h.ColorMode(obj.ShowImageSettings.ColorIndex).Profile1;
                if obj.ShowImageSettings.BothCrossSections && (obj.ShowImageSettings.HasChannel2 && h.hasChannel1)
                    if ~isempty(h.ChildLine)
                        if ~isstruct(h.ChildLine)
                            if ~isvalid(h.ChildLine)
                                h.ChildLine = [];
                            end
                        end
                    end
                    h.ChildLine.Visible = 'off';
                    
                    if h.B(31).Value && ...
                            h.Class{h.MainIndex}.OverlayGroup.hasOverlayGroup &&...
                            h.Class{h.ChildIndex}.OverlayGroup.hasOverlayGroup &&...
                            isequal(h.Class{h.MainIndex}.OverlayGroup.Names,h.Class{h.ChildIndex}.OverlayGroup.Names)
                        XDiff = h.Class{h.MainIndex}.Channel(1).OriginX - h.Class{h.ChildIndex}.Channel(1).OriginX;
                        SizePerPixelX = h.ScanSizeX(h.ChildIndex)./h.NumPixelsX(h.ChildIndex);
                        XDiff = XDiff/SizePerPixelX;
                        YDiff = h.Class{h.MainIndex}.Channel(1).OriginY - h.Class{h.ChildIndex}.Channel(1).OriginY;
                        SizePerPixelY = h.Class{h.ChildIndex}.Channel(1).ScanSizeY./h.Class{h.ChildIndex}.Channel(1).NumPixelsY;
                        YDiff = YDiff/SizePerPixelY;
                        AngleDiff = h.Class{h.MainIndex}.Channel(1).ScanAngle - h.Class{h.ChildIndex}.Channel(1).ScanAngle;
                        AngleDiff = deg2rad(-AngleDiff);
                        
                        % Resize the bigger Channel (in ScanSize-per-Pixel) so
                        % imagesizes correspond to ScanSizes. Do nothing, if they are
                        % exactly the same
                        SizePerPixel1 = h.ScanSizeX(h.MainIndex)/h.NumPixelsX(h.MainIndex);
                        SizePerPixel2 = h.ScanSizeX(h.ChildIndex)/h.NumPixelsX(h.ChildIndex);
                        if SizePerPixel1 ~= SizePerPixel2
                            ScaleMultiplier = SizePerPixel1/SizePerPixel2;
                        else
                            ScaleMultiplier = 1;
                        end
                        
                        XDiff = h.Class{h.MainIndex}.Channel(1).OriginX - h.Class{h.ChildIndex}.Channel(1).OriginX;
                        SizePerPixelX = h.ScanSizeX(h.ChildIndex)./h.NumPixelsX(h.ChildIndex);
                        XDiff = XDiff/SizePerPixelX;
                        YDiff = h.Class{h.MainIndex}.Channel(1).OriginY - h.Class{h.ChildIndex}.Channel(1).OriginY;
                        SizePerPixelY = h.Class{h.ChildIndex}.Channel(1).ScanSizeY./h.Class{h.ChildIndex}.Channel(1).NumPixelsY;
                        YDiff = YDiff/SizePerPixelY;
                        AngleDiff = h.Class{h.MainIndex}.Channel(1).ScanAngle - h.Class{h.ChildIndex}.Channel(1).ScanAngle;
                        AngleDiff = deg2rad(-AngleDiff);
                        
                        InitPos = h.MainLine.Position;
                        
                        CPos1 = ScaleMultiplier.*[InitPos(1,1) InitPos(1,2)];
                        CPos2 = ScaleMultiplier.*[InitPos(2,1) InitPos(2,2)];
                        
                        ImCenter = [h.NumPixelsX(h.ChildIndex)/2 h.Class{h.ChildIndex}.Channel(1).NumPixelsY/2];
                        
                        TempCP1 = CPos1 - ImCenter;
                        TempCP2 = CPos2 - ImCenter;
                        
                        RotationMatrix = [cos(AngleDiff) -sin(AngleDiff);sin(AngleDiff) cos(AngleDiff)];
                        
                        CPos1 = [RotationMatrix*TempCP1']' + ImCenter + [XDiff -YDiff];
                        CPos2 = [RotationMatrix*TempCP2']' + ImCenter + [XDiff -YDiff];
                        
                        InitPos = [CPos1(1) CPos1(2); CPos2(1) CPos2(2)];
                    else
                        InitPos = h.MainLine.Position;
                    end
                    if ~isvalid(h.ImAx(h.ChildIndex))
                        return
                    end
                    h.ChildLine = drawline('Position',InitPos,...
                        'Parent',h.ImAx(h.ChildIndex),'Color',h.ColorMode(obj.ShowImageSettings.ColorIndex).Profile2,...
                        'LineWidth',obj.ShowImageSettings.ProfileLineWidth);
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
                    h.ImAx(3).Color = h.ColorMode(obj.ShowImageSettings.ColorIndex).Background;
                    h.ImAx(3).LineWidth = 1;
                    h.ImAx(3).FontSize = round(obj.ShowImageSettings.ReferenceFontSize*(CurrentAxHeight/756));
                    h.ImAx(3).XColor = h.ColorMode(obj.ShowImageSettings.ColorIndex).Text;
                    h.ImAx(3).YColor = h.ColorMode(obj.ShowImageSettings.ColorIndex).Profile2;
                    h.ImAx(3).GridColor = h.ColorMode(obj.ShowImageSettings.ColorIndex).Text;
                    xlabel(sprintf('[%s]',UnitX))
                    ylabel(sprintf('%s [%s]',h.Channel{h.ChildIndex},UnitY))
                    xlim([0 ChildPoints(end).*MultiplierX])
                    h.CP.LineWidth = 2;
                    h.CP.Color = h.ColorMode(obj.ShowImageSettings.ColorIndex).Profile2;
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
%                         'FontSize',obj.ShowImageSettings.ReferenceFontSize);
                    
                end
                hold off
            end
            
            function checked_both_cross_sections(varargin)
                obj.ShowImageSettings.BothCrossSections = ~obj.ShowImageSettings.BothCrossSections;
                try
                    delete(h.ImAx(3));
                catch
                end
                draw_channel_1
                draw_channel_2
            end
            
            function get_and_draw_profile(varargin)
                if ~get(h.B(1),'Value')
                    return
                end
                if ~isempty(h.MainLine)
                    if ~isvalid(h.MainLine)
                        h.MainLine = [];
                    end
                end
                h.MainLine.Visible = 'off';
                if obj.ShowImageSettings.BothCrossSections && (obj.ShowImageSettings.HasChannel2 && h.hasChannel1)
                    if ~isempty(h.ChildLine)
                        if ~isvalid(h.ChildLine)
                            h.ChildLine = [];
                        end
                    end
                    h.ChildLine.Visible = 'off';
                end
                if isequal(varargin{1}.Parent,h.ImAx(1))
                    h.MainIndex = 1;
                    h.ChildIndex = 2;
                elseif isequal(varargin{1}.Parent,h.ImAx(2))
                    h.MainIndex = 2;
                    h.ChildIndex = 1;
                end
                try
                    delete(h.ImAx(3))
                catch
                end
                h.MainLine = drawline('Color',h.ColorMode(obj.ShowImageSettings.ColorIndex).Profile1,...
                    'Parent',varargin{1}.Parent,'LineWidth',obj.ShowImageSettings.ProfileLineWidth);
                addlistener(h.MainLine,'MovingROI',@moving_cross_section);
                addlistener(h.MainLine,'ROIMoved',@moving_cross_section);
                h.MainProfilePosition = h.MainLine.Position;
                h.MainProfileParent = h.MainLine.Parent;
                Pos1 = [h.MainLine.Position(1,1) h.MainLine.Position(1,2)];
                Pos2 = [h.MainLine.Position(2,1) h.MainLine.Position(2,2)];
                if norm(Pos1-Pos2)==0
                    get_and_draw_profile;
                    return
                end
                MainProfile = improfile(h.Image{h.MainIndex},[Pos1(1) Pos2(1)],[Pos1(2) Pos2(2)],'bicubic');
                Len = sqrt(((Pos1(1)-Pos2(1))*h.ScanSizeX(h.MainIndex)/h.NumPixelsY(h.MainIndex))^2 + ...
                ((Pos1(2)-Pos2(2))*h.ScanSizeY(h.MainIndex)/h.NumPixelsX(h.MainIndex))^2);
                Points = [0:1/(length(MainProfile)-1):1].*Len;
                [MainMultiplierY,UnitY,~] = AFMImage.parse_unit_scale(range(MainProfile),h.BaseUnit{h.MainIndex},1);
                [MultiplierX,UnitX,~] = AFMImage.parse_unit_scale(range(Points),'m',1);
                h.ImAx(3) = subplot(10,10,[71:78 81:88 91:98]);
                h.P = plot(Points.*MultiplierX,MainProfile.*MainMultiplierY);
                if obj.ShowImageSettings.BothCrossSections && (obj.ShowImageSettings.HasChannel2 && h.hasChannel1)
                    yyaxis left
                end
                hold on
                grid on
                CurrentAxHeight = round(h.Fig.Position(4)*h.ImAx(h.MainIndex).Position(4));
                h.ImAx(3).Color = h.ColorMode(obj.ShowImageSettings.ColorIndex).Background;
                h.ImAx(3).LineWidth = 1;
                h.ImAx(3).FontSize = round(obj.ShowImageSettings.ReferenceFontSize*(CurrentAxHeight/756));
                h.ImAx(3).XColor = h.ColorMode(obj.ShowImageSettings.ColorIndex).Text;
                h.ImAx(3).YColor = h.ColorMode(obj.ShowImageSettings.ColorIndex).Profile1;
                h.ImAx(3).GridColor = h.ColorMode(obj.ShowImageSettings.ColorIndex).Text;
                xlabel(sprintf('[%s]',UnitX))
                ylabel(sprintf('%s [%s]',h.Channel{h.MainIndex},UnitY))
                xlim([0 Points(end).*MultiplierX])
                h.P.LineWidth = 2;
                h.P.Color = h.ColorMode(obj.ShowImageSettings.ColorIndex).Profile1;
                if obj.ShowImageSettings.BothCrossSections && (obj.ShowImageSettings.HasChannel2 && h.hasChannel1)
                    
                    if h.B(31).Value && ...
                            h.Class{h.MainIndex}.OverlayGroup.hasOverlayGroup &&...
                            h.Class{h.ChildIndex}.OverlayGroup.hasOverlayGroup &&...
                            isequal(h.Class{h.MainIndex}.OverlayGroup.Names,h.Class{h.ChildIndex}.OverlayGroup.Names)
                        
                        
                        % Resize the bigger Channel (in ScanSize-per-Pixel) so
                        % imagesizes correspond to ScanSizes. Do nothing, if they are
                        % exactly the same
                        SizePerPixel1 = h.ScanSizeX(h.MainIndex)/h.NumPixelsX(h.MainIndex);
                        SizePerPixel2 = h.ScanSizeX(h.ChildIndex)/h.NumPixelsX(h.ChildIndex);
                        if SizePerPixel1 ~= SizePerPixel2
                            ScaleMultiplier = SizePerPixel1/SizePerPixel2;
                        else
                            ScaleMultiplier = 1;
                        end
                        
                        XDiff = h.Class{h.MainIndex}.Channel(1).OriginX - h.Class{h.ChildIndex}.Channel(1).OriginX;
                        SizePerPixelX = h.ScanSizeX(h.ChildIndex)./h.NumPixelsX(h.ChildIndex);
                        XDiff = XDiff/SizePerPixelX;
                        YDiff = h.Class{h.MainIndex}.Channel(1).OriginY - h.Class{h.ChildIndex}.Channel(1).OriginY;
                        SizePerPixelY = h.Class{h.ChildIndex}.Channel(1).ScanSizeY./h.Class{h.ChildIndex}.Channel(1).NumPixelsY;
                        YDiff = YDiff/SizePerPixelY;
                        AngleDiff = h.Class{h.MainIndex}.Channel(1).ScanAngle - h.Class{h.ChildIndex}.Channel(1).ScanAngle;
                        AngleDiff = deg2rad(-AngleDiff);
                        
                        InitPos = h.MainLine.Position;
                        
                        CPos1 = ScaleMultiplier.*[InitPos(1,1) InitPos(1,2)];
                        CPos2 = ScaleMultiplier.*[InitPos(2,1) InitPos(2,2)];
                        
                        ImCenter = [h.NumPixelsX(h.ChildIndex)/2 h.Class{h.ChildIndex}.Channel(1).NumPixelsY/2];
                        
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
                        'Parent',h.ImAx(h.ChildIndex),'Color',h.ColorMode(obj.ShowImageSettings.ColorIndex).Profile2,...
                        'LineWidth',obj.ShowImageSettings.ProfileLineWidth);
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
                    h.ImAx(3).Color = h.ColorMode(obj.ShowImageSettings.ColorIndex).Background;
                    h.ImAx(3).LineWidth = 1;
                    h.ImAx(3).FontSize = round(obj.ShowImageSettings.ReferenceFontSize*(CurrentAxHeight/756));
                    h.ImAx(3).XColor = h.ColorMode(obj.ShowImageSettings.ColorIndex).Text;
                    h.ImAx(3).YColor = h.ColorMode(obj.ShowImageSettings.ColorIndex).Profile2;
                    h.ImAx(3).GridColor = h.ColorMode(obj.ShowImageSettings.ColorIndex).Text;
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
                    h.CP.Color = h.ColorMode(obj.ShowImageSettings.ColorIndex).Profile2;
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
                    CurrentZoomX = h.ImAx(Index).XLim;
                    CurrentZoomY = h.ImAx(Index).YLim;
                    delete(h.ImAx(Index));
                    delete(h.I(Index));
                catch
                end
                if isequal(FullPart,'FullOne')
                    h.ImAx(Index) = axes(h.Fig,'Position',[0.1 0.1 .7 .8]);
                elseif isequal(FullPart,'FullTwo')
                    if isequal(LeftRight,'Left')
                    h.ImAx(Index) = axes(h.Fig,'Position',[0.12 0.1 .3 .8]);
                    else
                    h.ImAx(Index) = axes(h.Fig,'Position',[.47 0.1 .3 .8]);
                    end
                elseif isequal(FullPart,'PartOne')
                    h.ImAx(Index) = axes(h.Fig,'Position',[0.1 .35 .7 .6]);
                elseif isequal(FullPart,'PartTwo')
                    if isequal(LeftRight,'Left')
                    h.ImAx(Index) = axes(h.Fig,'Position',[0.12 .35 .3 .6]);
                    else
                    h.ImAx(Index) = axes(h.Fig,'Position',[0.47 .35 .3 .6]);
                    end
                end
                
                
                if obj.ShowImageSettings.LockChannels && Index==1
                    h.CurChannel1Idx = h.B(16).Value;
                    h.B(17).Value = mod(h.CurChannel1Idx + obj.ShowImageSettings.RelativeChannelIndex,h.NumClasses+1);
                    h.CurChannel2Idx = h.B(17).Value;
                    CurIndex = h.CurChannel1Idx;
                elseif obj.ShowImageSettings.LockChannels && Index==2
                    h.CurChannel2Idx = h.B(17).Value;
                    h.B(16).Value = mod(h.CurChannel2Idx - obj.ShowImageSettings.RelativeChannelIndex,h.NumClasses+1);
                    h.CurChannel1Idx = h.B(16).Value;
                    CurIndex = h.CurChannel2Idx;
                else
                    CurIndex = h.B(15+Index).Value;
                    h.CurChannel2Idx = h.B(17).Value;
                    h.CurChannel1Idx = h.B(16).Value;
                end
                
                obj.ShowImageSettings.RelativeChannelIndex = h.CurChannel2Idx - h.CurChannel1Idx;
                
                h.Class{Index} = obj.get_class_instance(h.ClassIndex(CurIndex,:));
                CurrentChannelName = h.B(1+Index).String{h.B(1+Index).Value};
                PopUp = h.Class{Index}.string_of_existing();
                set(h.B(1+Index),'String',PopUp)
                
                IdxOfSameChannel = find(strcmp(PopUp,CurrentChannelName));
                
                if isempty(IdxOfSameChannel)
                    set(h.B(1+Index),'Value',2);
                elseif sum(IdxOfSameChannel == h.B(1+Index).Value)
                    set(h.B(1+Index),'Value',IdxOfSameChannel(find(IdxOfSameChannel == h.B(1+Index).Value)));
                else
                    set(h.B(1+Index),'Value',IdxOfSameChannel(1));
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
                        obj.ShowImageSettings.HasChannel2 = 0;
                    end
                    return
                else
                    [Channel,ChannelIndex] = h.Class{Index}.get_channel(h.Channel{Index});
                    if obj.ShowImageSettings.IsUpscaled
                        Channel.Image = fillmissing(Channel.Image,'linear','EndValues','nearest');
                        Channel = AFMImage.resize_channel(Channel,1920,false);
                    end
                    h.Image{Index} = fillmissing(Channel.Image,'linear','EndValues','nearest');
                    h.BaseUnit{Index} = Channel.Unit;
                    h.ScanSizeX(Index) = Channel.ScanSizeX;
                    h.ScanSizeY(Index) = Channel.ScanSizeY;
                    if isfield(h,'NumPixelsX') && isfield(h,'NumPixelsY') &&...
                            length(h.NumPixelsX) == 2 &&...
                            length(h.NumPixelsY) == 2 &&...
                            ((h.NumPixelsX(Index) ~= Channel.NumPixelsX) ||...
                            (h.NumPixelsY(Index) ~= Channel.NumPixelsY))
                        h.ResetOnNextZoom = true;
                    end
                    h.NumPixelsX(Index) = Channel.NumPixelsX;
                    h.NumPixelsY(Index) = Channel.NumPixelsY;
                    ColorPattern = h.Class{Index}.CMap;
                end
                
                if h.B(21).Value && obj.ShowImageSettings.HasChannel2 && h.hasChannel1 && isequal(h.BaseUnit{1},h.BaseUnit{2})
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
                
                if ~isempty(h.SubChannelName{Index}) &&...
                        ~isequal(h.SubChannelName{Index},Channel.Name)
                    h.B(MinIndex).String = [];
                    h.B(MaxIndex).String = [];
                    h.OldMultiplie{Index} = [];
                end
                h.SubChannelName{Index} = Channel.Name;
                
                if ~isempty(get(h.B(MinIndex),'String')) && ~isempty(get(h.B(MaxIndex),'String'))
                    if ~isempty(h.OldMultiplier{Index}) && ...
                            h.OldMultiplier{Index} ~= h.Multiplier{Index}
                        MinValue = str2double(get(h.B(MinIndex),'String'))...
                            *(h.Multiplier{Index}/h.OldMultiplier{Index})/h.Multiplier{Index};
                        MaxValue = str2double(get(h.B(MaxIndex),'String'))...
                            *(h.Multiplier{Index}/h.OldMultiplier{Index})/h.Multiplier{Index};
                    else
                        MinValue = str2double(get(h.B(MinIndex),'String'))/h.Multiplier{Index};
                        MaxValue = str2double(get(h.B(MaxIndex),'String'))/h.Multiplier{Index};
                    end
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
                if h.MenuTabs(1).Value
                    h.I(Index).ButtonDownFcn = @get_and_draw_profile;
                elseif h.MenuTabs(2).Value
                    h.I(Index).ButtonDownFcn = @get_and_draw_volume_information;
                elseif h.MenuTabs(3).Value
                    h.I(Index).ButtonDownFcn = {@select_clicked_roi,Index};
                end
                hold on
                CurrentAxHeight = round(h.Fig.Position(4)*h.ImAx(Index).Position(4));
                CurrentAxWidth = round(h.Fig.Position(3)*h.ImAx(Index).Position(3));
                AFMImage.draw_scalebar_into_current_image(Channel.NumPixelsX,Channel.NumPixelsY,Channel.ScanSizeX,BarToImageRatio,CurrentAxHeight,CurrentAxWidth);
                c = colorbar('northoutside');
                c.FontSize = round(obj.ShowImageSettings.ReferenceFontSize*(CurrentAxHeight/756));
                c.Color = h.ColorMode(obj.ShowImageSettings.ColorIndex).Text;
                c.Label.String = sprintf('%s [%s]',h.Channel{Index},h.Unit{Index});
                c.Label.FontSize = round(obj.ShowImageSettings.ReferenceFontSize*(CurrentAxHeight/756));
                c.Label.Color = h.ColorMode(obj.ShowImageSettings.ColorIndex).Text;
                
                if isequal(h.CurrentClassName{Index},h.Class{Index}.Name) && ~h.ResetOnNextZoom
                    try
                        % Set Zoom region to previous
                        zoom reset
                        h.ImAx(Index).XLim = CurrentZoomX;
                        h.ImAx(Index).YLim = CurrentZoomY;
                    catch
                    end
                else
                    h.CurrentClassName{Index} = h.Class{Index}.Name;
                end
                
                set(h.B(MinIndex+4),'String',{'Min',sprintf('[%s]',h.Unit{Index})});
                set(h.B(MaxIndex+4),'String',{'Max',sprintf('[%s]',h.Unit{Index})});
                if ~isempty(get(h.B(MinIndex),'String')) && ~isempty(get(h.B(MaxIndex),'String'))
                    h.OldMultiplier{Index} = h.Multiplier{Index};
                    h.B(MinIndex).String = num2str(MinValue.*h.Multiplier{Index});
                    h.B(MaxIndex).String = num2str(MaxValue.*h.Multiplier{Index});
                else
                    h.OldMultiplier{Index} = [];
                end
                
                moving_cross_section
                draw_volume_information
                draw_existing_segments(Index)
%                 if obj.ShowImageSettings.ResultsMenuValue
%                     if h.hasChannel1
%                         draw_existing_segments(1)
%                     end
%                     if obj.ShowImageSettings.HasChannel2
%                         draw_existing_segments(2)
%                     end
%                 end
                try
                    if (h.ScanSizeX(1) == h.ScanSizeX(2)) &&...
                            (h.ScanSizeY(1) == h.ScanSizeY(2))
                        h.B(18).BackgroundColor = [.302 .6902 .302];
                    else
                        h.B(18).BackgroundColor = [.8392 .2706 .2706];
                    end
                catch
                end
                check_for_overlay_group
            end
            
            function ColorMode = set_default_color_options()
                
                ColorMode(1).Background = 'k';
            	ColorMode(1).Profile1 = [219 21 223]./255; %[189 0 96]./255; % 'b';
            	ColorMode(1).Profile2 = 'c';
            	ColorMode(1).Text = 'w';
            	ColorMode(1).Backdrop = [.05 .05 .05];
            	ColorMode(1).Buttons = [.1 .1 .1];
            	ColorMode(1).ButtonText = 'w';
                ColorMode(1).SpecialHighlight = 'w';
                ColorMode(1).Disabled = [.65 .65 .65];
            
            
            
            	ColorMode(2).Background = 'w';
            	ColorMode(2).Profile1 = [219 21 223]./255; %[94 170 170]./255; % [189 0 96]./255; %'b';
            	ColorMode(2).Profile2 = [0 181 26]./255; % alternatives %[80 200 204]./255;%[0,0.870588235294118,0.407843137254902];
            	ColorMode(2).Text = 'k';
            	ColorMode(2).Backdrop = [.95 .95 .95];
            	ColorMode(2).Buttons = [.9 .9 .9];
            	ColorMode(2).ButtonText = 'k';
                ColorMode(2).SpecialHighlight = 'k';
                ColorMode(2).Disabled = [.65 .65 .65];
            
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
                
                filter = {obj.ShowImageSettings.DefaultSaveType;'*.eps';'*.emf';'*.png';'*.tif';'*.jpg'};
                Name = split(h.Class{1}.Name,'.');
                DefName = Name{1};
                DefFull = fullfile(obj.ShowImageSettings.DefaultSavePath,DefName);
                [file, path] = uiputfile(filter,'Select Format, Name and Location of your figure',DefFull);
                obj.ShowImageSettings.DefaultSavePath = path;
                FullFile = fullfile(path,file);
                Split = split(FullFile,'.');
                obj.ShowImageSettings.DefaultSaveType = strcat('*.',Split{end});
                if isequal(Split{end},'eps') || isequal(Split{end},'emf')
                    exportgraphics(h.Fig,FullFile,'ContentType','vector','Resolution',300,'BackgroundColor','current')
                else
                    exportgraphics(h.Fig,FullFile,'Resolution',300,'BackgroundColor','current')
                end
                    
            end
            
            function changed_color(varargin)
                
                if ~h.B(7).Value
                    obj.ShowImageSettings.ColorIndex = 1;
                else
                    obj.ShowImageSettings.ColorIndex = 2;
                end
                
                h.Fig.Color = h.ColorMode(obj.ShowImageSettings.ColorIndex).Background;
                h.Backdrop.BackgroundColor = h.ColorMode(obj.ShowImageSettings.ColorIndex).Backdrop;
                for i=1:length(h.B)
                    h.B(i).ForegroundColor = h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText;
                    h.B(i).BackgroundColor = h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons;
                end
                for i=1:length(h.FVM)
                    h.FVM(i).ForegroundColor = h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText;
                    h.FVM(i).BackgroundColor = h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons;
                end
                for i=1:length(h.Res)
                    h.Res(i).ForegroundColor = h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText;
                    h.Res(i).BackgroundColor = h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons;
                end
                h.Res(8).BackgroundColor = h.ColorMode(obj.ShowImageSettings.ColorIndex).Backdrop;
                for i=1:length(h.MenuTabs)
                    h.MenuTabs(i).ForegroundColor = h.ColorMode(obj.ShowImageSettings.ColorIndex).ButtonText;
                    h.MenuTabs(i).BackgroundColor = h.ColorMode(obj.ShowImageSettings.ColorIndex).Buttons;
                end
                
                draw_channel_1
                draw_channel_2
                
            end
            
            function upscale_images(varargin)
                
                obj.ShowImageSettings.IsUpscaled = h.B(19).Value;
                h.ResetOnNextZoom = true;
                
                draw_channel_1
                draw_channel_2
                
                h.ResetOnNextZoom = false;
                
            end
            
            function lock_channels(varargin)
                
                obj.ShowImageSettings.LockChannels = h.B(20).Value;
                
                obj.ShowImageSettings.RelativeChannelIndex = h.CurChannel2Idx - h.CurChannel1Idx;
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
            
            function changed_use_overlay(varargin)
                
                draw_channel_1
                draw_channel_2
                
            end
            
            function check_for_overlay_group(varargin)
                
                
                if obj.ShowImageSettings.HasChannel2 &&...
                            ~isempty(h.Class{1}.OverlayGroup) &&...
                            ~isempty(h.Class{2}.OverlayGroup) &&...
                            h.Class{1}.OverlayGroup.hasOverlayGroup &&...
                            h.Class{2}.OverlayGroup.hasOverlayGroup &&...
                            isequal(h.Class{1}.OverlayGroup.Names,h.Class{2}.OverlayGroup.Names)
                    h.B(31).BackgroundColor = [.302 .6902 .302];
                else
                    h.B(31).BackgroundColor = [.8392 .2706 .2706];
                end
                
            end
            
            function switch_to_main_menu(varargin)
                
                h.MenuTabs(1).Value = 1;
                obj.ShowImageSettings.MainMenuValue = 1;
                h.MenuTabs(2).Value = 0;
                obj.ShowImageSettings.VolumeMenuValue = 0;
                h.MenuTabs(3).Value = 0;
                obj.ShowImageSettings.ResultsMenuValue = 0;
                
                update_interface
                
                try
                    delete(h.ImAx(3));
                catch
                end
                
                draw_channel_1
                draw_channel_2
            end
            
            function switch_to_volume_menu(varargin)
                
                h.MenuTabs(1).Value = 0;
                obj.ShowImageSettings.MainMenuValue = 0;
                h.MenuTabs(2).Value = 1;
                obj.ShowImageSettings.VolumeMenuValue = 1;
                h.MenuTabs(3).Value = 0;
                obj.ShowImageSettings.ResultsMenuValue = 0;
                h.hasCrossSection = 0;
                h.B(1).Value = 0;
                
                update_interface
                
                try
                    delete(h.ImAx(3));
                catch
                end
                
                
                
                draw_channel_1
                draw_channel_2
            end
            
            function switch_to_results_menu(varargin)
                
                h.MenuTabs(1).Value = 0;
                obj.ShowImageSettings.MainMenuValue = 0;
                h.MenuTabs(2).Value = 0;
                obj.ShowImageSettings.VolumeMenuValue = 0;
                h.MenuTabs(3).Value = 1;
                obj.ShowImageSettings.ResultsMenuValue = 1;
                h.hasCrossSection = 0;
                h.B(1).Value = 0;
                
                update_interface
                
                try
                    delete(h.ImAx(3));
                catch
                end
                
                draw_channel_1
                draw_channel_2
            end
            
            function update_interface(varargin)
                
                for i=1:length(h.ToggledMainMenuIndizes)
                    h.B(h.ToggledMainMenuIndizes(i)).Visible =...
                        obj.ShowImageSettings.MainMenuValue;
                end
                for i=1:length(h.ToggledVolumeMenuIndizes)
                    h.FVM(h.ToggledVolumeMenuIndizes(i)).Visible =...
                        obj.ShowImageSettings.VolumeMenuValue;
                end
                for i=1:length(h.ToggledResultMenuIndizes)
                    h.Res(h.ToggledResultMenuIndizes(i)).Visible =...
                        obj.ShowImageSettings.ResultsMenuValue;
                end
            end
            
            function get_and_draw_volume_information(varargin)
                
                if ~h.MenuTabs(2).Value
                    return
                end
                
                if isequal(varargin{1}.Parent,h.ImAx(1))
                    TempIndex = 1;
                elseif isequal(varargin{1}.Parent,h.ImAx(2))
                    TempIndex = 2;
                end
                h.VolumeStruct.ClassType = class(h.Class{TempIndex});
                if ~isequal(h.VolumeStruct.ClassType,'ForceMap')
                    warning('This Class Instance does not have Volume information!')
                    return
                end
                h.VolumeStruct.Index = TempIndex;
                
                if ~isempty(h.VolumeStruct.Point)
                    delete(h.VolumeStruct.Point);
                end
                
                h.VolumeStruct.Point = drawpoint('Color',...
                    h.ColorMode(obj.ShowImageSettings.ColorIndex).Profile1,...
                    'Parent',varargin{1}.Parent,...
                    'LineWidth',obj.ShowImageSettings.ProfileLineWidth);
                
                addlistener(h.VolumeStruct.Point,'MovingROI',@draw_volume_information);
                addlistener(h.VolumeStruct.Point,'ROIMoved',@draw_volume_information);
                
                draw_volume_information
                
            end
            
            function draw_volume_information(varargin)
                
                if ~h.MenuTabs(2).Value ||...
                        isempty(h.VolumeStruct.Point) ||...
                        (~isvalid(h.VolumeStruct.Point) && isempty(h.VolumeStruct.MapIndex))
                    return
                end
                h.VolumeStruct.ClassType = class(h.Class{h.VolumeStruct.Index});
                if ~isequal(h.VolumeStruct.ClassType,'ForceMap')
                    warning('This Class Instance does not have Volume information!')
                    try
                        delete(h.ImAx(3));
                    catch
                    end
                    return
                end
                
                
                IsLoaded = h.Class{h.VolumeStruct.Index}.check_if_zipped_file_is_loaded(1);
                
                if ~IsLoaded
                    h.VolumeStruct.Point.Visible = 'off';
                    h.VolumeStruct.Point = [];
                    return
                end
                
                % Snap Point to center of chosen Pixel and
                % draw new point on old coordinates, if image has been
                % deleted at some point
                if ~isempty(h.VolumeStruct.Point) &&...
                        ~isvalid(h.VolumeStruct.Point)
                    h.VolumeStruct.Point = drawpoint('Color',...
                        h.ColorMode(obj.ShowImageSettings.ColorIndex).Profile1,...
                        'Parent',h.ImAx(h.VolumeStruct.Index),...
                        'LineWidth',obj.ShowImageSettings.ProfileLineWidth,...
                        'Position',h.VolumeStruct.MapIndex);
                    addlistener(h.VolumeStruct.Point,'MovingROI',@draw_volume_information);
                    addlistener(h.VolumeStruct.Point,'ROIMoved',@draw_volume_information);
                end
                h.VolumeStruct.Point.Position = ...
                    round(h.VolumeStruct.Point.Position);
                h.VolumeStruct.Point.Position(h.VolumeStruct.Point.Position <= 0) = 1;
                if h.VolumeStruct.Point.Position(1) >= h.NumPixelsY(h.VolumeStruct.Index)
                    h.VolumeStruct.Point.Position(1) = h.NumPixelsY(h.VolumeStruct.Index);
                end
                if h.VolumeStruct.Point.Position(2) >= h.NumPixelsX(h.VolumeStruct.Index)
                    h.VolumeStruct.Point.Position(2) = h.NumPixelsX(h.VolumeStruct.Index);
                end
                h.VolumeStruct.MapIndex = h.VolumeStruct.Point.Position;
                h.VolumeStruct.ListIndex = ...
                    h.Class{h.VolumeStruct.Index}.Map2List(...
                    h.VolumeStruct.MapIndex(2),h.VolumeStruct.MapIndex(1));
                
                
                if h.FVM(7).Value
                    PointValue = h.Image{h.VolumeStruct.Index}...
                        (h.VolumeStruct.Point.Position(2),h.VolumeStruct.Point.Position(1));
                    [ScalingFactor,Unit] = AFMImage.parse_unit_scale(PointValue,h.BaseUnit{h.VolumeStruct.Index},1);
                    Scaled = ScalingFactor*PointValue;
                    PointLabel = [num2str(Scaled) ' ' Unit];
                    h.VolumeStruct.Point.Label = PointLabel;
                end
                
                try
                    delete(h.ImAx(3));
                catch
                end
                
                h.ImAx(3) = subplot(10,10,[71:78 81:88 91:98]);
                
                % get and process requested data according to chosen
                % options
                [RawApp,RawHHApp] = ...
                    h.Class{h.VolumeStruct.Index}.get_force_curve_data(...
                    h.VolumeStruct.ListIndex,'AppRetSwitch',0,...
                    'BaselineCorrection',0,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                [RawRet,RawHHRet] = ...
                    h.Class{h.VolumeStruct.Index}.get_force_curve_data(...
                    h.VolumeStruct.ListIndex,'AppRetSwitch',1,...
                    'BaselineCorrection',0,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                if obj.ShowImageSettings.UseCorrectedSensitivity &&...
                        ~isempty(h.Class{h.VolumeStruct.Index}.RefSlopeCorrectedSensitivity) &&...
                        ~isequal(obj.ForceMapAnalysisOptions.SensitivityCorrection.SensitivityCorrectionMethod,'Adaptive')
                    CorrectedSens = ...
                        h.Class{h.VolumeStruct.Index}.RefSlopeCorrectedSensitivity;
                elseif obj.ShowImageSettings.UseCorrectedSensitivity &&...
                        isequal(obj.ForceMapAnalysisOptions.SensitivityCorrection.SensitivityCorrectionMethod,'Adaptive')
                    CorrectedSens = h.Class{...
                        h.VolumeStruct.Index}.AdaptiveSensitivity(h.VolumeStruct.ListIndex);
                else
                    CorrectedSens = h.Class{h.VolumeStruct.Index}.Sensitivity;
                end
                Sens = h.Class{h.VolumeStruct.Index}.Sensitivity;
                SpringConstant = h.Class{h.VolumeStruct.Index}.SpringConstant;
                if obj.ShowImageSettings.BaselineCorrected &&...
                        h.Class{h.VolumeStruct.Index}.BaseAndTiltFlag
                    BaseParams = h.Class{h.VolumeStruct.Index}.Basefit{...
                        h.VolumeStruct.ListIndex};
                else
                    BaseParams = [0 0];
                end
                if obj.ShowImageSettings.TipHeight
                    TipHeightSwitch = 1;
                else
                    TipHeightSwitch = 0;
                end
                if obj.ShowImageSettings.ContactPointShifted
                    CP = h.Class{h.VolumeStruct.Index}.CP(h.VolumeStruct.ListIndex,:);
                else
                    CP = [0 0];
                end
                
                % Calculate the data
                App = ((RawApp - polyval(BaseParams,RawHHApp)) - CP(2)).*(CorrectedSens/Sens);
                Ret = ((RawRet - polyval(BaseParams,RawHHRet)) - CP(2)).*(CorrectedSens/Sens);
                HHApp = (RawHHApp - CP(1)) - TipHeightSwitch.*(App/SpringConstant);
                HHRet = (RawHHRet - CP(1)) - TipHeightSwitch.*(Ret/SpringConstant);
                
                if obj.ShowImageSettings.PlotTime
                    TimeApp = (RawHHApp - RawHHApp(1))./h.Class{h.VolumeStruct.Index}.ExtendVelocity;
                    TimeRet = (-RawHHRet + RawHHRet(1))./h.Class{h.VolumeStruct.Index}.RetractVelocity...
                        + TimeApp(end) + h.Class{h.VolumeStruct.Index}.HoldingTime;
                    
                    X1 = TimeApp;
                    X2 = TimeRet;
                else
                    X1 = HHApp;
                    X2 = HHRet;
                end
                
                if h.FVM(8).Value == 1
                    Y1 = App;
                    Y2 = Ret;
                elseif h.FVM(8).Value == 2
                    Y1 = App./SpringConstant.*(CorrectedSens/Sens);
                    Y2 = Ret./SpringConstant.*(CorrectedSens/Sens);
                elseif h.FVM(8).Value == 3
                    Y1 = App./(SpringConstant.*Sens);
                    Y2 = Ret./(SpringConstant.*Sens);
                end
                
                % Add HertzFit, if requested
                if obj.ShowImageSettings.ShowHertzFit && ...
                        obj.ShowImageSettings.TipHeight
                    FitFunction = fittype(h.Class{h.VolumeStruct.Index}.HertzFitType);
                    FitCoeffValues = ...
                        h.Class{h.VolumeStruct.Index}.HertzFitValues{h.VolumeStruct.ListIndex};
                    if length(FitCoeffValues) == 1
                        Fit = cfit(FitFunction,FitCoeffValues);
                    elseif length(FitCoeffValues) == 2
                        Fit = cfit(FitFunction,FitCoeffValues(1),FitCoeffValues(2));
                    end
                    if ~obj.ShowImageSettings.ContactPointShifted
                        CP = h.Class{h.VolumeStruct.Index}.CP(h.VolumeStruct.ListIndex,:);
                        X3 = HHApp(HHApp - CP(1) >= 0) - CP(1);
                        Y3 = feval(Fit,X3);
                        X3 = X3 + CP(1);
                        Y3 = Y3 + CP(2);
                    else
                        X3 = HHApp(HHApp >= 0);
                        Y3 = feval(Fit,X3);
                    end
                else
                    X3 = [];
                    Y3 = [];
                end
                
                % scale the data and determine units and prefixes
                if h.FVM(6).Value
                    BaseUnitX = 's';
                else
                    BaseUnitX = 'm';
                end
                BaseUnitY = h.FVM(8).String{h.FVM(8).Value};
                
                RangeX = max(range(X1),range(X2));
                RangeY = max(range(Y1),range(Y2));
                
                [h.VolumeStruct.MultiplierX,UnitX,~] = ...
                    AFMImage.parse_unit_scale(RangeX,BaseUnitX,10);
                [h.VolumeStruct.MultiplierY,UnitY,~] = ...
                    AFMImage.parse_unit_scale(RangeY,BaseUnitY,10);
                X1 = X1.*h.VolumeStruct.MultiplierX;
                X2 = X2.*h.VolumeStruct.MultiplierX;
                X3 = X3.*h.VolumeStruct.MultiplierX;
                Y1 = Y1.*h.VolumeStruct.MultiplierY;
                Y2 = Y2.*h.VolumeStruct.MultiplierY;
                Y3 = Y3.*h.VolumeStruct.MultiplierY;
                
                if obj.ShowImageSettings.PlotTime
                    XLabel = ['Time [' UnitX ']'];
                elseif obj.ShowImageSettings.TipHeight
                    XLabel = ['Tip Height [' UnitX ']'];
                else
                    XLabel = ['Head Height [' UnitX ']'];
                end
                YLabel = ['vDeflection [' UnitY ']'];
                
                % TODO needs a dedicated button
%                 % Sort Height
%                 X1 = sort(X1);
%                 if obj.ShowImageSettings.PlotTime
%                     X2 = sort(X2);
%                 else
%                     X2 = sort(X2,'descend');
%                 end
                
                % plot data
                h.VolumeStruct.Plot = plot(X1,Y1,X2,Y2,X3,Y3);
                hold on
                grid on
                CurrentAxHeight = ...
                    round(h.Fig.Position(4)*h.ImAx(1).Position(4));
                h.ImAx(3).Color = ...
                    h.ColorMode(obj.ShowImageSettings.ColorIndex).Background;
                h.ImAx(3).LineWidth = 1;
                h.ImAx(3).FontSize = ...
                    round(obj.ShowImageSettings.ReferenceFontSize*(CurrentAxHeight/756));
                h.ImAx(3).XColor = ...
                    h.ColorMode(obj.ShowImageSettings.ColorIndex).Text;
                h.ImAx(3).YColor = ...
                    h.ColorMode(obj.ShowImageSettings.ColorIndex).Text;
                h.ImAx(3).GridColor = ...
                    h.ColorMode(obj.ShowImageSettings.ColorIndex).Text;
                 xlabel(XLabel)
                 ylabel(YLabel)
                h.VolumeStruct.Plot(1).LineWidth = 2;
                h.VolumeStruct.Plot(2).LineWidth = 2;
                h.VolumeStruct.Plot(1).Color = ...
                    h.ColorMode(obj.ShowImageSettings.ColorIndex).Profile1;
                h.VolumeStruct.Plot(2).Color = ...
                    h.ColorMode(obj.ShowImageSettings.ColorIndex).Profile2;
                if obj.ShowImageSettings.ShowHertzFit && ...
                        obj.ShowImageSettings.TipHeight
                    h.VolumeStruct.Plot(3).LineWidth = 2;
                    h.VolumeStruct.Plot(3).Color = 'r';
                    Legend = legend({'Approach','Retract','Hertz Fit'});
                else
                    Legend = legend({'Approach','Retract'});
                end
                Legend.TextColor = ...
                    h.ColorMode(obj.ShowImageSettings.ColorIndex).Text;
                Legend.Color = ...
                    h.ColorMode(obj.ShowImageSettings.ColorIndex).Backdrop;
                Legend.Location = 'northwest';
            end
            
            function changed_baseline_corrected(varargin)
                
                obj.ShowImageSettings.BaselineCorrected = ...
                    h.FVM(1).Value;
                
                draw_volume_information
            end
            
            function changed_tip_height(varargin)
                
                obj.ShowImageSettings.TipHeight = ...
                    h.FVM(2).Value;
                
                draw_volume_information
            end
            
            function changed_contact_point_shifted(varargin)
                
                obj.ShowImageSettings.ContactPointShifted = ...
                    h.FVM(3).Value;
                
                draw_volume_information
            end
            
            function changed_show_hertz_fit(varargin)
                
                obj.ShowImageSettings.ShowHertzFit = ...
                    h.FVM(4).Value;
                
                draw_volume_information
            end
            
            function changed_use_corrected_sensitivity(varargin)
                
                obj.ShowImageSettings.UseCorrectedSensitivity = ...
                    h.FVM(5).Value;
                
                draw_volume_information
            end
            
            function changed_plot_time(varargin)
                
                obj.ShowImageSettings.PlotTime = ...
                    h.FVM(6).Value;
                
                draw_volume_information
            end
            
            function changed_vdeflection_unit(varargin)
                
                obj.ShowImageSettings.vDeflectionUnitIndex = ...
                    h.FVM(7).Value;
                
                draw_volume_information
            end
            
            function changed_extended_information(varargin)
                
                obj.ShowImageSettings.ExtendedInformation = ...
                    h.FVM(8).Value;
                
                draw_volume_information
            end
            
            function changed_gramm_x(varargin)
                
                obj.ShowImageSettings.GrammX = ...
                    h.Res(17).Value;
                
            end
            
            function changed_gramm_y(varargin)
                
                obj.ShowImageSettings.GrammY = ...
                    h.Res(19).Value;
                
            end
            
            function changed_gramm_group(varargin)
                
                obj.ShowImageSettings.GrammGroup = ...
                    h.Res(21).Value;
                
            end
            
            function changed_gramm_facets(varargin)
                
                obj.ShowImageSettings.GrammFacets = ...
                    h.Res(23).Value;
                
            end
            
            function draw_results(varargin)
                
                if ~obj.ShowImageSettings.ResultsMenuValue
                    return
                end
                
%                 draw_existing_segments(1)
%                 draw_existing_segments(2)
                draw_result_plots
                
            end
            
            function draw_existing_segments(Index,varargin)
                
                if ~obj.ShowImageSettings.ResultsMenuValue
                    return
                end
                if isempty(h.Class{Index}.Segment)
                    return
                end
                if isempty(h.Class{Index}.Segment(1).Name)
                    return
                end
                
                Names = {h.Class{Index}.Segment.Name};
                UniqueNames = unique(Names);
                NumSegments = length(UniqueNames);
                NumSubSegments = length(Names);
                
                
                
                % This guy figured out one simple trick to get all
                % the colors. Painters hate him
                FakePlots = zeros(2,NumSegments*2);
                FakeAx = plot(FakePlots);
                SegmentColors = get(FakeAx,'Color');
                delete(FakeAx)
                clear FakeAx FakePlots
                
                
                h.ResultStruct.ChosenSubsetIndizes{Index} = [];
                for i=1:NumSubSegments
                    if (obj.ShowImageSettings.UseSnapped &&...
                            sum(strfind(Names{i},'Snapped'))) ||...
                            (~obj.ShowImageSettings.UseSnapped &&...
                            ~sum(strfind(Names{i},'Snapped')))
                        h.ResultStruct.ChosenSubsetIndizes{Index}(end+1) = i;
                    else
                        continue
                    end
                    SegmentIndex = find(strcmp(UniqueNames,h.Class{Index}.Segment(i).Name));
                    if sum(h.ResultStruct.DisabledSegmentIndizes{Index} == i)
                        h.ROIObjects{Index}{i} = drawpolyline('Position',h.Class{Index}.Segment(i).ROIObject.Position,...
                            'Deletable',0,...
                            'InteractionsAllowed','none',...
                            'LineWidth',h.Class{Index}.Segment(i).ROIObject.LineWidth,...
                            'Color',h.ColorMode(obj.ShowImageSettings.ColorIndex).Disabled);
                    elseif ismember(i,h.ResultStruct.SelectedSegmentIndex{Index})
                        h.ROIObjects{Index}{i} = drawpolyline('Position',h.Class{Index}.Segment(i).ROIObject.Position,...
                            'Deletable',0,...
                            'InteractionsAllowed','none',...
                            'LineWidth',h.Class{Index}.Segment(i).ROIObject.LineWidth*2,...
                            'Label',sprintf('%s || %s',h.Class{Index}.Segment(i).Name,h.Class{Index}.Segment(i).SubSegmentName),...
                            'Color',h.ColorMode(obj.ShowImageSettings.ColorIndex).SpecialHighlight,...
                            'LabelAlpha',1);
                    else
                        h.ROIObjects{Index}{i} = drawpolyline('Position',h.Class{Index}.Segment(i).ROIObject.Position,...
                            'Deletable',0,...
                            'InteractionsAllowed','none',...
                            'LineWidth',h.Class{Index}.Segment(i).ROIObject.LineWidth,...
                            'Color',SegmentColors{2*SegmentIndex-1},...
                            'StripeColor',SegmentColors{2*SegmentIndex});
                    end
                end
                
                create_map_of_rois(Index);
                
            end
            
            function create_map_of_rois(Index,varargin)
                
                Size = size(h.I(1).CData);
                h.ROIImage{Index} = zeros(Size);
                
                for i=1:length(h.Class{Index}.Segment)
                    if ismember(i,h.ResultStruct.ChosenSubsetIndizes{Index})
                        TempMask = h.ROIObjects{Index}{i}.createMask;
                        TempMask = imdilate(TempMask,strel('diamond',ceil(Size(1)/128)));
                    else
                        continue
                    end
                    h.ROIImage{Index}(TempMask) = i;
                end
                
            end
            
            function select_clicked_roi(varargin)
                
                Index = varargin{3};
                
                modifiers = get(h.Fig,'currentModifier');
                ctrlIsPressed = ismember('control',modifiers);
                
                Point = varargin{2}.IntersectionPoint;
                
                X = round(Point(1));
                Y = round(Point(2));
                
                NewIndex = h.ROIImage{Index}(Y,X);
                
                if h.ResultStruct.SelectedSegmentIndex{Index}(1) == NewIndex
                    h.ResultStruct.DisabledSegmentIndizes{Index} =...
                        unique([h.ResultStruct.DisabledSegmentIndizes{Index} NewIndex]);
                    h.ResultStruct.SelectedSegmentIndex{Index} = 0;
                elseif sum(h.ResultStruct.DisabledSegmentIndizes{Index} == NewIndex)
                    h.ResultStruct.DisabledSegmentIndizes{Index}(h.ResultStruct.DisabledSegmentIndizes{Index} == NewIndex) = [];
                elseif ~ctrlIsPressed
                    h.ResultStruct.SelectedSegmentIndex{Index}(1) = NewIndex;
                    if numel(h.ResultStruct.SelectedSegmentIndex{Index}) > 1
                        h.ResultStruct.SelectedSegmentIndex{Index}(2:end) = [];
                    end
                elseif ctrlIsPressed
                    h.ResultStruct.SelectedSegmentIndex{Index}(end+1) = NewIndex;
                end
                
                draw_channel_1
            end
            
            function Results = read_out_results()
                Results = struct();
                
                % Determine how many Sets are to be created from which
                % class instances
                if obj.ShowImageSettings.PlotGroup
                    MotherClass = obj.get_afm_base_class_by_name(h.CurrentClassName{1});
                    ClassNames = MotherClass.OverlayGroup.Names;
                    NumResults = length(ClassNames);
                    if h.ResultStruct.SelectedSegmentIndex{1}(1) ~= 0 ||...
                            length(h.ResultStruct.SelectedSegmentIndex{1}) > 1
                        FullNameCell = MotherClass.get_full_names_from_segment_indizes(h.ResultStruct.SelectedSegmentIndex{1});
                    end
                    for i=1:NumResults
                        Class{i} = obj.get_afm_base_class_by_name(ClassNames{i});
                        ChannelName{i} = h.Channel{1};
                    end
                else
                    if obj.ShowImageSettings.JustChannel1
                        NumResults = 1;
                        ChannelName{1} = h.Channel{1};
                        Class{1} = obj.get_afm_base_class_by_name(h.CurrentClassName{1});
                    elseif obj.ShowImageSettings.JustChannel2
                        NumResults = 1;
                        ChannelName{1} = h.Channel{2};
                        Class{1} = obj.get_afm_base_class_by_name(h.CurrentClassName{2});
                    else
                        NumResults = 2;
                        Class{1} = obj.get_afm_base_class_by_name(h.CurrentClassName{1});
                        Class{2} = obj.get_afm_base_class_by_name(h.CurrentClassName{2});
                        ChannelName{1} = h.Channel{1};
                        ChannelName{2} = h.Channel{2};
                    end
                end
                
                if obj.ShowImageSettings.JustChannel1 || ...
                        obj.ShowImageSettings.PlotGroup
                    Ind = [1 1];
                elseif obj.ShowImageSettings.JustChannel2
                    Ind = [2 2];
                else
                    Ind = [1 2];
                end
                
                for i=NumResults:-1:1
                    Waitbar = waitbar(1/2,['Reading out results: ' num2str(NumResults+1 -i) ' of ' num2str(NumResults)]);
                    TempChannel = Class{i}.get_channel(ChannelName{i});
                    if isempty(TempChannel) &&...
                            NumResults > i &&...
                            i < length(Results)
                        Results(i) = [];
                        close(Waitbar)
                        continue
                    elseif isempty(TempChannel)
                        close(Waitbar)
                        continue
                    end
                    if h.ResultStruct.SelectedSegmentIndex{Ind(min(i,2))}(1) == 0 &&...
                        length(h.ResultStruct.SelectedSegmentIndex{Ind(min(i,2))}) == 1
                        IncludeIndexVector = [];
                    elseif ~obj.ShowImageSettings.PlotGroup
                        IncludeIndexVector = h.ResultStruct.SelectedSegmentIndex{Ind(min(i,2))};
                    else
                        IncludeIndexVector = Class{i}.get_indizes_of_matching_segments(FullNameCell);
                    end
                    [Results(i).Data,...
                        Results(i).SubSegmentCell,...
                        Results(i).SegmentCell,...
                        Results(i).Mask,...
                        Results(i).SegmentNames,...
                        Results(i).SubSegmentNames] = Class{i}.get_segment_data_from_channel(ChannelName{i},...
                        'PixelDilation',0,...
                        'JustSnapped',obj.ShowImageSettings.UseSnapped,...
                        'JustUnsnapped',~obj.ShowImageSettings.UseSnapped,...
                        'MatchString',[],...
                        'IncludeIndexVector',IncludeIndexVector,...
                        'ExcludeIndexVector',h.ResultStruct.DisabledSegmentIndizes{Ind(min(i,2))});
                    Results(i).Name = Class{i}.Name;
                    Results(i).ChannelType = TempChannel.Name;
                    Results(i).Unit = TempChannel.Unit;
                    close(Waitbar)
                end
                
            end
            
            function draw_result_plots()
                
                h.ResultStruct.Gramm = get_and_assign_result_data();
%                 
%                 T = obj.compile_experiment_content_to_table_sqldatabase;
%                 
%                 % Table filtering and sorting
%                 
%                 rows = contains(T.SegmentName,'Snapped');% & T.FiberSegment_SegmentLength > 8e-6;
%                 
%                 SubT = T(rows,:);
%                 
%                 SortBy = 'FiberSegment_RelativePosition';
%                 Order = 'ascend';
%                 XColumn = 'Processed';
%                 YColumn = 'Processed';
%                 CColumn = 'Name';
%                 FColumn = 'SubSegmentFullName';
%                 
%                 ST = sortrows(SubT,SortBy,Order);
%                 
%                 X = ST{:,XColumn};
%                 Y = ST{:,YColumn};
%                 C = ST{:,CColumn};
%                 Facets = ST{:,FColumn};
%                 
%                 XName = ST.Properties.VariableDescriptions{XColumn};
%                 YName = ST.Properties.VariableDescriptions{YColumn};
%                 CName = ST.Properties.VariableDescriptions{CColumn};
%                 FName = ST.Properties.VariableDescriptions{FColumn};
%                 XUnit = ST.Properties.VariableUnits {XColumn};
%                 YUnit = ST.Properties.VariableUnits{YColumn};
%                 CUnit = ST.Properties.VariableUnits{CColumn};
%                 FUnit = ST.Properties.VariableUnits{FColumn};
%                 XText = [XName ' [' XUnit ']'];
%                 YText = [YName ' [' YUnit ']'];
%                 CText = [CName ' [' CUnit ']'];
%                 FText = [FName ' [' FUnit ']'];
%                 
%                 if isstring(X) || iscell(X)
%                     X = categorical(X);
%                 end
%                 if isstring(Y) || iscell(Y)
%                     Y = categorical(Y);
%                 end
%                 if isstring(C) || iscell(C)
%                     C = categorical(C);
%                 end
%                 if isstring(Facets) || iscell(Facets)
%                     Facets = categorical(Facets);
%                 end
%                 
%                 GrammStruct.X = X;
%                 GrammStruct.Y = Y;
%                 GrammStruct.Group = C;
%                 GrammStruct.Facets = Facets;
%                 GrammStruct.UnitX = XText;
%                 GrammStruct.UnitY = YText;
%                 
%                 h.ResultStruct.Gramm = GrammStruct;
                % draw gramm plot
                h.ResultStruct.GrammFig = draw_gramm_plot(h.ResultStruct.Gramm);
            end
            
            function changed_use_snapped(varargin)
                
                obj.ShowImageSettings.UseSnapped = ...
                    h.Res(1).Value;
                
                draw_channel_1
                draw_channel_2
            end
            
            function changed_plot_group(varargin)
                
                obj.ShowImageSettings.PlotGroup = ...
                    h.Res(2).Value;
                
                draw_channel_1
                draw_channel_2
            end
            
            function changed_ignore_zeros(varargin)
                
                obj.ShowImageSettings.IgnoreZeros = ...
                    h.Res(1).Value;
                
                draw_channel_1
                draw_channel_2
            end
            
            function changed_just_channel_1(varargin)
                
                obj.ShowImageSettings.JustChannel1 = ...
                    h.Res(4).Value;
                
                obj.ShowImageSettings.JustChannel2 = 0;
                h.Res(5).Value = 0;
                
                draw_channel_1
                draw_channel_2
            end
            
            function changed_just_channel_2(varargin)
                
                obj.ShowImageSettings.JustChannel2 = ...
                    h.Res(5).Value;
                
                obj.ShowImageSettings.JustChannel1 = 0;
                h.Res(4).Value = 0;
                
                draw_channel_1
                draw_channel_2
            end
            
            function export_data(varargin)
                
                filter = {'*.mat'};
                if obj.ShowImageSettings.JustChannel2
                    Name = split(h.Class{2}.Name,'.');
                else
                    Name = split(h.Class{1}.Name,'.');
                end
                if obj.ShowImageSettings.PlotGroup
                    DefName = ['GroupData_' Name{1}];
                else
                    DefName = ['Data_' Name{1}];
                end
                DefFull = fullfile(obj.ShowImageSettings.DefaultSavePath,DefName);
                [file, path] = uiputfile(filter,'Select Name and Location of your .mat-file',DefFull);
                obj.ShowImageSettings.DefaultSavePath = path;
                FullFile = fullfile(path,file);
                Result = h.ResultStruct.Results;
                save(FullFile,'Result')
                % send to workspace
                assignin('base','Result',Result)
            end
            
            function set_gramm_options(varargin)
                
                obj.set_gramm_options;
                
            end
            
            function redraw_gramm_plot(varargin)
                
                h.ResultStruct.GrammFig = draw_gramm_plot(h.ResultStruct.Gramm);
                
            end
            
            function draw_over_current_gramm_plot(varargin)
                
                h.ResultStruct.GrammFig = draw_gramm_plot(h.ResultStruct.Gramm,h.ResultStruct.GrammFig);
                
            end
            
            function GrammStruct = get_and_assign_result_data()
                
                GrammStruct = struct(...
                    'X',[],...
                    'Y',[],...
                    'Group',[],...
                    'Facets',[],...
                    'UnitX',[],...
                    'UnitY',[]);
                
                h.ResultStruct.Results = read_out_results();
                
                if isempty(h.ResultStruct.Results)
                    return
                end
                
                % Set up data for gramm object. Need to unroll data into single
                % columns. Loop over everything in Data.
                Incr = 1;
                for i=1:length(h.ResultStruct.Results)
                    for j=1:length(h.ResultStruct.Results(i).SegmentNames)
                        for k=1:size(h.ResultStruct.Results(i).Data,1)
                            TempValue = h.ResultStruct.Results(i).Data(k,j);
                            if isnan(TempValue) ||...
                                    (obj.ShowImageSettings.IgnoreZeros&&TempValue==0)
                                continue
                            end
                            Value(Incr) = TempValue;
                            SegmentName{Incr} = h.ResultStruct.Results(i).SegmentNames{j};
                            Name{Incr} = h.ResultStruct.Results(i).Name;
                            ResultIndex(Incr) = i;
                            Incr = Incr + 1;
                        end
                    end
                end
                
                GrammStruct.Group = ResultIndex;
                [MultiplierY,GrammStruct.UnitY,~] = ...
                    AFMImage.parse_unit_scale(range(Value),h.ResultStruct.Results(1).Unit,1);
                ScaledValue = Value.*MultiplierY;
                GrammStruct.Y = ScaledValue;
                GrammStruct.X = SegmentName;
                GrammStruct.Facets = Name;
                
            end
            
            function g = draw_gramm_plot(GrammStruct,GrammObject)
                
                X = GrammStruct.X;
                Y = GrammStruct.Y;
                Group = GrammStruct.Group;
                Facets = GrammStruct.Facets;
                UnitX = GrammStruct.UnitX;
                UnitY = GrammStruct.UnitY;
                
                if nargin == 1
                    g = gramm('x',X,'y',Y,obj.GrammOptions.GroupingStyle,Group);
                elseif nargin == 2
                    g = GrammObject;
                    g.update();
                end
                
                if obj.GrammOptions.SubplotOptions.EnableFacets
                    g.facet_wrap(Facets,...
                        'ncols',obj.GrammOptions.SubplotOptions.ncols,...
                        'scale',obj.GrammOptions.SubplotOptions.scale,...
                        'column_labels',obj.GrammOptions.SubplotOptions.column_labels,...
                        'force_ticks',obj.GrammOptions.SubplotOptions.force_ticks);
                end
                
                if isequal(obj.GrammOptions.VisualizationType,'Statistical')
                    
                    switch obj.GrammOptions.StatVisOptions.Style
                        case 'Summary'
                            if obj.GrammOptions.StatVisOptions.black_errorbar
                                Sum_Geom = {obj.GrammOptions.StatVisOptions.Summary_geom,'black_errorbar'};
                            else
                                Sum_Geom = obj.GrammOptions.StatVisOptions.Summary_geom;
                            end
                            g.stat_summary('geom',Sum_Geom,...
                                'type',obj.GrammOptions.StatVisOptions.Summary_type,...
                                'setylim',obj.GrammOptions.StatVisOptions.Summary_setylim,...
                                'width',obj.GrammOptions.StatVisOptions.Summary_width,...
                                'dodge',obj.GrammOptions.StatVisOptions.Summary_dodge);
                        case 'Smooth'
                            g.stat_smooth('geom',obj.GrammOptions.StatVisOptions.Smooth_geom,...
                                'method',obj.GrammOptions.StatVisOptions.Smooth_method,...
                                'lambda',obj.GrammOptions.StatVisOptions.Smooth_lambda,...
                                'npoints',obj.GrammOptions.StatVisOptions.Smooth_npoints);
                        case 'glm'
                            g.stat_glm('geom',obj.GrammOptions.StatVisOptions.Glm_geom,...
                                'distribution',obj.GrammOptions.StatVisOptions.Glm_distribution,...
                                'fullrange',obj.GrammOptions.StatVisOptions.Glm_fullrange,...
                                'disp_fit',obj.GrammOptions.StatVisOptions.Glm_disp_fit);
                        case 'fit'
                            g.stat_fit('geom',obj.GrammOptions.StatVisOptions.Fit_geom,...
                                'fun',obj.GrammOptions.StatVisOptions.Fit_fun,...
                                'intopt',obj.GrammOptions.StatVisOptions.Fit_intopt,...
                                'fullrange',obj.GrammOptions.StatVisOptions.Fit_fullrange,...
                                'disp_fit',obj.GrammOptions.StatVisOptions.Fit_disp_fit);
                        case 'bin'
                            g.stat_bin('nbins',obj.GrammOptions.StatVisOptions.Bin_nbins,...
                                'geom',obj.GrammOptions.StatVisOptions.Bin_geom,...
                                'normalization',obj.GrammOptions.StatVisOptions.Bin_normalization,...
                                'fill',obj.GrammOptions.StatVisOptions.Bin_fill,...
                                'width',obj.GrammOptions.StatVisOptions.Bin_width,...
                                'dodge',obj.GrammOptions.StatVisOptions.Bin_dodge);
                        case 'cornerhist'
                            g.stat_cornerhist('location',obj.GrammOptions.StatVisOptions.Cornerhist_location,...
                                'aspect',obj.GrammOptions.StatVisOptions.Cornerhist_aspect,...
                                'nbins',obj.GrammOptions.StatVisOptions.Cornerhist_nbins,...
                                'geom',obj.GrammOptions.StatVisOptions.Cornerhist_geom,...
                                'normalization',obj.GrammOptions.StatVisOptions.Cornerhist_normalization,...
                                'fill',obj.GrammOptions.StatVisOptions.Cornerhist_fill,...
                                'width',obj.GrammOptions.StatVisOptions.Cornerhist_width,...
                                'dodge',obj.GrammOptions.StatVisOptions.Cornerhist_dodge);
                        case 'density'
                            g.stat_density('bandwidth',obj.GrammOptions.StatVisOptions.Density_bandwidth,...
                                'function',obj.GrammOptions.StatVisOptions.Density_function,...
                                'kernel',obj.GrammOptions.StatVisOptions.Density_kernel,...
                                'npoints',obj.GrammOptions.StatVisOptions.Density_npoints,...
                                'extra_x',obj.GrammOptions.StatVisOptions.Density_extra_x);
                        case 'bin2d'
                            g.stat_bin2d('geom',obj.GrammOptions.StatVisOptions.Bin2D_geom);
                        case 'ellipse'
                            g.stat_ellipse('type',obj.GrammOptions.StatVisOptions.Ellipse_type,...
                                'geom',obj.GrammOptions.StatVisOptions.Ellipse_geom);
                        case 'qq'
                            g.stat_qq('distribution',eval(obj.GrammOptions.StatVisOptions.QQ_distribution));
                        case 'boxplot'
                            g.stat_boxplot('width',obj.GrammOptions.StatVisOptions.Boxplot_width,...
                                'dodge',obj.GrammOptions.StatVisOptions.Boxplot_dodge,...
                                'notch',obj.GrammOptions.StatVisOptions.Boxplot_notch);
                        case 'violin'
                            g.stat_violin('normalization',obj.GrammOptions.StatVisOptions.Violin_normalization,...
                                'half',obj.GrammOptions.StatVisOptions.Violin_half,...
                                'width',obj.GrammOptions.StatVisOptions.Violin_width,...
                                'dodge',obj.GrammOptions.StatVisOptions.Violin_dodge,...
                                'fill',obj.GrammOptions.StatVisOptions.Violin_fill,...
                                'bandwidth',obj.GrammOptions.StatVisOptions.Violin_bandwidth,...
                                'kernel',obj.GrammOptions.StatVisOptions.Violin_kernel,...
                                'npoints',obj.GrammOptions.StatVisOptions.Violin_npoints,...
                                'extra_y',obj.GrammOptions.StatVisOptions.Violin_extra_y);
                    end
                elseif isequal(obj.GrammOptions.VisualizationType,'Direct')
                    switch obj.GrammOptions.DirVisOptions.Style
                        case 'Point'
                            g.geom_point('dodge',obj.GrammOptions.DirVisOptions.Dodge,...
                                'alpha',obj.GrammOptions.DirVisOptions.Alpha);
                        case 'Jitter'
                            g.geom_jitter('dodge',obj.GrammOptions.DirVisOptions.Dodge,...
                                'alpha',obj.GrammOptions.DirVisOptions.Alpha,...
                                'width',obj.GrammOptions.DirVisOptions.Width,...
                                'height',obj.GrammOptions.DirVisOptions.Height);
                        case 'Line'
                            g.geom_line('dodge',obj.GrammOptions.DirVisOptions.Dodge,...
                                'alpha',obj.GrammOptions.DirVisOptions.Alpha);
                        case 'Raster'
                            g.geom_raster();
                        case 'Bar'
                            g.geom_bar('dodge',obj.GrammOptions.DirVisOptions.Dodge,...
                                'width',obj.GrammOptions.DirVisOptions.Width,...
                                'stacked',obj.GrammOptions.DirVisOptions.Stacked);
                        case 'Interval'
                            g.geom_interval('dodge',obj.GrammOptions.DirVisOptions.Dodge,...
                                'width',obj.GrammOptions.DirVisOptions.Width);
                    end
                end
                
                if obj.GrammOptions.NoLegend
                    g.no_legend();
                end
                if obj.GrammOptions.FlipCoordinates
                    g.coord_flip();
                end
                
                g.set_color_options('map',obj.GrammOptions.Colorset,...
                    'lightness',obj.GrammOptions.lightness,...
                    'chroma',obj.GrammOptions.chroma,...
                    'legend',obj.GrammOptions.color_legend);
                
                if nargin == 1
                    YName = [h.ResultStruct.Results(1).ChannelType ' [' UnitY ']'];
                    g.set_title('Placeholder');
                    
                    g.set_names('x','Fibril','y',YName,'color','# Cyl');
                    g.set_title(obj.ExperimentName);
                    
                    
                    FigPos = h.Fig.InnerPosition;
                    FigPos(3:4) = FigPos(3:4).*.8;
                    g.set_parent(figure('Position',FigPos));
                    
                    g.draw();
                    
                elseif nargin == 2
                    
                    g.draw();
                    
                end
                
            end
            
            uiwait(h.Fig)
        end
        
        function choose_segments_manually(obj,SegmentType)
            
            FiberLikeMenu = true;
            GeneralMenu = false;
            
            h.ColorMode(1).Background = 'k';
            h.ColorMode(1).Profile1 = [219 21 223]./255; %[189 0 96]./255; % 'b';
            h.ColorMode(1).Profile2 = 'c';
            h.ColorMode(1).Text = 'w';
            
            
            h.ColorMode(2).Background = 'w';
            h.ColorMode(2).Profile1 = [219 21 223]./255; %[94 170 170]./255; % [189 0 96]./255; %'b';
            h.ColorMode(2).Profile2 = [0 181 26]./255; % alternatives %[80 200 204]./255;%[0,0.870588235294118,0.407843137254902];
            h.ColorMode(2).Text = 'k';
            
            h.Fig = figure('Name',sprintf('%s',obj.ExperimentName),...
                'Units','pixels',...
                'Position',[200 200 1024 512],...
                'Color',h.ColorMode(obj.ShowImageSettings.ColorIndex).Background);
            
            h.B(1) = uicontrol('style','togglebutton',...
                'String','Cross Section',...
                'units','normalized',...
                'position',[.85 .45 .1 .05],...
                'Callback',@cross_section_toggle,...
                'Visible','off');
            
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
                'position',[.85 .7 .15 .04],...
                'Visible','off');
            
            h.B(17) = uicontrol('style','popupmenu',...
                'String',ClassPopUp,...
                'units','normalized',...
                'position',[.85 .65 .15 .05],...
                'Callback',@draw_channel_2,...
                'Visible','off');
            
            h.B(3) = uicontrol('style','popupmenu',...
                'String',PopUp,...
                'units','normalized',...
                'position',[.85 .6 .15 .05],...
                'Callback',@draw_channel_2,...
                'Visible','off');
            
            h.B(6) = uicontrol('style','pushbutton',...
                'String','Save Figure',...
                'units','normalized',...
                'position',[.8 .05 .1 .04],...
                'Callback',@save_figure_to_file);
            
            h.B(7) = uicontrol('style','checkbox',...
                'String','...with white background',...
                'units','normalized',...
                'position',[.875 .01 .1 .04],...
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
                'Callback',@changed_slider,...
                'Visible','off');
            
            h.B(11) = uicontrol('style','slider',...
                'Value',0,...
                'Units','normalized',...
                'Position',[.85 .56 .1 .02],...
                'Callback',@changed_slider,...
                'Visible','off');
            
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
                'Position',[.95 .58 .03 .02],...
                'Visible','off');
            
            h.B(15) = uicontrol('style','text',...
                'String','Min',...
                'Units','normalized',...
                'Position',[.95 .56 .03 .02],...
                'Visible','off');
            
            h.B(18) = uicontrol('Style','checkbox',...
                'String','Both Channels',...
                'Value',0,...
                'Tooltip','Green, if both Channels have the same size scaling',...
                'Units','normalized',...
                'Position',[.85 .42 .1 .03],...
                'Callback',@checked_both_cross_sections,...
                'Visible','off');
            
            h.B(19) = uicontrol('style','checkbox',...
                'String','Upscale Images',...
                'units','normalized',...
                'position',[.85 .15 .1 .04],...
                'Callback',@upscale_images,...
                'Visible','off');
            
            h.B(20) = uicontrol('style','checkbox',...
                'String','Lock Channels',...
                'units','normalized',...
                'position',[.85 .3   .1 .04],...
                'Callback',@lock_channels,...
                'Visible','off');
            
            h.B(21) = uicontrol('style','checkbox',...
                'String','Lock Scalebars',...
                'units','normalized',...
                'position',[.85 .20 .1 .04],...
                'Callback',@lock_scalebars,...
                'Visible','off');
            
            h.B(22) = uicontrol('style','checkbox',...
                'String','Statistical CMapping',...
                'units','normalized',...
                'position',[.85 .25 .1 .04],...
                'Callback',@statistical_cmapping,...
                'Visible','off');
            
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
                'Callback',@set_scale,...
                'Visible','off');
            
            h.B(26) = uicontrol('style','edit',...
                'String','',...
                'units','normalized',...
                'position',[.925 .5 .035 .05],...
                'Callback',@set_scale,...
                'Visible','off');
            
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
                'Position',[.885 .5 .04 .05],...
                'Visible','off');
            
            h.B(30) = uicontrol('style','text',...
                'String',{'Max','[]'},...
                'Units','normalized',...
                'Position',[.96 .5 .04 .05],...
                'Visible','off');
            
            h.B(31) = uicontrol('Style','checkbox',...
                'String','Use Overlay',...
                'Value',0,...
                'Tooltip','Green, if both Channels share an overlay',...
                'Units','normalized',...
                'Position',[.85 .39 .1 .03],...
                'Visible','off');
            
            h.c(32) = uicontrol(h.Fig,'style','pushbutton','units','normalized',...
                'position',[.85 .15 .07 .04],'string','Delete Segment',...
                'Callback',@delete_selected_segment);
            
            h.c(33) = uicontrol(h.Fig,'style','pushbutton','units','normalized',...
                'position',[.925 .15 .07 .04],'string','Delete Subsegment',...
                'Callback',@delete_selected_subsegment);
            
            h.c(34) = uicontrol(h.Fig,'style','text','units','normalized',...
                'position',[.85 .7 .15 .04],'string','List of Segments');
            
            h.c(35) = uicontrol(h.Fig,'style','text','units','normalized',...
                'position',[.85 .55 .15 .04],'string','List of Subsegments');
            
            h.c(36) = uicontrol(h.Fig,'style','pushbutton','units','normalized',...
                'position',[.85 .4 .07 .04],'string','New Segment',...
                'Callback',@create_new_segment);
            
            h.c(37) = uicontrol(h.Fig,'style','radiobutton','units','normalized',...
                'position',[.925 .4 .07 .04],'string','Enable Drawing',...
                'Callback',@draw_channel_1);
            
            h.c(38) = uicontrol(h.Fig,'style','radiobutton','units','normalized',...
                'position',[.925 .1 .07 .04],'string','Show Labels',...
                'Callback',@show_labels);
            
            h.c(39) = uicontrol(h.Fig,'style','pushbutton','units','normalized',...
                'position',[.85 .35 .07 .04],'string','Snap Segment',...
                'Callback',@snap_segment);
            
            h.c(40) = uicontrol(h.Fig,'style','pushbutton','units','normalized',...
                'position',[.925 .35 .07 .04],'string','Snap All',...
                'Callback',@snap_all);
            
            h.c(41) = uicontrol('style','edit',...
                'String','10e-9',...
                'units','normalized',...
                'position',[.925 .31 .07 .03]);
            
            h.c(42) = uicontrol('style','edit',...
                'String','300e-9',...
                'units','normalized',...
                'position',[.925 .28 .07 .03]);
            
            h.c(43) = uicontrol('style','edit',...
                'String','21',...
                'units','normalized',...
                'position',[.925 .25 .07 .03]);
            
            h.c(44) = uicontrol('style','edit',...
                'String','1',...
                'units','normalized',...
                'position',[.925 .22 .07 .03]);
            
            h.c(45) = uicontrol('style','text',...
                'String','Distance between Sampling Points',...
                'Units','normalized',...
                'Position',[.85 .31 .07 .03],...
                'Visible','on');
            
            h.c(46) = uicontrol('style','text',...
                'String','Local Window Size',...
                'Tooltip','Width of local perpendicular profile in meters',...
                'Units','normalized',...
                'Position',[.85 .28 .07 .03],...
                'Visible','on');
            
            h.c(47) = uicontrol('style','text',...
                'String','#Points in Moving Window',...
                'Tooltip','Number of Points in moving window used in local regression',...
                'Units','normalized',...
                'Position',[.85 .25 .07 .03],...
                'Visible','on');
            
            h.c(48) = uicontrol('style','text',...
                'String','Degree of Polyfit',...
                'Tooltip','Degree of Polyfit',...
                'Units','normalized',...
                'Position',[.85 .22 .07 .03],...
                'Visible','on');
            
            h.c(49) = uicontrol(h.Fig,'style','pushbutton','units','normalized',...
                'position',[.85 .1 .07 .04],'string','Apply Segmentation to Overlay Group',...
                'Tooltip','Apply Segmentation to Overlay Group',...
                'Callback',@apply_segmentation_to_overlay_group);
            
            
            h.c(50) = uicontrol(h.Fig,'style','pushbutton','units','normalized',...
                'position',[.91 .05 .1 .04],'string','Switch Drawmode',...
                'Tooltip',{'Switch between general drawmode for all kinds of segment shapes (line, polyline, circle, freehand, etc...) and the polyline drawmode with additional tools for fiber like object snapping'},...
                'Callback',@switch_drawmode);
            
            DrawOptions = {'Line','Polyline','Freehand','Circle',...
                'Ellipse','Polygon','Rectangle','Assisted'};
            
            h.c(51) = uicontrol('style','popupmenu',...
                'String',DrawOptions,...
                'units','normalized',...
                'position',[.85 .36 .15 .05],...
                'Visible','off');
            
            % Tagsystem editbox
            h.c(52) = uicontrol('style','edit',...
                'String','',...
                'units','normalized',...
                'position',[.85 .33 .1 .03],...
                'Visible','off');
            
            h.c(53) = uicontrol('style','text',...
                'String','Tag',...
                'units','normalized',...
                'position',[.95 .33 .045 .03],...
                'Tooltip','Fill in a tag for the current segment and choose how to add it below',...
                'Visible','off');
            
            % Tagsystem add tag title
            h.c(54) = uicontrol('style','text',...
                'String','Add tag to | Delete tag from...',...
                'Tooltip','Choose to which segments the current tag in the edit field should be added to',...
                'Units','normalized',...
                'Position',[.85 .31 .15 .02],...
                'Visible','off');
            
            GroupTagOptions = {'Current Segment','Identical Segments in Ov. Group',...
                'Identical Segments in Experiment',...
                'All Segments in image','All Segments in Ov. Group','All Segments in Experiment'};
            
            % Tagsystem add to Segment button
            h.c(55) = uicontrol('style','popupmenu',...
                'String',GroupTagOptions,...
                'units','normalized',...
                'position',[.85 .28 .15 .03],...
                'Visible','off');
            
            h.c(56) = uicontrol('style','pushbutton',...
                'String','Add tag',...
                'units','normalized',...
                'position',[.85 .25 .07 .03],...
                'Tooltip','Add the tag currenty in the edit box to the segments chosen in the popup above',...
                'Visible','off',...
                'Callback',@add_tag);
            
            h.c(57) = uicontrol('style','pushbutton',...
                'String','Delete tag',...
                'Tooltip','Delete the tag currently highlighted in the tag box from the segments chosen in the popup above',...
                'units','normalized',...
                'position',[.93 .25 .07 .03],...
                'Visible','off',...
                'Callback',@delete_tag);
            
            h.c(58) = uicontrol(h.Fig,...
                'Style','listbox',...
                'Max',1000000,'Min',1,...
                'Units','normalized',...
                'Position',[.85  .18 .15 .07],...
                'Visible','off');
            
            
            h.SegmentBox = uicontrol(h.Fig,...
                'Style','listbox',...
                'Max',1000000,'Min',1,...
                'Units','normalized',...
                'Position',[.85  .6 .15 .1],...
                'Callback',@segment_changed);
            
            h.SubsegmentBox = uicontrol(h.Fig,...
                'Style','listbox',...
                'Max',1000000,'Min',1,...
                'Units','normalized',...
                'Position',[.85  .45 .15 .1],...
                'Callback',@subsegment_changed);
            
            h.Channel1Max = 1;
            h.Channel1Min = 0;
            h.Channel2Max = 1;
            h.Channel2Min = 0;
            
            h.CurChannel1Idx = h.B(2).Value;
            h.CurChannel2Idx = h.B(3).Value;
            obj.ShowImageSettings.RelativeChannelIndex = 0;
            
            h.MainLine = [];
            h.ChildLine = [];
            h.hasCrossSection = 0;
            obj.ShowImageSettings.BothCrossSections = 0;
            h.hasChannel2 = 0;
            obj.ShowImageSettings.IsUpscaled = false;
            obj.ShowImageSettings.LockChannels = h.B(20).Value;
            [~,DefIndex] = Class{1}.get_channel('Processed');
            
            % set initial listbox items
            if isempty(Class{1}.Segment(1).Name)
                set(h.SegmentBox,'String',[]);
                set(h.SubsegmentBox,'String',[]);
            else
                InitialNames = {Class{1}.Segment.Name};
                InitialSubSegmentNames = {Class{1}.Segment.SubSegmentName};
                InitialUniqueNames = unique(InitialNames);
                set(h.SegmentBox,'String',InitialUniqueNames);
                set(h.SubsegmentBox,'String',{InitialSubSegmentNames{strcmp({h.SegmentBox.String{h.SegmentBox.Value}},InitialNames)}});
            end
            h.CurrentClassName = Class{1}.Name;
            h.EnableEditing = 1;
            
            if isempty(DefIndex)
                DefIndex = 2;
            else
                DefIndex = DefIndex + 1;
            end
            h.B(2).Value = DefIndex;
            draw_channel_1
            
            function switch_drawmode(varargin)
                
                GeneralIndizes = [51:58];
                FiberLikeIndizes = [39:48];
                
                if FiberLikeMenu
                    for i=FiberLikeIndizes
                        h.c(i).Visible = 'Off';
                    end
                    for i=GeneralIndizes
                        h.c(i).Visible = 'On';
                    end
                    FiberLikeMenu = false;
                    GeneralMenu = true;
                else
                    for i=FiberLikeIndizes
                        h.c(i).Visible = 'On';
                    end
                    for i=GeneralIndizes
                        h.c(i).Visible = 'Off';
                    end
                    FiberLikeMenu = true;
                    GeneralMenu = false;
                end
                
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
            
            function draw_image(LeftRight,FullPart)
                if isequal(LeftRight,'Left')
                    Index = 1;
                elseif isequal(LeftRight,'Right')
                    Index = 2;
                end
                BarToImageRatio = 1/5;
                try
                    CurrentZoomX = h.ImAx(Index).XLim; 
                    CurrentZoomY = h.ImAx(Index).YLim;
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
                
                
                if obj.ShowImageSettings.LockChannels && Index==1
                    h.CurChannel1Idx = h.B(16).Value;
                    h.B(17).Value = mod(h.CurChannel1Idx + obj.ShowImageSettings.RelativeChannelIndex,h.NumClasses+1);
                    h.CurChannel2Idx = h.B(17).Value;
                    CurIndex = h.CurChannel1Idx;
                elseif obj.ShowImageSettings.LockChannels && Index==2
                    h.CurChannel2Idx = h.B(17).Value;
                    h.B(16).Value = mod(h.CurChannel2Idx - obj.ShowImageSettings.RelativeChannelIndex,h.NumClasses+1);
                    h.CurChannel1Idx = h.B(16).Value;
                    CurIndex = h.CurChannel2Idx;
                else
                    CurIndex = h.B(15+Index).Value;
                    h.CurChannel2Idx = h.B(17).Value;
                    h.CurChannel1Idx = h.B(16).Value;
                end
                
                obj.ShowImageSettings.RelativeChannelIndex = h.CurChannel2Idx - h.CurChannel1Idx;
                
                Class{Index} = obj.get_class_instance(ClassIndex(CurIndex,:));
                if ~isequal(h.CurrentClassName,Class{1}.Name)
                    h.SegmentBox.Value = 1;
                    h.SubsegmentBox.Value = 1;
                end
                CurrentChannelName = h.B(1+Index).String{h.B(1+Index).Value};
                PopUp = Class{Index}.string_of_existing();
                set(h.B(1+Index),'String',PopUp)
                
                IdxOfSameChannel = find(strcmp(PopUp,CurrentChannelName));
                
                if isempty(IdxOfSameChannel)
                    set(h.B(1+Index),'Value',2);
                elseif sum(IdxOfSameChannel == h.B(1+Index).Value)
                    set(h.B(1+Index),'Value',IdxOfSameChannel(find(IdxOfSameChannel == h.B(1+Index).Value)));
                else
                    set(h.B(1+Index),'Value',IdxOfSameChannel(1));
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
                    if obj.ShowImageSettings.IsUpscaled
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
                if h.c(37).Value
                    h.I(Index).ButtonDownFcn = @create_new_subsegment;
                else
                    h.I(Index).ButtonDownFcn = @select_clicked_roi;
                end
                hold on
                CurrentAxHeight = round(h.Fig.Position(4)*h.ImAx(Index).Position(4));
                CurrentAxWidth = round(h.Fig.Position(3)*h.ImAx(Index).Position(3));
                %AFMImage.draw_scalebar_into_current_image(Channel.NumPixelsX,Channel.NumPixelsY,Channel.ScanSizeX,BarToImageRatio,CurrentAxHeight,CurrentAxWidth);
                c = colorbar('northoutside');
                c.FontSize = round(obj.ShowImageSettings.ReferenceFontSize*(CurrentAxHeight/756));
                c.Color = h.ColorMode(obj.ShowImageSettings.ColorIndex).Text;
                c.Label.String = sprintf('%s [%s]',h.Channel{Index},h.Unit{Index});
                c.Label.FontSize = round(obj.ShowImageSettings.ReferenceFontSize*(CurrentAxHeight/756));
                c.Label.Color = h.ColorMode(obj.ShowImageSettings.ColorIndex).Text;
                
                set(h.B(MinIndex+4),'String',{'Min',sprintf('[%s]',h.Unit{Index})});
                set(h.B(MaxIndex+4),'String',{'Max',sprintf('[%s]',h.Unit{Index})});
                
                
                if isequal(h.CurrentClassName,Class{Index}.Name)
                    try
                        % Set Zoom region to previous
                        zoom reset
                        h.ImAx(Index).XLim = CurrentZoomX;
                        h.ImAx(Index).YLim = CurrentZoomY;
                    catch
                    end
                else
                    h.CurrentClassName = Class{Index}.Name;
                end
                
                try
                    if (h.ScanSizeX(1) == h.ScanSizeX(2)) &&...
                            (h.ScanSizeY(1) == h.ScanSizeY(2))
                        h.B(18).BackgroundColor = 'g';
                    else
                        h.B(18).BackgroundColor = 'r';
                    end
                catch
                end
                check_for_overlay_group
                
                % Segmentation related
                if isempty(Class{1}.Segment(1).Name)
                    set(h.SegmentBox,'String',[]);
                    set(h.SubsegmentBox,'String',[]);
                else
                    Names = {Class{1}.Segment.Name};
                    SubSegmentNames = {Class{1}.Segment.SubSegmentName};
                    UniqueNames = unique(Names);
                    set(h.SegmentBox,'String',UniqueNames);
                    set(h.SubsegmentBox,'String',{SubSegmentNames{strcmp({h.SegmentBox.String{h.SegmentBox.Value}},Names)}});
                end
                draw_existing_segments
                update_tag_box
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
                
                filter = {obj.ShowImageSettings.DefaultSaveType;'*.eps';'*.emf';'*.png';'*.tif';'*.jpg'};
                Name = split(Class{1}.Name,'.');
                DefName = Name{1};
                DefFull = fullfile(obj.ShowImageSettings.DefaultSavePath,DefName);
                [file, path] = uiputfile(filter,'Select Format, Name and Location of your figure',DefFull);
                obj.ShowImageSettings.DefaultSavePath = path;
                FullFile = fullfile(path,file);
                Split = split(FullFile,'.');
                obj.ShowImageSettings.DefaultSaveType = strcat('*.',Split{end});
                if isequal(Split{end},'eps') || isequal(Split{end},'emf')
                    exportgraphics(h.Fig,FullFile,'ContentType','vector','Resolution',300,'BackgroundColor','current')
                else
                    exportgraphics(h.Fig,FullFile,'Resolution',300,'BackgroundColor','current')
                end
                    
            end
            
            function changed_color(varargin)
                
                if ~h.B(7).Value
                    obj.ShowImageSettings.ColorIndex = 1;
                else
                    obj.ShowImageSettings.ColorIndex = 2;
                end
                
                h.Fig.Color = h.ColorMode(obj.ShowImageSettings.ColorIndex).Background;
                
                draw_channel_1
                draw_channel_2
                
            end
            
            function upscale_images(varargin)
                
                obj.ShowImageSettings.IsUpscaled = h.B(19).Value;
                
                draw_channel_1
                draw_channel_2
                
            end
            
            function lock_channels(varargin)
                
                obj.ShowImageSettings.LockChannels = h.B(20).Value;
                
                obj.ShowImageSettings.RelativeChannelIndex = h.CurChannel2Idx - h.CurChannel1Idx;
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
            
            function check_for_overlay_group(varargin)
                
                
                if h.hasChannel2 &&...
                            ~isempty(Class{1}.OverlayGroup) &&...
                            ~isempty(Class{2}.OverlayGroup) &&...
                            Class{1}.OverlayGroup.hasOverlayGroup &&...
                            Class{2}.OverlayGroup.hasOverlayGroup &&...
                            isequal(Class{1}.OverlayGroup.Names,Class{2}.OverlayGroup.Names)
                    h.B(31).BackgroundColor = 'g';
                else
                    h.B(31).BackgroundColor = 'r';
                end
                
            end
            
            function select_clicked_roi(varargin)
                
                Point = varargin{2}.IntersectionPoint;
                
                X = round(Point(1));
                Y = round(Point(2));
                
                NewIndex = h.ROIImage(Y,X);
                
                if ~NewIndex
                    return
                end
                
                Names = {Class{1}.Segment.Name};
                SubSegmentNames = {Class{1}.Segment.SubSegmentName};
                UniqueNames = unique(Names);
                
                h.SegmentBox.Value = find(strcmp({Class{1}.Segment(NewIndex).Name},UniqueNames));
                
                set(h.SegmentBox,'String',UniqueNames);
                SSBNames = {SubSegmentNames{strcmp({h.SegmentBox.String{h.SegmentBox.Value}},Names)}};
                set(h.SubsegmentBox,'String',SSBNames);
                
                h.SubsegmentBox.Value = find(strcmp({Class{1}.Segment(NewIndex).SubSegmentName},SSBNames));
                
                draw_channel_1
            end
            
            function draw_existing_segments(varargin)
                
                if isempty(Class{1}.Segment)
                    return
                end
                if isempty(Class{1}.Segment(1).Name)
                    return
                end
                
                Names = {Class{1}.Segment.Name};
                UniqueNames = unique(Names);
                NumSegments = length(UniqueNames);
                
                % This guy figured out one simple trick to get all
                % the colors. Painters hate him
                FakePlots = zeros(2,NumSegments*2);
                FakeAx = plot(FakePlots);
                SegmentColors = get(FakeAx,'Color');
                delete(FakeAx)
                clear FakeAx FakePlots
                
                NumSubSegments = length(Names);
                
                for i=1:NumSubSegments
                    SegmentIndex = find(strcmp(UniqueNames,Class{1}.Segment(i).Name));
                    if isempty(Class{1}.Segment(i).ROIObject) || ...
                            (isequal(Class{1}.Segment(i).Type,'polyline') && size(Class{1}.Segment(i).ROIObject.Position,1) < 2) ||...
                            (~isequal(Class{1}.Segment(i).Type,'polyline') && size(Class{1}.Segment(i).ROIObject.Position,2) < 2)
                        Class{1}.Segment(i) = [];
                        draw_channel_1
                    end
                    if isequal(Class{1}.Segment(i).Name,h.SegmentBox.String{h.SegmentBox.Value}) &&...
                            isequal(Class{1}.Segment(i).SubSegmentName,h.SubsegmentBox.String{h.SubsegmentBox.Value}) &&...
                            h.EnableEditing
                        h.CurrentIndex = i;
                        CurrentDrawMode = lower(Class{1}.Segment(i).Type);
                        switch CurrentDrawMode
                            case 'line'
                                h.ROIObjects{i} = drawline('Position',Class{1}.Segment(i).ROIObject.Position,...
                                    'Deletable',1,...
                                    'InteractionsAllowed','all',...
                                    'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                    'Label',sprintf('%s || %s',Class{1}.Segment(i).Name,Class{1}.Segment(i).SubSegmentName),...
                                    'LabelAlpha',0.6);
                            case 'freehand'
                                h.ROIObjects{i} = drawfreehand('Position',Class{1}.Segment(i).ROIObject.Position,...
                                    'Deletable',1,...
                                    'InteractionsAllowed','all',...
                                    'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                    'Label',sprintf('%s || %s',Class{1}.Segment(i).Name,Class{1}.Segment(i).SubSegmentName),...
                                    'LabelAlpha',0.6);
                            case 'square'
                                h.ROIObjects{i} = drawsquare('Position',Class{1}.Segment(i).ROIObject.Position,...
                                    'Deletable',1,...
                                    'InteractionsAllowed','all',...
                                    'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                    'Label',sprintf('%s || %s',Class{1}.Segment(i).Name,Class{1}.Segment(i).SubSegmentName),...
                                    'LabelAlpha',0.6);
                            case 'circle'
                                h.ROIObjects{i} = drawcircle('Center',Class{1}.Segment(i).ROIObject.Center,...
                                    'Radius',Class{1}.Segment(i).ROIObject.Radius,...
                                    'Deletable',1,...
                                    'InteractionsAllowed','all',...
                                    'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                    'Label',sprintf('%s || %s',Class{1}.Segment(i).Name,Class{1}.Segment(i).SubSegmentName),...
                                    'LabelAlpha',0.6);
                            case 'ellipse'
                                h.ROIObjects{i} = drawellipse('Center',Class{1}.Segment(i).ROIObject.Center,...
                                    'SemiAxes',Class{1}.Segment(i).ROIObject.SemiAxes,...
                                    'RotationAngle',Class{1}.Segment(i).ROIObject.RotationAngle,...
                                    'Deletable',1,...
                                    'InteractionsAllowed','all',...
                                    'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                    'Label',sprintf('%s || %s',Class{1}.Segment(i).Name,Class{1}.Segment(i).SubSegmentName),...
                                    'LabelAlpha',0.6);
                            case 'polygon'
                                h.ROIObjects{i} = drawpolygon('Position',Class{1}.Segment(i).ROIObject.Position,...
                                    'Deletable',1,...
                                    'InteractionsAllowed','all',...
                                    'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                    'Label',sprintf('%s || %s',Class{1}.Segment(i).Name,Class{1}.Segment(i).SubSegmentName),...
                                    'LabelAlpha',0.6);
                            case 'rectangle'
                                h.ROIObjects{i} = drawrectangle('Position',Class{1}.Segment(i).ROIObject.Position,...
                                    'RotationAngle',Class{1}.Segment(i).ROIObject.RotationAngle,...
                                    'Rotatable',true,...
                                    'Deletable',1,...
                                    'InteractionsAllowed','all',...
                                    'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                    'Label',sprintf('%s || %s',Class{1}.Segment(i).Name,Class{1}.Segment(i).SubSegmentName),...
                                    'LabelAlpha',0.6);
                            case 'crosshair'
                                h.ROIObjects{i} = drawcrosshair('Position',Class{1}.Segment(i).ROIObject.Position,...
                                    'Deletable',1,...
                                    'InteractionsAllowed','all',...
                                    'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                    'Label',sprintf('%s || %s',Class{1}.Segment(i).Name,Class{1}.Segment(i).SubSegmentName),...
                                    'LabelAlpha',0.6);
                            case 'point'
                                h.ROIObjects{i} = drawpoint('Position',Class{1}.Segment(i).ROIObject.Position,...
                                    'Deletable',1,...
                                    'InteractionsAllowed','all',...
                                    'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                    'Label',sprintf('%s || %s',Class{1}.Segment(i).Name,Class{1}.Segment(i).SubSegmentName),...
                                    'LabelAlpha',0.6);
                            case 'assisted'
                                h.ROIObjects{i} = drawassisted('Position',Class{1}.Segment(i).ROIObject.Position,...
                                    'Deletable',1,...
                                    'InteractionsAllowed','all',...
                                    'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                    'Label',sprintf('%s || %s',Class{1}.Segment(i).Name,Class{1}.Segment(i).SubSegmentName),...
                                    'LabelAlpha',0.6);
                            case 'polyline'
                                h.ROIObjects{i} = drawpolyline('Position',Class{1}.Segment(i).ROIObject.Position,...
                                    'Deletable',1,...
                                    'InteractionsAllowed','all',...
                                    'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                    'Label',sprintf('%s || %s',Class{1}.Segment(i).Name,Class{1}.Segment(i).SubSegmentName),...
                                    'LabelAlpha',0.6);
                        end
                            addlistener(h.ROIObjects{i},'ROIMoved',@moved_roi);
                        if isequal(CurrentDrawMode,'polyline')
                            addlistener(h.ROIObjects{i},'VertexAdded',@moved_roi);
                            addlistener(h.ROIObjects{i},'VertexDeleted',@moved_roi);
                        end
                    else
                        if h.c(38).Value
                            
                            CurrentDrawMode = lower(Class{1}.Segment(i).Type);
                            switch CurrentDrawMode
                                case 'line'
                                    h.ROIObjects{i} = drawline('Position',Class{1}.Segment(i).ROIObject.Position,...
                                        'Deletable',0,...
                                        'InteractionsAllowed','none',...
                                        'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                        'Label',sprintf('%s || %s',Class{1}.Segment(i).Name,Class{1}.Segment(i).SubSegmentName),...
                                        'LabelAlpha',0.6,...
                                        'Color',SegmentColors{2*SegmentIndex-1},...
                                        'StripeColor',SegmentColors{2*SegmentIndex});
                                case 'freehand'
                                    h.ROIObjects{i} = drawfreehand('Position',Class{1}.Segment(i).ROIObject.Position,...
                                        'Deletable',0,...
                                        'InteractionsAllowed','none',...
                                        'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                        'Label',sprintf('%s || %s',Class{1}.Segment(i).Name,Class{1}.Segment(i).SubSegmentName),...
                                        'LabelAlpha',0.6,...
                                        'Color',SegmentColors{2*SegmentIndex-1},...
                                        'StripeColor',SegmentColors{2*SegmentIndex});
                                case 'square'
                                    h.ROIObjects{i} = drawsquare('Position',Class{1}.Segment(i).ROIObject.Position,...
                                        'Deletable',0,...
                                        'InteractionsAllowed','none',...
                                        'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                        'Label',sprintf('%s || %s',Class{1}.Segment(i).Name,Class{1}.Segment(i).SubSegmentName),...
                                        'LabelAlpha',0.6,...
                                        'Color',SegmentColors{2*SegmentIndex-1},...
                                        'StripeColor',SegmentColors{2*SegmentIndex});
                                case 'circle'
                                    h.ROIObjects{i} = drawcircle('Center',Class{1}.Segment(i).ROIObject.Center,...
                                        'Radius',Class{1}.Segment(i).ROIObject.Radius,...
                                        'Deletable',0,...
                                        'InteractionsAllowed','none',...
                                        'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                        'Label',sprintf('%s || %s',Class{1}.Segment(i).Name,Class{1}.Segment(i).SubSegmentName),...
                                        'LabelAlpha',0.6,...
                                        'Color',SegmentColors{2*SegmentIndex-1},...
                                        'StripeColor',SegmentColors{2*SegmentIndex});
                                case 'ellipse'
                                    h.ROIObjects{i} = drawellipse('Center',Class{1}.Segment(i).ROIObject.Center,...
                                        'SemiAxes',Class{1}.Segment(i).ROIObject.SemiAxes,...
                                        'RotationAngle',Class{1}.Segment(i).ROIObject.RotationAngle,...
                                        'Deletable',0,...
                                        'InteractionsAllowed','none',...
                                        'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                        'Label',sprintf('%s || %s',Class{1}.Segment(i).Name,Class{1}.Segment(i).SubSegmentName),...
                                        'LabelAlpha',0.6,...
                                        'Color',SegmentColors{2*SegmentIndex-1},...
                                        'StripeColor',SegmentColors{2*SegmentIndex});
                                case 'polygon'
                                    h.ROIObjects{i} = drawpolygon('Position',Class{1}.Segment(i).ROIObject.Position,...
                                        'Deletable',0,...
                                        'InteractionsAllowed','none',...
                                        'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                        'Label',sprintf('%s || %s',Class{1}.Segment(i).Name,Class{1}.Segment(i).SubSegmentName),...
                                        'LabelAlpha',0.6,...
                                        'Color',SegmentColors{2*SegmentIndex-1},...
                                        'StripeColor',SegmentColors{2*SegmentIndex});
                                case 'rectangle'
                                    h.ROIObjects{i} = drawrectangle('Position',Class{1}.Segment(i).ROIObject.Position,...
                                        'RotationAngle',Class{1}.Segment(i).ROIObject.RotationAngle,...
                                        'Rotatable',true,...
                                        'Deletable',0,...
                                        'InteractionsAllowed','none',...
                                        'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                        'Label',sprintf('%s || %s',Class{1}.Segment(i).Name,Class{1}.Segment(i).SubSegmentName),...
                                        'LabelAlpha',0.6,...
                                        'Color',SegmentColors{2*SegmentIndex-1},...
                                        'StripeColor',SegmentColors{2*SegmentIndex});
                                case 'crosshair'
                                    h.ROIObjects{i} = drawcrosshair('Position',Class{1}.Segment(i).ROIObject.Position,...
                                        'Deletable',0,...
                                        'InteractionsAllowed','none',...
                                        'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                        'Label',sprintf('%s || %s',Class{1}.Segment(i).Name,Class{1}.Segment(i).SubSegmentName),...
                                        'LabelAlpha',0.6,...
                                        'Color',SegmentColors{2*SegmentIndex-1},...
                                        'StripeColor',SegmentColors{2*SegmentIndex});
                                case 'point'
                                    h.ROIObjects{i} = drawpoint('Position',Class{1}.Segment(i).ROIObject.Position,...
                                        'Deletable',0,...
                                        'InteractionsAllowed','none',...
                                        'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                        'Label',sprintf('%s || %s',Class{1}.Segment(i).Name,Class{1}.Segment(i).SubSegmentName),...
                                        'LabelAlpha',0.6,...
                                        'Color',SegmentColors{2*SegmentIndex-1},...
                                        'StripeColor',SegmentColors{2*SegmentIndex});
                                case 'assisted'
                                    h.ROIObjects{i} = drawassisted('Position',Class{1}.Segment(i).ROIObject.Position,...
                                        'Deletable',0,...
                                        'InteractionsAllowed','none',...
                                        'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                        'Label',sprintf('%s || %s',Class{1}.Segment(i).Name,Class{1}.Segment(i).SubSegmentName),...
                                        'LabelAlpha',0.6,...
                                        'Color',SegmentColors{2*SegmentIndex-1},...
                                        'StripeColor',SegmentColors{2*SegmentIndex});
                                case 'polyline'
                                    h.ROIObjects{i} = drawpolyline('Position',Class{1}.Segment(i).ROIObject.Position,...
                                        'Deletable',0,...
                                        'InteractionsAllowed','none',...
                                        'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                        'Label',sprintf('%s || %s',Class{1}.Segment(i).Name,Class{1}.Segment(i).SubSegmentName),...
                                        'LabelAlpha',0.6,...
                                        'Color',SegmentColors{2*SegmentIndex-1},...
                                        'StripeColor',SegmentColors{2*SegmentIndex});
                            end
                        else
                            CurrentDrawMode = lower(Class{1}.Segment(i).Type);
                            switch CurrentDrawMode
                                case 'line'
                                    h.ROIObjects{i} = drawline('Position',Class{1}.Segment(i).ROIObject.Position,...
                                        'Deletable',0,...
                                        'InteractionsAllowed','none',...
                                        'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                        'Color',SegmentColors{2*SegmentIndex-1},...
                                        'StripeColor',SegmentColors{2*SegmentIndex});
                                case 'freehand'
                                    h.ROIObjects{i} = drawfreehand('Position',Class{1}.Segment(i).ROIObject.Position,...
                                        'Deletable',0,...
                                        'InteractionsAllowed','none',...
                                        'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                        'Color',SegmentColors{2*SegmentIndex-1},...
                                        'StripeColor',SegmentColors{2*SegmentIndex});
                                case 'square'
                                    h.ROIObjects{i} = drawsquare('Position',Class{1}.Segment(i).ROIObject.Position,...
                                        'Deletable',0,...
                                        'InteractionsAllowed','none',...
                                        'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                        'Color',SegmentColors{2*SegmentIndex-1},...
                                        'StripeColor',SegmentColors{2*SegmentIndex});
                                case 'circle'
                                    h.ROIObjects{i} = drawcircle('Center',Class{1}.Segment(i).ROIObject.Center,...
                                        'Radius',Class{1}.Segment(i).ROIObject.Radius,...
                                        'Deletable',0,...
                                        'InteractionsAllowed','none',...
                                        'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                        'Color',SegmentColors{2*SegmentIndex-1},...
                                        'StripeColor',SegmentColors{2*SegmentIndex});
                                case 'ellipse'
                                    h.ROIObjects{i} = drawellipse('Center',Class{1}.Segment(i).ROIObject.Center,...
                                        'SemiAxes',Class{1}.Segment(i).ROIObject.SemiAxes,...
                                        'RotationAngle',Class{1}.Segment(i).ROIObject.RotationAngle,...
                                        'Deletable',0,...
                                        'InteractionsAllowed','none',...
                                        'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                        'Color',SegmentColors{2*SegmentIndex-1},...
                                        'StripeColor',SegmentColors{2*SegmentIndex});
                                case 'polygon'
                                    h.ROIObjects{i} = drawpolygon('Position',Class{1}.Segment(i).ROIObject.Position,...
                                        'Deletable',0,...
                                        'InteractionsAllowed','none',...
                                        'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                        'Color',SegmentColors{2*SegmentIndex-1},...
                                        'StripeColor',SegmentColors{2*SegmentIndex});
                                case 'rectangle'
                                    h.ROIObjects{i} = drawrectangle('Position',Class{1}.Segment(i).ROIObject.Position,...
                                        'RotationAngle',Class{1}.Segment(i).ROIObject.RotationAngle,...
                                        'Rotatable',true,...
                                        'Deletable',0,...
                                        'InteractionsAllowed','none',...
                                        'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                        'Color',SegmentColors{2*SegmentIndex-1},...
                                        'StripeColor',SegmentColors{2*SegmentIndex});
                                case 'crosshair'
                                    h.ROIObjects{i} = drawcrosshair('Position',Class{1}.Segment(i).ROIObject.Position,...
                                        'Deletable',0,...
                                        'InteractionsAllowed','none',...
                                        'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                        'Color',SegmentColors{2*SegmentIndex-1},...
                                        'StripeColor',SegmentColors{2*SegmentIndex});
                                case 'point'
                                    h.ROIObjects{i} = drawpoint('Position',Class{1}.Segment(i).ROIObject.Position,...
                                        'Deletable',0,...
                                        'InteractionsAllowed','none',...
                                        'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                        'Color',SegmentColors{2*SegmentIndex-1},...
                                        'StripeColor',SegmentColors{2*SegmentIndex});
                                case 'assisted'
                                    h.ROIObjects{i} = drawassisted('Position',Class{1}.Segment(i).ROIObject.Position,...
                                        'Deletable',0,...
                                        'InteractionsAllowed','none',...
                                        'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                        'Color',SegmentColors{2*SegmentIndex-1},...
                                        'StripeColor',SegmentColors{2*SegmentIndex});
                                case 'polyline'
                                    h.ROIObjects{i} = drawpolyline('Position',Class{1}.Segment(i).ROIObject.Position,...
                                        'Deletable',0,...
                                        'InteractionsAllowed','none',...
                                        'LineWidth',Class{1}.Segment(i).ROIObject.LineWidth,...
                                        'Color',SegmentColors{2*SegmentIndex-1},...
                                        'StripeColor',SegmentColors{2*SegmentIndex});
                            end
                        end
                    end
                end
                
                if ~h.c(37).Value
                    create_map_of_rois();
                end
                
            end
            
            function create_map_of_rois(varargin)
                
                Size = size(h.I(1).CData);
                h.ROIImage = zeros(Size);
                
                for i=1:length(Class{1}.Segment)
                    if ~sum(strfind(Class{1}.Segment(i).Name,'Snapped'))
                        TempMask = h.ROIObjects{i}.createMask;
                        TempMask = imdilate(TempMask,strel('diamond',ceil(Size(1)/128)));
                        h.ROIImage(TempMask) = i;
                    end
                end
                
            end
            
            function moved_roi(varargin)
                
                    CurrentDrawMode = Class{1}.Segment(h.CurrentIndex).Type;
                    switch CurrentDrawMode
                        case 'line'
                        case 'freehand'
                        case 'circle'
                            Class{1}.Segment(h.CurrentIndex).ROIObject.Center = [];
                            Class{1}.Segment(h.CurrentIndex).ROIObject.Radius = [];
                            Class{1}.Segment(h.CurrentIndex).ROIObject.Center = h.ROIObjects{h.CurrentIndex}.Center;
                            Class{1}.Segment(h.CurrentIndex).ROIObject.Radius = h.ROIObjects{h.CurrentIndex}.Radius;
                        case 'ellipse'
                            Class{1}.Segment(h.CurrentIndex).ROIObject.Center = [];
                            Class{1}.Segment(h.CurrentIndex).ROIObject.SemiAxes = [];
                            Class{1}.Segment(h.CurrentIndex).ROIObject.RotationAngle = [];
                            Class{1}.Segment(h.CurrentIndex).ROIObject.Center = h.ROIObjects{h.CurrentIndex}.Center;
                            Class{1}.Segment(h.CurrentIndex).ROIObject.SemiAxes = h.ROIObjects{h.CurrentIndex}.SemiAxes;
                            Class{1}.Segment(h.CurrentIndex).ROIObject.RotationAngle = h.ROIObjects{h.CurrentIndex}.RotationAngle;
                        case 'polygon'
                        case 'rectangle'
                            Class{1}.Segment(h.CurrentIndex).ROIObject.RotationAngle = [];
                            Class{1}.Segment(h.CurrentIndex).ROIObject.RotationAngle = h.ROIObjects{h.CurrentIndex}.RotationAngle;
                        case 'crosshair'
                        case 'point'
                        case 'assisted'
                    end
                Class{1}.Segment(h.CurrentIndex).ROIObject.Position = [];
                Class{1}.Segment(h.CurrentIndex).ROIObject.Position = h.ROIObjects{h.CurrentIndex}.Position;
                
            end
            
            function show_labels(varargin)
                
                draw_channel_1
                
            end
            
            function create_new_subsegment(varargin)
                
                h.EnableEditing = 0;
                draw_channel_1
                
                Class{1}.Segment(end+1).Name = h.SegmentBox.String{h.SegmentBox.Value};
                
                k = 1;
                SubSegmentName = sprintf('SubS-%02d',k);
                while sum(strcmp({SubSegmentName},h.SubsegmentBox.String))
                    k = k + 1;
                    SubSegmentName = sprintf('SubS-%02d',k);
                end
                Class{1}.Segment(end).SubSegmentName = SubSegmentName;
                
                draw_segment('SubSegment')
                
                Class{1}.sort_segments_by_name_and_subsegmentname;
                
                h.EnableEditing = 1;
                draw_channel_1
                h.SubsegmentBox.Value = find(strcmp(h.SubsegmentBox.String,SubSegmentName));
                draw_channel_1
            end
            
            function create_new_segment(varargin)
                
                h.EnableEditing = 0;
                draw_channel_1
                
                k = 1;
                SegmentName = sprintf('Seg-%02d',k);
                while sum(strcmp({SegmentName},h.SegmentBox.String))
                    k = k + 1;
                    SegmentName = sprintf('Seg-%02d',k);
                end
                
                if isempty(Class{1}.Segment(1).Name)
                    Class{1}.Segment(1).Name = SegmentName;
                else
                    Class{1}.Segment(end+1).Name = SegmentName;
                end
                
                SubSegmentName = sprintf('SubS-%02d',1);
                Class{1}.Segment(end).SubSegmentName = SubSegmentName;
                
                draw_segment('Segment')
                
                Class{1}.sort_segments_by_name_and_subsegmentname;
                
                if length(Class{1}.Segment) == 1
                    set(h.SegmentBox,'String',{SegmentName})
                    set(h.SubsegmentBox,'String',{SubSegmentName})
                end
                
                h.EnableEditing = 1;
                draw_channel_1
                
                h.SegmentBox.Value = find(strcmp(h.SegmentBox.String,SegmentName));
                segment_changed
                h.SubsegmentBox.Value = find(strcmp(h.SubsegmentBox.String,SubSegmentName));
                draw_channel_1
            end
            
            function segment_changed(varargin)
                % set initial listbox items
                Class{1}.sort_segments_by_name_and_subsegmentname;
                Names = {Class{1}.Segment.Name};
                SubSegmentNames = {Class{1}.Segment.SubSegmentName};
                UniqueNames = unique(Names);
                set(h.SegmentBox,'String',UniqueNames);
                h.SubsegmentBox.Value = 1;
                set(h.SubsegmentBox,'String',{SubSegmentNames{strcmp({h.SegmentBox.String{h.SegmentBox.Value}},Names)}});
                draw_channel_1
            end
            
            function draw_segment(Mode)
                
                if isequal(Mode,'SubSegment')
                    Addendum = 0;
                elseif isequal(Mode,'Segment')
                    Addendum = 0;
                end
                
                if FiberLikeMenu ||...
                        (GeneralMenu &&...
                        isequal(lower(h.c(51).String{h.c(51).Value}),'polyline'))
                    ROIObject = drawpolyline;
                    Class{1}.Segment(end + Addendum).ROIObject.Position = ROIObject.Position;
                    Class{1}.Segment(end + Addendum).ROIObject.LineWidth = ROIObject.LineWidth;
                    Class{1}.Segment(end + Addendum).Type = 'polyline';
                elseif GeneralMenu
                    CurrentDrawMode = lower(h.c(51).String{h.c(51).Value});
                    switch CurrentDrawMode
                        case 'line'
                            ROIObject = drawline;
                        case 'freehand'
                            ROIObject = drawfreehand;
                        case 'square'
                            ROIObject = drawsq;
                        case 'circle'
                            ROIObject = drawcircle;
                            Class{1}.Segment(end + Addendum).ROIObject.Center = ROIObject.Center;
                            Class{1}.Segment(end + Addendum).ROIObject.Radius = ROIObject.Radius;
                        case 'ellipse'
                            ROIObject = drawellipse;
                            Class{1}.Segment(end + Addendum).ROIObject.Center = ROIObject.Center;
                            Class{1}.Segment(end + Addendum).ROIObject.SemiAxes = ROIObject.SemiAxes;
                            Class{1}.Segment(end + Addendum).ROIObject.RotationAngle = ROIObject.RotationAngle;
                        case 'polygon'
                            ROIObject = drawpolygon;
                        case 'rectangle'
                            ROIObject = drawrectangle('Rotatable',true);
                            Class{1}.Segment(end + Addendum).ROIObject.RotationAngle = ROIObject.RotationAngle;
                        case 'crosshair'
                            ROIObject = drawcrosshair;
                        case 'point'
                            ROIObject = drawpoint;
                        case 'assisted'
                            ROIObject = drawassisted;
                    end
                    Class{1}.Segment(end + Addendum).ROIObject.Position = ROIObject.Position;
                    Class{1}.Segment(end + Addendum).ROIObject.LineWidth = ROIObject.LineWidth;
                    Class{1}.Segment(end + Addendum).Type = CurrentDrawMode;
                end
                
            end
            
            function subsegment_changed(varargin)
                Class{1}.sort_segments_by_name_and_subsegmentname;
                draw_channel_1
            end
            
            function snap_segment(varargin)
                
                Class{1}.sort_segments_by_name_and_subsegmentname;
                Indizes = find(strcmp({Class{1}.Segment.Name},{Class{1}.Segment(h.CurrentIndex).Name}));
                
                Method = ['localreg' num2str(round(str2num(h.c(44).String)))];
                
                Class{1}.snap_line_segments_to_local_perpendicular_maximum(str2double(h.c(41).String),...
                    str2double(h.c(42).String),...
                    str2double(h.c(43).String),...
                    Method,Indizes)
                
                draw_channel_1
            end
            
            function snap_all(varargin)
                
                Class{1}.sort_segments_by_name_and_subsegmentname;
                Method = ['localreg' num2str(round(str2num(h.c(44).String)))];
                
                Class{1}.snap_line_segments_to_local_perpendicular_maximum(str2double(h.c(41).String),...
                    str2double(h.c(42).String),...
                    str2double(h.c(43).String),...
                    Method,[])
                
                draw_channel_1
            end
            
            function apply_segmentation_to_overlay_group(varargin)
                
                Class{1}.sort_segments_by_name_and_subsegmentname;
                obj.apply_segmentation_to_overlay_group(Class{1});
                
            end
            
            function delete_selected_segment(varargin)
                
                OldString = get(h.SegmentBox,'String');
                DeleteIdx = get(h.SegmentBox,'Value');
                
                Class{1}.Segment(strcmp({OldString{DeleteIdx}},{Class{1}.Segment.Name})) = [];
                
                if length(Class{1}.Segment) == 0
                    Class{1}.Segment = struct('Name',[],...
                            'Type',[],...
                            'SubSegmentName',[],...
                            'ROIObject',[]);
                end
                
                set(h.SegmentBox,'Value',1);
                Class{1}.sort_segments_by_name_and_subsegmentname;
                draw_channel_1
            end
            
            function delete_selected_subsegment(varargin)
                
                if length(h.SubsegmentBox.String) == 1
                    delete_selected_segment
                    return
                end
                
                OldString = get(h.SegmentBox,'String');
                DeleteIdx = get(h.SegmentBox,'Value');
                SubOldString = get(h.SubsegmentBox,'String');
                SubDeleteIdx = get(h.SubsegmentBox,'Value');
                
                Class{1}.Segment(strcmp({OldString{DeleteIdx}},{Class{1}.Segment.Name})&...
                    strcmp({SubOldString{SubDeleteIdx}},{Class{1}.Segment.SubSegmentName})) = [];
                
                if length(Class{1}.Segment) == 0
                    Class{1}.Segment = struct('Name',[],...
                            'Type',[],...
                            'SubSegmentName',[],...
                            'ROIObject',[]);
                end
                
                set(h.SubsegmentBox,'Value',length(h.SubsegmentBox.String)-1);
                Class{1}.sort_segments_by_name_and_subsegmentname;
                draw_channel_1
            end
            
            function add_tag(varargin)
                
                CurrentTagfieldContent = h.c(52).String;
                
                TagValid = Experiment.check_segment_tag_validity(CurrentTagfieldContent);
                
                if ~TagValid
                    warning('Tag is not valid. Refrain from using "*" in the string');
                    return
                end
                
                TagGrouping = h.c(55).String{h.c(55).Value};
                
                switch TagGrouping
                    case 'Current Segment'
                        FocusClass = obj.get_class_instance(ClassIndex(h.CurChannel1Idx,:));
                        KeySegmentName = h.SegmentBox.String{h.SegmentBox.Value};
                        [~,CurrentSegment] = AFMBaseClass.find_matching_segment(KeySegmentName,FocusClass.Segment);
                        FocusClass.add_tag_to_segment(CurrentTagfieldContent,CurrentSegment);
                    case 'Identical Segments in Ov. Group'
                        FocusClass = obj.get_class_instance(ClassIndex(h.CurChannel1Idx,:));
                        KeySegmentName = h.SegmentBox.String{h.SegmentBox.Value};
                        if ~FocusClass.OverlayGroup.hasOverlayGroup
                            [~,CurrentSegment] = AFMBaseClass.find_matching_segment(KeySegmentName,FocusClass.Segment);
                            FocusClass.add_tag_to_segment(CurrentTagfieldContent,CurrentSegment);
                        else
                            [~,CurrentSegment] = AFMBaseClass.find_matching_segment(KeySegmentName,FocusClass.Segment);
                            FocusClass.add_tag_to_segment(CurrentTagfieldContent,CurrentSegment);
                            for i=1:FocusClass.OverlayGroup.Size
                                if strcmp(FocusClass.OverlayGroup.Names{i},FocusClass.Name)
                                    continue
                                end
                                Index = find(strcmp({FocusClass.OverlayGroup.Names{i}},obj.ForceMapNames));
                                if ~isempty(Index)
                                    [~,CurrentSegment] = AFMBaseClass.find_matching_segment(KeySegmentName,obj.FM{Index}.Segment);
                                    obj.FM{Index}.add_tag_to_segment(CurrentTagfieldContent,CurrentSegment);
                                end
                                Index = find(strcmp({FocusClass.OverlayGroup.Names{i}},obj.AFMImageNames));
                                if ~isempty(Index)
                                    [~,CurrentSegment] = AFMBaseClass.find_matching_segment(KeySegmentName,obj.I{Index}.Segment);
                                    obj.I{Index}.add_tag_to_segment(CurrentTagfieldContent,CurrentSegment);
                                end
                            end
                        end
                    case 'Identical Segments in Experiment'
                        KeySegmentName = h.SegmentBox.String{h.SegmentBox.Value};
                        for i=1:obj.NumForceMaps
                            [~,CurrentSegment] = AFMBaseClass.find_matching_segment(KeySegmentName,obj.FM{i}.Segment);
                            obj.FM{i}.add_tag_to_segment(CurrentTagfieldContent,CurrentSegment);
                        end
                        for i = 1:obj.NumAFMImages
                            [~,CurrentSegment] = AFMBaseClass.find_matching_segment(KeySegmentName,obj.I{i}.Segment);
                            obj.I{i}.add_tag_to_segment(CurrentTagfieldContent,CurrentSegment);
                        end
                    case 'All Segments in image'
                        FocusClass = obj.get_class_instance(ClassIndex(h.CurChannel1Idx,:));
                        FocusClass.add_tag_to_segment(CurrentTagfieldContent);
                    case 'All Segments in Ov. Group'
                        FocusClass = obj.get_class_instance(ClassIndex(h.CurChannel1Idx,:));
                        if ~FocusClass.OverlayGroup.hasOverlayGroup
                            FocusClass.add_tag_to_segment(CurrentTagfieldContent);
                        else
                            FocusClass.add_tag_to_segment(CurrentTagfieldContent);
                            for i=1:FocusClass.OverlayGroup.Size
                                if strcmp(FocusClass.OverlayGroup.Names{i},FocusClass.Name)
                                    continue
                                end
                                Index = find(strcmp({FocusClass.OverlayGroup.Names{i}},obj.ForceMapNames));
                                if ~isempty(Index)
                                    obj.FM{Index}.add_tag_to_segment(CurrentTagfieldContent);
                                end
                                Index = find(strcmp({FocusClass.OverlayGroup.Names{i}},obj.AFMImageNames));
                                if ~isempty(Index)
                                    obj.I{Index}.add_tag_to_segment(CurrentTagfieldContent);
                                end
                            end
                        end
                    case 'All Segments in Experiment'
                        for i=1:obj.NumForceMaps
                            obj.FM{i}.add_tag_to_segment(CurrentTagfieldContent);
                        end
                        for i = 1:obj.NumAFMImages
                            obj.I{i}.add_tag_to_segment(CurrentTagfieldContent);
                        end
                end
                
                update_tag_box
            end
            
            function delete_tag(varargin)
                
                CurrentTagfieldContent = h.c(58).String{h.c(58).Value};
                
                TagValid = Experiment.check_segment_tag_validity(CurrentTagfieldContent);
                
                if ~TagValid
                    warning('Tag is not valid. Refrain from using "*" in the string');
                    return
                end
                
                TagGrouping = h.c(55).String{h.c(55).Value};
                
                switch TagGrouping
                    case 'Current Segment'
                        FocusClass = obj.get_class_instance(ClassIndex(h.CurChannel1Idx,:));
                        KeySegmentName = h.SegmentBox.String{h.SegmentBox.Value};
                        [~,CurrentSegment] = AFMBaseClass.find_matching_segment(KeySegmentName,FocusClass.Segment);
                        FocusClass.remove_tag_from_segment(CurrentTagfieldContent,CurrentSegment);
                    case 'Identical Segments in Ov. Group'
                        FocusClass = obj.get_class_instance(ClassIndex(h.CurChannel1Idx,:));
                        KeySegmentName = h.SegmentBox.String{h.SegmentBox.Value};
                        if ~FocusClass.OverlayGroup.hasOverlayGroup
                            [~,CurrentSegment] = AFMBaseClass.find_matching_segment(KeySegmentName,FocusClass.Segment);
                            FocusClass.remove_tag_from_segment(CurrentTagfieldContent,CurrentSegment);
                        else
                            [~,CurrentSegment] = AFMBaseClass.find_matching_segment(KeySegmentName,FocusClass.Segment);
                            FocusClass.remove_tag_from_segment(CurrentTagfieldContent,CurrentSegment);
                            for i=1:FocusClass.OverlayGroup.Size
                                if strcmp(FocusClass.OverlayGroup.Names{i},FocusClass.Name)
                                    continue
                                end
                                Index = find(strcmp({FocusClass.OverlayGroup.Names{i}},obj.ForceMapNames));
                                if ~isempty(Index)
                                    [~,CurrentSegment] = AFMBaseClass.find_matching_segment(KeySegmentName,obj.FM{Index}.Segment);
                                    obj.FM{Index}.remove_tag_from_segment(CurrentTagfieldContent,CurrentSegment);
                                end
                                Index = find(strcmp({FocusClass.OverlayGroup.Names{i}},obj.AFMImageNames));
                                if ~isempty(Index)
                                    [~,CurrentSegment] = AFMBaseClass.find_matching_segment(KeySegmentName,obj.I{Index}.Segment);
                                    obj.I{Index}.remove_tag_from_segment(CurrentTagfieldContent,CurrentSegment);
                                end
                            end
                        end
                    case 'Identical Segments in Experiment'
                        KeySegmentName = h.SegmentBox.String{h.SegmentBox.Value};
                        for i=1:obj.NumForceMaps
                            [~,CurrentSegment] = AFMBaseClass.find_matching_segment(KeySegmentName,obj.FM{i}.Segment);
                            obj.FM{i}.remove_tag_from_segment(CurrentTagfieldContent,CurrentSegment);
                        end
                        for i=1:obj.NumAFMImages
                            [~,CurrentSegment] = AFMBaseClass.find_matching_segment(KeySegmentName,obj.I{i}.Segment);
                            obj.I{i}.remove_tag_from_segment(CurrentTagfieldContent,CurrentSegment);
                        end
                    case 'All Segments in image'
                        FocusClass = obj.get_class_instance(ClassIndex(h.CurChannel1Idx,:));
                        FocusClass.remove_tag_from_segment(CurrentTagfieldContent);
                    case 'All Segments in Ov. Group'
                        FocusClass = obj.get_class_instance(ClassIndex(h.CurChannel1Idx,:));
                        if ~FocusClass.OverlayGroup.hasOverlayGroup
                            FocusClass.remove_tag_from_segment(CurrentTagfieldContent);
                        else
                            FocusClass.remove_tag_from_segment(CurrentTagfieldContent);
                            for i=1:FocusClass.OverlayGroup.Size
                                if strcmp(FocusClass.OverlayGroup.Names{i},FocusClass.Name)
                                    continue
                                end
                                Index = find(strcmp({FocusClass.OverlayGroup.Names{i}},obj.ForceMapNames));
                                if ~isempty(Index)
                                    obj.FM{Index}.remove_tag_from_segment(CurrentTagfieldContent);
                                end
                                Index = find(strcmp({FocusClass.OverlayGroup.Names{i}},obj.AFMImageNames));
                                if ~isempty(Index)
                                    obj.I{Index}.remove_tag_from_segment(CurrentTagfieldContent);
                                end
                            end
                        end
                    case 'All Segments in Experiment'
                        for i=1:obj.NumForceMaps
                            obj.FM{i}.remove_tag_from_segment(CurrentTagfieldContent);
                        end
                        for i=1:obj.NumAFMImages
                            obj.I{i}.remove_tag_from_segment(CurrentTagfieldContent);
                        end
                end
                
                update_tag_box
                
            end
            
            function update_tag_box()
                
                if isempty(h.SegmentBox.String)
                    return
                end
                
                FocusClass = obj.get_class_instance(ClassIndex(h.CurChannel1Idx,:));
                KeySegmentName = h.SegmentBox.String{h.SegmentBox.Value};
                [~,CurrentSegment] = AFMBaseClass.find_matching_segment(KeySegmentName,FocusClass.Segment);
                TagList = FocusClass.read_out_tag_list(CurrentSegment);
                
                OldValue = h.c(58).Value;
                if numel(TagList) < OldValue
                    OldValue = 1;
                end
                h.c(58).String = TagList;
                h.c(58).Value = OldValue;
                
            end
            
            uiwait(h.Fig)
        end
        
        function [Fig,DataMat] = visualize_listed_data(obj,Property,YAxisName,BaseUnit,Method,ErrorBars,UseGrouping,ListOfIndizes)
            % function [Fig, DataMat] = visualize_listed_data(obj,
            % Property, YAxisName, BaseUnit, Method, ErrorBars,
            % UseGrouping, ListOfIndizes)
            %
            % This function visualizes the specified properties of ForceMap
            % objects in a variety of graphical formats, such as boxplots,
            % bar plots, or median bar plots. It can handle grouping and
            % different error bar types, such as standard deviation,
            % standard error, or confidence interval.
            %
            % Input:
            %
            % obj: An object representing the AFM experiment content. It
            % should have properties for force maps (FM) and their
            % corresponding names (ForceMapNames).
            %
            % Property: A string indicating the property to visualize,
            % e.g., 'EModOliverPharr' or 'IndentationDepth'.
            %
            % YAxisName: A string representing the label for the y-axis.
            %
            % BaseUnit: A string representing the base unit for the
            % property, e.g., 'N/m' for modulus.
            %
            % Method: A string indicating the visualization method. Options
            % are 'Boxplot', 'Barplot', or 'BarMedian'.
            %
            % ErrorBars: A string specifying the type of error bars to
            % display. Options are 'std' (standard deviation), 'ste'
            % (standard error), or 'ci' (confidence interval).
            %
            % UseGrouping: A boolean indicating whether to use grouping
            % (true) or not (false).
            %
            % ListOfIndizes: An array of indices specifying which ForceMaps
            % to visualize. Default is all.
            %
            % Output:
            %
            % Fig: A figure handle for the generated plot.
            %
            % DataMat: A matrix containing the data used for visualization,
            % with each column corresponding to a ForceMap instance and
            % each row corresponding to a data point.
            %
            % Usage:
            %
            % [Fig, DataMat] = visualize_listed_data(obj, Property,
            % YAxisName, BaseUnit, Method, ErrorBars, UseGrouping,
            % ListOfIndizes)
            %
            % Example:
            %
            % To visualize the 'EModOliverPharr' property using a boxplot
            % with standard deviation error bars and grouping, call the
            % function as: [Fig, DataMat] = visualize_listed_data(obj,
            % 'EModOliverPharr', 'Elastic Modulus', 'N/m', 'Boxplot',
            % 'std', true)
            
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
                    if length(Data(k).Values) == 1
                        Range = Data(k).Values;
                    else
                        Range = range(Data(k).Values);
                    end
                    [Data(k).Multiplier,Data(k).Unit] = AFMImage.parse_unit_scale(Range,BaseUnit,5);
                    Data(k).Name = obj.GroupFM(i).Name;
                    Data(k).Length = length(Data(k).Values);
                    if isequal(ErrorBars,'std')
                        Data(k).Bars = std(Data(k).Values,'omitnan');
                    elseif isequal(ErrorBars,'ste')
                        Data(k).Bars = std(Data(k).Values,'omitnan')/sqrt(length(Data(k).Values(~isnan(Data(k).Values))));
                    elseif isequal(ErrorBars,'ci')
                        [~,~,Data(k).Bars,~] = ttest(Data(k).Values(~isnan(Data(k).Values)),mean(Data(k).Values,'omitnan'));
                    end
                    k = k + 1;
                end
            else
                for i=ListOfIndizes
                    Data(k).Values = obj.FM{i}.get(Property);
                    if length(Data(k).Values) == 1
                        Range = Data(k).Values;
                    else
                        Range = range(Data(k).Values);
                    end
                    [Data(k).Multiplier,Data(k).Unit] = AFMImage.parse_unit_scale(Range,BaseUnit,5);
                    Data(k).Name = obj.FM{i}.Name;
                    Data(k).Length = length(Data(k).Values);
                    if isequal(ErrorBars,'std')
                        Data(k).Bars = std(Data(k).Values,'omitnan');
                    elseif isequal(ErrorBars,'ste')
                        Data(k).Bars = std(Data(k).Values,'omitnan')/sqrt(length(Data(k).Values(~isnan(Data(k).Values))));
                    elseif isequal(ErrorBars,'ci')
                        [~,~,Data(k).Bars,~] = ttest(Data(k).Values(~isnan(Data(k).Values)),mean(Data(k).Values,'omitnan'));
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
                bar(mean(DataMat,1,'omitnan'));
                hold on
                for i=1:length(Data)
                    errorbar(i,mean(DataMat(:,i),'omitnan'),-Data(i).Bars(1),Data(i).Bars(end),'Color','k','LineWidth',1.5);
                end
                title(sprintf('Barplot of %s',YAxisName))
            elseif isequal(Method,'BarMedian')
                bar(median(DataMat,1,'omitnan'));
                hold on
                for i=1:length(Data)
                    errorbar(i,median(DataMat(:,i),'omitnan'),-Data(i).Bars(1),Data(i).Bars(end),'Color','k','LineWidth',1.5);
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
%                 DataOP(i,:) = obj.FM{i}.EModOliverPharr(obj.FM{i}.RectApexIndex);
                DataHS(i,:) = obj.FM{i}.EModHertz(obj.FM{i}.RectApexIndex);
%                 OutliersOP = isoutlier(DataOP(i,:));
                OutliersHS = isoutlier(DataHS(i,:));
                for j=1:length(obj.FM{i}.RectApexIndex)
%                     if DataOP(i,j) > (nanmedian(DataOP(i,:))+2.5*iqr(DataOP(i,:))) || ...
%                             DataOP(i,j) < (nanmedian(DataOP(i,:))-2.5*iqr(DataOP(i,:))) || ...
%                             obj.FM{i}.ExclMask(obj.FM{i}.List2Map(obj.FM{i}.RectApexIndex(j),1),obj.FM{i}.List2Map(obj.FM{i}.RectApexIndex(j),2)) == 0 ||...
%                             OutliersOP(j) == 1
%                         DataOP(i,j) = NaN;
%                     elseif DataOP(i,j) < 0
%                         DataOP(i,j) = NaN;
%                     end
                    if DataHS(i,j) > (median(DataHS(i,:),'omitnan')+2.5*iqr(DataHS(i,:))) || ...
                            DataHS(i,j) < (median(DataHS(i,:),'omitnan')-2.5*iqr(DataHS(i,:))) || ...
                            obj.FM{i}.ExclMask(obj.FM{i}.List2Map(obj.FM{i}.RectApexIndex(j),1),obj.FM{i}.List2Map(obj.FM{i}.RectApexIndex(j),2)) == 0 ||...
                            OutliersHS(j) == 1
                        DataHS(i,j) = NaN;
                    elseif DataHS(i,j) < 0
                        DataHS(i,j) = NaN;
                    end
                end
            end
            
%             DataMeansOP = nanmean(DataOP,2);
            DataMeansHS = mean(DataHS,2,'omitnan');
            
            for i=1:obj.NumForceMaps
%                 obj.FM{i}.FibrilEModOliverPharr = DataMeansOP(i);
                obj.FM{i}.FibrilEModHertz = DataMeansHS(i);
            end
%             
%             figure('Name','OliverPharr vs HertzSneddon','Color','w');
%             plot(DataMeansHS,DataMeansOP,'bO')
%             legend(sprintf('E-Mod Hertz vs. Oliver-Pharr (R=%.3f)',corr(DataMeansHS,DataMeansOP)))
% %             xlim([0,N+1])
%             xlabel('E-Mod Hertz [Pa]')
%             ylabel('E-Mod Oliver-Pharr [Pa]')
            
            % loop over all rows of the Test Matrix, doing paired ttests
            for i=1:size(TestMat,1)
%                 % Statistics for Oliver-Pharr Method
%                 [hOP(i),pOP(i)] = ...
%                     ttest(DataMeansOP(obj.GroupFM(TestMat(i,2)).Indices),...
%                     DataMeansOP(obj.GroupFM(TestMat(i,1)).Indices),'Tail','right');
%                 
%                 figure('Name','Paired Right Tailed T-Test','Units','normalized','Position',[0.2 0.2 0.5 0.5],'Color','w')
%                 boxplot([DataMeansOP(obj.GroupFM(TestMat(i,1)).Indices) DataMeansOP(obj.GroupFM(TestMat(i,2)).Indices)]*1e-6)
%                 title('Paired Right Tailed T-Test for Oliver-Pharr Method')
%                 xticklabels({obj.GroupFM(TestMat(i,1)).Name,obj.GroupFM(TestMat(i,2)).Name})
%                 xlabel('Test Group')
%                 ylabel('Indentation Modulus [MPa]')
%                 DeltaMean = mean(DataMeansOP(obj.GroupFM(TestMat(i,2)).Indices)) - mean(DataMeansOP(obj.GroupFM(TestMat(i,1)).Indices));
%                 Sigma = std(DataMeansOP(obj.GroupFM(TestMat(i,2)).Indices) - DataMeansOP(obj.GroupFM(TestMat(i,1)).Indices));
%                 Beta = sampsizepwr('t',[0 Sigma],DeltaMean,[],length(obj.GroupFM(TestMat(i,2)).Indices));
%                 Stats = {sprintf('\\DeltaMean = %.2fMPa',DeltaMean*1e-6),...
%                     sprintf('P-Value = %.4f%',pOP(i)),...
%                     sprintf('Power \\beta = %.2f%%',Beta*100),...
%                     sprintf('Number of Specimen = %i',length(obj.GroupFM(TestMat(i,2)).Indices))};
%                 text(0.5,0.8,Stats,...
%                     'Units','normalized',...
%                     'FontSize',12,...
%                     'HorizontalAlignment','center')
                
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
%             DiffControlOP = DataMeansOP(11:20) - DataMeansOP(1:10);
%             DiffMGOOP = DataMeansOP(31:40) - DataMeansOP(21:30);
            DiffControlHS = DataMeansHS(11:20) - DataMeansHS(1:10);
            DiffMGOHS = DataMeansHS(31:40) - DataMeansHS(21:30);
            
%             % For Oliver-Pharr
%             [hOP,pOP,ciOP,statsOP] = ttest2(DiffMGOOP,DiffControlOP);
%             figure('Name','Two Sample T-Test','Units','normalized','Position',[0.2 0.2 0.5 0.5],'Color','w')
%             yyaxis left
%             boxplot([DiffControlOP DiffMGOOP])
%             ax = gca;
%             YLim = ax.YLim;
%             ylabel('Difference Before-After E-Mod [Pa]')
%             DeltaMean = mean(DiffMGOOP) - mean(DiffControlOP);
%             PooledSTD = statsOP.sd;
%             yyaxis right
%             errorbar(1.5,DeltaMean,ciOP(2)-DeltaMean,'O');
%             ylim(YLim)
%             xticks([1 1.5 2])
%             title('Two Sample T-Test for E-Mod Oliver-Pharr Method')
%             ax = gca;
%             ax.TickLabelInterpreter = 'tex';
%             xticklabels({sprintf('%s - %s',obj.GroupFM(2).Name,obj.GroupFM(1).Name),...
%                 '\DeltaMean with CI',...
%                 sprintf('%s - %s',obj.GroupFM(4).Name,obj.GroupFM(3).Name)})
%             ylabel('Difference of Differences [Pa]')
%             Beta = sampsizepwr('t2',[mean(DiffControlOP) PooledSTD],mean(DiffMGOOP),[],length(DiffControlOP),'Ratio',length(DiffMGOOP)/length(DiffControlOP));
%             Stats = {sprintf('\\DeltaMean = %.2f MPa',DeltaMean*1e-6),...
%                 sprintf('P-Value = %.4f%',pOP),...
%                 sprintf('Power \\beta = %.2f%%',Beta*100),...
%                 sprintf('Degrees of freedom df = %i',statsOP.df)};
%             text(0.5,0.8,Stats,...
%                 'Units','normalized',...
%                 'FontSize',12,...
%                 'HorizontalAlignment','center')
            
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
        
        function reset_segmentations(obj)
            % function reset_segmentations(obj)
            %
            % This function resets the segmentations for all AFMImage (I)
            % and ForceMap (FM) instances within the given object. It
            % clears the 'Name', 'Type', 'SubSegmentName', and 'ROIObject'
            % fields for each instance's 'Segment' structure.
            %
            % Input:
            %
            % obj: An object representing the AFM experiment content. It
            % should have properties for AFM images (I), force maps (FM),
            % and their corresponding counts (NumAFMImages and
            % NumForceMaps).
            %
            % Output: None
            %
            % Usage:
            %
            % reset_segmentations(obj)
            %
            % Example: To reset all segmentations for the given object,
            % simply call the function as: reset_segmentations(obj)

            for i=1:obj.NumAFMImages
                obj.I{i}.Segment = struct('Name',[],...
                            'Type',[],...
                            'SubSegmentName',[],...
                            'ROIObject',[]);
            end
            for i=1:obj.NumForceMaps
                obj.FM{i}.Segment = struct('Name',[],...
                            'Type',[],...
                            'SubSegmentName',[],...
                            'ROIObject',[]);
            end
            
        end
        
        function RadiusNM = calculate_tip_radius(obj,TipDepthNM,TipIndex)
            % function RadiusNM = calculate_tip_radius(obj, TipDepthNM,
            % TipIndex)
            %
            % This function calculates the tip radius of the AFM cantilever
            % using a specified tip depth and index. It iteratively fits
            % spheres to the tip surface and averages the radii of the
            % spheres to determine the tip radius.
            %
            % Input:
            %
            % obj (object): An object representing the AFM experiment
            % content. It should have a property
            %
            % 'CantileverTips' containing information about cantilever
            % tips.
            %
            % TipDepthNM (double, optional): The depth in nanometers used
            % for the tip surface fitting process. Default value is 20 nm.
            %
            % TipIndex (integer): Index of the cantilever tip to use for
            % the calculations.
            %
            % Output:
            %
            % RadiusNM (double): The calculated tip radius in nanometers.
            %
            % Usage:
            %
            % calculate_tip_radius(obj, TipDepthNM, TipIndex)
            %
            % Example: To calculate the tip radius for a specific
            % cantilever tip in the object, call the function as:
            % tip_radius = calculate_tip_radius(obj, 20, 1)
            
            
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
            % function check_for_new_host(obj)
            %
            % This function checks if the system environment has changed
            % (e.g., when moving from one computer to another). If the
            % environment has changed, it adjusts the object properties
            % accordingly.
            %
            % Input:
            %
            % obj (object): An object containing properties related to the
            % system environment, such as
            %
            % 'HostOS' (operating system) and 'HostName' (name of the host
            % computer).
            %
            % Output:
            %
            % None. This function modifies the input object properties in
            % place.
            %
            % Usage:
            %
            % check_for_new_host(obj)
            %
            % Example:
            %
            % To check for changes in the system environment and update the
            % object properties accordingly, call:
            % check_for_new_host(myObject)
            %
            % Notes: This function supports the following operating
            % systems: Windows ('PCW'), Linux ('GLN'), and macOS ('MAC').
            
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
        
        function reference_slope_parser(obj,DefaultOption,SkipGUIBool)
            % function reference_slope_parser(obj, DefaultOption,
            % SkipGUIBool)
            %
            % This function parses the method to determine the reference
            % slope for force maps and applies the selected method to the
            % object. It provides several methods for setting reference
            % slopes, such as using a user-defined value, user input, or
            % automatically calculating the reference slope from a specific
            % area.
            %
            % Input:
            %
            % obj (object): An object containing the force map properties
            % and reference slope methods.
            %
            % DefaultOption (integer, optional): The default reference
            % slope method index. Default is 1.
            %
            % SkipGUIBool (boolean, optional): A flag to skip the GUI and
            % use the default option. Default is false.
            %
            % Output:
            %
            % None. This function modifies the input object properties in
            % place.
            %
            % Usage:
            %
            % reference_slope_parser(obj) reference_slope_parser(obj,
            % DefaultOption) reference_slope_parser(obj, DefaultOption,
            % SkipGUIBool)
            %
            % Example:
            %
            % To parse and apply the reference slope method to an object
            % using default options, call: reference_slope_parser(myObject)
            %
            % Notes:
            %
            % The reference slope methods supported by this function are:
            %
            % 1. SetAllToValue: Set all reference slopes to a specific
            % value.
            %
            % 2. UserInput: Manually input the reference slopes for each
            % force map.
            %
            % 3. FromRefFM: Use a reference force map to calculate the
            % reference slopes.
            %
            % 4. FromArea: Calculate the reference slopes from a specific
            % area.
            %
            % 5. AutomaticFibril: Automatically calculate the reference
            % slopes based on fibril detection.
            %
            % 6. Automatic: Automatically calculate the reference slopes
            % without fibril detection.
            %
            % 7. KeepOriginal: Keep the original reference slopes.


            if nargin < 2
                DefaultOption = 1;
                SkipGUIBool = false;
            elseif nargin < 3
                SkipGUIBool = false;
            end
            
            Methods = false(7,1);
            if ~SkipGUIBool
                Methods(DefaultOption) = true;
                [ChosenMethod, AppRetSwitch] = obj.reference_slope_parser_gui(Methods);
                Methods = false(6,1);
                Methods(ChosenMethod) = true;
                obj.ReferenceSlopeFlag.AppRetSwitch = AppRetSwitch;
            end
            
            obj.ReferenceSlopeFlag.SetAllToValue = Methods(1);
            obj.ReferenceSlopeFlag.UserInput = Methods(2);
            obj.ReferenceSlopeFlag.FromRefFM = Methods(3);
            obj.ReferenceSlopeFlag.FromArea = Methods(4);
            obj.ReferenceSlopeFlag.AutomaticFibril = Methods(5);
            obj.ReferenceSlopeFlag.Automatic = Methods(6);
            obj.ReferenceSlopeFlag.KeepOriginal = Methods(7);
            obj.ReferenceSlopeFlag.Adaptive = false;
            
            if SkipGUIBool
                obj.ReferenceSlopeFlag.(...
                    obj.ForceMapAnalysisOptions.SensitivityCorrection.SensitivityCorrectionMethod) = true;
                obj.ReferenceSlopeFlag.AppRetSwitch = ...
                    obj.ForceMapAnalysisOptions.SensitivityCorrection.AppRetSwitch;
            end
            
            % Process all the RefSlope-Methods that can be done frontloaded
            if obj.ReferenceSlopeFlag.SetAllToValue
                GetValue = inputdlg('To which value should the reference slopes be set?','Choose a value',1);
                Value = str2double(GetValue{1});
                for i=1:obj.NumForceMaps
                    obj.FM{i}.set_reference_slope_to_value(Value)
                end
            elseif obj.ReferenceSlopeFlag.KeepOriginal
                Value = 1;
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
                    obj.FM{i}.RefSlopeMask = Mask;
                end
            end
        end
        
        function reference_slope_calculator(obj,Index)
            % function reference_slope_calculator(obj, Index)
            %
            % This function calculates the reference slope for a force map
            % using the selected method, such as from a reference force
            % map, automatically with or without fibril detection, or from
            % a specific area. The calculated reference slope is then
            % applied to the force map object.
            %
            % Input:
            %
            % obj (object): An object containing the force map properties
            % and reference slope methods.
            %
            % Index (integer): The index of the force map within the object
            % for which the reference slope should be calculated.
            %
            % Output:
            %
            % None. This function modifies the input object properties in
            % place.
            %
            % Usage:
            %
            % reference_slope_calculator(obj, Index)
            %
            % Example:
            %
            % To calculate and apply the reference slope for the first
            % force map in an object, call:
            % reference_slope_calculator(myObject, 1)
            %
            % Notes:
            %
            % The reference slope methods supported by this function are:
            %
            % 1. FromRefFM: Use a reference force map to calculate the
            % reference slopes.
            %
            % 2. AutomaticFibril: Automatically calculate the reference
            % slopes based on fibril detection.
            %
            % 3. Automatic: Automatically calculate the reference slopes
            % without fibril detection.
            %
            % 4. Adaptive: Calculate the reference slopes adaptively based
            % on an automatically created background mask.
            %
            % 5. FromArea: Calculate the reference slopes from a specific
            % area.
            %
            % In addition to the reference slope methods, this function
            % also accepts various options as part of the input object's
            % ReferenceSlopeFlag.Options property, such as AppRetSwitch,
            % FitRangeMode, FitRangeLowerFraction, FitRangeUpperFraction,
            % FitRangeLowerValue, FitRangeUpperValue, and MovingWindowSize.
            % These options control the behavior of the reference slope
            % calculation process.
            
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
                        obj.RefFM{RefIdx}.estimate_cp_old
                        Mask = ones(obj.RefFM{RefIdx}.NumPixelsX,obj.RefFM{RefIdx}.NumPixelsY);
                        obj.RefFM{RefIdx}.calculate_adaptive_sensitivity_from_area('Mask',Mask,...
                            'AppRetSwitch',obj.ReferenceSlopeFlag.Options.AppRetSwitch,...
                            'FitRangeMode',obj.ReferenceSlopeFlag.Options.FitRangeMode,...
                            'FitRangeLowerFraction',obj.ReferenceSlopeFlag.Options.FitRangeLowerFraction,...
                            'FitRangeUpperFraction',obj.ReferenceSlopeFlag.Options.FitRangeUpperFraction,...
                            'FitRangeLowerValue',obj.ReferenceSlopeFlag.Options.FitRangeLowerValue,...
                            'FitRangeUpperValue',obj.ReferenceSlopeFlag.Options.FitRangeUpperValue,...
                            'MovingWindowSize',obj.RefFM{RefIdx}.NCurves,...
                            'Verbose',false);
                    end
                    obj.FM{i}.RefSlope = obj.RefFM{RefIdx}.RefSlope;
                    obj.FM{i}.HasRefSlope = true;
                elseif length(obj.RefFM) == 1
                    if ~obj.RefFM{1}.HasRefSlope
                        if ~obj.RefFM{1}.BaseAndTiltFlag
                            obj.RefFM{1}.base_and_tilt
                        end
                        obj.RefFM{1}.estimate_cp_old
                        Mask = ones(obj.RefFM{1}.NumPixelsX,obj.RefFM{1}.NumPixelsY);
                        obj.RefFM{1}.calculate_adaptive_sensitivity_from_area('Mask',Mask,...
                            'AppRetSwitch',obj.ReferenceSlopeFlag.Options.AppRetSwitch,...
                            'FitRangeMode',obj.ReferenceSlopeFlag.Options.FitRangeMode,...
                            'FitRangeLowerFraction',obj.ReferenceSlopeFlag.Options.FitRangeLowerFraction,...
                            'FitRangeUpperFraction',obj.ReferenceSlopeFlag.Options.FitRangeUpperFraction,...
                            'FitRangeLowerValue',obj.ReferenceSlopeFlag.Options.FitRangeLowerValue,...
                            'FitRangeUpperValue',obj.ReferenceSlopeFlag.Options.FitRangeUpperValue,...
                            'MovingWindowSize',obj.RefFM{1}.NCurves,...
                            'Verbose',false);
                    end
                    obj.FM{i}.RefSlope = obj.RefFM{1}.RefSlope;
                    obj.FM{i}.HasRefSlope = true;
                end
            elseif obj.ReferenceSlopeFlag.AutomaticFibril
                Mask = obj.FM{i}.BackgroundMask;
                obj.FM{i}.calculate_adaptive_sensitivity_from_area('Mask',Mask,...
                    'AppRetSwitch',obj.ReferenceSlopeFlag.Options.AppRetSwitch,...
                    'FitRangeMode',obj.ReferenceSlopeFlag.Options.FitRangeMode,...
                    'FitRangeLowerFraction',obj.ReferenceSlopeFlag.Options.FitRangeLowerFraction,...
                    'FitRangeUpperFraction',obj.ReferenceSlopeFlag.Options.FitRangeUpperFraction,...
                    'FitRangeLowerValue',obj.ReferenceSlopeFlag.Options.FitRangeLowerValue,...
                    'FitRangeUpperValue',obj.ReferenceSlopeFlag.Options.FitRangeUpperValue,...
                    'MovingWindowSize',obj.FM{i}.NCurves,...
                    'Verbose',false);
            elseif obj.ReferenceSlopeFlag.Automatic
                obj.FM{i}.create_automatic_background_mask(.8)
                Mask = obj.FM{i}.BackgroundMask;
                obj.FM{i}.calculate_adaptive_sensitivity_from_area('Mask',Mask,...
                    'AppRetSwitch',obj.ReferenceSlopeFlag.Options.AppRetSwitch,...
                    'FitRangeMode',obj.ReferenceSlopeFlag.Options.FitRangeMode,...
                    'FitRangeLowerFraction',obj.ReferenceSlopeFlag.Options.FitRangeLowerFraction,...
                    'FitRangeUpperFraction',obj.ReferenceSlopeFlag.Options.FitRangeUpperFraction,...
                    'FitRangeLowerValue',obj.ReferenceSlopeFlag.Options.FitRangeLowerValue,...
                    'FitRangeUpperValue',obj.ReferenceSlopeFlag.Options.FitRangeUpperValue,...
                    'MovingWindowSize',obj.FM{i}.NCurves,...
                    'Verbose',false);
            elseif obj.ReferenceSlopeFlag.Adaptive
                obj.FM{i}.create_automatic_background_mask(.8)
                Mask = obj.FM{i}.BackgroundMask;
                obj.FM{i}.calculate_adaptive_sensitivity_from_area('Mask',Mask,...
                    'AppRetSwitch',obj.ReferenceSlopeFlag.Options.AppRetSwitch,...
                    'FitRangeMode',obj.ReferenceSlopeFlag.Options.FitRangeMode,...
                    'FitRangeLowerFraction',obj.ReferenceSlopeFlag.Options.FitRangeLowerFraction,...
                    'FitRangeUpperFraction',obj.ReferenceSlopeFlag.Options.FitRangeUpperFraction,...
                    'FitRangeLowerValue',obj.ReferenceSlopeFlag.Options.FitRangeLowerValue,...
                    'FitRangeUpperValue',obj.ReferenceSlopeFlag.Options.FitRangeUpperValue,...
                    'MovingWindowSize',obj.ReferenceSlopeFlag.Options.MovingWindowSize,...
                    'Verbose',false);
            elseif obj.ReferenceSlopeFlag.FromArea
                obj.FM{i}.calculate_adaptive_sensitivity_from_area('Mask',obj.FM{i}.RefSlopeMask,...
                    'AppRetSwitch',obj.ReferenceSlopeFlag.Options.AppRetSwitch,...
                    'FitRangeMode',obj.ReferenceSlopeFlag.Options.FitRangeMode,...
                    'FitRangeLowerFraction',obj.ReferenceSlopeFlag.Options.FitRangeLowerFraction,...
                    'FitRangeUpperFraction',obj.ReferenceSlopeFlag.Options.FitRangeUpperFraction,...
                    'FitRangeLowerValue',obj.ReferenceSlopeFlag.Options.FitRangeLowerValue,...
                    'FitRangeUpperValue',obj.ReferenceSlopeFlag.Options.FitRangeUpperValue,...
                    'MovingWindowSize',obj.FM{i}.NCurves,...
                    'Verbose',false);
            end
        end
        
        function assign_reference_force_map(obj,DefaultValues)
%             function assign_reference_force_map(obj, DefaultValues)
%             
%             This function assigns reference force maps to the force maps
%             in the object. The user provides a list of force map indices
%             that correspond to each reference force map. The function
%             ensures that each force map is assigned exactly one reference
%             force map.
%             
%             Input:
%             
%             obj (object): An object containing the force maps and their
%             properties.
%             
%             DefaultValues (cell array, optional): A cell array containing
%             default values for force map assignments. Each cell contains
%             a string with force map indices, e.g., '1 2 3 4 8 9 10' or
%             '1:4 8:10'. If not provided, the default value for each
%             reference force map will be set as 'e.g. 1 2 3 4 8 9 10 or
%             1:4 8:10'.
%             
%             Output:
%             
%             None. This function modifies the input object properties in
%             place.
%             
%             Usage:
%             
%             assign_reference_force_map(obj, DefaultValues)
%             
%             Example:
%             
%             To assign reference force maps to the force maps in an object
%             using default values, call:
%             assign_reference_force_map(myObject, DefaultValues)
%             
%             Notes:
%             
%             1. The function ensures that each force map is assigned
%             exactly one reference force map. If the assignment is
%             incorrect, a warning will be shown, and the user will be
%             prompted to reassign.
%             
%             2. A table with force map names and their corresponding
%             numbers is shown in the background while the user is prompted
%             to assign the reference force maps. This helps the user make
%             the correct assignments.
            
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
%             function assign_cantilever_tips(obj, DefaultValues)
%             
%             This function assigns cantilever tips to the force maps in
%             the object. The user provides a list of force map indices
%             that correspond to each cantilever tip. The function ensures
%             that each force map is assigned exactly one cantilever tip.
%             
%             Input:
%             
%             obj (object): An object containing the force maps and their
%             properties.
%             
%             DefaultValues (cell array, optional): A cell array containing
%             default values for force map assignments. Each cell contains
%             a string with force map indices, e.g., '1 2 3 4 8 9 10' or
%             '1:4 8:10'. If not provided, the default value for each
%             cantilever tip will be set as 'e.g. 1 2 3 4 8 9 10 or 1:4
%             8:10'.
%             
%             Output:
%             
%             None. This function modifies the input object properties in
%             place.
%             
%             Usage:
%             
%             assign_cantilever_tips(obj, DefaultValues)
%             
%             Example:
%             
%             To assign cantilever tips to the force maps in an object
%             using default values, call: assign_cantilever_tips(myObject,
%             DefaultValues)
%             
%             Notes:
%             
%             1. The function ensures that each force map is assigned
%             exactly one cantilever tip. If the assignment is incorrect, a
%             warning will be shown, and the user will be prompted to
%             reassign.
%             
%             2. A table with force map names and their corresponding
%             numbers is shown in the background while the user is prompted
%             to assign the cantilever tips. This helps the user make the
%             correct assignments.
            
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
            T = table(reshape(Names,[],1));
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
%             function load_in_reference_maps(obj)
%             
%             This function loads a specified number of reference force
%             maps into the object from .jpk-force-map or .jpk-qi-data
%             files. The user is prompted to select the reference force
%             maps to be loaded. The function ensures that the correct
%             number of reference force maps is loaded.
%             
%             Input:
%             
%             obj (object): An object containing the force maps and their
%             properties.
%             
%             Output:
%             
%             None. This function modifies the input object properties in
%             place by adding the loaded reference force maps.
%             
%             Usage:
%             
%             load_in_reference_maps(obj)
%             
%             Example:
%             
%             To load reference force maps into an object, call:
%             load_in_reference_maps(myObject)
%             
%             Notes:
%             
%             1. The user is prompted to enter the number of reference
%             force maps to be loaded. The function will ensure that the
%             correct number of reference force maps is loaded.
%             
%             2. The user can select one or more reference force maps at a
%             time. The function will keep prompting the user to select
%             reference force maps until the specified number is loaded.
%             
%             3. Reference force maps are loaded from .jpk-force-map or
%             .jpk-qi-data files.
%             
%             4. Each loaded reference force map is assigned an ID in the
%             format 'ReferenceMap-X', where X is the index of the
%             reference force map.
            
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
%             function initialize_flags(obj)
%             
%             This function initializes a set of flags for the input object
%             to track the status of various analyses and assignments. It
%             sets flags for fibril analysis, force map analysis,
%             preprocessing, grouping, cantilever tips, reference maps, and
%             reference slopes, among others.
%             
%             Input:
%             
%             obj (object): An object containing the force maps and their
%             properties.
%             
%             Output:
%             
%             None. This function modifies the input object properties in
%             place by initializing the flags.
%             
%             Usage:
%             
%             initialize_flags(obj)
%             
%             Example:
%             
%             To initialize the flags for an object, call:
%             initialize_flags(myObject)
%             
%             Notes:
%             
%             1. The function initializes flags for multiple aspects,
%             including: - Fibril analysis - Force map analysis -
%             Preprocessing - Grouping - Surface potential maps - SMFS
%             (single-molecule force spectroscopy) - Cantilever tips -
%             Reference maps - Reference slopes
%             
%             2. The function checks if the object already has flags and
%             updates them accordingly.
%             
%             3. The flags are used to track the progress and status of
%             various analyses and assignments within the object.
            
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
                obj.SMFSFlag.Min= false(NFM,1);
                obj.SMFSFlag.Length= false(NFM,1);
                obj.SMFSFlag.Fit= false(NFM,1);
                
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
                obj.ReferenceSlopeFlag.Adaptive = false;
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
                obj.ReferenceSlopeFlag.Adaptive = false;
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
%             function [PopUp, ClassIndex] =
%             string_of_existing_class_instances(obj)
%             
%             This function generates a list of strings and corresponding
%             class indices for existing class instances of AFM images,
%             force maps, reference force maps, and cantilever tip images
%             in the given object. The output strings are used for
%             populating a drop-down menu (or pop-up menu) in the user
%             interface.
%             
%             Input:
%             
%             obj (object): An object containing AFM images, force maps,
%             reference force maps, and cantilever tip images.
%             
%             Output:
%             
%             PopUp (cell array of strings): A list of strings representing
%             existing class instances in the input object. Each string is
%             formatted as "<Type>: <Instance Name>%", where <Type> is the
%             class type (Image, Force Map, Ref. Force Map, or Cant. Tip
%             Image), and <Instance Name> is the name of the instance.
%             
%             ClassIndex (Nx2 double array): A corresponding list of class
%             indices for the existing class instances, where the first
%             column represents the class type (1: Image, 2: Force Map, 3:
%             Ref. Force Map, 4: Cant. Tip Image), and the second column
%             represents the index of the instance within its respective
%             class.
%             
%             Usage:
%             
%             [PopUp, ClassIndex] = string_of_existing_class_instances(obj)
%             
%             Example:
%             
%             To generate a list of strings and corresponding class indices
%             for an object, call:
%             
%             [PopUp, ClassIndex] =
%             string_of_existing_class_instances(myObject)
%             
%             Notes:
%             
%             1. The function iterates through the object's properties (AFM
%             images, force maps, reference force maps, and cantilever tip
%             images) and generates the corresponding strings and class
%             indices.
%             
%             2. The output PopUp and ClassIndex can be used for creating a
%             drop-down menu in a user interface to allow the user to
%             select an existing class instance for further analysis or
%             visualization.
            
            k = 1;
            for i=1:obj.NumAFMImages
                PopUp{k} = sprintf('Image %i: %s%',i,obj.I{i}.Name);
                ClassIndex(k,:) = [1 i];
                k = k + 1;
            end
            for i=1:obj.NumForceMaps
                PopUp{k} = sprintf('Force Map %i: %s%',i,obj.FM{i}.Name);
                ClassIndex(k,:) = [2 i];
                k = k + 1;
            end
            for i=1:obj.NumReferenceForceMaps
                PopUp{k} = sprintf('Ref. Force Map %i: %s%',i,obj.RefFM{i}.Name);
                ClassIndex(k,:) = [3 i];
                k = k + 1;
            end
            for i=1:obj.NumCantileverTips
                PopUp{k} = sprintf('Cant. Tip Image %i: %s%',i,obj.CantileverTips{i}.Name);
                ClassIndex(k,:) = [4 i];
                k = k + 1;
            end
        end
        
        function Class = get_class_instance(obj,ClassIndex)
%             function Class = get_class_instance(obj, ClassIndex)
%             
%             This function retrieves a specific class instance from the
%             given object based on the provided class index. The class
%             index contains information about the class type (AFM images,
%             force maps, reference force maps, or cantilever tip images)
%             and the instance index within that class type.
%             
%             Input:
%             
%             obj (object): An object containing AFM images, force maps,
%             reference force maps, and cantilever tip images.
%             
%             ClassIndex (1x2 double array): The class index, where the
%             first element represents the class type (1: Image, 2: Force
%             Map, 3: Ref. Force Map, 4: Cant. Tip Image), and the second
%             element represents the index of the instance within its
%             respective class.
%             
%             Output:
%             
%             Class (object): The specific class instance retrieved from
%             the input object based on the provided class index.
%             
%             Usage:
%             
%             Class = get_class_instance(obj, ClassIndex)
%             
%             Example:
%             
%             To retrieve a specific class instance from an object based on
%             a class index, call: Class = get_class_instance(myObject, [2,
%             5])
%             
%             Notes:
%             
%             1. The function checks the class type specified in the
%             ClassIndex and retrieves the corresponding class instance
%             from the input object.
%             
%             2. This function is useful for retrieving a specific class
%             instance for further analysis or visualization when a user
%             selects an instance from a drop-down menu in a user
%             interface.
            
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
%             function cast_volume_data_to_single_precision(obj)
%             
%             This function casts the volume data of all force maps in the
%             given object to single precision. The volume data includes
%             approach and retraction curves, as well as height-height, and
%             based approach and retraction curves. The casting to single
%             precision reduces memory usage and improves performance,
%             especially when working with large datasets.
%             
%             Input:
%             
%             obj (object): An object containing force maps with volume
%             data. The volume data is originally stored as double
%             precision.
%             
%             Output:
%             
%             None. This function modifies the input object's force maps
%             volume data in-place by converting it to single precision.
%             
%             Usage:
%             
%             cast_volume_data_to_single_precision(obj)
%             
%             Example:
%             
%             To cast the volume data of all force maps in an object to
%             single precision, call:
%             cast_volume_data_to_single_precision(myObject)
%             
%             Notes:
%             
%             1. This function iterates over all force maps in the object,
%             and for each force map, it iterates over all curves
%             (approach, retraction, height-height, based approach, and
%             based retraction curves).
%             
%             2. It converts the curve data to single precision using the
%             'single()' function.
%             
%             3. The THApp and THRet fields of the force maps are set to
%             empty arrays, as they are no longer needed after casting the
%             data to single precision.
            
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
        
        function reset_selected_and_corrupted_flags(obj)
            % function reset_selected_and_corrupted_flags(obj)
            %
            % This function resets the selected and corrupted flags of all force maps within the given object.
            % These flags are used to identify specific force maps during analysis and processing, and resetting
            % them allows for a fresh start or re-analysis of the data.
            %
            % Input:
            % obj (object): An object containing force maps, each of which has selected and corrupted flags.
            %
            % Output:
            % None. This function modifies the input object's force maps' flags in-place, resetting the
            % selected and corrupted flags.
            %
            % Usage:
            % reset_selected_and_corrupted_flags(obj)
            %
            % Example:
            % To reset the selected and corrupted flags of all force maps in an object, call:
            % reset_selected_and_corrupted_flags(myObject)
            %
            % Notes:
            % 1. This function iterates over all force maps in the object and calls the
            % 'reset_selected_and_corrupted_flags()' method for each force map. This method is responsible
            % for resetting the flags within each individual force map.

            for i=1:obj.NumForceMaps
                obj.FM{i}.reset_selected_and_corrupted_flags()
            end
        end
        
        function load_python_files_to_memory(obj,IndexVector,RefIndexVector)
            % function load_python_files_to_memory(obj, IndexVector, RefIndexVector)
            %
            % This function loads the zipped force maps and reference force maps files into memory using
            % Python-based extraction. The function allows selective loading of force maps and reference force
            % maps through the use of index vectors.
            %
            % Input:
            % obj (object): An object containing force maps and reference force maps.
            % IndexVector (vector, optional): A vector containing the indices of the force maps to be loaded.
            % Default is to load all force maps.
            % RefIndexVector (vector, optional): A vector containing the indices of the reference force maps
            % to be loaded. Default is to load all reference force maps.
            %
            % Output:
            % None. This function modifies the input object in-place, loading the selected force maps and
            % reference force maps into memory.
            %
            % Usage:
            % load_python_files_to_memory(obj, IndexVector, RefIndexVector)
            %
            % Example:
            % To load all force maps and reference force maps into memory, call:
            % load_python_files_to_memory(myObject)
            %
            % To load selected force maps and reference force maps, call:
            % load_python_files_to_memory(myObject, [1, 3, 5], [2, 4])
            %
            % Notes:
            % 1. This function iterates over the specified force maps and reference force maps and calls the
            % 'load_zipped_files_with_python()' method for each. This method is responsible for extracting
            % and loading the data from the zipped files using Python.

            if nargin < 2
                IndexVector = 1:obj.NumForceMaps;
            end
            for i=IndexVector
                obj.FM{i}.load_zipped_files_with_python;
            end
            
            if nargin < 3
                RefIndexVector = 1:obj.NumReferenceForceMaps;
            end
            for i=RefIndexVector
                obj.RefFM{i}.load_zipped_files_with_python;
            end
        end
        
        function clear_python_files_from_memory(obj,IndexVector,RefIndexVector)
            % function clear_python_files_from_memory(obj, IndexVector, RefIndexVector)
            %
            % This function clears the force maps and reference force maps data from memory, which were
            % previously loaded using Python-based extraction. The function allows selective clearing of force
            % maps and reference force maps through the use of index vectors.
            %
            % Input:
            % obj (object): An object containing force maps and reference force maps.
            % IndexVector (vector, optional): A vector containing the indices of the force maps to be cleared
            % from memory. Default is to clear all force maps.
            % RefIndexVector (vector, optional): A vector containing the indices of the reference force maps
            % to be cleared from memory. Default is to clear all reference
            % force maps.
            %
            % Output:
            % None. This function modifies the input object in-place, clearing the selected force maps and
            % reference force maps from memory.
            %
            % Usage:
            % clear_python_files_from_memory(obj, IndexVector, RefIndexVector)
            %
            % Example:
            % To clear all force maps and reference force maps from memory, call:
            % clear_python_files_from_memory(myObject)
            %
            % To clear selected force maps and reference force maps, call:
            % clear_python_files_from_memory(myObject, [1, 3, 5], [2, 4])
            %
            % Notes:
            % 1. This function iterates over the specified force maps and reference force maps and calls the
            % 'clear_zipped_files_from_memory()' method for each. This method is responsible for clearing
            % the data loaded from the zipped files using Python.
            
            if nargin < 2
                IndexVector = 1:obj.NumForceMaps;
            end
            for i=IndexVector
                obj.FM{i}.clear_zipped_files_from_memory;
            end
            
            if nargin < 3
                RefIndexVector = 1:obj.NumReferenceForceMaps;
            end
            for i=RefIndexVector
                obj.RefFM{i}.clear_zipped_files_from_memory;
            end
        end
        
        function automatic_segmentation_on_singular_vertical_fiber_batch(obj,IndexVectorFM,IndexVectorImage)
            % function automatic_segmentation_on_singular_vertical_fiber_batch(obj, IndexVectorFM, IndexVectorImage)
            %
            % This function applies automatic segmentation on singular vertical fibers for a batch of force
            % maps and AFM images. It allows the user to specify the indices of the force maps and images to be
            % processed or to process all force maps and images by default.
            %
            % Input:
            % obj (object): An object containing force maps and AFM images.
            % IndexVectorFM (vector, optional): A vector containing the indices of the force maps to be
            % processed for automatic segmentation. Default is to process
            % all force maps.
            % IndexVectorImage (vector, optional): A vector containing the indices of the AFM images to be
            % processed for automatic segmentation. Default is to process
            % all AFM images if IndexVectorFM is not provided.
            %
            % Output:
            % None. This function modifies the input object in-place, applying automatic segmentation on the
            % specified force maps and AFM images.
            %
            % Usage:
            % automatic_segmentation_on_singular_vertical_fiber_batch(obj, IndexVectorFM, IndexVectorImage)
            %
            % Example:
            % To apply automatic segmentation on all force maps and AFM images, call:
            % automatic_segmentation_on_singular_vertical_fiber_batch(myObject)
            %
            % To apply automatic segmentation on selected force maps and AFM images, call:
            % automatic_segmentation_on_singular_vertical_fiber_batch(myObject, [1, 3, 5], [2, 4])
            %
            % Notes:
            % 1. This function iterates over the specified force maps and AFM images and calls the
            % 'automatic_segmentation_on_singular
            
            if nargin < 2
                for i=1:obj.NumForceMaps
                    obj.FM{i}.automatic_segmentation_on_singular_vertical_fiber
                end
                for i=1:obj.NumAFMImages
                    obj.I{i}.automatic_segmentation_on_singular_vertical_fiber
                end
            elseif nargin < 3
                for i=IndexVectorFM
                    obj.FM{i}.automatic_segmentation_on_singular_vertical_fiber
                end
            else
                for i=IndexVectorFM
                    obj.FM{i}.automatic_segmentation_on_singular_vertical_fiber
                end
                for i=IndexVectorImage
                    obj.I{i}.automatic_segmentation_on_singular_vertical_fiber
                end
            end
        end
        
        function OutClass = get_afm_base_class_by_name(obj,Name,AllCellStructs)
            % function OutClass = get_afm_base_class_by_name(obj, Name, AllCellStructs)
            %
            % This function searches for an AFM base class object (force map, AFM image, reference force map, or
            % cantilever tip) within the provided object based on its name. The user can specify whether to
            % search in all cell structures or only in force maps and AFM images.
            %
            % Input:
            % obj (object): An object containing force maps, AFM images, reference force maps, and cantilever
            % tips.
            % Name (string): The name of the AFM base class object to search for.
            % AllCellStructs (logical, optional): A flag indicating whether to search for the AFM base class
            % object in all cell structures (force maps, AFM images,
            % reference force maps, and cantilever tips) or only in force
            % maps and AFM images. Default is false (search only in force
            % maps and AFM images).
            %
            % Output:
            % OutClass (object): The AFM base class object found with the specified name. If no object with
            % the given name is found, the function returns an empty output.
            %
            % Usage:
            % OutClass = get_afm_base_class_by_name(obj, Name, AllCellStructs)
            %
            % Example:
            % To search for an AFM base class object named 'Sample1' in force maps and AFM images, call:
            % OutClass = get_afm_base_class_by_name(myObject, 'Sample1')
            %
            % To search for an AFM base class object named 'Sample1' in all cell structures, call:
            % OutClass

            if nargin < 3
                AllCellStructs = false;
            end
            
            for i=1:obj.NumForceMaps
                if isequal(Name,obj.FM{i}.Name)
                    OutClass = obj.FM{i};
                    return
                end
            end
            for i=1:obj.NumAFMImages
                if isequal(Name,obj.I{i}.Name)
                    OutClass = obj.I{i};
                    return
                end
            end
            if AllCellStructs
                for i=1:obj.NumReferenceForceMaps
                    if isequal(Name,obj.RefFM{i}.Name)
                        OutClass = obj.RefFM{i};
                        return
                    end
                end
                for i=1:obj.NumCantileverTips
                    if isequal(Name,obj.CantileverTips{i}.Name)
                        OutClass = obj.CantileverTips{i};
                        return
                    end
                end
            end
            
        end
        
        function set_force_map_analysis_options(obj)
            % function set_force_map_analysis_options(obj)
            %
            % This function sets the force map analysis options for the provided object. If no options are
            % specified, it sets the default force map analysis options. The function then prompts the user to
            % adjust the options using a user interface.
            %
            % Input:
            % obj (object): An object for which the force map analysis options need to be set.
            %
            % Output:
            % None. This function modifies the input object in-place, updating the force map analysis options.
            %
            % Usage:
            % set_force_map_analysis_options(obj)
            %
            % Example:
            % To set the force map analysis options for a given object, call:
            % set_force_map_analysis_options(myObject)
            %
            % Notes:
            % 1. The function first checks if the object's ForceMapAnalysisOptions field is empty. If it is
            % empty, the function sets the default force map analysis options using the
            % set_default_fma_options function.
            % 2. The function then calls ui_set_struct_fields to prompt the user to adjust the options using
            % a user interface.
            
            if isempty(obj.ForceMapAnalysisOptions)
                obj.ForceMapAnalysisOptions = set_default_fma_options;
            end
            
            obj.ForceMapAnalysisOptions = ui_set_struct_fields(obj.ForceMapAnalysisOptions);
            
        end
        
        function set_gramm_options(obj)
            % This function sets the gramm plot options for the Experiment object. The
            % function first checks whether the gramm options exist and are valid, and
            % if not, sets them to the default options using the set_default_gramm_options
            % function. The function then prompts the user to set the options using a
            % user interface and stores the updated options back in the object.
            
            % Input:
            % - obj: the Experiment object for which to set the gramm options
            %
            % Output: None. The function modifies the Experiment object by setting the
            % gramm options.
            
            % Valid data types:
            % - obj: Experiment object
            %
            % Valid values:
            % - None
            
            % Example:
            % 1. Set gramm options for an Experiment object:
            % exp = Experiment();
            % exp.set_gramm_options();
            
            if isempty(obj.GrammOptions) || ...
                    ~check_structs_for_same_fields_recursively(...
                    obj.GrammOptions,Experiment.set_default_gramm_options)
                obj.GrammOptions = Experiment.set_default_gramm_options;
            end
            
            obj.GrammOptions = ui_set_struct_fields(obj.GrammOptions);
            
        end
        
        function map_fiber_segment_properties_to_image_pixels(obj,PoolingMethod)
            % The function 'map_fiber_segment_properties_to_image_pixels'
            % is a method of the Experiment class.
            % It maps fiber segment properties such as position, angle, and length to
            % image pixels. The function takes two input arguments: 'obj' is the
            % Experiment object, and 'PoolingMethod' is an optional argument that
            % specifies the pooling method to use when mapping fiber segment properties
            % to image pixels. If 'PoolingMethod' is not specified, the default method
            % is used. The valid data types for 'obj' are Experiment objects, and the
            % valid data types for 'PoolingMethod' are strings or character arrays. The
            % possible values for 'PoolingMethod' are 'mean', 'median', 'min', 'max',
            % and 'none'. The default value is 'median'. The function loops over all
            % force maps and AFM images in the Experiment object and calls the
            % 'map_fiber_segment_properties_to_image_pixels' method for each of them
            % with the specified pooling method or the default one. The function does
            % not return any output arguments.
            
            for i=1:obj.NumForceMaps
                obj.FM{i}.map_fiber_segment_properties_to_image_pixels(PoolingMethod);
%                 obj.FM{i}.map_segments_to_image_pixels;
            end
            for i=1:obj.NumAFMImages
                obj.I{i}.map_fiber_segment_properties_to_image_pixels(PoolingMethod);
%                 obj.I{i}.map_segments_to_image_pixels;
            end
            
        end
        
        function map_segments_to_image_pixels(obj)
            % map_segments_to_image_pixels(obj)
            %
            % This function maps the fiber segments from force maps and AFM images to
            % the corresponding image pixels. The function applies the same mapping
            % transformation for each segment in the given image to transform the
            % segment into the image pixel space. The mapping can be affected by the
            % different pixel sizes in different AFM images. The function loops through
            % all force maps and AFM images in the Experiment object and maps their
            % fiber segments to image pixels.
            %
            % Inputs:
            % - obj: An Experiment object containing force maps and AFM images.
            %
            % Outputs: None. However, the function modifies the force maps and AFM
            % images in the Experiment object by mapping their fiber segments to image
            % pixels.
            
            for i=1:obj.NumForceMaps
                obj.FM{i}.map_segments_to_image_pixels;
            end
            for i=1:obj.NumAFMImages
                obj.I{i}.map_segments_to_image_pixels;
            end
            
        end
        
        function [DynPropNames,ChannelNames] = write_unrolled_channels_to_dynamic_properties(obj,ShowWaitbar)
            % The function write_unrolled_channels_to_dynamic_properties is a wrapper
            % function that calls the AFMBaseClass method of the same name. This method
            % writes unrolled channels to temporary properties of the AFMBaseClass
            % objects. The function can be used to extract dynamic property names and
            % channel names that are subsequently used for plotting or further
            % analysis. The function loops through all the Force Maps, AFM Images, and
            % CantileverTips stored in the object and calls the write_unrolled_channels_to_dynamic_properties
            % method for each of them. The outputs of the method, which are the dynamic
            % property names and channel names, are concatenated and returned as the
            % output of the wrapper function.
            %
            % INPUTS:
            % obj - an instance of the Experiment class.
            %
            % ShowWaitbar - a logical variable indicating whether or not to show the
            % waitbar during the loop through the object's force maps, AFM images, and
            % cantilever tips. Default value is true.
            %
            % OUTPUTS:
            % DynPropNames - a cell array containing the names of all the dynamic
            % properties across all the force maps, AFM images, and cantilever tips in
            % the Experiment class.
            %
            % ChannelNames - a cell array containing the names of all the channels
            % across all the force maps, AFM images, and cantilever tips in the
            % Experiment class. The channel names can be used to extract specific
            % channels for plotting or further analysis.
            
            if nargin < 3
                ShowWaitbar = true;
            end
            
            ExProps = {'FM','I','CantileverTips'};
            LoopNumber = [obj.NumForceMaps obj.NumAFMImages obj.NumCantileverTips];
            LoopSum = sum(LoopNumber);
            if ShowWaitbar
                h = waitbar(0,'Writing Channels to temporary properties...');
            end
            k = 1;
            for i=1:3
                if LoopNumber(i) == 0
                    continue
                end
                for j=1:LoopNumber(i)
                if ShowWaitbar
                    waitbar(k/LoopSum,h,sprintf('Writing Channels in Class %i/%i',k,LoopSum))
                end
                    [DynPropNames,ChannelNames] = obj.(ExProps{i}){j}.write_unrolled_channels_to_dynamic_properties;
                    k = k + 1;
                end
            end
            
            ChannelNames = unique(ChannelNames);
            
            if ShowWaitbar
                close(h)
            end
        end
        
        function delete_all_dynamic_properties(obj)
            % Wrapper function for AFMBaseClass Method with same name
            % The delete_all_dynamic_properties function is a wrapper function for the
            % corresponding method in the AFMBaseClass superclass. This function deletes
            % all the temporary dynamic properties created during analysis from all
            % ForceMap, AFMImage, and CantileverTip objects in the Experiment
            % instance.
            %
            % Inputs:
            % - obj: an instance of the Experiment class
            %
            % Outputs: none
            %
            % Note: this function will issue a warning to notify the user that all
            % temporary dynamic class properties have been deleted.
            
            ExProps = {'FM','I','CantileverTips'};
            LoopNumber = [obj.NumForceMaps obj.NumAFMImages obj.NumCantileverTips];
            LoopSum = sum(LoopNumber);
            for i=1:3
                if LoopNumber(i) == 0
                    continue
                end
                for j=1:LoopNumber(i)
                    obj.(ExProps{i}){j}.delete_all_dynamic_properties;
                end
            end
            
            warning('All temporary dynamic class properties deleted')
        end
        
        function broadcast_current_fma_id(obj,FMA_ID)
            
            obj.CurrentFMA_ID = FMA_ID;
            
            for i=1:obj.NumForceMaps
                obj.FM{i}.CurrentFMA_ID = FMA_ID;
            end
            for i=1:obj.NumAFMImages
                obj.I{i}.CurrentFMA_ID = FMA_ID;
            end
            for i=1:obj.NumReferenceForceMaps
                obj.RefFM{i}.CurrentFMA_ID = FMA_ID;
            end
            for i=1:obj.NumCantileverTips
                obj.CantileverTips{i}.CurrentFMA_ID = FMA_ID;
            end
            
        end
        
        function broadcast_CNN_optimization(obj,Index)
            
            OptimizedBatchSize = obj.FM{Index}.MiniBatchSize;
            
            for i=1:obj.NumForceMaps
                obj.FM{i}.MiniBatchSize = OptimizedBatchSize;
                obj.FM{i}.CPFlag.CNNopt = true;
                disp('Broadcasting optimal MiniBatchSize to other ForceMaps')
            end
            
        end
        
    end
    methods(Static)
        % Static auxilary methods mainly for tip deconvolution (code by Orestis Andriotis)
        % and Experiment() loading
        
        function wdata = getinfo(filename)
            %             %getinfo.m reads header information and extracts
            %             wave data from an .ibw file. The function takes
            %             in a file path/filename as input and returns the
            %             wave data as a 3D array. The data is read in
            %             chunks and preallocated for speed. The data type
            %             is also determined based on the type number in
            %             the header information, and the appropriate
            %             format string is selected for use with the fread
            %             function. If the file type is unknown, the
            %             function displays an error message. The output
            %             variable wdata is a 3D array where the first two
            %             dimensions are the pixel dimensions of the image
            %             and the third dimension corresponds to the number
            %             of images in the file.
            %
            % Input:
            %
            % filename: a string containing the file path and name of the
            % .ibw file to be read. Output:
            %
            % wdata: a 3D array containing the wave data from the .ibw
            % file. The first two dimensions represent the pixel dimensions
            % of the image, and the third dimension corresponds to the
            % number of images in the file. Valid data types:
            %
            % filename: a string containing the file path and name of the
            % .ibw file to be read. The file must be a valid .ibw file,
            % otherwise the function will not be able to read the header
            % information and extract the wave data. Possible values:
            %
            % filename: a string containing a valid file path and name of
            % an .ibw file on the system. type: integer, 2, 3, or 4,
            % indicating the type of data contained in the file. x, y:
            % integer values representing the dimensions of the image in
            % pixels. n: integer value indicating the number of images
            % contained in the file. d: a string indicating the format
            % string to be used with the fread function, based on the type
            % number in the header information. wdata: a 3D array
            % containing the wave data from the .ibw file.

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
            % THIS FUNCTION IS DEPRECATED
            %             This function is a deprecated graphical user
            %             interface (GUI) function that provides the user
            %                 with a list of options to choose from to
            %                 determine how to obtain reference slopes
            %                 (RefSlopes) from atomic force microscopy (AFM)
            %                 data. RefSlopes are used in the data analysis
            %                 process of AFM images to determine the height of
            %                 features on the sample surface. The options for
            %                 obtaining RefSlopes are:
            %
            % Set all RefSlopes to a chosen value Get RefSlopes from user
            % input for each Force Map individually Get RefSlopes from
            % Reference Force Maps Get RefSlopes from chosen area on Force
            % Map Get RefSlope from automatically identified glass
            % background around fibril Get RefSlope from automatically
            % identified glass background around object This function takes
            % the user input and returns the chosen option as a scalar
            % integer variable Out. The function also has an optional
            % return variable AppRetSwitch which returns 0 or 1, depending
            % on which radio button was selected by the user to specify
            % whether to obtain the RefSlope from the approach or retract
            % curve.
            %
            % Input:
            %
            % Methods: A scalar or vector of integers indicating which
            % options should be pre-selected in the GUI. Output:
            %
            % Out: An integer scalar indicating which option the user chose
            % to obtain RefSlopes from AFM data. AppRetSwitch: A scalar
            % integer equal to 0 or 1, depending on which radio button was
            % selected by the user to specify whether to obtain the
            % RefSlope from the approach or retract curve. Note that this
            % output variable is optional, and may not be returned if the
            % user did not select one of the options for obtaining
            % RefSlopes from the approach or retract curve.
            
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
        
        function [FileTypes,OutStruct,IsValid,BigData,PythonLoaderFlag,KeepPythonFilesOpen] = constructor_user_input_parser(ExperimentName,OS)
            
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
                'position',[.3 .025 .3 .085],'string','Big Data mode',...
                'tooltip','Recommended when processing >~100000 force curves',...
                'Value',true,'FontSize',18,'Callback',@pushed_big_data);
            c(17) = uicontrol(h.f,'style','checkbox','units','normalized',...
                'position',[.3 .025 .3 .025],'string','Python Loader Mode',...
                'tooltip','Recommended Mode: Avoids unpacking large folder structures',...
                'Value',true,'FontSize',18,'Callback',@pushed_python_loader,'enable','on');
            c(18) = uicontrol(h.f,'style','checkbox','units','normalized',...
                'position',[.6 .025 .3 .025],'string','Keep Files open in RAM',...
                'tooltip','Trades off burst execution speed for RAM. Recommended when running long analysis scripts. Not Recommended for data/results review',...
                'Enable','on',...
                'Value',false,...
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
%                 h.DragNDrop = uicontrol(h.f,'style','text','units','normalized',...
%                     'position',[.5 .025 .5 .05],...
%                     'string',sprintf('Drag and Drop files into boxes\n (excl. to Windows and Linux) or load from browser'),...
%                     'FontSize',16);
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
            
            function pushed_big_data(varargin)
                
                if c(16).Value
                    set(c(17),'Enable','on');
                else
                    set(c(17),'Enable','off');
                end
                
            end
            
            function pushed_python_loader(varargin)
                
                if c(17).Value
                    set(c(18),'Enable','on');
                else
                    set(c(18),'Enable','off');
                end
                
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
                PythonLoaderFlag = c(17).Value;
                KeepPythonFilesOpen = c(18).Value;
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
                try
                    if TempTempFile == 0
                        return
                    end
                catch
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
                OutStruct(Index).FullFile;
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
        
        function NewPath = switch_old_with_new_toplvl_path(OldPath,OldToplvl,NewToplvl)
            
            if isempty(OldPath)
                NewPath = OldPath;
                return
            end
            
            NewPath = replace(OldPath,OldToplvl,NewToplvl);
            
        end
        
        function FMAOptions = set_default_fma_options()
            
            
            
            CorrectedDZSlopesOptions = struct(...
                'SensitivityCorrectionMethod','corrected',...
                'AllowedSensitivityCorrectionMethod',{{'corrected',...
                'original','adaptive'}},...
                'FitRangeMode','SameAsSensitivityCorrection',...
                'AllowedFitRangeMode',{{'SameAsSensitivityCorrection','ValueHorizontal',...
                'ValueVertical','FractionHorizontal','FractionVertical'}},...
                'FitRangeLowerValue',0,...
                'FitRangeUpperValue',4e-9,...
                'FitRangeLowerFraction',.75,...
                'FitRangeUpperFraction',1,...
                'AppRetSwitch',false,...
                'DTypeAppRetSwitch','logical'...
                );
            
            OptionalMetrics = struct(...
                'AdhesionEnergy',true,...
                'DTypeAdhesionEnergy','logical',...
                'MaximumAdhesionForce',true,...
                'DTypeMaximumAdhesionForce','logical',...
                'DissipatedAndElasticEnergy',true,...
                'DTypeDissipatedAndElasticEnergy','logical',...
                'PeakIndentationAngle',true,...
                'DTypePeakIndentationAngle','logical',...
                'HeightMapByCurrentCP',true,...
                'DTypeHeightMapByCurrentCP','logical',...
                'CorrectedDZSlopes',true,...
                'DTypeCorrectedDZSlopes','logical',...
                'CorrectedDZSlopesOptions',CorrectedDZSlopesOptions...
                );
            
            SensitivityCorrection = struct(...
                'SensitivityCorrectionMethod','Automatic',...
                'AllowedSensitivityCorrectionMethod',{{'Automatic',...
                'KeepOriginal','Adaptive','AutomaticFibril','FromArea',...
                'SetAllToValue','UserInput','FromRefFM'}},...
                'FitRangeMode','ValueHorizontal',...
                'AllowedFitRangeMode',{{'ValueHorizontal',...
                'ValueVertical','FractionHorizontal','FractionVertical'}},...
                'FitRangeLowerValue',0,...
                'FitRangeUpperValue',4e-9,...
                'FitRangeLowerFraction',.75,...
                'FitRangeUpperFraction',1,...
                'MovingWindowSize',512,...
                'TooltipMovingWindowSize','Only for Sensitivity Correction Method "Adaptive"',...
                'AppRetSwitch',false,...
                'DTypeAppRetSwitch','logical',...
                'SetAllToValue',1 ...
                );
            
            OliverPharr = struct(...
                'CurvePercent',0.75);
            
            Hertz = struct(...
                'TipShape','parabolic',...
                'AllowedTipShape',{{'parabolic','spheric approx.','parabolic on thin film (not bonded)','parabolic on thin film (bonded)'}},...
                'TopographyHeightChannel','Processed',...
                'AllowedTopographyHeightChannel',{{'Processed','Contact Height','Deconvoluted','Height','Height (measured)'}},...
                'ThinFilmMode','From Topography',...
                'AllowedThinFilmMode',{{'From Topography','From Value'}},...
                'ThinFilmThickness',1e-6,...
                'FitDegreeForSneddonPolySurf',10,...
                'CorrectSensitivity',true,...
                'DTypeCorrectSensitivity','logical',...
                'AllowXShift',true,...
                'DTypeAllowXShift','logical',...
                'PredictiveRSquare',true,...
                'DTypePredictiveRSquare','logical',...
                'TooltipPredictiveRSquare','Calculate the RSquare of the fit on the whole set of data bigger than the CP, regardless of fitting range.',...
                'UpperForceCutOff',[],...
                'TooltipUpperForceCutOff','Set the upper absolute y-limit in Newtons [N] for the portion of the curve that should be considered for fitting.',...
                'LowerForceCutOff',0,...
                'TooltipLowerForceCutOff','Set the lower absolute y-limit in Newtons [N] for the portion of the curve that should be considered for fitting.',...
                'LowerCurveFraction',0,...
                'TooltipLowerCurveFraction','Set the lower fraction y-limit for the portion of the curve that should be considered for fitting. This applies to the curve AFTER the ForceCutOffs.',...
                'UpperCurveFraction',1,...
                'TooltipUpperCurveFraction','Set the upper fraction y-limit for the portion of the curve that should be considered for fitting. This applies to the curve AFTER the ForceCutOffs.',...
                'UseTipInHertz',true,...
                'DTypeTipInHertz','logical',...
                'UseTopography',false,...
                'DTypeUseTopography','logical',...
                'WeighPointsByInverseDistance',false,...
                'DTypeWeighPointsByInverseDistance','logical');
            
            EModOption = struct(...
                'Type','Hertz',...
                'AllowedType',{{'Hertz','OliverPharr','both'}},...
                'OliverPharr',OliverPharr,...
                'Hertz',Hertz);
            
            BaselineCorrection = struct(...
                'BaselineCorrectionMethod','Automatic',...
                'AllowedBaselineCorrectionMethod',{{'Automatic',...
                'FromFixedRange'}},...
                'TiltCorrection',true,...
                'DTypeTiltCorrection','logical',...
                'FitRangeMode','FractionHorizontal',...
                'AllowedFitRangeMode',{{'ValueHorizontal',...
                'ValueVertical','FractionHorizontal','FractionVertical'}},...
                'FitRangeLowerValue',0,...
                'FitRangeUpperValue',4e-9,...
                'FitRangeLowerFraction',0,...
                'FitRangeUpperFraction',.5...
                );
                
            
            FMAOptions = struct('ContactPointOption','Old',...
                'ContactPointJustXAxis',false,...
                'DTypeContactPointJustXAxis','logical',...
                'AllowedContactPointOption',{{'Old','ZoomSweep','Fast','Snap-In','RoV','GoF','Combo','manual','HardSurface','CurveOrigin','None'}},...
                'EModOption',EModOption,...
                'SensitivityCorrection',SensitivityCorrection,...
                'BaselineCorrection',BaselineCorrection,...
                'OptionalMetrics',OptionalMetrics,...
                'TemporaryLoadIn',true,...
                'DTypeTemporaryLoadIn','logical',...
                'ReferenceAppRetSwitch',0,...
                'KeepOldResults',false,...
                'DTypeKeepOldResults','logical',...
                'SkipAlreadyProcessedData',false,...
                'DTypeSkipAlreadyProcessedData','logical',...
                'SkipAreaExclusion',true,...
                'DTypeSkipAreaExclusion','logical',...
                'ResetFlagsSelectedAndCorrupted',true,...
                'DTypeResetFlagsSelectedAndCorrupted','logical',...
                'UnselectCurveFragmentsThreshold',1/2,...
                'UnselectCurveFragmentsAppRetSwitch',2,...
                'SaveWhenFinished',true,...
                'DTypeSaveWhenFinished','logical',...
                'SaveAfterEachMap',false,...
                'DTypeSaveAfterEachMap','logical',...
                'DTypeSortHeightDataForFit','logical',...
                'SortHeightDataForFit',true);
            
        end
        
        function FSPOptions = set_default_fsp_options()
            
            DBandingOptions = struct(...
                'A','A'...
                );
            
            FiberCharOptions = struct(...
                'A','A'...
                );
            
            FSPOptions = struct(...
                'CalculateFiberChar',true,...
                'TooltipCalculateFiberChar','Calculate certain characteristics (ellipse parameters, numerical height/area, length etc...) of fiber like polyline segments',...
                'FiberCharOptions',FiberCharOptions,...
                'TooltipFiberCharOptions','Calculate certain characteristics (ellipse parameters, numerical height/area, length etc...) of fiber like polyline segments',...
                'CalculateDBanding',true,...
                'TooltipCalculateDBanding','Calculate spatial frequency and amplitude of periodic pattern on fiber apex',...
                'DTypeCalculateDBanding','logical',...
                'DBandingOptions',DBandingOptions,...
                'TooltipDBandingOptions','Calculate spatial frequency and amplitude of periodic pattern on fiber apex',...
                'SaveWhenFinished',true,...
                'DTypeSaveWhenFinished','logical',...
                'SaveAfterEachMap',false,...
                'DTypeSaveAfterEachMap','logical'...
                );
            
        end
        
        function SISettings = set_default_show_image_settings()
            
            
            current = what();
            
            SISettings = struct(...
                'DefaultChannel1Index',1,...
                'DefaultChannel2Index',1,...
                'DefaultChannel1SubIndex',0,...
                'DefaultChannel2SubIndex',0,...
                'MainMenuValue',1,...
                'VolumeMenuValue',0,...
                'ResultsMenuValue',0,...
                'HasChannel2',0,...
                'ColorIndex',1,...
                'ReferenceFontSize',24,...
                'ProfileLineWidth',3,...
                'DefaultSavePath',current.path,...
                'DefaultSaveType','*.png',...
                'BothCrossSections',0,...
                'LockChannels',0,...
                'RelativeChannelIndex',0,...
                'LockScalebars',0,...
                'UseOverlay',0,...
                'StatisticalCMap',0,...
                'IsUpscaled',0,...
                'BaselineCorrected',1,...
                'TipHeight',1,...
                'ContactPointShifted',1,...
                'ShowHertzFit',1,...
                'PlotTime',0,...
                'UseCorrectedSensitivity',0,...
                'vDeflectionUnitIndex',1,...
                'ExtendedInformation',1,...
                'UseSnapped',1,...
                'PlotGroup',0,...
                'IgnoreZeros',1,...
                'JustChannel1',0,...
                'JustChannel2',0,...
                'GrammX',4,...
                'GrammY',2,...
                'GrammGroup',5,...
                'GrammFacets',1);
            
        end
        
        function GrammOptions = set_default_gramm_options()
            
            Subplots = struct(...
                'EnableFacets',false,...
                'DTypeEnableFacets','logical',...
                'TooltipEnableFacets','Enable the split of data according to Facets into subplots',...
                'ncols',4,...
                'Tooltipncols','After how many columns do we wrap and create a new row',...
                'scale','fixed',...
                'Allowedscale',{{'fixed','free_x','free_y','free','independent'}},...
                'Tooltipscale','Use to provide data that will determine separation between subblots rows and columns. First argument provided will separate along rows, second will separate along columns',...
                'space','fixed',...
                'Allowedspace',{{'fixed','free_x','free_y','free'}},...
                'Tooltipspace','Use to provide data that will determine separation between subblots rows and columns. First argument provided will separate along rows, second will separate along columns',...
                'column_labels',true,...
                'DTypecolumn_labels','logical',...
                'Tooltipcolumn_labels','Do we label subplot columns',...
                'row_labels',true,...
                'DTyperow_labels','logical',...
                'Tooltiprow_labels','Do we label subplot rows',...
                'force_ticks',true,...
                'DTypeforce_ticks','logical',...
                'Tooltipforce_ticks','Do we override defaults and force ticks on all subplots'...
                );
            
            Direct = struct(...
                'Style','Point',...
                'AllowedStyle',{{'Point','Jitter','Line','Raster','Bar',...
                'Interval'}},...
                'Dodge',.5,...
                'DTypeDodge','double',...
                'TooltipDodge','(For Point, Jitter, Line, Bar, Interval, Label) When using multiple colors, use to dodge graphical elements between colors with the same x value',...
                'Alpha',1,...
                'DTypeAlpha','double',...
                'TooltipAlpha','(For Point, Jitter, Line) Set the transparency of lines (0:fully transparent, 1: solid;',...
                'Width',0.3,...
                'DTypeWidth','double',...
                'TooltipWidth','(For Jitter, Bar, Interval) Provide to set the width of errorbars OR How much are the points jittered in horizontal direction (in data units)',...
                'Height',0,...
                'DTypeHeight','double',...
                'TooltipHeight','(For Jitter) How much are the points jittered in vertical direction (in data units)',...
                'Stacked',true,...
                'DTypeStacked','logical',...
                'TooltipStacked','(For Bar) Set to true to have bars placed at the same x stacked;'...
                );
            
            Statistical = struct(...
                'Style','Summary',...
                'AllowedStyle',{{'Summary','Smooth','glm','fit','bin',...
                'cornerhist','density','bin2d','ellipse','qq','boxplot',...
                'violin'}},...
                'TooltipStyle','Choose the type of summary-style plot',...
                'Summary_geom','bar',...
                'AllowedSummary_geom',{{'area','lines','line','solid_area',...
                'black_errorbar','errorbar','bar','point','area_only'}},...
                'TooltipSummary_geom','Determines how means and errorbars are represented in summary, smooth, glm and fit',...
                'black_errorbar',true,...
                'DTypeblack_errorbar','logical',...
                'Tooltipblack_errorbar','Do we want an additional black errorbar displayed on top of "area" option.',...
                'Summary_type','ci',...
                'AllowedSummary_type',{{'ci','bootci','sem','std','quartile',...
                '95percentile','fitnormalci','fitpoissonci','fitbinomialci'}},...
                'TooltipSummarytype','Method of Bar Height calculation',...
                'Summary_setylim',true,...
                'DTypeSummary_setylim','logical',...
                'TooltipSummary_setylim','Do we set the YLim for the subplot according to the summary or the data?',...
                'Summary_width',.6,...
                'DTypeSummary_width','double',...
                'TooltipSummary_width','Provide to set the width of bars and errorbars',...
                'Summary_dodge',.7,...
                'DTypeSummary_dodge','double',...
                'TooltipSummary_dodge','When using multiple colors, use to dodge graphical elements between colors with the same x value',...
                'Smooth_geom','area',...
                'AllowedSmooth_geom',{{'area','lines','line','solid_area',...
                'black_errorbar','errorbar','bar','point','area_only'}},...
                'TooltipSmooth_geom','Determines how means and errorbars are represented in summary, smooth, glm and fit',...
                'Smooth_method','eilers',...
                'AllowedSmooth_method',{{'eilers','smoothingspline','moving','lowess','sgolay'}},...
                'TooltipSmooth_method','Determines method of smoothing',...
                'Smooth_lambda',[],...
                'DTypeSmooth_lambda','double',...
                'TooltipSmooth_lambda','Smoothing parameter, depends on method, see documentation',...
                'Smooth_npoints',200,...
                'DTypeSmooth_npoints','double',...
                'TooltipSmooth_npoints','Number of points over which the smooth is evaluated',...
                'Glm_geom','area',...
                'AllowedGlm_geom',{{'area','lines','line','solid_area',...
                'black_errorbar','errorbar','bar','point','area_only'}},...
                'TooltipGlm_geom','Determines how means and errorbars are represented in summary, smooth, glm and fit',...
                'Glm_distribution','normal',...
                'AllowedGlm_distribution',{{'normal','binomial','poisson','gamma','inverse gaussian'}},...
                'TooltipGlm_distribution','Type of distribution to be fitted',...
                'Glm_fullrange',false,...
                'DTypeGlm_fullrange','logical',...
                'TooltipGlm_fullrange','Do we display the fit over the whole x axis, or just on the range of the value used for the fit',...
                'Glm_disp_fit',true,...
                'DTypeGlm_disp_fit','logical',...
                'TooltipGlm_disp_fit','Do we display the fitted equations',...
                'Fit_geom','area',...
                'AllowedFit_geom',{{'area','lines','line','solid_area',...
                'black_errorbar','errorbar','bar','point','area_only'}},...
                'TooltipFit_geom','Determines how means and errorbars are represented in summary, smooth, glm and fit',...
                'Fit_fun','@(param1,param2,x)x.^param1+param2',...
                'DTypeFit_fun','char',...
                'TooltipFit_fun','Anonymous function with parameters to fit as first arguments and x as last argument',...
                'Fit_intopt','observation',...
                'AllowedFit_intopt',{{'observation','functional'}},...
                'TooltipFit_intopt','"observation": 95% bounds on a new observation; "functional": 95% bounds for the fitted function',...
                'Fit_fullrange',false,...
                'DTypeFit_fullrange','logical',...
                'TooltipFit_fullrange','Do we display the fit over the whole x axis, or just on the range of the value used for the fit',...
                'Fit_disp_fit',true,...
                'DTypeFit_disp_fit','logical',...
                'TooltipFit_disp_fit','Do we display the fitted equations',...
                'Bin_nbins',30,...
                'DTypeBin_nbins','double',...
                'TooltipBin_nbnins','Number of bins',...
                'Bin_geom','bar',...
                'AllowedBin_geom',{{'bar','line','overlaid_bar','stacked_bars','stairs','point'}},...
                'TooltipBin_geom','Determines how means and errorbars are represented',...
                'Bin_normalization','count',...
                'AllowedBin_normalization',{{'count','countdensity',...
                'cumcount','probability','pdf','cdf'}},...
                'TooltipBin_normalization','Type of normalization',...
                'Bin_fill','face',...
                'AllowedBin_fill',{{'face','edge','all','transparent'}},...
                'TooltipBin_fill','How are the bars colored/filled',...
                'Bin_width',.6,...
                'DTypeBin_width','double',...
                'TooltipBin_width','Provide to specify width of bars',...
                'Bin_dodge',.7,...
                'DTypeBin_dodge','double',...
                'TooltipBin_dodge','Provide to specify dodging between elements',...
                'Cornerhist_location',.8,...
                'DTypeCornerhist_location','double',...
                'TooltipCornerhist_location','x (or y) location of the inset axis on the unity line of the parent: 0<location<1',...
                'Cornerhist_aspect',.3,...
                'DTypeCornerhist_aspect','double',...
                'TooltipCornerhist_aspect','Aspect ratio (y/x) of the inset axis',...
                'Cornerhist_nbins',30,...
                'DTypeCornerhist_nbins','double',...
                'TooltipCornerhist_nbnins','Number of bins',...
                'Cornerhist_geom','bar',...
                'AllowedCornerhist_geom',{{'bar','line','overlaid_bar','stairs','point'}},...
                'TooltipCornerhist_geom','Determines how means and errorbars are represented',...
                'Cornerhist_normalization','count',...
                'AllowedCornerhist_normalization',{{'count','countdensity',...
                'cumcount','probability','pdf','cdf'}},...
                'TooltipCornerhist_normalization','Type of normalization',...
                'Cornerhist_fill','face',...
                'AllowedCornerhist_fill',{{'face','edge','all','transparent'}},...
                'TooltipCornerhist_fill','How are the bars colored/filled',...
                'Cornerhist_width',.6,...
                'DTypeCornerhist_width','double',...
                'TooltipCornerhist_width','Provide to specify width of bars',...
                'Cornerhist_dodge',.7,...
                'DTypeCornerhist_dodge','double',...
                'TooltipCornerhist_dodge','Provide to specify dodging between elements',...
                'Density_bandwidth',.1,...
                'DTypeDensity_bandwidth','double',...
                'TooltipDensity_bandwidth','The bandwidth of the kernel-smoothing window, which is a function of the number of points in x',...
                'Density_function','pdf',...
                'AllowedDensity_function',{{'pdf','cdf','icdf','survivor','cumhazard'}},...
                'TooltipDensity_function','Function to estimate',...
                'Density_kernel','normal',...
                'AllowedDensity_kernel',{{'normal','box','triangle','epanechnikov'}},...
                'TooltipDensity_kernel','Type of kernel smoother',...
                'Density_npoints',100,...
                'DTypeDensity_npoints','double',...
                'TooltipDensity_npoints','How many points are used to plot the density',...
                'Density_extra_x',0,...
                'DTypeDensity_extra_x','double',...
                'TooltipDensity_extra_x','Extend the x value range over which the density is evaluated',...
                'Bin2D_geom','area',...
                'AllowedBin2D_geom',{{'image','contour'}},...
                'TooltipBin2D_geom','Determines how means and errorbars are represented',...
                'Ellipse_type','95percentile',...
                'AllowedEllipse_type',{{'95percentile','ci'}},...
                'TooltipEllipse_type','95percentile: Fit ellipse that contains 95% of the points (assuming bivariate normal); ci: Fit ellipse that contains 95% of the bootstrapped xy means',...
                'Ellipse_geom','area',...
                'AllowedEllipse_geom',{{'area','line'}},...
                'TooltipEllipse_geom','area: Plot the ellipse as a shaded area with outline; line: Just plot the outline of the ellipse',...
                'QQ_distribution',"makedist('Normal',mean(X,'omitnan'),std(X,'omitnan'))",...
                'DTypeQQ_distribution','char',...
                'TooltipQQ_distribution',"Provide a theoretical distribution to plot x against using Matlab's makedist() function. Set to 'y' to plot x against y densities.",...
                'Boxplot_width',0.6,...
                'DTypeBoxplot_width','double',...
                'TooltipBoxplot_width','Width of boxes',...
                'Boxplot_dodge',0.7,...
                'DTypeBoxplot_dodge','double',...
                'TooltipBoxplot_dodge','Dodging between boxes of different colors within unique x values',...
                'Boxplot_notch',false,...
                'DTypeBoxplot_notch','logical',...
                'TooltipBoxplot_notch','Add notches at median ± 1.58 IQR /sqrt(N) to the boxplot',...
                'Violin_normalization','area',...
                'AllowedViolin_normalization',{{'area','count','width'}},...
                'TooltipViolin_normalization','area: Equal violin areas; count: Areas proportional to point count; Equal violin widths',...
                'Violin_half',false,...
                'DTypeViolin_half','logical',...
                'TooltipViolin_half','Plot half the violin',...
                'Violin_width',0.6,...
                'DTypeViolin_width','double',...
                'TooltipViolin_width','Width of violins',...
                'Violin_dodge',0.7,...
                'DTypeViolin_dodge','double',...
                'TooltipViolin_dodge','Dodging between violins of different colors within unique x values',...
                'Violin_fill','face',...
                'AllowedViolin_fill',{{'face','edge','all','transparent'}},...
                'TooltipViolin_fill','How are the violins colored/filled',...
                'Violin_bandwidth',.1,...
                'DTypeViolin_bandwidth','double',...
                'TooltipViolin_bandwidth','The bandwidth of the kernel-smoothing window, which is a function of the number of points in x',...
                'Violin_kernel','normal',...
                'AllowedViolin_kernel',{{'normal','box','triangle','epanechnikov'}},...
                'TooltipViolin_kernel','Type of kernel smoother',...
                'Violin_npoints',100,...
                'DTypeViolin_npoints','double',...
                'TooltipViolin_npoints','How many points are used to plot the density',...
                'Violin_extra_y',0,...
                'DTypeViolin_extra_y','double',...
                'TooltipViolin_extra_y','Extend the y value range over which the density is evaluated');
            
            GrammOptions = struct(...
                'GroupingStyle','lightness',...
                'AllowedGroupingStyle',{{'lightness','color','linestyle',...
                'marker','size'}},...
                'Colorset','lch',...
                'AllowedColorset',{{'lch','matlab','brewer_1','brewer_2','brewer_3',...
                'brewer_pastel','brewer_dark','brewer_paired',...
                'd3_10','d3_20','d3_20b','d3_20c' }},...
                'color_legend','separate',...
                'Allowedcolor_legend',{{'separate-gray','separate','expand','merge'}},...
                'lightness',65,...
                'DTypelighness','double',...
                'Tooltiplightness','Only applies to grouping style color and colorset lch',...
                'chroma',75,...
                'DTypechroma','double',...
                'Tooltipchroma','Only applies to grouping style color and colorset lch',...
                'NoLegend',false,...
                'DTypeNoLegend','logical',...
                'TooltipNoLegend','Disables legends',...
                'FlipCoordinates',false,...
                'DTypeFlipCoordinates','logical',...
                'TooltipFlipCoordinates','Flips the x with the y axis. Also flips boxplots, bars, violinplots etc.',...
                'ymin',[],...
                'DTypeymin','double',...
                'ymax',[],...
                'DTypeymax','double',...
                'VisualizationType','Statistical',...
                'AllowedVisualizationType',{{'Direct','Statistical'}},...
                'DirVisOptions',Direct,...
                'StatVisOptions',Statistical,...
                'SubplotOptions',Subplots,...
                'TooltipSubplots','Subplots/Faceting and multiple figures'...
                );
                
            
        end
        
        function CustomCantileverTipOptions = set_custom_cantilever_tip_options()
            
            ParaboloidOptions = struct(...
                'Radius',40e-9,...
                'TooltipRadius','Radius R in meters relating to the parabola of revolution z=x^2*a with a=1/(2R)'...
                );
            
            SphereOptions = struct(...
                'Radius',10e-6...
                );
            
            PyramidOptions = struct(...
                'FrontHalfAngle',20,...
                'TooltipFrontHalfAngle','Half angle of the front facing pyramid side in degrees',...
                'BackHalfAngle',25,...
                'TooltipBackHalfAngle','Half angle of the back facing pyramid side in degrees',...
                'SideHalfAngle',25,...
                'TooltipSideHalfAngle','Half angle of the two side facing pyramid sides in degrees'...
                );
            
            CustomCantileverTipOptions = struct(...
                'TipType','Paraboloid',...
                'AllowedTipType',{{'Paraboloid','Sphere',...
                'Pyramid'}},...
                'TipRadius',10e-9,...
                'TooltipTipRadius','Determines realistic rounding of the tip by tapering the original form into a paraboloid of given radius in meters. A radius of 0 leaves the original shape unchanged',...
                'TipHeight',1000e-9,...
                'TooltipTipHeight','Determines the overall height of the tip',...
                'TipTilt',10,...
                'TooltipTipTilt','Overall tilt of the chosen shape in degrees. Tilt is usually determined by the angle the cantilever is mounted at',...
                'ImageResolution',512,...
                'TooltipImageResolution','Chosen width in pixels for resulting square image',...
                'ParaboloidOptions',ParaboloidOptions,...
                'SphereOptions',SphereOptions,...
                'PyramidOptions',PyramidOptions...
                );
            
        end
        
        function OutPath = replace_fileseps(InPath)
            
            OutPath = strrep(InPath,'/',filesep);
            OutPath = strrep(InPath,'\',filesep);
            
        end
        
        function load_empty_afmimage(SourceFile,TargetPath)
            
            FullFile = which()
            
        end
        
        function IsValid = check_segment_tag_validity(String)
            
            IsValid = true;
            
            if isempty(String) ||...
                    (~ischar(String) && ~isstring(String)) ||...
                    contains(String,'*')
                IsValid = false;
            end
            
        end
        
    end
    
    methods(Static)
        % General static auxilary methods
        
        function [UnitCell,FancyNameCell] = assign_units_and_fancy_names_to_properties(InCell)
            % UnitCell = assign_units_to_properties(InCell)
            %
            % The Unit Table has to be manually maintained
            
            UnitFancyTable = {...
                'AdaptiveSensitivity','m/V','Adaptive Sensitivity';...
                'AdhesionEnergy','J','Adhesion Energy';...
                'AdhesionLength','m','Adhesion Length';...
                'BaseAndTiltFlag','logical','Base and Tilt Flag';...
                'BigDataFlag','logical','Big Data Flag';...
                'CorruptedCurves','logical','Corrupted Curves';...
                'Date','yy-mm-dd','Date';...
                'DateTime','yy-mm-dd-hh-mm-ss','DateTime';...
                'DissipatedEnergy','J','Dissipated Energy';...
                'EModHertz','Pa','Indentation Modulus Hertz';...
                'EModHertz_Old','Pa','Indentation Modulus Hertz (Old-CP)';...
                'ElasticEnergy','J','Elastic Energy';...
                'ExtendTime','s','Extend Time';...
                'ExtendVelocity','m/s','Extend Velocity';...
                'ExtendZLength','m','Extend Z-Length';...
                'FileType','-','File Type';...
                'FractionOfMaxRAM','-','Fraction of Maximum Ram';...
                'HHType','-','Head Height Type';...
                'HasRefSlope','logical','Has Reference Slope';...
                'HertzFitType','-','Hertz Fit Type';...
                'HoldingTime','s','Holding Time';...
                'ID','-','ID';...
                'ImagingType','-','Imaging Type';...
                'IndentationDepth','m','Indentationm Depth';...
                'IndentationDepthHertz','m','Indentation Depth Hertz';...
                'KeepPythonFilesOpen','logical','Keep Python Files Open';...
                'MaxAdhesionForce','N','Maximum Adhesion Force';...
                'MaxPointsPerCurve','-','Maximum Points per Force Curve';...
                'NCurves','-','Number of Curves';...
                'Name','-','Name';...
                'NumChannels','-','Number of Channels';...
                'NumPixelsX','-','Number of Pixels in X';...
                'NumPixelsY','-','Number of Pixels in Y';...
                'NumSegments','-','Number of Segments';...
                'OriginX','m','Origin X';...
                'OriginY','m','Origin Y';...
                'PeakIndentationAngle','°','Peak Indentaion Angle';...
                'PoissonR','-','Poisson Ratio';...
                'PythonLoaderFlag','logical','Python Loader Flag';...
                'RetractTime','s','Retract Time';...
                'RetractVelocity','m/s','Retract Velocity';...
                'RetractZLength','m','Retract Length';...
                'ScanAngle','°','Scan Angle';...
                'ScanSizeX','m','Scan Size X';...
                'ScanSizeY','m','Scan Size Y';...
                'SelectedCurves','logical','Selected Curves';...
                'Sensitivity','m/V','Sensitvity';...
                'Setpoint','N','Setpoint';...
                'SpringConstant','N/m','Cantilever Spring Constant';...
                'Time','s','Time';...
                'TipRadius','m','Tip Radius';...
                'hasBackgroundMask','logical','Has Background Mask';...
                'hasDeconvolutedCantileverTip','logical','Has Deconvoluted Cantilever Tip';...
                'hasErrorSignal','logical','Has Error Signal';...
                'hasHeight','logical','Has Height Channel';...
                'hasHeightMeasured','logical','Has Height(measured) Channel';...
                'hasLateralDeflection','logical','Has Lateral Deflection Channel';...
                'hasLockInAmplitude','logical','Has Lock-In Amplitude Channel';...
                'hasLockInPhase','logical','Has Lock-In Phase Channel';...
                'hasOverlay','logical','Has Overlay Group';...
                'hasProcessed','logical','Has Processed Channel';...
                'hasSensitivity','logical','Has Sensitivity';...
                'hasSpringConstant','logical','Has Cantilever Spring Constant';...
                'hasVerticalDeflection','logical','Has Vertical Deflection Channel';...
                'FiberSegment_Area','m2','Fiber Numerical Area';...
                'FiberSegment_AreaDerivedDiameter','m','Fiber Area-Derived Diameter';...
                'FiberSegment_Height','m','Fiber Height';...
                'FiberSegment_WidthHalfHeight','m','Fiber Width at Half-Height';...
                'FiberSegment_Prominence','m','Fiber Prominence';...
                'FiberSegment_WidthHalfProminence','m','Fiber Width at Half-Prominence';...
                'FiberSegment_WidthBase','m','Fiber Width at Base';...
                'FiberSegment_AspectRatioHalfHeight','m/m','Fiber Aspect Ratio at Half-Height';...
                'FiberSegment_AspectRatioBaseHeight','m/m','Fiber Aspect Ratio at Base-Height';...
                'FiberSegment_Mean_Area','m2','Fiber Mean Numerical Area';...
                'FiberSegment_Mean_AreaDerivedDiameter','m','Fiber Mean Area-Derived Diameter';...
                'FiberSegment_Mean_Height','m','Fiber Mean Height';...
                'FiberSegment_Mean_WidthHalfHeight','m','Fiber Mean Width at Half-Height';...
                'FiberSegment_Mean_Prominence','m','Fiber Mean Prominence';...
                'FiberSegment_Mean_WidthHalfProminence','m','Fiber Mean Width at Half-Prominence';...
                'FiberSegment_Mean_WidthBase','m','Fiber Mean Width at Base';...
                'FiberSegment_Mean_AspectRatioHalfHeight','m/m','Fiber Mean Aspect Ratio at Half-Height';...
                'FiberSegment_Mean_AspectRatioBaseHeight','m/m','Fiber Mean Aspect Ratio at Base-Height';...
                'FiberSegment_Median_Area','m2','Fiber Median Numerical Area';...
                'FiberSegment_Median_AreaDerivedDiameter','m','Fiber Median Area-Derived Diameter';...
                'FiberSegment_Median_Height','m','Fiber Median Height';...
                'FiberSegment_Median_WidthHalfHeight','m','Fiber Median Width at Half-Height';...
                'FiberSegment_Median_Prominence','m','Fiber Median Prominence';...
                'FiberSegment_Median_WidthHalfProminence','m','Fiber Median Width at Half-Prominence';...
                'FiberSegment_Median_WidthBase','m','Fiber Median Width at Base';...
                'FiberSegment_Median_AspectRatioHalfHeight','m/m','Fiber Median Aspect Ratio at Half-Height';...
                'FiberSegment_Median_AspectRatioBaseHeight','m/m','Fiber Median Aspect Ratio at Base-Height';...
                'FiberSegment_RelativePixelPosition','-','Fiber Relative Pixel Position';...
                'FiberSegment_RelativePosition','m','Fiber Relative Position';...
                'FiberSegment_SegmentLength','m','Fiber Segment Length';...
                'FiberSegment_Ellipse_a','m','Fiber Ellipse a';...
                'FiberSegment_Ellipse_b','m','Fiber Ellipse b';...
                'FiberSegment_Ellipse_AspectRatio','m/m','Fiber Ellipse Aspect Ratio';...
                'FiberSegment_Ellipse_Area','m2','Fiber Ellipse Area';...
                'FiberSegment_Ellipse_Height','m','Fiber Ellipse Height';...
                'FiberSegment_Ellipse_WidthHalfHeight','m','Fiber Ellipse Width at Half-Height';...
                'FiberSegment_Mean_Ellipse_a','m','Fiber Mean Ellipse a';...
                'FiberSegment_Mean_Ellipse_b','m','Fiber Mean Ellipse b';...
                'FiberSegment_Mean_Ellipse_AspectRatio','m/m','Fiber Mean Ellipse Aspect Ratio';...
                'FiberSegment_Mean_Ellipse_Area','m2','Fiber Mean Ellipse Area';...
                'FiberSegment_Mean_Ellipse_Height','m','Fiber Mean Ellipse Height';...
                'FiberSegment_Mean_Ellipse_WidthHalfHeight','m','Fiber Mean Ellipse Width at Half-Height';...
                'FiberSegment_Median_Ellipse_a','m','Fiber Median Ellipse a';...
                'FiberSegment_Median_Ellipse_b','m','Fiber Median Ellipse b';...
                'FiberSegment_Median_Ellipse_AspectRatio','m/m','Fiber Median Ellipse Aspect Ratio';...
                'FiberSegment_Median_Ellipse_Area','m2','Fiber Median Ellipse Area';...
                'FiberSegment_Median_Ellipse_Height','m','Fiber Median Ellipse Height';...
                'FiberSegment_Median_Ellipse_WidthHalfHeight','m','Fiber Median Ellipse Width at Half-Height';...
                };
            
            UnitCell = cell(size(InCell));
            FancyNameCell = cell(size(InCell));
            
            for i=1:length(InCell)
                IndexArray = strcmp(InCell{i},UnitFancyTable);
                if any(IndexArray(:,1))
                    UnitCell{i} = UnitFancyTable{IndexArray(:,1),2};
                    FancyNameCell{i} = UnitFancyTable{IndexArray(:,1),3};
                else
                    UnitCell{i} = '-';
                    FancyNameCell{i} = InCell{i};
                end
            end
            
        end
        
        function TCombo = combine_data_tables(T1,T2)
            
            T1Names = T1.Properties.VariableNames;
            T2Names = T2.Properties.VariableNames;
            T1Types =  varfun(@class,T1,'OutputFormat','cell');
            T2Types =  varfun(@class,T2,'OutputFormat','cell');
            
            U = unique([T1Names T2Names]);
            
            for i=1:numel(U)
                % T1
                Idx1 = find(strcmp(T1Names, U{i}));
                Idx2 = find(strcmp(T2Names, U{i}));
                if isempty(Idx1)
                    T1{:,T2Names{Idx2}} = eval([T2Types{Idx2} '(missing)']);
                end
                if isempty(Idx2)
                    T2{:,T1Names{Idx1}} = eval([T1Types{Idx1} '(missing)']);
                end
            end
            
            TCombo = [T1 ; T2];
            
            TCombo = standardizeMissing(TCombo,{'',[],""});
            
        end
        
    end
    
    methods(Static)
        % static methods for meta operations on multiple experiments
        
        function SummaryStruct = multiexperiment_force_map_analysis(FMAOptions)
            
            if nargin < 1
                FMAOptions = Experiment.set_default_fma_options();
            end
            
            k = 1;
            LoadMore = 'Yes';
            while isequal(LoadMore,'Yes')
                [File,Path] = uigetfile('*.mat','Choose Experiment .mat from folder');
                FullFile{k} = fullfile(Path,File);
                LoadMore = questdlg('Do you want to add more Experiments?',...
                    'Multiexperiment loader',...
                    'Yes',...
                    'No',...
                    'No');
                if isfile(FullFile{k})
                    k = k + 1;
                else
                    FullFile{k} = [];
                end
            end
            
            NumExperiments = length(FullFile);
            
            if isstruct(FMAOptions)
                OptionCell = cell(NumExperiments,1);
                for i=1:NumExperiments
                    OptionCell{i} = FMAOptions;
                end
            elseif NumExperiments~=length(FMAOptions)
                warning(['Number of given Force Map Analysis options ' ...
                    'does not equal number of experiments to be ' ...
                    'processed. Processing everything with first ' ...
                    'given option instead!'])
                OptionCell = cell(NumExperiments,1);
                for i=1:NumExperiments
                    OptionCell{i} = FMAOptions{i};
                end
            else
                OptionCell = FMAOptions;
            end
            
            SummaryStruct(1:NumExperiments) = struct(...
                'AnalysisSuccessful',true,...
                'ErrorStack',[]);
            tic
            h = waitbar(0,'Setting up...');
            for i=1:NumExperiments
                waitbar((i)/NumExperiments,h,...
                    ['Loading Experiment ' num2str(i) ' of ' num2str(NumExperiments)])
                try
                    E = Experiment.load(FullFile{i});
                    waitbar((i)/NumExperiments,h,...
                        ['Processing Experiment ' num2str(i) ' of ' num2str(NumExperiments)])
                    E.ForceMapAnalysisOptions = OptionCell{i};
                    E.force_map_analysis_general()
                catch ME
                    SummaryStruct(i).ErrorStack = ME;
                    SummaryStruct(i).AnalysisSuccessful = false;
                end
            end
            close(h)
            toc
        end
        
    end
end