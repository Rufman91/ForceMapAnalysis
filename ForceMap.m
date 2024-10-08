classdef ForceMap < matlab.mixin.Copyable & matlab.mixin.SetGet & handle  & dynamicprops & AFMBaseClass
    % The force map class represents a single jpk force map file and
    % contains all necessary functions to process the forcecurves.
    % General naming convention for this class is:
    % -class properties get capital letters and generally no space between
    % words
    % -class methods are in lowercase letters and get underscores between
    % words
    %%%%%%%%%%%%%%%%%%%%%%%%%%DISCLAIMER%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This is a handle-class and as such the associated class methods
    % should be  called by using the class-INSTANCE and not the classname
    % 'ForceMap'.
    % So if you load a map under the name 'FM' with the classconstructor
    % e.g.>> FM = ForceMap();
    % you should from now on call methods on this specific instance.
    % A valid call could for example be:
    % >> FM.base_and_tilt();
    % to conduct an operation (a base fit in this case) on the force map or
    % >> FM.TipRadius;
    % to get a class parameter of this force map (the tip radius of the used cantilever)
    
    properties
        % Properties shared for the whole Force Map. All data is given SI
        % units otherwise it would be stated separately 
        
        Date = ''            % date when the force map was detected
        Time = ''            % time when the force map was detected
        FileVersion = ''     % Version of jpk-force-map file
        DataStoreFolder = '' % 
        RawDataFilePath = ''
        OpenZipFile = []
        BigDataFlag = true     % If true, unpack data container into Experiment subfolder
                        % and always load force volume data from there
        PythonLoaderFlag = true
        KeepPythonFilesOpen = false % Decides whether to preload all PythonLoader Files into memory
                            % all the time
        FractionOfMaxRAM = 1/5 % Specifies how much of MaxRAM space can be taken for certain partitioned calculations 
        NCurves = []         % number of curves on the force map
        MaxPointsPerCurve = []
        HoldingTime = 0
        NumSegments = []
        ExtendTime = []
        RetractTime = []
        ExtendZLength = []
        RetractZLength = []
        ExtendVelocity = []        % Approach  velocity as defined in the force map settings
        RetractVelocity = []
        Sensitivity = []
        RefSlopeCorrectedSensitivity = []
        AdaptiveSensitivity = []
        SpringConstant = []
        Setpoint = []
        RescalingConstants = []
        DBanding = []        % Fourieranalysis-based estimate of DBanding perdiod (only available with sufficient resolution)
        RefSlope = []        % Refernce slope as determined from the upper curve slope from data from very hard      
        % surfaces (mica,glass), either from glass parts beneath the specimen or from
        % separate reference force maps
        PixApp = []          % maximum number of measured points during approach
        PixRet = []          % maximum number of measured points during retraction
        SelectedCurves = []  % logical vector of length NCurves with 0s for excluded and 1s for included curves. gets initialized with ones
        StashedSelectedCurves = [] % In some calculations, SelectedCurves is changed temporarily inside a method.
                              % Actual Selction is stored and then restored from here
        CorruptedCurves = [] % Curves that cant be loaded
        TipRadius = 10e-9  % (nominal, if not otherwise calculated) tip radius in m for the chosen type of tip model
        PoissonR = 0.5  % standard Poisson ratio for most mechanical models
        Medium = ''
        FibrilFlag = []
        FibPot = []
    end
    properties
        % Curve data Properties
        
        App = {}        % approach force data in Newton
        Ret = {}        % retraction force data in Newton
        HHApp = {}      % capacitive-sensor-height approach data in meters
        HHRet = {}      % capacitive-sensor-height retract data in meters
        HHType = ''          % Type of Head Height.(Default is capacitiveSensorHeight; switches to measuredHeight, if default does'nt exist)
        THApp = {}      % vertical tip height approach data in meters
        THRet = {}      % vertical tip height retract data in meters
        BasedApp = {}   % approach force data with subtracted base line and tilt in Newton
        BasedRet = {}   % retraction force data with subtracted base line and tilt in Newton
        BaseAndTiltFlag = false
    end
    properties
        % Properties related to Contact Point (CP) estimation
        
        RoV = {}        % ratio of variance curve for CP estiamation
        GoF = {}        % goodness of fit curve for each selected curve
        CP = []              % chosen contact point (CP) for every selected curve, which is then used in E module fitting
        CP_SnapIn = []       % Prefered CP-method for data with snap-in-effect (most data from air)
        CP_HertzFitted = []  % Computed from the 0-Point of the fit calculated in Hertz-Sneddon model with AllowXShift = true
        CP_RoV = []          % CP estimated with the ratio of variance method
        CP_GoF = []          % CP estimated with the goodness of fit method
        CP_Combo = []        % CP estimated with a combination of the GoF and RoV methods
        CPComboCurve = []    % combination of the various metrics for contact point estimation
        CP_CNN = []          % CP estimated with a conv. neural network in a single pass
        CP_CNNZoom = []      % CP estimated with a conv. neural network, first estimate CP, then zoom into the curve for more accurate pred
        CP_Dropout = []      % CP estimated with a conv. neural network in multiple Monte Carlo Dropout passes
        CP_CNNZoomSweep = [] % CP estimated with a conv. neural network, sweep over several magnifications and take mean (or median?) estimate
        CP_MonteCarlo = []   % All predictions from the multiple inference steps done in
        % the monte carlo method
        CP_MonteCarlo_STD = [] % Standard deviation of CP_MonteCarlo
        MiniBatchSize = [] % Optimal MiniBatchSize for CNN-prediction for current system environment
        DeltaE = {}     %
        YDropPred = []       % Contains the Dropoutpredictions for every curve in the forcemap
        CP_Old = []          % contact point estimation from old script 'A_nIAFM_analysis_main'
        Man_CP = []          % manually chosen contact point
        CP_HardSurface = []  % Detract cantilever deflection for CP estimation
        CP_None = []         % Fills the CP and CP_None properties with zeros
        CP_CurveOrigin = []  % Assums the start of the force curve as cp
        CPFlag = []          % Struct containing booleans to indicate if a certain CP-type has been estimated
        
    end
    properties
        % Properties related to topological calculations, such as mapping and masking and visualisation
        
        HeightMap = []       % height profile map taken from the maximum head-height from approach max(hhapp)
        EModMapHertz = []    % E modulus profile map. same ordering as HeightMap
        EModMapOliverPharr = [] % """"
        FibDiam = []    % Estimated fibril diameter
        FibDiamSTD = []      % Estimated fibril diameter std
        FibMask = []         % Logical mask marking the whole fibril
        BackgroundMask = []  % Logical mask marking the glass/mica/even-substrate background
        RefSlopeMask = []    % Logical mask marking the areas to consider for reference slope calculation
        ExclMask = []        % Manually chosen areas that are to be excluded for calculations of fibril EMod, FibDiam, DBanding etc.
        Apex = []            % Value of highest pixel in each profile
        RectApex = []        % Value of rectified apex location in each profile
        ApexIndex = []       % Index of highest pixel in each profile (List indexing!)
        RectApexIndex = []   % Index of rectified apex location in each profile (List indexing!)
        hasOverlay = false
    end
    properties
        % Properties related to EMod calculations
        
        Basefit = {}    % fit model used for the baseline fit in the base_and_tilt method
        Baseline = []
        TrueZero = []
        EModHertz = []       % List of reduced smaple E-Modulus based on a the Hertz-Sneddon model
        EModOliverPharr = [] % List of reduced sample E-Modulus based on the Oliver-Pharr method
        FibrilEModOliverPharr = []
        FibrilEModHertz = []
        HertzFitStore = []
        HertzFit = []        % HertzFit model generated in the calculate_e_mod_hertz method
        HertzFitType = ''
        HertzFitCoeffNames = ''
        HertzFitValues = []
        HertzFitSSE = []
        HertzFitRSquare = []
        HertzFitDFE = []
        HertzFitAdjRSquare = []
        HertzFitRMSE = []
        HertzFitPredictiveRSquare = []
        SnapIn = []
        MaxAdhesionForce = []
        AdhesionEnergy = []
        AdhesionLength = []
        DissipatedEnergy = []
        ElasticEnergy = []
        PeakIndentationAngle = []
        IndentationDepth = []
        IndentationDepthHertz = []
        IndentationDepthHertzFitRange = []
        DZslope = []
        DZSlopeCorrected = []
        Stiffness = []
        IndentationDepthOliverPharr = []
        IndentArea = []
        ProjTipArea = []
        HasRefSlope = false
    end
    properties
        % auxiliary properties to facilitate comparing different methods of
        % CP estimation
        NeuralNetAccelerator = ''
        EModOliverPharr_CNN = []
        EModOliverPharr_Old = []
        EModOliverPharr_RoV = []
        EModHertz_CNN = []
        EModHertz_Old = []
        EModHertz_RoV = []
    end
    properties
        % SMFS related 
        
        Linker          % Used NHS-PEG-MI Linker (short or long) 
        Substrate       % Used substrate for the measurement 
        EnvCond         % Environmental condition during the experiment
        ChipCant        % AFM-Chip number and Cantilever label
        Chipbox         % AFM-Chipbox number (in Roman numerals)
        SMFSFlag        %
        BasedAppDataPts
        BasedRetDataPts
        BasedAppAlign={}   % capacitive-sensor-height approach data in meters with aligned datapoints to the corresponding retraction data
        BasedRetAlign={}   % capacitive-sensor-height retraction data in meters with aligned datapoints to the corresponding approach data
        xAppRetSubstracted={} 
        BasedRetCorr    % BasedRet data corrected based on a selection of the approach data
        BasedRetCorr2   % BasedRet data corrected based on a selection of the retraction data        
        CorrMeanApp     % Mean of a selection of the approach data for baseline correction
        CorrStdApp      % Corresponding standard deviation of the selection of the approach data for baseline correction
        CorrMeanRet     % Mean of a selection of the retraction data for baseline correction
        CorrStdRet      % Corresponding standard deviation of the selection of the retraction data for baseline correction
        EndIdx          % Index correspoding to a predefined value
        UnbindingBoundaryIdx % Unbinding boundary index
        yRetLim        % Modified retraction data to allow a valid integration
        yRetLim2        % Modified retraction data to allow a valid integration
        yAppLim         %  Modified approach data to allow a valid integration
        RetAdhEnergy_ThresholdMethod % Adhesion Energy of the retraction curve based on the threshold method
        AppAdhEnergy_IdxMethod % Adhesion Energy of the retraction curve based on the pulling length index method
        RetAdhEnergy_IdxMethod % Adhesion Energy of the retraction curve based on the pulling length index method
        FMAppAdhEnergyMean  % Mean of the adhesion energy of a force map
        FMAppAdhEnergyStd   % Standard deviation of the adhesion energy of a force map
        FMRetAdhEnergyMean  % Mean of the adhesion energy of a force map
        FMRetAdhEnergyStd   % Standard deviation of the adhesion energy of a force map
        SnapInIdx           % Index of the data point indentified as snap-in value 
        SnapInLength        % Snap-in length
        FMSnapInMean        % Mean snap-in length of a force map in meters (m) 
        FMSnapInMin         % Min snap-in length of a force map in meters (m) 
        FMSnapInMax         % Max snap-in length of a force map in meters (m)        
        PullingLengthIdx % Index of the data point identified as pulling length value
        PullingLength   % Pulling length in meters (m)
        FMPullingLengthMean  % Mean pulling length of a force map in meters (m)   
        FMPullingLengthMin   % Min pulling length of a force map in meters (m) 
        FMPullingLengthMax   % Max pulling length of a force map in meters (m)
        AdhForceMaxApp      % Maximum adhesion force of the approach data
        AdhForceMaxRet      % Maximum adhesion force of the retraction data
        AdhForceUnbinding 
        AdhForceMaxAppIdx   % Index of the maximum adhesion force of the approach data 
        AdhForceMaxRetIdx   % Index of the maximum adhesion force of the retraction data
        AdhForceUnbindingIdx
        Idx50nm             % Index closest to 50 nm distance 
        FitCoeffa1          % Sinus fit coefficient a1 derived from the fit function
        FitCoeffb1          % Sinus fit coefficient b1 derived from the fit function
        FitCoeffc1          % Sinus fit coefficient c1 derived from the fit function
        FitCoeffd1          % Sinus fit coefficient d1 derived from the fit function
        FitCoeffe1          % Sinus fit coefficient e1 derived from the fit function
        FitRSquare          % R sqare of the sinus fit
        FitSSE              % SSE of the sinus fit
        yAppFitCorr
        yRetFitCorr
        yRetFitMean
        yRetFitStd
        xDataToFit      % ONLY A TESTING PROPERTY
        yDataToFit      % ONLY A TESTING PROPERTY

    end
    
    methods
        % Main methods of the class
        
        function obj = ForceMap(MapFullFile,DataFolder,TempID,BigData,PythonLoaderFlag,KeepPythonFilesOpen,FakeOpt,NSynthCurves)
            %%% Constructor of the class
            
            % Specify the folder where the files live. And import them.
            % Also get curent folder and return to it after import of
            % files.
            % Assigns the properties that can be found in the jpk-file
            % already
            
            obj.ID = TempID;
            obj.BigDataFlag = BigData;
            obj.PythonLoaderFlag = PythonLoaderFlag;
            obj.KeepPythonFilesOpen = KeepPythonFilesOpen;
            
            if nargin >= 7 && isequal(FakeOpt,'Dummy')
                obj.create_dummy_force_map(NSynthCurves);
            end
            
            % get OS and use appropriate fitting system command
            obj.check_for_new_host
            
            % Unpack jpk-force-map with according to data loader choice
            obj.unpack_force_map_data(MapFullFile,DataFolder);
            
            Index = regexp(obj.ID,'(?<=\-).','all');
            LoadMessage = sprintf('loading data into ForceMap Nr.%s',obj.ID(Index(end):end));
            disp(LoadMessage)
            
            % reading header properties into object
            obj.read_in_header_properties;
            obj.SelectedCurves = true(obj.NCurves,1);
            obj.CorruptedCurves = false(obj.NCurves,1);
            
            obj.construct_list_to_map_relations
            
            try
                if isequal(obj.FileType,'nhf-spectroscopy')
                    obj.nhf_create_map_from_spectroscopy_data('Position Z','Height (measured)','m')
                    obj.preprocess_image;
                    Processed = obj.get_channel('Processed');
                    obj.HeightMap = Processed.Image;
                    obj.nhf_create_map_from_spectroscopy_data('Deflection','Maximum Force','N')
                    obj.nhf_create_map_from_spectroscopy_data('Lateral','Lateral Force','N')
                else
                    obj.read_jpk_images_from_files
                end
            catch ME
                warning(['Could not read standard Image Channels from file ' MapFullFile])
            end
            
            if ~obj.BigDataFlag
                %loading curve data into cell arrays
                obj.load_force_curves;
                
                %clean up unzipped jpk-force-map file
                rmdir(obj.DataStoreFolder,'s');
            end
            
            % intitialize masks
            obj.ExclMask = logical(ones(obj.NumPixelsX,obj.NumPixelsY));
            obj.FibMask = logical(zeros(obj.NumPixelsX,obj.NumPixelsY));
            
            try
                obj.create_pixel_difference_channel
                obj.set_channel_positions(obj.OriginX,obj.OriginY,obj.ScanAngle);
            catch ME
                warning(['Could not create Pixel Difference Channel for ' MapFullFile])
            end
            
            obj.initialize_flags();
            if ~obj.KeepPythonFilesOpen && obj.PythonLoaderFlag && obj.BigDataFlag
                obj.clear_zipped_files_from_memory
            end
        end
        
        function load_zipped_files_with_python(obj)
            % loads ForceMap source file (*.jpk-forrce-map,*.jpk-qi-data)
            % into memory and keeps it, open for the session
            
            if ~isempty(obj.OpenZipFile) && isequal(class(obj.OpenZipFile),'py.zipfile.ZipFile')
                disp(sprintf('File %s already loaded in',obj.Name))
                return
            end
            
            current = what;
            Split = split(obj.RawDataFilePath,filesep);
            Folder = fullfile(Split{1:end-1});
            File = Split{end};
            % If file path had a leading filesep, add it back in
            if isequal(obj.RawDataFilePath(1),filesep)
                if isequal(obj.RawDataFilePath(2),filesep)
                    Folder = strcat(filesep,Folder);
                    Folder = strcat(filesep,Folder);
                else
                    Folder = strcat(filesep,Folder);
                end
            end
            cd(Folder)
            obj.OpenZipFile = py.zipfile.ZipFile(File);
            
            cd(current.path)
        end
        
        function clear_zipped_files_from_memory(obj)
            
            if isempty(obj.OpenZipFile)
                return;
            end
            
            try
                obj.OpenZipFile.close();
                obj.OpenZipFile = [];
            catch
                obj.OpenZipFile = [];
            end
            
        end
        
        function choose_curves(obj)
            % Interactive dialogue for curve selection. Returns a logic 0/1
            % vector, that is later used and can also be changed by the
            % other force map methods.
            f = figure('Name','Curve Selection','Position',[10000 10000 1200 900]);
            movegui(f);
            for i=1:obj.NCurves
                [App,HHApp] = obj.get_force_curve_data(i,'AppRetSwitch',0,...
                    'BaselineCorrection',0,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                [Ret,HHRet] = obj.get_force_curve_data(i,'AppRetSwitch',1,...
                    'BaselineCorrection',0,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                dlgtitle = sprintf('Curve selection %i/%i',i,obj.NCurves);
                plottitle = sprintf('Curve Nr. %i',i);
                plot(HHApp,App,HHRet,Ret);
                title(plottitle);
                legend('Approach','Retract');
                answ = questdlg('Do you want to process this indentation curve?',...
                    dlgtitle,'Keep it','Exclude it','Abort Process',...
                    'Abort Process');
                if isequal(answ,'Keep it')
                elseif isequal(answ,'Abort Process')
                    close(gcf)
                    current = what();
                    cd(obj.Folder)
                    savename = sprintf('%s.mat',obj.Name);
                    save(savename,'obj')
                    cd(current.path)
                    return
                elseif isequal(answ,'Exclude it')
                    obj.SelectedCurves(i) = 0;
                end
            end
            close(gcf)
            %             current = what();
            %             cd(obj.Folder)
            %             savename = sprintf('%s.mat',obj.Name);
            %             save(savename,'obj')
            %             cd(current.path)
        end
        
        function unselect_curves_by_fraction_of_max_data_points(obj,Thresh,AppRetSwitch)
            % AppRetSwitch ... 0...App, 1...Ret, 2...Both
            
            Thresh = Thresh*obj.MaxPointsPerCurve;
            
            for i=1:obj.NCurves
                if AppRetSwitch==0 || AppRetSwitch==2
                    [App,HHApp] = obj.get_force_curve_data(i,'AppRetSwitch',0,...
                    'BaselineCorrection',0,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                    if (length(App) < Thresh) || (length(HHApp) < Thresh)
                        obj.SelectedCurves(i) = 0;
                    end
                end
                if AppRetSwitch==1 || AppRetSwitch==2
                    [Ret,HHRet] = obj.get_force_curve_data(i,'AppRetSwitch',1,...
                    'BaselineCorrection',0,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                    if (length(Ret) < Thresh) || (length(HHRet) < Thresh)
                        obj.SelectedCurves(i) = 0;
                    end
                end
                
            end
        end
        
        function base_and_tilt(obj,Options)
            % subtract baseline and tilt from the forcecurve by fitting a function to
            % the non contact domain the function tries to fit a non-affine-linear
            % function. If the linear fit is too bad, the function tries a 9th grade
            % polynomial fit instead
            
            if nargin<2
                DefOpts = Experiment.set_default_fma_options;
                Options = DefOpts.BaselineCorrection;
            end
            
            Range = find(obj.SelectedCurves);
            h = waitbar(0,'Setting up...');
            for i=Range'
                [vDef,Height] = obj.get_force_curve_data(i,'AppRetSwitch',0,...
                    'BaselineCorrection',0,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                
                prog = i/obj.NCurves;
                waitbar(prog,h,'processing baseline fits...');
                try
                    if isequal(Options.BaselineCorrectionMethod,'Automatic')
                        [ncd , ncidx] = ForceMap.no_contact_domain(vDef);
                        Params = polyfit(Height(1:length(ncd)),ncd,1);
                    elseif isequal(Options.BaselineCorrectionMethod,'FromFixedRange')
                        
                        switch Options.FitRangeMode
                            case 'ValueHorizontal'
                                UValue = max(Height) - Options.FitRangeLowerValue;
                                LValue = max(Height) - Options.FitRangeUpperValue;
                                FitHeight = Height(Height>=LValue & Height<=UValue);
                                FitvDef = vDef(Height>=LValue & Height<=UValue);
                            case 'ValueVertical'
                                FitHeight = Height(vDef>=Options.FitRangeLowerValue & vDef<=Options.FitRangeUpperValue);
                                FitvDef = vDef(vDef>=Options.FitRangeLowerValue & vDef<=Options.FitRangeUpperValue);
                            case 'FractionHorizontal'
                                RangeHeight = range(Height);
                                FitHeight = Height((Height-min(Height))>=Options.FitRangeLowerFraction*RangeHeight &...
                                    (Height-min(Height))<=Options.FitRangeUpperFraction*RangeHeight);
                                FitvDef = vDef((Height-min(Height))>=Options.FitRangeLowerFraction*RangeHeight &...
                                    (Height-min(Height))<=Options.FitRangeUpperFraction*RangeHeight);
                            case 'FractionVertical'
                                MaxvDef = max(vDef);
                                FitHeight = Height(vDef>=Options.FitRangeLowerFraction*MaxvDef & vDef<=Options.FitRangeUpperFraction*MaxvDef);
                                FitvDef = vDef(vDef>=Options.FitRangeLowerFraction*MaxvDef & vDef<=Options.FitRangeUpperFraction*MaxvDef);
                        end
                        Params = polyfit(FitHeight,FitvDef,1);
                    end
                    if Options.TiltCorrection
                        obj.Basefit{i} = Params;
                    else
                        obj.Basefit{i} = [0 Params(2)];
                    end
                catch
                    warning('Error in base and tilt. Skipping current force curve and marked as unselected')
                    obj.SelectedCurves(i) = false;
                end
            end
            
            close(h);
            obj.BaseAndTiltFlag = true;
        end
        
        function base_and_tilt_using_cp(obj,FractionBeforeCP)
            
            if ~sum(struct2array(obj.CPFlag))
                warning('No contact point data found. Run contact point estimation first!')
                return
            end
            
            if nargin < 2
                FractionBeforeCP = 1;
            end
            
            Range = find(obj.SelectedCurves);
            h = waitbar(0,'Refined Base and Tilt...');
            
            for i=Range'
                waitbar(i/obj.NCurves,h,'Refined Base and Tilt...');
                [App,HHApp] = obj.get_force_curve_data(i,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                [Ret,HHRet] = obj.get_force_curve_data(i,'AppRetSwitch',1,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                CP = obj.CP(i,:);
                
                Cutoff = CP(1) - (1-FractionBeforeCP)*(CP(1) - min(HHApp));
                
                HHAppFit = HHApp(HHApp<Cutoff);
                AppFit = App(1:length(HHAppFit));
                
                warning('off')
                Params = polyfit(HHAppFit,AppFit,1);
                obj.Basefit{i} = obj.Basefit{i} + Params;
                warning('on')
                
                FittedCP = polyval(Params,CP(1));
                obj.CP(i,2) = CP(2) - FittedCP;
                
                if obj.CPFlag.CNN
                    FittedCP = polyval(Params,obj.CP_CNN(i,1));
                    obj.CP_CNN(i,2) = obj.CP_CNN(i,2) - FittedCP;
                end
                if obj.CPFlag.CNNZoom
                    FittedCP = polyval(Params,obj.CP_CNNZoom(i,1));
                    obj.CP_CNNZoom(i,2) = obj.CP_CNNZoom(i,2) - FittedCP;
                end
                if obj.CPFlag.CNNZoomSweep
                    FittedCP = polyval(Params,obj.CP_CNNZoomSweep(i,1));
                    obj.CP_CNNZoomSweep(i,2) = obj.CP_CNNZoomSweep(i,2) - FittedCP;
                end
                if obj.CPFlag.Combo
                    FittedCP = polyval(Params,obj.CP_Combo(i,1));
                    obj.CP_Combo(i,2) = obj.CP_Combo(i,2) - FittedCP;
                end
                if obj.CPFlag.Dropout
                    FittedCP = polyval(Params,obj.CP_Dropout(i,1));
                    obj.CP_Dropout(i,2) = obj.CP_Dropout(i,2) - FittedCP;
                end
                if obj.CPFlag.GoF
                    FittedCP = polyval(Params,obj.CP_GoF(i,1));
                    obj.CP_GoF(i,2) = obj.CP_GoF(i,2) - FittedCP;
                end
                if obj.CPFlag.HardSurface
                    FittedCP = polyval(Params,obj.CP_HardSurface(i,1));
                    obj.CP_HardSurface(i,2) = obj.CP_HardSurface(i,2) - FittedCP;
                end
                if obj.CPFlag.HertzFitted
                    FittedCP = polyval(Params,obj.CP_HertzFitted(i,1));
                    obj.CP_HertzFitted(i,2) = obj.CP_HertzFitted(i,2) - FittedCP;
                end
                if obj.CPFlag.Old
                    FittedCP = polyval(Params,obj.CP_Old(i,1));
                    obj.CP_Old(i,2) = obj.CP_Old(i,2) - FittedCP;
                end
                if obj.CPFlag.SnapIn
                    FittedCP = polyval(Params,obj.CP_SnapIn(i,1));
                    obj.CP_SnapIn(i,2) = obj.CP_SnapIn(i,2) - FittedCP;
                end
                if obj.CPFlag.RoV
                    FittedCP = polyval(Params,obj.CP_RoV(i,1));
                    obj.CP_RoV(i,2) = obj.CP_RoV(i,2) - FittedCP;
                end
                if obj.CPFlag.Manual
                    FittedCP = polyval(Params,obj.Man_CP(i,1));
                    obj.Man_CP(i,2) = obj.Man_CP(i,2) - FittedCP;
                end
            end
            close(h)
        end
        
        function estimate_cp_snap_in(obj)
            % Takes the first point that crosses the approach baseline
            % starting from the max indent and sets the CP as the
            % interpolated point between the two nearest points to the
            % baseline (above and below)
            
            h = waitbar(0,'Setting up...','Name',obj.Name);
            Range = find(obj.SelectedCurves);
            SnapIn = zeros(obj.NCurves,1);
            for i=Range'
                waitbar(i/obj.NCurves,h,'Finding contact point for Snap-In curve...');
                [App,HHApp] = obj.get_force_curve_data(i,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                [Ret,HHRet] = obj.get_force_curve_data(i,'AppRetSwitch',1,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                AboveZeroBool = zeros(length(App),1);
                AboveZeroBool(find(App>0)) = 1;
                k = 0;
                if ~AboveZeroBool(end) || (sum(AboveZeroBool)==length(AboveZeroBool))
                    obj.SelectedCurves(i) = 0;
                    obj.CorruptedCurves(i) = 1;
                    warning(['Curve Nr. ' num2str(i) ' is consistently above 0 force [N]'])
                    continue
                end
                while AboveZeroBool(end-k)
                    k = k + 1;
                end
                if k==0
                    obj.CP(i,1) = HHApp(floor(.5*end));
                    obj.CP(i,2) = App(floor(.5*end));
                    obj.CP_SnapIn(i,:) = obj.CP(i,:);
                    obj.SelectedCurves(i) = 0;
                    continue
                end
                    
                AboveBase = [HHApp(end-(k-1)) App(end-(k-1))];
                BelowBase = [HHApp(end-k) App(end-k)];
                obj.CP(i,1) = mean([AboveBase(1) BelowBase(1)]);
                obj.CP(i,2) = 0;
                obj.CP_SnapIn(i,:) = obj.CP(i,:);
                clear AboveZeroBool
                
                try
                    Force = App - obj.CP(i,2);
                    SnapIn(i) = -min(Force);
                    
                catch
                    SnapIn(i) = nan;
                end
            end
            
            obj.SnapIn = SnapIn;
            
            for i=1:obj.NumPixelsX
                for j=1:obj.NumPixelsY
                    SnapInMap(i,j) = obj.SnapIn(obj.Map2List(i,j));
                end
            end
            
            % Write to Channel
            Channel = obj.create_standard_channel(SnapInMap,'Snap-In','N');
            [~,Index] = obj.get_channel('Snap-In');
            if isempty(Index)
                obj.Channel(end+1) = Channel;
            else
                obj.Channel(Index) = Channel;
            end
            
            close(h)
            obj.CPFlag.SnapIn = true;
            
        end
        
        function estimate_cp_rov(obj,WindowSize)
            % estimate_cp_rov(obj,WindowSize)
            % find contact point with the method of ratio of variances. The method
            % iterates though every point and builds the ratio of the variance of a
            % bunch of points before and after the current point. the point with the
            % biggest ratio is the returned contact point [Nuria Gavara, 2016]
            if nargin<2
                WindowSize = ceil(20*obj.MaxPointsPerCurve/300);
            end
            Range = find(obj.SelectedCurves);
            h = waitbar(0,'Setting up...','Name',obj.Name);
            CP_RoV = zeros(obj.NCurves,2);
            obj.CP_RoV = CP_RoV;
            for i=Range'
                prog = i/obj.NCurves;
                [App,HHApp] = obj.get_force_curve_data(i,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                waitbar(prog,h,'applying ratio of variances method...');
                obj.RoV{i} = zeros(length(App),1);
                SmoothedApp = smoothdata(App);
                % loop through points and calculate the RoV
                for j=(WindowSize+1):(length(App)-WindowSize)
                    obj.RoV{i}(j,1) = var(smoothdata(App((j+1):(j+WindowSize))))/...
                        var(SmoothedApp((j-WindowSize):(j-1)));
                end
                % normalize RoV-curve
                obj.RoV{i} = obj.RoV{i}/range(obj.RoV{i});
                minrov = min(obj.RoV{i}(WindowSize+1:length(obj.RoV{i})-WindowSize));
                obj.RoV{i}(obj.RoV{i}==0) = minrov;
                [~,CPidx] = max(obj.RoV{i});
                obj.CP_RoV(i,:) = [HHApp(CPidx) App(CPidx)];
                obj.CP = obj.CP_RoV;
            end
            close(h)
            obj.CPFlag.RoV = 1;
            %             current = what();
            %             cd(obj.Folder)
            %             savename = sprintf('%s.mat',obj.Name);
            %             save(savename,'obj')
            %             cd(current.path)
        end
        
        function estimate_cp_gof(obj)
            % estimate_cp_gof(obj)
            %
            % determine contact point with the goodness of fit method
            
            obj.CP_GoF = zeros(obj.NCurves,2);
            Range = find(obj.SelectedCurves);
            h = waitbar(0,'Setting up...','Name',obj.Name);
            TipRadius = obj.TipRadius;
            PoissonR = obj.PoissonR;
            for i=Range'
                [Based,THApp] = obj.get_force_curve_data(i,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',1,...
                    'Sensitivity','original','Unit','N');
                smoothx = smoothdata(THApp);
                smoothy = smoothdata(Based);
                Rsquare = 2*ones(length(Based),1);
                E = ones(length(Based),1);
                obj.DeltaE = ones(length(Based),1);
                testrange = floor(0.5*length(Based)):(length(Based)-5);
                msg = sprintf('applying goodness of fit method on curve Nr.%i/%i',i,obj.NCurves);
                parfor j=testrange
%                     prog = (j-testrange(1))/length(testrange);
%                     
                    [~,gof] = ForceMap.hertz_fit_gof(smoothx,smoothy,j,1,'parabolic',TipRadius,PoissonR);
                    Rsquare(j) = gof.rsquare;
                end
                waitbar(i/obj.NCurves,h,msg);
                Rsquare(Rsquare==2) = min(Rsquare);
                obj.GoF{i} = normalize(Rsquare,'range');
                [~,CPidx] = max(obj.GoF{i});
                obj.CP_GoF(i,:) = [THApp(CPidx) Based(CPidx)];
            end
            obj.CP = obj.CP_GoF;
            close(h)
            obj.CPFlag.GoF = 1;
            %             current = what();
            %             cd(obj.Folder)
            %             savename = sprintf('%s.mat',obj.Name);
            %             save(savename,'obj')
            %             cd(current.path)
        end
        
        function estimate_cp_combined(obj)
            Range = find(obj.SelectedCurves);
            obj.CP_Combo = zeros(obj.NCurves,2);
            obj.CPComboCurve = cell(obj.NCurves,1);
            for i=Range'
                [App,HHApp] = obj.get_force_curve_data(i,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',1,...
                    'Sensitivity','original','Unit','N');
                obj.CPComboCurve{i} = obj.RoV{i}.*obj.GoF{i};
                [~,CPidx] = max(obj.CPComboCurve{i});
                obj.CP_Combo(i,:) = [HHApp(CPidx) App(CPidx)];
                obj.CP(i,:) = obj.CP_Combo(i,:);
            end
            
            obj.CPFlag.Combo = 1;
            %             current = what();
            %             cd(obj.Folder)
            %             savename = sprintf('%s.mat',obj.Name);
            %             save(savename,'obj')
            %             cd(current.path)
        end
        
        function estimate_cp_cnn(obj,NeuralNet,RunMode,NumPasses)
            % [] = estimate_cp_oliver_pharr(obj,NeuralNet,RunMode,NumPasses)
            %
            % NeuralNet has to be passed down from the parent object
            % 'Experiment' to the ForceMap. For RunMode 'Fast' pass the
            % 'ExpInWorkspaceName.CP_CNN' property and for RunMode
            % 'Dropout' pass the 'ExpInWorkspaceName.Dropout' property
            %
            % RunMode = 'Fast' just predict with one forwardpass using
            % CP_CNN
            % RunMode = 'Dropout' predict through n=NumPasses forwardpasses
            % through DropoutNet and additionally gain a metric for
            % uncertainty of the CP estimate
            %
            % NumPasses is the positive integer defining how many passes
            % are to be done in RunMode 'Dropout'
            % Default NumPasses = 100
            % NumPasses >= 30 is recommended
            
            obj.check_for_cuda_capable_gpu_device()
            
            if nargin < 2
                runmode = 0;
            elseif isequal(lower(RunMode),'fast')
                runmode = 0;
            elseif isequal(lower(RunMode),'dropout')
                runmode = 1;
                if nargin < 3
                    NumPasses = 100; % if not specified in arguments, NumPasses defaults to 100
                end
                obj.CP_MonteCarlo = zeros(NumPasses,2,obj.NCurves);
            elseif isequal(lower(RunMode),'zoom')
                runmode = 2;
            elseif isequal(lower(RunMode),'zoomdropout')
                runmode = 3;
                if nargin < 3
                    NumPasses = 100; % if not specified in arguments, NumPasses defaults to 100
                end
                obj.CP_MonteCarlo = zeros(NumPasses,2,obj.NCurves);
            elseif isequal(lower(RunMode),'zoomsweep')
                runmode = 4;
                if nargin < 3
                    NumPasses = 5; % if not specified in arguments, NumPasses defaults to 5
                end
            end
            ImgSize = NeuralNet.Layers(1).InputSize;
            
            % Partition data using the SelectedCurves property to not
            % exceed memory capabilities of current system. Loop over all
            % partitions and in the end set SelectedCurves back to initial
            % values
            obj.StashedSelectedCurves = obj.SelectedCurves;
            try
                if isequal(obj.HostOS,'PCW')
                    Mem = memory;
                    MaxArraySize = Mem.MaxPossibleArrayBytes;
                else
                    try
                        [r,w] = unix('free | grep Mem');
                        stats = str2double(regexp(w, '[0-9]*', 'match'));
                        memsize = stats(1)/1e6;
                        freemem = (stats(3)+stats(end))/1e6;
                        MaxArraySize = freemem*10^9;
                    catch
                        MaxArraySize = 1e9;
                    end
                end
                if sum(runmode==[1 3 4],'all') >= 1
                    MaxPartitionSize = round(MaxArraySize/(ImgSize(1)*ImgSize(2)*ImgSize(3)*NumPasses));
                else
                    MaxPartitionSize = round(MaxArraySize/(ImgSize(1)*ImgSize(2)*ImgSize(3)));
                end
                PartitionSize = round(MaxPartitionSize*obj.FractionOfMaxRAM);
                NumPartitions = ceil(obj.NCurves/PartitionSize);
                
                for BigLoop=1:NumPartitions
                    TempSelectedCurves = zeros(obj.NCurves,1);
                    TempSelectedCurves((1+(BigLoop-1)*PartitionSize):min((BigLoop*PartitionSize),obj.NCurves)) = 1;
                    TempSelectedCurves = TempSelectedCurves & obj.StashedSelectedCurves;
                    obj.SelectedCurves = TempSelectedCurves;
                    objcell{1,1} = obj;
                    X = obj.CP_batchprep_3_channel(objcell,ImgSize(1),ImgSize(1),[0 0.3 0.7]);
                    h = waitbar(0,'Setting up and optimizing runtime...','Name',obj.Name);
                    len = size(X,4);
                    
                    obj.check_for_new_host();
                    obj.cnn_runtime_optimization(NeuralNet,X);
                    
                    CantHandle = true;
                    switch runmode
                        case 0
                            % Fast
                            waitbar(1/2,h,'Predicting CP');
                            while CantHandle == true
                                try
                                    Ypredicted = predict(NeuralNet,X,'MiniBatchSize',obj.MiniBatchSize,'Acceleration','auto',...
                                        'ExecutionEnvironment',obj.NeuralNetAccelerator);
                                    CantHandle = false;
                                catch
                                    obj.CPFlag.CNNopt = 0;
                                    obj.cnn_runtime_optimization(NeuralNet,X);
                                end
                            end
                            waitbar(1,h,'Wrapping up');
                            iRange = find(obj.SelectedCurves);
                            k = 1;
                            for i=iRange'
                                [App,HHApp] = obj.get_force_curve_data(i,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                                obj.CP_CNN(i,1) = Ypredicted(k,1)*range(HHApp)+min(HHApp);
                                obj.CP_CNN(i,2) = Ypredicted(k,2)*range(App)+min(App);
                                obj.CP(i,1) = obj.CP_CNN(i,1);
                                obj.CP(i,2) = obj.CP_CNN(i,2);
                                k = k + 1;
                            end
                            if BigLoop == NumPartitions
                                obj.CPFlag.CNN = 1;
                            end
                        case 1
                            % Dropout
                            obj.YDropPred = zeros(NumPasses,2,len);
                            for j=1:NumPasses
                                waitbar(j/NumPasses,h,sprintf('Predicting CP for %i curves. %i/%i passes done',len,j,NumPasses));
                                while CantHandle == true
                                    try
                                        Temp = predict(NeuralNet,X,'MiniBatchSize',obj.MiniBatchSize,'Acceleration','auto',...
                                        'ExecutionEnvironment',obj.NeuralNetAccelerator);
                                        CantHandle = false;
                                    catch
                                        CantHandle = true;
                                        obj.CPFlag.CNNopt = 0;
                                        obj.cnn_runtime_optimization(NeuralNet,X);
                                    end
                                end
                                obj.YDropPred(j,:,:) = Temp';
                                CantHandle = true;
                            end
                            waitbar(1,h,'Wrapping up');
                            iRange = find(obj.SelectedCurves);
                            k = 1;
                            for i=iRange'
                                [App,HHApp] = obj.get_force_curve_data(i,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                                obj.CP_MonteCarlo(:,1,i) = obj.YDropPred(:,1,k)*range(HHApp)+min(HHApp);
                                obj.CP_MonteCarlo(:,2,i) = obj.YDropPred(:,2,k)*range(App)+min(App);
                                obj.CP(i,1) = mean(obj.CP_MonteCarlo(:,1,i));
                                obj.CP(i,2) = mean(obj.CP_MonteCarlo(:,2,i));
                                obj.CP_Dropout(i,1) = obj.CP(i,1);
                                obj.CP_Dropout(i,2) = obj.CP(i,2);
                                obj.CP_MonteCarlo_STD(i) = norm([std(obj.CP_MonteCarlo(:,1,i)) std(obj.CP_MonteCarlo(:,2,i))]);
                                k = k + 1;
                            end
                            
                            if BigLoop == NumPartitions
                                obj.CPFlag.Dropout = 1;
                            end
                        case 2
                            % Zoom
                            waitbar(1/3,h,'Predicting CP, first guess...');
                            while CantHandle == true
                                try
                                    Ypredicted = predict(NeuralNet,X,'MiniBatchSize',obj.MiniBatchSize,'Acceleration','auto',...
                                        'ExecutionEnvironment',obj.NeuralNetAccelerator);
                                    CantHandle = false;
                                catch
                                    obj.CPFlag.CNNopt = 0;
                                    obj.cnn_runtime_optimization(NeuralNet,X);
                                end
                            end
                            iRange = find(obj.SelectedCurves);
                            k = 1;
                            for i=iRange'
                                [App,HHApp] = obj.get_force_curve_data(i,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                                obj.CP_CNNZoom(i,1) = Ypredicted(k,1)*range(HHApp)+min(HHApp);
                                obj.CP_CNNZoom(i,2) = Ypredicted(k,2)*range(App)+min(App);
                                obj.CP(i,1) = obj.CP_CNNZoom(i,1);
                                obj.CP(i,2) = obj.CP_CNNZoom(i,2);
                                k = k + 1;
                            end
                            
                            waitbar(2/3,h,'Predicting zoomed CP');
                            ZoomObj = obj.copy;
                            ZoomObj.cnn_zoom_in(0.5);
                            ZoomCell{1,1} = ZoomObj;
                            X = obj.CP_batchprep_3_channel(ZoomCell,ImgSize(1),ImgSize(1),[0 0.3 0.7]);
                            CantHandle = true;
                            while CantHandle == true
                                try
                                    Ypredicted = predict(NeuralNet,X,'MiniBatchSize',obj.MiniBatchSize,'Acceleration','auto',...
                                        'ExecutionEnvironment',obj.NeuralNetAccelerator);
                                    CantHandle = false;
                                catch
                                    obj.CPFlag.CNNopt = 0;
                                    obj.cnn_runtime_optimization(NeuralNet,X);
                                end
                            end
                            waitbar(1,h,'Wrapping up');
                            iRange = find(obj.SelectedCurves);
                            k = 1;
                            for i=iRange'
                                [App,HHApp] = ZoomObj.get_force_curve_data(i,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                                obj.CP_CNNZoom(i,1) = Ypredicted(k,1)*range(HHApp)+min(HHApp);
                                obj.CP_CNNZoom(i,2) = Ypredicted(k,2)*range(App)+min(App);
                                obj.CP(i,1) = obj.CP_CNNZoom(i,1);
                                obj.CP(i,2) = obj.CP_CNNZoom(i,2);
                                k = k + 1;
                            end
                            
                            if BigLoop == NumPartitions
                                obj.CPFlag.CNNZoom = 1;
                            end
                        case 3
                            % ZoomDropout
                            waitbar(1/3,h,'Predicting CP, first guess...');
                            while CantHandle == true
                                try
                                    Ypredicted = predict(NeuralNet,X,'MiniBatchSize',obj.MiniBatchSize,'Acceleration','auto',...
                                        'ExecutionEnvironment',obj.NeuralNetAccelerator);
                                    CantHandle = false;
                                catch
                                    obj.CPFlag.CNNopt = 0;
                                    obj.cnn_runtime_optimization(NeuralNet,X);
                                end
                            end
                            iRange = find(obj.SelectedCurves);
                            k = 1;
                            for i=iRange'
                                [App,HHApp] = obj.get_force_curve_data(i,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                                obj.CP_CNNZoom(i,1) = Ypredicted(k,1)*range(HHApp)+min(HHApp);
                                obj.CP_CNNZoom(i,2) = Ypredicted(k,2)*range(App)+min(App);
                                obj.CP(i,1) = obj.CP_CNNZoom(i,1);
                                obj.CP(i,2) = obj.CP_CNNZoom(i,2);
                                k = k + 1;
                            end
                            
                            waitbar(2/3,h,'Predicting zoomed CP');
                            ZoomObj = obj.copy;
                            ZoomObj.cnn_zoom_in();
                            ZoomCell{1,1} = ZoomObj;
                            X = obj.CP_batchprep_3_channel(ZoomCell,ImgSize(1),ImgSize(1),[0 0.3 0.7]);
                            CantHandle = true;
                            while CantHandle == true
                                try
                                    Ypredicted = predict(NeuralNet,X,'MiniBatchSize',obj.MiniBatchSize,'Acceleration','auto',...
                                        'ExecutionEnvironment',obj.NeuralNetAccelerator);
                                    CantHandle = false;
                                catch
                                    obj.CPFlag.CNNopt = 0;
                                    obj.cnn_runtime_optimization(NeuralNet,X);
                                end
                            end
                            waitbar(1,h,'Wrapping up');
                            iRange = find(obj.SelectedCurves);
                            k = 1;
                            for i=iRange'
                                [App,HHApp] = ZoomObj.get_force_curve_data(i,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                                obj.CP_CNNZoom(i,1) = Ypredicted(k,1)*range(HHApp)+min(HHApp{i});
                                obj.CP_CNNZoom(i,2) = Ypredicted(k,2)*range(App)+min(App);
                                obj.CP(i,1) = obj.CP_CNNZoom(i,1);
                                obj.CP(i,2) = obj.CP_CNNZoom(i,2);
                                k = k + 1;
                            end
                            if BigLoop == NumPartitions
                                obj.CPFlag.CNNZoomDropout = 1;
                            end
                        case 4
                            % ZoomSweep
                            waitbar(1/3,h,sprintf('Predicting CP, first guess...(Partition %i/%i)',BigLoop,NumPartitions));
                            while CantHandle == true
                                try
                                    Ypredicted = predict(NeuralNet,X,'MiniBatchSize',obj.MiniBatchSize,'Acceleration','auto',...
                                        'ExecutionEnvironment',obj.NeuralNetAccelerator);
                                    CantHandle = false;
                                catch
                                    obj.CPFlag.CNNopt = 0;
                                    obj.cnn_runtime_optimization(NeuralNet,X);
                                end
                            end
                            iRange = find(obj.SelectedCurves);
                            k = 1;
                            for i=iRange'
                                [App,HHApp] = obj.get_force_curve_data(i,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                                obj.CP_CNNZoom(i,1) = Ypredicted(k,1)*range(HHApp)+min(HHApp);
                                obj.CP_CNNZoom(i,2) = Ypredicted(k,2)*range(App)+min(App);
                                obj.CP(i,1) = obj.CP_CNNZoom(i,1);
                                obj.CP(i,2) = obj.CP_CNNZoom(i,2);
                                k = k + 1;
                            end
                            obj.CPFlag.CNNZoom = true;
                            
                            waitbar(2/3,h,sprintf('Predicting zoomed CP, sweeping over multiple zooms \n(Partition %i/%i)',BigLoop,NumPartitions));
                            MaxZoom = 0.8;
                            ZoomFactor = 0:MaxZoom/(NumPasses-1):MaxZoom;
                            for i=1:NumPasses
                                ZoomCell{i} = obj.copy;
                                ZoomCell{i}.cnn_zoom_in(ZoomFactor(i));
                            end
                            X = obj.CP_batchprep_3_channel(ZoomCell,ImgSize(1),ImgSize(1),[0 0.3 0.7]);
                            CantHandle = true;
                            while CantHandle == true
                                try
                                    Ypredicted = predict(NeuralNet,X,'MiniBatchSize',obj.MiniBatchSize,'Acceleration','auto',...
                                        'ExecutionEnvironment',obj.NeuralNetAccelerator);
                                    CantHandle = false;
                                catch
                                    obj.CPFlag.CNNopt = 0;
                                    obj.cnn_runtime_optimization(NeuralNet,X);
                                end
                            end
                            waitbar(1,h,'Wrapping up');
                            iRange = find(obj.SelectedCurves);
                            k = 1;
                            for i=iRange'
                                for j=1:NumPasses
                                    [App,HHApp] = ZoomCell{j}.get_force_curve_data(i,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                                    TempCP(i,1,j) = Ypredicted(k+length(iRange)*(j-1),1)*range(HHApp)+min(HHApp);
                                    TempCP(i,2,j) = Ypredicted(k+length(iRange)*(j-1),2)*range(App)+min(App);
                                end
                                obj.CP_CNNZoomSweep(i,1) = mean(TempCP(i,1,:),3);
                                obj.CP_CNNZoomSweep(i,2) = mean(TempCP(i,2,:),3);
                                obj.CP(i,1) = obj.CP_CNNZoomSweep(i,1);
                                obj.CP(i,2) = obj.CP_CNNZoomSweep(i,2);
                                k = k + 1;
                            end
                            if BigLoop == NumPartitions
                                obj.CPFlag.CNNZoomSweep = 1;
                            end
                    end
                    close(h)
                end
                obj.SelectedCurves = obj.StashedSelectedCurves;
            catch ME
                obj.SelectedCurves = obj.StashedSelectedCurves;
                rethrow(ME)
            end
            
            
            %             current = what();
            %             cd(obj.Folder)
            %             savename = sprintf('%s.mat',obj.Name);
            %             save(savename,'obj')
            %             cd(current.path)
        end
        
        function estimate_cp_old(obj)
            % CP estimation using the the approach from the old fibril
            % analysis script 'A_nIAFM_analysis_main.m' but without first
            % subtracting the deflection as it is not needed in
            % Oliver-Pharr analysis
            iRange = find(obj.SelectedCurves);
            for i=iRange'
                [App,HHApp] = obj.get_force_curve_data(i,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                [Ret,HHRet] = obj.get_force_curve_data(i,'AppRetSwitch',1,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                try
                load = zeros(length(App),2);
                unload = zeros(length(Ret),2);
                for j=1:length(App)
                    load(end-(j-1),1) = HHApp(j);
                    load(j,2) = App(j);
                end
                for j=1:length(Ret)
                    unload(end-(j-1),1) = HHRet(j);
                    unload(j,2) = Ret(j);
                end
                [LoadOld{i},UnloadOld{i},Position,vDef] = ContactPoint_sort(load,unload);
                obj.CP(i,2) = App(Position);
                obj.CP(i,1) = HHApp(Position);
                obj.CP_Old(i,1) =obj.CP(i,1);
                obj.CP_Old(i,2) =obj.CP(i,2);
                catch
                    disp(sprintf('Failed to find CP on Curve Nr.%i. \nReplacing with minimum values and unselecting curve',i))
                    obj.CP(i,2) = App(1);
                    obj.CP(i,1) = HHApp(1);
                    obj.CP_Old(i,1) =obj.CP(i,1);
                    obj.CP_Old(i,2) =obj.CP(i,2);
                    obj.SelectedCurves(i) = 0;
                end
            end
            obj.CPFlag.Old = 1;
            %             current = what();
            %             cd(obj.Folder)
            %             savename = sprintf('%s.mat',obj.Name);
            %             save(savename,'obj')
            %             cd(current.path)
        end
        
        function estimate_cp_manually(obj,Zoom,Idx)
            % Manual contact point selection for NN training on plotted force curves
            % returning a position vector in meters.
            % 0 <= Zoom <= 1 determines how much of the
            % rough-estimate-non-contact-domain is cut off from the
            % picture (e.g. Zoom = 1 just shows the rough-estimate-contact-domain)
            
            if nargin < 2
                Zoom = 0.2;
            end
            
            jRange = find(obj.SelectedCurves);
            fig = figure('Name',obj.Name,'Units','normalized','Position',[0.2 0.125 0.8 0.8]);
            k = 1;
            j = jRange(k);
            if nargin == 3
                j = Idx;
            end
            while sum(jRange==j)
                %                 fig.WindowState = 'fullscreen';
                [App,HHApp] = obj.get_force_curve_data(j,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                [Ret,HHRet] = obj.get_force_curve_data(j,'AppRetSwitch',1,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                plot(HHApp,App);
                plottitle = sprintf('Curve Nr.%i/%i\n Click or click and drag the point to the contact point\n Click the red area to exclude the current curve from further analysis\n Click the green area to go back one curve',j,obj.NCurves);
                title(plottitle);
                [~, domainidx] = ForceMap.no_contact_domain(App);
                axis([HHApp(floor(domainidx*Zoom)) inf -inf inf])
                XRange = range(HHApp(floor(domainidx*Zoom):end));
                YRange = range(App);
                BtnSemAxisX = XRange/8;
                BtnSemAxisY = YRange/8;
                XPosDel = max(HHApp) - 5/8*XRange;
                YPosDel = max(App) - 1/3*YRange;
                XPosBack = max(HHApp) - 7/8*XRange;
                YPosBack = max(App) - 1/3*YRange;
                DelBtn = drawellipse('Center',[XPosDel YPosDel],...
                    'SemiAxes',[BtnSemAxisX BtnSemAxisY],...
                    'Color','r',...
                    'Label','Exclude Curve');
                BackBtn = drawellipse('Center',[XPosBack YPosBack],...
                    'SemiAxes',[BtnSemAxisX BtnSemAxisY],...
                    'Color','g',...
                    'Label','Go Back one entry');
                try
                    drawpoint('Position',[obj.Man_CP(j,1) obj.Man_CP(j,2)],'Color','r','Label','Old Man_CP')
                catch
                end
                CP_point = drawpoint();
                
                if ((CP_point.Position(1)-DelBtn.Center(1))/DelBtn.SemiAxes(1))^2 +...
                        ((CP_point.Position(2)-DelBtn.Center(2))/DelBtn.SemiAxes(2))^2 <= 1
                    obj.SelectedCurves(j) = 0;
                    k = k + 1;
                    try
                        j = jRange(k);
                    catch
                        close(fig)
                        return
                    end
                    continue
                elseif ((CP_point.Position(1)-BackBtn.Center(1))/BackBtn.SemiAxes(1))^2 +...
                        ((CP_point.Position(2)-BackBtn.Center(2))/BackBtn.SemiAxes(2))^2 <= 1
                    CurrentIndex = find(jRange==j);
                    if CurrentIndex == 1
                        continue
                    else
                        j = jRange(CurrentIndex - 1);
                        continue
                    end
                end
                
                obj.Man_CP(j,1) = CP_point.Position(1);
                obj.Man_CP(j,2) = CP_point.Position(2);
                k = k + 1;
                try
                    j = jRange(k);
                catch
                    close(fig)
                    return
                end
            end
            close(fig);
            obj.CPFlag.Manual = 1;
            %             current = what();
            %             cd(obj.Folder)
            %             savename = sprintf('%s.mat',obj.Name);
            %             save(savename,'obj')
            %             cd(current.path)
        end
        
        function estimate_cp_hardsurface(obj)
            % contact point estimation for force curves detected on hard
            % surfaces
            
            for i=1:obj.NCurves
                [App,HHApp] = obj.get_force_curve_data(i,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                obj.CP_HardSurface(i,1) = HHApp(end) - App(end)/obj.SpringConstant;
                obj.CP_HardSurface(i,2) = 0;
                %% Debugging
                % plot(HHApp,App);
                % drawpoint('Position',[obj.CP_HardSurface(i,1) obj.CP_HardSurface(i,2)]);
            end
            obj.CP = obj.CP_HardSurface;
            obj.CPFlag.HardSurface = 1;
        end
        
        function estimate_cp_curve_origin(obj)
            % contact point estimation for force curves assuming the curve
            % starts at contact point (this is useful in very
            % range-restricted hertz fits)
            
            for i=1:obj.NCurves
                [App,HHApp] = obj.get_force_curve_data(i,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                obj.CP_CurveOrigin(i,1) = HHApp(1);
                obj.CP_CurveOrigin(i,2) = App(1);
                %% Debugging
                % plot(HHApp,App);
                % drawpoint('Position',[obj.CP_HardSurface(i,1) obj.CP_HardSurface(i,2)]);
            end
            obj.CP = obj.CP_CurveOrigin;
            obj.CPFlag.CurveOrigin = 1;
        end
        
        function E = calculate_e_mod_hertz(obj,CPType,TipShape,LowerCurveFraction,...
                UpperCurveFraction,AllowXShift,CorrectSensitivity,UseTipData,...
                UseTopography,TipObject,DataWeightByDistanceBool,...
                SortHeightDataForFit,FitDegreeForSneddonPolySurf,...
                LowerForceCutoff,UpperForceCutoff,...
                SensitivityCorrectionMethod,...
                KeepOldResults,TopographyHeightChannelName,...
                ThinFilmMode,ThinFilmThickness)
%                 E = calculate_e_mod_hertz(obj,CPType,TipShape,LowerCurveFraction,...
%                 UpperCurveFraction,AllowXShift,CorrectSensitivity,UseTipData,...
%                 UseTopography,TipObject,DataWeightByDistanceBool,...
%                 SortHeightDataForFit,FitDegreeForSneddonPolySurf,...
%                 LowerForceCutoff,UpperForceCutoff)

            if nargin < 5
                AllowXShift = false;
                CorrectSensitivity = false;
                UseTipData = false;
                KeepOldResults = false;
            end
            
            obj.HertzFitSSE = zeros(1,obj.NCurves);
            obj.HertzFitRSquare = zeros(1,obj.NCurves);
            obj.HertzFitDFE = zeros(1,obj.NCurves);
            obj.HertzFitAdjRSquare = zeros(1,obj.NCurves);
            obj.HertzFitRMSE = zeros(1,obj.NCurves);
            
            if contains(TipShape,'on thin film') && isequal(ThinFilmMode,'From Topography')
                TopographyHeightChannel = obj.get_channel(TopographyHeightChannelName);
                if isempty(TopographyHeightChannel)
                    error('Channel %s not found',TopographyHeightChannelName);
                end
                TopographyList = obj.convert_map_to_data_list(TopographyHeightChannel.Image);
            end
            
            if UseTopography
                [~, AllMatches] = obj.search_channel('Radius of Curvature');
                TempChan = [];
                for i = 1:length(AllMatches)
                    TempChan = obj.get_channel(AllMatches{i});
                    if isequal(TempChan.FMA_ID, obj.CurrentFMA_ID)
                        break;
                    end
                end
                
                if isempty(TempChan) || ~isequal(TempChan.FMA_ID, obj.CurrentFMA_ID)
                    warning(sprintf('No channel with FMA_ID matching CurrentFMA_ID found. Using the first match "%s" instead.',AllMatches{1}));
                    TempChan = obj.get_channel(AllMatches{1});
                end
                
                TopographyCurvatureChannel = TempChan;
                
                if isempty(TopographyCurvatureChannel)
                    error('Channel %s not found', 'Local Radius of Curvature');
                end
                
                LROCList = obj.convert_map_to_data_list(TopographyCurvatureChannel.Image);
            end

            
            if isequal(lower(TipShape),'parabolic') || isequal(lower(TipShape),'spheric approx.')
            elseif isequal(lower(TipShape),'parabolic on thin film (not bonded)')
                nu = obj.PoissonR;
                alpha = -0.347*(3-2*nu)/(1-nu);
                beta = 0.056*(5-2*nu)/(1-nu);
            elseif isequal(lower(TipShape),'parabolic on thin film (bonded)')
                nu = obj.PoissonR;
                alpha = -(1.2876-1.4678*nu+1.3442*nu^2)/(1-nu);
                beta = (0.6387-1.0277*nu+1.5164*nu^2)/(1-nu);
            elseif isequal(lower(TipShape,'sneddonpolysurf'))
                Area = TipObject.DepthDependendTipRadius;
                Radius = sqrt(Area/pi);
                Depth = [1:length(Radius)]'.*1e-9;
                RangeR = range(Radius);
                RangeD = range(Depth);
                X = Radius/RangeR;
                Y = Depth/RangeD;
                TypeString = '';
                for i=1:FitDegreeForSneddonPolySurf
                    TypeString = [TypeString '+c' num2str(i,'%02u') '*x^' num2str(i)];
                end
                FitType = fittype(TypeString);
                FitOpts = fitoptions('Method','LinearLeastSquares');
                FitObject = fit(X,Y,FitType,FitOpts);
            end
            
            iRange = find(obj.SelectedCurves);
            obj.EModHertz = zeros(obj.NCurves,1);
            obj.IndentationDepth = zeros(obj.NCurves,1);
            obj.IndentationDepthHertz = zeros(obj.NCurves,1);
            obj.IndentationDepthHertzFitRange = zeros(obj.NCurves,1);
            EffectiveRadius = ones(obj.NCurves,1).*NaN;
            
            c = parcluster();
            NumWorkers = c.NumWorkers;
            
            while ~isempty(iRange')
                BatchSize = min(NumWorkers,length(iRange));
                if isequal(lower(CPType),'cnn')
                    CP = obj.CP(iRange(1:BatchSize),:);
                elseif isequal(lower(CPType),'old')
                    CP = obj.CP_Old(iRange(1:BatchSize),:);
                elseif isequal(lower(CPType),'rov')
                    CP = obj.CP_RoV(iRange(1:BatchSize),:);
                elseif isequal(lower(CPType),'gof')
                    CP = obj.CP_GoF(iRange(1:BatchSize),:);
                elseif isequal(lower(CPType),'combo')
                    CP = obj.CP_Combo(iRange(1:BatchSize),:);
                elseif isequal(lower(CPType),'manual')
                    CP = obj.Man_CP(iRange(1:BatchSize),:);
                elseif isequal(lower(CPType),'snap-in')
                    CP = obj.CP_SnapIn(iRange(1:BatchSize),:);
                else
                    CP = obj.CP(iRange(1:BatchSize),:);
                end
                for i=1:BatchSize
                    
                    if CorrectSensitivity
                        if isequal(SensitivityCorrectionMethod,'Adaptive')
                            [App{i},HHApp{i}] = obj.get_force_curve_data(iRange(i),...
                                'AppRetSwitch',0,...
                                'BaselineCorrection',1,...
                                'TipHeightCorrection',0,...
                                'Sensitivity','adaptive',...
                                'Unit','N');
                        else
                            [App{i},HHApp{i}] = obj.get_force_curve_data(iRange(i),'AppRetSwitch',0,...
                                'BaselineCorrection',1,'TipHeightCorrection',0,...
                                'Sensitivity','corrected',...
                                'Unit','N');
                        end
                    else
                        [App{i},HHApp{i}] = obj.get_force_curve_data(iRange(i),'AppRetSwitch',0,...
                            'BaselineCorrection',1,'TipHeightCorrection',0,...
                            'Sensitivity','original',...
                            'Unit','N');
                    end
                    if SortHeightDataForFit
                        HHApp{i} = sort(HHApp{i},'ascend');
                    end
                    force{i} = App{i} - CP(i,2);
                    tip_h{i} = (HHApp{i} - CP(i,1)) - force{i}/obj.SpringConstant;
                    if ~isequal(CPType,'None')
                        tip_h{i}(tip_h{i} < 0) = [];
                    end
                    force{i}(1:(length(force{i})-length(tip_h{i}))) = [];
                    if length(tip_h{i}) < 2
                        continue
                    end
                    Max{i} = max(tip_h{i});
                    obj.IndentationDepth(iRange(i)) = Max{i}(1);
                    
                    % Apply absolute force cut offs
                    if ~isempty(UpperForceCutoff)
                        tip_h{i}(force{i}>UpperForceCutoff) = [];
                        force{i}(force{i}>UpperForceCutoff) = [];
                    end
                    if ~isempty(LowerForceCutoff)
                        tip_h{i}(force{i}<LowerForceCutoff) = [];
                        force{i}(force{i}<LowerForceCutoff) = [];
                    end
                    % delete everything below curve_percent of the maximum
                    % force
                    force{i}(force{i}<(LowerCurveFraction)*max(force{i})) = [];
                    force{i}(force{i}>(UpperCurveFraction)*max(force{i})) = [];
                    tip_h{i}(1:(length(tip_h{i})-length(force{i}))) = [];
                    
                    if length(tip_h{i}) < 2
                        continue
                    end
                    % Allocate the maximum indentation for later
                    % calculation of max indent. depth of fitrange
                    MaxFitRange{i} = max(tip_h{i});
                    RangeF{i} = range(force{i});
                    RangeTH{i} = range(tip_h{i});
                    force{i} = force{i}/RangeF{i};
                    tip_h{i} = tip_h{i}/RangeTH{i};
                    
                    if contains(TipShape,'on thin film')
                        if isequal(ThinFilmMode,'FromValue')
                            FilmHeight{i} = ThinFilmThickness;
                        elseif isequal(ThinFilmMode,'From Topography')
                            FilmHeight{i} = max(TopographyList(iRange(i)),MaxFitRange{i}(1));
                        end
                    end
                    
                    if isempty(obj.FibDiam) || UseTopography
                        if UseTipData
                            DepthIndex = floor(obj.IndentationDepth(iRange(i))*1e9);
                            DepthRemainder = obj.IndentationDepth(iRange(i))*1e9 - DepthIndex;
                            if DepthIndex >= length(TipObject.DepthDependendTipRadius)
                                DepthIndex = length(TipObject.DepthDependendTipRadius) - 1;
                                DepthRemainder = 0;
                            end
                            if DepthIndex == 0
                                TipRadius = TipObject.DepthDependendTipRadius(DepthIndex+1)*DepthRemainder;
                            else
                                TipRadius = TipObject.DepthDependendTipRadius(DepthIndex)*(1-DepthRemainder) + TipObject.DepthDependendTipRadius(DepthIndex+1)*DepthRemainder;
                            end
                            if UseTopography
                                RTopo = -LROCList(iRange(i));
                                
                                R_eff{i} =1/(1/TipRadius + 1/RTopo);
                            else
                                R_eff{i} = TipRadius;
                            end
                        else
                            if UseTopography
                                RTip = obj.TipRadius;
                                
                                % Negative LROCList because of how the
                                % ROC is defined for the local surfface
                                % fits
                                RTopo = -LROCList(iRange(i));
                                
                                R_eff{i} = 1/(1/RTip + 1/RTopo);
                            else
                                R_eff{i} = obj.TipRadius;
                            end
                        end
                    else
                        if UseTipData
                            DepthIndex = floor(obj.IndentationDepth(iRange(i))*1e9);
                            DepthRemainder = obj.IndentationDepth(iRange(i))*1e9 - DepthIndex;
                            if DepthIndex >= length(TipObject.DepthDependendTipRadius)
                                DepthIndex = length(TipObject.DepthDependendTipRadius) - 1;
                            end
                            if DepthIndex == 0
                                TipRadius = TipObject.DepthDependendTipRadius(DepthIndex+1)*DepthRemainder;
                            else
                                TipRadius = TipObject.DepthDependendTipRadius(DepthIndex)*(1-DepthRemainder) + TipObject.DepthDependendTipRadius(DepthIndex+1)*DepthRemainder;
                            end
                            R_eff{i} = 1/(1/TipRadius + 1/(obj.FibDiam/2));
                        else
                            R_eff{i} = 1/(1/(obj.TipRadius) + 1/(obj.FibDiam/2));
                        end
                    end
                    EffectiveRadius(iRange(i)) = R_eff{i};
                    if isequal(lower(TipShape),'parabolic') || isequal(lower(TipShape),'spheric approx.')
                        if AllowXShift
                            FitFunction{i} = 'a*(x+b)^(3/2)';
                        else
                            FitFunction{i} = 'a*(x)^(3/2)';
                        end
                    elseif contains(TipShape,'parabolic on thin film')
                        FitR_eff{i} = R_eff{i}/RangeTH{i};
                        FitFilmHeight{i} = FilmHeight{i}/RangeTH{i};
                        c1 = -2*alpha/pi*sqrt(FitR_eff{i})/FitFilmHeight{i};
                        c2 = 4*alpha^2/pi^2*(sqrt(FitR_eff{i})/FitFilmHeight{i})^2;
                        c3 = -8/pi^3*(alpha^3 + 4*pi^2/15*beta)*(sqrt(FitR_eff{i})/FitFilmHeight{i})^3;
                        c4 = 16*alpha/pi^4*(alpha^3 + 3*pi^2/5*beta)*(sqrt(FitR_eff{i})/FitFilmHeight{i})^4;
                        if AllowXShift
                            FitFunction{i} = sprintf(...
                                'a*(x-b)^(3/2)*(1+%6e*(x-b)^(1/2)+%6e*(x-b)^(2/2)+%6e*(x-b)^(3/2)+%6e*(x-b)^(4/2))',...
                                c1,c2,c3,c4);
                        else
                            FitFunction{i} = sprintf(...
                                'a*x^(3/2)*(1+%6e*x^(1/2)+%6e*x^(2/2)+%6e*x^(3/2)+%6e*x^(4/2))',...
                                c1,c2,c3,c4);
                        end
                    end
                    
                    if DataWeightByDistanceBool
                        Points = [tip_h{i}(1:end-1) force{i}(1:end-1)];
                        ShiftedPoints = [tip_h{i}(2:end) force{i}(2:end)];
                        TempPointWeight = vecnorm(Points-ShiftedPoints,2,2);
                        PointWeights{i} = [TempPointWeight(1) ; TempPointWeight];
                    else
                        PointWeights{i} = ones(length(force{i}),1);
                    end
                end
                parfor i=1:BatchSize
                    try
                        if AllowXShift
                            s = fitoptions('Method','NonlinearLeastSquares',...
                                'Lower',[10^(-5) -min(tip_h{i})],...
                                'Upper',[inf min(tip_h{i})],...
                                'MaxIter',100,...
                                'Startpoint',[1 0],...
                                'Weights',PointWeights{i});
                            f = fittype(FitFunction{i},'options',s);
                        else
                            s = fitoptions('Method','NonlinearLeastSquares',...
                                'Lower',10^(-5),...
                                'Upper',inf,...
                                'Startpoint',1,...
                                'Weights',PointWeights{i});
                            f = fittype(FitFunction{i},'options',s);
                        end
                        [Hertzfit{i},GoF{i}] = fit(tip_h{i},...
                            force{i},f);
                    catch
                        Hertzfit{i} = nan;
                    end
                end
                for i=1:BatchSize
                    try
                        % calculate E module based on the Hertz model. Be careful
                        % to convert to unnormalized data again
                        EMod{i} = 3*(Hertzfit{i}.a*RangeF{i}/RangeTH{i}^(3/2))/(4*sqrt(R_eff{i}))*(1-obj.PoissonR^2);
                    catch ME
                        EMod{i} = nan;
%                         Hertzfit{i}.a = 0;
%                         if AllowXShift
%                             Hertzfit{i}.b = 0;
%                         end
                    end
                    obj.EModHertz(iRange(i)) = EMod{i};
                    % Convert the model to the right scale so it can be plotted
                    % correctly later
                    warning('off','all');
                    try
                        Hertzfit{i}.a = Hertzfit{i}.a*RangeF{i}/RangeTH{i}.^(3/2);
                        if contains(TipShape,'parabolic on thin film')
                            % Regular expression pattern to match constants
                            if AllowXShift
                                pattern = '\+([-\d\.eE]+)\*\(x-b\)\^\([\d/]+\)';
                            else
                                pattern = '\+([-\d\.eE]+)\*x\^\([\d/]+\)';
                            end
                            func_str = formula(Hertzfit{i});
                            
                            % Extract constants using regexp
                            tokens = regexp(func_str, pattern, 'tokens');
                            
                            % Convert tokens to numerical values
                            constants = cellfun(@(c) str2double(c{1}), tokens);
                            
                            % Rescaling factors for each constant
                            x_range_powers = [0.5, 1, 1.5, 2]; % Corresponding powers of x_range
                            
                            % Rescale the constants
                            rescaled_constants = constants ./ (RangeTH{i}.^ x_range_powers);
                            
                            % Create a copy of the original function string to modify
                            new_func_str = func_str;
                            
                            % Prepare formatted strings of rescaled constants
                            formatted_constants = arrayfun(@(c) sprintf('%.8e', c), rescaled_constants, 'UniformOutput', false);
                            
                            % Regular expression pattern to match constants in the original string
                            if AllowXShift
                                pattern_replace = '(\+)([-\d\.eE]+)(\*\(x-b\)\^\([\d/]+\))';
                            else
                                pattern_replace = '(\+)([-\d\.eE]+)(\*x\^\([\d/]+\))';
                            end
                            % Use regexp to find positions of constants to replace
                            [match_starts, match_ends] = regexp(new_func_str, pattern_replace);
                            
                            % Loop over each match and replace the constant
                            for j = 1:length(match_starts)
                                % Extract the full matched substring
                                original_substring = new_func_str(match_starts(j):match_ends(j));
                                
                                % Extract the parts of the substring
                                parts = regexp(original_substring, pattern_replace, 'tokens');
                                parts = parts{1};
                                
                                % Reassemble the substring with the rescaled constant
                                new_substring = [parts{1}, formatted_constants{j}, parts{3}];
                                
                                % Replace the original substring with the new substring
                                new_func_str = [new_func_str(1:match_starts(j)-1), new_substring, new_func_str(match_ends(j)+1:end)];
                                
                                % Adjust positions for next iteration due to change in string length
                                offset = length(new_substring) - length(original_substring);
                                match_starts = match_starts + offset;
                                match_ends = match_ends + offset;
                            end
                        else
                            new_func_str = formula(Hertzfit{i});
                        end
                        if AllowXShift
                            Hertzfit{i}.b = Hertzfit{i}.b*RangeTH{i};
                            obj.CP_HertzFitted(iRange(i),1) = CP(i,1)-Hertzfit{i}.b;
                            obj.CP_HertzFitted(iRange(i),2) = CP(i,2);
                            % Not sure about this one
                            % obj.IndentationDepth(i) = obj.IndentationDepth(i) + Hertzfit.b;
                            obj.IndentationDepthHertz(iRange(i)) = Max{i}(1)+Hertzfit{i}.b;
                            obj.IndentationDepthHertzFitRange(iRange(i)) = MaxFitRange{i}(1)+Hertzfit{i}.b;
                        end
                        warning('on','all');
                        %                         obj.HertzFit{iRange(i)} = Hertzfit{i};
                        try
                            obj.HertzFitType{iRange(i)} = new_func_str;
                            obj.HertzFitCoeffNames = coeffnames(Hertzfit{i});
                        catch
                        end
                        obj.HertzFitValues{iRange(i)} = coeffvalues(Hertzfit{i});
                        obj.HertzFitSSE(iRange(i)) = GoF{i}.sse;
                        obj.HertzFitRSquare(iRange(i)) = GoF{i}.rsquare;
                        obj.HertzFitDFE(iRange(i)) = GoF{i}.dfe;
                        obj.HertzFitAdjRSquare(iRange(i)) = GoF{i}.adjrsquare;
                        obj.HertzFitRMSE(iRange(i)) = GoF{i}.rmse;
                    catch
                        obj.SelectedCurves(iRange(i)) = 0;
                        warning('on','all');
                    end
                end
                iRange(1:BatchSize) = [];
            end
            E = obj.EModHertz;
            if isequal(lower(CPType),'cnn')
                obj.EModHertz_CNN = E;
            elseif isequal(lower(CPType),'old')
                obj.EModHertz_Old = E;
            elseif isequal(lower(CPType),'rov')
                obj.EModHertz_RoV = E;
            end
            for i=1:obj.NumPixelsX
                for j=1:obj.NumPixelsY
                    obj.EModMapHertz(i,j,1) = obj.EModHertz(obj.Map2List(i,j));
                    IndDepMap(i,j) = obj.IndentationDepth(obj.Map2List(i,j));
                    IndDepMapHertz(i,j) = obj.IndentationDepthHertz(obj.Map2List(i,j));
                    IndDepMapHertzFitRange(i,j) = obj.IndentationDepthHertzFitRange(obj.Map2List(i,j));
                    ERMap(i,j) = EffectiveRadius(obj.Map2List(i,j)); 
                end
            end
            
            % Write to Channel
            EModMap = obj.create_standard_channel(obj.EModMapHertz(:,:,1),'Indentation Modulus Hertz','Pa');
            obj.add_channel(EModMap,~KeepOldResults)
            
            EModLog = obj.create_standard_channel(log(obj.EModMapHertz(:,:,1)),'Indentation Modulus Hertz (log)','Pa');
            obj.add_channel(EModLog,~KeepOldResults)
            
            Channel = obj.create_standard_channel(IndDepMap,'Indentation Depth','m');
            obj.add_channel(Channel,~KeepOldResults)
            
            Channel = obj.create_standard_channel(IndDepMapHertz,'Indentation Depth Hertz','m');
            obj.add_channel(Channel,~KeepOldResults)
            
            Channel = obj.create_standard_channel(IndDepMapHertzFitRange,'Indentation Depth Hertz Fit Range','m');
            obj.add_channel(Channel,~KeepOldResults)
            
            RSquareMap = obj.convert_data_list_to_map(obj.HertzFitRSquare);
            RSquare = obj.create_standard_channel(RSquareMap,'Hertz Fit RSquare','');
            obj.add_channel(RSquare,~KeepOldResults)
            
            Channel = obj.create_standard_channel(ERMap,'Effective Radius Hertz','m');
            obj.add_channel(Channel,~KeepOldResults)
            
            if AllowXShift
                obj.CPFlag.HertzFitted = 1;
            end
            
            % Write current fit parameters to HertzFitStore
            if isempty(obj.HertzFitStore)
                obj.HertzFitStore(1).FMA_ID = obj.CurrentFMA_ID;
            else
                obj.HertzFitStore(end + 1).FMA_ID = obj.CurrentFMA_ID;
            end
            obj.HertzFitStore(end).HertzFitType = obj.HertzFitType;
            obj.HertzFitStore(end).HertzFitCoeffNames = obj.HertzFitCoeffNames;
            obj.HertzFitStore(end).HertzFitValues = obj.HertzFitValues;
            
        end
        
        function calculate_predictive_rsquare_hertz(obj,CorrectSensitivity,...
                SensitivityCorrectionMethod,AllowXShift,SortHeightData,KeepOldResults)
            
            Range = find(obj.SelectedCurves);
            PredictiveRSquare = ones(obj.NCurves,1).*NaN;
            
            if ~iscell(obj.HertzFitType)
                FitFunction = fittype(obj.HertzFitType);
            end
            
            for i=Range'
                if iscell(obj.HertzFitType)
                    FitFunction = fittype(obj.HertzFitType{i});
                end
                if CorrectSensitivity
                    if isequal(SensitivityCorrectionMethod,'Adaptive')
                        [App,HHApp] = obj.get_force_curve_data(i,...
                            'AppRetSwitch',0,...
                            'BaselineCorrection',1,...
                            'TipHeightCorrection',1,...
                            'Sensitivity','adaptive',...
                            'Unit','N');
                    else
                        [App,HHApp] = obj.get_force_curve_data(i,'AppRetSwitch',0,...
                            'BaselineCorrection',1,'TipHeightCorrection',1,...
                            'Sensitivity','corrected',...
                            'Unit','N');
                    end
                else
                    [App,HHApp] = obj.get_force_curve_data(i,'AppRetSwitch',0,...
                        'BaselineCorrection',1,'TipHeightCorrection',1,...
                        'Sensitivity','original',...
                        'Unit','N');
                end
                try
                    % Trim everything below CP_Hertz
                    FitCoeffValues = obj.HertzFitValues{i};
                    if length(FitCoeffValues) == 1
                        Fit = cfit(FitFunction,FitCoeffValues);
                        FitCoeffValues(2) = 0;
                    elseif length(FitCoeffValues) == 2
                        Fit = cfit(FitFunction,FitCoeffValues(1),FitCoeffValues(2));
                    end
                    if SortHeightData
                        HHApp = sort(HHApp,'ascend');
                    end
                    X = HHApp - obj.CP(i,1);
                    Y = App - obj.CP(i,2);
                    X(X<0-FitCoeffValues(2)) = [];
                    Y(1:end-length(X)) = [];
                    YFitted = feval(Fit,X);
                    % Now calculate the RSquare
                    YMean = mean(Y);
                    L = length(Y);
                    SSTot = sum((Y - YMean).^2);
                    SSRes = sum((Y - YFitted).^2);
                    PredictiveRSquare(i) = 1 - SSRes./SSTot;
                catch
                    PredictiveRSquare(i) = NaN;
                end
                %Debug
%                 plot(X,Y,'bx',X,YFitted,'r-')
%                 title(sprintf('R-Square = %f   SSRes = %e   SSTot = %e',PredictiveRSquare(i),SSRes,SSTot))
%                 drawnow
            end
            
            obj.HertzFitPredictiveRSquare = PredictiveRSquare;
            
            PredictiveRSquareMap = obj.convert_data_list_to_map(obj.HertzFitPredictiveRSquare);
            PredictiveRSquareChannel = obj.create_standard_channel(PredictiveRSquareMap,'Hertz Fit Predictive RSquare','-');
            obj.add_channel(PredictiveRSquareChannel,~KeepOldResults)
        end
        
        function calculate_indentation_depth_from_chosen_cp(obj,CPVector)
            
            
            iRange = find(obj.SelectedCurves);
            IndentationDepth = zeros(obj.NCurves,1);
            for i=iRange'
                [App{i},HHApp{i}] = obj.get_force_curve_data(iRange(i),'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                force{i} = App{i} - CP(i,2);
                if CorrectSensitivity
                    force{i} = force{i}.*obj.RefSlopeCorrectedSensitivity/obj.Sensitivity;
                end
                tip_h{i} = (HHApp{i} - CP(i,1)) - force{i}/obj.SpringConstant;
                tip_h{i}(tip_h{i} < 0) = [];
                if length(tip_h{i}) < 2
                    continue
                end
                Max{i} = max(tip_h{i});
                IndentationDepth(iRange(i)) = Max{i}(1);
            end
                    
            obj.IndentationDepth = IndentationDepth;
        end
        
        function EMod = calculate_e_mod_oliverpharr(obj,TipProjArea,CurvePercent,KeepOldResults)
            
            if nargin < 3
                CurvePercent = 0.75;
            end
            
            Range = find(obj.SelectedCurves);
            Epsilon = 0.73; % Correction constant from Oliver Pharr Method (1992)
            Beta = 1.0226; %Correction constant from Oliver Pharr Method (1992)
            obj.EModOliverPharr = zeros(obj.NCurves,1);
            EMod = zeros(obj.NCurves,1);
            obj.ProjTipArea = TipProjArea;
            obj.DZslope = zeros(obj.NCurves,1);
            obj.Stiffness = zeros(obj.NCurves,1);
            obj.IndentationDepth = zeros(obj.NCurves,1);
            obj.IndentationDepthOliverPharr = zeros(obj.NCurves,1);
            obj.IndentArea = zeros(obj.NCurves,1);
            for i=Range'
                [Ret,HHRet] = obj.get_force_curve_data(i,'AppRetSwitch',1,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                Z = HHRet - obj.CP(i,1);
                D = (Ret - obj.CP(i,2))/obj.SpringConstant;
                Zmax(i) = max(Z);
                Dmax(i) = max(D);
                DCurvePercent = D(D>=(1-CurvePercent)*Dmax(i));
                ZCurvePercent = Z(1:length(DCurvePercent));
                LineFit = polyfit(ZCurvePercent,DCurvePercent,1);
                obj.DZslope(i) = LineFit(1);
                Hmax(i) = Zmax(i) - Dmax(i);
                dD = Dmax(i) - (1-CurvePercent)*Dmax(i);
                dh = dD*(1./obj.DZslope(i) - 1/(obj.RefSlope));
                df = dD*obj.SpringConstant;
                obj.Stiffness(i) = df/dh;
                Fmax(i) = Dmax(i).*obj.SpringConstant;
                obj.IndentationDepth(i) = Hmax(i);
                obj.IndentationDepthOliverPharr(i) = Hmax(i) - Epsilon.*Fmax(i)/obj.Stiffness(i);
                % IndentArea is taken as the linear interpolation between
                % the two numeric values of TipProjArea the Hc(i) falls
                % inbetween
                try
                    obj.IndentArea(i) = ((1-(obj.IndentationDepthOliverPharr(i)*1e9-floor(obj.IndentationDepthOliverPharr(i)*1e9)))*TipProjArea(floor(obj.IndentationDepthOliverPharr(i)*1e9))...
                        + (obj.IndentationDepthOliverPharr(i)*1e9-floor(obj.IndentationDepthOliverPharr(i)*1e9))*TipProjArea(ceil(obj.IndentationDepthOliverPharr(i)*1e9)));
                    EMod(i) = sqrt(pi/obj.IndentArea(i))*1/Beta*...
                        obj.Stiffness(i)/2*(1-obj.PoissonR^2);
                    if EMod(i) <= 0
                        EMod(i) = NaN;
                    end
                catch
                    EMod(i) = NaN;
                end
            end
            
            obj.EModOliverPharr = EMod;
            % Write values into EModMapOliverPharr
            
            for i=1:obj.NumPixelsX
                for j=1:obj.NumPixelsY
                    obj.EModMapOliverPharr(i,j,1) = obj.EModOliverPharr(obj.Map2List(i,j));
                    IndDepthMapOP(i,j) = obj.IndentationDepthOliverPharr(obj.Map2List(i,j));
                    DZMap(i,j) = obj.DZslope(obj.Map2List(i,j));
                    IndDepthMap(i,j) = obj.IndentationDepth(obj.Map2List(i,j));
                end
            end
            
            % Write to Channel
            EModMap = obj.create_standard_channel(obj.EModMapOliverPharr(:,:,1),'Indentation Modulus OliverPharr','Pa');
            obj.add_channel(EModMap,~KeepOldResults)
            
            EModLog = obj.create_standard_channel(log(obj.EModMapOliverPharr(:,:,1)),'Indentation Modulus OliverPharr (log)','Pa');
            obj.add_channel(EModLog,~KeepOldResults)
            
            Channel = obj.create_standard_channel(IndDepthMapOP,'Indentation Depth Oliver-Pharr','m');
            obj.add_channel(Channel,~KeepOldResults)
            
            Channel = obj.create_standard_channel(IndDepthMap,'Indentation Depth','m');
            obj.add_channel(Channel,~KeepOldResults)
            
            Channel = obj.create_standard_channel(DZMap,'DZ-Slope','m/m');
            obj.add_channel(Channel,~KeepOldResults)
        end
        
        function MaxAdhesionForce = calculate_adhesion_force(obj,KeepOldResults)
            
            Range = find(obj.SelectedCurves);
            MaxAdhesionForce = zeros(obj.NCurves,1);
            for i=Range'
                [Ret,HHRet] = obj.get_force_curve_data(i,'AppRetSwitch',1,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                Force = Ret - obj.CP(i,2);
                TipHeight = (HHRet - obj.CP(i,1)) - Force/obj.SpringConstant;
                MaxAdhesionForce(i) = -min(Force);
            end
            
            obj.MaxAdhesionForce = MaxAdhesionForce;
            
            
            % Write to Channel
            MaxAdhesionForceMap = obj.convert_data_list_to_map(obj.MaxAdhesionForce);
            Channel = obj.create_standard_channel(MaxAdhesionForceMap,'Maximum Adhesion Force','N');
            obj.add_channel(Channel,~KeepOldResults)
        end
        
        function [AdhesionEnergy,AdhesionLength] = calculate_adhesion_energy_and_length(obj,ThresholdMult,KeepOldResults)
            
            if nargin < 2
                ThresholdMult = 2;
            end
            
            Range = find(obj.SelectedCurves);
            AdhesionEnergy = zeros(obj.NCurves,1);
            AdhesionLength = zeros(obj.NCurves,1);
            for i=Range'
                [App,HHApp] = obj.get_force_curve_data(i,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                [Ret,HHRet] = obj.get_force_curve_data(i,'AppRetSwitch',1,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                AppForce = App - obj.CP(i,2);
                AppTipHeight = (HHApp - obj.CP(i,1)) - AppForce/obj.SpringConstant;
                RetForce = flipud(Ret) - obj.CP(i,2);
                RetTipHeight = (flipud(HHRet) - obj.CP(i,1)) - RetForce/obj.SpringConstant;
                
                try
                    % Baseline the retraction curve
                    NoContact = ForceMap.no_contact_domain(RetForce);
                    RetFit = polyfit(1:length(NoContact),NoContact,1);
                    RetForce = RetForce - polyval(RetFit,1:length(RetForce))';
                    
                    NoContactSTD = std(NoContact);
                    AdhesionX = RetTipHeight(RetForce<-ThresholdMult*NoContactSTD);
                    AdhesionY = RetTipHeight(RetForce<-ThresholdMult*NoContactSTD);
                    
                    AdhesionLength(i) = -min(AdhesionX);
                    
                    SummedEnergy = 0;
                    AdhesionBool = RetForce<-ThresholdMult*NoContactSTD;
                    k = 1;
                    while k <= length(AdhesionBool)
                        if AdhesionBool(k)
                            StartIndex = k;
                            while AdhesionBool(k)
                                k = k + 1;
                            end
                            EndIndex = k - 1;
                            if EndIndex-StartIndex < 2
                                continue
                            end
                            SummedEnergy = SummedEnergy + abs(trapz(RetTipHeight(StartIndex:EndIndex),RetForce(StartIndex:EndIndex)));
                        else
                            k = k + 1;
                        end
                    end
                    AdhesionEnergy(i) = SummedEnergy;
                catch
                    AdhesionEnergy(i) = nan;
                    AdhesionLength(i) = nan;
                end
                
            end
            
            obj.AdhesionEnergy = AdhesionEnergy;
            obj.AdhesionLength = AdhesionLength;
            
            
            % Write to Channel
            AEMap = obj.convert_data_list_to_map(obj.AdhesionEnergy);
            Channel1 = obj.create_standard_channel(AEMap,'Adhesion Energy','J');
            obj.add_channel(Channel1,~KeepOldResults)
            
            CoulombConstant = 1.602176634e-19;
            Channel1 = obj.create_standard_channel(AEMap./CoulombConstant,'Adhesion eV-Energy','eV');
            obj.add_channel(Channel1,~KeepOldResults)
            
            ALMap = obj.convert_data_list_to_map(obj.AdhesionLength);
            Channel2 = obj.create_standard_channel(ALMap,'Adhesion Length','m');
            obj.add_channel(Channel2,~KeepOldResults)
        end
        
        function PeakIndentationAngle = calculate_peak_indentation_angle(obj,FitPortion,KeepOldResults)
            
            if nargin < 2
                FitPortion = .5;
            end
            
            Range = find(obj.SelectedCurves);
            PeakIndentationAngle = zeros(obj.NCurves,1);
            for i=Range'
                [App,HHApp] = obj.get_force_curve_data(i,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                [Ret,HHRet] = obj.get_force_curve_data(i,'AppRetSwitch',1,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                AppForce = App - obj.CP(i,2);
                AppHeadHeight = (HHApp - obj.CP(i,1));% - AppForce/obj.SpringConstant;
                RetForce = Ret - obj.CP(i,2);
                RetHeadHeight = (HHRet - obj.CP(i,1));% - RetForce/obj.SpringConstant;
                
                % Convert from force to vDef
                AppForce = AppForce/obj.SpringConstant;
                RetForce = RetForce/obj.SpringConstant;
                try
                    CutOffForceApp = (1-FitPortion)*AppForce(end);
                    CutOffForceRet = (1-FitPortion)*RetForce(1);
                    k = 0;
                    while AppForce(end-k) > CutOffForceApp
                        k = k + 1;
                    end
                    AppY = AppForce(end-(k+1):end);
                    AppX = AppHeadHeight(end-(k+1):end);
                    k = 1;
                    while RetForce(k) > CutOffForceRet
                        k = k + 1;
                    end
                    RetY = RetForce(1:k-1);
                    RetX = RetHeadHeight(1:k-1);
                    
                    if length(AppY) < 2
                        AppY = AppForce(end-1:end);
                        AppX = AppHeadHeight(end-1:end);
                    end
                    if length(RetY) < 2
                        RetY = RetForce(1:2);
                        RetX = RetHeadHeight(1:2);
                    end
                    Fit1 = polyfit(AppX,AppY,1);
                    Fit2 = polyfit(RetX,RetY,1);
                    Slope1 = 1/(1/Fit1(1)-1/obj.RefSlope);
                    Slope2 = 1/(1/Fit2(1)-1/obj.RefSlope);
                    
                    PeakIndentationAngle(i) = rad2deg(atan(abs((Slope1 - Slope2)/(1+Slope1*Slope2))));
                catch
                    PeakIndentationAngle(i) = nan;
                end
                %                 % DebugSection
                %                 plot(AppTipHeight,AppForce,RetTipHeight,RetForce,...
                %                     AppX,AppY,RetX,RetY,...
                %                     AppTipHeight,polyval(Fit1,AppTipHeight),RetTipHeight,polyval(Fit2,RetTipHeight))
                %                 title(sprintf('%.4f OnFib: %i',PeakIndentationAngle(i),obj.FibMask(obj.List2Map(i,1),obj.List2Map(i,2))))
                %                 legend({'App','Ret','FitApp','FitRet','Fit1','Fit2'})
                %                 drawnow
            end
            
            obj.PeakIndentationAngle = PeakIndentationAngle;
            
            
            % Write to Channel
            PIAMap = obj.convert_data_list_to_map(obj.PeakIndentationAngle);
            Channel = obj.create_standard_channel(PIAMap,'Peak Indentation Angle','Degrees');
            obj.add_channel(Channel,~KeepOldResults)
        end
        
        function [DissipatedEnergy,ElasticEnergy] = calculate_dissipated_and_elastic_energy(obj,KeepOldResults)
            
            Range = find(obj.SelectedCurves);
            DissipatedEnergy = zeros(obj.NCurves,1);
            ElasticEnergy = zeros(obj.NCurves,1);
            for i=Range'
                [App,HHApp] = obj.get_force_curve_data(i,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                [Ret,HHRet] = obj.get_force_curve_data(i,'AppRetSwitch',1,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                AppForce = App - obj.CP(i,2);
                AppTipHeight = (HHApp - obj.CP(i,1)) - AppForce/obj.SpringConstant;
                RetForce = Ret - obj.CP(i,2);
                RetTipHeight = (HHRet - obj.CP(i,1)) - RetForce/obj.SpringConstant;
                try
                    AppX = AppTipHeight(AppTipHeight > 0);
                    AppY = AppForce(AppTipHeight > 0);
                    RetX = RetTipHeight(RetTipHeight > 0);
                    RetY = RetForce(RetTipHeight > 0);
                    
                    AppEnergy = trapz(AppX,AppY);
                    RetEnergy = trapz(RetX,RetY);
                    
                    DissipatedEnergy(i) = abs(AppEnergy - RetEnergy);
                    ElasticEnergy(i) = abs(RetEnergy);
                catch
                    DissipatedEnergy(i) = nan;
                    ElasticEnergy(i) = nan;
                end
            end
            
            CoulombConstant = 1.602176634e-19;
            
            obj.DissipatedEnergy = DissipatedEnergy;
            
            DEMap = obj.convert_data_list_to_map(obj.DissipatedEnergy);
            
            % Write to Channel
            Channel = obj.create_standard_channel(DEMap,'Dissipated Energy','J');
            obj.add_channel(Channel,~KeepOldResults)
            
            DEMap = DEMap./CoulombConstant;
            Channel = obj.create_standard_channel(DEMap,'Dissipated eV-Energy','eV');
            obj.add_channel(Channel,~KeepOldResults)
            
            obj.ElasticEnergy = ElasticEnergy;
            
            EEMap = obj.convert_data_list_to_map(obj.ElasticEnergy);
            Channel = obj.create_standard_channel(EEMap,'Elastic Energy','J');
            obj.add_channel(Channel,~KeepOldResults)
            
            EEMap = EEMap./CoulombConstant;
            Channel = obj.create_standard_channel(EEMap,'Elastic eV-Energy','eV');
            obj.add_channel(Channel,~KeepOldResults)
            
            Channel = obj.create_standard_channel(100.*EEMap./(DEMap+EEMap),'Elastic Fraction','%');
            obj.add_channel(Channel,~KeepOldResults)
            
            Channel = obj.create_standard_channel(100.*DEMap./(DEMap+EEMap),'Inelastic Fraction','%');
            obj.add_channel(Channel,~KeepOldResults)
        end
        
        function DZslopeCorrected = calculate_corrected_dzslopes(obj,varargin)
            % function DZslopeCorrected = calculate_corrected_dzslopes(obj,varargin)
            %
            % Calculates DZ-Slopes based on SensitivityCorrectionMethod's
            % corrected sensitivity
            %
            %
            % Required inputs
            % obj ... ForceMap instance
            %
            % Name-Value pairs
            % "FitRangeMode" ... <NAMEVALUE DESCRIPTION>
            % "FitRangeLowerFraction" ... <NAMEVALUE DESCRIPTION>
            % "FitRangeUpperFraction" ... <NAMEVALUE DESCRIPTION>
            % "FitRangeLowerValue" ... <NAMEVALUE DESCRIPTION>
            % "FitRangeUpperValue" ... <NAMEVALUE DESCRIPTION>
            % "AppRetSwitch" ... <NAMEVALUE DESCRIPTION>
            % "SensitivityCorrectionMethod" ... <NAMEVALUE DESCRIPTION>
            % "KeepOldResults"
            
            p = inputParser;
            p.FunctionName = "calculate_corrected_dzslopes";
            p.CaseSensitive = false;
            p.PartialMatching = true;
            
            % Required inputs
            validobj = @(x)true;
            addRequired(p,"obj",validobj);
            
            % NameValue inputs
            defaultAppRetSwitch = 0;
            defaultFitRangeMode = 'ValueHorizontal';
            defaultFitRangeLowerFraction = 0.75;
            defaultFitRangeUpperFraction = 1;
            defaultFitRangeLowerValue = 0;
            defaultFitRangeUpperValue = 4e-9;
            defaultKeepOldResults = false;
            defaultSensitivityCorrectionMethod = 'corrected';
            validAppRetSwitch = @(x)x==0|x==1|islogical(x);
            validFitRangeMode = @(x)any(validatestring(x,{'ValueHorizontal','FractionHorizontal','ValueVertical','FractionVertical'}));
            validSensitivityCorrectionMethod = @(x)any(validatestring(x,{'corrected','adaptive','original'}));
            validFitRangeLowerFraction = @(x)(0<=x)&&(x<=1);
            validFitRangeUpperFraction = @(x)(0<=x)&&(x<=1);
            validFitRangeLowerValue = @(x)isscalar(x);
            validFitRangeUpperValue = @(x)isscalar(x);
            validKeepOldResults = @(x)true;
            addParameter(p,"FitRangeMode",defaultFitRangeMode,validFitRangeMode);
            addParameter(p,"FitRangeLowerFraction",defaultFitRangeLowerFraction,validFitRangeLowerFraction);
            addParameter(p,"FitRangeUpperFraction",defaultFitRangeUpperFraction,validFitRangeUpperFraction);
            addParameter(p,"FitRangeLowerValue",defaultFitRangeLowerValue,validFitRangeLowerValue);
            addParameter(p,"FitRangeUpperValue",defaultFitRangeUpperValue,validFitRangeUpperValue);
            addParameter(p,"AppRetSwitch",defaultAppRetSwitch,validAppRetSwitch);
            addParameter(p,"SensitivityCorrectionMethod",defaultSensitivityCorrectionMethod,validSensitivityCorrectionMethod);
            addParameter(p,"KeepOldResults",defaultKeepOldResults,validKeepOldResults);
            parse(p,obj,varargin{:});
            
            % Assign parsing results to named variables
            obj = p.Results.obj;
            FitRangeMode = p.Results.FitRangeMode;
            FitRangeLowerFraction = p.Results.FitRangeLowerFraction;
            FitRangeUpperFraction = p.Results.FitRangeUpperFraction;
            FitRangeLowerValue = p.Results.FitRangeLowerValue;
            FitRangeUpperValue = p.Results.FitRangeUpperValue;
            AppRetSwitch = p.Results.AppRetSwitch;
            SensitivityCorrectionMethod = p.Results.SensitivityCorrectionMethod;
            KeepOldResults = p.Results.KeepOldResults;
            
            
            DZslopeCorrected = ones(obj.NCurves,1)*NaN;
            
            Range = find(obj.SelectedCurves);
            
            % first, calculate all the DZ-Slopes
            for i=Range'
                [vDef,Height] = obj.get_force_curve_data(i,'AppRetSwitch',AppRetSwitch,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity',SensitivityCorrectionMethod,'Unit','m');
                switch FitRangeMode
                    case 'ValueHorizontal'
                        UValue = max(Height) - FitRangeLowerValue;
                        LValue = max(Height) - FitRangeUpperValue;
                        FitHeight = Height(Height>=LValue & Height<=UValue);
                        FitvDef = vDef(Height>=LValue & Height<=UValue);
                    case 'ValueVertical'
                        FitHeight = Height(vDef>=FitRangeLowerValue & vDef<=FitRangeUpperValue);
                        FitvDef = vDef(vDef>=FitRangeLowerValue & vDef<=FitRangeUpperValue);
                    case 'FractionHorizontal'
                        RangeHeight = range(Height);
                        FitHeight = Height((Height-min(Height))>=FitRangeLowerFraction*RangeHeight &...
                            (Height-min(Height))<=FitRangeUpperFraction*RangeHeight);
                        FitvDef = vDef((Height-min(Height))>=FitRangeLowerFraction*RangeHeight &...
                            (Height-min(Height))<=FitRangeUpperFraction*RangeHeight);
                    case 'FractionVertical'
                        MaxvDef = max(vDef);
                        FitHeight = Height(vDef>=FitRangeLowerFraction*MaxvDef & vDef<=FitRangeUpperFraction*MaxvDef);
                        FitvDef = vDef(vDef>=FitRangeLowerFraction*MaxvDef & vDef<=FitRangeUpperFraction*MaxvDef);
                end
                Params(i,:) = polyfit(FitHeight,FitvDef,1);
                DZslopeCorrected(i) = Params(i,1);
            end
            
            % Convert to map and write to Channel
            Map = obj.convert_data_list_to_map(DZslopeCorrected);
            if isequal(lower(SensitivityCorrectionMethod),'corrected')
                NameModifier = 'global';
            else
                NameModifier = SensitivityCorrectionMethod;
            end
            Channel = obj.create_standard_channel(Map,...
                sprintf('Corrected %s DZ-Slope',NameModifier),...
                'm/m');
            obj.add_channel(Channel,~KeepOldResults);
            
            % Assign property
            obj.DZSlopeCorrected = DZslopeCorrected;            
        end
        
        function [NCMEnergy,PreCPNCMEnergy] = calculate_non_contact_model_energy(obj,varargin)
            % function [NCMEnergy,PreCPNCMEnergy] = calculate_non_contact_model_energy(obj,varargin)
            %
            % <FUNCTION DESCRIPTION HERE>
            %
            %
            % Required inputs
            % obj ... <VARIABLE DESCRIPTION>
            %
            % Name-Value pairs
            % "FMA_ID" ... <NAMEVALUE DESCRIPTION>
            % "SensitivityCorrectionMethod" ... <NAMEVALUE DESCRIPTION>
            % "KeepOldResults" ... <NAMEVALUE DESCRIPTION>
            
            p = inputParser;
            p.FunctionName = "calculate_non_contact_model_energy";
            p.CaseSensitive = false;
            p.PartialMatching = true;
            
            % Required inputs
            validobj = @(x)true;
            addRequired(p,"obj",validobj);
            
            % NameValue inputs
            defaultFMA_ID = obj.CurrentFMA_ID;
            defaultSensitivityCorrectionMethod = 'original';
            defaultKeepOldResults = true;
            validFMA_ID = @(x)true;
            validSensitivityCorrectionMethod = @(x)true;
            validKeepOldResults = @(x)true;
            addParameter(p,"FMA_ID",defaultFMA_ID,validFMA_ID);
            addParameter(p,"SensitivityCorrectionMethod",defaultSensitivityCorrectionMethod,validSensitivityCorrectionMethod);
            addParameter(p,"KeepOldResults",defaultKeepOldResults,validKeepOldResults);
            
            parse(p,obj,varargin{:});
            
            % Assign parsing results to named variables
            obj = p.Results.obj;
            FMA_ID = p.Results.FMA_ID;
            SensitivityCorrectionMethod = p.Results.SensitivityCorrectionMethod;
            KeepOldResults = p.Results.KeepOldResults;
            
            
            NCMEnergy = ones(obj.NCurves,1)*NaN;
            PreCPNCMEnergy = ones(obj.NCurves,1)*NaN;
            HertzEnergy = ones(obj.NCurves,1)*NaN;
            
            Range = find(obj.SelectedCurves);
            
            
            ID_Index = find(contains(...
                {obj.HertzFitStore.FMA_ID},...
                FMA_ID));
            if isempty(ID_Index)
                warning(sprintf('Found no Hertz fit data with ID: %s', FMA_ID))
                return
            end
            
            % first, calculate all the DZ-Slopes
            for i=Range'
                [Force,Height] = obj.get_force_curve_data(i,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',1,...
                    'Sensitivity',SensitivityCorrectionMethod);
                
                if iscell(obj.HertzFitStore(ID_Index).HertzFitType)
                    FitFunction = ...
                        fittype(obj.HertzFitStore(ID_Index).HertzFitType{i});
                else
                    FitFunction = ...
                        fittype(obj.HertzFitStore(ID_Index).HertzFitType);
                end
                FitCoeffValues = ...
                    obj.HertzFitStore(ID_Index).HertzFitValues{i};
                if length(FitCoeffValues) == 1
                    Fit = cfit(FitFunction,FitCoeffValues);
                elseif length(FitCoeffValues) == 2
                    Fit = cfit(FitFunction,FitCoeffValues(1),FitCoeffValues(2));
                end
                CP = obj.CP(i,:);
                Force = Force(Height - CP(1) >= 0);
                X3 = Height(Height - CP(1) >= 0) - CP(1);
                Y3 = feval(Fit,X3);
%                 Y3 = Y3 - CP(2);
                
                % Compute Energies
                FullEnergy = trapz(X3,Force);
                
                HertzEnergy(i) = trapz(X3(~isnan(Y3)), Y3(~isnan(Y3)));
                
                NCMEnergy(i) = FullEnergy - HertzEnergy(i);
                if length(FitCoeffValues) == 2
                    try
                        PreCPNCMEnergy(i) = trapz(X3(X3 >= 0 & X3 < -FitCoeffValues(2)),Force(X3 >= 0 & X3 < -FitCoeffValues(2)));
                    catch
                    end
                end
                
                clear FullEnergy
            end
            
            HertzEnergyMap = obj.convert_data_list_to_map(HertzEnergy);
            Channel = obj.create_standard_channel(HertzEnergyMap,'Contact Model Energy','J');
            obj.add_channel(Channel,~KeepOldResults)
            
            NCMEnergyMap = obj.convert_data_list_to_map(NCMEnergy);
            Channel = obj.create_standard_channel(NCMEnergyMap,'Non-Contact-Model Energy','J');
            obj.add_channel(Channel,~KeepOldResults)
            
            if length(FitCoeffValues) == 2
                PreCPNCMEnergyMap = obj.convert_data_list_to_map(PreCPNCMEnergy);
                Channel = obj.create_standard_channel(PreCPNCMEnergyMap,'Pre-Non-Contact-Model Energy','J');
                obj.add_channel(Channel,~KeepOldResults)
            end
            
        end
        
        function create_and_level_height_map_by_current_cp(obj,KeepOldResults)
            
            Map = obj.convert_data_list_to_map(-obj.CP(:,1));
            
            for i=1:obj.NumPixelsX
                Map(i,:) = AFMImage.replace_points_of_certain_value_in_line(Map(i,:),0);
            end
            
            for i=1:5
                Map = AFMImage.subtract_line_fit_vertical_rov(Map,.2,0);
            end
            
            Map = AFMImage.find_and_replace_outlier_lines(Map,10);
            
            % Write to Channel
            Channel = obj.create_standard_channel(Map,'Contact Height','m');
            obj.add_channel(Channel,~KeepOldResults)
        end
        
        function create_height_map_by_current_cp_and_level_by_other_channel(obj,OtherChannelName,KeepOldResults)
            
            Map = obj.convert_data_list_to_map(-obj.CP(:,1));
            
            for i=1:obj.NumPixelsX
                Map(i,:) = AFMImage.replace_points_of_certain_value_in_line(Map(i,:),0);
            end
            
            Target = obj.create_standard_channel(Map,'Raw Contact Height','m');
            
            Source = obj.get_channel(OtherChannelName);
            
            Map = obj.flatten_image_by_other_channels_fitparams(Source,Target);
            
            Map = AFMImage.find_and_replace_outlier_lines(Map,10);
            
            % Write to Channel
            Channel = obj.create_standard_channel(Map,'Contact Height','m');
            obj.add_channel(Channel,~KeepOldResults)
        end
        
        function manual_exclusion(obj)
            
            CheckSum = 100;
            f = figure('Name','Choose areas to be excluded');
            f.WindowState = 'maximized';
            
            if ~isempty(obj.ExclMask)
                obj.ExclMask = logical(ones(obj.NumPixelsX,obj.NumPixelsY));
            else
                TempMap = obj.HeightMap(:,:,1);
                TempMap(obj.ExclMask == 0) = max(TempMap,[],'all');
                imshow(TempMap,[])
                
                answer = questdlg('There already exists an exclusion mask!', ...
                    'Keep it!', ...
                    'Draw new!','Draw new!');
                % Handle response
                switch answer
                    case 'Keep it!'
                        close(f)
                        return
                    case 'Draw new!'
                        obj.ExclMask = logical(ones(obj.NumPixelsX,obj.NumPixelsY));
                        close(f)
                        f = figure('Name','Choose areas to be excluded');
                        f.WindowState = 'maximized';
                end
            end
            
            while CheckSum > 1
                subplot(2,1,2)
                surf(imresize(imrotate(obj.HeightMap(:,:,1)',90),[1024 1024]),'LineStyle','none','FaceLighting','gouraud','FaceColor','interp')
                light('Style','local')
                subplot(2,1,1)
                TempMap = obj.HeightMap(:,:,1);
                TempMap(obj.ExclMask == 0) = max(TempMap,[],'all');
                imshow(TempMap,[])
                title(sprintf('%s: Draw Freehand ROI around areas, that are to be excluded\nThe area will be taken out and the same map redrawn \n If there is nothing to do just click on the image once without dragging the cursor',obj.Name))
                ROI = drawfreehand;
                CheckSum = length(ROI.Waypoints);
                Mask = ~createMask(ROI);
                obj.ExclMask = obj.ExclMask.*Mask;
                drawnow
            end
            
            close(f)
            %             current = what();
            %             cd(obj.Folder)
            %             savename = sprintf('%s.mat',obj.Name);
            %             save(savename,'obj')
            %             cd(current.path)
        end
        
        function calculate_fib_diam(obj)
            
            Chan = obj.get_channel('Processed');
            if isempty(Chan)
                Chan = obj.get_channel('Height');
                if isempty(Chan)
                    HeightMap = obj.HeightMap(:,:,1);
                else
                    HeightMap = Chan.Image;
                end
            else
                HeightMap = Chan.Image;
            end
            
            [obj.Apex,obj.ApexIndex] = max(HeightMap.*obj.FibMask,[],2);
            
            if obj.FibrilFlag.Straight == 1
                obj.RectApex = zeros(obj.NumPixelsX,1);
                obj.RectApexIndex = zeros(obj.NumPixelsX,1);
                obj.RectApexIndex = round(predictGP_mean([1:obj.NumPixelsX]',[1:obj.NumPixelsX]',1,5*obj.NumPixelsX,obj.ApexIndex,1));
                for i=1:obj.NumPixelsX
                    obj.RectApex(i) = obj.HeightMap(i,obj.RectApexIndex(i),1);
                end
            else
                obj.RectApex = obj.Apex;
                obj.RectApexIndex = obj.ApexIndex;
            end
            
            k = 1;
            for i=1:obj.NumPixelsX
                if obj.ExclMask(i,obj.RectApexIndex(i)) == 1
                    FibHeight(k) = obj.RectApex(i);
                    k = k + 1;
                end
            end
            obj.FibDiam = mean(FibHeight);
            obj.FibDiamSTD = std(FibHeight,'omitnan');
            
            % Convert RectApexIndex and ApexIndex, so that they are
            % consistent with the Map2List-List2Map format (its entrances
            % being the corresponding list values)
            for i=1:obj.NumPixelsX
                obj.RectApexIndex(i) = obj.Map2List(i,obj.RectApexIndex(i));
                obj.ApexIndex(i) = obj.Map2List(i,obj.ApexIndex(i));
            end
            
            %             current = what();
            %             cd(obj.Folder)
            %             savename = sprintf('%s.mat',obj.Name);
            %             save(savename,'obj')
            %             cd(current.path)
        end
        
        function [img,imgorsize] = force2img(obj,ImgSize)
            % Generate a 2D image of the forcecurve and scale to a
            % ImgSizexImgSize grayscale image. Returns a
            % sum(obj.selcted_curves)x1 cell array of binary images
            if nargin<2
                ImgSize = 128;
            end
            Range = find(obj.SelectedCurves);
            imres = obj.Header{2,2};
            k = 1;
            img = cell(sum(obj.SelectedCurves),1);
            imgorsize = cell(sum(obj.SelectedCurves),1);
            h = waitbar(0,'Setting up...');
            for i=Range'
                [App,~] = obj.get_force_curve_data(i,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                prog = i/obj.NCurves;
                waitbar(prog,h,'Converting force curves to images...');
                norm =  round(normalize(App,'range',[1 imres]));
                mat = zeros(imres,imres);
                for m=1:length(norm)
                    mat(m,norm(m)) = 1;
                end
                imgorsize{k} = imrotate(mat,90);
                img{k} = imresize(imgorsize{k},[ImgSize ImgSize],'bilinear');
                k = k + 1;
            end
            close(h)
        end
        
        function choose_fibril(obj,pth_quantile)
            % sets every forcecurve below a certain threshold to zero in
            % obj.SelectedCurves. Note that curves on the fibril already
            % chosen to be disregarded will keep that status!.This function
            % is a crude way to find fibril areas and should not be used, if
            % its important to select exclusively indentations on the fibril
            % or ,in the more general case,other elevated region of interest
            if nargin < 2
                pth_quantile = 0.8;
            end
            q = quantile(obj.HeightMap(:,:,1),pth_quantile,'All');
            mask = ones(size(obj.HeightMap));
            for i=1:obj.NumPixelsX
                for j=1:obj.NumPixelsY
                    if obj.HeightMap(i,j,1) < q
                        obj.SelectedCurves(obj.HeightMap(i,j,2)) = 0;
                        mask(i,j,1) = 0;
                    end
                end
            end
            f = figure();
            imshowpair(imresize(mat2gray(obj.HeightMap(:,:,1)),[1024 1024]),imresize(mask(:,:,1),[1024 1024]),'montage')
            %            pause(5)
            close(f)
        end
        
        function determine_relevant_radius_of_indentation(obj,varargin)
            % function determine_relevant_radius_of_indentation(obj,varargin)
            %
            % <FUNCTION DESCRIPTION HERE>
            %
            %
            % Required inputs
            % obj ... <VARIABLE DESCRIPTION>
            %
            % Name-Value pairs
            % "Method" ... <NAMEVALUE DESCRIPTION>
            % "CantileverTipInstance" ... <NAMEVALUE DESCRIPTION>
            % "CantileverTipChannelName" ... <NAMEVALUE DESCRIPTION>
            % "ChannelName" ... <NAMEVALUE DESCRIPTION>
            % "IndentDepthChannelName" ... <NAMEVALUE DESCRIPTION>
            % "KeepOldResults" ... <NAMEVALUE DESCRIPTION>
            % "Verbose" ... <NAMEVALUE DESCRIPTION>
            % "KeyFrame" ... <NAMEVALUE DESCRIPTION>
            
            p = inputParser;
            p.FunctionName = "determine_relevant_radius_of_indentation";
            p.CaseSensitive = false;
            p.PartialMatching = true;
            
            % Required inputs
            validobj = @(x)true;
            addRequired(p,"obj",validobj);
            
            % NameValue inputs
            defaultMethod = [];
            defaultCantileverTipInstance = [];
            defaultCantileverTipChannelName = 'Eroded Tip';
            defaultChannelName = 'Contact Height';
            defaultIndentDepthChannelName = 'Indentation Depth';
            defaultKeepOldResults = true;
            defaultVerbose = false;
            defaultKeyFrame = 23;
            validMethod = @(x)true;
            validCantileverTipInstance = @(x)true;
            validCantileverTipChannelName = @(x)true;
            validChannelName = @(x)true;
            validIndentDepthChannelName = @(x)true;
            validKeepOldResults = @(x)true;
            validVerbose = @(x)true;
            validKeyFrame = @(x)true;
            addParameter(p,"Method",defaultMethod,validMethod);
            addParameter(p,"CantileverTipInstance",defaultCantileverTipInstance,validCantileverTipInstance);
            addParameter(p,"CantileverTipChannelName",defaultCantileverTipChannelName,validCantileverTipChannelName);
            addParameter(p,"ChannelName",defaultChannelName,validChannelName);
            addParameter(p,"IndentDepthChannelName",defaultIndentDepthChannelName,validIndentDepthChannelName);
            addParameter(p,"KeepOldResults",defaultKeepOldResults,validKeepOldResults);
            addParameter(p,"Verbose",defaultVerbose,validVerbose);
            addParameter(p,"KeyFrame",defaultKeyFrame,validKeyFrame);
            
            parse(p,obj,varargin{:});
            
            % Assign parsing results to named variables
            obj = p.Results.obj;
            Method = p.Results.Method;
            CantileverTipInstance = p.Results.CantileverTipInstance;
            CantileverTipChannelName = p.Results.CantileverTipChannelName;
            ChannelName = p.Results.ChannelName;
            IndentDepthChannelName = p.Results.IndentDepthChannelName;
            KeepOldResults = p.Results.KeepOldResults;
            Verbose = p.Results.Verbose;
            KeyFrame = p.Results.KeyFrame;
           
            
            TopoChannel = obj.get_channel(ChannelName);
            IndentDepthChannel = obj.get_channel(IndentDepthChannelName);
            CTChannel = CantileverTipInstance.get_channel(CantileverTipChannelName);
            
            if isempty(TopoChannel)
                error(sprintf('Channel not found: %s',ChannelName));
            end
            if isempty(IndentDepthChannel)
                error(sprintf('Channel not found: %s',IndentDepthChannelName));
            end
            if isempty(CTChannel)
                error(sprintf('Channel not found: %s', CantileverTipChannelName));
            end
            
            % Resize CTChannel to obj resolution
            CTChannel = obj.resize_channel(CTChannel,[TopoChannel.NumPixelsX TopoChannel.NumPixelsY]);
            
            
            
            IndentDepthList = obj.convert_map_to_data_list(IndentDepthChannel.Image);
            TopoList = obj.convert_map_to_data_list(TopoChannel.Image);
            
            PixelArea = TopoChannel.ScanSizeX/TopoChannel.NumPixelsY * TopoChannel.ScanSizeY/TopoChannel.NumPixelsX;
            ContactRadius = zeros(size(TopoList));
            MaxCTHeight = max(CTChannel.Image,[],'all');
            [MaxCTIndexX,MaxCTIndexY] = find(CTChannel.Image == MaxCTHeight);
            
            h = waitbar(0,'Setting up...');
            for i=1:obj.NCurves
                waitbar(i/obj.NCurves,h,'Calculating contact radius...');
                TempCTMask = CTChannel.Image >= (MaxCTHeight - IndentDepthList(i));
                CurPos = obj.List2Map(i,:);
                PixShiftX = CurPos(1) - MaxCTIndexX ;
                PixShiftY = CurPos(2) - MaxCTIndexY;
                CTMask = shiftLogicalArray(TempCTMask,PixShiftX,PixShiftY);
                TopoMask = TopoChannel.Image >= (TopoList(i) - IndentDepthList(i));
                AreaMask = CTMask & TopoMask;
                SumPixels = sum(AreaMask,'all');
                ContactRadius(i) = sqrt((SumPixels*PixelArea)/pi);
                if Verbose && ~mod(i,KeyFrame)
                    subplot(2,1,1)
                    CurPosMap = zeros(size(TempCTMask));
                    CurPosMap(CurPos(1),CurPos(2)) = true;
                    imshowpair(CTChannel.Image >= (MaxCTHeight - IndentDepthList(i)),TopoChannel.Image >= (TopoList(i) - IndentDepthList(i)),'montage');
                    subplot(2,1,2)
                    imshowpair(CurPosMap,AreaMask,'montage')
                    title(sprintf('Contact Radius: %.2e',ContactRadius(i)))
                    drawnow
                end
            end
            
            CRMap = obj.convert_data_list_to_map(ContactRadius);
            
            CRChannel = obj.create_standard_channel(CRMap,'Contact Radius','m');
            obj.add_channel(CRChannel,~KeepOldResults);
            
            close(h)
        end

        
    %% SMFS related
    
        function fc_flag_workaround(obj)
    % fc_flag_workaround: Function to fix some differences in the
    % prerequirements between the SMFS batch and imaging parts of the class
      
             for ii=1:obj.NCurves
             obj.BaseAndTiltFlag=1;
         %    obj.BigDataFlag=1;
             end
        end
        
        function fc_estimate_cp_hardsurface(obj)
            % contact point estimation for force curves detected on hard
            % surfaces
            
            for ii=1:obj.NCurves 
                obj.CP_HardSurface(ii,1) = obj.HHApp{ii}(end) - obj.App{ii}(end)/obj.SpringConstant; % Determine the contact point by simply substracting the last entries of height and force and correct with the spring constant
                obj.CP_HardSurface(ii,2) = 0;
                %% Debugging
                % plot(HHApp,App);
                % drawpoint('Position',[obj.CP_HardSurface(i,1) obj.CP_HardSurface(i,2)]);
            end
            obj.CPFlag.HardSurface = 1;
        end
        
        
        function fc_selection_threshold(obj,ThresholdDist,ThreshValue)
            % fc_selection_threshold: A function to distinguish between
            % force curves that fulfil or do not fulfil the logical
            % statement
            % SMFSFlag.Min = 0:  Force curves indicate a naked tip
            % SMFSFlag.Min = 1:  Force curves indicate a functionalized tip
            if nargin <2
                ThresholdDist1=50e-9;   % 50 nm
                ThresholdDist2=75e-9;   % 75 nm
                ThresholdDist3=100e-9;  % 100 nm
                ThresholdDist4=150e-9;  % 150 nm
                ThreshValue=25e-12;     % 25 pN
            elseif nargin<3
                ThreshValue=25e-12;    % 25 pN
            end
            % loop over all force curves
            for kk=1:100
            % Determine the index corresponding to the threshold distance
            ThreshDist1=abs(obj.HHRet{kk}-obj.CP_HardSurface(kk,1)+ThresholdDist1);
            [~, ThreshIdx1]=min(ThreshDist1);
            ThreshDist2=abs(obj.HHRet{kk}-obj.CP_HardSurface(kk,1)+ThresholdDist2);
            [~, ThreshIdx2]=min(ThreshDist2);
            ThreshDist3=abs(obj.HHRet{kk}-obj.CP_HardSurface(kk,1)+ThresholdDist3);
            [~, ThreshIdx3]=min(ThreshDist3);
            ThreshDist4=abs(obj.HHRet{kk}-obj.CP_HardSurface(kk,1)+ThresholdDist4);
            [~, ThreshIdx4]=min(ThreshDist4);
            % Check if the force curve is selected 
               if obj.BasedRet{kk}(ThreshIdx1)>ThreshValue
                   obj.SMFSFlag.Min(kk)=1;
               elseif obj.BasedRet{kk}(ThreshIdx1)<ThreshValue && obj.BasedRet{kk}(ThreshIdx2)>ThreshValue
                   obj.SMFSFlag.Min(kk)=1;
               elseif obj.BasedRet{kk}(ThreshIdx1)<ThreshValue && obj.BasedRet{kk}(ThreshIdx2)<ThreshValue && obj.BasedRet{kk}(ThreshIdx3)>ThreshValue
                   obj.SMFSFlag.Min(kk)=1;
               elseif obj.BasedRet{kk}(ThreshIdx1)<ThreshValue && obj.BasedRet{kk}(ThreshIdx2)<ThreshValue && obj.BasedRet{kk}(ThreshIdx3)<ThreshValue && obj.BasedRet{kk}(ThreshIdx4)>ThreshValue
                   obj.SMFSFlag.Min(kk)=1;
               else
                   obj.SMFSFlag.Min(kk)=0;
               end
            end           
%             %% Appendix
%             close all
%             % Define variables
%             kk=1
%             x100=-100e-9; % Defines 100nm
%             x500=-500e-9; % Defines 500nm
%             % Graphical preview
%             fig=gcf;
%             fig.Units='normalized'; % changes to normalized unit settings, necessary to receive the full screen size in the next line
%             fig.Color='white'; % changes the background color of the figure
%             fig.OuterPosition=[0.5 0 0.5 1];% changes the size of the figure to half screen
%             fig.PaperOrientation='landscape';
%             grid on
%             hold on
%             plot(obj.FM{1}.THApp{kk}-obj.FM{1}.CP_HardSurface(kk,1),obj.FM{1}.BasedApp{kk},'b');
%             plot(obj.FM{1}.THRet{kk}-obj.FM{1}.CP_HardSurface(kk,1),obj.FM{1}.BasedRet{kk},'r'); 
%             plot(obj.FM{1}.THRet{kk}(ThreshIdx,1)-obj.FM{1}.CP_HardSurface(kk,1),obj.FM{1}.BasedRet{kk}(ThreshIdx,1),'kx','MarkerSize',20);
%             line([x100 x100], ylim,'Color','k'); % Draws a vertical line                  
%             line([x500 x500], ylim,'Color','k'); % Draws a vertical line      
        end
                
        function fc_measurement_prop(obj)
               % fc_chipprop: A fct to read out properties about the SMFS measurements from the name of the jpk data file (i.e. "JPK-FORCE-MAP").  
                % Chip number and Cantilever
                exp15='(\d+\D{1}\>)'; % Finds the chip number and the cantilever   
                obj.ChipCant = regexp(obj.Name, exp15, 'match','once');                 
                % Chip box
                exp16 = '(?!L)[CLXVI]+'; % Finds the chip number given in roman numerals 
                obj.Chipbox = regexp(obj.Name, exp16, 'match','once');                 
                % Glass
                exp21='glass';
                pat=regexpPattern(exp21,"IgnoreCase",true);
                ext21=extract(obj.Name,pat);
                % Mica
                exp22='mica';
                pat=regexpPattern(exp22,"IgnoreCase",true);
                ext22=extract(obj.Name,pat);
                % inorganic bone (Hydroxyapatite)
                exp23='inorganic';
                pat=regexpPattern(exp23,"IgnoreCase",true);
                ext23=extract(obj.Name,pat);
                % organic bone 
                exp24='organic';
                pat=regexpPattern(exp24,"IgnoreCase",true);
                ext24=extract(obj.Name,pat);
                % poly Lysine 
                exp25='Lysine';
                pat=regexpPattern(exp25,"IgnoreCase",true);
                ext25=extract(obj.Name,pat);
                % Substrate
                if isempty(ext21)==0
                    obj.Substrate='glass';
                elseif isempty(ext22)==0
                    obj.Substrate='mica';
                elseif isempty(ext23)==0
                    obj.Substrate='hydAp'; % Hydroxyapatite
                elseif isempty(ext24)==0
                    obj.Substrate='orgBone'; % Organic Bone
                elseif isempty(ext25)==0
                    obj.Substrate='lysine'; % poly-Lysine
                else
                    obj.Substrate='glass';
                end
                % Milli-Q water
                exp31='mil';
                pat=regexpPattern(exp31,"IgnoreCase",true);
                ext31=extract(obj.Name,pat);
                % Acetic acid (HAc)
                exp32='HAc';
                pat=regexpPattern(exp32,"IgnoreCase",true);
                ext32=extract(obj.Name,pat);
                % CAPS buffer
                exp33='CAPS';
                pat=regexpPattern(exp33,"IgnoreCase",true);
                ext33=extract(obj.Name,pat);
                % Dulbecco's phospate buffered saline (DPBS) buffer
                exp34='DPBS';
                pat=regexpPattern(exp34,"IgnoreCase",true);
                ext34=extract(obj.Name,pat);
                % Environmental conditions
                if isempty(ext31)==0
                    obj.EnvCond='Water'; % Milli-Q water
                elseif isempty(ext32)==0
                    obj.EnvCond='HAc'; % Acetic acid
                elseif isempty(ext33)==0
                    obj.EnvCond='CAPS'; % CAPS
                elseif isempty(ext34)==0
                    obj.EnvCond='DPBS'; % Dulbecco's phospate buffered saline
                else
                    obj.EnvCond='';
                end
                % Long linker
                exp41='long';
                pat=regexpPattern(exp41,"IgnoreCase",true);
                ext41=extract(obj.Name,pat);
                % Short linker
                exp42='short';
                pat=regexpPattern(exp42,"IgnoreCase",true);
                ext42=extract(obj.Name,pat);              
                % Linker
                if isempty(ext41)==0
                    obj.Linker='long'; % long linker
                elseif isempty(ext42)==0
                    obj.Linker='short'; % short linker
                else
                  %  obj.Linker='long';
                    obj.Linker='short';
                end
                
        end
        
        function [fitresult, gof] = fc_sinoidal_fit(obj)
        %CREATEFIT1(X,Y)
        %  Create a fit.
        %
        %  Data for 'untitled fit 1' fit:
        %      X Input : x
        %      Y Output: y
        %  Output:
        %      fitresult : a fit object representing the fit.
        %      gof : structure with goodness-of fit info.
        %
        %  See also FIT, CFIT, SFIT.

        %  Auto-generated by MATLAB on 02-Sep-2021 17:40:40

            % Define variables
            RGB3=[80 220 100]./255; % Emerald
            RGB4=[200 81 160]./255; % Compl to RGB3
            mega=10e6;
            giga=10e9;
            SelectedPercentage=0.7;
            NFigures=4;
            NLoop=25;
            ExtendVelocityConvert=num2str(obj.ExtendVelocity*1e9);
            RetractVelocityConvert=num2str(obj.RetractVelocity*1e9);
            % Classification criteria
            figname=strcat(obj.Date,{'_'},obj.Time,{'_'},obj.ID,{'_'},obj.Substrate,{'_'},obj.EnvCond,{'_'},obj.Linker,{'_'},obj.Chipbox,{'_'},obj.ChipCant,{'_'},ExtendVelocityConvert,{'_'},RetractVelocityConvert);
            figname=char(figname);    
                % Figure loop
                for jj=1:NFigures
                % figure properties
                h_fig=figure(jj);
                h_fig.Color='white'; % changes the background color of the figure
                h_fig.Units='normalized'; % Defines the units 
                h_fig.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
                h_fig.PaperOrientation='landscape';
                h_fig.Name=figname; 
                % Tile properties
                t = tiledlayout(5,5);
                t.TileSpacing = 'none'; % To reduce the spacing between the tiles
                t.Padding = 'none'; % To reduce the padding of perimeter of a tile
                %% Tile loop
                    for kk=1:NLoop
                    % Tile
                    oo=kk+25*(jj-1);    
                    % Allocate data
                    x=obj.HHApp{oo}*mega; % Transform head height to nm, needed for the fit
                    y=obj.App{oo}*giga; % Transform approach force data to pN, needed for the fit
                    % Select only a percentage of the data so that a potential
                    % snap-in is excluded from the fit
                    x(ceil(end*SelectedPercentage):end)=[];
                    y(ceil(end*SelectedPercentage):end)=[];                    
                    % Prepare data
                    [xData, yData] = prepareCurveData( x, y );                  
                    % Set up fittype and options.
                    % Fit 1
                    ft1 = fittype( 'a1*sin(b1*x+c1)+d1', 'independent', 'x', 'dependent', 'y' );
                    opts1= fitoptions( 'Method', 'NonlinearLeastSquares' );
                    opts1.Display = 'Off';
                    opts1.Robust = 'Bisquare';
                    a1Range=range(yData)/2;
                    %opts1.StartPoint = [a1Range -1.8 x(1,1) y(1,1)];
                    opts1.StartPoint = [-0.2 -1.8 x(1,1) y(1,1)];      
                    [fitresult1, gof1] = fit( xData, yData, ft1, opts1 ); % Fit                    
                    % Fit 2
                    % Set up fittype and options.
                    ft2 = fittype( 'a1*sin(b1*x+c1)+d1+e1*x', 'independent', 'x', 'dependent', 'y' );
                    opts2= fitoptions( 'Method', 'NonlinearLeastSquares' );
                    opts2.Display = 'Off';                   
                    opts2.StartPoint = [fitresult1.a1 fitresult1.b1 fitresult1.c1 fitresult1.d1 0];                  
                    [fitresult2, gof2] = fit( xData, yData, ft2, opts2 ); % Fit                    
                    % Allocate variables
                    obj.FitCoeffa1(oo)=fitresult2.a1;
                    obj.FitCoeffb1(oo)=fitresult2.b1;
                    obj.FitCoeffc1(oo)=fitresult2.c1;
                    obj.FitCoeffd1(oo)=fitresult2.d1;
                    obj.FitCoeffe1(oo)=fitresult2.e1;
                    obj.FitRSquare(oo)=gof2.rsquare;
                    obj.FitSSE(oo)=gof2.sse;                     
                    % Flag fitting results basaed on the received R-squared
                    % value
                    RSquareLimit=0.5;
                    if obj.FitRSquare(oo)<RSquareLimit
                         obj.SMFSFlag.Fit(oo)=0;
                    else
                         obj.SMFSFlag.Fit(oo)=1;
                    end
                    
                    %% Plotting the tiles
                    nexttile;     
                    hold on
                    grid on
                    plot(xData,yData,'.','Color',RGB4) 
                    plot(xData,obj.FitCoeffa1(oo)*sin(obj.FitCoeffb1(oo)*xData+obj.FitCoeffc1(oo))+obj.FitCoeffd1(oo)+obj.FitCoeffe1(oo)*xData,'LineWidth',2,'Color',RGB3)
                    % Plotted Text  
                    NE = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.01; % Define the position in the plot  
                    partstrA1='R^2=';
                    partstrA2=num2str(round(gof2.rsquare,2));
                    fullstrA=strcat(partstrA1,partstrA2); % Define the string that shall be shown in the plot
                    te1=text(NE(1), NE(2),fullstrA, 'VerticalAlignment','top', 'HorizontalAlignment','right');
                    te1.FontSize = 14;
                    % Title for each Subplot                    
                    ti=title(sprintf('%i',oo),'Color','k');                    
                    ti.Units='normalized'; % Set units to 'normalized'  
                    ti.Position=[0.5,0.95]; % Position the subplot title within the subplot 
                    end                              
                 %% Save figures
                 %%% Define the name for the figure title   
                 partname=sprintf('-part%d',jj);  
                 fullname=sprintf('%s%s',figname,partname);
                 %%% Save the current figure in the current folder
                 print(gcf,fullname,'-dpng');             
                end
                close Figure 1 Figure 2 Figure 3 Figure 4
        end
                
        function fc_fit_based_yData(obj)
            % Based on the fitting equation this function calculates the
            % y-data of the corresponding retraction x-data
            % f(x) = a1*sin(b1*x+c1)+d1
            
            for jj=1:obj.NCurves
            %for jj=30
                % Define variables
                mega=10e6;
                giga=10e9;
                % Allocate data
                xApp=obj.HHApp{jj}*mega;
                xRet=obj.HHRet{jj}*mega; % Transform head height to nm
                yApp=obj.App{jj}*giga; % Transform approach force data to pN
                yRet=obj.Ret{jj}*giga; % Transform approach force data to pN                
                % Determine the corresponding y-data of the x retraction
                % data by using the determined fitting coefficiantsb (a1,b1,c1,d1) and the
                % fitting equation: f(x) = a1*sin(b1*x+c1)+d1+e1*x
                yAppFitbased=obj.FitCoeffa1(jj)*sin(obj.FitCoeffb1(jj)*xApp+obj.FitCoeffc1(jj))+obj.FitCoeffd1(jj)+obj.FitCoeffe1(jj)*xApp;
                yRetFitbased=obj.FitCoeffa1(jj)*sin(obj.FitCoeffb1(jj)*xRet+obj.FitCoeffc1(jj))+obj.FitCoeffd1(jj)+obj.FitCoeffe1(jj)*xRet;                
                % Substract the determined y-data from the original data
                % and allocate
                obj.BasedApp{jj}=(yApp-yAppFitbased)/giga;
                obj.BasedRet{jj}=(yRet-yRetFitbased)/giga;                
            end            
            %% Appendix
%             close all
%             % Plot  
%             hold on
%             plot(xApp,yApp,'b')
%             plot(xRet,yRet,'r')
%             plot(xRet,yRetFitbased,'g')
%             plot(xApp,obj.BasedApp(jj),'k')
%             plot(xRet,obj.BasedRet(jj),'m')
           
        end
                
        function fc_TipHeight_calculation(obj)
            % Function to determine the TipHeight based on the Head Height
            % and the fit based y data
            
            for jj=1:obj.NCurves
                obj.THApp{jj} = obj.HHApp{jj} - obj.BasedApp{jj}/obj.SpringConstant;
                obj.THRet{jj} = obj.HHRet{jj} - obj.BasedRet{jj}/obj.SpringConstant;
            end
        end

        
        function fc_visual_selection(obj,XMin,XMax,YMin,YMax) % fc ... force curve
            % fc_visual_selection: function plots all force curves of a force map
                        
            if nargin < 2
                XMin= -inf;
                XMax= inf;
                YMin= -inf;
                YMax= inf;
            end
            % Define remainder situation
            Remainder=mod(obj.NCurves,25);
            NFigures=floor(obj.NCurves./25);
            if Remainder ~= 0
                NFigures=NFigures+1;
            end    
            %% Define some variables
            RGB1=[0 26 255]./255;  % Blue 
            RGB2=[255 119 0]./255; % Orange
            RGB10=[69 22 113]./255; % Violet
            % Define variables for the plotted tiles 
                    x50=-50e-9; % Defines 50nm
                    x150=-150e-9; % Defines 150nm
                    x500=-500e-9; % Defines 500nm
            % Define variables for the figure name
            ExtendVelocityConvert=num2str(obj.ExtendVelocity*1e9);
            RetractVelocityConvert=num2str(obj.RetractVelocity*1e9);
            % Classification criteria
            figname=strcat(obj.Date,{'_'},obj.Time,{'_'},obj.ID,{'_'},obj.Substrate,{'_'},obj.EnvCond,{'_'},obj.Linker,{'_'},obj.Chipbox,{'_'},obj.ChipCant,{'_'},ExtendVelocityConvert,{'_'},RetractVelocityConvert);
            figname=char(figname);
            %% Figure loop   
            for ii=1:NFigures           
            % Figure    
            h_fig=figure(ii);
            h_fig.Color='white'; % changes the background color of the figure
            h_fig.Units='normalized'; % Defines the units 
            h_fig.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
            h_fig.PaperOrientation='landscape';
            h_fig.Name=figname;         
            % Plotting the tiles
            t = tiledlayout(5,5);
            %t.TileSpacing = 'compact';
            %t.Padding = 'compact';
            t.TileSpacing = 'none'; % To reduce the spacing between the tiles
            t.Padding = 'none'; % To reduce the padding of perimeter of a tile           
            % Defining variables
            if ii==NFigures && Remainder~=0
                NLoop=Remainder;
            else
                NLoop=25;
            end
                %% Plot loop    
                for jj=1:NLoop
                    % Tile jj
                    kk=jj+25*(ii-1);                    
                    % Plot tile
                    ax=nexttile;      
                    ax.XLim = [XMin XMax];
                    ax.YLim = [YMin YMax];
                    hold on
                    grid on
                    % For sinoidal fitted data
                    plot(obj.HHApp{kk}-obj.CP_HardSurface(kk,1),obj.BasedApp{kk},'Color',RGB1);
                    plot(obj.HHRet{kk}-obj.CP_HardSurface(kk,1),obj.BasedRet{kk},'Color',RGB2);
                    % Non fitted data
                    %plot(obj.THApp{kk}-obj.CP_HardSurface(kk,1),obj.BasedApp{kk});
                    %plot(obj.THRet{kk}-obj.CP_HardSurface(kk,1),obj.BasedRet{kk});
                    xline(x50,'Color',RGB10); % Draws a vertical line
                    xline(x150,'Color',RGB10); % Draws a vertical line   
                    xline(x500,'Color',RGB10); % Draws a vertical line
                    % Title for each Subplot
                    if obj.SMFSFlag.Uncorrupt(kk)==0
                        ti=title(sprintf('%i',kk),'Color','r');
                    elseif obj.SMFSFlag.Uncorrupt(kk)==1
                        ti=title(sprintf('%i',kk),'Color','b');
                    end
                    ti.Units='normalized'; % Set units to 'normalized'  
                    ti.Position=[0.5,0.95]; % Position the subplot title within the subplot                 
                end
                
                %% Dialog boxes
                % Function 'bttnChoiseDialog.m' is needed to excute this section               
                inputOptions={'Select all', 'Select none', 'Select all - except of', 'Select none - except of'}; % Define the input arguments
                % 'Select all' = 1
                % 'Select none' = 2
                % 'Select all - except of' = 3
                % 'Select none - except of' = 4
                defSelection=inputOptions{1}; % Default selection; Defined selection if the window is closed without choosing a selection possibility
                SelectBttns=bttnChoiseDialog(inputOptions, 'Force curve selection', defSelection,...
                'Please choose the appropriate button ...'); % Stores the selected button number per figure
                % Case 1: Select all
                if SelectBttns == 1
                    obj.SMFSFlag.Uncorrupt(kk-24:kk)=1;
                end
                % Case 2: Select none
                if SelectBttns == 2
                    obj.SMFSFlag.Uncorrupt(kk-24:kk)=0;
                end
                % Case 3: Select all - except of
                if SelectBttns == 3
                    obj.SMFSFlag.Uncorrupt(kk-24:kk)=1;
                    prompt = {'Enter the force curve number you do not want to keep for analysis (For multiple selections just use the space key to separeat entries)'};
                    definput = {''};
                    opts.Interpreter = 'tex';
                    IdxExc=inputdlg(prompt,'Select all - except of ...',[1 150],definput,opts); % Stores the individual selected fc as a cell array of character vectors 
                    IdxExc=str2num(IdxExc{1}); % Convert the cell array to numerals
                    obj.SMFSFlag.Uncorrupt(IdxExc)=0;
                end
                if obj.SMFSFlag.Uncorrupt(kk-24:kk)==0
                    title(sprintf('%i',kk),'Color','r');
                elseif obj.SMFSFlag.Uncorrupt(kk-24:kk)==1    
                    title(sprintf('%i',kk),'Color','b');
                end        
                % Case 4: Select none - except of
                if SelectBttns == 4
                    obj.SMFSFlag.Uncorrupt(kk-24:kk)=0;
                    prompt = {'Enter the force curve number you want want to keep for analysis (For multiple selections just use the space key to separeat entries)'};
                    definput = {''};
                    opts.Interpreter = 'tex';
                    IdxExc=inputdlg(prompt,'Select none - except of ...',[1 150],definput,opts); % Stores the individual selected fc as a cell array of character vectors 
                    IdxExc=str2num(IdxExc{1}); % Convert the cell array to numerals
                    obj.SMFSFlag.Uncorrupt(IdxExc)=1;
                end
            end
            close all
                
            %%% Colour highlighting of the force curves regarding the choosen answer and storage in a structure
            %% Figure loop
            figname=strcat(obj.ID,{'-'},obj.Name);
            figname=char(figname);
            for ii=1:NFigures  
            h_fig=figure(ii);
            h_fig.Color='white'; % changes the background color of the figure
            h_fig.Units='normalized'; % Defines the units 
            h_fig.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
            h_fig.PaperOrientation='landscape';
            h_fig.Name=figname; 
                %% Plotting the tiles
                t = tiledlayout(5,5);
                %t.TileSpacing = 'compact';
                %t.Padding = 'compact';
                t.TileSpacing = 'none'; % To reduce the spacing between the tiles
                t.Padding = 'none'; % To reduce the padding of perimeter of a tile
                % Defining variables
                if ii==NFigures && Remainder~=0
                    NLoop=Remainder;
                else
                    NLoop=25;
                end               
                %% Title loop
                for jj=1:NLoop
                    % Tile jj
                    kk=jj+25*(ii-1);
                    ax=nexttile;      
                    ax.XLim = [XMin XMax];
                    ax.YLim = [YMin YMax];
                    hold on
                    grid on
                    plot(obj.HHApp{kk}-obj.CP_HardSurface(kk,1),obj.BasedApp{kk},'Color',RGB1);
                    plot(obj.HHRet{kk}-obj.CP_HardSurface(kk,1),obj.BasedRet{kk},'Color',RGB2);
                    %plot(obj.THApp{kk}-obj.CP_HardSurface(kk,1),obj.BasedApp{kk});
                    %plot(obj.THRet{kk}-obj.CP_HardSurface(kk,1),obj.BasedRet{kk});
                    xline(x50,'Color',RGB10); % Draws a vertical line
                    xline(x150,'Color',RGB10); % Draws a vertical line   
                    xline(x500,'Color',RGB10); % Draws a vertical line
                    % Title for each Subplot
                    if obj.SMFSFlag.Uncorrupt(kk)==0
                        ti=title(sprintf('%i',kk),'Color','r');
                    elseif obj.SMFSFlag.Uncorrupt(kk)==1
                        ti=title(sprintf('%i',kk),'Color','b');
                    end
                    ti.Units='normalized'; % Set units to 'normalized'  
                    ti.Position=[0.5,0.95]; % Position the subplot title within the subplot                                  
                end

            %% Save figures
            %%% Define the name for the figure title    
            partname=sprintf('-p%d',ii);        
            fullname=sprintf('%s%s',figname,partname);
            %%% Save the current figure in the current folder
            print(gcf,fullname,'-dpng'); 
            end
        close Figure 1 Figure 2 Figure 3 Figure 4
        end
        
        
        function fc_based_ret_correction(obj,DataShareStartApp,DataShareEndApp,DataShareStartRet,DataShareEndRet)
            % fc_based_ret_correction: A function to correct for an AFM
            % based baseline deviation between the approach and retraction
            % data
        if nargin <2
            DataShareStartApp=0.05; % 5%
            DataShareEndApp=0.15; % 15%
            DataShareStartRet=0.01; % 1%
            DataShareEndRet=0.11; % 5%
        end
        % Loop over all force curves  
        for ii=1:obj.NCurves
            
            % Define limits
            DataPtsApp=size(obj.BasedApp{ii}); % Determine the amount of data points in the force curve 
            LimitIdxApp1=round(DataPtsApp(1)*DataShareStartApp); % Determine the corresponidng index
            LimitIdxApp2=round(DataPtsApp(1)*DataShareEndApp); % Determine the corresponidng index
            DataPtsRet=size(obj.BasedRet{ii}); % Determine the amount of data points in the force curve 
            LimitIdxRet1=round(DataPtsRet(1)*DataShareStartRet); % Determine the corresponidng index
            LimitIdxRet2=round(DataPtsRet(1)*DataShareEndRet); % Determine the corresponidng index
            
            if obj.SMFSFlag.Fit(ii)==1
            % Mean and standard deviation of a selection of the fitted y retantion data           
            obj.yRetFitMean(ii)=mean(obj.BasedRet{ii}(DataPtsRet(1)-LimitIdxRet2:DataPtsRet(1)-LimitIdxRet1,1)); % Calculate the mean of the difference data
            obj.yRetFitStd(ii)=std(obj.BasedRet{ii}(DataPtsRet(1)-LimitIdxRet2:DataPtsRet(1)-LimitIdxRet1,1));   % Calculate the standard deviation of the difference data
     
            else
            % Correction based on the approach data            
            obj.CorrMeanApp(ii)=mean(abs(obj.BasedApp{ii}(LimitIdxApp1:LimitIdxApp2,1))-abs(obj.BasedRet{ii}(DataPtsApp(1)-LimitIdxApp2:DataPtsApp(1)-LimitIdxApp1,1))); % Calculate the mean of the difference data
            obj.CorrStdApp(ii)=std(obj.BasedApp{ii}(DataPtsApp(1)-LimitIdxApp2:DataPtsApp(1)-LimitIdxApp1,1));             
            obj.BasedRetCorr{ii}=obj.BasedRet{ii}-obj.CorrMeanApp(ii); % Correct the BasedRet data with the mean of the correction data
            
            % Correction based on the retraction data                     
            obj.CorrMeanRet(ii)=mean(obj.BasedRet{ii}(DataPtsRet(1)-LimitIdxRet2:DataPtsRet(1)-LimitIdxRet1,1)); % Calculate the mean of the difference data
            obj.CorrStdRet(ii)=std(obj.BasedRet{ii}(DataPtsRet(1)-LimitIdxRet2:DataPtsRet(1)-LimitIdxRet1,1));   % Calculate the standard deviation of the difference data            
            obj.BasedRetCorr2{ii}=obj.BasedRet{ii}-obj.CorrMeanRet(ii); % Correct the BasedRet data with the mean of the correction data          
                        
            end
        end   
              
        %% Appendix
%         close all
%         % Define variables
%         kk=1
%         x100=-100e-9; % Defines 100nm
%         x500=-500e-9; % Defines 500nm
%         % Graphical preview
%         fig=gcf;
%         fig.Units='normalized'; % changes to normalized unit settings, necessary to receive the full screen size in the next line
%         fig.Color='white'; % changes the background color of the figure
%         fig.OuterPosition=[0.5 0 0.5 1];% changes the size of the figure to half screen
%         fig.PaperOrientation='landscape';
%         grid on
%         hold on
%         % "Origin" data
%         plot(a.FM{1}.THApp{kk}-a.FM{1}.CP_HardSurface(kk,1),a.FM{1}.BasedApp{kk},'b');
%         plot(a.FM{1}.THRet{kk}-a.FM{1}.CP_HardSurface(kk,1),a.FM{1}.BasedRet{kk},'r');
%         % Retention data corrected
%         plot(a.FM{1}.THRet{kk}-a.FM{1}.CP_HardSurface(kk,1),a.FM{ii}.BasedRetCorr{kk},'g');
%         % DataShare part of the data
%         plot(a.FM{1}.THApp{kk}(LimitIdx1:LimitIdx2,1)-a.FM{1}.CP_HardSurface(kk,1),a.FM{1}.BasedApp{kk}(LimitIdx1:LimitIdx2,1),'y');
%         plot(a.FM{1}.THRet{kk}(DataPts(1)-LimitIdx2:DataPts(1)-LimitIdx1,1)-a.FM{1}.CP_HardSurface(kk,1),a.FM{1}.BasedRet{kk}(DataPts(1)-LimitIdx2:DataPts(1)-LimitIdx1,1),'m');
%         % Markers
%         plot(a.FM{1}.THRet{kk}(ThreshIdx,1)-a.FM{1}.CP_HardSurface(kk,1),a.FM{1}.BasedRet{kk}(ThreshIdx,1),'kx','MarkerSize',20);
%         plot(a.FM{1}.THApp{kk}(ThreshIdx,1)-a.FM{1}.CP_HardSurface(kk,1),a.FM{1}.BasedApp{kk}(ThreshIdx,1),'mx','MarkerSize',20);
%         Distances lines
%         line([x100 x100], ylim,'Color','k'); % Draws a vertical line
%         line([x500 x500], ylim,'Color','k'); % Draws a vertical line
        end
                
        function fc_pulling_length_sigma(obj)
            sigma=1;
            for ii=1:obj.NCurves
            %for ii=55 % debugging
                if ~obj.SMFSFlag.Uncorrupt(jj) || ~obj.SMFSFlag.Min(jj)     % Exclude corrupted force curves from the analysis     
                continue
                end
                % Allocate data
                xApp=obj.THApp{ii}-obj.CP_HardSurface(ii);
                xRet=obj.THRet{ii}-obj.CP_HardSurface(ii);
                if obj.SMFSFlag.Fit(ii)==1
                yApp=obj.BasedApp{ii};
                yRet=obj.BasedRet{ii};
                else
                yApp=obj.BasedApp{ii};
                yRet=obj.BasedRetCorr2{ii};
                end

                % Find the index and determine the pulling length
                if obj.SMFSFlag.Fit(ii)==1
                obj.PullingLengthIdx(ii)=find(yRet<obj.yRetFitMean(ii)-obj.yRetFitStd(ii)*sigma,1,'last'); % Finds the index of the value that fulfils the condition         
                obj.PullingLength(ii)=abs(xRet(obj.PullingLengthIdx(ii))); % Corresponding x-value of the index
                else
                obj.PullingLengthIdx(ii)=find(yRet<obj.CorrMeanRet(ii)-obj.CorrStdApp(ii)*sigma,1,'last'); % Finds the index of the value that fulfils the condition         
                obj.PullingLength(ii)=abs(xRet(obj.PullingLengthIdx(ii))); % Corresponding x-value of the index
                end                 
            end           
            obj.FMPullingLengthMean=mean(obj.PullingLength);
            obj.FMPullingLengthMin=min(obj.PullingLength);
            obj.FMPullingLengthMax=max(obj.PullingLength);
            
%             %% Appendix
%             close all
%             % Define variables
%             qq=1;
%             % Allocate data
%             xApp=obj.THApp{qq}-obj.CP_HardSurface(qq);
%             xRet=obj.THRet{qq}-obj.CP_HardSurface(qq);
%             yApp=obj.BasedApp{qq};
%             yRet=obj.BasedRetCorr2{qq};    
%             % Graphical preview
%             h_fig=figure(1);
%             h_fig.Color='white'; % changes the background color of the figure
%             h_fig.Units='normalized'; % Defines the units
%             h_fig.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
%             h_fig.PaperOrientation='landscape';
%             %% Plotting the tiles
%             t = tiledlayout(2,2);
%             %t.TileSpacing = 'compact';
%             %t.Padding = 'compact';
%             t.TileSpacing = 'none'; % To reduce the spacing between the tiles
%             t.Padding = 'none'; % To reduce the padding of perimeter of a tile
%             nexttile
%             hold on
%             grid on
%             plot(xApp,yApp,'b');
%             plot(xRet,yRet,'r');
%             nexttile
%             hold on
%             grid on
%             plot(xApp,yApp,'b');
%             plot(xRet,yRet,'r');
%             plot(obj.THRet{qq}(obj.PullingLengthIdx(qq))-obj.CP_HardSurface(qq),obj.BasedRetCorr2{qq}(obj.PullingLengthIdx(qq)),'*','MarkerSize',10,'MarkerEdgeColor','g')              
        end        
        
        function fc_pulling_length_MAD(obj)
                      
             % Define variables
             beforePts=5;
             afterPts=10;
             movWindow=[beforePts afterPts]; % moving data point window for the moving median absolute deviation         
             yLimit1=0.045; 
             
             for jj=1:obj.NCurves
                if ~obj.SMFSFlag.Uncorrupt(jj) || ~obj.SMFSFlag.Min(jj)     % Exclude corrupted force curves from the analysis     
                continue
                end
                % Allocate data
                xRet{jj}=obj.HHRet{jj}-obj.CP_HardSurface(jj);
                yRet{jj}=obj.BasedRet{jj};
                % Determine the pulling length index 
                normRet{jj}=normalize(flip(yRet{jj}));  % Normalize and flip the data             
                normRetMAD{jj} = movmad(normRet{jj},movWindow); % Apply the moving median absolute deviation                 
                Peak1{jj}=find(normRetMAD{jj}>yLimit1,1,'first'); % Find the index in the data fulfilling the condition
                obj.PullingLengthIdx(jj)=length(yRet{jj})-Peak1{jj}; % Correct for the flipped data the peak index is based on by substracting from the number of data points
                obj.PullingLength(jj)=abs(xRet{jj}(obj.PullingLengthIdx(jj))); % Corresponding x-value of the index
             end           
            % Allocate data  
            obj.FMPullingLengthMean=mean(obj.PullingLength);
            obj.FMPullingLengthMin=min(obj.PullingLength);
            obj.FMPullingLengthMax=max(obj.PullingLength);
        
             %% Appendix
%             close all
%              % Graphical preview
%              for ii=1:obj.NCurves
%                 h_fig=figure(1);
%                 h_fig.Color='white'; % changes the background color of the figure
%                 h_fig.Units='normalized'; % Defines the units
%                 h_fig.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
%                 h_fig.PaperOrientation='landscape';
%                 subplot(3,1,1)
%                 hold on
%                 plot(flip(yRet{ii}))
%                 plot(obj.PullingLengthIdx(ii),flip(yRet{ii}(obj.PullingLengthIdx(ii))),'kx','MarkerSize',20)
%                 plot(Peak1{ii},flip(yRet{ii}(Peak1{ii})),'r*','MarkerSize',20)
%                 hold off
%                 ax21=gca;
%                 ax21.XLim = [100 inf];
%                 ax21.YLim = [-inf inf];
%                 subplot(3,1,2)
%                 hold on
%                 plot(xRet{ii},yRet{ii})
%                 plot(xRet{ii}(length(yRet{ii})-obj.PullingLengthIdx(ii)),yRet{ii}(length(yRet{ii})-obj.PullingLengthIdx(ii)),'kx','MarkerSize',20)
%                 plot(xRet{ii}(length(yRet{ii})-Peak1{ii}),yRet{ii}(length(yRet{ii})-Peak1{ii}),'r*','MarkerSize',20)
%                 hold off
%                 ax22=gca;
%                 ax22.XLim = [-inf inf];
%                 ax22.YLim = [-inf inf];
%                 subplot(3,1,3)
%                 hold on
%                 plot(normRetMAD{jj})
%                 plot(Peak1{ii},normRetMAD{ii}(Peak1{ii}),'kx','MarkerSize',20)
%                 hold off
%                 ax23=gca;
%                 ax23.XLim = [100 inf];
%                 ax23.YLim = [0 0.2];
%                 drawnow     
%                 pause(3)
%                 close figure 1
%              end
        end
             
        function fc_snap_in_length_MAD(obj)
                      
             % Define variables
             beforePts=75;
             afterPts=75;
             movWindow=[beforePts afterPts]; % moving data point window for the moving median absolute deviation        
             yLimit1=0.1; 
             
             for jj=1:obj.NCurves
                if ~obj.SMFSFlag.Uncorrupt(jj) || ~obj.SMFSFlag.Min(jj)     % Exclude corrupted force curves from the analysis     
                continue
                end
                % Allocate data               
                xApp{jj}=obj.HHApp{jj}-obj.CP_HardSurface(jj); 
                yApp{jj}=obj.BasedApp{jj};
                % Determine the pulling length index 
                normApp{jj}=normalize(flip(yApp{jj}));  % Normalize and flip the data             
                normAppMAD{jj} = movmad(normApp{jj},movWindow); % Apply the moving median absolute deviation                           
                Peak1{jj}=find(normAppMAD{jj}>yLimit1,1,'last'); % Find the index in the data fulfilling the condition
                obj.SnapInIdx(jj)=length(yApp{jj})-Peak1{jj}; % Correct for the data the peak index is based on by substracting from the number of data points               
                obj.SnapInLength(jj)=abs(xApp{jj}(obj.SnapInIdx(jj))); % Corresponding x-value of the index
             end
            % Allocate data  
            obj.FMSnapInMean=mean(obj.PullingLength);
            obj.FMSnapInMin=min(obj.PullingLength);
            obj.FMSnapInMax=max(obj.PullingLength);
        
             %% Appendix
%            close all
%             % Graphical preview
%             % Figure 1
%              for ii=1:obj.NCurves
%                 h_fig=figure(1);
%                 h_fig.Color='white'; % changes the background color of the figure
%                 h_fig.Units='normalized'; % Defines the units
%                 h_fig.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
%                 h_fig.PaperOrientation='landscape';
%                 subplot(3,1,1)
%                 hold on
%                 plot(yApp{ii})
%                 plot(obj.SnapInIdx(ii),yApp{ii}(obj.SnapInIdx(ii)),'kx','MarkerSize',20)
%                 plot(Peak1{ii},yApp{ii}(Peak1{ii}),'r*','MarkerSize',20)
%                 hold off
%                 ax21=gca;
%                 ax21.XLim = [100 inf];
%                 ax21.YLim = [-inf inf];
%                 subplot(3,1,2)
%                 hold on
%                 plot(xApp{ii},yApp{ii})
%                 plot(xApp{ii}(obj.SnapInIdx(ii)),yApp{ii}(obj.SnapInIdx(ii)),'kx','MarkerSize',20)
%                 plot(xApp{ii}(Peak1{ii}),yApp{ii}(Peak1{ii}),'r*','MarkerSize',20)
%                 hold off
%                 ax22=gca;
%                 ax22.XLim = [-inf inf];
%                 ax22.YLim = [0 2];
%                 subplot(3,1,3)
%                 hold on
%                 plot(normAppMAD{ii})
%                 plot(Peak2{ii},normAppMAD{ii}(Peak1{ii}),'kx','MarkerSize',20)
%                 hold off
%                 ax23=gca;
%                 ax23.XLim = [-inf inf];
%                 ax23.YLim = [0 2];
%                 drawnow     
%                 pause(3)
%                 close figure 1
%              end
            % Figure 2
            %  Test different windows
%                for jj=1:obj.NCurves
%                   movWindow1=[50 50]; % moving data point window for the moving median absolute deviation      
%                   normAppMAD1{jj} = movmad(normApp{jj},movWindow1); % Apply the moving median absolute deviation   
%                   Peak2{jj}=find(normAppMAD1{jj}>yLimit1,1,'last'); % Find the index in the data fulfilling the condition
%                   SnapInIdx(jj)=length(yApp{jj})-Peak2{jj}; % Correct for the data the peak index is based on by substracting from the number of data points                  
%                end
%               for ii=1:obj.NCurves
%                 h_fig=figure(2);
%                 h_fig.Color='white'; % changes the background color of the figure
%                 h_fig.Units='normalized'; % Defines the units
%                 h_fig.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
%                 h_fig.PaperOrientation='landscape';
%                 subplot(3,1,1)
%                 hold on
%                 plot(flip(yApp{ii}))              
%              %   plot(obj.SnapInIdx(ii),flip(yApp{ii}(obj.SnapInIdx(ii))),'kx','MarkerSize',20)
%              %   plot(SnapInIdx(ii),flip(yApp{ii}(SnapInIdx(ii))),'r*','MarkerSize',20)
%                 hold off
%                 ax21=gca;
%                 ax21.XLim = [-inf inf];
%                 ax21.YLim = [-inf inf];
%                 subplot(3,1,2)
%                 hold on
%                 plot(normAppMAD{ii})
%                 plot(Peak1{ii},normAppMAD{ii}(Peak1{ii}),'kx','MarkerSize',20)
%                 hold off
%                 ax22=gca;
%                 ax22.XLim = [-inf inf];
%                 ax22.YLim = [0 0.7];
%                 subplot(3,1,3)
%                 hold on
%                 plot(normAppMAD1{ii})
%                 plot(Peak2{ii},normAppMAD1{ii}(Peak2{ii}),'kx','MarkerSize',20)
%                 hold off
%                 ax23=gca;
%                 ax23.XLim = [-inf inf];
%                 ax23.YLim = [0 0.7];
%                 drawnow     
%                 pause(3)
%                 close figure 2
%              end
            
        end
        
        function fc_find_idx(obj)
            % Finds the corresponding indices of an array of data points
            
            % Define the corresponding data points            
            xLimit=-50e-9;  % Value on the x-axis
            xDistance=10e-9;
            % Loop
            for ii=1:obj.NCurves
            if ~obj.SMFSFlag.Uncorrupt(ii) || ~obj.SMFSFlag.Min(ii)     % Exclude corrupted force curves from the analysis     
                continue
            end
            % Allocate data
            xRet=obj.HHRet{ii}-obj.CP_HardSurface(ii);
            if obj.SMFSFlag.Fit(ii)==1
            yApp=obj.BasedApp{ii};
            yRet=obj.BasedRet{ii};
            else
            yApp=obj.BasedApp{ii};
            yRet=obj.BasedRetCorr2{ii};
            end      
            % 50 nm index                          
            EndLogical=xRet<xLimit; % Determine the elements that fulfil the logical argument           
            EndIdx(ii)=find(EndLogical,1,'first'); % Read out the index of the first cell that fulfil the argument
            % Allocate data
            obj.Idx50nm(ii)=EndIdx(ii);
            % 10 nm before pulling length index
            xUnbdingBoundary=obj.PullingLength(ii)-xDistance;
            [~,obj.UnbindingBoundaryIdx(ii)]=min(abs(abs(xRet)-xUnbdingBoundary));         
            end            
        end
        
        function fc_adh_force_max(obj)
            % Function to find the maximum adhesion force value in the
            % approach and retraction data of a force curve
            
            for ii=1:100
            % for ii=55:100 % for debugging
            if ~obj.SMFSFlag.Uncorrupt(ii) || ~obj.SMFSFlag.Min(ii)     % Exclude corrupted force curves from the analysis     
                continue
            end
            % Allocate data
            %xApp=obj.THApp{ii}-obj.CP_HardSurface(ii);
            %xRet=obj.THRet{ii}-obj.CP_HardSurface(ii);            
            if obj.SMFSFlag.Fit(ii)==1
            yApp=obj.BasedApp{ii};
            yRet=obj.BasedRet{ii};
            else
            yApp=obj.BasedApp{ii};
            yRet=obj.BasedRetCorr2{ii};
            end           
            % Determine the maximum adhesion force 
            obj.AdhForceMaxApp(ii)=min(yApp(1:end)); % Determine maximum adhesion forces from the pulling length index to the last data point
            obj.AdhForceMaxRet(ii)=min(yRet(obj.Idx50nm(ii):obj.PullingLengthIdx(ii))); % Determine maximum adhesion forces from the 50nm index to the pulling length index
            obj.AdhForceUnbinding(ii)=min(yRet(obj.UnbindingBoundaryIdx(ii):obj.PullingLengthIdx(ii))); % Determine maximum adhesion forces close to the pulling length
            % Maximum adhesion force in a predefined range before the
            % pulling length
            % Determine the corresponding indices
            obj.AdhForceMaxAppIdx(ii)=find(yApp==obj.AdhForceMaxApp(ii)); % Finds the index of the value that fulfils the condition         
            obj.AdhForceMaxRetIdx(ii)=find(yRet==obj.AdhForceMaxRet(ii)); % Finds the index of the value that fulfils the condition          
            UnbindingBoundaryIdx=find(yRet(obj.UnbindingBoundaryIdx(ii):obj.PullingLengthIdx(ii))==obj.AdhForceUnbinding(ii)); % Finds the index of the value that fulfils the condition
            obj.AdhForceUnbindingIdx(ii)=obj.UnbindingBoundaryIdx(ii)+UnbindingBoundaryIdx; % Correct for the shifted starting position
            end
           
%             %% Appendix
%             close all
%             % Define variables
%             qq=1;
%             % Allocate data
%             xApp=obj.THApp{qq}-obj.CP_HardSurface(qq);
%             xRet=obj.THRet{qq}-obj.CP_HardSurface(qq);
%             yApp=obj.BasedApp{qq};
%             yRet=obj.BasedRetCorr2{qq};    
%             % Graphical preview
%             h_fig=figure(1);
%             h_fig.Color='white'; % changes the background color of the figure
%             h_fig.Units='normalized'; % Defines the units
%             h_fig.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
%             h_fig.PaperOrientation='landscape';
%             %% Plotting the tiles
%             t = tiledlayout(2,2);
%             %t.TileSpacing = 'compact';
%             %t.Padding = 'compact';
%             t.TileSpacing = 'none'; % To reduce the spacing between the tiles
%             t.Padding = 'none'; % To reduce the padding of perimeter of a tile
%             nexttile
%             hold on
%             grid on
%             plot(xApp,yApp,'b');
%             plot(xRet,yRet,'r');
%             nexttile
%             hold on
%             grid on
%             plot(xApp,yApp,'b');
%             plot(xRet,yRet,'r');
%             plot(obj.THRet{qq}(obj.AdhForceMaxRetIdx(qq))-obj.CP_HardSurface(qq),obj.BasedRetCorr2{qq}(obj.AdhForceMaxRetIdx(qq)),'*','MarkerSize',10,'MarkerEdgeColor','g')              
%             plot(obj.THApp{qq}(obj.AdhForceMaxAppIdx(qq))-obj.CP_HardSurface(qq),obj.BasedApp{qq}(obj.AdhForceMaxAppIdx(qq)),'*','MarkerSize',10,'MarkerEdgeColor','m')              
%             nexttile
%             hold on
%             grid on
%             plot(xApp,yApp,'b');
%             plot(xRet,yRet,'r');
%             plot(xRet(obj.AdhForceMaxRetIdx(qq)),yRet(obj.AdhForceMaxRetIdx(qq)),'*','MarkerSize',10,'MarkerEdgeColor','g')              
%             plot(xApp(obj.AdhForceMaxAppIdx(qq)),yApp(obj.AdhForceMaxAppIdx(qq)),'*','MarkerSize',10,'MarkerEdgeColor','m')              
       
        end
        
        function fc_adhesion_energy_threshold(obj)
            % Determine the adhesion energy using a predefined force
            % threshold to distinguish interactions from background noise
            
            %% Loop over all force curves
            for ii=1:obj.NCurves
            %for ii=97 % For debugging and testing 
                if ~obj.SMFSFlag.Uncorrupt(ii) || ~obj.SMFSFlag.Min(ii)     % Exclude corrupted force curves from the analysis     
                continue
                end            
                % Allocate data
                xRet=obj.HHRet{ii}-obj.CP_HardSurface(ii); % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method 
                xApp=obj.HHApp{ii}-obj.CP_HardSurface(ii); % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method 
                if obj.SMFSFlag.Fit(ii)==1
                yApp=obj.BasedApp{ii};
                yRet=obj.BasedRet{ii};
                else
                yApp=obj.BasedApp{ii};
                yRet=obj.BasedRetCorr2{ii};
                end
                % Define variables
                limit1=0;   % Define the limit
                % Apply the limit
                % Approach
                yApp(yApp>Limit1)=0;
                yApp(obj.PullingLengthIdx(ii):end)=0; % Set all data points with a higher index (surface distance is higher) than the pulling length index to 0
                % Retention 
                yRet(yRet>limit1)=0;  % Set all values above the zero line of the x-axis 0
                yRet(obj.PullingLengthIdx(ii):end)=0; % Set all data points with a higher index (surface distance is higher) than the pulling length index to 0
                % Allocate data
                obj.yRetLim{ii}=yRet; 
                obj.yAppLim{ii}=yApp; 
                % Determine the adhesion energy
                IntApp(ii)=trapz(yApp,xApp); % Integrates over the modified y-retraction data with respect to the corresponding x-retraction data 
                obj.AppAdhEnergy_IdxMethod(ii)=IntApp(ii);
                IntRet(ii)=trapz(yRet,xRet); % Integrates over the modified y-retraction data with respect to the corresponding x-retraction data 
                obj.RetAdhEnergy_IdxMethod(ii)=IntRet(ii);
            end
                obj.FMAppAdhEnergyMean=mean(obj.AppAdhEnergy_IdxMethod);
                obj.FMAppAdhEnergyStd=std(obj.AppAdhEnergy_IdxMethod); 
                obj.FMRetAdhEnergyMean=mean(obj.RetAdhEnergy_IdxMethod);
                obj.FMRetAdhEnergyStd=std(obj.RetAdhEnergy_IdxMethod);      
                      
%             % %% Appendix
%             close all
%             % Testing
%             % Allocate data
%             xApp=obj.THApp{ii}-obj.CP_HardSurface(ii); % Approach x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method 
%             yApp=obj.BasedApp{ii};             
%             % Define variables
%             limit2=-20e-12;  % Define the limit
%             yRetLimTest=obj.yRetLim{ii};           
%             yRetLimTest(yRetLimTest>limit2)=0;
%             yRetLimTest(1:obj.PullingLengthIdx(ii))=-100e-12; % Set all data points with a lower index (surface distance is lower and lies within origin and pulling length) than the pulling length index to a constant value
%             % Determine the adhesion energy
%             AdhEne2=trapz(yRetLimTest,xRet); % Integrate the modified retraction data
%             AdhEneCum2=cumtrapz(yRetLimTest,xRet);
%             % Graphical preview
%             h_fig=figure(1);
%             h_fig.Color='white'; % changes the background color of the figure
%             h_fig.Units='normalized'; % Defines the units
%             h_fig.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
%             h_fig.PaperOrientation='landscape';
%             %% Plotting the tiles
%             t = tiledlayout(2,2);
%             %t.TileSpacing = 'compact';
%             %t.Padding = 'compact';
%             t.TileSpacing = 'none'; % To reduce the spacing between the tiles
%             t.Padding = 'none'; % To reduce the padding of perimeter of a tile
%             % Tile 1
%             nexttile
%             hold on
%             grid on
%             plot(xApp,yApp,'g')
%             plot(xRet,obj.BasedRet{ii},'m')
%             plot(xRet,yRetLim,'b')
%             plot(xRet,yRetLimTest,'k')
%             plot(xRet(obj.PullingLengthIdx(ii)),yRetLim(obj.PullingLengthIdx(ii)),'*','MarkerSize',10,'MarkerEdgeColor','r')
%             legend('Approach data','Retraction data, y-data corrected','Retraction data, y-data corrected, limits included','Test Retraction data, y-data corrected, limits included, levelled','Pulling length position')
%             % Tile 2
%             nexttile;
%             plot(AdhEneCum2)
%             title('Cumulative adhesion energy of "yRetLimTest"')
%             % Tile 3
%             nexttile;
%             hold on
%             plot(xApp,yApp,'g');
%             plot(xRet,yRetLimTest,'k');
%             plot(xRet(obj.PullingLengthIdx(ii)),yRetLim(obj.PullingLengthIdx(ii)),'*','MarkerSize',10,'MarkerEdgeColor','r') 
%             legend('Approach x-data','Test Retraction x-data corrected, limits included, levelled','Pulling length position')     
%             % Tile 4
%             nexttile
%             hold on
%             grid on
%             xlim([-inf 10e-9]) % Define the x axes limits [min max]
%             ylim([-inf 200e-12]) % Define the y axes limits [min max]
%             plot(xApp,yApp,'g');
%             plot(xRet,yRetLim,'b');
%             area(xRet,obj.BasedRetCorr2{ii},'FaceColor','y')
%             plot(xRet(obj.PullingLengthIdx(ii)),yRetLim(obj.PullingLengthIdx(ii)),'*','MarkerSize',10,'MarkerEdgeColor','r')
%             area(xRet(1:obj.PullingLengthIdx(ii)),yRetLim(1:obj.PullingLengthIdx(ii)),'FaceColor','c') % Highlights the area with starting and end points: Pulling length and origin
%             legend('Approach data','Retraction data, y-data corrected, limits included','Adhesion force based on Retraction data, y-data corrected','Pulling length position','Adhesion force based on Retraction data, y-data corrected, limits included')    
           end
        
        function fc_adhesion_energy_idxpulllength(obj)
            % Determine the adhesion energy using the previous defined pulling length index
            
            %% Loop over all force curves
                for ii=1:obj.NCurves
            %for ii=97 % % For debugging and testing
                if ~obj.SMFSFlag.Uncorrupt(ii) || ~obj.SMFSFlag.Min(ii)     % Exclude corrupted force curves from the analysis     
                continue
                end               
                % Allocate data
                xRet=obj.HHRet{ii}-obj.CP_HardSurface(ii); % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method 
                xApp=obj.HHApp{ii}-obj.CP_HardSurface(ii); % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method 
                if obj.SMFSFlag.Fit(ii)==1
                yApp=obj.BasedApp{ii};
                yRet=obj.BasedRet{ii};
                else
                yApp=obj.BasedApp{ii};
                yRet=obj.BasedRetCorr2{ii};
                end
                % Define variables
                limit1=0;   % Define the limit
                % Apply the limit
                % Approach
                yApp(yApp>limit1)=0; % Set all data points representing surface contact (set point related, above the y-axis 0 line) to 0
                yApp(1:obj.PullingLengthIdx(ii))=0; % Set all data points with a higher index (surface distance is higher) than the pulling length index to 0
                % Retention 
                yRet(yRet>limit1)=0;  % Set all data points representing surface contact (set point related, above the y-axis 0 line) to 0
                yRet(obj.PullingLengthIdx(ii):end)=0; % Set all data points with a higher index (surface distance is higher) than the pulling length index to 0
                % Allocate data
                obj.yRetLim{ii}=yRet; 
                obj.yAppLim{ii}=yApp; 
                % Determine the adhesion energy
                IntApp(ii)=trapz(yApp,xApp); % Integrates over the modified y-retraction data with respect to the corresponding x-retraction data 
                obj.AppAdhEnergy_IdxMethod(ii)=IntApp(ii);
                IntRet(ii)=trapz(yRet,xRet); % Integrates over the modified y-retraction data with respect to the corresponding x-retraction data 
                obj.RetAdhEnergy_IdxMethod(ii)=IntRet(ii);
                end
                obj.FMAppAdhEnergyMean=mean(obj.AppAdhEnergy_IdxMethod);
                obj.FMAppAdhEnergyStd=std(obj.AppAdhEnergy_IdxMethod); 
                obj.FMRetAdhEnergyMean=mean(obj.RetAdhEnergy_IdxMethod);
                obj.FMRetAdhEnergyStd=std(obj.RetAdhEnergy_IdxMethod);       
            
%             % %% Appendix
%             close all
%             % Testing
%             % Allocate data
%             xApp=obj.THApp{ii}-obj.CP_HardSurface(ii); % Approach x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method 
%             yApp=obj.BasedApp{ii};             
%             % Define variables
%             limit2=-20e-12;  % Define the limit
%             yRetLimTest=obj.yRetLim{ii};           
%             yRetLimTest(yRetLimTest>limit2)=0;
%             yRetLimTest(1:obj.PullingLengthIdx(ii))=-100e-12; % Set all data points with a lower index (surface distance is lower and lies within origin and pulling length) than the pulling length index to a constant value
%             % Determine the adhesion energy
%             AdhEne2=trapz(yRetLimTest,xRet); % Integrate the modified retraction data
%             AdhEneCum2=cumtrapz(yRetLimTest,xRet);
%             % Graphical preview
%             h_fig=figure(1);
%             h_fig.Color='white'; % changes the background color of the figure
%             h_fig.Units='normalized'; % Defines the units
%             h_fig.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
%             h_fig.PaperOrientation='landscape';
%             %% Plotting the tiles
%             t = tiledlayout(2,2);
%             %t.TileSpacing = 'compact';
%             %t.Padding = 'compact';
%             t.TileSpacing = 'none'; % To reduce the spacing between the tiles
%             t.Padding = 'none'; % To reduce the padding of perimeter of a tile
%             % Tile 1
%             nexttile
%             hold on
%             grid on
%             plot(xApp,yApp,'g')
%             plot(xRet,obj.BasedRet{ii},'m')
%             plot(xRet,yRetLim,'b')
%             plot(xRet,yRetLimTest,'k')
%             plot(xRet(obj.PullingLengthIdx(ii)),yRetLim(obj.PullingLengthIdx(ii)),'*','MarkerSize',10,'MarkerEdgeColor','r')
%             legend('Approach data','Retraction data, y-data corrected','Retraction data, y-data corrected, limits included','Test Retraction data, y-data corrected, limits included, levelled','Pulling length position')
%             % Tile 2
%             nexttile;
%             plot(AdhEneCum2)
%             title('Cumulative adhesion energy of "yRetLimTest"')
%             % Tile 3
%             nexttile;
%             hold on
%             plot(xApp,yApp,'g');
%             plot(xRet,yRetLimTest,'k');
%             plot(xRet(obj.PullingLengthIdx(ii)),yRetLim(obj.PullingLengthIdx(ii)),'*','MarkerSize',10,'MarkerEdgeColor','r') 
%             legend('Approach x-data','Test Retraction x-data corrected, limits included, levelled','Pulling length position')     
%             % Tile 4
%             nexttile
%             hold on
%             grid on
%             xlim([-inf 10e-9]) % Define the x axes limits [min max]
%             ylim([-inf 200e-12]) % Define the y axes limits [min max]
%             plot(xApp,yApp,'g');
%             plot(xRet,yRetLim,'b');
%             area(xRet,obj.BasedRetCorr2{ii},'FaceColor','y')
%             plot(xRet(obj.PullingLengthIdx(ii)),yRetLim(obj.PullingLengthIdx(ii)),'*','MarkerSize',10,'MarkerEdgeColor','r')
%             area(xRet(1:obj.PullingLengthIdx(ii)),yRetLim(1:obj.PullingLengthIdx(ii)),'FaceColor','c') % Highlights the area with starting and end points: Pulling length and origin
%             legend('Approach data','Retraction data, y-data corrected, limits included','Adhesion force based on Retraction data, y-data corrected','Pulling length position','Adhesion force based on Retraction data, y-data corrected, limits included')    
        end
                 
        function fc_print_properties(obj,XMin,XMax,YMin,YMax,NumFcMax,NumFcUncorrupt,hh) % fc ... force curve
            % fc_print_adhenergy_pulllength: A function to plot all selected force curves of a
            % force map including adhesion energy and pulling length in
            % each force curve
            if nargin < 2
                XMin= -inf;
                XMax= inf;
                YMin= -inf;
                YMax= inf;
            end
            % Define variables
            RGB1=[0 26 255]./255;  % Blue 
            RGB2=[255 119 0]./255; % Orange
            RGB7=[255 230 0]./255; % Yellow
            RGB8=[80 200 204]./255; % Turquoise
            RGB10=[200 0 255]./255; % Violet
            RGB11=[200 255 150]./255; % Light Green
            RGB12=[185 230 254]./255; % Light Blue
            RGB13=[200 0 0]./255; % Red
            % Define variables for the figure name
            ExtendVelocityConvert=num2str(obj.ExtendVelocity*1e9);
            RetractVelocityConvert=num2str(obj.RetractVelocity*1e9);
            % Classification criteria
            figname=strcat(obj.Date,{'_'},obj.Time,{'_'},obj.ID,{'_'},obj.Substrate,{'_'},obj.EnvCond,{'_'},{'_'},obj.Chipbox,{'_'},obj.ChipCant,{'_'},ExtendVelocityConvert,{'_'},RetractVelocityConvert);
            figname=char(figname); 
            % Define variables for the plot loop
            mm=ceil(sqrt(NumFcMax));
            nn=mm;
            ww=1; % "flag while loop" variable
            DiffFc=0;
            NumFcUncorrupt=nnz(obj.SMFSFlag.Uncorrupt.*obj.SMFSFlag.Min); % Determine the amount of uncorrupted force curves            
            NumFigures=ceil(NumFcUncorrupt/NumFcMax);
            if NumFigures==0     % If condition is fulfilled stop function and return to calling function     
                return              
            end 
            RemainderMax=mod(NumFcUncorrupt,NumFcMax); % Check for remainder           
            if RemainderMax ~= 0
                oo=round(sqrt(RemainderMax)); % Determine the number of rows in the figure
                pp=ceil(sqrt(RemainderMax)); % Determine the number of columns in the figure
                RemainderReal=mod(NumFcUncorrupt,oo*pp); % Correct the remainder based on the determined rows times columns
            end
            %% figure loop
            for ii=1:NumFigures               
                % Define variables
                jj=1; % "force curve plotted per figure" variable        
                % Figure
                h_fig=figure(ii);
                h_fig.Color='white'; % changes the background color of the figure
                h_fig.Units='normalized'; % Defines the units
                h_fig.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
                h_fig.PaperOrientation='landscape';
                h_fig.Name=figname;
                %% Plotting the tiles
                if RemainderMax == 0
                    t = tiledlayout(mm,nn);
                    t.TileSpacing = 'none'; % To reduce the spacing between the tiles
                    t.Padding = 'none'; % To reduce the padding of perimeter of a tile
                    if ii==1
                        NumFcCorSelec(ii)=nnz(~obj.SMFSFlag.Uncorrupt(1:NumFcMax));
                    elseif ii==2
                        NumFcCorSelec(ii)=nnz(obj.SMFSFlag.Uncorrupt((NumFcMax+1):(NumFcMax*(ii))));
                    else
                        NumFcCorSelec(ii)=nnz(obj.SMFSFlag.Uncorrupt((NumFcMax*(ii-1)+1):(NumFcMax*(ii))));
                    end
                    
                    if ii==1
                        kk=jj;
                    else
                        kk=jj+mm*nn*(ii-1);
                    end
                    %% Plot loop
                    for qq=kk:obj.NCurves % Loop over all force curves in the force map
                        if ww<qq+DiffFc
                            ww=qq+DiffFc;
                        end
                        while ~obj.SMFSFlag.Uncorrupt(ww)     % Stay in the while loop as long as the entry is zero
                            ww=ww+1;
                            if ww>qq
                                DiffFc=ww-qq;
                            end
                        end
                        % if condition
                        if ww>qq
                            ax=nexttile;
                            ax.XLim = [XMin XMax];
                            ax.YLim = [YMin YMax];                                     
                            if obj.SMFSFlag.Fit(qq+DiffFc)==1                            
                            grid on
                            hold on
                            area(obj.HHRet{qq+DiffFc}(1:obj.PullingLengthIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.yRetLim{qq+DiffFc}(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                            area(obj.HHApp{qq+DiffFc}(obj.SnapInIdx(qq+DiffFc):end)-obj.CP_HardSurface(qq+DiffFc),obj.yAppLim{qq+DiffFc}(obj.SnapInIdx(qq+DiffFc):end),'FaceColor',RGB7)
                            plot(obj.HHApp{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc),obj.BasedApp{qq+DiffFc},'Color',RGB1);
                            plot(obj.HHRet{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc),obj.BasedRet{qq+DiffFc},'Color',RGB2);
                            plot(obj.HHRet{qq+DiffFc}(obj.PullingLengthIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedRet{qq+DiffFc}(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)                          
                            plot(obj.HHApp{qq+DiffFc}(obj.SnapInIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedApp{qq+DiffFc}(obj.SnapInIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB10,'MarkerEdgeColor',RGB10)                                           
                            plot(obj.HHRet{qq+DiffFc}(obj.AdhForceMaxRetIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedRet{qq+DiffFc}(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)              
                            plot(obj.HHApp{qq+DiffFc}(obj.AdhForceMaxAppIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedApp{qq+DiffFc}(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)                                                                   
                            plot(obj.HHRet{qq+DiffFc}(obj.AdhForceUnbindingIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedRet{qq+DiffFc}(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)                                                                     
                            else
                            area(obj.HHRet{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc),obj.yRetLim{qq+DiffFc},'FaceColor',RGB12)
                            plot(obj.HHApp{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc),obj.BasedApp{qq+DiffFc},'Color',RGB1);
                            plot(obj.HHRet{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc),obj.BasedRetCorr2{qq+DiffFc},'Color',RGB2);
                            plot(obj.HHRet{qq+DiffFc}(obj.PullingLengthIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedRetCorr2{qq+DiffFc}(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','b')                          
                            plot(obj.HHRet{qq+DiffFc}(obj.AdhForceMaxRetIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedRetCorr2{qq+DiffFc}(obj.AdhForceMaxRetIdx(qq+DiffFc)),'p','MarkerSize',10,'MarkerFaceColor','g','MarkerEdgeColor','g')              
                            plot(obj.HHApp{qq+DiffFc}(obj.AdhForceMaxAppIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedApp{qq+DiffFc}(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB7,'MarkerEdgeColor',RGB7)                                                                                            
                            end
                            
                            % Title for each Subplot
                            ti=title(sprintf('%i',qq+DiffFc),'Color','k');
                            ti.Units='normalized'; % Set units to 'normalized'
                            ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                        else
                            ax=nexttile;
                            ax.XLim = [XMin XMax];
                            ax.YLim = [YMin YMax];
                            if obj.SMFSFlag.Fit(qq+DiffFc)==1                            
                            grid on
                            hold on
                            area(obj.HHRet{qq+DiffFc}(1:obj.PullingLengthIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.yRetLim{qq+DiffFc}(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                            area(obj.HHApp{qq+DiffFc}(obj.SnapInIdx(qq+DiffFc):end)-obj.CP_HardSurface(qq+DiffFc),obj.yAppLim{qq+DiffFc}(obj.SnapInIdx(qq+DiffFc):end),'FaceColor',RGB7)
                            plot(obj.HHApp{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc),obj.BasedApp{qq+DiffFc},'Color',RGB1);
                            plot(obj.HHRet{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc),obj.BasedRet{qq+DiffFc},'Color',RGB2);
                            plot(obj.HHRet{qq+DiffFc}(obj.PullingLengthIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedRet{qq+DiffFc}(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)                          
                            plot(obj.HHApp{qq+DiffFc}(obj.SnapInIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedApp{qq+DiffFc}(obj.SnapInIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB10,'MarkerEdgeColor',RGB10)                                           
                            plot(obj.HHRet{qq+DiffFc}(obj.AdhForceMaxRetIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedRet{qq+DiffFc}(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)              
                            plot(obj.HHApp{qq+DiffFc}(obj.AdhForceMaxAppIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedApp{qq+DiffFc}(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)                                                                   
                            plot(obj.HHRet{qq+DiffFc}(obj.AdhForceUnbindingIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedRet{qq+DiffFc}(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)                                                                     
                            else
                            area(obj.HHRet{qq}-obj.CP_HardSurface(qq),obj.yRetLim{qq},'FaceColor',RGB12)
                            plot(obj.HHApp{qq}-obj.CP_HardSurface(qq),obj.BasedApp{qq},'Color',RGB1);
                            plot(obj.HHRet{qq}-obj.CP_HardSurface(qq),obj.BasedRetCorr2{qq},'Color',RGB2);
                            plot(obj.HHRet{qq}(obj.PullingLengthIdx(qq))-obj.CP_HardSurface(qq),obj.BasedRetCorr2{qq}(obj.PullingLengthIdx(qq)),'d','MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','b')                          
                            plot(obj.HHRet{qq}(obj.AdhForceMaxRetIdx(qq))-obj.CP_HardSurface(qq),obj.BasedRetCorr2{qq}(obj.AdhForceMaxRetIdx(qq)),'p','MarkerSize',10,'MarkerFaceColor','g','MarkerEdgeColor','g')              
                            plot(obj.HHApp{qq}(obj.AdhForceMaxAppIdx(qq))-obj.CP_HardSurface(qq),obj.BasedApp{qq}(obj.AdhForceMaxAppIdx(qq)),'h','MarkerSize',10,'MarkerFaceColor',RGB7,'MarkerEdgeColor',RGB7)                                                                                            
                            end
                            % Title for each Subplot
                            ti=title(sprintf('%i',qq),'Color','k');
                            ti.Units='normalized'; % Set units to 'normalized'
                            ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                        end
                        if jj == NumFcMax
                            break
                        end
                        jj=jj+1;
                    end
                    
                else
                    if ii~=NumFigures
                        t = tiledlayout(mm,nn);
                        t.TileSpacing = 'none'; % To reduce the spacing between the tiles
                        t.Padding = 'none'; % To reduce the padding of perimeter of a tile
                        if ii==1
                            NumFcCorSelec(ii)=nnz(~obj.SMFSFlag.Uncorrupt(1:NumFcMax));
                        elseif ii==2
                            NumFcCorSelec(ii)=nnz(obj.SMFSFlag.Uncorrupt((NumFcMax+1):(NumFcMax*(ii))));
                        else
                            NumFcCorSelec(ii)=nnz(obj.SMFSFlag.Uncorrupt((NumFcMax*(ii-1)+1):(NumFcMax*(ii))));
                        end
                        
                        if ii==1
                            kk=jj;
                        else
                            kk=jj+mm*nn*(ii-1);
                        end
                        %% Plot loop
                        for qq=kk:obj.NCurves % Loop over all force curves in the force map
                            if ww<qq+DiffFc
                                ww=qq+DiffFc;
                            end
                            while ~obj.SMFSFlag.Uncorrupt(ww)     % Stay in the while loop as long as the entry is zero
                                ww=ww+1;
                                if ww>qq
                                    DiffFc=ww-qq;
                                end
                            end
                            if ww>qq
                                ax=nexttile;
                                ax.XLim = [XMin XMax];
                                ax.YLim = [YMin YMax];
                                if obj.SMFSFlag.Fit(qq+DiffFc)==1                            
                                grid on
                                hold on
                                area(obj.HHRet{qq+DiffFc}(1:obj.PullingLengthIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.yRetLim{qq+DiffFc}(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                area(obj.HHApp{qq+DiffFc}(obj.SnapInIdx(qq+DiffFc):end)-obj.CP_HardSurface(qq+DiffFc),obj.yAppLim{qq+DiffFc}(obj.SnapInIdx(qq+DiffFc):end),'FaceColor',RGB7)
                                plot(obj.HHApp{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc),obj.BasedApp{qq+DiffFc},'Color',RGB1);
                                plot(obj.HHRet{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc),obj.BasedRet{qq+DiffFc},'Color',RGB2);
                                plot(obj.HHRet{qq+DiffFc}(obj.PullingLengthIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedRet{qq+DiffFc}(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)
                                plot(obj.HHApp{qq+DiffFc}(obj.SnapInIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedApp{qq+DiffFc}(obj.SnapInIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB10,'MarkerEdgeColor',RGB10)
                                plot(obj.HHRet{qq+DiffFc}(obj.AdhForceMaxRetIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedRet{qq+DiffFc}(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                plot(obj.HHApp{qq+DiffFc}(obj.AdhForceMaxAppIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedApp{qq+DiffFc}(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)
                                plot(obj.HHRet{qq+DiffFc}(obj.AdhForceUnbindingIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedRet{qq+DiffFc}(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)                                                                     
                                else
                                area(obj.HHRet{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc),obj.yRetLim{qq+DiffFc},'FaceColor',RGB12)
                                plot(obj.HHApp{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc),obj.BasedApp{qq+DiffFc},'Color',RGB1);
                                plot(obj.HHRet{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc),obj.BasedRetCorr2{qq+DiffFc},'Color',RGB2);
                                plot(obj.HHRet{qq+DiffFc}(obj.PullingLengthIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedRetCorr2{qq+DiffFc}(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','b')                          
                                plot(obj.HHRet{qq+DiffFc}(obj.AdhForceMaxRetIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedRetCorr2{qq+DiffFc}(obj.AdhForceMaxRetIdx(qq+DiffFc)),'p','MarkerSize',10,'MarkerFaceColor','g','MarkerEdgeColor','g')              
                                plot(obj.HHApp{qq+DiffFc}(obj.AdhForceMaxAppIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedApp{qq+DiffFc}(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB7,'MarkerEdgeColor',RGB7)                                                                                            
                                end
                                % Title for each Subplot
                                ti=title(sprintf('%i',qq+DiffFc),'Color','k');
                                ti.Units='normalized'; % Set units to 'normalized'
                                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                            else
                                ax=nexttile;
                                ax.XLim = [XMin XMax];
                                ax.YLim = [YMin YMax];
                                if obj.SMFSFlag.Fit(qq+DiffFc)==1                            
                                grid on
                                hold on
                                area(obj.HHRet{qq+DiffFc}(1:obj.PullingLengthIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.yRetLim{qq+DiffFc}(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                area(obj.HHApp{qq+DiffFc}(obj.SnapInIdx(qq+DiffFc):end)-obj.CP_HardSurface(qq+DiffFc),obj.yAppLim{qq+DiffFc}(obj.SnapInIdx(qq+DiffFc):end),'FaceColor',RGB7)
                                plot(obj.HHApp{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc),obj.BasedApp{qq+DiffFc},'Color',RGB1);
                                plot(obj.HHRet{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc),obj.BasedRet{qq+DiffFc},'Color',RGB2);
                                plot(obj.HHRet{qq+DiffFc}(obj.PullingLengthIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedRet{qq+DiffFc}(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)                          
                                plot(obj.HHApp{qq+DiffFc}(obj.SnapInIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedApp{qq+DiffFc}(obj.SnapInIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB10,'MarkerEdgeColor',RGB10)                                           
                                plot(obj.HHRet{qq+DiffFc}(obj.AdhForceMaxRetIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedRet{qq+DiffFc}(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)              
                                plot(obj.HHApp{qq+DiffFc}(obj.AdhForceMaxAppIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedApp{qq+DiffFc}(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)                                                                   
                                plot(obj.HHRet{qq+DiffFc}(obj.AdhForceUnbindingIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedRet{qq+DiffFc}(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)                                                                     
                                else
                                area(obj.HHRet{qq}-obj.CP_HardSurface(qq),obj.yRetLim{qq},'FaceColor',RGB12)
                                plot(obj.HHApp{qq}-obj.CP_HardSurface(qq),obj.BasedApp{qq},'Color',RGB1);
                                plot(obj.HHRet{qq}-obj.CP_HardSurface(qq),obj.BasedRetCorr2{qq},'Color',RGB2);
                                plot(obj.HHRet{qq}(obj.PullingLengthIdx(qq))-obj.CP_HardSurface(qq),obj.BasedRetCorr2{qq}(obj.PullingLengthIdx(qq)),'d','MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','b')                          
                                plot(obj.HHRet{qq}(obj.AdhForceMaxRetIdx(qq))-obj.CP_HardSurface(qq),obj.BasedRetCorr2{qq}(obj.AdhForceMaxRetIdx(qq)),'p','MarkerSize',10,'MarkerFaceColor','g','MarkerEdgeColor','g')              
                                plot(obj.HHApp{qq}(obj.AdhForceMaxAppIdx(qq))-obj.CP_HardSurface(qq),obj.BasedApp{qq}(obj.AdhForceMaxAppIdx(qq)),'h','MarkerSize',10,'MarkerFaceColor',RGB7,'MarkerEdgeColor',RGB7)                                                                                            
                                end     
                                % Title for each Subplot
                                ti=title(sprintf('%i',qq),'Color','k');
                                ti.Units='normalized'; % Set units to 'normalized'
                                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                            end
                            if jj == NumFcMax
                                break
                            end
                            jj=jj+1;
                        end
                    else % corresponds to the last figure plotted
                        t = tiledlayout(oo,pp);
                        t.TileSpacing = 'none'; % To reduce the spacing between the tiles
                        t.Padding = 'none'; % To reduce the padding of perimeter of a tile
                        NumFcPlot=oo*pp;
                        if ii==1
                            NumFcCorSelec(ii)=nnz(~obj.SMFSFlag.Uncorrupt(1:NumFcPlot));
                        elseif ii==2
                            NumFcCorSelec(ii)=nnz(obj.SMFSFlag.Uncorrupt((NumFcPlot+1):(NumFcPlot*(ii))));
                        else
                            NumFcCorSelec(ii)=nnz(obj.SMFSFlag.Uncorrupt((NumFcPlot*(ii-1)+1):(NumFcPlot*(ii))));
                        end
                        
                        if ii==1
                            kk=jj;
                        else
                            kk=jj+mm*nn*(ii-1);
                        end
                        %% Plot loop
                        for qq=kk:obj.NCurves % Loop over all force curves in the force map
                            if ww<qq+DiffFc
                                ww=qq+DiffFc;
                            end
                            while ~obj.SMFSFlag.Uncorrupt(ww)     % Stay in the while loop as long as the entry is zero
                                ww=ww+1;
                                if ww>qq
                                    DiffFc=ww-qq;
                                end
                            end
                            if ww>qq
                                ax=nexttile;
                                ax.XLim = [XMin XMax];
                                ax.YLim = [YMin YMax];
                                if obj.SMFSFlag.Fit(qq+DiffFc)==1                            
                                grid on
                                hold on
                                area(obj.HHRet{qq+DiffFc}(1:obj.PullingLengthIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.yRetLim{qq+DiffFc}(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                area(obj.HHApp{qq+DiffFc}(obj.SnapInIdx(qq+DiffFc):end)-obj.CP_HardSurface(qq+DiffFc),obj.yAppLim{qq+DiffFc}(obj.SnapInIdx(qq+DiffFc):end),'FaceColor',RGB7)
                                plot(obj.HHApp{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc),obj.BasedApp{qq+DiffFc},'Color',RGB1);
                                plot(obj.HHRet{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc),obj.BasedRet{qq+DiffFc},'Color',RGB2);
                                plot(obj.HHRet{qq+DiffFc}(obj.PullingLengthIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedRet{qq+DiffFc}(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)                          
                                plot(obj.HHApp{qq+DiffFc}(obj.SnapInIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedApp{qq+DiffFc}(obj.SnapInIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB10,'MarkerEdgeColor',RGB10)                                           
                                plot(obj.HHRet{qq+DiffFc}(obj.AdhForceMaxRetIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedRet{qq+DiffFc}(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)              
                                plot(obj.HHApp{qq+DiffFc}(obj.AdhForceMaxAppIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedApp{qq+DiffFc}(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)                                                                   
                                plot(obj.HHRet{qq+DiffFc}(obj.AdhForceUnbindingIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedRet{qq+DiffFc}(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)                                                                     
                                else
                                area(obj.HHRet{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc),obj.yRetLim{qq+DiffFc},'FaceColor',RGB12)
                                plot(obj.HHApp{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc),obj.BasedApp{qq+DiffFc},'Color',RGB1);
                                plot(obj.HHRet{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc),obj.BasedRetCorr2{qq+DiffFc},'Color',RGB2);
                                plot(obj.HHRet{qq+DiffFc}(obj.PullingLengthIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedRetCorr2{qq+DiffFc}(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','b')                          
                                plot(obj.HHRet{qq+DiffFc}(obj.AdhForceMaxRetIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedRetCorr2{qq+DiffFc}(obj.AdhForceMaxRetIdx(qq+DiffFc)),'p','MarkerSize',10,'MarkerFaceColor','g','MarkerEdgeColor','g')              
                                plot(obj.HHApp{qq+DiffFc}(obj.AdhForceMaxAppIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedApp{qq+DiffFc}(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor','c','MarkerEdgeColor','c')                                                                                            
                                end
                                % Title for each Subplot
                                ti=title(sprintf('%i',qq+DiffFc),'Color','k');
                                %ti=title(sprintf('%i',(kk+ww)/2),'Color','k');
                                ti.Units='normalized'; % Set units to 'normalized'
                                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                            else
                                ax=nexttile;
                                ax.XLim = [XMin XMax];
                                ax.YLim = [YMin YMax];
                                if obj.SMFSFlag.Fit(qq+DiffFc)==1                            
                                grid on
                                hold on
                                area(obj.HHRet{qq+DiffFc}(1:obj.PullingLengthIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.yRetLim{qq+DiffFc}(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                area(obj.HHApp{qq+DiffFc}(obj.SnapInIdx(qq+DiffFc):end)-obj.CP_HardSurface(qq+DiffFc),obj.yAppLim{qq+DiffFc}(obj.SnapInIdx(qq+DiffFc):end),'FaceColor',RGB7)
                                plot(obj.HHApp{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc),obj.BasedApp{qq+DiffFc},'Color',RGB1);
                                plot(obj.HHRet{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc),obj.BasedRet{qq+DiffFc},'Color',RGB2);
                                plot(obj.HHRet{qq+DiffFc}(obj.PullingLengthIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedRet{qq+DiffFc}(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)                          
                                plot(obj.HHApp{qq+DiffFc}(obj.SnapInIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedApp{qq+DiffFc}(obj.SnapInIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB10,'MarkerEdgeColor',RGB10)                                           
                                plot(obj.HHRet{qq+DiffFc}(obj.AdhForceMaxRetIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedRet{qq+DiffFc}(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)              
                                plot(obj.HHApp{qq+DiffFc}(obj.AdhForceMaxAppIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedApp{qq+DiffFc}(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)                                                                   
                                plot(obj.HHRet{qq+DiffFc}(obj.AdhForceUnbindingIdx(qq+DiffFc))-obj.CP_HardSurface(qq+DiffFc),obj.BasedRet{qq+DiffFc}(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)                                                                     
                                else
                                area(obj.HHRet{qq}-obj.CP_HardSurface(qq),obj.yRetLim{qq},'FaceColor',RGB12)
                                plot(obj.HHApp{qq}-obj.CP_HardSurface(qq),obj.BasedApp{qq},'Color',RGB1);
                                plot(obj.HHRet{qq}-obj.CP_HardSurface(qq),obj.BasedRetCorr2{qq},'Color',RGB2);
                                plot(obj.HHRet{qq}(obj.PullingLengthIdx(qq))-obj.CP_HardSurface(qq),obj.BasedRetCorr2{qq}(obj.PullingLengthIdx(qq)),'d','MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','b')                          
                                plot(obj.HHRet{qq}(obj.AdhForceMaxRetIdx(qq))-obj.CP_HardSurface(qq),obj.BasedRetCorr2{qq}(obj.AdhForceMaxRetIdx(qq)),'p','MarkerSize',10,'MarkerFaceColor','g','MarkerEdgeColor','g')              
                                plot(obj.HHApp{qq}(obj.AdhForceMaxAppIdx(qq))-obj.CP_HardSurface(qq),obj.BasedApp{qq}(obj.AdhForceMaxAppIdx(qq)),'h','MarkerSize',10,'MarkerFaceColor','c','MarkerEdgeColor','c')                                                                                            
                                end
                                % Title for each Subplot
                                ti=title(sprintf('%i',qq),'Color','k');
                                %ti=title(sprintf('%i',(kk+ww)/2),'Color','k');
                                ti.Units='normalized'; % Set units to 'normalized'
                                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                            end
                            if jj == RemainderReal
                                break
                            end
                            jj=jj+1;
                        end
                    end
                end
                %% Save figures
                %%% Define the name for the figure title
                partname=sprintf('-p%d',ii);
                % fullname=sprintf('%s%s',figname,partname);
                fullname=sprintf('%s%s',figname,partname);
                %%% Save the current figure in the current folder
                print(gcf,fullname,'-dpng');
            end
            close all
        end
         
                   
        function fc_print_raw(obj,XMin,XMax,YMin,YMax) % fc ... force curve
            % fc_print_raw: A function to simply plot all force curves of a
            % force map without any selection taking place
       
            if nargin < 2
                XMin= -inf;
                XMax= inf;
                YMin= -inf;
                YMax= inf;
            end
            % Define remainder situation
            Remainder=mod(obj.NCurves,25);
            NFigures=floor(obj.NCurves./25);
            if Remainder ~= 0
                NFigures=NFigures+1;
            end 
            % Define variables
            RGB1=[0 26 255]./255;  % Blue 
            RGB2=[255 119 0]./255; % Orange
            RGB10=[69 22 113]./255; % Violet
            % Define variables for the figure name
            ExtendVelocityConvert=num2str(obj.ExtendVelocity*1e9);
            RetractVelocityConvert=num2str(obj.RetractVelocity*1e9);
            % Classification criteria
            figname=strcat(obj.Date,{'_'},obj.Time,{'_'},obj.ID,{'_'},obj.Substrate,{'_'},obj.EnvCond,{'_'},obj.Linker,{'_'},obj.Chipbox,{'_'},obj.ChipCant,{'_'},ExtendVelocityConvert,{'_'},RetractVelocityConvert);
            figname=char(figname);
            %% Figure loop
            for ii=1:NFigures           
            % Figure    
            h_fig=figure(ii);
            h_fig.Color='white'; % changes the background color of the figure
            h_fig.Units='normalized'; % Defines the units 
            h_fig.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
            h_fig.PaperOrientation='landscape';
            h_fig.Name=figname;         
            %% Plotting the tiles
            t = tiledlayout(5,5);
            %t.TileSpacing = 'compact';
            %t.Padding = 'compact';
            t.TileSpacing = 'none'; % To reduce the spacing between the tiles
            t.Padding = 'none'; % To reduce the padding of perimeter of a tile
            
            % Defining variables
            if ii==NFigures && Remainder~=0
                NLoop=Remainder;
            else
                NLoop=25;
            end
                %% Plot loop    
                for jj=1:NLoop
                    % Tile jj
                    kk=jj+25*(ii-1);
                    %%% Define some variables
                    x50=-50e-9; % Defines 50nm
                    x150=-150e-9; % Defines 150nm
                    x500=-500e-9; % Defines 500nm
                    % Plot tile
                    ax=nexttile;      
                    ax.XLim = [XMin XMax];
                    ax.YLim = [YMin YMax];
                    hold on
                    grid on
                    %plot(obj.HHApp{kk}-obj.CP_HardSurface(kk,1),obj.BasedApp{kk},'Color',RGB1);
                    plot(obj.HHApp{kk},obj.App{kk},'Color',RGB1);
                    %plot(obj.HHRet{kk}-obj.CP_HardSurface(kk,1),obj.BasedRetCorr{kk});
                    %plot(obj.HHRet{kk}-obj.CP_HardSurface(kk,1),obj.BasedRet{kk},'Color',RGB2); 
                    plot(obj.HHRet{kk},obj.Ret{kk},'Color',RGB2); 
                    %xline(x50,'Color',RGB10); % Draws a vertical line
                    %xline(x150,'Color',RGB10); % Draws a vertical line   
                    %xline(x500,'Color',RGB10); % Draws a vertical line
                    % Title for each Subplot
                    ti=title(sprintf('%i',kk),'Color','k');                                     
                    ti.Units='normalized'; % Set units to 'normalized'  
                    ti.Position=[0.5,0.95]; % Position the subplot title within the subplot                   
                    % Legend, x- and y-labels
                    %legend('Approach','Retraction','Location','best')
                    %xlabel('Tip-sample seperation  (nm)','FontSize',11,'Interpreter','latex');
                    %ylabel('Force (nN)','FontSize',11,'Interpreter','latex');                  
                end

            %% Save figures
            %%% Define the name for the figure title    
            partname=sprintf('-p%d',ii);        
            % fullname=sprintf('%s%s',figname,partname);
            fullname=sprintf('%s%s',figname,partname);
            %%% Save the current figure in the current folder
            print(gcf,fullname,'-dpng'); 
            end
        close Figure 1 Figure 2 Figure 3 Figure 4
        end
    
        %%
        %%
        function fc_print_coloured_section(obj,XMin,XMax,YMin,YMax) % fc ... force curve
            % fc_print_raw: A function to simply plot all force curves of a
            % force map without any selection taking place
       
            if nargin < 2
                XMin= -inf;
                XMax= inf;
                YMin= -inf;
                YMax= inf;
            end
            % Define remainder situation
            Remainder=mod(obj.NCurves,25);
            NFigures=floor(obj.NCurves./25);
            if Remainder ~= 0
                NFigures=NFigures+1;
            end 
                      
            % Define variables for the figure name
            ExtendVelocityConvert=num2str(obj.ExtendVelocity*1e9);
            RetractVelocityConvert=num2str(obj.RetractVelocity*1e9);
            % Classification criteria
            figname=strcat(obj.Date,{'_'},obj.Time,{'_'},obj.ID,{'_'},obj.Substrate,{'_'},obj.EnvCond,{'_'},obj.Linker,{'_'},obj.Chipbox,{'_'},obj.ChipCant,{'_'},ExtendVelocityConvert,{'_'},RetractVelocityConvert);
            figname=char(figname);
            %% Figure loop
            for ii=1:NFigures           
            % Figure    
            h_fig=figure(ii);
            h_fig.Color='white'; % changes the background color of the figure
            h_fig.Units='normalized'; % Defines the units 
            h_fig.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
            h_fig.PaperOrientation='landscape';
            h_fig.Name=figname;         
            %% Plotting the tiles
            t = tiledlayout(5,5);
            %t.TileSpacing = 'compact';
            %t.Padding = 'compact';
            t.TileSpacing = 'none'; % To reduce the spacing between the tiles
            t.Padding = 'none'; % To reduce the padding of perimeter of a tile
            
            % Defining variables
            if ii==NFigures && Remainder~=0
                NLoop=Remainder;
            else
                NLoop=25;
            end
                %% Plot loop    
                for jj=1:NLoop
                    % Tile jj
                    kk=jj+25*(ii-1);
                    %%% Define some variables
                    x50=-50e-9; % Defines 50nm
                    x150=-150e-9; % Defines 150nm
                    x500=-500e-9; % Defines 500nm
                    % Plot tile
                    ax=nexttile;      
                    ax.XLim = [XMin XMax];
                    ax.YLim = [YMin YMax];
                    hold on
                    grid on
                    plot(obj.THApp{kk}-obj.CP_HardSurface(kk,1),obj.BasedApp{kk});
                    %plot(obj.THRet{kk}-obj.CP_HardSurface(kk,1),obj.BasedRetCorr{kk});
                    plot(obj.THRet{kk}-obj.CP_HardSurface(kk,1),obj.BasedRet{kk});                    
                    xline(x50,'Color','r'); % Draws a vertical line
                    xline(x150,'Color','r'); % Draws a vertical line   
                    xline(x500,'Color','r'); % Draws a vertical line
                    % Title for each Subplot
                    ti=title(sprintf('%i',kk),'Color','k');                                     
                    ti.Units='normalized'; % Set units to 'normalized'  
                    ti.Position=[0.5,0.95]; % Position the subplot title within the subplot                   
                    % Legend, x- and y-labels
                    %legend('Approach','Retraction','Location','best')
                    %xlabel('Tip-sample seperation  (nm)','FontSize',11,'Interpreter','latex');
                    %ylabel('Force (nN)','FontSize',11,'Interpreter','latex');                  
                end

            %% Save figures
            %%% Define the name for the figure title    
            partname=sprintf('-p%d',ii);        
            % fullname=sprintf('%s%s',figname,partname);
            fullname=sprintf('%s%s',figname,partname);
            %%% Save the current figure in the current folder
            print(gcf,fullname,'-dpng'); 
            end
        close Figure 1 Figure 2 Figure 3 Figure 4
        end
          
       
        function fc_fine_figure(obj,XMin,XMax,YMin,YMax,ii)
             
            if nargin < 2
                XMin= -inf;
                XMax= inf;
                YMin= -inf;
                YMax= inf;
            end
            
            % Define RGB colours
            RGB1=[0 26 255]./255;  % Blue 
            RGB2=[255 119 0]./255; % Orange
            RGB3=[80 220 100]./255; % Emerald
            RGB10=[205 207 208]./255; % Grey
              
            [Xmultiplier,Xunit,~] = AFMImage.parse_unit_scale(1e+9,'nm',1);
            [Ymultiplier,Yunit,~] = AFMImage.parse_unit_scale(1e+9,'nN',1);
            
            % Define some variables
            kk=21;
            VelocityConvert=num2str(obj.Velocity/Xmultiplier); % Convert into nm
            figname=strcat(obj.Date,{'_'},obj.Time,{'_'},obj.ID,{'_'},obj.Substrate,{'_'},obj.EnvCond,{'_'},VelocityConvert,{'_'},obj.Chipbox,{'_'},obj.ChipCant);
            figname=char(figname);
          
            % Allocate data
           xApp1=(obj.THApp{kk}-obj.CP_HardSurface(kk))/-Xmultiplier;
          %  xApp1b=(obj.THApp{kk})*-1e9;
            yApp1=(obj.BasedApp{kk})/Ymultiplier;
            xRet1=(obj.THRet{kk}-obj.CP_HardSurface(kk))/-Xmultiplier;
            xRet1b=(obj.THRet{kk})/-Xmultiplier;
            yRet1=(obj.BasedRetCorr2{kk})/Ymultiplier;
            
            yRet_2=(obj.yRetLim2{kk})*1e9;
            yRet_1=(obj.yRetLim{kk})*1e9;
        
            % Figure
            h_fig=figure(ii);
            h_fig.Color='white'; % changes the background color of the figure
            h_fig.Units='normalized'; % Defines the units
            h_fig.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
            h_fig.PaperOrientation='landscape';
            h_fig.Name=figname;
            % Plot
            hold on
            grid on
           % area(xRet1,yRet_1,'FaceColor',RGB10);
            area(xRet1,yRet1,'FaceColor',RGB10);
            plot(xApp1,yApp1,'Color',RGB1,'LineWidth',6);
            plot(xRet1,yRet1,'Color',RGB2,'LineWidth',6);
          %  plot(xApp1,yApp1,'Color',RGB1,'LineWidth',4);
        %    plot(xRet1(obj.PullingLengthIdx:end),yRet1(obj.PullingLengthIdx:end),'Color',RGB3,'LineWidth',4);
            % Title for each Subplot
           % ti=title(sprintf('Typical force-distance curve',kk),'Color','k');
           % ti.Units='normalized'; % Set units to 'normalized'
           % ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
            % Legend
            le=legend(' Adhesion energy',' Approach data',' Retraction data','Location','best');
            le.FontSize = 48;      
            le.EdgeColor='w';
            %le.Box = 'off';
            %%% Axes
            ax = gca; % current axes
            ax.FontSize = 48;
            ax.LineWidth = 5;
            ax.XTick=0:100:400;
           % ax.XTickLabel=[];
            ax.YTick=-0.3:0.1:0.2;
           % ax.YTickLabel=[];
            %ax.XLabel.String = sprintf('Tip-Substrate separation (%s)',Xunit);
            ax.XLabel.String = 'Tip-Substrate separation (nm)';
            ax.XLabel.FontSize = 52;
            %ax.YLabel.String = sprintf('Force (%s)',Yunit);
            ax.YLabel.String = 'Force (nN)';
            ax.YLabel.FontSize = 52;
            ax.XLim = [XMin XMax];
            ax.YLim = [YMin YMax];
            %% Save figures
            %%% Define the name for the figure title
            partname=sprintf('-ForceCurve%d',kk);
            % fullname=sprintf('%s%s',figname,partname);
            fullname=sprintf('%s%s',figname,partname);
            %%% Save the current figure in the current folder
            print(gcf,fullname,'-dpng');
        end
                
        function test(obj) 
            % TEST FUNCTION
            mega=10e6;
                giga=10e9;
                tera=10e12;
            jj=5;
           obj.xDataToFit=obj.HHApp{jj}*mega
           obj.yDataToFit=obj.App{jj}*giga
           

            obj.yDataToFit(ceil(end*0.9):end)=[]
            obj.xDataToFit(ceil(end*0.9):end)=[]
          
%               xApp=obj.HHApp{ii}-obj.CP_HardSurface(ii);
%                 xRet=obj.HHRet{ii}-obj.CP_HardSurface(ii);
              
            
        end

      
    end    
    methods (Static)
        % Auxiliary methods
        
        function [nocontvDef,domainidx] = no_contact_domain(vDef)
            % finds the index of a point in the curve well before contact point and
            % returns a vector with the curve up to that point and the assotiated
            % index.
            
            diffvDef = diff(smoothdata(vDef));
            leng = length(diffvDef);
            stddiff = std(diffvDef(1:floor(0.5*leng)));
            domainidx = 1;
            while diffvDef(domainidx) < 6*stddiff
                domainidx = domainidx + 1;
                if domainidx == leng
                    domainidx = floor(0.6*leng);
                    break
                end
            end
            domainidx = domainidx - floor(0.05*leng);
            if domainidx < floor(1/3*length(vDef))
                domainidx = floor(1/3*length(vDef));
            end
            nocontvDef = vDef(1:domainidx);
        end
        
        function [OutForce,DidProcess] = remove_sinusoidal_approach(InForce,PowerCutoff)
            
            if nargin < 2
                PowerCutoff = .2;
            end
            
            LossBefore = ForceMap.low_frequency_loss(InForce);
            
            FFT = fft(InForce);
            MidIndex = floor(length(FFT)/2+1);
            FFT([1 end]) = 0; % Removes constant component of signal
            PowerSpectrum = abs(FFT).^2/length(FFT)/range(InForce)^2;
            
%             % find peaks in the spectrum
%             [Peaks,PeakPos,PeakWidths] = findpeaks(PowerSpectrum(1:MidIndex),...
%                 'MinPeakDistance',10,'MinPeakHeight',PowerCutoff,'SortStr','descend');
%             if isempty(Peaks) || PeakPos(1) < 3
%                 OutForce = InForce;
%                 DidProcess = false;
%                 return
%             end
            BandWidth = 400;
            
            FiltFFT = FFT;
            FiltFFT([1:BandWidth (end-(BandWidth)-1:end)]) = 0;
            OutForce = real(ifft((FiltFFT)));
            DidProcess = true;
            
            LossAfter = ForceMap.low_frequency_loss(OutForce);
            
            % Debug Section
            subplot(4,1,1)
            plot(InForce)
            subplot(4,1,2)
            plot(PowerSpectrum)
            subplot(4,1,3)
            plot(1:10000,sin([1:10000]./8))
            subplot(4,1,4)
            plot(real(OutForce))
            drawnow
            
        end
        
        function [Loss] = low_frequency_loss(InSignal)
            Loss = sum(abs(InSignal));
        end
        
        function [n,V,p] = affine_fit(X)
            %Computes the plane that fits best (lest square of the normal distance
            %to the plane) a set of sample points.
            %INPUTS:
            %
            %X: a N by 3 matrix where each line is a sample point
            %
            %OUTPUTS:
            %
            %n : a unit (column) vector normal to the plane
            %V : a 3 by 2 matrix. The columns of V form an orthonormal basis of the
            %plane
            %p : a point belonging to the plane
            %
            %NB: this code actually works in any dimension (2,3,4,...)
            %Author: Adrien Leygue
            %Date: August 30 2013
            
            %the mean of the samples belongs to the plane
            p = mean(X,1);
            
            %The samples are reduced:
            R = bsxfun(@minus,X,p);
            %Computation of the principal directions if the samples cloud
            [V,~] = eig(R'*R);
            %Extract the output from the eigenvectors
            n = V(:,1);
            V = V(:,2:end);
        end
        
        function [HeadHeight ,Force,spring_constant,sensitivity]=...
                writedata(HeaderFileDirectory,SegmentHeaderFileDirectory,...
                HeighDataDirectory,vDelfDataDirectory,HHType)
            % Author: Orestis Andriotis (changed and adapted by Manuel Rufin)
            % HeaderFileDirectory: file directory of the main header properties. There
            % is the information about the scaling factors to convert the raw data.
            % SegmentHeaderFileDirectory: file directory of the header of each segment.
            % Each segment means the loading (#0) and unloading (#1).
            
            % HeighDataDirectory: file directory of capacitiveSensorHeight.dat file
            % vDelfDataDirectory: file directory of vDeflection.dat file
            %
            
            % HeaderFileDirectory = headerDir;
            % SegmentHeaderFileDirectory=segHeaderDir;
            % HeighDataDirectory=heightdir;
            % vDelfDataDirectory=vDefldir;
            % find the number of total points
            A=fileread(SegmentHeaderFileDirectory);
            FileID=fopen(SegmentHeaderFileDirectory,'rt','n','UTF-8');
            % FileID = fopen(filename,permission,machinefmt,encodingIn)
            
            B=strfind(A,'force-segment-header.num-points=');
            % strfind(file,string) is looking for a specific string in the file.
            fseek(FileID,B,'cof');
            % moves at the location where specific string is located
            
            tline = fgetl(FileID);
            where=strfind(tline,'=');
            % "where" is the position of the "=" symbol in the tline string
            n = str2double(tline(where+1:end)); % number of points
            fclose(FileID);
            
            % Read the height measured data and the vertical deflection.
            % Reading the raw height data into a column of length n
            FileID = fopen(HeighDataDirectory);
            % fread(FileID,sizeA,precision,skip,machinefmt)
            RawHeight = fread(FileID,n,'int32',0,'b'); %raw data
            fclose(FileID);
            
            
            % Reading the deflection height data into a column of length n
            FileID = fopen(vDelfDataDirectory);
            RawvDeflection = fread(FileID,n,'int32',0,'b'); %raw data
            
            
            % get scaling factors from header file
            % calls getheaderinfo function to extract the scaling factors.
            % All the scaling factors are included in the header.properties file and
            % not the segment-header.properties.
            [mult_height_meters1, offset_height_meters1,...
                mult_height_meters2, offset_height_meters2,...
                mult_vDefl_volts, offset_vDefl_volts,...
                sensitivity, spring_constant] =...
                ForceMap.getheaderinfo(HeaderFileDirectory,HHType);
            
            % converts raw into volts by calling the getheaderinfo function
            heightvolts = RawHeight.*mult_height_meters1 + offset_height_meters1;
            
            deflvolts=RawvDeflection.*mult_vDefl_volts + offset_vDefl_volts;
            
            % converts volts into meters
            height_meters = heightvolts.*mult_height_meters2 + offset_height_meters2;
            % height_meters = heightvolts.*mult_height_meters + offset_height_meters;
            deflmeter=deflvolts.*sensitivity;
            
            HeadHeight = height_meters;
            
            Force = deflmeter;
            
            fclose(FileID);
        end
        
        function [mult_height_meters1, offset_height_meters1,...
                mult_height_meters2, offset_height_meters2,...
                mult_vDefl_volts, offset_vDefl_volts,...
                sensitivity, spring_constant] = getheaderinfo(filedirectory,HHType)
            % Author: Orestis Andriotis
            % getheaderinfo
            % reads header from a force scan series and extracts the scaling
            % parameters.
            % filedirectory=HeaderFileDirectory;
            FileID=fopen(filedirectory,'rt','n','UTF-8'); % FileID = fopen(filename,permission,machinefmt,encodingIn)
            A=fileread(filedirectory);
            % Height: 1. CONVERSION raw-meters & 2. SCALE meters
            % Conversion RAW -> VOLTS
            fseek(FileID,1,'cof'); % goes at the first position in the file
            
            exp1='\d{1}';
            
            frewind(FileID);
            B=strfind(A,HHType);
            tline = A(B:end);
            HHNum = regexp(tline,exp1,'match','once');
            
            clear tline;
            frewind(FileID);
            B=regexp(A,'lcd-info\.\d+\.channel\.name=vDeflection');
            tline = A(B:end);
            vDefNum = regexp(tline,exp1,'match','once');
            
            %   Multiplier
            clear tline;
            frewind(FileID);
            B=strfind(A,strcat('lcd-info.',HHNum,'.encoder.scaling.multiplier='));
            % strfind(file,string) is looking for a specific string in the file.
            fseek(FileID,B,'cof');
            % moves at the location where specific string is located
            tline = fgetl(FileID);
            % stores that string in a character
            where=strfind(tline,'=');
            % "where" is the position of the "=" symbol in the tline string
            mult_height_meters1 = str2double(... % convert the string to number
                tline(where+1:end)... % this is the number
                );
            
            %   Offset
            clear tline where;
            frewind(FileID);
            B=strfind(A,strcat('lcd-info.',HHNum,'.encoder.scaling.offset='));
            fseek(FileID,B,'cof');
            tline = fgetl(FileID);
            where=strfind(tline,'=');
            offset_height_meters1 = str2double(tline(where+1:end));
            
            %   Multiplier
            clear tline where;
            frewind(FileID);
            B=strfind(A,strcat('lcd-info.',HHNum,'.conversion-set.conversion.nominal.scaling.multiplier='));
            fseek(FileID,B,'cof');
            tline = fgetl(FileID);
            where=strfind(tline,'=');
            mult_height_meters2 = str2double(tline(where+1:end));
            
            
            %   Offset
            clear tline where;
            frewind(FileID);
            B=strfind(A,strcat('lcd-info.',HHNum,'.conversion-set.conversion.nominal.scaling.offset='));
            fseek(FileID,B,'cof');
            tline = fgetl(FileID);
            where=strfind(tline,'=');
            offset_height_meters2 = str2double(tline(where+1:end));
            
            
            % vDeflection: 1. CONVERSION raw-volts & 2. volts to meters
            % Conversion RAW -> VOLTS
            
            %   Multiplier
            clear tline where;
            frewind(FileID);
            B=strfind(A,strcat('lcd-info.',vDefNum,'.encoder.scaling.multiplier='));
            fseek(FileID,B,'cof');
            tline = fgetl(FileID);
            where=strfind(tline,'=');
            mult_vDefl_volts = str2double(tline(where+1:end)); % multiplier for scaling the raw height data and convert to volts
            
            %   Offset
            clear tline where;
            frewind(FileID);
            B=strfind(A,strcat('lcd-info.',vDefNum,'.encoder.scaling.offset='));
            fseek(FileID,B,'cof');
            tline = fgetl(FileID);
            where=strfind(tline,'=');
            offset_vDefl_volts = str2double(tline(where+1:end));
            
            
            % Conversion VOLTS -> METERS
            
            %   Multiplier (that is the sensitivity measured in meters per Volts)
            clear tline where;
            frewind(FileID);
            B=strfind(A,strcat('lcd-info.',vDefNum,'.conversion-set.conversion.distance.scaling.multiplier='));
            fseek(FileID,B,'cof');
            tline = fgetl(FileID);
            where=strfind(tline,'=');
            sensitivity = str2double(tline(where+1:end));
            
            % Spring constant
            
            clear tline where;
            frewind(FileID);
            B=strfind(A,strcat('lcd-info.',vDefNum,'.conversion-set.conversion.force.scaling.multiplier='));
            fseek(FileID,B,'cof');
            tline = fgetl(FileID);
            where=strfind(tline,'=');
            spring_constant = str2double(tline(where+1:end));
            
            fclose(FileID);
        end
        
        function X = CP_batchprep_3_channel(objcell,ImgSizeFinal,ImgSize,CutPercent)
            % ImgSizeFinal ...Size of the Output (X(:,:,-,-))
            % ImgSize ...Size of the Image the final Image is resized from
            % CutPercent ...1x3 Vector with the percentage that is to be
            % cut off from the left end of the curve image; one percentage
            % for each channel
            
            if nargin<2
                ImgSize = 478;
                ImgSizeFinal = 128;
                CutPercent = [0 0.3 0.7];
            elseif nargin < 3
                ImgSize = 478;
                CutPercent = [0 0.3 0.7];
            elseif nargin<4
                CutPercent = [0 0.3 0.7];
            end
            Nmaps = length(objcell);
            
            % Preallocate sizes of the outputvariables. The trainNetwork() function
            % requires the batch X of predictor varables (in this case the force-images)
            % to be a height-by-width-by-channelnumber-by-Numberofimages array and the
            % prelabeled regression responses Y to be a
            % Numberofimages-by-Numberofresponses array
            Nimgs = 0;
            for i=1:Nmaps
                Nimgs = Nimgs + sum(objcell{i}.SelectedCurves);
            end
            if ImgSize == ImgSizeFinal
                X = zeros(ImgSizeFinal,ImgSizeFinal,3,Nimgs,'logical');
            else
                X = zeros(ImgSizeFinal,ImgSizeFinal,3,Nimgs,'single');
            end
            
            
            k = 1;
            for i=1:Nmaps
                jRange = find(objcell{i}.SelectedCurves);
                for j=jRange'
                    o = 1;
                    for n=CutPercent
                        if ImgSize == ImgSizeFinal
                            Image = zeros(ImgSizeFinal,ImgSizeFinal,'logical');
                        else
                            Image = zeros(ImgSizeFinal,ImgSizeFinal,'single');
                        end
                        
                        [App,HHApp] = objcell{i}.get_force_curve_data(j,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                        Points(:,1) = HHApp(floor(n*length(HHApp))+1:end);
                        Points(:,2) = App(floor(n*length(App))+1:end);
                        Points(:,1) = (Points(:,1)-min(Points(:,1)))/range(Points(:,1))*(ImgSize-1);
                        Points(:,2) = (Points(:,2)-min(Points(:,2)))/range(Points(:,2))*(ImgSize-1);
                        
                        L = length(Points);
                        for m=1:L
                            
                            if m==1
                                Position = [floor(Points(m,1))+1,floor(Points(m,2))+1];
                            end
                            if m<L
                                NextPosition = [floor(Points(m+1,1))+1,floor(Points(m+1,2))+1];
                            end
                            Image((ImgSize+1)-Position(2),Position(1)) = 1;
                            
                            % fill out points between actual data points
                            if m<L
                                L1Norm = norm(Points(m+1,:)-Points(m,:));
                                IncrX = (Points(m+1,1) - Points(m,1))/L1Norm;
                                IncrY = (Points(m+1,2) - Points(m,2))/L1Norm;
                                FillerPos = Points(m,:);
                                while norm((FillerPos + [IncrX IncrY]) - Points(m,:)) < norm(Points(m+1,:) - Points(m,:))
                                    FillerPos(1) = FillerPos(1) + IncrX;
                                    FillerPos(2) = FillerPos(2) + IncrY;
                                    Image((ImgSize+1)-(floor(FillerPos(2))+1),floor(FillerPos(1))+1) = 1;
                                end
                            end
                            Position = NextPosition;
                        end
                        if ImgSize ~= ImgSizeFinal
                            Image = imresize(Image,[ImgSizeFinal ImgSizeFinal]);
                        end
                        X(:,:,o,k) = Image;
                        clear Points
                        o = o + 1;
                    end
                    k = k + 1;
                end
            end
        end
        
        function [E_mod,GoF,Hertzfit] = hertz_fit_gof(tip_h,force,CP,curve_percent,shape,TipRadius,PoissonR)
            
            if TipRadius == -1
                prompt = {'What is the nominal tip radius of the used cantilever in m?'};
                dlgtitle = 'Cantilever tip';
                dims = [1 35];
                definput = {'10'};
                TipRadius = str2double(inputdlg(prompt,dlgtitle,dims,definput));
            end
            
            if nargin < 6
                shape = 'parabolic';
            end
            % define fitrange as the lower 'curve_percent' of the curve and
            % normalize the data for optimal fitting.
            fitrange = [CP:(CP+floor(curve_percent*(length(tip_h)-CP)))];
            CPforce = force(fitrange) - force(CP);
            CPtip_h = tip_h(fitrange) - tip_h(CP);
            ranf = range(CPforce);
            rant = range(CPtip_h);
            CPforce = CPforce/ranf;
            CPtip_h = CPtip_h/rant;
            
            if isequal(shape,'parabolic')
                s = fitoptions('Method','NonlinearLeastSquares',...
                    'Lower',10^(-30),...
                    'Upper',inf,...
                    'Startpoint',1);
                f = fittype('a*(x)^(3/2)','options',s);
                [Hertzfit,GoF] = fit(CPtip_h,...
                    CPforce,f);
                % calculate E module based on the Hertz model. Be careful
                % to convert to unnormalized data again
                E_mod = 3*(Hertzfit.a*ranf/rant^(3/2))/(4*sqrt(TipRadius))*(1-PoissonR^2);
            elseif isequal(shape,'spherical')
            elseif isequal(shape,'conical')
            elseif isequal(shape,'pyramid')
            end
        end
        
    end
    
    methods
        % non-static auxiliary methods
        
        function calculate_reference_slope_from_area(obj,Mask,AppRetSwitch)
            % Calculates the distribution of DZslopes on the curves
            % that that are neihter on the fibril nor the excluded zones.
            %  the upper 25% of the curve are considered for the
            %  calculation
            
            Range = find(obj.SelectedCurves & ~obj.CorruptedCurves);
            CurvePercent = 0.25;
            % Calculate the DZslopes
            k = 1;
            for i=Range'
                if (Mask(obj.List2Map(i,1),obj.List2Map(i,2)) == 1) &&...
                        (obj.ExclMask(obj.List2Map(i,1),obj.List2Map(i,2)) == 1)
                    [Force,Height] = obj.get_force_curve_data(i,'AppRetSwitch',AppRetSwitch,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                    Z = Height - obj.CP(i,1);
                    D = (Force - obj.CP(i,2))/obj.SpringConstant;
                    Dmax = max(D);
                    DCurvePercent = D(D>=(1-CurvePercent)*Dmax);
                    if length(DCurvePercent) < 2
                        DCurvePercent = D(end-1:end);
                    end
                    ZCurvePercent = Z(1:length(DCurvePercent));    
                    LineFit = polyfit(ZCurvePercent,DCurvePercent,1);
                    DZslope(k) = LineFit(1);
%                     % Debug
%                     plot(Z,D,Z,polyval(LineFit,Z))
                    k = k + 1;
                end
            end
            % Fit the log-Gaussian
            DZslopeCleaned = DZslope;
            DZslopeCleaned(DZslopeCleaned <= 0) = [];
            DZslopeCleaned = exp(rmoutliers(log(DZslopeCleaned),'quartiles'));
            Gaussian = fitdist(DZslopeCleaned','Lognormal');
            % Expected value of a lognormal is exp(mu + sigma^2/2)
            obj.RefSlope = exp(Gaussian.mu+Gaussian.sigma^2/2);
            obj.RefSlopeCorrectedSensitivity = obj.Sensitivity/obj.RefSlope;
            obj.HasRefSlope = true;
        end
        
        function calculate_adaptive_sensitivity_from_area(obj,varargin)
            % function calculate_adaptive_sensitivity_from_area(obj,varargin)
            %
            % <FUNCTION DESCRIPTION HERE>
            %
            %
            % Required inputs
            % obj ... <VARIABLE DESCRIPTION>
            %
            % Name-Value pairs
            % "Mask" ... <NAMEVALUE DESCRIPTION>
            % "AppRetSwitch" ... <NAMEVALUE DESCRIPTION>
            % "MovingWindowSize" ... <NAMEVALUE DESCRIPTION>
            % "FitRangeMode" ... <NAMEVALUE DESCRIPTION>
            % "FitRangeLowerFraction" ... <NAMEVALUE DESCRIPTION>
            % "FitRangeUpperFraction" ... <NAMEVALUE DESCRIPTION>
            % "FitRangeLowerValue" ... <NAMEVALUE DESCRIPTION>
            % "FitRangeUpperValue" ... <NAMEVALUE DESCRIPTION>
            % "Verbose" ... <NAMEVALUE DESCRIPTION>
            
            p = inputParser;
            p.FunctionName = "calculate_adaptive_sensitivity_from_area";
            p.CaseSensitive = false;
            p.PartialMatching = true;
            
            % Required inputs
            validobj = @(x)true;
            addRequired(p,"obj",validobj);
            
            % NameValue inputs
            defaultMask = ones(size(obj.Channel(1).Image));
            defaultAppRetSwitch = 0;
            defaultMovingWindowSize = obj.NumPixelsX;
            defaultFitRangeMode = 'ValueHorizontal';
            defaultFitRangeLowerFraction = 0.75;
            defaultFitRangeUpperFraction = 1;
            defaultFitRangeLowerValue = 0;
            defaultFitRangeUpperValue = 4e-9;
            defaultVerbose = false;
            validMask = @(x)isequal(size(x),size(obj.Channel(1).Image));
            validAppRetSwitch = @(x)x==0|x==1|islogical(x);
            validMovingWindowSize = @(x)x<=obj.NCurves;
            validFitRangeMode = @(x)any(validatestring(x,{'ValueHorizontal','FractionHorizontal','ValueVertical','FractionVertical'}));
            validFitRangeLowerFraction = @(x)(0<=x)&&(x<=1);
            validFitRangeUpperFraction = @(x)(0<=x)&&(x<=1);
            validFitRangeLowerValue = @(x)isscalar(x);
            validFitRangeUpperValue = @(x)isscalar(x);
            validVerbose = @(x)islogical(x);
            addParameter(p,"Mask",defaultMask,validMask);
            addParameter(p,"AppRetSwitch",defaultAppRetSwitch,validAppRetSwitch);
            addParameter(p,"MovingWindowSize",defaultMovingWindowSize,validMovingWindowSize);
            addParameter(p,"FitRangeMode",defaultFitRangeMode,validFitRangeMode);
            addParameter(p,"FitRangeLowerFraction",defaultFitRangeLowerFraction,validFitRangeLowerFraction);
            addParameter(p,"FitRangeUpperFraction",defaultFitRangeUpperFraction,validFitRangeUpperFraction);
            addParameter(p,"FitRangeLowerValue",defaultFitRangeLowerValue,validFitRangeLowerValue);
            addParameter(p,"FitRangeUpperValue",defaultFitRangeUpperValue,validFitRangeUpperValue);
            addParameter(p,"Verbose",defaultVerbose,validVerbose);
            
            parse(p,obj,varargin{:});
            
            % Assign parsing results to named variables
            obj = p.Results.obj;
            Mask = p.Results.Mask;
            AppRetSwitch = p.Results.AppRetSwitch;
            MovingWindowSize = p.Results.MovingWindowSize;
            FitRangeMode = p.Results.FitRangeMode;
            FitRangeLowerFraction = p.Results.FitRangeLowerFraction;
            FitRangeUpperFraction = p.Results.FitRangeUpperFraction;
            FitRangeLowerValue = p.Results.FitRangeLowerValue;
            FitRangeUpperValue = p.Results.FitRangeUpperValue;
            Verbose = p.Results.Verbose;
            
            
            Range = find(obj.SelectedCurves & ~obj.CorruptedCurves);
            isCalibrationCurve = false(obj.NCurves,1);
            DZSlopes = ones(obj.NCurves,1);
            
            % first, calculate all the DZ-Slopes
            for i=Range'
                [vDef,Height] = obj.get_force_curve_data(i,'AppRetSwitch',AppRetSwitch,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','m');
                switch FitRangeMode
                    case 'ValueHorizontal'
                        UValue = max(Height) - FitRangeLowerValue;
                        LValue = max(Height) - FitRangeUpperValue;
                        FitHeight = Height(Height>=LValue & Height<=UValue);
                        FitvDef = vDef(Height>=LValue & Height<=UValue);
                    case 'ValueVertical'
                        FitHeight = Height(vDef>=FitRangeLowerValue & vDef<=FitRangeUpperValue);
                        FitvDef = vDef(vDef>=FitRangeLowerValue & vDef<=FitRangeUpperValue);
                    case 'FractionHorizontal'
                        RangeHeight = range(Height);
                        FitHeight = Height((Height-min(Height))>=FitRangeLowerFraction*RangeHeight &...
                            (Height-min(Height))<=FitRangeUpperFraction*RangeHeight);
                        FitvDef = vDef((Height-min(Height))>=FitRangeLowerFraction*RangeHeight &...
                            (Height-min(Height))<=FitRangeUpperFraction*RangeHeight);
                    case 'FractionVertical'
                        MaxvDef = max(vDef);
                        FitHeight = Height(vDef>=FitRangeLowerFraction*MaxvDef & vDef<=FitRangeUpperFraction*MaxvDef);
                        FitvDef = vDef(vDef>=FitRangeLowerFraction*MaxvDef & vDef<=FitRangeUpperFraction*MaxvDef);
                end
                Params(i,:) = polyfit(FitHeight,FitvDef,1);
                DZSlopes(i) = Params(i,1);
                
                if (Mask(obj.List2Map(i,1),obj.List2Map(i,2)) == 1) &&...
                        (obj.ExclMask(obj.List2Map(i,1),obj.List2Map(i,2)) == 1)
                    isCalibrationCurve(i) = true;
                end
            end
            CalibrationCurveIndizes = find(isCalibrationCurve);
            
            if Verbose 
                BackgroundImage = obj.convert_data_list_to_map(isCalibrationCurve).*.2;
                F = figure('Color','w');
            end
            
            if MovingWindowSize == obj.NCurves
            % Fit the log-Gaussian
            DZslopeCleaned = DZSlopes(CalibrationCurveIndizes);
            DZslopeCleaned(DZslopeCleaned <= 0) = [];
            DZslopeCleaned = exp(rmoutliers(log(DZslopeCleaned),'quartiles'));
            Gaussian = fitdist(DZslopeCleaned,'Lognormal');
            % Expected value of a lognormal is exp(mu + sigma^2/2)
            obj.RefSlope = exp(Gaussian.mu+Gaussian.sigma^2/2);
            obj.RefSlopeCorrectedSensitivity = obj.Sensitivity/obj.RefSlope;
            obj.HasRefSlope = true;
            return
            end
            
            % Calculate median slope inside sliding window around current
            % point. The sliding window does not shrink below
            % MovingWindowSize/2 around the boundaries.
            LowerWindow = floor((MovingWindowSize-1)/2)+1;
            UpperWindow = ceil((MovingWindowSize-1)/2);
            for i=Range'
                TempLowerIndex = CalibrationCurveIndizes(CalibrationCurveIndizes <= i);
                TempUpperIndex = CalibrationCurveIndizes(CalibrationCurveIndizes > i);
                DiffLow = max(LowerWindow - length(TempLowerIndex),0);
                DiffUp = max(UpperWindow - length(TempUpperIndex),0);
                LowerIndizes = TempLowerIndex(max(1,end+1-(DiffUp + LowerWindow - DiffLow)):end);
                UpperIndizes = TempUpperIndex(1:min(end,(DiffLow + UpperWindow - DiffUp)));
                WindowIndizes = [LowerIndizes' UpperIndizes'];
                WindowSlopes = DZSlopes(WindowIndizes);
                AdaptiveMedianRefSlope(i) = median(WindowSlopes);
                if Verbose
                    CurrentWindowList = zeros(obj.NCurves,1);
                    CurrentWindowList(WindowIndizes) = .8;
                    CurrentWindowImage = obj.convert_data_list_to_map(CurrentWindowList);
                    CurrentCombinedImage = BackgroundImage + CurrentWindowImage;
                    gcf = F;
                    subplot(2,2,1)
                    imshow(CurrentCombinedImage,[]);
                    subplot(2,2,2)
                    histogram(WindowSlopes)
                    hold on
                    plot(AdaptiveMedianRefSlope(i),MovingWindowSize)
                    hold off
                    subplot(2,2,3)
                    title('Adaptive Median Slope')
                    plot(AdaptiveMedianRefSlope(1:i))
                    ylim([0 2])
                    subplot(2,2,4)
                    [vDef,Height] = obj.get_force_curve_data(i,'AppRetSwitch',AppRetSwitch,...
                        'BaselineCorrection',1,'TipHeightCorrection',0,...
                        'Sensitivity','original','Unit','m');
                    switch FitRangeMode
                        case 'ValueHorizontal'
                            UValue = max(Height) - FitRangeLowerValue;
                            LValue = max(Height) - FitRangeUpperValue;
                            FitHeight = Height(Height>=LValue & Height<=UValue);
                            FitvDef = vDef(Height>=LValue & Height<=UValue);
                        case 'ValueVertical'
                            FitHeight = Height(vDef>=FitRangeLowerValue & vDef<=FitRangeUpperValue);
                            FitvDef = vDef(vDef>=FitRangeLowerValue & vDef<=FitRangeUpperValue);
                        case 'FractionHorizontal'
                            RangeHeight = range(Height);
                            FitHeight = Height(Height>=FitRangeLowerValue*RangeHeight & Height<=FitRangeUpperValue*RangeHeight);
                            FitvDef = vDef(Height>=FitRangeLowerValue*RangeHeight & Height<=FitRangeUpperValue*RangeHeight);
                        case 'FractionVertical'
                            MaxvDef = range(vDef);
                            FitHeight = Height(vDef>=FitRangeLowerValue*MaxvDef & vDef<=FitRangeUpperValue*MaxvDef);
                            FitvDef = vDef(vDef>=FitRangeLowerValue*MaxvDef & vDef<=FitRangeUpperValue*MaxvDef);
                    end
                    plot(Height,vDef,FitHeight,FitvDef,FitHeight,polyval(Params(i,:),FitHeight))
                    xlim([max(FitHeight)-range(FitHeight)*2 max(FitHeight)])
                    legend({'Data','Fit Data','Fitted Slope'})
                    title(['Current Point: ' num2str(i) '    CurrentSlope: ' num2str(AdaptiveMedianRefSlope(i))])
                    drawnow
                end
            end
            
            obj.AdaptiveSensitivity = obj.Sensitivity./AdaptiveMedianRefSlope;
            
        end
        
        function set_reference_slope_to_user_input(obj)
            IsValidInput = false;
            while ~IsValidInput
                UsrInput = inputdlg(sprintf('Reference Slope of %s',obj.Name),'Manual Reference Slope');
                if ~isempty(regexp(UsrInput{1},'[^0-9.]', 'once')) || isempty(UsrInput{1})
                    IsValidInput = false;
                    Warn = warndlg('Input can only contain numbers and periods. Please try again');
                    uiwait(Warn)
                    continue
                end
                obj.RefSlope = str2double(UsrInput{1});
                obj.RefSlopeCorrectedSensitivity = obj.Sensitivity/obj.RefSlope;
                obj.HasRefSlope = true;
                IsValidInput = true;
            end
        end
        
        function set_reference_slope_to_value(obj,Value)
            obj.RefSlope = Value; % BEST.FUNCTION.EVER.WRITTEN.
            obj.RefSlopeCorrectedSensitivity = obj.Sensitivity/obj.RefSlope;
            obj.HasRefSlope = true; % so true
        end
        
        function Mask = create_mask_general(obj)
            
            Mask = logical(zeros(obj.NumPixelsX,obj.NumPixelsY));
            CheckSum = 100;
            while CheckSum > 1
                f = figure('Name','Choose areas');
                f.WindowState = 'maximized';
                subplot(2,1,2)
                surf(imresize(imrotate(obj.HeightMap(:,:,1)',90),[1024 1024]),'LineStyle','none','FaceLighting','gouraud','FaceColor','interp')
                light('Style','local')
                subplot(2,1,1)
                imshow(obj.HeightMap(:,:,1).*~Mask,[min(obj.HeightMap(:,:,1),[],'all') max(obj.HeightMap(:,:,1),[],'all')])
                title(sprintf('%s: Draw Freehand ROI around areas, that are to be included in the mask\nThe area will be marked and the same map redrawn \n If you are done masking just click on the image once without dragging the cursor',obj.Name))
                ROI = drawfreehand;
                CheckSum = length(ROI.Waypoints);
                TempMask = createMask(ROI);
                Mask = Mask | TempMask;
                close(f)
            end
        end
        
        function read_jpk_images_from_files(obj)
            
            if obj.PythonLoaderFlag
                current = what;
                TempFolder = fullfile(obj.DataStoreFolder,obj.ID);
                warning('off', 'all');
                mkdir(TempFolder)
                warning('on', 'all');
                cd(TempFolder)
                if isequal(obj.FileType,'quantitative-imaging-map')
                    obj.OpenZipFile.extract('data-image.jpk-qi-image');
                elseif isequal(obj.FileType,'force-scan-map')
                    obj.OpenZipFile.extract('data-image.force');
                end
                cd(current.path)
                Folder = TempFolder;
            else
                Folder = obj.DataStoreFolder;
            end
            
            if isequal(obj.FileType,'quantitative-imaging-map')
                I = AFMImage(fullfile(Folder,'data-image.jpk-qi-image'));
            elseif isequal(obj.FileType,'force-scan-map')
                I = AFMImage(fullfile(Folder,'data-image.force'));
                for i=1:length(I.Channel)
                    I.Channel(i).Image = fliplr(I.Channel(i).Image);
                end
            end
            
            I.preprocess_image;
            
            obj.Channel = I.Channel;
            
            Processed = obj.get_channel('Processed');
            
            obj.HeightMap = Processed.Image;
            
            if obj.PythonLoaderFlag
                rmdir(Folder,'s');
            end
            
        end
        
        function create_and_level_height_map(obj, ForceFraction)
            % first set Height Map to default values for reproducable
            % results
            
            obj.construct_list_to_map_relations();
            
            if nargin < 2
                ForceFraction = 1;
            end
            
            Max = zeros(obj.NCurves,1);
            for i=1:obj.NCurves
                if ForceFraction == 1
                [Force,HHApp] = obj.get_force_curve_data(i,'AppRetSwitch',0,...
                    'BaselineCorrection',0,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                    Max(i) = -max(HHApp);
                else
                [Force,HHApp] = obj.get_force_curve_data(i,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                    interp_Force = linspace(min(Force), obj.Setpoint, numel(Force));
                    [~, ForceIdx] = min(abs(interp_Force - ForceFraction * obj.Setpoint));
                    Max(i) = -HHApp(ForceIdx);
                end
            end
            TempHeightMap = obj.convert_data_list_to_map(Max);
            
            if ForceFraction == 1
                Channel1Name = 'Indented Height';
                Channel2Name = 'Processed Indented Height';
            else
                Percent = ForceFraction*100;
                Channel1Name = ['Indented Height ForceFraction ' num2str(Percent) '%' ];
                Channel2Name = ['Processed Indented Height ForceFraction ' num2str(Percent) '%' ];
            end
            
            % write to Channel
            Height = obj.create_standard_channel(TempHeightMap,Channel1Name,'m');
            
            obj.add_channel(Height,true);
            
            if size(obj.HeightMap,1) < 128
                Map = imresize(TempHeightMap,[256 256],'nearest');
            elseif size(obj.HeightMap,1) < 512
                Map = imresize(TempHeightMap,[512 512],'nearest');
            else
                Map = imresize(TempHeightMap,[1024 1024],'nearest');
            end
            for i=1:5
                Map = AFMImage.subtract_line_fit_vertical_rov(Map,.2,0);
            end
            Map = imresize(Map,[obj.NumPixelsX obj.NumPixelsY],'nearest');
            
            try
                Map = AFMImage.find_and_replace_outlier_lines(Map,10);
            catch
                warning('Could not replace outlier lines')
            end
            
            % write to Channel
            Processed = obj.create_standard_channel(Map,Channel2Name,'m');
            
            obj.add_channel(Processed,true);
            
        end
        
        function create_fibril_mask(obj,MaskParam)
            
            if nargin < 2
                MaskParam = 0.5;
            end
            HghtRange = range(obj.HeightMap,'all');
            HghtMin = min(obj.HeightMap,[],'all');
            HghtNorm = (obj.HeightMap-HghtMin)/HghtRange;
            mask = zeros(size(HghtNorm));
            STDLine = zeros(obj.NumPixelsX,obj.NumPixelsY);
            [HghtSorted,Idx] = sort(HghtNorm,[2],'descend');
            for j=1:obj.NumPixelsY
                STDLine(:,j) = std(HghtSorted(:,1:j),0,[2]);
            end
            [~,MaxIdx] = max(STDLine,[],[2]);
            for i=1:obj.NumPixelsX
                mask(i,Idx(i,1:floor(MaxIdx*MaskParam))) = 1;
            end
            mask = logical(mask);
            mask = bwareafilt(mask,1,4);
            
            % determine linear orientation angle of the fibril and pad
            % accordingly
            k = 1;
            for i=1:obj.NumPixelsX
                for j = 1:obj.NumPixelsY
                    if mask(i,j) == 1
                        xy(k,:) = [i,j];
                        k = k + 1;
                    end
                end
            end
            FitLine = fit(xy(:,1).*obj.ScanSizeX/obj.NumPixelsX,xy(:,2).*obj.ScanSizeY/obj.NumPixelsY,'poly1');
            obj.FibrilFlag.Straight = 0;
            if abs(FitLine.p1) <= 0.05
                obj.FibrilFlag.Straight = 1;
                for i=1:10
                    mask = medfilt2(mask,[5 3],'symmetric');
%                     imshow(mask)
                end
                mask = bwareafilt(mask,1,4);
            end
            obj.FibMask = mask;
            
            %             current = what();
            %             cd(obj.Folder)
            %             savename = sprintf('%s.mat',obj.Name);
            %             save(savename,'obj')
            %             cd(current.path)
        end
        
        function create_automatic_background_mask(obj,PercentOfRange)
            
            Channel = obj.get_channel('Processed');
            
            obj.BackgroundMask = AFMImage.mask_background_by_threshold(Channel.Image,PercentOfRange);
        end
        
        function check_for_new_host(obj)
            % check_for_new_host(obj)
            %
            % Checks, if the system environment has changed and, if so,
            % resets CPFlag.CNNopt
            
            FullOS = computer;
            OS = FullOS(1:3);
            if isequal(OS,'PCW')
                Host = getenv('COMPUTERNAME');
            elseif isequal(OS,'GLN')
                Host = getenv('HOSTNAME');
            elseif isequal(OS,'MAC')
                Host = getenv('HOSTNAME');
            end
            
            if isequal(obj.HostOS,OS) && isequal(obj.HostName,Host)
                return
            else
                obj.CPFlag.CNNopt = 0;
                obj.HostOS = OS;
                obj.HostName = Host;
            end
        end
        
        function cnn_runtime_optimization(obj,NeuralNet,X)
            % cnn_runtime_optimization(obj,NeuralNet,X)
            %
            % optimizes runtime for neural network prediction by adjusting
            % MiniBatchSize
            
            len = size(X,4);
            if obj.CPFlag.CNNopt == 0
                CantHandle = true;
                obj.MiniBatchSize = min(1024,len);
                DynMBSdone = false;
                HasFailed = false;
                NumIter = 0;
                while CantHandle == true
                    if NumIter >= 20
                        error('Neural Net keeps failing despite adjusting MiniBatchSize. Call your programm administrator and an exorcist')
                    end
                    try
                        predict(NeuralNet,X,'MiniBatchSize',obj.MiniBatchSize,'Acceleration','auto',...
                                        'ExecutionEnvironment',obj.NeuralNetAccelerator);
                        CantHandle = false;
                        if DynMBSdone == false
                            if HasFailed == true
                                obj.MiniBatchSize = ceil(obj.MiniBatchSize/3);
                            end
                            fprintf('Dynamically adjusting MiniBatchSize = %i\n',obj.MiniBatchSize)
                            DynMBSdone = true;
                        end
                    catch ME
                        switch ME.identifier
                            case 'parallel:gpu:array:OOM'
                                obj.MiniBatchSize = ceil(obj.MiniBatchSize*3/4);
                                fprintf('Dynamically adjusting MiniBatchSize = %i\n',obj.MiniBatchSize)
                                HasFailed = true;
                            case 'nnet_cnn:internal:cnn:layer:CustomLayer:PredictErrored'
                                obj.MiniBatchSize = ceil(obj.MiniBatchSize*3/4);
                                fprintf('Dynamically adjusting MiniBatchSize = %i\n',obj.MiniBatchSize)
                                HasFailed = true;
                            case 'nnet_cnn:dlAccel:MEXCallFailed'
                                obj.MiniBatchSize = ceil(obj.MiniBatchSize*3/4);
                                fprintf('Dynamically adjusting MiniBatchSize = %i\n',obj.MiniBatchSize)
                                HasFailed = true;
                            otherwise
                                obj.MiniBatchSize = ceil(obj.MiniBatchSize*3/4);
                                fprintf('Dynamically adjusting MiniBatchSize = %i\n',obj.MiniBatchSize)
                                HasFailed = true;
                        end
                    end
                    NumIter = NumIter + 1;
                end
                obj.CPFlag.CNNopt = 1;
            end
            
        end
        
        function cnn_zoom_in(obj,ZoomFactor)
            % cnn_zoom_in(obj)
            %
            % Caution! do not use this function on your ForceMaps! This is
            % an aux method for the Zoom mode in CNN CP estimation and is
            % only run on throwaway ForceMap copys
            
            iRange = find(obj.SelectedCurves)';
            for i=iRange
                [App,HHApp] = obj.get_force_curve_data(i,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                CutOff = min(HHApp) + range(HHApp(HHApp<=obj.CP_CNNZoom(i,1)))*ZoomFactor;
                HHApp(HHApp<CutOff) = [];
                App(1:(end-length(HHApp))) = [];
            end
            
            
        end
        
        function unpack_force_map_data(obj,MapFullFile,DataFolder)
            
            FileExtension = regexp(MapFullFile, '\.[^.]+$', 'match', 'once');
            
            if obj.BigDataFlag
                if obj.PythonLoaderFlag
                    TempFolderName = 'DataStore';
                    % Check if the directory already exists
                    if ~exist(fullfile(DataFolder, TempFolderName), 'dir')
                        % Create the directory if it does not exist
                        mkdir(DataFolder, TempFolderName);
                    end
                    TempFolder = fullfile(DataFolder,TempFolderName);
                else
                    TempFolderName = sprintf('FM_DataStore_%s',obj.ID);
                end
            else
                TempFolderName = sprintf('Temp%s',obj.ID);
            end
            
            if (obj.BigDataFlag && ~obj.PythonLoaderFlag) || ~obj.BigDataFlag
                if isequal('PCW',obj.HostOS)
                    % unpack jpk-file into temporary folder to read out data
                    cmd1 = '"C:\Program Files\7-Zip\7z.exe" x ';
                    cmd2 = '"';
                    cmd3 = MapFullFile;
                    cmd4 = '"';
                    cmd5 = ' -o';
                    mkdir(DataFolder,TempFolderName)
                    cmd6 = '"';
                    TempFolder = fullfile(DataFolder,TempFolderName,filesep);
                    cmd8 = '"';
                    CMD = append(cmd1,cmd2,cmd3,cmd4,cmd5,cmd6,TempFolder,cmd8);
                    disp('extracting file...')
                    h = system(CMD);
                elseif isequal('GLN',obj.HostOS)
                    % unpack jpk-file into temporary folder to read out data
                    cmd1 = 'unzip -o ';
                    cmd2 = MapFullFile;
                    cmd3 = ' -d ';
                    mkdir(DataFolder,TempFolderName)
                    TempFolder = fullfile(DataFolder,TempFolderName,filesep);
                    CMD = append(cmd1,cmd2,cmd3,TempFolder);
                    h = system(CMD);
                    disp('extracting file...')
                elseif isequal('MAC',obj.HostOS)
                    % unpack jpk-file into temporary folder to read out data
                    cmd1 = 'unzip -o ';
                    cmd2 = MapFullFile;
                    cmd3 = ' -d ';
                    mkdir(DataFolder,TempFolderName)
                    TempFolder = fullfile(DataFolder,TempFolderName,filesep);
                    CMD = append(cmd1,cmd2,cmd3,TempFolder);
                    h = system(CMD);
                    disp('extracting file...')
                end
                if h==0
                    disp('unzipping successfull')
                elseif h==1
                    disp('unzipping failed')
                end
            elseif obj.BigDataFlag && obj.PythonLoaderFlag
                % Copy the jpk file to the DataStore folder and open it
                % with the python loader.
                copyfile(MapFullFile,TempFolder)
                Split = split(MapFullFile,filesep);
                obj.RawDataFilePath = fullfile(TempFolder,Split{end});
                switch FileExtension
                    case '.nhf'
                        obj.NHF_FileInfo = h5info(MapFullFile);
                    case {'.jpk-qi-data','.jpk-force-map'}
                        obj.load_zipped_files_with_python
                end
            end
            
            
            Strings = split(MapFullFile,filesep);
            switch FileExtension
                case '.nhf'
                    exp1 = '.*'; % Finds the typed-in name of the user during the AFM experiment
                    NameCell = regexp(Strings(end,1), exp1, 'match');
                    obj.Name = NameCell{1}{1};
                    DateTime = obj.nhf_get_group_attribute_value('spectroscopy','created');
                    parts = regexp(DateTime, '([0-9-]+)T([0-9:.]+)Z', 'tokens');
                    obj.Date = parts{1}{1};
                    obj.Time = parts{1}{2};
                    obj.Time=strrep(obj.Time,'.','-'); % Remove dots in obj.Time
                    obj.Time=strrep(obj.Time,':','-'); % Remove dots in obj.Time
                case {'.jpk-force-map','.jpk-qi-data'}
                    %%% Define the search expressions
                    % Comment: JPK includes a few attributes of a measurement into the name:
                    % The name attachement is the same for all experiments:
                    % 'Typedname'+'-'+'data'+'-'+'year'+'.'+'months'+'.'+'day'+'-hour'+'.'+'minute'+'.'+'second'+'.'+'thousandths'+'.'+'jpk file extension'
                    exp1 = '.*(?=\-data)'; % Finds the typed-in name of the user during the AFM experiment
                    exp2 = '\d{4}\.\d{2}\.\d{2}'; % Finds the date
                    exp3 = '\d{2}\.\d{2}\.\d{2}\.\d{3}'; % Finds the time
                    obj.Name = regexp(Strings(end,1), exp1, 'match');
                    obj.Name = char(obj.Name{1});
                    if isequal(obj.Name,'')
                        exp4 = '.*(?=.jpk)';
                        obj.Name = regexp(Strings(end,1), exp4, 'match');
                        obj.Name = char(obj.Name{1});
                    end
                    obj.Date = regexp(Strings(end,1), exp2, 'match');
                    obj.Date = char(obj.Date{1});
                    obj.Time = regexp(Strings(end,1), exp3, 'match');
                    obj.Time = char(obj.Time{1});
                    obj.Date=strrep(obj.Date,'.','-'); % Remove dots in obj.Date
                    obj.Time=strrep(obj.Time,'.','-'); % Remove dots in obj.Time
            end
            
            obj.Folder = [DataFolder replace(obj.Name,'.','') '-' obj.ID];
            obj.DataStoreFolder = TempFolder;
            
        end
        
        function read_in_header_properties(obj)
            % Check for file type and get important ForceMap
            % properties
            
            FileExtension = regexp(obj.RawDataFilePath, '\.[^.]+$', 'match', 'once');
            
            switch FileExtension
                case {'.jpk-force-map','.jpk-qi-data'}
                obj.jpk_read_in_header_properties;
                case '.nhf'
                obj.nhf_read_in_header_properties;
            end
            
        end
        
        function nhf_read_in_header_properties(obj)
            
            obj.FileType = 'nhf-spectroscopy';
            
            obj.FileVersion = obj.nhf_get_group_attribute_value('spectroscopy','software_version');
            RectAxisSize = obj.nhf_get_group_attribute_value('spectroscopy','rect_axis_size');
            obj.NCurves = double(RectAxisSize(1)*RectAxisSize(2));
            obj.NumPixelsX = double(RectAxisSize(1));
            obj.NumPixelsY = double(RectAxisSize(2));
            [RectAxisRange,~,GroupIndex] = obj.nhf_get_group_attribute_value('spectroscopy','rect_axis_range');
            obj.ScanSizeX = double(RectAxisRange(1));
            obj.ScanSizeY = double(RectAxisRange(2));
            obj.SpringConstant = double(obj.nhf_get_group_attribute_value('spectroscopy','spm_probe_calibration_spring_constant'));
            obj.Sensitivity = double(obj.nhf_get_group_attribute_value('spectroscopy','spm_probe_calibration_deflection_sensitivity'));
            obj.TipRadius = double(obj.nhf_get_group_attribute_value('spectroscopy','spm_probe_nominal_tip_radius'));
            obj.Medium = obj.nhf_get_group_attribute_value('spectroscopy','measurement_environment');
            obj.ScanAngle = double(obj.nhf_get_group_attribute_value('spectroscopy','rect_rotation'));
            
            JSONConfigApp = AFMBaseClass.nhf_get_segment_attribute_value(...
                obj.NHF_FileInfo.Groups(GroupIndex),'Advance to Setpoint 1','segment_configuration');
            JSONConfigRet = AFMBaseClass.nhf_get_segment_attribute_value(...
                obj.NHF_FileInfo.Groups(GroupIndex),'Retract 1','segment_configuration');
            JSONConfigHold = AFMBaseClass.nhf_get_segment_attribute_value(...
                obj.NHF_FileInfo.Groups(GroupIndex),'Wait 1','segment_configuration');
            
            obj.MaxPointsPerCurve =  double(string(AFMBaseClass.nhf_get_json_name_value(JSONConfigApp,'datapoints')));
            obj.NumSegments = numel(obj.NHF_FileInfo.Groups(GroupIndex).Groups);
            obj.HoldingTime = double(string(AFMBaseClass.nhf_get_json_name_value(JSONConfigHold,'duration')));
            obj.ExtendTime = double(string(AFMBaseClass.nhf_get_json_name_value(JSONConfigApp,'duration')));
            obj.ExtendZLength =  double(string(AFMBaseClass.nhf_get_json_name_value(JSONConfigApp,'length')));
            obj.RetractTime = double(string(AFMBaseClass.nhf_get_json_name_value(JSONConfigRet,'duration')));
            obj.RetractZLength =  double(string(AFMBaseClass.nhf_get_json_name_value(JSONConfigRet,'length')));
            obj.ExtendVelocity =  double(string(AFMBaseClass.nhf_get_json_name_value(JSONConfigApp,'speed')));
            obj.RetractVelocity =  double(string(AFMBaseClass.nhf_get_json_name_value(JSONConfigRet,'speed')));
            obj.HHType = 'Position Z';
            try
                obj.Setpoint =  double(string(AFMBaseClass.nhf_get_json_name_value(JSONConfigApp,'setpoint')))...
                    *obj.Sensitivity*obj.SpringConstant;
            catch
                obj.Setpoint = [];
            end
            if isempty(obj.Setpoint)
                try
                    obj.Setpoint =  double(string(AFMBaseClass.nhf_get_json_name_value(JSONConfigApp,'set_point')))...
                        *obj.Sensitivity*obj.SpringConstant;
                catch
                    obj.Setpoint = [];
                end
            end
            
        end
        
        function jpk_read_in_header_properties(obj)
            
            if obj.PythonLoaderFlag
                current = what;
                TempFolder = fullfile(obj.DataStoreFolder,obj.ID);
                TempSharedFolder = fullfile(obj.DataStoreFolder,obj.ID,'shared-data');
                mkdir(TempFolder)
                cd(TempFolder)
                obj.OpenZipFile.extract('header.properties');
                obj.OpenZipFile.extract('shared-data/header.properties');
                cd(current.path)
                
                filedirectory = fullfile(TempFolder,'header.properties');
                FileDirectoryShared = fullfile(TempSharedFolder,'header.properties');
            else
                filedirectory = fullfile(obj.DataStoreFolder,'header.properties');
                FileDirectoryShared = fullfile(obj.DataStoreFolder,'shared-data','header.properties');
            end
            
            FileID=fopen(filedirectory,'rt','n','UTF-8'); % FileID = fopen(filename,permission,machinefmt,encodingIn)
            FileCont=fileread(filedirectory);
            % Shared-data file directory
            FileIDShared=fopen(FileDirectoryShared,'rt','n','UTF-8'); % FileID = fopen(filename,permission,machinefmt,encodingIn)
            FileContShared=fileread(FileDirectoryShared);
            % Height: 1. CONVERSION raw-meters & 2. SCALE meters
            % Conversion RAW -> VOLTS
            fseek(FileID,1,'cof'); % goes at the first position in the file
            
            %   Check for file type (.jpk-force-map, .jpk-qi-data)
            frewind(FileID);
            StrPos=strfind(FileCont,'jpk-data-file=');
            fseek(FileID,StrPos,'cof');
            tline = fgetl(FileID);
            where=strfind(tline,'=');
            TempType = tline(where+1:end);
            if isequal(TempType, 'spm-quantitative-image-data-file')   % Valid for software versions 6.1.158  
                obj.FileType = 'quantitative-imaging-map';
            elseif isequal(TempType,'spm-force-scan-map-file')         % Valid for software versions 6.1.158 
                obj.FileType = 'force-scan-map';
            end            
            % Restore initial conditions
            frewind(FileID); % Move file position indicator back to beginning of the open file
                       
            %%   Check for the file version            
            StrPos=strfind(FileCont,strcat(obj.FileType,'.description.source-software=')); % finds the position of the string searched for
            fseek(FileID,StrPos,'cof'); % Move to specified position in file
            % 'cof': Current position in file
            TextLine = fgetl(FileID); % reads the rest of the line starting from the defined string until the end of the line
            LinePos=strfind(TextLine,'='); % finds the position of the '='-sign in the line
            obj.FileVersion = TextLine(LinePos+1:end); % allocates all found next to the '='-sign until the end of the line
            % Restore initial conditions
            frewind(FileID); % Move file position indicator back to beginning of the open file

            %   Number of curves (COMMENT: JPK starts from 0 not from 1)
            StrPos=strfind(FileCont,strcat(obj.FileType,'.indexes.max=')); % strfind(file,string) is looking for a specific string in the file.
            fseek(FileID,StrPos,'cof'); % moves at the location where specific string is located
            TextLine = fgetl(FileID); % reads the rest of the line starting from the defined string until the end of the line
            % stores that string in a character
            LinePos=strfind(TextLine,'='); % finds the position of the '='-sign in the line
            obj.NCurves = 1 + str2double(... % convert the string to number and add '1' to start from 1 not from 0
                TextLine(LinePos+1:end)... % read out the number
                );
            
            %   NumPixelsX
            clear tline where;
            frewind(FileID);
            StrPos=strfind(FileCont,strcat(obj.FileType,'.position-pattern.grid.jlength='));
            fseek(FileID,StrPos,'cof');
            tline = fgetl(FileID);
            where=strfind(tline,'=');
            obj.NumPixelsX = str2double(tline(where+1:end));
            
            %   NumPixelsY
            clear tline where;
            frewind(FileID);
            StrPos=strfind(FileCont,strcat(obj.FileType,'.position-pattern.grid.ilength='));
            fseek(FileID,StrPos,'cof');
            tline = fgetl(FileID);
            where=strfind(tline,'=');
            obj.NumPixelsY = str2double(tline(where+1:end));
            
            %   ScanSizeX
            clear tline where;
            frewind(FileID);
            StrPos=strfind(FileCont,strcat(obj.FileType,'.position-pattern.grid.ulength='));
            fseek(FileID,StrPos,'cof');
            tline = fgetl(FileID);
            where=strfind(tline,'=');
            obj.ScanSizeX = str2double(tline(where+1:end));
            
            %   ScanSizeY
            clear tline where;
            frewind(FileID);
            StrPos=strfind(FileCont,strcat(obj.FileType,'.position-pattern.grid.vlength='));
            fseek(FileID,StrPos,'cof');
            tline = fgetl(FileID);
            where=strfind(tline,'=');
            obj.ScanSizeY = str2double(tline(where+1:end));
            
            %   MaxPointsPerCurve
            clear tline where;
            frewind(FileID);
            if isequal(obj.FileType,'force-scan-map')
                StrPos=strfind(FileCont,strcat(obj.FileType,'.settings.force-settings.type='));
                fseek(FileID,StrPos,'cof');
                TextLine = fgetl(FileID);
                LinePos=strfind(TextLine,'=');
                TempType=TextLine(LinePos+1:end);
                % Restore initial conditions
                frewind(FileID); % Move file position indicator back to beginning of the open file
                if isequal(TempType,'relative-force-settings')
                    StrPos1=strfind(FileCont,strcat(obj.FileType,'.settings.force-settings.extend-k-length='));
                    fseek(FileID,StrPos1,'cof');
                    TextLine = fgetl(FileID);
                    LinePos=strfind(TextLine,'=');
                    obj.MaxPointsPerCurve=str2double(TextLine(LinePos+1:end));
                    % Restore initial conditions
                    frewind(FileID); % Move file position indicator back to beginning of the open file
                elseif isequal(TempType,'segmented-force-settings')                   
                    % Number of data points
                    StrPos2=strfind(FileCont,strcat(obj.FileType,'.settings.force-settings.segment.0.num-points='));          
                    fseek(FileID,StrPos2,'cof');
                    TextLine = fgetl(FileID);
                    LinePos=strfind(TextLine,'=');
                    obj.MaxPointsPerCurve=str2double(TextLine(LinePos+1:end));
                    % Restore initial conditions
                    frewind(FileID); % Move file position indicator back to beginning of the open file
                end                    
            elseif isequal(obj.FileType,'quantitative-imaging-map')
                StrPos=strfind(FileCont,strcat(obj.FileType,'.settings.force-settings.extend.num-points='));
                fseek(FileID,StrPos,'cof');
                TextLine = fgetl(FileID);
                LinePos=strfind(TextLine,'=');
                obj.MaxPointsPerCurve=str2double(TextLine(LinePos+1:end));
                % Restore initial conditions
                frewind(FileID); % Move file position indicator back to beginning of the open file
            end            

            % Number of Segments
            clear tline where;
            frewind(FileIDShared);
                StrPos=strfind(FileContShared,'force-segment-header-infos.count=');
                fseek(FileIDShared,StrPos,'cof');
                TextLine = fgetl(FileIDShared);
                LinePos=strfind(TextLine,'=');
                obj.NumSegments=str2double(TextLine(LinePos+1:end)); % jpk assigns approach and retraction as well as possible holding times into segments
                % Restore initial conditions
                frewind(FileIDShared); % Move file position indicator back to beginning of the open file                   
                      
            %% Holding time
            if isequal(TempType,'relative-force-settings')
                    % Holding time
                    StrPos=strfind(FileCont,strcat(obj.FileType,'.settings.force-settings.extended-pause-time'));
                    fseek(FileID,StrPos,'cof');
                    TextLine = fgetl(FileID);
                    LinePos=strfind(TextLine,'=');
                    obj.HoldingTime = str2double(TextLine(LinePos+1:end));
                    % Restore initial conditions
                    frewind(FileID); % Move file position indicator back to beginning of the open file  
            elseif isequal(TempType,'segmented-force-settings') && obj.NumSegments == 3
                    % Holding time
                    StrPos=strfind(FileCont,strcat(obj.FileType,'.settings.force-settings.segment.1.duration'));
                    fseek(FileID,StrPos,'cof');
                    TextLine = fgetl(FileID);
                    LinePos=strfind(TextLine,'=');
                    obj.HoldingTime = str2double(TextLine(LinePos+1:end));
                    % Restore initial conditions
                    frewind(FileID); % Move file position indicator back to beginning of the open file                    
            end
            
            %%   Velocity
            if isequal(obj.FileType,'force-scan-map')
                frewind(FileID); % Move file position indicator back to beginning of the open file
                if isequal(TempType,'relative-force-settings')
                    % Extend
                    % Scan time
                    StrPos=strfind(FileCont,strcat(obj.FileType,'.settings.force-settings.extend-scan-time='));
                    fseek(FileID,StrPos,'cof');
                    TextLine = fgetl(FileID);
                    LinePos=strfind(TextLine,'=');
                    obj.ExtendTime = str2double(TextLine(LinePos+1:end));
                    % Restore initial conditions
                    frewind(FileID); % Move file position indicator back to beginning of the open file
                    % z-height
                    StrPos=strfind(FileCont,strcat(obj.FileType,'.settings.force-settings.relative-z-start='));
                    fseek(FileID,StrPos,'cof');
                    TextLine = fgetl(FileID);
                    LinePos=strfind(TextLine,'=');
                    obj.ExtendZLength = str2double(TextLine(LinePos+1:end));
                    % Restore initial conditions
                    frewind(FileID); % Move file position indicator back to beginning of the open file
                    % Retract
                    % Scan time
                    StrPos=strfind(FileCont,strcat(obj.FileType,'.settings.force-settings.retract-scan-time='));
                    fseek(FileID,StrPos,'cof');
                    TextLine = fgetl(FileID);
                    LinePos=strfind(TextLine,'=');
                    obj.RetractTime = str2double(TextLine(LinePos+1:end));
                    % Restore initial conditions
                    frewind(FileID); % Move file position indicator back to beginning of the open file
                    % z-height
                    StrPos=strfind(FileCont,strcat(obj.FileType,'.settings.force-settings.relative-z-start='));
                    fseek(FileID,StrPos,'cof');
                    TextLine = fgetl(FileID);
                    LinePos=strfind(TextLine,'=');
                    obj.RetractZLength = str2double(TextLine(LinePos+1:end));
                    % Restore initial conditions
                    frewind(FileID); % Move file position indicator back to beginning of the open file
                elseif isequal(TempType,'segmented-force-settings')
                    % Extend
                    % Scan time
                    StrPos=strfind(FileCont,strcat(obj.FileType,'.settings.force-settings.segment.0.duration='));
                    fseek(FileID,StrPos,'cof');
                    TextLine = fgetl(FileID);
                    LinePos=strfind(TextLine,'=');
                    obj.ExtendTime = str2double(TextLine(LinePos+1:end));
                    % Restore initial conditions
                    frewind(FileID); % Move file position indicator back to beginning of the open file
                    % z-height
                    StrPos=strfind(FileCont,strcat(obj.FileType,'.settings.force-settings.segment.0.z-start'));
                    fseek(FileID,StrPos,'cof');
                    TextLine = fgetl(FileID);
                    LinePos=strfind(TextLine,'=');
                    obj.ExtendZLength = str2double(TextLine(LinePos+1:end));
                    % Restore initial conditions
                    frewind(FileID); % Move file position indicator back to beginning of the open file
                    % Retract
                    % Scan time
                    StrPos=strfind(FileCont,strcat(obj.FileType,sprintf('.settings.force-settings.segment.%d.duration=',obj.NumSegments-1)));
                    fseek(FileID,StrPos,'cof');
                    TextLine = fgetl(FileID);
                    LinePos=strfind(TextLine,'=');
                    obj.RetractTime = str2double(TextLine(LinePos+1:end));
                    % Restore initial conditions
                    frewind(FileID); % Move file position indicator back to beginning of the open file
                    % z-height
                    StrPos=strfind(FileCont,strcat(obj.FileType,sprintf('.settings.force-settings.segment.%d.z-start',obj.NumSegments-1)));                  
                    fseek(FileID,StrPos,'cof');
                    TextLine = fgetl(FileID);
                    LinePos=strfind(TextLine,'=');
                    obj.RetractZLength = str2double(TextLine(LinePos+1:end));
                    % Restore initial conditions
                    frewind(FileID); % Move file position indicator back to beginning of the open file   
                end     
            elseif isequal(obj.FileType,'quantitative-imaging-map')
                clear tline where;
                frewind(FileID);
                StrPos=strfind(FileCont,'quantitative-imaging-map.settings.force-settings.extend.duration=');
                fseek(FileID,StrPos,'cof');
                tline = fgetl(FileID);
                where=strfind(tline,'=');
                obj.ExtendTime = str2double(tline(where+1:end));
                
                clear tline where;
                frewind(FileID);
                StrPos=strfind(FileCont,'quantitative-imaging-map.settings.force-settings.extend.duration=');
                fseek(FileID,StrPos,'cof');
                tline = fgetl(FileID);
                where=strfind(tline,'=');
                obj.RetractTime = str2double(tline(where+1:end));
                
                clear tline where;
                frewind(FileID);
                StrPos=strfind(FileCont,'quantitative-imaging-map.settings.force-settings.extend.z-start=');
                fseek(FileID,StrPos,'cof');
                tline = fgetl(FileID);
                where=strfind(tline,'=');
                obj.ExtendZLength = str2double(tline(where+1:end));
                obj.RetractZLength = obj.ExtendZLength;
            end
            
            obj.ExtendVelocity = obj.ExtendZLength/obj.ExtendTime;
            obj.RetractVelocity = obj.RetractZLength/obj.RetractTime;
            
            %   ScanAngle
            clear tline where;
            frewind(FileID);
            StrPos=strfind(FileCont,strcat(obj.FileType,'.position-pattern.grid.theta='));
            fseek(FileID,StrPos,'cof');
            tline = fgetl(FileID);
            where=strfind(tline,'=');
            obj.ScanAngle = str2double(tline(where+1:end));
            obj.ScanAngle = obj.ScanAngle*180/pi;
            
            obj.HHType = 'capacitiveSensorHeight';  
            obj.read_in_rescaling_constants(FileDirectoryShared)
            
            %   Setpoint
            clear tline where;
            frewind(FileID);
            if isequal(obj.FileType,'quantitative-imaging-map')
                StrPos=strfind(FileCont,strcat(obj.FileType,'.settings.force-settings.extend.setpoint='));
                fseek(FileID,StrPos,'cof');
                tline = fgetl(FileID);
                where=strfind(tline,'=');
                obj.Setpoint = str2double(tline(where+1:end)).*obj.Sensitivity.*obj.SpringConstant;
                % Restore initial conditions
                frewind(FileID); % Move file position indicator back to beginning of the open file
            elseif isequal(obj.FileType,'force-scan-map')
                if isequal(TempType,'relative-force-settings')
                    % Extend                   
                    StrPos=strfind(FileCont,strcat(obj.FileType,'.settings.force-settings.relative-setpoint='));
                    fseek(FileID,StrPos,'cof');
                    tline = fgetl(FileID);
                    where=strfind(tline,'=');
                    obj.Setpoint = str2double(tline(where+1:end)).*obj.Sensitivity.*obj.SpringConstant;
                    % Restore initial conditions
                    frewind(FileID); % Move file position indicator back to beginning of the open file                
                elseif isequal(TempType,'segmented-force-settings')
                    % Extend                   
                    StrPos=strfind(FileCont,strcat(obj.FileType,'.settings.force-settings.segment.0.duration='));
                    fseek(FileID,StrPos,'cof');
                    tline = fgetl(FileID);
                    where=strfind(tline,'=');
                    obj.Setpoint = str2double(tline(where+1:end)).*obj.Sensitivity.*obj.SpringConstant;
                    % Restore initial conditions
                    frewind(FileID); % Move file position indicator back to beginning of the open file
                end                
            end                                 
            
            fclose(FileIDShared);
            fclose(FileID);
            
            
%             if obj.PythonLoaderFlag
%                 rmdir(TempFolder,'s');
%             end
        end
        
        function read_in_rescaling_constants(obj,SharedPropertiesFullFile)
            
            obj.RescalingConstants = struct();
            
            ChannelNames = {'capacitiveSensorHeight',...
                                'measuredHeight',...
                                'head-height',...
                                'height'};
            for i=1:4
                obj.RescalingConstants(i).ChannelName = ChannelNames{i};
                obj.RescalingConstants(i).exists = false;
                try
                    [mult_height_meters1, offset_height_meters1,...
                        mult_height_meters2, offset_height_meters2,...
                        mult_vDefl_volts, offset_vDefl_volts,...
                        sensitivity, spring_constant] = ForceMap.getheaderinfo(SharedPropertiesFullFile,ChannelNames{i});
                    
                    obj.RescalingConstants(i).Mult1 = mult_height_meters1;
                    obj.RescalingConstants(i).Off1 = offset_height_meters1;
                    obj.RescalingConstants(i).Mult2 = mult_height_meters2;
                    obj.RescalingConstants(i).Off2 = offset_height_meters2;
                    obj.RescalingConstants(i).exists = true;
                catch
                end
            end
            
            % Special case vDef
            obj.RescalingConstants(5).ChannelName = 'vDeflection';
            obj.RescalingConstants(5).Mult1 = mult_vDefl_volts;
            obj.RescalingConstants(5).Off1 = offset_vDefl_volts;
            obj.RescalingConstants(5).Mult2 = 1;
            obj.RescalingConstants(5).Off2 = 0;
            obj.RescalingConstants(5).exists = true;
            obj.SpringConstant = spring_constant;
            obj.Sensitivity = sensitivity;
            
            if obj.RescalingConstants(...
                    reshape(strcmp({obj.RescalingConstants.ChannelName}, obj.HHType),...
                    size(obj.RescalingConstants))).exists
            elseif obj.RescalingConstants(...
                    reshape(strcmp({obj.RescalingConstants.ChannelName}, 'measuredHeight'),...
                    size(obj.RescalingConstants))).exists
                obj.HHType = 'measuredHeight';
            elseif obj.RescalingConstants(...
                    reshape(strcmp({obj.RescalingConstants.ChannelName}, 'head-height'),...
                    size(obj.RescalingConstants))).exists
                obj.HHType = 'head-height';
            elseif obj.RescalingConstants(...
                    reshape(strcmp({obj.RescalingConstants.ChannelName}, 'height'),...
                    size(obj.RescalingConstants))).exists
                obj.HHType = 'height';
            end
        end
        
        function load_force_curves(obj)
            TempFolder = obj.DataStoreFolder;
            obj.HHType = 'capacitiveSensorHeight';
            for i=1:obj.NCurves
                HeaderFileDirectory = fullfile(TempFolder,'shared-data','header.properties');
                % 0 segment = approach part
                SegmentHeaderFileDirectory = fullfile(TempFolder,'index',string((i-1)),'segments','0','segment-header.properties');
                HeightDataDirectory = fullfile(TempFolder,'index',string((i-1)),'segments','0','channels','capacitiveSensorHeight.dat');
                vDefDataDirectory = fullfile(TempFolder,'index',string((i-1)),'segments','0','channels','vDeflection.dat');
                
                if ~isfile(HeightDataDirectory) || isequal(obj.HHType,'measuredHeight')
                    HeightDataDirectory = fullfile(TempFolder,'index',string((i-1)),'segments','0','channels','measuredHeight.dat');
                    obj.HHType = 'measuredHeight';
                end
                if ~isfile(HeightDataDirectory) || isequal(obj.HHType,'Height')
                    HeightDataDirectory = fullfile(TempFolder,'index',string((i-1)),'segments','0','channels','Height.dat');
                    obj.HHType = 'Height';
                end
                
                try
                    [TempHHApp,obj.App{i},obj.SpringConstant,obj.Sensitivity]=...
                        obj.writedata(HeaderFileDirectory,SegmentHeaderFileDirectory,...
                        HeightDataDirectory,vDefDataDirectory,obj.HHType);
                    
                    obj.HHApp{i} = -TempHHApp;
                    obj.App{i} = obj.App{i}.*obj.SpringConstant;
                    clear TempHHApp
                catch
                    if i ~= 1
                        disp(sprintf('Curve Nr. %i seems to be corrupted. Replacing with previous value instead',i))
                        obj.HHApp{i} = obj.HHApp{i-1};
                        obj.App{i} = obj.App{i-1};
                    else
                        disp(sprintf('Curve Nr. %i seems to be corrupted. Replacing with zeros instead',i))
                        obj.HHApp{i} = zeros(obj.MaxPointsPerCurve,1);
                        obj.App{i} = zeros(obj.MaxPointsPerCurve,1);
                        obj.SelectedCurves(i) = 0;
                    end
                end
                
                % Following segments
                
                SegmentHeaderFileDirectory = fullfile(TempFolder,'index',string((i-1)),'segments',string(obj.NumSegments-1),'segment-header.properties');
                HeightDataDirectory = fullfile(TempFolder,'index',string((i-1)),'segments',string(obj.NumSegments-1),'channels','capacitiveSensorHeight.dat');
                vDefDataDirectory = fullfile(TempFolder,'index',string((i-1)),'segments',string(obj.NumSegments-1),'channels','vDeflection.dat');
                
                if ~isfile(HeightDataDirectory) || isequal(obj.HHType,'measuredHeight')
                    HeightDataDirectory = fullfile(TempFolder,'index',string((i-1)),'segments',string(obj.NumSegments-1),'channels','measuredHeight.dat');
                end
                if ~isfile(HeightDataDirectory) || isequal(obj.HHType,'Height')
                    HeightDataDirectory = fullfile(TempFolder,'index',string((i-1)),'segments',string(obj.NumSegments-1),'channels','Height.dat');
                    obj.HHType = 'Height';
                end
                
                try
                    [TempHHRet,obj.Ret{i}]=...
                        obj.writedata(HeaderFileDirectory,SegmentHeaderFileDirectory,...
                        HeightDataDirectory,vDefDataDirectory,obj.HHType);
                    
                    obj.HHRet{i} = -TempHHRet;
                    obj.Ret{i} = obj.Ret{i}.*obj.SpringConstant;
                catch
                    if i ~= 1
                        disp(sprintf('Curve Nr. %i seems to be corrupted. Replacing with previous value instead',i))
                        obj.HHRet{i} = obj.HHRet{i-1};
                        obj.Ret{i} = obj.Ret{i-1};
                        clear TempHHRet
                    else
                        disp(sprintf('Curve Nr. %i seems to be corrupted. Replacing with zeros instead',i))
                        obj.HHRet{i} = zeros(obj.MaxPointsPerCurve,1);
                        obj.Ret{i} = zeros(obj.MaxPointsPerCurve,1);
                    end
                end
            end
        end
        
        function [OutvDef,OutHeight] = get_force_curve_data(obj,CurveNumber,varargin)
            % function [OutvDef,OutHeight] = get_force_curve_data(obj,CurveNumber,varargin)
            %
            % Get channel data from connected file container or class property and apply
            % several corrections/operations to it. 
            %
            %
            % Required inputs
            % obj ... ForceMap instance
            % CurveNumber ... Index of curve to be loaded (counting from 1 to obj.NCurves)
            %
            % Name-Value pairs
            % "AppRetSwitch" ... 0 for approach curve, 1 for retract curve
            % "BaselineCorrection" ... logical for base and tilt subtraction
            % "TipHeightCorrection" ... logical to decide if OutHeight is
            %                           HeadHeight (=0) or TipHeight (=1)
            % "Sensitivity" ... type of sensitivity value:
            %                   'original','corrected','adaptive'
            % "Unit" ... conversion unit of OutvDef: 'N','m','V'
            
            p = inputParser;
            p.FunctionName = "get_force_curve_data";
            p.CaseSensitive = false;
            p.PartialMatching = true;
            
            % Required inputs
            validobj = @(x)true;
            validCurveNumber = @(x)true;
            addRequired(p,"obj",validobj);
            addRequired(p,"CurveNumber",validCurveNumber);
            
            % NameValue inputs
            defaultAppRetSwitch = 0;
            defaultBaselineCorrection = false;
            defaultTipHeightCorrection = false;
            defaultSensitivity = 'original';
            defaultUnit = 'N';
            validAppRetSwitch = @(x)x==0|x==1|islogical(x);
            validBaselineCorrection = @(x)x==0|x==1|islogical(x);
            validTipHeightCorrection = @(x)x==0|x==1|islogical(x);
            validSensitivity = @(x)any(validatestring(x,{'original','corrected','adaptive'}));
            validUnit = @(x)any(validatestring(x,{'N','m','V'}));
            addParameter(p,"AppRetSwitch",defaultAppRetSwitch,validAppRetSwitch);
            addParameter(p,"BaselineCorrection",defaultBaselineCorrection,validBaselineCorrection);
            addParameter(p,"TipHeightCorrection",defaultTipHeightCorrection,validTipHeightCorrection);
            addParameter(p,"Sensitivity",defaultSensitivity,validSensitivity);
            addParameter(p,"Unit",defaultUnit,validUnit);
            
            parse(p,obj,CurveNumber,varargin{:});
            
            % Assign parsing results to named variables
            obj = p.Results.obj;
            CurveNumber = p.Results.CurveNumber;
            AppRetSwitch = p.Results.AppRetSwitch;
            BaselineCorrection = p.Results.BaselineCorrection;
            TipHeightCorrection = p.Results.TipHeightCorrection;
            Sensitivity = p.Results.Sensitivity;
            Unit = p.Results.Unit;
            
            
            if CurveNumber > obj.NCurves
                error('Requested curve index (%i) is bigger than maximum number of curves (%i)',CurveNumber,obj.NCurves);
            end
            if (BaselineCorrection || TipHeightCorrection) && ~obj.BaseAndTiltFlag
                error('No Base and Tilt data found');
            end
            if isempty(obj.BigDataFlag) || ~obj.BigDataFlag
                if AppRetSwitch==0 && BaselineCorrection==0 && TipHeightCorrection==0
                    OutvDef = obj.App{CurveNumber};
                    OutHeight = obj.HHApp{CurveNumber};
                elseif AppRetSwitch==1 && BaselineCorrection==0 && TipHeightCorrection==0
                    OutvDef = obj.Ret{CurveNumber};
                    OutHeight = obj.HHRet{CurveNumber};
                elseif AppRetSwitch==0 && BaselineCorrection==1 && TipHeightCorrection==0
                    OutvDef = obj.App{CurveNumber};
                    OutHeight = obj.HHApp{CurveNumber};
                    FittedApp = polyval(obj.Basefit{CurveNumber},OutHeight);
                    OutvDef = OutvDef - FittedApp;
                elseif AppRetSwitch==1 && BaselineCorrection==1 && TipHeightCorrection==0
                    OutvDef = obj.Ret{CurveNumber};
                    OutHeight = obj.HHRet{CurveNumber};
                    FittedRet = polyval(obj.Basefit{CurveNumber},OutHeight);
                    OutvDef = OutvDef - FittedRet;
                elseif AppRetSwitch==0 && BaselineCorrection==0 && TipHeightCorrection==1
                    OutvDef = obj.App{CurveNumber};
                    OutHeight = obj.HHApp{CurveNumber};
                    FittedApp = polyval(obj.Basefit{CurveNumber},OutHeight);
                    OutvDef = OutvDef - FittedApp;
                    OutHeight = OutHeight - OutvDef./obj.SpringConstant;
                    OutvDef = obj.App{CurveNumber};
                elseif AppRetSwitch==1 && BaselineCorrection==0 && TipHeightCorrection==1
                    OutvDef = obj.Ret{CurveNumber};
                    OutHeight = obj.HHRet{CurveNumber};
                    FittedRet = polyval(obj.Basefit{CurveNumber},OutHeight);
                    OutvDef = OutvDef - FittedRet;
                    OutHeight = OutHeight - OutvDef./obj.SpringConstant;
                    OutvDef = obj.Ret{CurveNumber};
                elseif AppRetSwitch==0 && BaselineCorrection==1 && TipHeightCorrection==1
                    OutvDef = obj.App{CurveNumber};
                    OutHeight = obj.HHApp{CurveNumber};
                    FittedApp = polyval(obj.Basefit{CurveNumber},OutHeight);
                    OutvDef = OutvDef - FittedApp;
                    OutHeight = OutHeight - OutvDef./obj.SpringConstant;
                elseif AppRetSwitch==1 && BaselineCorrection==1 && TipHeightCorrection==1
                    OutvDef = obj.Ret{CurveNumber};
                    OutHeight = obj.HHRet{CurveNumber};
                    FittedRet = polyval(obj.Basefit{CurveNumber},OutHeight);
                    OutvDef = OutvDef - FittedRet;
                    OutHeight = OutHeight - OutvDef./obj.SpringConstant;
                end
                if ~isequal(Unit,'V')
                    switch Sensitivity
                        case 'original'
                        case 'corrected'
                            OutvDef = OutvDef./obj.Sensitivity.*obj.RefSlopeCorrectedSensitivity;
                        case 'adaptive'
                            OutvDef = OutvDef./obj.Sensitivity.*obj.AdaptiveSensitivity(CurveNumber);
                    end
                end
                switch Unit
                    case 'N'
                    case 'm'
                        OutvDef = OutvDef./obj.SpringConstant;
                    case 'V'
                        OutvDef = OutvDef./(obj.SpringConstant*obj.Sensitivity);
                end
                return
            end
            
            TempFolder = obj.DataStoreFolder;
            
            % Correct for potential holding segment
            if AppRetSwitch
                AppRetSwitch=obj.NumSegments-1;
            end
            
            AppRetSwitch = int8(AppRetSwitch);
            
            if ~obj.PythonLoaderFlag
                HeaderFileDirectory = fullfile(TempFolder,'shared-data','header.properties');
                SegmentHeaderFileDirectory = fullfile(TempFolder,'index',string((CurveNumber-1)),'segments',string(AppRetSwitch),'segment-header.properties');
                HeightDataDirectory = fullfile(TempFolder,'index',string((CurveNumber-1)),'segments',string(AppRetSwitch),'channels','capacitiveSensorHeight.dat');
                vDefDataDirectory = fullfile(TempFolder,'index',string((CurveNumber-1)),'segments',string(AppRetSwitch),'channels','vDeflection.dat');
                
                if ~isfile(HeightDataDirectory) || isequal(obj.HHType,'measuredHeight')
                    HeightDataDirectory = fullfile(TempFolder,'index',string((CurveNumber-1)),'segments',string(AppRetSwitch),'channels','measuredHeight.dat');
                    obj.HHType = 'measuredHeight';
                end
                if ~isfile(HeightDataDirectory) || isequal(obj.HHType,'height')
                    HeightDataDirectory = fullfile(TempFolder,'index',string((CurveNumber-1)),'segments',string(AppRetSwitch),'channels','height.dat');
                    obj.HHType = 'height';
                end
            end
            
            try
                if obj.CorruptedCurves(CurveNumber)
                    error('skip this one')
                end
                
                if isequal(obj.FileType,'nhf-spectroscopy')
                    
                    TempHeight = obj.nhf_load_single_curve_channel_data(CurveNumber,AppRetSwitch,obj.HHType);
                    OutHeight = -TempHeight;
                    OutvDef = obj.nhf_load_single_curve_channel_data(CurveNumber,AppRetSwitch,'Deflection');
                    
                elseif isequal(obj.FileType,'jpk-force-map') ||...
                        isequal(obj.FileType,'jpk-qi-data') ||...
                        isequal(obj.FileType,'force-scan-map') ||...
                        isequal(obj.FileType,'quantitative-imaging-map')
                    if obj.PythonLoaderFlag
                        TempHeight = obj.load_single_curve_channel_data_with_python(string(CurveNumber-1),string(AppRetSwitch),obj.HHType);
                        OutvDef = obj.load_single_curve_channel_data_with_python(string(CurveNumber-1),string(AppRetSwitch),'vDeflection');
                    else
                        [TempHeight,OutvDef,SC,Sens]=...
                            obj.writedata(HeaderFileDirectory,SegmentHeaderFileDirectory,...
                            HeightDataDirectory,vDefDataDirectory,obj.HHType);
                    end
                    if isempty(obj.Sensitivity)
                        obj.Sensitivity = Sens;
                    end
                    if isempty(obj.SpringConstant)
                        obj.SpringConstant = SC;
                    end
                    
                    OutHeight = -TempHeight;
                    OutvDef = OutvDef.*obj.SpringConstant;
                end
                
                if AppRetSwitch==0 && BaselineCorrection==0 && TipHeightCorrection==0
                elseif AppRetSwitch==1 && BaselineCorrection==0 && TipHeightCorrection==0
                elseif AppRetSwitch==0 && BaselineCorrection==1 && TipHeightCorrection==0
                    FittedApp = polyval(obj.Basefit{CurveNumber},OutHeight);
                    OutvDef = OutvDef - FittedApp;
                elseif AppRetSwitch==1 && BaselineCorrection==1 && TipHeightCorrection==0
                    FittedRet = polyval(obj.Basefit{CurveNumber},OutHeight);
                    OutvDef = OutvDef - FittedRet;
                elseif AppRetSwitch==0 && BaselineCorrection==0 && TipHeightCorrection==1
                    TempApp = OutvDef;
                    FittedApp = polyval(obj.Basefit{CurveNumber},OutHeight);
                    OutvDef = OutvDef - FittedApp;
                    OutHeight = OutHeight - OutvDef./obj.SpringConstant;
                    OutvDef = TempApp;
                elseif AppRetSwitch==1 && BaselineCorrection==0 && TipHeightCorrection==1
                    TempRet = OutvDef;
                    FittedRet = polyval(obj.Basefit{CurveNumber},OutHeight);
                    OutvDef = OutvDef - FittedRet;
                    OutHeight = OutHeight - OutvDef./obj.SpringConstant;
                    OutvDef = TempRet;
                elseif AppRetSwitch==0 && BaselineCorrection==1 && TipHeightCorrection==1
                    FittedApp = polyval(obj.Basefit{CurveNumber},OutHeight);
                    OutvDef = OutvDef - FittedApp;
                    OutHeight = OutHeight - OutvDef./obj.SpringConstant;
                elseif AppRetSwitch==1 && BaselineCorrection==1 && TipHeightCorrection==1
                    FittedRet = polyval(obj.Basefit{CurveNumber},OutHeight);
                    OutvDef = OutvDef - FittedRet;
                    OutHeight = OutHeight - OutvDef./obj.SpringConstant;
                end
                
                if ~isequal(Unit,'V')
                    switch Sensitivity
                        case 'original'
                        case 'corrected'
                            OutvDef = OutvDef./obj.Sensitivity.*obj.RefSlopeCorrectedSensitivity;
                        case 'adaptive'
                            OutvDef = OutvDef./obj.Sensitivity.*obj.AdaptiveSensitivity(CurveNumber);
                    end
                end
                switch Unit
                    case 'N'
                    case 'm'
                        OutvDef = OutvDef./obj.SpringConstant;
                    case 'V'
                        OutvDef = OutvDef./(obj.SpringConstant*obj.Sensitivity);
                end
                
                
            catch
                disp(sprintf('Curve Nr. %i seems to be corrupted. Replacing with next non corrupted curve instead',CurveNumber))
                obj.CorruptedCurves(CurveNumber) = true;
                obj.SelectedCurves(CurveNumber) = false;
                k = 1;
                while obj.CorruptedCurves(mod(CurveNumber+k-1,obj.NCurves)+1)
                    k = k + 1;
                    if k > obj.NCurves
                        error('All force curves seem to be corrupted');
                    end
                end
                [OutvDef,OutHeight] = obj.get_force_curve_data(mod(CurveNumber+k-1,obj.NCurves)+1,...
                    'AppRetSwitch',AppRetSwitch,'BaselineCorrection',BaselineCorrection,...
                    'TipHeightCorrection',TipHeightCorrection,'Sensitivity',Sensitivity,...
                    'Unit',Unit);
            end
        end
        
        function OutVector = load_single_curve_channel_data_with_python(obj,Index,Segment,ChannelName)
            
            % First check if zip file is even loaded.If not, load it
            Meta = metaclass(obj.OpenZipFile);
            if ~isequal(Meta.Name,'py.zipfile.ZipFile')
                warning('Zip file not loaded. Loading... You may clear the files again using .clear_zipped_files_from_memory')
                obj.load_zipped_files_with_python;
            end
            
            TargetFile = strcat('index/',Index,'/segments/',Segment,'/channels/',strcat(ChannelName,'.dat'));
            
            File = obj.OpenZipFile.open(TargetFile);
            
            Data = uint8(File.read(-1));
            OutVector = double(swapbytes(typecast(Data,'int32')))';
            File.close()
            
            % Rescale stored data to real world units
            OutVector = obj.rescale_channel_data(OutVector,ChannelName);
            
        end
        
        function OutVector = rescale_channel_data(obj,InVector,ChannelName)
            
            for i=1:5
                if isequal(ChannelName,obj.RescalingConstants(i).ChannelName)
                    TempVector = InVector.*obj.RescalingConstants(i).Mult1 + obj.RescalingConstants(i).Off1;
                    OutVector = TempVector.*obj.RescalingConstants(i).Mult2 + obj.RescalingConstants(i).Off2;
                    if isequal(ChannelName,'vDeflection')
                        OutVector = OutVector.*obj.Sensitivity;
                    end
                end
            end
        end
        
        function OutVector = nhf_load_single_curve_channel_data(obj,Index,Segment,ChannelName)
            
            if ~Segment
                SegmentName = 'Advance to Setpoint 1';
            else
                SegmentName = 'Retract 1';
            end
            
            [~,~,GroupIndex] = obj.nhf_get_group_attribute_value('spectroscopy','');
            [~,SegmentIndex] = AFMBaseClass.nhf_get_segment_attribute_value(...
                obj.NHF_FileInfo.Groups(GroupIndex),SegmentName,'segment_configuration');
            
            ChannelIndex = obj.nhf_search_dataset_with_matching_attribute_name_value(...
                obj.NHF_FileInfo.Groups(GroupIndex).Groups(SegmentIndex).Datasets,'signal_name',ChannelName);
            
            BlockSizeID = obj.NHF_FileInfo.Groups(GroupIndex).Groups(SegmentIndex).Datasets(ChannelIndex).Attributes(...
                find(contains({obj.NHF_FileInfo.Groups(GroupIndex).Groups(SegmentIndex).Datasets(ChannelIndex).Attributes.Name},'dataset_block_size_source'))).Value;
            
            ChunkingIndex = obj.nhf_search_dataset_with_matching_attribute_name_value(...
                obj.NHF_FileInfo.Groups(GroupIndex).Groups(SegmentIndex).Datasets,'dataset_block_size_id',BlockSizeID);
            
            ChunkingData = h5read(obj.RawDataFilePath,...
                ['/group_' sprintf('%04d',(GroupIndex - 1)) '/subgroup_' sprintf('%04d',(SegmentIndex - 1)) '/dataset_' sprintf('%04d',(ChunkingIndex - 1))]);
            
            ChunkIndizes = [1 ; cumsum(ChunkingData)]+1;
            
            Start = double(ChunkIndizes(Index));
            Count = double(ChunkingData(Index));
            
            Data = h5read(obj.RawDataFilePath,...
                ['/group_' sprintf('%04d',(GroupIndex - 1)) '/subgroup_' sprintf('%04d',(SegmentIndex - 1)) '/dataset_' sprintf('%04d',(ChannelIndex - 1))],...
                Start,Count);
            
            % Calibrate Data
            Data = obj.nhf_calibrate_channel_data(Data,GroupIndex,SegmentIndex,ChannelIndex);
            
            OutVector = Data;
            %             OutVector = Data(ChunkIndizes(Index)+1:ChunkIndizes(Index+1));
            %             OutVector = Data(ChunkIndizes(2:end)-Index);
        end
        
        function [OutVector,ChunkIndizes,ChunkingData] = nhf_load_full_channel_data(obj,Segment,ChannelName)
            
            if ~Segment
                SegmentName = 'Advance to Setpoint 1';
            else
                SegmentName = 'Retract 1';
            end
            
            [~,~,GroupIndex] = obj.nhf_get_group_attribute_value('spectroscopy','');
            [~,SegmentIndex] = AFMBaseClass.nhf_get_segment_attribute_value(...
                obj.NHF_FileInfo.Groups(GroupIndex),SegmentName,'segment_configuration');
            
            ChannelIndex = obj.nhf_search_dataset_with_matching_attribute_name_value(...
                obj.NHF_FileInfo.Groups(GroupIndex).Groups(SegmentIndex).Datasets,'signal_name',ChannelName);
            
            BlockSizeID = obj.NHF_FileInfo.Groups(GroupIndex).Groups(SegmentIndex).Datasets(ChannelIndex).Attributes(...
                find(contains({obj.NHF_FileInfo.Groups(GroupIndex).Groups(SegmentIndex).Datasets(ChannelIndex).Attributes.Name},'dataset_block_size_source'))).Value;
            
            ChunkingIndex = obj.nhf_search_dataset_with_matching_attribute_name_value(...
                obj.NHF_FileInfo.Groups(GroupIndex).Groups(SegmentIndex).Datasets,'dataset_block_size_id',BlockSizeID);
            
            ChunkingData = h5read(obj.RawDataFilePath,...
                ['/group_' sprintf('%04d',(GroupIndex - 1)) '/subgroup_' sprintf('%04d',(SegmentIndex - 1)) '/dataset_' sprintf('%04d',(ChunkingIndex - 1))]);
            
            ChunkIndizes = [1 ; cumsum(ChunkingData)];
            
            Data = h5read(obj.RawDataFilePath,...
                ['/group_' sprintf('%04d',(GroupIndex - 1)) '/subgroup_' sprintf('%04d',(SegmentIndex - 1)) '/dataset_' sprintf('%04d',(ChannelIndex - 1))]);
            
            % Calibrate Data
            Data = obj.nhf_calibrate_channel_data(Data,GroupIndex,SegmentIndex,ChannelIndex);
            
            OutVector = Data;
%             OutVector = Data(ChunkIndizes(Index)+1:ChunkIndizes(Index+1));
%             OutVector = Data(ChunkIndizes(2:end)-Index);
        end
        
        function nhf_create_map_from_spectroscopy_data(obj,ChannelName,FancyName,Unit)
            
            [Data,ChunkIndizes] = obj.nhf_load_full_channel_data(0,ChannelName);
            
            MapData = Data(ChunkIndizes(2:end));
            
            Image = obj.convert_data_list_to_map(MapData);
            
            ImageChannel = obj.create_standard_channel(Image,FancyName,Unit);
            
            obj.add_channel(ImageChannel);
        end
        
        function Data = nhf_calibrate_channel_data(obj,Data,GroupIndex,SegmentIndex,ChannelIndex)
            
            [~,~,CalibrationIndex] = obj.nhf_get_group_attribute_value('calibration','');
            CalibrationSource = 'signal_calibration_source';
            
            % How many Calibration steps?
            NumCalibs = sum(contains({obj.NHF_FileInfo.Groups(GroupIndex).Groups(SegmentIndex).Datasets(ChannelIndex).Attributes.Name},CalibrationSource));
            
            if NumCalibs == 1
                CalibSourceNames = {CalibrationSource};
            else
                for i=1:NumCalibs
                    CalibSourceNames{i} = [CalibrationSource sprintf('_%d',i-1)];
                end
            end
            
            Data = double(Data);
            for i=1:NumCalibs
                CalibrationUUID = obj.NHF_FileInfo.Groups(GroupIndex).Groups(SegmentIndex).Datasets(ChannelIndex).Attributes(...
                    find(contains({obj.NHF_FileInfo.Groups(GroupIndex).Groups(SegmentIndex).Datasets(ChannelIndex).Attributes.Name},CalibSourceNames{i}))).Value;
                
                % Find Calibration source
                CalibDataIndex = obj.nhf_search_dataset_with_matching_attribute_name_value(obj.NHF_FileInfo.Groups(CalibrationIndex).Datasets,'signal_calibration_id',CalibrationUUID);
                
                MappingType = obj.NHF_FileInfo.Groups(CalibrationIndex).Datasets(CalibDataIndex).Attributes(...
                    find(contains({obj.NHF_FileInfo.Groups(CalibrationIndex).Datasets(CalibDataIndex).Attributes.Name},'signal_mapping_type'))).Value;
                
                MappingParameters = obj.NHF_FileInfo.Groups(CalibrationIndex).Datasets(CalibDataIndex).Attributes(...
                    find(contains({obj.NHF_FileInfo.Groups(CalibrationIndex).Datasets(CalibDataIndex).Attributes.Name},'signal_mapping_parameters'))).Value;
                
                switch MappingType
                    case {'linear','linear_non_invertible'}
                        Data = Data.*MappingParameters(1) + MappingParameters(2);
                    case 'exponential'
                        Data = MappingParameters(2).*exp((Data - MappingParameters(3))./MappingParameters(1));
                    case 'logarithmic'
                        Data = MappingParameters(1).*log(Data./MappingParameters(2)) + MappingParameters(3);
                end
                
            end
            
        end
        
        function initialize_flags(obj)
            % initialize all flags related to the ForceMap class
            
            obj.BaseAndTiltFlag = false;
            obj.CPFlag.SnapIn = 0;
            obj.CPFlag.HertzFitted = 0;
            obj.CPFlag.RoV = 0;
            obj.CPFlag.GoF = 0;
            obj.CPFlag.Combo = 0;
            obj.CPFlag.CNN = 0;
            obj.CPFlag.CNNZoom = 0;
            obj.CPFlag.CNNZoomDropout = 0;
            obj.CPFlag.CNNZoomSweep = 0;
            obj.CPFlag.Dropout = 0;
            obj.CPFlag.Manual = 0;
            obj.CPFlag.Old = 0;
            obj.CPFlag.CNNopt = 0;
            obj.CPFlag.HardSurface = 0;
            obj.HasRefSlope = false;
            % SMFS
            obj.SMFSFlag.SelectFM=ones(1,obj.NCurves);
            obj.SMFSFlag.Preprocessed=ones(1,obj.NCurves);
            obj.SMFSFlag.Presorted=ones(1,obj.NCurves);
            obj.SMFSFlag.Uncorrupt=ones(1,obj.NCurves);
            obj.SMFSFlag.Min=zeros(1,obj.NCurves);
            obj.SMFSFlag.Length=zeros(1,obj.NCurves);
            obj.SMFSFlag.Fit=zeros(1,obj.NCurves);      
        end
        
        function create_dummy_force_map(obj,NSynthCurves)
            % creates blank bare miinimum ForceMap for synthetic force
            % curves needed in certain CP scripts
            
            obj.Name = 'DummyForceMap';
            obj.NCurves = NSynthCurves;
            obj.SelectedCurves = ones(obj.NCurves,1);
            obj.CPFlag.RoV = 0;
            obj.CPFlag.GoF = 0;
            obj.CPFlag.Combo = 0;
            obj.CPFlag.CNN = 0;
            obj.CPFlag.CNNZoom = 0;
            obj.CPFlag.CNNZoomDropout = 0;
            obj.CPFlag.CNNZoomSweep = 0;
            obj.CPFlag.Dropout = 0;
            obj.CPFlag.Manual = 0;
            obj.CPFlag.Old = 0;
            obj.CPFlag.CNNopt = 0;
            obj.Ret = cell(1,NSynthCurves);
            obj.BasedRet = cell(1,NSynthCurves);
            obj.create_and_level_height_map();
            return
        end
        
        function temporary_data_load_in(obj,OnOffBool)
            
            if OnOffBool && obj.BigDataFlag
                if isequal(obj.FileType,'nhf-spectroscopy')
                    [DefAppData,DefAppChunkIndizes,DefAppChunkSizes] = ...
                        obj.nhf_load_full_channel_data(0,'Deflection');
                    [HHAppData,HHAppChunkIndizes,HHAppChunkSizes] = ...
                        obj.nhf_load_full_channel_data(0,'Position Z');
                    [DefRetData,DefRetChunkIndizes,DefRetChunkSizes] = ...
                        obj.nhf_load_full_channel_data(1,'Deflection');
                    [HHRetData,HHRetChunkIndizes,HHRetChunkSizes] = ...
                        obj.nhf_load_full_channel_data(1,'Position Z');
                    for i=1:obj.NCurves
                        obj.App{i} = DefAppData(DefAppChunkIndizes(i)+1:DefAppChunkIndizes(i+1));
                        obj.HHApp{i} = -HHAppData(HHAppChunkIndizes(i)+1:HHAppChunkIndizes(i+1));
                        obj.Ret{i} = DefRetData(DefRetChunkIndizes(i)+1:DefRetChunkIndizes(i+1));
                        obj.HHRet{i} = -HHRetData(HHRetChunkIndizes(i)+1:HHRetChunkIndizes(i+1));
                    end
                else
                    for i=1:obj.NCurves
                        [obj.App{i},obj.HHApp{i}] = obj.get_force_curve_data(i,'AppRetSwitch',0,...
                            'BaselineCorrection',0,'TipHeightCorrection',0,...
                            'Sensitivity','original','Unit','N');
                        [obj.Ret{i},obj.HHRet{i}] = obj.get_force_curve_data(i,'AppRetSwitch',1,...
                            'BaselineCorrection',0,'TipHeightCorrection',0,...
                            'Sensitivity','original','Unit','N');
                    end
                end
                obj.BigDataFlag = 0;
            elseif ~OnOffBool
                obj.BigDataFlag = 1;
                obj.App = cell(0,0);
                obj.HHApp = cell(0,0);
                obj.Ret = cell(0,0);
                obj.HHRet = cell(0,0);
            end
        end
        
        function check_for_cuda_capable_gpu_device(obj)
            % Check if current machine is capable of gpu CUDA accel.
            try
                GPU = gpuDevice;
                if GPU.DeviceAvailable && GPU.DeviceSupported
                    obj.NeuralNetAccelerator = 'gpu'; % Drake approves
                else
                    obj.NeuralNetAccelerator = 'cpu'; % YIKES, get a new system
                end
            catch ME
                obj.NeuralNetAccelerator = 'cpu'; % YIKES, get a new system
            end
            
        end
        
        function reset_selected_and_corrupted_flags(obj)
            
            obj.SelectedCurves(:) = 1;
            obj.CorruptedCurves(:) = 0;
            
        end
        
        function convert_hertzfit_property_to_fractioned_properties(obj)
            
            for i=1:length(obj.HertzFit)
                try
                    Formula = formula(obj.HertzFit{i});
                    Coeffnames = coeffnames(obj.HertzFit{i});
                    Coeffs = coeffvalues(obj.HertzFit{i});
                    
                    obj.HertzFitType = Formula;
                    obj.HertzFitCoeffNames = Coeffnames;
                    obj.HertzFitValues{i} = Coeffs;
                catch
                end
            end
            
            obj.HertzFit = [];
            
        end
        
        function IsLoaded = check_if_zipped_file_is_loaded(obj,AskIfLoadDialogueBool)
            
            if nargin < 2
                AskIfLoadDialogueBool = false;
            end
            
            if isequal(obj.FileType,'nhf-spectroscopy')
                IsLoaded = true;
                return
            end
            
            if obj.PythonLoaderFlag && obj.BigDataFlag
                Meta = metaclass(obj.OpenZipFile);
                if ~isequal(Meta.Name,'py.zipfile.ZipFile')
                    IsLoaded = false;
                    if AskIfLoadDialogueBool
                        answer = questdlg('File container is not loaded. Load to memory?',...
                            'Load Zip file',...
                            'Yes','No','No');
                        if isequal(answer,'Yes')
                            h = waitbar(1/3,'Loading can take a few minutes',...
                                'Name','Loading file container');
                            obj.load_zipped_files_with_python;
                            IsLoaded = true;
                            close(h)
                        end 
                    end
                else
                    IsLoaded = true;
                end
            else
                IsLoaded = true;
                return
            end
            
        end
        
    end
    
    methods
        % methods for visualization, plotting, statistics and quality control
        
        function plotting_show_force_curve(obj,ZoomMult,k,fig)
            if nargin < 2
                jRange = find(obj.SelectedCurves);
                k = jRange(randi(length(jRange)));
                ZoomMult = 0;
            end
            if nargin < 3
                jRange = find(obj.SelectedCurves);
                k = jRange(randi(length(jRange)));
            end
            if nargin < 4
                fig = figure('Units','normalized',...
                    'Position',[0.6 0.1 0.4 0.8],...
                    'Color','white',...
                    'Name',obj.Name);
            end
            
            try
                [AppY,AppX] = obj.get_force_curve_data(k,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                [RetY,RetX] = obj.get_force_curve_data(k,'AppRetSwitch',1,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
            catch
                [AppY,AppX] = obj.get_force_curve_data(k,'AppRetSwitch',0,...
                    'BaselineCorrection',0,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                [RetY,RetX] = obj.get_force_curve_data(k,'AppRetSwitch',1,...
                    'BaselineCorrection',0,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
            end
            
            
            subplot(2,1,1)
            title(sprintf('Curve Nr.%i of %s',k,obj.Name))
            hold on
            [MultiplierX,UnitX,~] = AFMImage.parse_unit_scale(range(RetX),'m',10);
            [MultiplierY,UnitY,~] = AFMImage.parse_unit_scale(range(RetY),'N',5);
            plot(AppX*MultiplierX,AppY*MultiplierY,RetX*MultiplierX,RetY*MultiplierY,'LineWidth',1.5);
            Legends = {'Approach','Retract'};
            
            if obj.CPFlag.HertzFitted == 1 && ~isempty(obj.CP_HertzFitted)
                plot(obj.CP_HertzFitted(k,1)*MultiplierX, obj.CP_HertzFitted(k,2)*MultiplierY,'O',...
                    'LineWidth',1.5,...
                    'MarkerSize',7,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor',[0.4940 0.1840 0.5560]);
                Legends{end+1} = 'Origin of Hertz-Sneddon fit';
                xlim([(AppX(1)+ZoomMult*(obj.CP_HertzFitted(k,1) - AppX(1)))*MultiplierX inf]);
            end
            if obj.CPFlag.SnapIn == 1
                plot(obj.CP_SnapIn(k,1)*MultiplierX, obj.CP_SnapIn(k,2)*MultiplierY,'O',...
                    'LineWidth',1.5,...
                    'MarkerSize',7,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','r');
                Legends{end+1} = 'SnapIn';
                xlim([(AppX(1)+ZoomMult*(obj.CP_SnapIn(k,1) - AppX(1)))*MultiplierX inf]);
            end
            if obj.CPFlag.Manual == 1
                plot(obj.Man_CP(k,1)*MultiplierX, obj.Man_CP(k,2)*MultiplierY,'O',...
                    'LineWidth',1.5,...
                    'MarkerSize',7,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','r');
                Legends{end+1} = 'Manual CP';
                xlim([(AppX(1)+ZoomMult*(obj.Man_CP(k,1) - AppX(1)))*MultiplierX inf]);
            end
            if obj.CPFlag.RoV == 1
                plot(obj.CP_RoV(k,1)*MultiplierX, obj.CP_RoV(k,2)*MultiplierY,'O',...
                    'LineWidth',1.5,...
                    'MarkerSize',7,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','b');
                Legends{end+1} = 'CP RoV';
                xlim([(AppX(1)+ZoomMult*(obj.CP_RoV(k,1) - AppX(1)))*MultiplierX inf]);
            end
            if obj.CPFlag.GoF == 1
                plot(obj.CP_GoF(k,1)*MultiplierX, obj.CP_GoF(k,2)*MultiplierY,'O',...
                    'LineWidth',1.5,...
                    'MarkerSize',7,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','c');
                Legends{end+1} = 'CP GoF';
                xlim([(AppX(1)+ZoomMult*(obj.CP_GoF(k,1) - AppX(1)))*MultiplierX inf]);
            end
            if obj.CPFlag.Combo == 1
                plot(obj.CP_Combo(k,1)*MultiplierX, obj.CP_Combo(k,2)*MultiplierY,'O',...
                    'LineWidth',1.5,...
                    'MarkerSize',7,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','m');
                Legends{end+1} = 'CP Combo';
                xlim([(AppX(1)+ZoomMult*(obj.CP_Combo(k,1) - AppX(1)))*MultiplierX inf]);
            end
            if obj.CPFlag.CNN == 1
                plot(obj.CP_CNN(k,1)*MultiplierX, obj.CP_CNN(k,2)*MultiplierY,'O',...
                    'LineWidth',1.5,...
                    'MarkerSize',7,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','y');
                Legends{end+1} = 'CP CNN';
                xlim([(AppX(1)+ZoomMult*(obj.CP_CNN(k,1) - AppX(1)))*MultiplierX inf]);
            end
            if obj.CPFlag.CNNZoom == 1
                plot(obj.CP_CNNZoom(k,1)*MultiplierX, obj.CP_CNNZoom(k,2)*MultiplierY,'O',...
                    'LineWidth',1.5,...
                    'MarkerSize',7,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor',[0.8500 0.3250 0.0980]);
                Legends{end+1} = 'CP CNNZoom';
                xlim([(AppX(1)+ZoomMult*(obj.CP_CNNZoom(k,1) - AppX(1)))*MultiplierX inf]);
            end
            if obj.CPFlag.CNNZoomSweep == 1
                plot(obj.CP_CNNZoomSweep(k,1)*MultiplierX, obj.CP_CNNZoomSweep(k,2)*MultiplierY,'gs',...
                    'LineWidth',1.5,...
                    'MarkerSize',7,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor',[0.8500 0.3250 0.0980]);
                Legends{end+1} = 'CP CNNZoomSweep';
                xlim([(AppX(1)+ZoomMult*(obj.CP_CNNZoomSweep(k,1) - AppX(1)))*MultiplierX inf]);
            end
            if obj.CPFlag.Old == 1
                plot(obj.CP_Old(k,1)*MultiplierX, obj.CP_Old(k,2)*MultiplierY,'O',...
                    'LineWidth',1.5,...
                    'MarkerSize',7,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor',[0.9290 0.6940 0.1250]);
                Legends{end+1} = 'CP SD6';
                xlim([(AppX(1)+ZoomMult*(obj.CP_Old(k,1) - AppX(1)))*MultiplierX inf]);
            end
            if obj.CPFlag.Dropout == 1
                X = obj.CP_Dropout(k,1)*MultiplierX;
                Y = obj.CP_Dropout(k,2)*MultiplierY;
                XSTD = std(obj.CP_MonteCarlo(:,1,k))*MultiplierX;
                YSTD = std(obj.CP_MonteCarlo(:,2,k))*MultiplierY;
                plot(X, Y,'O',...
                    'LineWidth',1.5,...
                    'MarkerSize',7,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor',[0.4940 0.1840 0.5560]);
                plot([X-XSTD X+XSTD],[Y Y],'-+',[X X],[Y-YSTD Y+YSTD],'-+',...
                    'LineWidth',1.5,...
                    'Color',[0.4940 0.1840 0.5560]);
                Legends{end+1} = 'CP Dropout';
                Legends{end+1} = sprintf('Dropout Uncertainty %.2f nm',obj.CP_MonteCarlo_STD(k)*MultiplierX);
                xlim([(AppX(1)+ZoomMult*(obj.CP_Dropout(k,1) - AppX(1)))*MultiplierX inf]);
            end
            
            legend(Legends,...
                'Location','northwest');
            xlabel(sprintf('Cantilever Head Height [%s]',UnitX));
            ylabel(sprintf('vDeflection-Force [%s]',UnitY));
            grid on
            grid minor
            %             dim = [0.4 0.3 0.6 0.6];
            %             str = {'     Network Uncertainty',sprintf('Standard Deviation = %.2f nm',obj.CP_MonteCarlo_STD(k)*1e9)};
            %             annotation('textbox',dim,'String',str,'FitBoxToText','on');
            
            subplot(2,1,2)
            I = obj.get_channel('Processed');
            if isempty(I)
                I = obj.HeightMap;
            end
            imshow(I.Image,[],'Colormap',AFMImage.define_afm_color_map)
            axis on
            hold on;
            try
                for i=1:obj.NumPixelsX
                    plot((obj.List2Map(obj.RectApexIndex(i),2)),...
                        (obj.List2Map(obj.RectApexIndex(i),1)),...
                        'g+', 'MarkerSize', 10, 'LineWidth', 2);
                    %                 plot((obj.List2Map(obj.ApexIndex(i),2)-1/2)*1024/obj.NumPixelsY,...
                    %                     (obj.List2Map(obj.ApexIndex(i),1)-1/2)*1024/obj.NumPixelsX,...
                    %                     'g+', 'MarkerSize', 10, 'LineWidth', 1);
                end
            catch
            end
            plot(obj.List2Map(k,2),obj.List2Map(k,1), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
            title('Choose a different point...')
            NextK = drawpoint();
            MapPos1 = round(NextK.Position(2));
            MapPos2 = round(NextK.Position(1));
            if MapPos1 == 0
                MapPos1 = 1;
            elseif MapPos1 > size(obj.Map2List,1)
                MapPos1 = size(obj.Map2List,1);
            end
            if MapPos2 == 0
                MapPos2 = 1;
            elseif MapPos2 > size(obj.Map2List,2)
                MapPos2 = size(obj.Map2List,2);
            end
            k = obj.Map2List(MapPos1,MapPos2);
            try
                hold off
                cla(fig.Children(1))
                cla(fig.Children(2))
                cla(fig.Children(3))
                obj.show_force_curve(ZoomMult,k,fig);
            catch
            end
        end
        
        function Fig = plotting_show_analyzed_fibril(obj)
            T = sprintf('Height Map of %s\nwith chosen indentation points',obj.Name);
            Fig = figure('Name',T,'Units','normalized','Color','w','Position',[0.5 0.1 0.5 0.8]);
            
            subplot(2,2,1)
            I = imresize(obj.HeightMap(:,:,1).*1e9,[1024 1024]);
            %             I = (I*range(obj.HeightMap(:,:,1),'all') + min(obj.HeightMap(:,:,1),[],'all'))*1e9;
            imshow(I,[min(I,[],'all') max(I,[],'all')],'Colormap',hot)
            title(sprintf('%s Plane Fitted Height',obj.Name))
            c1 = colorbar;
            c1.Label.String = 'Height [nm]';
            hold on;
            try
                for i=1:obj.NumPixelsX
                    plot((obj.List2Map(obj.RectApexIndex(i),2)-1/2)*1024/obj.NumPixelsY,...
                        (obj.List2Map(obj.RectApexIndex(i),1)-1/2)*1024/obj.NumPixelsX,...
                        'g+', 'MarkerSize', 10, 'LineWidth', 2);
                    %                 plot((obj.List2Map(obj.ApexIndex(i),2)-1/2)*1024/obj.NumPixelsY,...
                    %                     (obj.List2Map(obj.ApexIndex(i),1)-1/2)*1024/obj.NumPixelsX,...
                    %                     'g+', 'MarkerSize', 10, 'LineWidth', 1);
                end
            catch
            end
            title(T);
            
            subplot(2,2,3)
            Masks = imresize(obj.FibMask.*obj.ExclMask,[1024 1024]);
            imshow(Masks)
            title('Fibril mask with excluded areas')
            
            subplot(2,2,2)
            surf(imrotate(obj.HeightMap',90),'LineStyle','none','FaceLighting','gouraud','FaceColor','interp')
            axis manual
            light('Style','local')
            title(sprintf('Fibril Diameter = %.2fnm',obj.FibDiam*1e9))
            
            subplot(2,2,4)
            try
                EM = imresize(mat2gray(log(obj.EModMapOliverPharr(:,:,1))),[1024 1024]);
            catch
                try
                    EM = imresize(mat2gray(log(obj.EModMapHertz(:,:,1))),[1024 1024]);
                catch
                    EM = zeros(1024,1024);
                end
            end
            imshow(EM)
            title('Indentation Modulus Map')
        end
        
        function plotting_show_height_map(obj)
            T = sprintf('Height Map of %s',obj.Name);
            Fig = figure('Name',T,'Units','normalized','Color','w','Position',[0.5 0.1 0.5 0.8]);
            
            subplot(2,1,1)
            I = obj.HeightMap(:,:,1).*1e9;
            %             I = (I*range(obj.HeightMap(:,:,1),'all') + min(obj.HeightMap(:,:,1),[],'all'))*1e9;
            imshow(I,[min(I,[],'all') max(I,[],'all')],'Colormap',hot)
            title(sprintf('%s Plane Fitted Height',obj.Name))
            c1 = colorbar;
            c1.Label.String = 'Height [nm]';
            
            subplot(2,1,2)
            I = imresize(obj.HeightMap(:,:,1).*1e9,[1024 1024]);
            %             I = (I*range(obj.HeightMap(:,:,1),'all') + min(obj.HeightMap(:,:,1),[],'all'))*1e9;
            imshow(I,[min(I,[],'all') max(I,[],'all')],'Colormap',hot)
            title(sprintf('%s Plane Fitted Height Resized',obj.Name))
            c1 = colorbar;
            c1.Label.String = 'Height [nm]';
            
        end
        
        function plotting_show_e_mod_map(obj)
            Title = sprintf('Indentation Modulus Map of %s',obj.Name);
            figure('Name',Title);
            subplot(2,2,1)
            Lower = mean(obj.EModMapHertz,'all','omitnan')-1.5*std(obj.EModMapHertz,0,'all','omitnan');
            Upper = mean(obj.EModMapHertz,'all','omitnan')+1.5*std(obj.EModMapHertz,0,'all','omitnan');
            I = imresize(obj.EModMapHertz,[1024 1024]);
            imshow(I,[Lower Upper],'Colormap',copper);
            colorbar
            title('EMod - Hertz-Sneddon Method')
            subplot(2,2,2)
            surf(obj.EModMapHertz(:,:,1),'LineStyle','none','FaceLighting','gouraud','FaceColor','interp')
            light('Style','local')
            subplot(2,2,3)
            Lower = mean(obj.EModMapOliverPharr,'all','omitnan')-1.5*std(obj.EModMapOliverPharr,0,'all','omitnan');
            Upper = mean(obj.EModMapOliverPharr,'all','omitnan')+1.5*std(obj.EModMapOliverPharr,0,'all','omitnan');
            I = imresize(obj.EModMapOliverPharr,[1024 1024]);
            imshow(I,[Lower Upper]);
            colorbar
            title('Emod Oliver-Pharr Method')
            subplot(2,2,4)
            surf(obj.EModMapOliverPharr,'LineStyle','none','FaceLighting','gouraud','FaceColor','interp')
            light('Style','local')
        end
        
        function plotting_quality_control_oliver_pharr_fibril(obj,PauseTime)
            % shows some relevant plots for the E-Mod calculation
            % rectified apex force curves
            
            if nargin < 2
                PauseTime = 1.5;
            end
            
            T = sprintf('Quality control of fibril %s',obj.Name);
            Fig = figure('Name',T,'Units','normalized','Color','w','Position',[0.1 0.1 0.8 0.8]);
            m = 1;
            while 1==1
                try
                    figure(Fig);
                catch
                    return
                end
                k = obj.RectApexIndex(m);
                subplot(2,3,1)
                I = obj.get_channel('Processed');
                if isempty(I)
                    I = obj.HeightMap;
                end
                imshow(imresize(I.Image,[1024 1024]),[],'Colormap',AFMImage.define_afm_color_map)
                hold on;
                for i=1:obj.NumPixelsX
                    if obj.RectApexIndex(i)==k
                        plot((obj.List2Map(obj.RectApexIndex(i),2)-1/2)*1024/obj.NumPixelsY,...
                            (obj.List2Map(obj.RectApexIndex(i),1)-1/2)*1024/obj.NumPixelsX,...
                            'g*', 'MarkerSize', 10, 'LineWidth', 2);
                    else
                        plot((obj.List2Map(obj.RectApexIndex(i),2)-1/2)*1024/obj.NumPixelsY,...
                            (obj.List2Map(obj.RectApexIndex(i),1)-1/2)*1024/obj.NumPixelsX,...
                            'r+', 'MarkerSize', 5, 'LineWidth', 1);
                    end
                end
                title('Height Map with Apex Points');
                
                [AppY,AppX] = obj.get_force_curve_data(k,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                [RetY,RetX] = obj.get_force_curve_data(k,'AppRetSwitch',1,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                
                [MultiplierX,UnitX,~] = AFMImage.parse_unit_scale(range(RetX),'m',10);
                [MultiplierY,UnitY,~] = AFMImage.parse_unit_scale(range(RetY./obj.SpringConstant),'m',5);
                [MultiplierPaY,UnitPaY,~] = AFMImage.parse_unit_scale(obj.EModOliverPharr(k),'Pa',5);
                HHApp = AppX.*MultiplierX;
                App = AppY/obj.SpringConstant.*MultiplierY;
                HHRet = RetX.*MultiplierX;
                Ret = RetY/obj.SpringConstant.*MultiplierY;
                subplot(2,3,2)
                plot(HHApp,App,...
                    HHRet,Ret)
                xlim([min(HHApp)+range(HHApp)/2 ...
                    max(HHApp)+range(HHApp)*0.1])
                title(sprintf('Apparent Indentation Modulus = %.2f %s',obj.EModOliverPharr(k).*MultiplierPaY,UnitPaY))
                legend('Approach','Retract','Location','northwest')
                xlabel(sprintf('Head Height [%s]',UnitX))
                ylabel(sprintf('vDeflection [%s]',UnitY))
                drawpoint('Position',[obj.CP(k,1).*MultiplierX obj.CP(k,2).*MultiplierY]);
                
                if m == 1
                    subplot(2,3,3)
                    boxplot(obj.EModOliverPharr(obj.RectApexIndex));
                    xticklabels(obj.Name)
                    title(sprintf('mean = %.2f MPa\nmedian = %.2f MPa\nstd = %.3f MPa',...
                        mean(obj.EModOliverPharr(obj.RectApexIndex),'omitnan')*1e-6,...
                        median(obj.EModOliverPharr(obj.RectApexIndex),'omitnan')*1e-6,...
                        std(obj.EModOliverPharr(obj.RectApexIndex),'omitnan')*1e-6));
                end
                
                subplot(2,3,4)
                plot(0:obj.NumPixelsX+1,obj.RefSlope*ones(obj.NumPixelsX+2,1))
                ylim([0 1.3])
                xlim([0 obj.NumPixelsX+1])
                hold on
                plot(1:obj.NumPixelsX,obj.DZslope(obj.RectApexIndex),'bO')
                plot(m,obj.DZslope(obj.RectApexIndex(m)),'rO','MarkerFaceColor','r')
                xlabel('Index')
                ylabel('DZ-Slope')
                legend(sprintf('Glass Reference Slope = %.3f',obj.RefSlope))
                hold off
                
                subplot(2,3,5)
                plot(obj.IndentationDepthOliverPharr(obj.RectApexIndex)*1e9,...
                    obj.EModOliverPharr(obj.RectApexIndex)*1e-6,'bO')
                hold on
                plot(obj.IndentationDepthOliverPharr(obj.RectApexIndex(m))*1e9,...
                    obj.EModOliverPharr(obj.RectApexIndex(m))*1e-6,'rO','MarkerFaceColor','r')
                xlabel('Indentation Depth [nm]')
                ylabel('Apparent Indentation Modulus [MPa]')
                hold off
                
                subplot(2,3,6)
                plot(obj.IndentationDepthOliverPharr(obj.RectApexIndex)*1e9,...
                    obj.IndentArea(obj.RectApexIndex),'bO','MarkerSize',10,'MarkerFaceColor','b')
                hold on
                plot(obj.IndentationDepthOliverPharr(obj.RectApexIndex(m))*1e9,...
                    obj.IndentArea(obj.RectApexIndex(m)),'rO','MarkerSize',10,'MarkerFaceColor','r')
                Xmax = round(max(obj.IndentationDepthOliverPharr(obj.RectApexIndex))*1e9+5);
                xlim([0 Xmax])
                plot(1:Xmax,obj.ProjTipArea(1:Xmax),'Color','black')
                xlabel('Indentation Depth [nm]')
                ylabel('Projected Area [m^2]')
                hold off
                
                pause(PauseTime)
                if m<obj.NumPixelsX
                    m = m + 1;
                else
                    m = 1;
                end
            end
        end
        
        function plotting_quality_control_hertz_sneddon_fibril(obj,PauseTime)
            % shows some relevant plots for the E-Mod calculation
            % rectified apex force curves
            
            if nargin < 2
                PauseTime = 1.5;
            end
            
            T = sprintf('Quality control of fibril %s',obj.Name);
            Fig = figure('Name',T,'Units','normalized','Color','w','Position',[0.1 0.1 0.8 0.8]);
            m = 1;
            while 1==1
                try
                    figure(Fig);
                catch
                    return
                end
                k = obj.RectApexIndex(m);
                
                subplot(2,2,1)
                I = obj.get_channel('Processed');
                if isempty(I)
                    I = obj.HeightMap;
                end
                imshow(imresize(I.Image,[1024 1024]),[],'Colormap',AFMImage.define_afm_color_map)
                hold on;
                for i=1:obj.NumPixelsX
                    if obj.RectApexIndex(i)==k
                        plot((obj.List2Map(obj.RectApexIndex(i),2)-1/2)*1024/obj.NumPixelsY,...
                            (obj.List2Map(obj.RectApexIndex(i),1)-1/2)*1024/obj.NumPixelsX,...
                            'g*', 'MarkerSize', 10, 'LineWidth', 2);
                    else
                        plot((obj.List2Map(obj.RectApexIndex(i),2)-1/2)*1024/obj.NumPixelsY,...
                            (obj.List2Map(obj.RectApexIndex(i),1)-1/2)*1024/obj.NumPixelsX,...
                            'r+', 'MarkerSize', 5, 'LineWidth', 1);
                    end
                end
                title('Height Map with Apex Points');
                
                [AppY,AppX] = obj.get_force_curve_data(k,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',1,...
                    'Sensitivity','original','Unit','N');
                [RetY,RetX] = obj.get_force_curve_data(k,'AppRetSwitch',1,...
                    'BaselineCorrection',1,'TipHeightCorrection',1,...
                    'Sensitivity','original','Unit','N');
                
                subplot(2,2,2)
                
                % Determine X-Range for HertzModel
                X = AppX - obj.CP(k,1);
                X(X<0) = [];
                HertzModelX = 0:range(X)/100:2*max(X);
                
                try
                    FitModel = obj.HertzFit{m};
                catch
                    ft = fittype(obj.HertzFitType);
                    if length(obj.HertzFitValues{k}) == 1
                        FitModel = cfit(ft,obj.HertzFitValues{k}(1));
                    elseif length(obj.HertzFitValues{k}) == 2
                        FitModel = cfit(ft,obj.HertzFitValues{k}(1),obj.HertzFitValues{k}(2));
                    elseif length(obj.HertzFitValues{k}) == 3
                        FitModel = cfit(ft,obj.HertzFitValues{k}(1),obj.HertzFitValues{k}(2),...
                            obj.HertzFitValues{k}(3));
                    elseif length(obj.HertzFitValues{k}) == 4
                        FitModel = cfit(ft,obj.HertzFitValues{k}(1),obj.HertzFitValues{k}(2),...
                            obj.HertzFitValues{k}(3),obj.HertzFitValues{k}(4));
                    elseif length(obj.HertzFitValues{k}) == 5
                        FitModel = cfit(ft,obj.HertzFitValues{k}(1),obj.HertzFitValues{k}(2),...
                            obj.HertzFitValues{k}(3),obj.HertzFitValues{k}(4),...
                            obj.HertzFitValues{k}(5));
                    elseif length(obj.HertzFitValues{k}) == 6
                        FitModel = cfit(ft,obj.HertzFitValues{k}(1),obj.HertzFitValues{k}(2),...
                            obj.HertzFitValues{k}(3),obj.HertzFitValues{k}(4),...
                            obj.HertzFitValues{k}(5),obj.HertzFitValues{k}(6));
                    end
                end
                
                HertzModelY = feval(FitModel,HertzModelX);
                
                
                plot(HertzModelX(HertzModelY<=max(AppY - obj.CP(k,2))),HertzModelY(HertzModelY<=max(AppY - obj.CP(k,2))),...
                    AppX - obj.CP(k,1),AppY - obj.CP(k,2),...
                    RetX - obj.CP(k,1),RetY - obj.CP(k,2))
                xlim([min(AppX - obj.CP(k,1))+range(AppX - obj.CP(k,1))/2 ...
                    max(AppX - obj.CP(k,1))+range(AppX - obj.CP(k,1))*0.1])
                title(sprintf('Apparent Indentation Modulus = %.2f MPa',obj.EModHertz(k)*1e-6))
                legend('Hertz Fit','Approach','Retract','Location','northwest')
                xlabel('Cantilever Tip Height [m]')
                ylabel('Force [N]')
                drawpoint('Position',[0 0]);
                
                if m == 1
                    subplot(2,2,3)
                    boxplot(obj.EModHertz(obj.RectApexIndex));
                    xticklabels(obj.Name)
                    title(sprintf('mean = %.2f MPa\nmedian = %.2f MPa\nstd = %.3f MPa',...
                        mean(obj.EModHertz(obj.RectApexIndex),'omitnan')*1e-6,...
                        median(obj.EModHertz(obj.RectApexIndex),'omitnan')*1e-6,...
                        std(obj.EModHertz(obj.RectApexIndex),'omitnan')*1e-6));
                end
                
                
                subplot(2,2,4)
                plot((obj.IndentationDepth(obj.RectApexIndex))*1e9,...
                    obj.EModHertz(obj.RectApexIndex)*1e-6,'bO')
                hold on
                plot((obj.IndentationDepth(obj.RectApexIndex(m)))*1e9,...
                    obj.EModHertz(obj.RectApexIndex(m))*1e-6,'rO','MarkerFaceColor','r')
                xlabel('Indentation Depth [nm]')
                ylabel('Apparent Indentation Modulus [MPa]')
                hold off
                
                pause(PauseTime)
                if m<obj.NumPixelsX
                    m = m + 1;
                else
                    m = 1;
                end
            end
        end
        
        function plotting_quality_control_oliver_pharr(obj,PauseTime)
            % shows some relevant plots for the E-Mod calculation
            
            if nargin < 2
                PauseTime = 1.5;
            end
            
            T = sprintf('Quality control of force map %s',obj.Name);
            Fig = figure('Name',T,'Units','normalized','Color','w','Position',[0.1 0.1 0.8 0.8]);
            m = 1;
            while 1==1
                try
                    figure(Fig);
                catch
                    return
                end
                subplot(2,3,1)
                I = obj.get_channel('Processed');
                if isempty(I)
                    I = obj.HeightMap;
                end
                imshow(imresize(I.Image,[1024 1024]),[],'Colormap',AFMImage.define_afm_color_map)
                hold on;
                
                plot((obj.List2Map(m,2)-1/2)*1024/obj.NumPixelsY,...
                    (obj.List2Map(m,1)-1/2)*1024/obj.NumPixelsX,...
                    'g*', 'MarkerSize', 10, 'LineWidth', 2);
                
                title('Height Map');
                
                subplot(2,3,2)
                
                [App,HHApp] = obj.get_force_curve_data(m,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                [Ret,HHRet] = obj.get_force_curve_data(m,'AppRetSwitch',1,...
                    'BaselineCorrection',1,'TipHeightCorrection',0,...
                    'Sensitivity','original','Unit','N');
                
                [MultiplierX,UnitX,~] = AFMImage.parse_unit_scale(range(HHRet),'m',10);
                [MultiplierY,UnitY,~] = AFMImage.parse_unit_scale(range(Ret./obj.SpringConstant),'m',5);
                [MultiplierPaY,UnitPaY,~] = AFMImage.parse_unit_scale(obj.EModOliverPharr(m),'Pa',5);
                HHApp = HHApp.*MultiplierX;
                App = App/obj.SpringConstant.*MultiplierY;
                HHRet = HHRet.*MultiplierX;
                Ret = Ret/obj.SpringConstant.*MultiplierY;
                subplot(2,3,2)
                plot(HHApp,App,...
                    HHRet,Ret)
                xlim([min(HHApp)+range(HHApp)/2 ...
                    max(HHApp)+range(HHApp)*0.1])
                title(sprintf('Apparent Indentation Modulus = %.2f %s',obj.EModOliverPharr(m).*MultiplierPaY,UnitPaY))
                legend('Approach','Retract','Location','northwest')
                xlabel(sprintf('Head Height [%s]',UnitX))
                ylabel(sprintf('vDeflection [%s]',UnitY))
                drawpoint('Position',[obj.CP(m,1).*MultiplierX obj.CP(m,2).*MultiplierY]);
                
                if m == 1
                    subplot(2,3,3)
                    boxplot(obj.EModOliverPharr);
                    xticklabels(obj.Name)
                    title(sprintf('mean = %.2f MPa\nmedian = %.2f MPa\nstd = %.3f MPa',...
                        mean(obj.EModOliverPharr,'omitnan')*1e-6,...
                        median(obj.EModOliverPharr,'omitnan')*1e-6,...
                        std(obj.EModOliverPharr,'omitnan')*1e-6));
                end
                
                subplot(2,3,4)
                plot(0:obj.NCurves+1,obj.RefSlope*ones(obj.NCurves+2,1))
                ylim([0 1.3])
                xlim([0 obj.NCurves+1])
                hold on
                plot(1:obj.NCurves,obj.DZslope,'bO')
                plot(m,obj.DZslope(m),'rO','MarkerFaceColor','r')
                xlabel('Index')
                ylabel('DZ-Slope')
                legend(sprintf('Glass Reference Slope = %.3f',obj.RefSlope))
                hold off
                
                subplot(2,3,5)
                plot(obj.IndentationDepthOliverPharr*1e9,...
                    obj.EModOliverPharr*1e-6,'bO')
                hold on
                plot(obj.IndentationDepthOliverPharr(m)*1e9,...
                    obj.EModOliverPharr(m)*1e-6,'rO','MarkerFaceColor','r')
                xlabel('Indentation Depth [nm]')
                ylabel('Apparent Indentation Modulus [MPa]')
                hold off
                
                subplot(2,3,6)
                plot(obj.IndentationDepthOliverPharr*1e9,...
                    obj.IndentArea,'bO','MarkerSize',10,'MarkerFaceColor','b')
                hold on
                plot(obj.IndentationDepthOliverPharr(m)*1e9,...
                    obj.IndentArea(m),'rO','MarkerSize',10,'MarkerFaceColor','r')
                Xmax = round(max(obj.IndentationDepthOliverPharr)*1e9+5);
                XmaxAlt = length(obj.ProjTipArea);
                Xmax = min([Xmax XmaxAlt]);
                xlim([0 Xmax])
                plot(1:Xmax,obj.ProjTipArea(1:Xmax),'Color','black')
                xlabel('Indentation Depth [nm]')
                ylabel('Projected Area [m^2]')
                hold off
                
                pause(PauseTime)
                if m<obj.NCurves
                    m = m + 1;
                else
                    m = 1;
                end
            end
        end
        
        function plotting_quality_control_hertz_sneddon(obj,PauseTime)
            % shows some relevant plots for the E-Mod calculation
            % rectified apex force curves
            
            if nargin < 2
                PauseTime = 1.5;
            end
            
            T = sprintf('Quality control of force map %s',obj.Name);
            Fig = figure('Name',T,'Units','normalized','Color','w','Position',[0.1 0.1 0.8 0.8]);
            m = 1;
            while 1==1
                try
                    figure(Fig);
                catch
                    return
                end
                subplot(2,3,1)
                I = obj.get_channel('Processed');
                if isempty(I)
                    I = obj.HeightMap;
                end
                imshow(imresize(I.Image,[1024 1024]),[],'Colormap',AFMImage.define_afm_color_map)
                hold on;
                
                plot((obj.List2Map(m,2)-1/2)*1024/obj.NumPixelsY,...
                    (obj.List2Map(m,1)-1/2)*1024/obj.NumPixelsX,...
                    'g*', 'MarkerSize', 10, 'LineWidth', 2);
                
                title('Height Map');
                
                subplot(2,2,2)
                
                [AppY,AppX] = obj.get_force_curve_data(m,'AppRetSwitch',0,...
                    'BaselineCorrection',1,'TipHeightCorrection',1,...
                    'Sensitivity','original','Unit','N');
                [RetY,RetX] = obj.get_force_curve_data(m,'AppRetSwitch',1,...
                    'BaselineCorrection',1,'TipHeightCorrection',1,...
                    'Sensitivity','original','Unit','N');
                
                % Determine X-Range for HertzModel
                X = AppX - obj.CP(m,1);
                X(X<0) = [];
                HertzModelX = 0:range(X)/100:2*max(X);
                try
                    FitModel = obj.HertzFit{m};
                catch
                    ft = fittype(obj.HertzFitType);
                    if length(obj.HertzFitValues{m}) == 1
                        FitModel = cfit(ft,obj.HertzFitValues{m}(1));
                    elseif length(obj.HertzFitValues{m}) == 2
                        FitModel = cfit(ft,obj.HertzFitValues{m}(1),obj.HertzFitValues{m}(2));
                    elseif length(obj.HertzFitValues{m}) == 3
                        FitModel = cfit(ft,obj.HertzFitValues{m}(1),obj.HertzFitValues{m}(2),...
                            obj.HertzFitValues{m}(3));
                    elseif length(obj.HertzFitValues{m}) == 4
                        FitModel = cfit(ft,obj.HertzFitValues{m}(1),obj.HertzFitValues{m}(2),...
                            obj.HertzFitValues{m}(3),obj.HertzFitValues{m}(4));
                    elseif length(obj.HertzFitValues{m}) == 5
                        FitModel = cfit(ft,obj.HertzFitValues{m}(1),obj.HertzFitValues{m}(2),...
                            obj.HertzFitValues{m}(3),obj.HertzFitValues{m}(4),...
                            obj.HertzFitValues{m}(5));
                    elseif length(obj.HertzFitValues{m}) == 6
                        FitModel = cfit(ft,obj.HertzFitValues{m}(1),obj.HertzFitValues{m}(2),...
                            obj.HertzFitValues{m}(3),obj.HertzFitValues{m}(4),...
                            obj.HertzFitValues{m}(5),obj.HertzFitValues{m}(6));
                    end
                end
                
                HertzModelY = feval(FitModel,HertzModelX);
                
                plot(HertzModelX(HertzModelY<=max(AppY - obj.CP(m,2))),HertzModelY(HertzModelY<=max(AppY - obj.CP(m,2))),...
                    AppX - obj.CP(m,1),AppY - obj.CP(m,2),...
                    RetX - obj.CP(m,1),RetY - obj.CP(m,2))
                xlim([min(AppX - obj.CP(m,1))+range(AppX - obj.CP(m,1))/2 ...
                    max(AppX - obj.CP(m,1))+range(AppX - obj.CP(m,1))*0.1])
                title(sprintf('Apparent Indentation Modulus = %.2f MPa',obj.EModHertz(m)*1e-6))
                legend('Hertz Fit','Approach','Retract','Location','northwest')
                xlabel('Cantilever Tip Height [m]')
                ylabel('Force [N]')
                drawpoint('Position',[0 0]);
                
                if m == 1
                    subplot(2,2,3)
                    boxplot(obj.EModHertz);
                    xticklabels(obj.Name)
                    title(sprintf('mean = %.2f MPa\nmedian = %.2f MPa\nstd = %.3f MPa',...
                        mean(obj.EModHertz,'omitnan')*1e-6,...
                        median(obj.EModHertz,'omitnan')*1e-6,...
                        std(obj.EModHertz,'omitnan')*1e-6));
                end
                
                
                subplot(2,2,4)
                if obj.NCurves <= 512
                    plot((obj.IndentationDepth)*1e9,...
                        obj.EModHertz*1e-6,'bO')
                    hold on
                    plot((obj.IndentationDepth(m))*1e9,...
                        obj.EModHertz(m)*1e-6,'rO','MarkerFaceColor','r')
                    xlabel('Indentation Depth [nm]')
                    ylabel('Apparent Indentation Modulus [MPa]')
                else
                    histogram((obj.IndentationDepth)*1e9)
                    xlabel('Indentation Depth [nm]')
                end    
                hold off
                    
                pause(PauseTime)
                if m<obj.NCurves
                    m = m + 1;
                else
                    m = 1;
                end
            end
        end
        
        function plotting_compare_hertz_oliver_fibril(obj)
            % Boxplots comparing the two methods of elastic modulus
            % calculation
            figure('Name',obj.Name)
            % Oliver Pharr Method
            subplot(2,2,1)
            boxplot(obj.EModOliverPharr(obj.RectApexIndex));
            xticklabels('Oliver-Pharr Method')
            title(sprintf('mean = %.2f MPa\nmedian = %.2f MPa\nstd = %.3f MPa',...
                mean(obj.EModOliverPharr(obj.RectApexIndex))*1e-6,...
                median(obj.EModOliverPharr(obj.RectApexIndex))*1e-6,...
                std(obj.EModOliverPharr(obj.RectApexIndex))*1e-6));
            
            % Hertz Sneddon Method
            subplot(2,2,2)
            boxplot(obj.EModHertz(obj.RectApexIndex));
            xticklabels('Hertz-Sneddon Method')
            title(sprintf('mean = %.2f MPa\nmedian = %.2f MPa\nstd = %.3f MPa',...
                mean(obj.EModHertz(obj.RectApexIndex))*1e-6,...
                median(obj.EModHertz(obj.RectApexIndex))*1e-6,...
                std(obj.EModHertz(obj.RectApexIndex))*1e-6));
            
            % Comparing all apex force curves
            subplot(2,2,[3 4])
            plot(1:obj.NumPixelsX,obj.EModOliverPharr(obj.RectApexIndex)*1e-6,'bO')
            hold on
            plot(1:obj.NumPixelsX,obj.EModHertz(obj.RectApexIndex)*1e-6,'rO')
            xlabel('Index')
            ylabel('Apparent Indentation Modulus [MPa]')
        end
        
        function plotting_compare_loglog_with_shifted_cp(obj,CurveI)
            
            if nargin < 2
                CurveI = randi(obj.NCurves);
            end
            
            figure('Name',['compare loglog with shifted cp. curve nr.: ' num2str(CurveI)],...
                'Color','w',...
                'Position',[100 100 1480 850]);
            
            [Force,TH] = obj.get_force_curve_data(CurveI,'AppRetSwitch',0,...
                'BaselineCorrection',1,'TipHeightCorrection',1,...
                'Sensitivity','original','Unit','N');
            TH = TH - obj.CP_CNN(CurveI,1);
            TH = TH.*1e6;
            Force = Force.*1e9;
            VarVec = [-1:.2:1].*max(TH)*0.5;
            for i=1:length(VarVec)
                ax1 = subplot(2,1,1);
                plot(TH - VarVec(i),Force)
                hold on
                xlabel('Tip Height [\mum]')
                ylabel('Force [nN]')
                ax2 = subplot(2,1,2);
                loglog(TH - VarVec(i),Force)
                hold on
                if i==length(VarVec)
                    loglog(TH,10*TH.^(3/2),'Color','k','LineWidth',1.5)
                    hold on
                    TextBoxI = find(TH>0);
                    text(TH(TextBoxI(1)),10*TH(TextBoxI(1)).^(3/2),'Hertzian contact slope: 3/2','FontSize',22,'EdgeColor','k')
                end
                xlabel('Tip Height [\mum]')
                ylabel('Force [nN]')
            end
            FS = 18;
            ax1.FontSize = FS;
            ax2.FontSize = FS;
        end
        
    end
end