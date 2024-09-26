classdef ForceMap < matlab.mixin.Copyable & matlab.mixin.SetGet & handle & AFMBaseClass
    % The force map class represents a single jpk force map file and
    % contains all necessary functions to process the forcecurves.
    % General naming convention for this class is:
    % -class properties get capital letters and generally no space between
    % words
    % -class methods are in lowercase letters and get underscores between
    % words
    %%%%%%%%%%%%%%%%%%%%%%%%%%DISCLAIMER%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        
        Date            % date when the force map was detected
        Time            % time when the force map was detected
        FileVersion     % Version of jpk-force-map file
        Folder          % location of the .csv files of the force map
        DataStoreFolder % 
        BigDataFlag     % If true, unpack data container into Experiment subfolder
                        % and always load force volume data from there
        FractionOfMaxRAM = 1/5 % Specifies how much of MaxRAM space can be taken for certain partitioned calculations 
        NCurves         % number of curves on the force map
        MaxPointsPerCurve
        HoldingTime = 0
        NumSegments
        ExtendTime
        RetractTime
        ExtendZLength
        RetractZLength
        ExtendVelocity        % Approach  velocity as defined in the force map settings
        RetractVelocity
        Sensitivity
        RefSlopeCorrectedSensitivity
        SpringConstant
        Setpoint
        DBanding        % Fourieranalysis-based estimate of DBanding perdiod (only available with sufficient resolution)
        RefSlope        % Refernce slope as determined from the upper curve slope from data from very hard      
        % surfaces (mica,glass), either from glass parts beneath the specimen or from
        % separate reference force maps
        PixApp          % maximum number of measured points during approach
        PixRet          % maximum number of measured points during retraction
        SelectedCurves  % logical vector of length NCurves with 0s for excluded and 1s for included curves. gets initialized with ones
        StashedSelectedCurves % In some calculations, SelectedCurves is changed temporarily inside a method.
                              % Actual Selction is stored and then restored from here
        CorruptedCurves % Curves that cant be loaded
        TipRadius = 2  % (nominal, if not otherwise calculated) tip radius in nm for the chosen type of tip model
        PoissonR = 0.5  % standard Poisson ratio for most mechanical models
        Medium
        FibrilFlag
        FibPot
    end
    properties
        % Curve data Properties
        
        App = {}        % approach force data in Newton
        Ret = {}        % retraction force data in Newton
        HHApp = {}      % capacitive-sensor-height approach data in meters
        HHRet = {}      % capacitive-sensor-height retract data in meters
        HHType          % Type of Head Height.(Default is capacitiveSensorHeight; switches to measuredHeight, if default does'nt exist)
        THApp = {}      % vertical tip height approach data in meters
        THRet = {}      % vertical tip height retract data in meters
        BasedApp = {}   % approach force data with subtracted base line and tilt in Newton
        BasedRet = {}   % retraction force data with subtracted base line and tilt in Newton
        BasedRetFit = {}  
        BaseAndTiltFlag
    end
    properties
        % Properties related to Contact Point (CP) estimation
        
        RoV = {}        % ratio of variance curve for CP estiamation
        GoF = {}        % goodness of fit curve for each selected curve
        CP              % chosen contact point (CP) for every selected curve, which is then used in E module fitting
        CP_SnapIn       % Prefered CP-method for data with snap-in-effect (most data from air)
        CP_HertzFitted  % Computed from the 0-Point of the fit calculated in Hertz-Sneddon model with AllowXShift = true
        CP_RoV          % CP estimated with the ratio of variance method
        CP_GoF          % CP estimated with the goodness of fit method
        CP_Combo        % CP estimated with a combination of the GoF and RoV methods
        CPComboCurve    % combination of the various metrics for contact point estimation
        CP_CNN          % CP estimated with a conv. neural network in a single pass
        CP_CNNZoom      % CP estimated with a conv. neural network, first estimate CP, then zoom into the curve for more accurate pred
        CP_Dropout      % CP estimated with a conv. neural network in multiple Monte Carlo Dropout passes
        CP_CNNZoomSweep % CP estimated with a conv. neural network, sweep over several magnifications and take mean (or median?) estimate
        CP_MonteCarlo   % All predictions from the multiple inference steps done in
        % the monte carlo method
        CP_MonteCarlo_STD % Standard deviation of CP_MonteCarlo
        MiniBatchSize % Optimal MiniBatchSize for CNN-prediction for current system environment
        DeltaE = {}     %
        YDropPred       % Contains the Dropoutpredictions for every curve in the forcemap
        CP_Old          % contact point estimation from old script 'A_nIAFM_analysis_main'
        LoadOld         % comes from same script as CP_old
        UnloadOld       % comes from same script as CP_old
        Man_CP          % manually chosen contact point
        CP_HardSurface  % Detract cantilever deflection for CP estimation
        CPFlag          % Struct containing booleans to indicate if a certain CP-type has been estimated
        
    end
    properties
        % Properties related to topological calculations, such as mapping and masking and visualisation
        
        HeightMap       % height profile map taken from the maximum head-height from approach max(hhapp)
        EModMapHertz    % E modulus profile map. same ordering as HeightMap
        EModMapOliverPharr % """"
        FibDiam = []    % Estimated fibril diameter
        FibDiamSTD      % Estimated fibril diameter std
        FibMask         % Logical mask marking the whole fibril
        BackgroundMask  % Logical mask marking the glass/mica/even-substrate background
        ExclMask        % Manually chosen areas that are to be excluded for calculations of fibril EMod, FibDiam, DBanding etc.
        Apex            % Value of highest pixel in each profile
        RectApex        % Value of rectified apex location in each profile
        ApexIndex       % Index of highest pixel in each profile (List indexing!)
        RectApexIndex   % Index of rectified apex location in each profile (List indexing!)
        hasOverlay = 0
    end
    properties
        % Properties related to EMod calculations
        
        Basefit = {}    % fit model used for the baseline fit in the base_and_tilt method
        Baseline
        TrueZero
        EModHertz       % List of reduced smaple E-Modulus based on a the Hertz-Sneddon model
        EModOliverPharr % List of reduced sample E-Modulus based on the Oliver-Pharr method
        FibrilEModOliverPharr
        FibrilEModHertz
        HertzFit        % HertzFit model generated in the calculate_e_mod_hertz method
        SnapIn
        MaxAdhesionForce
        AdhesionEnergy
        AdhesionLength
        DissipatedEnergy
        ElasticEnergy
        PeakIndentationAngle
        IndentationDepth
        DZslope
        Stiffness
        IndentationDepthOliverPharr
        IndentArea
        ProjTipArea
        HasRefSlope
    end
    properties
        % auxiliary properties to facilitate comparing different methods of
        % CP estimation
        EModOliverPharr_CNN
        EModOliverPharr_Old
        EModOliverPharr_RoV
        EModHertz_CNN
        EModHertz_Old
        EModHertz_RoV
    end

    % SM related
    properties                 
        Linker          % Used NHS-PEG-MI Linker (short or long) 
        Substrate       % Used substrate for the measurement 
        EnvCond         % Environmental condition during the experiment
        ChipCant        % AFM-Chip number and Cantilever label
        Chipbox         % AFM-Chipbox number (in Roman numerals)
        BasedAppDataPts
        BasedRetDataPts
        BasedAppAlign={}   % capacitive-sensor-height approach data in meters with aligned datapoints to the corresponding retraction data
        BasedRetAlign={}   % capacitive-sensor-height retraction data in meters with aligned datapoints to the corresponding approach data
        xAppRetSubstracted={} 
        xOriRet         % x Origin data
        yOriRet         % y Origin data
        BasedRetCorr2   % BasedRet data corrected based on a selection of the retraction data     
        CorrBasedRet    % BasedRet data corrected based on a selection of the approach data
        CorrBasedRetFit    % BasedRet data corrected based on a selection of the approach data
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
        AdhForceMaxAppIdx   % Index of the maximum adhesion force of the approach data Ba
        AdhForceMaxRetIdx   % Index of the maximum adhesion force of the retraction data
        AdhForceUnbindingIdx
        AppIdx1             % Index closest to 200 nm distance of the approach dat
        RetIdx1             % Index closest to 50 nm distance of the retraction data
        RetIdx2             % Index closest to 300 nm distance of the retraction data
        SinFit1Coeffa         % Sinus fit coefficient a derived from the first sinoidal fit
        SinFit1Coeffb         % Sinus fit coefficient b derived from the first sinoidal fit
        SinFit1Coeffc         % Sinus fit coefficient c derived from the first sinoidal fit
        SinFit1Coeffd         % Sinus fit coefficient d derived from the first sinoidal fit
        SinFit1RSquare          % R sqare of the first sinoidal fit
        SinFit1SSE              % SSE of the first sinoidal fit
        FitxData        % x data used in the sinoidal fit
        FityData        % y data used in the sinoidal fit
        SinFit2Coeffa         % Sinus fit coefficient a derived from the second sinoidal fit
        SinFit2Coeffb         % Sinus fit coefficient b derived from the second sinoidal fit
        SinFit2Coeffc         % Sinus fit coefficient c derived from the second sinoidal fit
        SinFit2Coeffd         % Sinus fit coefficient d derived from the second sinoidal fit
        SinFit2Coeffe         % Sinus fit coefficient e derived from the second sinoidal fit
        SinFit2RSquare          % R sqare of the second sinoidal fit
        SinFit2SSE              % SSE of the second sinoidal fit
        LinFitxData         % x data used in the linear fit
        LinFityData         % y data used in the linear fit
        LinFitCoeffa        % Linear fit coefficient a derived from the linear fit
        LinFitCoeffb        % Linear fit coefficient b derived from the linear fit
        LinFitRSquare       % R sqare of the linear fit
        LinFitSSE           % SSE of the linear fit
        LinFitRetxData         % x data used in the linear fit
        LinFitRetyData         % y data used in the linear fit
        LinFitRetCoeffa        % Linear fit coefficient a derived from the linear fit
        LinFitRetCoeffb        % Linear fit coefficient b derived from the linear fit
        LinFitRetRSquare       % R sqare of the linear fit
        LinFitRetSSE           % SSE of the linear fit
        yRetFitCorr
        yAppFitMean
        yRetFitMean
        yRet2FitMean
        yRetFitStd
        SensitivityFcSlope
        SensitivityFc
        SensitivityFcLinFitRSquare
        SensitivityFcLinFitSSE
        AppCorr
        RetCorr
        BasedAppCorr     % Based data calculated using the sensitivity of the force curve
        BasedRetCorr     % Based data calculated using the sensitivity of the force curve
        THAppCorr
        THRetCorr
        CP_HardSurfaceCorr
        xDataToFit      % ONLY A TESTING PROPERTY
        yDataToFit      % ONLY A TESTING PROPERTY

    end
    % Flags
    properties
        SMFSFlag
        SMFSFlagDown
        DebugFlag
    end
    methods
        % Main methods of the class
        
        function obj = ForceMap(MapFullFile,DataFolder,TempID,BigData,FakeOpt,NSynthCurves)
            %%% Constructor of the class
            
            % Specify the folder where the files live. And import them.
            % Also get curent folder and return to it after import of
            % files.
            % Assigns the properties that can be found in the jpk-file
            % already
            
            current = what();
            
            obj.ID = TempID;
            obj.BigDataFlag = BigData;
            
            if nargin >= 5 && isequal(FakeOpt,'Dummy')
                obj.create_dummy_force_map(NSynthCurves);
            end
            
            % get OS and use appropriate fitting system command
            obj.check_for_new_host
            
            % Unpack jpk-force-map with 7zip call to the terminal
            obj.unpack_jpk_force_map(MapFullFile,DataFolder);
            
            Index = regexp(obj.ID,'(?<=\-).','all');
            LoadMessage = sprintf('loading data into ForceMap Nr.%s',obj.ID(Index(end):end));
            disp(LoadMessage)
            
            % reading header properties into object
            obj.read_in_header_properties;
            obj.SelectedCurves = true(obj.NCurves,1);
            obj.CorruptedCurves = false(obj.NCurves,1);
            
            try
                obj.read_jpk_images_from_files
            catch
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
            
            obj.construct_list_to_map_relations
            obj.create_pixel_difference_channel
            
            obj.initialize_flags();
        end
        
        function choose_curves(obj)
            % Interactive dialogue for curve selection. Returns a logic 0/1
            % vector, that is later used and can also be changed by the
            % other force map methods.
            f = figure('Name','Curve Selection','Position',[10000 10000 1200 900]);
            movegui(f);
            for i=1:obj.NCurves
                [App,HHApp] = obj.get_force_curve_data(i,0,0,0);
                [Ret,HHRet] = obj.get_force_curve_data(i,1,0,0);
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
                    [App,HHApp] = obj.get_force_curve_data(i,0,0,0);
                    if (length(App) < Thresh) || (length(HHApp) < Thresh)
                        obj.SelectedCurves(i) = 0;
                    end
                end
                if AppRetSwitch==1 || AppRetSwitch==2
                    [Ret,HHRet] = obj.get_force_curve_data(i,1,0,0);
                    if (length(Ret) < Thresh) || (length(HHRet) < Thresh)
                        obj.SelectedCurves(i) = 0;
                    end
                end
                
            end
        end
        
        function base_and_tilt(obj,RunMode,TiltCorrectionBool)
            % subtract baseline and tilt from the forcecurve by fitting a function to
            % the non contact domain the function tries to fit a non-affine-linear
            % function. If the linear fit is too bad, the function tries a 9th grade
            % polynomial fit instead
            
            if nargin < 3
                TiltCorrectionBool = true;
            end
            
            Range = find(obj.SelectedCurves);
            h = waitbar(0,'Setting up...');
            for i=Range'
                [AppForce,HHApp] = obj.get_force_curve_data(i,0,0,0);
                
                prog = i/obj.NCurves;
                waitbar(prog,h,'processing baseline fits...');
                try
                    [ncd , ncidx] = ForceMap.no_contact_domain(AppForce);
                    Params = polyfit(HHApp(1:length(ncd)),ncd,1);
                    if TiltCorrectionBool
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
                [App,HHApp] = obj.get_force_curve_data(i,0,1,0);
                [Ret,HHRet] = obj.get_force_curve_data(i,1,1,0);
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
                [App,HHApp] = obj.get_force_curve_data(i,0,1,0);
                [Ret,HHRet] = obj.get_force_curve_data(i,1,1,0);
                AboveZeroBool = zeros(length(App),1);
                AboveZeroBool(find(App>0)) = 1;
                k = 0;
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
                WindowSize = 20;
            end
            Range = find(obj.SelectedCurves);
            h = waitbar(0,'Setting up...','Name',obj.Name);
            CP_RoV = zeros(obj.NCurves,2);
            obj.CP_RoV = CP_RoV;
            for i=Range'
                prog = i/obj.NCurves;
                [App,HHApp] = obj.get_force_curve_data(i,0,1,0);
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
                [Based,THApp] = obj.get_force_curve_data(i,0,1,1);
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
                [App,HHApp] = obj.get_force_curve_data(i,0,1,1);
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
                    NumPasses = 20; % if not specified in arguments, NumPasses defaults to 20
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
                    [r,w] = unix('free | grep Mem');
                    stats = str2double(regexp(w, '[0-9]*', 'match'));
                    memsize = stats(1)/1e6;
                    freemem = (stats(3)+stats(end))/1e6;
                    MaxArraySize = freemem*10^9;
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
                                    Ypredicted = predict(NeuralNet,X,'MiniBatchSize',obj.MiniBatchSize,'Acceleration','auto');
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
                                [App,HHApp] = obj.get_force_curve_data(i,0,1,0);
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
                                        Temp = predict(NeuralNet,X,'MiniBatchSize',obj.MiniBatchSize,'Acceleration','auto');
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
                                [App,HHApp] = obj.get_force_curve_data(i,0,1,0);
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
                                    Ypredicted = predict(NeuralNet,X,'MiniBatchSize',obj.MiniBatchSize,'Acceleration','auto');
                                    CantHandle = false;
                                catch
                                    obj.CPFlag.CNNopt = 0;
                                    obj.cnn_runtime_optimization(NeuralNet,X);
                                end
                            end
                            iRange = find(obj.SelectedCurves);
                            k = 1;
                            for i=iRange'
                                [App,HHApp] = obj.get_force_curve_data(i,0,1,0);
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
                                    Ypredicted = predict(NeuralNet,X,'MiniBatchSize',obj.MiniBatchSize,'Acceleration','auto');
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
                                [App,HHApp] = ZoomObj.get_force_curve_data(i,0,1,0);
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
                                    Ypredicted = predict(NeuralNet,X,'MiniBatchSize',obj.MiniBatchSize,'Acceleration','auto');
                                    CantHandle = false;
                                catch
                                    obj.CPFlag.CNNopt = 0;
                                    obj.cnn_runtime_optimization(NeuralNet,X);
                                end
                            end
                            iRange = find(obj.SelectedCurves);
                            k = 1;
                            for i=iRange'
                                [App,HHApp] = obj.get_force_curve_data(i,0,1,0);
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
                                    Ypredicted = predict(NeuralNet,X,'MiniBatchSize',obj.MiniBatchSize,'Acceleration','auto');
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
                                [App,HHApp] = ZoomObj.get_force_curve_data(i,0,1,0);
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
                                    Ypredicted = predict(NeuralNet,X,'MiniBatchSize',obj.MiniBatchSize,'Acceleration','auto');
                                    CantHandle = false;
                                catch
                                    obj.CPFlag.CNNopt = 0;
                                    obj.cnn_runtime_optimization(NeuralNet,X);
                                end
                            end
                            iRange = find(obj.SelectedCurves);
                            k = 1;
                            for i=iRange'
                                [App,HHApp] = obj.get_force_curve_data(i,0,1,0);
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
                                    Ypredicted = predict(NeuralNet,X,'MiniBatchSize',obj.MiniBatchSize,'Acceleration','auto');
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
                                    [App,HHApp] = ZoomCell{j}.get_force_curve_data(i,0,1,0);
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
                [App,HHApp] = obj.get_force_curve_data(i,0,1,0);
                [Ret,HHRet] = obj.get_force_curve_data(i,1,1,0);
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
                [obj.LoadOld{i},obj.UnloadOld{i},Position,vDef] = ContactPoint_sort(load,unload);
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
                [App,HHApp] = obj.get_force_curve_data(j,0,1,0);
                [Ret,HHRet] = obj.get_force_curve_data(j,1,1,0);
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
                [App,HHApp] = obj.get_force_curve_data(i,0,1,0);
                obj.CP_HardSurface(i,1) = HHApp(end) - App(end)/obj.SpringConstant;
                obj.CP_HardSurface(i,2) = 0;
                %% Debugging
                % plot(HHApp,App);
                % drawpoint('Position',[obj.CP_HardSurface(i,1) obj.CP_HardSurface(i,2)]);
            end
            obj.CPFlag.HardSurface = 1;
        end
        
        function [E,HertzFit] = calculate_e_mod_hertz(obj,CPType,TipShape,curve_percent,AllowXShift,CorrectSensitivity,UseTipData,UseTopology,TipObject)
            % [E,HertzFit] = calculate_e_mod_hertz(obj,CPType,TipShape,curve_percent)
            %
            % calculate the E modulus of the chosen curves using the CP
            % type chosen in the arguments fitting the upper curve_percent
            % part of the curves
            if ~exist('curve_percent','var') && ischar('curve_percent')
                curve_percent = 0.75;
            end
            if ~exist('TipShape','var') && ~ischar(TipShape)
                TipShape = 'parabolic';
            end
            if nargin < 5
                AllowXShift = false;
                CorrectSensitivity = false;
                UseTipData = false;
            end
            iRange = find(obj.SelectedCurves);
            obj.EModHertz = zeros(obj.NCurves,1);
            obj.IndentationDepth = zeros(obj.NCurves,1);
            while ~isempty(iRange')
                NumWorkers = 8;
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
                    [App{i},HHApp{i}] = obj.get_force_curve_data(iRange(i),0,1,0);
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
                    obj.IndentationDepth(iRange(i)) = Max{i}(1);
                    % delete everything below curve_percent of the maximum
                    % force
                    force{i}(1:(length(force{i})-length(tip_h{i}))) = [];
                    force{i}(force{i}<(1-curve_percent)*max(force{i})) = [];
                    tip_h{i}(1:(length(tip_h{i})-length(force{i}))) = [];
                    RangeF{i} = range(force{i});
                    RangeTH{i} = range(tip_h{i});
                    force{i} = force{i}/RangeF{i};
                    tip_h{i} = tip_h{i}/RangeTH{i};
                end
                parfor i=1:BatchSize
                    if AllowXShift
                        s = fitoptions('Method','NonlinearLeastSquares',...
                            'Lower',[10^(-5) -min(tip_h{i})],...
                            'Upper',[inf min(tip_h{i})],...
                            'MaxIter',100,...
                            'Startpoint',[1 0]);
                        f = fittype('a*(x+b)^(3/2)','options',s);
                    else
                        s = fitoptions('Method','NonlinearLeastSquares',...
                            'Lower',10^(-5),...
                            'Upper',inf,...
                            'Startpoint',1);
                        f = fittype('a*(x)^(3/2)','options',s);
                    end
                    try
                        Hertzfit{i} = fit(tip_h{i},...
                            force{i},f);
                    catch
                        Hertzfit{i} = nan;
                    end
                end
                for i=1:BatchSize
                    try
                        % calculate E module based on the Hertz model. Be careful
                        % to convert to unnormalized data again
                        if isempty(obj.FibDiam) || UseTopology
                            if UseTipData
                                if UseTopology
                                    
                                else
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
                                    R_eff = TipRadius;
                                end
                            else
                                R_eff = obj.TipRadius*1e-9;
                            end
                        else
                            if UseTipData
                                DepthIndex = floor(obj.IndentationDepth(iRange(i))*1e9);
                                DepthRemainder = obj.IndentationDepth(iRange(i)) - DepthIndex;
                                if DepthIndex >= length(TipObject.DepthDependendTipRadius)
                                    DepthIndex = length(TipObject.DepthDependendTipRadius) - 1;
                                end
                                if DepthIndex == 0
                                    TipRadius = TipObject.DepthDependendTipRadius(DepthIndex+1)*DepthRemainder;
                                else
                                    TipRadius = TipObject.DepthDependendTipRadius(DepthIndex)*(1-DepthRemainder) + TipObject.DepthDependendTipRadius(DepthIndex+1)*DepthRemainder;
                                end
                                R_eff = 1/(1/TipRadius + 1/(obj.FibDiam/2));
                            else
                                R_eff = 1/(1/(obj.TipRadius*1e-9) + 1/(obj.FibDiam/2));
                            end
                        end
                        EMod{i} = 3*(Hertzfit{i}.a*RangeF{i}/RangeTH{i}^(3/2))/(4*sqrt(R_eff))*(1-obj.PoissonR^2);
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
                        if AllowXShift
                            Hertzfit{i}.b = Hertzfit{i}.b*RangeTH{i};
                            obj.CP_HertzFitted(iRange(i),1) = CP(i,1)-Hertzfit{i}.b;
                            obj.CP_HertzFitted(iRange(i),2) = CP(i,2);
                            % Not sure about this one
                            % obj.IndentationDepth(i) = obj.IndentationDepth(i) + Hertzfit.b;
                        end
                        warning('on','all');
                        obj.HertzFit{iRange(i)} = Hertzfit{i};
                    catch
                        obj.SelectedCurves(iRange(i)) = 0;
                        warning('on','all');
                    end
                end
                iRange(1:BatchSize) = [];
            end
            E = obj.EModHertz;
            HertzFit = obj.HertzFit;
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
                end
            end
            
            % Write to Channel
            EModMap = obj.create_standard_channel(obj.EModMapHertz(:,:,1),'E-Mod Hertz','Pa');
            [~,Index] = obj.get_channel('E-Mod Hertz');
            if isempty(Index)
                obj.Channel(end+1) = EModMap;
            else
                obj.Channel(Index) = EModMap;
            end
            
            EModLog = obj.create_standard_channel(log(obj.EModMapHertz(:,:,1)),'E-Mod Hertz (log)','Pa');
            [~,Index] = obj.get_channel('E-Mod Hertz (log)');
            if isempty(Index)
                obj.Channel(end+1) = EModLog;
            else
                obj.Channel(Index) = EModLog;
            end
            
            % Write to Channel
            Channel = obj.create_standard_channel(IndDepMap,'Indentation Depth','m');
            [~,Index] = obj.get_channel('Indentation Depth');
            if isempty(Index)
                obj.Channel(end+1) = Channel;
            else
                obj.Channel(Index) = Channel;
            end
            
            if AllowXShift
                obj.CPFlag.HertzFitted = 1;
            end
            
        end
        
        function EMod = calculate_e_mod_oliverpharr(obj,TipProjArea,CurvePercent)
            
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
                [Ret,HHRet] = obj.get_force_curve_data(i,1,1,0);
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
            EModMap = obj.create_standard_channel(obj.EModMapOliverPharr(:,:,1),'E-Mod OliverPharr','Pa');
            [~,Index] = obj.get_channel('E-Mod OliverPharr');
            if isempty(Index)
                obj.Channel(end+1) = EModMap;
            else
                obj.Channel(Index) = EModMap;
            end
            
            EModLog = obj.create_standard_channel(log(obj.EModMapOliverPharr(:,:,1)),'E-Mod OliverPharr (log)','Pa');
            [~,Index] = obj.get_channel('E-Mod OliverPharr (log)');
            if isempty(Index)
                obj.Channel(end+1) = EModLog;
            else
                obj.Channel(Index) = EModLog;
            end
            
            % Write to Channel
            Channel = obj.create_standard_channel(IndDepthMapOP,'Indentation Depth Oliver-Pharr','m');
            [~,Index] = obj.get_channel('Indentation Depth Oliver-Pharr');
            if isempty(Index)
                obj.Channel(end+1) = Channel;
            else
                obj.Channel(Index) = Channel;
            end
            
            % Write to Channel
            Channel = obj.create_standard_channel(IndDepthMap,'Indentation Depth','m');
            [~,Index] = obj.get_channel('Indentation Depth');
            if isempty(Index)
                obj.Channel(end+1) = Channel;
            else
                obj.Channel(Index) = Channel;
            end
            
            % Write to Channel
            Channel = obj.create_standard_channel(DZMap,'DZ-Slope','m/m');
            [~,Index] = obj.get_channel('DZ-Slope');
            if isempty(Index)
                obj.Channel(end+1) = Channel;
            else
                obj.Channel(Index) = Channel;
            end
        end
        
        function MaxAdhesionForce = calculate_adhesion_force(obj)
            
            Range = find(obj.SelectedCurves);
            MaxAdhesionForce = zeros(obj.NCurves,1);
            for i=Range'
                [Ret,HHRet] = obj.get_force_curve_data(i,1,1,0);
                Force = Ret - obj.CP(i,2);
                TipHeight = (HHRet - obj.CP(i,1)) - Force/obj.SpringConstant;
                MaxAdhesionForce(i) = -min(Force);
            end
            
            obj.MaxAdhesionForce = MaxAdhesionForce;
            
            for i=1:obj.NumPixelsX
                for j=1:obj.NumPixelsY
                    MaxAdhesionForceMap(i,j) = obj.MaxAdhesionForce(obj.Map2List(i,j));
                end
            end
            
            % Write to Channel
            Channel = obj.create_standard_channel(MaxAdhesionForceMap,'Maximum Adhesion Force','N');
            [~,Index] = obj.get_channel('Maximum Adhesion Force');
            if isempty(Index)
                obj.Channel(end+1) = Channel;
            else
                obj.Channel(Index) = Channel;
            end
        end
        
        function [AdhesionEnergy,AdhesionLength] = calculate_adhesion_energy_and_length(obj,ThresholdMult)
            
            if nargin < 2
                ThresholdMult = 2;
            end
            
            Range = find(obj.SelectedCurves);
            AdhesionEnergy = zeros(obj.NCurves,1);
            AdhesionLength = zeros(obj.NCurves,1);
            for i=Range'
                [App,HHApp] = obj.get_force_curve_data(i,0,1,0);
                [Ret,HHRet] = obj.get_force_curve_data(i,1,1,0);
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
            
            for i=1:obj.NumPixelsX
                for j=1:obj.NumPixelsY
                    AEMap(i,j) = obj.AdhesionEnergy(obj.Map2List(i,j));
                    ALMap(i,j) = obj.AdhesionLength(obj.Map2List(i,j));
                end
            end
            
            % Write to Channel
            Channel1 = obj.create_standard_channel(AEMap,'Adhesion Energy','J');
            [~,Index] = obj.get_channel('Adhesion Energy');
            if isempty(Index)
                obj.Channel(end+1) = Channel1;
            else
                obj.Channel(Index) = Channel1;
            end
            
            
            CoulombConstant = 1.602176634e-19;
            % Write to Channel
            Channel1 = obj.create_standard_channel(AEMap./CoulombConstant,'Adhesion eV-Energy','eV');
            [~,Index] = obj.get_channel('Adhesion eV-Energy');
            if isempty(Index)
                obj.Channel(end+1) = Channel1;
            else
                obj.Channel(Index) = Channel1;
            end
            
            Channel2 = obj.create_standard_channel(ALMap,'Adhesion Length','m');
            [~,Index] = obj.get_channel('Adhesion Length');
            if isempty(Index)
                obj.Channel(end+1) = Channel2;
            else
                obj.Channel(Index) = Channel2;
            end
        end
        
        function PeakIndentationAngle = calculate_peak_indentation_angle(obj,FitPortion)
            
            if nargin < 2
                FitPortion = .5;
            end
            
            Range = find(obj.SelectedCurves);
            PeakIndentationAngle = zeros(obj.NCurves,1);
            for i=Range'
                [App,HHApp] = obj.get_force_curve_data(i,0,1,0);
                [Ret,HHRet] = obj.get_force_curve_data(i,1,1,0);
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
            
            for i=1:obj.NumPixelsX
                for j=1:obj.NumPixelsY
                    PIAMap(i,j) = obj.PeakIndentationAngle(obj.Map2List(i,j));
                end
            end
            
            % Write to Channel
            Channel = obj.create_standard_channel(PIAMap,'Peak Indentation Angle','Degrees');
            [~,Index] = obj.get_channel('Peak Indentation Angle');
            if isempty(Index)
                obj.Channel(end+1) = Channel;
            else
                obj.Channel(Index) = Channel;
            end
        end
        
        function [DissipatedEnergy,ElasticEnergy] = calculate_dissipated_and_elastic_energy(obj)
            
            Range = find(obj.SelectedCurves);
            DissipatedEnergy = zeros(obj.NCurves,1);
            ElasticEnergy = zeros(obj.NCurves,1);
            for i=Range'
                [App,HHApp] = obj.get_force_curve_data(i,0,1,0);
                [Ret,HHRet] = obj.get_force_curve_data(i,1,1,0);
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
            
            for i=1:obj.NumPixelsX
                for j=1:obj.NumPixelsY
                    DEMap(i,j) = obj.DissipatedEnergy(obj.Map2List(i,j));
                end
            end
            
            % Write to Channel
            Channel = obj.create_standard_channel(DEMap,'Dissipated Energy','J');
            [~,Index] = obj.get_channel('Dissipated Energy');
            if isempty(Index)
                obj.Channel(end+1) = Channel;
            else
                obj.Channel(Index) = Channel;
            end
            
            DEMap = DEMap./CoulombConstant;
            
            % Write to Channel
            Channel = obj.create_standard_channel(DEMap,'Dissipated eV-Energy','eV');
            [~,Index] = obj.get_channel('Dissipated eV-Energy');
            if isempty(Index)
                obj.Channel(end+1) = Channel;
            else
                obj.Channel(Index) = Channel;
            end
            
            obj.ElasticEnergy = ElasticEnergy;
            
            for i=1:obj.NumPixelsX
                for j=1:obj.NumPixelsY
                    EEMap(i,j) = obj.ElasticEnergy(obj.Map2List(i,j));
                end
            end
            
            % Write to Channel
            Channel = obj.create_standard_channel(EEMap,'Elastic Energy','J');
            [~,Index] = obj.get_channel('Elastic Energy');
            if isempty(Index)
                obj.Channel(end+1) = Channel;
            else
                obj.Channel(Index) = Channel;
            end
            
            
            EEMap = EEMap./CoulombConstant;
            
            % Write to Channel
            Channel = obj.create_standard_channel(EEMap,'Elastic eV-Energy','eV');
            [~,Index] = obj.get_channel('Elastic eV-Energy');
            if isempty(Index)
                obj.Channel(end+1) = Channel;
            else
                obj.Channel(Index) = Channel;
            end
            
            % Write to Channel
            Channel = obj.create_standard_channel(100.*EEMap./(DEMap+EEMap),'Elastic Fraction','%');
            [~,Index] = obj.get_channel('Elastic Fraction');
            if isempty(Index)
                obj.Channel(end+1) = Channel;
            else
                obj.Channel(Index) = Channel;
            end
            
            
            % Write to Channel
            Channel = obj.create_standard_channel(100.*DEMap./(DEMap+EEMap),'Inelastic Fraction','%');
            [~,Index] = obj.get_channel('Inelastic Fraction');
            if isempty(Index)
                obj.Channel(end+1) = Channel;
            else
                obj.Channel(Index) = Channel;
            end
        end
        
        function manual_exclusion(obj)
            
            obj.ExclMask = logical(ones(obj.NumPixelsX,obj.NumPixelsY));
            CheckSum = 100;
            while CheckSum > 1
                f = figure('Name','Choose areas to be excluded');
                f.WindowState = 'maximized';
                subplot(2,1,2)
                surf(imresize(imrotate(obj.HeightMap(:,:,1)',90),[1024 1024]),'LineStyle','none','FaceLighting','gouraud','FaceColor','interp')
                light('Style','local')
                subplot(2,1,1)
                imshow(obj.HeightMap(:,:,1).*obj.ExclMask,[min(obj.HeightMap(:,:,1),[],'all') max(obj.HeightMap(:,:,1),[],'all')])
                title(sprintf('%s: Draw Freehand ROI around areas, that are to be excluded\nThe area will be taken out and the same map redrawn \n If there is nothing to do just click on the image once without dragging the cursor',obj.Name))
                ROI = drawfreehand;
                CheckSum = length(ROI.Waypoints);
                Mask = ~createMask(ROI);
                obj.ExclMask = obj.ExclMask.*Mask;
                close(f)
            end
            
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
                obj.RectApexIndex = round(predictGP_mean([1:obj.NumPixelsX],[1:obj.NumPixelsX],1,5*obj.NumPixelsX,obj.ApexIndex,1));
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
            obj.FibDiamSTD = nanstd(FibHeight);
            
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
                [App,~] = obj.get_force_curve_data(i,0,1,0);
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
    end    
    %% Fc related
    methods
       
        function fc_measurement_prop(obj)
            % A fct to read out measurement properties written in the name of the jpk data file (i.e. "JPK-FORCE-MAP").
            
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
                obj.Substrate='';
            end
            % Milli-Q water
            exp31a='mil';
            pata=regexpPattern(exp31a,"IgnoreCase",true);
            ext31a=extract(obj.Name,pata);
            exp31b='HOH';
            patb=regexpPattern(exp31b,"IgnoreCase",true);
            ext31b=extract(obj.Name,patb);
            % Acetic acid (HAc)
            exp32='HAc';
            pat=regexpPattern(exp32,"IgnoreCase",true);
            ext32=extract(obj.Name,pat);
            % CAPS buffer
            exp33='CAPS';
            pat=regexpPattern(exp33,"IgnoreCase",true);
            ext33=extract(obj.Name,pat);
            % Dulbecco's phospate buffered saline (DPBS) buffer
            exp34='PBS';
            pat=regexpPattern(exp34,"IgnoreCase",true);
            ext34=extract(obj.Name,pat);
            % Dulbecco's phospate buffered saline (DPBS) buffer
            exp35='DPBS';
            pat=regexpPattern(exp35,"IgnoreCase",true);
            ext35=extract(obj.Name,pat);
            % CaCl2
            exp36='CaCl2';
            pat=regexpPattern(exp36,"IgnoreCase",true);
            ext36=extract(obj.Name,pat);            
            % Environmental conditions
            if ~isempty(ext31a) || ~isempty(ext31b)
                obj.EnvCond='Water'; % Milli-Q water
            elseif ~isempty(ext32) && ~isempty(ext36)
                obj.EnvCond='HAcCaCl2'; % Acetic acid
            elseif ~isempty(ext32)
                obj.EnvCond='HAc'; % Acetic acid
            elseif ~isempty(ext33)
                obj.EnvCond='CAPS'; % CAPS
            elseif ~isempty(ext34)
                obj.EnvCond='PBS'; % Phospate buffered saline
            elseif ~isempty(ext35)
                obj.EnvCond='DPBS'; % Dulbecco's phospate buffered saline
            elseif ~isempty(ext36)
                obj.EnvCond='CaCl2'; % Calcium chloride solution
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
                obj.Linker='';

                %obj.Linker='amino'; 
                %obj.Linker='naked';                 
                %obj.Linker='long';
                %obj.Linker='short';
            end

        end
                        
        function fc_sensitivity_correction(obj)
            % Function to correct the measured sensitivity before the
            % AFM SMFS measurment with the sensitivity data from each
            % individual force curve

        for Fc=1:obj.NCurves
            %for Fc=1
            % Allocate data
            xApp=obj.HHApp{Fc};
            % Convert y data from Force to distance:
            % Force(N)/SpringConstant(N/m)
            yApp=obj.App{Fc}/obj.SpringConstant; % Convert y data into a distance
            yRet=obj.Ret{Fc}/obj.SpringConstant; % Convert y data into a distance
            % 1. Determine maximum value of y-data and select data
            [yAppMax,yAppMaxIdx]=max(yApp); % Read out the maximum value and corresponding index
            [yAppMin,yAppMinIdx]=min(yApp); % Read out the maximum value and corresponding index
            AppSelIdx1=find(yApp>yAppMin+(yAppMax-yAppMin)/3,1,'first'); % Read out the index of the first cell that fulfils the argument
            yAppSel=yApp(AppSelIdx1:yAppMaxIdx); % Select the data used for the fit
            xAppSel=xApp(AppSelIdx1:yAppMaxIdx); % Select corresponding x data
            % 2. Fit the selected data
            %  Create a fit.
            %
            %  Data for 'untitled fit 1' fit:
            %      X Input : xAppSel
            %      Y Output: yAppSel
            %  Output:
            %      fitresult : a fit object representing the fit.
            %      gof : structure with goodness-of fit info.
            %
            %  See also FIT, CFIT, SFIT.
            %  Auto-generated by MATLAB on 25-May-2022 10:17:39
            % Fittype: linear
            % y=a*x+b
            % Fit
            [xData, yData] = prepareCurveData( xAppSel, yAppSel );
            % Set up fittype and options.
            ft = fittype( {'x', '1'}, 'independent', 'x', 'dependent', 'y', 'coefficients', {'a', 'b'} );
            % Fit model to data.
            [fitresult, gof] = fit( xData, yData, ft );
            % 3. Determine sensitivity and allocate data
            obj.SensitivityFcSlope(Fc)=fitresult.a; % Reference slope (1)
            obj.SensitivityFc(Fc)=obj.Sensitivity/fitresult.a; % Old Sensitivity (m/V) divided by the reference slope (1) = New Sensitivity (m/V)
            obj.SensitivityFcLinFitRSquare(Fc)=gof.rsquare;
            obj.SensitivityFcLinFitSSE(Fc)=gof.sse;
            % 4. Recalculate corrected y data 
            % Force (N) = Force (N) / old Sensitivity (m/V) * new Sensitivity (m/V)
            obj.AppCorr{Fc}=obj.App{Fc}/obj.Sensitivity*obj.SensitivityFc(Fc); % Convert y data into Force
            obj.RetCorr{Fc}=obj.Ret{Fc}/obj.Sensitivity*obj.SensitivityFc(Fc); % Convert y data into Force
            % 5. Set flag
            obj.SMFSFlag.FitSensitivity(Fc)=1;
         end

    end
    
        function fc_sinoidal_fit(obj)
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
            RGB7=[255 230 0]./255; % Yellow
            RGB10=[200 0 255]./255; % Violet
            mega=10e6;
            giga=10e9;
            DataShareStartApp=0.15; % 15% 
            DataShareEndApp=0.65; % 65%
            NFigures=4;
            NLoop=25;
            ExtendVelocityConvert=num2str(round(obj.ExtendVelocity*1e9));
            RetractVelocityConvert=num2str(round(obj.RetractVelocity*1e9));
            HoldingTimeConvert=num2str(round(obj.HoldingTime));        
            Res=[1 1 2560 1440]; % Define the figure resolution
            % Classification criteria
            figname=strcat(obj.Date,{'_'},obj.Time,{'_'},obj.ID,{'_'},obj.Substrate,{'_'},obj.EnvCond,{'_'},obj.Linker,{'_'},obj.Chipbox,{'_'},obj.ChipCant,{'_'},ExtendVelocityConvert,{'_'},RetractVelocityConvert,{'_'},HoldingTimeConvert);
            figname=char(figname);
            % Fit loop
             for Fc=1:obj.NCurves
                %% Debugging
             %for kk=88 % for debugging
                % sprintf('Force curve No. %d',kk) % Gives current Force curve
                % for debugging
                % Allocate data
                xApp=obj.HHApp{Fc}*mega;
                xRet=obj.HHRet{Fc}*mega; % Transform head height to nm
                yApp=obj.AppCorr{Fc}*giga; % Transform approach force data to pN
                yRet=obj.RetCorr{Fc}*giga; % Transform approach force data to pN
                xAppSel=obj.HHApp{Fc}*mega;
                yAppSel=obj.AppCorr{Fc}*giga;
                % Select only a percentage of the data so that a potential
                % snap-in is excluded from the fit
                % Define limits
                DataPtsApp=size(yAppSel); % Determine the amount of data points in the force curve
                LimitIdxApp1=round(DataPtsApp(1)*DataShareStartApp); % Determine the corresponding index
                LimitIdxApp2=round(DataPtsApp(1)*DataShareEndApp); % Determine the corresponding index
                % Select data points for the fit
                xAppSel([1:LimitIdxApp1, LimitIdxApp2:end])=[]; % Remove all entries except of the selections made using the limits
                yAppSel([1:LimitIdxApp1, LimitIdxApp2:end])=[]; % Remove all entries except of the selections made using the limits 
                % Prepare data
                [xData, yData] = prepareCurveData( xAppSel, yAppSel );
                % Set up fittype and options.
                % Fit 1
                ft1 = fittype( 'a1*sin(b1*x+c1)+d1', 'independent', 'x', 'dependent', 'y' );
                opts1= fitoptions( 'Method', 'NonlinearLeastSquares' );
                opts1.Display = 'Off';
                opts1.Robust = 'Bisquare';
                a1Range=range(yData)/2;
                %opts1.StartPoint = [a1Range -1.8 x(1,1) y(1,1)];
                opts1.StartPoint = [-0.2 -1.8 xAppSel(1,1) yAppSel(1,1)];
                [fitresult1, gof1] = fit( xData, yData, ft1, opts1 ); % Fit
                % Allocate variables
                obj.FitxData{Fc}=xData;
                obj.FityData{Fc}=yData;
                obj.SinFit1Coeffa(Fc)=fitresult1.a1;
                obj.SinFit1Coeffb(Fc)=fitresult1.b1;
                obj.SinFit1Coeffc(Fc)=fitresult1.c1;
                obj.SinFit1Coeffd(Fc)=fitresult1.d1;
                obj.SinFit1RSquare(Fc)=gof1.rsquare;
                obj.SinFit1SSE(Fc)=gof1.sse;
                % Fit 2 - Sinoidal fit
                % Set up fittype and options.
                ft2 = fittype( 'a1*sin(b1*x+c1)+d1+e1*x', 'independent', 'x', 'dependent', 'y' );
                opts2= fitoptions( 'Method', 'NonlinearLeastSquares' );
                opts2.Display = 'Off';
                opts2.StartPoint = [fitresult1.a1 fitresult1.b1 fitresult1.c1 fitresult1.d1 0];
                [fitresult2, gof2] = fit( xData, yData, ft2, opts2 ); % Fit                
                % Allocate variables                
                obj.SinFit2Coeffa(Fc)=fitresult2.a1;
                obj.SinFit2Coeffb(Fc)=fitresult2.b1;
                obj.SinFit2Coeffc(Fc)=fitresult2.c1;
                obj.SinFit2Coeffd(Fc)=fitresult2.d1;
                obj.SinFit2Coeffe(Fc)=fitresult2.e1;
                obj.SinFit2RSquare(Fc)=gof2.rsquare;
                obj.SinFit2SSE(Fc)=gof2.sse;
                % Fit 3 - Linear fit
                % Allocate data
                % Set up fittype and options.
                ft3 = fittype( {'x', '1'}, 'independent', 'x', 'dependent', 'y', 'coefficients', {'a', 'b'} );
                [fitresult3, gof3] = fit( xData, yData, ft3 );          
                % Allocate variables
                obj.LinFitCoeffa(Fc)=fitresult3.a;
                obj.LinFitCoeffb(Fc)=fitresult3.b;
                obj.LinFitRSquare(Fc)=gof3.rsquare;
                obj.LinFitSSE(Fc)=gof3.sse;
                %% Compare fits                
                % Choose better fit
                if obj.SinFit2RSquare(Fc)>obj.LinFitRSquare(Fc) % Apply sinoidal fit
                obj.BasedApp{Fc}=(yApp-feval(fitresult2,xApp))/giga; % Determine the corresponding y-data of the x retraction data by using the determined fitting coefficients
                obj.BasedRet{Fc}=(yRet-feval(fitresult2,xRet))/giga; % Determine the corresponding y-data of the x retraction data by using the determined fitting coefficients    
                obj.BasedRetFit{Fc}=(yRet-feval(fitresult2,xRet))/giga; % Determine the corresponding y-data of the x retraction data by using the determined fitting coefficients      
                obj.SMFSFlag.FitSinoidal(Fc)=1;
                else % Apply linear fit
                obj.BasedApp{Fc}=(yApp-feval(fitresult3,xApp))/giga; % Determine the corresponding y-data of the x retraction data by using the determined fitting coefficients
                obj.BasedRet{Fc}=(yRet-feval(fitresult3,xRet))/giga; % Determine the corresponding y-data of the x retraction data by using the determined fitting coefficients    
                obj.BasedRetFit{Fc}=(yRet-feval(fitresult3,xRet))/giga; % Determine the corresponding y-data of the x retraction data by using the determined fitting coefficients            
                obj.SMFSFlag.FitLinear(Fc)=1;
                end                
            end
            % Figure loop
            for jj=1:NFigures
                %% Plot condition                
                %if  ~obj.DebugFlag.Plot % Suppress plotting
                if  obj.DebugFlag.Plot % Allow plotting
                    continue
                end
                % Figure properties
                h_fig=figure(jj);
                h_fig.Color='white'; % changes the background color of the figure
                h_fig.Units='pixel'; % Defines the units
                h_fig.OuterPosition=Res;
                h_fig.PaperOrientation='landscape';
                h_fig.Name=figname;
                % Tile properties
                t = tiledlayout(5,5);
                t.TileSpacing = 'none'; % To reduce the spacing between the tiles
                t.Padding = 'none'; % To reduce the padding of perimeter of a tile
                %% Tile loop
                for Fc=1:NLoop
                    oo=Fc+25*(jj-1); % Tile variable
                    %% Debugging
                    % oo=45 % for debugging
                    % sprintf('Force curve No. %d',oo) % Gives current Force curve
                    % for debugging
                    % Allocate data
                    xData=obj.FitxData{oo};
                    yData=obj.FityData{oo};
                    %% Plotting the tiles
                    nexttile;
                    hold on
                    grid on
                    plot(xData,yData,'.','Color',RGB4)
                    % plot(xData,obj.SinFit1Coeffa(oo)*sin(obj.SinFit1Coeffb(oo)*xData+obj.SinFit1Coeffc(oo))+obj.SinFit1Coeffd(oo),'LineWidth',2,'Color',RGB10)
                    if obj.SMFSFlag.FitSinoidal(oo)
                    plot(xData,obj.SinFit2Coeffa(oo)*sin(obj.SinFit2Coeffb(oo)*xData+obj.SinFit2Coeffc(oo))+obj.SinFit2Coeffd(oo)+obj.SinFit2Coeffe(oo)*xData,'LineWidth',2,'Color',RGB3) % The plot shows the second sinoidal fit including a linear fit element        
                    elseif  obj.SMFSFlag.FitLinear(oo)
                    plot(xData,obj.LinFitCoeffa(oo)*xData+obj.LinFitCoeffb(oo),'LineWidth',2,'Color',RGB3) % The plot shows the linear fit
                    end
                    % Plotted Text
                    NE = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.01; % Define the position in the plot
                    partstrA1='R^2=';
                    if obj.SMFSFlag.FitSinoidal(oo)
                    partstrA2=num2str(round(gof2.rsquare,2));                    
                    elseif obj.SMFSFlag.FitLinear(oo)
                    partstrA2=num2str(round(gof3.rsquare,2));     
                    end
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
                partname=sprintf('-fitfct1-part%d',jj);
                fullname=sprintf('%s%s',figname,partname);
                %%% Save the current figure in the current folder
                print(gcf,fullname,'-dpng');
            end           
            close all
        end
        
        function fc_linear_fit(obj)
            %CREATEFIT1(XCAT,YCAT)
            %  Create a fit.
            %
            %  Data for 'Linear fit' fit:
            %      X Input : xCat
            %      Y Output: yCat
            %  Output:
            %      fitresult : a fit object representing the fit.
            %      gof : structure with goodness-of fit info.
            %
            %  See also FIT, CFIT, SFIT.

            %  Auto-generated by MATLAB on 01-Oct-2021 15:32:09

            % Fittype: linear
            % y=a*x+b

            % Define variables
            mega=10e6;
            giga=10e9;
            RGB1=[0 26 255]./255;  % Blue
            RGB2=[255 119 0]./255; % Orange
            RGB3=[80 220 100]./255; % Emerald
            RGB4=[200 81 160]./255; % Compl to RGB3
            RGB7=[255 230 0]./255; % Yellow
            DataShareStartApp=0.15; % 15% 
            DataShareEndApp=0.65; % 65%
            DataShareStartRet=0.07; % 7 % 
            DataShareEndRet=0.02; % 2 %
            NFigures=4;
            NLoop=25;
            Res=[1 1 2560 1440]; % Define the figure resolution
            ExtendVelocityConvert=num2str(round(obj.ExtendVelocity*1e9));
            RetractVelocityConvert=num2str(round(obj.RetractVelocity*1e9));
            HoldingTimeConvert=num2str(round(obj.HoldingTime));
            % Classification criteria
            figname=strcat(obj.Date,{'_'},obj.Time,{'_'},obj.ID,{'_'},obj.Substrate,{'_'},obj.EnvCond,{'_'},obj.Linker,{'_'},obj.Chipbox,{'_'},obj.ChipCant,{'_'},ExtendVelocityConvert,{'_'},RetractVelocityConvert,{'_'},HoldingTimeConvert);
            figname=char(figname);
            % Fit loop
            for Fc=1:obj.NCurves
                %% Debugging
                % for kk=23 % for debugging
                % sprintf('Force curve No. %d',kk) % Gives current Force curve
                % for debugging
                xApp=obj.HHApp{Fc}*mega;
                yApp=obj.BasedApp{Fc}*giga;
                xRet=obj.HHRet{Fc}*mega;
                yRet=obj.BasedRet{Fc}*giga;
                xRetSel=obj.HHRet{Fc}*mega;
                yRetSel=obj.BasedRet{Fc}*giga;
                % Define limits
                DataPtsRet=size(yRetSel); % Determine the amount of data points in the force curve
                LimitIdxRet1=round(DataPtsRet(1)-DataPtsRet(1)*DataShareStartRet); % Determine the corresponding index
                LimitIdxRet2=round(DataPtsRet(1)-DataPtsRet(1)*DataShareEndRet); % Determine the corresponding index
                % Select data points for the fit
                xRetSel([1:LimitIdxRet1, LimitIdxRet2:end])=[]; % Remove all entries except of the selections made using the limits
                xRetSelMean=mean(xRetSel);
                yRetSel([1:LimitIdxRet1, LimitIdxRet2:end])=[]; % Remove all entries except of the selections made using the limits 
                yRetSelMean=mean(yRetSel);
                % Define limits
                DataPtsApp=size(yApp); % Determine the amount of data points in the force curve
                LimitIdxApp1=round(DataPtsApp(1)*DataShareStartApp); % Determine the corresponidng index
                LimitIdxApp2=round(DataPtsApp(1)*DataShareEndApp); % Determine the corresponidng index
                % Find retention index corresponding to the approach baseline
                yAppMean=mean(yApp(LimitIdxApp1:LimitIdxApp2,1)); % Calculate the mean of the difference data
                RetIdx=find(yRet<yAppMean,1,'first'); % Find the closest data point smaller than the computed mean
                obj.xOriRet(Fc)=mean(xRet(RetIdx-1:RetIdx)); % Approximate corresponding x-value to the computed mean
                obj.yOriRet(Fc)=yAppMean;
                % Concanate generated data with zero as first data value
                xCat=cat(1,obj.xOriRet(Fc),xRetSelMean);
                yCat=cat(1,obj.yOriRet(Fc),yRetSelMean);
                % Prepare data
                [xData, yData] = prepareCurveData(xCat,yCat);
                % Set up fittype and options.
                ft = fittype( {'x', '1'}, 'independent', 'x', 'dependent', 'y', 'coefficients', {'a', 'b'} );
                [fitresult, gof] = fit( xData, yData, ft );
                % Allocate variables
                obj.LinFitRetxData{Fc}=xData;
                obj.LinFitRetyData{Fc}=yData;
                obj.LinFitRetCoeffa(Fc)=fitresult.a;
                obj.LinFitRetCoeffb(Fc)=fitresult.b;
                obj.LinFitRetRSquare(Fc)=gof.rsquare;
                obj.LinFitRetSSE(Fc)=gof.sse;
                obj.BasedRet{Fc}=(yRet-feval(fitresult,xRet))/giga;  % Determine the corresponding y-data of the x retraction data by using the determined fitting coefficients               
                % Flag
                obj.SMFSFlag.FitLinearRet(Fc)=1;
                if obj.SMFSFlag.FitSinoidal(Fc)
                obj.SMFSFlag.Fit(Fc)=obj.SMFSFlag.FitSinoidal(Fc)*obj.SMFSFlag.FitLinearRet(Fc)*obj.SMFSFlag.FitSensitivity(Fc);
                elseif obj.SMFSFlag.FitLinear(Fc)
                obj.SMFSFlag.Fit(Fc)=obj.SMFSFlag.FitLinear(Fc)*obj.SMFSFlag.FitLinearRet(Fc)*obj.SMFSFlag.FitSensitivity(Fc);
                end
            end
            % Figure loop
            for jj=1:NFigures
                %if  ~obj.DebugFlag.Plot % Suppress plotting
                if  obj.DebugFlag.Plot % Allow plotting
                    continue
                end
                % Figure properties
                h_fig=figure(jj);
                h_fig.Color='white'; % changes the background color of the figure h_fig.Units='pixel'; % Defines the units
                h_fig.Units='pixel'; % Defines the units
                h_fig.OuterPosition=Res;
                h_fig.PaperOrientation='landscape';
                h_fig.Name=figname;
                % Tile properties
                t = tiledlayout(5,5);
                t.TileSpacing = 'none'; % To reduce the spacing between the tiles
                t.Padding = 'none'; % To reduce the padding of perimeter of a tile
                % Tile loop
                for Fc=1:NLoop
                    % Tile
                    oo=Fc+25*(jj-1);
                    % Allocate data
                    xApp=obj.HHApp{oo}*mega;
                    yApp=obj.BasedApp{oo}*giga;
                    xRet=obj.HHRet{oo}*mega;
                    yRet=obj.BasedRetFit{oo}*giga; % Allocate the fit function 1 fitted y data
                    xData=obj.LinFitRetxData{oo};
                    yData=obj.LinFitRetyData{oo};
                    %% Plot condition                   
                    %if  ~obj.DebugFlag.Plot % Suppress plotting
                    if  obj.DebugFlag.Plot % Allow plotting
                        continue
                    end
                    % Preparation for the plot
                    AppIdx=find(yApp<obj.yOriRet(oo),1,'last');
                    yApp(AppIdx:end)=[]; % Remove "high" values
                    xApp(AppIdx:end)=[]; % Remove "high" values
                    xRet(1:RetIdx)=[]; % Remove "high" values
                    yRet(1:RetIdx)=[]; % Remove "high" values
                    % Tiles
                    nexttile;
                    hold on
                    grid on
                    plot(xApp,yApp,'Color',RGB1)
                    plot(xRet,yRet,'Color',RGB2)
                    plot(xData,yData,'*','MarkerSize',4,'Color',RGB4)
                    plot(xData,obj.LinFitRetCoeffa(oo)*xData+obj.LinFitRetCoeffb(oo),'LineWidth',2,'Color',RGB3)
                    % Plotted Text
                    SW = [min(xlim) min(ylim)]+[diff(xlim)*2 diff(ylim)*2]*0.1;
                    partstr1='slope=';
                    partstr2=num2str(round(obj.LinFitRetCoeffa(Fc),2));                 
                    fullstr=strcat(partstr1,partstr2); % Define the string that shall be shown in the plot     
                    te1=text(SW(1),SW(2),fullstr, 'VerticalAlignment','top', 'HorizontalAlignment','right');
                    te1.FontSize = 14;
                    % Title for each Subplot
                    ti=title(sprintf('%i',oo),'Color','k');
                    ti.Units='normalized'; % Set units to 'normalized'
                    %% Save figures
                    %%% Define the name for the figure title
                    partname=sprintf('-fitfct2-part%d',jj);
                    fullname=sprintf('%s%s',figname,partname);
                    %%% Save the current figure in the current folder
                    print(gcf,fullname,'-dpng');
                end
            end
           close all
        end

        function fc_TipHeight_calculation(obj)
            % Function to determine the TipHeight based on the Head Height
            % and the fit based y data
            
            for Fc=1:obj.NCurves
            % for Fc=1
                obj.THApp{Fc} = obj.HHApp{Fc} - obj.BasedApp{Fc}/obj.SpringConstant;
                obj.THRet{Fc} = obj.HHRet{Fc} - obj.BasedRet{Fc}/obj.SpringConstant;
                %% Appendix
                % Define variables
%                 RGB3=[80 220 100]./255; % Emerald
%                 RGB4=[200 81 160]./255; % Compl to RGB3
%                 ExtendVelocityConvert=num2str(round(obj.ExtendVelocity*1e9));
%                 RetractVelocityConvert=num2str(round(obj.RetractVelocity*1e9));
%                 HoldingTimeConvert=num2str(round(obj.HoldingTime));
%                 Res=[1 1 2560 1440]; % Define the figure resolution
%                 % Classification criteria
%                 figname=strcat(obj.Date,{'_'},obj.Time,{'_'},obj.ID,{'_'},obj.Substrate,{'_'},obj.EnvCond,{'_'},obj.Linker,{'_'},obj.Chipbox,{'_'},obj.ChipCant,{'_'},ExtendVelocityConvert,{'_'},RetractVelocityConvert,{'_'},HoldingTimeConvert);
%                 figname=char(figname);
%                 % Figure properties               
%                 h_fig=figure(Fc);
%                 h_fig.Color='white'; % changes the background color of the figure h_fig.Units='pixel'; % Defines the units
%                 h_fig.Units='pixel'; % Defines the units
%                 h_fig.OuterPosition=Res;
%                 h_fig.PaperOrientation='landscape';
%                 h_fig.Name=figname;
%                 hold on
%                 grid on
%                 plot(obj.HHApp{Fc},obj.BasedApp{Fc},'Color',RGB3,'LineWidth',3);
%                 plot(obj.THApp{Fc},obj.BasedApp{Fc},'Color',RGB4,'LineWidth',3);
%                 % Axes
%                 ax=gca;
%                 ax.LineWidth=2;
%                 ax.FontSize = 16;
%                 ax.XLabel.String = 'Height (m)';
%                 ax.XLabel.FontSize = 20;
%                 ax.YLabel.String = 'Vertical deflection (N)';
%                 ax.YLabel.FontSize = 20;
%                 %% Save figures
%                 %%% Define the name for the figure title
%                 partname=sprintf('-THvsHH-Fc%d',Fc);
%                 fullname=sprintf('%s%s',figname,partname);
%                 %%% Save the current figure in the current folder
%                 print(gcf,fullname,'-dpng');
%                 close all
            end
        end
                                        
        function fc_estimate_cp_hardsurface(obj)
            % contact point estimation for force curves detected on hard
            % surfaces
            
            for Fc=1:obj.NCurves 
           %  for Fc=1
           %    obj.CP_HardSurface(Fc,1) = obj.HHApp{Fc}(end) - obj.BasedApp{Fc}(end)/obj.SpringConstant; % = THApp of the last entry point corresponding to the applied set point
              obj.CP_HardSurface(Fc,1) = obj.THApp{Fc}(end); % = THApp of the last entry point corresponding to the applied set point
              obj.CP_HardSurface(Fc,2) = 0;
                %% Appendix
                % Define variables
%                 RGB3=[80 220 100]./255; % Emerald
%                 RGB4=[200 81 160]./255; % Compl to RGB3
%                 ExtendVelocityConvert=num2str(round(obj.ExtendVelocity*1e9));
%                 RetractVelocityConvert=num2str(round(obj.RetractVelocity*1e9));
%                 HoldingTimeConvert=num2str(round(obj.HoldingTime));
%                 Res=[1 1 2560 1440]; % Define the figure resolution
%                 % Classification criteria
%                 figname=strcat(obj.Date,{'_'},obj.Time,{'_'},obj.ID,{'_'},obj.Substrate,{'_'},obj.EnvCond,{'_'},obj.Linker,{'_'},obj.Chipbox,{'_'},obj.ChipCant,{'_'},ExtendVelocityConvert,{'_'},RetractVelocityConvert,{'_'},HoldingTimeConvert);
%                 figname=char(figname);
%                 % Figure properties               
%                 h_fig=figure(Fc);
%                 h_fig.Color='white'; % changes the background color of the figure h_fig.Units='pixel'; % Defines the units
%                 h_fig.Units='pixel'; % Defines the units
%                 h_fig.OuterPosition=Res;
%                 h_fig.PaperOrientation='landscape';
%                 h_fig.Name=figname;
%                 hold on
%                 grid on
%                 plot(obj.THApp{Fc},obj.App{Fc},'Color',RGB3,'LineWidth',5);
%                 plot(obj.THApp{Fc}-obj.CP_HardSurface(Fc,1),obj.App{Fc},'Color',RGB4,'LineWidth',5);
%                 % Axes
%                 ax=gca;
%                 ax.LineWidth=2;
%                 ax.FontSize = 16;
%                 ax.XLabel.String = 'Height (m)';
%                 ax.XLabel.FontSize = 20;
%                 ax.YLabel.String = 'Vertical deflection (N)';
%                 ax.YLabel.FontSize = 20;
%                 %% Save figures
%                 %%% Define the name for the figure title
%                 partname=sprintf('-cp-Fc%d',Fc);
%                 fullname=sprintf('%s%s',figname,partname);
%                 %%% Save the current figure in the current folder
%                 print(gcf,fullname,'-dpng');
%                 close all
            end
            obj.CPFlag.HardSurface = 1;
        end
                        
        function fc_xLimit_idx(obj,xLimit1,xLimit2,xLimit3)
            % Finds the corresponding indices of an array of data points
            
            % Limits used for the analysis of TC            
            %xLimit1=20e-9; Retraction limit 1
            %xLimit2=200e-9; Approach limit
            %xLimit3=300e-9; Retraction limit 2
            % Loop
            for Fc=1:obj.NCurves
             %for Fc=20 % For Debugging
             %  sprintf('Force curve No. %d',Fc) % Gives current Force
             %  curve for debugging
                if ~obj.SMFSFlag.Uncorrupt(Fc)     % Exclude corrupted force curves from the analysis     
                    continue
                end
                % Allocate data
                xApp=obj.THApp{Fc}-obj.CP_HardSurface(Fc);
                xRet=obj.THRet{Fc}-obj.CP_HardSurface(Fc);
                % Approach index 
                AppCond1=abs(xApp)<xLimit2; % Determine the elements that fulfil the logical argument           
                AppIdx1=find(AppCond1,1,'first'); % Read out the index of the first cell that fulfil the argument
                % Retraction indices  
                RetCond=abs(xRet)>xLimit1; % Determine the elements that fulfil the logical argument           
                RetIdx1=find(RetCond,1,'first'); % Read out the index of the first cell that fulfil the argument
                RetCond2=abs(xRet)>xLimit3; % Determine the elements that fulfil the logical argument           
                RetIdx2=find(RetCond2,1,'first'); % Read out the index of the first cell that fulfil the argument
                % Allocate data
                obj.AppIdx1(Fc)=AppIdx1;    % 200 nm index of approach data
                obj.RetIdx1(Fc)=RetIdx1;    % 20 nm index of retraction data
                obj.RetIdx2(Fc)=RetIdx2;    % 300 nm index of retraction data
            end     
              %% Appendix
%             close all
%             % Define variables
%             kk=1   
%             % Graphical preview
%             fig=gcf;
%             fig.Units='normalized'; % changes to normalized unit settings, necessary to receive the full screen size in the next line
%             fig.Color='white'; % changes the background color of the figure
%             fig.OuterPosition=[0.5 0 0.5 1];% changes the size of the figure to half screen
%             fig.PaperOrientation='landscape';
%             grid on
%             hold on
%             plot(obj.HHApp{kk}-obj.CP_HardSurface(kk,1),obj.BasedApp{kk},'b');
%             plot(obj.HHRet{kk}-obj.CP_HardSurface(kk,1),obj.BasedRet{kk},'r');
%             plot(obj.HHApp{kk}(AppIdx1,1)-obj.CP_HardSurface(kk,1),obj.BasedApp{kk}(AppIdx1,1),'gd','MarkerFaceColor','g','MarkerSize',20);
%             plot(obj.HHRet{kk}(RetIdx1,1)-obj.CP_HardSurface(kk,1),obj.BasedRet{kk}(RetIdx1,1),'kp','MarkerFaceColor','k','MarkerSize',20);
%             plot(obj.HHRet{kk}(RetIdx2,1)-obj.CP_HardSurface(kk,1),obj.BasedRet{kk}(RetIdx2,1),'kp','MarkerFaceColor','k','MarkerSize',20);
        end
        
        function fc_selection_threshold(obj,AppThreshValue,RetThreshValue)
            % fc_selection_threshold: A function to distinguish between
            % force curves based on prerequirements
            % 
            % SMFSFlag.AppMinCrit = 0:  Flags that force curves show snap-in interaction
            % SMFSFlag.AppMinCrit = 1:  Flags that force curves show no snap-in interaction
            % SMFSFlag.RetMinCrit = 0:  Flags that force curves shows interaction
            % SMFSFlag.RetMinCrit = 1:  Flags that force curves shows no or only tip-substrate interaction
%             if nargin <2
%                 AppThreshValue=15e-12;     % 15 pN 
%                 RetThreshValue=25e-12;     % 25 pN 
%             end
            % Force curves loop
            for Fc=1:obj.NCurves
             %for kk=20 % For Debugging
               % sprintf('Force curve No. %d',kk) % Gives current Force curve for debugging 
                % Allocate data
                yApp=obj.BasedApp{Fc};
                yRet=obj.BasedRet{Fc};
                % Determine the minimum value within the selected range
                AppMin=abs(min(yApp(obj.AppIdx1(Fc):end))); % Finds the max adhesion value within the defined range
                RetMin=abs(min(yRet(obj.RetIdx1(Fc):obj.RetIdx2(Fc)))); % Finds the max adhesion value within the defined range
                % Flags
                if AppMin>AppThreshValue
                    obj.SMFSFlag.AppMinCrit(Fc)=1; % Flags that force curves show snap-in interaction
                else
                    obj.SMFSFlag.AppMinCrit(Fc)=0; % Flags that force curves show no snap-in interaction
                end
                if RetMin>RetThreshValue
                    obj.SMFSFlag.RetMinCrit(Fc)=1; % Flags force curves that show interaction
                else
                    obj.SMFSFlag.RetMinCrit(Fc)=0; % Flags force curves that show no or only tip-substrate interaction
                end
            end
         
        end
        
     
        function fc_snap_in_length_MAD(obj,SelectedFlagStatus)
            % Fct determines the snap-in length using a moving median absolute deviation

            % Preallocate
            obj.SnapInIdx=zeros(1,obj.NCurves);
            obj.SnapInLength=zeros(1,obj.NCurves);
            % Define variables
            WindowBeforePercentage=0.1;
            WindowAfterPercentage=0.1;
            DataShareStartApp=0.65;
            DataShareEndApp=0.15;
            LimitFactor=1.5;
            %% Loop over all force curves
            for Fc=1:obj.NCurves
                % Debugging
                %for Fc=74
                sprintf('Force curve No. %d',Fc) % Gives current
                % Force curve for debugging
                %% Force curve selection criteria
                switch SelectedFlagStatus
                case 'Pre'
                if ~obj.SMFSFlag.Uncorrupt(Fc) || ~obj.SMFSFlag.AppMinCrit(Fc)     % Condition before the first analysis: Exclude corrupted force curves or force curves showing no snap-in from the analysis
                continue
                end
                case 'Post'
                if ~obj.SMFSFlag.Uncorrupt(Fc) || ~obj.SMFSFlag.Selected(Fc) || ~obj.SMFSFlag.AppMinCrit(Fc)     % Condition if FM Flag Selected has been set: Exclude corrupted force curves or force curves showing no snap-in from the analysis
                continue
                end
                end                   
                % Allocate data
                xApp=obj.THApp{Fc}-obj.CP_HardSurface(Fc);
                yApp=obj.BasedApp{Fc};
                % Determine the snap-in index
                DataPtsApp=size(yApp); % Determine the amount of data points in the force curve
                LimitIdxApp1=round(DataPtsApp(1)-DataPtsApp(1)*DataShareStartApp); % Determine the corresponding index
                LimitIdxApp2=round(DataPtsApp(1)-DataPtsApp(1)*DataShareEndApp); % Determine the corresponding index
                WindowBeforePts=round(length(yApp)*WindowBeforePercentage);
                WindowAfterPts=round(length(yApp)*WindowAfterPercentage);
                MovWindow=[WindowBeforePts WindowAfterPts]; % moving data point window for the moving median absolute deviation
                NormApp=normalize(flip(yApp));  % Normalize and flip the data
                NormAppMAD = movmad(NormApp,MovWindow); % Apply the moving median absolute deviation
                NormAppSelMAD=NormAppMAD;
                NormAppSelMAD([1:LimitIdxApp1, LimitIdxApp2:end])=[]; % Select defined range of data points
                MaxSelMAD=max(NormAppSelMAD);  % Find the maximun value of the selection
                PeakSelIdx=find(NormAppSelMAD==MaxSelMAD,1,'first'); % Find the index of the maximum value in the selection
                PeakSelIdx=PeakSelIdx+LimitIdxApp1; % Correct for the selection to determine the corresponding index in the whole data set
                yLimit=MaxSelMAD*LimitFactor; % Define selection limit based on the maximum peak of the selection - representing the background
                PeakIdx=find(NormAppMAD>yLimit,1,'last'); % Find the index in the data fulfilling the condition
                if ~isempty(PeakIdx)  
                    % Flag
                    obj.SMFSFlag.SnapIn(Fc)=1;
                    % Allocate data
                    obj.SnapInIdx(Fc)=length(yApp)-PeakIdx; % Correct for the data the peak index is based on by substracting from the number of data points
                    obj.SnapInLength(Fc)=abs(xApp(obj.SnapInIdx(Fc))); % Corresponding x-value of the index
                    %% Plot condition
                    if  ~obj.DebugFlag.Plot % Suppress plotting
                        %if  obj.DebugFlag.Plot % Allow plotting
                        continue
                    end
                    % Define variables for the figure name
                    ExtendVelocityConvert=num2str(obj.ExtendVelocity*1e9);
                    RetractVelocityConvert=num2str(obj.RetractVelocity*1e9);
                    HoldingTimeConvert=num2str(obj.HoldingTime);
                    % Classification criteria
                    figname=strcat(obj.Date,{'_'},obj.Time,{'_'},obj.ID,{'_'},obj.Substrate,{'_'},obj.EnvCond,{'_'},obj.Linker,{'_'},obj.Chipbox,{'_'},obj.ChipCant,{'_'},ExtendVelocityConvert,{'_'},RetractVelocityConvert,{'_'},HoldingTimeConvert);
                    figname=char(figname);
                    %% Figure
                    h_fig=figure(Fc);
                    h_fig.Color='white'; % changes the background color of the figure
                    h_fig.Units='normalized'; % Defines the units
                    h_fig.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
                    h_fig.PaperOrientation='landscape';
                    subplot(3,1,1)
                    hold on
                    plot(xApp,yApp)
                    plot(xApp(obj.SnapInIdx(Fc)),yApp(obj.SnapInIdx(Fc)),'rx','MarkerSize',20)
                    hold off
                    title('Force-distance curve')
                    ax1=gca;
                    ax1.XLim = [-inf inf];
                    ax1.YLim = [-inf inf];
                    ax1.XLabel.String = 'Tip-sample seperation  (m)';
                    ax1.YLabel.String = 'Force (N)';
                    subplot(3,1,2)
                    hold on
                    plot(NormApp)
                    plot(PeakIdx,NormApp(PeakIdx),'rx','MarkerSize',20)
                    plot(PeakSelIdx,NormApp(PeakSelIdx),'g*','MarkerSize',20)
                    hold off
                    title('flipped normalized approach data')
                    ax2=gca;
                    ax2.XLim = [0 inf];
                    ax2.YLim = [-inf 1];
                    ax2.XLabel.String = 'Index (1)';
                    ax2.YLabel.String = '';
                    subplot(3,1,3)
                    hold on
                    plot(NormAppMAD)
                    plot(PeakIdx,NormAppMAD(PeakIdx),'rx','MarkerSize',20)
                    plot(PeakSelIdx,NormAppMAD(PeakSelIdx),'g*','MarkerSize',20)
                    hold off
                    title('MAD data')
                    ax3=gca;
                    ax3.XLim = [0 inf];
                    ax3.YLim = [0 0.5];
                    ax3.XLabel.String = 'Index (1)';
                    ax3.YLabel.String = 'Median deviation';
                    %% Save figures
                    %%% Define the name for the figure title
                    partname=sprintf('-Fc %d',Fc);
                    % fullname=sprintf('%s%s',figname,partname);
                    fullname=sprintf('%s%s',figname,partname);
                    %%% Save the current figure in the current folder
                    print(gcf,fullname,'-dpng');
                    continue
                end
            end
            close all
            % Allocate data
            obj.FMSnapInMean=mean(obj.SnapInLength);
            obj.FMSnapInMin=min(obj.SnapInLength);
            obj.FMSnapInMax=max(obj.SnapInLength);
        end

        function fc_pulling_length_MAD(obj,SelectedFlagStatus)
            
            % Preallocate
            obj.PullingLengthIdx=zeros(1,obj.NCurves);
            obj.PullingLength=zeros(1,obj.NCurves);
            % Define variables
            DataShareStartRet=0.02; %
            DataShareEndRet=0.07; %
            WindowBeforePercentage=0.01;
            WindowAfterPercentage=0.01;
            LimitFactor=2;
            %% Loop over all force curves
            for Fc=1:obj.NCurves
            % Debugging
            %for fc=7
                % sprintf('Force curve No. %d',jj) % Gives current Force curve
                %% Force curve selection criteria
                switch SelectedFlagStatus
                case 'Pre'
                if ~obj.SMFSFlag.Uncorrupt(Fc) || ~obj.SMFSFlag.RetMinCrit(Fc)     % Condition before the first analysis: Exclude corrupted force curves or force curves showing no snap-in from the analysis
                continue
                end
                case 'Post'
                if ~obj.SMFSFlag.Uncorrupt(Fc) || ~obj.SMFSFlag.Selected(Fc) || ~obj.SMFSFlag.RetMinCrit(Fc)     % Condition if FM Flag Selected has been set: Exclude corrupted force curves or force curves showing no snap-in from the analysis
                continue
                end
                end
                % Allocate data
                xRet=obj.THRet{Fc}-obj.CP_HardSurface(Fc);
                yRet=obj.BasedRet{Fc};                
                % Determine the pulling length index                
                DataPtsRet=size(yRet); % Determine the amount of data points in the force curve
                LimitIdxRet1=round(DataPtsRet(1)*DataShareStartRet); % Determine the corresponding index
                LimitIdxRet2=round(DataPtsRet(1)*DataShareEndRet); % Determine the corresponding index
                BeforePts=round(length(yRet)*WindowBeforePercentage); % Before data points of the moving window
                AfterPts=round(length(yRet)*WindowAfterPercentage); % After data points of the moving window
                MovWindow=[BeforePts AfterPts]; % moving data point window for the moving median absolute deviation
                NormRet=normalize(flip(yRet));  % Normalize and flip the data
                NormRetMAD=movmad(NormRet,MovWindow); % Apply the moving median absolute deviation
                NormRetSelMAD=NormRetMAD;      
                NormRetSelMAD([1:LimitIdxRet1, LimitIdxRet2:end])=[]; % Select defined range of data points
                MaxSelMAD=max(NormRetSelMAD);  % Find the maximun value of the selection
                PeakSelIdx=find(NormRetSelMAD==MaxSelMAD,1,'first'); % Find the index of the maximum value in the selection
                PeakSelIdx=PeakSelIdx+LimitIdxRet1; % Correct for the selection to determine the corresponding index in the whole data set
                yLimit=MaxSelMAD*LimitFactor; % Define selection limit based on the maximum peak of the selection - representing the background
                PeakIdx=find(NormRetMAD>yLimit,1,'first'); % Find the pulling length index
                % Allocate data
                obj.PullingLengthIdx(Fc)=length(yRet)-PeakIdx; % Correct for the flipped data the peak index is based on by substracting from the number of data points
                obj.PullingLength(Fc)=abs(xRet(obj.PullingLengthIdx(Fc))); % Corresponding x-value of the index
                % Flags
                if obj.RetIdx1(Fc)>obj.PullingLengthIdx(Fc)
                    obj.SMFSFlag.LengthRequisite(Fc)=0;
                else
                    obj.SMFSFlag.LengthRequisite(Fc)=1;
                end
                obj.SMFSFlag.PullingLength(Fc)=1;                
                %% Plot condition
                if  ~obj.DebugFlag.Plot % Suppress plotting
                %if  obj.DebugFlag.Plot % Allow plotting
                    continue
                end
                % Define variables for the figure name
                ExtendVelocityConvert=num2str(obj.ExtendVelocity*1e9);
                RetractVelocityConvert=num2str(obj.RetractVelocity*1e9);
                HoldingTimeConvert=num2str(obj.HoldingTime);
                % Classification criteria
                figname=strcat(obj.Date,{'_'},obj.Time,{'_'},obj.ID,{'_'},obj.Substrate,{'_'},obj.EnvCond,{'_'},obj.Linker,{'_'},obj.Chipbox,{'_'},obj.ChipCant,{'_'},ExtendVelocityConvert,{'_'},RetractVelocityConvert,{'_'},HoldingTimeConvert);
                figname=char(figname);
                %% Figure
                h_fig=figure(Fc);
                h_fig.Color='white'; % changes the background color of the figure
                h_fig.Units='normalized'; % Defines the units
                h_fig.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
                h_fig.PaperOrientation='landscape';
                subplot(3,1,1)
                hold on
                plot(xRet,yRet)
                plot(xRet(obj.PullingLengthIdx(Fc)),yRet(obj.PullingLengthIdx(Fc)),'rx','MarkerSize',20)
                hold off
                title('Force-distance curve')
                ax1=gca;
                ax1.XLim = [-inf inf];
                ax1.YLim = [-inf inf];
                ax1.XLabel.String = 'Tip-sample seperation  (m)';
                ax1.YLabel.String = 'Force (N)';
                subplot(3,1,2)
                hold on
                plot(NormRet)
                plot(PeakIdx,NormRet(PeakIdx),'rx','MarkerSize',20)
                plot(PeakSelIdx,NormRet(PeakSelIdx),'g*','MarkerSize',20)
                hold off
                title('flipped normalized retention data')
                ax2=gca;
                ax2.XLim = [0 inf];
                ax2.YLim = [-inf 1];
                ax2.XLabel.String = 'Index (1)';
                ax2.YLabel.String = '';
                subplot(3,1,3)
                hold on
                plot(NormRetMAD)
                plot(PeakIdx,NormRetMAD(PeakIdx),'rx','MarkerSize',20)
                plot(PeakSelIdx,NormRetMAD(PeakSelIdx),'g*','MarkerSize',20)
                hold off
                title('MAD data')
                ax3=gca;
                ax3.XLim = [0 inf];
                ax3.YLim = [0 0.15];
                ax3.XLabel.String = 'Index (1)';
                ax3.YLabel.String = 'Median deviation';
                %% Save figures
                %%% Define the name for the figure title
                partname=sprintf('-Fc %d',Fc);
                % fullname=sprintf('%s%s',figname,partname);
                fullname=sprintf('%s%s',figname,partname);
                %%% Save the current figure in the current folder
                print(gcf,fullname,'-dpng');             
            end
            close all
            % Allocate data
            obj.FMPullingLengthMean=mean(obj.PullingLength);
            obj.FMPullingLengthMin=min(obj.PullingLength);
            obj.FMPullingLengthMax=max(obj.PullingLength);
        end
                
        function fc_adh_force_max(obj,SelectedFlagStatus)
            % Function to find the maximum adhesion force value in the
            % approach and retraction data of a force curve
            
            % Preallocate
            obj.AdhForceMaxApp=zeros(1,obj.NCurves);
            obj.AdhForceMaxAppIdx=zeros(1,obj.NCurves);
            obj.AdhForceMaxRet=zeros(1,obj.NCurves);
            obj.AdhForceMaxRetIdx=zeros(1,obj.NCurves);
            obj.AdhForceUnbinding=zeros(1,obj.NCurves);
            obj.AdhForceUnbindingIdx=zeros(1,obj.NCurves);
            % Define variables            
            RGB_A1=[0 25 255]./255;  % Blue 
            RGB_A3=[255 102 0]./255; % Orange
            RGB_A5=[255 0 26]./255; % Red
            xDistance=15e-9;
            %% %% Loop over all force curves
            for Fc=1:obj.NCurves
            % Debugging
            %for Fc=7 % for debugging
             sprintf('Force curve No. %d',Fc) % Gives current Force curve
            % for debugging
            %% Force curve selection criteria
                switch SelectedFlagStatus
                case 'Pre'
                if ~obj.SMFSFlag.Uncorrupt(Fc) || ~obj.SMFSFlag.RetMinCrit(Fc) || ~obj.SMFSFlag.LengthRequisite(Fc)    % Condition before the first analysis: Exclude corrupted force curves or force curves showing no snap-in from the analysis
                continue
                end
                case 'Post'
                if ~obj.SMFSFlag.Uncorrupt(Fc) || ~obj.SMFSFlag.Selected(Fc) || ~obj.SMFSFlag.RetMinCrit(Fc) || ~obj.SMFSFlag.LengthRequisite(Fc)    % Condition if FM Flag Selected has been set: Exclude corrupted force curves or force curves showing no snap-in from the analysis
                continue
                end
                end                   
            % Allocate data
            xApp=obj.THApp{Fc}-obj.CP_HardSurface(Fc);
            xRet=obj.THRet{Fc}-obj.CP_HardSurface(Fc);            
            yApp=obj.BasedApp{Fc};
            yRet=obj.BasedRet{Fc};         
            % Define start index for the unbinding event 
            xUnbdingBoundary=obj.PullingLength(Fc)-xDistance;
            [~,obj.UnbindingBoundaryIdx(Fc)]=min(abs(abs(xRet)-xUnbdingBoundary));     
            %% Determine the maximum adhesion force
            % Approach
            if obj.SMFSFlag.SnapIn(Fc)==1 && obj.SMFSFlag.AppMinCrit(Fc)% Select only force curves that show a snap-in behavior
                obj.AdhForceMaxApp(Fc)=min(yApp(obj.SnapInIdx(Fc):end)); % Determine maximum adhesion forces from the pulling length index to the last data point
                obj.AdhForceMaxAppIdx(Fc)=find(yApp==obj.AdhForceMaxApp(Fc)); % Finds the corresponding index of the value
            end 
            % Retract
            obj.AdhForceMaxRet(Fc)=min(yRet(obj.RetIdx1(Fc):obj.PullingLengthIdx(Fc))); % Determine maximum adhesion forces from the 50nm index to the pulling length index
            obj.AdhForceMaxRetIdx(Fc)=find(yRet==obj.AdhForceMaxRet(Fc)); % Finds the index of the value that fulfils the condition               
            % Last unbinding event on the retract curve
            obj.AdhForceUnbinding(Fc)=min(yRet(obj.UnbindingBoundaryIdx(Fc):obj.PullingLengthIdx(Fc))); % Determine maximum adhesion forces close to the pulling length                         
            AdhForceUnbindingSelIdx=find(yRet(obj.UnbindingBoundaryIdx(Fc):obj.PullingLengthIdx(Fc))==obj.AdhForceUnbinding(Fc)); % Finds the index of the value that fulfils the condition
            obj.AdhForceUnbindingIdx(Fc)=obj.UnbindingBoundaryIdx(Fc)+AdhForceUnbindingSelIdx; % Correct for the shifted starting position
            %% Plot condition
                if  ~obj.DebugFlag.Plot % Suppress plotting
                %if  obj.DebugFlag.Plot % Allow plotting
                    continue
                end
                % Define variables for the figure name
                ExtendVelocityConvert=num2str(obj.ExtendVelocity*1e9);
                RetractVelocityConvert=num2str(obj.RetractVelocity*1e9);
                HoldingTimeConvert=num2str(obj.HoldingTime);
                % Classification criteria
                figname=strcat(obj.Date,{'_'},obj.Time,{'_'},obj.ID,{'_'},obj.Substrate,{'_'},obj.EnvCond,{'_'},obj.Linker,{'_'},obj.Chipbox,{'_'},obj.ChipCant,{'_'},ExtendVelocityConvert,{'_'},RetractVelocityConvert,{'_'},HoldingTimeConvert);
                figname=char(figname);
                % Figure
                h_fig=figure(Fc);
                h_fig.Color='white'; % changes the background color of the figure
                h_fig.Units='normalized'; % Defines the units
                h_fig.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
                h_fig.PaperOrientation='landscape';
                hold on
                plot(xApp,yApp,'Color',RGB_A1);
                plot(xRet,yRet,'Color',RGB_A3);
                plot(xApp(obj.AdhForceMaxAppIdx(Fc)),yApp(obj.AdhForceMaxAppIdx(Fc)),'h','MarkerSize',10,'MarkerFaceColor',RGB_A1,'MarkerEdgeColor',RGB_A1)
                plot(xRet(obj.AdhForceMaxRetIdx(Fc)),yRet(obj.AdhForceMaxRetIdx(Fc)),'h','MarkerSize',10,'MarkerFaceColor', RGB_A3,'MarkerEdgeColor',RGB_A3)
                plot(xRet(obj.AdhForceUnbindingIdx(Fc)),yRet(obj.AdhForceUnbindingIdx(Fc)),'p','MarkerSize',12,'MarkerFaceColor',RGB_A5,'MarkerEdgeColor',RGB_A5)
                hold off
                title('Force-distance curve')
                ax1=gca;
                ax1.XLim = [-inf inf];
                ax1.YLim = [-inf inf];
                ax1.XLabel.String = 'Tip-sample seperation  (m)';
                ax1.YLabel.String = 'Force (N)';
                %% Save figures
                %%% Define the name for the figure title
                partname=sprintf('-Fc %d',Fc);
                % fullname=sprintf('%s%s',figname,partname);
                fullname=sprintf('%s%s',figname,partname);
                %%% Save the current figure in the current folder
                print(gcf,fullname,'-dpng');                 
            end
        close all
        end
       
        function fc_adhesion_energy_idxlength(obj,SelectedFlagStatus)
            % Determine the adhesion energy using the previous defined pulling length index
                % Preallocate
                obj.yAppLim=cell(1,obj.NCurves);
                obj.AppAdhEnergy_IdxMethod=zeros(1,obj.NCurves);
                obj.yRetLim=cell(1,obj.NCurves);
                obj.RetAdhEnergy_IdxMethod=zeros(1,obj.NCurves);
                IntApp=zeros(1,obj.NCurves);
                IntRet=zeros(1,obj.NCurves);
                % Define variables
                limit1=0;   % Define the limit
            %% Loop over all force curves
                for Fc=1:obj.NCurves
                %for Fc=97 % % For debugging and testing
                %% Force curve selection criteria
                switch SelectedFlagStatus
                case 'Pre'
                if ~obj.SMFSFlag.Uncorrupt(Fc) || ~obj.SMFSFlag.RetMinCrit(Fc) || ~obj.SMFSFlag.LengthRequisite(Fc)    % Condition before the first analysis: Exclude corrupted force curves or force curves showing no snap-in from the analysis
                continue
                end
                case 'Post'
                if ~obj.SMFSFlag.Uncorrupt(Fc) || ~obj.SMFSFlag.Selected(Fc) || ~obj.SMFSFlag.RetMinCrit(Fc) || ~obj.SMFSFlag.LengthRequisite(Fc)    % Condition if FM Flag Selected has been set: Exclude corrupted force curves or force curves showing no snap-in from the analysis
                continue
                end
                end         
                % Allocate data
                xRet=obj.THRet{Fc}-obj.CP_HardSurface(Fc); % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method 
                xApp=obj.THApp{Fc}-obj.CP_HardSurface(Fc); % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method 
                yApp=obj.BasedApp{Fc};
                yRet=obj.BasedRet{Fc};                
                % Apply the limit                
                    if obj.SMFSFlag.SnapIn(Fc)==1                                    
                    % Approach
                    yApp(yApp>limit1)=0; % Set all data points representing surface contact (set point related, above the y-axis 0 line) to 0
                    yApp(1:obj.SnapInIdx(Fc))=0; % Set all data points with a higher index (surface distance is higher) than the snap-in index to 0
                    % Allocate data
                    obj.yAppLim{Fc}=yApp; 
                    % Determine the adhesion energy
                    IntApp(Fc)=trapz(yApp,xApp); % Integrates over the modified y-retraction data with respect to the corresponding x-retraction data 
                    obj.AppAdhEnergy_IdxMethod(Fc)=IntApp(Fc);
                    end
                    % Retention 
                    if obj.SMFSFlag.PullingLength(Fc)==1                       
                    yRet(yRet>limit1)=0;  % Set all data points representing surface contact (set point related, above the y-axis 0 line) to 0
                    yRet(obj.PullingLengthIdx(Fc):end)=0; % Set all data points with a higher index (surface distance is higher) than the pulling length index to 0
                    % Allocate data
                    obj.yRetLim{Fc}=yRet; 
                    % Determine the adhesion energy
                    IntRet(Fc)=trapz(yRet,xRet); % Integrates over the modified y-retraction data with respect to the corresponding x-retraction data 
                    obj.RetAdhEnergy_IdxMethod(Fc)=IntRet(Fc);
                    end
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
        
        function fc_print_analysed(obj,XMin,XMax,YMin,YMax,NumFcMax) % fc ... force curve
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
            Res=[1 1 2560 1440]; % Define the figure resolution
            % Define variables for the figure name
            ExtendVelocityConvert=num2str(round(obj.ExtendVelocity*1e9));
            RetractVelocityConvert=num2str(round(obj.RetractVelocity*1e9));
            HoldingTimeConvert=num2str(round(obj.HoldingTime));        
            % Classification criteria
            figname=strcat(obj.Date,{'_'},obj.Time,{'_'},obj.ID,{'_'},obj.Substrate,{'_'},obj.EnvCond,{'_'},obj.Linker,{'_'},obj.Chipbox,{'_'},obj.ChipCant,{'_'},ExtendVelocityConvert,{'_'},RetractVelocityConvert,{'_'},HoldingTimeConvert);
            figname=char(figname);
            % Parse unit scale function
            [Xmultiplier,Xunit,~] = AFMImage.parse_unit_scale(1e+9,'nm',1);
            [Ymultiplier,Yunit,~] = AFMImage.parse_unit_scale(1e+9,'nN',1);
            %% Define variables for the plot loop
            mm=ceil(sqrt(NumFcMax)); % Variable defines the tiledlayout and thereby the tile arrangement
            nn=mm; % Variable defines the tiledlayout and thereby the tile arrangement
            ww=1; % Variable used to loop through the flag based while loop
            DiffFc=0; % Variable to correct for differences of the variables between the plot loop and flag based while loop
            NumFcUncorrupt=nnz(obj.SMFSFlag.Selected.*obj.SMFSFlag.Uncorrupt.*obj.SMFSFlag.RetMinCrit.*obj.SMFSFlag.LengthRequisite); % Determine the number of force curves that could been analysed
            NumFigures=ceil(NumFcUncorrupt/NumFcMax);
            if NumFigures==0     % If condition is fulfilled stop function and return to calling function
                return
            end
            RemainderMax=mod(NumFcUncorrupt,NumFcMax); % Check for remainder
            if RemainderMax ~= 0
                % Additional variables if there is a remainder
                oo=round(sqrt(RemainderMax)); % Determine the number of rows in the figure
                pp=ceil(sqrt(RemainderMax)); % Determine the number of columns in the figure
                RemainderReal=mod(NumFcUncorrupt,oo*pp); % Correct the remainder based on the determined rows times columns
            end
            %% figure loop
            for kk=1:NumFigures
                % if  ~obj.DebugFlag.Plot % Suppress plotting
                if  obj.DebugFlag.Plot % Allow plotting
                    continue
                end
                % Figure
                h_fig=figure(kk);
                h_fig.Color='white'; % changes the background color of the figure
                h_fig.Units='pixel'; % Defines the units
                h_fig.OuterPosition=Res;
                h_fig.PaperOrientation='landscape';
                h_fig.Name=figname;
                %% Verify the remainder condition
                if RemainderMax == 0 % No remainder existing
                    t = tiledlayout(mm,nn);
                    t.TileSpacing = 'none'; % To reduce the spacing between the tiles
                    t.Padding = 'none'; % To reduce the padding of perimeter of a tile
                    if kk==1
                        NumFcCorSelec(kk)=nnz(~obj.SMFSFlag.Uncorrupt(1:NumFcMax));
                    elseif kk==2
                        NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt((NumFcMax+1):(NumFcMax*(kk))));
                    else
                        NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt((NumFcMax*(kk-1)+1):(NumFcMax*(kk))));
                    end
                    
                    if kk==1
                        qq=1;
                    else
                        qq=(kk-1)*mm*nn+1; % Determines the index of the next force curve to be plotted on the next figure
                    end
                    %% Plot loop
                    for qq=qq:obj.NCurves % Loop over all force curves in the force map                       
                        if ww<qq+DiffFc
                            ww=qq+DiffFc;
                        end
                        % Flag based while loop
                        while ~obj.SMFSFlag.Selected(ww) || ~obj.SMFSFlag.Uncorrupt(ww) || ~obj.SMFSFlag.RetMinCrit(ww) || ~obj.SMFSFlag.LengthRequisite(ww) || ~obj.SMFSFlag.Selected(ww)  % Stay in the while loop as long as the entry is zero
                            ww=ww+1;                            
                            if ww>qq
                                DiffFc=ww-qq;
                            end
                            if ww>obj.NCurves
                                break
                            end
                        end
                        if ww>obj.NCurves
                                break
                        end
                        %% Allocate data
                        xApp=(obj.THApp{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                        xRet=(obj.THRet{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                        yApp=obj.BasedApp{qq+DiffFc}/Ymultiplier;
                        yRet=obj.BasedRet{qq+DiffFc}/Ymultiplier;    
                        yRetLim=obj.yRetLim{qq+DiffFc}/Ymultiplier;
                        yAppLim=obj.yAppLim{qq+DiffFc}/Ymultiplier;                 
                        % if condition with flag based while loop variable
                        % ww > plot loop variable qq
                        if ww>qq
                            ax=nexttile;
                            ax.XLim = [XMin XMax];
                            ax.YLim = [YMin YMax];
                                if obj.SMFSFlag.PullingLength(qq+DiffFc)==1 && obj.SMFSFlag.SnapIn(qq+DiffFc)==1 % Force curve posses pulling length and snap-in variables
                                    grid on
                                    hold on
                                    area(xRet(1:obj.PullingLengthIdx(qq+DiffFc)),yRetLim(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                    area(xApp(obj.SnapInIdx(qq+DiffFc):end),yAppLim(obj.SnapInIdx(qq+DiffFc):end),'FaceColor',RGB7)
                                    plot(xApp,yApp,'Color',RGB1);
                                    plot(xRet,yRet,'Color',RGB2);
                                    plot(xRet(obj.PullingLengthIdx(qq+DiffFc)),yRet(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)
                                    plot(xApp(obj.SnapInIdx(qq+DiffFc)),yApp(obj.SnapInIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB10,'MarkerEdgeColor',RGB10)
                                    plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                    plot(xApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),yApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)
                                    plot(xRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),yRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)
                                elseif obj.SMFSFlag.PullingLength(qq+DiffFc)==1 && obj.SMFSFlag.SnapIn(qq+DiffFc)==0 % Force curve posses pulling length variables but no snap-in variables
                                    grid on
                                    hold on
                                    area(xRet(1:obj.PullingLengthIdx(qq+DiffFc)),yRetLim(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                    plot(xApp,yApp,'Color',RGB1);
                                    plot(xRet,yRet,'Color',RGB2);
                                    plot(xRet(obj.PullingLengthIdx(qq+DiffFc)),yRet(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)
                                    plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                    plot(xRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),yRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)
                                elseif obj.SMFSFlag.PullingLength(qq+DiffFc)==0 && obj.SMFSFlag.SnapIn(qq+DiffFc)==0 % Force curve posses neither pulling length variables nor snap-in variables
                                    grid on
                                    hold on
                                    plot(xApp,yApp,'Color',RGB1);
                                    plot(xRet,yRet,'Color',RGB2);
                                    plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                    plot(xApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),yApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)
                                end                 
                            % Title for each Subplot
                            ti=title(sprintf('%i',qq+DiffFc),'Color','k');
                            ti.Units='normalized'; % Set units to 'normalized'
                            ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                        else
                            ax=nexttile;
                            ax.XLim = [XMin XMax];
                            ax.YLim = [YMin YMax];
                                if obj.SMFSFlag.PullingLength(qq+DiffFc)==1 && obj.SMFSFlag.SnapIn(qq+DiffFc)==1 % Force curve posses pulling length and snap-in variables
                                    grid on
                                    hold on
                                    area(xRet(1:obj.PullingLengthIdx(qq+DiffFc)),yRetLim(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                    area(xApp(obj.SnapInIdx(qq+DiffFc):end),yAppLim(obj.SnapInIdx(qq+DiffFc):end),'FaceColor',RGB7)
                                    plot(xApp,yApp,'Color',RGB1);
                                    plot(xRet,yRet,'Color',RGB2);
                                    plot(xRet(obj.PullingLengthIdx(qq+DiffFc)),yRet(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)
                                    plot(xApp(obj.SnapInIdx(qq+DiffFc)),yApp(obj.SnapInIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB10,'MarkerEdgeColor',RGB10)
                                    plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                    plot(xApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),yApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)
                                    plot(xRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),yRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)
                                elseif obj.SMFSFlag.PullingLength(qq+DiffFc)==1 && obj.SMFSFlag.SnapIn(qq+DiffFc)==0 % Force curve posses pulling length variables but no snap-in variables
                                    grid on
                                    hold on
                                    area(xRet(1:obj.PullingLengthIdx(qq+DiffFc)),yRetLim(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                    plot(xApp,yApp,'Color',RGB1);
                                    plot(xRet,yRet,'Color',RGB2);
                                    plot(xRet(obj.PullingLengthIdx(qq+DiffFc)),yRet(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)
                                    plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                    plot(xRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),yRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)
                                elseif obj.SMFSFlag.PullingLength(qq+DiffFc)==0 && obj.SMFSFlag.SnapIn(qq+DiffFc)==0 % Force curve posses neither pulling length variables nor snap-in variables
                                    grid on
                                    hold on
                                    plot(xApp,yApp,'Color',RGB1);
                                    plot(xRet,yRet,'Color',RGB2);
                                    plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                    plot(xApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),yApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)
                                end
                            % Title for each Subplot
                            ti=title(sprintf('%i',qq+DiffFc),'Color','k');
                            ti.Units='normalized'; % Set units to 'normalized'
                            ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                        end
                        if qq == NumFcMax*kk
                            break
                        end
                        if ww==obj.NCurves
                                break
                        end
                    end
                    
                else % Remainder exisiting
                    %% Last figure conditions
                    if kk~=NumFigures
                        t = tiledlayout(mm,nn);
                        t.TileSpacing = 'none'; % To reduce the spacing between the tiles
                        t.Padding = 'none'; % To reduce the padding of perimeter of a tile
                        if kk==1
                            NumFcCorSelec(kk)=nnz(~obj.SMFSFlag.Uncorrupt(1:NumFcMax));
                        elseif kk==2
                            NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt((NumFcMax+1):(NumFcMax*(kk))));
                        else
                            NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt((NumFcMax*(kk-1)+1):(NumFcMax*(kk))));
                        end
                        
                        if kk==1
                            qq=1;
                        else
                            qq=(kk-1)*mm*nn+1; % Determines the index of the next force curve to be plotted on the next figure
                        end
                        %% Plot loop
                        for qq=qq:obj.NCurves % Loop over all force curves in the force map
                            if ww<qq+DiffFc
                                ww=qq+DiffFc;
                            end
                            while ~obj.SMFSFlag.Selected(ww) || ~obj.SMFSFlag.Uncorrupt(ww) || ~obj.SMFSFlag.RetMinCrit(ww) || ~obj.SMFSFlag.LengthRequisite(ww)    % Stay in the while loop as long as the entry is zero
                                ww=ww+1;
                                if ww>qq
                                    DiffFc=ww-qq;
                                end
                            end
                            %% Allocate data
                            xApp=(obj.THApp{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                            xRet=(obj.THRet{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                            yApp=obj.BasedApp{qq+DiffFc}/Ymultiplier;
                            yRet=obj.BasedRet{qq+DiffFc}/Ymultiplier;    
                            yRetLim=obj.yRetLim{qq+DiffFc}/Ymultiplier;
                            yAppLim=obj.yAppLim{qq+DiffFc}/Ymultiplier;     
                            if ww>qq
                                ax=nexttile;
                                ax.XLim = [XMin XMax];
                                ax.YLim = [YMin YMax];
                                    if obj.SMFSFlag.PullingLength(qq+DiffFc)==1 && obj.SMFSFlag.SnapIn(qq+DiffFc)==1 % Force curve posses pulling length and snap-in variables
                                        grid on
                                        hold on
                                        area(xRet(1:obj.PullingLengthIdx(qq+DiffFc)),yRetLim(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                        area(xApp(obj.SnapInIdx(qq+DiffFc):end),yAppLim(obj.SnapInIdx(qq+DiffFc):end),'FaceColor',RGB7)
                                        plot(xApp,yApp,'Color',RGB1);
                                        plot(xRet,yRet,'Color',RGB2);
                                        plot(xRet(obj.PullingLengthIdx(qq+DiffFc)),yRet(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)
                                        plot(xApp(obj.SnapInIdx(qq+DiffFc)),yApp(obj.SnapInIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB10,'MarkerEdgeColor',RGB10)
                                        plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                        plot(xApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),yApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)
                                        plot(xRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),yRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)
                                    elseif obj.SMFSFlag.PullingLength(qq+DiffFc)==1 && obj.SMFSFlag.SnapIn(qq+DiffFc)==0 % Force curve posses pulling length variables but no snap-in variables
                                        grid on
                                        hold on
                                        area(xRet(1:obj.PullingLengthIdx(qq+DiffFc)),yRetLim(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                        plot(xApp,yApp,'Color',RGB1);
                                        plot(xRet,yRet,'Color',RGB2);
                                        plot(xRet(obj.PullingLengthIdx(qq+DiffFc)),yRet(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)
                                        plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                        plot(xRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),yRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)
                                    elseif obj.SMFSFlag.PullingLength(qq+DiffFc)==0 && obj.SMFSFlag.SnapIn(qq+DiffFc)==0 % Force curve posses neither pulling length variables nor snap-in variables
                                        grid on
                                        hold on
                                        plot(xApp,yApp,'Color',RGB1);
                                        plot(xRet,yRet,'Color',RGB2);
                                        plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                        plot(xApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),yApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)
                                    end
                                % Title for each Subplot
                                ti=title(sprintf('%i',qq+DiffFc),'Color','k');
                                ti.Units='normalized'; % Set units to 'normalized'
                                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                            else
                                ax=nexttile;
                                ax.XLim = [XMin XMax];
                                ax.YLim = [YMin YMax];
                                    if obj.SMFSFlag.PullingLength(qq+DiffFc)==1 && obj.SMFSFlag.SnapIn(qq+DiffFc)==1 % Force curve posses pulling length and snap-in variables
                                        grid on
                                        hold on
                                        area(xRet(1:obj.PullingLengthIdx(qq+DiffFc)),yRetLim(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                        area(xApp(obj.SnapInIdx(qq+DiffFc):end),yAppLim(obj.SnapInIdx(qq+DiffFc):end),'FaceColor',RGB7)
                                        plot(xApp,yApp,'Color',RGB1);
                                        plot(xRet,yRet,'Color',RGB2);
                                        plot(xRet(obj.PullingLengthIdx(qq+DiffFc)),yRet(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)
                                        plot(xApp(obj.SnapInIdx(qq+DiffFc)),yApp(obj.SnapInIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB10,'MarkerEdgeColor',RGB10)
                                        plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                        plot(xApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),yApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)
                                        plot(xRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),yRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)
                                    elseif obj.SMFSFlag.PullingLength(qq+DiffFc)==1 && obj.SMFSFlag.SnapIn(qq+DiffFc)==0 % Force curve posses pulling length variables but no snap-in variables
                                        grid on
                                        hold on
                                        area(xRet(1:obj.PullingLengthIdx(qq+DiffFc)),yRetLim(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                        plot(xApp,yApp,'Color',RGB1);
                                        plot(xRet,yRet,'Color',RGB2);
                                        plot(xRet(obj.PullingLengthIdx(qq+DiffFc)),yRet(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)
                                        plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                        plot(xRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),yRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)
                                    elseif obj.SMFSFlag.PullingLength(qq+DiffFc)==0 && obj.SMFSFlag.SnapIn(qq+DiffFc)==0 % Force curve posses neither pulling length variables nor snap-in variables
                                        grid on
                                        hold on
                                        plot(xApp,yApp,'Color',RGB1);
                                        plot(xRet,yRet,'Color',RGB2);
                                        plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                        plot(xApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),yApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)
                                    end
                                % Title for each Subplot
                                ti=title(sprintf('%i',qq+DiffFc),'Color','k');
                                ti.Units='normalized'; % Set units to 'normalized'
                                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                            end
                            if qq == NumFcMax*kk
                                break
                            end
                        end
                    else % kk==NumFigures % corresponds to the last figure plotted
                        t = tiledlayout(oo,pp);
                        t.TileSpacing = 'none'; % To reduce the spacing between the tiles
                        t.Padding = 'none'; % To reduce the padding of perimeter of a tile
                        NumFcPlot=oo*pp;
                        if kk==1
                            NumFcCorSelec(kk)=nnz(~obj.SMFSFlag.Uncorrupt(1:NumFcPlot));
                        elseif kk==2
                            NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt((NumFcPlot+1):(NumFcPlot*(kk))));
                        else
                            NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt((NumFcPlot*(kk-1)+1):(NumFcPlot*(kk))));
                        end
                        
                        if kk==1
                            qq=1;
                        else
                            qq=(kk-1)*mm*nn+1; % Determines the index of the next force curve to be plotted on the next figure
                        end
                        %% Plot loop
                        for qq=qq:obj.NCurves % Loop over all force curves in the force map
                            if ww<qq+DiffFc
                                ww=qq+DiffFc;
                            end
                            while ~obj.SMFSFlag.Selected(ww) || ~obj.SMFSFlag.Uncorrupt(ww) || ~obj.SMFSFlag.RetMinCrit(ww) || ~obj.SMFSFlag.LengthRequisite(ww)    % Stay in the while loop as long as the entry is zero
                                ww=ww+1;
                                if ww>qq
                                    DiffFc=ww-qq;
                                end
                            end
                            %% Allocate data
                        xApp=(obj.THApp{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                        xRet=(obj.THRet{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                        yApp=obj.BasedApp{qq+DiffFc}/Ymultiplier;
                        yRet=obj.BasedRet{qq+DiffFc}/Ymultiplier;    
                        yRetLim=obj.yRetLim{qq+DiffFc}/Ymultiplier;
                        yAppLim=obj.yAppLim{qq+DiffFc}/Ymultiplier;    
                            if ww>qq
                                ax=nexttile;
                                ax.XLim = [XMin XMax];
                                ax.YLim = [YMin YMax];
                                    if obj.SMFSFlag.PullingLength(qq+DiffFc)==1 && obj.SMFSFlag.SnapIn(qq+DiffFc)==1 % Force curve posses pulling length and snap-in variables
                                        grid on
                                        hold on
                                        area(xRet(1:obj.PullingLengthIdx(qq+DiffFc)),yRetLim(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                        area(xApp(obj.SnapInIdx(qq+DiffFc):end),yAppLim(obj.SnapInIdx(qq+DiffFc):end),'FaceColor',RGB7)
                                        plot(xApp,yApp,'Color',RGB1);
                                        plot(xRet,yRet,'Color',RGB2);
                                        plot(xRet(obj.PullingLengthIdx(qq+DiffFc)),yRet(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)
                                        plot(xApp(obj.SnapInIdx(qq+DiffFc)),yApp(obj.SnapInIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB10,'MarkerEdgeColor',RGB10)
                                        plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                        plot(xApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),yApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)
                                        plot(xRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),yRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)
                                    elseif obj.SMFSFlag.PullingLength(qq+DiffFc)==1 && obj.SMFSFlag.SnapIn(qq+DiffFc)==0 % Force curve posses pulling length variables but no snap-in variables
                                        grid on
                                        hold on
                                        area(xRet(1:obj.PullingLengthIdx(qq+DiffFc)),yRetLim(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                        plot(xApp,yApp,'Color',RGB1);
                                        plot(xRet,yRet,'Color',RGB2);
                                        plot(xRet(obj.PullingLengthIdx(qq+DiffFc)),yRet(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)
                                        plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                        plot(xRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),yRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)
                                    elseif obj.SMFSFlag.PullingLength(qq+DiffFc)==0 && obj.SMFSFlag.SnapIn(qq+DiffFc)==0 % Force curve posses neither pulling length variables nor snap-in variables
                                        grid on
                                        hold on
                                        plot(xApp,yApp,'Color',RGB1);
                                        plot(xRet,yRet,'Color',RGB2);
                                        plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                        plot(xApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),yApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)
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
                                    if obj.SMFSFlag.PullingLength(qq+DiffFc)==1 && obj.SMFSFlag.SnapIn(qq+DiffFc)==1 % Force curve posses pulling length and snap-in variables
                                        grid on
                                        hold on
                                        area(xRet(1:obj.PullingLengthIdx(qq+DiffFc)),yRetLim(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                        area(xApp(obj.SnapInIdx(qq+DiffFc):end),yAppLim(obj.SnapInIdx(qq+DiffFc):end),'FaceColor',RGB7)
                                        plot(xApp,yApp,'Color',RGB1);
                                        plot(xRet,yRet,'Color',RGB2);
                                        plot(xRet(obj.PullingLengthIdx(qq+DiffFc)),yRet(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)
                                        plot(xApp(obj.SnapInIdx(qq+DiffFc)),yApp(obj.SnapInIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB10,'MarkerEdgeColor',RGB10)
                                        plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                        plot(xApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),yApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)
                                        plot(xRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),yRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)
                                    elseif obj.SMFSFlag.PullingLength(qq+DiffFc)==1 && obj.SMFSFlag.SnapIn(qq+DiffFc)==0 % Force curve posses pulling length variables but no snap-in variables
                                        grid on
                                        hold on
                                        area(xRet(1:obj.PullingLengthIdx(qq+DiffFc)),yRetLim(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                        plot(xApp,yApp,'Color',RGB1);
                                        plot(xRet,yRet,'Color',RGB2);
                                        plot(xRet(obj.PullingLengthIdx(qq+DiffFc)),yRet(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)
                                        plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                        plot(xRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),yRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)
                                    elseif obj.SMFSFlag.PullingLength(qq+DiffFc)==0 && obj.SMFSFlag.SnapIn(qq+DiffFc)==0 % Force curve posses neither pulling length variables nor snap-in variables
                                        grid on
                                        hold on
                                        plot(xApp,yApp,'Color',RGB1);
                                        plot(xRet,yRet,'Color',RGB2);
                                        plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                        plot(xApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),yApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)
                                    end                               
                                % Title for each Subplot
                                ti=title(sprintf('%i',qq+DiffFc),'Color','k');
                                %ti=title(sprintf('%i',(kk+ww)/2),'Color','k');
                                ti.Units='normalized'; % Set units to 'normalized'
                                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                            end
                            if qq == NumFcUncorrupt
                                break
                            end
                        end
                    end
                end
                %% Save figures
                %%% Define the name for the figure title
                partname=sprintf('-p%d',kk);
                % fullname=sprintf('%s%s',figname,partname);
                fullname=sprintf('%s%s',figname,partname);
                %%% Save the current figure in the current folder
                print(gcf,fullname,'-dpng');
            end
            close all
        end
        
       
        function fc_visual_selection_flag_pre(obj,XMin,XMax,YMin,YMax,NumFcMax) % fc ... force curve

            %% Define variables
            RGB1=[0 26 255]./255;  % Blue
            RGB2=[255 119 0]./255; % Orange
            RGB7=[255 230 0]./255; % Yellow
            RGB8=[80 200 204]./255; % Turquoise
            RGB10=[200 0 255]./255; % Violet
            RGB11=[200 255 150]./255; % Light Green
            RGB12=[185 230 254]./255; % Light Blue
            RGB13=[200 0 0]./255; % Red
            x50=-50e-9; % Defines 50nm
            x150=-150e-9; % Defines 150nm
            x500=-500e-9; % Defines 500nm
            DataShareStartApp=0.15; % 15%
            DataShareEndApp=0.65; % 65%
            DataShareStartRet=0.07; % 7 %
            DataShareEndRet=0.02; % 2 %
            % Parse unit scale function
            [Xmultiplier,Xunit,~] = AFMImage.parse_unit_scale(1e+9,'nm',1);
            [Ymultiplier,Yunit,~] = AFMImage.parse_unit_scale(1e+9,'nN',1);
            % Define variables for the figure name
            ExtendVelocityConvert=num2str(round(obj.ExtendVelocity*1e9));
            RetractVelocityConvert=num2str(round(obj.RetractVelocity*1e9));
            HoldingTimeConvert=num2str(round(obj.HoldingTime));        
            % Classification criteria
            figname=strcat(obj.Date,{'_'},obj.Time,{'_'},obj.ID,{'_'},obj.Substrate,{'_'},obj.EnvCond,{'_'},obj.Linker,{'_'},obj.Chipbox,{'_'},obj.ChipCant,{'_'},ExtendVelocityConvert,{'_'},RetractVelocityConvert,{'_'},HoldingTimeConvert);
            figname=char(figname);
            %% Define variables for the plot loop
            mm=ceil(sqrt(NumFcMax)); % Variable defines the tiledlayout and thereby the tile arrangement
            nn=mm; % Variable defines the tiledlayout and thereby the tile arrangement
            ww=1; % Variable used to loop through the flag based while loop
            DiffFc=0; % Variable to correct for differences of the variables between the plot loop and flag based while loop
            NumFcUncorrupt=nnz(obj.SMFSFlag.Uncorrupt.*obj.SMFSFlag.RetMinCrit); % Determine the number of force curves that could been analysed
            NumFigures=ceil(NumFcUncorrupt/NumFcMax);
%             if NumFigures==0     % If condition is fulfilled stop function and return to calling function
%                 return
%             end
            RemainderMax=mod(NumFcUncorrupt,NumFcMax); % Check for remainder
            if RemainderMax ~= 0
                % Additional variables if there is a remainder
                oo=round(sqrt(RemainderMax)); % Determine the number of rows in the figure
                pp=ceil(sqrt(RemainderMax)); % Determine the number of columns in the figure
                RemainderReal=mod(NumFcUncorrupt,oo*pp); % Correct the remainder based on the determined rows times columns
            end
            %% figure loop 
            for kk=1:NumFigures
                % Figure
                h_fig=figure(kk);
                h_fig.Color='white'; % changes the background color of the figure
                h_fig.Units='normalized'; % Defines the units
                %h_fig.OuterPosition=Res;% changes the size to the defined
                % resolution
                h_fig.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
                h_fig.PaperOrientation='landscape';
                h_fig.Name=figname;
                %% Verify the remainder condition
                if RemainderMax == 0 % No remainder existing
                    t = tiledlayout(mm,nn);
                    t.TileSpacing = 'none'; % To reduce the spacing between the tiles
                    t.Padding = 'none'; % To reduce the padding of perimeter of a tile
                    if kk==1
                        NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt(1:NumFcMax));
                    elseif kk==2
                        NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt((NumFcMax+1):(NumFcMax*(kk))));
                    else
                        NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt((NumFcMax*(kk-1)+1):(NumFcMax*(kk))));
                    end
                    if kk==1
                        qq=1;
                    else
                        qq=(kk-1)*mm*nn+1; % Determines the index of the next force curve to be plotted on the next figure
                    end
                    %% Plot loop
                    for qq=qq:obj.NCurves % Loop over all force curves in the force map
                        if ww<qq+DiffFc
                            ww=qq+DiffFc;
                        end
                        % Flag based while loop
                        while ~obj.SMFSFlag.Uncorrupt(ww) || ~obj.SMFSFlag.RetMinCrit(ww)   % Stay in the while loop as long as the entry is zero
                            ww=ww+1;
                            if ww>qq
                                DiffFc=ww-qq;
                            end
                            %NumFcStartIdx(kk)
                        end
                        % Save the indeces of plotted force curves
                        IdxFcPlot(ww)=ww;
                        %% Allocate data
                        xApp=(obj.THApp{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                        xRet=(obj.THRet{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                        yApp=obj.BasedApp{qq+DiffFc}/Ymultiplier;
                        yRet=obj.BasedRet{qq+DiffFc}/Ymultiplier;
                        % Determine indices
                        IdxStartApp=round(length(yApp)*DataShareStartApp);
                        IdxEndApp=round(length(yApp)*DataShareEndApp);
                        IdxStartRet=round(length(yRet)-length(yRet)*DataShareStartRet);
                        IdxEndRet=round(length(yRet)-length(yRet)*DataShareEndRet);
                        % if condition with flag based while loop variable
                        % ww > plot loop variable qq
                        if ww>qq
                            ax=nexttile;
                            ax.XLim = [XMin XMax];
                            ax.YLim = [YMin YMax];
                            grid on
                            hold on
                            plot(xApp,yApp,'Color',RGB1);
                            plot(xRet,yRet,'Color',RGB2);
                            plot(xApp(IdxStartApp:IdxEndApp),yApp(IdxStartApp:IdxEndApp),'Color',RGB13);
                            plot(xRet(IdxStartRet:IdxEndRet),yRet(IdxStartRet:IdxEndRet),'Color','k');
                            xline(x50,'Color',RGB10); % Draws a vertical line
                            xline(x150,'Color',RGB10); % Draws a vertical line
                            xline(x500,'Color',RGB10); % Draws a vertical line
                            % Title for each Subplot
                            ti=title(sprintf('%i',qq+DiffFc),'Color','k');
                            ti.Units='normalized'; % Set units to 'normalized'
                            ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                        else
                            ax=nexttile;
                            ax.XLim = [XMin XMax];
                            ax.YLim = [YMin YMax];
                            grid on
                            hold on
                            plot(xApp,yApp,'Color',RGB1);
                            plot(xRet,yRet,'Color',RGB2);
                            plot(xApp(IdxStartApp:IdxEndApp),yApp(IdxStartApp:IdxEndApp),'Color',RGB13);
                            plot(xRet(IdxStartRet:IdxEndRet),yRet(IdxStartRet:IdxEndRet),'Color','k');
                            xline(x50,'Color',RGB10); % Draws a vertical line
                            xline(x150,'Color',RGB10); % Draws a vertical line
                            xline(x500,'Color',RGB10); % Draws a vertical line
                            % Title for each Subplot
                            ti=title(sprintf('%i',qq+DiffFc),'Color','k');
                            ti.Units='normalized'; % Set units to 'normalized'
                            ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                        end
                        if qq == NumFcMax*kk
                            break
                        end
                    end
                    
                else % Remainder exisiting
                    %% Last figure conditions
                    if kk~=NumFigures
                        t = tiledlayout(mm,nn);
                        t.TileSpacing = 'none'; % To reduce the spacing between the tiles
                        t.Padding = 'none'; % To reduce the padding of perimeter of a tile
                        if kk==1
                            NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt(1:NumFcMax));
                        elseif kk==2
                            NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt((NumFcMax+1):(NumFcMax*(kk))));
                        else
                            NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt((NumFcMax*(kk-1)+1):(NumFcMax*(kk))));
                        end
                        
                        if kk==1
                            qq=1;
                        else
                            qq=(kk-1)*mm*nn+1; % Determines the index of the next force curve to be plotted on the next figure
                        end
                        %% Plot loop
                        for qq=qq:obj.NCurves % Loop over all force curves in the force map
                            if ww<qq+DiffFc
                                ww=qq+DiffFc;
                            end
                            while ~obj.SMFSFlag.Uncorrupt(ww) || ~obj.SMFSFlag.RetMinCrit(ww)    % Stay in the while loop as long as the entry is zero
                                ww=ww+1;
                                if ww>qq
                                    DiffFc=ww-qq;
                                end
                            end
                            % Save the indeces of plotted force curves
                            IdxFcPlot(ww)=ww;
                            %% Allocate data
                        xApp=(obj.THApp{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                        xRet=(obj.THRet{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                        yApp=obj.BasedApp{qq+DiffFc}/Ymultiplier;
                        yRet=obj.BasedRet{qq+DiffFc}/Ymultiplier;
                            % Determine indices
                            IdxStartApp=round(length(yApp)*DataShareStartApp);
                            IdxEndApp=round(length(yApp)*DataShareEndApp);
                            IdxStartRet=round(length(yRet)-length(yRet)*DataShareStartRet);
                            IdxEndRet=round(length(yRet)-length(yRet)*DataShareEndRet);
                            if ww>qq
                                ax=nexttile;
                                ax.XLim = [XMin XMax];
                                ax.YLim = [YMin YMax];
                                grid on
                                hold on
                                plot(xApp,yApp,'Color',RGB1);
                                plot(xRet,yRet,'Color',RGB2);
                                plot(xApp(IdxStartApp:IdxEndApp),yApp(IdxStartApp:IdxEndApp),'Color',RGB13);
                                plot(xRet(IdxStartRet:IdxEndRet),yRet(IdxStartRet:IdxEndRet),'Color','k');
                                xline(x50,'Color',RGB10); % Draws a vertical line
                                xline(x150,'Color',RGB10); % Draws a vertical line
                                xline(x500,'Color',RGB10); % Draws a vertical line
                                % Title for each Subplot
                                ti=title(sprintf('%i',qq+DiffFc),'Color','k');
                                ti.Units='normalized'; % Set units to 'normalized'
                                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                            else
                                ax=nexttile;
                                ax.XLim = [XMin XMax];
                                ax.YLim = [YMin YMax];
                                grid on
                                hold on
                                plot(xApp,yApp,'Color',RGB1);
                                plot(xRet,yRet,'Color',RGB2);
                                plot(xApp(IdxStartApp:IdxEndApp),yApp(IdxStartApp:IdxEndApp),'Color',RGB13);
                                plot(xRet(IdxStartRet:IdxEndRet),yRet(IdxStartRet:IdxEndRet),'Color','k');
                                xline(x50,'Color',RGB10); % Draws a vertical line
                                xline(x150,'Color',RGB10); % Draws a vertical line
                                xline(x500,'Color',RGB10); % Draws a vertical line
                                % Title for each Subplot
                                ti=title(sprintf('%i',qq+DiffFc),'Color','k');
                                ti.Units='normalized'; % Set units to 'normalized'
                                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                            end
                            if qq == NumFcMax*kk
                                break
                            end
                        end
                    else % kk==NumFigures % corresponds to the last figure plotted
                        t = tiledlayout(oo,pp);
                        t.TileSpacing = 'none'; % To reduce the spacing between the tiles
                        t.Padding = 'none'; % To reduce the padding of perimeter of a tile
                        NumFcPlot=oo*pp;
                        if kk==1
                            NumFcCorSelec(kk)=nnz(~obj.SMFSFlag.Uncorrupt(1:NumFcPlot));
                        elseif kk==2
                            NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt((NumFcPlot+1):(NumFcPlot*(kk))));
                        else
                            NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt((NumFcPlot*(kk-1)+1):(NumFcPlot*(kk))));
                        end
                        
                        if kk==1
                            qq=1;
                        else
                            qq=(kk-1)*mm*nn+1; % Determines the index of the next force curve to be plotted on the next figure
                        end
                        %% Plot loop
                        for qq=qq:obj.NCurves % Loop over all force curves in the force map
                            if ww<qq+DiffFc
                                ww=qq+DiffFc;
                            end
                            while ~obj.SMFSFlag.Uncorrupt(ww) || ~obj.SMFSFlag.RetMinCrit(ww)   % Stay in the while loop as long as the entry is zero
                                ww=ww+1;
                                if ww>qq
                                    DiffFc=ww-qq;
                                end
                            end
                            % Save the indeces of plotted force curves
                            IdxFcPlot(ww)=ww;
                            %% Allocate data
                        xApp=(obj.THApp{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                        xRet=(obj.THRet{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                        yApp=obj.BasedApp{qq+DiffFc}/Ymultiplier;
                        yRet=obj.BasedRet{qq+DiffFc}/Ymultiplier;
                            % Determine indices
                            IdxStartApp=round(length(yApp)*DataShareStartApp);
                            IdxEndApp=round(length(yApp)*DataShareEndApp);
                            IdxStartRet=round(length(yRet)-length(yRet)*DataShareStartRet);
                            IdxEndRet=round(length(yRet)-length(yRet)*DataShareEndRet);
                            if ww>qq
                                ax=nexttile;
                                ax.XLim = [XMin XMax];
                                ax.YLim = [YMin YMax];
                                grid on
                                hold on
                                plot(xApp,yApp,'Color',RGB1);
                                plot(xRet,yRet,'Color',RGB2);
                                plot(xApp(IdxStartApp:IdxEndApp),yApp(IdxStartApp:IdxEndApp),'Color',RGB13);
                                plot(xRet(IdxStartRet:IdxEndRet),yRet(IdxStartRet:IdxEndRet),'Color','k');
                                xline(x50,'Color',RGB10); % Draws a vertical line
                                xline(x150,'Color',RGB10); % Draws a vertical line
                                xline(x500,'Color',RGB10); % Draws a vertical line
                                % Title for each Subplot
                                ti=title(sprintf('%i',qq+DiffFc),'Color','k');
                                %ti=title(sprintf('%i',(kk+ww)/2),'Color','k');
                                ti.Units='normalized'; % Set units to 'normalized'
                                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                            else
                                ax=nexttile;
                                ax.XLim = [XMin XMax];
                                ax.YLim = [YMin YMax];
                                grid on
                                hold on
                                plot(xApp,yApp,'Color',RGB1);
                                plot(xRet,yRet,'Color',RGB2);
                                plot(xApp(IdxStartApp:IdxEndApp),yApp(IdxStartApp:IdxEndApp),'Color',RGB13);
                                plot(xRet(IdxStartRet:IdxEndRet),yRet(IdxStartRet:IdxEndRet),'Color','k');
                                xline(x50,'Color',RGB10); % Draws a vertical line
                                xline(x150,'Color',RGB10); % Draws a vertical line
                                xline(x500,'Color',RGB10); % Draws a vertical line
                                % Title for each Subplot
                                ti=title(sprintf('%i',qq+DiffFc),'Color','k');
                                %ti=title(sprintf('%i',(kk+ww)/2),'Color','k');
                                ti.Units='normalized'; % Set units to 'normalized'
                                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                            end
                            if qq == NumFcUncorrupt
                                break
                            end
                        end
                    end
                end
                % Remove zeros from the plotted force curves index array
                IdxFcPlot=nonzeros(IdxFcPlot);
                FcPerFig(1)=1;
                FcPerFig(kk+1)=length(IdxFcPlot);
                IdxFcCorr=IdxFcPlot; % Allocate array which will correct the plotted fc array by the selected fc
                if kk>1
                IdxFcCorr(FcPerFig(1):FcPerFig(kk)) = []; % Remove force curves of previous figure
                else
                    
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
                    obj.SMFSFlag.Uncorrupt(IdxFcCorr)=1;
                    % 'Blue' loop - Repaint tile titles for force curves
                    % kept for the analysis
                    for ii=1:length(IdxFcCorr)
                        nexttile(ii) % Access the tile corresponding to the selected force curve
                        ti=title(sprintf('%i',IdxFcCorr(ii)),'Color','b'); % Change title color
                        ti.Units='normalized'; % Set units to 'normalized'
                        ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                    end                    
                end
                % Case 2: Select none
                 if SelectBttns == 2
                    obj.SMFSFlag.Uncorrupt(IdxFcCorr)=0;
                    % 'Red' loop - Repaint tile titles for force curves
                    % removed from the analysis
                    for jj=1:length(IdxFcCorr)
                        nexttile(jj) % Access the tile corresponding to the selected force curve
                        ti=title(sprintf('%i',IdxFcCorr(jj)),'Color','r'); % Change title color
                        ti.Units='normalized'; % Set units to 'normalized'
                        ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                    end                    
                 end                   
                % Case 3: Select all - except of
                if SelectBttns == 3                  
                    obj.SMFSFlag.Uncorrupt(IdxFcCorr)=1;
                    prompt = {'Enter the force curve number you do not want to keep for analysis (For multiple selections just use the space key to separate entries)'};
                    definput = {''};
                    opts.Interpreter = 'tex';
                    IdxFcSel=inputdlg(prompt,'Select all - except of ...',[1 150],definput,opts); % Stores the individual selected fc as a cell array of character vectors
                    IdxFcSel=str2num(IdxFcSel{1}); % Convert the cell array to numerals
                    obj.SMFSFlag.Uncorrupt(IdxFcSel)=0; % Set flag for the selected force curves                  
                    % 'Blue' loop - Repaint tile titles for force curves
                    % kept for the analysis
                    for ii=1:length(IdxFcCorr)
                        nexttile(ii); % Access the tile corresponding to the selected force curve
                        ti=title(sprintf('%i',IdxFcCorr(ii)),'Color','b'); % Change title color
                        ti.Units='normalized'; % Set units to 'normalized'
                        ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                    end                    
                    % 'Red' loop - Repaint tile titles for force curves
                    % removed from the analysis
                    for jj=1:length(IdxFcSel)
                        TileIdx(jj)=find(IdxFcCorr==IdxFcSel(jj));
                        nexttile(TileIdx(jj)) % Access the tile corresponding to the selected force curve
                        ti=title(sprintf('%i',IdxFcSel(jj)),'Color','r'); % Change title color
                        ti.Units='normalized'; % Set units to 'normalized'
                        ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                    end
                end
                
                % Case 4: Select none - except of
                if SelectBttns == 4
                    obj.SMFSFlag.Uncorrupt(IdxFcCorr)=0;
                    prompt = {'Enter the force curve number you want want to keep for analysis (For multiple selections just use the space key to separate entries)'};
                    definput = {''};
                    opts.Interpreter = 'tex';
                    IdxFcSel=inputdlg(prompt,'Select none - except of ...',[1 150],definput,opts); % Stores the individual selected fc as a cell array of character vectors
                    IdxFcSel=str2num(IdxFcSel{1}); % Convert the cell array to numerals
                    obj.SMFSFlag.Uncorrupt(IdxFcSel)=1; % Set flag for the selected force curves                    
                    % 'Red' loop - Repaint tile titles for force curves
                    % removed from the analysis
                    for ii=1:length(IdxFcCorr)
                        nexttile(ii) % Access the tile corresponding to the selected force curve
                        ti=title(sprintf('%i',IdxFcCorr(ii)),'Color','r'); % Change title color
                        ti.Units='normalized'; % Set units to 'normalized'
                        ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                    end                    
                    % 'Blue' loop - Repaint tile titles for force curves
                    % kept for the analysis                    
                    for jj=1:length(IdxFcSel)
                        TileIdx(jj)=find(IdxFcCorr==IdxFcSel(jj));
                        nexttile(TileIdx(jj)) % Access the tile corresponding to the selected force curve
                        ti=title(sprintf('%i',IdxFcSel(jj)),'Color','b'); % Change title color
                        ti.Units='normalized'; % Set units to 'normalized'
                        ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                    end                    
                end
                %% Save figures
                %%% Define the name for the figure title
                partname=sprintf('-p%d',kk);
                fullname=sprintf('%s%s',figname,partname);
                %%% Save the current figure in the current folder
                print(gcf,fullname,'-dpng');  
            end                         
            close all        
        end
        
        function fc_visual_selection_post(obj,XMin,XMax,YMin,YMax,NumFcMax) % fc ... force curve
            % fc_print_adhenergy_pulllength: A function to plot all selected force curves of a
            % force map including adhesion energy and pulling length in
            % each force curve
            if nargin < 2
                XMin= -inf;
                XMax= inf;
                YMin= -inf;
                YMax= inf;
                NumFcMax=25;
            end
            % Define variables
            RGB1=[0 26 255]./255;  % Blue
            RGB2=[255 119 0]./255; % Orange
            RGB7=[255 230 0]./255; % Yellow
            RGB8=[80 200 204]./255; % Turquoise
            RGB10=[200 0 255]./255; % Violet
            RGB11=[200 255 150]./255; % Light Green
            RGB13=[200 0 0]./255; % Red     
            Res=[1 1 2560 1440]; % Define the figure resolution
            % NumFcMax=25;
            % Parse unit scale function
            [Xmultiplier,Xunit,~] = AFMImage.parse_unit_scale(1e+9,'nm',1);
            [Ymultiplier,Yunit,~] = AFMImage.parse_unit_scale(1e+9,'nN',1);
            % Define variables for the figure name
            ExtendVelocityConvert=num2str(round(obj.ExtendVelocity*1e9));
            RetractVelocityConvert=num2str(round(obj.RetractVelocity*1e9));
            HoldingTimeConvert=num2str(round(obj.HoldingTime));        
            % Classification criteria
            figname=strcat(obj.Date,{'_'},obj.Time,{'_'},obj.ID,{'_'},obj.Substrate,{'_'},obj.EnvCond,{'_'},obj.Linker,{'_'},obj.Chipbox,{'_'},obj.ChipCant,{'_'},ExtendVelocityConvert,{'_'},RetractVelocityConvert,{'_'},HoldingTimeConvert);
            figname=char(figname);
            %% Define variables for the plot loop
            mm=ceil(sqrt(NumFcMax)); % Variable defines the tiledlayout and thereby the tile arrangement
            nn=mm; % Variable defines the tiledlayout and thereby the tile arrangement
            ww=1; % Variable used to loop through the flag based while loop
            DiffFc=0; % Variable to correct for differences of the variables between the plot loop and flag based while loop       
            NumFcUncorrupt=nnz(obj.SMFSFlag.Uncorrupt.*obj.SMFSFlag.RetMinCrit.*obj.SMFSFlag.LengthRequisite); % Determine the number of force curves that could been analysed
            NumFigures=ceil(NumFcUncorrupt/NumFcMax);
            if NumFigures==0     % If condition is fulfilled stop function and return to calling function
                return
            end
            RemainderMax=mod(NumFcUncorrupt,NumFcMax); % Check for remainder
            if RemainderMax ~= 0
                % Additional variables if there is a remainder
                oo=round(sqrt(RemainderMax)); % Determine the number of rows in the figure
                pp=ceil(sqrt(RemainderMax)); % Determine the number of columns in the figure
                RemainderReal=mod(NumFcUncorrupt,oo*pp); % Correct the remainder based on the determined rows times columns
            end
            %% figure loop
            for kk=1:NumFigures
               % if  ~obj.DebugFlag.Plot % Suppress plotting
                if  obj.DebugFlag.Plot % Allow plotting
                    continue
                end
                % Figure
                h_fig=figure(kk);
                h_fig.Color='white'; % changes the background color of the figure
                h_fig.Units='pixel'; % Defines the units
                h_fig.OuterPosition=Res;
                h_fig.PaperOrientation='landscape';
                h_fig.Name=figname;
                %% Verify the remainder condition
                if RemainderMax == 0 % No remainder existing
                    t = tiledlayout(mm,nn);
                    t.TileSpacing = 'none'; % To reduce the spacing between the tiles
                    t.Padding = 'none'; % To reduce the padding of perimeter of a tile
                    if kk==1
                        NumFcCorSelec(kk)=nnz(~obj.SMFSFlag.Uncorrupt(1:NumFcMax));
                    elseif kk==2
                        NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt((NumFcMax+1):(NumFcMax*(kk))));
                    else
                        NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt((NumFcMax*(kk-1)+1):(NumFcMax*(kk))));
                    end
                    
                    if kk==1
                        qq=1;
                    else
                        qq=(kk-1)*mm*nn+1; % Determines the index of the next force curve to be plotted on the next figure
                    end
                    %% Plot loop
                    for qq=qq:obj.NCurves % Loop over all force curves in the force map
                        % Debugging
                    %for qq=5
                        if ww<qq+DiffFc
                            ww=qq+DiffFc;
                        end
                        % Flag based while loop
                        while ~obj.SMFSFlag.Uncorrupt(ww) || ~obj.SMFSFlag.RetMinCrit(ww) || ~obj.SMFSFlag.LengthRequisite(ww)   % Stay in the while loop as long as the entry is zero
                            ww=ww+1;
                            if ww>qq
                                DiffFc=ww-qq;
                            end
                        end
                        % Save the indeces of plotted force curves
                        IdxFcPlot(ww)=ww;
                        %% Allocate data
                        xApp=(obj.THApp{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                        xRet=(obj.THRet{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                        yApp=obj.BasedApp{qq+DiffFc}/Ymultiplier;
                        yRet=obj.BasedRet{qq+DiffFc}/Ymultiplier;    
                        yRetLim=obj.yRetLim{qq+DiffFc}/Ymultiplier;
                        yAppLim=obj.yAppLim{qq+DiffFc}/Ymultiplier;
                        % if condition with flag based while loop variable
                        % ww > plot loop variable qq
                        if ww>qq
                            ax=nexttile;
                            ax.XLim = [XMin XMax];
                            ax.YLim = [YMin YMax];
                                if obj.SMFSFlag.PullingLength(qq+DiffFc)==1 && obj.SMFSFlag.SnapIn(qq+DiffFc)==1 % Force curve posses pulling length and snap-in variables
                                    grid on
                                    hold on
                                    area(xRet(1:obj.PullingLengthIdx(qq+DiffFc)),yRetLim(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                    area(xApp(obj.SnapInIdx(qq+DiffFc):end),yAppLim(obj.SnapInIdx(qq+DiffFc):end),'FaceColor',RGB7)
                                    plot(xApp,yApp,'Color',RGB1);
                                    plot(xRet,yRet,'Color',RGB2);
                                    plot(xRet(obj.PullingLengthIdx(qq+DiffFc)),yRet(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)
                                    plot(xApp(obj.SnapInIdx(qq+DiffFc)),yApp(obj.SnapInIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB10,'MarkerEdgeColor',RGB10)
                                    plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                    plot(xApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),yApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)
                                    plot(xRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),yRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)
                                elseif obj.SMFSFlag.PullingLength(qq+DiffFc)==1 && obj.SMFSFlag.SnapIn(qq+DiffFc)==0 % Force curve posses pulling length variables but no snap-in variables
                                    grid on
                                    hold on
                                    area(xRet(1:obj.PullingLengthIdx(qq+DiffFc)),yRetLim(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                    plot(xApp,yApp,'Color',RGB1);
                                    plot(xRet,yRet,'Color',RGB2);
                                    plot(xRet(obj.PullingLengthIdx(qq+DiffFc)),yRet(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)
                                    plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                    plot(xRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),yRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)
                                elseif obj.SMFSFlag.PullingLength(qq+DiffFc)==0 && obj.SMFSFlag.SnapIn(qq+DiffFc)==0 % Force curve posses neither pulling length variables nor snap-in variables
                                    grid on
                                    hold on
                                    plot(xApp,yApp,'Color',RGB1);
                                    plot(xRet,yRet,'Color',RGB2);
                                    plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                    plot(xApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),yApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)
                                end                 
                            % Title for each Subplot
                            ti=title(sprintf('%i',qq+DiffFc),'Color','k');
                            ti.Units='normalized'; % Set units to 'normalized'
                            ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                        else
                            ax=nexttile;
                            ax.XLim = [XMin XMax];
                            ax.YLim = [YMin YMax];
                                if obj.SMFSFlag.PullingLength(qq+DiffFc)==1 && obj.SMFSFlag.SnapIn(qq+DiffFc)==1 % Force curve posses pulling length and snap-in variables
                                    grid on
                                    hold on
                                    area(xRet(1:obj.PullingLengthIdx(qq+DiffFc)),yRetLim(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                    area(xApp(obj.SnapInIdx(qq+DiffFc):end),yAppLim(obj.SnapInIdx(qq+DiffFc):end),'FaceColor',RGB7)
                                    plot(xApp,yApp,'Color',RGB1);
                                    plot(xRet,yRet,'Color',RGB2);
                                    plot(xRet(obj.PullingLengthIdx(qq+DiffFc)),yRet(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)
                                    plot(xApp(obj.SnapInIdx(qq+DiffFc)),yApp(obj.SnapInIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB10,'MarkerEdgeColor',RGB10)
                                    plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                    plot(xApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),yApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)
                                    plot(xRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),yRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)
                                elseif obj.SMFSFlag.PullingLength(qq+DiffFc)==1 && obj.SMFSFlag.SnapIn(qq+DiffFc)==0 % Force curve posses pulling length variables but no snap-in variables
                                    grid on
                                    hold on
                                    area(xRet(1:obj.PullingLengthIdx(qq+DiffFc)),yRetLim(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                    plot(xApp,yApp,'Color',RGB1);
                                    plot(xRet,yRet,'Color',RGB2);
                                    plot(xRet(obj.PullingLengthIdx(qq+DiffFc)),yRet(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)
                                    plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                    plot(xRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),yRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)
                                elseif obj.SMFSFlag.PullingLength(qq+DiffFc)==0 && obj.SMFSFlag.SnapIn(qq+DiffFc)==0 % Force curve posses neither pulling length variables nor snap-in variables
                                    grid on
                                    hold on
                                    plot(xApp,yApp,'Color',RGB1);
                                    plot(xRet,yRet,'Color',RGB2);
                                    plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                    plot(xApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),yApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)
                                end
                            % Title for each Subplot
                            ti=title(sprintf('%i',qq+DiffFc),'Color','k');
                            ti.Units='normalized'; % Set units to 'normalized'
                            ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                        end
                        if qq == NumFcMax*kk
                            break
                        end
                    end
                    
                else % Remainder exisiting
                    %% Last figure conditions
                    if kk~=NumFigures
                        t = tiledlayout(mm,nn);
                        t.TileSpacing = 'none'; % To reduce the spacing between the tiles
                        t.Padding = 'none'; % To reduce the padding of perimeter of a tile
                        if kk==1
                            NumFcCorSelec(kk)=nnz(~obj.SMFSFlag.Uncorrupt(1:NumFcMax));
                        elseif kk==2
                            NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt((NumFcMax+1):(NumFcMax*(kk))));
                        else
                            NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt((NumFcMax*(kk-1)+1):(NumFcMax*(kk))));
                        end
                        
                        if kk==1
                            qq=1;
                        else
                            qq=(kk-1)*mm*nn+1; % Determines the index of the next force curve to be plotted on the next figure
                        end
                        %% Plot loop
                        for qq=qq:obj.NCurves % Loop over all force curves in the force map
                            if ww<qq+DiffFc
                                ww=qq+DiffFc;
                            end
                            while ~obj.SMFSFlag.Uncorrupt(ww) || ~obj.SMFSFlag.RetMinCrit(ww) || ~obj.SMFSFlag.LengthRequisite(ww)    % Stay in the while loop as long as the entry is zero
                                ww=ww+1;
                                if ww>qq
                                    DiffFc=ww-qq;
                                end
                            end
                            % Save the indeces of plotted force curves
                            IdxFcPlot(ww)=ww;
                            %% Allocate data
                        xApp=(obj.THApp{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                        xRet=(obj.THRet{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                        yApp=obj.BasedApp{qq+DiffFc}/Ymultiplier;
                        yRet=obj.BasedRet{qq+DiffFc}/Ymultiplier;   
                        yRetLim=obj.yRetLim{qq+DiffFc}/Ymultiplier;
                        yAppLim=obj.yAppLim{qq+DiffFc}/Ymultiplier;
                            if ww>qq
                                ax=nexttile;
                                ax.XLim = [XMin XMax];
                                ax.YLim = [YMin YMax];
                                    if obj.SMFSFlag.PullingLength(qq+DiffFc)==1 && obj.SMFSFlag.SnapIn(qq+DiffFc)==1 % Force curve posses pulling length and snap-in variables
                                        grid on
                                        hold on
                                        area(xRet(1:obj.PullingLengthIdx(qq+DiffFc)),yRetLim(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                        area(xApp(obj.SnapInIdx(qq+DiffFc):end),yAppLim(obj.SnapInIdx(qq+DiffFc):end),'FaceColor',RGB7)
                                        plot(xApp,yApp,'Color',RGB1);
                                        plot(xRet,yRet,'Color',RGB2);
                                        plot(xRet(obj.PullingLengthIdx(qq+DiffFc)),yRet(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)
                                        plot(xApp(obj.SnapInIdx(qq+DiffFc)),yApp(obj.SnapInIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB10,'MarkerEdgeColor',RGB10)
                                        plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                        plot(xApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),yApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)
                                        plot(xRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),yRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)
                                    elseif obj.SMFSFlag.PullingLength(qq+DiffFc)==1 && obj.SMFSFlag.SnapIn(qq+DiffFc)==0 % Force curve posses pulling length variables but no snap-in variables
                                        grid on
                                        hold on
                                        area(xRet(1:obj.PullingLengthIdx(qq+DiffFc)),yRetLim(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                        plot(xApp,yApp,'Color',RGB1);
                                        plot(xRet,yRet,'Color',RGB2);
                                        plot(xRet(obj.PullingLengthIdx(qq+DiffFc)),yRet(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)
                                        plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                        plot(xRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),yRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)
                                    elseif obj.SMFSFlag.PullingLength(qq+DiffFc)==0 && obj.SMFSFlag.SnapIn(qq+DiffFc)==0 % Force curve posses neither pulling length variables nor snap-in variables
                                        grid on
                                        hold on
                                        plot(xApp,yApp,'Color',RGB1);
                                        plot(xRet,yRet,'Color',RGB2);
                                        plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                        plot(xApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),yApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)
                                    end
                                % Title for each Subplot
                                ti=title(sprintf('%i',qq+DiffFc),'Color','k');
                                ti.Units='normalized'; % Set units to 'normalized'
                                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                            else
                                ax=nexttile;
                                ax.XLim = [XMin XMax];
                                ax.YLim = [YMin YMax];
                                    if obj.SMFSFlag.PullingLength(qq+DiffFc)==1 && obj.SMFSFlag.SnapIn(qq+DiffFc)==1 % Force curve posses pulling length and snap-in variables
                                        grid on
                                        hold on
                                        area(xRet(1:obj.PullingLengthIdx(qq+DiffFc)),yRetLim(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                        area(xApp(obj.SnapInIdx(qq+DiffFc):end),yAppLim(obj.SnapInIdx(qq+DiffFc):end),'FaceColor',RGB7)
                                        plot(xApp,yApp,'Color',RGB1);
                                        plot(xRet,yRet,'Color',RGB2);
                                        plot(xRet(obj.PullingLengthIdx(qq+DiffFc)),yRet(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)
                                        plot(xApp(obj.SnapInIdx(qq+DiffFc)),yApp(obj.SnapInIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB10,'MarkerEdgeColor',RGB10)
                                        plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                        plot(xApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),yApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)
                                        plot(xRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),yRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)
                                    elseif obj.SMFSFlag.PullingLength(qq+DiffFc)==1 && obj.SMFSFlag.SnapIn(qq+DiffFc)==0 % Force curve posses pulling length variables but no snap-in variables
                                        grid on
                                        hold on
                                        area(xRet(1:obj.PullingLengthIdx(qq+DiffFc)),yRetLim(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                        plot(xApp,yApp,'Color',RGB1);
                                        plot(xRet,yRet,'Color',RGB2);
                                        plot(xRet(obj.PullingLengthIdx(qq+DiffFc)),yRet(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)
                                        plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                        plot(xRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),yRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)
                                    elseif obj.SMFSFlag.PullingLength(qq+DiffFc)==0 && obj.SMFSFlag.SnapIn(qq+DiffFc)==0 % Force curve posses neither pulling length variables nor snap-in variables
                                        grid on
                                        hold on
                                        plot(xApp,yApp,'Color',RGB1);
                                        plot(xRet,yRet,'Color',RGB2);
                                        plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                        plot(xApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),yApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)
                                    end
                                % Title for each Subplot
                                ti=title(sprintf('%i',qq+DiffFc),'Color','k');
                                ti.Units='normalized'; % Set units to 'normalized'
                                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                            end
                            if qq == NumFcMax*kk
                                break
                            end
                        end
                    else % kk==NumFigures % corresponds to the last figure plotted
                        t = tiledlayout(oo,pp);
                        t.TileSpacing = 'none'; % To reduce the spacing between the tiles
                        t.Padding = 'none'; % To reduce the padding of perimeter of a tile
                        NumFcPlot=oo*pp;
                        if kk==1
                            NumFcCorSelec(kk)=nnz(~obj.SMFSFlag.Uncorrupt(1:NumFcPlot));
                        elseif kk==2
                            NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt((NumFcPlot+1):(NumFcPlot*(kk))));
                        else
                            NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt((NumFcPlot*(kk-1)+1):(NumFcPlot*(kk))));
                        end
                        
                        if kk==1
                            qq=1;
                        else
                            qq=(kk-1)*mm*nn+1; % Determines the index of the next force curve to be plotted on the next figure
                        end
                        %% Plot loop
                        for qq=qq:obj.NCurves % Loop over all force curves in the force map
                            if ww<qq+DiffFc
                                ww=qq+DiffFc;
                            end
                            while ~obj.SMFSFlag.Uncorrupt(ww) || ~obj.SMFSFlag.RetMinCrit(ww) || ~obj.SMFSFlag.LengthRequisite(ww)    % Stay in the while loop as long as the entry is zero
                                ww=ww+1;
                                if ww>qq
                                    DiffFc=ww-qq;
                                end
                            end
                            % Save the indeces of plotted force curves
                            IdxFcPlot(ww)=ww;
                            %% Allocate data
                        xApp=(obj.THApp{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                        xRet=(obj.THRet{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                        yApp=obj.BasedApp{qq+DiffFc}/Ymultiplier;
                        yRet=obj.BasedRet{qq+DiffFc}/Ymultiplier; 
                        yRetLim=obj.yRetLim{qq+DiffFc}/Ymultiplier;
                        yAppLim=obj.yAppLim{qq+DiffFc}/Ymultiplier;
                            if ww>qq
                                ax=nexttile;
                                ax.XLim = [XMin XMax];
                                ax.YLim = [YMin YMax];
                                    if obj.SMFSFlag.PullingLength(qq+DiffFc)==1 && obj.SMFSFlag.SnapIn(qq+DiffFc)==1 % Force curve posses pulling length and snap-in variables
                                        grid on
                                        hold on
                                        area(xRet(1:obj.PullingLengthIdx(qq+DiffFc)),yRetLim(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                        area(xApp(obj.SnapInIdx(qq+DiffFc):end),yAppLim(obj.SnapInIdx(qq+DiffFc):end),'FaceColor',RGB7)
                                        plot(xApp,yApp,'Color',RGB1);
                                        plot(xRet,yRet,'Color',RGB2);
                                        plot(xRet(obj.PullingLengthIdx(qq+DiffFc)),yRet(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)
                                        plot(xApp(obj.SnapInIdx(qq+DiffFc)),yApp(obj.SnapInIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB10,'MarkerEdgeColor',RGB10)
                                        plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                        plot(xApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),yApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)
                                        plot(xRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),yRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)
                                    elseif obj.SMFSFlag.PullingLength(qq+DiffFc)==1 && obj.SMFSFlag.SnapIn(qq+DiffFc)==0 % Force curve posses pulling length variables but no snap-in variables
                                        grid on
                                        hold on
                                        area(xRet(1:obj.PullingLengthIdx(qq+DiffFc)),yRetLim(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                        plot(xApp,yApp,'Color',RGB1);
                                        plot(xRet,yRet,'Color',RGB2);
                                        plot(xRet(obj.PullingLengthIdx(qq+DiffFc)),yRet(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)
                                        plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                        plot(xRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),yRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)
                                    elseif obj.SMFSFlag.PullingLength(qq+DiffFc)==0 && obj.SMFSFlag.SnapIn(qq+DiffFc)==0 % Force curve posses neither pulling length variables nor snap-in variables
                                        grid on
                                        hold on
                                        plot(xApp,yApp,'Color',RGB1);
                                        plot(xRet,yRet,'Color',RGB2);
                                        plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                        plot(xApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),yApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)
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
                                    if obj.SMFSFlag.PullingLength(qq+DiffFc)==1 && obj.SMFSFlag.SnapIn(qq+DiffFc)==1 % Force curve posses pulling length and snap-in variables
                                        grid on
                                        hold on
                                        area(xRet(1:obj.PullingLengthIdx(qq+DiffFc)),yRetLim(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                        area(xApp(obj.SnapInIdx(qq+DiffFc):end),yAppLim(obj.SnapInIdx(qq+DiffFc):end),'FaceColor',RGB7)
                                        plot(xApp,yApp,'Color',RGB1);
                                        plot(xRet,yRet,'Color',RGB2);
                                        plot(xRet(obj.PullingLengthIdx(qq+DiffFc)),yRet(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)
                                        plot(xApp(obj.SnapInIdx(qq+DiffFc)),yApp(obj.SnapInIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB10,'MarkerEdgeColor',RGB10)
                                        plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                        plot(xApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),yApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)
                                        plot(xRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),yRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)
                                    elseif obj.SMFSFlag.PullingLength(qq+DiffFc)==1 && obj.SMFSFlag.SnapIn(qq+DiffFc)==0 % Force curve posses pulling length variables but no snap-in variables
                                        grid on
                                        hold on
                                        area(xRet(1:obj.PullingLengthIdx(qq+DiffFc)),yRetLim(1:obj.PullingLengthIdx(qq+DiffFc)),'FaceColor',RGB11)
                                        plot(xApp,yApp,'Color',RGB1);
                                        plot(xRet,yRet,'Color',RGB2);
                                        plot(xRet(obj.PullingLengthIdx(qq+DiffFc)),yRet(obj.PullingLengthIdx(qq+DiffFc)),'d','MarkerSize',14,'MarkerFaceColor',RGB8,'MarkerEdgeColor',RGB8)
                                        plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                        plot(xRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),yRet(obj.AdhForceUnbindingIdx(qq+DiffFc)),'p','MarkerSize',12,'MarkerFaceColor',RGB13,'MarkerEdgeColor',RGB13)
                                    elseif obj.SMFSFlag.PullingLength(qq+DiffFc)==0 && obj.SMFSFlag.SnapIn(qq+DiffFc)==0 % Force curve posses neither pulling length variables nor snap-in variables
                                        grid on
                                        hold on
                                        plot(xApp,yApp,'Color',RGB1);
                                        plot(xRet,yRet,'Color',RGB2);
                                        plot(xRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),yRet(obj.AdhForceMaxRetIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB2,'MarkerEdgeColor',RGB2)
                                        plot(xApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),yApp(obj.AdhForceMaxAppIdx(qq+DiffFc)),'h','MarkerSize',10,'MarkerFaceColor',RGB1,'MarkerEdgeColor',RGB1)
                                    end                               
                                % Title for each Subplot
                                ti=title(sprintf('%i',qq+DiffFc),'Color','k');
                                %ti=title(sprintf('%i',(kk+ww)/2),'Color','k');
                                ti.Units='normalized'; % Set units to 'normalized'
                                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                            end
                            if qq == NumFcUncorrupt
                                break
                            end
                        end               
                    end             
                end
                                % Remove zeros from the plotted force curves index array
                IdxFcPlot=nonzeros(IdxFcPlot);
                FcPerFig(1)=1;
                FcPerFig(kk+1)=length(IdxFcPlot);
                IdxFcCorr=IdxFcPlot; % Allocate array which will correct the plotted fc array by the selected fc
                if kk>1
                IdxFcCorr(FcPerFig(1):FcPerFig(kk)) = []; % Remove force curves of previous figure
                else
                    
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
                    obj.SMFSFlag.Selected(IdxFcCorr)=1;
                    % 'Blue' loop - Repaint tile titles for force curves
                    % kept for the analysis
                    for ii=1:length(IdxFcCorr)
                        nexttile(ii) % Access the tile corresponding to the selected force curve
                        ti=title(sprintf('%i',IdxFcCorr(ii)),'Color','b'); % Change title color
                        ti.Units='normalized'; % Set units to 'normalized'
                        ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                    end                    
                end
                % Case 2: Select none
                 if SelectBttns == 2
                    obj.SMFSFlag.Selected(IdxFcCorr)=0;
                    %'Red' loop - Repaint tile titles for force curves
                    % removed from the analysis
                    for jj=1:length(IdxFcCorr)
                        nexttile(jj) % Access the tile corresponding to the selected force curve
                        ti=title(sprintf('%i',IdxFcCorr(jj)),'Color','r'); % Change title color
                        ti.Units='normalized'; % Set units to 'normalized'
                        ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                    end                    
                 end                   
                % Case 3: Select all - except of
                if SelectBttns == 3                  
                    obj.SMFSFlag.Selected(IdxFcCorr)=1;
                    prompt = {'Enter the force curve number you do not want to keep for analysis (For multiple selections just use the space key to separate entries)'};
                    definput = {''};
                    opts.Interpreter = 'tex';
                    IdxFcSel=inputdlg(prompt,'Select all - except of ...',[1 150],definput,opts); % Stores the individual selected fc as a cell array of character vectors
                    IdxFcSel=str2num(IdxFcSel{1}); % Convert the cell array to numerals
                    obj.SMFSFlag.Selected(IdxFcSel)=0; % Set flag for the selected force curves                  
                    % 'Blue' loop - Repaint tile titles for force curves
                    % kept for the analysis
                    for ii=1:length(IdxFcCorr)
                        nexttile(ii); % Access the tile corresponding to the selected force curve
                        ti=title(sprintf('%i',IdxFcCorr(ii)),'Color','b'); % Change title color
                        ti.Units='normalized'; % Set units to 'normalized'
                        ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                    end                    
                    % 'Red' loop - Repaint tile titles for force curves
                    % removed from the analysis
                    for jj=1:length(IdxFcSel)
                        TileIdx(jj)=find(IdxFcCorr==IdxFcSel(jj));
                        nexttile(TileIdx(jj)) % Access the tile corresponding to the selected force curve
                        ti=title(sprintf('%i',IdxFcSel(jj)),'Color','r'); % Change title color
                        ti.Units='normalized'; % Set units to 'normalized'
                        ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                    end
                end
                
                % Case 4: Select none - except of
                if SelectBttns == 4
                    obj.SMFSFlag.Selected(IdxFcCorr)=0;
                    prompt = {'Enter the force curve number you want want to keep for analysis (For multiple selections just use the space key to separate entries)'};
                    definput = {''};
                    opts.Interpreter = 'tex';
                    IdxFcSel=inputdlg(prompt,'Select none - except of ...',[1 150],definput,opts); % Stores the individual selected fc as a cell array of character vectors
                    IdxFcSel=str2num(IdxFcSel{1}); % Convert the cell array to numerals
                    obj.SMFSFlag.Selected(IdxFcSel)=1; % Set flag for the selected force curves                    
                    % 'Red' loop - Repaint tile titles for force curves
                    % removed from the analysis
                    for ii=1:length(IdxFcCorr)
                        nexttile(ii) % Access the tile corresponding to the selected force curve
                        ti=title(sprintf('%i',IdxFcCorr(ii)),'Color','r'); % Change title color
                        ti.Units='normalized'; % Set units to 'normalized'
                        ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                    end                    
                    % 'Blue' loop - Repaint tile titles for force curves
                    % kept for the analysis                    
                    for jj=1:length(IdxFcSel)
                        TileIdx(jj)=find(IdxFcCorr==IdxFcSel(jj));
                        nexttile(TileIdx(jj)) % Access the tile corresponding to the selected force curve
                        ti=title(sprintf('%i',IdxFcSel(jj)),'Color','b'); % Change title color
                        ti.Units='normalized'; % Set units to 'normalized'
                        ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                    end                    
                end
                %% Save figures
                %%% Define the name for the figure title
                partname=sprintf('-p%d',kk);
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
            Res=[1 1 2560 1440]; % Define the figure resolution
            % Define variables for the figure name
            ExtendVelocityConvert=num2str(obj.ExtendVelocity*1e9);
            RetractVelocityConvert=num2str(obj.RetractVelocity*1e9);
            HoldingTime=num2str(obj.HoldingTime);
            % Classification criteria
            figname=strcat(obj.Date,{'_'},obj.Time,{'_'},obj.ID,{'_'},obj.Substrate,{'_'},obj.EnvCond,{'_'},obj.Linker,{'_'},obj.Chipbox,{'_'},obj.ChipCant,{'_'},ExtendVelocityConvert,{'_'},RetractVelocityConvert,{'_'},HoldingTime);
            figname=char(figname);
            %% Figure loop
            for kk=1:NFigures
                % Figure
                h_fig=figure(kk);
                h_fig.Color='white'; % changes the background color of the figure
                h_fig.Units='pixel'; % Defines the units
                h_fig.OuterPosition=Res;
                h_fig.PaperOrientation='landscape';
                h_fig.Name=figname;
                %% Plotting the tiles
                t = tiledlayout(5,5);
                %t.TileSpacing = 'compact';
                %t.Padding = 'compact';
                t.TileSpacing = 'none'; % To reduce the spacing between the tiles
                t.Padding = 'none'; % To reduce the padding of perimeter of a tile
                
                % Defining variables
                if kk==NFigures && Remainder~=0
                    NLoop=Remainder;
                else
                    NLoop=25;
                end
                %% Plot loop
                for jj=1:NLoop
                    %% Debugging
                    %for jj=4 % for debugging
                    % sprintf('Force curve No. %d',jj) % Gives current Force curve for debugging
                    % Tile jj
                    qq=jj+25*(kk-1);
                    % Allocate data
                    xApp=obj.HHApp{qq}; 
                    xRet=obj.HHRet{qq}; 
                    yApp=obj.App{qq};
                    yRet=obj.Ret{qq}; 
                    % Plot tile
                    nexttile;                 
                    hold on
                    grid on
                    plot(xApp,yApp,'Color',RGB1);
                    plot(xRet,yRet,'Color',RGB2);
                    % Title for each Subplot
                    ti=title(sprintf('%i',qq),'Color','k');
                    ti.Units='normalized'; % Set units to 'normalized'
                    ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                end
                %% Save figures
                %%% Define the name for the figure title
                partname=sprintf('-p%d',kk);
                % fullname=sprintf('%s%s',figname,partname);
                fullname=sprintf('%s%s',figname,partname);
                %%% Save the current figure in the current folder
                print(gcf,fullname,'-dpng');
                
            end
            close all
        end
        
        
        function fc_print_fitted(obj,XMin,XMax,YMin,YMax,Flags) % fc ... force curve
            % fc_print_raw: A function to simply plot all force curves of a
            % force map without any selection taking place
            
          
            % Define variables
            RGB1=[0 26 255]./255;  % Blue
            RGB2=[255 119 0]./255; % Orange
            RGB3=[80 220 100]./255; % Emerald
            RGB4=[200 81 160]./255; % Compl to emerald
            RGB7=[255 230 0]./255; % Yellow
            RGB9=[0 181 26]./255; % RAL6038 Green
            RGB10=[69 22 113]./255; % Violet
            RGB13=[200 8 0]./255; % Red
            DataShareStartApp=0.15; % 15%
            DataShareEndApp=0.65; % 65%
            DataShareStartRet=0.07; % 7 %
            DataShareEndRet=0.02; % 2 %
            NumFcMax=25;
            % Parse unit scale function
            [Xmultiplier,Xunit,~] = AFMImage.parse_unit_scale(1e+9,'nm',1);
            [Ymultiplier,Yunit,~] = AFMImage.parse_unit_scale(1e+9,'nN',1);
            % Define variables for the figure name
            ExtendVelocityConvert=num2str(obj.ExtendVelocity*1e9);
            RetractVelocityConvert=num2str(obj.RetractVelocity*1e9);
            HoldingTimeConvert=num2str(obj.HoldingTime);
            % Classification criteria
            figname=strcat(obj.Date,{'_'},obj.Time,{'_'},obj.ID,{'_'},obj.Substrate,{'_'},obj.EnvCond,{'_'},obj.Linker,{'_'},obj.Chipbox,{'_'},obj.ChipCant,{'_'},ExtendVelocityConvert,{'_'},RetractVelocityConvert,{'_'},HoldingTimeConvert);
            figname=char(figname);
            %% Define variables for the plot loop
            mm=ceil(sqrt(NumFcMax)); % Variable defines the tiledlayout and thereby the tile arrangement
            nn=mm; % Variable defines the tiledlayout and thereby the tile arrangement
            ww=1; % Variable used to loop through the flag based while loop
            DiffFc=0; % Variable to correct for differences of the variables between the plot loop and flag based while loop
            NumFcUncorrupt=nnz(obj.SMFSFlag.Uncorrupt.*obj.SMFSFlag.RetMinCrit.*obj.SMFSFlag.LengthRequisite); % Determine the number of force curves that could been analysed
            if strcmpi(Flags,'Presorted')
            NumFigures=ceil(NumFcUncorrupt/NumFcMax);
            elseif strcmpi(Flags,'Fit')
            NumFigures=4;
            end
            if NumFigures==0     % If condition is fulfilled stop function and return to calling function
                return
            end
            RemainderMax=mod(NumFcUncorrupt,NumFcMax); % Check for remainder
            if RemainderMax ~= 0
                % Additional variables if there is a remainder
                oo=round(sqrt(RemainderMax)); % Determine the number of rows in the figure
                pp=ceil(sqrt(RemainderMax)); % Determine the number of columns in the figure
                RemainderReal=mod(NumFcUncorrupt,oo*pp); % Correct the remainder based on the determined rows times columns
            end
            %% Figure loop
            for kk=1:NumFigures
                % Figure
                h_fig=figure(kk);
                h_fig.Color='white'; % changes the background color of the figure
                h_fig.Units='normalized'; % Defines the units
                h_fig.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
                h_fig.PaperOrientation='landscape';
                h_fig.Name=figname;
                %% Verify the remainder condition
                if RemainderMax == 0 % No remainder existing
                    t = tiledlayout(mm,nn);
                    %t.TileSpacing = 'compact';
                    %t.Padding = 'compact';
                    t.TileSpacing = 'none'; % To reduce the spacing between the tiles
                    t.Padding = 'none'; % To reduce the padding of perimeter of a tile
                    if kk==1
                        qq=1;
                    else
                        qq=(kk-1)*mm*nn+1; % Determines the index of the next force curve to be plotted on the next figure
                    end
                    %% Plot loop
                    for qq=qq:obj.NCurves % Loop over all force curves in the force map
                        if ww<qq+DiffFc
                            ww=qq+DiffFc;
                        end
                        % Flag based while loop
                        if strcmpi(Flags,'Presorted')
                        while ~obj.SMFSFlag.Uncorrupt(ww) || ~obj.SMFSFlag.RetMinCrit(ww) || ~obj.SMFSFlag.LengthRequisite(ww)   % Stay in the while loop as long as the entry is zero
                            ww=ww+1;
                            if ww>qq
                                DiffFc=ww-qq;
                            end
                        end
                        % Allocate data
                        xApp=(obj.THApp{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                        xRet=(obj.THRet{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                        yApp=obj.BasedApp{qq+DiffFc}/Ymultiplier;
                        yRet=obj.BasedRet{qq+DiffFc}/Ymultiplier;
                        % Determine indices
                        IdxStartApp=round(length(yApp)*DataShareStartApp);
                        IdxEndApp=round(length(yApp)*DataShareEndApp);
                        IdxStartRet=round(length(yRet)-length(yRet)*DataShareStartRet);
                        IdxEndRet=round(length(yRet)-length(yRet)*DataShareEndRet);
                        % Plot tile
                        ax=nexttile;
                        ax.XLim = [XMin XMax];
                        ax.YLim = [YMin YMax];
                        hold on
                        grid on
                        plot(xApp,yApp,'Color',RGB1);
                        plot(xRet,yRet,'Color',RGB2);
                        plot(xApp(IdxStartApp:IdxEndApp),yApp(IdxStartApp:IdxEndApp),'Color',RGB13);
                        plot(xRet(IdxStartRet:IdxEndRet),yRet(IdxStartRet:IdxEndRet),'Color','k');
                        % Title for each Subplot
                        ti=title(sprintf('%i',qq+DiffFc),'Color','k');
                        ti.Units='normalized'; % Set units to 'normalized'
                        ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                        if qq == NumFcMax*kk
                            break
                        end
                        elseif strcmpi(Flags,'Fit')
                        while ~obj.SMFSFlag.Fit(ww)   % Stay in the while loop as long as the entry is zero
                            ww=ww+1;
                            if ww>qq
                                DiffFc=ww-qq;
                            end
                        end
                        % Allocate data
                        xApp=(obj.THApp{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                        xRet=(obj.THRet{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                        yApp=obj.BasedApp{qq+DiffFc}/Ymultiplier;
                        yRet=obj.BasedRet{qq+DiffFc}/Ymultiplier;
                        % Determine indices
                        IdxStartApp=round(length(yApp)*DataShareStartApp);
                        IdxEndApp=round(length(yApp)*DataShareEndApp);
                        IdxStartRet=round(length(yRet)-length(yRet)*DataShareStartRet);
                        IdxEndRet=round(length(yRet)-length(yRet)*DataShareEndRet);
                        % Plot tile
                        ax=nexttile;
                        ax.XLim = [XMin XMax];
                        ax.YLim = [YMin YMax];
                        hold on
                        grid on
                        plot(xApp,yApp,'Color',RGB1);
                        plot(xRet,yRet,'Color',RGB2);
                        plot(xApp(IdxStartApp:IdxEndApp),yApp(IdxStartApp:IdxEndApp),'Color',RGB13);
                        plot(xRet(IdxStartRet:IdxEndRet),yRet(IdxStartRet:IdxEndRet),'Color','k');
                        % Title for each Subplot
                        ti=title(sprintf('%i',qq+DiffFc),'Color','k');
                        ti.Units='normalized'; % Set units to 'normalized'
                        ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                        if qq == NumFcMax*kk
                            break
                        end
                        end
                    end
                else % Remainder exisiting
                    %% Last figure conditions
                    if kk~=NumFigures
                        t = tiledlayout(mm,nn);
                        t.TileSpacing = 'none'; % To reduce the spacing between the tiles
                        t.Padding = 'none'; % To reduce the padding of perimeter of a tile
                        if kk==1
                            qq=1;
                        else
                            qq=(kk-1)*mm*nn+1; % Determines the index of the next force curve to be plotted on the next figure
                        end
                        %% Plot loop
                        for qq=qq:obj.NCurves % Loop over all force curves in the force map
                            if ww<qq+DiffFc
                                ww=qq+DiffFc;
                            end
                            % Flag based while loop
                            while ~obj.SMFSFlag.Uncorrupt(ww) || ~obj.SMFSFlag.RetMinCrit(ww) || ~obj.SMFSFlag.LengthRequisite(ww)   % Stay in the while loop as long as the entry is zero
                                ww=ww+1;
                                if ww>qq
                                    DiffFc=ww-qq;
                                end
                            end
                            % Allocate data
                            xApp=obj.HHApp{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc); % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                            xRet=obj.HHRet{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc); % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                            yApp=obj.BasedApp{qq+DiffFc};
                            yRet=obj.BasedRet{qq+DiffFc};
                            % Determine indices
                            IdxStartApp=round(length(yApp)*DataShareStartApp);
                            IdxEndApp=round(length(yApp)*DataShareEndApp);
                            IdxStartRet=round(length(yRet)-length(yRet)*DataShareStartRet);
                            IdxEndRet=round(length(yRet)-length(yRet)*DataShareEndRet);
                            % Plot tile
                            ax=nexttile;
                            ax.XLim = [XMin XMax];
                            ax.YLim = [YMin YMax];
                            hold on
                            grid on
                            plot(xApp,yApp,'Color',RGB1);
                            plot(xRet,yRet,'Color',RGB2);
                            plot(xApp(IdxStartApp:IdxEndApp),yApp(IdxStartApp:IdxEndApp),'Color',RGB13);
                            plot(xRet(IdxStartRet:IdxEndRet),yRet(IdxStartRet:IdxEndRet),'Color','k');
                            % Title for each Subplot
                            ti=title(sprintf('%i',qq+DiffFc),'Color','k');
                            ti.Units='normalized'; % Set units to 'normalized'
                            ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                            if qq == NumFcMax*kk
                                break
                            end
                        end
                    else % kk==NumFigures % corresponds to the last figure plotted
                        t = tiledlayout(oo,pp);
                        t.TileSpacing = 'none'; % To reduce the spacing between the tiles
                        t.Padding = 'none'; % To reduce the padding of perimeter of a tile
                        NumFcPlot=oo*pp;
                        if kk==1
                            qq=1;
                        else
                            qq=(kk-1)*mm*nn+1; % Determines the index of the next force curve to be plotted on the next figure
                        end
                        %% Plot loop
                        for qq=qq:obj.NCurves % Loop over all force curves in the force map
                            if ww<qq+DiffFc
                                ww=qq+DiffFc;
                            end
                            if strcmpi(Flags,'Presorted')
                            while ~obj.SMFSFlag.Uncorrupt(ww) || ~obj.SMFSFlag.RetMinCrit(ww) || ~obj.SMFSFlag.LengthRequisite(ww)    % Stay in the while loop as long as the entry is zero
                                ww=ww+1;
                                if ww>qq
                                    DiffFc=ww-qq;
                                end
                            end
                            % Allocate data
                        xApp=(obj.THApp{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                        xRet=(obj.THRet{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                        yApp=obj.BasedApp{qq+DiffFc}/Ymultiplier;
                        yRet=obj.BasedRet{qq+DiffFc}/Ymultiplier;
                            % Determine indices
                            IdxStartApp=round(length(yApp)*DataShareStartApp);
                            IdxEndApp=round(length(yApp)*DataShareEndApp);
                            IdxStartRet=round(length(yRet)-length(yRet)*DataShareStartRet);
                            IdxEndRet=round(length(yRet)-length(yRet)*DataShareEndRet);
                            % Plot tile
                            ax=nexttile;
                            ax.XLim = [XMin XMax];
                            ax.YLim = [YMin YMax];
                            hold on
                            grid on
                            plot(xApp,yApp,'Color',RGB1);
                            plot(xRet,yRet,'Color',RGB2);
                            plot(xApp(IdxStartApp:IdxEndApp),yApp(IdxStartApp:IdxEndApp),'Color',RGB13);
                            plot(xRet(IdxStartRet:IdxEndRet),yRet(IdxStartRet:IdxEndRet),'Color','k');
                            % Title for each Subplot
                            ti=title(sprintf('%i',qq+DiffFc),'Color','k');
                            ti.Units='normalized'; % Set units to 'normalized'
                            ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                            if qq == NumFcUncorrupt
                                break
                            end
                            elseif strcmpi(Flags,'Fit')
                            while ~obj.SMFSFlag.Fit(ww)    % Stay in the while loop as long as the entry is zero
                                ww=ww+1;
                                if ww>qq
                                    DiffFc=ww-qq;
                                end
                            end
                            % Allocate data
                        xApp=(obj.THApp{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                        xRet=(obj.THRet{qq+DiffFc}-obj.CP_HardSurface(qq+DiffFc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
                        yApp=obj.BasedApp{qq+DiffFc}/Ymultiplier;
                        yRet=obj.BasedRet{qq+DiffFc}/Ymultiplier;
                            % Determine indices
                            IdxStartApp=round(length(yApp)*DataShareStartApp);
                            IdxEndApp=round(length(yApp)*DataShareEndApp);
                            IdxStartRet=round(length(yRet)-length(yRet)*DataShareStartRet);
                            IdxEndRet=round(length(yRet)-length(yRet)*DataShareEndRet);
                            % Plot tile
                            ax=nexttile;
                            ax.XLim = [XMin XMax];
                            ax.YLim = [YMin YMax];
                            hold on
                            grid on
                            plot(xApp,yApp,'Color',RGB1);
                            plot(xRet,yRet,'Color',RGB2);
                            plot(xApp(IdxStartApp:IdxEndApp),yApp(IdxStartApp:IdxEndApp),'Color',RGB13);
                            plot(xRet(IdxStartRet:IdxEndRet),yRet(IdxStartRet:IdxEndRet),'Color','k');
                            % Title for each Subplot
                            ti=title(sprintf('%i',qq+DiffFc),'Color','k');
                            ti.Units='normalized'; % Set units to 'normalized'
                            ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                            if qq == NumFcUncorrupt
                                break
                            end
                            end
                         end
                    end                   
                end
                %% Save figures
                    %%% Define the name for the figure title
                    partname=sprintf('-p%d',kk);
                    % fullname=sprintf('%s%s',figname,partname);
                    fullname=sprintf('%s%s',figname,partname);
                    %%% Save the current figure in the current folder
                    print(gcf,fullname,'-dpng');               
            end
            close all
        end
        
        function fc_flag_status(obj)
            % Find not processed force curves

            obj.SMFSFlagDown.FitSensitivity=find(~obj.SMFSFlag.FitSensitivity);
            obj.SMFSFlagDown.FitSinoidal=find(~obj.SMFSFlag.FitSinoidal);
            obj.SMFSFlagDown.FitLinear=find(~obj.SMFSFlag.FitLinear);
            obj.SMFSFlagDown.FitLinearRet=find(~obj.SMFSFlag.FitLinearRet);
            obj.SMFSFlagDown.Fit=find(~obj.SMFSFlag.Fit);
            obj.SMFSFlagDown.Uncorrupt=find(~obj.SMFSFlag.Uncorrupt);
            obj.SMFSFlagDown.AppMinCrit=find(~obj.SMFSFlag.AppMinCrit);
            obj.SMFSFlagDown.RetMinCrit=find(~obj.SMFSFlag.RetMinCrit);
            obj.SMFSFlagDown.LengthRequisite=find(~obj.SMFSFlag.LengthRequisite);
            obj.SMFSFlagDown.SnapIn=find(~obj.SMFSFlag.SnapIn);
            obj.SMFSFlagDown.PullingLength=find(~obj.SMFSFlag.PullingLength);
            obj.SMFSFlagDown.Selected=find(~obj.SMFSFlag.Selected);

        end

        function [CS1,CS2,CS3,CS4,CS5,CS6,CS7,CS8,CS9,CS10]=fc_ColorScheme(obj)
            % Color scheme from colorbrewer2.org
            % https://colorbrewer2.org/#type=diverging&scheme=RdYlBu&n=10
            % Properties: colorblind safe
            % Number of data classes (colors): 10
            
            CS1=[165 0 38]./255; % Dark reddish
            CS2=[215 48 39]./255; % Light reddish
            CS3=[244 109 67]./255; % Orangish
            CS4=[253 174 97]./255; % Ochreish
            CS5=[254 224 144]./255; % Yellowish
            CS6=[224 243 248]./255; % Pastel blueish
            CS7=[171 217 233]./255; % Light blueish
            CS8=[116 173 209]./255; % Steel blueish
            CS9=[69 117 180]./255; % Distant blueish
            CS10=[49 54 149]./255; % Pale ultramarineish
            
        end
        
        function fc_fine_figure(obj,XMin,XMax,YMin,YMax,Fm,Fc,Linker)
           % function fc_fine_figure(obj,XMin,XMax,YMin,YMax,ii)

            if nargin < 2
                XMin= -inf;
                XMax= inf;
                YMin= -inf;
                YMax= inf;
            end

            % Define RGB colours
            % Define variables
            HEX8DB600=[141 182 0]./255; % Apple green
            HEXCE1620=[206 22 32]./255; % Fire Engine Red
            HEX8F00FF=[143 0 255]./255; % Violet
            RGB1=[0 26 255]./255;  % Blue
            RGB2=[255 119 0]./255; % Orange
            RGB7=[255 230 0]./255; % Yellow
            RGB8=[80 200 204]./255; % Turquoise
            RGB10=[200 0 255]./255; % Violet
            RGB11=[200 255 150]./255; % Light Green
            RGB12=[185 230 254]./255; % Light Blue
            RGB13=[200 0 0]./255; % Red
            RGB_A4=[66 255 0]./255; % Green
            CS1=[165 0 38]./255; % Dark reddish
            CS2=[215 48 39]./255; % Light reddish
            CS3=[244 109 67]./255; % Orangish
            CS4=[253 174 97]./255; % Ochreish
            CS5=[254 224 144]./255; % Yellowish
            CS6=[224 243 248]./255; % Pastel blueish
            CS7=[171 217 233]./255; % Light blueish
            CS8=[116 173 209]./255; % Steel blueish
            CS9=[69 117 180]./255; % Distant blueish
            CS10=[49 54 149]./255; % Pale ultramarineish
            %% General variables 1

            TheoretCollLength1=310;
            TheoretCollLength2=463;
            if strcmpi(Linker,'Long')
            LimitLengthRet1=378; 
            LimitLengthRet2=522;
            LimitLengthRetLabel1='378';
            LimitLengthRetLabel2='522';
            LimitLengthApp1=50;
            LimitLengthApp2=120;
            LimitLengthAppLabel1='50';
            LimitLengthAppLabel2='120';
            elseif strcmpi(Linker,'Short')
            LimitLengthRet1=333; 
            LimitLengthRet2=463;
            LimitLengthRetLabel1='333';
            LimitLengthRetLabel2='463';
            LimitLengthApp1=50;
            LimitLengthApp2=120;
            LimitLengthAppLabel1='50';
            LimitLengthAppLabel2='120';
            end
            Res=[1 1 2560 1440]; % Define the figure resolution
            % Parse unit scale function
            [Xmultiplier,Xunit,~] = AFMImage.parse_unit_scale(1e+9,'nm',1);
            [Ymultiplier,Yunit,~] = AFMImage.parse_unit_scale(1e+9,'nN',1);
            % Define variables for the figure name         
            ExtendVelocityConvert=num2str(round(obj.ExtendVelocity,8)*1e9);
            RetractVelocityConvert=num2str(round(obj.RetractVelocity,8)*1e9);
            HoldingTimeConvert=num2str(round(obj.HoldingTime,2));
            FcNumConvert=num2str(Fc);
            % Classification criteria
            figname=strcat(obj.Date,{'_'},obj.Time,{'_'},obj.ID,{'_'},'Fc',FcNumConvert,{'_'},obj.Substrate,{'_'},obj.EnvCond,{'_'},obj.Linker,{'_'},obj.Chipbox,{'_'},obj.ChipCant,{'_'},ExtendVelocityConvert,{'_'},RetractVelocityConvert,{'_'},HoldingTimeConvert);
            figname=char(figname);
            %% Allocate data       
            xApp=(obj.THApp{Fc}-obj.CP_HardSurface(Fc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
            xRet=(obj.THRet{Fc}-obj.CP_HardSurface(Fc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
            yApp=obj.BasedApp{Fc}/Ymultiplier;
            yRet=obj.BasedRet{Fc}/Ymultiplier;
  %          yAppLim=obj.yAppLim{Fc}/Ymultiplier;
  %          yRetLim=obj.yRetLim{Fc}/Ymultiplier;
            %      xPolygon1=[-10 1000 1000 -10];
            %      yPolygon1=[0 0 -14e-3 -14e-3];
            %      yPolygon2=[-14e-3 -14e-3 -0.9 -0.9];
            %      PolygonShape1=polyshape(xPolygon1,yPolygon1);
            %      PolygonShape2=polyshape(xPolygon1,yPolygon2);
            %
            %% Figure
            h_fig=figure(Fc);
            h_fig.Color='white'; % changes the background color of the figure
            h_fig.Units='normalized'; % Defines the units
            h_fig.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
            %h_fig.Units='pixel'; % Defines the units
            %h_fig.OuterPosition=res;
            h_fig.PaperOrientation='landscape';
            h_fig.Name=figname;
            % Plot
            hold on
             grid on
%              if obj.PullingLengthIdx(Fc)
%              area(xRet(1:obj.PullingLengthIdx(Fc)),yRetLim(1:obj.PullingLengthIdx(Fc)),'FaceColor',CS5)
%              else 
%              end
%              if obj.SnapInIdx(Fc)
%              area(xApp(obj.SnapInIdx(Fc):end),yAppLim(obj.SnapInIdx(Fc):end),'FaceColor',CS7)
%              else
%             end    
            plot(xApp,yApp,'Color',RGB1,'LineWidth',6);
            plot(xRet,yRet,'Color',RGB2,'LineWidth',6);
            if strcmpi(Linker,'Long')
            xlineLabel={LimitLengthRetLabel1,LimitLengthRetLabel2};
            xline([LimitLengthRet1,LimitLengthRet2],'--',xlineLabel,'LabelOrientation','horizontal','FontSize',32,'LineWidth',6);
            elseif strcmpi(Linker,'short')
            xlineLabel={LimitLengthRetLabel1,LimitLengthRetLabel2};
            xline([LimitLengthRet1,LimitLengthRet2],'--',xlineLabel,'LabelOrientation','horizontal','FontSize',32,'LineWidth',6);
            end
            xlineLabel={LimitLengthAppLabel1,LimitLengthAppLabel2};
            xline([LimitLengthApp1,LimitLengthApp2],'--',xlineLabel,'LabelOrientation','horizontal','FontSize',32,'LineWidth',6);
            %            Polygon1=plot(PolygonShape1);
            %            Polygon1.FaceColor=CS10;
            %            Polygon1.EdgeColor='none';            
            % Legend
            % le=legend(' Adhesion energy',' Approach data',' Retraction data','Location','best');
            % le.FontSize = 48;
            % le.EdgeColor='w';
            %le.Box = 'off';
            %%% Axes
            ax = gca; % current axes
            ax.FontSize = 46;
            ax.LineWidth = 5;
            %   ax.XTick=0:100:400;
            % ax.XTickLabel=[];
            %    ax.YTick=-0.3:0.1:0.2;
            % ax.YTickLabel=[];
            ax.XLabel.String = 'Tip-Substrate separation (nm)';
            ax.XLabel.FontSize = 46;
            ax.YLabel.String = 'Force (nN)';
            ax.YLabel.FontSize = 46;
            ax.YLimMode='manual';
            ax.XLim = [XMin XMax];
            ax.YLim = [YMin YMax];
            % Title
            ax.Title.String=sprintf('Force map %d Force-distance curve %d',Fm,Fc);
            ax.Subtitle.String=strcat(obj.Date,{' '},obj.Time,{' '},obj.ID);
            ax.Subtitle.FontSize = 32;
%            ax.Title.Color=HEX8DB600;
%             ti=title(sprintf('Force map %d Force-distance curve %d',Fm,Fc));
%             ti.Units='normalized'; % Set units to 'normalized'
%             ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
            %% Save figures
            %%% Define the name for the figure title
            partname=sprintf('-ForceCurve%d',Fc);
            % fullname=sprintf('%s%s',figname,partname);
            fullname=sprintf('%s%s',figname,partname);
            %%% Save the current figure in the current folder
            print(gcf,fullname,'-dpng');
            % House keeping
            close all
        end

        function fc_fine_figure_publication(obj,XMin,XMax,YMin,YMax,Fm,Fc,Linker,CArea)
           % function fc_fine_figure(obj,XMin,XMax,YMin,YMax,ii)

            if nargin < 2
                XMin= -inf;
                XMax= inf;
                YMin= -inf;
                YMax= inf;
            end

            % Define RGB colours
            % Define variables
            ColorBrewerMap1=[[253 174 97]./255; % Ochreish
               [116 173 209]./255]; % Steel blueish
            HEXfdae61=[253 174 97]./255; % Ochreish, ColorBrewerMap1 color
            HEX74add1=[116 173 209]./255; % Steel blueish, ColorBrewerMap2 color
            HEX8DB600=[141 182 0]./255; % Apple green
            HEXCE1620=[206 22 32]./255; % Fire Engine Red
            HEX8F00FF=[143 0 255]./255; % Violet
            RGB1=[0 26 255]./255;  % Blue
            RGB2=[255 119 0]./255; % Orange
            RGB7=[255 230 0]./255; % Yellow
            RGB8=[80 200 204]./255; % Turquoise
            RGB10=[200 0 255]./255; % Violet
            RGB11=[200 255 150]./255; % Light Green
            RGB12=[185 230 254]./255; % Light Blue
            RGB13=[200 0 0]./255; % Red
            RGB_A4=[66 255 0]./255; % Green
            CS1=[165 0 38]./255; % Dark reddish
            CS2=[215 48 39]./255; % Light reddish
            CS3=[244 109 67]./255; % Orangish
            CS4=[253 174 97]./255; % Ochreish
            CS5=[254 224 144]./255; % Yellowish
            CS6=[224 243 248]./255; % Pastel blueish
            CS7=[171 217 233]./255; % Light blueish
            CS8=[116 173 209]./255; % Steel blueish
            CS9=[69 117 180]./255; % Distant blueish
            CS10=[49 54 149]./255; % Pale ultramarineish
            Res=[1 1 2560 1440]; % Define the figure resolution
            FontName='Helvetica';
            FontSize=46;
            % Parse unit scale function
            [Xmultiplier,Xunit,~] = AFMImage.parse_unit_scale(1e+9,'nm',1);
            [Ymultiplier,Yunit,~] = AFMImage.parse_unit_scale(1e+9,'nN',1);
            % Define variables for the figure name         
            ExtendVelocityConvert=num2str(round(obj.ExtendVelocity,8)*1e9);
            RetractVelocityConvert=num2str(round(obj.RetractVelocity,8)*1e9);
            HoldingTimeConvert=num2str(round(obj.HoldingTime,2));
            FcNumConvert=num2str(Fc);
            % Classification criteria
            figname=strcat(obj.Date,{'_'},obj.Time,{'_'},obj.ID,{'_'},'Fc',FcNumConvert,{'_'},obj.Substrate,{'_'},obj.EnvCond,{'_'},obj.Linker,{'_'},obj.Chipbox,{'_'},obj.ChipCant,{'_'},ExtendVelocityConvert,{'_'},RetractVelocityConvert,{'_'},HoldingTimeConvert);
            figname=char(figname);
            %% Allocate data       
            xApp=(obj.THApp{Fc}-obj.CP_HardSurface(Fc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
            xRet=(obj.THRet{Fc}-obj.CP_HardSurface(Fc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
            yApp=obj.BasedApp{Fc}/Ymultiplier;
            yRet=obj.BasedRet{Fc}/Ymultiplier;
            if nnz(obj.PullingLengthIdx(Fc)) && strcmpi(CArea,'yes')
            yRetLim=obj.yRetLim{Fc}/Ymultiplier;
            else
            end
            if nnz(obj.SnapInIdx(Fc)) && strcmpi(CArea,'yes')
            yAppLim=obj.yAppLim{Fc}/Ymultiplier;
            else
            end
            LengthLong1=378;
            LengthLong2=522;
            LengthShort1=333;
            LengthShort2=463;
            if strcmpi(Linker,'Long')
            xPolygon1=[XMin XMin LengthLong1 LengthLong1];
            xPolygon2=[LengthLong1 LengthLong1 LengthLong2 LengthLong2];
            yPolygon=[YMin YMax YMax YMin];
            PolygonShape1=polyshape(xPolygon1,yPolygon);
            PolygonShape2=polyshape(xPolygon2,yPolygon);
            elseif strcmpi(Linker,'Short')
            xPolygon1=[XMin XMin LengthShort1 LengthShort1];
            xPolygon2=[LengthShort1 LengthShort1 LengthShort2 LengthShort2];
            yPolygon=[YMin YMax YMax YMin];
            PolygonShape1=polyshape(xPolygon1,yPolygon);
            PolygonShape2=polyshape(xPolygon2,yPolygon);
            end
            %% Figure
            % h_fig=figure(ii);
            h_fig=figure(1);
            set(h_fig,'DefaultAxesFontName',FontName);
            set(h_fig,'DefaultTextFontName',FontName);
            set(h_fig,'DefaultAxesFontSize',FontSize);
            set(h_fig,'DefaultTextFontSize',FontSize);
            h_fig.Color='white'; % changes the background color of the figure
            h_fig.Units='normalized'; % Defines the units
            h_fig.OuterPosition=[0 0 2 2];% changes the size of the to the whole screen
            h_fig.Units='pixel'; % Defines the units
    %        h_fig.OuterPosition=Res;
            h_fig.PaperType='a5';
            h_fig.PaperOrientation='landscape';
            h_fig.Name=figname;
            % Plot
            hold on
             grid on
             if nnz(obj.PullingLengthIdx(Fc)) && strcmpi(CArea,'yes')
             area(xRet(1:obj.PullingLengthIdx(Fc)),yRetLim(1:obj.PullingLengthIdx(Fc)),'FaceColor',CS5)
             else 
             end
             if nnz(obj.SnapInIdx(Fc)) && strcmpi(CArea,'yes')
             area(xApp(obj.SnapInIdx(Fc):end),yAppLim(obj.SnapInIdx(Fc):end),'FaceColor',CS7)
             else
             end    
            plot(xApp,yApp,'Color',RGB1,'LineWidth',6);
            plot(xRet,yRet,'Color',RGB2,'LineWidth',6);
            if strcmpi(Linker,'Short') || strcmpi(Linker,'Long')
              Polygon1=plot(PolygonShape1);
              Polygon1.FaceColor=HEXfdae61;
              Polygon1.EdgeColor='none';       
               Polygon2=plot(PolygonShape2);
               Polygon2.FaceColor=HEX74add1;
               Polygon2.EdgeColor='none';  
            elseif strcmpi(Linker,'No')   
            end
            % Legend
            % le=legend(' Adhesion energy',' Approach data',' Retraction data','Location','best');
            % le.FontSize = 48;
            % le.EdgeColor='w';
            %le.Box = 'off';
            %%% Axes
            ax = gca; % current axes                
            %ax.Visible='off';
        %    ax.TickLabelInterpreter='latex';
        %    ax.FontSize = 46;
            ax.LineWidth = 5;
            % ax.XTick=0:100:400;
            % ax.XTickLabel=[];
            %    ax.YTick=-0.3:0.1:0.2;
            % ax.YTickLabel=[];
            ax.XLabel.String = 'Tip-surface separation distance (nm)';
        %    ax.XLabel.Interpreter='latex';
         %   ax.XLabel.FontSize = 46;
       %     ax.XLabel.FontName = 'Helvetica';
            ax.YLabel.String = 'Force (nN)';
       %     ax.YLabel.Interpreter='latex';
          %  ax.YLabel.FontSize = 46;
        %    ax.YLabel.FontName = 'Helvetica';
            ax.YLimMode='manual';
            ax.XLim = [XMin XMax];
            ax.YLim = [YMin YMax];
            ax.XTick=0:100:900;
          %  ax.Position=[0 0 1 1]; % Sets the Position to the same values as OuterPosition
            % Title
%             ax.Title.String=sprintf('Force map %d Force-distance curve %d',Fm,Fc);
%             ax.Subtitle.String=strcat(obj.Date,{' '},obj.Time,{' '},obj.ID);
%             ax.Subtitle.FontSize = 32;
%            ax.Title.Color=HEX8DB600;
%             ti=title(sprintf('Force map %d Force-distance curve %d',Fm,Fc));
%             ti.Units='normalized'; % Set units to 'normalized'
%             ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
            %% Save figures
            %%% Define the name for the figure title
            partname=sprintf('-ForceCurve%d',Fc);
            % fullname=sprintf('%s%s',figname,partname);
            fullname=sprintf('%s%s',figname,partname);
            %%% Save the current figure in the current folder
            print(gcf,fullname,'-dpng');
            print(gcf,fullname,'-depsc');
            % House keeping
            close all
        end

       %% SMFS Old functions 
          
        function fc_adhesion_energy_threshold(obj)
            % Determine the adhesion energy using a predefined force
            % threshold to distinguish interactions from background noise
            
            %% Loop over all force curves
            for kk=1:obj.NCurves
            %for ii=97 % For debugging and testing 
                if ~obj.SMFSFlag.Uncorrupt(kk) || ~obj.SMFSFlag.RetMinCrit(kk) || ~obj.SMFSFlag.LengthRequisite(kk)    % Exclude corrupted force curves from the analysis     
                continue
                end            
                % Allocate data
                xRet=obj.HHRet{kk}-obj.CP_HardSurface(kk); % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method 
                xApp=obj.HHApp{kk}-obj.CP_HardSurface(kk); % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method 
                if obj.SMFSFlag.Fit(kk)==1
                yApp=obj.BasedApp{kk};
                yRet=obj.BasedRet{kk};
                else
                yApp=obj.BasedApp{kk};
                yRet=obj.BasedRetCorr2{kk};
                end
                % Define variables
                limit1=0;   % Define the limit
                % Apply the limit
                % Approach
                yApp(yApp>Limit1)=0;
                yApp(obj.PullingLengthIdx(kk):end)=0; % Set all data points with a higher index (surface distance is higher) than the pulling length index to 0
                % Retention 
                yRet(yRet>limit1)=0;  % Set all values above the zero line of the x-axis 0
                yRet(obj.PullingLengthIdx(kk):end)=0; % Set all data points with a higher index (surface distance is higher) than the pulling length index to 0
                % Allocate data
                obj.yRetLim{kk}=yRet; 
                obj.yAppLim{kk}=yApp; 
                % Determine the adhesion energy
                IntApp(kk)=trapz(yApp,xApp); % Integrates over the modified y-retraction data with respect to the corresponding x-retraction data 
                obj.AppAdhEnergy_IdxMethod(kk)=IntApp(kk);
                IntRet(kk)=trapz(yRet,xRet); % Integrates over the modified y-retraction data with respect to the corresponding x-retraction data 
                obj.RetAdhEnergy_IdxMethod(kk)=IntRet(kk);
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
                   
        function fc_print_scatter(obj,XMin,XMax,YMin,YMax,NumFcMax,hh) % fc ... force curve
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
            NumFcUncorrupt=nnz(obj.SMFSFlag.Uncorrupt.*obj.SMFSFlag.RetMinCrit); % Determine the amount of uncorrupted force curves            
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
            for kk=1:NumFigures               
                % Define variables
                jj=1; % "force curve plotted per figure" variable        
                % Figure
                h_fig=figure(kk);
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
                    if kk==1
                        NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt(1:NumFcMax));
                    elseif kk==2
                        NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt((NumFcMax+1):(NumFcMax*(kk))));
                    else
                        NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt((NumFcMax*(kk-1)+1):(NumFcMax*(kk))));
                    end
                    
                    if kk==1
                        kk=jj;
                    else
                        kk=jj+mm*nn*(kk-1);
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
                           
                            
                            % Title for each Subplot
                            ti=title(sprintf('%i',qq+DiffFc),'Color','k');
                            ti.Units='normalized'; % Set units to 'normalized'
                            ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                        else
                            ax=nexttile;
                            ax.XLim = [XMin XMax];
                            ax.YLim = [YMin YMax];                 
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
                    if kk~=NumFigures
                        t = tiledlayout(mm,nn);
                        t.TileSpacing = 'none'; % To reduce the spacing between the tiles
                        t.Padding = 'none'; % To reduce the padding of perimeter of a tile
                        if kk==1
                            NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt(1:NumFcMax));
                        elseif kk==2
                            NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt((NumFcMax+1):(NumFcMax*(kk))));
                        else
                            NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt((NumFcMax*(kk-1)+1):(NumFcMax*(kk))));
                        end
                        
                        if kk==1
                            kk=jj;
                        else
                            kk=jj+mm*nn*(kk-1);
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
                               
                                % Title for each Subplot
                                ti=title(sprintf('%i',qq+DiffFc),'Color','k');
                                ti.Units='normalized'; % Set units to 'normalized'
                                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                            else
                                ax=nexttile;
                                ax.XLim = [XMin XMax];
                                ax.YLim = [YMin YMax];
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
                        if kk==1
                            NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt(1:NumFcPlot));
                        elseif kk==2
                            NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt((NumFcPlot+1):(NumFcPlot*(kk))));
                        else
                            NumFcCorSelec(kk)=nnz(obj.SMFSFlag.Uncorrupt((NumFcPlot*(kk-1)+1):(NumFcPlot*(kk))));
                        end
                        
                        if kk==1
                            kk=jj;
                        else
                            kk=jj+mm*nn*(kk-1);
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
                                % Title for each Subplot
                                ti=title(sprintf('%i',qq+DiffFc),'Color','k');
                                %ti=title(sprintf('%i',(kk+ww)/2),'Color','k');
                                ti.Units='normalized'; % Set units to 'normalized'
                                ti.Position=[0.5,0.95]; % Position the subplot title within the subplot
                            else
                                ax=nexttile;
                                ax.XLim = [XMin XMax];
                                ax.YLim = [YMin YMax];
                                
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
                partname=sprintf('-p%d',kk);
                % fullname=sprintf('%s%s',figname,partname);
                fullname=sprintf('%s%s',figname,partname);
                %%% Save the current figure in the current folder
                print(gcf,fullname,'-dpng');
            end
            close all
        end
        
        %% SMFS Testing functions
        
           function fc_WLC_fitting(obj,Fc,Fm)



            RGB1=[0 26 255]./255;  % Blue
            RGB2=[255 119 0]./255; % Orange
            RGB8=[80 200 204]./255; % Turquoise
            % Parse unit scale function
            [Xmultiplier,Xunit,~] = AFMImage.parse_unit_scale(1e+9,'nm',1);
            [Ymultiplier,Yunit,~] = AFMImage.parse_unit_scale(1e+9,'nN',1);
            % Define variables for the figure name         
%              ExtendVelocityConvert=num2str(round(.ExtendVelocity,8)*1e9);
%              RetractVelocityConvert=num2str(round(Trial15.RetractVelocity,8)*1e9);
%              HoldingTimeConvert=num2str(round(Trial15.HoldingTime,2));
%              FcNumConvert=num2str(Fc);
            % Define limits
            XMin=-10;
            XMax=500;            
            YMin=-0.5;
            YMax=0.1;
            % Classification criteria
%              figname=strcat(Trial15.Date,{'_'},Trial15.Time,{'_'},Trial15.ID,{'_'},'Fc',FcNumConvert,{'_'},Trial15.Substrate,{'_'},Trial15.EnvCond,{'_'},Trial15.Linker,{'_'},Trial15.Chipbox,{'_'},Trial15.ChipCant,{'_'},ExtendVelocityConvert,{'_'},RetractVelocityConvert,{'_'},HoldingTimeConvert);
%              figname=char(figname);
            %% Allocate data       
            xRet=(obj.THRet{Fc}-obj.CP_HardSurface(Fc))/-Xmultiplier; % Retraction x-data (m): Vertical tip height data corrected by the determined contact point using the hard surface method
            yRet=obj.BasedRet{Fc}/Ymultiplier;

            %% Figure
            % h_fig=figure(ii);
            h_fig=figure(1);
            h_fig.Color='white'; % changes the background color of the figure
            h_fig.Units='normalized'; % Defines the units
            h_fig.OuterPosition=[0 0 1 1];% changes the size of the to the whole screen
            %h_fig.Units='pixel'; % Defines the units
            %h_fig.OuterPosition=res;
            h_fig.PaperOrientation='landscape';
       %     h_fig.Name=figname;
            % Plot
            hold on
            grid on
            plot(xRet,yRet,'Color',RGB2,'LineWidth',6);
            % Axes
            ax = gca; % current axes
            ax.FontSize = 46;
            ax.LineWidth = 5;
            ax.XLabel.String = 'Tip-Substrate separation (nm)';
            ax.XLabel.FontSize = 46;
            ax.YLabel.String = 'Force (nN)';
            ax.YLabel.FontSize = 46;
            ax.XLim = [XMin XMax];
            ax.YLim = [YMin YMax];
            % Draw points
            roi1 = drawpoint(ax); % Draw point in the figure
            xPos1=roi1.Position(1); % Read out the x-value of the drawn point
            roi2 = drawpoint(ax); % Draw point in the figure
            xPos2=roi2.Position(1); % Read out the x-value of the drawn point

            % Determine the lower and higher x-value
            if xPos1<xPos2
                xPtLow=xPos1;
                xPtHigh=xPos2;
            elseif xPos1>xPos2
                xPtLow=xPos2;
                xPtHigh=xPos1;
            end

            % Find indices
            [~,PtLowIdx] = min(abs(xRet-xPtLow));
            [~,PtHighIdx] = min(abs(xRet-xPtHigh));

            % Select Data within the two indices
            xSel=xRet(PtLowIdx:PtHighIdx);
            ySel=yRet(PtLowIdx:PtHighIdx);

            % Plot the selected data points
            plot(xSel,ySel,'Color',RGB8,'LineWidth',6)
            
            close all
        end

        function fc_testing(obj)
                           
            % For loop
            %for jj=1:obj.NCurves
            %% Debugging
            %for jj=33 % for debugging
    
            obj.SMFSFlag.SnapIn=zeros(1,obj.NCurves);
            %end
            
        end  
        
        function [x,y]=test(obj)
            % TEST FUNCTION
            mega=10e6;
                giga=10e9;
                tera=10e12;
            jj=37;
           % App
        %   obj.xDataToFit=obj.HHApp{jj}*mega
         %  obj.yDataToFit=obj.App{jj}*giga
           % Ret
           obj.xDataToFit=obj.HHRet{jj}*mega
           obj.yDataToFit=obj.Ret{jj}*giga
           
           % App
          %  obj.yDataToFit(ceil(end*0.7):end)=[]
          %  obj.xDataToFit(ceil(end*0.7):end)=[]
           % Ret
            obj.yDataToFit(1:ceil(end*0.35))=[]
            obj.xDataToFit(1:ceil(end*0.35))=[]
            
            
    %        x=obj.xDataToFit  
    %        y=obj.yDataToFit
            
        x=e.FM{135}.xDataToFit  
       y=e.FM{135}.yDataToFit
            
            
            
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
            % Author: Orestis Andriotis (slightly changed and adapted by Manuel Rufin)
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
            B=strfind(A,'vDeflection');
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
                        
                        [App,HHApp] = objcell{i}.get_force_curve_data(j,0,1,0);
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
                prompt = {'What is the nominal tip radius of the used cantilever in nm?'};
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
                E_mod = 3*(Hertzfit.a*ranf/rant^(3/2))/(4*sqrt(TipRadius*10^(-9)))*(1-PoissonR^2);
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
                [Force,Height] = obj.get_force_curve_data(i,AppRetSwitch,1,0);
                if (Mask(obj.List2Map(i,1),obj.List2Map(i,2)) == 1) &&...
                        (obj.ExclMask(obj.List2Map(i,1),obj.List2Map(i,2)) == 1)
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
                    k = k + 1;
                end
            end
            % Fit the Gaussian
            Gaussian = fitdist(DZslope','Normal');
            obj.RefSlope = Gaussian.mean;
            obj.RefSlopeCorrectedSensitivity = obj.Sensitivity/obj.RefSlope;
            obj.HasRefSlope = true;
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
            
            if isequal(obj.FileType,'quantitative-imaging-map')
                I = AFMImage(fullfile(obj.DataStoreFolder,'data-image.jpk-qi-image'));
            elseif isequal(obj.FileType,'force-scan-map')
                I = AFMImage(fullfile(obj.DataStoreFolder,'data-image.force'));
                for i=1:length(I.Channel)
                    I.Channel(i).Image = fliplr(I.Channel(i).Image);
                end
            end
            
            obj.Channel = I.Channel;
            
            Height = obj.get_channel('Height');
            obj.HeightMap = Height.Image;
            
            if size(obj.HeightMap,1) < 128
                Map = imresize(obj.HeightMap,[256 256],'nearest');
            elseif size(obj.HeightMap,1) < 512
                Map = imresize(obj.HeightMap,[512 512],'nearest');
            else
                Map = imresize(obj.HeightMap,[1024 1024],'nearest');
            end
            for i=1:5
                Map = AFMImage.subtract_line_fit_vertical_rov(Map,.2,0);
            end
            Map = imresize(Map,[obj.NumPixelsX obj.NumPixelsY],'nearest');
            
            Map = AFMImage.find_and_replace_outlier_lines(Map,10);
            
            % write to Channel
            obj.delete_channel('Processed')
            Processed = obj.create_standard_channel(Map,'Processed','m');
            
            [Channel,Index] = obj.get_channel('Processed');
            if isempty(Channel)
                Len = length(obj.Channel);
                if ~Len
                    obj.Channel = Processed;
                else
                    obj.Channel(Len+1) = Processed;
                end
            else
                obj.Channel(Index) = Processed;
            end
            
            obj.HeightMap = Processed.Image;
            
        end
        
        function create_and_level_height_map(obj)
            % first set Height Map to default values for reproducable
            % results
            
            obj.construct_list_to_map_relations()
            
            Max = zeros(obj.NCurves,1);
            for i=1:obj.NCurves
                [~,HHApp] = obj.get_force_curve_data(i,0,0,0);
                Max(i) = -max(HHApp);
            end
            TempHeightMap = obj.convert_data_list_to_map(Max);
            
            % write to Channel
            obj.delete_channel('Indented Height')
            Height = obj.create_standard_channel(TempHeightMap,'Indented Height','m');
            
            [Channel,Index] = obj.get_channel('Indented Height');
            if isempty(Channel)
                Len = length(obj.Channel);
                if ~Len
                    obj.Channel = Height;
                else
                    obj.Channel(Len+1) = Height;
                end
            else
                obj.Channel(Index) = Height;
            end
            
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
            
            Map = AFMImage.find_and_replace_outlier_lines(Map,10);
            
            % write to Channel
            obj.delete_channel('Processed Indented Height')
            Processed = obj.create_standard_channel(Map,'Processed Indented Height','m');
            
            [Channel,Index] = obj.get_channel('Processed Indented Height');
            if isempty(Channel)
                Len = length(obj.Channel);
                if ~Len
                    obj.Channel = Processed;
                else
                    obj.Channel(Len+1) = Processed;
                end
            else
                obj.Channel(Index) = Processed;
            end
            
        end
        
        function create_and_level_height_map_by_current_cp(obj)
            
            Map = obj.convert_data_list_to_map(-obj.CP(:,1));
            
            for i=1:obj.NumPixelsX
                Map(i,:) = AFMImage.replace_points_of_certain_value_in_line(Map(i,:),0);
            end
            
            for i=1:5
                Map = AFMImage.subtract_line_fit_vertical_rov(Map,.2,0);
            end
            
            Map = AFMImage.find_and_replace_outlier_lines(Map,10);
            
            % Write to Channel
            obj.delete_channel('Contact Height')
            Chan = obj.create_standard_channel(Map,'Contact Height','m');
            
            obj.Channel(end+1) = Chan;
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
                obj.MiniBatchSize = 1024;
                DynMBSdone = false;
                HasFailed = false;
                while CantHandle == true
                    try
                        predict(NeuralNet,X,'MiniBatchSize',obj.MiniBatchSize,'Acceleration','auto');
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
                                rethrow(ME)
                        end
                    end
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
                [App,HHApp] = obj.get_force_curve_data(i,0,1,0);
                CutOff = min(HHApp) + range(HHApp(HHApp<=obj.CP_CNNZoom(i,1)))*ZoomFactor;
                HHApp(HHApp<CutOff) = [];
                App(1:(end-length(HHApp))) = [];
            end
            
            
        end
        
        function unpack_jpk_force_map(obj,MapFullFile,DataFolder)
            
            if obj.BigDataFlag
                TempFolderName = sprintf('FM_DataStore_%s',obj.ID);
            else
                TempFolderName = sprintf('Temp%s',obj.ID);
            end
            
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
                if h==0
                    disp('unzipping successfull')
                elseif h==1
                    disp('unzipping failed')
                end
                Strings = split(MapFullFile,filesep);
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
                %%% Create a data folder to store the force data
                mkdir(DataFolder,'ForceData')
                obj.Folder = fullfile(DataFolder,'ForceData',filesep);
                
                %             system(['unzip -o ', fullfile(datadir,fnamemap), ' ''*shared-data/header.properties'' -d ', tempdir{fib,1}]);
                %
            elseif isequal('GLN',obj.HostOS)
                % unpack jpk-file into temporary folder to read out data
                cmd1 = 'unzip -o ';
                cmd2 = MapFullFile;
                cmd3 = ' -d ';
                mkdir(DataFolder,TempFolderName)
                TempFolder = fullfile(DataFolder,TempFolderName,filesep);
                CMD = append(cmd1,cmd2,cmd3,TempFolder);
                system(CMD);
                disp('extracting file...')
                Strings = split(MapFullFile,filesep);
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
                %%% Create a data folder to store the force data
                mkdir(DataFolder,'ForceData')
                obj.Folder = fullfile(DataFolder,'ForceData',filesep);
                
            elseif isequal('MAC',obj.HostOS)
                % unpack jpk-file into temporary folder to read out data
                cmd1 = 'unzip -o ';
                cmd2 = MapFullFile;
                cmd3 = ' -d ';
                mkdir(DataFolder,TempFolderName)
                TempFolder = fullfile(DataFolder,TempFolderName,filesep);
                CMD = append(cmd1,cmd2,cmd3,TempFolder);
                system(CMD);
                disp('extracting file...')
                Strings = split(MapFullFile,filesep);
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
                %%% Create a data folder to store the force data
                mkdir(DataFolder,'ForceData')
                obj.Folder = fullfile(DataFolder,'ForceData',filesep);
            end
            
            obj.DataStoreFolder = TempFolder;
            
        end
        
        function read_in_header_properties(obj)
            % Check for jpk-software version and get important ForceMap
            % properties
            filedirectory = fullfile(obj.DataStoreFolder,'header.properties');
            FileID=fopen(filedirectory,'rt','n','UTF-8'); % FileID = fopen(filename,permission,machinefmt,encodingIn)
            FileCont=fileread(filedirectory);
            % Shared-data file directory 
            FileDirectoryShared = fullfile(obj.DataStoreFolder,'shared-data','header.properties');
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
            obj.HHType = 'capacitiveSensorHeight';           
            fclose(FileID);
            fclose(FileIDShared);
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
        
        function [OutForce,OutHeight] = get_force_curve_data(obj,CurveNumber,AppRetSwitch,isBased,isTipHeightCorrected)
            % [OutForce,OutHeight] = get_force_curve_data(obj,CurveNumber,AppRetSwitch,isBased,isTipHeightCorrected)
            % CurveNumber ... Index of curve to be loaded (counting from 1 to obj.NCurves)
            % AppRetSwitch ... 0 for approach curve, 1 for retract curve
            % isBased ... bool for base and tilt subtraction
            % isTipHeightCorrected ... bool to decide if OutHeight is
            %                           HeadHeight (=0) or TipHeight (=1)
            
            if CurveNumber > obj.NCurves
                error('Requested curve index (%i) is bigger than maximum number of curves (%i)',CurveNumber,obj.NCurves);
            end
            if (isBased || isTipHeightCorrected) && ~obj.BaseAndTiltFlag
                error('No Base and Tilt data found');
            end
            if isempty(obj.BigDataFlag) || ~obj.BigDataFlag
                if AppRetSwitch==0 && isBased==0 && isTipHeightCorrected==0
                    OutForce = obj.App{CurveNumber};
                    OutHeight = obj.HHApp{CurveNumber};
                elseif AppRetSwitch==1 && isBased==0 && isTipHeightCorrected==0
                    OutForce = obj.Ret{CurveNumber};
                    OutHeight = obj.HHRet{CurveNumber};
                elseif AppRetSwitch==0 && isBased==1 && isTipHeightCorrected==0
                    OutForce = obj.App{CurveNumber};
                    OutHeight = obj.HHApp{CurveNumber};
                    FittedApp = polyval(obj.Basefit{CurveNumber},OutHeight);
                    OutForce = OutForce - FittedApp;
                elseif AppRetSwitch==1 && isBased==1 && isTipHeightCorrected==0
                    OutForce = obj.Ret{CurveNumber};
                    OutHeight = obj.HHRet{CurveNumber};
                    FittedRet = polyval(obj.Basefit{CurveNumber},OutHeight);
                    OutForce = OutForce - FittedRet;
                elseif AppRetSwitch==0 && isBased==0 && isTipHeightCorrected==1
                    OutForce = obj.App{CurveNumber};
                    OutHeight = obj.HHApp{CurveNumber};
                    FittedApp = polyval(obj.Basefit{CurveNumber},OutHeight);
                    OutForce = OutForce - FittedApp;
                    OutHeight = OutHeight - OutForce./obj.SpringConstant;
                    OutForce = obj.App{CurveNumber};
                elseif AppRetSwitch==1 && isBased==0 && isTipHeightCorrected==1
                    OutForce = obj.Ret{CurveNumber};
                    OutHeight = obj.HHRet{CurveNumber};
                    FittedRet = polyval(obj.Basefit{CurveNumber},OutHeight);
                    OutForce = OutForce - FittedRet;
                    OutHeight = OutHeight - OutForce./obj.SpringConstant;
                    OutForce = obj.Ret{CurveNumber};
                elseif AppRetSwitch==0 && isBased==1 && isTipHeightCorrected==1
                    OutForce = obj.App{CurveNumber};
                    OutHeight = obj.HHApp{CurveNumber};
                    FittedApp = polyval(obj.Basefit{CurveNumber},OutHeight);
                    OutForce = OutForce - FittedApp;
                    OutHeight = OutHeight - OutForce./obj.SpringConstant;
                elseif AppRetSwitch==1 && isBased==1 && isTipHeightCorrected==1
                    OutForce = obj.Ret{CurveNumber};
                    OutHeight = obj.HHRet{CurveNumber};
                    FittedRet = polyval(obj.Basefit{CurveNumber},OutHeight);
                    OutForce = OutForce - FittedRet;
                    OutHeight = OutHeight - OutForce./obj.SpringConstant;
                end
                return
            end
            
            TempFolder = obj.DataStoreFolder;
            
            % Correct for potential holding segment
            if AppRetSwitch
                AppRetSwitch=obj.NumSegments-1;
            end
            
            HeaderFileDirectory = fullfile(TempFolder,'shared-data','header.properties');
            SegmentHeaderFileDirectory = fullfile(TempFolder,'index',string((CurveNumber-1)),'segments',string(AppRetSwitch),'segment-header.properties');
            HeightDataDirectory = fullfile(TempFolder,'index',string((CurveNumber-1)),'segments',string(AppRetSwitch),'channels','capacitiveSensorHeight.dat');
            vDefDataDirectory = fullfile(TempFolder,'index',string((CurveNumber-1)),'segments',string(AppRetSwitch),'channels','vDeflection.dat');
            
            if ~isfile(HeightDataDirectory) || isequal(obj.HHType,'measuredHeight')
                HeightDataDirectory = fullfile(TempFolder,'index',string((CurveNumber-1)),'segments',string(AppRetSwitch),'channels','measuredHeight.dat');
                obj.HHType = 'measuredHeight';
            end
            if ~isfile(HeightDataDirectory) || isequal(obj.HHType,'Height')
                HeightDataDirectory = fullfile(TempFolder,'index',string((CurveNumber-1)),'segments',string(AppRetSwitch),'channels','Height.dat');
                obj.HHType = 'Height';
            end
            
            try
                if obj.CorruptedCurves(CurveNumber)
                    error('skip this one')
                end
                [TempHeight,OutForce,SC,Sens]=...
                    obj.writedata(HeaderFileDirectory,SegmentHeaderFileDirectory,...
                    HeightDataDirectory,vDefDataDirectory,obj.HHType);
                
                if isempty(obj.Sensitivity)
                    obj.Sensitivity = Sens;
                end
                if isempty(obj.SpringConstant)
                    obj.SpringConstant = SC;
                end
                
                OutHeight = -TempHeight;
                OutForce = OutForce.*obj.SpringConstant;
                
                if AppRetSwitch==0 && isBased==0 && isTipHeightCorrected==0
                elseif AppRetSwitch==1 && isBased==0 && isTipHeightCorrected==0
                elseif AppRetSwitch==0 && isBased==1 && isTipHeightCorrected==0
                    FittedApp = polyval(obj.Basefit{CurveNumber},OutHeight);
                    OutForce = OutForce - FittedApp;
                elseif AppRetSwitch==1 && isBased==1 && isTipHeightCorrected==0
                    FittedRet = polyval(obj.Basefit{CurveNumber},OutHeight);
                    OutForce = OutForce - FittedRet;
                elseif AppRetSwitch==0 && isBased==0 && isTipHeightCorrected==1
                    TempApp = OutForce;
                    FittedApp = polyval(obj.Basefit{CurveNumber},OutHeight);
                    OutForce = OutForce - FittedApp;
                    OutHeight = OutHeight - OutForce./obj.SpringConstant;
                    OutForce = TempApp;
                elseif AppRetSwitch==1 && isBased==0 && isTipHeightCorrected==1
                    TempRet = OutForce;
                    FittedRet = polyval(obj.Basefit{CurveNumber},OutHeight);
                    OutForce = OutForce - FittedRet;
                    OutHeight = OutHeight - OutForce./obj.SpringConstant;
                    OutForce = TempRet;
                elseif AppRetSwitch==0 && isBased==1 && isTipHeightCorrected==1
                    FittedApp = polyval(obj.Basefit{CurveNumber},OutHeight);
                    OutForce = OutForce - FittedApp;
                    OutHeight = OutHeight - OutForce./obj.SpringConstant;
                elseif AppRetSwitch==1 && isBased==1 && isTipHeightCorrected==1
                    FittedRet = polyval(obj.Basefit{CurveNumber},OutHeight);
                    OutForce = OutForce - FittedRet;
                    OutHeight = OutHeight - OutForce./obj.SpringConstant;
                end
            catch
                disp(sprintf('Curve Nr. %i seems to be corrupted. Replacing with next non corrupted curve instead',CurveNumber))
                obj.CorruptedCurves(CurveNumber) = true;
                obj.SelectedCurves(CurveNumber) = false;
                k = 1;
                while obj.CorruptedCurves(mod(CurveNumber+k-1,obj.NCurves)+1)
                    k = k + 1;
                end
                [OutForce,OutHeight] = obj.get_force_curve_data(mod(CurveNumber+k-1,obj.NCurves)+1,AppRetSwitch,isBased,isTipHeightCorrected);
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
            obj.SMFSFlag.FitSensitivity=zeros(1,obj.NCurves);
            obj.SMFSFlag.FitSinoidal=zeros(1,obj.NCurves);
            obj.SMFSFlag.FitLinear=zeros(1,obj.NCurves);
            obj.SMFSFlag.FitLinearRet=zeros(1,obj.NCurves);
            obj.SMFSFlag.Fit=zeros(1,obj.NCurves);
            obj.SMFSFlag.Uncorrupt=ones(1,obj.NCurves);
            obj.SMFSFlag.AppMinCrit=zeros(1,obj.NCurves);
            obj.SMFSFlag.RetMinCrit=zeros(1,obj.NCurves);
            obj.SMFSFlag.LengthRequisite=zeros(1,obj.NCurves);
            obj.SMFSFlag.SnapIn=zeros(1,obj.NCurves);
            obj.SMFSFlag.PullingLength=zeros(1,obj.NCurves);
            obj.SMFSFlag.Selected=zeros(1,obj.NCurves);
            % SMFS Downflags
            obj.SMFSFlagDown.FitSensitivity=zeros(1,obj.NCurves);
            obj.SMFSFlagDown.FitSinoidal=zeros(1,obj.NCurves);
            obj.SMFSFlagDown.FitLinear=zeros(1,obj.NCurves);
            obj.SMFSFlagDown.FitLinearRet=zeros(1,obj.NCurves);
            obj.SMFSFlagDown.Fit=zeros(1,obj.NCurves);
            obj.SMFSFlagDown.Uncorrupt=ones(1,obj.NCurves);
            obj.SMFSFlagDown.AppMinCrit=zeros(1,obj.NCurves);
            obj.SMFSFlagDown.RetMinCrit=zeros(1,obj.NCurves);
            obj.SMFSFlagDown.LengthRequisite=zeros(1,obj.NCurves);
            obj.SMFSFlagDown.SnapIn=zeros(1,obj.NCurves);
            obj.SMFSFlagDown.PullingLength=zeros(1,obj.NCurves);
            obj.SMFSFlagDown.Selected=zeros(1,obj.NCurves);
            % Debugging
            obj.DebugFlag.Plot=zeros(1,obj.NCurves);
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
                for i=1:obj.NCurves
                    [obj.App{i},obj.HHApp{i}] = obj.get_force_curve_data(i,0,0,0);
                    [obj.Ret{i},obj.HHRet{i}] = obj.get_force_curve_data(i,1,0,0);
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
        
    end
    
    methods
        % methods for visualization, plotting, statistics and quality control
        
        
        function show_force_curve(obj,ZoomMult,k,fig)
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
                [AppY,AppX] = obj.get_force_curve_data(k,0,1,0);
                [RetY,RetX] = obj.get_force_curve_data(k,1,1,0);
            catch
                [AppY,AppX] = obj.get_force_curve_data(k,0,0,0);
                [RetY,RetX] = obj.get_force_curve_data(k,1,0,0);
            end
            
            
            subplot(2,1,1)
            title(sprintf('Curve Nr.%i of %s',k,obj.Name))
            hold on
            [MultiplierX,UnitX,~] = AFMImage.parse_unit_scale(range(RetX),'m',10);
            [MultiplierY,UnitY,~] = AFMImage.parse_unit_scale(range(RetY),'N',5);
            plot(AppX*MultiplierX,AppY*MultiplierY,RetX*MultiplierX,RetY*MultiplierY,'LineWidth',1.5);
            Legends = {'Approach','Retract'};
            
            if obj.CPFlag.HertzFitted == 1
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
        
        function Fig = show_analyzed_fibril(obj)
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
            title('E-Modulus Map')
        end
        
        function show_height_map(obj)
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
        
        function show_e_mod_map(obj)
            Title = sprintf('E-Mod Map of %s',obj.Name);
            figure('Name',Title);
            subplot(2,2,1)
            Lower = nanmean(obj.EModMapHertz,'all')-1.5*nanstd(obj.EModMapHertz,0,'all');
            Upper = nanmean(obj.EModMapHertz,'all')+1.5*nanstd(obj.EModMapHertz,0,'all');
            I = imresize(obj.EModMapHertz,[1024 1024]);
            imshow(I,[Lower Upper],'Colormap',copper);
            colorbar
            title('EMod - Hertz-Sneddon Method')
            subplot(2,2,2)
            surf(obj.EModMapHertz(:,:,1),'LineStyle','none','FaceLighting','gouraud','FaceColor','interp')
            light('Style','local')
            subplot(2,2,3)
            Lower = nanmean(obj.EModMapOliverPharr,'all')-1.5*nanstd(obj.EModMapOliverPharr,0,'all');
            Upper = nanmean(obj.EModMapOliverPharr,'all')+1.5*nanstd(obj.EModMapOliverPharr,0,'all');
            I = imresize(obj.EModMapOliverPharr,[1024 1024]);
            imshow(I,[Lower Upper]);
            colorbar
            title('Emod Oliver-Pharr Method')
            subplot(2,2,4)
            surf(obj.EModMapOliverPharr,'LineStyle','none','FaceLighting','gouraud','FaceColor','interp')
            light('Style','local')
        end
        
        function quality_control_oliver_pharr_fibril(obj,PauseTime)
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
                
                [AppY,AppX] = obj.get_force_curve_data(k,0,1,0);
                [RetY,RetX] = obj.get_force_curve_data(k,1,1,0);
                
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
                        nanmean(obj.EModOliverPharr(obj.RectApexIndex))*1e-6,...
                        nanmedian(obj.EModOliverPharr(obj.RectApexIndex))*1e-6,...
                        nanstd(obj.EModOliverPharr(obj.RectApexIndex))*1e-6));
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
        
        function quality_control_hertz_sneddon_fibril(obj,PauseTime)
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
                
                [AppY,AppX] = obj.get_force_curve_data(k,0,1,1);
                [RetY,RetX] = obj.get_force_curve_data(k,1,1,1);
                
                subplot(2,2,2)
                
                % Determine X-Range for HertzModel
                X = AppX - obj.CP(k,1);
                X(X<0) = [];
                HertzModelX = 0:range(X)/100:2*max(X);
                HertzModelY = feval(obj.HertzFit{k},HertzModelX);
                
                
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
                        nanmean(obj.EModHertz(obj.RectApexIndex))*1e-6,...
                        nanmedian(obj.EModHertz(obj.RectApexIndex))*1e-6,...
                        nanstd(obj.EModHertz(obj.RectApexIndex))*1e-6));
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
        
        function quality_control_oliver_pharr(obj,PauseTime)
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
                
                [App,HHApp] = obj.get_force_curve_data(m,0,1,0);
                [Ret,HHRet] = obj.get_force_curve_data(m,1,1,0);
                
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
                        nanmean(obj.EModOliverPharr)*1e-6,...
                        nanmedian(obj.EModOliverPharr)*1e-6,...
                        nanstd(obj.EModOliverPharr)*1e-6));
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
        
        function quality_control_hertz_sneddon(obj,PauseTime)
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
                
                [AppY,AppX] = obj.get_force_curve_data(m,0,1,1);
                [RetY,RetX] = obj.get_force_curve_data(m,1,1,1);
                
                % Determine X-Range for HertzModel
                X = AppX - obj.CP(m,1);
                X(X<0) = [];
                HertzModelX = 0:range(X)/100:2*max(X);
                FitModel = obj.HertzFit{m};
                FitModel.b = 0;
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
                        nanmean(obj.EModHertz)*1e-6,...
                        nanmedian(obj.EModHertz)*1e-6,...
                        nanstd(obj.EModHertz)*1e-6));
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
        
        function compare_hertz_oliver_fibril(obj)
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
        
    end
end