classdef ForceMap < matlab.mixin.Copyable
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
        
        Name            % name of the force map. taken as the name of the folder, containing the .csv files
        Date            % date when the force map was detected
        Time            % time when the force map was detected
        ID              % Identifier for relation to Experiment
        FileVersion     % Version of jpk-force-map file
        FileType        % File Type: force map or qi-map
        Folder          % location of the .csv files of the force map
        HostOS          % Operating System
        HostName        % Name of hosting system
        NCurves         % number of curves on the force map
        NumProfiles     % number of scanned profiles along the YSize of the force map
        NumPoints       % number of scanned points per profile along the XSize of the force map
        NumSegments     % number of segments
        SeriesTime      % ongoing time from start of experiment to end
        SegDuration     % duration of one segment
        SegFrequency    % frequency of one segment
        SegNumPoints    % number of points in every segment
        SecPerPoint     % seconds per point
        TStart          % starting time for time vector
        TEnd            % ending time for time vector
        SegTime         % time vector
        InterpTimeF      % time vector for force data interpolation
        InterpTimeH      % time vector for height data interpolation
        MaxPointsPerCurve
        XSize           % Size of imaged window in X-direction
        YSize           % Size of imaged window in Y-direction
        Velocity        % Approach and retraction velocity as defined in the force map settings
        GridAngle       % in degrees (°)
        Sensitivity
        SpringConstant
        DBanding        % Fourieranalysis-based estimate of DBanding perdiod (only available with sufficient resolution)
        RefSlope        % Refernce slope as determined from the upper curve slope from data from very hard
        % surfaces (mica,glass), either from glass parts beneath the specimen or from
        % separate reference force maps
        PixApp          % maximum number of measured points during approach
        PixRet          % maximum number of measured points during retraction
        SelectedCurves  % logical vector of length NCurves with 0s for excluded and 1s for included curves. gets initialized with ones
        TipRadius = 8  % (nominal, if not otherwise calculated) tip radius in nm for the chosen type of tip model
        PoissonR = 0.5  % standard Poisson ratio for most mechanical models
        Medium
        FibrilFlag
        DeltaPhi        % phase shift between indentation and force
        LossTangent     % Loss tangent 
        SineVarsF       % Variables of Force fitted Sine (Amplitude, Phaseshift)
        SineVarsH       % Variables of Height fitted Sine (..)
        FilterF         % filtered force data
        FilterH         % filtered height data
        FZShift         % original force data shifted to the zero line
        HZShift         % original height data shifted to the zero line
        psF             % phase shift of forve
        psH             % phase shift of indentation
        slopeF          % slope of linear force data
        slopeH          % slope of linear indentation data
        interceptF      % intercept of linear force data
        interceptH      % intercept of linear indentation data

    end
    properties
        % Curve data Properties
        
        App = {}        % approach force data in Newton
        Force = {}      % all force data in Newton
        Ret = {}        % retraction force data in Newton
        Height = {}     % all height data in meters
        Indentation = {}% all indentation data in meters
        HHApp = {}      % capacitive-sensor-height approach data in meters
        HHRet = {}      % capacitive-sensor-height retract data in meters
        HHType          % Type of Head Height.(Default is capacitiveSensorHeight; switches to measuredHeight, if default does'nt exist)
        THApp = {}      % vertical tip height approach data in meters
        THRet = {}      % vertical tip height retract data in meters
        BasedApp = {}   % approach force data with subtracted base line and tilt in Newton
        BasedRet = {}   % retraction force data with subtracted base line and tilt in Newton
        BaseAndTiltFlag
        BasedForce = {}
        THeight = {}
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
        EModMapMicro1 % """"
        EModMapMicro2 % """"
        List2Map        % An R->RxR ((k)->(i,j)) mapping of indices to switch between the two representations
        Map2List        % An RxR->R ((i,j)->(k))mapping of indices to switch between the two representations
        FibDiam = []    % Estimated fibril diameter
        FibDiamSTD      % Estimated fibril diameter std
        FibMask         % Logical mask marking the whole fibril
        BackgroundMask  % Logical mask marking the glass/mica/even-substrate background
        ExclMask        % Manually chosen areas that are to be excluded for calculations of fibril EMod, FibDiam, DBanding etc.
        Apex            % Value of highest pixel in each profile
        RectApex        % Value of rectified apex location in each profile
        ApexIndex       % Index of highest pixel in each profile (List indexing!)
        RectApexIndex   % Index of rectified apex location in each profile (List indexing!)
    end
    properties
        % Properties related to EMod calculations
        
        Basefit = {}    % fit model used for the baseline fit in the base_and_tilt method
        EModHertz       % List of reduced smaple E-Modulus based on a the Hertz-Sneddon model
        EModOliverPharr % List of reduced sample E-Modulus based on the Oliver-Pharr method
        EModMicro1      % List of sample E-Modulus based on the for microrheology adapted Oliver-Pharr method
        EModMicro2      % List of sample E-Modulus based on the for microrheology adapted Oliver-Pharr method
        HertzFit        % HertzFit model generated in the calculate_e_mod_hertz method
        IndDepthHertz
        DZslope
        DZAppslope
        Stiffness
        IndDepth
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
    properties
        % SMFS related 
        MinRet % Minimum value of the
        Substrate       % Used substrate for the measurement 
        EnvCond         % Environmental condition during the experiment
        ChipCant        % AFM-Chip number and Cantilever label
        Chipbox         % AFM-Chipbox number (in Roman numerals)
        ModDate         % Modified Date is a modification of the poperty Date. Dots are removed
        SMFSFlag        %
        BasedRetCorr    % BasedRet data corrected
    end
    
    methods
        % Main methods of the class
        
        function obj = ForceMap(MapFullFile,DataFolder,TempID,FakeOpt,NSynthCurves)
            %%% Constructor of the class
            
            % Specify the folder where the files live. And import them.
            % Also get curent folder and return to it after import of
            % files.
            % Assigns the properties that can be found in the jpk-file
            % already
            
            current = what();
            
            obj.ID = TempID;
            
            if nargin >= 4 && isequal(FakeOpt,'Dummy')
                obj.create_dummy_force_map(NSynthCurves);
            end
            
            % get OS and use appropriate fitting system command
            FullOS = computer;
            OS = FullOS(1:3);
            obj.HostOS = OS;
            if isequal(OS,'PCW')
                obj.HostName = getenv('COMPUTERNAME');
            elseif isequal(OS,'GLN')
                obj.HostName = getenv('HOSTNAME');
            elseif isequal(OS,'MAC')
                obj.HostName = getenv('HOSTNAME');
            end
            
            % Unpack jpk-force-map with 7zip call to the terminal
            TempFolder = obj.unpack_jpk_force_map(MapFullFile,DataFolder);
            
            Index = regexp(obj.ID,'(?<=\-).','all');
            LoadMessage = sprintf('loading data into ForceMap Nr.%s',obj.ID(Index(end):end));
            disp(LoadMessage)
            
            % reading header properties into object
            obj.read_in_header_properties(TempFolder);
            
            %loading curve data into cell arrays
            obj.load_force_curves(TempFolder);
            
            %clean up unzipped jpk-force-map file
            rmdir(TempFolder,'s');
            
            % intitialize masks
            obj.ExclMask = logical(ones(obj.NumProfiles,obj.NumPoints));
            obj.FibMask = logical(zeros(obj.NumProfiles,obj.NumPoints));
            
            obj.create_and_level_height_map();
            
            obj.SelectedCurves = ones(obj.NCurves,1);
            
            obj.initialize_flags();
            
            % Save ForceMap and then change back into original folder
            cd(current.path);
            current = what();
            cd(obj.Folder)
            savename = sprintf('%s.mat',obj.Name);
            save(savename,'obj','-v7.3')
            cd(current.path)
            disp('loading successfull. object saved in objects folder')
        end
        
        function choose_curves(obj)
            % Interactive dialogue for curve selection. Returns a logic 0/1
            % vector, that is later used and can also be changed by the
            % other force map methods.
            f = figure('Name','Curve Selection','Position',[10000 10000 1200 900]);
            movegui(f);
            for i=1:obj.NCurves
                dlgtitle = sprintf('Curve selection %i/%i',i,obj.NCurves);
                plottitle = sprintf('Curve Nr. %i',i);
                plot(obj.HHApp{i},obj.App{i},obj.HHRet{i},obj.Ret{i});
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
        
        function base_and_tilt(obj,RunMode)
            % subtract baseline and tilt from the forcecurve by fitting a function to
            % the non contact domain the function tries to fit a non-affine-linear
            % function. If the linear fit is too bad, the function tries a 9th grade
            % polynomial fit instead
            if nargin<2
                RunMode = 'poly9';
            end
            Range = find(obj.SelectedCurves);
            h = waitbar(0,'Setting up...');
            for i=Range'
                prog = i/obj.NCurves;
                waitbar(prog,h,'processing baseline fits...');
                [ncd , ncidx] = ForceMap.no_contact_domain(obj.App{i});
                gof.rsquare = 0;
                for j=1:1
                    [testfit,goftest] = fit(obj.HHApp{i}(1:ncidx),smooth(ncd),'poly1','Normalize','on');
                    if goftest.rsquare>gof.rsquare
                        gof = goftest;
                        obj.Basefit{i} = testfit;
                    end
                end
                
                if gof.rsquare < 0.9 && isequal(RunMode,'poly9')
                    %ncd(ncidx+1:length(obj.App{i})) = feval(obj.Basefit{i},obj.HHApp{i}(ncidx+1:length(obj.HHApp{i})));
                    ext_range = obj.HHApp{i};
                    ext_l = length(obj.HHApp{i});
                    rangefit = fit([1:length(obj.HHApp{i})]',obj.HHApp{i},'poly1','Normalize','on');
                    ext_range(ext_l+1:floor(1.01*ext_l)) = feval(rangefit,[ext_l+1:floor(1.01*ext_l)]');
                    ncd(ncidx+1:length(ext_range)) = feval(obj.Basefit{i},ext_range(ncidx+1:length(ext_range)));
                    [obj.Basefit{i},gof] =  fit(ext_range,smooth(ncd),'poly9','Normalize','on');
                end
                %                 if gof.rsquare < 0.5
                %                     f = figure('Name','Base and tilt','Position',[10000 10000 800 600]);
                %                     movegui(f);
                %                     plot(obj.HHApp{i},feval(obj.Basefit{i},obj.HHApp{i}),obj.HHApp{i},obj.App{i});
                %                     legend('Fit','Data');
                %                     plottitle = sprintf('Curve Nr. %i',i);
                %                     title(plottitle);
                %                     answ = questdlg('Quality of base line fit is very low. consider excluding this force curve',...
                %                         'Base line and tilt','Exclude it',...
                %                         'Keep current fit','Keep current fit');
                %                     if isequal(answ,'Exclude it')
                %                         obj.SelectedCurves(i) = 0;
                %                         continue
                %                     else
                %                         return
                %                     end
                %                     close(gcf)
                %                 end
                obj.BasedApp{i} = (obj.App{i}-feval(obj.Basefit{i},obj.HHApp{i}));
                obj.BasedRet{i} = (obj.Ret{i}-feval(obj.Basefit{i},obj.HHRet{i}));
                
                if ~isempty(obj.NumSegments)  && (obj.NumSegments > 2)
                    
                    %for i=1:obj.NCurves
                        for j=1:obj.NumSegments
                            obj.BasedForce{i,j} = (obj.Force{i,j}-feval(obj.Basefit{i},obj.Height{i,j}));
                        end
                    %end
                end
            end
            % calculate vertical tip position by subtracting vertical tip deflection from head height
            iRange = find(obj.SelectedCurves);
            for i=iRange'
                obj.THApp{i} = obj.HHApp{i} - obj.BasedApp{i}/obj.SpringConstant;
                obj.THRet{i} = obj.HHRet{i} - obj.BasedRet{i}/obj.SpringConstant;
                
                if ~isempty(obj.NumSegments)  && (obj.NumSegments > 2)
                    
                    for j=1:obj.NumSegments
                        obj.THeight{i,j} = obj.Height{i,j} - obj.BasedForce{i,j}./obj.SpringConstant;
                    end
                end
            end
            close(h);
            obj.BaseAndTiltFlag = true;
            %             current = what();
            %             cd(obj.Folder)
            %             savename = sprintf('%s.mat',obj.Name);
            %             save(savename,'obj')
            %             cd(current.path)
        end
        
        function estimate_cp_snap_in(obj)
            % Takes the first point that crosses the approach baseline
            % starting from the max indent and sets the CP as the
            % interpolated point between the two nearest points to the
            % baseline (above and below)
            
            h = waitbar(0,'Setting up...','Name',obj.Name);
            Range = find(obj.SelectedCurves);
            for i=Range'
                waitbar(i/obj.NCurves,h,'Finding contact point for Snap-In curve...');
                AboveZeroBool = zeros(length(obj.BasedApp{i}),1);
                AboveZeroBool(find(obj.BasedApp{i}>0)) = 1;
                k = 0;
                while AboveZeroBool(end-k)
                    k = k + 1;
                end
                if k<5
                    obj.CP(i,1) = obj.HHApp{i}(floor(.5*end));
                    obj.CP(i,2) = obj.BasedApp{i}(floor(.5*end));
                    obj.CP_SnapIn(i,:) = obj.CP(i,:);
                    obj.SelectedCurves(i) = 0;
                    continue
                end
                    
                AboveBase = [obj.HHApp{i}(end-(k-1)) obj.BasedApp{i}(end-(k-1))];
                BelowBase = [obj.HHApp{i}(end-k) obj.BasedApp{i}(end-k)];
                obj.CP(i,1) = mean([AboveBase(1) BelowBase(1)]);
                obj.CP(i,2) = 0;
                obj.CP_SnapIn(i,:) = obj.CP(i,:);
                clear AboveZeroBool
            end
            close(h)
            obj.CPFlag.SnapIn = true;
            
        end
        
        function estimate_cp_rov(obj,batchsize)
            % find contact point with the method of ratio of variances. The method
            % iterates though every point and builds the ratio of the variance of a
            % bunch of points before and after the current point. the point with the
            % biggest ratio is the returned contact point [Nuria Gavara, 2016]
            if nargin<2
                batchsize = 20;
            end
            Range = find(obj.SelectedCurves);
            h = waitbar(0,'Setting up...','Name',obj.Name);
            CP_RoV = zeros(obj.NCurves,2);
            obj.CP_RoV = CP_RoV;
            for i=Range'
                prog = i/obj.NCurves;
                waitbar(prog,h,'applying ratio of variances method...');
                obj.RoV{i} = zeros(length(obj.BasedApp{i}),1);
                SmoothedApp = smoothdata(obj.BasedApp{i});
                % loop through points and calculate the RoV
                for j=(batchsize+1):(length(obj.BasedApp{i})-batchsize)
                    obj.RoV{i}(j,1) = var(smoothdata(obj.BasedApp{i}((j+1):(j+batchsize))))/...
                        var(SmoothedApp((j-batchsize):(j-1)));
                end
                % normalize RoV-curve
                obj.RoV{i} = obj.RoV{i}/range(obj.RoV{i});
                minrov = min(obj.RoV{i}(batchsize+1:length(obj.RoV{i})-batchsize));
                obj.RoV{i}(obj.RoV{i}==0) = minrov;
                [~,CPidx] = max(obj.RoV{i});
                obj.CP_RoV(i,:) = [obj.HHApp{i}(CPidx) obj.BasedApp{i}(CPidx)];
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
            for i=Range'
                Based = obj.BasedApp{i};
                THApp = obj.HHApp{i} - obj.BasedApp{i}/obj.SpringConstant;
                smoothx = smoothdata(THApp);
                smoothy = smoothdata(Based);
                Rsquare = 2*ones(length(Based),1);
                E = ones(length(Based),1);
                obj.DeltaE = ones(length(Based),1);
                testrange = floor(0.5*length(Based)):(length(Based)-5);
                msg = sprintf('applying goodness of fit method on curve Nr.%i/%i',i,obj.NCurves);
                for j=testrange
                    prog = (j-testrange(1))/length(testrange);
                    waitbar(prog,h,msg);
                    [emod,gof] = obj.hertz_fit_gof(smoothx,smoothy,j,1,'parabolic');
                    Rsquare(j) = gof.rsquare;
                    E(j) = emod;
                end
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
                obj.CPComboCurve{i} = obj.RoV{i}.*obj.GoF{i};
                [~,CPidx] = max(obj.CPComboCurve{i});
                obj.CP_Combo(i,:) = [obj.HHApp{i}(CPidx) obj.BasedApp{i}(CPidx)];
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
            elseif isequal(lower(RunMode),'zoom')
                runmode = 2;
            elseif isequal(lower(RunMode),'zoomdropout')
                runmode = 3;
                if nargin < 3
                    NumPasses = 100; % if not specified in arguments, NumPasses defaults to 100
                end
            elseif isequal(lower(RunMode),'zoomsweep')
                runmode = 4;
                if nargin < 3
                    NumPasses = 20; % if not specified in arguments, NumPasses defaults to 20
                end
            end
            ImgSize = NeuralNet.Layers(1).InputSize;
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
                        obj.CP_CNN(i,1) = Ypredicted(k,1)*range(obj.HHApp{i})+min(obj.HHApp{i});
                        obj.CP_CNN(i,2) = Ypredicted(k,2)*range(obj.BasedApp{i})+min(obj.BasedApp{i});
                        obj.CP(i,1) = obj.CP_CNN(i,1);
                        obj.CP(i,2) = obj.CP_CNN(i,2);
                        k = k + 1;
                    end
                    obj.CPFlag.CNN = 1;
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
                    obj.CP_MonteCarlo = zeros(NumPasses,2,obj.NCurves);
                    for i=iRange'
                        obj.CP_MonteCarlo(:,1,i) = obj.YDropPred(:,1,k)*range(obj.HHApp{i})+min(obj.HHApp{i});
                        obj.CP_MonteCarlo(:,2,i) = obj.YDropPred(:,2,k)*range(obj.BasedApp{i})+min(obj.BasedApp{i});
                        obj.CP(i,1) = mean(obj.CP_MonteCarlo(:,1,i));
                        obj.CP(i,2) = mean(obj.CP_MonteCarlo(:,2,i));
                        obj.CP_Dropout(i,1) = obj.CP(i,1);
                        obj.CP_Dropout(i,2) = obj.CP(i,2);
                        obj.CP_MonteCarlo_STD(i) = norm([std(obj.CP_MonteCarlo(:,1,i)) std(obj.CP_MonteCarlo(:,2,i))]);
                        k = k + 1;
                    end
                    obj.CPFlag.Dropout = 1;
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
                        obj.CP_CNNZoom(i,1) = Ypredicted(k,1)*range(obj.HHApp{i})+min(obj.HHApp{i});
                        obj.CP_CNNZoom(i,2) = Ypredicted(k,2)*range(obj.BasedApp{i})+min(obj.BasedApp{i});
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
                        obj.CP_CNNZoom(i,1) = Ypredicted(k,1)*range(ZoomObj.HHApp{i})+min(ZoomObj.HHApp{i});
                        obj.CP_CNNZoom(i,2) = Ypredicted(k,2)*range(ZoomObj.BasedApp{i})+min(ZoomObj.BasedApp{i});
                        obj.CP(i,1) = obj.CP_CNNZoom(i,1);
                        obj.CP(i,2) = obj.CP_CNNZoom(i,2);
                        k = k + 1;
                    end
                    
                    obj.CPFlag.CNNZoom = 1;
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
                        obj.CP_CNNZoom(i,1) = Ypredicted(k,1)*range(obj.HHApp{i})+min(obj.HHApp{i});
                        obj.CP_CNNZoom(i,2) = Ypredicted(k,2)*range(obj.BasedApp{i})+min(obj.BasedApp{i});
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
                        obj.CP_CNNZoom(i,1) = Ypredicted(k,1)*range(ZoomObj.HHApp{i})+min(ZoomObj.HHApp{i});
                        obj.CP_CNNZoom(i,2) = Ypredicted(k,2)*range(ZoomObj.BasedApp{i})+min(ZoomObj.BasedApp{i});
                        obj.CP(i,1) = obj.CP_CNNZoom(i,1);
                        obj.CP(i,2) = obj.CP_CNNZoom(i,2);
                        k = k + 1;
                    end
                    
                    obj.CPFlag.CNNZoomDropout = 1;
                case 4
                    % ZoomSweep
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
                        obj.CP_CNNZoom(i,1) = Ypredicted(k,1)*range(obj.HHApp{i})+min(obj.HHApp{i});
                        obj.CP_CNNZoom(i,2) = Ypredicted(k,2)*range(obj.BasedApp{i})+min(obj.BasedApp{i});
                        obj.CP(i,1) = obj.CP_CNNZoom(i,1);
                        obj.CP(i,2) = obj.CP_CNNZoom(i,2);
                        k = k + 1;
                    end
                    
                    waitbar(2/3,h,'Predicting zoomed CP, sweeping over multiple zooms');
                    MaxZoom = 0.7;
                    ZoomFactor = (1-MaxZoom):MaxZoom/(NumPasses-1):1;
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
                            TempCP(i,1,j) = Ypredicted(k+length(iRange)*(j-1),1)*range(ZoomCell{j}.HHApp{i})+min(ZoomCell{j}.HHApp{i});
                            TempCP(i,2,j) = Ypredicted(k+length(iRange)*(j-1),2)*range(ZoomCell{j}.BasedApp{i})+min(ZoomCell{j}.BasedApp{i});
                        end
                        obj.CP_CNNZoomSweep(i,1) = mean(TempCP(i,1,:),3);
                        obj.CP_CNNZoomSweep(i,2) = mean(TempCP(i,2,:),3);
                        obj.CP(i,1) = obj.CP_CNNZoomSweep(i,1);
                        obj.CP(i,2) = obj.CP_CNNZoomSweep(i,2);
                        k = k + 1;
                    end
                    
                    obj.CPFlag.CNNZoomSweep = 1;
            end
            close(h)
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
                load = zeros(length(obj.BasedApp{i}),2);
                unload = zeros(length(obj.Ret{i}),2);
                for j=1:length(obj.BasedApp{i})
                    load(end-(j-1),1) = obj.HHApp{i}(j);
                    load(j,2) = obj.BasedApp{i}(j);
                end
                for j=1:length(obj.Ret{i})
                    unload(end-(j-1),1) = obj.HHRet{i}(j);
                    unload(j,2) = obj.Ret{i}(j);
                end
                [obj.LoadOld{i},obj.UnloadOld{i},Position,vDef] = ContactPoint_sort(load,unload);
                obj.CP(i,2) = obj.BasedApp{i}(Position);
                obj.CP(i,1) = obj.HHApp{i}(Position);
                obj.CP_Old(i,1) =obj.CP(i,1);
                obj.CP_Old(i,2) =obj.CP(i,2);
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
                plot(obj.HHApp{j},obj.BasedApp{j});
                plottitle = sprintf('Curve Nr.%i/%i\n Click or click and drag the point to the contact point\n Click the red area to exclude the current curve from further analysis\n Click the green area to go back one curve',j,obj.NCurves);
                title(plottitle);
                [~, domainidx] = ForceMap.no_contact_domain(obj.App{j});
                axis([obj.HHApp{j}(floor(domainidx*Zoom)) inf -inf inf])
                XRange = range(obj.HHApp{j}(floor(domainidx*Zoom):end));
                YRange = range(obj.BasedApp{j});
                BtnSemAxisX = XRange/8;
                BtnSemAxisY = YRange/8;
                XPosDel = max(obj.HHApp{j}) - 5/8*XRange;
                YPosDel = max(obj.BasedApp{j}) - 1/3*YRange;
                XPosBack = max(obj.HHApp{j}) - 7/8*XRange;
                YPosBack = max(obj.BasedApp{j}) - 1/3*YRange;
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
                obj.CP_HardSurface(i,1) = obj.HHApp{i}(end) - obj.BasedApp{i}(end)/obj.SpringConstant;
                obj.CP_HardSurface(i,2) = 0;
                %% Debugging
                % plot(obj.HHApp{i},obj.BasedApp{i});
                % drawpoint('Position',[obj.CP_HardSurface(i,1) obj.CP_HardSurface(i,2)]);
            end
            obj.CPFlag.HardSurface = 1;
        end
        
        function [E,HertzFit] = calculate_e_mod_hertz(obj,CPType,TipShape,curve_percent,AllowXShift)
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
            end
            iRange = find(obj.SelectedCurves);
            obj.EModHertz = zeros(obj.NCurves,1);
            obj.IndDepthHertz = zeros(obj.NCurves,1);
            for i=iRange'
                if isequal(lower(CPType),'cnn')
                    CP = obj.CP(i,:);
                elseif isequal(lower(CPType),'old')
                    CP = obj.CP_Old(i,:);
                elseif isequal(lower(CPType),'rov')
                    CP = obj.CP_RoV(i,:);
                elseif isequal(lower(CPType),'gof')
                    CP = obj.CP_GoF(i,:);
                elseif isequal(lower(CPType),'combo')
                    CP = obj.CP_Combo(i,:);
                elseif isequal(lower(CPType),'manual')
                    CP = obj.Man_CP(i,:);
                elseif isequal(lower(CPType),'snap-in')
                    CP = obj.CP_SnapIn(i,:);
                else
                    CP = obj.CP(i,:);
                end
                force = obj.BasedApp{i} - CP(2);
                tip_h = (obj.HHApp{i} - CP(1)) - force/obj.SpringConstant;
                tip_h(tip_h < 0) = [];
                Max = max(tip_h);
                obj.IndDepthHertz(i) = Max(1);
                % delete everything below curve_percent of the maximum
                % force
                force(1:(length(force)-length(tip_h))) = [];
                force(force<(1-curve_percent)*max(force)) = [];
                tip_h(1:(length(tip_h)-length(force))) = [];
                RangeF = range(force);
                RangeTH = range(tip_h);
                force = force/RangeF;
                tip_h = tip_h/RangeTH;
                if isequal(TipShape,'parabolic')
                    if AllowXShift
                        s = fitoptions('Method','NonlinearLeastSquares',...
                            'Lower',[10^(-5) 0],...
                            'Upper',[inf inf],...
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
                        Hertzfit = fit(tip_h,...
                            force,f);
                        % calculate E module based on the Hertz model. Be careful
                        % to convert to unnormalized data again
                        if isempty(obj.FibDiam)
                            R_eff = obj.TipRadius*1e-9;
                        else
                            R_eff = 1/(1/(obj.TipRadius*1e-9) + 1/(obj.FibDiam/2));
                        end
                        EMod = 3*(Hertzfit.a*RangeF/RangeTH^(3/2))/(4*sqrt(R_eff))*(1-obj.PoissonR^2);
                    catch
                        EMod = nan;
                        Hertzfit.a = 0;
                        if AllowXShift
                            Hertzfit.b = 0;
                        end
                    end
                elseif isequal(shape,'spherical')
                elseif isequal(shape,'conical')
                elseif isequal(shape,'pyramid')
                end
                obj.EModHertz(i) = EMod;
                % Convert the model to the right scale so it can be plotted
                % correctly later
                warning('off','all');
                Hertzfit.a = Hertzfit.a*RangeF/RangeTH^(3/2);
                if AllowXShift
                Hertzfit.b = Hertzfit.b*RangeTH;
                obj.CP_HertzFitted(i,1) = CP(1)-Hertzfit.b;
                obj.CP_HertzFitted(i,2) = CP(2);
                obj.CP(i,:) = obj.CP_HertzFitted(i,:);
                % Not sure about this one
                % obj.IndDepthHertz(i) = obj.IndDepthHertz(i) + Hertzfit.b;
                end
                warning('on','all');
                obj.HertzFit{i} = Hertzfit;
                
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
            for i=1:obj.NumProfiles
                for j=1:obj.NumPoints
                    obj.EModMapHertz(i,j,1) = obj.EModHertz(obj.Map2List(i,j));
                end
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
            obj.DZAppslope = zeros(obj.NCurves,1);
            obj.Stiffness = zeros(obj.NCurves,1);
            obj.IndDepth = zeros(obj.NCurves,1);
            obj.IndentArea = zeros(obj.NCurves,1);
            for i=Range'
                Z = obj.HHRet{i} - obj.CP(i,1);
                D = (obj.BasedRet{i} - obj.CP(i,2))/obj.SpringConstant;
                Zmax(i) = max(Z);
                Dmax(i) = max(D);
                DCurvePercent = D(D>=(1-CurvePercent)*Dmax(i));
                ZCurvePercent = Z(1:length(DCurvePercent));
                LineFit = polyfit(ZCurvePercent,DCurvePercent,1);
                obj.DZslope(i) = LineFit(1);
                
                ZApp = obj.HHApp{i} - obj.CP(i,1);
                DApp = (obj.BasedApp{i} - obj.CP(i,2))/obj.SpringConstant;
                ZAppmax(i) = max(ZApp);
                DAppmax(i) = max(DApp);
                DAppCurvePercent = DApp(DApp>=(1-CurvePercent)*DAppmax(i));
                ZAppCurvePercent = ZApp(1:length(DAppCurvePercent));
                LineFitApp = polyfit(ZAppCurvePercent,DAppCurvePercent,1);
                obj.DZAppslope(i) = LineFitApp(1);
                    
                Hmax(i) = Zmax(i) - Dmax(i);
                dD = Dmax(i) - (1-CurvePercent)*Dmax(i);
                dh = dD*(1./obj.DZslope(i) - 1/(obj.RefSlope));
                df = dD*obj.SpringConstant;
                obj.Stiffness(i) = df/dh;
                Fmax(i) = Dmax(i).*obj.SpringConstant;
                obj.IndDepth(i) = Hmax(i) - Epsilon.*Fmax(i)/obj.Stiffness(i);
                % IndentArea is taken as the linear interpolation between
                % the two numeric values of TipProjArea the Hc(i) falls
                % inbetween
                try
                    obj.IndentArea(i) = ((1-(obj.IndDepth(i)*1e9-floor(obj.IndDepth(i)*1e9)))*TipProjArea(floor(obj.IndDepth(i)*1e9))...
                        + (obj.IndDepth(i)*1e9-floor(obj.IndDepth(i)*1e9))*TipProjArea(ceil(obj.IndDepth(i)*1e9)));
                      
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
            
            for i=1:obj.NumProfiles
                for j=1:obj.NumPoints
                    obj.EModMapOliverPharr(i,j,1) = obj.EModOliverPharr(obj.Map2List(i,j));
                end
            end
        end
        
        function [EModMicro1, EModMicro2] = calculate_e_mod_microrheology(obj,TipProjArea,CurvePercent)
            
            if nargin < 3
                CurvePercent = 0.75;
            end
            
            Range = find(obj.SelectedCurves);
            Epsilon = 0.73; % Correction constant from Oliver Pharr Method (1992)
            Beta = 1.0226; %Correction constant from Oliver Pharr Method (1992)
            obj.EModMicro1 = zeros(obj.NCurves,1);
            obj.EModMicro2 = zeros(obj.NCurves,1);
            EModMicro1 = zeros(obj.NCurves,1);
            EModMicro2 = zeros(obj.NCurves,1);
            obj.ProjTipArea = TipProjArea;
            obj.DZslope = zeros(obj.NCurves,1);
            obj.DZAppslope = zeros(obj.NCurves,1);
            obj.Stiffness = zeros(obj.NCurves,1);
            obj.IndDepth = zeros(obj.NCurves,1);
            obj.IndentArea = zeros(obj.NCurves,1);
            for i=Range'
                for j=1:obj.NumSegments
                    Z = obj.HHRet{i} - obj.CP(i,1);
                    D = (obj.BasedRet{i} - obj.CP(i,2))/obj.SpringConstant;
                    Zmax(i) = max(Z);
                    Dmax(i) = max(D);
                    DCurvePercent = D(D>=(1-CurvePercent)*Dmax(i));
                    ZCurvePercent = Z(1:length(DCurvePercent));
                    LineFit = polyfit(ZCurvePercent,DCurvePercent,1);
                    obj.DZslope(i) = LineFit(1);
                    
                    ZApp = obj.HHApp{i} - obj.CP(i,1);
                    DApp = (obj.BasedApp{i} - obj.CP(i,2))/obj.SpringConstant;
                    ZAppmax(i) = max(ZApp);
                    DAppmax(i) = max(DApp);
                    DAppCurvePercent = DApp(DApp>=(1-CurvePercent)*DAppmax(i));
                    ZAppCurvePercent = ZApp(1:length(DAppCurvePercent));
                    LineFitApp = polyfit(ZAppCurvePercent,DAppCurvePercent,1);
                    obj.DZAppslope(i) = LineFitApp(1);
                    
                    Hmax(i) = Zmax(i) - Dmax(i);
                    dD = Dmax(i) - (1-CurvePercent)*Dmax(i);
                    dh = dD*(1./obj.DZslope(i) - 1/(obj.RefSlope));
                    df = dD*obj.SpringConstant;
                    obj.Stiffness(i) = df/dh;
                    Fmax(i) = Dmax(i).*obj.SpringConstant;
                    obj.IndDepth(i) = Hmax(i) - Epsilon.*Fmax(i)/obj.Stiffness(i);
                    % IndentArea is taken as the linear interpolation between
                    % the two numeric values of TipProjArea the Hc(i) falls
                    % inbetween
                    
                    if obj.SegFrequency{j} > 0
                        
                        x = obj.SegTime{j};
                        
                        %Y-values fitted sine of indentation and force:
                        %obj.SineFunctionF = obj.SineVarsF{i,j}(1)*(sin((2*pi*x)/obj.SineVarsF{i,j}(2) + 2*pi/obj.SineVarsF{i,j}(3)));
                        %obj.SineFunctionH = obj.SineVarsH{i,j}(1)*(sin((2*pi*x)/obj.SineVarsH{i,j}(2) + 2*pi/obj.SineVarsH{i,j}(3)));


                        try
                             obj.IndentArea(i) = ((1-(obj.IndDepth(i)*1e9-floor(obj.IndDepth(i)*1e9)))*TipProjArea(floor(obj.IndDepth(i)*1e9))...
                                + (obj.IndDepth(i)*1e9-floor(obj.IndDepth(i)*1e9))*TipProjArea(ceil(obj.IndDepth(i)*1e9)));

                            EModMicro1(i,j) = sqrt(pi/obj.IndentArea(i))*1/2*...
                                (1-obj.PoissonR^2)*cosd(obj.DeltaPhi{i,j})*...
                                (obj.SineVarsF{i,j}(1)/obj.SineVarsH{i,j}(1));

                            EModMicro2(i,j) = sqrt(pi/obj.IndentArea(i))*1/2*...
                                (1-obj.PoissonR^2)*sind(obj.DeltaPhi{i,j})*...
                                (obj.SineVarsF{i,j}(1)/obj.SineVarsH{i,j}(1));

                            if EModMicro1(i,j) <= 0
                                EModMicro1(i,j) = NaN;
                            end

                            if EModMicro2(i,j) <= 0
                                EModMicro2(i,j) = NaN;
                            end

                        catch
                            EModMicro1(i,j) = NaN;
                            EModMicro2(i,j) = NaN;
                        end
                    end
                end
            end
            
            obj.EModMicro1 = EModMicro1;
            obj.EModMicro2 = EModMicro2;
            % Write values into EModMapMicrorheology
            
            for i=1:obj.NumProfiles
                for j=1:obj.NumPoints
                    obj.EModMapMicro1(i,j,1) = obj.EModMicro1(obj.Map2List(i,j));
                    obj.EModMapMicro2(i,j,1) = obj.EModMicro2(obj.Map2List(i,j));
                end
            end
        end
        
        function manual_exclusion(obj)
            
            obj.ExclMask = logical(ones(obj.NumProfiles,obj.NumPoints));
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
            [obj.Apex,obj.ApexIndex] = max(obj.HeightMap(:,:,1).*obj.FibMask,[],2);
            
            if obj.FibrilFlag.Straight == 1
                obj.RectApex = zeros(obj.NumProfiles,1);
                obj.RectApexIndex = zeros(obj.NumProfiles,1);
                obj.RectApexIndex = round(predictGP_mean([1:obj.NumProfiles],[1:obj.NumProfiles],1,5*obj.NumProfiles,obj.ApexIndex,1));
                for i=1:obj.NumProfiles
                    obj.RectApex(i) = obj.HeightMap(i,obj.RectApexIndex(i),1);
                end
            else
                obj.RectApex = obj.Apex;
                obj.RectApexIndex = obj.ApexIndex;
            end
            
            k = 1;
            for i=1:obj.NumProfiles
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
            for i=1:obj.NumProfiles
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
                prog = i/obj.NCurves;
                waitbar(prog,h,'Converting force curves to images...');
                norm =  round(normalize(obj.BasedApp{i},'range',[1 imres]));
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
            for i=1:obj.NumProfiles
                for j=1:obj.NumPoints
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

        function [MinApp] = min_force(obj)           
            for ii=1:obj.NCurves
            MinApp(ii)=min(obj.BasedApp{ii});
            MinRet(ii)=min(obj.BasedRet{ii}); 
            end
            obj.MinRet = MinRet;
        end
        
        function fc_based_ret_correction(obj,DataShareStart,DataShareEnd)  
            % fc_based_ret_correction: A function to correct for an AFM
            % based baseline deviation between the approach and retraction
            % data
        if nargin <2
            DataShareStart=0.05; % 5%
            DataShareEnd=0.1; % 1%
        end
        % loop over all force curves  
        for kk=1:100
            DataPts=size(obj.BasedApp{kk}); % Determine the quantity of data points in the force curve 
            LimitIdx1=round(DataPts(1)*DataShareStart); % Determine the corresponidng index
            LimitIdx2=round(DataPts(1)*DataShareEnd);
            CorrMean=mean(abs(obj.BasedApp{kk}(LimitIdx1:LimitIdx2,1))-abs(obj.BasedRet{kk}(DataPts(1)-LimitIdx2:DataPts(1)-LimitIdx1,1))); % Calculate the mean of the difference data
            obj.BasedRetCorr{kk}=obj.BasedRet{kk}-CorrMean; % Correct the BasedRet data with the mean of the correction data
        end        
        % %% Appendix
        % close all
        % % Define variables
        % kk=1
        % x100=-100e-9; % Defines 100nm
        % x500=-500e-9; % Defines 500nm
        % % Graphical preview
        % fig=gcf;
        % fig.Units='normalized'; % changes to normalized unit settings, necessary to receive the full screen size in the next line
        % fig.Color='white'; % changes the background color of the figure
        % fig.OuterPosition=[0.5 0 0.5 1];% changes the size of the figure to half screen
        % fig.PaperOrientation='landscape';
        % grid on
        % hold on
        % % "Origin" data
        % plot(a.FM{1}.THApp{kk}-a.FM{1}.CP_HardSurface(kk,1),a.FM{1}.BasedApp{kk},'b');
        % plot(a.FM{1}.THRet{kk}-a.FM{1}.CP_HardSurface(kk,1),a.FM{1}.BasedRet{kk},'r');
        % % Retention data corrected
        % plot(a.FM{1}.THRet{kk}-a.FM{1}.CP_HardSurface(kk,1),a.FM{ii}.BasedRetCorr{kk},'g');
        % % DataShare part of the data
        % plot(a.FM{1}.THApp{kk}(LimitIdx1:LimitIdx2,1)-a.FM{1}.CP_HardSurface(kk,1),a.FM{1}.BasedApp{kk}(LimitIdx1:LimitIdx2,1),'y');
        % plot(a.FM{1}.THRet{kk}(DataPts(1)-LimitIdx2:DataPts(1)-LimitIdx1,1)-a.FM{1}.CP_HardSurface(kk,1),a.FM{1}.BasedRet{kk}(DataPts(1)-LimitIdx2:DataPts(1)-LimitIdx1,1),'m');
        % % Markers
        % plot(a.FM{1}.THRet{kk}(ThreshIdx,1)-a.FM{1}.CP_HardSurface(kk,1),a.FM{1}.BasedRet{kk}(ThreshIdx,1),'kx','MarkerSize',20);
        % plot(a.FM{1}.THApp{kk}(ThreshIdx,1)-a.FM{1}.CP_HardSurface(kk,1),a.FM{1}.BasedApp{kk}(ThreshIdx,1),'mx','MarkerSize',20);
        % Distances lines
        % line([x100 x100], ylim,'Color','k'); % Draws a vertical line
        % line([x500 x500], ylim,'Color','k'); % Draws a vertical line
        end
        
        function fc_chipprop(obj)
                 
                % Chip number and Cantilever
                exp15='(?!10)\d+\w{1}'; % Finds the chip number and the cantilever   
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
                % Environmental conditions
                if isempty(ext31)==0
                    obj.EnvCond='Water'; % Milli-Q water
                elseif isempty(ext32)==0
                    obj.EnvCond='HAc'; % Acetic acid
                else
                    obj.EnvCond='PBS'; % Phosphate buffered saline
                end
        end
        
        function fc_print(obj) % fc ... force curve
            % fc_print: A function to simply plot all force curves of a
            % force map without any selection taking place
            
            % Remove dots in obj.Date
            obj.ModDate=strrep(obj.Date,'.','');
            
            % Define remainder situation
            Remainder=mod(obj.NCurves,25);
            NFigures=floor(obj.NCurves./25);
            if Remainder ~= 0
                NFigures=NFigures+1;
            end    
            % Define variables for the figure name
            VelocityConvert=num2str(obj.Velocity*1e+9); % Convert into nm
            %figname=strcat(obj.ID,{'-'},obj.ModDate,{'-'},obj.Velocity,{'-'},obj.Name);
            figname=strcat(obj.ID,{'-'},obj.ModDate,{'-'},VelocityConvert,{'-'},obj.Substrate,{'-'},obj.EnvCond,{'-'},obj.Chipbox,{'-'},obj.ChipCant);
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
                    x50=-50e-9; % Defines 100nm
                    x500=-500e-9; % Defines 500nm
                    % Plot tile
                    nexttile
                    hold on
                    grid on
                    plot(obj.THApp{kk}-obj.CP_HardSurface(kk,1),obj.BasedApp{kk});
                    plot(obj.THRet{kk}-obj.CP_HardSurface(kk,1),obj.BasedRetCorr{kk});
                    line([x50 x50], ylim,'Color','k'); % Draws a vertical line                  
                    line([x500 x500], ylim,'Color','k'); % Draws a vertical line
                    % Title for each Subplot
                    ti=title(sprintf('%i',kk),'Color','k');                                     
                    ti.Units='normalized'; % Set units to 'normalized'  
                    ti.Position=[0.5,0.9]; % Position the subplot title within the subplot
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
        
        function fc_selection(obj) % fc ... force curve
            
            % Define remainder situation
            Remainder=mod(obj.NCurves,25);
            NFigures=floor(obj.NCurves./25);
            if Remainder ~= 0
                NFigures=NFigures+1;
            end    
            %% Figure loop          
            %figname=strcat(obj.ID,{'-'},obj.ModDate,{'-'},obj.Velocity,{'-'},obj.Name);
            figname=strcat(obj.ID,{'-'},obj.ModDate,{'-'},obj.Velocity,{'-'},obj.Substrate,{'-'},obj.EnvCond,{'-'},obj.Chipbox,{'-'},obj.ChipCant);
            figname=char(figname);
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
                    x100=-100e-9; % Defines 100nm
                    x500=-500e-9; % Defines 100nm
                    % Plot tile
                    nexttile
                    hold on
                    grid on
                    plot(obj.THApp{kk}-obj.CP_HardSurface(kk,1),obj.BasedApp{kk});
                    plot(obj.THRet{kk}-obj.CP_HardSurface(kk,1),obj.BasedRet{kk});
                    line([x100 x100], ylim,'Color','k'); % Draws a vertical line                  
                    line([x500 x500], ylim,'Color','k'); % Draws a vertical line
                    % Title for each Subplot
                    if obj.SelectedCurves(kk) == 0
                        ti=title(sprintf('%i',kk),'Color','r');
                    elseif obj.SelectedCurves(kk) == 1
                        ti=title(sprintf('%i',kk),'Color','b');
                    end
                    ti.Units='normalized'; % Set units to 'normalized'  
                    ti.Position=[0.5,0.9]; % Position the subplot title within the subplot
                    % Legend, x- and y-labels
                    %legend('Approach','Retraction','Location','best')
                    %xlabel('Tip-sample seperation  (nm)','FontSize',11,'Interpreter','latex');
                    %ylabel('Force (nN)','FontSize',11,'Interpreter','latex');                  
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
                    obj.SelectedCurves(kk-24:kk) = 1;
                end

                % Case 2: Select none
                if SelectBttns == 2
                    obj.SelectedCurves(kk-24:kk) = 0;
                end

                % Case 3: Select all - except of
                if SelectBttns == 3
                    obj.SelectedCurves(kk-24:kk) = 1;
                    prompt = {'Enter the force curve number you do not want to keep for analysis (For multiple selections just use the space key to separeat entries)'};
                    definput = {''};
                    opts.Interpreter = 'tex';
                    IndSelec=inputdlg(prompt,'Select all except of ...',[1 150],definput,opts); % Stores the individual selected fc as a cell array of character vectors 
                    IndSelec=str2num(IndSelec{1}); % Convert the cell array to numerals
                    obj.SelectedCurves(IndSelec) = 0;
                end

                if obj.SelectedCurves(kk-24:kk) == 0
                    t=title(sprintf('%i',kk),'Color','r');
                    elseif obj.SelectedCurves(kk-24:kk) == 1
                    t=title(sprintf('%i',kk),'Color','b');
                end
        
                % Case 4: Select none - except of
                if SelectBttns == 4
                    obj.SelectedCurves(kk-24:kk) = 0;
                    prompt = {'Enter the force curve number you want want to keep for analysis (For multiple selections just use the space key to separeat entries)'};
                    definput = {''};
                    opts.Interpreter = 'tex';
                    IndSelec=inputdlg(prompt,'Select all except of ...',[1 150],definput,opts); % Stores the individual selected fc as a cell array of character vectors 
                    IndSelec=str2num(IndSelec{1}); % Convert the cell array to numerals
                    obj.SelectedCurves(IndSelec) = 1;
                end
            end
                % Housekeeping
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
                    nexttile
                    hold on
                    grid on
                    plot(obj.THApp{kk}-obj.CP_HardSurface(kk,1),obj.BasedApp{kk});
                    plot(obj.THRet{kk}-obj.CP_HardSurface(kk,1),obj.BasedRet{kk});                    
                    % Title for each Subplot
                    if obj.SelectedCurves(kk) == 0
                        ti=title(sprintf('%i',kk),'Color','r');
                    elseif obj.SelectedCurves(kk) == 1
                        ti=title(sprintf('%i',kk),'Color','b');
                    end
                    ti.Units='normalized'; % Set units to 'normalized'  
                    ti.Position=[0.5,0.9]; % Position the subplot title within the subplot
                    % Legend, x- and y-labels
                    %legend('Approach','Retraction','Location','best')
                    %xlabel('Tip-sample seperation  (nm)','FontSize',11,'Interpreter','latex');
                    %ylabel('Force (nN)','FontSize',11,'Interpreter','latex');                  
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
        
        function fc_selection_procedure(obj,ThresholdDist,ThreshValue)
            % fc_selection_procedure: A function to distinguish between
            % force curves that fulfil or do not fulfil the logical
            % statement
            % SMFSFlag.Min = 0:  Force curves indicate a naked tip
            % SMFSFlag.Min = 1:  Force curves indicate a functionalized tip
            if nargin <2
                ThresholdDist=50e-9;  % 50 nm
                ThreshValue=50e-12;    % 50 pN
            elseif nargin<3
                ThreshValue=50e-12;    % 50 pN
            end
            % loop over all force curves
            for kk=1:100
            % Determine the index corresponding to the threshold distance
            ThreshDist=abs(obj.THRet{kk}-obj.CP_HardSurface(kk,1)+ThresholdDist);
            [~, ThreshIdx]=min(ThreshDist);
            % Check if the force curve is selected 
                if (obj.BasedApp{kk}(ThreshIdx)-obj.BasedRetCorr{kk}(ThreshIdx))>ThreshValue
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
        
        function sine_fit_for_microrheology(obj)

            
            obj.slopeF = zeros(obj.NCurves,obj.NumSegments);
            obj.slopeH = zeros(obj.NCurves,obj.NumSegments);
            obj.interceptF = zeros(obj.NCurves,obj.NumSegments);
            obj.interceptH = zeros(obj.NCurves,obj.NumSegments);
                 
            for i=1:obj.NCurves
                lastseg = obj.NumSegments-1;
                 
                for j=2:lastseg
                    
                    if obj.SegFrequency{j} > 0
                        

                        % Divide data through their range
                        rangeF = range(obj.BasedForce{i,j});
                        rangeH = range(obj.Indentation{i,j});
                        
                        obj.BasedForce{i,j} = obj.BasedForce{i,j}/rangeF;
                        obj.Indentation{i,j} = obj.Indentation{i,j}/rangeH;
                        
                        % Max values of Force and Indentation
                         %maxF1 = max(obj.BasedForce{i,j});
                         %maxH1 = max(obj.Indentation{i,j});

                         % Min values of Force and Indentation
                         %minF1 = min(obj.BasedForce{i,j});
                         %minH1 = min(obj.Indentation{i,j});

                         % Difference max min
                         %DiffF1 = maxF1 - minF1;
                         %DiffH1 = maxH1 - minH1;
                         

                         % Shift to Zero Line
                         %obj.FZShift{i,j} = obj.BasedForce{i,j}-maxF1+(DiffF1/2);
                         %obj.HZShift{i,j} = obj.Indentation{i,j}-maxH1+(DiffH1/2);
                         
                         

                         % Calculation of Sampling rate and Invariance to
                         % be able to subsequently choose the right filter
                         Samplingrate = 1/obj.SecPerPoint{j};
                         Invariance = Samplingrate/obj.SegFrequency{j};

                         % Force-Data are only placed horizontally in 
                         % the zero line with detrend
                         %ForceTrend{i,j} = detrend(obj.FZShift{i,j});
                         linearfitF = polyfit(obj.SegTime{j},obj.BasedForce{i,j},1);
                         Fvalueslinfit = polyval(linearfitF,obj.SegTime{j});
                         linearfitH = polyfit(obj.SegTime{j},obj.Indentation{i,j},1);
                         Hvalueslinfit = polyval(linearfitH,obj.SegTime{j});
                         
                         %y =k*x+d;
                         obj.slopeF(i,j) = linearfitF(1);
                         obj.slopeH(i,j) = linearfitH(1);
                         
                         obj.interceptF(i,j) = linearfitF(2);
                         obj.interceptH(i,j) = linearfitH(2);

                         
                         %manually detrended data
                         ForceTrend{i,j} = obj.BasedForce{i,j} - Fvalueslinfit;
                         IndentTrend{i,j} = obj.Indentation{i,j} - Hvalueslinfit;
                         
                         
                         %1D-digital filter
                         iN = Invariance/100;
                         d = ones(1,fix(iN))/iN;
                         obj.FilterF{i,j} = filtfilt(d,1,ForceTrend{i,j});
                         
                         % Modification of Indentation-Data
                         %IndentTrend{i,j} = detrend(obj.HZShift{i,j});
                         iN = Invariance/100;
                         d = ones(1,fix(iN))/iN;
                         obj.FilterH{i,j} = filtfilt(d,1,IndentTrend{i,j});
                         
                         
                         % Amplitude Correction
                         AmplFiltF = trapz(abs(obj.FilterF{i,j}));
                         AmplOrigF = trapz(abs(ForceTrend{i,j}));
                         AmplCorrectionF = AmplOrigF/AmplFiltF;
                         obj.FilterF{i,j} = obj.FilterF{i,j}*AmplCorrectionF;
                         
                         AmplFiltH = trapz(abs(obj.FilterH{i,j}));
                         AmplOrigH = trapz(abs(IndentTrend{i,j}));
                         AmplCorrectionH = AmplOrigH/AmplFiltH;
                         obj.FilterH{i,j} = obj.FilterH{i,j}*AmplCorrectionH;
                         

                         % Max values of Force and Indentation
                         maxF = max(obj.FilterF{i,j});
                         maxH = max(obj.FilterH{i,j});

                         % Min values of Force and Indentation
                         minF = min(obj.FilterF{i,j});
                         minH = min(obj.FilterH{i,j});

                         % Difference max min
                         DiffF = maxF - minF;
                         DiffH = maxH - minH;

                         % Amplitude
                         AmplitudeF=(DiffF/2);
                         AmplitudeH=(DiffH/2);
                         

                         % Estimate offset
                         %meanF = mean(obj.FilterF{i,j});
                         %meanH = mean(obj.FilterH{i,j});

                         x = obj.SegTime{j};

                         % Function to fit force data 
                         %b(1) (max-min)/2 b(2) FFT b(3) first sign change b(4) mean
                         fit = @(b,x)  b(1).*(sin(2*pi*x.*obj.SegFrequency{j} + b(3)));    
                         % Least-Squares cost function:
                         fcn = @(b) sum((fit(b,x) - obj.FilterF{i,j}).^2);       
                         % Minimise Least-Squares with estimated start values:

                         options = optimset('FunValCheck','off');
                         x0 = [AmplitudeF,obj.SegFrequency{j},0];
                         lb = [0,-Inf,-pi];
                         ub = [Inf,Inf,pi];
                         obj.SineVarsF{i,j} = fmincon(fcn,x0,[],[],[],[],lb,ub); 
                         % Spacing of time vector:
                         %xpF = linspace(min(obj.InterpTimeF{j}),max(obj.InterpTimeF{j}),100000);
                         %obj.SineVarsF{i,j}(1)= obj.SineVarsF{i,j}(1)*rangeF;
                         obj.SineVarsF{i,j}(2)= obj.SegFrequency{j};
                         %obj.SineVarsF{i,j}(3)= firstsignchangeF;

                         
                         % Function to fit indentation data 
                         %b(1) (max-min)/2 b(2) FFT b(3) first sign change b(4) mean
                         fit = @(a,x)  a(1).*(sin(2*pi*x.*obj.SegFrequency{j} + a(3)));    
                         % Least-Squares cost function:
                         fcn = @(a) sum((fit(a,x) - obj.FilterH{i,j}).^2);       
                         % Minimise Least-Squares with estimated start values:
                         options = optimset('FunValCheck','off');
                         lb = [0,-Inf,-pi];
                         ub = [Inf,Inf,pi];
                         obj.SineVarsH{i,j} = fmincon(fcn, [AmplitudeH; obj.SegFrequency{j}; 0],[],[],[],[],lb,ub,[],options); 
                         % Spacing of time vector:
                         %xpH = linspace(min(obj.InterpTimeH{j}),max(obj.InterpTimeH{j}),100000);
                         %obj.SineVarsH{i,j}(1)= obj.SineVarsH{i,j}(1)*rangeH;
                         obj.SineVarsH{i,j}(2)= obj.SegFrequency{j};
                         %obj.SineVarsH{i,j}(3)= firstsignchangeH;
                         
                         
                         %multiply range back 
                         obj.SineVarsF{i,j}(1) = obj.SineVarsF{i,j}(1)*rangeF;
                         obj.SineVarsH{i,j}(1) = obj.SineVarsH{i,j}(1)*rangeH;
                         
                         
                        % phase shift of force and indentation
                         obj.psF{i,j} = obj.SineVarsF{i,j}(3);
                         obj.psH{i,j} = obj.SineVarsH{i,j}(3);
                        
                        % phase shift between indentation and force in degrees:
                        obj.DeltaPhi{i,j} = (obj.psF{i,j} - obj.psH{i,j})*(180/pi);
                        
                        % loss tangent:
                        obj.LossTangent{i,j} = tand(obj.DeltaPhi{i,j});
                        
                        % approximated dynamic modulus:
%                         EModMicro1(i,j) = sqrt(pi/obj.IndentArea(i))*1/2*...
%                                 (1-obj.PoissonR^2)*cos(obj.DeltaPhi{i,j})*...
%                                 (obj.SineVarsF{i,j}(1)/obj.SineVarsH{i,j}(1));
                       
                        
                        % turn range back to normal
                        obj.interceptF(i,j) = obj.interceptF(i,j)*rangeF;
                        obj.interceptH(i,j) = obj.interceptH(i,j)*rangeH;
                        obj.slopeF(i,j) = obj.slopeF(i,j)*rangeF;
                        obj.slopeH(i,j) = obj.slopeH(i,j)*rangeH;
                        obj.BasedForce{i,j} = obj.BasedForce{i,j}*rangeF;
                        obj.Indentation{i,j} = obj.Indentation{i,j}*rangeH;
                        
                        

                        
                    end
                end
            end
           
        end  
        
        function indentation(obj)
            
            
            for i=1:obj.NCurves
                %lastseg = obj.NumSegments-1;
                for j=1:obj.NumSegments
                    Z{i,j} = obj.Height{i,j} - obj.CP(i,1);
                    D{i,j} = (obj.BasedForce{i,j} - obj.CP(i,2))./obj.SpringConstant;
                    obj.Indentation{i,j} = Z{i,j} - D{i,j};
                end
                
            end
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
            [V,D] = eig(R'*R);
            %Extract the output from the eigenvectors
            n = V(:,1);
            V = V(:,2:end);
        end
        
        function [HeadHeight ,Force,spring_constant,sensitivity]=...
                writedata(HeaderFileDirectory,SegmentHeaderFileDirectory,...
                HeightDataDirectory,vDelfDataDirectory,HHType)
            % Author: Orestis Andriotis (slightly changed and adapted by Manuel Rufin)
            % HeaderFileDirectory: file directory of the main header properties. There
            % is the information about the scaling factors to convert the raw data.
            % SegmentHeaderFileDirectory: file directory of the header of each segment.
            % Each segment means the loading (#0) and unloading (#1).
            
            % HeightDataDirectory: file directory of capacitiveSensorHeight.dat file
            % vDelfDataDirectory: file directory of vDeflection.dat file
            %
            
            % HeaderFileDirectory = headerDir;
            % SegmentHeaderFileDirectory=segHeaderDir;
            % HeightDataDirectory=heightdir;
            % vDelfDataDirectory=vDefldir;
            % find the number of total points
            A=fileread(SegmentHeaderFileDirectory);
            fileID=fopen(SegmentHeaderFileDirectory,'rt','n','UTF-8');
            % fileID = fopen(filename,permission,machinefmt,encodingIn)
            
            B=strfind(A,'force-segment-header.num-points=');
            % strfind(file,string) is looking for a specific string in the file.
            fseek(fileID,B,'cof');
            % moves at the location where specific string is located
            
            tline = fgetl(fileID);
            where=strfind(tline,'=');
            % "where" is the position of the "=" symbol in the tline string
            n = str2double(tline(where+1:end)); % number of points
            fclose(fileID);
            
            % Read the height measured data and the vertical deflection.
            % Reading the raw height data into a column of length n
            fileID = fopen(HeightDataDirectory);
            % fread(fileID,sizeA,precision,skip,machinefmt)
            RawHeight = fread(fileID,n,'int32',0,'b'); %raw data
            fclose(fileID);
            
            
            % Reading the deflection height data into a column of length n
            fileID = fopen(vDelfDataDirectory);
            RawvDeflection = fread(fileID,n,'int32',0,'b'); %raw data
            
            
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
            
            clear tline SegmentHeaderFileDirectory HeightDataDirectory...
                vDelfDataDirectory mult_height_volts offset_height_volts...
                mult_height_meters offset_height_meters...
                mult_vDefl_volts offset_vDefl_volts...
                mult_vDefl_meters ;
            fclose(fileID);
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
            fileID=fopen(filedirectory,'rt','n','UTF-8'); % fileID = fopen(filename,permission,machinefmt,encodingIn)
            A=fileread(filedirectory);
            % Height: 1. CONVERSION raw-meters & 2. SCALE meters
            % Conversion RAW -> VOLTS
            fseek(fileID,1,'cof'); % goes at the first position in the file
            
            exp1='\d{1}';
            
            frewind(fileID);
            B=strfind(A,HHType);
            tline = A(B:end);
            HHNum = regexp(tline,exp1,'match','once');
            
            clear tline;
            frewind(fileID);
            B=strfind(A,'vDeflection');
            tline = A(B:end);
            vDefNum = regexp(tline,exp1,'match','once');
            
            %   Multiplier
            clear tline;
            frewind(fileID);
            B=strfind(A,strcat('lcd-info.',HHNum,'.encoder.scaling.multiplier='));
            % strfind(file,string) is looking for a specific string in the file.
            fseek(fileID,B,'cof');
            % moves at the location where specific string is located
            tline = fgetl(fileID);
            % stores that string in a character
            where=strfind(tline,'=');
            % "where" is the position of the "=" symbol in the tline string
            mult_height_meters1 = str2double(... % convert the string to number
                tline(where+1:end)... % this is the number
                );
            
            %   Offset
            clear tline where;
            frewind(fileID);
            B=strfind(A,strcat('lcd-info.',HHNum,'.encoder.scaling.offset='));
            fseek(fileID,B,'cof');
            tline = fgetl(fileID);
            where=strfind(tline,'=');
            offset_height_meters1 = str2double(tline(where+1:end));
            
            %   Multiplier
            clear tline where;
            frewind(fileID);
            B=strfind(A,strcat('lcd-info.',HHNum,'.conversion-set.conversion.nominal.scaling.multiplier='));
            fseek(fileID,B,'cof');
            tline = fgetl(fileID);
            where=strfind(tline,'=');
            mult_height_meters2 = str2double(tline(where+1:end));
            
            
            %   Offset
            clear tline where;
            frewind(fileID);
            B=strfind(A,strcat('lcd-info.',HHNum,'.conversion-set.conversion.nominal.scaling.offset='));
            fseek(fileID,B,'cof');
            tline = fgetl(fileID);
            where=strfind(tline,'=');
            offset_height_meters2 = str2double(tline(where+1:end));
            
            
            % vDeflection: 1. CONVERSION raw-volts & 2. volts to meters
            % Conversion RAW -> VOLTS
            
            %   Multiplier
            clear tline where;
            frewind(fileID);
            B=strfind(A,strcat('lcd-info.',vDefNum,'.encoder.scaling.multiplier='));
            fseek(fileID,B,'cof');
            tline = fgetl(fileID);
            where=strfind(tline,'=');
            mult_vDefl_volts = str2double(tline(where+1:end)); % multiplier for scaling the raw height data and convert to volts
            
            %   Offset
            clear tline where;
            frewind(fileID);
            B=strfind(A,strcat('lcd-info.',vDefNum,'.encoder.scaling.offset='));
            fseek(fileID,B,'cof');
            tline = fgetl(fileID);
            where=strfind(tline,'=');
            offset_vDefl_volts = str2double(tline(where+1:end));
            
            
            % Conversion VOLTS -> METERS
            
            %   Multiplier (that is the sensitivity measured in meters per Volts)
            clear tline where;
            frewind(fileID);
            B=strfind(A,strcat('lcd-info.',vDefNum,'.conversion-set.conversion.distance.scaling.multiplier='));
            fseek(fileID,B,'cof');
            tline = fgetl(fileID);
            where=strfind(tline,'=');
            sensitivity = str2double(tline(where+1:end));
            
            % Spring constant
            
            clear tline where;
            frewind(fileID);
            B=strfind(A,strcat('lcd-info.',vDefNum,'.conversion-set.conversion.force.scaling.multiplier='));
            fseek(fileID,B,'cof');
            tline = fgetl(fileID);
            where=strfind(tline,'=');
            spring_constant = str2double(tline(where+1:end));
            
            clear tline A B where
            
            fclose(fileID);
        end
        
        function X = CP_CNN_batchprep(objcell,ImgSize)
            % This function takes as input a 1xN objcell of N objects of the class
            % 'ForceMap' and combines them into the matrizes needed for NN-training
            if nargin < 2
                ImgSize = 128;
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
            X = zeros(ImgSize,ImgSize,1,Nimgs);
            
            fig = figure('Color','w');
            k = 1;
            for i=1:Nmaps
                jRange = find(objcell{i}.SelectedCurves);
                for j=jRange'
                    % Save the plots as images and convert them into cropped [0
                    % 1]-range grayscale images
                    
                    figure(fig)
                    fig.Name = sprintf('%s curve nr.%i',objcell{i}.Name,j);
                    plot(objcell{i}.THApp{j},objcell{i}.BasedApp{j},'color','black');
                    axis([min(objcell{i}.THApp{j}) max(objcell{i}.THApp{j}) min(objcell{i}.BasedApp{j}) max(objcell{i}.BasedApp{j})])
                    axis off
                    graytest = imcomplement(rgb2gray(frame2im(getframe(fig))));
                    grayscaled = double(graytest)/255;
                    % Crop off the rows and columns of the image only containing zeros
                    crop = [1 size(grayscaled,2) 1 size(grayscaled,1)];
                    while sum(grayscaled(:,crop(1)))== 0
                        crop(1) = crop(1) + 1;
                    end
                    crop(1) = crop(1) - 1;
                    while sum(grayscaled(:,crop(2)))== 0
                        crop(2) = crop(2) - 1;
                    end
                    crop(2) = crop(2) + 1;
                    while sum(grayscaled(crop(3),:))== 0
                        crop(3) = crop(3) + 1;
                    end
                    crop(3) = crop(3) - 1;
                    while sum(grayscaled(crop(4),:))== 0
                        crop(4) = crop(4) - 1;
                    end
                    crop(4) = crop(4) + 1;
                    grayscaled(:,[1:crop(1) crop(2):size(grayscaled,2)]) = [];
                    grayscaled([1:crop(3) crop(4):size(grayscaled,1)],:) = [];
                    grayfinal = imresize(grayscaled,[ImgSize ImgSize],'bilinear');
                    % Fill the output variable X
                    X(:,:,1,k) = grayfinal;
                    k = k + 1;
                end
            end
            close(fig);
            
        end
        
        function X = CP_oliver_pharr_batchprep(objcell,ImgSize)
            % This function takes as input a 1xN objcell of N objects of the class
            % 'ForceMap' and converts them into 128x128 images needed for
            % NN-inference
            if nargin < 2
                ImgSize = 128;
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
            X = zeros(ImgSize,ImgSize,1,Nimgs);
            
            fig = figure('Color','w');
            k = 1;
            for i=1:Nmaps
                jRange = find(objcell{i}.SelectedCurves);
                for j=jRange'
                    % Save the plots as images and convert them into cropped [0
                    % 1]-range grayscale images
                    
                    figure(fig)
                    fig.Name = sprintf('%s curve nr.%i',objcell{i}.Name,j);
                    plot(objcell{i}.HHApp{j},objcell{i}.BasedApp{j},'color','black');
                    axis([min(objcell{i}.HHApp{j}) max(objcell{i}.HHApp{j}) min(objcell{i}.BasedApp{j}) max(objcell{i}.BasedApp{j})])
                    axis off
                    graytest = imcomplement(rgb2gray(frame2im(getframe(fig))));
                    grayscaled = double(graytest)/255;
                    % Crop off the rows and columns of the image only containing zeros
                    crop = [1 size(grayscaled,2) 1 size(grayscaled,1)];
                    while sum(grayscaled(:,crop(1)))== 0
                        crop(1) = crop(1) + 1;
                    end
                    crop(1) = crop(1) - 1;
                    while sum(grayscaled(:,crop(2)))== 0
                        crop(2) = crop(2) - 1;
                    end
                    crop(2) = crop(2) + 1;
                    while sum(grayscaled(crop(3),:))== 0
                        crop(3) = crop(3) + 1;
                    end
                    crop(3) = crop(3) - 1;
                    while sum(grayscaled(crop(4),:))== 0
                        crop(4) = crop(4) - 1;
                    end
                    crop(4) = crop(4) + 1;
                    grayscaled(:,[1:crop(1) crop(2):size(grayscaled,2)]) = [];
                    grayscaled([1:crop(3) crop(4):size(grayscaled,1)],:) = [];
                    grayfinal = imresize(grayscaled,[ImgSize ImgSize],'bilinear');
                    % Fill the output variable X
                    X(:,:,1,k) = grayfinal;
                    k = k + 1;
                end
            end
            close(fig);
            
        end
        
        function X = CP_batchprep_new(objcell,ImgSizeFinal,ImgSize)
            
            if nargin<2
                ImgSize = 478;
                ImgSizeFinal = 128;
            elseif nargin < 3
                ImgSize = 478;
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
            X = zeros(ImgSizeFinal,ImgSizeFinal,1,Nimgs);
            
            k = 1;
            for i=1:Nmaps
                jRange = find(objcell{i}.SelectedCurves);
                for j=jRange'
                    
                    Image = zeros(ImgSize,ImgSize);
                    
                    Points(:,1) = objcell{i}.HHApp{j};
                    Points(:,2) = objcell{i}.BasedApp{j};
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
                    X(:,:,1,k) = Image;
                    k = k + 1;
                    clear Points
                end
            end
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
                CutPercent = [0 0.5 0.75];
            elseif nargin < 3
                ImgSize = 478;
                CutPercent = [0 0.5 0.75];
            elseif nargin<4
                CutPercent = [0 0.5 0.75];
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
                        
                        Points(:,1) = objcell{i}.HHApp{j}(floor(n*length(objcell{i}.HHApp{j}))+1:end);
                        Points(:,2) = objcell{i}.BasedApp{j}(floor(n*length(objcell{i}.BasedApp{j}))+1:end);
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
        
    end
    
    methods
        % non-static auxiliary methods
        
        
        function save(obj)
            current = what();
            cd(obj.Folder)
            savename = sprintf('%s.mat',obj.Name);
            save(savename,'obj','-v7.3')
            cd(current.path)
            savemsg = sprintf('Changes to ForceMap %s saved to %s',obj.Name,obj.Folder);
            disp(savemsg);
        end
        
        function [E_mod,GoF,Hertzfit] = hertz_fit_gof(obj,tip_h,force,CP,curve_percent,shape)
            
            if obj.TipRadius == -1
                prompt = {'What is the nominal tip radius of the used cantilever in nm?'};
                dlgtitle = 'Cantilever tip';
                dims = [1 35];
                definput = {'10'};
                obj.TipRadius = str2double(inputdlg(prompt,dlgtitle,dims,definput));
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
                E_mod = 3*(Hertzfit.a*ranf/rant^(3/2))/(4*sqrt(obj.TipRadius*10^(-9)))*(1-obj.PoissonR^2);
            elseif isequal(shape,'spherical')
            elseif isequal(shape,'conical')
            elseif isequal(shape,'pyramid')
            end
        end
        
        function calculate_reference_slope_from_area(obj,Mask)
            % Calculates the distribution of DZslopes on the curves
            % that that are neihter on the fibril nor the excluded zones.
            %  the upper 25% of the curve are considered for the
            %  calculation
            
            Range = find(obj.SelectedCurves);
            CurvePercent = 0.25;
            % Calculate the DZslopes
            k = 1;
            for i=Range'
                if (Mask(obj.List2Map(i,1),obj.List2Map(i,2)) == 1) &&...
                        (obj.ExclMask(obj.List2Map(i,1),obj.List2Map(i,2)) == 1)
                    Z(:,i) = obj.HHRet{i} - obj.CP(i,1);
                    D(:,i) = (obj.BasedRet{i} - obj.CP(i,2))/obj.SpringConstant;
                    Dmax = max(D(:,i));
                    DCurvePercent = D(D(:,i)>=(1-CurvePercent)*Dmax,i);
                    ZCurvePercent = Z(1:length(DCurvePercent),i);
                    LineFit = polyfit(ZCurvePercent,DCurvePercent,1);
                    DZslope(k) = LineFit(1);
                    k = k + 1;
                end
            end
            % Fit the Gaussian
            Gaussian = fitdist(DZslope','Normal');
            obj.RefSlope = Gaussian.mean;
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
                obj.HasRefSlope = true;
                IsValidInput = true;
            end
        end
        
        function set_reference_slope_to_value(obj,Value)
            obj.RefSlope = Value; % BEST.FUNCTION.EVER.WRITTEN.
            obj.HasRefSlope = true;
        end
        
        function Mask = create_mask_general(obj)
            
            Mask = logical(zeros(obj.NumProfiles,obj.NumPoints));
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
        
        function create_and_level_height_map(obj)
            % first set Height Map to default values for reproducable
            % results
            k = 1;
            obj.List2Map = zeros(obj.NCurves,2);
            if isequal(obj.FileType,'quantitative-imaging-map')
                for i=1:obj.NumProfiles
                    for j=1:obj.NumPoints
                        obj.Map2List(i,j) = k;
                        obj.List2Map(k,:) = [i j];
                        k = k + 1;
                    end
                end
            elseif isequal(obj.FileType,'force-scan-map')
                for i=1:obj.NumProfiles
                    if ~mod(i,2)
                        for j=1:obj.NumPoints
                            obj.Map2List(i,j) = k;
                            obj.List2Map(k,:) = [i j];
                            k = k + 1;
                        end
                    else
                        for j=1:obj.NumPoints
                            obj.Map2List(i,obj.NumPoints-j+1) = k;
                            obj.List2Map(k,:) = [i obj.NumPoints-j+1];
                            k = k + 1;
                        end
                    end
                end
            end
            Max = zeros(obj.NCurves,1);
            for i=1:obj.NCurves
                Max(i) = -max(obj.HHApp{i});
            end
            obj.HeightMap = obj.convert_data_list_to_map(Max);
            % Create an Exclusion Mask with the standard deviation method
            MaskParam = 0.65;
            HghtRange = range(obj.HeightMap,'all');
            HghtMin = min(obj.HeightMap,[],'all');
            HghtNorm = (obj.HeightMap-HghtMin)/HghtRange;
            mask = zeros(size(HghtNorm));
            STDLine = zeros(obj.NumProfiles,obj.NumPoints);
            [HghtSorted,Idx] = sort(HghtNorm,[2],'descend');
            for j=1:obj.NumPoints
                STDLine(:,j) = std(HghtSorted(:,1:j),0,[2]);
            end
            [~,MaxIdx] = max(STDLine,[],[2]);
            for i=1:obj.NumProfiles
                mask(i,Idx(i,1:floor(MaxIdx*MaskParam))) = 1;
            end
            mask = logical(~mask);
            
            masked_map = mask(:,:,1).*obj.HeightMap;
            % create a N-by-3 matrix with each column being a point in
            % 3D-space
            k = 1;
            HghtVctr = zeros(sum(mask,'all'),3);
            for i=1:size(masked_map,1)
                for j=1:size(masked_map,2)
                    if mask(i,j) == 0
                    else
                        HghtVctr(k,:) = [size(masked_map,2)/size(masked_map,1)*i,size(masked_map,1)/size(masked_map,2)*j,masked_map(i,j)];
                        k = k + 1 ;
                    end
                end
            end
            % fit a plane to the glass part
            [Norm,~,Point] = obj.affine_fit(HghtVctr);
            Plane = zeros(size(masked_map,1),size(masked_map,2));
            % Create the plane that can then be subtracted from the
            % complete height data to generate the leveled height data.
            for i=1:size(masked_map,1)
                for j=1:size(masked_map,2)
                    Plane(i,j) = (Point(3)-Norm(1)/Norm(3)*(size(masked_map,2)/size(masked_map,1)*i)-Norm(2)/Norm(3)*(size(masked_map,1)/size(masked_map,2)*j));
                    Plane(isnan(Plane)) = Point(3);
                end
            end
            obj.HeightMap = obj.HeightMap - Plane;
        end
        
        function create_fibril_mask(obj,MaskParam)
            
            if nargin < 2
                MaskParam = 0.5;
            end
            HghtRange = range(obj.HeightMap,'all');
            HghtMin = min(obj.HeightMap,[],'all');
            HghtNorm = (obj.HeightMap-HghtMin)/HghtRange;
            mask = zeros(size(HghtNorm));
            STDLine = zeros(obj.NumProfiles,obj.NumPoints);
            [HghtSorted,Idx] = sort(HghtNorm,[2],'descend');
            for j=1:obj.NumPoints
                STDLine(:,j) = std(HghtSorted(:,1:j),0,[2]);
            end
            [~,MaxIdx] = max(STDLine,[],[2]);
            for i=1:obj.NumProfiles
                mask(i,Idx(i,1:floor(MaxIdx*MaskParam))) = 1;
            end
            mask = logical(mask);
            mask = bwareafilt(mask,1,4);
            
            % determine linear orientation angle of the fibril and pad
            % accordingly
            k = 1;
            for i=1:obj.NumProfiles
                for j = 1:obj.NumPoints
                    if mask(i,j) == 1
                        xy(k,:) = [i,j];
                        k = k + 1;
                    end
                end
            end
            FitLine = fit(xy(:,1).*obj.XSize/obj.NumProfiles,xy(:,2).*obj.YSize/obj.NumPoints,'poly1');
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
        
        function create_automatic_background_mask(obj,MaskParam)
            
            if nargin < 2
                MaskParam = 0.5;
            end
            mask = zeros(size(obj.HeightMap));
            HeightList = zeros(obj.NCurves,1);
            for i=1:obj.NCurves
                HeightList(i) = -obj.HHApp{i}(end);
            end
            STDLine = zeros(obj.NCurves,1);
            [HeightSorted,Idx] = sort(HeightList,'descend');
            for i=1:obj.NCurves
                STDLine(i) = std(HeightSorted(1:i));
            end
            [~,MaxIdx] = max(STDLine);
            MapIndex = obj.List2Map(Idx(1:floor(MaxIdx*MaskParam)),:);
            for i=1:length(MapIndex)
                mask(MapIndex(i,1),MapIndex(i,2)) = 1;
            end
            mask = logical(mask);
            mask = bwareafilt(mask,1,4);
            
            obj.BackgroundMask = ~mask;
            
            %             current = what();
            %             cd(obj.Folder)
            %             savename = sprintf('%s.mat',obj.Name);
            %             save(savename,'obj')
            %             cd(current.path)
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
                obj.MiniBatchSize = len;
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
                Lb = abs(min(obj.HHApp{i})-obj.CP_CNNZoom(i,1));
                La = abs(max(obj.HHApp{i})-obj.CP_CNNZoom(i,1));
                L = Lb + La;
                if Lb/L <= ZoomFactor
                    continue
                end
                Lsub = (Lb - ZoomFactor*L)/(1-ZoomFactor);
                CutOff = min(obj.HHApp{i}) + Lsub;
                obj.HHApp{i}(obj.HHApp{i}<CutOff) = [];
                obj.BasedApp{i}(1:(end-length(obj.HHApp{i}))) = [];
            end
            
            
        end
        
        function TempFolder = unpack_jpk_force_map(obj,MapFullFile,DataFolder)
            
            if isequal('PCW',obj.HostOS)
                % unpack jpk-file into temporary folder to read out data
                cmd1 = '"C:\Program Files\7-Zip\7z.exe" x ';
                cmd2 = '"';
                cmd3 = MapFullFile;
                cmd4 = '"';
                cmd5 = ' -o';
                TempFolderName = sprintf('Temp%s',obj.ID);
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
                TempFolderName = sprintf('Temp%s',obj.ID);
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
                %%% Create a data folder to store the force data
                mkdir(DataFolder,'ForceData')
                obj.Folder = fullfile(DataFolder,'ForceData',filesep);
                
            elseif isequal('MAC',obj.HostOS)
                % unpack jpk-file into temporary folder to read out data
                cmd1 = 'unzip -o ';
                cmd2 = MapFullFile;
                cmd3 = ' -d ';
                TempFolderName = sprintf('Temp%s',obj.ID);
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
                %%% Create a data folder to store the force data
                mkdir(DataFolder,'ForceData')
                obj.Folder = fullfile(DataFolder,'ForceData',filesep);
            end
        end
        
        function read_in_header_properties(obj,TempFolder)
            % Check for jpk-software version and get important ForceMap
            % properties
            
            filedirectory = fullfile(TempFolder,'header.properties');
            fileID=fopen(filedirectory,'rt','n','UTF-8'); % fileID = fopen(filename,permission,machinefmt,encodingIn)
            A=fileread(filedirectory);
            % Height: 1. CONVERSION raw-meters & 2. SCALE meters
            % Conversion RAW -> VOLTS
            fseek(fileID,1,'cof'); % goes at the first position in the file
            
            %   Check for file type (.jpk-force-map, .jpk-qi-data)
            frewind(fileID);
            B=strfind(A,'jpk-data-file=');
            fseek(fileID,B,'cof');
            tline = fgetl(fileID);
            where=strfind(tline,'=');
            TempType = tline(where+1:end);
            if isequal(TempType, 'spm-quantitative-image-data-file')   % Valid for software versions 6.1.158  
                obj.FileType = 'quantitative-imaging-map';
            elseif isequal(TempType,'spm-force-scan-map-file')         % Valid for software versions 6.1.158 
                obj.FileType = 'force-scan-map';
            end
            
            %   Check for file version
            clear tline where;
            frewind(fileID);
            B=strfind(A,strcat(obj.FileType,'.description.source-software='));
            fseek(fileID,B,'cof');
            tline = fgetl(fileID);
            where=strfind(tline,'=');
            obj.FileVersion = tline(where+1:end);
            
            %   NCurves
            clear tline where;
            frewind(fileID);
            B=strfind(A,strcat(obj.FileType,'.indexes.max='));
            % strfind(file,string) is looking for a specific string in the file.
            fseek(fileID,B,'cof');
            % moves at the location where specific string is located
            tline = fgetl(fileID);
            % stores that string in a character
            where=strfind(tline,'=');
            % "where" is the position of the "=" symbol in the tline string
            obj.NCurves = 1 + str2double(... % convert the string to number
                tline(where+1:end)... % this is the number
                );
            
            % NumSegments    
            clear tline where;
            frewind(fileID);
            B=strfind(A,strcat(obj.FileType,'.settings.force-settings.segments.size='));
            if isempty(B)
                clear tline where;
                frewind(fileID);
            else
                fseek(fileID,B,'cof');
                tline = fgetl(fileID);
                where=strfind(tline,'=');
                obj.NumSegments = str2double(tline(where+1:end)); %not str2num?
                
                for i=1:obj.NumSegments
                    clear tline where;
                    frewind(fileID);
                    B=strfind(A,strcat(obj.FileType,'.settings.force-settings.segment.',string((i-1)),'.duration='));
                    fseek(fileID,B,'cof');
                    tline = fgetl(fileID);
                    where=strfind(tline,'=');
                    obj.SegDuration{i} = str2double(tline(where+1:end)); %not str2num?
                    
                    if i == 1
                        obj.SeriesTime{i} = obj.SegDuration{i};
                    else
                        obj.SeriesTime{i} = obj.SeriesTime{i-1}+obj.SegDuration{i};
                    end
                    
                    clear tline where;
                    frewind(fileID);
                    B=strfind(A,strcat(obj.FileType,'.settings.force-settings.segment.',string((i-1)),'.num-points='));
                    fseek(fileID,B,'cof');
                    tline = fgetl(fileID);
                    where=strfind(tline,'=');
                    obj.SegNumPoints{i} = str2double(tline(where+1:end));
                    
                    obj.SecPerPoint{i} = obj.SegDuration{i}/obj.SegNumPoints{i};
                    
                    if i == 1
                        obj.TStart{i} = obj.SecPerPoint{i}/2;
                    else 
                        obj.TStart{i} = obj.SeriesTime{i-1}+(obj.SecPerPoint{i}/2);
                    end
                    
                    obj.TEnd{i} = obj.SeriesTime{i};
                    obj.SegTime{i} = obj.TStart{i}:obj.SecPerPoint{i}:obj.TEnd{i};
                    obj.InterpTimeF{i} = obj.TStart{i}:0.000001:obj.TEnd{i};
                    obj.InterpTimeH{i} = obj.TStart{i}:0.000001:obj.TEnd{i};
                    obj.InterpTimeF{i} = obj.InterpTimeF{i}.';
                    obj.InterpTimeH{i} = obj.InterpTimeH{i}.';
                    obj.SegTime{i} = obj.SegTime{i}.';
                        
                        
                    clear tline where;
                    frewind(fileID);
                    B=strfind(A,strcat(obj.FileType,'.settings.force-settings.segment.',string((i-1)),'.frequency='));
                    if isempty(B)
                        obj.SegFrequency{i} = 0.0;
                    else
                        fseek(fileID,B,'cof');
                        tline = fgetl(fileID);
                        where=strfind(tline,'=');
                        obj.SegFrequency{i} = str2double(tline(where+1:end));
                    end
                end
            end

                
                
                
            
            %   NumProfiles
            clear tline where;
            frewind(fileID);
            B=strfind(A,strcat(obj.FileType,'.position-pattern.grid.jlength='));
            fseek(fileID,B,'cof');
            tline = fgetl(fileID);
            where=strfind(tline,'=');
            obj.NumProfiles = str2double(tline(where+1:end));
            
            %   NumPoints
            clear tline where;
            frewind(fileID);
            B=strfind(A,strcat(obj.FileType,'.position-pattern.grid.ilength='));
            fseek(fileID,B,'cof');
            tline = fgetl(fileID);
            where=strfind(tline,'=');
            obj.NumPoints = str2double(tline(where+1:end));
            
            %   XSize
            clear tline where;
            frewind(fileID);
            B=strfind(A,strcat(obj.FileType,'.position-pattern.grid.ulength='));
            fseek(fileID,B,'cof');
            tline = fgetl(fileID);
            where=strfind(tline,'=');
            obj.XSize = str2double(tline(where+1:end));
            
            %   YSize
            clear tline where;
            frewind(fileID);
            B=strfind(A,strcat(obj.FileType,'.position-pattern.grid.vlength='));
            fseek(fileID,B,'cof');
            tline = fgetl(fileID);
            where=strfind(tline,'=');
            obj.YSize = str2double(tline(where+1:end));
            
            %   MaxPointsPerCurve
            clear tline where;
            frewind(fileID);
            if isequal(obj.FileType,'force-scan-map') && isempty(obj.NumSegments)
                B=strfind(A,strcat(obj.FileType,'.settings.force-settings.extend-k-length='));
            elseif isequal(obj.FileType,'quantitative-imaging-map')
                B=strfind(A,strcat(obj.FileType,'.settings.force-settings.extend.num-points='));
            end
            fseek(fileID,B,'cof');
            tline = fgetl(fileID);
            where=strfind(tline,'=');
            obj.MaxPointsPerCurve = str2double(tline(where+1:end));
            
            %   Velocity
            
            if isequal(obj.FileType,'force-scan-map') && isempty(obj.NumSegments)
                clear tline where;
                frewind(fileID);
                B=strfind(A,strcat(obj.FileType,'.settings.force-settings.extend-scan-time='));
                fseek(fileID,B,'cof');
                tline = fgetl(fileID);
                where=strfind(tline,'=');
                ExtendTime = str2double(tline(where+1:end));
                
                clear tline where;
                frewind(fileID);
                B=strfind(A,strcat(obj.FileType,'.settings.force-settings.relative-z-start='));
                fseek(fileID,B,'cof');
                tline = fgetl(fileID);
                where=strfind(tline,'=');
                ZLength = str2double(tline(where+1:end));
                obj.Velocity = ZLength/ExtendTime;
            elseif isequal(obj.FileType,'quantitative-imaging-map')
                clear tline where;
                frewind(fileID);
                B=strfind(A,'quantitative-imaging-map.settings.force-settings.extend.duration=');
                fseek(fileID,B,'cof');
                tline = fgetl(fileID);
                where=strfind(tline,'=');
                ExtendTime = str2double(tline(where+1:end));
                
                clear tline where;
                frewind(fileID);
                B=strfind(A,'quantitative-imaging-map.settings.force-settings.extend.z-start=');
                fseek(fileID,B,'cof');
                tline = fgetl(fileID);
                where=strfind(tline,'=');
                ZLength = str2double(tline(where+1:end));
                obj.Velocity = ZLength/ExtendTime;
            end
            
            
            
            %   GridAngle
            clear tline where;
            frewind(fileID);
            B=strfind(A,strcat(obj.FileType,'.position-pattern.grid.theta='));
            fseek(fileID,B,'cof');
            tline = fgetl(fileID);
            where=strfind(tline,'=');
            obj.GridAngle = str2double(tline(where+1:end));
            obj.GridAngle = obj.GridAngle*180/pi;
            
            clear tline A B where
            
            fclose(fileID);
        end
        
        function load_force_curves(obj,TempFolder)
            
          if ~isempty(obj.NumSegments)  && (obj.NumSegments > 2)
            obj.HHType = 'capacitiveSensorHeight';
            for i=1:obj.NCurves
                for j=1:obj.NumSegments
                    HeaderFileDirectory = fullfile(TempFolder,'shared-data','header.properties');
                    SegmentHeaderFileDirectory = fullfile(TempFolder,'index',string((i-1)),'segments',string((j-1)),'segment-header.properties');
                    HeightDataDirectory = fullfile(TempFolder,'index',string((i-1)),'segments',string((j-1)),'channels','capacitiveSensorHeight.dat');
                    vDefDataDirectory = fullfile(TempFolder,'index',string((i-1)),'segments',string((j-1)),'channels','vDeflection.dat');

                    if ~isfile(HeightDataDirectory) || isequal(obj.HHType,'measuredHeight')
                        HeightDataDirectory = fullfile(TempFolder,'index',string((i-1)),'segments',string((j-1)),'channels','measuredHeight.dat');
                        obj.HHType = 'measuredHeight';
                    end
                    if ~isfile(HeightDataDirectory) || isequal(obj.HHType,'Height')
                        HeightDataDirectory = fullfile(TempFolder,'index',string((i-1)),'segments',string((j-1)),'channels','Height.dat');
                        obj.HHType = 'Height';
                    end
                    
                    
                    [TempHHApp,obj.Force{i,j},obj.SpringConstant,obj.Sensitivity]=...
                        obj.writedata(HeaderFileDirectory,SegmentHeaderFileDirectory,...
                        HeightDataDirectory,vDefDataDirectory,obj.HHType);

                    obj.Height{i,j} = -TempHHApp;
                    % force = deflection * spring constant
                    obj.Force{i,j} = obj.Force{i,j}.*obj.SpringConstant;
                    clear TempHHApp
                    

                end
                    
                 obj.HHApp{i} = obj.Height{i,1};
                 obj.App{i} = obj.Force{i,1};
                 
 
                 %lastseg = obj.NumSegments - 1;
                 obj.HHRet{i} = obj.Height{i,obj.NumSegments};
                 obj.Ret{i} = obj.Force{i,obj.NumSegments};

            end
         
          else
            
            obj.HHType = 'capacitiveSensorHeight';
            for i=1:obj.NCurves
                HeaderFileDirectory = fullfile(TempFolder,'shared-data','header.properties');
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
                
                [TempHHApp,obj.App{i},obj.SpringConstant,obj.Sensitivity]=...
                    obj.writedata(HeaderFileDirectory,SegmentHeaderFileDirectory,...
                    HeightDataDirectory,vDefDataDirectory,obj.HHType);
                
                obj.HHApp{i} = -TempHHApp;
                obj.App{i} = obj.App{i}.*obj.SpringConstant;
                clear TempHHApp
                
                % Below there is a workaround for jpk-force-map files,
                % where the capacitiveSensorHeight
                % written into an additional segment folder '2' instead of
                % '1'. This occurs, when there is a nonzero holding time,
                % which gets writtem into the '1'-folder instead
                
                % To be implemented: read holding segment into ForceMap
                % class
                
                SegmentHeaderFileDirectory = fullfile(TempFolder,'index',string((i-1)),'segments','1','segment-header.properties');
                HeightDataDirectory = fullfile(TempFolder,'index',string((i-1)),'segments','1','channels','capacitiveSensorHeight.dat');
                vDefDataDirectory = fullfile(TempFolder,'index',string((i-1)),'segments','1','channels','vDeflection.dat');
                
                if isfolder(fullfile(TempFolder,'index',string((i-1)),'segments','2'))
                    SegmentHeaderFileDirectory = fullfile(TempFolder,'index',string((i-1)),'segments','2','segment-header.properties');
                    HeightDataDirectory = fullfile(TempFolder,'index',string((i-1)),'segments','2','channels','capacitiveSensorHeight.dat');
                    vDefDataDirectory = fullfile(TempFolder,'index',string((i-1)),'segments','2','channels','vDeflection.dat');
                    if ~isfile(HeightDataDirectory)  || isequal(obj.HHType,'measuredHeight')
                        HeightDataDirectory = fullfile(TempFolder,'index',string((i-1)),'segments','2','channels','measuredHeight.dat');
                    end
                end
                
                if ~isfile(HeightDataDirectory) || isequal(obj.HHType,'measuredHeight')
                    HeightDataDirectory = fullfile(TempFolder,'index',string((i-1)),'segments','1','channels','measuredHeight.dat');
                end
                
                [TempHHRet,obj.Ret{i}]=...
                    obj.writedata(HeaderFileDirectory,SegmentHeaderFileDirectory,...
                    HeightDataDirectory,vDefDataDirectory,obj.HHType);
                
                obj.HHRet{i} = -TempHHRet;
                obj.Ret{i} = obj.Ret{i}.*obj.SpringConstant;
                clear TempHHRet
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
            obj.SMFSFlag.Min=zeros(1,obj.NCurves);
            obj.SMFSFlag.Length=zeros(1,obj.NCurves);
            
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
            subplot(2,1,1)
            title(sprintf('Curve Nr.%i of %s',k,obj.Name))
            hold on
            [MultiplierX,UnitX,~] = AFMImage.parse_unit_scale(range(obj.HHRet{k}),'m',10);
            [MultiplierY,UnitY,~] = AFMImage.parse_unit_scale(range(obj.BasedRet{k}),'N',10);
            plot(obj.HHApp{k}*MultiplierX,obj.BasedApp{k}*MultiplierY,obj.HHRet{k}*MultiplierX,obj.BasedRet{k}*MultiplierY,'LineWidth',1.5);
            Legends = {'Approach','Retract'};
            
            if obj.CPFlag.HertzFitted == 1
                plot(obj.CP_HertzFitted(k,1)*1e9, obj.CP_HertzFitted(k,2)*1e9,'O',...
                    'LineWidth',1.5,...
                    'MarkerSize',7,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor',[0.4940 0.1840 0.5560]);
                Legends{end+1} = 'Origin of Hertz-Sneddon fit';
                xlim([(obj.HHApp{k}(1)+ZoomMult*(obj.CP_SnapIn(k,1) - obj.HHApp{k}(1)))*1e9 inf]);
            end
            if obj.CPFlag.SnapIn == 1
                plot(obj.CP_SnapIn(k,1)*1e9, obj.CP_SnapIn(k,2)*1e9,'O',...
                    'LineWidth',1.5,...
                    'MarkerSize',7,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','r');
                Legends{end+1} = 'SnapIn';
                xlim([(obj.HHApp{k}(1)+ZoomMult*(obj.CP_SnapIn(k,1) - obj.HHApp{k}(1)))*1e9 inf]);
            end
            if obj.CPFlag.Manual == 1
                plot(obj.Man_CP(k,1)*1e9, obj.Man_CP(k,2)*1e9,'O',...
                    'LineWidth',1.5,...
                    'MarkerSize',7,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','r');
                Legends{end+1} = 'Manual CP';
                xlim([(obj.HHApp{k}(1)+ZoomMult*(obj.Man_CP(k,1) - obj.HHApp{k}(1)))*1e9 inf]);
            end
            if obj.CPFlag.RoV == 1
                plot(obj.CP_RoV(k,1)*1e9, obj.CP_RoV(k,2)*1e9,'O',...
                    'LineWidth',1.5,...
                    'MarkerSize',7,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','b');
                Legends{end+1} = 'CP RoV';
                xlim([(obj.HHApp{k}(1)+ZoomMult*(obj.CP_RoV(k,1) - obj.HHApp{k}(1)))*1e9 inf]);
            end
            if obj.CPFlag.GoF == 1
                plot(obj.CP_GoF(k,1)*1e9, obj.CP_GoF(k,2)*1e9,'O',...
                    'LineWidth',1.5,...
                    'MarkerSize',7,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','c');
                Legends{end+1} = 'CP GoF';
                xlim([(obj.HHApp{k}(1)+ZoomMult*(obj.CP_GoF(k,1) - obj.HHApp{k}(1)))*1e9 inf]);
            end
            if obj.CPFlag.Combo == 1
                plot(obj.CP_Combo(k,1)*1e9, obj.CP_Combo(k,2)*1e9,'O',...
                    'LineWidth',1.5,...
                    'MarkerSize',7,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','m');
                Legends{end+1} = 'CP Combo';
                xlim([(obj.HHApp{k}(1)+ZoomMult*(obj.CP_Combo(k,1) - obj.HHApp{k}(1)))*1e9 inf]);
            end
            if obj.CPFlag.CNN == 1
                plot(obj.CP_CNN(k,1)*1e9, obj.CP_CNN(k,2)*1e9,'O',...
                    'LineWidth',1.5,...
                    'MarkerSize',7,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','y');
                Legends{end+1} = 'CP CNN';
                xlim([(obj.HHApp{k}(1)+ZoomMult*(obj.CP_CNN(k,1) - obj.HHApp{k}(1)))*1e9 inf]);
            end
            if obj.CPFlag.CNNZoom == 1
                plot(obj.CP_CNNZoom(k,1)*1e9, obj.CP_CNNZoom(k,2)*1e9,'O',...
                    'LineWidth',1.5,...
                    'MarkerSize',7,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor',[0.8500 0.3250 0.0980]);
                Legends{end+1} = 'CP CNNZoom';
                xlim([(obj.HHApp{k}(1)+ZoomMult*(obj.CP_CNNZoom(k,1) - obj.HHApp{k}(1)))*1e9 inf]);
            end
            if obj.CPFlag.CNNZoomSweep == 1
                plot(obj.CP_CNNZoomSweep(k,1)*1e9, obj.CP_CNNZoomSweep(k,2)*1e9,'gs',...
                    'LineWidth',1.5,...
                    'MarkerSize',7,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor',[0.8500 0.3250 0.0980]);
                Legends{end+1} = 'CP CNNZoomSweep';
                xlim([(obj.HHApp{k}(1)+ZoomMult*(obj.CP_CNNZoomSweep(k,1) - obj.HHApp{k}(1)))*1e9 inf]);
            end
            if obj.CPFlag.Old == 1
                plot(obj.CP_Old(k,1)*1e9, obj.CP_Old(k,2)*1e9,'O',...
                    'LineWidth',1.5,...
                    'MarkerSize',7,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor',[0.9290 0.6940 0.1250]);
                Legends{end+1} = 'CP SD6';
                xlim([(obj.HHApp{k}(1)+ZoomMult*(obj.CP_Old(k,1) - obj.HHApp{k}(1)))*1e9 inf]);
            end
            if obj.CPFlag.Dropout == 1
                X = obj.CP_Dropout(k,1)*1e9;
                Y = obj.CP_Dropout(k,2)*1e9;
                XSTD = std(obj.CP_MonteCarlo(:,1,k))*1e9;
                YSTD = std(obj.CP_MonteCarlo(:,2,k))*1e9;
                plot(X, Y,'O',...
                    'LineWidth',1.5,...
                    'MarkerSize',7,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor',[0.4940 0.1840 0.5560]);
                plot([X-XSTD X+XSTD],[Y Y],'-+',[X X],[Y-YSTD Y+YSTD],'-+',...
                    'LineWidth',1.5,...
                    'Color',[0.4940 0.1840 0.5560]);
                Legends{end+1} = 'CP Dropout';
                Legends{end+1} = sprintf('Dropout Uncertainty %.2f nm',obj.CP_MonteCarlo_STD(k)*1e9);
                xlim([(obj.HHApp{k}(1)+ZoomMult*(obj.CP_Dropout(k,1) - obj.HHApp{k}(1)))*1e9 inf]);
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
            I = obj.HeightMap;
            I = mat2gray(I);
            imshow(I)
            axis on
            hold on;
            try
                for i=1:obj.NumProfiles
                    plot((obj.List2Map(obj.RectApexIndex(i),2)),...
                        (obj.List2Map(obj.RectApexIndex(i),1)),...
                        'g+', 'MarkerSize', 10, 'LineWidth', 2);
                    %                 plot((obj.List2Map(obj.ApexIndex(i),2)-1/2)*1024/obj.NumPoints,...
                    %                     (obj.List2Map(obj.ApexIndex(i),1)-1/2)*1024/obj.NumProfiles,...
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
        
        function Map = convert_data_list_to_map(obj,List)
            
            Map = zeros(obj.NumProfiles,obj.NumPoints);
            for i=1:obj.NumProfiles
                for j=1:obj.NumPoints
                    Map(i,j) = List(obj.Map2List(i,j));
                end
            end
            
        end
        
        function Fig = show_analyzed_fibril(obj)
            T = sprintf('Height Map of %s\nwith chosen indentation points',obj.Name);
            Fig = figure('Name',T,'Units','normalized','Position',[0.5 0.1 0.5 0.8]);
            
            subplot(2,2,1)
            I = imresize(obj.HeightMap(:,:,1).*1e9,[1024 1024]);
            %             I = (I*range(obj.HeightMap(:,:,1),'all') + min(obj.HeightMap(:,:,1),[],'all'))*1e9;
            imshow(I,[min(I,[],'all') max(I,[],'all')],'Colormap',hot)
            title(sprintf('%s Plane Fitted Height',obj.Name))
            c1 = colorbar;
            c1.Label.String = 'Height [nm]';
            hold on;
            try
                for i=1:obj.NumProfiles
                    plot((obj.List2Map(obj.RectApexIndex(i),2)-1/2)*1024/obj.NumPoints,...
                        (obj.List2Map(obj.RectApexIndex(i),1)-1/2)*1024/obj.NumProfiles,...
                        'g+', 'MarkerSize', 10, 'LineWidth', 2);
                    %                 plot((obj.List2Map(obj.ApexIndex(i),2)-1/2)*1024/obj.NumPoints,...
                    %                     (obj.List2Map(obj.ApexIndex(i),1)-1/2)*1024/obj.NumProfiles,...
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
            EM = imresize(mat2gray(log(obj.EModMapOliverPharr(:,:,1))),[1024 1024]);
            imshow(EM)
            title('E-Modulus Map')
        end
        
        function Fig = show_analyzed_fibril_micro(obj)
            T = sprintf('Height Map of %s\nwith chosen indentation points',obj.Name);
            Fig = figure('Name',T,'Units','normalized','Position',[0.5 0.1 0.5 0.8]);
            
            subplot(2,2,1)
            I = imresize(obj.HeightMap(:,:,1).*1e9,[1024 1024]);
            %             I = (I*range(obj.HeightMap(:,:,1),'all') + min(obj.HeightMap(:,:,1),[],'all'))*1e9;
            imshow(I,[min(I,[],'all') max(I,[],'all')],'Colormap',hot)
            title(sprintf('%s Plane Fitted Height',obj.Name))
            c1 = colorbar;
            c1.Label.String = 'Height [nm]';
            hold on;
            try
                for i=1:obj.NumProfiles
                    plot((obj.List2Map(obj.RectApexIndex(i),2)-1/2)*1024/obj.NumPoints,...
                        (obj.List2Map(obj.RectApexIndex(i),1)-1/2)*1024/obj.NumProfiles,...
                        'g+', 'MarkerSize', 10, 'LineWidth', 2);
                    %                 plot((obj.List2Map(obj.ApexIndex(i),2)-1/2)*1024/obj.NumPoints,...
                    %                     (obj.List2Map(obj.ApexIndex(i),1)-1/2)*1024/obj.NumProfiles,...
                    %                     'g+', 'MarkerSize', 10, 'LineWidth', 1);
                end
            catch
            end
            title(T);
            
            subplot(2,2,3)
            Masks = imresize(obj.FibMask.*obj.ExclMask,[1024 1024]);
            imshow(Masks)
            title('Fibril mask with excluded areas')
            
%             subplot(2,2,2)
%             surf(imrotate(obj.HeightMap',90),'LineStyle','none','FaceLighting','gouraud','FaceColor','interp')
%             axis manual
%             light('Style','local')
%             title(sprintf('Fibril Diameter = %.2fnm',obj.FibDiam*1e9))
%             
%             subplot(2,2,4)
%             EM = imresize(mat2gray(log(obj.EModMapOliverPharr(:,:,1))),[1024 1024]);
%             imshow(EM)
%             title('E-Modulus Map')
        end
        
        function show_height_map(obj)
            T = sprintf('Height Map of %s',obj.Name);
            Fig = figure('Name',T,'Units','normalized','Position',[0.5 0.1 0.5 0.8]);
            
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
        
        function show_stephold(obj)
            
            close all
            
            
            for i=1:obj.NCurves
                
                
                [MultiplierF,UnitF,~] = AFMImage.parse_unit_scale(range(obj.Force{i,2}),'N',10);
                [MultiplierF2,UnitF2,~] = AFMImage.parse_unit_scale(range(obj.BasedForce{i,2}),'N',10);
                [MultiplierH,UnitH,~] = AFMImage.parse_unit_scale(range(obj.Height{i,2}),'m',10);
                
                % find min/max of indentation and force modulation
                Hmin = min(obj.Height{i,2});
                Hmax = max(obj.Height{i,2});
                Fmin = min(obj.Force{i,2});
                Fmax = max(obj.Force{i,2});
                Fmin2 = min(obj.BasedForce{i,2});
                Fmax2 = max(obj.BasedForce{i,2});
                
                yHmin = min(Hmin)*MultiplierH;
                yHmax = max(Hmax)*MultiplierH;
                yFmin = min(Fmin)*MultiplierF;
                yFmax = max(Fmax)*MultiplierF;
                yFmin2 = min(Fmin2)*MultiplierF2;
                yFmax2 = max(Fmax2)*MultiplierF2;
                
                if yHmin > 0
                    yHmin = yHmin * 0.995;
                    
                else
                    yHmin = yHmin * 1.005;
                end
                
                if yHmax > 0
                    yHmax = yHmax * 1.005;
                    
                else
                    yHmax = yHmax * 0.995;
                end
                
                if yFmin > 0
                    yFmin = yFmin * 0.99;
                    
                else
                    yFmin = yFmin * 1.01;
                end
                
                if yFmax > 0
                    yFmax = yFmax * 1.01;
                    
                else
                    yFmax = yFmax * 0.99;
                end
                
                if yFmin2 > 0
                    yFmin2 = yFmin2 * 0.99;
                    
                else
                    yFmin2 = yFmin2 * 1.01;
                end
                
                if yFmax2 > 0
                    yFmax2 = yFmax2 * 1.01;
                    
                else
                    yFmax2 = yFmax2 * 0.99;
                end
                
                
                
                %Colours 
                lila = [0.368, 0.058, 0.721];
                lightblue = [0.101, 0.701, 0.976];
                darkblue = [0.109, 0.078, 0.941];
                
                
                for j=1:obj.NumSegments
                    
                       lengthApp = length(obj.App{i});
                       obj.SecPerPoint{1} = obj.SegDuration{1}/lengthApp;
                       obj.TStart{1} = obj.SecPerPoint{1}/2;
                       obj.TEnd{1} = obj.SeriesTime{1};
                       obj.SegTime{1} = obj.TStart{1}:obj.SecPerPoint{1}:obj.TEnd{1};
                       obj.SegTime{1} = obj.SegTime{1}.';
                       
                       lastseg = obj.NumSegments - 1;
                       lengthRet = length(obj.Ret{i});
                       obj.SecPerPoint{obj.NumSegments} = obj.SegDuration{obj.NumSegments}/lengthRet;
                       obj.TStart{obj.NumSegments} = obj.SeriesTime{lastseg}+(obj.SecPerPoint{obj.NumSegments}/2);
                       obj.TEnd{obj.NumSegments} = obj.SeriesTime{obj.NumSegments};
                       obj.SegTime{obj.NumSegments} = obj.TStart{obj.NumSegments}:obj.SecPerPoint{obj.NumSegments}:obj.TEnd{obj.NumSegments};
                       obj.SegTime{obj.NumSegments} = obj.SegTime{obj.NumSegments}.';
%                        
                end
                
                 % Plot
                figure('Name',sprintf('Stephold Curves %i',i))
                set(gcf,'units','normalized','outerposition',[0 0 1 1])
                hold on
                for j=1:obj.NumSegments
                    
                    subplot(2,2,1)
                    hold on
                    plot(obj.SegTime{j},obj.Height{i,j}*MultiplierH,'-','color',darkblue)
                    ylim([yHmin yHmax])
                    title(sprintf('Displacement Time Curve %i',i),'FontSize', 18)
                    xlabel('time in s','FontSize', 12)
                    ylabel(sprintf('Displacement [%s]',UnitH),'FontSize', 12)
                    grid on
                    grid minor
                    
                    subplot(2,2,2)
                    hold on
                    plot(obj.SegTime{j},obj.Force{i,j}*MultiplierF,'-r')
                    ylim([yFmin yFmax])
                    title(sprintf('Force Time Curve %i',i),'FontSize', 18)
                    xlabel('time in s','FontSize', 12)
                    ylabel(sprintf('vDeflection-Force [%s]',UnitF),'FontSize', 12)
                    grid on
                    grid minor
                    
                    subplot(2,2,3)
                    hold on
                    plot(obj.SegTime{j},obj.BasedForce{i,j}*MultiplierF2,'-m')
                    ylim([yFmin2 yFmax2])
                    title(sprintf('Based Force Time Curve %i',i),'FontSize', 18)
                    xlabel('time in s','FontSize', 12)
                    ylabel(sprintf('vDeflection-Force [%s]',UnitF2),'FontSize', 12)
                    grid on
                    grid minor
                    
                end
            end
        end
        
        function show_height(obj)
            
            close all
            
            DirectoryPath = uigetdir();
            
            for i=1:obj.NCurves
                
                frequencies = cell2mat(obj.SegFrequency);
                FirstFreq = find(frequencies,1,'first');
                
                [MultiplierF,UnitF,~] = AFMImage.parse_unit_scale(range(obj.Force{i,FirstFreq}),'N',10);
                [MultiplierF2,UnitF2,~] = AFMImage.parse_unit_scale(range(obj.BasedForce{i,FirstFreq}),'N',10);
                [MultiplierH,UnitH,~] = AFMImage.parse_unit_scale(range(obj.Height{i,FirstFreq}),'m',10);
                
                % find min/max of indentation and force modulation
                Hmin = zeros(obj.NumSegments,1);
                Hmax = zeros(obj.NumSegments,1);
                Fmin = zeros(obj.NumSegments,1);
                Fmax = zeros(obj.NumSegments,1);
                for j=1:obj.NumSegments
                    if obj.SegFrequency{j} > 0

                       Hmin(j,:) = min(obj.Height{i,j});
                       Hmax(j,:) = max(obj.Height{i,j});
                       Fmin(j,:) = min(obj.Force{i,j});
                       Fmax(j,:) = max(obj.Force{i,j});
                    else
                       Hmin(j,:) = NaN;
                       Hmax(j,:) = NaN;
                       Fmin(j,:) = NaN;
                       Fmax(j,:) = NaN;

                    end
                end
                
                yHmin = min(Hmin)*MultiplierH;
                yHmax = max(Hmax)*MultiplierH;
                yFmin = min(Fmin)*MultiplierF;
                yFmax = max(Fmax)*MultiplierF;
                
                if yHmin > 0
                    yHmin = yHmin * 0.99;
                    
                else
                    yHmin = yHmin * 1.005;
                end
                
                if yHmax > 0
                    yHmax = yHmax * 1.01;
                    
                else
                    yHmax = yHmax * 0.995;
                end
                
                if yFmin > 0
                    yFmin = yFmin * 0.7;
                    
                else
                    yFmin = yFmin * 1.2;
                end
                
                if yFmax > 0
                    yFmax = yFmax * 1.2;
                    
                else
                    yFmax = yFmax * 0.7;
                end
                
                
                
                %Colours 
                lila = [0.368, 0.058, 0.721];
                lightblue = [0.101, 0.701, 0.976];
                darkblue = [0.109, 0.078, 0.941];
                
                
                 for j=1:obj.NumSegments
                    
                       lengthApp = length(obj.App{i});
                       obj.SecPerPoint{1} = obj.SegDuration{1}/lengthApp;
                       obj.TStart{1} = obj.SecPerPoint{1}/2;
                       obj.TEnd{1} = obj.SeriesTime{1};
                       obj.SegTime{1} = obj.TStart{1}:obj.SecPerPoint{1}:obj.TEnd{1};
                       obj.SegTime{1} = obj.SegTime{1}.';
                       
                       lastseg = obj.NumSegments - 1;
                       lengthRet = length(obj.Ret{i});
                       obj.SecPerPoint{obj.NumSegments} = obj.SegDuration{obj.NumSegments}/lengthRet;
                       obj.TStart{obj.NumSegments} = obj.SeriesTime{lastseg}+(obj.SecPerPoint{obj.NumSegments}/2);
                       obj.TEnd{obj.NumSegments} = obj.SeriesTime{obj.NumSegments};
                       obj.SegTime{obj.NumSegments} = obj.TStart{obj.NumSegments}:obj.SecPerPoint{obj.NumSegments}:obj.TEnd{obj.NumSegments};
                       obj.SegTime{obj.NumSegments} = obj.SegTime{obj.NumSegments}.';
%                        
                 end
                 
                 % Plot
                figure('Name',sprintf('Microrheology Curves %i',i))
                set(gcf,'units','normalized','outerposition',[0 0 1 1])
                hold on
                for j=1:obj.NumSegments
                    
                    subplot(2,2,1)
                    hold on
                    plot(obj.SegTime{j},obj.Height{i,j}*MultiplierH,'-','color',darkblue)
                    ylim([yHmin yHmax])
                    title(sprintf('Displacement Time Curve %i',i),'FontSize', 18)
                    xlabel('time in s','FontSize', 12)
                    ylabel(sprintf('Displacement [%s]',UnitH),'FontSize', 12)
                    grid on
                    grid minor
                    
                    subplot(2,2,2)
                    hold on
                    plot(obj.SegTime{j},obj.Force{i,j}*MultiplierF,'-r')
                    ylim([yFmin yFmax])
                    title(sprintf('Force Time Curve %i',i),'FontSize', 18)
                    xlabel('time in s','FontSize', 12)
                    ylabel(sprintf('vDeflection-Force [%s]',UnitF),'FontSize', 12)
                    grid on
                    grid minor
                    
                    % subplot 4: force vs indentation
                   subplot(2,2,3)
                   plot(obj.Height{i,1}*MultiplierH,obj.Force{i,1}*MultiplierF,'-r',obj.Height{i,obj.NumSegments}*MultiplierH,obj.Force{i,obj.NumSegments}*MultiplierF,'-b')
                   hold on
                   plot(obj.Height{i,j}*MultiplierH,obj.Force{i,j}*MultiplierF,':m')
                   xlim([yHmin-100 yHmax+50])
                   title(sprintf('Force Displacement Curve %i',i),'FontSize', 18)
                   xlabel(sprintf('Displacement [%s]',UnitH),'FontSize', 12);
                   ylabel(sprintf('vDeflection-Force [%s]',UnitF),'FontSize', 12);
                   grid on
                   grid minor
                   
                   l1 = plot(nan, nan, '-r');
                   hold on
                   l2 = plot(nan, nan, '-b');
                   l3 = plot(nan, nan, ':m');
                   l1.LineWidth = 2;
                   l2.LineWidth = 2;
                   l3.LineWidth = 2;
                   legend([l1, l2, l3], {'approach', 'retract','modulation'}, 'FontSize', 8)
                   
                   
                   subplot(2,2,4)
                   plot(obj.Height{i,1}*MultiplierH,obj.BasedForce{i,1}*MultiplierF2,'-r',obj.Height{i,obj.NumSegments}*MultiplierH,obj.BasedForce{i,obj.NumSegments}*MultiplierF2,'-b')
                   hold on
                   plot(obj.Height{i,j}*MultiplierH,obj.BasedForce{i,j}*MultiplierF2,':m')
                   xlim([yHmin-100 yHmax+50])
                   title(sprintf('Based Force Displacement Curve %i',i),'FontSize', 18)
                   xlabel(sprintf('Displacement [%s]',UnitH),'FontSize', 12);
                   ylabel(sprintf('vDeflection-Force [%s]',UnitF2),'FontSize', 12);
                   grid on
                   grid minor
                   
                   l1 = plot(nan, nan, '-r');
                   hold on
                   l2 = plot(nan, nan, '-b');
                   l3 = plot(nan, nan, ':m');
                   l1.LineWidth = 2;
                   l2.LineWidth = 2;
                   l3.LineWidth = 2;
                   legend([l1, l2, l3], {'approach', 'retract','modulation'}, 'FontSize', 8)
                    
                   
                   k = obj.RectApexIndex;
                   if DirectoryPath~=0
                       whereToStore=fullfile(DirectoryPath,['height_curves_' num2str(i) '.svg']);
                       saveas(gcf, whereToStore);
                   end
               
                end
            end
        end
        
        function show_modulation(obj)
            
            close all
            
            DirectoryPath = uigetdir();
            
            for i=1:obj.NCurves
                
                %Identify position of first modulation for Multiplier later
                frequencies = zeros(obj.NumSegments,1);
                for j=1:obj.NumSegments
                    frequencies(j,:) = obj.SegFrequency{j};
                end
                FirstFreq = find(frequencies,1,'first');
                
                
                % Multiplier
                [MultiplierF,UnitF,~] = AFMImage.parse_unit_scale(range(obj.BasedForce{i,FirstFreq}),'N',10);
                [MultiplierI,UnitI,~] = AFMImage.parse_unit_scale(range(obj.Indentation{i,FirstFreq}),'m',10);
                
                
                % find min/max of indentation and force modulation
                Hmin = zeros(obj.NumSegments,1);
                Hmax = zeros(obj.NumSegments,1);
                Fmin = zeros(obj.NumSegments,1);
                Fmax = zeros(obj.NumSegments,1);
                for j=1:obj.NumSegments
                    if obj.SegFrequency{j} > 0

                       Hmin(j,:) = min(obj.Indentation{i,j});
                       Hmax(j,:) = max(obj.Indentation{i,j});
                       Fmin(j,:) = min(obj.BasedForce{i,j});
                       Fmax(j,:) = max(obj.BasedForce{i,j});
                    else
                       Hmin(j,:) = NaN;
                       Hmax(j,:) = NaN;
                       Fmin(j,:) = NaN;
                       Fmax(j,:) = NaN;

                    end
                end
                
                yHmin = min(Hmin)*MultiplierI;
                yHmax = max(Hmax)*MultiplierI;
                yFmin = min(Fmin)*MultiplierF;
                yFmax = max(Fmax)*MultiplierF;
                
                if yHmin > 0
                    yHmin = yHmin * 0.7;
                    
                else
                    yHmin = yHmin * 1.2;
                end
                
                if yHmax > 0
                    yHmax = yHmax * 1.2;
                    
                else
                    yHmax = yHmax * 0.7;
                end
                
                if yFmin > 0
                    yFmin = yFmin * 0.7;
                    
                else
                    yFmin = yFmin * 1.2;
                end
                
                if yFmax > 0
                    yFmax = yFmax * 1.2;
                    
                else
                    yFmax = yFmax * 0.7;
                end
                
                
                %Colours 
                lila = [0.368, 0.058, 0.721];
                lightblue = [0.101, 0.701, 0.976];
                darkblue = [0.109, 0.078, 0.941];
                
                for j=1:obj.NumSegments
                    
                       lengthApp = length(obj.App{i});
                       obj.SecPerPoint{1} = obj.SegDuration{1}/lengthApp;
                       obj.TStart{1} = obj.SecPerPoint{1}/2;
                       obj.TEnd{1} = obj.SeriesTime{1};
                       obj.SegTime{1} = obj.TStart{1}:obj.SecPerPoint{1}:obj.TEnd{1};
                       obj.SegTime{1} = obj.SegTime{1}.';
                       
                       lastseg = obj.NumSegments - 1;
                       lengthRet = length(obj.Ret{i});
                       obj.SecPerPoint{obj.NumSegments} = obj.SegDuration{obj.NumSegments}/lengthRet;
                       obj.TStart{obj.NumSegments} = obj.SeriesTime{lastseg}+(obj.SecPerPoint{obj.NumSegments}/2);
                       obj.TEnd{obj.NumSegments} = obj.SeriesTime{obj.NumSegments};
                       obj.SegTime{obj.NumSegments} = obj.TStart{obj.NumSegments}:obj.SecPerPoint{obj.NumSegments}:obj.TEnd{obj.NumSegments};
                       obj.SegTime{obj.NumSegments} = obj.SegTime{obj.NumSegments}.';
%                        
%                        lengthHHApp = length(obj.HHApp{i});
%                        obj.SecPerPoint{1} = obj.SegDuration{1}/lengthHHApp;
%                        obj.TStart{1} = obj.SecPerPoint{1}/2;
%                        obj.TEnd{1} = obj.SeriesTime{1};
%                        obj.SegTime{1} = obj.TStart{1}:obj.SecPerPoint{1}:obj.TEnd{1};
%                        obj.SegTime{1} = obj.SegTime{1}.';
%                        
%                        lengthHHRet = length(obj.HHRet{i});
%                        obj.SecPerPoint{obj.NumSegments} = obj.SegDuration{obj.NumSegments}/lengthHHRet;
%                        obj.TStart{obj.NumSegments} = obj.SeriesTime{obj.NumSegments - 1}+(obj.SecPerPoint{obj.NumSegments}/2);
%                        obj.TEnd{obj.NumSegments} = obj.SeriesTime{obj.NumSegments};
%                        obj.SegTime{obj.NumSegments} = obj.TStart{obj.NumSegments}:obj.SecPerPoint{obj.NumSegments}:obj.TEnd{obj.NumSegments};
%                        obj.SegTime{obj.NumSegments} = obj.SegTime{obj.NumSegments}.';
                       
                
                end
                
                % Plot
                figure('Name',sprintf('Microrheology Curves %i',i))
                set(gcf,'units','normalized','outerposition',[0 0 1 1])
                %set(gcf, 'Position', [100, 100, 800, 600]);
                hold on
                for j=1:obj.NumSegments
                    
                    %subplot 1: force time 
                    subplot(2,2,1)
                    hold on
                    plot(obj.SegTime{j},obj.BasedForce{i,j}*MultiplierF,'m')
                    title(sprintf('Force Time Curve %i',i),'FontSize', 18)
                    xlabel('time in s','FontSize', 12)
                    ylabel(sprintf('vDeflection-Force [%s]',UnitF),'FontSize', 12)
                    grid on
                    grid minor

                    %subplot 2: indentation time
                    subplot(2,2,2)
                    hold on
                    plot(obj.SegTime{j},obj.Indentation{i,j}*MultiplierI,'-','color',lightblue)
                    ylim([yHmin yHmax])
                    title(sprintf('Indentation Time Curve %i',i),'FontSize', 18)
                    xlabel('time in s','FontSize', 12)
                    ylabel(sprintf('Indentation [%s]',UnitI),'FontSize', 12);
                    grid on
                    grid minor

                    %subplot 3: force and indentation over time
                    subplot(2,2,3)
                    hold on

                    yyaxis left
                    ylim([yFmin yFmax])
                    plot(obj.SegTime{j},obj.BasedForce{i,j}*MultiplierF,'-m')
                    set(gca, 'YColor', 'm')
                    xlabel('time in s','FontSize', 12)
                    ylabel(sprintf('vDeflection-Force [%s]',UnitF),'FontSize', 12)


                    yyaxis right
                    ylim([yHmin yHmax])
                    plot(obj.SegTime{j},obj.Indentation{i,j}*MultiplierI,'-','color',lightblue)
                    set(gca, 'YColor', lightblue)
                    title(sprintf('Force and Indentation over Time Curve %i',i),'FontSize', 18)
                    ylabel(sprintf('Indentation [%s]',UnitI),'FontSize', 12);
                    grid on
                    grid minor


                    % subplot 4: force vs indentation
                    subplot(2,2,4)
                    plot(obj.Indentation{i,1}*MultiplierI,obj.BasedForce{i,1}*MultiplierF,'-r',obj.Indentation{i,obj.NumSegments}*MultiplierI,obj.BasedForce{i,obj.NumSegments}*MultiplierF,'-b')
                    hold on
                    plot(obj.Indentation{i,j}*MultiplierI,obj.BasedForce{i,j}*MultiplierF,':m')
                    xlim([yHmin-40 yHmax+20])
                    title(sprintf('Force Indentation Curve %i',i),'FontSize', 18)
                    xlabel(sprintf('Indentation [%s]',UnitI),'FontSize', 12);
                    ylabel(sprintf('vDeflection-Force [%s]',UnitF),'FontSize', 12);
                    grid on
                    grid minor
                   
                    l1 = plot(nan, nan, '-r');
                    hold on
                    l2 = plot(nan, nan, '-b');
                    l3 = plot(nan, nan, ':m');
                    l1.LineWidth = 2;
                    l2.LineWidth = 2;
                    l3.LineWidth = 2;
                    legend([l1, l2, l3], {'approach', 'retract','modulation'}, 'FontSize', 8)

               
                end
               
               k = obj.RectApexIndex;
               if DirectoryPath~=0
                   whereToStore=fullfile(DirectoryPath,['microrheology_curve_' num2str(i) '.svg']);
                   saveas(gcf, whereToStore);
               end
                       
          
                       
                
                %force indentation only segments
%                 s=2*obj.NCurves + 10*i;
%                 for j=1:obj.NumSegments
%                     
%                     
%                     if obj.SegFrequency{j} > 0
%                         
%                         Freq = obj.SegFrequency{j};
%                         
%                         figure(s)
%                         plot(obj.Indentation{i,j}*MultiplierI,obj.BasedForce{i,j}*MultiplierF,'b')
%                         title(sprintf('Force Indentation Curve %i, %.1f [Hz]',i,Freq),'FontSize', 18)
%                         xlabel(sprintf('Indentation [%s]',UnitI),'FontSize', 16);
%                         ylabel(sprintf('vDeflection-Force [%s]',UnitF),'FontSize', 16);
%                         grid on
%                         grid minor
%                         
%                         if DirectoryPath~=0
%                            whereToStore=fullfile(DirectoryPath,['force_indentation_curve_' num2str(i) '_segment_' num2str(j) '.svg']);
%                            saveas(gcf, whereToStore);
%                         end
%                        
%                         
%                         
%                     end
%                     s=s+1;
%                        
%                 end
                
                

            end
            
        end
        
        function show_microrheology_results(obj)
            % showing the results (phaseshift, losstangent, storage and
            % loss modulus) over the different frequencies
            
            close all
            DirectoryPath = uigetdir();
            k=1;
            
            figure('Name',sprintf('Results'))
            set(gcf,'units','normalized','outerposition',[0 0 1 1])
            hold on
            
            Dphi = zeros(obj.NumSegments,1);
            frequencies = zeros(obj.NumSegments,1);
            losstangent = zeros(obj.NumSegments,1);
            emodmicro1 = zeros(obj.NumSegments,1);
            emodmicro2 = zeros(obj.NumSegments,1);
            
            for j=1:obj.NumSegments
                    frequencies(j,:) = obj.SegFrequency{j};
            end
            frequencies = frequencies(frequencies ~= 0);
            lf = length(frequencies);
            
            Dphi = cell2mat(obj.DeltaPhi);
            losstangent = cell2mat(obj.LossTangent);
            emodmicro1 = obj.EModMicro1;
            emodmicro2 = obj.EModMicro2;
            emodmicro1(:,all(emodmicro1 == 0))=[];
            emodmicro2(:,all(emodmicro2 == 0))=[];
            %emodmicro1(isnan(emodmicro1))=0;
            %emodmicro2(isnan(emodmicro2))=0;
            
            
            indices = find(obj.DZslope > 0.95);
            
            Dphi(indices,:) = [];
            losstangent(indices,:) = [];
            emodmicro1(indices,:) = [];
            emodmicro2(indices,:) = [];
                
            %figure('Name',sprintf('Results'))
            hold on

            subplot(2,2,1)
            %plot(xq,vqDphi,'-','DisplayName',sprintf('Curve %i',i))
            %hold on
            boxplot(Dphi,frequencies)
            %ylim([-270 270])
            title('Phaseshift of all curves','FontSize', 18)
            xlabel('frequency [Hz]','FontSize', 16)
            ylabel('phaseshift [°]','FontSize', 16)
            %legend show
            %drawnow
            grid on
            grid minor

            subplot(2,2,2)
            %plot(xq,vqLosstangent,'-','DisplayName',sprintf('Curve %i',i))
            %hold on
            boxplot(losstangent,frequencies)
            title('Loss Tangent of all curves','FontSize', 18)
            xlabel('frequency [Hz]','FontSize', 16)
            ylabel('losstangent','FontSize', 16)
            %legend show
            %drawnow
            grid on
            grid minor

            subplot(2,2,3)
            %plot(xq,vqEmodmicro1,'-','DisplayName',sprintf('Curve %i',i))
            %,'HandleVisibility','off'
            %hold on
            boxplot(emodmicro1*1e-6, frequencies)
            title('Storage modulus','FontSize', 18)
            xlabel('frequency [Hz]','FontSize', 16)
            ylabel('elastic modulus [MPa]','FontSize', 16)
            %legend show
            %drawnow
            grid on
            grid minor

            subplot(2,2,4)
            %plot(xq,vqEmodmicro2,'-','DisplayName',sprintf('Curve %i',i))
            %hold on
            boxplot(emodmicro2*1e-6,frequencies)
            title('Loss modulus','FontSize', 18)
            xlabel('frequency [Hz]','FontSize', 16)
            ylabel('viscous modulus [MPa]','FontSize', 16)
            %legend show
            %drawnow
            grid on
            grid minor

            if DirectoryPath~=0
               whereToStore=fullfile(DirectoryPath,['results.svg']);
               saveas(gcf, whereToStore);
            end
            
 
        end
        
        function show_sine(obj)
            % this function is to display the filtered force and
            % indentation data and the respective fits to see the quality
            % of the fit
            
            close all
            DirectoryPath = uigetdir();
            k=1;
            
            
            for i=1:obj.NCurves
                
                %Identify position of first modulation for Multiplier later
                 frequencies = zeros(obj.NumSegments,1);
                 time = obj.SegTime;
                 for j=1:obj.NumSegments
                     frequencies(j,:) = obj.SegFrequency{j};
                     
                 end
                %frequencies = cell2mat(obj.SegFrequency);
                FF = obj.FilterF;
                FH = obj.FilterH;
                SVF = obj.SineVarsF;
                SVH = obj.SineVarsH;
                BF = obj.BasedForce;
                BH = obj.Indentation;
                empty_fields = cellfun(@isempty,FF);
                FF(:,empty_fields(1,:))=[];
                FH(:,empty_fields(1,:))=[];
                SVF(:,empty_fields(1,:))=[];
                SVH(:,empty_fields(1,:))=[];

                time = time(frequencies ~= 0);
                BF = BF(:,frequencies ~= 0);
                BH = BH(:,frequencies ~= 0);
                frequencies = frequencies(frequencies ~= 0);
                lf = length(frequencies);
                
                maxfreq = max(frequencies);
                minfreq = min(frequencies);
                faktor = maxfreq/minfreq;
                
%                 if faktor > 10
%                             
%                 end
                
                %Colours 
                lila = [0.368, 0.058, 0.721];
                lightblue = [0.101, 0.701, 0.976];
                darkblue = [0.109, 0.078, 0.941];
                
                
                %Plot
                figure('Name',sprintf('Normalized curves with Fit %i',i))
                set(gcf,'units','normalized','outerposition',[0 0 1 1])
                
                %t = tiledlayout(1,lf);
                %t.TileSpacing = 'none';
                
                %n = 1*lf;
                %h = gobjects(2,n);
                
                hold on
                
                %nexttile
                %semilogx(time{1},FF{i,1},'-m')
                ax1 = subplot(121);
                x = time(1);
                x =cell2mat(x);
                e = 1.05*time{1}(end);
                a = 0.95*time{1}(1);
                %x = exp(x);
                semilogx(x,FF{i,1},'-m');

                hold on
                ax2 = subplot(122);
                for j=2:lf
                    

                        % Divide data through their range
                        rangeF = range(BF{i,j});
                        rangeH = range(BH{i,j});
                        
                        %obj.SineVarsF{i,j}(1) = obj.SineVarsF{i,j}(1)/rangeF;
                        %obj.SineVarsH{i,j}(1) = obj.SineVarsH{i,j}(1)/rangeH;
                        %obj.slopeF = obj.slopeF/rangeF;
                        %obj.slopeH = obj.slopeH/rangeH;
                        %obj.interceptF = obj.interceptF/rangeF;
                        %obj.interceptH = obj.interceptH/rangeH;
                        x = time(j);
                        x = cell2mat(x);
                        
                        x1 = 0.95*time{2}(1);
                        x2 = 1.15*time{2}(end);
                        
                        
                        
                        %Y-values fitted sine of indentation and force:
                         ypF = SVF{i,j}(1)*(sin(2*pi*x.*SVF{i,j}(2) + SVF{i,j}(3)));
                         ypH = SVH{i,j}(1)*(sin(2*pi*x.*SVH{i,j}(2) + SVH{i,j}(3)));
                           
                        

                        
                            
                            
                        hold on


                        
                        %subplot(2,1,1)
                        %nexttile
                        plot(x,FF{i,j}*rangeF,'-m');
                        hold on
                        %[MultiplierF,UnitF,~] = AFMImage.parse_unit_scale(range(obj.FilterF{i,FirstFreq}),'N',10);
%                         ax2 = subplot(122);
%                         for j=2:lf
%                             plot(time{j},FF{i,j},'-m');
%                         end
                        %hold on
                        %plot(x,ypF,'-','color',lila)
                        %title(sprintf('Normalized Force over Time incl. Fit Curve %i',i),'FontSize', 18)
                        %xlabel('time [s]','FontSize', 16)
                        %ylabel(sprintf('vDeflection-Force'),'FontSize', 16)
                        %xticks(10:0.1:150)
                        %xticklabels({'10','100','150'})
                        %set(gca,'xscale','log')
                        %grid on
                        %grid minor
%                         
%                         subplot(2,1,2)
%                         %[MultiplierI,UnitI,~] = AFMImage.parse_unit_scale(range(obj.FilterH{i,FirstFreq}),'m',10);
%                         semilogx(x,FH{i,j},'-','color', lightblue)
%                         hold on
%                         semilogx(x,ypH,'-','color',darkblue)
%                         title(sprintf('Normalized Indentation over Time incl. Fit Curve %i',i),'FontSize', 18)
%                         xlabel('time [s]','FontSize', 16)
%                         ylabel(sprintf('Indentation'),'FontSize', 16);
%                         grid on
%                         grid minor
                        
                        
%                         l1 = plot(nan, nan, 'm-');
%                         hold on
%                         l2 = plot(nan, nan, '-','color', lila);
%                         l3 = plot(nan, nan, '-', 'color', lightblue);
%                         l4 = plot(nan, nan, '-','color',darkblue);
%                         l1.LineWidth = 3;
%                         l2.LineWidth = 3;
%                         l3.LineWidth = 3;
%                         l4.LineWidth = 3;
%                         legend([l1, l2, l3, l4], {'normalized force data', 'force fit','normalized indentation data', 'indentation fit'}, 'Location', 'southoutside','FontSize', 14)

        %                         obj.SineVarsF{i,j}(1) = obj.SineVarsF{i,j}(1)*rangeF;
        %                         obj.SineVarsH{i,j}(1) = obj.SineVarsH{i,j}(1)*rangeH;
        %                         obj.slopeF = obj.slopeF*rangeF;
        %                         obj.slopeH = obj.slopeH*rangeH;
        %                         obj.interceptF = obj.interceptF*rangeF;
        %                         obj.interceptH = obj.interceptH*rangeH;


                        
                      
                        
                end
                
                set(ax1,'units','normalized','position',[0.1 0.4 0.3 0.3]);
                set(ax2,'units','normalized','position',[0.4 0.4 0.7 0.3]);
                set(ax1,'xscale','log','xlim',[a e]);
                set(ax2,'xlim',[x1 x2],'ytick',[],'yticklabel','');

                set([ax1 ax2],'box','off');
                set(ax2,'yticklabel','');
                uistack(ax2,'top');
                title(sprintf('Normalized Force over Time incl. Fit Curve %i',i),'FontSize', 18)
                xlabel('time [s]','FontSize', 16)
                ylabel(sprintf('vDeflection-Force'),'FontSize', 16)
                grid(ax1,'on');
                grid(ax2,'on');
                
                k = obj.RectApexIndex;
                if DirectoryPath~=0
                   whereToStore=fullfile(DirectoryPath,['filtered_force_indentation_fit_curve_' num2str(i) '.svg']);
                   saveas(gcf, whereToStore);
                end
            end
        end
        
        function show_sine_2(obj)
            % showing the original and raw force and indentation data with
            % the respective fits
            
            close all
            DirectoryPath = uigetdir();
            k=1;
            
            
            for i=1:obj.NCurves
                
                %Identify position of first modulation for Multiplier later
                frequencies = zeros(obj.NumSegments,1);
                for j=1:obj.NumSegments
                    frequencies(j,:) = obj.SegFrequency{j};
                end
                FirstFreq = find(frequencies,1,'first');
                
                %Colours 
                lila = [0.368, 0.058, 0.721];
                lightblue = [0.101, 0.701, 0.976];
                darkblue = [0.109, 0.078, 0.941];
                
                
                %Plot
                figure('Name',sprintf('Raw curves with Fit %i',i))
                set(gcf,'units','normalized','outerposition',[0 0 1 1])
                lastseg = obj.NumSegments - 2;
                hold on
                for j=FirstFreq:lastseg
                    
                    
                    %if obj.SegFrequency{j} > 0
                        
                        x= obj.SegTime{j};
                        
                        
                        % Estimate offset
                         meanF = mean(obj.BasedForce{i,j});
                         meanH = mean(obj.Indentation{i,j});
                        
                        %Y-values fitted sine of indentation and force:
                        try
                            ypF = obj.SineVarsF{i,j}(1)*(sin(2*pi*x.*obj.SineVarsF{i,j}(2) + obj.SineVarsF{i,j}(3)));
                            ypH = obj.SineVarsH{i,j}(1)*(sin(2*pi*x.*obj.SineVarsH{i,j}(2) + obj.SineVarsH{i,j}(3)));
                            
                            ypFtrend = ypF + obj.slopeF(i,j).*x + obj.interceptF(i,j);
                            ypHtrend = ypH + obj.slopeH(i,j).*x + obj.interceptH(i,j);
                        catch
                            ypF = zeros(length(x),1);
                            ypH = zeros(length(x),1);
                            
                            ypFtrend = zeros(length(x),1);
                            ypHtrend = zeros(length(x),1);
                        end
                        
                        % calculate linear fit
                        lF = obj.slopeF(i,j).*x;
                        lH = obj.slopeH(i,j).*x;
                        

                        
                        hold on

                        yyaxis left
                        [MultiplierF,UnitF,~] = AFMImage.parse_unit_scale(range(obj.BasedForce{i,FirstFreq}),'N',10);
                        semilogx(x,obj.BasedForce{i,j}*MultiplierF,'-m')
                        hold on
                        if obj.SegFrequency{j} > 0
                            semilogx(x,ypFtrend*MultiplierF,'-','color',lila)
                        end
                        set(gca, 'YColor', 'm')
                        xlabel('time [s]','FontSize', 16)
                        ylabel(sprintf('vDeflection-Force [%s]',UnitF),'FontSize', 16)
                        
                        
                        yyaxis right
                        [MultiplierI,UnitI,~] = AFMImage.parse_unit_scale(range(obj.Indentation{i,FirstFreq}),'m',10);
                        semilogx(x,obj.Indentation{i,j}*MultiplierI,'-','color', lightblue)
                        hold on
                        if obj.SegFrequency{j} > 0
                            semilogx(x,ypHtrend*MultiplierI,'-','color',darkblue)
                        end
                        set(gca, 'YColor', lightblue)
                        title(sprintf('Force and Indentation over Time incl. Fit Curve %i',i),'FontSize', 18)
                        ylabel(sprintf('Indentation [%s]',UnitI),'FontSize', 16);
                        grid on
                        grid minor


                        l1 = plot(nan, nan, 'm-');
                        hold on
                        l2 = plot(nan, nan, '-','color', lila);
                        l3 = plot(nan, nan, '-', 'color', lightblue);
                        l4 = plot(nan, nan, '-','color',darkblue);
                        l1.LineWidth = 3;
                        l2.LineWidth = 3;
                        l3.LineWidth = 3;
                        l4.LineWidth = 3;
                        legend([l1, l2, l3, l4], {'raw force data', 'force fit','raw indentation data', 'indentation fit'}, 'Location', 'southoutside','FontSize', 14)
                       
                        k = obj.RectApexIndex;
                        if DirectoryPath~=0
                           whereToStore=fullfile(DirectoryPath,['raw_force_indentation_fit_curve_' num2str(i) '.svg']);
                           saveas(gcf, whereToStore);
                       end
                        
                        
                       
                    %end
                        
                end
            end
 
        end
        
        function show_phaseshift(obj)
            close all
            DirectoryPath = uigetdir();
            k=1;
            for i=1:obj.NCurves
                
                %Identify position of first modulation for Multiplier later
                frequencies = zeros(obj.NumSegments,1);
                for j=1:obj.NumSegments
                    frequencies(j,:) = obj.SegFrequency{j};
                end
                FirstFreq = find(frequencies,1,'first');
                frequencies = frequencies(frequencies ~= 0);
                lf = length(frequencies);
                k = obj.NumSegments/2;
                
                        
               % Plot
                figure('Name',sprintf('Phaseshift Curve %i',i))
                set(gcf,'units','normalized','outerposition',[0 0 1 1])
                hold on
                for j=1:obj.NumSegments
                    
                     if obj.SegFrequency{j} > 0
                         
                        Freq = obj.SegFrequency{j};
                        
                        %Period = 2*pi/obj.SegFrequency{j};
                        End = obj.TStart{j} + 2/obj.SegFrequency{j};
                        x = linspace(obj.TStart{j},End);
                        %End = 2/obj.SegFrequency{j};
                        %x = linspace(0,End);
                        
                       
                        
                        %Y-values fitted sine of indentation and force:
                        try
                            ypF = obj.SineVarsF{i,j}(1)*(sin(2*pi*x.*obj.SineVarsF{i,j}(2) + obj.SineVarsF{i,j}(3)));
                            ypH = obj.SineVarsH{i,j}(1)*(sin(2*pi*x.*obj.SineVarsH{i,j}(2) + obj.SineVarsH{i,j}(3)));
                        catch
                            ypF = zeros(length(x),1);
                            ypH = zeros(length(x),1);
                        end
                        
                        

                        hold on
                        
                        [MultiplierI,UnitI,~] = AFMImage.parse_unit_scale(range(obj.Indentation{i,FirstFreq}),'m',10);
                        [MultiplierF,UnitF,~] = AFMImage.parse_unit_scale(range(obj.BasedForce{i,FirstFreq}),'N',10);
                        subplot(1,obj.NumSegments,j)
                        findpeaks(ypF*MultiplierF,x)
                        hold on
                        findpeaks(ypH*MultiplierI,x)
                        title(sprintf('Curve %i, %.1f [Hz]',i,Freq),'FontSize', 16)
                        xlabel('time [s]','FontSize', 14)
                        ylabel(sprintf('fit of vDeflection-Force [%s] and Indentation [%s]',UnitF, UnitI),'FontSize', 14)
                        legend({'force fit','force fit peak','indentation fit','indentation fit peak'},'Location','southoutside','FontSize', 12)
                        %drawnow
                        grid on
                        grid minor
                        %Nx = 2; % in x-Richtung verdoppeln
                        %Ny = 1; % in y-Richtung unverändert
                        %pos = get(gcf,'Position');
                        %set(gcf,'Position',[pos(1:2) Nx*pos(3) Ny*pos(4)]); 


                        k = obj.RectApexIndex;
                        if DirectoryPath~=0
                           whereToStore=fullfile(DirectoryPath,['force_indentation_fit_curve_' num2str(i) '.svg']);
                           saveas(gcf, whereToStore);
                        end
                       
                     end
                        
                end
            end
        end
        
        function quality_control_oliver_pharr_fibril(obj,PauseTime)
            % shows some relevant plots for the E-Mod calculation
            % rectified apex force curves
            
            if nargin < 2
                PauseTime = 1.5;
            end
            
            T = sprintf('Quality control of fibril %s',obj.Name);
            Fig = figure('Name',T,'Units','normalized','Position',[0.1 0.1 0.8 0.8]);
            m = 1;
            while 1==1
                try
                    figure(Fig);
                catch
                    return
                end
                k = obj.RectApexIndex(m);
                subplot(2,3,1)
                I = imresize(mat2gray(obj.HeightMap(:,:,1)),[1024 1024]);
                imshow(I);
                hold on;
                for i=1:obj.NumProfiles
                    if obj.RectApexIndex(i)==k
                        plot((obj.List2Map(obj.RectApexIndex(i),2)-1/2)*1024/obj.NumPoints,...
                            (obj.List2Map(obj.RectApexIndex(i),1)-1/2)*1024/obj.NumProfiles,...
                            'g*', 'MarkerSize', 10, 'LineWidth', 2);
                    else
                        plot((obj.List2Map(obj.RectApexIndex(i),2)-1/2)*1024/obj.NumPoints,...
                            (obj.List2Map(obj.RectApexIndex(i),1)-1/2)*1024/obj.NumProfiles,...
                            'r+', 'MarkerSize', 5, 'LineWidth', 1);
                    end
                end
                title('Height Map with Apex Points');
                
                subplot(2,3,2)
                plot(obj.HHApp{k},obj.BasedApp{k}/obj.SpringConstant,...
                    obj.HHRet{k},obj.BasedRet{k}/obj.SpringConstant)
                xlim([min(obj.HHApp{k})+range(obj.HHApp{k})/2 ...
                    max(obj.HHApp{k})+range(obj.HHApp{k})*0.1])
                title(sprintf('Elastic Modulus = %.2f MPa',obj.EModOliverPharr(k)*1e-6))
                legend('Approach','Retract','Location','northwest')
                xlabel('Head Height [m]')
                ylabel('vDeflection [m]')
                drawpoint('Position',[obj.CP(k,1) obj.CP(k,2)]);
                
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
                plot(0:obj.NumProfiles+1,obj.RefSlope*ones(obj.NumProfiles+2,1))
                ylim([0 1.3])
                xlim([0 obj.NumProfiles+1])
                hold on
                plot(1:obj.NumProfiles,obj.DZslope(obj.RectApexIndex),'bO')
                plot(m,obj.DZslope(obj.RectApexIndex(m)),'rO','MarkerFaceColor','r')
                xlabel('Index')
                ylabel('DZ-Slope')
                legend(sprintf('Glass Reference Slope = %.3f',obj.RefSlope))
                hold off
                
                subplot(2,3,5)
                plot(obj.IndDepth(obj.RectApexIndex)*1e9,...
                    obj.EModOliverPharr(obj.RectApexIndex)*1e-6,'bO')
                hold on
                plot(obj.IndDepth(obj.RectApexIndex(m))*1e9,...
                    obj.EModOliverPharr(obj.RectApexIndex(m))*1e-6,'rO','MarkerFaceColor','r')
                xlabel('Indentation Depth [nm]')
                ylabel('Elastic Modulus [MPa]')
                hold off
                
                subplot(2,3,6)
                plot(obj.IndDepth(obj.RectApexIndex)*1e9,...
                    obj.IndentArea(obj.RectApexIndex),'bO','MarkerSize',10,'MarkerFaceColor','b')
                hold on
                plot(obj.IndDepth(obj.RectApexIndex(m))*1e9,...
                    obj.IndentArea(obj.RectApexIndex(m)),'rO','MarkerSize',10,'MarkerFaceColor','r')
                Xmax = round(max(obj.IndDepth(obj.RectApexIndex))*1e9+5);
                xlim([0 Xmax])
                plot(1:Xmax,obj.ProjTipArea(1:Xmax),'Color','black')
                xlabel('Indentation Depth [nm]')
                ylabel('Projected Area [m^2]')
                hold off
                
                pause(PauseTime)
                if m<obj.NumProfiles
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
            Fig = figure('Name',T,'Units','normalized','Position',[0.1 0.1 0.8 0.8]);
            m = 1;
            while 1==1
                try
                    figure(Fig);
                catch
                    return
                end
                k = obj.RectApexIndex(m);
                
                subplot(2,2,1)
                I = imresize(mat2gray(obj.HeightMap(:,:,1)),[1024 1024]);
                imshow(I);
                hold on;
                for i=1:obj.NumProfiles
                    if obj.RectApexIndex(i)==k
                        plot((obj.List2Map(obj.RectApexIndex(i),2)-1/2)*1024/obj.NumPoints,...
                            (obj.List2Map(obj.RectApexIndex(i),1)-1/2)*1024/obj.NumProfiles,...
                            'g*', 'MarkerSize', 10, 'LineWidth', 2);
                    else
                        plot((obj.List2Map(obj.RectApexIndex(i),2)-1/2)*1024/obj.NumPoints,...
                            (obj.List2Map(obj.RectApexIndex(i),1)-1/2)*1024/obj.NumProfiles,...
                            'r+', 'MarkerSize', 5, 'LineWidth', 1);
                    end
                end
                title('Height Map with Apex Points');
                
                subplot(2,2,2)
                
                % Determine X-Range for HertzModel
                X = obj.THApp{m} - obj.CP(m,1);
                X(X<0) = [];
                HertzModelX = 0:range(X)/100:2*max(X);
                HertzModelY = feval(obj.HertzFit{m},HertzModelX);
                
                plot(HertzModelX(HertzModelY<=max(obj.BasedApp{m} - obj.CP(m,2))),HertzModelY(HertzModelY<=max(obj.BasedApp{m} - obj.CP(m,2))),...
                    obj.THApp{m} - obj.CP(m,1),obj.BasedApp{m} - obj.CP(m,2),...
                    obj.THRet{m} - obj.CP(m,1),obj.BasedRet{m} - obj.CP(m,2))
                xlim([min(obj.THApp{k} - obj.CP(k,1))+range(obj.THApp{k} - obj.CP(k,1))/2 ...
                    max(obj.THApp{k} - obj.CP(k,1))+range(obj.THApp{k} - obj.CP(k,1))*0.1])
                title(sprintf('Elastic Modulus = %.2f MPa',obj.EModHertz(k)*1e-6))
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
                plot((obj.IndDepthHertz(obj.RectApexIndex))*1e9,...
                    obj.EModHertz(obj.RectApexIndex)*1e-6,'bO')
                hold on
                plot((obj.IndDepthHertz(obj.RectApexIndex(m)))*1e9,...
                    obj.EModHertz(obj.RectApexIndex(m))*1e-6,'rO','MarkerFaceColor','r')
                xlabel('Indentation Depth [nm]')
                ylabel('Elastic Modulus [MPa]')
                hold off
                
                pause(PauseTime)
                if m<obj.NumProfiles
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
            Fig = figure('Name',T,'Units','normalized','Position',[0.1 0.1 0.8 0.8]);
            m = 1;
            while 1==1
                try
                    figure(Fig);
                catch
                    return
                end
                subplot(2,3,1)
                I = imresize(mat2gray(obj.HeightMap(:,:,1)),[1024 1024]);
                imshow(I);
                hold on;
                
                plot((obj.List2Map(m,2)-1/2)*1024/obj.NumPoints,...
                    (obj.List2Map(m,1)-1/2)*1024/obj.NumProfiles,...
                    'g*', 'MarkerSize', 10, 'LineWidth', 2);
                
                title('Height Map');
                
                subplot(2,3,2)
                plot(obj.HHApp{m},obj.BasedApp{m}/obj.SpringConstant,...
                    obj.HHRet{m},obj.BasedRet{m}/obj.SpringConstant)
                xlim([min(obj.HHApp{m})+range(obj.HHApp{m})/2 ...
                    max(obj.HHApp{m})+range(obj.HHApp{m})*0.1])
                title(sprintf('Elastic Modulus = %.2f MPa',obj.EModOliverPharr(m)*1e-6))
                legend('Approach','Retract','Location','northwest')
                xlabel('Head Height [m]')
                ylabel('vDeflection [m]')
                drawpoint('Position',[obj.CP(m,1) obj.CP(m,2)]);
                
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
                plot(obj.IndDepth*1e9,...
                    obj.EModOliverPharr*1e-6,'bO')
                hold on
                plot(obj.IndDepth(m)*1e9,...
                    obj.EModOliverPharr(m)*1e-6,'rO','MarkerFaceColor','r')
                xlabel('Indentation Depth [nm]')
                ylabel('Elastic Modulus [MPa]')
                hold off
                
                subplot(2,3,6)
                plot(obj.IndDepth*1e9,...
                    obj.IndentArea,'bO','MarkerSize',10,'MarkerFaceColor','b')
                hold on
                plot(obj.IndDepth(m)*1e9,...
                    obj.IndentArea(m),'rO','MarkerSize',10,'MarkerFaceColor','r')
                Xmax = round(max(obj.IndDepth)*1e9+5);
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
            Fig = figure('Name',T,'Units','normalized','Position',[0.1 0.1 0.8 0.8]);
            m = 1;
            while 1==1
                try
                    figure(Fig);
                catch
                    return
                end
                subplot(2,3,1)
                I = imresize(mat2gray(obj.HeightMap(:,:,1)),[1024 1024]);
                imshow(I);
                hold on;
                
                plot((obj.List2Map(m,2)-1/2)*1024/obj.NumPoints,...
                    (obj.List2Map(m,1)-1/2)*1024/obj.NumProfiles,...
                    'g*', 'MarkerSize', 10, 'LineWidth', 2);
                
                title('Height Map');
                
                subplot(2,2,2)
                
                % Determine X-Range for HertzModel
                X = obj.THApp{m} - obj.CP(m,1);
                X(X<0) = [];
                HertzModelX = 0:range(X)/100:2*max(X);
                FitModel = obj.HertzFit{m};
                FitModel.b = 0;
                HertzModelY = feval(FitModel,HertzModelX);
                
                plot(HertzModelX(HertzModelY<=max(obj.BasedApp{m} - obj.CP(m,2))),HertzModelY(HertzModelY<=max(obj.BasedApp{m} - obj.CP(m,2))),...
                    obj.THApp{m} - obj.CP(m,1),obj.BasedApp{m} - obj.CP(m,2),...
                    obj.THRet{m} - obj.CP(m,1),obj.BasedRet{m} - obj.CP(m,2))
                xlim([min(obj.THApp{m} - obj.CP(m,1))+range(obj.THApp{m} - obj.CP(m,1))/2 ...
                    max(obj.THApp{m} - obj.CP(m,1))+range(obj.THApp{m} - obj.CP(m,1))*0.1])
                title(sprintf('Elastic Modulus = %.2f MPa',obj.EModHertz(m)*1e-6))
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
                    plot((obj.IndDepthHertz)*1e9,...
                        obj.EModHertz*1e-6,'bO')
                    hold on
                    plot((obj.IndDepthHertz(m))*1e9,...
                        obj.EModHertz(m)*1e-6,'rO','MarkerFaceColor','r')
                    xlabel('Indentation Depth [nm]')
                    ylabel('Elastic Modulus [MPa]')
                else
                    histogram((obj.IndDepthHertz)*1e9)
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
            plot(1:obj.NumProfiles,obj.EModOliverPharr(obj.RectApexIndex)*1e-6,'bO')
            hold on
            plot(1:obj.NumProfiles,obj.EModHertz(obj.RectApexIndex)*1e-6,'rO')
            xlabel('Index')
            ylabel('Elastic Modulus [MPa]')
        end
        
    end
    
end