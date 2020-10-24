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
    % e.g. FM = ForceMap();
    % you should from now on call methods on this specific instance.
    % A valid call could for example be:
    % FM.base_and_tilt();
    % to conduct an operation (a base fit in this case) on the force map or
    % FM.TipRadius;
    % to get a class parameter of this force map (the tip radius of the used cantilever)
    
    properties
        % Properties shared for the whole Force Map
        
        Name            % name of the force map. taken as the name of the folder, containing the .csv files
        Folder          % location of the .csv files of the force map
        Header          % header properties taken from the JPK force map container file
        NCurves         % number of curves on the force map
        NumProfiles
        NumPoints
        Sensitivity
        SpringConstant
        DBanding        % Fourieranalysis-based estimate of DBanding perdiod (only available with sufficient resolution)
        RefSlope
        PixApp          % maximum number of measured points during approach
        PixRet          % maximum number of measured points during retraction
        SelectedCurves  % logical vector of length NCurves with 0s for excluded and 1s for included curves. gets initialized with ones
        TipRadius = 8  % (nominal, if not otherwise calculated) tip radius in nm for the chosen type of tip model
        PoissonR = 0.5  % standard Poisson ratio for most mechanical models
    end
    properties
       % Curve data Properties 
        
        App = {}        % approach force data in Newton
        Ret = {}        % retraction force data in Newton
        HHApp = {}      % head-height approach data in meters
        HHRet = {}      % head-height retract data in meters
        THApp = {}      % vertical tip height approach data in meters
        THRet = {}      % vertical tip height retract data in meters
        BasedApp = {}   % approach force data with subtracted base line and tilt in Newton
        BasedRet = {}   % retraction force data with subtracted base line and tilt in Newton
    end
    properties
        % Properties related to Contact Point (CP) estimation
        
        RoV = {}        % ratio of variance curve for CP estiamation
        GoF = {}        % goodness of fit curve for each selected curve
        CPComboCurve    % combination of the various metrics for contact point estimation
        CP              % chosen contact point (CP) for every selected curve, which is then used in E module fitting
        CP_RoV          % CP estimated with the ratio of variance method
        CP_GoF          % CP estimated with the goodness of fit method
        CP_Combo        % CP estimated with a combination of the GoF and RoV methods
        DeltaE = {}     %
        CP_MonteCarlo   % Contains all Monte Carlo Predictions of the CP
        YDropPred       % Contains the Dropoutpredictions for every curve in the forcemap
        CP_CNN_Error    % Contains the estimated uncertainty of the CP_CNN
        Man_CP          % manually chosen contact point
        CP_old          % contact point estimation from old script 'A_nIAFM_analysis_main'
        LoadOld         % comes from same script as CP_old
        UnloadOld       % comes from same script as CP_old
        CP_OliverPharr  % CP needed in Oliver Pharr type analysis, done on HeadHeight vs Displacement curves
        CP_OliverPharr_MonteCarlo
        CP_OliverPharr_MonteCarlo_STD
        
    end
    properties
       % Properties related to topological calculations, such as mapping and masking and visualisation 
       
        HeightMap       % height profile map taken from the maximum head-height from approach max(hhapp)
        EModMapHertz    % E modulus profile map. same ordering as HeightMap
        EModMapOliverPharr % """"
        List2Map        % An R->RxR ((k)->(i,j)) mapping of indices to switch between the two representations
        Map2List        % An RxR->R ((i,j)->(k))mapping of indices to switch between the two representations
        FibDiam = []    % Estimated fibril diameter
        FibMask         % Logical mask marking the whole fibril
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
        HertzFit        % HertzFit model generated in the calculate_e_mod_hertz method
        IndDepthHertz
        DZslope
        Stiffness 
        IndDepth
        IndentArea
        ProjTipArea
    end
    properties
        % auxiliary properties to facilitate comparing different methods of
        % CP estimation
        CP_OliverPharr_CNN
        CP_OliverPharr_Old
        CP_OliverPharr_RoV
        EModOliverPharr_CNN
        EModOliverPharr_Old
        EModOliverPharr_RoV
        EModHertz_CNN
        EModHertz_Old
        EModHertz_RoV
    end
    methods
        % Main methods of the class

        %%%%% this constructor method version goes together with Martin
        %%%%% Handelhausers .py-script for .jpk-force-map -> .cvs
        %%%%% conversion
%         function obj = ForceMap(mapfilepath,mapname)
%             %%% Constructor of the class
%             
%             % Specify the folder where the files live. And import them.
%             % Also get curent folder and return to it after import of
%             % files.
%             % Assigns the properties that can be found in the jpk-file
%             % already
%             
%             current = what();
%             
%             % determine if ForceMap is given a loadpath for existing .mat
%             if nargin > 0
%                 cd(mapfilepath);
%                 file = dir(sprintf('%s.mat',mapname));
%                 load(file.name);
%                 cd(current.path);
%                 msg = sprintf('loading %s',obj.Name);
%                 disp(msg);
%                 disp('loading successfull')
%                 return
%             end
%             
%             quest = 'Do you want to load an already existing .mat file of the force map?';
%             answer = questdlg(quest,'Load map...','...from .mat-file','...from .cvs-folder','...from .cvs-folder');
%             if isequal(answer, '...from .mat-file')
%                 [mapfile, mapfilepath] = uigetfile();
%                 cd(mapfilepath);
%                 load(mapfile,'-mat');
%                 cd(current.path);
%                 msg = sprintf('loading %s',obj.Name);
%                 disp(msg);
%                 disp('loading successfull')
%                 obj.Folder = mapfilepath;
%                 return
%             else
%             end
%             
%             
%             obj.Folder = uigetdir;
%             cd(obj.Folder);
%             splitfolder = strsplit(obj.Folder,'\');
%             obj.Name = string(splitfolder(end));
%             msg = sprintf('loading %s',obj.Name);
%             disp(msg);
%             
%             % Get a list of all files in the folder with the desired file name pattern.
%             % Change to whatever pattern you used for the names of your force map files.
%             filePattern = fullfile(obj.Folder, '*.csv');
%             theFiles = dir(filePattern);
%             % Check to make sure files of specified patterns have been found
%             if length(theFiles) < 1
%                 errorMessage = 'No files of the specified name pattern have been found in the current directory';
%                 uiwait(warndlg(errorMessage));
%                 return;
%             end
%             
%             % Import data from csv-files
%             import_header = importdata(theFiles(2).name);
%             k = 1;
%             for i=1:length(import_header)
%                 if isempty(strfind(import_header{i},'start-option.'))
%                     obj.Header{k,1} = import_header{i};
%                     k = k + 1;
%                 end
%             end
%             for i=1:length(obj.Header)
%                 split = strsplit(obj.Header{i,1},",");
%                 if length(split) > 1
%                     obj.Header{i,2} = str2double(split{2});
%                     obj.Header{i,1} = split{1};
%                 else
%                 end
%             end
%             % Apply scaling multipliers from header to curve data and
%             % delete placeholder data points from curves where piezo range
%             % was too small
%             TempApp = readmatrix(theFiles(1).name);
%             TempHHApp = readmatrix(theFiles(3).name)*(-obj.Header{9,2}) - obj.Header{8,2};
%             obj.NCurves = obj.Header{5,2}*obj.Header{6,2};
%             
%             for i=1:obj.NCurves
%                 DelNApp = round(TempApp(end,i),6,'significant');
%                 if DelNApp == round(TempApp(end-1,i),6,'significant')
%                     DelMap = DelNApp;
%                     k = 1;
%                     while DelNApp == round(TempApp(end-k,i),6,'significant')
%                         k = k + 1;
%                     end
%                     obj.App{i} = TempApp(1:(end-k),i);
%                 else
%                     obj.App{i} = TempApp(:,i);
%                 end
%                 obj.HHApp{i} = TempHHApp(1:length(obj.App{i}),i);
%             end
%             
%             for i=1:obj.NCurves
%                 if round(obj.App{i}(end),6,'significant') == DelMap
%                     obj.App{i}(end) = [];
%                     obj.HHApp{i}(end) = [];
%                 end
%             end
%             
%             obj.HHRet = mat2cell(readmatrix(theFiles(4).name)*obj.Header{11,2},...
%                 [obj.Header{3,2}],ones(obj.Header{5,2}*obj.Header{6,2},1));
%             obj.Ret = mat2cell(readmatrix(theFiles(5).name),...
%                 [obj.Header{3,2}],ones(obj.Header{5,2}*obj.Header{6,2},1));
%             obj.PixApp = obj.Header{2,2};
%             obj.PixRet = obj.Header{3,2};
%             obj.Sensitivity = obj.Header{16,2};
%             obj.SelectedCurves = ones(obj.NCurves,1);
%             obj.Man_CP = zeros(obj.NCurves,2);            
%             obj.NumProfiles = obj.Header{6,2};
%             obj.NumPoints = obj.Header{5,2};
%             obj.SpringConstant = obj.Header{17,2};
%             
%             obj.create_and_level_height_map();
%             
%             cd(current.path);
%             current = what();
%             cd(obj.Folder)
%             savename = sprintf('%s.mat',obj.Name);
%             save(savename,'obj')
%             cd(current.path)
%             disp('loading successfull. object saved in objects folder')
%         end
        
        function obj = ForceMap(MapFullFile,DataFolder)
            %%% Constructor of the class
            
            % Specify the folder where the files live. And import them.
            % Also get curent folder and return to it after import of
            % files.
            % Assigns the properties that can be found in the jpk-file
            % already
            
            current = what();
            
            if nargin < 2
                DataFolder = current.path;
            end
            
            % get OS and use appropriate fitting system command
            
            FullOS = computer;
            OS = FullOS(1:3);
            
            if isequal('PCW',OS)
                % unpack jpk-file into temporary folder to read out data
                cmd1 = '"C:\Program Files\7-Zip\7z.exe" x ';
                cmd2 = '"';
                cmd3 = MapFullFile;
                cmd4 = '"';
                cmd5 = ' -o';
                mkdir(DataFolder,'Temp')
                cmd6 = '"';
                TempFolder = fullfile(DataFolder,'Temp',filesep);
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
                CutExtension = split(Strings{end},'.');
                obj.Name = CutExtension{1};
                mkdir(DataFolder,'ForceData')
                obj.Folder = fullfile(DataFolder,'ForceData',filesep);
                
                
%             system(['unzip -o ', fullfile(datadir,fnamemap), ' ''*shared-data/header.properties'' -d ', tempdir{fib,1}]);
%                 
            elseif isequal('GLN',OS)
                % unpack jpk-file into temporary folder to read out data
                cmd1 = 'unzip -o ';
                cmd2 = MapFullFile;
                cmd3 = ' -d ';
                mkdir(DataFolder,'Temp')
                TempFolder = fullfile(DataFolder,'Temp',filesep);
                CMD = append(cmd1,cmd2,cmd3,TempFolder);
                system(CMD);
                disp('extracting file...')
                Strings = split(MapFullFile,filesep);
                CutExtension = split(Strings{end},'.');
                obj.Name = CutExtension{1};
                mkdir(DataFolder,'ForceData')
                obj.Folder = fullfile(DataFolder,'ForceData',filesep);
                
            elseif isequal('MAC',OS)
                % unpack jpk-file into temporary folder to read out data
                cmd1 = 'unzip -o ';
                cmd2 = MapFullFile;
                cmd3 = ' -d ';
                mkdir(DataFolder,'Temp')
                TempFolder = fullfile(DataFolder,'Temp',filesep);
                CMD = append(cmd1,cmd2,cmd3,TempFolder);
                system(CMD);
                disp('extracting file...')
                Strings = split(MapFullFile,filesep);
                CutExtension = split(Strings{end},'.');
                obj.Name = CutExtension{1};
                mkdir(DataFolder,'ForceData')
                obj.Folder = fullfile(DataFolder,'ForceData',filesep);
            end
            
            
            % reading header properties into object
            filedirectory = fullfile(TempFolder,'header.properties');
            fileID=fopen(filedirectory,'rt','n','UTF-8'); % fileID = fopen(filename,permission,machinefmt,encodingIn)
            A=fileread(filedirectory);
            % Height: 1. CONVERSION raw-meters & 2. SCALE meters
            % Conversion RAW -> VOLTS
            fseek(fileID,1,'cof'); % goes at the first position in the file
            
            %   NCurves
            B=strfind(A,'force-scan-map.indexes.max=');
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
            
            %   NumPoints
            clear tline where;
            frewind(fileID);
            B=strfind(A,'force-scan-map.position-pattern.grid.ilength=');
            fseek(fileID,B,'cof');
            tline = fgetl(fileID);
            where=strfind(tline,'=');
            obj.NumPoints = str2double(tline(where+1:end));
            
            %   NumProfiles
            clear tline where;
            frewind(fileID);
            B=strfind(A,'force-scan-map.position-pattern.grid.jlength=');
            fseek(fileID,B,'cof');
            tline = fgetl(fileID);
            where=strfind(tline,'=');
            obj.NumProfiles = str2double(tline(where+1:end));
            
            
            clear tline A B where
            
            fclose(fileID);
            
            %loading curve data into cell arrays
            for i=1:obj.NCurves
                HeaderFileDirectory = fullfile(TempFolder,'shared-data','header.properties');
                SegmentHeaderFileDirectory = fullfile(TempFolder,'index',string((i-1)),'segments','0','segment-header.properties');
                HeightDataDirectory = fullfile(TempFolder,'index',string((i-1)),'segments','0','channels','capacitiveSensorHeight.dat');
                vDefDataDirectory = fullfile(TempFolder,'index',string((i-1)),'segments','0','channels','vDeflection.dat');
                
                [TempHHApp,obj.App{i},obj.SpringConstant,obj.Sensitivity]=...
                obj.writedata(HeaderFileDirectory,SegmentHeaderFileDirectory,...
                HeightDataDirectory,vDefDataDirectory);
                
                obj.HHApp{i} = -TempHHApp;
                obj.App{i} = obj.App{i}.*obj.SpringConstant;
                clear TempHHApp
            
                SegmentHeaderFileDirectory = fullfile(TempFolder,'index',string((i-1)),'segments','1','segment-header.properties');
                HeightDataDirectory = fullfile(TempFolder,'index',string((i-1)),'segments','1','channels','capacitiveSensorHeight.dat');
                vDefDataDirectory = fullfile(TempFolder,'index',string((i-1)),'segments','1','channels','vDeflection.dat');
            
                [TempHHRet,obj.Ret{i}]=...
                obj.writedata(HeaderFileDirectory,SegmentHeaderFileDirectory,...
                HeightDataDirectory,vDefDataDirectory);
            
                obj.HHRet{i} = -TempHHRet;
                obj.Ret{i} = obj.Ret{i}.*obj.SpringConstant;
                clear TempHHRet
            end
            rmdir(TempFolder,'s');
            
            obj.ExclMask = logical(ones(obj.NumProfiles,obj.NumPoints));
            
            obj.create_and_level_height_map();
            
            obj.SelectedCurves = ones(obj.NCurves,1);
            
            cd(current.path);
            current = what();
            cd(obj.Folder)
            savename = sprintf('%s.mat',obj.Name);
            save(savename,'obj')
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
            end
            % calculate vertical tip position by subtracting vertical tip deflection from head height
            iRange = find(obj.SelectedCurves);
            for i=iRange'
                obj.THApp{i} = obj.HHApp{i} - obj.BasedApp{i}/obj.SpringConstant;
                obj.THRet{i} = obj.HHRet{i} - obj.BasedRet{i}/obj.SpringConstant;
            end
            close(h);
%             current = what();
%             cd(obj.Folder)
%             savename = sprintf('%s.mat',obj.Name);
%             save(savename,'obj')
%             cd(current.path)
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
                obj.CP_RoV(i,:) = [obj.THApp{i}(CPidx) obj.BasedApp{i}(CPidx)];
                CP_RoV = obj.CP_RoV;
            end
            close(h)
%             current = what();
%             cd(obj.Folder)
%             savename = sprintf('%s.mat',obj.Name);
%             save(savename,'obj')
%             cd(current.path)
        end
        
        function estimate_cp_rov_oliver_pharr(obj,batchsize)
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
            obj.CP_OliverPharr_RoV = CP_RoV;
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
                obj.CP_OliverPharr_RoV(i,:) = [obj.HHApp{i}(CPidx) obj.BasedApp{i}(CPidx)];
                obj.CP_OliverPharr = obj.CP_OliverPharr_RoV;
            end
            close(h)
%             current = what();
%             cd(obj.Folder)
%             savename = sprintf('%s.mat',obj.Name);
%             save(savename,'obj')
%             cd(current.path)
        end
        
        function estimate_cp_gof(obj)
            obj.CP_GoF = zeros(obj.NCurves,2);
            Range = find(obj.SelectedCurves);
            h = waitbar(0,'Setting up...','Name',obj.Name);
            for i=Range'
                smoothx = smoothdata(obj.THApp{i});
                smoothy = smoothdata(obj.BasedApp{i});
                Rsquare = 2*ones(length(obj.BasedApp{i}),1);
                E = ones(length(obj.BasedApp{i}),1);
                obj.DeltaE = ones(length(obj.BasedApp{i}),1);
                testrange = floor(0.5*length(obj.BasedApp{i})):(length(obj.BasedApp{i})-5);
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
                obj.CP_GoF(i,:) = [obj.THApp{i}(CPidx) obj.BasedApp{i}(CPidx)];
                %deltaE(floor(0.5*length(obj.BasedApp{i})):(length(obj.BasedApp{i})-6)) = -diff(log(E));
            end
            close(h)
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
                obj.CP_Combo(i,:) = [obj.THApp{i}(CPidx) obj.BasedApp{i}(CPidx)];
            end
%             current = what();
%             cd(obj.Folder)
%             savename = sprintf('%s.mat',obj.Name);
%             save(savename,'obj')
%             cd(current.path)
        end
                
        function estimate_cp_cnn(obj,NeuralNet,RunMode,NumPasses)
            % RunMode = 'Fast' just predict with one forwardpass using
            % CP_CNN
            % RunMode = 'Dropout' predict through n=NumPasses forwardpasses
            % through DropoutNet and additionally gain a metric for
            % uncertainty of the CP estimate
            if nargin < 2
                runmode = 0;
            elseif isequal(RunMode,'Fast')
                runmode = 0;
            elseif isequal(RunMode,'Dropout')
                runmode = 1;
                if nargin < 3
                    NumPasses = 100; % if not specified in arguments, NumPasses defaults to 100
                end
            end
            ImgSize = NeuralNet.Layers(1).InputSize;
            objcell{1,1} = obj;
            X = obj.CP_CNN_batchprep(objcell,ImgSize(1));
            len = size(X,4);
            h = waitbar(0,'Setting up','Name',obj.Name);
            switch runmode
                case 0
                    waitbar(1/2,h,'Predicting CP');
                    Ypredicted = predict(NeuralNet,X);
                    iRange = find(obj.SelectedCurves);
                    k = 1;
                    for i=iRange'
                        obj.CP(i,1) = Ypredicted(k,1)*range(obj.THApp{i})+min(obj.THApp{i});
                        obj.CP(i,2) = Ypredicted(k,2)*range(obj.BasedApp{i})+min(obj.BasedApp{i});
                        k = k + 1;
                    end
                case 1
                    obj.YDropPred = zeros(NumPasses,2,len);
                    CantHandle = false;
                    for j=1:NumPasses
                        waitbar(j/NumPasses,h,sprintf('Predicting CP for %i curves. %i/%i passes done',len,j,NumPasses));
                        if CantHandle == false
                            try
                                Temp = predict(NeuralNet,X,'MiniBatchSize',len);
                            catch
                                CantHandle = true;
                                warning("Can't compute with optimal MiniBatchSize. Choosing a lower one")
                                try
                                    Temp = predict(NeuralNet,X);
                                catch
                                    warning("Can't compute with lowered MiniBatchSize. Choosing an even lower one. Computation time is gonna suffer")
                                    Temp = predict(NeuralNet,X,'MiniBatchSize',10);
                                end
                            end
                        else
                            try
                                Temp = predict(NeuralNet,X);
                            catch
                                warning("Can't compute with lowered MiniBatchSize. Choosing an even lower one \nComputation time is gonna suffer")
                                Temp = predict(NeuralNet,X,'MiniBatchSize',10);
                            end
                        end
                        obj.YDropPred(j,:,:) = Temp';
                    end
                    waitbar(1,h,'Wrapping up');
                    iRange = find(obj.SelectedCurves);
                    k = 1;
                    obj.CP_MonteCarlo = zeros(NumPasses,2,obj.NCurves);
                    for i=iRange'
                        obj.CP_MonteCarlo(:,1,i) = obj.YDropPred(:,1,k)*range(obj.THApp{i})+min(obj.THApp{i});
                        obj.CP_MonteCarlo(:,2,i) = obj.YDropPred(:,2,k)*range(obj.BasedApp{i})+min(obj.BasedApp{i});
                        obj.CP(i,1) = mean(obj.CP_MonteCarlo(:,1,i));
                        obj.CP(i,2) = mean(obj.CP_MonteCarlo(:,2,i));
                        k = k + 1;
                    end
            end
            close(h)
%             current = what();
%             cd(obj.Folder)
%             savename = sprintf('%s.mat',obj.Name);
%             save(savename,'obj')
%             cd(current.path)
        end
        
        function estimate_cp_oliver_pharr(obj,NeuralNet,RunMode,NumPasses)
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
            % NumPasses <= 30 is recommended
            
            if nargin < 2
                runmode = 0;
            elseif isequal(RunMode,'Fast')
                runmode = 0;
            elseif isequal(RunMode,'Dropout')
                runmode = 1;
                if nargin < 3
                    NumPasses = 100; % if not specified in arguments, NumPasses defaults to 100
                end
            end
            ImgSize = NeuralNet.Layers(1).InputSize;
            objcell{1,1} = obj;
            X = obj.CP_oliver_pharr_batchprep(objcell,ImgSize(1));
            len = size(X,4);
            h = waitbar(0,'Setting up','Name',obj.Name);
            switch runmode
                case 0
                    waitbar(1/2,h,'Predicting CP');
                    Ypredicted = predict(NeuralNet,X);waitbar(1,h,'Wrapping up');
                    iRange = find(obj.SelectedCurves);
                    k = 1;
                    obj.CP_OliverPharr_MonteCarlo = zeros(NumPasses,2,obj.NCurves);
                    for i=iRange'
                        obj.CP_OliverPharr(i,1) = Ypredicted(k,1)*range(obj.HHApp{i})+min(obj.HHApp{i});
                        obj.CP_OliverPharr(i,2) = Ypredicted(k,2)*range(obj.BasedApp{i})+min(obj.BasedApp{i});
                        k = k + 1;
                    end
                case 1
                    obj.YDropPred = zeros(NumPasses,2,len);
                    CantHandle = false;
                    for j=1:NumPasses
                        waitbar(j/NumPasses,h,sprintf('Predicting CP for %i curves. %i/%i passes done',len,j,NumPasses));
                        if CantHandle == false
                            try
                                Temp = predict(NeuralNet,X,'MiniBatchSize',len);
                            catch
                                CantHandle = true;
                                warning("Can't compute with optimal MiniBatchSize. Choosing a lower one")
                                try
                                    Temp = predict(NeuralNet,X);
                                catch
                                    warning("Can't compute with lowered MiniBatchSize. Choosing an even lower one. Computation time is gonna suffer")
                                    Temp = predict(NeuralNet,X,'MiniBatchSize',10);
                                end
                            end
                        else
                            try
                                Temp = predict(NeuralNet,X);
                            catch
                                warning("Can't compute with lowered MiniBatchSize. Choosing an even lower one \nComputation time is gonna suffer")
                                Temp = predict(NeuralNet,X,'MiniBatchSize',10);
                            end
                        end
                        obj.YDropPred(j,:,:) = Temp';
                    end
                    waitbar(1,h,'Wrapping up');
                    iRange = find(obj.SelectedCurves);
                    k = 1;
                    obj.CP_OliverPharr_MonteCarlo = zeros(NumPasses,2,obj.NCurves);
                    for i=iRange'
                        obj.CP_OliverPharr_MonteCarlo(:,1,i) = obj.YDropPred(:,1,k)*range(obj.HHApp{i})+min(obj.HHApp{i});
                        obj.CP_OliverPharr_MonteCarlo(:,2,i) = obj.YDropPred(:,2,k)*range(obj.BasedApp{i})+min(obj.BasedApp{i});
                        obj.CP_OliverPharr(i,1) = mean(obj.CP_OliverPharr_MonteCarlo(:,1,i));
                        obj.CP_OliverPharr(i,2) = mean(obj.CP_OliverPharr_MonteCarlo(:,2,i));
                        obj.CP_OliverPharr_CNN(i,1) = obj.CP_OliverPharr(i,1);
                        obj.CP_OliverPharr_CNN(i,2) = obj.CP_OliverPharr(i,2);
                        obj.CP_OliverPharr_MonteCarlo_STD(i) = norm([std(obj.CP_OliverPharr_MonteCarlo(:,1,i)) std(obj.CP_OliverPharr_MonteCarlo(:,2,i))]);
                        k = k + 1;
                    end
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
            % analysis script 'A_nIAFM_analysis_main.m'
            iRange = find(obj.SelectedCurves);
            for i=iRange'
                load = zeros(length(obj.BasedApp{i}),2);
                unload = zeros(length(obj.Ret{i}),2);
                for j=1:length(obj.BasedApp{i})
                    load(end-(j-1),1) = obj.THApp{i}(j);
                    load(j,2) = obj.BasedApp{i}(j);
                end
                for j=1:length(obj.Ret{i})
                    unload(end-(j-1),1) = obj.HHRet{i}(j);
                    unload(j,2) = obj.Ret{i}(j);
                end
                [obj.LoadOld{i},obj.UnloadOld{i},Position,obj.CP_old(i,2)] = ContactPoint_sort_mr(load,unload);
                obj.CP_old(i,1) = obj.THApp{i}(Position);
            end
%             current = what();
%             cd(obj.Folder)
%             savename = sprintf('%s.mat',obj.Name);
%             save(savename,'obj')
%             cd(current.path)
        end
        
        function estimate_cp_old_oliver_pharr(obj)
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
                [obj.LoadOld{i},obj.UnloadOld{i},Position,obj.CP_OliverPharr(i,2)] = ContactPoint_sort_mr(load,unload);
                obj.CP_OliverPharr(i,1) = obj.HHApp{i}(Position);
                obj.CP_OliverPharr_Old(i,1) =obj.CP_OliverPharr(i,1);
                obj.CP_OliverPharr_Old(i,2) =obj.CP_OliverPharr(i,2);
            end
%             current = what();
%             cd(obj.Folder)
%             savename = sprintf('%s.mat',obj.Name);
%             save(savename,'obj')
%             cd(current.path)
        end
        
        function estimate_cp_manually(obj)
            % Manual contact point selection for NN training on plotted force curves
            % returning a position vector in meters.
            jRange = find(obj.SelectedCurves);
            fig = figure('Name',obj.Name);
            for j=jRange'
                fig.WindowState = 'fullscreen';
                plot(obj.THApp{j},obj.BasedApp{j});
                plottitle = sprintf('Curve Nr.%i/%i\n Click and drag the point to the contact point\n Confirm with any key press',j,obj.NCurves);
                
                plottitle = sprintf('Curve Nr.%i/%i\n Click and drag the point to the contact point\n Confirm with any key press',j,obj.NCurves);
                title(plottitle);
                [~, domainidx] = ForceMap.no_contact_domain(obj.App{j});
                axis([obj.THApp{j}(floor(domainidx*0.2)) inf -inf inf])
                CP_point = drawpoint();
                
                obj.Man_CP(j,1) = CP_point.Position(1);
                obj.Man_CP(j,2) = CP_point.Position(2);
            end
            close(fig);
%             current = what();
%             cd(obj.Folder)
%             savename = sprintf('%s.mat',obj.Name);
%             save(savename,'obj')
%             cd(current.path)
        end
        
        function [E,HertzFit] = calculate_e_mod_hertz(obj,CPType,TipShape,curve_percent)
            % calculate the E modulus of the chosen curves using the CP
            % type chosen in the arguments fitting the lower curve_percent
            % part of the curves
            if ~exist('curve_percent','var') && ischar('curve_percent')
                curve_percent = 0.75;
            end
            if ~exist('TipShape','var') && ~ischar(TipShape)
                TipShape = 'parabolic';
            end
            iRange = find(obj.SelectedCurves);
            obj.EModHertz = zeros(obj.NCurves,1);
            for i=iRange'
                if isequal(lower(CPType),'cnn')
                    CP = obj.CP(i,:);
                elseif isequal(lower(CPType),'old')
                    CP = obj.CP_old(i,:);
                elseif isequal(lower(CPType),'rov')
                    CP = obj.CP_RoV(i,:);
                elseif isequal(lower(CPType),'gof')
                    CP = obj.CP_GoF(i,:);
                elseif isequal(lower(CPType),'combo')
                    CP = obj.CP_Combo(i,:);
                elseif isequal(lower(CPType),'manual')
                    CP = obj.Man_CP(i,:);
                else
                    CP = obj.CP(i,:);
                end
                tip_h = obj.THApp{i} - CP(1);
                tip_h(tip_h < 0) = [];
                obj.IndDepthHertz(i) = max(tip_h);
                force = obj.BasedApp{i} - CP(2);
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
                    s = fitoptions('Method','NonlinearLeastSquares',...
                        'Lower',10^(-30),...
                        'Upper',inf,...
                        'Startpoint',1);
                    f = fittype('a*(x)^(3/2)','options',s);
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
                elseif isequal(shape,'spherical')
                elseif isequal(shape,'conical')
                elseif isequal(shape,'pyramid')
                end
                obj.EModHertz(i) = EMod;
                % Convert the model to the right scale so it can be plotted
                % correctly later
                warning('off','all');
                Hertzfit.a = Hertzfit.a*RangeF/RangeTH^(3/2);
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
                    obj.EModMapHertz(i,j,1) = (obj.EModHertz(mod(i,2)*j+(1-mod(i,2))*(obj.NumPoints-(j-1))+(obj.NumProfiles-i)*obj.NumPoints));
                end
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
            obj.IndDepth = zeros(obj.NCurves,1);
            obj.IndentArea = zeros(obj.NCurves,1);
            
            obj.calculate_reference_slope();

            for i=Range'
                Z(:,i) = obj.HHRet{i} - obj.CP_OliverPharr(i,1);
                D(:,i) = (obj.BasedRet{i} - obj.CP_OliverPharr(i,2))/obj.SpringConstant;
                Zmax(i) = max(Z(:,i));
                Dmax(i) = max(D(:,i));
                DCurvePercent = D(D(:,i)>=(1-CurvePercent)*Dmax(i),i);
                ZCurvePercent = Z(1:length(DCurvePercent),i);
                LineFit = polyfit(ZCurvePercent,DCurvePercent,1);
                obj.DZslope(i) = LineFit(1);
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
                    obj.EModMapOliverPharr(i,j,1) = (obj.EModOliverPharr(mod(i,2)*j+(1-mod(i,2))*(obj.NumPoints-(j-1))+(obj.NumProfiles-i)*obj.NumPoints));
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
            obj.RectApex = zeros(obj.NumProfiles,1);
            obj.RectApexIndex = zeros(obj.NumProfiles,1);
            obj.RectApexIndex = round(predictGP_mean([1:obj.NumProfiles],[1:obj.NumProfiles],1,5*obj.NumProfiles,obj.ApexIndex,1));
            for i=1:obj.NumProfiles
                obj.RectApex(i) = obj.HeightMap(i,obj.RectApexIndex(i),1);
            end
            k = 1;
            for i=1:obj.NumProfiles
                if obj.ExclMask(i,obj.RectApexIndex(i)) == 1
                    FibHeight(k) = obj.RectApex(i);
                    k = k + 1;
                end
            end
            obj.FibDiam = mean(FibHeight);
            
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
                HeighDataDirectory,vDelfDataDirectory)
            % Author: Orestis Andriotis (slightly changed by Manuel Rufin)
            % HeaderFileDirectory: file directory of the main header properties. There
            % is the information about the scaling factors to convert the raw data.
            % SegmentHeaderFileDirectory: file directory of the header of each segment.
            % Each segment means the loading (#0) and unloading (#1).
            
            % HeighDataDirectory: file directory of height.dat file
            % vDelfDataDirectory: file directory of vDeflection.dat file
            %
            
            % HeaderFileDirectory = headerDir;
            % SegmentHeaderFileDirectory=segHeaderDir;
            % HeighDataDirectory=heightdir;
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
            fileID = fopen(HeighDataDirectory);
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
                ForceMap.getheaderinfo(HeaderFileDirectory);
            % converts raw data and writes them into one varialble 'fdata'
            
            % converts raw into volts by calling the getheaderinfo function
            heightvolts = RawHeight.*mult_height_meters1 + offset_height_meters1;
            
            deflvolts=RawvDeflection.*mult_vDefl_volts + offset_vDefl_volts;
            
            % converts volts into meters
            height_meters = heightvolts.*mult_height_meters2 + offset_height_meters2;
            % height_meters = heightvolts.*mult_height_meters + offset_height_meters;
            deflmeter=deflvolts.*sensitivity;
            
            HeadHeight = height_meters;
            
            Force = deflmeter;
            
            clear tline SegmentHeaderFileDirectory HeighDataDirectory...
                vDelfDataDirectory mult_height_volts offset_height_volts...
                mult_height_meters offset_height_meters...
                mult_vDefl_volts offset_vDefl_volts...
                mult_vDefl_meters ;
            fclose(fileID);
        end
        
        function [mult_height_meters1, offset_height_meters1,...
                mult_height_meters2, offset_height_meters2,...
                mult_vDefl_volts, offset_vDefl_volts,...
                sensitivity, spring_constant] = getheaderinfo(filedirectory)
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
            
            %   Multiplier
            B=strfind(A,'lcd-info.3.encoder.scaling.multiplier=');
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
            B=strfind(A,'lcd-info.3.encoder.scaling.offset=');
            fseek(fileID,B,'cof');
            tline = fgetl(fileID);
            where=strfind(tline,'=');
            offset_height_meters1 = str2double(tline(where+1:end));
            
            
            % SCALING METERS -> METERS
            
            %   Multiplier
            clear tline where;
            frewind(fileID);
            B=strfind(A,'lcd-info.3.conversion-set.conversion.nominal.scaling.multiplier=');
            fseek(fileID,B,'cof');
            tline = fgetl(fileID);
            where=strfind(tline,'=');
            mult_height_meters2 = str2double(tline(where+1:end));
            
            
            %   Offset
            clear tline where;
            frewind(fileID);
            B=strfind(A,'lcd-info.3.conversion-set.conversion.nominal.scaling.offset=');
            fseek(fileID,B,'cof');
            tline = fgetl(fileID);
            where=strfind(tline,'=');
            offset_height_meters2 = str2double(tline(where+1:end));
            
            
            % vDeflection: 1. CONVERSION raw-volts & 2. volts to meters
            % Conversion RAW -> VOLTS
            
            %   Multiplier
            clear tline where;
            frewind(fileID);
            B=strfind(A,'lcd-info.2.encoder.scaling.multiplier=');
            fseek(fileID,B,'cof');
            tline = fgetl(fileID);
            where=strfind(tline,'=');
            mult_vDefl_volts = str2double(tline(where+1:end)); % multiplier for scaling the raw height data and convert to volts
            
            %   Offset
            clear tline where;
            frewind(fileID);
            B=strfind(A,'lcd-info.2.encoder.scaling.offset=');
            fseek(fileID,B,'cof');
            tline = fgetl(fileID);
            where=strfind(tline,'=');
            offset_vDefl_volts = str2double(tline(where+1:end));
            
            
            % Conversion VOLTS -> METERS
            
            %   Multiplier (that is the sensitivity measured in meters per Volts)
            clear tline where;
            frewind(fileID);
            B=strfind(A,'lcd-info.2.conversion-set.conversion.distance.scaling.multiplier=');
            fseek(fileID,B,'cof');
            tline = fgetl(fileID);
            where=strfind(tline,'=');
            sensitivity = str2double(tline(where+1:end));
            
            % Spring constant
            
            clear tline where;
            frewind(fileID);
            B=strfind(A,'lcd-info.2.conversion-set.conversion.force.scaling.multiplier=');
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
                    Image = imresize(Image,[ImgSizeFinal ImgSizeFinal]);
                    X(:,:,1,k) = Image;
                    k = k + 1;
                    clear Points
                end
            end
        end
        
        function X = CP_batchprep_new_fast(objcell,ImgSize)
            
            if nargin<2
                ImgSize = 128;
                ImgSizeFinal = 128;
            else
                ImgSizeFinal = ImgSize;
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
                    X(:,:,1,k) = Image;
                    k = k + 1;
                    clear Points
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
            save(savename,'obj')
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
        
        function calculate_reference_slope(obj)
            % Calculates the distribution of DZslopes on the curves 
            % that that are neihter on the fibril nor the excluded zones.
            %  the upper 25% of the curve are considered for the
            %  calculation
            
            Range = find(obj.SelectedCurves);
            CurvePercent = 0.25;
            % Calculate the DZslopes
            k = 1;
            for i=Range'
                if (obj.FibMask(obj.List2Map(i,1),obj.List2Map(i,2)) == 0) &&...
                        (obj.ExclMask(obj.List2Map(i,1),obj.List2Map(i,2)) == 1)
                    Z(:,i) = obj.HHRet{i} - obj.CP_OliverPharr(i,1);
                    D(:,i) = (obj.BasedRet{i} - obj.CP_OliverPharr(i,2))/obj.SpringConstant;
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
        end
        
        function create_and_level_height_map(obj)
            % first set Height Map to default values for reproducable
            % results
            k = 1;
            obj.List2Map = zeros(obj.NCurves,2);
            for i=1:obj.NumProfiles
                for j=1:obj.NumPoints
                    obj.HeightMap(i,j) = -max(obj.HHApp{mod(i,2)*j+(1-mod(i,2))*(obj.NumPoints-(j-1))+(obj.NumProfiles-i)*obj.NumPoints});
                    obj.Map2List(i,j) = k;
                    obj.List2Map(k,:) = [i j];
                    k = k + 1;
                end
            end
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
            for i=1:10
                mask = medfilt2(mask,[5 3],'symmetric');
            end
            mask = bwareafilt(mask,1,4);
            obj.FibMask = mask;
            
%             current = what();
%             cd(obj.Folder)
%             savename = sprintf('%s.mat',obj.Name);
%             save(savename,'obj')
%             cd(current.path)
        end
        
    end
    
    methods
       % methods for visualization, plotting, statistics and quality control 
        
        function show_random_curve_THvsForce(obj)
            k = randi(obj.NCurves);
            figure('Units','normalized','Position',[0.6 0.1 0.4 0.8])
            subplot(2,1,1)
            plot(obj.THApp{k},obj.BasedApp{k},obj.THRet{k},obj.BasedRet{k})
            title(sprintf('Curve Nr.%i of %s',k,obj.Name))
            legend('Based Approach','Based Retract','Location','northwest')
            xlabel('Tip Height [m]')
            ylabel('Force [N]')
            drawpoint('Position',[obj.CP(k,1) obj.CP(k,2)]);
            subplot(2,1,2)
            imshow(mat2gray(obj.HeightMap))
            axis on
            hold on;
            plot(obj.List2Map(k,2),obj.List2Map(k,1), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
            
        end
        
        function show_random_curve_HHvsDefl(obj,k)
            if nargin < 2
                k = randi(obj.NCurves);
            end
            figure('Units','normalized','Position',[0.6 0.1 0.4 0.8],'Color','white')
            subplot(2,1,1)
            plot(obj.HHApp{k},obj.BasedApp{k}/obj.SpringConstant,...
                obj.HHRet{k},obj.BasedRet{k}/obj.SpringConstant)
            title(sprintf('Curve Nr.%i of %s',k,obj.Name))
            legend('Based Approach','Based Retract','Location','northwest')
            xlabel('Head Height [m]')
            ylabel('vDeflection [m]')
            drawpoint('Position',[obj.CP_OliverPharr(k,1) obj.CP_OliverPharr(k,2)/obj.SpringConstant]);
            subplot(2,1,2)
            I = obj.HeightMap;
            I = mat2gray(I);
            imshow(I)
            axis on
            hold on;
            plot(obj.List2Map(k,2),obj.List2Map(k,1), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
            
        end
        
        function show_height_map(obj)
            T = sprintf('Height Map of %s\nwith chosen indentation points',obj.Name);
            figure('Name',T,'Units','normalized','Position',[0.5 0.1 0.5 0.8]);
            
            subplot(2,2,1)
            I = imresize(mat2gray(obj.HeightMap(:,:,1)),[1024 1024]);
            I = (I*range(obj.HeightMap(:,:,1),'all') + min(obj.HeightMap(:,:,1),[],'all'))*1e9;
            imshow(I,'Colormap',hot)
            title(sprintf('%s Plane Fitted Height',obj.Name))
            c1 = colorbar;
            c1.Label.String = 'Height [nm]';
            hold on;
            for i=1:obj.NumProfiles
                plot((obj.List2Map(obj.RectApexIndex(i),2)-1/2)*1024/obj.NumPoints,...
                    (obj.List2Map(obj.RectApexIndex(i),1)-1/2)*1024/obj.NumProfiles,...
                    'r+', 'MarkerSize', 10, 'LineWidth', 1);
                plot((obj.List2Map(obj.ApexIndex(i),2)-1/2)*1024/obj.NumPoints,...
                    (obj.List2Map(obj.ApexIndex(i),1)-1/2)*1024/obj.NumProfiles,...
                    'g+', 'MarkerSize', 10, 'LineWidth', 1);
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
        
        function quality_control_oliver_pharr(obj,PauseTime)
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
                drawpoint('Position',[obj.CP_OliverPharr(k,1) obj.CP_OliverPharr(k,2)]);
                
                if m == 1
                    subplot(2,3,3)
                    boxplot(obj.EModOliverPharr(obj.RectApexIndex));
                    xticklabels(obj.Name)
                    title(sprintf('mean = %.2f MPa\nmedian = %.2f MPa\nstd = %.3f MPa',...
                        mean(obj.EModOliverPharr(obj.RectApexIndex))*1e-6,...
                        median(obj.EModOliverPharr(obj.RectApexIndex))*1e-6,...
                        std(obj.EModOliverPharr(obj.RectApexIndex))*1e-6));
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
        
        function quality_control_hertz_sneddon(obj,PauseTime)
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
                plot(obj.THApp{k} - obj.CP(k,1),feval(obj.HertzFit{k},obj.THApp{k} - obj.CP(k,1)),...
                    obj.THApp{k} - obj.CP(k,1),obj.BasedApp{k} - obj.CP(k,2),...
                    obj.THRet{k} - obj.CP(k,1),obj.BasedRet{k} - obj.CP(k,2))
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
                        mean(obj.EModHertz(obj.RectApexIndex))*1e-6,...
                        median(obj.EModHertz(obj.RectApexIndex))*1e-6,...
                        std(obj.EModHertz(obj.RectApexIndex))*1e-6));
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
        
        function compare_hertz_oliver(obj)
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