classdef ForceMap < handle
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
    % FM.tip_radius;
    % to get a class parameter of this force map (the tip radius of the used cantilever)
    
    properties
        Name            % name of the force map. taken as the name of the folder, containing the .csv files
        Folder          % location of the .csv files of the force map
        Header          % header properties taken from the JPK force map container file
        App = {}        % approach force data in Newton
        Ret = {}        % retraction force data in Newton
        HHApp = {}      % head-height approach data in meters
        HHRet = {}      % head-height retract data in meters
        THApp = {}      % vertical tip height approach data in meters
        THRet = {}      % vertical tip height retract data in meters
        BasedApp = {}   % approach force data with subtracted base line and tilt in Newton
        BasedRet = {}   % retraction force data with subtracted base line and tilt in Newton
        Basefit = {}    % fit model used for the baseline fit
        HeightMap       % height profile map taken from the inverted maximum head-height from approach max(hhapp)
                        % the additional dimension in the height map holds
                        % the index from the corresponding 1-D lists of
                        % force curves, so one can easily identify the
                        % right curves
        Sensitivity
        NCurves         % number of curves on the force map
        PixApp 
        PixRet
        SelectedCurves  % logical vector of length n_curves with 0s for excluded and 1s for included curves. gets initialized with ones
        RoV = {}        % ratio of variance curve for CP estiamation
        RoVTH = {}      % ratio of variance curve, based on the vertical tip position. [does not work well/experimental development stage]
        GoF = {}        % goodness of fit curve for each selected curve
        CPCombo = {}    % combination of the various metrics for contact point estimation
        CP              % chosen contact point (CP) for every selected curve, which is then used in E module fitting
        CP_CNN          % Convolutional Neural Network for CP-estimation used in CP_CNN_predict
        DropoutNet      % Convolutional Neural Network for uncertainty estimation
        YDropPred       % Contains the Dropoutpredictions for every curve in the forcemap
        Man_CP          % manually chosen contact point
        CP_old          % contact point estimation from old script 'A_nIAFM_analysis_main'
        LoadOld         % comes from same script as CP_old
        UnloadOld       % comes from same script as CP_old
        DeltaE = {}     %
        TipRadius = -1  % (nominal) tip radius in nm for the chosen type of tip model
        PoissonR = 0.5  % standard Poisson ratio for most mechanical models
        
    end
    methods
        
        function obj = ForceMap(mapfilepath)
            %%% Constructor of the class
            
            % Specify the folder where the files live. And import them.
            % Also get curent folder and return to it after import of
            % files.
            % Assigns the properties: name, folder, header, app, ret,
            % hhapp, hhret, height_map, sensitivity, n_curves, pix_app,
            % pix_ret and initiates selected_curves with ones (all curves chosen)
            
            current = what();
            
            % determine if ForceMap is given a loadpath for existing .mat
            if nargin==1
                cd(mapfilepath);
                file = dir('*.mat');
                load(file.name);
                cd(current.path);
                msg = sprintf('loading %s',obj.Name);
                disp(msg);
                disp('loading successfull')
                return
            end
            
            quest = 'Do you want to load an already existing .mat file of the force map?';
            answer = questdlg(quest,'Load map...','...from .mat-file','...from .cvs-folder','...from .cvs-folder');
            if isequal(answer, '...from .mat-file')
                [mapfile, mapfilepath] = uigetfile();
                cd(mapfilepath);
                load(mapfile,'-mat');
                cd(current.path);
                msg = sprintf('loading %s',obj.Name);
                disp(msg);
                disp('loading successfull')
                obj.Folder = mapfilepath;
                return
            else
            end
            
            
            obj.Folder = uigetdir;
            cd(obj.Folder);
            splitfolder = strsplit(obj.Folder,'\');
            obj.Name = string(splitfolder(end));
            msg = sprintf('loading %s',obj.Name);
            disp(msg);
            
            % Get a list of all files in the folder with the desired file name pattern.
            % Change to whatever pattern you used for the names of your force map files.
            filePattern = fullfile(obj.Folder, '*.csv');
            theFiles = dir(filePattern);
            % Check to make sure files of specified patterns have been found
            if length(theFiles) < 1
                errorMessage = 'No files of the specified name pattern have been found in the current directory';
                uiwait(warndlg(errorMessage));
                return;
            end
            
            % Import data from csv-files
            import_header = importdata(theFiles(2).name);
            k = 1;
            for i=1:length(import_header)
                if isempty(strfind(import_header{i},'start-option.'))
                    obj.Header{k,1} = import_header{i};
                    k = k + 1;
                end
            end
            for i=1:length(obj.Header)
                split = strsplit(obj.Header{i,1},",");
                if length(split) > 1
                    obj.Header{i,2} = str2double(split{2});
                    obj.Header{i,1} = split{1};
                else
                end
            end
            % Apply scaling multipliers from header to curve data and
            % delete placeholder data points from curves where piezo range
            % was too small
            TempApp = readmatrix(theFiles(1).name);
            TempHHApp = readmatrix(theFiles(3).name)*obj.Header{11,2};
            obj.NCurves = obj.Header{5,2}*obj.Header{6,2};
            
            for i=1:obj.NCurves
                DelNApp = round(TempApp(end,i),6,'significant');
                if DelNApp == round(TempApp(end-1,i),6,'significant')
                    DelMap = DelNApp;
                    k = 1;
                    while DelNApp == round(TempApp(end-k,i),6,'significant')
                        k = k + 1;
                    end
                    obj.App{i} = TempApp(1:(end-k),i);
                else
                    obj.App{i} = TempApp(:,i);
                end
                obj.HHApp{i} = TempHHApp(1:length(obj.App{i}),i);
            end
            
            for i=1:obj.NCurves
                if round(obj.App{i}(end),6,'significant') == DelMap
                    obj.App{i}(end) = [];
                    obj.HHApp{i}(end) = [];
                end
            end
            
            obj.HHRet = mat2cell(readmatrix(theFiles(4).name)*obj.Header{11,2},...
                [obj.Header{3,2}],ones(obj.Header{5,2}*obj.Header{6,2},1));
            obj.Ret = mat2cell(readmatrix(theFiles(5).name),...
                [obj.Header{3,2}],ones(obj.Header{5,2}*obj.Header{6,2},1));
            obj.PixApp = obj.Header{2,2};
            obj.PixRet = obj.Header{3,2};
            obj.Sensitivity = obj.Header{16,2};
            obj.SelectedCurves = ones(obj.NCurves,1);
            for i=1:obj.Header{6,2}
                for j=1:obj.Header{5,2}
                    obj.HeightMap(i,j,1) = -max(obj.HHApp{mod(i,2)*j+(1-mod(i,2))*(obj.Header{5,2}-(j-1))+(obj.Header{6,2}-i)*obj.Header{5,2}});
                    obj.HeightMap(i,j,2) = (mod(i,2)*j+(1-mod(i,2))*(obj.Header{5,2}-(j-1))+(obj.Header{6,2}-i)*obj.Header{5,2});
                end
            end
            obj.Man_CP = zeros(obj.NCurves,2);
            try
                temp = load('CP_CNN');
                obj.CP_CNN = temp.CP_CNN;
            catch ME1
                disp("Can't find Neural Network named 'CP_CNN' in Folder. Trying to load from Workspace instead...")
                try
                    obj.CP_CNN = CP_CNN;
                catch ME2
                    disp("Can't find Neural Network named 'CP_CNN' in Workspace. Loading ForceMap without it")
                end
            end
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
            current = what();
            cd(obj.Folder)
            savename = sprintf('%s.mat',obj.Name);
            save(savename,'obj')
            cd(current.path)
        end
        
        function base_and_tilt(obj)
            % subtract baseline and tilt from the forcecurve by fitting a function to
            % the non contact domain the function tries to fit a non-affine-linear
            % function. If the linear fit is too bad, the function tries a 9th grade
            % polynomial fit instead
            Range = find(obj.SelectedCurves);
            h = waitbar(0,'Setting up...');
            for i=Range'
                prog = i/obj.NCurves;
                waitbar(prog,h,'processing baseline fits...');
                [ncd , ncidx] = ForceMap.no_contact_domain(obj.App{i});
                gof.rsquare = 0;
                for j=1:10
                    [testfit,goftest] = fit(obj.HHApp{i}(1:ncidx),smooth(ncd),'poly1','Normalize','on');
                    if goftest.rsquare>gof.rsquare
                        gof = goftest;
                        obj.Basefit{i} = testfit;
                    end
                end
                
                if gof.rsquare < 0.9
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
                obj.THApp{i} = obj.HHApp{i} - obj.BasedApp{i}/obj.Header{17,2};
                obj.THRet{i} = obj.HHRet{i} - obj.BasedRet{i}/obj.Header{17,2};
            end
            close(h);
            current = what();
            cd(obj.Folder)
            savename = sprintf('%s.mat',obj.Name);
            save(savename,'obj')
            cd(current.path)
        end
        
        function cp_rov(obj,batchsize)
            % find contact point with the method of ratio of variances. The method
            % iterates though every point and builds the ratio of the variance of a
            % bunch of points before and after the current point. the point with the
            % biggest ratio is the returned contact point [Nuria Gavara, 2016]
            if nargin<2
                batchsize = 10;
            end
            Range = find(obj.SelectedCurves);
            h = waitbar(0,'Setting up...','Name',obj.Name);
            for i=Range'
                prog = i/obj.NCurves;
                waitbar(prog,h,'applying ratio of variances method...');
                obj.RoV{i} = zeros(length(obj.BasedApp{i}),1);
                
                % loop through points and calculate the RoV
                for j=(batchsize+1):(length(obj.BasedApp{i})-batchsize)
                    obj.RoV{i}(j,1) = var(smoothdata(obj.BasedApp{i}((j+1):(j+batchsize))))/...
                        var(smoothdata(obj.BasedApp{i}((j-batchsize):(j-1))));
                end
                % normalize RoV-curve
                obj.RoV{i} = obj.RoV{i}/range(obj.RoV{i});
                minrov = min(obj.RoV{i}(batchsize+1:length(obj.RoV{i})-batchsize));
                obj.RoV{i}(obj.RoV{i}==0) = minrov;
            end
            close(h)
            current = what();
            cd(obj.Folder)
            savename = sprintf('%s.mat',obj.Name);
            save(savename,'obj')
            cd(current.path)
        end
        
        function [E_mod,GoF,Hertzfit] = hertz_fit(obj,tip_h,force,CP,curve_percent,shape)
            
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
            current = what();
            cd(obj.Folder)
            savename = sprintf('%s.mat',obj.Name);
            save(savename,'obj')
            cd(current.path)
        end
        
        function cp_gof(obj)
            Range = find(obj.SelectedCurves);
            for i=Range'
                smoothx =obj.THApp{i};
                smoothy =obj.BasedApp{i};
                Rsquare = 2*ones(length(obj.BasedApp{i}),1);
                E = ones(length(obj.BasedApp{i}),1);
                deltaE = ones(length(obj.BasedApp{i}),1);
                testrange = floor(0.5*length(obj.BasedApp{i})):(length(obj.BasedApp{i})-5);
                h = waitbar(0,'Setting up...');
                msg = sprintf('applying goodness of fit method on curve Nr.%i/%i',i,obj.NCurves);
                for j=testrange
                    prog = (j-testrange(1))/length(testrange);
                    waitbar(prog,h,msg);
                    [emod,gof] = obj.hertz_fit(smoothx,smoothy,j,1,'parabolic');
                    Rsquare(j) = gof.rsquare;
                    E(j) = emod;
                end
                Rsquare(Rsquare==2) = min(Rsquare);
                obj.GoF{i} = normalize(Rsquare,'range');
                %deltaE(floor(0.5*length(obj.BasedApp{i})):(length(obj.BasedApp{i})-6)) = -diff(log(E));
                close(h)
            end
            current = what();
            cd(obj.Folder)
            savename = sprintf('%s.mat',obj.Name);
            save(savename,'obj')
            cd(current.path)
        end
        
        function cp_combine(obj)
            Range = find(obj.SelectedCurves);
            for i=Range'
                obj.CPCombo{i} = obj.RoV{i}.*obj.GoF{i};
                [~,obj.CP(i)] = max(obj.CPCombo{i});
            end
            current = what();
            cd(obj.Folder)
            savename = sprintf('%s.mat',obj.Name);
            save(savename,'obj')
            cd(current.path)
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
        
        function manual_cp(obj)
            % Manual contact point selection for NN training on plotted force curves
            % returning a position vector in meters.
            jRange = find(obj.SelectedCurves);
            fig = figure('Name',obj.Name);
            for j=jRange'
                fig.WindowState = 'fullscreen';
                %                 Btn_disc=uicontrol('Parent',fig,...
                %                     'Style','pushbutton',...
                %                     'String','Discard Curve',...
                %                     'Units','normalized',...
                %                     'Position',[0.9 0.95 0.14 0.05],...
                %                     'Visible','on',...
                %                     'Callback',{@obj.disc,j,obj});
                plot(obj.THApp{j},obj.BasedApp{j});
                plottitle = sprintf('Curve Nr.%i/%i\n Click and drag the point to the contact point\n Confirm with any key press',j,obj.NCurves);
                title(plottitle);
                [~, domainidx] = ForceMap.no_contact_domain(obj.App{j});
                axis([obj.THApp{j}(floor(domainidx*0.2)) inf -inf inf])
                CP_point = drawpoint();
                %                 w = 0;
                %                 while w == 0
                %                     w = waitforbuttonpress;
                %                 end
                obj.Man_CP(j,1) = CP_point.Position(1);
                obj.Man_CP(j,2) = CP_point.Position(2);
            end
            close(fig);
            current = what();
            cd(obj.Folder)
            savename = sprintf('%s.mat',obj.Name);
            save(savename,'obj')
            cd(current.path)
        end
        
        function old_cp(obj)
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
                [obj.LoadOld{i},obj.UnloadOld{i},Position,obj.CP_old(i,2)] = ContactPoint_sort(load,unload);
                obj.CP_old(i,1) = obj.THApp{i}(Position);
            end
            current = what();
            cd(obj.Folder)
            savename = sprintf('%s.mat',obj.Name);
            save(savename,'obj')
            cd(current.path)
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
            for i=1:obj.Header{6,2}
                for j=1:obj.Header{5,2}
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
        
        function level_height_map(obj)
            % Do a plane fit on the elevationwise lower 40%-quantile of the height_map
            % and subtract the plane from the map. Only works properly on
            % force maps, where 40% or more of the curves are from the glass slide
            % and the plane tilt isn't too big to begin with.
            % (the latter conditions will in practice almost always be fulfilled)
            
            % create a mask that contains the lower 40%-quantile
            pth_quantile = 0.5;
            q = quantile(obj.HeightMap(:,:,1),pth_quantile,'All');
            mask = zeros(size(obj.HeightMap,[1,2]));
            for i=1:obj.Header{6,2}
                for j=1:obj.Header{5,2}
                    if obj.HeightMap(i,j,1) < q
                        mask(i,j,1) = 1;
                    end
                end
            end
            masked_map = mask(:,:,1).*obj.HeightMap(:,:,1);
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
            [Norm,~,Point] = affine_fit(HghtVctr);
            Plane = zeros(size(masked_map,1),size(masked_map,2));
            % Create the plane that can then be subtracted from the
            % complete height data to generate the leveled height data.
            for i=1:size(masked_map,1)
                for j=1:size(masked_map,2)
                    Plane(i,j) = (Point(3)-Norm(1)/Norm(3)*(size(masked_map,2)/size(masked_map,1)*i)-Norm(2)/Norm(3)*(size(masked_map,1)/size(masked_map,2)*j));
                end
            end
            obj.HeightMap(:,:,1) = obj.HeightMap(:,:,1) - Plane;
        end
        
        function save(obj)
            current = what();
            cd(obj.Folder)
            savename = sprintf('%s.mat',obj.Name);
            save(savename,'obj')
            cd(current.path)
            savemsg = sprintf('Changes to ForceMap:%s saved to %s',obj.Name,obj.Folder);
            disp(savemsg);
        end
        
        function cp_cnn_predict(obj,String,NumPasses)
            % String = 'Fast' just predict with one forwardpass using
            % CP_CNN
            % String = 'Dropout' predict through multiple forwardpasses
            % through DropoutNet and additionally gain a metric for
            % uncertainty of the CP estimate
            if nargin < 2
                runmode = 0;
            elseif isequal(String,'Fast')
                runmode = 0;
            elseif isequal(String,'Dropout')
                runmode = 1;
                if nargin < 3
                    NumPasses = 100; % if not specified in arguments, NumPasses defaults to 100
                end
            end
            if isempty(obj.CP_CNN)
                try
                    temp = load('CP_CNN');
                    obj.CP_CNN = temp.CP_CNN;
                catch ME1
                    disp("Can't find Neural Network named 'CP_CNN' in Folder. Trying to load from Workspace instead...")
                    try
                        obj.CP_CNN = CP_CNN;
                    catch ME2
                        disp("Can't find Neural Network named 'CP_CNN' in Workspace. Can't compute Contact Points")
                        rethrow(ME2);
                    end
                end
            end
            if isempty(obj.DropoutNet) && runmode==1
                try
                    temp = load('DropoutNet.mat');
                    obj.DropoutNet = temp.DropoutNet;
                catch ME1
                    disp("Can't find Neural Network named 'CP_CNN' in Folder. Trying to load from Workspace instead...")
                    try
                        obj.CP_CNN = CP_CNN;
                    catch ME2
                        disp("Can't find Neural Network named 'CP_CNN' in Workspace. Can't compute Contact Points")
                        rethrow(ME2);
                    end
                end
            end
            ImgSize = obj.CP_CNN.Layers(1).InputSize;
            objcell{1,1} = obj;
            X = CP_CNN_batchprep_alt(objcell,ImgSize(1));
            len = length(X);
            h = waitbar(0,'Setting up','Name',obj.Name);
            switch runmode
                case 0
                    waitbar(1/2,h,'Predicting CP');
                    Ypredicted = predict(obj.CP_CNN,X);
                case 1
                    obj.YDropPred = zeros(NumPasses,2,len);
                    %                    Ypredicted = zeros(len,2);
                    %                     for i=1:len
                    %                         for j=1:NumPasses
                    %                             waitbar(j/NumPasses,h,sprintf('Predicting CP for curve %i/%i',i,len));
                    %                             obj.YDropPred(j,:,i) = predict(obj.DropoutNet,X(:,:,:,i));
                    %                         end
                    %                         Ypredicted(i,:) = [mean(obj.YDropPred(:,1,i)) mean(obj.YDropPred(:,2,i))];
                    %                     end
                    Ypredicted = zeros(len,2);
                    CantHandle = false;
                    for j=1:NumPasses
                        waitbar(j/NumPasses,h,sprintf('Predicting CP for %i curves. %i/%i passes done',len,j,NumPasses));
                        if CantHandle == false
                            try
                                Temp = predict(obj.DropoutNet,X,'MiniBatchSize',len);
                            catch
                                CantHandle = true;
                                warning("Can't compute with optimal MiniBatchSize. Choosing a lower one")
                                try
                                    Temp = predict(obj.DropoutNet,X);
                                catch
                                    warning("Can't compute with lowered MiniBatchSize. Choosing an even lower one \nComputation time is gonna suffer")
                                    Temp = predict(obj.DropoutNet,X,'MiniBatchSize',10);
                                end
                            end
                        else
                            try
                                Temp = predict(obj.DropoutNet,X);
                            catch
                                warning("Can't compute with lowered MiniBatchSize. Choosing an even lower one \nComputation time is gonna suffer")
                                Temp = predict(obj.DropoutNet,X,'MiniBatchSize',10);
                            end
                        end
                        obj.YDropPred(j,:,:) = Temp';
                    end
                    for i=1:len
                        Ypredicted(i,:) = [mean(obj.YDropPred(:,1,i)) mean(obj.YDropPred(:,2,i))];
                    end
            end
            waitbar(1,h,'Wrapping up');
            iRange = find(obj.SelectedCurves);
            k = 1;
            for i=iRange'
                obj.CP(i,1) = Ypredicted(k,1)*range(obj.THApp{i})+min(obj.THApp{i});
                obj.CP(i,2) = Ypredicted(k,2)*range(obj.BasedApp{i})+min(obj.BasedApp{i});
                k = k + 1;
            end
            close(h)
            current = what();
            cd(obj.Folder)
            savename = sprintf('%s.mat',obj.Name);
            save(savename,'obj')
            cd(current.path)
        end
        
    end
    
    
    methods (Static)
        % Auxilliary functions
        
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
        
        function disc(src,evt,j,obj)
            % This function is under arrested development
            obj.SelectedCurves(j) = 0;
            disp(sprintf('force curve %i discarded from %s',j,obj.Name));
            obj.btnprss.wait = 1;
        end
        
    end
    
    methods
        % non-static auxilliary methods
        
        function show_height_map(obj)
            title = sprintf('Height Map of %s',obj.Name);
            figure('Name',title);
            imshow(imresize(mat2gray(obj.HeightMap(:,:,1)),[1024 1024]));
        end
    end
    
end