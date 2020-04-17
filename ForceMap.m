classdef ForceMap < handle
    % The force map class represents a single jpk force map file and
    % contains all necessary functions to process the forcecurves
    properties
        name            % name of the force map. taken as the name of the folder, containing the .csv files
        folder          % location of the .csv files of the force map
        header = {}     % header properties taken from the JPK force map container file
        app = {}        % approach force data in Newton
        ret = {}        % retraction force data in Newton
        hhapp = {}      % head-height approach data in meters
        hhret = {}      % head-height retract data in meters
        thapp = {}      % vertical tip height approach data in meters
        thret = {}      % vertical tip height retract data in meters
        basedapp = {}   % approach force data with subtracted base line and tilt in Newton
        basedret = {}   % retraction force data with subtracted base line and tilt in Newton
        basefit = {}    % fit model used for the baseline fit
        height_map      % height profile map taken from the inverted maximum head-height from approach max(hhapp)
                        % the additional dimension in the height map holds
                        % the index from the corresponding 1-D lists of
                        % force curves, so one can easily identify the
                        % right curves
        sensitivity
        n_curves        % number of curves on the force map
        pix_app
        pix_ret
        selected_curves % logical vector of length n_curves with 0s for excluded and 1s for included curves. gets initialized with ones
        RoV = {}        % ratio of variance curve for CP estiamation
        RoVTH = {}      % ratio of variance curve, based on the vertical tip position. [does not work well/experimental development stage]
        GoF = {}        % goodness of fit curve for each selected curve
        CPCombo = {}    % combination of the various metrics for contact point estimation
        CP              % chosen contact point (CP) for every selected curve, which is then used in E module fitting
        Man_CP          % manually chosen contact point
        deltaE = {}     % 
        tip_radius = -1 % (nomianl) tip radius for the chosen type of tip model
        poisson_r = 0.5 % standard Poisson ratio for most mechanical models
        btnprss = struct('disc',0,...    % struct to handle buttonpresses that need to communicate outside of their own callback functions
                        'prev',0,...
                        'save',0,...
                        'wait',0)
    end
    methods
        
        function obj = ForceMap()
            %%% Constructor of the class
            
            % Specify the folder where the files live. And import them.
            % Also get curent folder and return to it after import of
            % files.
            % Assigns the properties: name, folder, header, app, ret,
            % hhapp, hhret, height_map, sensitivity, n_curves, pix_app,
            % pix_ret and initiates selected_curves with ones (all curves chosen)
           
            current = what();
            
            quest = 'Do you want to load an already existing .mat file of the force map?';
            answer = questdlg(quest,'Load map...','...from .mat-file','...from .cvs-folder','...from .cvs-folder');
            if isequal(answer, 'Yes')
                [mapfile, mapfilepath] = uigetfile();
                cd(mapfilepath);
                load(mapfile,'-mat');
                cd(current.path);
                msg = sprintf('loading %s',obj.name);
                display(msg);
                disp('loading successfull')
                return
            else
            end
            
            obj.folder = uigetdir;
            cd(obj.folder);
            splitfolder = strsplit(obj.folder,'\');
            obj.name = string(splitfolder(end));
            msg = sprintf('loading %s',obj.name);
            disp(msg);
            
            % Get a list of all files in the folder with the desired file name pattern.
            % Change to whatever pattern you used for the names of your force map files.
            filePattern = fullfile(obj.folder, '*.csv');
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
                    obj.header{k,1} = import_header{i};
                    k = k + 1;
                end
            end
            for i=1:length(obj.header)
                split = strsplit(obj.header{i,1},",");
                if length(split) > 1
                    obj.header{i,2} = str2double(split{2});
                    obj.header{i,1} = split{1};
                else
                end
            end
            % Apply scaling multipliers from header to curve data and
            % delete placeholder data points from curves where piezo range
            % was too small
            TempApp = readmatrix(theFiles(1).name);
            TempHHApp = readmatrix(theFiles(3).name)*obj.header{11,2};
            obj.n_curves = obj.header{5,2}*obj.header{6,2};
            
            for i=1:obj.n_curves
                DelNApp = round(TempApp(end,i),6,'significant');
                if DelNApp == round(TempApp(end-1,i),6,'significant')
                    DelMap = DelNApp;
                    k = 1;
                    while DelNApp == round(TempApp(end-k,i),6,'significant')
                        k = k + 1;
                    end
                    obj.app{i} = TempApp(1:(end-k),i);
                else
                    obj.app{i} = TempApp(:,i);
                end
                obj.hhapp{i} = TempHHApp(1:length(obj.app{i}),i);
            end
            
            for i=1:obj.n_curves
                if round(obj.app{i}(end),6,'significant') == DelMap
                    obj.app{i}(end) = [];
                    obj.hhapp{i}(end) = [];
                end
            end
            
            obj.hhret = mat2cell(readmatrix(theFiles(4).name)*obj.header{11,2},...
                [obj.header{3,2}],ones(obj.header{5,2}*obj.header{6,2},1));
            obj.ret = mat2cell(readmatrix(theFiles(5).name),...
                [obj.header{3,2}],ones(obj.header{5,2}*obj.header{6,2},1));
            obj.pix_app = obj.header{2,2};
            obj.pix_ret = obj.header{3,2};
            obj.sensitivity = obj.header{16,2};
            obj.selected_curves = ones(obj.n_curves,1);
            for i=1:obj.header{6,2}
                for j=1:obj.header{5,2}
                    obj.height_map(i,j,1) = -max(obj.hhapp{mod(i,2)*j+(1-mod(i,2))*(obj.header{5,2}-(j-1))+(obj.header{6,2}-i)*obj.header{5,2}});
                    obj.height_map(i,j,2) = (mod(i,2)*j+(1-mod(i,2))*(obj.header{5,2}-(j-1))+(obj.header{6,2}-i)*obj.header{5,2});
                end
            end
            obj.Man_CP = zeros(obj.n_curves,2);
            cd(current.path);
            current = what();
            cd(obj.folder)
            savename = sprintf('%s.mat',obj.name);
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
            for i=1:obj.n_curves
                dlgtitle = sprintf('Curve selection %i/%i',i,obj.n_curves);
                plottitle = sprintf('Curve Nr. %i',i);
                plot(obj.hhapp{i},obj.app{i},obj.hhret{i},obj.ret{i});
                title(plottitle);
                legend('Approach','Retract');
                answ = questdlg('Do you want to process this indentation curve?',...
                    dlgtitle,'Keep it','Exclude it','Abort Process',...
                    'Abort Process');
                if isequal(answ,'Keep it')
                elseif isequal(answ,'Abort Process')
                    close(gcf)
                    current = what();
                    cd(obj.folder)
                    savename = sprintf('%s.mat',obj.name);
                    save(savename,'obj')
                    cd(current.path)
                    return
                elseif isequal(answ,'Exclude it')
                    obj.selected_curves(i) = 0;
                end
            end
            close(gcf)
            current = what();
            cd(obj.folder)
            savename = sprintf('%s.mat',obj.name);
            save(savename,'obj')
            cd(current.path)
        end
        
        function base_and_tilt(obj)
            % subtract baseline and tilt from the forcecurve by fitting a function to
            % the non contact domain the function tries to fit a non-affine-linear
            % function. If the linear fit is too bad, the function tries a 9th grade
            % polynomial fit instead
            Range = find(obj.selected_curves);
            h = waitbar(0,'Setting up...');
            for i=Range'
                prog = i/obj.n_curves;
                waitbar(prog,h,'processing baseline fits...');
                [ncd , ncidx] = ForceMap.no_contact_domain(obj.app{i});
                gof.rsquare = 0;
                for j=1:10
                    [testfit,goftest] = fit(obj.hhapp{i}(1:ncidx),smooth(ncd),'poly1','Normalize','on');
                    if goftest.rsquare>gof.rsquare
                        gof = goftest;
                        obj.basefit{i} = testfit;
                    end
                end
                
                if gof.rsquare < 0.9
                    %ncd(ncidx+1:length(obj.app{i})) = feval(obj.basefit{i},obj.hhapp{i}(ncidx+1:length(obj.hhapp{i})));
                    ext_range = obj.hhapp{i};
                    ext_l = length(obj.hhapp{i});
                    rangefit = fit([1:length(obj.hhapp{i})]',obj.hhapp{i},'poly1','Normalize','on');
                    ext_range(ext_l+1:floor(1.01*ext_l)) = feval(rangefit,[ext_l+1:floor(1.01*ext_l)]');
                    ncd(ncidx+1:length(ext_range)) = feval(obj.basefit{i},ext_range(ncidx+1:length(ext_range)));
                    [obj.basefit{i},gof] =  fit(ext_range,smooth(ncd),'poly9','Normalize','on');
                end
%                 if gof.rsquare < 0.5
%                     f = figure('Name','Base and tilt','Position',[10000 10000 800 600]);
%                     movegui(f);
%                     plot(obj.hhapp{i},feval(obj.basefit{i},obj.hhapp{i}),obj.hhapp{i},obj.app{i});
%                     legend('Fit','Data');
%                     plottitle = sprintf('Curve Nr. %i',i);
%                     title(plottitle);
%                     answ = questdlg('Quality of base line fit is very low. consider excluding this force curve',...
%                         'Base line and tilt','Exclude it',...
%                         'Keep current fit','Keep current fit');
%                     if isequal(answ,'Exclude it')
%                         obj.selected_curves(i) = 0;
%                         continue
%                     else
%                         return
%                     end
%                     close(gcf)
%                 end
                obj.basedapp{i} = (obj.app{i}-feval(obj.basefit{i},obj.hhapp{i}));
                obj.basedret{i} = (obj.ret{i}-feval(obj.basefit{i},obj.hhret{i}));
            end
            % calculate vertical tip position by subtracting vertical tip deflection from head height
            iRange = find(obj.selected_curves);
            for i=iRange'
                obj.thapp{i} = obj.hhapp{i} - obj.basedapp{i}/obj.header{17,2};
                obj.thret{i} = obj.hhret{i} - obj.basedret{i}/obj.header{17,2};
            end
            close(h);
            current = what();
            cd(obj.folder)
            savename = sprintf('%s.mat',obj.name);
            save(savename,'obj')
            cd(current.path)
        end
        
        function CP_RoV(obj,batchsize)
            % find contact point with the method of ratio of variances. The method
            % iterates though every point and builds the ratio of the variance of a
            % bunch of points before and after the current point. the point with the
            % biggest ratio is the returned contact point [Nuria Gavara, 2016]
            if nargin<2
                batchsize = 10;
            end
            Range = find(obj.selected_curves);
            h = waitbar(0,'Setting up...');
            for i=Range'
                prog = i/obj.n_curves;
                waitbar(prog,h,'applying ratio of variances method...');
                obj.RoV{i} = zeros(length(obj.basedapp{i}),1);
                
                % loop through points and calculate the RoV
                for j=(batchsize+1):(length(obj.basedapp{i})-batchsize)
                    obj.RoV{i}(j,1) = var(smoothdata(obj.basedapp{i}((j+1):(j+batchsize))))/...
                        var(smoothdata(obj.basedapp{i}((j-batchsize):(j-1))));
                end
                % normalize RoV-curve
                obj.RoV{i} = obj.RoV{i}/range(obj.RoV{i});
                minrov = min(obj.RoV{i}(batchsize+1:length(obj.RoV{i})-batchsize));
                obj.RoV{i}(obj.RoV{i}==0) = minrov;
            end
            close(h)
            current = what();
            cd(obj.folder)
            savename = sprintf('%s.mat',obj.name);
            save(savename,'obj')
            cd(current.path)
        end
        
        function [E_mod,GoF,Hertzfit] = Hertz_fit(obj,tip_h,force,CP,curve_percent,shape)
            
            if obj.tip_radius == -1
                prompt = {'What is the nominal tip radius of the used cantilever in nm?'};
                dlgtitle = 'Cantilever tip';
                dims = [1 35];
                definput = {'10'};
                obj.tip_radius = str2double(inputdlg(prompt,dlgtitle,dims,definput));
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
                E_mod = 3*(Hertzfit.a*ranf/rant^(3/2))/(4*sqrt(obj.tip_radius*10^(-9)))*(1-obj.poisson_r^2);
            elseif isequal(shape,'spherical')
            elseif isequal(shape,'conical')
            elseif isequal(shape,'pyramid')
            end
            current = what();
            cd(obj.folder)
            savename = sprintf('%s.mat',obj.name);
            save(savename,'obj')
            cd(current.path)
        end
        
        function CP_GoF_deltaE(obj)
            Range = find(obj.selected_curves);
            for i=Range'
                smoothx =obj.thapp{i};
                smoothy =obj.basedapp{i};
                Rsquare = 2*ones(length(obj.basedapp{i}),1);
                E = ones(length(obj.basedapp{i}),1);
                deltaE = ones(length(obj.basedapp{i}),1);
                testrange = floor(0.5*length(obj.basedapp{i})):(length(obj.basedapp{i})-5);
                h = waitbar(0,'Setting up...');
                msg = sprintf('applying goodness of fit method on curve Nr.%i/%i',i,obj.n_curves);
                for j=testrange
                    prog = (j-testrange(1))/length(testrange);
                    waitbar(prog,h,msg);
                    [emod,gof] = obj.Hertz_fit(smoothx,smoothy,j,1,'parabolic');
                    Rsquare(j) = gof.rsquare;
                    E(j) = emod;
                end
                Rsquare(Rsquare==2) = min(Rsquare);
                obj.GoF{i} = normalize(Rsquare,'range');
                %deltaE(floor(0.5*length(obj.basedapp{i})):(length(obj.basedapp{i})-6)) = -diff(log(E));
                close(h)
            end
            current = what();
            cd(obj.folder)
            savename = sprintf('%s.mat',obj.name);
            save(savename,'obj')
            cd(current.path)
        end
        
        function CP_combined(obj)
            Range = find(obj.selected_curves);
            for i=Range'
                obj.CPCombo{i} = obj.RoV{i}.*obj.GoF{i};
                [~,obj.CP(i)] = max(obj.CPCombo{i});
            end
            current = what();
            cd(obj.folder)
            savename = sprintf('%s.mat',obj.name);
            save(savename,'obj')
            cd(current.path)
        end
        
        function CP_RoVTH(obj,batchsize)
            % find contact point with the method of ratio of variances. The method
            % iterates though every point and builds the ratio of the variance of a
            % bunch of points before and after the current point. the point with the
            % biggest ratio is the returned contact point [Nuria Gavara, 2016]
            if nargin<2
                batchsize = 10;
            end
            Range = find(obj.selected_curves);
            h = waitbar(0,'Setting up...');
            for i=Range'
                prog = i/obj.n_curves;
                waitbar(prog,h,'applying ratio of variances method...');
                obj.RoVTH{i} = zeros(length(obj.thapp{i}),1);
                
                % loop through points and calculate the RoV
                for j=(batchsize+1):(length(obj.thapp{i})-batchsize)
                    obj.RoVTH{i}(j,1) = var(smoothdata(obj.thapp{i}((j+1):(j+batchsize))))/...
                        var(smoothdata(obj.thapp{i}((j-batchsize):(j-1))));
                end
                % normalize RoVTH-curve
                obj.RoVTH{i} = obj.RoVTH{i}/range(obj.RoVTH{i});
                minrov = min(obj.RoVTH{i}(batchsize+1:length(obj.RoVTH{i})-batchsize));
                obj.RoVTH{i}(obj.RoVTH{i}==0) = minrov;
            end
            close(h)
            current = what();
            cd(obj.folder)
            savename = sprintf('%s.mat',obj.name);
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
            Range = find(obj.selected_curves);
            imres = obj.header{2,2};
            k = 1;
            img = cell(sum(obj.selected_curves),1);
            imgorsize = cell(sum(obj.selected_curves),1);
            h = waitbar(0,'Setting up...');
            for i=Range'
                prog = i/obj.n_curves;
                waitbar(prog,h,'Converting force curves to images...');
                norm =  round(normalize(obj.basedapp{i},'range',[1 imres]));
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
        
        function manual_CP(obj)
            % Manual contact point selection for NN training on plotted force curves
            % returning a position vector in meters.
            jRange = find(obj.selected_curves);
            fig = figure('Name',obj.name);
            for j=jRange'
                fig.WindowState = 'fullscreen';
%                 Btn_disc=uicontrol('Parent',fig,...
%                     'Style','pushbutton',...
%                     'String','Discard Curve',...
%                     'Units','normalized',...
%                     'Position',[0.9 0.95 0.14 0.05],...
%                     'Visible','on',...
%                     'Callback',{@obj.disc,j,obj});
                plot(obj.thapp{j},obj.basedapp{j});
                plottitle = sprintf('Curve Nr.%i/%i\n Click and drag the point to the contact point\n Confirm with any key press',j,obj.n_curves);
                title(plottitle);
                [~, domainidx] = ForceMap.no_contact_domain(obj.app{j});
                axis([obj.thapp{j}(floor(domainidx*0.2)) inf -inf inf])
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
            cd(obj.folder)
            savename = sprintf('%s.mat',obj.name);
            save(savename,'obj')
            cd(current.path)
        end
        
        function choose_fibril(obj,pth_quantile)
           % sets every forcecurve below a certain threshold to zero in
           % obj.selected_curves. Note that curves on the fibril already
           % chosen to be disregarded will keep that status!.This function
           % is a crude way to find fibril areas and should not be used, if
           % its important to select exclusively indentations on the fibril
           % or ,in the more general case,other elevated region of interest
           if nargin < 2
               pth_quantile = 0.8;
           end
           q = quantile(obj.height_map(:,:,1),pth_quantile,'All');
           mask = ones(size(obj.height_map));
           for i=1:obj.header{6,2}
               for j=1:obj.header{5,2}
                   if obj.height_map(i,j,1) < q
                       obj.selected_curves(obj.height_map(i,j,2)) = 0;
                       mask(i,j,1) = 0;
                   end
               end
           end
           f = figure();
           imshowpair(imresize(mat2gray(obj.height_map(:,:,1)),[1024 1024]),imresize(mask(:,:,1),[1024 1024]),'montage')
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
            q = quantile(obj.height_map(:,:,1),pth_quantile,'All');
            mask = zeros(size(obj.height_map,[1,2]));
            for i=1:obj.header{6,2}
                for j=1:obj.header{5,2}
                    if obj.height_map(i,j,1) < q
                        mask(i,j,1) = 1;
                    end
                end
            end
            masked_map = mask(:,:,1).*obj.height_map(:,:,1);
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
            obj.height_map(:,:,1) = obj.height_map(:,:,1) - Plane;
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
            obj.selected_curves(j) = 0;
            disp(sprintf('force curve %i discarded from %s',j,obj.name));
            obj.btnprss.wait = 1;
        end
        
    end
    
end