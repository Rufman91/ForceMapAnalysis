classdef SurfacePotentialMap < matlab.mixin.Copyable
    
    properties
        Name
        Folder
        NumPoints
        NumProfiles
        XScale
        YScale
        Height
        Potential
        HeightBased
        FibMask
        FPotMask
        GPotMask
        PotExclMask
        Apex
        ApexIndex
        RectApex
        RectApexIndex
        FibDiam
        DBanding
        FibPot
        FibPotSTD
    end
    
    methods
        
        function obj = SurfacePotentialMap(mapfilepath,mapname)
            
            %%% Constructor of the class
            
            % Specify the folder where the files live. And import them.
            % Also get curent folder and return to it after import of
            % files.
            
            current = what();
            
            % determine if SurfacePotentialMap is given a loadpath for existing .mat
            if nargin > 0
                cd(mapfilepath);
                file = dir(sprintf('%s.mat',mapname));
                load(file.name);
                cd(current.path);
                msg = sprintf('loading %s',obj.Name);
                disp(msg);
                disp('loading successfull')
                return
            end
            
            quest = 'Do you want to load an already existing .mat file of the surface potential map?';
            answer = questdlg(quest,'Load map...','...from .mat-file','...from .sdf-file','...from .sdf-file');
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
            
            
            [TempName , obj.Folder] = uigetfile('*.sdf','Select 2 files: Height and Potential','MultiSelect','on');
            obj.Name = TempName{1}(1:(end-6));
            cd(obj.Folder);
            msg = sprintf('loading %s',obj.Name);
            disp(msg);
            
            TempP = importdata(TempName{2},'',16);
            TempH = importdata(TempName{1},'',16);
            
            % Create an array the size of the image for the potential data
            obj.NumPoints = sscanf(TempP.textdata{5,1},'NumPoints   = %d');
            obj.NumProfiles = sscanf(TempP.textdata{6,1},'NumProfiles   = %d');
            obj.YScale = sscanf(TempP.textdata{8,1},'Yscale      = %e');
            obj.XScale = sscanf(TempP.textdata{7,1},'Xscale      = %e');
            obj.Potential = zeros(obj.NumProfiles,obj.NumPoints);
            for i = 1:obj.NumProfiles
                for j = 1:obj.NumPoints
                    obj.Potential(i,j) = TempP.data(obj.NumPoints*(i-1)+j);
                end
            end
            
            % Create an array the size of the image for the height data
            obj.Height = zeros(obj.NumProfiles,obj.NumPoints);
            for i = 1:obj.NumProfiles
                for j = 1:obj.NumPoints
                    obj.Height(i,j) = TempH.data(obj.NumPoints*(i-1)+j);
                end
            end
            
            cd(current.path)
            cd(obj.Folder)
            savename = sprintf('%s.mat',obj.Name);
            save(savename,'obj')
            cd(current.path)
            disp('loading successfull. object saved in objects folder')
        end
        
        function level_height(obj,MaskParam)
            if nargin < 2
                MaskParam = 0.5;
            end
            HghtRange = range(obj.Height,'all');
            HghtMin = min(obj.Height,[],'all');
            HghtNorm = (obj.Height-HghtMin)/HghtRange;
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
            
%             subplot(2,2,1)
%             imshow(mask)
%             subplot(2,2,2)
%             imshow(HghtNorm)
%             subplot(2,2,3)
%             imshow(bwareafilt(mask,1,4))
%             subplot(2,2,4)
%             surf(HghtNorm.*mask)
            X=zeros(sum(mask,'All'),3);
            k = 1;
            for i=1:obj.NumProfiles
                for j=1:obj.NumPoints
                    if mask(i,j) == 0
                        X(k,:)=[obj.XScale*i obj.YScale*j obj.Height(i,j)];
                        k = k + 1;
                    end
                end
            end
            [Norm,~,Point] = obj.affine_fit(X);
            Plane = zeros(obj.NumProfiles,obj.NumPoints);
            % Create the plane that can then be subtracted from the
            % complete height data to generate the leveled height data.
            for i=1:obj.NumProfiles
                for j=1:obj.NumPoints
                    Plane(i,j) = (Point(3)-Norm(1)/Norm(3)*(obj.XScale*i)-Norm(2)/Norm(3)*(obj.YScale*j));
                end
            end
            obj.HeightBased = obj.Height - Plane;
            
            current = what();
            cd(obj.Folder)
            savename = sprintf('%s.mat',obj.Name);
            save(savename,'obj')
            cd(current.path)
        end
        
        function create_fibril_mask(obj,MaskParam)
            
            if nargin < 2
                MaskParam = 0.5;
            end
            HghtRange = range(obj.HeightBased,'all');
            HghtMin = min(obj.HeightBased,[],'all');
            HghtNorm = (obj.HeightBased-HghtMin)/HghtRange;
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
            for i=1:100
                mask = medfilt2(mask,[11 5],'symmetric');
            end
            mask = bwareafilt(mask,1,4);
            obj.FibMask = mask;
            
            current = what();
            cd(obj.Folder)
            savename = sprintf('%s.mat',obj.Name);
            save(savename,'obj')
            cd(current.path)
        end
        
        function calc_fib_diam_dband(obj)
            [obj.Apex,obj.ApexIndex] = max(obj.HeightBased.*obj.FibMask,[],2);
            obj.RectApex = zeros(obj.NumProfiles,1);
            obj.RectApexIndex = zeros(obj.NumProfiles,1);
            obj.RectApexIndex = round(predictGP_mean([1:obj.NumProfiles],[1:obj.NumProfiles],1,8*obj.NumProfiles,obj.ApexIndex,1));
            for i=1:obj.NumProfiles
                obj.RectApex(i) = obj.HeightBased(i,obj.RectApexIndex(i));
            end
            X = [-(obj.NumProfiles/2):(obj.NumProfiles/2-1)]/(obj.NumProfiles*obj.YScale);
            Y = (obj.RectApex - mean(obj.RectApex))/range(obj.RectApex);
            FourierTrafo = abs(fftshift(fft(Y)));
%             plot(X,FourierTrafo)
            XPos = X(X>=0);
            YPos = FourierTrafo((length(XPos)+1):end);
            [~,Idx] = sort(YPos,'descend');
            k = XPos(Idx);
            j = 1;
            for i=1:length(k)
                if (k(i) < 1/(57e-9)) && (k(i) > 1/(77e-9))
                    Tempk(j) = k(i);
                    j = j + 1;
                end
            end
            obj.DBanding = 1/Tempk(1);
            k = 1;
            for i=1:obj.NumProfiles
                if obj.PotExclMask(i,obj.RectApexIndex(i)) == 1
                    FibHeight(k) = obj.RectApex(i);
                    k = k + 1;
                end
            end
            obj.FibDiam = mean(FibHeight);
            current = what();
            cd(obj.Folder)
            savename = sprintf('%s.mat',obj.Name);
            save(savename,'obj')
            cd(current.path)
        end
        
        function manual_exclusion(obj)
            
            obj.PotExclMask = logical(ones(obj.NumProfiles,obj.NumPoints));
            CheckSum = 100;
            while CheckSum > 5
                f = figure('Name','Choose areas to be excluded');
                f.WindowState = 'maximized';
                subplot(2,1,2)
                imshow(obj.HeightBased.*obj.PotExclMask,[min(obj.HeightBased,[],'all') max(obj.HeightBased,[],'all')])
                subplot(2,1,1)
                imshow(obj.Potential.*obj.PotExclMask,[min(obj.Potential,[],'all') max(obj.Potential,[],'all')])
                title(sprintf('%s: Draw Freehand ROI around areas, that are to be excluded\nThe area will be taken out and the same map redrawn \n If there is nothing to do just click on the image once without dragging the cursor',obj.Name))
                ROI = drawfreehand;
                CheckSum = length(ROI.Waypoints);
                Mask = ~createMask(ROI);
                obj.PotExclMask = obj.PotExclMask.*Mask;
                close(f)
            end
            
            current = what();
            cd(obj.Folder)
            savename = sprintf('%s.mat',obj.Name);
            save(savename,'obj')
            cd(current.path)
        end
        
        function calc_glass_fib_potdiff(obj,Percent)
            % Fibril potential is everything thats +-Percent*FibDiam from the
            % Apex minus points that have been excluded under the method
            % manual_exclusion. Glass potential is everything starting
            % 3*Percent*obj.FibDiam away from the border of obj.FibMask minus
            % points that have been excluded under method manual_exclusion.
            if nargin < 2
                Percent = 1;
            end
            
            % Create Mask for Fibril Potential
            obj.FPotMask = zeros(obj.NumProfiles,obj.NumPoints);
                Range = floor(Percent*obj.FibDiam/obj.YScale);
            for i=1:obj.NumProfiles
                if (obj.RectApexIndex(i) - Range) < 1
                    obj.FPotMask(i,[obj.RectApexIndex(i):...
                         (obj.RectApexIndex(i) + Range)]) = 1;
                elseif (obj.RectApexIndex(i) + Range) > obj.NumPoints
                    obj.FPotMask(i,[(obj.RectApexIndex(i) - Range):...
                         obj.RectApexIndex(i)]) = 1;
                else
                    obj.FPotMask(i,[(obj.RectApexIndex(i) - Range):...
                         (obj.RectApexIndex(i) + Range)]) = 1;
                end
            end
            obj.FPotMask = obj.FPotMask.*obj.PotExclMask;
            
            % Create Mask for Glass Potential
            obj.GPotMask = zeros(obj.NumProfiles,obj.NumPoints);
            for i=1:obj.NumProfiles
                j = 1;
                k = obj.NumPoints;
                while obj.FibMask(i,j+3*Range) == 0
                    obj.GPotMask(i,j) = 1;
                    j = j + 1;
                end
                while obj.FibMask(i,k-3*Range) == 0
                    obj.GPotMask(i,k) = 1;
                    k = k - 1;
                end
            end
            obj.GPotMask = obj.GPotMask.*obj.PotExclMask;
            
            %Calculate Fibril Surface Potential
            Fib = obj.Potential.*obj.FPotMask;
            Fib = Fib(Fib~=0);
            Glass = obj.Potential.*obj.GPotMask;
            Glass = Glass(Glass~=0);
            obj.FibPot = mean(Fib,'all') - mean(Glass,'all');
            % According to Gaussian propagation of uncertainty
            obj.FibPotSTD = sqrt(std(Fib,0,'all')^2 + std(Glass,0,'all')^2);
        end
        
    end
    
    methods
        % Auxiliary methods
        
        function save(obj)
            current = what();
            cd(obj.Folder)
            savename = sprintf('%s.mat',obj.Name);
            save(savename,'obj')
            cd(current.path)
            savemsg = sprintf('Changes to SurfacePotentialMap %s saved to %s',obj.Name,obj.Folder);
            disp(savemsg);
        end
        
        function show_height(obj,ResMult)
            if nargin < 2
                ResMult = 1;
            end
            subplot(2,1,1)
            I = imresize(obj.Height,ResMult*[obj.NumProfiles*obj.YScale/obj.XScale obj.NumPoints]);
            imshow(I,[min(I,[],'all') max(I,[],'all')])
            title(sprintf('%s Height',obj.Name))
            subplot(2,1,2)
            surf(I,'LineStyle','none','FaceLighting','gouraud','FaceColor','interp')
            light('Style','local')
        end
        
        function show_heightbased(obj,ResMult)
            if nargin < 2
                ResMult = 1;
            end
            I = imresize(obj.HeightBased ,ResMult*[obj.NumProfiles*obj.YScale/obj.XScale obj.NumPoints]);
            subplot(2,1,1)
            imshow(I,[min(I,[],'all') max(I,[],'all')])
            title(sprintf('%s Plane Fitted Height',obj.Name))
            subplot(2,1,2)
            surf(I,'LineStyle','none','FaceLighting','gouraud','FaceColor','interp')
            light('Style','local')
        end
        
        function show_potential(obj,ResMult)
            if nargin < 2
                ResMult = 1;
            end
            I = imresize(obj.Potential ,ResMult*[obj.NumProfiles*obj.YScale/obj.XScale obj.NumPoints]);
            subplot(2,1,1)
            imshow(I,[min(I,[],'all') max(I,[],'all')])
            title(sprintf('%s Surface Potential',obj.Name))
            subplot(2,1,2)
            surf(I,'LineStyle','none','FaceLighting','gouraud','FaceColor','interp')
            light('Style','local')
        end
        
        function show_summary(obj,ResMult)
            if nargin < 2
                ResMult = 1;
            end
            
            figure('Name',sprintf('Summarized Properties of %s',obj.Name),'Units','normalized','Position',[0.15 0.15 0.7 0.7])
            
            I = imresize(obj.HeightBased ,ResMult*[obj.NumProfiles*obj.YScale/obj.XScale obj.NumPoints]);
            I = I*1e9;
            
            subplot(2,2,1)
            imshow(I,[min(I,[],'all') 1.2*max(I,[],'all')],'Colormap',hot)
            title(sprintf('%s Plane Fitted Height',obj.Name))
            c1 = colorbar;
            c1.Label.String = 'Height [nm]';
            
            subplot(2,2,2)
            surf(imrotate(I',90),'LineStyle','none','FaceLighting','gouraud','FaceColor','interp')
            light('Style','local')
            title({sprintf('Fibril Diameter = %.1fnm',obj.FibDiam*1e9),...
                sprintf('D-Banding = %.1fnm',obj.DBanding*1e9)})
            
            subplot(2,2,3)
            Potential = imresize(obj.Potential,[size(I,1) size(I,2)]);
            imshow(Potential,[min(Potential,[],'all') max(Potential,[],'all')],'Colormap',winter)
            c3 = colorbar;
            c3.Label.String = 'Potential [V]';
            title({sprintf('%s Surface Potential',obj.Name),...
                sprintf('Surface Potential relative to Glass = %.2fmV',obj.FibPot*1e3)});
            
            subplot(2,2,4)
            CombinedMask = obj.FPotMask+obj.GPotMask-~obj.PotExclMask;
            CombinedMask = imresize(CombinedMask,[size(I,1) size(I,2)]);
            imshow(CombinedMask);
            title('Combined Mask of Fibril Potential-, Glass Potential and Excluded Region')
        end
    end
    
    methods(Static)
        % Auxiliary static methods
        
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
        
    end
    
end