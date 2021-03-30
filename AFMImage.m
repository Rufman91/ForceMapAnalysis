classdef AFMImage < matlab.mixin.Copyable
    % This is supposed to be a class for analysis and processing of AFM
    % data at the image level in general and isn't restricted to a specific
    % mode of image acquisition. As such, it should be able to load, access
    % and process images from all AC-modes, contact mode, QI-images, image
    % projections of ForceMaps and so on (maybe subclass SurfacePotentialMap
    % should eventually be incorporated into this).
   
    properties
        % File and Header properties
        Name
        ID
        HostOS
        HostName
        ImagingType
        FileVersion
        DateTime
        NumChannels
    end
    properties
        % Image Data Properties. All dimensions in SI-units, Angles in
        % degrees
        OriginX
        OriginY
        ScanSizeX
        ScanSizeY
        ScanAngle
        NumPixelsX
        NumPixelsY
        IGain
        PGain
        RelativeSetpoint
        DriveAmplitude
        DriveFrequency
        PhaseShift
        LineRate
        TipVelocity
        Sensitivity
        SpringConstant
    end
    properties
        % Additional properties needed during (ErrorSignal)channel readout
        Baseline_Raw
        Bline_adjust
        SetP_V
        Raw
        SetP_m
        SetP_N
        Baseline_N
    end
    properties
        % All possible image channels. The Channels are themselves structs
        % containing some infos and the imagedata in their raw and
        % processed form.
        Height
        HeightMeasured
        ErrorSignal
        LockInAmplitude
        LockInPhase
        LateralDeflection
        VerticalDeflection
        Processed
    end
    properties
        % Properties related to Image processing/segmenting/classification
        CMap
        MaskBackground
    end
    properties
        ErodedTip
        DepthDependendTipRadius
        DepthDependendTipShape
        ProjectedTipArea
    end
    properties
        % All the Flags
        HasHeight
        HasHeightMeasured
        HasErrorSignal
        HasLateralDeflection
        HasLockInAmplitude
        HasLockInPhase
        HasVerticalDeflection
        HasProcessed
        HasDeconvolutedCantileverTip
    end
    
    methods
        % Main methods of the class
        
        function obj = AFMImage(ImageFullFile,DataFolder,TempID)
            % Constructor of the class. Extracts Header properties as well
            % as all available channel-data
            
            if nargin == 0
                [File, Path] = uigetfile('*.jpk','Choose a jpk image file');
                ImageFullFile = fullfile(Path, File);
                TempID = 'AFMImage detached from Experiment-class 1';
            end
            
            obj. CMap = obj.define_afm_color_map(0);
            
            obj.initialize_flags
            
            obj.Name = obj.parse_file_name(ImageFullFile);
            obj.ID = TempID;
            
            % get OS and use appropriate fitting system command
            obj.check_for_new_host
            
            Index = regexp(obj.ID,'(?<=\-).','all');
            LoadMessage = sprintf('loading data into AFMImage Nr.%s',obj.ID(Index(end):end));
            disp(LoadMessage)
            
            obj.read_in_header_properties(ImageFullFile);
            
            obj.load_image_channels(ImageFullFile);
            
        end
        
        function deconvolute_cantilever_tip(obj)
            
            Based = imgaussfilt(obj.subtract_line_fit_hist(obj.HeightMeasured.Trace,obj.NumPixelsY,0.4));
            obj.MaskBackground = obj.mask_background_by_threshold(Based,1);
            Based = obj.masked_plane_fit(Based,obj.MaskBackground);
            ConeHeight = range(Based,'all');
            Cone = obj.cone(obj.NumPixelsX,obj.NumPixelsY,ConeHeight,obj.ScanSizeX,obj.ScanSizeY,10e-9);
            obj.ErodedTip = obj.deconvolute_by_mathematical_morphology(Based,Cone);
            
            StepSize = 1e-9;
            
            PixelSizeX = obj.ScanSizeX/obj.NumPixelsX;
            PixelSizeY = obj.ScanSizeY/obj.NumPixelsY;
            
            obj.ProjectedTipArea = obj.calculate_projected_area(obj.ErodedTip,PixelSizeX,PixelSizeY,StepSize);
            
            [obj.DepthDependendTipRadius,obj.DepthDependendTipShape] = obj.calculate_depth_dependend_tip_data(obj.ProjectedTipArea,RangePercent);
            
            obj.hasDeconvolutedCantileverTip = true;
        end
        
        function base_image()
            
        end
        
    end
    
    methods(Static)
        % Static Main Methods
        
        function OutImage = subtract_line_fit_hist(InImage,CutOff)
            
            NumProfiles = size(InImage,1);
            NumPoints = size(InImage,2);
            CutOff = ceil(CutOff*NumPoints);
            
            for i=1:NumProfiles
                Line = InImage(i,:)';
                [~, SortedIndex] = sort(Line,'ascend');
                LineFit = polyfit(SortedIndex(1:CutOff),Line(SortedIndex(1:CutOff)),1);
                LineEval = [1:NumPoints]'*LineFit(1) + LineFit(2);
                Line = Line - LineEval;
                InImage(i,:) = Line;
            end
            OutImage = InImage;
        end
        
        function subtract_line_fit_automatically(InImage)
            if nargin<2
                Mask = ones(size(InImage));
            end
            
            for i=1:NumProfiles
                Line = InImage(i,:)';
                [~, SortedIndex] = sort(Line,'ascend');
                LineFit = polyfit(SortedIndex(1:CutOff),Line(SortedIndex(1:CutOff)),1);
                LineEval = [1:NumPoints]'*LineFit(1) + LineFit(2);
                Line = Line - LineEval;
                InImage(i,:) = Line;
            end
            OutImage = InImage;
        end
        
        function OutMask = mask_background_by_threshold(Image,PercentOfRange)
            
            if nargin < 2
                PercentOfRange = 5;
            end
            
            Thresh = range(Image,'all')*PercentOfRange/100;
            [Row,Col] = find((abs(Image)<=Thresh));
            OutMask = zeros(size(Image));
            for i=1:length(Row)
                OutMask(Row(i),Col(i)) = 1;
            end
            
        end
        
        function OutImage = masked_plane_fit(Image,Mask)
            
            % Convert Image to Point Cloud for plane fit
            [X,Y,Z] = AFMImage.convert_masked_to_point_cloud(Image,Mask);
            
            [Norm,~,Point] = AFMImage.affine_fit([X Y Z]);
            Plane = zeros(size(Image));
            % Create the plane that can then be subtracted from the
            % complete height data to generate the leveled height data.
            for i=1:size(Image,1)
                for j=1:size(Image,2)
                    Plane(i,j) = (Point(3)-Norm(1)/Norm(3)*i-Norm(2)/Norm(3)*j);
                end
            end
            
            OutImage = Image - Plane;
        end
        
        function ProjectedArea = calculate_projected_area(Image,PixelSizeX,PixelSizeY,StepSize)
            
            MinImage = min(Image,[],'all');
            MaxImage = max(Image,[],'all');
            
            Thresh = MaxImage;
            k = 1;
            while Thresh >= MinImage
                NumPointsInArea = length(find(Image >= Thresh));
                ProjectedArea(k) = NumPointsInArea*PixelSizeX*PixelSizeY;
                Thresh = Thresh -StepSize;
                k = k + 1;
            end
            ProjectedArea = ProjectedArea';
        end
        
        function CMap = define_afm_color_map(PlusBrightness)
            if nargin == 0
                PlusBrightness = 0;
            end
            CMap(:,1) = (0:1/255:1).*2 + PlusBrightness;
            CMap(:,2) = (0:1/255:1).*2 - 0.5 + PlusBrightness;
            CMap(:,3) = (0:1/255:1).*2 - 1 + PlusBrightness;
            CMap(CMap < 0) = 0;
            CMap(CMap > 1) = 1;
        end
        
        function OutImage = deconvolute_by_mathematical_morphology(InImage,ErodingGeometry)
            
            if ~isequal(size(InImage),size(ErodingGeometry))
                error('The Image and the Eroding Geometry need to have the same pixel-dimensions')
            end
            MaxPeakValue = max(ErodingGeometry,[],'all');
            [PeakIndexX,PeakIndexY] = find(ErodingGeometry==MaxPeakValue);
            EGPixelsX = size(ErodingGeometry,1);
            EGPixelsY = size(ErodingGeometry,2);
            InPixelsX = size(InImage,1);
            InPixelsY = size(InImage,2);
            OutImage=ones(InPixelsX,InPixelsY); %creates the empty image array
            
            
            h = waitbar(0,'Please wait, processing deconvolution...');
            
            %loops over all the elements and find the minimum value of w and allocate it
            for j=1:InPixelsY %loops over points in image output
                waitbar(j/InPixelsY);
                s_ymin=max(-j,-PeakIndexY); %determines the allowed range for tip scanning
                s_ymax=min(EGPixelsY-PeakIndexY,InPixelsY-j)-1; %idem
                for i=1:InPixelsX %loops over points in other direction in image output
                    s_xmin=max(-PeakIndexX,-i); %determines allowed range for tip scanning
                    s_xmax=min(EGPixelsX-PeakIndexX,InPixelsX-i)-1; %idem
                    TempMat = InImage((i+1+s_xmin ):(i+1+s_xmax),(j+1+s_ymin):(j+1+s_ymax))-...
                        ErodingGeometry((s_xmin+PeakIndexX+1):(s_xmax+PeakIndexX+1),...
                        +(1+PeakIndexY+s_ymin):(1+PeakIndexY+s_ymax));
                    minimum=min(TempMat,[],'all'); %checks if w is minimum
                    OutImage(i,j)=minimum; %allocates the minimum value
                end
            end
            close(h);
        end
        
        function [DepthDependendTipRadius,DepthDependendTipShape] = calculate_depth_dependend_tip_data(ProjectedTipArea,RangePercent)
            
            if nargin < 2
                RangePercent = 100;
            end
            
            MaxIdx = floor(RangePercent/100*length(ProjectedTipArea));
            ProjectedTipArea = ProjectedTipArea*(1e9)^2;
            DepthDependendTipRadius = zeros(MaxIdx,1);
            DepthDependendTipShape = cell(MaxIdx,1);
            
            % Fit a sphere and a parabola for every depthstep and choose
            % the one with better fit. Start at 5nm ind. depth 
            for i=5:MaxIdx
                % fit projected area of a parabolic tip
                SphOpt = fitoptions('Method','NonlinearLeastSquares',...
                    'Lower',0,...
                    'Upper',inf,...
                    'MaxIter',4000,...
                    'StartPoint',1e-6,...
                    'Normalize','off',...
                    'DiffMaxChange',1e20,...
                    'DiffMinChange',1e-20,...
                    'MaxFunEvals',4000,...
                    'TolFun',1e-20,...
                    'TolX',1e-20);
                ProjAParabola = fittype('pi*a*x',...
                    'dependent',{'y'},'independent',{'x'},...
                    'coefficients',{'a'},...
                    'options',SphOpt);
                Depth = [1:i]';
                [ParabolaFit,GoFParabola] = fit(Depth,...
                    ProjectedTipArea(1:i),...
                    ProjAParabola);
                RParabola = 1/(2*ParabolaFit.a);
                % fit projected area of a spherical tip
                ProjASphere = fittype('(R^2-(R-x)^2)',...
                    'dependent',{'y'},'independent',{'x'},...
                    'coefficients',{'R'},...
                    'options',SphOpt);
                [SphereFit,GoFSphere] = fit(Depth,...
                    ProjectedTipArea(1:i),...
                    ProjASphere);
                RSphere = SphereFit.R;
                if GoFParabola.rmse <= GoFSphere.rmse
                    DepthDependendTipRadius(i) = RParabola;
                    DepthDependendTipShape{i} = 'parabolic';
                else
                    DepthDependendTipRadius(i) = RSphere;
                    DepthDependendTipShape{i} = 'spherical';
                end
                plot(Depth,ProjectedTipArea(1:i),'rO',Depth,feval(SphereFit,Depth),'b',Depth,feval(ParabolaFit,Depth),'g')
                legend({'Proj. A. from Eroded Tip','Proj. A. Spherical Fit','Proj. A. Parabolic Fit'})
                title({'Spherical Fit',sprintf('Radius:%d nm  Depth:%i nm GoF.rmse: %d',SphereFit.R,i,GoFSphere.rmse),...
                    'Parabolic Fit',sprintf('Radius:%d nm  Depth:%i nm GoF.rmse: %d',RParabola,i,GoFParabola.rmse)})
                % choose the better fit and fill Output
            end
            % Fill the first 4 nm with the data from the 5th nm
            for i=1:4
                DepthDependendTipRadius(i) = DepthDependendTipRadius(5);
                DepthDependendTipShape(i) = DepthDependendTipShape(5);
            end
            
        end
        
    end
    
    methods
        % Methods for image visualisation and output
        function show_image(obj)
            % TODO: implement ui elements for customization
            
            h.Fig = figure('Name',sprintf('[Processed] %s',obj.Name),...
                'Units','pixels',...
                'Position',[200 200 1024 512],...
                'Color','k');
            
            h.B(1) = uicontrol('style','togglebutton',...
                'String','Cross Section',...
                'units','normalized',...
                'position',[.85 .5 .1 .05],...
                'Callback',@cross_section_toggle);
            
            PopUp = obj.string_of_existing();
            
            h.B(4) = uicontrol('style','text',...
                'String','Channel 1',...
                'units','normalized',...
                'position',[.85 .85 .1 .05]);
            
            h.B(2) = uicontrol('style','popupmenu',...
                'String',PopUp,...
                'units','normalized',...
                'position',[.85 .8 .1 .05],...
                'Callback',@draw_channel_1);
            
            h.B(5) = uicontrol('style','text',...
                'String','Channel 2',...
                'units','normalized',...
                'position',[.85 .7 .1 .05]);
            
            h.B(3) = uicontrol('style','popupmenu',...
                'String',PopUp,...
                'units','normalized',...
                'position',[.85 .65 .1 .05],...
                'Callback',@draw_channel_2);
            
            h.B(6) = uicontrol('style','pushbutton',...
                'String','Save Figure',...
                'units','normalized',...
                'position',[.85 .1 .1 .05],...
                'Callback',@save_figure_to_file);
            
            h.Line = [];
            h.hasCrossSection = 0;
            h.hasChannel2 = 0;
            h.B(2).Value = 2;
            draw_channel_1
            
            function cross_section_toggle(varargin)
                h.hasCrossSection = ~h.hasCrossSection;
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
                Pos1 = [h.Line.Position(1,1) h.Line.Position(1,2)];
                Pos2 = [h.Line.Position(2,1) h.Line.Position(2,2)];
                if norm(Pos1-Pos2)==0
                    get_and_draw_profile;
                    return
                end
                Profile = improfile(h.Image{1},[Pos1(1) Pos2(1)],[Pos1(2) Pos2(2)]);
                Len = norm(Pos1-Pos2)/obj.NumPixelsX*obj.ScanSizeX;
                Points = [0:1/(length(Profile)-1):1].*Len;
                [MultiplierY,UnitY,~] = AFMImage.parse_unit_scale(range(Profile),h.BaseUnit{1},1);
                [MultiplierX,UnitX,~] = AFMImage.parse_unit_scale(range(Points),'m',1);
                Ax = subplot(10,10,[71:78 81:88 91:98]);
                P = plot(Points.*MultiplierX,Profile.*MultiplierY);
                grid on
                Ax.Color = 'k';
                Ax.LineWidth = 1;
                Ax.FontSize = 18;
                Ax.XColor = 'w';
                Ax.YColor = 'w';
                Ax.GridColor = 'w';
                xlabel(sprintf('[%s]',UnitX))
                ylabel(sprintf('%s [%s]',h.Channel{1},UnitY))
                xlim([0 Points(end).*MultiplierX])
                P.LineWidth = 2;
                P.Color = 'b';
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
                Len = norm(Pos1-Pos2)/obj.NumPixelsX*obj.ScanSizeX;
                Points = [0:1/(length(Profile)-1):1].*Len;
                [MultiplierY,UnitY,~] = AFMImage.parse_unit_scale(range(Profile),h.BaseUnit{2},1);
                [MultiplierX,UnitX,~] = AFMImage.parse_unit_scale(range(Points),'m',1);
                Ax = subplot(10,10,[71:78 81:88 91:98]);
                P = plot(Points.*MultiplierX,Profile.*MultiplierY);
                grid on
                Ax.Color = 'k';
                Ax.LineWidth = 1;
                Ax.FontSize = 18;
                Ax.XColor = 'w';
                Ax.YColor = 'w';
                Ax.GridColor = 'w';
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
                Len = norm(Pos1-Pos2)/obj.NumPixelsX*obj.ScanSizeX;
                Points = [0:1/(length(Profile)-1):1].*Len;
                [MultiplierY,UnitY,~] = AFMImage.parse_unit_scale(range(Profile),h.BaseUnit{1},1);
                [MultiplierX,UnitX,~] = AFMImage.parse_unit_scale(range(Points),'m',1);
                Ax = subplot(10,10,[71:78 81:88 91:98]);
                P = plot(Points.*MultiplierX,Profile.*MultiplierY);
                grid on
                Ax.Color = 'k';
                Ax.LineWidth = 1;
                Ax.FontSize = 18;
                Ax.XColor = 'w';
                Ax.YColor = 'w';
                Ax.GridColor = 'w';
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
                h.Line.Visible = 'off';
                h.Line = drawline('Color','b','Parent',h.ImAx(2));
                addlistener(h.Line,'MovingROI',@moving_cross_section_channel_2);
                addlistener(h.Line,'ROIMoved',@moving_cross_section_channel_2);
                Pos1 = [h.Line.Position(1,1) h.Line.Position(1,2)];
                Pos2 = [h.Line.Position(2,1) h.Line.Position(2,2)];
                if norm(Pos1-Pos2)==0
                    get_and_draw_profile_channel_2;
                    return
                end
                Profile = improfile(h.Image{2},[Pos1(1) Pos2(1)],[Pos1(2) Pos2(2)]);
                Len = norm(Pos1-Pos2)/obj.NumPixelsX*obj.ScanSizeX;
                Points = [0:1/(length(Profile)-1):1].*Len;
                [MultiplierY,UnitY,~] = AFMImage.parse_unit_scale(range(Profile),h.BaseUnit{2},1);
                [MultiplierX,UnitX,~] = AFMImage.parse_unit_scale(range(Points),'m',1);
                Ax = subplot(10,10,[71:78 81:88 91:98]);
                P = plot(Points.*MultiplierX,Profile.*MultiplierY);
                grid on
                Ax.Color = 'k';
                Ax.LineWidth = 1;
                Ax.FontSize = 18;
                Ax.XColor = 'w';
                Ax.YColor = 'w';
                Ax.GridColor = 'w';
                xlabel(sprintf('[%s]',UnitX))
                ylabel(sprintf('%s [%s]',h.Channel{2},UnitY))
                xlim([0 Points(end).*MultiplierX])
                P.LineWidth = 2;
                P.Color = 'b';
            end
            
            function draw_image(LeftRight,FullPart)
                Index = 1;
                BarToImageRatio = 1/5;
                PlusBrightness = 0;
                if isequal(FullPart,'FullOne')
                    Domain = [1:8 11:18 21:28 31:38 41:48 51:58 61:68 71:78 81:88 91:98];
                elseif isequal(FullPart,'FullTwo')
                    if isequal(LeftRight,'Left')
                    Domain = [1:4 11:14 21:24 31:34 41:44 51:54 61:64 71:74 81:84 91:94];
                    else
                    Domain = [5:8 15:18 25:28 35:38 45:48 55:58 65:68 75:78 85:88 95:98];
                    Index = 2;
                    end
                elseif isequal(FullPart,'PartOne')
                    Domain = [1:8 11:18 21:28 31:38 41:48 51:58 61:68];
                elseif isequal(FullPart,'PartTwo')
                    if isequal(LeftRight,'Left')
                    Domain = [1:4 11:14 21:24 31:34 41:44 51:54 61:64];
                    Index = 1;
                    else
                    Domain = [5:8 15:18 25:28 35:38 45:48 55:58 65:68];
                    Index = 2;
                    end
                end
                
                h.Channel{Index} = h.B(1+Index).String{h.B(1+Index).Value};
                switch h.Channel{Index}
                    case 'Height(Trace)'
                        h.BaseUnit{Index} = 'm';
                        h.Image{Index} = obj.Height.Trace;
                        ColorPattern = obj.CMap;
                    case 'Height(Retrace)'
                        h.BaseUnit{Index} = 'm';
                        h.Image{Index} = obj.Height.ReTrace;
                        ColorPattern = obj.CMap;
                    case 'Height-Measured(Trace)'
                        h.BaseUnit{Index} = 'm';
                        h.Image{Index} = obj.HeightMeasured.Trace;
                        ColorPattern = obj.CMap;
                    case 'Height-Measured(Retrace)'
                        h.BaseUnit{Index} = 'm';
                        h.Image{Index} = obj.HeightMeasured.ReTrace;
                        ColorPattern = obj.CMap;
                    case 'ErrorSignal(Trace)'
                        h.BaseUnit{Index} = 'V';
                        h.Image{Index} = obj.ErrorSignal.Trace;
                        ColorPattern = obj.CMap;
                    case 'ErrorSignal(Retrace)'
                        h.BaseUnit{Index} = 'V';
                        h.Image{Index} = obj.ErrorSignal.ReTrace;
                        ColorPattern = obj.CMap;
                    case 'LateralDeflection(Trace)'
                        h.BaseUnit{Index} = 'm';
                        h.Image{Index} = obj.LateralDeflection.Trace;
                        ColorPattern = obj.CMap;
                    case 'LateralDeflection(Retrace)'
                        h.BaseUnit{Index} = 'm';
                        h.Image{Index} = obj.LateralDeflection.ReTrace;
                        ColorPattern = obj.CMap;
                    case 'LockInAmplitude(Trace)'
                        h.BaseUnit{Index} = 'm';
                        h.Image{Index} = obj.LockInAmplitude.Trace;
                        ColorPattern = obj.CMap;
                    case 'LockInAmplitude(Retrace)'
                        h.BaseUnit{Index} = 'm';
                        h.Image{Index} = obj.LockInAmplitude.ReTrace;
                        ColorPattern = obj.CMap;
                    case 'LockInPhase(Trace)'
                        h.BaseUnit{Index} = 'deg';
                        h.Image{Index} = obj.LockInPhase.Trace;
                        ColorPattern = obj.CMap;
                    case 'LockInPhase(Retrace)'
                        h.BaseUnit{Index} = 'deg';
                        h.Image{Index} = obj.LockInPhase.ReTrace;
                        ColorPattern = obj.CMap;
                    case 'VerticalDeflection(Trace)'
                        h.BaseUnit{Index} = 'm';
                        h.Image{Index} = obj.Height.Trace;
                        ColorPattern = obj.CMap;
                    case 'VerticalDeflection(Retrace)'
                        h.BaseUnit{Index} = 'm';
                        h.Image{Index} = obj.Height.Trace;
                        ColorPattern = obj.CMap;
                    case 'Processed'
                        h.BaseUnit{Index} = 'm';
                        h.Image{Index} = obj.Processed;
                        ColorPattern = obj.CMap;
                    case 'none'
                        if Index == 1
                            h.hasChannel1 = 0;
                        elseif Index == 2
                            h.hasChannel2 = 0;
                        end
                        return
                end
                
                [Multiplier,Unit,~] = AFMImage.parse_unit_scale(range(h.Image{Index},'all'),h.BaseUnit{Index},1);
                h.ImAx(Index) = subplot(10,10,Domain);
                set(gca,'LooseInset',get(gca,'TightInset'));
%                 if Index == 1
%                     Shift = 0;
%                 else
%                     Shift = 1;
%                 end
%                 h.ImAx(Index).InnerPosition(1) = 0.06+Shift*h.ImAx(Index).InnerPosition(3);
                h.I(Index) = imshow(h.Image{Index}*Multiplier,[],'Colormap',ColorPattern);
                if Index == 1
                    h.I(Index).ButtonDownFcn = @get_and_draw_profile_channel_1;
                else
                    h.I(Index).ButtonDownFcn = @get_and_draw_profile_channel_2;
                end
                hold on
                AFMImage.draw_scalebar_into_current_image(obj.NumPixelsX,obj.ScanSizeX,BarToImageRatio,h.ImAx(Index).Position(3));
                c = colorbar;
                c.FontSize = round(18*(obj.NumPixelsX/1024));
                c.Color = 'w';
                c.Label.String = sprintf('%s [%s]',h.Channel{Index},Unit);
                c.Label.FontSize = round(22*(obj.NumPixelsX/1024));
                c.Label.Color = 'w';
            end
            
            function save_figure_to_file(varargin)
                Pos = h.Fig.InnerPosition;
                Frame = getframe(h.Fig);
                
                filter = {'*.png';'*.tif'};
                [file, path] = uiputfile(filter);
                
                FullFile = fullfile(path,file);
                save(FullFile,'Frame.cdata');
            end
            
            uiwait(h.Fig)
        end
        
    end
    
    methods
        % Auxiliary methods
        
        function read_in_header_properties(obj,ImageFullFile)
            % determines imaging type (contact,AC) and reads out general image
            % properties
            
            FileInfo = imfinfo(ImageFullFile);
            
            obj.NumChannels = numel(FileInfo) - 1;
            
            obj.ImagingType = upper(FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32816)).Value);
            
            try
                obj.FileVersion = FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32768)).Value;
            catch
                warning("Couldn't determine FileVersion")
                obj.FileVersion = nan;
            end
            try
                obj.DateTime = FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32771)).Value;
            catch
                warning("Couldn't determine Date and Time of file creation")
                obj.DateTime = nan;
            end
            try
                obj.OriginX = FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32832)).Value;
            catch
                warning("Couldn't determine OriginX")
                obj.OriginX = nan;
            end
            try
                obj.OriginY = FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32833)).Value;
            catch
                warning("Couldn't determine OriginY")
                obj.OriginY = nan;
            end
            try
                obj.ScanSizeX = FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32834)).Value;
            catch
                warning("Couldn't determine ScanSizeX")
                obj.ScansizeX = nan;
            end
            try
                obj.ScanSizeY = FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32835)).Value;
            catch
                warning("Couldn't determine ScanSizeY")
                obj.ScanSizeY = nan;
            end
            try
                obj.NumPixelsX = FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32838)).Value;
            catch
                warning("Couldn't determine NumPixelsX")
                obj.NumPixelsX = nan;
            end
            try
                obj.NumPixelsY = FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32839)).Value;
            catch
                warning("Couldn't determine NumPixelsY")
                obj.NumPixelsY = nan;
            end
            
            % Now read out ImagingType-specific properties
            
            if isequal(obj.ImagingType,'AC') || isequal(obj.ImagingType,'AC-DIRECT-DRIVE')
                obj.read_in_ac_header(FileInfo)
            elseif isequal(obj.ImagingType,'CONTACT')
                obj.read_in_contact_header(FileInfo)
            else
                error('Unknown imaging mode')
            end
            
            
            
        end
        
        function read_in_ac_header(obj,FileInfo)
            
            try
                obj.IGain = abs(FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32818)).Value);
            catch
                obj.IGain = nan;
                warning ("Couldn't determine integral Gain")
            end
            try
                obj.PGain = abs(FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32817)).Value);
            catch
                obj.PGain = nan;
                warning ("Couldn't determine proportional Gain")
            end
            try
                obj.RelativeSetpoint = FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32819)).Value;
            catch
                obj.RelativeSetpoint = nan;
                warning ("Couldn't determine relative setpoint")
            end
            try
                obj.DriveAmplitude = FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32822)).Value;
            catch
                obj.DriveAmplitude = nan;
                warning ("Couldn't determine drive amplitude")
            end
            try
                obj.DriveFrequency = FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32823)).Value;
            catch
                obj.DriveFrequency = nan;
                warning ("Couldn't determine drive frequency")
            end
            try
                obj.PhaseShift = rad2deg(FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32824)).Value);
            catch
                obj.PhaseShift = nan;
                warning ("Couldn't determine phase shift")
            end
            try
                obj.LineRate = FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32841)).Value;
                obj.TipVelocity = obj.LineRate*obj.ScanSizeX;
            catch
                obj.LineRate = nan;
                obj.TipVelocity = nan;
                warning ("Couldn't determine line rate and tip velocity")
            end
            try
                obj.Sensitivity = FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==33028)).Value/...
                    FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32980)).Value;
            catch
                obj.Sensitivity = nan;
                warning ("Couldn't determine sensitivity. Image isn't calibrated")
            end
            try
                obj.SpringConstant = FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==33076)).Value/...
                    FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==33028)).Value;
            catch
                obj.SpringConstant = nan;
                warning ("Couldn't determine spring constant. Image isn't calibrated")
            end
            
        end
        
        function read_in_contact_header(obj,FileInfo)
            
            flag_data=strsplit(FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32830)).Value)';
            
            flag=find(~cellfun(@isempty,strfind(flag_data,'setpoint-feedback-settings.i-gain')));
            if(~isempty(flag))
                obj.IGain=cellfun(@str2double, flag_data(flag+2,1));
            else
                obj.IGain=nan;
            end
            flag=find(~cellfun(@isempty,strfind(flag_data,'setpoint-feedback-settings.p-gain')));
            if(~isempty(flag))
                obj.PGain=cellfun(@str2double, flag_data(flag+2,1));
            else
                obj.PGain=nan;
            end
            
            if(~isempty(find([FileInfo(1).UnknownTags.ID]==33028)))&&(~isempty(find([FileInfo(1).UnknownTags.ID]==32980)))
                obj.Sensitivity=(FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==33028)).Value)/(FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32980)).Value);
            else
                warning('No Sensitivity calibration slot found, image is uncalibarated')
                Uncalibrated=1;
                obj.Sensitivity=1;
            end
            
            if(~isempty(find([FileInfo(1).UnknownTags.ID]==33076)))&&(~isempty(find([FileInfo(1).UnknownTags.ID]==32980)))
                obj.SpringConstant=(FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==33076)).Value)/(FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==33028)).Value);
            else
                warning('No Spring Constant calibration slot found, image is uncalibarated')
                Uncalibrated=1;
                obj.SpringConstant=1;
            end
            
            if(~isempty(find([FileInfo(1).UnknownTags.ID]==32836)))
                obj.ScanAngle=rad2deg(FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32836)).Value);
            else
                warning('No Scan Angle slot found')
                obj.ScanAngle=nan;
            end
            
            if(FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32820)).Value==1)
                obj.Baseline_Raw=((FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32819)).Value-FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32821)).Value)-FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32980)).Value);
                obj.Bline_adjust='Yes';
            else
                obj.Baseline_Raw=((FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32819)).Value)-FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32980)).Value);
                obj.Bline_adjust='No';
            end
            
            obj.SetP_V=obj.Baseline_Raw;
            
            if(~isempty(find([FileInfo(1).UnknownTags.ID]==32981)))&&(~isempty(find([FileInfo(1).UnknownTags.ID]==32980)))
                obj.Raw=(obj.Baseline_Raw-(FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32981)).Value))/(FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32980)).Value); % Setpoint in Volts [V]
            else
                warning('No Baseline correction found')
                obj.Raw=obj.SetP_V;
                
            end
            
            if(~isempty(find([FileInfo(1).UnknownTags.ID]==33028)))&&(~isempty(find([FileInfo(1).UnknownTags.ID]==33029)))
                obj.SetP_m=(obj.Raw)*(FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==33028)).Value)+(FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==33029)).Value); % Setpoint in meters [m]
            else
                warning('No Setpoint (in meters) calibration slot found, image might be uncalibarated')
                obj.SetP_m=(obj.Raw);
            end
            
            if(~isempty(find([FileInfo(1).UnknownTags.ID]==33076)))&&(~isempty(find([FileInfo(1).UnknownTags.ID]==33077)))
                obj.SetP_N=(obj.Raw)*(FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==33076)).Value)+(FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==33077)).Value); % Setpoint in force [N]
            else
                warning('No Setpoint (in Newton) calibration slot found, image might be uncalibarated')
                obj.SetP_N=(obj.Raw);
            end
            
        end
        
        function load_image_channels(obj,ImageFullFile)
            % load in existing image channels and related properties
            % Most parts of the following code have been copied from:
            % Ortuso, Roberto D., Kaori Sugihara.
            % "Detailed Study on the Failure of the Wedge Calibration Method at Nanonewton Setpoints for Friction Force Microscopy."
            % The Journal of Physical Chemistry C 122.21 (2018): 11464-11474.
            
            FileInfo = imfinfo(ImageFullFile);
            
            for i=2:obj.NumChannels+1
                Channel_Name=(FileInfo(i).UnknownTags(find([FileInfo(i).UnknownTags.ID]==32850)).Value);
                
                strsp=(strsplit((FileInfo(i).UnknownTags(find([FileInfo(i).UnknownTags.ID]==32851)).Value)))';
                
                for k=1:size(strsp,1)
                    if(strcmp(strsp{k,1},'retrace')==1)
                        if(strcmp(strsp{k+2,1},'true'))
                            trace_type_flag='ReTrace';
                        else
                            trace_type_flag='Trace';
                        end
                        break
                    end
                end
                
                typpe_of_ch=FileInfo(i).UnknownTags(find([FileInfo(i).UnknownTags.ID]==32897)).Value;
                
                if(strcmp(typpe_of_ch,'nominal')||(strcmp(typpe_of_ch,'voltsamplitude')))
                    m_ID=33028;
                    off_ID=33029;
                elseif ((strcmp(typpe_of_ch,'force'))||(strcmp(typpe_of_ch,'calibrated'))||(strcmp(typpe_of_ch,'distanceamplitude')))
                    m_ID=33076;
                    off_ID=33077;
                elseif(strcmp(typpe_of_ch,'volts'))
                    typpe_of_ch_det=FileInfo(i).UnknownTags(find([FileInfo(i).UnknownTags.ID]==32848)).Value;
                    if(strcmp(typpe_of_ch_det,'capacitiveSensorXPosition'))||(strcmp(typpe_of_ch_det,'servoDacY'))||(strcmp(typpe_of_ch_det,'servoDacX'))||(strcmp(typpe_of_ch_det,'capacitiveSensorYPosition'))
                        m_ID=33028;
                        off_ID=33029;
                    else
                        m_ID=32980;
                        off_ID=32981;
                    end
                else
                    m_ID=32980;
                    off_ID=32981;
                end
                
                if(~isempty(find([FileInfo(i).UnknownTags.ID]==m_ID)))
                    multiplyer=FileInfo(i).UnknownTags(find([FileInfo(i).UnknownTags.ID]==m_ID)).Value;
                else
                    warning('No multiplyer slot found, image might be uncalibrated uncalibarated')
                    multiplyer=1;
                end
                
                if(~isempty(find([FileInfo(i).UnknownTags.ID]==off_ID)))
                    offset=FileInfo(i).UnknownTags(find([FileInfo(i).UnknownTags.ID]==off_ID)).Value;
                else
                    warning('No offset slot found, image might be uncalibrated uncalibarated')
                    offset=0;
                end
                
                
                
                if(~strcmp(Channel_Name,'Vertical Deflection'))
                    afm_image=((double(imread(ImageFullFile,i))*multiplyer))+offset;
                else
                    if(strcmp(obj.Bline_adjust,'No'))
                        afm_image=((double(imread(ImageFullFile,i))*multiplyer))+offset;
                    else
                        obj.Baseline_N=(obj.Baseline_Raw*multiplyer)+offset;
                        afm_image=((double(imread(ImageFullFile,i))*multiplyer))+offset;
                    end
                end
                
                afm_image = afm_image(end:-1:1,:); % mirror Y-pixels to flip image to same orientation as in jpk data processing
                
                if isequal(Channel_Name,'Height') && isequal(trace_type_flag,'Trace')
                    obj.Height.Trace = afm_image;
                    obj.Height.Multiplier = multiplyer;
                    obj.Height.Offset = offset;
                    obj.HasHeight = true;
                elseif isequal(Channel_Name,'Height (measured)') && isequal(trace_type_flag,'Trace')
                    obj.HeightMeasured.Trace = afm_image;
                    obj.HeightMeasured.Multiplier = multiplyer;
                    obj.HeightMeasured.Offset = offset;
                    obj.HasHeightMeasured = true;
                elseif isequal(Channel_Name,'Lock-In Amplitude') && isequal(trace_type_flag,'Trace')
                    obj.LockInAmplitude.Trace = afm_image;
                    obj.LockInAmplitude.Multiplier = multiplyer;
                    obj.LockInAmplitude.Offset = offset;
                    obj.HasLockInAmplitude = true;
                elseif isequal(Channel_Name,'Lock-In Phase') && isequal(trace_type_flag,'Trace')
                    obj.LockInPhase.Trace = afm_image;
                    obj.LockInPhase.Multiplier = multiplyer;
                    obj.LockInPhase.Offset = offset;
                    obj.HasLockInPhase = true;
                elseif isequal(Channel_Name,'Error Signal') && isequal(trace_type_flag,'Trace')
                    obj.ErrorSignal.Trace = afm_image;
                    obj.ErrorSignal.Multiplier = multiplyer;
                    obj.ErrorSignal.Offset = offset;
                    obj.HasErrorSignal = true;
                elseif isequal(Channel_Name,'Lateral Deflection') && isequal(trace_type_flag,'Trace')
                    obj.LateralDeflection.Trace = afm_image;
                    obj.LateralDeflection.Multiplier = multiplyer;
                    obj.LateralDeflection.Offset = offset;
                    obj.HasLateralDeflection = true;
                elseif isequal(Channel_Name,'Vertical Deflection') && isequal(trace_type_flag,'Trace')
                    obj.VerticalDeflection.Trace = afm_image;
                    obj.VerticalDeflection.Multiplier = multiplyer;
                    obj.VerticalDeflection.Offset = offset;
                    obj.HasVerticalDeflection = true;
                elseif isequal(Channel_Name,'Height') && isequal(trace_type_flag,'ReTrace')
                    obj.Height.ReTrace = afm_image;
                    obj.Height.Multiplier = multiplyer;
                    obj.Height.Offset = offset;
                    obj.HasHeight = true;
                elseif isequal(Channel_Name,'Height (measured)') && isequal(trace_type_flag,'ReTrace')
                    obj.HeightMeasured.ReTrace = afm_image;
                    obj.HeightMeasured.Multiplier = multiplyer;
                    obj.HeightMeasured.Offset = offset;
                    obj.HasHeightMeasured = true;
                elseif isequal(Channel_Name,'Lock-In Amplitude') && isequal(trace_type_flag,'ReTrace')
                    obj.LockInAmplitude.ReTrace = afm_image;
                    obj.LockInAmplitude.Multiplier = multiplyer;
                    obj.LockInAmplitude.Offset = offset;
                    obj.HasLockInAmplitude = true;
                elseif isequal(Channel_Name,'Lock-In Phase') && isequal(trace_type_flag,'ReTrace')
                    obj.LockInPhase.ReTrace = afm_image;
                    obj.LockInPhase.Multiplier = multiplyer;
                    obj.LockInPhase.Offset = offset;
                    obj.HasLockInPhase = true;
                elseif isequal(Channel_Name,'Error Signal') && isequal(trace_type_flag,'ReTrace')
                    obj.ErrorSignal.ReTrace = afm_image;
                    obj.ErrorSignal.Multiplier = multiplyer;
                    obj.ErrorSignal.Offset = offset;
                    obj.HasErrorSignal = true;
                elseif isequal(Channel_Name,'Lateral Deflection') && isequal(trace_type_flag,'ReTrace')
                    obj.LateralDeflection.ReTrace = afm_image;
                    obj.LateralDeflection.Multiplier = multiplyer;
                    obj.LateralDeflection.Offset = offset;
                    obj.HasLateralDeflection = true;
                elseif isequal(Channel_Name,'Vertical Deflection') && isequal(trace_type_flag,'ReTrace')
                    obj.VerticalDeflection.ReTrace = afm_image;
                    obj.VerticalDeflection.Multiplier = multiplyer;
                    obj.VerticalDeflection.Offset = offset;
                    obj.HasVerticalDeflection = true;
                end
                
            end
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
                obj.HostOS = OS;
                obj.HostName = Host;
            end
        end
        
        function initialize_flags(obj)
            
            obj.HasHeight = false;
            obj.HasHeightMeasured = false;
            obj.HasLateralDeflection = false;
            obj.HasLockInAmplitude = false;
            obj.HasLockInPhase = false;
            obj.HasVerticalDeflection = false;
            obj.HasProcessed = false;
            
            obj.hasDeconvolutedCantileverTip = false;
        end
        
        function PopUp = string_of_existing(obj)
            PopUp{1} = 'none';
            k = 2;
            if obj.HasProcessed
                PopUp{k} = 'Processed';
                k = k + 1;
            end
            if obj.HasHeight
                PopUp{k} = 'Height(Trace)';
                PopUp{k+1} = 'Height(Retrace)';
                k = k + 2;
            end
            if obj.HasHeightMeasured
                PopUp{k} = 'Height-Measured(Trace)';
                PopUp{k+1} = 'Height-Measured(Retrace)';
                k = k + 2;
            end
            if obj.HasErrorSignal
                PopUp{k} = 'ErrorSignal(Trace)';
                PopUp{k+1} = 'ErrorSignal(Retrace)';
                k = k + 2;
            end
            if obj.HasLateralDeflection
                PopUp{k} = 'LateralDeflection(Trace)';
                PopUp{k+1} = 'LateralDeflection(Retrace)';
                k = k + 2;
            end
            if obj.HasLockInAmplitude
                PopUp{k} = 'LockInAmplitude(Trace)';
                PopUp{k+1} = 'LockInAmplitude(Retrace)';
                k = k + 2;
            end
            if obj.HasLockInPhase
                PopUp{k} = 'LockInPhase(Trace)';
                PopUp{k+1} = 'LockInPhase(Retrace)';
                k = k + 2;
            end
            if obj.HasVerticalDeflection
                PopUp{k} = 'VerticalDeflection(Trace)';
                PopUp{k+1} = 'VerticalDeflection(Retrace)';
                k = k + 2;
            end
        end
        
    end
    
    methods(Static)
        % Static auxiliary methods.
        function Name = parse_file_name(FullFile)
            % Parses filename to fill property Name. Static method to
            % possibly be used in other constructors aswell
            FileSepPos = strfind(FullFile,filesep);
            FileExtPos = strfind(FullFile,'.');
            File = FullFile(FileSepPos(end)+1:FileExtPos(end)-1);
            File = replace(File,'_','-');
            
            Name = File;
        end
        
        function [Cone,SizeOfPixelX,SizeOfPixelY] = cone(ConePixelX,ConePixelY,ConeHeigth,ScanSizeX,ScanSizeY,ConeTipRadius)
            %cone.m, version 1.1
            %The file creates a TGT1 grating surface with a cone peaking at the centre
            %of the sample. The angle of the is assumed to be 50 degrees and the height
            %is taken to be equal to the scan height. The radius of the cone tip can be
            %changed so that different levels of accuracy may be achieved.
            %Ask the user for the size of the scan.
            %Code adapted from original code by O. Andriotis(2014)
            
            %Variables
            height_cone=ConeHeigth;
            angle_cone=50;
            height_loss=(ConeTipRadius*cosd(angle_cone/2))/(tand(angle_cone/2))+(ConeTipRadius*sind(angle_cone/2))-ConeTipRadius;
            %Calculates the height that is lost from a perfect tip when a rounded tip
            %is used instead.
            height_tip=height_cone+height_loss;%The ideal tip is derived from the required
            %height of the cone (from the scan height) and the amount of height loss
            %experienced for the desired grating radius. This value is then used to
            %generate the ideal cone with a perfect tip such that when the curved tip
            %is added the height of the cone is equal to the scan height.
            
            radius_cone=tand(angle_cone/2)*height_tip; %Calculates the cone radius.
            Cone=zeros(ConePixelX,ConePixelY); %initiates a flat surface of size equal to scan.
            SizeOfPixelX=ScanSizeX/ConePixelX;
            SizeOfPixelY=ScanSizeY/ConePixelY;
            Cone(floor(ConePixelX/2), floor(ConePixelY/2))=height_tip; % Positions the cone.
            max_pixel_movement_x=floor(radius_cone/SizeOfPixelX); %Radius of cone divided
            %by size of a pixel to find the maximal number of pixels in line in the
            %cone radius.
            max_pixel_movement_y=floor(radius_cone/SizeOfPixelY);
            
            %Determine the limits of the cone whether it fits completely in the image
            %or not. Done for each dimension and limit using the centre of the image as
            %a reference.
            if ConePixelX/2-max_pixel_movement_x>=1
                limit_x_1=ConePixelX/2-max_pixel_movement_x;
            else limit_x_1=1;
            end
            
            if ConePixelX/2+max_pixel_movement_x<=ConePixelX
                limit_x_2=ConePixelX/2+max_pixel_movement_x;
            else limit_x_2=ConePixelX;
            end
            
            if ConePixelY/2-max_pixel_movement_y>=1
                limit_y_1=ConePixelY/2-max_pixel_movement_y;
            else limit_y_1=1;
            end
            
            if ConePixelY/2+max_pixel_movement_y<=ConePixelY
                limit_y_2=ConePixelY/2+max_pixel_movement_y;
            else limit_y_2=ConePixelY;
            end
            
            
            curve_start_height = height_tip - ((ConeTipRadius*cosd(angle_cone/2))/(tand(angle_cone/2)));
            %Calculates the hieght at which the cone leaves from a constant gradient
            %into the curved profile.
            
            %Generates the cone
            for i=limit_x_1:limit_x_2
                for j=limit_y_1:limit_y_2
                    distance=sqrt(((ConePixelX/2-i)*SizeOfPixelX)^2+((ConePixelY/2-j)*SizeOfPixelY)^2);
                    %Distance of point i,j with reference from the centre of the image.
                    
                    curve_height=sqrt(((ConeTipRadius)^2)-(distance^2))-(ConeTipRadius*sind(angle_cone/2));
                    %The absolute hieght of each point of the curved tip.
                    
                    if distance<=radius_cone;
                        Cone(i,j)=(radius_cone-distance)*height_tip/radius_cone;
                        %If the dstance is smaller than than the radius then the
                        %constant slope of the cone is generated.
                    end
                    
                    if distance<=ConeTipRadius*cosd(angle_cone/2);
                        Cone(i,j) = curve_height + curve_start_height;
                        %If the distance is smaller than the radius of the curved peak
                        %radius then the curved peak is assumed.
                    end
                    
                end
            end
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
        
        function [X,Y,Z] = convert_masked_to_point_cloud(Image,Mask,ScanSizeX,ScanSizeY,OriginX,OriginY)
            
            if nargin<4
                ScanSizeX = 1;
                ScanSizeY = 1;
                OriginX = 0;
                OriginY = 0;
            end
            if nargin<6
                OriginX = 0;
                OriginY = 0;
            end
            
            
            [X,Y] = find(Mask);
            Z = zeros(length(X),1);
            for i=1:length(X)
                Z(i) = Image(X(i),Y(i));
            end
            
            X = X.*ScanSizeX - OriginX;
            Y = Y.*ScanSizeY - OriginY;
        end
        
        function R = draw_scalebar_into_current_image(NumPixelsX,ScanSizeX,BarToImageRatio,FontSizeMult)
            
            if nargin < 5
                BarToImageRatio = 1/5;
            end
            
            ScalebarThickness = 1/40;
            DistFromBorder = 0.1;
            
            [Multiplier,Unit,SnapTo] = AFMImage.parse_unit_scale(ScanSizeX,'m',BarToImageRatio);
            
            Width = (SnapTo)/(ScanSizeX*Multiplier);
            Height = ScalebarThickness;
            Left = 1-Width*(1+DistFromBorder);
            Bottom = 1-(DistFromBorder-Height);
            
            R = rectangle('Position',[Left Bottom Width Height].*NumPixelsX);
            R.FaceColor = 'w';
            R.EdgeColor = 'k';
            R.LineWidth = 2;
            
            
            A = text((Left+Width/2)*NumPixelsX,(Bottom-Height)*NumPixelsX,sprintf('%i %s',SnapTo,Unit),'HorizontalAlignment','center');
            A.Color = 'w';
            A.FontSize = round(32*(NumPixelsX/1024*FontSizeMult));
            A.FontWeight = 'bold';
            
            
        end
        
        function [Multiplier,Unit,SnapTo] = parse_unit_scale(ScanSize,Unit,Mult)
            
            ScaleSize = ScanSize*Mult;
            
            if (ScaleSize < 1e-11)
                Multiplier = 1e12;
                Prefix = 'p';
                SnapTo = round(ScaleSize*Multiplier,0);
            elseif (ScaleSize < 1e-10) && (ScaleSize >= 1e-11)
                Multiplier = 1e12;
                Prefix = 'p';
                SnapTo = round(ScaleSize*Multiplier,-1);
            elseif (ScaleSize < 1e-9) && (ScaleSize >= 1e-10)
                Multiplier = 1e12;
                Prefix = 'p';
                SnapTo = round(ScaleSize*Multiplier,-2);
            elseif (ScaleSize < 1e-8) && (ScaleSize >= 1e-9)
                Multiplier = 1e9;
                Prefix = 'n';
                SnapTo = round(ScaleSize*Multiplier,0);
            elseif (ScaleSize < 1e-7) && (ScaleSize >= 1e-8)
                Multiplier = 1e9;
                Prefix = 'n';
                SnapTo = round(ScaleSize*Multiplier,-1);
            elseif (ScaleSize < 1e-6) && (ScaleSize >= 1e-7)
                Multiplier = 1e9;
                Prefix = 'n';
                SnapTo = round(ScaleSize*Multiplier,-2);
            elseif (ScaleSize < 1e-5) && (ScaleSize >= 1e-6)
                Multiplier = 1e6;
                Prefix = '\mu';
                SnapTo = round(ScaleSize*Multiplier,0);
            elseif (ScaleSize < 1e-4) && (ScaleSize >= 1e-5)
                Multiplier = 1e6;
                Prefix = '\mu';
                SnapTo = round(ScaleSize*Multiplier,-1);
            elseif (ScaleSize < 1e-3) && (ScaleSize >= 1e-4)
                Multiplier = 1e6;
                Prefix = '\mu';
                SnapTo = round(ScaleSize*Multiplier,-2);
            elseif (ScaleSize < 1e-4) && (ScaleSize >= 1e-3)
                Multiplier = 1e3;
                Prefix = 'm';
                SnapTo = round(ScaleSize*Multiplier,0);
            elseif (ScaleSize < 1e-3) && (ScaleSize >= 1e-2)
                Multiplier = 1e3;
                Prefix = 'm';
                SnapTo = round(ScaleSize*Multiplier,-1);
            elseif (ScaleSize < 1e-2) && (ScaleSize >= 1e-1)
                Multiplier = 1e3;
                Prefix = 'm';
                SnapTo = round(ScaleSize*Multiplier,-2);
            elseif (ScaleSize < 1e1) && (ScaleSize >= 1e0)
                Multiplier = 1e0;
                Prefix = '';
                SnapTo = round(ScaleSize*Multiplier,0);
            elseif (ScaleSize < 1e2) && (ScaleSize >= 1e1)
                Multiplier = 1e0;
                Prefix = '';
                SnapTo = round(ScaleSize*Multiplier,-1);
            elseif (ScaleSize < 1e3) && (ScaleSize >= 1e2)
                Multiplier = 1e0;
                Prefix = '';
                SnapTo = round(ScaleSize*Multiplier,-2);
            elseif (ScaleSize < 1e4) && (ScaleSize >= 1e3)
                Multiplier = 1e-3;
                Prefix = 'k';
                SnapTo = round(ScaleSize*Multiplier,0);
            elseif (ScaleSize < 1e5) && (ScaleSize >= 1e4)
                Multiplier = 1e-3;
                Prefix = 'k';
                SnapTo = round(ScaleSize*Multiplier,-1);
            elseif (ScaleSize < 1e6) && (ScaleSize >= 1e5)
                Multiplier = 1e-3;
                Prefix = 'k';
                SnapTo = round(ScaleSize*Multiplier,-2);
            elseif (ScaleSize < 1e7) && (ScaleSize >= 1e6)
                Multiplier = 1e-6;
                Prefix = 'M';
                SnapTo = round(ScaleSize*Multiplier,0);
            elseif (ScaleSize < 1e8) && (ScaleSize >= 1e7)
                Multiplier = 1e-6;
                Prefix = 'M';
                SnapTo = round(ScaleSize*Multiplier,-1);
            elseif (ScaleSize < 1e9) && (ScaleSize >= 1e8)
                Multiplier = 1e-6;
                Prefix = 'M';
                SnapTo = round(ScaleSize*Multiplier,-2);
            elseif (ScaleSize < 1e10) && (ScaleSize >= 1e9)
                Multiplier = 1e-9;
                Prefix = 'G';
                SnapTo = round(ScaleSize*Multiplier,0);
            elseif (ScaleSize < 1e11) && (ScaleSize >= 1e10)
                Multiplier = 1e-9;
                Prefix = 'G';
                SnapTo = round(ScaleSize*Multiplier,-1);
            elseif (ScaleSize < 1e12) && (ScaleSize >= 1e11)
                Multiplier = 1e-9;
                Prefix = 'G';
                SnapTo = round(ScaleSize*Multiplier,-2);
            elseif (ScaleSize < 1e13) && (ScaleSize >= 1e12)
                Multiplier = 1e-12;
                Prefix = 'T';
                SnapTo = round(ScaleSize*Multiplier,0);
            elseif (ScaleSize < 1e14) && (ScaleSize >= 1e13)
                Multiplier = 1e-12;
                Prefix = 'T';
                SnapTo = round(ScaleSize*Multiplier,-1);
            elseif (ScaleSize < 1e15) && (ScaleSize >= 1e14)
                Multiplier = 1e-12;
                Prefix = 'T';
                SnapTo = round(ScaleSize*Multiplier,-2);
            else
                Multiplier = 1;
                Prefix = '';
                SnapTo = round(ScaleSize*Multiplier);
            end
            
            Unit = sprintf('%s%s',Prefix,Unit);
            
        end
        
    end
end
