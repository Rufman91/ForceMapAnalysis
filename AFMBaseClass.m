classdef AFMBaseClass < matlab.mixin.Copyable & matlab.mixin.SetGet & handle
    % This a baseclass for the classes AFMImage and ForceMap to inherit
    % shared methods and properties from
    
    properties
        Name
        Folder
        ID
        HostOS          % Operating System
        HostName        % Name of hosting system
        FileType = 'Image'
        ScanSizeX           % Size of imaged window in X-direction
        ScanSizeY           % Size of imaged window in Y-direction
        ScanAngle = 0   % in degrees (°)
        NumPixelsX
        NumPixelsY
        OriginX = 0
        OriginY = 0
        List2Map        % An R->RxR ((k)->(i,j)) mapping of indices to switch between the two representations
        Map2List        % An RxR->R ((i,j)->(k))mapping of indices to switch between the two representations
    end
    properties
        % All possible image channels. The Channels are all part of the
        % struct Channel and should each contain the properties Image,
        % Unit, Name, ScanSizeX, ScanSizeY, NumPixelsX, NumPixelsY,
        % ScanAngle, OriginX and OriginY. This might seem redundant but
        % allows for image cropping, image-overlays and easy addition of
        % other kinds of image-data (e.g. AdhesionMaps, EModMaps)
        Channel
        CMap = AFMImage.define_afm_color_map
    end
    properties
        % Properties for handling of image segmentation, fibril
        % segmentation, specifically
        Segment = struct('Name',[],...
                            'Type',[],...
                            'SubSegmentName',[],...
                            'ROIObject',[],...
                            'ProximityMap',[])
        OverlayGroup
    end
    
    methods
        % Main Methods
        
        function obj = AFMBaseClass()
            
            obj.OverlayGroup.hasOverlayGroup = false;
        end
        
        function save_afm_class(obj,DataFolder)
            
            current = what;
            
            if nargin > 1
                obj.Folder = [DataFolder regexprep(obj.Name,'[.]','') obj.ID];
            end
            
            warning('off')
            mkdir(obj.Folder);
            warning('on')
            cd(obj.Folder);
            PropertyStruct = obj.get();
            PropertyNames = fieldnames(PropertyStruct);
            
            for i=1:length(PropertyNames)
                TempProp = getfield(PropertyStruct,PropertyNames{i});
                save([PropertyNames{i} '.mat'],'TempProp','-v7');
            end
            
            cd(current.path)
        end
        
        function load_afm_class_properties(obj,Folder)
            
            
            current = what;
            
            cd(Folder);
            
            MetaClass = metaclass(obj);
            MetaProperties = MetaClass.PropertyList;
            
            PropertyStruct = obj.get();
            PropertyNames = fieldnames(PropertyStruct);
            
            for i=1:length(PropertyNames)
                try
                    LoadStruct = load([PropertyNames{i} '.mat']);
                    set(obj,PropertyNames{i},LoadStruct.TempProp);
                catch
                    if MetaProperties(i).HasDefault
                        set(obj,PropertyNames{i},MetaProperties.DefaultValue);
                    else
                        set(obj,PropertyNames{i},[]);
                    end
                end
            end
            
            cd(current.path)
            
        end
        
        function deconvolute_image(obj,CTClassInstance)
            % Comment
            
            ErodedTip = CTClassInstance.get_channel('Eroded Tip');
            if isempty(ErodedTip)
                CTClassInstance.deconvolute_cantilever_tip;
                ErodedTip = CTClassInstance.get_channel('Eroded Tip');
            end
            Processed = obj.get_channel('Processed');
            if isempty('Processed')
                warning('Image needs to be flattened first')
                return
            end
            
            % Resize and pad images to same resolution at correct spacial
            % ratio
            TempMultiplier = (ErodedTip.NumPixelsX/ErodedTip.ScanSizeX)/(Processed.NumPixelsX/Processed.ScanSizeX);
            Multiplier = TempMultiplier;
            if TempMultiplier*Processed.NumPixelsX > 1024
                ReductionFactor = 1024/(TempMultiplier*Processed.NumPixelsX);
                Multiplier = TempMultiplier*ReductionFactor;
                ReducedRes = round(ErodedTip.NumPixelsX*ReductionFactor);
                ErodedTip.Image = imresize(ErodedTip.Image,[ReducedRes nan]);
            end
            if TempMultiplier >= 1
                Smaller = ErodedTip.Image;
                Bigger = Processed.Image;
            else
                Smaller = Processed.Image;
                Bigger = ErodedTip.Image;
            end
            Bigger = imresize(Bigger,[round(Multiplier*size(Bigger,1)) nan]);
            Smaller = padarray(Smaller,size(Bigger)-size(Smaller),...
                min(Smaller,[],'all'),'post');
            
            if TempMultiplier >= 1
                In1 = Bigger;
                In2 = Smaller;
            else
                In1 = Smaller;
                In2 = Bigger;
            end
            
            OutImage = obj.deconvolute_by_mathematical_morphology(In1,In2);
            OutImage = imresize(OutImage,[Processed.NumPixelsX Processed.NumPixelsY]);
            
            Deconvoluted = Processed;
            Deconvoluted.Name = 'Deconvoluted';
            
            MinIn = min(Processed.Image,[],'all');
            MinOut = min(OutImage,[],'all');
            
            OutImage = OutImage + (MinIn - MinOut);
            
            Deconvoluted.Image = OutImage;
            
            [~,Index] = obj.get_channel('Deconvoluted');
            if isempty(Index)
                obj.Channel(end+1) = Deconvoluted;
            else
                obj.Channel(Index) = Deconvoluted;
            end
        end
        
        function clear_all_properties(obj)
            
            PropertyStruct = obj.get();
            PropertyNames = fieldnames(PropertyStruct);
            
            for i=1:length(PropertyNames)
                set(obj,PropertyNames{i},[])
            end
        end
        
        function construct_list_to_map_relations(obj)
            k = 1;
            obj.List2Map = zeros(obj.NumPixelsX*obj.NumPixelsY,2);
            if isequal(obj.FileType,'quantitative-imaging-map')
                for i=1:obj.NumPixelsX
                    for j=1:obj.NumPixelsY
                        obj.Map2List(i,j) = k;
                        obj.List2Map(k,:) = [i j];
                        k = k + 1;
                    end
                end
            elseif isequal(obj.FileType,'force-scan-map')
                for i=1:obj.NumPixelsX
                    if ~mod(i,2)
                        for j=1:obj.NumPixelsY
                            obj.Map2List(i,j) = k;
                            obj.List2Map(k,:) = [i j];
                            k = k + 1;
                        end
                    else
                        for j=1:obj.NumPixelsY
                            obj.Map2List(i,obj.NumPixelsY-j+1) = k;
                            obj.List2Map(k,:) = [i obj.NumPixelsY-j+1];
                            k = k + 1;
                        end
                    end
                end
            else
                for i=1:obj.NumPixelsX
                    for j=1:obj.NumPixelsY
                        obj.Map2List(i,j) = k;
                        obj.List2Map(k,:) = [i j];
                        k = k + 1;
                    end
                end
            end
        end
        
        function Map = convert_data_list_to_map(obj,List)
            
            Map = zeros(obj.NumPixelsX,obj.NumPixelsY);
            for i=1:obj.NumPixelsX
                for j=1:obj.NumPixelsY
                    Map(i,j) = List(obj.Map2List(i,j));
                end
            end
            
        end
        
        function [ChannelStruct,Index] = get_channel(obj,ChannelName)
            k = 0;
            for i=1:length(obj.Channel)
                if isequal(obj.Channel(i).Name,ChannelName)
                    ChannelStruct = obj.Channel(i);
                    Index = i;
                    k = k+1;
                end
            end
            if k > 1
                warning(sprintf('Caution! There are more than one channels named %s (%i)',ChannelName,k))
            end
            if k == 0
                ChannelStruct = [];
                Index = [];
            end
        end
        
        function [Channel,Index,FoundRequested] = get_unprocessed_height_channel(obj,ChannelName)
            % Goes over all possible variations of height channels in a
            % certain order and spits out an alternative one should the
            % specified one not exist
            
            HeightChannelList = {'Height (measured) (Trace)',...
                                    'Height (Trace)',...
                                    'Height (measured)',...
                                    'Height',...
                                    'Height (measured) (Retrace)',...
                                    'Height (Retrace)'};
            
            [Channel,Index] = obj.get_channel(ChannelName);
            if ~isempty(Channel)
                FoundRequested = true;
                return
            end
            
            for i=1:length(HeightChannelList)
                [Channel,Index] = obj.get_channel(HeightChannelList{i});
                if ~isempty(Channel)
                    warning(sprintf('Could not find a Channel named "%s", loaded Channel "%s" instead',ChannelName,HeightChannelList{i}))
                    FoundRequested = false;
                    return
                end
            end
            
            if isempty(Channel)
                error('No Channel found')
            end
        end
        
        function delete_channel(obj,ChannelName)
            k = 0;
            Index = [];
            for i=1:length(obj.Channel)
                if isequal(obj.Channel(i).Name,ChannelName)
                    ChannelStruct = obj.Channel(i);
                    Index(k+1) = i;
                    k = k+1;
                end
            end
            if isempty(Index)
                return
            end
            obj.Channel(Index) = [];
            if k > 1
                warning(sprintf('Caution! There are more than one channels named %s (%i). Deleted all of them',ChannelName,k))
            end
        end
        
        function OutChannel = create_standard_channel(obj,Image,Name,Unit)
            
            OutChannel.Image = Image;
            OutChannel.Name = Name;
            OutChannel.Unit = Unit;
            OutChannel.ScanSizeX = obj.ScanSizeX;
            OutChannel.ScanSizeY = obj.ScanSizeY;
            OutChannel.ScanAngle = obj.ScanAngle;
            OutChannel.NumPixelsX = obj.NumPixelsX;
            OutChannel.NumPixelsY = obj.NumPixelsY;
            OutChannel.OriginX = obj.OriginX;
            OutChannel.OriginY = obj.OriginY;
            
        end
        
        function set_channel_positions(obj,OriginX,OriginY,ScanAngle)
            
            NumChan = length(obj.Channel);
            
            for i=1:NumChan
                obj.Channel(i).OriginX = OriginX;
                obj.Channel(i).OriginY = OriginY;
                obj.Channel(i).ScanAngle = ScanAngle;
            end
        end
        
        function create_pixel_difference_channel(obj,useSlowScanDirBool)
            
            obj.delete_channel('Pixel Difference');
            
            InChannel = obj.get_unprocessed_height_channel('Processed');
            
            if nargin < 2
                useSlowScanDirBool = false;
            end
            
            InImage = InChannel.Image;
            
            if useSlowScanDirBool
                InImage = imrotate(InImage,90);
            end
            
            OutImage = AFMImage.create_pixel_difference_map(InImage);
            
            OutChannel = obj.create_standard_channel(OutImage,'Pixel Difference','m');
            
            obj.Channel(end+1) = OutChannel;
        end
        
        function PopUp = string_of_existing(obj)
            PopUp{1} = 'none';
            for i=1:length(obj.Channel)
                PopUp{i+1} = obj.Channel(i).Name;
            end
        end
        
        function snap_line_segments_to_local_perpendicular_maximum(obj,SampleDistanceMeters,WidthLocalWindowMeters,SmoothingWindowSize,SmoothingWindowWeighting,Indizes)
            % snap_line_segments_to_local_perpendicular_maximum(obj,SampleDistanceMeters,WidthLocalWindowMeters,SmoothingWindowSize,SmoothingWindowWeighting)
            %
            % SampleDistanceMeters ...(def=50e-9) determines the distance
            %   between polyline vertices created for snap-in
            % WidthLocalWindowMeters ...(def=300e-9) determines the full
            %   Length of the perpendicular line from which the local
            %   profile is drawn
            % SmoothingWindowSize ...(def=1) determines how many vertices
            %   are used for several smoothing steps during the process
            % SmoothingWindowWeighting ...(def='flat') determines smoothing
            %   window weighting; possible options:
            %         'flat','linear','gaussian','localregN' with N being a
            %         number between 1 and 9, representing an Nth grade
            %         polynomial fit
            % Indizes ... Indizes of Segments to be snapped. Snap all if
            %               empty
            
            if nargin < 6
                Indizes = [];
            end
            if nargin < 5
                Indizes = [];
                SmoothingWindowWeighting = 'flat';
            end
            if nargin < 4
                Indizes = [];
                SmoothingWindowSize = 21;
                SmoothingWindowWeighting = 'flat';
            end
            if nargin < 3
                Indizes = [];
                SampleDistanceMeters = 50*1e-9;
                SmoothingWindowSize = 21;
                SmoothingWindowWeighting = 'flat';
            end
            if nargin < 2
                Indizes = [];
                SampleDistanceMeters = 50*1e-9;
                WidthLocalWindowMeters = 300*1e-9;
                SmoothingWindowSize = 21;
                SmoothingWindowWeighting = 'flat';
            end
            
            if isempty(Indizes)
                for i=1:length(obj.Segment)
                    if isempty(strfind(obj.Segment(i).Name,'Snapped'))
                        Indizes(end+1) = i;
                    end
                end
            end
            
            NumSnaps = length(Indizes);
            h = waitbar(0,sprintf('Preparing to snap %i Segments',NumSnaps));
            
            Channel = obj.get_channel('Processed');
            if isempty(Channel)
                warning('No Channel "Processed" found. Need flattened height channel for snap to local maximum')
                return
            end
            
            % If Channel has an unequal SizePerPixel, resize it for
            % following operations
            [Channel,XMult,YMult] = obj.resize_channel_to_same_size_per_pixel(Channel);
            
            CenterIndex = floor(SmoothingWindowSize/2+1);
            
            % Convert inputs from meters to pixels
            SizePerPixel = Channel.ScanSizeX/Channel.NumPixelsX;
            WidthLocalWindowPixels = WidthLocalWindowMeters/SizePerPixel;
            SampleDistancePixels = SampleDistanceMeters/SizePerPixel;
            
            % Loop over all Segments and create 'Snapped' segments for all
            % Polylines
            k = 1;
            for i=Indizes
                if isequal(obj.Segment(i).Type,'polyline') && isempty(strfind(obj.Segment(i).Name,'Snapped'))
                    Snapped(k) = obj.Segment(i);
                    Snapped(k).Name = ['Snapped - ' obj.Segment(i).Name];
                    NewVertices = [];
                    OriginalVertices = obj.Segment(i).ROIObject.Position;
                    OriginalVertices = [OriginalVertices(:,1).*YMult OriginalVertices(:,2).*XMult];
                    while size(OriginalVertices,1) > 1
                        Length = norm(OriginalVertices(1,:) - OriginalVertices(2,:));
                        N = round(Length/SampleDistancePixels);
                        if N <= 2 
                            N = 2;
                        end
                        t = linspace(0,1,N)';
                        TempNewVertices = (1-t)*OriginalVertices(1,:) + t*OriginalVertices(2,:);
                        if length(OriginalVertices(:,1)) > 2
                            TempNewVertices(end,:) = [];
                        end
                        NewVertices = [NewVertices; TempNewVertices];
                        OriginalVertices(1,:) = [];
                    end
                    OutOfBounds = or(or(NewVertices(:,1) < 1,NewVertices(:,2) < 1),...
                        or(NewVertices(:,1) > Channel.NumPixelsX,NewVertices(:,2) > Channel.NumPixelsX));
                    NewVertices(OutOfBounds,:) = [];
                    if size(NewVertices,1) < 3
                        continue
                    end
                    Snapped(k).ROIObject.Position = NewVertices;
                    k = k + 1;
                end
            end
            if k == 1
                disp('No polyline segments found')
                return
            end
            
            for i=1:length(Snapped)
                waitbar(i/NumSnaps,h,sprintf('Snapping %s %s',obj.Segment(i).Name,obj.Segment(i).SubSegmentName))
                SnappedPos = [];
                SmoothedSnappedPos = [];
                SnappedToOriginalDistance = [];
                PerpendicularVector = [];
                OriginalPos = [];
                FinalSnappedPos = [];
                [LocalDirectionVector,WeightingVector] = AFMBaseClass.find_local_direction_vector_in_ordered_vector_list(...
                    Snapped(i).ROIObject.Position,SmoothingWindowSize,'flat');
                for j=1:size(Snapped(i).ROIObject.Position,1)
                    OriginalPos(j,:) = Snapped(i).ROIObject.Position(j,:);
                    PerpendicularVector(j,:) = [LocalDirectionVector(j,2) -LocalDirectionVector(j,1)]/norm(LocalDirectionVector(j,:));
                    WindowStart = Snapped(i).ROIObject.Position(j,:) + PerpendicularVector(j,:).*WidthLocalWindowPixels/2;
                    WindowEnd = Snapped(i).ROIObject.Position(j,:) - PerpendicularVector(j,:).*WidthLocalWindowPixels/2;
                    [LocalX,LocalY,LocalProfile] = improfile(Channel.Image,[WindowStart(1) WindowEnd(1)],[WindowStart(2) WindowEnd(2)]);
                    if length(LocalProfile) >= 4 && ~sum(isnan(LocalProfile))
                        LocalProfile(1) = min(LocalProfile);
                        LocalProfile(end) = min(LocalProfile);
                        LocalX = interp1(LocalX,linspace(1,length(LocalX),100)','spline');
                        LocalY = interp1(LocalY,linspace(1,length(LocalY),100)','spline');
                        LocalProfile = interp1(LocalProfile,linspace(1,length(LocalProfile),100)','spline');
                    end
                    [~,MaxIndex] = max(LocalProfile);
                    SnappedPos(j,:) = [LocalX(MaxIndex) LocalY(MaxIndex)];
                    DisplacementVector(j,:) = (SnappedPos(j,:) - OriginalPos(j,:));
                    SnappedToOriginalDistance(j) = PerpendicularVector(j,:)*DisplacementVector(j,:)';
%                     %                     Debug Start
%                     if j==1
%                         f = figure;
%                     end
%                     if j>1
%                         subplot(2,1,1)
%                         scatter([WindowStart(1) WindowEnd(1)],[WindowStart(2) WindowEnd(2)])
%                         Ax = gca;
%                         hold on
%                         plot(Snapped(i).ROIObject.Position(:,1),Snapped(i).ROIObject.Position(:,2))
%                         AutoX = Ax.XLim;
%                         AutoY = Ax.YLim;
%                         MaxDiff = max(diff(AutoX),diff(AutoY));
%                         xlim([AutoX(1) AutoX(1)+MaxDiff])
%                         ylim([AutoY(1) AutoY(1)+MaxDiff])
%                         plot(SnappedPos(:,1),SnappedPos(:,2))
%                         subplot(2,1,2)
%                         plot(LocalProfile)
%                         hold on
%                         plot(MaxIndex,LocalProfile(MaxIndex),'gv','MarkerEdgeColor','g','MarkerFaceColor','r','MarkerSize',12)
%                         drawnow
%                         hold off
%                     end
%                     %                     Debug end
                end
                SnappedToOriginalDistance = SnappedToOriginalDistance';
                % Apply smoothing to SnappedPos
                if ~strfind(SmoothingWindowWeighting,'localreg')
                    SmoothedSnappedDistance = zeros(size(SnappedToOriginalDistance));
                    for j=1:size(SnappedToOriginalDistance,1)
                        LowerIndex = max(j - (CenterIndex-1),1);
                        UpperIndex = min(j + (CenterIndex-1),size(SnappedToOriginalDistance,1));
                        SmoothedSnappedDistance(j) = sum(SnappedToOriginalDistance(LowerIndex:UpperIndex).*WeightingVector(CenterIndex - (j-LowerIndex):CenterIndex + (UpperIndex-j)),1);
                        FinalSnappedPos(j,:) = OriginalPos(j,:) + PerpendicularVector(j,:).*SmoothedSnappedDistance(j);
                    end
                else
                    for j=1:size(SnappedToOriginalDistance,1)
                        LowerIndex = max(j - (CenterIndex-1),1);
                        UpperIndex = min(j + (CenterIndex-1),size(SnappedToOriginalDistance,1));
                        IndexOI = (j+1) - LowerIndex;
                        LocalPoints = SnappedPos(LowerIndex:UpperIndex,:);
                        Means = mean(LocalPoints,1);
                        CenteredLP = LocalPoints - Means;
                        % replace outliers for robustness
%                         CenteredLP = filloutliers(CenteredLP,'linear',1);
                        [PCACoeff,TransformedSnappedPos,Latent] = pca(CenteredLP);
                        DegreeOfPolyfit = str2num(SmoothingWindowWeighting(end));
                        if isempty(DegreeOfPolyfit)
                            DegreeOfPolyfit = 2;
                        end
                        LocalPFit = polyfit(TransformedSnappedPos(:,1),TransformedSnappedPos(:,2),DegreeOfPolyfit);
                        SmoothedPos = [TransformedSnappedPos(IndexOI,1) polyval(LocalPFit,TransformedSnappedPos(IndexOI,1))];
                        FinalSnappedPos(j,:) = ((PCACoeff')\SmoothedPos')' + Means;
                    end
                end
                
                FinalSnappedPos(:,1) = FinalSnappedPos(:,1)./YMult;
                FinalSnappedPos(:,2) = FinalSnappedPos(:,2)./XMult;
                Snapped(i).ROIObject.Position = FinalSnappedPos;
            end
            
            for i=1:length(Snapped)
                ToReplaceIndex = find(strcmp({Snapped(i).Name},{obj.Segment.Name}) & strcmp({Snapped(i).SubSegmentName},{obj.Segment.SubSegmentName}));
                if ~isempty(ToReplaceIndex)
                    obj.Segment(ToReplaceIndex) = Snapped(i);
                else
                    obj.Segment(end+1) = Snapped(i);
                end
            end
            
            close(h)
        end
        
        function apply_segmentation_to_other_baseclass(DonorClass,AcceptorClass)
            
            AcceptorClass.Segment = DonorClass.Segment;
            
            % Loop over Segments
            for i=1:length(DonorClass.Segment)
                for j=1:size(DonorClass.Segment(i).ROIObject.Position,1)
                    [Xout,Yout] = DonorClass.transform_pixels_to_other_coordinate_basis(AcceptorClass,...
                        DonorClass.Segment(i).ROIObject.Position(j,1),...
                        DonorClass.Segment(i).ROIObject.Position(j,2));
                    AcceptorClass.Segment(i).ROIObject.Position(j,1) = Xout;
                    AcceptorClass.Segment(i).ROIObject.Position(j,2) = Yout;
                end
            end
            
        end
        
        function [Xout,Yout] = transform_pixels_to_other_coordinate_basis(DonorClass,AcceptorClass,X,Y)
            
            XDiff = DonorClass.Channel(1).OriginX - AcceptorClass.Channel(1).OriginX;
            SizePerPixelX = AcceptorClass.Channel(1).ScanSizeX./AcceptorClass.Channel(1).NumPixelsX;
            XDiff = XDiff/SizePerPixelX;
            YDiff = DonorClass.Channel(1).OriginY - AcceptorClass.Channel(1).OriginY;
            SizePerPixelY = AcceptorClass.Channel(1).ScanSizeY./AcceptorClass.Channel(1).NumPixelsY;
            YDiff = YDiff/SizePerPixelY;
            AngleDiff = DonorClass.Channel(1).ScanAngle - AcceptorClass.Channel(1).ScanAngle;
            AngleDiff = deg2rad(-AngleDiff);
            
            Vector = [X Y];
            
            ImCenter = [AcceptorClass.Channel(1).NumPixelsX/2 AcceptorClass.Channel(1).NumPixelsY/2];
            
            TempVector = Vector - ImCenter;
            
            RotationMatrix = [cos(AngleDiff) -sin(AngleDiff);sin(AngleDiff) cos(AngleDiff)];
            
            Vector = [RotationMatrix*TempVector']' + ImCenter + [XDiff -YDiff];
            
            Xout = Vector(1);
            Yout = Vector(2);
            
        end
        
        function [OutMat,OutConcatCell,OutCell,OutMask,SegmentNameList,FullSegmentNameList]...
                = get_segment_data_from_channel(obj,ChannelName,varargin)
                % function [OutMat,OutConcatCell,OutCell,OutMask,SegmentNameList,FullSegmentNameList]...
                % = get_segment_data_from_channel(obj,ChannelName,varargin)
                %
                % Reads out data points from under an objects Segments from
                % any channel existing in obj.Channel
                %
                %
                % Required inputs
                % obj ... AFMBaseClass Object with valid nonempty Segment
                %          properties
                % ChannelName ... Name of ChannelStruct in property Channel
                %          of obj
                %
                % Name-Value pairs
                % "PixelDilation" ... (def=0) Integer>=0. Not only output pixels from
                %               directly under the segment line but widened
                %               by PixelDilation pixels.
                % "JustSnapped" ...(def=false) Only use Data from Segments with the
                %               matched string 'Snapped'. Has priority over
                %               JustUnsnapped
                % "JustUnsnapped" ...(def=false) Only use Data from Segments without
                %               the match string 'Snapped'
                % "MatchString" ...(def=[]) Only take Data from Segments matching
                %               the string MatchString. Multiple matches
                %               still return True
                % "IncludeIndexVector" ...(def=[]) Include segments with
                %                   Index contained in IncludeIndexVector.
                %                   Gets overwritten by any of the
                %                   discriminators above.
                % "ExcludeIndexVector" ...(def=[]) Exclude segments with
                %                   Index contained in ExcludeIndexVector.
                %                   Gets overwritten by IncludeIndexVector.
                
                p = inputParser;
                p.FunctionName = "get_segment_data_from_channel";
                p.CaseSensitive = false;
                p.PartialMatching = true;
                
                % Required inputs
                validobj = @(x)true;
                validChannelName = @(x)ischar(x);
                addRequired(p,"obj",validobj);
                addRequired(p,"ChannelName",validChannelName);
                
                % NameValue inputs
                defaultPixelDilation = 0;
                defaultJustSnapped = false;
                defaultJustUnsnapped = false;
                defaultMatchString = [];
                defaultIncludeIndexVector = [];
                defaultExcludeIndexVector = [];
                validPixelDilation = @(x)isscalar(x)&&x>=0;
                validJustSnapped = @(x)islogical(x)||x==0||x==1;
                validJustUnsnapped = @(x)islogical(x)||x==0||x==1;
                validMatchString = @(x)ischar(x)||isempty(x);
                validIncludeIndexVector = @(x)isnumeric(x)&&(size(x,1)==1||size(x,2)==1)||isempty(x);
                validExcludeIndexVector = @(x)isnumeric(x)&&(size(x,1)==1||size(x,2)==1)||isempty(x);
                addParameter(p,"PixelDilation",defaultPixelDilation,validPixelDilation);
                addParameter(p,"JustSnapped",defaultJustSnapped,validJustSnapped);
                addParameter(p,"JustUnsnapped",defaultJustUnsnapped,validJustUnsnapped);
                addParameter(p,"MatchString",defaultMatchString,validMatchString);
                addParameter(p,"IncludeIndexVector",defaultIncludeIndexVector,validIncludeIndexVector);
                addParameter(p,"ExcludeIndexVector",defaultExcludeIndexVector,validExcludeIndexVector);
                
                parse(p,obj,ChannelName,varargin{:});
                
                % Assign parsing results to named variables
                obj = p.Results.obj;
                ChannelName = p.Results.ChannelName;
                PixelDilation = round(p.Results.PixelDilation);
                JustSnapped = p.Results.JustSnapped;
                JustUnsnapped = p.Results.JustUnsnapped&&~JustSnapped;
                MatchString = p.Results.MatchString;
                IncludeIndexVector = p.Results.IncludeIndexVector;
                ExcludeIndexVector = setdiff(p.Results.ExcludeIndexVector,IncludeIndexVector);
                
                
                
            if isempty(obj.Segment) || isempty(obj.Segment(1).Name)
                warning(sprintf('%s has no segementation data',obj.Name))
            end
            
            OutMat = [];
            OutConcatCell = {};
            OutCell = {};
            
            Channel = obj.get_channel(ChannelName);
            Data = Channel.Image;
            OutMask = zeros(size(Data));
            
            SegmentNames = unique({obj.Segment.Name});
            set(0,'DefaultFigureVisible','off');
            f = figure;
            I = imshow(Data);
            Parent = I.Parent;
            m = 1;
            for i=1:length(SegmentNames)
                k = 1;
                for j=1:length(obj.Segment)
                    if isequal(SegmentNames{i},obj.Segment(j).Name)
                        if any(find(ExcludeIndexVector==j))
                            continue
                        end
                        if ~isempty(IncludeIndexVector)&&~any(find(IncludeIndexVector==j))
                            continue
                        end
                        if ~isempty(MatchString)&&~any(strfind(obj.Segment(j).Name,MatchString))
                            continue
                        end
                        if JustSnapped&&~any(strfind(obj.Segment(j).Name,'Snapped'))
                            continue
                        end
                        if JustUnsnapped&&any(strfind(obj.Segment(j).Name,'Snapped'))
                            continue
                        end
                        gca = Parent;
                        Line = drawpolyline('Position',obj.Segment(j).ROIObject.Position,'Parent',Parent);
                        Mask = Line.createMask;
                        if PixelDilation
                            StrEl = strel('Disk',PixelDilation,0);
                            Mask = imdilate(Mask,StrEl);
                        end
                        OutMask = OutMask | Mask;
                        SubSegmentPoints = Data(Mask == 1);
                        OutCell{i}{k} = SubSegmentPoints;
                        SegmentNameList{m} = obj.Segment(j).Name;
                        FullSegmentNameList{m} = [obj.Segment(j).Name ' || ' obj.Segment(j).SubSegmentName];
                        m = m + 1;
                        k = k + 1;
                    end
                end
            end
            close(f)
            set(0,'DefaultFigureVisible','on');
            
            SegmentNameList = unique(SegmentNameList);
            FullSegmentNameList = unique(FullSegmentNameList);
            
            for i=length(OutCell):-1:1
                if isempty(OutCell{i})
                    OutCell(i) = [];
                end
            end
            if isempty(OutCell)
                return
            end
            
            for i=1:length(OutCell)
                ConcVec = [];
                for j=1:length(OutCell{i})
                    Temp = cell2mat(OutCell{i}(j));
                    ConcVec = cat(1,reshape(Temp,[],1),reshape(ConcVec,[],1));
                end
                OutConcatCell{i} = ConcVec;
                SegmentLength(i) = length(ConcVec);
            end
            
            OutMat = nan(max(SegmentLength),length(OutCell));
            for i=1:length(OutCell)
                OutMat(1:SegmentLength(i),i) = OutConcatCell{i};
            end
            
        end
        
        function create_proximity_mapping_for_segments(obj,TargetObject)
            
            NumSegments = length(obj.Segment);
            
            for i=1:NumSegments
                OriginalSeg = obj.Segment(i);
                [MatchedSeg,MatchedIndex] = AFMBaseClass.find_matching_segment(OriginalSeg,TargetObject.Segment);
                if isempty(MatchedIndex)
                    continue
                end
                % Coordinate Transform MatchedSeg
                for j=1:length(MatchedSeg.ROIObject.Position(:,1))
                    [MatchedSeg.ROIObject.Position(j,1),MatchedSeg.ROIObject.Position(j,2)] = ...
                        TargetObject.transform_pixels_to_other_coordinate_basis(obj,...
                        MatchedSeg.ROIObject.Position(j,1),MatchedSeg.ROIObject.Position(j,2));
                end
                [Seg1,Seg2] = AFMBaseClass.proximity_map_two_segmentation_elements(OriginalSeg,MatchedSeg);
                Seg1.ProximityMap.Names = {obj.Name TargetObject.Name};
                Seg2.ProximityMap.Names = {obj.Name TargetObject.Name};
                if isfield(obj.Segment(i),'ProximityMap')
                    obj.Segment(i).ProximityMap(end+1) = Seg1.ProximityMap;
                else
                    obj.Segment(i).ProximityMap = Seg1.ProximityMap;
                end
                
                if isfield(TargetObject.Segment(i),'ProximityMap')
                    TargetObject.Segment(MatchedIndex).ProximityMap(end+1) = Seg2.ProximityMap;
                else
                    TargetObject.Segment(MatchedIndex).ProximityMap = Seg2.ProximityMap;
                end
            end
            
        end
        
        function reset_proximity_mapping(obj)
            
            for i=1:length(obj.Segment)
                obj.Segment(i).ProximityMap = [];
            end
            
        end
        
        function reset_property_to_default(obj,PropertyName)
            
            Meta = metaclass(obj);
            Prop = Meta.PropertyList;
            PropNames = {Prop.Name};
            Index = 0;
            for i=1:length(PropNames)
                if isequal(PropNames{i},PropertyName)
                    Index = i;
                end
            end
            if ~Index
                warning([PropertyName ' is not a Property'])
                return
            end
            if Prop(Index).HasDefault
                set(obj,PropertyName,Prop(Index).DefaultValue);
            else
                warning(['Property ' PropertyName ' does not have a default value!'])
            end
            
        end
        
        function automatic_segmentation_on_singular_vertical_fiber(obj,...
                SampleDistanceMeters,WidthLocalWindowMeters,...
                SmoothingWindowSize,SmoothingWindowWeighting,Indizes)
            
            
            HeightChannel = obj.get_unprocessed_height_channel('Processed');
            
            if nargin < 2
                SampleDistanceMeters = 20e-9;
                WidthLocalWindowMeters = 200e-9;
                SmoothingWindowSize = 41;
                SmoothingWindowWeighting = 'localreg1';
            end
            
            HeightMap = HeightChannel.Image;
            
            [~,MaxIdx] = max(filloutliers(HeightMap,'linear','movmedian',ceil(HeightChannel.NumPixelsY/5),2),[],2);
            Points = [MaxIdx [1:size(HeightMap,1)]'];
            
            
            obj.Segment(1).Name = 'SingleFiber';
            obj.Segment(1).SubSegmentName = 'SubS-01';
            obj.Segment(1).Type = 'polyline';
            obj.Segment(1).ROIObject.Position = Points;
            obj.Segment(1).ROIObject.LineWidth = 2.25;
            
            obj.snap_line_segments_to_local_perpendicular_maximum(SampleDistanceMeters,WidthLocalWindowMeters,SmoothingWindowSize,SmoothingWindowWeighting)
            
        end
        
        function [OutArrayStruct,OutStruct,OutStructAll] = characterize_fiber_like_polyline_segments(obj,varargin)
            % function [OutArrayStruct,OutStruct,OutStructAll] = characterize_fiber_like_polyline_segments(obj,varargin)
            %
            % Takes in SNAPPED polyline-segmented fiber-like structures and
            % determines several characterstics such as Height, FWHM, Area.
            %
            % The outputs are as follows:
            %       OutArrayStruct ... All nx1 Segment characteristics put
            %           together in an nxm array with m the number of unique
            %           segments (subsegments are concatenated).
            %       OutStruct ... All available information, barring the
            %           distance vectors, but only for unique segments
            %           (subsegments are concatenated).
            %       OutStructAll ... All available information on every
            %           subsegment separately this is also the information
            %           that is saved as part of the Experiment as a
            %           substruct of AFMBaseClass.Segment
            %
            %
            % Required inputs
            % obj ... AFMBaseClass object that already has polyline
            %         segments and a 'Processed' channel
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
            defaultWidthLocalWindowMeters = 800e-9;
            defaultSmoothingWindowSize = 41;
            defaultMinPeakDistanceMeters = 50e-9;
            defaultLowerEndThreshold = .1;
            defaultThresholdType = 'Fraction';
            defaultVerbose = false;
            defaultRecordMovieBool = false;
            defaultKeyFrames = 3;
            validWidthLocalWindowMeters = @(x)isnumeric(x)&&isscalar(x);
            validSmoothingWindowSize = @(x)isnumeric(x)&&mod(x,1)==0;
            validMinPeakDistanceMeters = @(x)isnumeric(x)&&isscalar(x);
            validLowerEndThreshold = @(x)isnumeric(x)&&isscalar(x);
            validThresholdType = @(x)any(validatestring(x,{'Fraction','Meters'}));
            validVerbose = @(x)islogical(x);
            validRecordMovieBool = @(x)islogical(x);
            validKeyFrames = @(x)isnumeric(x)&&mod(x,1)==0;
            addParameter(p,"WidthLocalWindowMeters",defaultWidthLocalWindowMeters,validWidthLocalWindowMeters);
            addParameter(p,"SmoothingWindowSize",defaultSmoothingWindowSize,validSmoothingWindowSize);
            addParameter(p,"MinPeakDistanceMeters",defaultMinPeakDistanceMeters,validMinPeakDistanceMeters);
            addParameter(p,"LowerEndThreshold",defaultLowerEndThreshold,validLowerEndThreshold);
            addParameter(p,"ThresholdType",defaultThresholdType,validThresholdType);
            addParameter(p,"Verbose",defaultVerbose,validVerbose);
            addParameter(p,"RecordMovieBool",defaultRecordMovieBool,validRecordMovieBool);
            addParameter(p,"KeyFrames",defaultKeyFrames,validKeyFrames);
            
            parse(p,obj,varargin{:});
            
            % Assign parsing results to named variables
            obj = p.Results.obj;
            WidthLocalWindowMeters = p.Results.WidthLocalWindowMeters;
            SmoothingWindowSize = p.Results.SmoothingWindowSize;
            MinPeakDistanceMeters = p.Results.MinPeakDistanceMeters;
            LowerEndThreshold = p.Results.LowerEndThreshold;
            ThresholdType = p.Results.ThresholdType;
            Verbose = p.Results.Verbose;
            RecordMovieBool = p.Results.RecordMovieBool;
            KeyFrames = p.Results.KeyFrames;
            
            OutArrayStruct = struct();
            OutStruct = struct();
            OutStructAll = struct();
            Struct = struct([]);
            
            Channel = obj.get_channel('Processed');
            if isempty(Channel)
                warning('No Channel "Processed" found. Need flattened height channel for snap to local maximum')
                return
            end
            
            % Convert inputs from meters to pixels
            SizePerPixel = Channel.ScanSizeX/Channel.NumPixelsX;
            WidthLocalWindowPixels = WidthLocalWindowMeters/SizePerPixel;
            MinPeakDistancePixels = MinPeakDistanceMeters/SizePerPixel;
            
            h = waitbar(0,sprintf('Preparing to analyze %i Segments',length(obj.Segment)));
            
            Size = size(Channel.Image);
            SegmentMasks = zeros(Size(1),Size(2),length(obj.Segment));
            
            f = figure('Color','w','Units','normalized','Position',[.2 .1 .6 .8]);
            if RecordMovieBool
                v = VideoWriter([obj.Folder,'FibrilCharacterization.avi'],'Motion JPEG AVI');
                v.Quality = 95;
                v.FrameRate = 60;
                v.open
            end
            subplot(1,2,1)
            Ax = imshow(Channel.Image,[]);
            for i=1:length(obj.Segment)
                if ~isequal(obj.Segment(i).Type,'polyline') || isempty(strfind(obj.Segment(i).Name,'Snapped'))
                    continue
                end
                
                Polyline = drawpolyline(Ax.Parent,'Position',obj.Segment(i).ROIObject.Position);
                FirstMask = Polyline.createMask;
                TempMask = FirstMask;
                PixelDil = 1;
                while PixelDil <= floor(WidthLocalWindowPixels)
                    Dilated = imdilate(FirstMask,strel('disk',PixelDil));
                    TempMask = TempMask + Dilated;
                    PixelDil = PixelDil + 1;
                end
                SegmentMasks(:,:,i) = TempMask;
                if Verbose
                    subplot(1,2,2)
                    imshowpair(Channel.Image,max(SegmentMasks,[],3))
                    drawnow
                    if RecordMovieBool
                        Frame = getframe(f);
                        v.writeVideo(Frame)
                    end
                end
            end
            
            if ~Verbose
                close(f)
            end
            
            k = 1;
            for i=1:length(obj.Segment)
                waitbar(i/length(obj.Segment),h,sprintf('Analyzing %s %s',obj.Segment(i).Name,obj.Segment(i).SubSegmentName))
                if ~isequal(obj.Segment(i).Type,'polyline') || isempty(strfind(obj.Segment(i).Name,'Snapped'))
                    continue
                end
                IgnoreMask = zeros(Size);
                for j=1:length(obj.Segment)
                    if j==i
                        continue
                    end
                    IgnoreMask = or(IgnoreMask,(SegmentMasks(:,:,j) > SegmentMasks(:,:,i)));
                end
                SegmentPositions = obj.Segment(i).ROIObject.Position;
                LocalDirectionVector = AFMBaseClass.find_local_direction_vector_in_ordered_vector_list(obj.Segment(i).ROIObject.Position,SmoothingWindowSize,'flat');
                LocalHeight = zeros(length(SegmentPositions(:,1)),1);
                LocalProminence = zeros(length(SegmentPositions(:,1)),1);
                LocalWidthHalfHeight = zeros(length(SegmentPositions(:,1)),1);
                LocalWidthHalfProminence = zeros(length(SegmentPositions(:,1)),1);
                LocalArea = zeros(length(SegmentPositions(:,1)),1);
                LocalWidthBase = zeros(length(SegmentPositions(:,1)),1);
                for j=1:length(SegmentPositions(:,1))
                    PerpendicularVector(j,:) = [LocalDirectionVector(j,2) -LocalDirectionVector(j,1)]/norm(LocalDirectionVector(j,:));
                    WindowStart =SegmentPositions(j,:) + PerpendicularVector(j,:).*WidthLocalWindowPixels/2;
                    WindowEnd = SegmentPositions(j,:) - PerpendicularVector(j,:).*WidthLocalWindowPixels/2;
                    [LocalX,LocalY,LocalProfile] = improfile(Channel.Image,[WindowStart(1) WindowEnd(1)],[WindowStart(2) WindowEnd(2)]);
                    if length(LocalProfile) >= 4 && ~sum(isnan(LocalProfile))
                        LocalX = interp1(LocalX,linspace(0,1,100)'.*length(LocalX),'spline');
                        LocalY = interp1(LocalY,linspace(0,1,100)'.*length(LocalY),'spline');
                        LocalProfile = interp1(LocalProfile,linspace(0,1,100)'.*length(LocalProfile),'spline');
                    end
                    LocalDistance = sign(LocalX - SegmentPositions(j,1)).*vecnorm([LocalX LocalY] - SegmentPositions(j,:),2,2);
                    LocalDistance = sign(LocalDistance(end) - LocalDistance(1)).*LocalDistance;
                    LocalDistance = LocalDistance.*WidthLocalWindowMeters/range(LocalDistance);
                    ForbiddenX = LocalX;
                    ForbiddenY = LocalY;
                    ForbiddenProfile = LocalProfile;
                    ForbiddenLocalDistance = LocalDistance;
                    isForbidden = true(size(LocalX));
                    for jj=length(ForbiddenX):-1:1
                        if (round(ForbiddenX(jj))<1 || round(ForbiddenX(jj))>Size(1)) ||...
                                (round(ForbiddenY(jj))<1 || round(ForbiddenY(jj))>Size(2)) 
                            continue
                        end
                        if ~IgnoreMask(round(ForbiddenY(jj)),round(ForbiddenX(jj)))
                            ForbiddenX(jj) = [];
                            ForbiddenY(jj) = [];
                            ForbiddenProfile(jj) = [];
                            ForbiddenLocalDistance(jj) = [];
                            isForbidden(jj) = false;
                        end
                    end
                    [Heights,Locations,TempWidthHalfHeight,Prom] = findpeaks(LocalProfile,LocalDistance,...
                        'WidthReference','halfheight','MinPeakDistance',MinPeakDistanceMeters);
                    [~,~,TempWidthHalfProm] = findpeaks(LocalProfile,LocalDistance,...
                        'WidthReference','halfprom','MinPeakDistance',MinPeakDistanceMeters);
                    % Assign a score to the quality of the peak and 
                    TempScore = ~ismember(Locations,ForbiddenLocalDistance).*Heights.*Prom./abs(Locations).^2;
                    [~,WinnerIndex] = max(TempScore);
                    
                    if Verbose && mod(j,KeyFrames)==0
                        subplot('Position',[0.05 0.05 .5 .25])
                        imshow(IgnoreMask,[])
                        title('Forbidden Areas')
                        drawline('Position',[WindowStart; WindowEnd]);
                        subplot('Position',[0.05 0.35 .5 .5])
                        imshow(rescale(Channel.Image),'Colormap',obj.CMap);
                        drawline('Position',[WindowStart; WindowEnd]);
                        Ax1 = subplot('Position',[.6 .1 .3 .4]);
                        findpeaks(LocalProfile,LocalDistance,'Annotate','extents',...
                            'WidthReference','halfheight','MinPeakDistance',MinPeakDistanceMeters);
                        hold on
                        plot(ForbiddenLocalDistance,ForbiddenProfile,'rX',...
                            Locations(WinnerIndex),1.05*Heights(WinnerIndex),'rv',...
                            'MarkerSize',16,'MarkerFaceColor','g')
                        hold off
                        Ax1.Legend.String{end-1} = 'Forbidden Points';
                        Ax1.Legend.String{end} = 'Chosen Point';
                        Ax2 = subplot('Position',[.6 .55 .3 .4]);
                        Ax1.Legend.Location = 'northeast';
                        findpeaks(LocalProfile,LocalDistance,'Annotate','extents',...
                            'WidthReference','halfprom','MinPeakDistance',MinPeakDistanceMeters);
                        hold on
                        plot(ForbiddenLocalDistance,ForbiddenProfile,'rX',...
                            Locations(WinnerIndex),1.05*Heights(WinnerIndex),'rv',...
                            'MarkerSize',16,'MarkerFaceColor','g')
                        Ax2.Legend.String{end-1} = 'Forbidden Points';
                        Ax2.Legend.String{end} = 'Chosen Point';
                        Ax2.Legend.Location = 'northeast';
                        hold off
                        drawnow
                        if RecordMovieBool
                            Frame = getframe(f);
                            v.writeVideo(Frame)
                        end
                    end
                    if isempty(Heights)
                        LocalHeight(j) = nan;
                        LocalProminence(j) = nan;
                        LocalWidthHalfHeight(j) = nan;
                        LocalWidthHalfProminence(j) = nan;
                        LocalArea(j) = nan;
                        LocalWidthBase(j) = nan;
                        continue
                    end
                    % find out local area
                    PeakIndex = find(LocalDistance==Locations(WinnerIndex));
                    if isequal(ThresholdType,'Fraction')
                        Thresh = LowerEndThreshold.*LocalProfile(PeakIndex);
                    else
                        Thresh = LowerEndThreshold;
                    end
                    % To the left
                    kk = 1;
                    LeftBoundIndex = [];
                    while (PeakIndex-kk)>0 && isempty(LeftBoundIndex)
                        if LocalProfile(PeakIndex-kk) < Thresh
                            LeftBoundIndex = PeakIndex - kk;
                        end
                        kk = kk + 1;
                    end
                    RightBoundIndex = [];
                    % To the right
                    kk = 1;
                    while (PeakIndex+kk)<=length(LocalProfile) && isempty(RightBoundIndex)
                        if LocalProfile(PeakIndex+kk) < Thresh
                            RightBoundIndex = PeakIndex + kk;
                        end
                        kk = kk + 1;
                    end
                    % If info is sufficient, calculate area and base width
                    if ~isempty(LeftBoundIndex) && ~isempty(RightBoundIndex)
                        LocalArea(j) = trapz(...
                            LocalDistance(LeftBoundIndex:RightBoundIndex),...
                            LocalProfile(LeftBoundIndex:RightBoundIndex));
                        LocalWidthBase(j) = abs(LocalDistance(LeftBoundIndex)...
                             - LocalDistance(RightBoundIndex));
                    else
                        LocalArea(j) = nan;
                        LocalWidthBase(j) = nan;
                    end
                    
                    % Write down results
                    LocalHeight(j) = Heights(WinnerIndex);
                    LocalProminence(j) = Prom(WinnerIndex);
                    LocalWidthHalfHeight(j) = TempWidthHalfHeight(WinnerIndex);
                    LocalWidthHalfProminence(j) = TempWidthHalfProm(WinnerIndex);
                end
                % find out overall segment pathlength
                for m=size(SegmentPositions,1):-1:1
                    if SegmentPositions(m,1)>Channel.NumPixelsX ||...
                            SegmentPositions(m,1)<0
                        SegmentPositions(m,:) = [];
                        continue
                    end
                    if SegmentPositions(m,2)>Channel.NumPixelsY ||...
                            SegmentPositions(m,2)<0
                        SegmentPositions(m,:) = [];
                    end
                end
                SegmentLength = sum(vecnorm(...
                    (SegmentPositions(1:end-1,:) - SegmentPositions(2:end,:))...
                    *SizePerPixel,2,2),'all');
                
                % write down results
                Struct(k).Name = obj.Segment(i).Name;
                Struct(k).Height = LocalHeight;
                obj.Segment(i).Height = LocalHeight;
                Struct(k).WidthHalfHeight = LocalWidthHalfHeight;
                obj.Segment(i).WidthHalfHeight = LocalWidthHalfHeight;
                Struct(k).Prominence = LocalProminence;
                obj.Segment(i).Prominence = LocalProminence;
                Struct(k).WidthHalfProminence = LocalWidthHalfProminence;
                obj.Segment(i).WidthHalfProminence = LocalWidthHalfProminence;
                Struct(k).DirectionVector = LocalDirectionVector;
                obj.Segment(i).DirectionVector = LocalDirectionVector;
                Struct(k).SegmentLength = SegmentLength;
                obj.Segment(i).SegmentLength = SegmentLength;
                Struct(k).Mean_WidthHalfHeight = nanmean(LocalWidthHalfHeight);
                obj.Segment(i).Mean_WidthHalfHeight = nanmean(LocalWidthHalfHeight);
                Struct(k).Median_WidthHalfHeight = nanmedian(LocalWidthHalfHeight);
                obj.Segment(i).Median_WidthHalfHeight = nanmedian(LocalWidthHalfHeight);
                Struct(k).Mean_Height = nanmean(LocalHeight);
                obj.Segment(i).Mean_Height = nanmean(LocalHeight);
                Struct(k).Median_Height = nanmedian(LocalHeight);
                obj.Segment(i).Median_Height = nanmedian(LocalHeight);
                Struct(k).Mean_WidthHalfProminence = nanmean(LocalWidthHalfProminence);
                obj.Segment(i).Mean_WidthHalfProminence = nanmean(LocalWidthHalfProminence);
                Struct(k).Median_WidthHalfProminence = nanmedian(LocalWidthHalfProminence);
                obj.Segment(i).Median_WidthHalfProminence = nanmedian(LocalWidthHalfProminence);
                Struct(k).Mean_Prominence = nanmean(LocalProminence);
                obj.Segment(i).Mean_Prominence = nanmean(LocalProminence);
                Struct(k).Median_Prominence = nanmedian(LocalProminence);
                obj.Segment(i).Median_Prominence = nanmedian(LocalProminence);
                Struct(k).Area = LocalArea;
                obj.Segment(i).Area = LocalArea;
                Struct(k).Mean_Area = nanmean(LocalArea);
                obj.Segment(i).Mean_Area = nanmean(LocalArea);
                Struct(k).Median_Area = nanmedian(LocalArea);
                obj.Segment(i).Median_Area = nanmedian(LocalArea);
                Struct(k).WidthBase = LocalWidthBase;
                obj.Segment(i).WidthBase = LocalWidthBase;
                Struct(k).Mean_WidthBase = nanmean(LocalWidthBase);
                obj.Segment(i).Mean_WidthBase = nanmean(LocalWidthBase);
                Struct(k).Median_WidthBase = nanmedian(LocalWidthBase);
                obj.Segment(i).Median_WidthBase = nanmedian(LocalWidthBase);
                Struct(k).AspectRatioHalfHeight = LocalHeight./LocalWidthHalfProminence;
                obj.Segment(i).AspectRatioHalfHeight = LocalHeight./LocalWidthHalfProminence;
                Struct(k).Mean_AspectRatioHalfHeight = nanmean(LocalHeight./LocalWidthHalfProminence);
                obj.Segment(i).Mean_AspectRatioHalfHeight = nanmean(LocalHeight./LocalWidthHalfProminence);
                Struct(k).Median_AspectRatioHalfHeight = nanmedian(LocalHeight./LocalWidthHalfProminence);
                obj.Segment(i).Median_AspectRatioHalfHeight = nanmedian(LocalHeight./LocalWidthHalfProminence);
                Struct(k).AspectRatioBaseHeight = LocalHeight./LocalWidthBase;
                obj.Segment(i).AspectRatioBaseHeight = LocalHeight./LocalWidthBase;
                Struct(k).Mean_AspectRatioBaseHeight = nanmean(LocalHeight./LocalWidthBase);
                obj.Segment(i).Mean_AspectRatioBaseHeight = nanmean(LocalHeight./LocalWidthBase);
                Struct(k).Median_AspectRatioBaseHeight = nanmedian(LocalHeight./LocalWidthBase);
                obj.Segment(i).Median_AspectRatioBaseHeight = nanmedian(LocalHeight./LocalWidthBase);
                Struct(k).AreaDerivedDiameter = sqrt(LocalArea/pi).*2;
                obj.Segment(i).AreaDerivedDiameter = sqrt(LocalArea/pi).*2;
                Struct(k).Mean_AreaDerivedDiameter = nanmean(sqrt(LocalArea/pi).*2);
                obj.Segment(i).Mean_AreaDerivedDiameter = nanmean(sqrt(LocalArea/pi).*2);
                Struct(k).Median_AreaDerivedDiameter = nanmedian(sqrt(LocalArea/pi).*2);
                obj.Segment(i).Median_AreaDerivedDiameter = nanmedian(sqrt(LocalArea/pi).*2);
                k = k + 1;
            end
            
            % Now for some tedious data handling. Life is pain            
            if isempty(Struct)
                close(h)
                return
            end
            
            NameList = {Struct.Name};
            Unique = unique(NameList);
            
            for i=1:length(Unique)
                OutStruct(i).Name = Unique{i};
                OutStruct(i).Height = [];
                OutStruct(i).WidthHalfHeight = [];
                OutStruct(i).Prominence = [];
                OutStruct(i).WidthHalfProminence = [];
                OutStruct(i).Area = [];
                OutStruct(i).WidthBase = [];
                OutStruct(i).AspectRatioHalfHeight = [];
                OutStruct(i).AspectRatioBaseHeight = [];
                OutStruct(i).AreaDerivedDiameter = [];
                TempSegmentLength = 0;
                for j=1:length(Struct)
                    if isequal(Unique{i},Struct(j).Name)
                        OutStruct(i).Height = vertcat(OutStruct(i).Height,Struct(j).Height);
                        OutStruct(i).WidthHalfHeight = vertcat(OutStruct(i).WidthHalfHeight,Struct(j).WidthHalfHeight);
                        OutStruct(i).Prominence = vertcat(OutStruct(i).Prominence,Struct(j).Prominence);
                        OutStruct(i).WidthHalfProminence = vertcat(OutStruct(i).WidthHalfProminence,Struct(j).WidthHalfProminence);
                        OutStruct(i).Area = vertcat(OutStruct(i).Area,Struct(j).Area);
                        OutStruct(i).WidthBase = vertcat(OutStruct(i).WidthBase,Struct(j).WidthBase);
                        OutStruct(i).AspectRatioHalfHeight = vertcat(OutStruct(i).AspectRatioHalfHeight,Struct(j).AspectRatioHalfHeight);
                        OutStruct(i).AspectRatioBaseHeight = vertcat(OutStruct(i).AspectRatioBaseHeight,Struct(j).AspectRatioBaseHeight);
                        OutStruct(i).AreaDerivedDiameter = vertcat(OutStruct(i).AreaDerivedDiameter,Struct(j).AreaDerivedDiameter);
                        TempSegmentLength = TempSegmentLength + Struct(j).SegmentLength;
                    end
                end
                OutStruct(i).SegmentLength = TempSegmentLength;
                OutStruct(i).Mean_Height = nanmean(OutStruct(i).Height);
                OutStruct(i).Median_Height = nanmedian(OutStruct(i).Height);
                OutStruct(i).Mean_WidthHalfHeight = nanmean(OutStruct(i).WidthHalfHeight);
                OutStruct(i).Median_WidthHalfHeight = nanmedian(OutStruct(i).WidthHalfHeight);
                OutStruct(i).Mean_Prominence = nanmean(OutStruct(i).Prominence);
                OutStruct(i).Median_Prominence = nanmedian(OutStruct(i).Prominence);
                OutStruct(i).Mean_WidthHalfProminence = nanmean(OutStruct(i).WidthHalfProminence);
                OutStruct(i).Median_WidthHalfProminence = nanmedian(OutStruct(i).WidthHalfProminence);
                OutStruct(i).Mean_Area = nanmean(OutStruct(i).Area);
                OutStruct(i).Median_Area = nanmedian(OutStruct(i).Area);
                OutStruct(i).Mean_WidthBase = nanmean(OutStruct(i).WidthBase);
                OutStruct(i).Median_WidthBase = nanmedian(OutStruct(i).WidthBase);
                OutStruct(i).Mean_AspectRatioHalfHeight = nanmean(OutStruct(i).AspectRatioHalfHeight);
                OutStruct(i).Median_AspectRatioHalfHeight = nanmedian(OutStruct(i).AspectRatioHalfHeight);
                OutStruct(i).Mean_AspectRatioBaseHeight = nanmean(OutStruct(i).AspectRatioBaseHeight);
                OutStruct(i).Median_AspectRatioBaseHeight = nanmedian(OutStruct(i).AspectRatioBaseHeight);
                OutStruct(i).Mean_AreaDerivedDiameter = nanmean(OutStruct(i).AreaDerivedDiameter);
                OutStruct(i).Median_AreaDerivedDiameter = nanmedian(OutStruct(i).AreaDerivedDiameter);
            end
            
            for i=1:length(OutStruct)
                Lengths(i) = length(OutStruct(i).Height);
            end
            MaxLength = max(Lengths);
            OutArrayStruct.Height = nan.*ones(MaxLength,length(OutStruct));
            OutArrayStruct.WidthHalfHeight = nan.*ones(MaxLength,length(OutStruct));
            OutArrayStruct.Prominence = nan.*ones(MaxLength,length(OutStruct));
            OutArrayStruct.Area = nan.*ones(MaxLength,length(OutStruct));
            OutArrayStruct.WidthBase = nan.*ones(MaxLength,length(OutStruct));
            OutArrayStruct.AspectRatioHalfHeight = nan.*ones(MaxLength,length(OutStruct));
            OutArrayStruct.AspectRatioBaseHeight = nan.*ones(MaxLength,length(OutStruct));
            OutArrayStruct.AreaDerivedDiameter = nan.*ones(MaxLength,length(OutStruct));
            
            for i=1:length(OutStruct)
                OutArrayStruct.Height(1:length(OutStruct(i).Height),i) = OutStruct(i).Height;
                OutArrayStruct.WidthHalfHeight(1:length(OutStruct(i).Height),i) = OutStruct(i).WidthHalfHeight;
                OutArrayStruct.Prominence(1:length(OutStruct(i).Height),i) = OutStruct(i).Prominence;
                OutArrayStruct.WidthHalfProminence(1:length(OutStruct(i).Height),i) = OutStruct(i).WidthHalfProminence;
                OutArrayStruct.Area(1:length(OutStruct(i).Height),i) = OutStruct(i).Area;
                OutArrayStruct.WidthBase(1:length(OutStruct(i).Height),i) = OutStruct(i).WidthBase;
                OutArrayStruct.AspectRatioHalfHeight(1:length(OutStruct(i).Height),i) = OutStruct(i).AspectRatioHalfHeight;
                OutArrayStruct.AspectRatioBaseHeight(1:length(OutStruct(i).Height),i) = OutStruct(i).AspectRatioBaseHeight;
                OutArrayStruct.AreaDerivedDiameter(1:length(OutStruct(i).Height),i) = OutStruct(i).AreaDerivedDiameter;
            end
            
            
            if RecordMovieBool
                v.close
            end
            
            OutStructAll = Struct;
            close(h)
        end
        
        function sort_segments_by_name_and_subsegmentname(obj)
            
            if isempty(obj.Segment(1).ROIObject)
                return
            end
            
            for i =1:length(obj.Segment)
                NameList{i} = [obj.Segment(i).Name obj.Segment(i).SubSegmentName ];
            end
            [~,NewOrder] = sort(NameList);
            
            NewSegment = obj.Segment;
            for i=1:length(obj.Segment)
                NewSegment(i) = obj.Segment(NewOrder(i));
            end
            
            obj.Segment = NewSegment;
        end
        
    end
    methods (Static)
        % Static main methods
        
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
        
        function OutChannel = resize_channel(InChannel,Multiplicator,TargetRes,TransformToSquare)
            
            if nargin < 4
                TransformToSquare = false;
            end
            if nargin == 3
                Multiplicator = TargetRes/InChannel.NumPixelsX;
            end
            OutChannel = InChannel;
            
            if ~TransformToSquare
                OutChannel.Image = imresize(InChannel.Image,TargetRes/InChannel.NumPixelsX,'bilinear');
                OutChannel.NumPixelsX = size(OutChannel.Image,1);
                OutChannel.NumPixelsY = size(OutChannel.Image,2);
            elseif TransformToSquare
                OutChannel.Image = imresize(InChannel.Image,[TargetRes TargetRes],'bilinear');
                OutChannel.NumPixelsX = size(OutChannel.Image,1);
                OutChannel.NumPixelsY = size(OutChannel.Image,2);
            end
        end
        
        function [OutChannel,XMultiplier,YMultiplier] = resize_channel_to_same_size_per_pixel(InChannel)
            % OutChannel = resize_channel_to_same_size_per_pixel(InChannel)
            % 
            % Resizes the image dimension with bigger size per pixel so
            % the OutChannel always has more pixels.
            % New Image will have aspect ratio according to the ratio of
            % ScanSizeX and ScanSizeY
            
            SizePerPixelX = InChannel.ScanSizeX./InChannel.NumPixelsX;
            SizePerPixelY = InChannel.ScanSizeY./InChannel.NumPixelsY;
            
            OutChannel = InChannel;
            
            if SizePerPixelX == SizePerPixelY
                XMultiplier = 1;
                YMultiplier = 1;
                return
            elseif SizePerPixelX > SizePerPixelY
                NewNumPixels = round(InChannel.ScanSizeX*InChannel.NumPixelsY/InChannel.ScanSizeY);
                OutChannel.Image = imresize(InChannel.Image,[InChannel.NumPixelsY NewNumPixels]);
                OutChannel.NumPixelsX = NewNumPixels;
                XMultiplier = OutChannel.NumPixelsX/InChannel.NumPixelsX;
                YMultiplier = 1;
            elseif SizePerPixelY > SizePerPixelX
                NewNumPixels = round(InChannel.ScanSizeY*InChannel.NumPixelsX/InChannel.ScanSizeX);
                OutChannel.Image = imresize(InChannel.Image,[NewNumPixels InChannel.NumPixelsX]);
                OutChannel.NumPixelsY = NewNumPixels;
                YMultiplier = OutChannel.NumPixelsY/InChannel.NumPixelsY;
                XMultiplier = 1;
            end
        end
        
        function OutChannel = resize_channel_to_padded_same_size_per_pixel_square_image(InChannel,varargin)
            % function OutChannel = resize_channel_to_padded_same_size_per_pixel_square_image(InChannel,varargin)
            %
            % <FUNCTION DESCRIPTION HERE>
            %
            %
            % Required inputs
            % InChannel ... <VARIABLE DESCRIPTION>
            %
            % Name-Value pairs
            % "PaddingType" ... <NAMEVALUE DESCRIPTION>
            % "TargetResolution" ... <NAMEVALUE DESCRIPTION>
            
            p = inputParser;
            p.FunctionName = "resize_channel_to_padded_same_size_per_pixel_square_image";
            p.CaseSensitive = false;
            p.PartialMatching = true;
            
            % Required inputs
            validInChannel = @(x)true;
            addRequired(p,"InChannel",validInChannel);
            
            % NameValue inputs
            defaultPaddingType = 'Min';
            defaultTargetResolution = max(InChannel.NumPixelsX,InChannel.NumPixelsY);
            validPaddingType = @(x)any(validatestring(x,{'Min','Max','Zero'}));
            validTargetResolution = @(x)true;
            addParameter(p,"PaddingType",defaultPaddingType,validPaddingType);
            addParameter(p,"TargetResolution",defaultTargetResolution,validTargetResolution);
            
            parse(p,InChannel,varargin{:});
            
            % Assign parsing results to named variables
            InChannel = p.Results.InChannel;
            PaddingType = p.Results.PaddingType;
            TargetResolution = p.Results.TargetResolution;
            
            % First, equalize the size-per-pixel relation
            OutChannel = AFMBaseClass.resize_channel_to_same_size_per_pixel(InChannel);
            
            if isequal(lower(PaddingType),'min')
                PaddingValue = min(OutChannel.Image,[],'all');
            elseif isequal(lower(PaddingType),'max')
                PaddingValue = max(OutChannel.Image,[],'all');
            elseif isequal(lower(PaddingType),'zero')
                PaddingValue = 0;
            end
            
            %Padd the side wih less pixels
            if OutChannel.NumPixelsX < OutChannel.NumPixelsY
                PixelDiff = OutChannel.NumPixelsY - OutChannel.NumPixelsX;
                OutChannel.Image(:,end+1:end+PixelDiff) = PaddingValue;
                OutChannel.NumPixelsX = OutChannel.NumPixelsY;
                OutChannel.ScanSizeX = OutChannel.ScanSizeY;
            else
                PixelDiff = OutChannel.NumPixelsX - OutChannel.NumPixelsY;
                OutChannel.Image(end+1:end+PixelDiff,:) = PaddingValue;
                OutChannel.NumPixelsY = OutChannel.NumPixelsX;
                OutChannel.ScanSizeY = OutChannel.ScanSizeX;
            end
            
            OutChannel = AFMBaseClass.resize_channel(OutChannel,[],TargetResolution,true);
            
        end
        
    end
    methods (Static)
        % Static auxiliary methods
        
        function [LocalDirectionVector,WeightingVector] = find_local_direction_vector_in_ordered_vector_list(VectorList,SmoothingWindowSize,SmoothingWindowWeighting)
            
            if ~mod(SmoothingWindowSize,2)
                SmoothingWindowSize = SmoothingWindowSize + 1;
            end
            CenterIndex = SmoothingWindowSize/2 + .5;
            WeightingVector = zeros(SmoothingWindowSize,1);
            
            if isequal(lower(SmoothingWindowWeighting),'flat')
                WeightingVector = ones(SmoothingWindowSize,1).*1/SmoothingWindowSize;
            end
            if isequal(lower(SmoothingWindowWeighting),'linear')
                Lin = linspace(0,1,CenterIndex)';
                WeightingVector(1:CenterIndex) = Lin;
                WeightingVector(CenterIndex:end) = Lin(end:-1:1);
                WeightingVector = WeightingVector/sum(WeightingVector);
            end
            if isequal(lower(SmoothingWindowWeighting),'gaussian')
                Gauss = makedist('Normal','mu',CenterIndex,'sigma',sqrt(SmoothingWindowSize));
                WeightingVector = pdf(Gauss,linspace(1,SmoothingWindowSize,SmoothingWindowSize)');
                WeightingVector = WeightingVector/sum(WeightingVector);
            end
            
            LocalDirectionVector = zeros(size(VectorList));
            DiffVectorList = zeros(size(VectorList));
            
            TempDiff = diff(VectorList,1);
            DiffVectorList(2:end,:) = TempDiff;
            DiffVectorList(1,:) = TempDiff(1,:);
            
            for i=1:size(DiffVectorList,1)
                LowerIndex = max(i - (CenterIndex-1),1);
                UpperIndex = min(i + (CenterIndex-1),size(DiffVectorList,1));
                if UpperIndex - LowerIndex < length(WeightingVector)
                    NormalizingFactor = sum(WeightingVector(CenterIndex - (i-LowerIndex):CenterIndex + (UpperIndex-i)));
                else
                    NormalizingFactor = 1;
                end
                SumVector = sum(DiffVectorList(LowerIndex:UpperIndex,:).*WeightingVector(CenterIndex - (i-LowerIndex):CenterIndex + (UpperIndex-i)),1)./NormalizingFactor;
                LocalDirectionVector(i,:) = SumVector/norm(SumVector);
            end
            
        end
        
        function [Seg1,Seg2] = proximity_map_two_segmentation_elements(Seg1,Seg2)
            
            Seg1.ProximityMap = [];
            Seg2.ProximityMap = [];
            
            NumPointsSeg1 = length(Seg1.ROIObject.Position(:,1));
            NumPointsSeg2 = length(Seg2.ROIObject.Position(:,1));
            IndexListSeg1 = (1:NumPointsSeg1)';
            IndexListSeg2 = (1:NumPointsSeg2)';
            
            NumMappings = min(NumPointsSeg1,NumPointsSeg2);
            IndexPairs = zeros(NumMappings,2);
            
            SearchedPoints1 = Seg1.ROIObject.Position;
            SearchedPoints2 = Seg2.ROIObject.Position;
            
            DistanceMatrix = zeros(NumPointsSeg1,NumPointsSeg2);
            for i=1:NumPointsSeg1
                for j=1:NumPointsSeg2
                    DistanceMatrix(i,j) = norm(SearchedPoints1(i,:) - SearchedPoints2(j,:));
                end
            end
                
            k = 1;
            while ~isempty(DistanceMatrix)
                [~,LinIdx] = min(DistanceMatrix,[],'all','linear');
                [row,col] = ind2sub(size(DistanceMatrix),LinIdx);
                IndexPairs(k,:) = [IndexListSeg1(row) ; IndexListSeg2(col)];
%                 % Debug
%                 plot(SearchedPoints1(IndexListSeg1,1),SearchedPoints1(IndexListSeg1,2),'bO',...
%                     SearchedPoints2(IndexListSeg2,1),SearchedPoints2(IndexListSeg2,2),'rO',...
%                     SearchedPoints1(row,1),SearchedPoints1(row,2),'bX',SearchedPoints2(col,1),SearchedPoints2(col,2),'rX')
%                 drawnow
                IndexListSeg1(row) = [];
                IndexListSeg2(col) = [];
                DistanceMatrix(row,:) = [];
                if isempty(DistanceMatrix)
                    break
                end
                DistanceMatrix(:,col) = [];
                k= k + 1;
            end
            
            Seg1.ProximityMap.IndexPairs = IndexPairs;
            Seg2.ProximityMap.IndexPairs = IndexPairs;
            
        end
        
        function [OutSegment,Index] = find_matching_segment(SegmentElement,SegmentStruct)
            
            TargetName = SegmentElement.Name;
            TargetSubName = SegmentElement.SubSegmentName;
            
            NameList = {SegmentStruct.Name};
            SubNameList = {SegmentStruct.SubSegmentName};
            
            Index = find(strcmp({TargetName},NameList) & strcmp(TargetSubName,SubNameList));
            
            if isempty(Index)
                Index = [];
                OutSegment = [];
                return
            end
            
            OutSegment = SegmentStruct(Index);
        end
        
    end
    methods
        % auxiliary methods
        
        function IndexVector = get_indizes_of_matching_segments(obj,FullNameCell)
            % IndexVector = get_indizes_of_matching_segments(obj,FullNameCell)
            
            IndexVector = 0;
            
            for i=1:length(obj.Segment)
                TempName = [obj.Segment(i).Name ' || ' obj.Segment(i).SubSegmentName];
                if any(matches(TempName,FullNameCell))
                    IndexVector(end+1) = i;
                end
            end
            
        end
        
        function FullNameCell = get_full_names_from_segment_indizes(obj,IndexVector)
            % FullNameCell = get_full_names_from_segment_indizes(obj,IndexVector)
            
            IndexVector(IndexVector == 0) = [];
            
            FullNameCell = {};
            
            k = 1;
            for i=IndexVector
                FullNameCell{k} = [obj.Segment(i).Name ' || ' obj.Segment(i).SubSegmentName];
                k=k+1;
            end
            
        end
        
    end
end