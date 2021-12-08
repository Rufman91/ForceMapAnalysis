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
                obj.Folder = [DataFolder regexprep(obj.Name,'[.]','')];
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
                        NewVertices = [NewVertices; TempNewVertices];
                        OriginalVertices(1,:) = [];
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
                [LocalDirectionVector,WeightingVector] = AFMBaseClass.find_local_direction_vector_in_ordered_vector_list(Snapped(i).ROIObject.Position,4*SmoothingWindowSize,'flat');
                for j=1:size(Snapped(i).ROIObject.Position,1)
                    OriginalPos(j,:) = Snapped(i).ROIObject.Position(j,:);
                    PerpendicularVector(j,:) = [LocalDirectionVector(j,2) -LocalDirectionVector(j,1)]/norm(LocalDirectionVector(j,:));
                    WindowStart = Snapped(i).ROIObject.Position(j,:) + PerpendicularVector(j,:).*WidthLocalWindowPixels/2;
                    WindowEnd = Snapped(i).ROIObject.Position(j,:) - PerpendicularVector(j,:).*WidthLocalWindowPixels/2;
                    %Debug
%                     scatter([WindowStart(1) WindowEnd(1)],[WindowStart(2) WindowEnd(2)])
%                     xlim([0 512])
%                     ylim([0 512])
%                     drawnow
                    [LocalX,LocalY,LocalProfile] = improfile(Channel.Image,[WindowStart(1) WindowEnd(1)],[WindowStart(2) WindowEnd(2)]);
                    if length(LocalProfile) >= 4 && ~sum(isnan(LocalProfile))
                        LocalX = interp1(LocalX,linspace(0,1,100)'.*length(LocalX),'spline');
                        LocalY = interp1(LocalY,linspace(0,1,100)'.*length(LocalY),'spline');
                        LocalProfile = interp1(LocalProfile,linspace(0,1,100)'.*length(LocalProfile),'spline');
                    end
                    [~,MaxIndex] = max(LocalProfile);
                    SnappedPos(j,:) = [LocalX(MaxIndex) LocalY(MaxIndex)];
                    DisplacementVector(j,:) = (SnappedPos(j,:) - OriginalPos(j,:));
                    SnappedToOriginalDistance(j) = PerpendicularVector(j,:)*DisplacementVector(j,:)';
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
                        CenteredLP = filloutliers(CenteredLP,'linear',1);
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
        
        function [OutMat,OutConcatCell,OutCell,OutMask] = get_segment_data_from_channel(obj,ChannelName,PixelDilation)
            
            if isempty(obj.Segment) || isempty(obj.Segment(1).Name)
                warning(sprintf('%s has no segementation data',obj.Name))
            end
            if nargin < 3
                PixelDilation = 0;
            end
            
            Channel = obj.get_channel(ChannelName);
            Data = Channel.Image;
            OutMask = zeros(size(Data));
            
            SegmentNames = unique({obj.Segment.Name});
            set(0,'DefaultFigureVisible','off');
            
            for i=1:length(SegmentNames)
                k = 1;
                for j=1:length(obj.Segment)
                    if isequal(SegmentNames{i},obj.Segment(j).Name)
                        imshow(Data)
                        Line = drawpolyline('Position',obj.Segment(j).ROIObject.Position);
                        Mask = Line.createMask;
                        if PixelDilation
                            StrEl = strel('Disk',PixelDilation,0);
                            Mask = imdilate(Mask,StrEl);
                        end
                        OutMask = OutMask | Mask;
                        SubSegmentPoints = Data(Mask == 1);
                        OutCell{i}{k} = SubSegmentPoints;
                        k = k + 1;
                    end
                end
            end
            close gcf
            set(0,'DefaultFigureVisible','on');
            
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
        
        function automatic_segmentation_on_singular_vertical_fiber(obj,SampleDistanceMeters,WidthLocalWindowMeters,SmoothingWindowSize,SmoothingWindowWeighting,Indizes)
            
            
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
        
        function characterize_fiber_like_polyline_segments(obj,SampleDistanceMeters,WidthLocalWindowMeters,SmoothingWindowSize,SmoothingWindowWeighting)
            
            
            
        end
        
    end
    methods (Static)
        % Static main methods
        
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
    
end