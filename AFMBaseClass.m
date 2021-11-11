classdef AFMBaseClass < matlab.mixin.Copyable & matlab.mixin.SetGet & handle
    % This a baseclass for the classes AFMImage and ForceMap to inherit
    % shared methods and properties from
    
    properties
        Name
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
                            'ROIObject',[])
        OverlayGroup
    end
    
    methods
        % Main Methods
        
        function obj = AFMBaseClass()
            
            obj.OverlayGroup.hasOverlayGroup = false;
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
        
        function snap_line_segments_to_local_perpendicular_maximum(obj,SampleDistanceMeters,WidthLocalWindowMeters,SmoothingWindowSize,SmoothingWindowWeighting)
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
            %   window weighting; possibler options:
            %         'flat','linear','gaussian'
            
            if nargin < 5
                SmoothingWindowWeighting = 'flat';
            end
            if nargin < 4
                SmoothingWindowSize = 1;
                SmoothingWindowWeighting = 'flat';
            end
            if nargin < 3
                SampleDistanceMeters = 50*1e-9;
                SmoothingWindowSize = 1;
                SmoothingWindowWeighting = 'flat';
            end
            if nargin < 2
                SampleDistanceMeters = 50*1e-9;
                WidthLocalWindowMeters = 300*1e-9;
                SmoothingWindowSize = 1;
                SmoothingWindowWeighting = 'flat';
            end
            
            for i=length(obj.Segment):-1:1
                if ~isempty(strfind(obj.Segment(i).Name,'Snapped'))
                    obj.Segment(i) = [];
                end
            end
            
            Channel = obj.get_channel('Processed');
            if isempty(Channel)
                warning('No Channel "Processed" found. Need flattened height channel for snap to local maximum')
                return
            end
            
            % If Channel has an unequal SizePerPixel, resize it for
            % following operations
            [Channel,XMult,YMult] = obj.resize_channel_to_same_size_per_pixel(Channel);
            
            % Convert inputs from meters to pixels
            SizePerPixel = Channel.ScanSizeX/Channel.NumPixelsX;
            WidthLocalWindowPixels = WidthLocalWindowMeters/SizePerPixel;
            SampleDistancePixels = SampleDistanceMeters/SizePerPixel;
            
            % Loop over all Segments and create 'Snapped' segments for all
            % Polylines
            k = 1;
            for i=1:length(obj.Segment)
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
                SnappedPos = [];
                SmoothedSnappedPos = [];
                SnappedToOriginalDistance = [];
                PerpendicularVector = [];
                OriginalPos = [];
                FinalSnappedPos = [];
                [LocalDirectionVector,WeightingVector] = AFMBaseClass.find_local_direction_vector_in_ordered_vector_list(Snapped(i).ROIObject.Position,SmoothingWindowSize,SmoothingWindowWeighting);
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
                % Apply smoothing to SnappedPos
                SnappedToOriginalDistance = SnappedToOriginalDistance';
                SmoothedSnappedDistance = zeros(size(SnappedToOriginalDistance));
                CenterIndex = floor(SmoothingWindowSize/2+1);
                for j=1:size(SnappedToOriginalDistance,1)
                    LowerIndex = max(j - (CenterIndex-1),1);
                    UpperIndex = min(j + (CenterIndex-1),size(SnappedToOriginalDistance,1));
                    SmoothedSnappedDistance(j) = sum(SnappedToOriginalDistance(LowerIndex:UpperIndex).*WeightingVector(CenterIndex - (j-LowerIndex):CenterIndex + (UpperIndex-j)),1);
                    FinalSnappedPos(j,:) = OriginalPos(j,:) + PerpendicularVector(j,:).*SmoothedSnappedDistance(j);
                end
% 
%                 SmoothedSnappedPos = zeros(size(SnappedPos));
%                 CenterIndex = floor(SmoothingWindowSize/2+1);
%                 for j=1:size(SnappedPos,1)
%                     LowerIndex = max(j - (CenterIndex-1),1);
%                     UpperIndex = min(j + (CenterIndex-1),size(SnappedPos,1));
%                     if UpperIndex - LowerIndex < length(WeightingVector)
%                         NormalizingFactor = sum(WeightingVector(CenterIndex - (j-LowerIndex):CenterIndex + (UpperIndex-j)));
%                     else
%                         NormalizingFactor = 1;
%                     end
%                     SegmentMean = mean(SnappedPos(LowerIndex:UpperIndex,:),1);
%                     SmoothedSnappedPos(j,:) = sum((SnappedPos(LowerIndex:UpperIndex,:) - SegmentMean).*WeightingVector(CenterIndex - (j-LowerIndex):CenterIndex + (UpperIndex-j))./NormalizingFactor,1);
%                     FinalSnappedPos(j,:) = SmoothedSnappedPos(j,:) + SegmentMean;
%                 end
% 
%                 SmoothedSnappedPos = zeros(size(SnappedPos));
%                 
%                 SegmentMean = mean(SnappedPos,1);
%                     
%                 SmoothedSnappedPos(:,1) = conv(SnappedPos(:,1) - SegmentMean(1),WeightingVector,'same');
%                 SmoothedSnappedPos(:,2) = conv(SnappedPos(:,2) - SegmentMean(2),WeightingVector,'same');
%                 FinalSnappedPos = SmoothedSnappedPos + SegmentMean;
                
                FinalSnappedPos(:,1) = FinalSnappedPos(:,1)./YMult;
                FinalSnappedPos(:,2) = FinalSnappedPos(:,2)./XMult;
                Snapped(i).ROIObject.Position = FinalSnappedPos;
            end
            
            obj.Segment(end+1:end+length(Snapped)) = Snapped;
            
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
                OutChannel.Image = imresize(InChannel.Image,Multiplicator);
                OutChannel.NumPixelsX = size(OutChannel.Image,1);
                OutChannel.NumPixelsY = size(OutChannel.Image,2);
            elseif TransformToSquare
                OutChannel.Image = imresize(InChannel.Image,[TargetRes TargetRes]);
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
        
    end
    
end