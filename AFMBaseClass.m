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
        
    end
    methods (Static)
        % Static auxiliary methods
        
        
    end
    
end