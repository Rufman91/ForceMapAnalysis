classdef AFMBaseClass < matlab.mixin.Copyable & matlab.mixin.SetGet & handle & dynamicprops
    % This a baseclass for the classes AFMImage and ForceMap to inherit
    % shared methods and properties from
    
    properties
        Name = ''
        Folder = ''
        ID = ''
        HostOS = ''         % Operating System
        HostName = ''       % Name of hosting system
        FileType = 'Image'
        ScanSizeX = []          % Size of imaged window in X-direction
        ScanSizeY = []           % Size of imaged window in Y-direction
        ScanAngle = 0   % in degrees (°)
        NumPixelsX = []
        NumPixelsY = []
        OriginX = 0
        OriginY = 0
        List2Map = []        % An R->RxR ((k)->(i,j)) mapping of indices to switch between the two representations
        Map2List = []      % An RxR->R ((i,j)->(k))mapping of indices to switch between the two representations
        Metadata = []
        CurrentFMA_ID = 'none'
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
        SegmentName = ''
        SubSegmentName = ''
        SubSegmentFullName = ''
        FiberSegment_Area = []
        FiberSegment_AreaDerivedDiameter = []
        FiberSegment_Height = []
        FiberSegment_WidthHalfHeight = []
        FiberSegment_Prominence = []
        FiberSegment_WidthHalfProminence = []
        FiberSegment_WidthBase = []
        FiberSegment_AspectRatioHalfHeight = []
        FiberSegment_AspectRatioBaseHeight = []
        FiberSegment_Mean_Area = []
        FiberSegment_Mean_AreaDerivedDiameter = []
        FiberSegment_Mean_Height = []
        FiberSegment_Mean_WidthHalfHeight = []
        FiberSegment_Mean_Prominence = []
        FiberSegment_Mean_WidthHalfProminence = []
        FiberSegment_Mean_WidthBase = []
        FiberSegment_Mean_AspectRatioHalfHeight = []
        FiberSegment_Mean_AspectRatioBaseHeight = []
        FiberSegment_Median_Area = []
        FiberSegment_Median_AreaDerivedDiameter = []
        FiberSegment_Median_Height = []
        FiberSegment_Median_WidthHalfHeight = []
        FiberSegment_Median_Prominence = []
        FiberSegment_Median_WidthHalfProminence = []
        FiberSegment_Median_WidthBase = []
        FiberSegment_Median_AspectRatioHalfHeight = []
        FiberSegment_Median_AspectRatioBaseHeight = []
        FiberSegment_RelativePixelPosition = []
        FiberSegment_RelativePosition = []
        FiberSegment_SegmentLength = []
        FiberSegment_Ellipse_a = []
        FiberSegment_Ellipse_b = []
        FiberSegment_Ellipse_AspectRatio = []
        FiberSegment_Ellipse_Area = []
        FiberSegment_Ellipse_Height = []
        FiberSegment_Ellipse_WidthHalfHeight = []
        FiberSegment_Mean_Ellipse_a = []
        FiberSegment_Mean_Ellipse_b = []
        FiberSegment_Mean_Ellipse_AspectRatio = []
        FiberSegment_Mean_Ellipse_Area = []
        FiberSegment_Mean_Ellipse_Height = []
        FiberSegment_Mean_Ellipse_WidthHalfHeight = []
        FiberSegment_Median_Ellipse_a = []
        FiberSegment_Median_Ellipse_b = []
        FiberSegment_Median_Ellipse_AspectRatio = []
        FiberSegment_Median_Ellipse_Area = []
        FiberSegment_Median_Ellipse_Height = []
        FiberSegment_Median_Ellipse_WidthHalfHeight = []
        OverlayGroup = ''
        OverlayGroupName = ''
        OverlayGroupIndex = []
        hasOverlayGroup = false
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
                        set(obj,PropertyNames{i},MetaProperties(i).DefaultValue);
                    else
                        set(obj,PropertyNames{i},[]);
                    end
                end
            end
            
            cd(current.path)
            
        end
        
        function deconvolute_image(obj,CTClassInstance,ChannelName,MaxResolution,KeepOldResults)
            % Deconvolutes the specified image channel by removing the convolution effects of the AFM tip.
            %
            % This function applies mathematical morphology to deconvolute the image of the AFM tip from
            % the sample image, aiming to enhance the image resolution and quality. It adjusts the image based
            % on the tip's known shape and size, potentially improving the accuracy of subsequent analyses.
            %
            % Parameters:
            %   obj (AFMBaseClass): The instance of the class containing the method, which should include
            %                       image data and metadata.
            %   CTClassInstance (CantileverTipClass): An instance of the CantileverTipClass containing information
            %                                         about the AFM cantilever tip, including its shape and size.
            %   ChannelName (string): The name of the image channel to be deconvoluted. Default is 'Processed'.
            %   MaxResolution (integer): The maximum resolution for the deconvolution process, specified as the
            %                            length of one side of the image in pixels. This parameter controls the
            %                            trade-off between processing time and deconvolution quality. Higher values
            %                            lead to finer deconvolution at the cost of increased processing time.
            %                            Default is 1024.
            %   KeepOldResults (logical): A flag indicating whether to keep the original, unprocessed image data
            %                             alongside the deconvoluted results. Setting this to true allows for
            %                             comparison between original and processed images but requires additional
            %                             memory. Default is true.
            %
            % Returns:
            %   This function does not return a value. Instead, it modifies the instance of AFMBaseClass in place,
            %   updating the specified image channel with the deconvoluted image.
            %
            % Examples:
            %   afmImageInstance.deconvolute_image(cantileverTipInstance, 'Height', 1024, false);
            %   This example deconvolutes the 'Height' channel of the afmImageInstance using the specified
            %   cantileverTipInstance. The maximum resolution is set to 1024 pixels, and the original image data
            %   is not retained.
            %
            % See also:
            %   CantileverTipClass, AFMBaseClass
            
            if nargin < 3
                ChannelName = 'Processed';
                MaxResolution = 1024;
                KeepOldResults = true;
            end
            if nargin < 4
                MaxResolution = 1024;
                KeepOldResults = true;
            end
            if nargin < 5
                KeepOldResults = true;
            end
            
            ErodedTip = CTClassInstance.get_channel('Eroded Tip');
            if isempty(ErodedTip)
                CTClassInstance.deconvolute_cantilever_tip;
                ErodedTip = CTClassInstance.get_channel('Eroded Tip');
            end
            Processed = obj.get_channel(ChannelName);
            if isempty(Processed)
                warning(sprintf('Channel "%s" not found. Cannot deconvolute!',ChannelName));
                return
            end
            
            PPHasUnequalSPP = false;
            [PP,XMult,YMult] = AFMBaseClass.resize_channel_to_padded_same_size_per_pixel_square_image(Processed);
            ET = AFMBaseClass.resize_channel_to_padded_same_size_per_pixel_square_image(ErodedTip);
            
            if XMult ~= 1 || YMult ~= 1
                PPHasUnequalSPP = true;
            end
            
            imshowpair(Processed.Image,PP.Image,'montage')
            % Resize and pad images to same resolution at correct spacial
            % ratio
            Multiplier = (ET.NumPixelsY/ET.ScanSizeX)/(PP.NumPixelsY/PP.ScanSizeX);
            
            if Multiplier >= 1
                PP.Image = imresize(PP.Image,[round(Multiplier*PP.NumPixelsX) nan]);
                PPResized = true;
                ETResized = false;
            else
                ET.Image = imresize(ET.Image,[round(1/Multiplier*ET.NumPixelsX) nan]);
                PPResized = false;
                ETResized = true;
            end
            
            PaddingSize = abs(size(PP.Image) - size(ET.Image));
            PaddingSizePre = floor(PaddingSize./2);
            PaddingSizePost = ceil(PaddingSize./2);
            if size(PP.Image,1) >= size(ET.Image,1)
                ET.Image = padarray(ET.Image,PaddingSizePre,...
                    min(ET.Image,[],'all'),'pre');
                ET.Image = padarray(ET.Image,PaddingSizePost,...
                    min(ET.Image,[],'all'),'post');
                ETPadded = true;
                PPPadded = false;
            else
                PP.Image = padarray(PP.Image,PaddingSizePre,...
                    min(PP.Image,[],'all'),'pre');
                PP.Image = padarray(PP.Image,PaddingSizePost,...
                    min(PP.Image,[],'all'),'post');
                ETPadded = false;
                PPPadded = true;
            end
            
            % if images ended up being too big, scale them down for
            % deconvolution
            isDownscaled = false; 
            if size(ET.Image,1) > MaxResolution
                DownscaleFactor = MaxResolution/size(ET.Image,1);
                ET.Image = imresize(ET.Image,[MaxResolution nan]);
                PP.Image = imresize(PP.Image,[MaxResolution nan]);
                isDownscaled = true;
            end
%             subplot(2,1,1)
%             imshowpair(PP.Image,ET.Image,'montage')
%             subplot(2,1,2)
%             imshowpair(Processed.Image,PP.Image,'montage')
%             drawnow
            
            % Make sure the Tip is below the convoluted image EVERYWHERE
            Min = min(PP.Image,[],'all');
            PP.Image  = PP.Image - Min;

            Out = obj.deconvolute_by_mathematical_morphology(PP.Image,ET.Image);
            
            % put it back
            Out = Out + Min;
            
            if isDownscaled
                Out = imresize(Out,[round(MaxResolution/DownscaleFactor) nan]);
            end
            if PPPadded
                Out = Out((PaddingSizePre(1) + 1):(end - PaddingSizePost(1)),...
                    (PaddingSizePre(1) + 1):(end - PaddingSizePost(2)));
            end
            if PPResized || PPHasUnequalSPP
                Out = imresize(Out,[round(PP.NumPixelsX*YMult) round(PP.NumPixelsY*XMult)]);
            end
            
            OutImage = Out(1:Processed.NumPixelsX,1:Processed.NumPixelsY);
            
            Deconvoluted = Processed;
            Deconvoluted.Name = 'Deconvoluted';
            
            MinIn = min(Processed.Image,[],'all');
            MinOut = min(OutImage,[],'all');
            
            OutImage = OutImage + (MinIn - MinOut);
            
            Deconvoluted.Image = OutImage;
            
            obj.add_channel(Deconvoluted,~KeepOldResults)
        end
        
        function clear_all_properties(obj)
            
            PropertyStruct = obj.get();
            PropertyNames = fieldnames(PropertyStruct);
            
            for i=1:length(PropertyNames)
                set(obj,PropertyNames{i},[])
            end
        end
        
        function construct_list_to_map_relations(obj)
            obj.Map2List = zeros(obj.NumPixelsX, obj.NumPixelsY);
            obj.List2Map = zeros(obj.NumPixelsX * obj.NumPixelsY, 2);
            
            [i_mat, j_mat] = meshgrid(1:obj.NumPixelsY, 1:obj.NumPixelsX);
            i_mat = i_mat';
            j_mat = j_mat';
            
            k = 1;
            k_mat = reshape(1:(obj.NumPixelsX * obj.NumPixelsY), obj.NumPixelsY, obj.NumPixelsX);
            
            if isequal(obj.FileType, 'quantitative-imaging-map')
                obj.Map2List = k_mat';
                obj.List2Map = [j_mat(:), i_mat(:)];
            elseif isequal(obj.FileType, 'force-scan-map')
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
                obj.Map2List = k_mat';
                obj.List2Map = [j_mat(:), i_mat(:)];
            end
        end
                
        function construct_list_to_map_relations_old(obj)
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
        
        function List = convert_map_to_data_list(obj,Map)
            
            NPixels = obj.NumPixelsX*obj.NumPixelsY;
            List = zeros(NPixels,1);
            
            for i=1:NPixels
                List(i) = Map(obj.List2Map(i,1),obj.List2Map(i,2));
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
        
        function [ChannelNameResult, AllMatches] = search_channel(obj, ChannelName)
            k = 0;
            ChannelNameResult = '';
            AllMatches = {};
            for i = 1:length(obj.Channel)
                if contains(obj.Channel(i).Name, ChannelName, 'IgnoreCase', true)
                    if k == 0
                        ChannelNameResult = obj.Channel(i).Name;
                    end
                    k = k + 1;
                    AllMatches{end+1} = obj.Channel(i).Name; %#ok<AGROW>
                end
            end
            if k == 0
                warning('No channel name matches the searched one: %s', ChannelName);
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
        
        function add_channel(obj,Channel,ReplaceSameNamed)
            
            if nargin < 3
                ReplaceSameNamed = false;
            end
            
            if ~isfield(obj.Channel,'FMA_ID')
                for i=1:length(obj.Channel)
                    obj.Channel(i).FMA_ID = 'none';
                end
            end
            Channel.FMA_ID = obj.CurrentFMA_ID;
            
            % Write to Channel
            [~,Index] = obj.get_channel(Channel.Name);
            if isempty(Index)
                obj.Channel(end+1) = Channel;
            else
                if ~ReplaceSameNamed
                    if ~isempty(regexp(Channel.Name(end-4:end),'.\(\d\d\)','once'))
                        OriginalName = Channel.Name(1:end-5);
                        Channel.Name = sprintf('%s (%s)',OriginalName,sprintf('%02i',str2num(Channel.Name(end-2:end-1))+1));
                    else
                        Channel.Name = sprintf('%s (%s)',Channel.Name,sprintf('%02i',1));
                    end
                    obj.add_channel(Channel,false);
                else
                    obj.Channel(Index) = Channel;
                end
            end
            
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
            
            obj.add_channel(OutChannel,true)
            
            obj.assert_channel_number;
        end
        
        function assert_channel_number(obj)
            
            obj.NumChannels = numel(obj.Channel);
            
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
            
            SizePerPixelX1 = DonorClass.Channel(1).ScanSizeX/DonorClass.Channel(1).NumPixelsX;
            SizePerPixelX2 = AcceptorClass.Channel(1).ScanSizeX/AcceptorClass.Channel(1).NumPixelsX;
            if SizePerPixelX1 ~= SizePerPixelX2
                ScaleMultiplierX = SizePerPixelX1/SizePerPixelX2;
            else
                ScaleMultiplierX = 1;
            end
            
            SizePerPixelY1 = DonorClass.Channel(1).ScanSizeY/DonorClass.Channel(1).NumPixelsY;
            SizePerPixelY2 = AcceptorClass.Channel(1).ScanSizeY/AcceptorClass.Channel(1).NumPixelsY;
            if SizePerPixelY1 ~= SizePerPixelY2
                ScaleMultiplierY = SizePerPixelY1/SizePerPixelY2;
            else
                ScaleMultiplierY = 1;
            end
            
            XDiff = DonorClass.Channel(1).OriginX - AcceptorClass.Channel(1).OriginX;
            SizePerPixelX = AcceptorClass.Channel(1).ScanSizeX./AcceptorClass.Channel(1).NumPixelsX;
            XDiff = XDiff/SizePerPixelX;
            YDiff = DonorClass.Channel(1).OriginY - AcceptorClass.Channel(1).OriginY;
            SizePerPixelY = AcceptorClass.Channel(1).ScanSizeY./AcceptorClass.Channel(1).NumPixelsY;
            YDiff = YDiff/SizePerPixelY;
            AngleDiff = DonorClass.Channel(1).ScanAngle - AcceptorClass.Channel(1).ScanAngle;
            AngleDiff = deg2rad(-AngleDiff);
            
            Vector = [ScaleMultiplierX*X ScaleMultiplierY*Y];
            
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
                        %%%%%%%
                        
                        CurrentDrawMode = lower(obj.Segment(j).Type);
                        switch CurrentDrawMode
                            case 'line'
                                DrawObject = drawline('Position',obj.Segment(j).ROIObject.Position,...
                                    'Deletable',1,...
                                    'InteractionsAllowed','all',...
                                    'LineWidth',obj.Segment(j).ROIObject.LineWidth,...
                                    'Label',sprintf('%s || %s',obj.Segment(j).Name,obj.Segment(j).SubSegmentName),...
                                    'LabelAlpha',0.6,...
                                    'Parent',Parent);
                            case 'freehand'
                                DrawObject = drawfreehand('Position',obj.Segment(j).ROIObject.Position,...
                                    'Deletable',1,...
                                    'InteractionsAllowed','all',...
                                    'LineWidth',obj.Segment(j).ROIObject.LineWidth,...
                                    'Label',sprintf('%s || %s',obj.Segment(j).Name,obj.Segment(j).SubSegmentName),...
                                    'LabelAlpha',0.6,...
                                    'Parent',Parent);
                            case 'square'
                                DrawObject = drawsquare('Position',obj.Segment(j).ROIObject.Position,...
                                    'Deletable',1,...
                                    'InteractionsAllowed','all',...
                                    'LineWidth',obj.Segment(j).ROIObject.LineWidth,...
                                    'Label',sprintf('%s || %s',obj.Segment(j).Name,obj.Segment(j).SubSegmentName),...
                                    'LabelAlpha',0.6,...
                                    'Parent',Parent);
                            case 'circle'
                                DrawObject = drawcircle('Center',obj.Segment(j).ROIObject.Center,...
                                    'Radius',obj.Segment(j).ROIObject.Radius,...
                                    'Deletable',1,...
                                    'InteractionsAllowed','all',...
                                    'LineWidth',obj.Segment(j).ROIObject.LineWidth,...
                                    'Label',sprintf('%s || %s',obj.Segment(j).Name,obj.Segment(j).SubSegmentName),...
                                    'LabelAlpha',0.6,...
                                    'Parent',Parent);
                            case 'ellipse'
                                DrawObject = drawellipse('Center',obj.Segment(j).ROIObject.Center,...
                                    'SemiAxes',obj.Segment(j).ROIObject.SemiAxes,...
                                    'RotationAngle',obj.Segment(j).ROIObject.RotationAngle,...
                                    'Deletable',1,...
                                    'InteractionsAllowed','all',...
                                    'LineWidth',obj.Segment(j).ROIObject.LineWidth,...
                                    'Label',sprintf('%s || %s',obj.Segment(j).Name,obj.Segment(j).SubSegmentName),...
                                    'LabelAlpha',0.6,...
                                    'Parent',Parent);
                            case 'polygon'
                                DrawObject = drawpolygon('Position',obj.Segment(j).ROIObject.Position,...
                                    'Deletable',1,...
                                    'InteractionsAllowed','all',...
                                    'LineWidth',obj.Segment(j).ROIObject.LineWidth,...
                                    'Label',sprintf('%s || %s',obj.Segment(j).Name,obj.Segment(j).SubSegmentName),...
                                    'LabelAlpha',0.6,...
                                    'Parent',Parent);
                            case 'rectangle'
                                DrawObject = drawrectangle('Position',obj.Segment(j).ROIObject.Position,...
                                    'RotationAngle',obj.Segment(j).ROIObject.RotationAngle,...
                                    'Rotatable',true,...
                                    'Deletable',1,...
                                    'InteractionsAllowed','all',...
                                    'LineWidth',obj.Segment(j).ROIObject.LineWidth,...
                                    'Label',sprintf('%s || %s',obj.Segment(j).Name,obj.Segment(j).SubSegmentName),...
                                    'LabelAlpha',0.6,...
                                    'Parent',Parent);
                            case 'crosshair'
                                DrawObject = drawcrosshair('Position',obj.Segment(j).ROIObject.Position,...
                                    'Deletable',1,...
                                    'InteractionsAllowed','all',...
                                    'LineWidth',obj.Segment(j).ROIObject.LineWidth,...
                                    'Label',sprintf('%s || %s',obj.Segment(j).Name,obj.Segment(j).SubSegmentName),...
                                    'LabelAlpha',0.6,...
                                    'Parent',Parent);
                            case 'point'
                                DrawObject = drawpoint('Position',obj.Segment(j).ROIObject.Position,...
                                    'Deletable',1,...
                                    'InteractionsAllowed','all',...
                                    'LineWidth',obj.Segment(j).ROIObject.LineWidth,...
                                    'Label',sprintf('%s || %s',obj.Segment(j).Name,obj.Segment(j).SubSegmentName),...
                                    'LabelAlpha',0.6,...
                                    'Parent',Parent);
                            case 'assisted'
                                DrawObject = drawassisted('Position',obj.Segment(j).ROIObject.Position,...
                                    'Deletable',1,...
                                    'InteractionsAllowed','all',...
                                    'LineWidth',obj.Segment(j).ROIObject.LineWidth,...
                                    'Label',sprintf('%s || %s',obj.Segment(j).Name,obj.Segment(j).SubSegmentName),...
                                    'LabelAlpha',0.6,...
                                    'Parent',Parent);
                            case 'polyline'
                                DrawObject = drawpolyline('Position',obj.Segment(j).ROIObject.Position,...
                                    'Deletable',1,...
                                    'InteractionsAllowed','all',...
                                    'LineWidth',obj.Segment(j).ROIObject.LineWidth,...
                                    'Label',sprintf('%s || %s',obj.Segment(j).Name,obj.Segment(j).SubSegmentName),...
                                    'LabelAlpha',0.6,...
                                    'Parent',Parent);
                        end
                        %%%%%%%
                        Mask = DrawObject.createMask;
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
                if isfield(obj.Segment(i),'ProximityMap') && ~isempty(obj.Segment(i).ProximityMap)
                    obj.Segment(i).ProximityMap(end+1) = Seg1.ProximityMap;
                else
                    obj.Segment(i).ProximityMap = Seg1.ProximityMap;
                end
                
                if isfield(TargetObject.Segment(i),'ProximityMap') && ~isempty(TargetObject.Segment(i).ProximityMap)
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
            
            [~,MaxIdx] = max(filloutliers(HeightMap,'linear','movmedian',max(3,ceil(HeightChannel.NumPixelsY/100)),2),[],2);
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
            % "EllipseFitThreshold" ... Determines the cutoff above
            %                       baseline in terms of ThresholdType
            %                       above which profile points are to be
            %                       taken into account for the ellipse fit.
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
            defaultHeightChannelName = 'Processed';
            defaultWidthLocalWindowMeters = 800e-9;
            defaultSmoothingWindowSize = 41;
            defaultMinPeakDistanceMeters = 50e-9;
            defaultLowerEndThreshold = .1;
            defaultEllipseFitThreshold = .75;
            defaultThresholdType = 'Fraction';
            defaultVerbose = false;
            defaultRecordMovieBool = false;
            defaultKeyFrames = 3;
            validHeightChannelName = @(x)ischar(x);
            validWidthLocalWindowMeters = @(x)isnumeric(x)&&isscalar(x);
            validSmoothingWindowSize = @(x)isnumeric(x)&&mod(x,1)==0;
            validMinPeakDistanceMeters = @(x)isnumeric(x)&&isscalar(x);
            validLowerEndThreshold = @(x)isnumeric(x)&&isscalar(x);
            validEllipseFitThreshold = @(x)isnumeric(x)&&isscalar(x); 
            validThresholdType = @(x)any(validatestring(x,{'Fraction','Meters'}));
            validVerbose = @(x)islogical(x);
            validRecordMovieBool = @(x)islogical(x);
            validKeyFrames = @(x)isnumeric(x)&&mod(x,1)==0;
            addParameter(p,"HeightChannelName",defaultHeightChannelName,validHeightChannelName);
            addParameter(p,"WidthLocalWindowMeters",defaultWidthLocalWindowMeters,validWidthLocalWindowMeters);
            addParameter(p,"SmoothingWindowSize",defaultSmoothingWindowSize,validSmoothingWindowSize);
            addParameter(p,"MinPeakDistanceMeters",defaultMinPeakDistanceMeters,validMinPeakDistanceMeters);
            addParameter(p,"LowerEndThreshold",defaultLowerEndThreshold,validLowerEndThreshold);
            addParameter(p,"EllipseFitThreshold",defaultEllipseFitThreshold,validEllipseFitThreshold);
            addParameter(p,"ThresholdType",defaultThresholdType,validThresholdType);
            addParameter(p,"Verbose",defaultVerbose,validVerbose);
            addParameter(p,"RecordMovieBool",defaultRecordMovieBool,validRecordMovieBool);
            addParameter(p,"KeyFrames",defaultKeyFrames,validKeyFrames);
            
            parse(p,obj,varargin{:});
            
            % Assign parsing results to named variables
            obj = p.Results.obj;
            HeightChannelName = p.Results.HeightChannelName;
            WidthLocalWindowMeters = p.Results.WidthLocalWindowMeters;
            SmoothingWindowSize = p.Results.SmoothingWindowSize;
            MinPeakDistanceMeters = p.Results.MinPeakDistanceMeters;
            LowerEndThreshold = p.Results.LowerEndThreshold;
            EllipseFitThreshold = p.Results.EllipseFitThreshold;
            ThresholdType = p.Results.ThresholdType;
            Verbose = p.Results.Verbose;
            RecordMovieBool = p.Results.RecordMovieBool;
            KeyFrames = p.Results.KeyFrames;
            
            if isequal(ThresholdType,'Fraction') && EllipseFitThreshold < 0.5
                error('EllipseFitThreshold cannot be smaller than 0.5, if the ThresholdType is "Fraction"!')
            end
            
            OutArrayStruct = struct();
            OutStruct = struct();
            OutStructAll = struct();
            Struct = struct([]);
            
            Channel = obj.get_channel(HeightChannelName);
            if isempty(Channel)
                warning('No Channel "Processed" found. Need flattened height channel for snap to local maximum')
                return
            end
            
            % All calculations assume same size per pixel so convert
            % Channel to same physical size per pixel and use multipliers
            % to convert back afterwards
            [Channel,XMult,YMult] = AFMBaseClass.resize_channel_to_same_size_per_pixel(Channel);
            
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
                % Adapt Segment positions to new image dimensions
                ConvSegPos = obj.Segment(i).ROIObject.Position;
                ConvSegPos = [ConvSegPos(:,1)*YMult ConvSegPos(:,2)*XMult];
                
                Polyline = drawpolyline(Ax.Parent,'Position',ConvSegPos);
                FirstMask = Polyline.createMask;
                TempMask = FirstMask;
                [TempCol,TempRow] = find(FirstMask == 1);
                SegmentPixelIndizes{i} = [TempRow TempCol];
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
                
                % Adapt Segment positions to new image dimensions
                ConvSegPos = obj.Segment(i).ROIObject.Position;
                ConvSegPos = [ConvSegPos(:,1)*YMult ConvSegPos(:,2)*XMult];
                SegmentPositions = ConvSegPos;
                
                LocalDirectionVector = AFMBaseClass.find_local_direction_vector_in_ordered_vector_list(SegmentPositions,SmoothingWindowSize,'flat');
                LocalHeight = zeros(length(SegmentPositions(:,1)),1);
                LocalProminence = zeros(length(SegmentPositions(:,1)),1);
                LocalWidthHalfHeight = zeros(length(SegmentPositions(:,1)),1);
                LocalWidthHalfProminence = zeros(length(SegmentPositions(:,1)),1);
                LocalArea = zeros(length(SegmentPositions(:,1)),1);
                LocalWidthBase = zeros(length(SegmentPositions(:,1)),1);
                LocalPosition = zeros(length(SegmentPositions(:,1)),2);
                LocalEllipse_a = zeros(length(SegmentPositions(:,1)),1);
                LocalEllipse_b = zeros(length(SegmentPositions(:,1)),1);
                LocalEllipse_AspectRatio = zeros(length(SegmentPositions(:,1)),1);
                LocalEllipse_Area = zeros(length(SegmentPositions(:,1)),1);
                LocalEllipse_Height = zeros(length(SegmentPositions(:,1)),1);
                LocalEllipse_WidthHalfHeight = zeros(length(SegmentPositions(:,1)),1);
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
                        Ax1.Legend.Location = 'northeast';
                        Ax2 = subplot('Position',[.6 .55 .3 .4]);
%                         findpeaks(LocalProfile,LocalDistance,'Annotate','extents',...
%                             'WidthReference','halfprom','MinPeakDistance',MinPeakDistanceMeters);
%                         hold on
%                         plot(ForbiddenLocalDistance,ForbiddenProfile,'rX',...
%                             Locations(WinnerIndex),1.05*Heights(WinnerIndex),'rv',...
%                             'MarkerSize',16,'MarkerFaceColor','g')
%                         Ax2.Legend.String{end-1} = 'Forbidden Points';
%                         Ax2.Legend.String{end} = 'Chosen Point';
%                         Ax2.Legend.Location = 'northeast';
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
                        LocalPosition(j,:) = [nan nan];
                        continue
                    end
                    % find out local area
                    PeakIndex = find(LocalDistance==Locations(WinnerIndex));
                    if isequal(ThresholdType,'Fraction')
                        Thresh = LowerEndThreshold.*LocalProfile(PeakIndex);
                        EllThresh = EllipseFitThreshold.*LocalProfile(PeakIndex);
                    else
                        Thresh = LowerEndThreshold;
                        EllThresh = EllipseFitThreshold;
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
                    
                    
                    % To the left
                    kk = 1;
                    EllLeftBoundIndex = [];
                    while (PeakIndex-kk)>0 && isempty(EllLeftBoundIndex)
                        if LocalProfile(PeakIndex-kk) < EllThresh
                            EllLeftBoundIndex = PeakIndex - kk;
                        end
                        kk = kk + 1;
                    end
                    EllRightBoundIndex = [];
                    % To the right
                    kk = 1;
                    while (PeakIndex+kk)<=length(LocalProfile) && isempty(EllRightBoundIndex)
                        if LocalProfile(PeakIndex+kk) < EllThresh
                            EllRightBoundIndex = PeakIndex + kk;
                        end
                        kk = kk + 1;
                    end
                    % If info is sufficient, fit ellipse
                    if ~isempty(EllLeftBoundIndex) && ~isempty(EllRightBoundIndex)
                        X = LocalDistance(EllLeftBoundIndex:EllRightBoundIndex);
                        Y = LocalProfile(EllLeftBoundIndex:EllRightBoundIndex);
                        if Verbose && mod(j,KeyFrames)==0
                            EllResultStruct = EllipseFit_fit_ellipse(X,Y,Ax2);
                            hold off
                            drawnow
                        else
                            EllResultStruct = EllipseFit_fit_ellipse(X,Y);
                        end
                    else
                        EllResultStruct.status = 'AllWentWrong';
                    end
                    
                    % Write down results
                    LocalHeight(j) = Heights(WinnerIndex);
                    LocalProminence(j) = Prom(WinnerIndex);
                    LocalWidthHalfHeight(j) = TempWidthHalfHeight(WinnerIndex);
                    LocalWidthHalfProminence(j) = TempWidthHalfProm(WinnerIndex);
                    LocalPosition(j,1) = LocalX(WinnerIndex);
                    LocalPosition(j,2) = LocalY(WinnerIndex);
                    if isequal(EllResultStruct.status,'')
                        LocalEllipse_a(j) = EllResultStruct.a;
                        LocalEllipse_b(j) = EllResultStruct.b;
                        LocalEllipse_AspectRatio(j) = EllResultStruct.b./EllResultStruct.a;
                        LocalEllipse_Area(j) = EllResultStruct.a.*EllResultStruct.b.*pi;
                        LocalEllipse_Height(j) = 2*EllResultStruct.b;
                        LocalEllipse_WidthHalfHeight(j) = 2*EllResultStruct.a;
                    else
                        LocalEllipse_a(j) = nan;
                        LocalEllipse_b(j) = nan;
                        LocalEllipse_AspectRatio(j) = nan;
                        LocalEllipse_Area(j) = nan;
                        LocalEllipse_Height(j) = nan;
                        LocalEllipse_WidthHalfHeight(j) = nan;
                    end
                    
                    
                    % Find corresponding Pixel by distance
%                     DistanceVector = SegmentPixelIndizes{i} - [LocalPosition(j,1) LocalPosition(j,2)];
                    DistanceVector = SegmentPixelIndizes{i} - [SegmentPositions(j,1) SegmentPositions(j,2)];
                    Distances = vecnorm(DistanceVector,2,2);
                    [~,IndexOfClosest] = min(Distances);
                    try
                        TempCorrespondingPixelIndex(j,:) = SegmentPixelIndizes{i}(IndexOfClosest,:);
                    catch ME
                        warning(ME.message)
                    end
                    
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
                
                % Now rescale back all relevant pixel positions so they
                % make sense in the context of the original channel
                % dimensions. Skip if XMult==YMult==1
                if XMult == 1 && YMult==1
                    % do nothing
                else
                    SegmentPositions = [SegmentPositions(:,1)/YMult SegmentPositions(:,2)/XMult];
                    TempCorrespondingPixelIndex = [TempCorrespondingPixelIndex(:,1)/YMult TempCorrespondingPixelIndex(:,2)/XMult];
                    
                    % Bin integer pixel positions
                    integerListX = 1:round(Channel.NumPixelsX/XMult);
                    integerListY = 1:round(Channel.NumPixelsY/YMult);
                    TempCorrespondingPixelIndex(:,2) = binToNearest(TempCorrespondingPixelIndex(:,2),integerListX);
                    TempCorrespondingPixelIndex(:,1) = binToNearest(TempCorrespondingPixelIndex(:,1),integerListY);
                    SegmentPixelIndizes{i}(:,2) = binToNearest(SegmentPixelIndizes{i}(:,2)/XMult,integerListX);
                    SegmentPixelIndizes{i}(:,1) = binToNearest(SegmentPixelIndizes{i}(:,1)/YMult,integerListY);
                end
                
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
                Struct(k).Mean_WidthHalfHeight = mean(LocalWidthHalfHeight,'omitnan');
                obj.Segment(i).Mean_WidthHalfHeight = mean(LocalWidthHalfHeight,'omitnan');
                Struct(k).Median_WidthHalfHeight = median(LocalWidthHalfHeight,'omitnan');
                obj.Segment(i).Median_WidthHalfHeight = median(LocalWidthHalfHeight,'omitnan');
                Struct(k).Mean_Height = mean(LocalHeight,'omitnan');
                obj.Segment(i).Mean_Height = mean(LocalHeight,'omitnan');
                Struct(k).Median_Height = median(LocalHeight,'omitnan');
                obj.Segment(i).Median_Height = median(LocalHeight,'omitnan');
                Struct(k).Mean_WidthHalfProminence = mean(LocalWidthHalfProminence,'omitnan');
                obj.Segment(i).Mean_WidthHalfProminence = mean(LocalWidthHalfProminence,'omitnan');
                Struct(k).Median_WidthHalfProminence = median(LocalWidthHalfProminence,'omitnan');
                obj.Segment(i).Median_WidthHalfProminence = median(LocalWidthHalfProminence,'omitnan');
                Struct(k).Mean_Prominence = mean(LocalProminence,'omitnan');
                obj.Segment(i).Mean_Prominence = mean(LocalProminence,'omitnan');
                Struct(k).Median_Prominence = median(LocalProminence,'omitnan');
                obj.Segment(i).Median_Prominence = median(LocalProminence,'omitnan');
                Struct(k).Area = LocalArea;
                obj.Segment(i).Area = LocalArea;
                Struct(k).Mean_Area = mean(LocalArea,'omitnan');
                obj.Segment(i).Mean_Area = mean(LocalArea,'omitnan');
                Struct(k).Median_Area = median(LocalArea,'omitnan');
                obj.Segment(i).Median_Area = median(LocalArea,'omitnan');
                Struct(k).WidthBase = LocalWidthBase;
                obj.Segment(i).WidthBase = LocalWidthBase;
                Struct(k).Mean_WidthBase = mean(LocalWidthBase,'omitnan');
                obj.Segment(i).Mean_WidthBase = mean(LocalWidthBase,'omitnan');
                Struct(k).Median_WidthBase = median(LocalWidthBase,'omitnan');
                obj.Segment(i).Median_WidthBase = median(LocalWidthBase,'omitnan');
                Struct(k).AspectRatioHalfHeight = LocalHeight./LocalWidthHalfProminence;
                obj.Segment(i).AspectRatioHalfHeight = LocalHeight./LocalWidthHalfProminence;
                Struct(k).Mean_AspectRatioHalfHeight = mean(LocalHeight./LocalWidthHalfProminence,'omitnan');
                obj.Segment(i).Mean_AspectRatioHalfHeight = mean(LocalHeight./LocalWidthHalfProminence,'omitnan');
                Struct(k).Median_AspectRatioHalfHeight = median(LocalHeight./LocalWidthHalfProminence,'omitnan');
                obj.Segment(i).Median_AspectRatioHalfHeight = median(LocalHeight./LocalWidthHalfProminence,'omitnan');
                Struct(k).AspectRatioBaseHeight = LocalHeight./LocalWidthBase;
                obj.Segment(i).AspectRatioBaseHeight = LocalHeight./LocalWidthBase;
                Struct(k).Mean_AspectRatioBaseHeight = mean(LocalHeight./LocalWidthBase,'omitnan');
                obj.Segment(i).Mean_AspectRatioBaseHeight = mean(LocalHeight./LocalWidthBase,'omitnan');
                Struct(k).Median_AspectRatioBaseHeight = median(LocalHeight./LocalWidthBase,'omitnan');
                obj.Segment(i).Median_AspectRatioBaseHeight = median(LocalHeight./LocalWidthBase,'omitnan');
                Struct(k).AreaDerivedDiameter = real(sqrt(LocalArea/pi).*2);
                obj.Segment(i).AreaDerivedDiameter = real(sqrt(LocalArea/pi).*2);
                Struct(k).Mean_AreaDerivedDiameter = mean(real(sqrt(LocalArea/pi).*2),'omitnan');
                obj.Segment(i).Mean_AreaDerivedDiameter = mean(real(sqrt(LocalArea/pi).*2),'omitnan');
                Struct(k).Median_AreaDerivedDiameter = median(real(sqrt(LocalArea/pi).*2),'omitnan');
                obj.Segment(i).Median_AreaDerivedDiameter = median(real(sqrt(LocalArea/pi).*2),'omitnan');
                Struct(k).Ellipse_a = LocalEllipse_a;
                obj.Segment(i).Ellipse_a = LocalEllipse_a;
                Struct(k).Mean_Ellipse_a = mean(LocalEllipse_a,'omitnan');
                obj.Segment(i).Mean_Ellipse_a = mean(LocalEllipse_a,'omitnan');
                Struct(k).Median_Ellipse_a = median(LocalEllipse_a,'omitnan');
                obj.Segment(i).Median_Ellipse_a = median(LocalEllipse_a,'omitnan');
                Struct(k).Ellipse_b = LocalEllipse_b;
                obj.Segment(i).Ellipse_b = LocalEllipse_b;
                Struct(k).Mean_Ellipse_b = mean(LocalEllipse_b,'omitnan');
                obj.Segment(i).Mean_Ellipse_b = mean(LocalEllipse_b,'omitnan');
                Struct(k).Median_Ellipse_b = median(LocalEllipse_b,'omitnan');
                obj.Segment(i).Median_Ellipse_b = median(LocalEllipse_b,'omitnan');
                Struct(k).Ellipse_AspectRatio = LocalEllipse_AspectRatio;
                obj.Segment(i).Ellipse_AspectRatio = LocalEllipse_AspectRatio;
                Struct(k).Mean_Ellipse_AspectRatio = mean(LocalEllipse_AspectRatio,'omitnan');
                obj.Segment(i).Mean_Ellipse_AspectRatio = mean(LocalEllipse_AspectRatio,'omitnan');
                Struct(k).Median_Ellipse_AspectRatio = median(LocalEllipse_AspectRatio,'omitnan');
                obj.Segment(i).Median_Ellipse_AspectRatio = median(LocalEllipse_AspectRatio,'omitnan');
                Struct(k).Ellipse_Area = LocalEllipse_Area;
                obj.Segment(i).Ellipse_Area = LocalEllipse_Area;
                Struct(k).Mean_Ellipse_Area = mean(LocalEllipse_Area,'omitnan');
                obj.Segment(i).Mean_Ellipse_Area = mean(LocalEllipse_Area,'omitnan');
                Struct(k).Median_Ellipse_Area = median(LocalEllipse_Area,'omitnan');
                obj.Segment(i).Median_Ellipse_Area = median(LocalEllipse_Area,'omitnan');
                Struct(k).Ellipse_Height = LocalEllipse_Height;
                obj.Segment(i).Ellipse_Height = LocalEllipse_Height;
                Struct(k).Mean_Ellipse_Height = mean(LocalEllipse_Height,'omitnan');
                obj.Segment(i).Mean_Ellipse_Height = mean(LocalEllipse_Height,'omitnan');
                Struct(k).Median_Ellipse_Height = median(LocalEllipse_Height,'omitnan');
                obj.Segment(i).Median_Ellipse_Height = median(LocalEllipse_Height,'omitnan');
                Struct(k).Ellipse_WidthHalfHeight = LocalEllipse_WidthHalfHeight;
                obj.Segment(i).Ellipse_WidthHalfHeight = LocalEllipse_WidthHalfHeight;
                Struct(k).Mean_Ellipse_WidthHalfHeight = mean(LocalEllipse_WidthHalfHeight,'omitnan');
                obj.Segment(i).Mean_Ellipse_WidthHalfHeight = mean(LocalEllipse_WidthHalfHeight,'omitnan');
                Struct(k).Median_Ellipse_WidthHalfHeight = median(LocalEllipse_WidthHalfHeight,'omitnan');
                obj.Segment(i).Median_Ellipse_WidthHalfHeight = median(LocalEllipse_WidthHalfHeight,'omitnan');
                
                Struct(k).SegmentPixelIndex = SegmentPixelIndizes{i};
                obj.Segment(i).SegmentPixelIndex = SegmentPixelIndizes{i};
                Struct(k).CorrespondingPixelIndex = TempCorrespondingPixelIndex;
                obj.Segment(i).CorrespondingPixelIndex = TempCorrespondingPixelIndex;
                TempCorrespondingPixelIndex = [];
                
                LocalPositionRealWorld = [SegmentPositions(:,1).*Channel.ScanSizeX/(Channel.NumPixelsX/YMult) ...
                    SegmentPositions(:,2).*Channel.ScanSizeY/(Channel.NumPixelsY/XMult)];
                RelativePixel(1) = 0;
                RelativeRealWorld(1) = 0;
                RelativePixel(2:length(SegmentPositions)) = cumsum(vecnorm(...
                    SegmentPositions(2:end,:) - SegmentPositions(1:end-1,:),2,2));
                RelativeRealWorld(2:length(SegmentPositions)) = cumsum(vecnorm(...
                    LocalPositionRealWorld(2:end,:) - LocalPositionRealWorld(1:end-1,:),2,2));
                obj.Segment(i).RelativePixelPosition = reshape(RelativePixel,[],1);
                Struct(k).RelativePixelPosition = reshape(RelativePixel,[],1);
                obj.Segment(i).RelativePosition = reshape(RelativeRealWorld,[],1);
                Struct(k).RelativePosition = reshape(RelativeRealWorld,[],1);
                
                RelativePixel = [];
                RelativeRealWorld = [];
                
                k = k + 1;
            end
            
            % Now for some tedious data shoveling. Life is pain            
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
                OutStruct(i).Mean_Height = mean(OutStruct(i).Height,'omitnan');
                OutStruct(i).Median_Height = median(OutStruct(i).Height,'omitnan');
                OutStruct(i).Mean_WidthHalfHeight = mean(OutStruct(i).WidthHalfHeight,'omitnan');
                OutStruct(i).Median_WidthHalfHeight = median(OutStruct(i).WidthHalfHeight,'omitnan');
                OutStruct(i).Mean_Prominence = mean(OutStruct(i).Prominence,'omitnan');
                OutStruct(i).Median_Prominence = median(OutStruct(i).Prominence,'omitnan');
                OutStruct(i).Mean_WidthHalfProminence = mean(OutStruct(i).WidthHalfProminence,'omitnan');
                OutStruct(i).Median_WidthHalfProminence = median(OutStruct(i).WidthHalfProminence,'omitnan');
                OutStruct(i).Mean_Area = mean(OutStruct(i).Area,'omitnan');
                OutStruct(i).Median_Area = median(OutStruct(i).Area,'omitnan');
                OutStruct(i).Mean_WidthBase = mean(OutStruct(i).WidthBase,'omitnan');
                OutStruct(i).Median_WidthBase = median(OutStruct(i).WidthBase,'omitnan');
                OutStruct(i).Mean_AspectRatioHalfHeight = mean(OutStruct(i).AspectRatioHalfHeight,'omitnan');
                OutStruct(i).Median_AspectRatioHalfHeight = median(OutStruct(i).AspectRatioHalfHeight,'omitnan');
                OutStruct(i).Mean_AspectRatioBaseHeight = mean(OutStruct(i).AspectRatioBaseHeight,'omitnan');
                OutStruct(i).Median_AspectRatioBaseHeight = median(OutStruct(i).AspectRatioBaseHeight,'omitnan');
                OutStruct(i).Mean_AreaDerivedDiameter = mean(OutStruct(i).AreaDerivedDiameter,'omitnan');
                OutStruct(i).Median_AreaDerivedDiameter = median(OutStruct(i).AreaDerivedDiameter,'omitnan');
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
            
             obj.map_fiber_segment_properties_to_image_pixels('Median');
            
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
        
        function TagList = read_out_tag_list(obj,Indizes)
            
            if nargin < 2
                NumSegments = numel(obj.Segment); 
                Indizes = [1:NumSegments];
            end
            
            LoopIndizes = reshape(Indizes,1,[]);
            TagString = '';
            for i=LoopIndizes
                if ~isfield(obj.Segment(i),'Tags')
                    continue
                end
                TagString = [TagString obj.Segment(i).Tags];
            end
            
            TagList = AFMBaseClass.split_tags_string(TagString);
            
        end
        
        function add_tag_to_segment(obj,TagString,Indizes)
            
            if nargin < 3
                NumSegments = numel(obj.Segment);
                Indizes = 1:NumSegments;
            end
            
            IsValid = Experiment.check_segment_tag_validity(TagString);
            
            if ~IsValid
                error([TagString ' is not a valid tag string'])
            end
            
            LoopIndizes = reshape(Indizes, 1, []);
            Default = '';
            TagField = 'Tags';
            
            % Check if 'Tags' field exists and add it if it doesn't
            if ~isfield(obj.Segment, TagField)
                [obj.Segment(1:numel(obj.Segment)).(TagField)] = deal(Default);
            end
            
            % Add TagString to the Tags field of specified segments
            for i = LoopIndizes
                TagList = obj.read_out_tag_list(i);
                if sum(contains(TagList,TagString))
                    continue
                end
                obj.Segment(i).Tags = append(obj.Segment(i).Tags, [TagString '**']);
            end
        end
        
        function remove_tag_from_segment(obj, TagString, Indizes)
            
            if nargin < 3
                NumSegments = numel(obj.Segment);
                Indizes = 1:NumSegments;
            end
            
            IsValid = Experiment.check_segment_tag_validity(TagString);
            
            if ~IsValid
                error([TagString ' is not a valid tag string'])
            end
            
            LoopIndizes = reshape(Indizes, 1, []);
            TagField = 'Tags';
            
            % Remove TagString from the Tags field of specified segments
            for i = LoopIndizes
                if ~isfield(obj.Segment(i), TagField)
                    continue
                end
                
                TagList = obj.read_out_tag_list(i);
                TagList = TagList(cellfun(@(tag) ~isequal(tag, TagString), TagList)); % Remove the specified tag with exact one-to-one match
                
                % Reconstruct the Tags field without the specified tag
                obj.Segment(i).Tags = [strjoin(TagList, '**') '**'];
            end
        end
        
        function [slope, radius_of_curvature, slope_direction, CorrectedRadius] = localSurfaceFit(obj, HeightChannelName, Radius)
            %[slope, radius_of_curvature, slope_direction, CorrectedRadius] = localSurfaceFit(obj, HeightChannelName, Radius)
            %
            %LOCALSURFACEFIT Calculates the local slope and curvature of an AFMBaseClass object's height image
            %   This method takes an AFMBaseClass object instance (obj) and a HeightChannelName as input,
            %   reads the height image using the get_channel method, and computes the local slope and
            %   curvature within a given Radius. It also creates a map of the slope direction.
            %
            %   Input:
            %       obj - An AFMBaseClass object instance
            %       HeightChannelName - The name of the height channel
            %       Radius - The radius within which the surface fit is determined for each area
            %
            %   Output:
            %       slope - The local slope as a first-order surface fit
            %       curvature - The local curvature as a second-order surface fit
            %       slope_direction - A map of the direction of the slope
            
            % Read the height image from the object using the get_channel method
            height_channel = obj.get_channel(HeightChannelName);
            
            orig_num_pixels_x = height_channel.NumPixelsX;
            orig_num_pixels_y = height_channel.NumPixelsY;
            
            height_channel = AFMImage.resize_channel_to_padded_same_size_per_pixel_square_image(height_channel);
            
            height_image = height_channel.Image;
            scan_size_x = height_channel.ScanSizeX;
            scan_size_y = height_channel.ScanSizeY;
            num_pixels_x = height_channel.NumPixelsX;
            num_pixels_y = height_channel.NumPixelsY;
            
            required_data_points = 6;
            min_radius = AFMBaseClass.calculateMinimumRadius(height_channel, required_data_points);
            if isscalar(Radius)
                if min_radius > Radius
                    Radius = min_radius;
                    warning(['The ROI-Radius you chose for curvature fitting is too small. Replaced with the minimum radius of ' num2str(min_radius)]);
                end
            else
                CompRadius = Radius;
                MinRadVector = min_radius.*ones(length(CompRadius),1);
                Radius = max([MinRadVector CompRadius],[],2);
                if ~isequal(Radius,CompRadius)
                    warning(['Some or all of the ROI-Radii you chose for curvature fitting are too small. Replaced with the minimum radius of ' num2str(min_radius)]);
                end
                TempRadiusMap = obj.convert_data_list_to_map(Radius);
                RadiusMap = ones(size(height_image)).*min_radius;
                RadiusMap(1:size(TempRadiusMap,1),1:size(TempRadiusMap,2)) = TempRadiusMap;
            end
            
            % Calculate the pixel size in x and y directions
            pixel_size_x = scan_size_x / (num_pixels_x - 1);
            pixel_size_y = scan_size_y / (num_pixels_y - 1);
            
            % Initialize the output matrices
            slope = zeros(num_pixels_y, num_pixels_x);
            radius_of_curvature = zeros(num_pixels_y, num_pixels_x);
            slope_direction = zeros(num_pixels_y, num_pixels_x);
            
            PixRadiusX = round(Radius/pixel_size_x);
            PixRadiusY = round(Radius/pixel_size_y);
            
            
            % Initialize the warning count
            warningCount = struct('fitPlaneWarnings', 0, 'fitParaboloidWarnings', 0);
            
            % Suppress warnings for fitting functions
            plane_warning_state = warning('off', 'MATLAB:rankDeficientMatrix');
            singular_matrix_warning_state = warning('off', 'MATLAB:singularMatrix');
            ill_conditioned_matrix_warning_state = warning('off', 'MATLAB:nearlySingularMatrix');
            
            
            
            % Loop through each pixel in the height image
            for i = 1:num_pixels_y
                for j = 1:num_pixels_x
                    if ~isscalar(Radius)
                        PixRadiusX = round(RadiusMap(i,j)/pixel_size_x);
                        PixRadiusY = round(RadiusMap(i,j)/pixel_size_y);
                    end
                    % Determine the region of interest (ROI) around the current pixel
                    [x_roi, y_roi] = meshgrid(max(1, j-PixRadiusX):min(num_pixels_x, j+PixRadiusX), ...
                        max(1, i-PixRadiusY):min(num_pixels_y, i+PixRadiusY));
                    
                    % Extract the height values within the ROI
                    z_roi = height_image(sub2ind(size(height_image), y_roi(:), x_roi(:)));
                    
                    % Fit a first-order surface (plane) to the ROI data
                    [p1, ~, ~, ~, plane_warnings] = AFMBaseClass.fitPlane(x_roi(:) * pixel_size_x, y_roi(:) * pixel_size_y, z_roi);
                    warningCount.fitPlaneWarnings = warningCount.fitPlaneWarnings + plane_warnings;
                    
                    % Fit a second-order surface (paraboloid) to the ROI data
                    [p2, ~, ~, ~, paraboloid_warnings] = AFMBaseClass.fitParaboloid(x_roi(:) * pixel_size_x, y_roi(:) * pixel_size_y, z_roi);
                    warningCount.fitParaboloidWarnings = warningCount.fitParaboloidWarnings + paraboloid_warnings;
                    
                    % Calculate the local slope and curvature at the current pixel
                    temp_slope(i, j) = norm(p1(1:2));
                    temp_radius_of_curvature(i, j) = sign(p2(1) + p2(2))/sqrt(p2(1)^2 + p2(2)^2);
                    
                    % Calculate the slope direction and convert it to degrees
                    temp_slope_direction(i, j) = atan2d(p1(2), p1(1));
                end
            end
            % Restore the original warning state
            warning(plane_warning_state);
            warning(singular_matrix_warning_state);
            warning(ill_conditioned_matrix_warning_state);
            
            slope = temp_slope(1:orig_num_pixels_x,1:orig_num_pixels_y);
            slope_direction = temp_slope_direction(1:orig_num_pixels_x,1:orig_num_pixels_y);
            radius_of_curvature = temp_radius_of_curvature(1:orig_num_pixels_x,1:orig_num_pixels_y);
            
            CorrectedRadius = Radius;
            
        end
        
        function localSurfaceFit_ClassWrapper(obj, HeightChannelName, Radius, KeepOldResults)
            %localSurfaceFit_ClassWrapper(obj, HeightChannelName, Radius, KeepOldResults)
            
            % Calculate the metrics
            [slope, radius_of_curvature, slope_direction, CorrectedRadius] = localSurfaceFit(obj, HeightChannelName, Radius);
            
            % Determine Channel Names
            if isscalar(CorrectedRadius)
                LRoC_Name = sprintf('Local Radius of Curvature Kernel-R. = %.2e',CorrectedRadius);
                S_Name = sprintf('Local Slope Kernel-R. = %.2e',CorrectedRadius);
                SD_Name = sprintf('Local Slope Direction Kernel-R. = %.2e',CorrectedRadius);
            else
                LRoC_Name = 'Dynamic Kernel Radius of Curvature';
                S_Name = 'Dynamic Kernel Slope';
                SD_Name = 'Dynamic Kernel Slope Direction';
            end
            
            % Write new Channels
            obj.add_channel(...
                obj.create_standard_channel(radius_of_curvature,LRoC_Name,'m'), ~KeepOldResults);
            obj.add_channel(...
                obj.create_standard_channel(slope,S_Name,'m/m'), ~KeepOldResults);
            obj.add_channel(...
                obj.create_standard_channel(slope_direction,SD_Name,'°'), ~KeepOldResults);
            
        end
        
        function [ContactRadii,ContactArea] = local_contact_area_ClassWrapper(...
                obj, TipHeightChannel, HeightChannel, IndendationDepthChannel,KeepOldResults)
            
            THChan = TipHeightChannel;
            SHChan = HeightChannel;
            SIDChan = IndendationDepthChannel;
            
            ContactArea = AFMBaseClass.calculateContactArea(THChan,SHChan,SIDChan,'Verbose',false);
            ContactRadii = sqrt(ContactArea./pi);
            
            % Write new Channels
            obj.add_channel(...
                obj.create_standard_channel(ContactArea,'Dynamic Contact Area','m^2'), ~KeepOldResults);
            obj.add_channel(...
                obj.create_standard_channel(ContactRadii,'Dynamic Contact Radius','m'), ~KeepOldResults);
            
        end
        
        function create_logical_mask_channel_from_segment(obj,SegmentName,MaskName,DilationMeters)
            
            if nargin < 4
                DilationMeters = 0;
            end
            
            Indices = find(contains({obj.Segment.Name},SegmentName));
            
            Mask = zeros(obj.NumPixelsX,obj.NumPixelsY);
            
            F = figure;
            imshow(Mask)
            
            for i=Indices
                CurrentDrawMode = lower(obj.Segment(i).Type);
                switch CurrentDrawMode
                    case 'line'
                        ROIObjects{i} = drawline('Position',obj.Segment(i).ROIObject.Position);
                    case 'freehand'
                        ROIObjects{i} = drawfreehand('Position',obj.Segment(i).ROIObject.Position);
                    case 'square'
                        ROIObjects{i} = drawsquare('Position',obj.Segment(i).ROIObject.Position);
                    case 'circle'
                        ROIObjects{i} = drawcircle('Center',obj.Segment(i).ROIObject.Center,...
                            'Radius',obj.Segment(i).ROIObject.Radius);
                    case 'ellipse'
                        ROIObjects{i} = drawellipse('Center',obj.Segment(i).ROIObject.Center,...
                            'SemiAxes',obj.Segment(i).ROIObject.SemiAxes,...
                            'RotationAngle',obj.Segment(i).ROIObject.RotationAngle);
                    case 'polygon'
                        ROIObjects{i} = drawpolygon('Position',obj.Segment(i).ROIObject.Position);
                    case 'rectangle'
                        ROIObjects{i} = drawrectangle('Position',obj.Segment(i).ROIObject.Position);
                    case 'crosshair'
                        ROIObjects{i} = drawcrosshair('Position',obj.Segment(i).ROIObject.Position);
                    case 'point'
                        ROIObjects{i} = drawpoint('Position',obj.Segment(i).ROIObject.Position);
                    case 'assisted'
                        ROIObjects{i} = drawassisted('Position',obj.Segment(i).ROIObject.Position);
                    case 'polyline'
                        ROIObjects{i} = drawpolyline('Position',obj.Segment(i).ROIObject.Position);
                end
                TempMask = ROIObjects{i}.createMask;
                
                if DilationMeters ~= 0
                    %TODO Implement size specific pixel dilation
                    TempMask = TempMask;
                end
                
                Mask = Mask | TempMask;
            end
            
            MaskChan = obj.create_standard_channel(Mask,MaskName,'logical');
            
            obj.add_channel(MaskChan,true);
            
            close(F)
        end
        
        function [Map,FitParams] = flatten_image_by_vertical_rov(obj,SourceChannelName)
            
            if nargin < 2
                SourceChannelName = 'NoChannelNameGiven';
            end
            
            if nargin == 2
                Height = obj.get_channel(SourceChannelName);
                if isempty(Height)
                    warning(sprintf('Channel %s not found, getting unprocessed height channel instead',SourceChannelName))
                    Height = obj.get_unprocessed_height_channel('Height (Trace)');
                end
            else
                Height = obj.get_unprocessed_height_channel('Height (Trace)');
            end
            
            if isempty(Height)
                warning("Could not find height channel needed for preprocessing")
                return
            end
            
            if size(Height.Image,1) < 128
                Map = imresize(Height.Image,[256 256],'nearest');
            elseif size(Height.Image,1) < 512
                Map = imresize(Height.Image,[512 512],'nearest');
            else
                Map = imresize(Height.Image,[1024 1024],'nearest');
            end
            FitParams = zeros(size(Map,1),2);
            for i=1:5
                [Map,TempFitParams] = AFMImage.subtract_line_fit_vertical_rov(Map,.2,0);
                FitParams = FitParams + TempFitParams;
            end
            Map = imresize(Map,[obj.NumPixelsX obj.NumPixelsY],'nearest');
            
        end
        
        function Map = flatten_image_by_other_channels_fitparams(obj,Source,Target)
            
            % Chack equal size
            if (size(Source.Image,1) ~= size(Target.Image,1)) ||...
                    (size(Source.Image,2) ~= size(Target.Image,2))
                
                error(sprintf('Channels %s and %s must be of the same size',SourceChannelName,TargetChannelName))
                
            end
            
            [~,FitParams] = obj.flatten_image_by_vertical_rov(Source.Name);
            
            if size(Target.Image,1) < 128
                Map = imresize(Target.Image,[256 256],'nearest');
            elseif size(Target.Image,1) < 512
                Map = imresize(Target.Image,[512 512],'nearest');
            else
                Map = imresize(Target.Image,[1024 1024],'nearest');
            end
            NumProfiles = size(Map,1);
            NumPoints = size(Map,2);
            InImage = Map;
            for i=1:NumProfiles
                LineEval = [1:NumPoints]'*FitParams(i,1) + FitParams(i,2);
                Line = InImage(i,:)';
                Line = Line - LineEval;
                InImage(i,:) = Line;
            end
            
            Map = imresize(InImage,[Target.NumPixelsX Target.NumPixelsY],'nearest');
            
        end
        
    end
    methods (Static)
        % Static main methods
        
        function OutImage = deconvolute_by_mathematical_morphology_python(InImage, ErodingGeometry)
            % Make sure Python is properly set up in MATLAB by running "pyversion" command in MATLAB
            % Save the Python function in a file named "run_deconvolution.py"
            
            % Find the full path of the Python file on the MATLAB search path
            python_file_path = which('run_deconvolution.py');
            
            if isempty(python_file_path)
                error('run_deconvolution.py not found on the MATLAB search path.');
            end
            
            % Extract the directory containing the Python file
            python_file_dir = fileparts(python_file_path);
            
            % Convert MATLAB arrays to Python numpy arrays
            InImage_np = py.numpy.array(InImage);
            ErodingGeometry_np = py.numpy.array(ErodingGeometry);
            
            % Call the Python function
            OutImage_np = py.run_deconvolution.run_deconvolution(python_file_dir, InImage_np, ErodingGeometry_np);
            
            % Convert the result back to a MATLAB array
            OutImage = double(OutImage_np);
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
        
        function OutChannel = resize_channel(InChannel,TargetRes,TransformToSquare)
            
            if nargin < 3
                TransformToSquare = false;
            end
            
            OutChannel = InChannel;
            
            if length(TargetRes) == 2
                OutChannel.Image = imresize(InChannel.Image,[TargetRes(1) TargetRes(2)],'bilinear');
                OutChannel.NumPixelsX = size(OutChannel.Image,1);
                OutChannel.NumPixelsY = size(OutChannel.Image,2);
                return
            end
            
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
            % [OutChannel,XMultiplier,YMultiplier] = resize_channel_to_same_size_per_pixel(InChannel)
            % 
            % Resizes the image dimension with bigger size per pixel so
            % the OutChannel always has more pixels.
            % New Image will have aspect ratio according to the ratio of
            % ScanSizeX and ScanSizeY
            
            SizePerPixelX = InChannel.ScanSizeY./InChannel.NumPixelsX;
            SizePerPixelY = InChannel.ScanSizeX./InChannel.NumPixelsY;
            
            OutChannel = InChannel;
            
            if SizePerPixelX == SizePerPixelY
                XMultiplier = 1;
                YMultiplier = 1;
                return
            elseif SizePerPixelX > SizePerPixelY
                NewNumPixels = round(InChannel.ScanSizeY*InChannel.NumPixelsY/InChannel.ScanSizeX);
                OutChannel.Image = imresize(InChannel.Image,[NewNumPixels InChannel.NumPixelsY],'bilinear');
                OutChannel.NumPixelsX = NewNumPixels;
                XMultiplier = OutChannel.NumPixelsX/InChannel.NumPixelsX;
                YMultiplier = 1;
            elseif SizePerPixelY > SizePerPixelX
                NewNumPixels = round(InChannel.ScanSizeX*InChannel.NumPixelsX/InChannel.ScanSizeY);
                OutChannel.Image = imresize(InChannel.Image,[InChannel.NumPixelsX NewNumPixels],'bilinear');
                OutChannel.NumPixelsY = NewNumPixels;
                YMultiplier = OutChannel.NumPixelsY/InChannel.NumPixelsY;
                XMultiplier = 1;
            end
        end
        
        function [OutChannel,XMultiplier,YMultiplier] = resize_channel_to_padded_same_size_per_pixel_square_image(InChannel,varargin)
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
            [OutChannel,XMultiplier,YMultiplier] = AFMBaseClass.resize_channel_to_same_size_per_pixel(InChannel);
            
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
                OutChannel.Image(end+1:end+PixelDiff,:) = PaddingValue;
                OutChannel.NumPixelsX = OutChannel.NumPixelsY;
                OutChannel.ScanSizeY = OutChannel.ScanSizeX;
            else
                PixelDiff = OutChannel.NumPixelsX - OutChannel.NumPixelsY;
                OutChannel.Image(:,end+1:end+PixelDiff) = PaddingValue;
                OutChannel.NumPixelsY = OutChannel.NumPixelsX;
                OutChannel.ScanSizeX = OutChannel.ScanSizeY;
            end
            
            OutChannel = AFMBaseClass.resize_channel(OutChannel,TargetResolution,true);
            
        end
        
        function [OutChannel1,OutChannel2] = resize_channels_to_same_physical_size_per_pixel(InChannel1,InChannel2,varargin)
            % function [OutChannel1,OutChannel2] = resize_channels_to_same_physical_size_per_pixel(InChannel1,InChannel2,varargin)
            %
            % <FUNCTION DESCRIPTION HERE>
            %
            %
            % Required inputs
            % InChannel1 ... <VARIABLE DESCRIPTION>
            % InChannel2 ... <VARIABLE DESCRIPTION>
            %
            % Name-Value pairs
            % "KeepResolution" ... <NAMEVALUE DESCRIPTION>
            % "PaddingType" ... <NAMEVALUE DESCRIPTION>
            % "EquateResolution" ... <NAMEVALUE DESCRIPTION>
            
            p = inputParser;
            p.FunctionName = "resize_channels_to_same_physical_size_per_pixel";
            p.CaseSensitive = false;
            p.PartialMatching = true;
            
            % Required inputs
            validInChannel1 = @(x)true;
            validInChannel2 = @(x)true;
            addRequired(p,"InChannel1",validInChannel1);
            addRequired(p,"InChannel2",validInChannel2);
            
            % NameValue inputs
            defaultKeepResolution = true;
            defaultPaddingType = 'min';
            defaultEquateResolution = false;
            validKeepResolution = @(x)true;
            validPaddingType = @(x)any(validatestring(x,{'Min','Max','Zero'}));
            validEquateResolution = @(x)true;
            addParameter(p,"KeepResolution",defaultKeepResolution,validKeepResolution);
            addParameter(p,"PaddingType",defaultPaddingType,validPaddingType);
            addParameter(p,"EquateResolution",defaultEquateResolution,validEquateResolution);
            
            parse(p,InChannel1,InChannel2,varargin{:});
            
            % Assign parsing results to named variables
            InChannel1 = p.Results.InChannel1;
            InChannel2 = p.Results.InChannel2;
            KeepResolution = p.Results.KeepResolution;
            PaddingType = p.Results.PaddingType;
            EquateResolution = p.Results.EquateResolution;
            
            
            InChannel1 = AFMBaseClass.resize_channel_to_padded_same_size_per_pixel_square_image(InChannel1);
            InChannel2 = AFMBaseClass.resize_channel_to_padded_same_size_per_pixel_square_image(InChannel2);
            
            if EquateResolution
                if InChannel1.NumPixelsX >= InChannel2.NumPixelsX
                    InChannel2 = AFMBaseClass.resize_channel(InChannel2,InChannel1.NumPixelsX);
                else
                    InChannel1 = AFMBaseClass.resize_channel(InChannel1,InChannel2.NumPixelsX);
                end
            end
            
            SizePerPixelX1 = InChannel1.ScanSizeX/InChannel1.NumPixelsX;
            SizePerPixelY1 = InChannel1.ScanSizeY/InChannel1.NumPixelsY;
            SizePerPixelX2 = InChannel2.ScanSizeX/InChannel2.NumPixelsX;
            SizePerPixelY2 = InChannel2.ScanSizeY/InChannel2.NumPixelsY;
            
            MultiplicatorX = SizePerPixelX2/SizePerPixelX1;
            MultiplicatorY = SizePerPixelY2/SizePerPixelY1;
            
            if isequal(lower(PaddingType),'min')
                PaddingValue1 = min(InChannel1.Image,[],'all');
                PaddingValue2 = min(InChannel2.Image,[],'all');
            elseif isequal(lower(PaddingType),'max')
                PaddingValue1 = max(InChannel1.Image,[],'all');
                PaddingValue2 = max(InChannel2.Image,[],'all');
            elseif isequal(lower(PaddingType),'zero')
                PaddingValue1 = 0;
                PaddingValue2 = 0;
            end
            
            OutChannel1 = InChannel1;
            OutChannel2 = InChannel2;
            
            % X Direction
            if MultiplicatorX==1
                % Just take it as it is
            elseif MultiplicatorX<1
                NewX = round(OutChannel2.NumPixelsX.*MultiplicatorX);
                NewX = NewX + rem(NewX,2);
                OutChannel2.Image = imresize(OutChannel2.Image,[OutChannel2.NumPixelsY NewX],'bilinear');
                OutChannel2.NumPixelsX = NewX;
                PadX = (InChannel2.NumPixelsX - NewX)/2;
                OutChannel2.Image = padarray(OutChannel2.Image,[0 PadX],PaddingValue2,'both');
                OutChannel2.NumPixelsX = NewX + 2*PadX;
                OutChannel2.ScanSizeX = OutChannel2.ScanSizeX*OutChannel2.NumPixelsX/NewX;
                % Y Direction
                NewY = round(OutChannel2.NumPixelsY.*MultiplicatorY);
                NewY = NewY + rem(NewY,2);
                OutChannel2.Image = imresize(OutChannel2.Image,[NewY OutChannel2.NumPixelsX],'bilinear');
                OutChannel2.NumPixelsY = NewY;
                PadY = (InChannel2.NumPixelsY - NewY)/2;
                OutChannel2.Image = padarray(OutChannel2.Image,[PadY 0],PaddingValue2,'both');
                OutChannel2.NumPixelsY = NewY + 2*PadY;
                OutChannel2.ScanSizeY = OutChannel2.ScanSizeY*OutChannel2.NumPixelsY/NewY;
            elseif MultiplicatorX>1
                [OutChannel2,OutChannel1] = AFMBaseClass.resize_channels_to_same_physical_size_per_pixel(...
                    InChannel2,InChannel1,...
                    'PaddingType',PaddingType,...
                    'KeepResolution',KeepResolution);
            end
            
            
            if ~KeepResolution
                InChannel
            end
            
        end
        
        function OutChannel = crop_channel_by_value_range(InChannel,MinValue,MaxValue)
            
            if nargin < 3
                MaxValue = max(InChannel.Image,[],'all');
            end
            if nargin < 2
                OutChannel = InChannel;
                return
            end
            
            Mask = (InChannel.Image <= MaxValue) & (InChannel.Image >= MinValue);
            
            XCumSum = cumsum(Mask,2);
            YCumSum = cumsum(Mask,1);
            
            MinX = 0;
            MinY = 0;
            MaxX = 0;
            MaxY = 0;
            
            
            Positions = [MinX MinY MaxX MaxY];
            
            OutChannel = AFMBaseClass.crop_channel(InChannel,Positions);
            
        end
        
        function OutChannel = crop_channel(InChannel,Positions)
            
            OutChannel = InChannel;
            
        end
        
        function OutChannel = rotate_channel(InChannel,RotDegrees)
            % OutChannel = rotate_channel(InChannel,RotDegrees)
            
            OutChannel = InChannel;
            
            if rem(abs(RotDegrees),90) ~= 0
                error('At this point only multiples of 90° rotations have been implemented')
            end
            
            OutChannel.Image = imrotate(InChannel.Image,RotDegrees);
            
        end
        
        function OutChannel = set_freehand_area_in_channel_to_value(InChannel, Value)
            
            F = figure('Color','w');
            imshow(InChannel.Image,[])
            colormap(AFMImage.define_afm_color_map)
            
            ROI = drawfreehand;
            
            Mask = ROI.createMask;
            
            OutChannel = InChannel;
            OutChannel.Image(Mask) = Value;
            
            imshowpair(InChannel.Image, OutChannel.Image, 'Montage')
            colormap(AFMImage.define_afm_color_map)
            
        end
        
        function TagList = split_tags_string(TagString)
            
            TagList = split(TagString,'**');
            
            TagList = unique(TagList);
            
            % Remove empty cells from the TagList
            emptyCells = cellfun(@isempty, TagList);
            TagList = TagList(~emptyCells);
        end
        
        function [coeffs, x_fit, z_fit, y_fit, warnings] = fitPlane(x, y, z)
            %FITPLANE Fits a first-order surface (plane) to the input data
            %   Input:
            %       x, y, z - The input data points
            %
            %   Output:
            %       coeffs
            %       coeffs - The coefficients of the fitted plane
            %       x_fit, z_fit - The fitted x and z values
            
            % Standardize x and y values
            x_mean = mean(x);
            x_std = std(x);
            y_mean = mean(y);
            y_std = std(y);
            x_stdized = (x - x_mean) / x_std;
            y_stdized = (y - y_mean) / y_std;
            
            % Formulate the design matrix A and the observation vector b
            A = [x_stdized, y_stdized, ones(size(x))];
            b = z;
            
            % Solve the least squares problem to find the coefficients of the fitted plane
            coeffs_std = A \ b;
            
            % Back-transform the coefficients to the original scale
            coeffs = zeros(3, 1);
            coeffs(1) = coeffs_std(1) / x_std;
            coeffs(2) = coeffs_std(2) / y_std;
            coeffs(3) = coeffs_std(3) - coeffs_std(1) * x_mean - coeffs_std(2) * y_mean;
            
            % Calculate the fitted x and z values in the original scale
            x_fit = x;
            y_fit = y;
            z_fit = coeffs(1) * x + coeffs(2) * y + coeffs(3);
            
            
            % Count rank deficient warnings
            [~, warnings] = lastwarn;
            warnings = double(contains(warnings, 'rankDeficientMatrix')) + double(contains(warnings, 'singularMatrix'));
        end
        
        function [coeffs, X_fit, Y_fit, Z_fit, warnings] = fitParaboloid(x, y, z)
            %FITPARABOLOID Fits a second-order surface (paraboloid) to the input data
            % Input:
            % x, y, z - The input data points
            %
            % Output:
            % coeffs - The coefficients of the fitted paraboloid
            % X_fit, Y_fit, Z_fit - The fitted x, y, and z values
            
            % Standardize x and y values
            x_mean = mean(x);
            x_std = std(x);
            y_mean = mean(y);
            y_std = std(y);
            x_stdized = (x - x_mean) / x_std;
            y_stdized = (y - y_mean) / y_std;
            
            % Formulate the design matrix A and the observation vector b
            A = [x_stdized.^2, y_stdized.^2, x_stdized.*y_stdized, x_stdized, y_stdized, ones(size(x))];
            b = z;
            
            % Solve the least squares problem to find the coefficients of the fitted paraboloid
            coeffs_std = A \ b;
            
            % Back-transform the coefficients to the original scale
            coeffs = zeros(6, 1);
            coeffs(1) = coeffs_std(1) / (x_std^2);
            coeffs(2) = coeffs_std(2) / (y_std^2);
            coeffs(3) = coeffs_std(3) / (x_std * y_std);
            coeffs(4) = coeffs_std(4) / x_std - 2 * coeffs_std(1) * x_mean / x_std + coeffs_std(3) * y_mean / y_std;
            coeffs(5) = coeffs_std(5) / y_std - 2 * coeffs_std(2) * y_mean / y_std + coeffs_std(3) * x_mean / x_std;
            coeffs(6) = coeffs_std(6) - coeffs_std(1) * x_mean^2 - coeffs_std(2) * y_mean^2 - coeffs_std(3) * x_mean * y_mean + coeffs_std(4) * x_mean + coeffs_std(5) * y_mean;
            
            % Create a meshgrid for X and Y coordinates
            [X_fit, Y_fit] = meshgrid(linspace(min(x), max(x), 100), linspace(min(y), max(y), 100));
            
            % Calculate the fitted Z values using the fit coefficients and meshgrid
            Z_fit = coeffs(1) * X_fit.^2 + coeffs(2) * Y_fit.^2 + coeffs(3) * X_fit .* Y_fit + coeffs(4) * X_fit + coeffs(5) * Y_fit + coeffs(6);
            
            
            % Count rank deficient warnings
            [~, warnings] = lastwarn;
            warnings = double(contains(warnings, 'rankDeficientMatrix')) + double(contains(warnings, 'singularMatrix'));
        end
        
        function min_radius = calculateMinimumRadius(Channel, required_data_points)
            %CALCULATEMINIMUMRADIUS Calculates the minimum radius required for fitting a second-order surface
            %   This function takes an AFMBaseClass object instance (obj) and the required number of data points
            %   for fitting a second-order surface as input and calculates the minimum radius based on the grid
            %   spacing in the x and y directions.
            %
            %   Input:
            %       Channel - An AFMBaseClass channel
            %       required_data_points - The required number of data points for fitting a second-order surface
            %
            %   Output:
            %       min_radius - The minimum radius required for fitting a second-order surface
            
            % Get the grid spacing in the x and y directions
            num_pixels_y = Channel.NumPixelsX;
            num_pixels_x = Channel.NumPixelsY;
            scan_size_x = Channel.ScanSizeX;
            scan_size_y = Channel.ScanSizeY;
            
            grid_spacing_x = scan_size_x / (num_pixels_x - 1);
            grid_spacing_y = scan_size_y / (num_pixels_y - 1);
            
            % Calculate the minimum radius in the x and y directions
            min_radius_x = grid_spacing_x*ceil(sqrt(required_data_points / (num_pixels_x - 1)));
            min_radius_y = grid_spacing_y*ceil(sqrt(required_data_points / (num_pixels_y - 1)));
            
            % Take the maximum of the two minimum radii to ensure enough data points
            min_radius = max(min_radius_x, min_radius_y);
        end
        
        function contactAreas = calculateContactArea(tipStruct, sampleStruct, indentationStruct, varargin)
            % contactAreas = calculateContactArea(tipStruct, sampleStruct, indentationStruct, varargin)
            % Calculate contact area between AFM tip and sample based on indentation depth
            %
            % Inputs:
            % tipStruct - struct representing the AFM tip height map
            % sampleStruct - struct representing the sample height map
            % indentationStruct - struct representing the indentation depth map
            % 'Verbose', false, @islogical
            %
            % Outputs:
            % contactAreas - contact area for each pixel position
            
            % Parse input arguments
            p = inputParser;
            addParameter(p, 'Verbose', false, @islogical);
            parse(p, varargin{:});
            verbose = p.Results.Verbose;
            
            if verbose
                figure
            end
            
            % Extract image data
            tipHeightMap = tipStruct.Image;
            sampleHeightMap = sampleStruct.Image;
            indentationDepthMap = indentationStruct.Image;
            
            % Extract scan sizes and resolutions
            tipScanSizeX = tipStruct.ScanSizeX;
            tipScanSizeY = tipStruct.ScanSizeY;
            tipNumPixelsX = tipStruct.NumPixelsY;
            tipNumPixelsY = tipStruct.NumPixelsX;
            
            sampleScanSizeX = sampleStruct.ScanSizeX;
            sampleScanSizeY = sampleStruct.ScanSizeY;
            sampleNumPixelsX = sampleStruct.NumPixelsY;
            sampleNumPixelsY = sampleStruct.NumPixelsX;
            
            % Determine the physical sizes of both images
            tipPhysicalSizeX = tipScanSizeX;
            tipPhysicalSizeY = tipScanSizeY;
            
            samplePhysicalSizeX = sampleScanSizeX;
            samplePhysicalSizeY = sampleScanSizeY;
            
            % Equalize the AFM images
            
            % Equalize the AFM images
            [equalizedTipStruct, equalizedSampleStruct, indicesTip, indicesSample] =...
                AFMBaseClass.equalizeAFMImages(tipStruct, sampleStruct,...
                'PreserveResolution',true,'SquarePhysicalPixel',false);
            [~, equalizedIndentationStruct, ~, indicesIndentation] =...
                AFMBaseClass.equalizeAFMImages(equalizedSampleStruct, indentationStruct,...
                'PreserveResolution',true,'SquarePhysicalPixel',false);
            
            if equalizedTipStruct.NumPixelsX*equalizedTipStruct.NumPixelsY >= 1024^2
                warning(sprintf('Trying to compute with maximal resolution results in a %i sized image. Downscaling instead...',sqrt(equalizedTipStruct.NumPixelsX*equalizedTipStruct.NumPixelsY)))
                [equalizedTipStruct, equalizedSampleStruct, indicesTip, indicesSample] =...
                    AFMBaseClass.equalizeAFMImages(tipStruct, sampleStruct,...
                    'PreserveResolution',false,'SquarePhysicalPixel',false);
                [~, equalizedIndentationStruct, ~, indicesIndentation] =...
                    AFMBaseClass.equalizeAFMImages(equalizedSampleStruct, indentationStruct,...
                    'PreserveResolution',false,'SquarePhysicalPixel',false);
            end
            
            tipHeightMapRescaled = equalizedTipStruct.Image;
            sampleHeightMapRescaled = equalizedSampleStruct.Image;
            indentationDepthMapRescaled = equalizedIndentationStruct.Image;
            
            % Prune tip height map background down to not affect high
            % sample map regions
            LowerCutFraction = 0.2;
            tipHeightMapRescaled(tipHeightMapRescaled < min(tipHeightMap,[],'all')...
                + LowerCutFraction*range(tipHeightMap,'all'))...
                = min(tipHeightMapRescaled - range(sampleHeightMapRescaled,'all'),[],'all');
            
            % Calculate padding sizes
            padX = floor(size(tipHeightMapRescaled, 2));
            padY = floor(size(tipHeightMapRescaled, 1));
            
            % Pad the sample height map and indentation depth map to accommodate shifts
            paddedSampleHeightMap = padarray(sampleHeightMapRescaled, [padY, padX], 'replicate', 'both');
            paddedIndentationDepthMap = padarray(indentationDepthMapRescaled, [padY, padX], 'replicate', 'both');
            
            % Initialize results
            contactAreas = zeros(sampleNumPixelsY, sampleNumPixelsX);
            
            % Calculate pixel area
            pixelArea = (equalizedSampleStruct.ScanSizeX / equalizedSampleStruct.NumPixelsX) * ...
                (equalizedSampleStruct.ScanSizeY / equalizedSampleStruct.NumPixelsY);
            
            % Find the position of the highest point in the tip height map
            [maxTipHeight, maxTipIndex] = max(tipHeightMapRescaled(:));
            [maxTipY, maxTipX] = ind2sub(size(tipHeightMapRescaled), maxTipIndex);
            
            % Invert Tip Data so, tip apex is at zero
            tipHeightMapRescaled = -tipHeightMapRescaled;
            tipHeightMapRescaled = tipHeightMapRescaled + maxTipHeight;
            
            % Iterate over each pixel in the sample height map
            for i = 1:sampleNumPixelsY
                for j = 1:sampleNumPixelsX
                    % Current position in the padded maps
                    currentX = round(j/sampleNumPixelsX*(indicesSample.endX - indicesSample.startX + 1)) + indicesSample.startX + padX;
                    currentY = round(i/sampleNumPixelsY*(indicesSample.endY - indicesSample.startY + 1)) + indicesSample.startY + padY;
                    
                    % Calculate the shift required to align the highest point of the tip with the current pixel
                    shiftY = currentY - maxTipY;
                    shiftX = currentX - maxTipX;
                    
                    % Determine the bounds of the region to be copied
                    startY = shiftY;
                    endY = shiftY + size(tipHeightMapRescaled, 1) - 1;
                    startX = shiftX;
                    endX = shiftX + size(tipHeightMapRescaled, 2) - 1;
                    
                    % Check if the start and end indices are within the padded map boundaries
                    if startY > 0 && endY <= size(paddedSampleHeightMap, 1) && startX > 0 && endX <= size(paddedSampleHeightMap, 2)
                        % Extract the region of interest from the padded maps
                        overlapRegionSample = paddedSampleHeightMap(startY:endY, startX:endX);
                        overlapRegionIndentation = paddedIndentationDepthMap(startY:endY, startX:endX);
                        
                        % Ensure the sizes match for shifting operation
                        if all(size(overlapRegionSample) == size(tipHeightMapRescaled))
                            % Calculate overlap
                            overlapMask = ((overlapRegionSample - sampleHeightMap(i,j)) >= (tipHeightMapRescaled - indentationDepthMap(i,j)));
                            
                            % Calculate contact area
                            contactAreas(i, j) = sum(overlapMask(:)) * pixelArea;
                            
                            if verbose && mod(j,256)==0
                                subplot(2,3,1)
                                imshow(overlapMask,[])
                                subplot(2,3,2)
                                imshow(contactAreas,[])
                                
                                subplot(2,3,4)
                                imshow(overlapRegionSample,[])
                                subplot(2,3,5)
                                imshowpair(overlapRegionSample,-tipHeightMapRescaled)
                                
                                subplot(2,3,[3 6])
                                mesh((overlapRegionSample - sampleHeightMap(i,j)))
                                hold on
                                mesh((tipHeightMapRescaled - indentationDepthMap(i,j)))
                                hold off
                                zlim([-3e-7 4e-7])
                                drawnow
                            end
                        else
                            fprintf('somethings wrong %i,%i \n',i,j)
                        end
                    end
                end
            end
            
            % Debugging and visualization
            if verbose
                figure;
                subplot(1, 3, 1);
                imagesc(sampleHeightMapRescaled);
                title('Sample Height Map');
                colorbar;
                
                subplot(1, 3, 2);
                imagesc(tipHeightMapRescaled);
                title('Rescaled Tip Height Map');
                colorbar;
                
                subplot(1, 3, 3);
                imagesc(contactAreas);
                title('Contact Areas');
                colorbar;
            end
        end
        
        function [outStruct1, outStruct2, indices1, indices2] = equalizeAFMImages(imageStruct1, imageStruct2, varargin)
            % Parse optional input arguments
            p = inputParser;
            addParameter(p, 'PreserveResolution', false, @(x) islogical(x) || x == 0 || x == 1);
            addParameter(p, 'SquarePhysicalPixel', false, @(x) islogical(x) || x == 0 || x == 1);
            parse(p, varargin{:});
            preserveResolution = p.Results.PreserveResolution;
            squarePhysicalPixel = p.Results.SquarePhysicalPixel;
            
            % Extract fields from input structs
            image1 = imageStruct1.Image;
            numPixelsX1 = imageStruct1.NumPixelsY;
            numPixelsY1 = imageStruct1.NumPixelsX;
            scanSizeX1 = imageStruct1.ScanSizeX;
            scanSizeY1 = imageStruct1.ScanSizeY;
            
            image2 = imageStruct2.Image;
            numPixelsX2 = imageStruct2.NumPixelsY;
            numPixelsY2 = imageStruct2.NumPixelsX;
            scanSizeX2 = imageStruct2.ScanSizeX;
            scanSizeY2 = imageStruct2.ScanSizeY;
            
            % Calculate physical sizes per pixel
            physSizePerPixelX1 = scanSizeX1 / numPixelsX1;
            physSizePerPixelY1 = scanSizeY1 / numPixelsY1;
            physSizePerPixelX2 = scanSizeX2 / numPixelsX2;
            physSizePerPixelY2 = scanSizeY2 / numPixelsY2;
            
            % Determine the maximum scan size
            [maxScanSizeX, ~] = max([scanSizeX1, scanSizeX2]);
            [maxScanSizeY, ~] = max([scanSizeY1, scanSizeY2]);
            
            % Determine the minimum physical size per pixel for upscaling if needed
            if preserveResolution
                [minPhysSizePerPixelX, ~] = min([physSizePerPixelX1, physSizePerPixelX2]);
                [minPhysSizePerPixelY, ~] = min([physSizePerPixelY1, physSizePerPixelY2]);
            else
                [minPhysSizePerPixelX, ~] = max([physSizePerPixelX1, physSizePerPixelX2]);
                [minPhysSizePerPixelY, ~] = max([physSizePerPixelY1, physSizePerPixelY2]);
            end
            
            % Ensure square physical pixels if requested
            if squarePhysicalPixel
                minPhysSizePerPixelX = min(minPhysSizePerPixelX, minPhysSizePerPixelY);
                minPhysSizePerPixelY = minPhysSizePerPixelX;
            end
            
            % Calculate new number of pixels
            newNumPixelsX = ceil(maxScanSizeX / minPhysSizePerPixelX);
            newNumPixelsY = ceil(maxScanSizeY / minPhysSizePerPixelY);
            
            % Resize and pad the images
            [resizedImage1, indices1] = AFMBaseClass.resizeAndPadImage(image1, ...
                numPixelsX1, numPixelsY1, newNumPixelsX, newNumPixelsY, ...
                scanSizeX1, scanSizeY1, minPhysSizePerPixelX, minPhysSizePerPixelY);
            [resizedImage2, indices2] = AFMBaseClass.resizeAndPadImage(image2, ...
                numPixelsX2, numPixelsY2, newNumPixelsX, newNumPixelsY, ...
                scanSizeX2, scanSizeY2, minPhysSizePerPixelX, minPhysSizePerPixelY);
            
            % Create output structs
            outStruct1 = imageStruct1;
            outStruct2 = imageStruct2;
            
            outStruct1.Image = resizedImage1;
            outStruct1.NumPixelsX = newNumPixelsY;
            outStruct1.NumPixelsY = newNumPixelsX;
            outStruct1.ScanSizeX = maxScanSizeX; 
            outStruct1.ScanSizeY = maxScanSizeY;
            
            outStruct2.Image = resizedImage2;
            outStruct2.NumPixelsX = newNumPixelsY;
            outStruct2.NumPixelsY = newNumPixelsX;
            outStruct2.ScanSizeX = maxScanSizeX; 
            outStruct2.ScanSizeY = maxScanSizeY;
        end
        
        function [resizedImage, indices] = resizeAndPadImage(image, ...
                numPixelsX, numPixelsY, newNumPixelsX, newNumPixelsY, ...
                scanSizeX, scanSizeY, physSizePerPixelX, physSizePerPixelY)
            % Calculate scale factors
            scaleFactorX = (scanSizeX / numPixelsX) / physSizePerPixelX;
            scaleFactorY = (scanSizeY / numPixelsY) / physSizePerPixelY;
            
            % Resize image
            resizedImage = imresize(image, [round(numPixelsY * scaleFactorY), round(numPixelsX * scaleFactorX)]);
            
            % Determine padding
            padX = newNumPixelsX - size(resizedImage, 2);
            padY = newNumPixelsY - size(resizedImage, 1);
            
            % Pad image
            padLeft = floor(padX / 2);
            padRight = padX - padLeft;
            padTop = floor(padY / 2);
            padBottom = padY - padTop;
            resizedImage = padarray(resizedImage, [padTop, padLeft], min(resizedImage, [], 'all'), 'pre');
            resizedImage = padarray(resizedImage, [padBottom, padRight], min(resizedImage, [], 'all'), 'post');
            
            % Calculate indices
            startIndexX = padLeft + 1;
            endIndexX = startIndexX + size(image, 2) * scaleFactorX - 1;
            startIndexY = padTop + 1;
            endIndexY = startIndexY + size(image, 1) * scaleFactorY - 1;
            indices = struct('startX', round(startIndexX), 'endX', round(endIndexX), 'startY', round(startIndexY), 'endY', round(endIndexY));
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
        
        function [OutSegment, Index] = find_matching_segment(SegmentElement, SegmentStruct)
            
            if isstruct(SegmentElement)
                TargetName = SegmentElement.Name;
                TargetSubName = SegmentElement.SubSegmentName;
                
                NameList = {SegmentStruct.Name};
                SubNameList = {SegmentStruct.SubSegmentName};
                
                Index = find(strcmp(TargetName, NameList) & strcmp(TargetSubName, SubNameList));
                
            elseif ischar(SegmentElement) || isstring(SegmentElement)
                TargetName = char(SegmentElement); % Convert to char if input is a string
                
                NameList = {SegmentStruct.Name};
                
                Index = find(strcmp(TargetName, NameList));
            else
                error('Invalid input: SegmentElement must be a struct, character, or string.');
            end
            
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
        
        function map_fiber_segment_properties_to_image_pixels(obj,PoolingMethod)
            
            if nargin < 2
                PoolingMethod = 'Median';
            end
            
            NumSegments = length(obj.Segment);
            
            FieldsList = {'Height','WidthHalfHeight',...
                'Prominence','WidthHalfProminence',...
                'SegmentLength','Mean_WidthHalfHeight',...
                'Median_WidthHalfHeight','Mean_Height',...
                'Median_Height','Mean_WidthHalfProminence',...
                'Median_WidthHalfProminence','Mean_Prominence',...
                'Median_Prominence','Area','Mean_Area','Median_Area',...
                'WidthBase','Mean_WidthBase','Median_WidthBase',...
                'AspectRatioHalfHeight','Mean_AspectRatioHalfHeight',...
                'Median_AspectRatioHalfHeight','AspectRatioBaseHeight',...
                'Mean_AspectRatioBaseHeight','Median_AspectRatioBaseHeight',...
                'AreaDerivedDiameter','Mean_AreaDerivedDiameter',...
                'Median_AreaDerivedDiameter',...
                'Ellipse_a','Mean_Ellipse_a','Median_Ellipse_a',...
                'Ellipse_b','Mean_Ellipse_b','Median_Ellipse_b',...
                'Ellipse_AspectRatio','Mean_Ellipse_AspectRatio',...
                'Median_Ellipse_AspectRatio',...
                'Ellipse_Area','Mean_Ellipse_Area','Median_Ellipse_Area',...
                'Ellipse_Height','Mean_Ellipse_Height',...
                'Median_Ellipse_Height',...
                'Ellipse_WidthHalfHeight','Mean_Ellipse_WidthHalfHeight',...
                'Median_Ellipse_WidthHalfHeight',...
                'RelativePixelPosition','RelativePosition'};
            
            switch PoolingMethod
                case 'Mean'
                    PoolingFcn = @(x)mean(x,'omitnan');
                case 'Median'
                    PoolingFcn = @(x)median(x,'omitnan');
                case 'Max'
                    PoolingFcn = @(x)max(x,'omitnan');
                case 'Min'
                    PoolingFcn = @(x)min(x,'omitnan');
            end
            
            h = waitbar(0,'setting up...',...
                    'Name',obj.Name);
            % initialize object properties
            obj.SegmentName = cell(obj.NumPixelsX*obj.NumPixelsY,1);
            obj.SubSegmentName = cell(obj.NumPixelsX*obj.NumPixelsY,1);
            obj.SubSegmentFullName = cell(obj.NumPixelsX*obj.NumPixelsY,1);
            [obj.SegmentName(:)] = {''};
            [obj.SubSegmentName(:)] = {''};
            [obj.SubSegmentFullName(:)] = {''};
            for i=1:length(FieldsList)
                for j=1:length(obj.Segment)
                    if ~isfield(obj.Segment(j),FieldsList{i}) || isempty(obj.Segment(j).(FieldsList{i}))
                        continue
                    elseif isnumeric(obj.Segment(j).(FieldsList{i}))
                        obj.(['FiberSegment_' FieldsList{i}]) = ones(obj.NumPixelsX*obj.NumPixelsY,1).*NaN;
                    else
                        obj.(['FiberSegment_' FieldsList{i}]) = cell(obj.NumPixelsX*obj.NumPixelsY,1);
                        [obj.(['FiberSegment_' FieldsList{i}])(:)] = {''};
                    end
                end
            end
            
            % skip non polyline segments
            for i=1:NumSegments
                waitbar(i/NumSegments,h,obj.Segment(i).Name,...
                    'Name',obj.Name)
                if ~isfield(obj.Segment(i),'SegmentPixelIndex') ||...
                        ~isfield(obj.Segment(i),'ROIObject') ||...
                        ~isfield(obj.Segment(i),'CorrespondingPixelIndex') ||...
                        isempty(obj.Segment(i).SegmentPixelIndex) ||...
                        isempty(obj.Segment(i).ROIObject) ||...
                        isempty(obj.Segment(i).CorrespondingPixelIndex) ||...
                        ~isequal(obj.Segment(i).Type,'polyline')
                    continue
                end
                CorrespondingPixelIndex = obj.Segment(i).CorrespondingPixelIndex;
                SegLength = size(CorrespondingPixelIndex,1);
                Processed = false(SegLength,1);
                while ~all(Processed)
                    Unprocessed = find(Processed == 0);
                    FirstUnprocessed = Unprocessed(1);
                    CurrentPixel = CorrespondingPixelIndex(FirstUnprocessed,:);
                    Indizes = find(CorrespondingPixelIndex(:,1) == CurrentPixel(1) &...
                        CorrespondingPixelIndex(:,2) == CurrentPixel(2));
                    if CurrentPixel(1) < 1 ||...
                            CurrentPixel(1) > obj.NumPixelsY ||...
                            CurrentPixel(2) < 1 ||...
                            CurrentPixel(2) > obj.NumPixelsX
                        Processed(Indizes) = true;
                        continue
                    end
                    % Pool data by PoolingMethod and assign to
                    % corresponding image pixel
                    for j=1:length(FieldsList)
                        if length(obj.Segment(i).(FieldsList{j})) == SegLength
                            Values = obj.Segment(i).(FieldsList{j})(Indizes);
                            PooledValue = PoolingFcn(Values);
                        elseif size(obj.Segment(i).(FieldsList{j})) == size([1])
                            PooledValue = obj.Segment(i).(FieldsList{j});
                        else
                            PooledValue = NaN;
                        end
                        try
                        obj.(['FiberSegment_' FieldsList{j}])...
                            (obj.Map2List(CurrentPixel(2),CurrentPixel(1))) = PooledValue;
                        catch
                            warning("Couldn't assign property")
                            i
                            j
                        end
                    end
                    obj.SegmentName{obj.Map2List(CurrentPixel(2),CurrentPixel(1))} =...
                        obj.Segment(i).Name;
                    obj.SubSegmentName{obj.Map2List(CurrentPixel(2),CurrentPixel(1))} =...
                        obj.Segment(i).SubSegmentName;
                    obj.SubSegmentFullName{obj.Map2List(CurrentPixel(2),CurrentPixel(1))} =...
                        [obj.Segment(i).Name '-' obj.Segment(i).SubSegmentName];
                    Processed(Indizes) = true;
                end
            end
            close(h)
        end
        
        function map_segments_to_image_pixels(obj)
            
            h = waitbar(0,'setting up...',...
                'Name',obj.Name);
            
            
            F = figure;
            
            ROIObjects = obj.draw_segments_to_figure(F);
            
            % initialize object properties
            obj.SegmentName = cell(obj.NumPixelsX*obj.NumPixelsY,1);
            obj.SubSegmentName = cell(obj.NumPixelsX*obj.NumPixelsY,1);
            obj.SubSegmentFullName = cell(obj.NumPixelsX*obj.NumPixelsY,1);
            [obj.SegmentName(:)] = {''};
            [obj.SubSegmentName(:)] = {''};
            [obj.SubSegmentFullName(:)] = {''};
            
            NumSegments = numel(obj.Segment);
            
            UniqueMask = zeros(obj.NumPixelsX,obj.NumPixelsY);
            UnionMask = zeros(obj.NumPixelsX,obj.NumPixelsY);
            
            for i=1:NumSegments
                waitbar(i/NumSegments,h,obj.Segment(i).Name,...
                    'Name',obj.Name)
                if ~isfield(obj.Segment(i),'ROIObject') ||...
                        isempty(obj.Segment(i).ROIObject)
                    continue
                end
                CurMask = ROIObjects{i}.createMask;
%                 CurList = obj.convert_map_to_data_list(CurMask);
%                 [obj.UnionSegmentName(CurList)] = {obj.Segment(i).Name};
                UniqueMask = xor(CurMask,UniqueMask);
%                 [obj.UniqueSegmentName(CurList)] = {obj.Segment(i).Name};
            end
            for i=1:NumSegments
                CurMask = ROIObjects{i}.createMask;
                UniqueList = obj.convert_map_to_data_list(UniqueMask & CurMask);
                [obj.UniqueSegmentName(CurList)] = {obj.Segment(i).Name};
            end
            close(h)
            
        end
        
        function [DynPropNames,ChannelNames] = write_unrolled_channels_to_dynamic_properties(obj)
            % unrolls all images in Channel struct into nx1 lists and
            % assigns them to dynamic properties derived from the Channel
            % name
            
            ChannelNames = {obj.Channel.Name};
            DynPropNames = genvarname(ChannelNames);
            
            for i=1:length(DynPropNames)
                if obj.NumPixelsX*obj.NumPixelsY ~=...
                        obj.Channel(i).NumPixelsX*obj.Channel(i).NumPixelsY
                    continue
                end
                if ~isprop(obj,DynPropNames{i})
                    mp = obj.addprop(DynPropNames{i});
                    mp.Description = 'DynamicProperty';
                else
                   mp = findprop(obj,DynPropNames{i});
                   if ~isequal(mp.Description,'DynamicProperty')
                       warning("You're trying to overwrite a static class property. Try renaming your Channels. Skipping this one...")
                       continue
                   end
                end
                List = obj.convert_map_to_data_list(obj.Channel(i).Image);
                obj.(DynPropNames{i}) = List;
            end
            
        end
        
        function delete_all_dynamic_properties(obj)
            
            PropertyNames = fieldnames(obj);
            for i=1:length(PropertyNames)
                mp = findprop(obj,PropertyNames{i});
                if isequal(mp.Description,'DynamicProperty')
                    delete(mp);
                end
            end
        end
        
        function [ROIObjects,F] = draw_segments_to_figure(obj,F)
            
            if nargin < 2
                F = figure;
            end
            
            NumSegments = numel(obj.Segment);
            ROIObjects = cell(NumSegments,1);
            
            Image = zeros(obj.NumPixelsX,obj.NumPixelsY);
            
            gcf = F;
            imshow(Image)
            
            for i=1:NumSegments
                if ~isfield(obj.Segment(i),'ROIObject') ||...
                        isempty(obj.Segment(i).ROIObject)
                    continue
                end
                CurrentDrawMode = obj.Segment(i).Type;
                switch CurrentDrawMode
                    case 'line'
                        ROIObjects{i} = drawline('Position',obj.Segment(i).ROIObject.Position,...
                            'Deletable',1,...
                            'InteractionsAllowed','all',...
                            'LineWidth',obj.Segment(i).ROIObject.LineWidth,...
                            'Label',sprintf('%s || %s',obj.Segment(i).Name,obj.Segment(i).SubSegmentName),...
                            'LabelAlpha',0.6);
                    case 'freehand'
                        ROIObjects{i} = drawfreehand('Position',obj.Segment(i).ROIObject.Position,...
                            'Deletable',1,...
                            'InteractionsAllowed','all',...
                            'LineWidth',obj.Segment(i).ROIObject.LineWidth,...
                            'Label',sprintf('%s || %s',obj.Segment(i).Name,obj.Segment(i).SubSegmentName),...
                            'LabelAlpha',0.6);
                    case 'square'
                        ROIObjects{i} = drawsquare('Position',obj.Segment(i).ROIObject.Position,...
                            'Deletable',1,...
                            'InteractionsAllowed','all',...
                            'LineWidth',obj.Segment(i).ROIObject.LineWidth,...
                            'Label',sprintf('%s || %s',obj.Segment(i).Name,obj.Segment(i).SubSegmentName),...
                            'LabelAlpha',0.6);
                    case 'circle'
                        ROIObjects{i} = drawcircle('Center',obj.Segment(i).ROIObject.Center,...
                            'Radius',obj.Segment(i).ROIObject.Radius,...
                            'Deletable',1,...
                            'InteractionsAllowed','all',...
                            'LineWidth',obj.Segment(i).ROIObject.LineWidth,...
                            'Label',sprintf('%s || %s',obj.Segment(i).Name,obj.Segment(i).SubSegmentName),...
                            'LabelAlpha',0.6);
                    case 'ellipse'
                        ROIObjects{i} = drawellipse('Center',obj.Segment(i).ROIObject.Center,...
                            'SemiAxes',obj.Segment(i).ROIObject.SemiAxes,...
                            'RotationAngle',obj.Segment(i).ROIObject.RotationAngle,...
                            'Deletable',1,...
                            'InteractionsAllowed','all',...
                            'LineWidth',obj.Segment(i).ROIObject.LineWidth,...
                            'Label',sprintf('%s || %s',obj.Segment(i).Name,obj.Segment(i).SubSegmentName),...
                            'LabelAlpha',0.6);
                    case 'polygon'
                        ROIObjects{i} = drawpolygon('Position',obj.Segment(i).ROIObject.Position,...
                            'Deletable',1,...
                            'InteractionsAllowed','all',...
                            'LineWidth',obj.Segment(i).ROIObject.LineWidth,...
                            'Label',sprintf('%s || %s',obj.Segment(i).Name,obj.Segment(i).SubSegmentName),...
                            'LabelAlpha',0.6);
                    case 'rectangle'
                        ROIObjects{i} = drawrectangle('Position',obj.Segment(i).ROIObject.Position,...
                            'RotationAngle',obj.Segment(i).ROIObject.RotationAngle,...
                            'Rotatable',true,...
                            'Deletable',1,...
                            'InteractionsAllowed','all',...
                            'LineWidth',obj.Segment(i).ROIObject.LineWidth,...
                            'Label',sprintf('%s || %s',obj.Segment(i).Name,obj.Segment(i).SubSegmentName),...
                            'LabelAlpha',0.6);
                    case 'crosshair'
                        ROIObjects{i} = drawcrosshair('Position',obj.Segment(i).ROIObject.Position,...
                            'Deletable',1,...
                            'InteractionsAllowed','all',...
                            'LineWidth',obj.Segment(i).ROIObject.LineWidth,...
                            'Label',sprintf('%s || %s',obj.Segment(i).Name,obj.Segment(i).SubSegmentName),...
                            'LabelAlpha',0.6);
                    case 'point'
                        ROIObjects{i} = drawpoint('Position',obj.Segment(i).ROIObject.Position,...
                            'Deletable',1,...
                            'InteractionsAllowed','all',...
                            'LineWidth',obj.Segment(i).ROIObject.LineWidth,...
                            'Label',sprintf('%s || %s',obj.Segment(i).Name,obj.Segment(i).SubSegmentName),...
                            'LabelAlpha',0.6);
                    case 'assisted'
                        ROIObjects{i} = drawassisted('Position',obj.Segment(i).ROIObject.Position,...
                            'Deletable',1,...
                            'InteractionsAllowed','all',...
                            'LineWidth',obj.Segment(i).ROIObject.LineWidth,...
                            'Label',sprintf('%s || %s',obj.Segment(i).Name,obj.Segment(i).SubSegmentName),...
                            'LabelAlpha',0.6);
                    case 'polyline'
                        ROIObjects{i} = drawpolyline('Position',obj.Segment(i).ROIObject.Position,...
                            'Deletable',1,...
                            'InteractionsAllowed','all',...
                            'LineWidth',obj.Segment(i).ROIObject.LineWidth,...
                            'Label',sprintf('%s || %s',obj.Segment(i).Name,obj.Segment(i).SubSegmentName),...
                            'LabelAlpha',0.6);
                end
            end
            
        end
        
        function visualize_curvature(obj,varargin)
            % function visualize_curvature(obj,varargin)
            %
            % <FUNCTION DESCRIPTION HERE>
            %
            %
            % Required inputs
            % obj ... <VARIABLE DESCRIPTION>
            %
            % Optional inputs
            % ZSpacingFactor ... <OPTIONAL POSITIONAL VARIABLE DESCRIPTION>
            % HeightChannelName ... <OPTIONAL POSITIONAL VARIABLE DESCRIPTION>
            % SurfFitChannelKeyphrase ... <OPTIONAL POSITIONAL VARIABLE DESCRIPTION>
            
            p = inputParser;
            p.FunctionName = "visualize_curvature";
            p.CaseSensitive = false;
            p.PartialMatching = true;
            
            % Required inputs
            validobj = @(x)true;
            addRequired(p,"obj",validobj);
            
            % Optional positional inputs
            defaultZSpacingFactor = 5;
            defaultHeightChannelName = 'Processed';
            defaultSurfFitChannelKeyphrase = '';
            validZSpacingFactor = @(x)true;
            validHeightChannelName = @(x)true;
            validSurfFitChannelKeyphrase = @(x)true;
            addParameter(p,"ZSpacingFactor",defaultZSpacingFactor,validZSpacingFactor);
            addParameter(p,"HeightChannelName",defaultHeightChannelName,validHeightChannelName);
            addParameter(p,"SurfFitChannelKeyphrase",defaultSurfFitChannelKeyphrase,validSurfFitChannelKeyphrase);
            
            parse(p,obj);
            
            % Assign parsing results to named variables
            obj = p.Results.obj;
            ZSpacingFactor = p.Results.ZSpacingFactor;
            HeightChannelName = p.Results.HeightChannelName;
            SurfFitChannelKeyphrase = p.Results.SurfFitChannelKeyphrase;
            
            
            
            heightChannel = obj.get_channel(HeightChannelName);
            
            % Determine which SurfFit channels should be displayed if there
            % are multiple.
            PossibleChannels = contains({obj.Channel.Name},'Radius of Curvature') & contains({obj.Channel.Name},SurfFitChannelKeyphrase);
            if ~any(PossibleChannels)
                PossibleChannels = contains({obj.Channel.Name},'Radius of Curvature');
                if ~any(PossibleChannels)
                    error('No surface topography channels found');
                end
            end
            ChannelIndex = find(PossibleChannels);
            radiusOfCurvatureChannel = obj.Channel(ChannelIndex(1));
            
            PossibleChannels = contains({obj.Channel.Name},'Slope Kernel') | contains({obj.Channel.Name},'Kernel Slope')...
                & contains({obj.Channel.Name},SurfFitChannelKeyphrase);
            if ~any(PossibleChannels)
                PossibleChannels =  contains({obj.Channel.Name},'Slope Kernel') | contains({obj.Channel.Name},'Kernel Slope');
                if ~any(PossibleChannels)
                    error('No surface topography channels found');
                end
            end
            ChannelIndex = find(PossibleChannels);
            slopeChannel = obj.Channel(ChannelIndex(1));
            
            PossibleChannels = contains({obj.Channel.Name},'Slope Direction') & contains({obj.Channel.Name},SurfFitChannelKeyphrase);
            if ~any(PossibleChannels)
                PossibleChannels = contains({obj.Channel.Name},'Slope Direction');
                if ~any(PossibleChannels)
                    error('No surface topography channels found');
                end
            end
            ChannelIndex = find(PossibleChannels);
            slopeDirectionChannel = obj.Channel(ChannelIndex(1));
            
            
            
            if nargin < 5
                ZSpacingFactor = 1;
            end
            % Create the figure and set up the layout
            fig = figure;
            ax1 = subplot(1, 2, 1);
            ax2 = subplot(1, 2, 2);
            
            % Show the height image in the bottom part of the figure
            imshow(flipud(heightChannel.Image),[], 'Parent', ax2);
            
            colormap(ax2,AFMImage.define_afm_color_map);
            
            % Initialize the drawpoint object
            hPoint = drawpoint(ax2, 'Position', [heightChannel.NumPixelsX/2, heightChannel.NumPixelsY/2], ...
                'Color', 'r', 'Selected', true);
            
            % Set up the callback function for the drawpoint object
            addlistener(hPoint, 'ROIMoved', @pointMoved);
            
            infoTextBox = uicontrol('Style', 'text', 'Position', [10, 10, 300, 80], 'BackgroundColor', 'white', 'HorizontalAlignment', 'left');
            
            % Function to handle the drawpoint object's position change
            function pointMoved(src, eventData)
                
                position = hPoint.Position;
                
                % Get the position in the height image
                xIdx = round(position(1));
                yIdx = heightChannel.NumPixelsY - round(position(2)) + 1;
                
                % Get the values from the channels
                height = heightChannel.Image(yIdx, xIdx);
                radius = radiusOfCurvatureChannel.Image(yIdx, xIdx);
                slope = slopeChannel.Image(yIdx, xIdx);
                slopeDirection = slopeDirectionChannel.Image(yIdx, xIdx);
                [heightMultiplier, heightUnit, ~] = AFMImage.parse_unit_scale(height, 'm', 1);
                [radiusMultiplier, radiusUnit, ~] = AFMImage.parse_unit_scale(radius, 'm', 1);
                [slopeMultiplier, slopeUnit, ~] = AFMImage.parse_unit_scale(slope, 'm/m', 1);
                [slopeDirectionMultiplier, slopeDirectionUnit, ~] = AFMImage.parse_unit_scale(slopeDirection, '°', 1);
                
                
                % Calculate the tangent plane
                [x, y] = meshgrid(linspace(0, heightChannel.ScanSizeX, heightChannel.NumPixelsX), linspace(0, heightChannel.ScanSizeY, heightChannel.NumPixelsY));
                z = height + slope * (cosd(slopeDirection) * (x - x(yIdx, xIdx)) + sind(slopeDirection) * (y - y(yIdx, xIdx)));
                
                % Calculate the sphere center
                sphereCenter = [x(yIdx, xIdx), y(yIdx, xIdx), height + sign(radius) * abs(radius)];
                
                % Set up the surface plot in the top part of the figure
                cla(ax1);
                hold(ax1, 'on');
                surf(ax1, x, y, heightChannel.Image, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
                colormap(ax1, AFMImage.define_afm_color_map);
                
                % Add this line to set the color scaling of the height map plot
                caxis(ax1, [min(heightChannel.Image,[],'all'), max(heightChannel.Image,[],'all')]);
                
                num_points = 8;
                [x_down, y_down] = meshgrid(linspace(0, heightChannel.ScanSizeX, num_points), linspace(0, heightChannel.ScanSizeY, num_points));
                z_down = height + slope * (cosd(slopeDirection) * (x_down - x(yIdx, xIdx)) + sind(slopeDirection) * (y_down - y(yIdx, xIdx)));
                surf(ax1, x_down, y_down, z_down, 'EdgeColor', 'k', 'FaceColor', 'g', 'FaceAlpha', 0.3);
                
                
                % Plot the sphere
                [sx, sy, sz] = sphere(50);
                sx = sx * abs(radius) + sphereCenter(1);
                sy = sy * abs(radius) + sphereCenter(2);
                sz = sz * abs(radius) + sphereCenter(3);
                
                % Rotate the sphere about the chosen point
                angleToRotate = atand(slope);
                rotationAxis = [-sind(slopeDirection), cosd(slopeDirection), 0];
                rotationCenter = [x(yIdx, xIdx), y(yIdx, xIdx), heightChannel.Image(yIdx, xIdx)];
                [sx, sy, sz] = rotateSphere(sx, sy, sz, rotationCenter, rotationAxis, angleToRotate);
                
                surf(ax1, sx, sy, sz, 'EdgeColor', 'k', 'FaceColor', 'r', 'FaceAlpha', 0.3);
                
                % Adjust the axis settings
                axis(ax1, 'equal');
                xlim(ax1, [0 heightChannel.ScanSizeX]);
                ylim(ax1, [0 heightChannel.ScanSizeY]);
                Spacer = range(heightChannel.Image,'all').*ZSpacingFactor;
                zlim(ax1, [min(heightChannel.Image,[],'all')-Spacer max(heightChannel.Image,[],'all')+Spacer])
                hold(ax1, 'off');
                view(ax1, 3);
                xlabel(ax1, 'x');
                ylabel(ax1, 'y');
                zlabel(ax1, 'z');
                title(ax1, 'Height map, tangent plane, and sphere');
                infoText = sprintf('Height: %.3f %s\nRadius: %.3f %s\nSlope: %.3f %s\nSlope Direction: %.3f %s', ...
                    height * heightMultiplier, heightUnit, ...
                    radius * radiusMultiplier, radiusUnit, ...
                    slope * slopeMultiplier, slopeUnit, ...
                    slopeDirection * slopeDirectionMultiplier, slopeDirectionUnit);
                
                set(infoTextBox, 'String', infoText);
                
            end
            % Call the pointMoved function to update the visualization initially
            pointMoved();
            
        end
        
        
    end
end