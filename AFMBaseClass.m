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
        Segment
        OverlayGroup
    end
    
    methods
        % Main Methods
        
        function obj = AFMBaseClass()
            
            obj.OverlayGroup.hasOverlayGroup = false;
        end
        
        function choose_segments_manually(obj,SegmentType)
            
            if nargin < 2
                SegmentType = 'Apexline';
            end
            
            h.Fig = figure('Color','w');
            h.Fig.WindowState = 'fullscreen';
            
            h.c(1) = uicontrol(h.Fig,'style','pushbutton','units','normalized',...
                'position',[.8 .25 .14 .05],'string','Accept Current Segment', 'FontSize',22,...
                'Callback',@accept_current_segment);
            
            h.c(2) = uicontrol(h.Fig,'style','pushbutton','units','normalized',...
                'position',[.8 .2 .14 .05],'string','Delete Selected Segment', 'FontSize',22,...
                'Callback',@delete_selected);
            
            h.c(3) = uicontrol(h.Fig,'style','text','units','normalized',...
                'position',[.8 .9 .15 .05],'string','List of Segments',...
                'FontSize',22);
            
            h.c(4) = uicontrol(h.Fig,'style','pushbutton','units','normalized',...
                'position',[.8 .2 .14 .05],'string','New Segment', 'FontSize',22,...
                'Callback',@delete_selected);
            
            h.ListBox = uicontrol(h.Fig,...
                'Style','listbox',...
                'Max',1000000,'Min',1,...
                'Units','normalized',...
                'Position',[.8  .5 .15 .4]);
            
            Channel = obj.get_channel('Processed');
            axes(h.Fig,'Position',[0.1 0.1 .6 .8]);
            [h.Multiplier,h.Unit,~] = AFMImage.parse_unit_scale(range(Channel.Image,'all'),Channel.Unit,1);
            h.I = imshow(Channel.Image*h.Multiplier,[],'Colormap',obj.CMap);
            h.I.ButtonDownFcn = @get_and_draw_profile;
            hold on
            c = colorbar;
            c.FontSize = round(24);
            c.Label.String = sprintf('%s [%s]',Channel.Name,h.Unit);
            c.Label.FontSize = round(24);
            
            
            function delete_selected(varargin)
                
                OldString = get(h.ListBox,'String');
                DeleteIdx = get(h.ListBox,'Value');
                OldString(DeleteIdx) = [];
                OutStruct(Index).FullFile(DeleteIdx) = [];
                
                set(h.ListBox,'Value',1); 
                set(h.ListBox,'String',OldString);
            end
            
            uiwait(h.Fig)
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
    
end