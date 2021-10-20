classdef AFMImage < matlab.mixin.Copyable & matlab.mixin.SetGet & handle & AFMBaseClass
    % This is supposed to be a class for analysis and processing of AFM
    % data at the image level in general and isn't restricted to a specific
    % mode of image acquisition. As such, it should be able to load, access
    % and process images from all AC-modes, contact mode, QI-images, image
    % projections of ForceMaps and so on (maybe subclass SurfacePotentialMap
    % should eventually be incorporated into this).
   
    properties
        % File and Header properties
        ImagingType
        FileVersion
        DateTime
        NumChannels
    end
    properties
        % Image Data Properties. All dimensions in SI-units, Angles in
        % degrees
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
        % Properties related to Image processing/segmenting/classification
        MaskBackground
    end
    properties
        ErodedTip
        DepthDependendTipRadius
        ProjectedTipArea
    end
    properties
        % All the Flags
        hasSensitivity
        hasSpringConstant
        hasHeight
        hasHeightMeasured
        hasErrorSignal
        hasLateralDeflection
        hasLockInAmplitude
        hasLockInPhase
        hasVerticalDeflection
        hasProcessed
        hasBackgroundMask
        hasOverlay
        hasDeconvolutedCantileverTip
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
            if nargin == 1
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
            
            obj.construct_list_to_map_relations();
            
            obj.preprocess_image
            
        end
        
        function preprocess_image(obj)
            
            Height = obj.get_unprocessed_height_channel('Height (Trace)');
            
            if isempty(Height)
                return
            end
            
            if size(Height.Image,1) < 128
                Map = imresize(Height.Image,[256 256],'nearest');
            elseif size(Height.Image,1) < 512
                Map = imresize(Height.Image,[512 512],'nearest');
            else
                Map = imresize(Height.Image,[1024 1024],'nearest');
            end
            for i=1:5
                Map = AFMImage.subtract_line_fit_vertical_rov(Map,.2,0);
            end
            Map = imresize(Map,[obj.NumPixelsY obj.NumPixelsX],'nearest');
            
            Map = AFMImage.find_and_replace_outlier_lines(Map,10);
            
            % write to Channel
            obj.delete_channel('Processed')
            Processed = obj.create_standard_channel(Map,'Processed','m');
            
            [Channel,Index] = obj.get_channel('Processed');
            if isempty(Channel)
                Len = length(obj.Channel);
                if ~Len
                    obj.Channel = Processed;
                else
                    obj.Channel(Len+1) = Processed;
                end
            else
                obj.Channel(Index) = Processed;
            end
            
            obj.create_pixel_difference_channel;
            
        end
        
        function deconvolute_cantilever_tip(obj,StepSize)
            
            Channel = obj.get_unprocessed_height_channel('Height (measured) (Trace)');
            Based = imgaussfilt(AFMImage.subtract_line_fit_hist(Channel.Image,0.4));
            Channel.Image = Based;
            obj.Channel(end+1) = Channel;
            obj.Channel(end).Name = 'Background Mask';
            obj.Channel(end).Unit = 'Logical';
            obj.Channel(end).Image = obj.mask_background_by_threshold(Channel.Image,1);
            Channel.Image = obj.masked_plane_fit(Channel,obj.Channel(end).Image);
            ConeHeight = range(Channel.Image,'all');
            Cone = obj.cone(obj.NumPixelsX,obj.NumPixelsY,ConeHeight,obj.ScanSizeX,obj.ScanSizeY,10e-9);
            obj.Channel(end+1) = Channel;
            obj.Channel(end).Name = 'Eroded Tip';
            obj.Channel(end).Unit = 'm';
            obj.Channel(end).Image = obj.deconvolute_by_mathematical_morphology(Channel.Image,Cone);
            
            if nargin < 2
                StepSize = 1e-9;
            end
            
            PixelSizeX = obj.ScanSizeX/obj.NumPixelsX;
            PixelSizeY = obj.ScanSizeY/obj.NumPixelsY;
            
            h = waitbar(1/3,'Calculating projected tip area');
            
            obj.ProjectedTipArea = obj.calculate_projected_area(obj.Channel(end).Image,PixelSizeX,PixelSizeY,StepSize);
            
            waitbar(2/3,h,'Calculating depth dependent tip radius');
            
            obj.DepthDependendTipRadius = obj.calculate_depth_dependend_tip_data(obj.ProjectedTipArea,75);
            
            obj.hasDeconvolutedCantileverTip = true;
            close(h); 
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
        
        function base_image()
            
        end
        
        function subtract_overlayed_image(obj,OverlayedImageClassInstance)
            
            if ~obj.hasOverlay
                warning('Overlay missing from your AFMImage. Determine Overlay pairings and parameters first')
                return
            end
            
            Channel1 = obj.get_channel('Overlay');
            Channel2 = OverlayedImageClassInstance.get_channel('Overlay');
            
            OutChannel = Channel1;
            OutChannel.Name = 'Overlay Difference';
            OutChannel.Image = Channel1.Image - Channel2.Image;
            [Old,OldIndex] = obj.get_channel('Overlay Difference');
            if isempty(Old)
                obj.Channel(end+1) = OutChannel;
            else
                obj.Channel(OldIndex) = OutChannel;
            end
        end
        
        function find_and_classify_fibrils(obj,RadiusMeters,NumRotations,MinKernelSize)
            
            Channel = obj.get_channel('Processed');
            
            [ActivationImage,AngleImage] = AFMImage.mask_apex_points_by_convolution(Channel,RadiusMeters,NumRotations,MinKernelSize);
            
            
            
        end
        
        function find_and_classify_fibrils_old(obj, KernelWindowX,...
                KernelWindowY, KernelStepSizeX, KernelStepSizeY,ComboSumThresh,DebugBool)
           
            [Processed,~,isProcessed] = obj.get_unprocessed_height_channel('Processed');
            if ~isProcessed
                Processed.Name = 'Processed';
                for i=1:5
                    Processed.Image = AFMImage.subtract_line_fit_vertical_rov(Processed.Image,.2,0);
                end
                obj.Channel(end+1) = Processed;
            end
            BackgroundMask = obj.get_channel('Background Mask');
            if isempty(BackgroundMask)
                BackgroundMask = Processed;
                BackgroundMask.Image = AFMImage.mask_background_by_threshold(Processed.Image,5);
                BackgroundMask.Name = 'Background Mask';
                BackgroundMask.Unit = 'Boolean';
                obj.Channel(end+1) = BackgroundMask;
            end
            ApexMap = obj.get_channel('Apex Map');
            if isempty(ApexMap)
                ApexMap = Processed;
                ApexMap.Image = AFMImage.mask_apex_points_by_rotation(Processed.Image,10,BackgroundMask.Image);
                ApexMap.Name = 'Apex Map';
                ApexMap.Unit = 'Apex Confidence';
                obj.Channel(end+1) = ApexMap;
            end
            
            % Convert window- and stepsizes to pixels
            KernelWindowX = round(KernelWindowX/Processed.ScanSizeX*Processed.NumPixelsX);
            KernelWindowY = round(KernelWindowY/Processed.ScanSizeY*Processed.NumPixelsY);
            KernelStepSizeX = round(KernelStepSizeX/Processed.ScanSizeX*Processed.NumPixelsX);
            KernelStepSizeY = round(KernelStepSizeY/Processed.ScanSizeY*Processed.NumPixelsY);
            if ~mod(KernelWindowX,2)
                KernelWindowX = KernelWindowX + 1;
            end
            if ~mod(KernelWindowY,2)
                KernelWindowY = KernelWindowY + 1;
            end
            
            % scan over the image with the feature kernel, that snaps onto
            % fibrils, if it finds regions with enough cumsum of
            % apexpoints. Once the kernel has scanned over all of the
            % image, the outer while loop terminates
            ObjectIndex = 1;
            CurOuterPosX = (KernelWindowX+1)/2;
            CurOuterPosY = (KernelWindowY+1)/2;
            OuterCondition = true;
            if DebugBool
                Image = imshow(Processed.Image,[]);
                Rect = drawrectangle('Position',[(CurOuterPosX - (KernelWindowX-1)/2) ...
                    (CurOuterPosY - (KernelWindowY-1)/2) ...
                    KernelWindowX, KernelWindowY],...
                    'RotationAngle',0);
            end
            while OuterCondition
                [ComboSum,AspectRatio,MeanX,MeanY,Angle] = ...
                    AFMImage.fibril_feature_kernel(...
                    Processed,ApexMap,CurOuterPosX,CurOuterPosY,...
                    KernelWindowX,KernelWindowY);
                if ComboSum > ComboSumThresh
                    ReachedImageBorder = false;
                    ObjectEnded = false;
                    CurInnerPosX = CurOuterPosX;
                    CurInnerPosY = CurOuterPosY;
                    Angle = 0;
                    AdaptiveWindowX = KernelWindowX;
                    AdaptiveWindowY = KernelWindowY;
                    while ~ReachedImageBorder || ObjectEnded
                        [ComboSum,AspectRatio,AdaptiveWindowX,AdaptiveWindowY,Angle] = ...
                            AFMImage.fibril_feature_kernel(...
                            Processed,ApexMap,CurInnerPosX,CurInnerPosY,Angle,...
                            AdaptiveWindowX,AdaptiveWindowY);
                    end
                    ObjectIndex = ObjectIndex + 1;
                end
                if CurOuterPosX+(KernelWindowX+1)/2 == Processed.NumPixelsX
                    if CurOuterPosY+(KernelWindowY+1)/2 == Processed.NumPixelsY
                        OuterCondition = false;
                    elseif CurOuterPosY+KernelStepSizeY+(KernelWindowY-1)/2 > Processed.NumPixelsY
                        CurOuterPosX = (KernelWindowX+1)/2;
                        CurOuterPosY = Processed.NumPixelsY - (KernelWindowY+1)/2;
                    else
                        CurOuterPosX = (KernelWindowX+1)/2;
                        CurOuterPosY = CurOuterPosY + KernelStepSizeY;
                    end
                elseif CurOuterPosX+KernelStepSizeX+(KernelWindowX-1)/2 > Processed.NumPixelsX
                    CurOuterPosX = Processed.NumPixelsX - (KernelWindowX+1)/2;
                else
                    CurOuterPosX = CurOuterPosX + KernelStepSizeX;
                end
                if DebugBool
                    Rect.Position = [(CurOuterPosX - (KernelWindowX-1)/2) ...
                        (CurOuterPosY - (KernelWindowY-1)/2) ...
                        KernelWindowX, KernelWindowY];
                    drawnow
                end
            end
            
        end
        
    end
    
    methods(Static)
        % Static Main Methods
        
        function ChannelOut = overlay_parameters_by_bayesopt(Channel1,Channel2,BackgroundPercent,MinOverlap,AngleRange,UseParallel,MaxFunEval,PreMaxFunEval,NumPreSearches,NClusters)
            % ChannelOut = overlay_parameters_by_bayesopt(Channel1,Channel2,BackgroundPercent,MinOverlap,AngleRange,UseParallel,MaxFunEval,PreMaxFunEval,NumPreSearches,NClusters)
            
            % Assume Channel2 has same angle and origin as Channel1
            Channel2.OriginX = Channel1.OriginX;
            Channel2.OriginY = Channel1.OriginY;
            Channel2.ScanAngle = Channel1.ScanAngle;
            
            % Resize the bigger Channel (in ScanSize-per-Pixel) so
            % imagesizes correspond to ScanSizes. Do nothing, if they are
            % exactly the same
            SizePerPixel1 = Channel1.ScanSizeX/Channel1.NumPixelsX;
            SizePerPixel2 = Channel2.ScanSizeX/Channel2.NumPixelsX;
            if SizePerPixel1 > SizePerPixel2
                ScaleMultiplier = SizePerPixel1/SizePerPixel2;
                Channel1.Image = imresize(Channel1.Image,ScaleMultiplier);
                Channel1.NumPixelsX = size(Channel1.Image,1);
                Channel1.NumPixelsY = size(Channel1.Image,2);
            elseif SizePerPixel1 > SizePerPixel2
                ScaleMultiplier = SizePerPixel2/SizePerPixel1;
                Channel2.Image = imresize(Channel2.Image,ScaleMultiplier);
                Channel2.NumPixelsX = size(Channel2.Image,1);
                Channel2.NumPixelsY = size(Channel2.Image,2);
            end
            % Create Background masks
            Mask1 = AFMImage.mask_background_by_threshold(Channel1.Image,BackgroundPercent,'on');
            Mask2 = AFMImage.mask_background_by_threshold(Channel2.Image,BackgroundPercent,'on');
            
            
            ShiftPixX = optimizableVariable(...
                'ShiftPixX',[-floor(Channel1.NumPixelsX/2 + Channel2.NumPixelsX*(1/2-MinOverlap)),...
                floor(Channel1.NumPixelsX/2 + Channel2.NumPixelsX*(1/2-MinOverlap))],...
                'Type','integer');
            ShiftPixY = optimizableVariable(...
                'ShiftPixY',[-floor(Channel1.NumPixelsY/2 + Channel2.NumPixelsY*(1/2-MinOverlap)),...
                floor(Channel1.NumPixelsY/2 + Channel2.NumPixelsY*(1/2-MinOverlap))],...
                'Type','integer');
            Angle = optimizableVariable('Angle',[-AngleRange/2,AngleRange/2],'Type','real');
            fun = @(x)AFMImage.overlay_loss(Channel1,Channel2,~Mask1,~Mask2,...
                x.ShiftPixX,x.ShiftPixY,x.Angle);
            
            if contains(struct2array(ver), 'Parallel Computing Toolbox') && UseParallel
                parfor i=1:NumPreSearches
                    TempResults{i} = bayesopt(fun,[ShiftPixX,ShiftPixY,Angle],...
                        'MaxObjectiveEvaluations',PreMaxFunEval,...
                        'AcquisitionFunctionName','expected-improvement-plus',...
                        'UseParallel',false);
                    MinObj(i) = TempResults{i}.MinObjective;
                    MinVal(i,:) = TempResults{i}.XAtMinObjective.Variables;
                end
                Ranges = [range(ShiftPixX.Range),range(ShiftPixY.Range),range(Angle.Range)];
                [BestParam,HalfRange] = AFMImage.choose_best_from_cluster(MinVal,MinObj,NClusters,Ranges);
                parfor i=1:NClusters
                    ShiftPixX = optimizableVariable(...
                        'ShiftPixX',[BestParam(i,1)-HalfRange(i,1),...
                        BestParam(i,1)+HalfRange(i,1)],...
                        'Type','integer');
                    ShiftPixY = optimizableVariable(...
                        'ShiftPixY',[BestParam(i,2)-HalfRange(i,2),...
                        BestParam(i,2)+HalfRange(i,2)],...
                        'Type','integer');
                    Angle = optimizableVariable(...
                        'Angle',[BestParam(i,3)-HalfRange(i,3),...
                        BestParam(i,3)+HalfRange(i,3)],...
                        'Type','real');
                    ClusterResults{i} = bayesopt(fun,[ShiftPixX,ShiftPixY,Angle],...
                        'MaxObjectiveEvaluations',MaxFunEval,...
                        'AcquisitionFunctionName','lower-confidence-bound',...
                        'UseParallel',false);
                    ClusterMin(i) = ClusterResults{i}.MinObjective;
                end
                
                [~,Best] = min(ClusterMin);
                Results = ClusterResults{Best};
            else
                Results = bayesopt(fun,[ShiftPixX,ShiftPixY,Angle],...
                        'MaxObjectiveEvaluations',MaxFunEval,...
                        'AcquisitionFunctionName','lower-confidence-bound',...
                        'UseParallel',UseParallel);
            end
            
            Fig = figure;
            
            MinObj = Results.XAtMinObjective.Variables;
            MinShiftPixX = MinObj(1);
            MinShiftPixY = MinObj(2);
            MinAngle = MinObj(3);
            [Overlay1,Overlay2,OverlayMask1,OverlayMask2] = AFMImage.overlay_two_images(Channel1,Channel2,MinShiftPixX,MinShiftPixY,MinAngle,Mask1,Mask2);
            imshowpair(OverlayMask1,OverlayMask2)
            title(sprintf('%i,%i,%.3f',MinShiftPixX,MinShiftPixY,MinAngle));
            
            Fig = figure;
            
            subplot(2,2,1)
            MinObj = ClusterResults{1}.XAtMinObjective.Variables;
            ShiftPixX = MinObj(1);
            ShiftPixY = MinObj(2);
            Angle = MinObj(3);
            [Overlay1,Overlay2,OverlayMask1,OverlayMask2] = AFMImage.overlay_two_images(Channel1,Channel2,ShiftPixX,ShiftPixY,Angle,Mask1,Mask2);
            imshowpair(OverlayMask1,OverlayMask2)
            title(sprintf('%i,%i,%.3f',ShiftPixX,ShiftPixY,Angle));
            subplot(2,2,2)
            MinObj = ClusterResults{2}.XAtMinObjective.Variables;
            ShiftPixX = MinObj(1);
            ShiftPixY = MinObj(2);
            Angle = MinObj(3);
            [Overlay1,Overlay2,OverlayMask1,OverlayMask2] = AFMImage.overlay_two_images(Channel1,Channel2,ShiftPixX,ShiftPixY,Angle,Mask1,Mask2);
            imshowpair(OverlayMask1,OverlayMask2)
            title(sprintf('%i,%i,%.3f',ShiftPixX,ShiftPixY,Angle));
            subplot(2,2,3)
            MinObj = ClusterResults{3}.XAtMinObjective.Variables;
            ShiftPixX = MinObj(1);
            ShiftPixY = MinObj(2);
            Angle = MinObj(3);
            [Overlay1,Overlay2,OverlayMask1,OverlayMask2] = AFMImage.overlay_two_images(Channel1,Channel2,ShiftPixX,ShiftPixY,Angle,Mask1,Mask2);
            imshowpair(OverlayMask1,OverlayMask2)
            title(sprintf('%i,%i,%.3f',ShiftPixX,ShiftPixY,Angle));
            subplot(2,2,4)
            MinObj = ClusterResults{4}.XAtMinObjective.Variables;
            ShiftPixX = MinObj(1);
            ShiftPixY = MinObj(2);
            Angle = MinObj(3);
            [Overlay1,Overlay2,OverlayMask1,OverlayMask2] = AFMImage.overlay_two_images(Channel1,Channel2,ShiftPixX,ShiftPixY,Angle,Mask1,Mask2);
            imshowpair(OverlayMask1,OverlayMask2)
            title(sprintf('%i,%i,%.3f',ShiftPixX,ShiftPixY,Angle));
            
            drawnow
            
            % Create standard AFMImage-ChannelStruct with new Origin and
            % ScanAngle. Compute real-world-shift from ShiftPixX/Y
            ShiftX = MinShiftPixX*Channel2.ScanSizeX/Channel2.NumPixelsX;
            ShiftY = MinShiftPixY*Channel2.ScanSizeY/Channel2.NumPixelsY;
            ChannelOut = Channel2;
            ChannelOut.OriginX = Channel2.OriginX + ShiftX;
            ChannelOut.OriginY = Channel2.OriginY + ShiftY;
            ChannelOut.ScanAngle = mod(Channel2.ScanAngle + MinAngle,360);
            ChannelOut.Name = 'Overlay';
        end
        
        function Loss = overlay_loss(Channel1,Channel2,Mask1,Mask2,ShiftPixX,ShiftPixY,Angle)
            
            [Overlay1,Overlay2,OverlayMask1,OverlayMask2] = AFMImage.overlay_two_images(Channel1,Channel2,ShiftPixX,ShiftPixY,Angle,Mask1,Mask2);
            NullMask1 = zeros(size(OverlayMask1));
            NullMask2 = zeros(size(OverlayMask1));
            NullMask1(Overlay1 == 0) = 1;
            NullMask2(Overlay2 == 0) = 1;
            NullCombo = ~NullMask1 & ~NullMask2;
            
            XorCombo = ~xor(OverlayMask1,OverlayMask2) & NullCombo;
            OverlapIndizes = find(XorCombo);
            
%             AndCombo = and(OverlayMask1,OverlayMask2) & NullCombo;
%             Loss = -sum(AndCombo,'all');
            
            WeightedArea = sum(abs(Overlay1(OverlapIndizes))) + sum(abs(Overlay2(OverlapIndizes)));
            WeightedDifference = sum(abs(Overlay1(OverlapIndizes)-Overlay2(OverlapIndizes)));
            
            NormFactor = sum(abs(Overlay1),'all');
            
            Loss = (WeightedDifference*.5 - WeightedArea)/NormFactor;
        end
        
        function [Overlay1,Overlay2,OverlayMask1,OverlayMask2] = overlay_two_images(Channel1,Channel2,ShiftPixX,ShiftPixY,Angle,Mask1,Mask2)
            
            % Resize the bigger Channel (in ScanSize-per-Pixel) so
            % imagesizes correspond to ScanSizes. Do nothing, if they are
            % exactly the same
            SizePerPixel1 = Channel1.ScanSizeX/Channel1.NumPixelsX;
            SizePerPixel2 = Channel2.ScanSizeX/Channel2.NumPixelsX;
            if SizePerPixel1 > SizePerPixel2
                ScaleMultiplier = SizePerPixel1/SizePerPixel2;
                Channel1.Image = imresize(Channel1.Image,ScaleMultiplier);
                Channel1.NumPixelsX = size(Channel1.Image,1);
                Channel1.NumPixelsY = size(Channel1.Image,2);
            elseif SizePerPixel1 > SizePerPixel2
                ScaleMultiplier = SizePerPixel2/SizePerPixel1;
                Channel2.Image = imresize(Channel2.Image,ScaleMultiplier);
                Channel2.NumPixelsX = size(Channel2.Image,1);
                Channel2.NumPixelsY = size(Channel2.Image,2);
            end
            
            % If no shifts and angles are provided, they will be calculated
            % from the positional data of the Channel-struct (Operation mode
            % in standard overlay)
            if nargin == 2
                ShiftPixX = (Channel2.OriginX - Channel1.OriginX)*Channel2.NumPixelsX/Channel2.ScanSizeX;
                ShiftPixY = (Channel2.OriginY - Channel1.OriginY)*Channel2.NumPixelsY/Channel2.ScanSizeY;
                Angle = Channel2.ScanAngle - Channel1.ScanAngle;
            end
            
            C2NumX = Channel2.NumPixelsX;
            C2NumY = Channel2.NumPixelsY;
            Channel2.Image = imrotate(Channel2.Image,Angle,'bicubic','loose');
            Channel2.NumPixelsX = size(Channel2.Image,1);
            Channel2.NumPixelsY = size(Channel2.Image,2);
            DiffX = Channel2.NumPixelsX - C2NumX;
            DiffY = Channel2.NumPixelsY - C2NumY;
            ShiftPixX = (ShiftPixX+1) - (floor(DiffX/2)+1);
            ShiftPixY = (ShiftPixY+1) + (ceil(DiffY/2)-1);
            NumPixOverlayX = 2*Channel1.NumPixelsX + 2*Channel2.NumPixelsX;
            NumPixOverlayY = 2*Channel1.NumPixelsY + 2*Channel2.NumPixelsY;
            TempOverlay1 = zeros(NumPixOverlayX,NumPixOverlayY);
            TempOverlay2 = TempOverlay1;
            TempOverlayMask1 = TempOverlay1;
            TempOverlayMask2 = TempOverlay1;
            X1 = (floor(NumPixOverlayX/2)+1):(floor(NumPixOverlayX/2)+Channel1.NumPixelsX);
            Y1 = (floor(NumPixOverlayY/2)+1):(floor(NumPixOverlayY/2)+Channel1.NumPixelsY);
            X2 = (floor(NumPixOverlayX/2)+1-ShiftPixY):(floor(NumPixOverlayX/2)+Channel2.NumPixelsX-ShiftPixY);
            Y2 =  (floor(NumPixOverlayY/2+ShiftPixX)+1):(floor(NumPixOverlayY/2)+Channel2.NumPixelsY+ShiftPixX);
            TempOverlay1(X1,Y1) = Channel1.Image;
            TempOverlay2(X2,Y2) = Channel2.Image;
            if nargin == 7
                Mask2 = imrotate(Mask2,Angle,'nearest','loose');
                TempOverlayMask1(X1,Y1) = Mask1;
                TempOverlayMask2(X2,Y2) = Mask2;
            end
            
            MinX = min([X1(1) X2(1)]);
            MaxX = max([X1(end) X2(end)]);
            MinY = min([Y1(1) Y2(1)]);
            MaxY = max([Y1(end) Y2(end)]);
            Overlay1 = TempOverlay1(MinX:MaxX,MinY:MaxY);
            Overlay2 = TempOverlay2(MinX:MaxX,MinY:MaxY);
            if nargin == 7
                OverlayMask1 = TempOverlayMask1(MinX:MaxX,MinY:MaxY);
                OverlayMask2 = TempOverlayMask2(MinX:MaxX,MinY:MaxY);
            end
        end
        
        function overlay_wrapper_old(FirstClassInstance,OverlayedClassInstance,...
                BackgroundPercent,MinOverlap,AngleRange,UseParallel,...
                MaxFunEval,PreMaxFunEval,NumPreSearches,NClusters)
            % overlay_wrapper_old(obj,OverlayedImageClassInstance,BackgroundPercent,MinOverlap,AngleRange,UseParallel,MaxFunEval,PreMaxFunEval,NumPreSearches,NClusters)
            
            if nargin < 3
                BackgroundPercent = 0;
                MinOverlap = .6;
                AngleRange = 30;
                UseParallel = true;
                MaxFunEval = 200;
                PreMaxFunEval = 30;
                NumPreSearches = 84;
                NClusters = 4;
            end
            
            Channel1 = FirstClassInstance.get_channel('Processed');
            Channel2 = OverlayedClassInstance.get_channel('Processed');
            
            [OutChannel,ScaleMultiplier,WhoScaled] = AFMImage.overlay_parameters_by_bayesopt(Channel1,Channel2,...
                BackgroundPercent,MinOverlap,AngleRange,UseParallel,MaxFunEval,PreMaxFunEval,NumPreSearches,NClusters);
            
            [Overlay1,Overlay2] = AFMImage.overlay_two_images(Channel1,OutChannel);
            
            % Create Overlay Background masks
            OverlayMask1 = AFMImage.mask_background_by_threshold(Overlay1,BackgroundPercent,'on');
            OverlayMask2 = AFMImage.mask_background_by_threshold(Overlay2,BackgroundPercent,'on');
            
            if isequal(WhoScaled,'Channel1')
                Scale1 = ScaleMultiplier;
                Scale2 = 1;
            elseif isequal(WhoScaled,'Channel2')
                Scale1 = 1;
                Scale2 = ScaleMultiplier;
            elseif isequal(WhoScaled,'No one')
                Scale1 = 1;
                Scale2 = 1;
            end
            
            % write the overlays into the corresponding Class instances
            SecondOut = OutChannel;
            SecondMaskOut = OutChannel;
            SecondOut.Image = Overlay2;
            SecondMaskOut.Image = OverlayMask2;
            SecondMaskOut.Unit = 'Logical';
            SecondMaskOut.Name = 'Overlay Mask';
            SecondOut.NumPixelsX = size(SecondOut.Image,1);
            SecondOut.NumPixelsY = size(SecondOut.Image,2);
            SecondOut.ScanSizeX = Channel2.ScanSizeX*(SecondOut.NumPixelsX/(Channel2.NumPixelsX*Scale2));
            SecondOut.ScanSizeY = Channel2.ScanSizeY*(SecondOut.NumPixelsY/(Channel2.NumPixelsY*Scale2));
            SecondMaskOut.NumPixelsX = SecondOut.NumPixelsX;
            SecondMaskOut.NumPixelsY = SecondOut.NumPixelsY;
            SecondMaskOut.ScanSizeX = SecondOut.ScanSizeX;
            SecondMaskOut.ScanSizeY = SecondOut.ScanSizeY;
            
            FirstOut = Channel1;
            FirstMaskOut = Channel1;
            FirstOut.Image = Overlay1;
            FirstOut.Name = 'Overlay';
            FirstMaskOut.Image = OverlayMask1;
            FirstMaskOut.Unit = 'Logical';
            FirstMaskOut.Name = 'Overlay Mask';
            FirstOut.NumPixelsX = size(FirstOut.Image,1);
            FirstOut.NumPixelsY = size(FirstOut.Image,2);
            FirstOut.ScanSizeX = Channel2.ScanSizeX*(FirstOut.NumPixelsX/(Channel1.NumPixelsX*Scale1));
            FirstOut.ScanSizeY = Channel2.ScanSizeY*(FirstOut.NumPixelsY/(Channel1.NumPixelsY*Scale1));
            FirstMaskOut.NumPixelsX = FirstOut.NumPixelsX;
            FirstMaskOut.NumPixelsY = FirstOut.NumPixelsY;
            FirstMaskOut.ScanSizeX = FirstOut.ScanSizeX;
            FirstMaskOut.ScanSizeY = FirstOut.ScanSizeY;
            
            FirstClassInstance.Channel(end+1) = FirstOut;
            FirstClassInstance.Channel(end+1) = FirstMaskOut;
            
            OverlayedClassInstance.Channel(end+1) = SecondOut;
            OverlayedClassInstance.Channel(end+1) = SecondMaskOut;
            
            FirstClassInstance.hasOverlay = true;
            OverlayedClassInstance.hasOverlay = true;
        end
        
        function overlay_wrapper(FirstClassInstance,OverlayedClassInstance,...
                BackgroundPercent,MinOverlap,AngleRange,UseParallel,...
                MaxFunEval,PreMaxFunEval,NumPreSearches,NClusters)
            % overlay_wrapper(FirstClassInstance,OverlayedClassInstance,...
            %    BackgroundPercent,MinOverlap,AngleRange,UseParallel,...
            %    MaxFunEval,PreMaxFunEval,NumPreSearches,NClusters)
            
            if nargin < 3
                BackgroundPercent = 0;
                MinOverlap = .6;
                AngleRange = 30;
                UseParallel = true;
                MaxFunEval = 200;
                PreMaxFunEval = 30;
                NumPreSearches = 84;
                NClusters = 4;
            end
            
            Channel1 = FirstClassInstance.get_channel('Processed');
            Channel2 = OverlayedClassInstance.get_channel('Processed');
            
            OutChannel = AFMImage.overlay_parameters_by_bayesopt(Channel1,Channel2,...
                BackgroundPercent,MinOverlap,AngleRange,UseParallel,MaxFunEval,PreMaxFunEval,NumPreSearches,NClusters);
            
            OverlayedClassInstance.set_channel_positions(OutChannel.OriginX,OutChannel.OriginY,OutChannel.ScanAngle);
            
        end
        
        function create_overlay_group(AllToOneBool,varargin)
            % create_overlay_group(AllToOneBool,varargin)
            %
            % AllToOneBool... if true, all overlay-relations are tied to
            % the first Class-instance that is input to the function
            %
            % varargin... Input all the class-instances that should be
            % overlayed (Works for all classes inheriting from
            % AFMBaseClass).
            %
            % Example: E.create_overlay_group(false,E.FM{1},E.FM{2},E.FM{3},...
            %                      E.I{1},E.I{2})
            % This will overlay FM{1}->FM{2}->FM{3}->I{1}->I{2} in this
            % exact order
            
            for i=1:length(varargin)
                Names{i} = varargin{i}.Name;
            end
            
            for i=1:(length(varargin)-1)
                if AllToOneBool
                    AFMImage.overlay_wrapper(varargin{1},varargin{i+1})
                else
                    AFMImage.overlay_wrapper(varargin{i},varargin{i+1})
                end
            end
            
            for i=1:length(varargin)
                varargin{i}.OverlayGroup.hasOverlayGroup = true;
                varargin{i}.OverlayGroup.Size = length(varargin);
                varargin{i}.OverlayGroup.Names = Names;
            end
            
        end
        
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
        
        function OutImage = subtract_line_fit_diffhist_method(InImage)
            
            NumProfiles = size(InImage,1);
            NumPoints = size(InImage,2);
            
            for i=1:NumProfiles
                Line = InImage(i,:)';
                CutOff = AFMImage.automatic_cutoff(Line);
                Indizes = find(Line<=CutOff);
                LineFit = polyfit(Indizes,Line(Indizes),1);
                LineEval = [1:NumPoints]'*LineFit(1) + LineFit(2);
                Line = Line - LineEval;
                InImage(i,:) = Line;
            end
            OutImage = InImage;
        end
        
        function OutImage = subtract_line_fit_std_method(InImage)
            NumProfiles = size(InImage,1);
            NumPoints = size(InImage,2);
            for i=1:NumProfiles
                Line = InImage(i,:);
                Line = reshape(Line,[],1);
                [Sorted,Indizes] = sort(Line,'descend');
                InvSampleRate = ceil(NumPoints/32^2);
                k = 1;
                for i=1:InvSampleRate:NumPoints
                    STDLine(k) =  std(Sorted(1:i));
                    k = k + 1;
                end
                [~,PeakIdx] = findpeaks(STDLine);
                if isempty(PeakIdx)
                    ThreshIndex = 1;
                else
                    ThreshIndex = PeakIdx(end);
                end
                LineFit = polyfit(Indizes(ThreshIndex:end),Line(Indizes(ThreshIndex:end)),1);
                LineEval = [1:NumPoints]'*LineFit(1) + LineFit(2);
                Line = Line - LineEval;
                InImage(i,:) = Line;
                Line = [];
            end
            OutImage = InImage;
        end
        
        function OutImage = subtract_line_fit_rov(InImage,WindowSize)
            
            NumProfiles = size(InImage,1);
            NumPoints = size(InImage,2);
            ForwardRoV = zeros(NumProfiles,NumPoints);
            BackwardRoV = zeros(NumProfiles,NumPoints);
            for i=1:10
                InImage = smoothdata(InImage,2,'movmedian',5);
                imshow(InImage,[])
                drawnow
            end
            for i=1:NumProfiles
                for j=(1+WindowSize):(NumPoints-WindowSize)
                    ForwardRoV(i,j) = std(InImage(i,(j+1):(j+WindowSize)))/std(InImage(i,(j-WindowSize):(j-1)));
                    BackwardRoV(i,j) = 1/ForwardRoV(i,j);
                end
                ForwardRoV(i,:) = normalize(smoothdata(ForwardRoV(i,:),'movmean',5),'range');
                BackwardRoV(i,:) = normalize(smoothdata(BackwardRoV(i,:),'movmean',5),'range');
                yyaxis left
                plot(InImage(i,:))
                yyaxis right
                plot(1:NumPoints,ForwardRoV(i,:),1:NumPoints,BackwardRoV(i,:))
                drawnow
                
            end
            % Dont really know where to go with this atm
            OutImage = 1;
        end
        
        function OutImage = subtract_line_fit_vertical_rov(InImage,WindowSize,DebugBool)
            
            NumProfiles = size(InImage,1);
            NumPoints = size(InImage,2);
            if nargin < 2
                DebugBool = false;
                WindowSize = .2;
            elseif nargin < 3
                DebugBool = false;
            end
            
            WindowSize = round(WindowSize*NumPoints);
            
            RoV = zeros(NumProfiles,NumPoints);
            MaxIdx = zeros(NumProfiles,1);
            Smoothed = smoothdata(InImage,2,'gaussian',5);
            Smoothed = smoothdata(Smoothed,2,'gaussian',5);
            [Sorted, SortedIndex] = sort(Smoothed,2,'descend');
            
            % Pre debug block    
            if DebugBool
                f = figure('Units','normalized','Position',[0 0 .9 .9]);
            end
            
            for i=1:NumProfiles
                for j=(1+WindowSize):(NumPoints-WindowSize)
                    RoV(i,j) = std(Sorted(i,(j-WindowSize):(j-1)))/std(Sorted(i,(j+1):(j+WindowSize)));
                end
                [~,MaxIdx(i)] = max(RoV(i,:));
                
                % Debug Block
                if DebugBool
                    yyaxis left
                    plot(SortedIndex(i,1:MaxIdx(i)), Sorted(i,1:MaxIdx(i)),'ro',...
                        SortedIndex(i,MaxIdx(i)+1:end),Sorted(i,MaxIdx(i)+1:end),'go');
                    yyaxis right
                    plot(1:NumPoints,RoV(i,:))
                    drawnow
                end
            end
            % Smooth out the hist cutoff by fitting a gp and using its mean
            KernelLambda = 5000;
            KernelSigma = 1;
            GPNoise = 0;
            RectMaxIdx = round(predictGP_mean(1:NumProfiles,1:NumProfiles,KernelSigma,KernelLambda,MaxIdx,GPNoise));
            for i=1:NumProfiles
                LineFit = polyfit(SortedIndex(i,RectMaxIdx(i)+round(.25*(NumPoints-RectMaxIdx(i))):end),...
                    Sorted(i,RectMaxIdx(i)+round(.25*(NumPoints-RectMaxIdx(i))):end),1);
                LineEval = [1:NumPoints]'*LineFit(1) + LineFit(2);
                Line = InImage(i,:)';
                Line = Line - LineEval;
                InImage(i,:) = Line;
            end
            
            OutImage = InImage;
            
            % Second debug block
            if DebugBool
                close(f)
                f = figure('Units','normalized','Position',[0 0 .9 .9]);
                subplot(2,1,1)
                plot(1:NumProfiles,MaxIdx,'rx',1:NumProfiles,RectMaxIdx,'g-')
                subplot(2,1,2)
                imshowpair(Smoothed,OutImage,'montage')
                drawnow
                pause(5)
                close(f)
            end
        end
        
        function OutMask = mask_background_by_threshold(Image,PercentOfRange,AutoMode)
            
            if nargin < 2
                PercentOfRange = 5;
                AutoMode = 'off';
            end
            
            if nargin < 3
                AutoMode = 'off';
            end
            
            if isequal(lower(AutoMode),'on')
                try
                    VecImage = reshape(Image,1,[]);
                    Sorted = sort(VecImage,'descend');
                    InvSampleRate = ceil(length(VecImage)/32^2);
                    k = 1;
                    for i=1:InvSampleRate:length(VecImage)
                        STDLine(k) =  std(Sorted(1:i));
                        k = k + 1;
                    end
                    Peaks = findpeaks(STDLine);
                    Thresh = Peaks(end);
                catch
                    PercentOfRange = 5;
                    Thresh = range(Image,'all')*PercentOfRange/100;
                end
            else
                Thresh = range(Image,'all')*PercentOfRange/100;
            end
            
            [Row,Col] = find((abs(Image)<=Thresh));
            OutMask = zeros(size(Image));
            for i=1:length(Row)
                OutMask(Row(i),Col(i)) = 1;
            end
            
        end
        
        function OutMask = mask_apex_points_by_rotation(InImage,NumRotations,BackgroundMask)
            
            RotAngles = 0:90/NumRotations:90;
            RotAngles(end) = [];
            OutMask = zeros(size(InImage));
            if (nargin > 2) && ~isequal(size(InImage),size(BackgroundMask))
                BackgroundMask = ones(size(InImage));
            elseif nargin < 2
                BackgroundMask = ones(size(InImage));
            end
            h = waitbar(0,'Lets start rotating');
            for i=1:NumRotations
                waitbar(i/NumRotations,h,sprintf('%i/%i Rotation-scans done',i,NumRotations))
                TempImage = imrotate(InImage,RotAngles(i));
                TempMask = zeros(size(TempImage));
                NumProfiles = size(TempImage,1);
                for j=1:NumProfiles
                    [Peaks,PeakIndices] = findpeaks(TempImage(j,:));
                    TempMask(j,PeakIndices) = 1;
                end
                TempMask = imrotate(TempMask,(-RotAngles(i)),'crop');
                TempImage = imrotate(TempImage,(-RotAngles(i)),'crop');
                while ~sum(find(TempImage(1,:)))
                    TempImage(1,:) = [];
                    TempMask(1,:) = [];
                end
                while ~sum(find(TempImage(end,:)))
                    TempImage(end,:) = [];
                    TempMask(end,:) = [];
                end
                while ~sum(find(TempImage(:,1)))
                    TempImage(:,1) = [];
                    TempMask(:,1) = [];
                end
                while ~sum(find(TempImage(:,end)))
                    TempImage(:,end) = [];
                    TempMask(:,end) = [];
                end
                TempMask = imresize(TempMask,size(InImage));
                OutMask = OutMask + TempMask;
            end
            close(h)
            
            OutMask(BackgroundMask == 1) = 0;
            
            OutMask = OutMask./NumRotations;
            
        end
        
        function [OutMask,AngleMask] = mask_apex_points_by_convolution(InChannel,RadiusMeters,NumRotations,MinKernelSize,BackgroundMask)
            
            OutMask = zeros(size(InChannel.Image));
            RotAngles = 0:180/NumRotations:180;
            RotAngles(end) = [];
            RotAngles = deg2rad(RotAngles);
            MinKernelSize = MinKernelSize + (1-mod(MinKernelSize,2));
            
            MaxSizePerPixel = RadiusMeters/MinKernelSize;
            SizePerPixel = InChannel.ScanSizeX/InChannel.NumPixelsX;
            KernelSize = round(RadiusMeters/SizePerPixel);
            KernelSize = KernelSize + (1-mod(KernelSize,2));
            if SizePerPixel > MaxSizePerPixel
                Image = imresize(InChannel.Image,SizePerPixel/MaxSizePerPixel,'bicubic');
                KernelSize = MinKernelSize;
                SizePerPixel = MaxSizePerPixel;
            else
                Image = InChannel.Image;
            end
            
            Kernel = zeros(KernelSize,KernelSize,NumRotations);
            FlattenKernel = zeros(KernelSize,KernelSize);
            KernelCenter = (KernelSize+1)/2;
            Radius = KernelCenter;
            % Create Kernel
            for k=1:NumRotations
                for i=1:KernelSize
                    for j=1:KernelSize
                        if (i-KernelCenter)^2 + (j-KernelCenter)^2 <= Radius^2
                        a = [(i-KernelCenter) (j-KernelCenter)];
                        b = [cos(RotAngles(k)) sin(RotAngles(k))];
                        Dist = dot(a,b)/norm(b);
                        Kernel(i,j,k) = Radius^3/8 - abs(Dist^3);
                        if k == 1
                            FlattenKernel(i,j) = 1;
                        end
                        end
                    end
                end
            end
%             FlattenKernel = imresize(FlattenKernel,[10*KernelSize+1 10*KernelSize+1],'nearest');
%             FlattenKernel = FlattenKernel/(sum(FlattenKernel,'all'));
            
            FlattenedImage = Image - conv2(Image,FlattenKernel,'same');
            
%             % Debug section
%             tiledlayout(5,5)
%             for i = 1:NumRotations
%                 nexttile
%                 imshow(Kernel(:,:,i),[])
%             end
            
            TempMask = zeros([size(Image) NumRotations]);
            TempMask = gpuArray(TempMask);
            FlattenedImage = gpuArray(FlattenedImage);
            Kernel = gpuArray(Kernel);
            TempMask = convn(FlattenedImage,Kernel);
            AngleMask = gpuArray(zeros(size(TempMask,1),size(TempMask,2)));
            
            [OutMask,AngleIndex] = max(TempMask,[],3);
            
            for i=1:NumRotations
                AngleMask(AngleIndex == i) = RotAngles(i);
            end
            
            if nargin == 5
                OutMask(imresize(BackgroundMask,size(TempMask(:,:,1)),'nearest')==1) = 0;
            end
            
            OutMask(:,end-(KernelCenter-2):end) = [];
            OutMask(end-(KernelCenter-2):end,:) = [];
            OutMask(1:(KernelCenter-1),:) = [];
            OutMask(:,1:(KernelCenter-1)) = [];
            
            AngleMask(:,end-(KernelCenter-2):end) = [];
            AngleMask(end-(KernelCenter-2):end,:) = [];
            AngleMask(1:(KernelCenter-1),:) = [];
            AngleMask(:,1:(KernelCenter-1)) = [];
        end
        
        function OutImage = masked_plane_fit(Channel,Mask)
            
            Image = Channel.Image;
            
            % Convert Image to Point Cloud for plane fit
            [X,Y,Z] = AFMImage.convert_masked_to_point_cloud(Channel,Mask);
            
            [Norm,~,Point] = AFMImage.affine_fit([X Y Z]);
            Plane = zeros(size(Image));
            % Create the plane that can then be subtracted from the
            % complete height data to generate the leveled height data.
            for i=1:size(Image,1)
                for j=1:size(Image,2)
                    Plane(i,j) = (Point(3)-Norm(1)/Norm(3)*i*Channel.ScanSizeX/Channel.NumPixelsX-Norm(2)/Norm(3)*j*Channel.ScanSizeY/Channel.NumPixelsY);
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
        
        function DepthDependendTipRadius = calculate_depth_dependend_tip_data(ProjectedTipArea,RangePercent)
            
            if nargin < 2
                RangePercent = 100;
            end
            k = 2;
            while ProjectedTipArea(1) == ProjectedTipArea(k)
                k=k+1;
            end
            MinIdx = k;
            MaxIdx = floor(RangePercent/100*length(ProjectedTipArea));
            DepthDependendTipRadius = zeros(MaxIdx,1);
            
            for i=MinIdx:MaxIdx
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
                ProjAParabola = fittype('pi*x/a',...
                    'dependent',{'y'},'independent',{'x'},...
                    'coefficients',{'a'},...
                    'options',SphOpt);
                Depth = [1:i]'.*1e-9;
                X = Depth/range(Depth);
                Y = ProjectedTipArea(1:i)/range(ProjectedTipArea(1:i));
                ParabolaFit = fit(X,...
                    Y,...
                    ProjAParabola,...
                    'Weights',[1:i]'.^2);
                warning('off')
                ParabolaFit.a = ParabolaFit.a*range(Depth)/range(ProjectedTipArea(1:i));
                warning('on')
                RParabola = 1/(2*ParabolaFit.a);
                DepthDependendTipRadius(i) = RParabola;
%                 plot(Depth,ProjectedTipArea(1:i),'rO',Depth,feval(ParabolaFit,Depth),'g')
            end
            % Fill the first 4 nm with the data from the 5th nm
            for i=1:MinIdx
                DepthDependendTipRadius(i) = DepthDependendTipRadius(MinIdx);
            end
        end
        
        function [ComboSum,AspectRatio,MeanX,MeanY,Angle] = fibril_feature_kernel(Height, ApexMap, PosX, PosY,...
                AdaptiveWindowX,AdaptiveWindowY)
            
            % In order to avoid directional biases, points should be taken
            % from a circle, not a rectangle
            Radius = (min([AdaptiveWindowX AdaptiveWindowY])-1)/2;
            InIndices(1,:) = [PosX PosY];
            i = 0;
            j = 0;
            k = 1;
            while i^2+j^2 <= Radius^2
                while i^2+j^2 <= Radius^2
                    InIndices(k,:) = [PosX+i PosY+j];
                    k = k + 1;
                    if j
                    InIndices(k,:) = [PosX+i PosY-j];
                    k = k + 1;
                    end
                    if i
                    InIndices(k,:) = [PosX-i PosY+j];
                    k = k + 1;
                    end
                    if i && j
                    InIndices(k,:) = [PosX-i PosY-j];
                    k = k + 1;
                    end
                    j = j + 1;
                end
                i = i + 1;
                j = 0;
            end
            Weights = zeros(length(InIndices),1);
            for i=1:length(InIndices)
                Weights(i) = abs(Height.Image(InIndices(i,1),InIndices(i,2))...
                    *ApexMap.Image(InIndices(i,1),InIndices(i,2)));
            end
            [EigenVectors,~,EigenValues,~,~,Means] = pca(InIndices,'Weights',Weights);
            MeanX = Means(1);
            MeanY = Means(2);
            Angle = rad2deg(atan(EigenVectors(1,1)/EigenVectors(2,1)));
            ComboSum = sum(Weights,'all');
            AspectRatio = EigenValues(1)/EigenValues(2);
        end
        
        function OutImage = find_and_replace_outlier_lines(InImage,IQRMult)
            % creates line to line cahnge statistics and then detects and
            % replaces outlies
            
            PixelsX = size(InImage,2);
            PixelsY = size(InImage,1);
            
            DiffPreviousLine = zeros(PixelsY,1);
            DiffNextLine = zeros(PixelsY,1);
            
            for i=2:(PixelsY-1)
                DiffPreviousLine(i) = sum(abs(InImage(i,:) - InImage(i-1,:)));
                DiffNextLine(i) = sum(abs(InImage(i,:) - InImage(i+1,:)));
            end
            
            DiffSum = DiffPreviousLine + DiffNextLine;
            
            Median = nanmedian(DiffSum);
            IQR = iqr(DiffSum);
            
            OutlierLines = find(DiffSum >= Median + IQRMult*IQR);
            
            OutImage = AFMImage.replace_outlier_lines(InImage,OutlierLines);
        end
        
        function OutImage = replace_outlier_lines(InImage,Indizes)
            % replaces lines in an image interpolating between the
            % neighbouring Lines. If Indizes contains neighbouring numbers
            % the function will replace with the non-outliers before and
            % after the outlierblock
            
            OutImage = InImage;
            
            while length(Indizes) > 0
                SpanOfBlock = 1;
                while (SpanOfBlock ~= length(Indizes)) && (Indizes(SpanOfBlock)+1 == Indizes(SpanOfBlock+1))
                    SpanOfBlock = SpanOfBlock + 1;
                end
                BeforeLine = InImage(Indizes(1)-1,:);
                AfterLine = InImage(Indizes(1)+(SpanOfBlock),:);
                NewLine = zeros(SpanOfBlock,size(InImage,2));
                for i=1:SpanOfBlock
                    NewLine(i,:) = (BeforeLine*(SpanOfBlock+1-i)/(SpanOfBlock+1) + AfterLine*i/(SpanOfBlock+1));
                end
                OutImage(Indizes(1):Indizes(1)+SpanOfBlock-1,:) = NewLine;
                Indizes(1:SpanOfBlock) = [];
            end
            
        end
        
        function OutLine = replace_points_of_certain_value_in_line(InLine,Value)
            
            OutLine = InLine;
            k = 1;
            while k < length(InLine)
                SpanOfBlock = 0;
                Indizes = [];
                while (k ~= length(InLine)) && (InLine(k) == Value)
                    SpanOfBlock = SpanOfBlock + 1;
                    Indizes(SpanOfBlock) = k;
                    k = k + 1;
                end
                if isempty(Indizes)
                    k = k + 1;
                end
                if ~isempty(Indizes)
                    if length(Indizes) >= length(InLine) - 2
                        OutLine = InLine;
                        return
                    end
                    if Indizes(1) == 1
                        OutLine(Indizes) = OutLine(Indizes(end)+1);
                        continue
                    else
                        BeforePoint = OutLine(Indizes(1)-1);
                        BeforeIndex = Indizes(1)-1;
                    end
                    if Indizes(end) == length(InLine)
                        OutLine(Indizes) = OutLine(Indizes(1)-1);
                        continue
                    else
                        AfterPoint = OutLine(Indizes(end)+1);
                        AfterIndex = Indizes(end)+1;
                    end
                    OutLine(Indizes) = interp1([BeforeIndex AfterIndex],[BeforePoint AfterPoint],Indizes);
                end
            end
        end
        
        function [OutChannel,GridX,GridY,ProjLength] = project_height_image_to_tilted_surface(InChannel,PolarAngle,AzimuthalAngle,ResMultiplier,UpscaleMult)
            % Creates an alternative projection of a list of points in
            % space X,Y,Z, quantizes them onto a grid, choosing the closest
            % point to the surface, should multiple points fall into a
            % pixel.
            
            % Upscale to UpscaleMult*ResMultiplier in order to get enough points to
            % fill all pixels in the endresult
            
            InChannel.Image = imresize(InChannel.Image,UpscaleMult*ResMultiplier);
            InChannel.NumPixelsX = size(InChannel.Image,1);
            InChannel.NumPixelsY = size(InChannel.Image,2);
            
            [X,Y,Z] = AFMImage.convert_masked_to_point_cloud(InChannel,AFMImage.mask_background_by_threshold(InChannel.Image,50,0));
            
            u1 = [ cos(PolarAngle)*cos(AzimuthalAngle) ; cos(PolarAngle)*sin(AzimuthalAngle) ; -sin(PolarAngle) ];
            u2 = [ -sin(AzimuthalAngle) ; cos(AzimuthalAngle) ; 0 ];
            
            % The pointcloud needs to be shifted perpendicular to the
            % surface to be projected on to avoid numerical instabilities
            % as well as deformations caused by zero-crossings ignored by the
            % norm-operation. Add a vector to X,Y,Z a few ranges in length
            % and in direction u1 x u2 (cross-operation) = u3, with u3
            % being the radial base vector in spherical coordinates;
            
            u3 = [ sin(PolarAngle)*cos(AzimuthalAngle) ; sin(PolarAngle)*sin(AzimuthalAngle) ; cos(PolarAngle) ];
            
            X = X - mean(X);
            Y = Y - mean(Y);
            Z = Z - mean(Z);
            
            Range = max([range(X) range(Y) range(Z)]);
            
            X = X + 10*Range*u3(1);
            Y = Y + 10*Range*u3(2);
            Z = Z + 10*Range*u3(3);
            
            % Define the projection matrix. It is calculated from P=A*A^T,
            % where A=[u1,u2] and u1,u2 are the base vectors of the plane the
            % point cloud needs to be projected to
            P = [cos(PolarAngle)^2*cos(AzimuthalAngle)^2 + sin(AzimuthalAngle)^2 ...
                -cos(AzimuthalAngle)*sin(AzimuthalAngle) + cos(PolarAngle)^2*cos(AzimuthalAngle)*sin(AzimuthalAngle) ...
                -cos(PolarAngle)*cos(AzimuthalAngle)*sin(PolarAngle) ;...
                -cos(AzimuthalAngle)*sin(AzimuthalAngle) + cos(PolarAngle)^2*cos(AzimuthalAngle)*sin(AzimuthalAngle) ...
                cos(AzimuthalAngle)^2 + cos(PolarAngle)^2*sin(AzimuthalAngle)^2 ...
                -cos(PolarAngle)*sin(PolarAngle)*sin(AzimuthalAngle) ;...
                -cos(PolarAngle)*cos(AzimuthalAngle)*sin(PolarAngle) ...
                -cos(PolarAngle)*sin(PolarAngle)*sin(AzimuthalAngle) ...
                sin(PolarAngle)^2];
            
            
            XHat = zeros(size(X));
            YHat = zeros(size(X));
            ZHat = zeros(size(X));
            ProjLength = zeros(size(X));
            GridX = zeros(size(X));
            GridY = zeros(size(X));
            
            L = length(X);
            for i=1:L
                V = [X(i); Y(i); Z(i)];
                VHat = P*V;
                XHat(i) = VHat(1);
                YHat(i) = VHat(2);
                ZHat(i) = VHat(3);
                ProjLength(i) = norm(VHat - V);
                GridX(i) = VHat'*u1;
                GridY(i) = VHat'*u2;
            end
            
            OutChannel = AFMImage.convert_point_cloud_to_image(GridX,GridY,ProjLength,InChannel,1/UpscaleMult);
            
        end
        
        function OutChannel = convert_point_cloud_to_image(X,Y,Z,InChannel,ResMultiplier)
            % Embedds X and Y into an Image grid with Pixel values Z. If
            % multiple pointas fall into one pixel, the one with higher Z
            % value is chosen. If a pixel is empty, it is interpolated
            % from neighboring points.
            
            if nargin < 5
                ResMultiplier = 1;
            end
            
            DummyScaled = imresize(InChannel.Image,ResMultiplier);
            OutChannel = InChannel;
            OutChannel.NumPixelsX = size(DummyScaled,1);
            OutChannel.NumPixelsY = size(DummyScaled,2);
            OutChannel.ScanSizeX = range(X);
            OutChannel.ScanSizeY = range(Y);
            
            % Quantize the X and Y coordinates and sort out multiples with
            % lower z-values.
            XMult = (OutChannel.NumPixelsX-1)/OutChannel.ScanSizeX;
            YMult = (OutChannel.NumPixelsY-1)/OutChannel.ScanSizeY;
            XQ = floor((X-min(X)).*XMult) + 1;
            YQ = floor((Y-min(Y)).*YMult) + 1;
            
            I = zeros(OutChannel.NumPixelsX,OutChannel.NumPixelsY);
            for i=1:length(XQ)
                if I(XQ(i),YQ(i)) < Z(i)
                    I(XQ(i),YQ(i)) = Z(i);
                end
            end
            
            I(I==0) = min(Z);
            
            OutChannel.Image = I;
        end
        
        function OutImage = create_pixel_difference_map(InImage)
            
            OutImage = diff(InImage,1,2);
            OutImage = imresize(OutImage,size(InImage));
            
        end
    end
    
    methods
        % Auxiliary methods
        
        function read_in_header_properties(obj,ImageFullFile)
            % determines imaging type (contact,AC) and reads out general image
            % properties
            
            FileInfo = imfinfo(ImageFullFile);
            
            obj.NumChannels = numel(FileInfo) - 1;
            
            try
                obj.ImagingType = upper(FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==32816)).Value);
            catch
                obj.ImagingType = 'UNKNOWN';
            end
            
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
            elseif isequal(ImageFullFile(end-5:end),'.force') || isequal(ImageFullFile(end-12:end),'.jpk-qi-image')
            elseif isequal(obj.ImagingType,'UNKNOWN')
                error('Unknown imaging mode')
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
                obj.hasSensitivity = true;
            catch
                obj.Sensitivity = nan;
                warning ("Couldn't determine sensitivity. Image isn't calibrated")
                obj.hasSensitivity =false;
            end
            try
                obj.SpringConstant = FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==33076)).Value/...
                    FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==33028)).Value;
                obj.hasSpringConstant = true;
            catch
                obj.SpringConstant = nan;
                warning ("Couldn't determine spring constant. Image isn't calibrated")
                obj.hasSpringConstant = false;
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
                obj.hasSensitivity = true;
            else
                warning('No Sensitivity calibration slot found, image is uncalibarated')
                Uncalibrated=1;
                obj.hasSensitivity = false;
                obj.Sensitivity=1;
            end
            
            if(~isempty(find([FileInfo(1).UnknownTags.ID]==33076)))&&(~isempty(find([FileInfo(1).UnknownTags.ID]==32980)))
                obj.SpringConstant=(FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==33076)).Value)/(FileInfo(1).UnknownTags(find([FileInfo(1).UnknownTags.ID]==33028)).Value);
                obj.hasSpringConstant = true;
            else
                warning('No Spring Constant calibration slot found, image is uncalibarated')
                Uncalibrated=1;
                obj.SpringConstant=1;
                obj.hasSpringConstant = false;
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
                
                trace_type_flag = [];
                for k=1:size(strsp,1)
                    if(strcmp(strsp{k,1},'retrace')==1)
                        if(strcmp(strsp{k+2,1},'true'))
                            trace_type_flag=' (Retrace)';
                        else
                            trace_type_flag=' (Trace)';
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
                
%                 afm_image = afm_image(end:-1:1,:); % mirror Y-pixels to flip image to same orientation as in jpk data processing
                
                obj.Channel(i-1).Image = afm_image;
                if ~isempty(trace_type_flag)
                    obj.Channel(i-1).Name = strcat(Channel_Name,' ',trace_type_flag);
                else
                    obj.Channel(i-1).Name = Channel_Name;
                end
                if isequal(Channel_Name,'Height') ||...
                        isequal(Channel_Name,'Height (measured)') ||...
                        (isequal(Channel_Name,'Lock-In Amplitude') && obj.hasSensitivity) ||...
                        (isequal(Channel_Name,'Vertical Deflection') && obj.hasSensitivity && ~obj.hasSpringConstant) ||...
                        (isequal(Channel_Name,'Lateral Deflection') && obj.hasSensitivity && ~obj.hasSpringConstant)
                    obj.Channel(i-1).Unit = 'm';
                elseif isequal(Channel_Name,'Lock-In Phase')
                    obj.Channel(i-1).Unit = 'deg';
                elseif isequal(Channel_Name,'Vertical Deflection') && obj.hasSensitivity && obj.hasSpringConstant ||...
                        isequal(Channel_Name,'Lateral Deflection') && obj.hasSensitivity && obj.hasSpringConstant
                    obj.Channel(i-1).Unit = 'N';
                elseif isequal(Channel_Name,'Error Signal') ||...
                        (isequal(Channel_Name,'Lock-In Amplitude') && ~obj.hasSensitivity && ~obj.hasSpringConstant) ||...
                        isequal(Channel_Name,'Lateral Deflection')
                    obj.Channel(i-1).Unit = 'V';
                elseif isequal(Channel_Name,'Slope')
                    obj.Channel(i-1).Unit = 'm/m';
                elseif isequal(Channel_Name,'Adhesion')
                    obj.Channel(i-1).Unit = 'N';
                end
                obj.Channel(i-1).NumPixelsX = obj.NumPixelsX;
                obj.Channel(i-1).NumPixelsY = obj.NumPixelsY;
                obj.Channel(i-1).ScanSizeX = obj.ScanSizeX;
                obj.Channel(i-1).ScanSizeY = obj.ScanSizeY;
                obj.Channel(i-1).OriginX = obj.OriginX;
                obj.Channel(i-1).OriginY = obj.OriginY;
                if ~isempty(obj.ScanAngle)
                    obj.Channel(i-1).ScanAngle = obj.ScanAngle;
                else
                    obj.Channel(i-1).ScanAngle = 0;
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
            
            obj.hasHeight = false;
            obj.hasHeightMeasured = false;
            obj.hasErrorSignal = false;
            obj.hasLateralDeflection = false;
            obj.hasLockInAmplitude = false;
            obj.hasLockInPhase = false;
            obj.hasVerticalDeflection = false;
            obj.hasProcessed = false;
            obj.hasSensitivity = false;
            obj.hasSpringConstant = false;
            obj.hasOverlay = false;
            obj.hasBackgroundMask = false;
            obj.hasDeconvolutedCantileverTip = false;
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
        
        function [X,Y,Z] = convert_masked_to_point_cloud(InChannel,Mask)
            
            if nargin < 2
                Mask = ones(size(InChannel.Image));
            end
            
            Image = InChannel.Image;
            
            [X,Y] = find(Mask);
            Z = zeros(length(X),1);
            for i=1:length(X)
                Z(i) = Image(X(i),Y(i));
            end
            
            % the rotation is done around the center of the IMAGE and not
            % the Origin so first shift the origin to mid image, then rotate and finally scale and shift to OriginX/Y;
            X1 = X - (size(Image,1)/2 + 1);
            Y1 = Y - (size(Image,2)/2 + 1);
            X = X1.*cosd(InChannel.ScanAngle) + Y1.*sind(InChannel.ScanAngle);
            Y = -X1.*sind(InChannel.ScanAngle) + Y1.*cosd(InChannel.ScanAngle);
            
            X = X.*InChannel.ScanSizeX/InChannel.NumPixelsX - InChannel.OriginX;
            Y = Y.*InChannel.ScanSizeY/InChannel.NumPixelsY - InChannel.OriginY;
        end
        
        function [XOut,YOut] = rotate_and_shift_point_cloud(X,Y,ShiftX,ShiftY,Angle)
            % Center PC, rotate it and shift it back to original positions
            CenterX = min(X) + range(X)/2;
            CenterY = min(Y) + range(Y)/2;
            X1 = X - CenterX;
            Y1 = Y - CenterY;
            X = X1*cosd(Angle) + Y1*sind(Angle);
            Y = -X1*sind(Angle) + Y1*cosd(Angle);
            X = X + CenterX;
            Y = Y + CenterY;
            
            % Shift it by InputShift and return
            XOut = X - ShiftX;
            YOut = Y - ShiftY;
        end
        
        function R = draw_scalebar_into_current_image(NumPixelsX,NumPixelsY,ScanSizeX,BarToImageRatio,CurAxHeight,CurAxWidth)
            
            if nargin < 4
                BarToImageRatio = 1/5;
                CurAxHeight = NumPixelsY;
                CurAxWidth = NumPixelsX;
            end
            
            ImageRatio = NumPixelsY/NumPixelsX;
            
            if CurAxHeight/CurAxWidth >= ImageRatio
                FontSizeMult = CurAxWidth/1000;
            else
                FontSizeMult = CurAxHeight/1000;
            end
            
            ScalebarThickness = 1/40;
            DistFromBorder = 0.08;
            
            [Multiplier,Unit,SnapTo] = AFMImage.parse_unit_scale(ScanSizeX,'m',BarToImageRatio);
            
            Width = (SnapTo)/(ScanSizeX*Multiplier);
            Height = ScalebarThickness;
            Left = 1-Width*(1+DistFromBorder*1.8);
            Bottom = 1-(DistFromBorder-Height);
            
            BoxDeltaX = Width/4;
            BoxDeltaY = Height*3;
            
            BackR = rectangle('Position',[Left-BoxDeltaX/2 Bottom-5*BoxDeltaY/6 Width+BoxDeltaX Height+BoxDeltaY].*NumPixelsX);
            BackR.FaceColor = 'k';
            BackR.EdgeColor = 'w';
            BackR.LineWidth = 1;
            
            R = rectangle('Position',[Left Bottom Width Height].*NumPixelsX);
            R.FaceColor = 'w';
            R.EdgeColor = 'k';
            R.LineWidth = 2;
            
            FontSize = round(42*FontSizeMult);
            A = text((Left*1.13)*NumPixelsX,(1.005*Bottom-3*BoxDeltaY/6)*NumPixelsY,sprintf('%i %s',SnapTo,Unit),'HorizontalAlignment','center');
            A.Color = 'w';
            A.FontSize = FontSize;
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
        
        function [pBest,IndizesBest,GoFBest,SweepGoF] = line_fit_sweep(Line,StepSizePercent)
            % AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA Don't
            % use this!......
            X = 1:length(Line);
            Indizes = X;
            StepSize = StepSizePercent/100;
            k = 1;
            for Percentage=1:-StepSize:0
            [p,S] = polyfit(X(Indizes),Line(Indizes),1);
            SweepGoF(k) = S.normr;
            Based = Line - polyval(p,X);
%             plot(X,Line,'b.',X,polyval(p,X),'r')
            Indizes = find(Based <= max(Based)*Percentage);
            k = k + 1;
            end
            
        end
        
        function [BestParam,HalfRange] = choose_best_from_cluster(MinVal,MinObj,NClusters,Ranges)
            % [BestParam,HalfRange] = AFMImage.choose_best_from_cluster(MinVal,MinObj,NClusters)
            
            NormedVals(:,1) = MinVal(:,1)/Ranges(1);
            NormedVals(:,2) = MinVal(:,2)/Ranges(2);
            NormedVals(:,3) = MinVal(:,3)/Ranges(3);
            
            % Cluster the Points in parameter-space
            [Indizes,~] = kmeans(NormedVals,NClusters,...
                'Distance','sqeuclidean',...
                'Display','off',...
                'Replicates',5);
            
            for i=1:NClusters
                ClusterVal{i} = MinVal(Indizes==i,:);
                ClusterObj{i} = MinObj(Indizes==i);
                if round(std(ClusterVal{i}(:,1))) > 0
                    HalfRange(i,1) = round(std(ClusterVal{i}(:,1)));
                else
                    HalfRange(i,1) = 50;
                end
                if round(std(ClusterVal{i}(:,2))) > 0
                    HalfRange(i,2) = round(std(ClusterVal{i}(:,2)));
                else
                    HalfRange(i,2) = 50;
                end
                if std(ClusterVal{i}(:,3)) > 0
                    HalfRange(i,3) = std(ClusterVal{i}(:,3));
                else
                    HalfRange(i,3) = 5;
                end
                [~,Best] = min(ClusterObj{i});
                BestParam(i,:) = ClusterVal{i}(Best,:);
            end
            
            
            
        end
        
        function CutOff = automatic_cutoff(Input)
            
            if sum(Input,'all')==0
                CutOff = 0;
                return
            end
            
            NumSamples = ceil(numel(Input)/10);
            Input = reshape(Input,[],1);
            Thresholds = [max(Input):-range(Input)/(NumSamples-1):min(Input)];
            InRangeLine = zeros(NumSamples,1);
            for i=1:NumSamples
                InRangeLine(i) = length(Input(Input >= Thresholds(i)));
            end
            
            DiffIRL = diff(InRangeLine);
            [~,MaxIndex] = max(DiffIRL);
            
            CutOff = Thresholds(MaxIndex);
            
        end
    end
end