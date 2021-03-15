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
    end
    properties
        % Properties related to Image processing/segmenting/classification
        CMap
        MaskBackground
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
            
            obj. CMap = obj.define_afm_color_map(0.35)
            
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
        
    end
    
    methods(Static)
        % Static Main Methods
        
        function [OutImage] = subtract_line_fit_hist(InImage,NumProfiles,CutOff)
            
            NumPoints = length(InImage(1,:));
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
        
        function [OutMask] = mask_background(Image)
            
        end
        
        function CMap = define_afm_color_map(PlusBrightness)
            if nargin == 0
                PlusBrightness = .35;
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
            
            
            h = waitbar(0,'Please wait, Processing...');
            
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
            
            % TODO
            
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
                    if(strcmp(Bline_adjust,'No'))
                        afm_image=((double(imread(ImageFullFile,i))*multiplyer))+offset;
                    else
                        Details_Img.Baseline_N=(Baseline_Raw*multiplyer)+offset;
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
        
    end
end
