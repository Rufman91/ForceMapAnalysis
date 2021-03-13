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
            
            obj.define_afm_color_map
            
            obj.initialize_flags
            
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
        
        function define_afm_color_map(obj)
            CMap(:,1) = (0:1/255:1).*2;
            CMap(:,2) = (0:1/255:1).*2 - 0.5;
            CMap(:,3) = (0:1/255:1).*2 - 1;
            CMap(CMap < 0) = 0;
            CMap(CMap > 1) = 1;
            
            obj.CMap = CMap;
            
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
end
