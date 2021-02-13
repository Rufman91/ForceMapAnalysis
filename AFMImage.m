classdef AFMImage < matlab.mixin.Copyable
   
    properties
        % File and Header properties
        Name
        ID
        HostOS
        HostName
        ImagingType
        
    end
    properties
        % Image Data Properties
    
    end
    properties
        % Properties related to Image processing/segmenting/classification
        
    end
    
    methods
        % Main methods of the class
        
        function obj = AFMImage(MapFullFile,DataFolder,TempID)
            % Constructor of the class. Extracts Header properties as well
            % as all available channel-data
            
            obj.ID = TempID;
            
            % get OS and use appropriate fitting system command
            FullOS = computer;
            OS = FullOS(1:3);
            obj.HostOS = OS;
            if isequal(OS,'PCW')
                obj.HostName = getenv('COMPUTERNAME');
            elseif isequal(OS,'GLN')
                obj.HostName = getenv('HOSTNAME');
            elseif isequal(OS,'MAC')
                obj.HostName = getenv('HOSTNAME');
            end
            
            %Open File
            
            
            obj.read_in_header_properties(File);
            
            obj.load_image_channels(File);
            
            
        end
        
    end
    
    methods
        % Auxiliary methods
        
        function read_in_header_properties(File)
            % determines imaging type (contact,AC) and reads out general image
            % properties
            
            
            
        end
        
        function load_image_channels(obj)
            % load in existing image channels and related properties
            
            
            
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
        
    end
end
