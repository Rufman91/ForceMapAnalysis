classdef AFMImage < matlab.mixin.Copyable
   
    properties
        % File and Header properties
        Name
        ID
        HostOS
        HostName
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
            
        end
        
    end
end
