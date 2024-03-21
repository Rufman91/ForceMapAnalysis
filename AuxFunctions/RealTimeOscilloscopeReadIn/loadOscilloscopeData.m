function dataTable = loadOscilloscopeData(filePath)
    % Open the file
    fileId = fopen(filePath, 'r');
    
    % Read the header lines and extract relevant information
    sampleFrequencyLine = fgetl(fileId);
    columnsLine = fgetl(fileId);
    fancyNamesLine = fgetl(fileId);
    calibrationSlotsLine = fgetl(fileId);
    unitsLine = fgetl(fileId);
    
    % Initialize a variable to read lines
    tline = fgetl(fileId); % Prime reading with the next line
    
    % Skip any additional header lines
    while ischar(tline)
        if startsWith(tline, '#')
            tline = fgetl(fileId); % Read the next line
        else
            fseek(fileId, -length(tline)-2, 'cof'); % Rewind to start of the current line
            break; % Break the loop once data lines start
        end
    end
    
    % Prepare column headers from the columns line
    columnsHeader = strsplit(columnsLine);
    columnsHeader = strrep(columnsHeader(3:end), ':', ''); % Remove the "# columns:" part and trailing colon
    
    % Read data into a table
    dataTable = readtable(filePath, 'FileType', 'text', 'HeaderLines', 5, ...
                          'Delimiter', ' ', 'ReadVariableNames', false);
                      
    % Assign the column headers
    dataTable.Properties.VariableNames = columnsHeader;
    
    % Close the file
    fclose(fileId);
end
