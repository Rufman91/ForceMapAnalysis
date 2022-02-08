function FunctionHeader = set_up_function_header(OutputNames,FunctionName,Required,NameValue,OptionalPositional)
% Sets up a standardized function header to use with matlab inputParser and
% optional conditions, expected values and defaults. Function output is a
% string that can be copied to the editor from the workspace or variables
% tab but it also copies the text to the clipboard for easy insert via
% ctrl+v.
%
% If an input type is not used, replace it with empty brackets.
% Valid example:
%       >> function_header_setup_helper({'Out1','Out2'},'this_is_a_test',{'ReqVar1','ReqVar2','ReqVar3'},'NV1',[]);
%
% Default defaults are [] and default conditions are true (accepts every input).
% See inputParser documentation for details.
%
% OutputNames ... Names of outputs (input as string or cell of strings)
% FunctionName ... Type in function name as a string
% Required ... Names of positional required inputs (input as string or cell of strings)
% NameValue ... Names of optional NameValue inputs (input as string or cell of strings)
% OptionalPositional ... Names of positional optional inputs (input as string or cell of strings)

% Set up FunctionHeader and Code string
FHFirstLine = 'function ';
FHDescription = ['% <FUNCTION DESCRIPTION HERE>' newline '%'];
FHCode = ['p = inputParser;' newline...
    'p.FunctionName = "' FunctionName '";' newline...
    'p.CaseSensitive = false;' newline...
    'p.PartialMatching = true;'];
FHCodeReassign = '% Assign parsing results to named variables';
ReqNameCommaSep = '';

% Outputs
if ~isempty(OutputNames)
    if ~iscell(OutputNames)
        Temp = OutputNames;
        OutputNames = [];
        OutputNames{1} = Temp;
    end
    FHFirstLine = [FHFirstLine '['];
    for i=1:length(OutputNames)
        if i~=1
            FHFirstLine = [FHFirstLine ','];
        end
        FHFirstLine = [FHFirstLine OutputNames{i}];
    end
    FHFirstLine = [FHFirstLine '] = '];
    if length(OutputNames)==1
        FHFirstLine = strrep(FHFirstLine,'[','');
        FHFirstLine = strrep(FHFirstLine,']','');
    end
end

% Name
FHFirstLine = [FHFirstLine FunctionName '('];

% Required positional inputs
if ~isempty(Required)
    if ~iscell(Required)
        Temp = Required;
        Required = [];
        Required{1} = Temp;
    end
    DefaultsBlock = '';
    ValidationBlock = '';
    ParserBlock = '';
    FHDescription = [FHDescription newline '%' ...
        newline '% Required inputs'];
    for i=1:length(Required)
        if i~=1
            FHFirstLine = [FHFirstLine ','];
            ReqNameCommaSep = [ReqNameCommaSep ','];
        end
        FHFirstLine = [FHFirstLine Required{i}];
        ReqNameCommaSep = [ReqNameCommaSep Required{i}];
        FHDescription = [FHDescription newline...
            '% ' Required{i} ' ... <VARIABLE DESCRIPTION>'];
        ValidationBlock = [ValidationBlock newline...
            'valid' Required{i} ' = @(x)true;'];
        ParserBlock = [ParserBlock newline...
            'addRequired(p,"' Required{i} '",valid' Required{i} ');'];
        FHCodeReassign = [FHCodeReassign newline...
            Required{i} ' = p.Results.' Required{i} ';'];
    end
    FHCode = [FHCode newline newline...
        '% Required inputs'...
        ValidationBlock...
        ParserBlock];
    ReqNameCommaSep = [ReqNameCommaSep ','];
end

% Name Value arguments
if ~isempty(NameValue)
    if ~iscell(NameValue)
        Temp = NameValue;
        NameValue = [];
        NameValue{1} = Temp;
    end
    DefaultsBlock = '';
    ValidationBlock = '';
    ParserBlock = '';
    FHDescription = [FHDescription newline '%' ...
        newline '% Name-Value pairs'];
    for i=1:length(NameValue)
        FHDescription = [FHDescription newline...
            '% "' NameValue{i} '" ... <NAMEVALUE DESCRIPTION>'];
        DefaultsBlock = [DefaultsBlock newline...
            'default' NameValue{i} ' = [];'];
        ValidationBlock = [ValidationBlock newline...
            'valid' NameValue{i} ' = @(x)true;'];
        ParserBlock = [ParserBlock newline...
            'addParameter(p,"' NameValue{i} '",default' NameValue{i} ...
            ',valid' NameValue{i} ');'];
        FHCodeReassign = [FHCodeReassign newline...
            NameValue{i} ' = p.Results.' NameValue{i} ';'];
    end
    FHCode = [FHCode newline newline...
        '% NameValue inputs'...
        DefaultsBlock...
        ValidationBlock...
        ParserBlock];
end

% Optional positional arguments
if ~isempty(OptionalPositional)
    if ~iscell(OptionalPositional)
        Temp = OptionalPositional;
        OptionalPositional = [];
        OptionalPositional{1} = Temp;
    end
    DefaultsBlock = '';
    ValidationBlock = '';
    ParserBlock = '';
    FHDescription = [FHDescription newline '%' ...
        newline '% Optional inputs'];
    for i=1:length(OptionalPositional)
        FHDescription = [FHDescription newline...
            '% ' OptionalPositional{i} ' ... <OPTIONAL POSITIONAL VARIABLE DESCRIPTION>'];
        DefaultsBlock = [DefaultsBlock newline...
            'default' OptionalPositional{i} ' = [];'];
        ValidationBlock = [ValidationBlock newline...
            'valid' OptionalPositional{i} ' = @(x)true;'];
        ParserBlock = [ParserBlock newline...
            'addOptional(p,"' OptionalPositional{i} '",default' OptionalPositional{i} ...
            ',valid' OptionalPositional{i} ');'];
        FHCodeReassign = [FHCodeReassign newline...
            OptionalPositional{i} ' = p.Results.' OptionalPositional{i} ';'];
    end
    FHCode = [FHCode newline newline...
        '% Optional positional inputs'...
        DefaultsBlock...
        ValidationBlock...
        ParserBlock];
end

if ~(isempty(OptionalPositional) && isempty(NameValue))
    if isempty(Required)
        FHFirstLine = [FHFirstLine 'varargin)'];
    else
        FHFirstLine = [FHFirstLine ',varargin)'];
    end
else
    FHFirstLine = [FHFirstLine ')'];
end
FHDescription = ['% ' FHFirstLine newline '%' newline FHDescription];
FHCode = [FHCode newline newline...
    'parse(p,' ReqNameCommaSep 'varargin{:});'];
FunctionHeader = [FHFirstLine newline... 
    FHDescription newline newline...
    FHCode newline newline...
    FHCodeReassign...
    newline newline];
try
    clipboard('copy',FunctionHeader)
    disp([newline 'Function header copied to clipboard!'...
        newline 'Insert into your editor with ctrl+v' newline])
catch
    warning(['Failed to add text to the clipboard!' newline 'Use direct output string instead'])
end

end