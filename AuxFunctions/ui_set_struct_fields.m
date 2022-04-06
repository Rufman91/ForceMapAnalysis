function OutStruct = ui_set_struct_fields(InStruct,varargin)
% function OutStruct = ui_set_struct_fields(InStruct,varargin)
%
% Open UI window to set struct fields. What type of inpuit choice is given
% to the user depends on the data type of the fields values or, if given,
% by the contents of an associated field of the naming convention:
% '<FIELDNAME>' ---> 'DType<FIELDNAME>', e.g. InStruct.OptionFlag --->
% InStruct.DTypeOptionFlag Value='logical'
% If additionally there is a Field named 'Allowed<FIELDNAME>' the UI will
% show a popuplist to choose from the allowed values, e.g.:
% InStruct.Method = 'Default', InStruct.DTypeMethod = 'char',
% InStruct.AllowedMethod = {'Default','Other1','Other2'} will give a
% popupmenu with the 3 possibilities.
% If InStruct has fields that are themselves structs a pushbutton will be
% displayed that calls this very function on the substruct.
%
% Caution! This function does at the moment not allow for anything other
% than 1x1 dimensional inputs. Be it the structs themselves or the values
% of fields.
%
% Required inputs
% InStruct ... struct with fields
%
% Name-Value pairs
% "Title" ... character string, Title of the ui-option window
% "Message" ... character string, Message displayed in the top row of the ui-option window
% "NOptsPerPanel" ... <NAMEVALUE DESCRIPTION>
% "PixelsPerLine" ... <NAMEVALUE DESCRIPTION>
% "PixelsPerPanelWidth" ... <NAMEVALUE DESCRIPTION>

p = inputParser;
p.FunctionName = "ui_set_struct_fields";
p.CaseSensitive = false;
p.PartialMatching = true;

% Required inputs
validInStruct = @(x)true;
addRequired(p,"InStruct",validInStruct);

% NameValue inputs
defaultTitle = 'Options';
defaultMessage = '';
defaultNOptsPerPanel = 5;
defaultPixelsPerLine = 40;
defaultPixelsPerPanelWidth = 400;
validTitle = @(x)ischar(x);
validMessage = @(x)ischar(x);
validNOptsPerPanel = @(x)true;
validPixelsPerLine = @(x)true;
validPixelsPerPanelWidth = @(x)true;
addParameter(p,"Title",defaultTitle,validTitle);
addParameter(p,"Message",defaultMessage,validMessage);
addParameter(p,"NOptsPerPanel",defaultNOptsPerPanel,validNOptsPerPanel);
addParameter(p,"PixelsPerLine",defaultPixelsPerLine,validPixelsPerLine);
addParameter(p,"PixelsPerPanelWidth",defaultPixelsPerPanelWidth,validPixelsPerPanelWidth);

parse(p,InStruct,varargin{:});

% Assign parsing results to named variables
InStruct = p.Results.InStruct;
Title = p.Results.Title;
Message = p.Results.Message;
NOptsPerPanel = p.Results.NOptsPerPanel;
PixelsPerLine = p.Results.PixelsPerLine;
PixelsPerPanelWidth = p.Results.PixelsPerPanelWidth;

OutStruct = InStruct;
AcceptChanges = false;
FontSize = round(PixelsPerLine*.3);

% First, determine the needed size of the options window. Each individual
% field takes up a fixed amount of space and options are grouped in modular
% panels that will fill up to max 4x4 panels filling first left to right
% then up to down. Panels contain up to NOptsPerPanel options
FieldNames = fieldnames(InStruct);
FieldNames((contains(FieldNames,'DType'))|...
    (contains(FieldNames,'Allowed'))) = [];
NFields = length(FieldNames);
NPanels = ceil(NFields/NOptsPerPanel);
NRows = ceil(NPanels/4);
NColumns = min(NPanels,4);

% set up figure
PixelsPerPanelHeight = PixelsPerLine*NOptsPerPanel;
ScreenSize = get(0,'ScreenSize');
ScreenSize(1:2) = [];
WindowWidth = NColumns*PixelsPerPanelWidth;
WindowHeight = NRows*PixelsPerPanelHeight + 2*PixelsPerLine;
NLines = NRows*NOptsPerPanel+2;
Left = (ScreenSize(1)-WindowWidth)/2;
Bottom = (ScreenSize(2)-WindowHeight)/2;
Position = [Left Bottom WindowWidth WindowHeight];
F = figure('Name',Title,...
    'Units','pixels',...
    'Position',Position...
    );

DefUI(1) = uicontrol('Style','text',...
    'String',Message,...
    'Units','normalized',...
    'Position',[0 (NLines-1)/NLines 1 1/NLines],...
    'FontSize',FontSize);
DefUI(2) = uicontrol('Style','pushbutton',...
    'String','Accept',...
    'Units','normalized',...
    'Position',[(1/3-.15/NColumns) .05/NRows .3/NColumns .9/NLines],...
    'FontSize',FontSize,...
    'Callback',@accepted_changes);
DefUI(3) = uicontrol('Style','pushbutton',...
    'String','Cancel',...
    'Units','normalized',...
    'Position',[(2/3-.15/NColumns) .05/NRows .3/NColumns .9/NLines],...
    'FontSize',FontSize,...
    'Callback',@canceled_options);

% Loop over all fields and create their edit lines
m = 1;
for i=1:NRows
    for j=1:NColumns
        for k=1:NOptsPerPanel
            if m > NFields
                m = m+1;
                continue
            end
            if isfield(InStruct,['DType' FieldNames{m}])
                CurrentDType = InStruct.(['DType' FieldNames{m}]);
            elseif ~isempty(InStruct.(FieldNames{m}))
                CurrentDType = class(InStruct.(FieldNames{m}));
            else
                CurrentDType = 'double';
            end
            % Determine the position of the ui element
            TwoElementSpacer = .45/NColumns;
            CurrentLeft = (j-1+.1)/NColumns;
            CurrentBottom = 1 - (i+.25)/NLines - (k-1)*NRows/NLines;
            CurrentWidth = .35/NColumns;
            CurrentHeight = .5/NLines;
            CurrentTextPosition = [CurrentLeft CurrentBottom CurrentWidth CurrentHeight];
            CurrentUIPosition = [CurrentLeft+TwoElementSpacer CurrentBottom CurrentWidth CurrentHeight];
            % If there are allowed values, create a popup list
            if isfield(InStruct,['Allowed' FieldNames{m}])
                UI(2*(m-1)+1) = uicontrol('Style','text',...
                    'String',FieldNames{m},...
                    'Tooltip',FieldNames{m},...
                    'Units','normalized',...
                    'Position',CurrentTextPosition,...
                    'FontSize',FontSize,...
                    'HorizontalAlignment','right');
                UI(2*m) = uicontrol('Style','popupmenu',...
                    'String',InStruct.(['Allowed' FieldNames{m}]),...
                    'Units','normalized',...
                    'Position',CurrentUIPosition,...
                    'FontSize',FontSize);
            elseif any(matches({'logical'},CurrentDType))
                UI(2*(m-1)+1) = uicontrol('Style','text',...
                    'String',FieldNames{m},...
                    'Tooltip',FieldNames{m},...
                    'Units','normalized',...
                    'Position',CurrentTextPosition,...
                    'FontSize',FontSize,...
                    'HorizontalAlignment','right');
                UI(2*m) = uicontrol('Style','radiobutton',...
                    'String','',...
                    'Units','normalized',...
                    'Position',CurrentUIPosition,...
                    'FontSize',FontSize);
            else
                UI(2*(m-1)+1) = uicontrol('Style','text',...
                    'String',FieldNames{m},...
                    'Tooltip',FieldNames{m},...
                    'Units','normalized',...
                    'Position',CurrentTextPosition,...
                    'FontSize',FontSize,...
                    'HorizontalAlignment','right');
                UI(2*m) = uicontrol('Style','edit',...
                    'String',InStruct.(FieldNames{m}),...
                    'Units','normalized',...
                    'Position',CurrentUIPosition,...
                    'FontSize',FontSize);
            end
            m = m + 1;
        end
    end
end

    function canceled_options(varargin)
        AcceptChanges = false;
        close(F)
    end
    function accepted_changes(varargin)
        AcceptChanges = true;
        close(F)
    end
uiwait(F)

if AcceptChanges
    OutStruct = InStruct;
end

end
