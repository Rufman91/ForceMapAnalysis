function [Editbox,Listbox] = create_autocomplete_edit_field(ACDictionary,varargin)
% function [Editbox,Listbox] = create_autocomplete_edit_field(ACDictionary,varargin)
%
% <FUNCTION DESCRIPTION HERE>
%
%
% Required inputs
% ACDictionary ... <VARIABLE DESCRIPTION>
%
% Name-Value pairs
% "Parent" ... <NAMEVALUE DESCRIPTION>
% "Arrangement" ... <NAMEVALUE DESCRIPTION>
% "Units" ... <NAMEVALUE DESCRIPTION>
% "Position" ... <NAMEVALUE DESCRIPTION>
% "MatchCase" ... <NAMEVALUE DESCRIPTION>
% "RelativeEditPosition" ... <NAMEVALUE DESCRIPTION>

p = inputParser;
p.FunctionName = "create_autocomplete_edit_field";
p.CaseSensitive = false;
p.PartialMatching = true;

% Required inputs
validACDictionary = @(x)true;
addRequired(p,"ACDictionary",validACDictionary);

% NameValue inputs
defaultParent = gcf;
defaultArrangement = 'Vertical';
defaultUnits = 'normalized';
defaultPosition = [.2 .2 .4 .4];
defaultMatchCase = false;
defaultRelativeEditPosition = 0;
validParent = @(x)true;
validArrangement = @(x)any(validatestring(x,{'Vertical','Horizontal'}));
validUnits = @(x)any(validatestring(x,{'Fraction','Meters'}));
validPosition = @(x)(size(x,1)==1 & size(x,2)==4 & isnumeric(x));
validMatchCase = @(x)(islogical(x) | (isscalar(x) & isinteger(x)));
validRelativeEditPosition = @(x)true;
addParameter(p,"Parent",defaultParent,validParent);
addParameter(p,"Arrangement",defaultArrangement,validArrangement);
addParameter(p,"Units",defaultUnits,validUnits);
addParameter(p,"Position",defaultPosition,validPosition);
addParameter(p,"MatchCase",defaultMatchCase,validMatchCase);
addParameter(p,"RelativeEditPosition",defaultRelativeEditPosition,validRelativeEditPosition);

parse(p,ACDictionary,varargin{:});

% Assign parsing results to named variables
ACDictionary = p.Results.ACDictionary;
Parent = p.Results.Parent;
Arrangement = p.Results.Arrangement;
Units = p.Results.Units;
Position = p.Results.Position;
MatchCase = p.Results.MatchCase;
RelativeEditPosition = min(p.Results.RelativeEditPosition,1);
RelativeEditPosition = max(RelativeEditPosition,0);

ListPos = Position;
EditPos = Position;
Height = Position(4);
Width = Position(3);

if isequal(Arrangement,'Vertical')
    ListPos(4) = Height*9/10;
    EditPos(4) = Height/10;
    EditPos(2) = Position(2) + ListPos(4);
elseif isequal(Arrangement,'Horizontal')
    ListPos(4) = Height*9/10;
    EditPos(4) = Height/10;
    ListPos(3) = Width/2;
    EditPos(3) = Width/2;
    ListPos(1) = Position(1) + EditPos(3);
    EditPos(2) = Position(2) + RelativeEditPosition*(ListPos(4) - EditPos(4));
end

Editbox = uicontrol(Parent,...
    'style','edit',...
    'String','',...
    'units',Units,...
    'position',EditPos);

Listbox = uicontrol(Parent,...
    'Style','listbox',...
    'String',ACDictionary,...
    'Units',Units,...
    'Position',ListPos);


set(Editbox,'KeyPressFcn',{@update_suggestion_list,MatchCase,Listbox,ACDictionary})
set(Editbox,'ButtonDownFcn',{@activate_edit_list_pair,Listbox})
set(Editbox,'Enable', 'Inactive');
set(Listbox,'Callback',{@autofill_edit_field_from_list,Editbox})

end