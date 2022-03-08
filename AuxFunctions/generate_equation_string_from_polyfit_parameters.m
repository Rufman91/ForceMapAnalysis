function EquationString = generate_equation_string_from_polyfit_parameters(Parameters,varargin)
% function EquationString = generate_equation_string_from_polyfit_parameters(Parameters,varargin)
%
% Generates an equation or expression from polyfit parameters. The result
% is fit to be used as expression or equation with the symbolic toolbox 
% via the str2sym(EquationString) function.
%
%
% Required inputs
% Parameters ... nx1 or 1xn numerical list of polynomial coefficients.
%
% Name-Value pairs
% "Order" ... 'Descending'(def) Parameters go from high to low power terms
%             'Ascending' Parameters go from low to high power terms
% "Type" ... 'Expression' e.g. EquationString = '2+x+3*x'
%            'Equation'(def) e.g. EquationString = 'y == 2+x+3*x'
% "DependendVarName" ... Name of the dependend variable (def='y')
% "IndependentVarName" ... Name of the independend variable (def='x')

p = inputParser;
p.FunctionName = "generate_equation_string_from_polyfit_parameters";
p.CaseSensitive = false;
p.PartialMatching = true;

% Required inputs
validParameters = @(x)isnumeric(x)&&(size(x,1)==1||size(x,2)==1);
addRequired(p,"Parameters",validParameters);

% NameValue inputs
defaultOrder = 'Descending';
defaultType = 'Equation';
defaultDependendVarName = 'y';
defaultIndependentVarName = 'x';
validOrder = @(x)any(validatestring(x,{'Ascending','Descending'}));
validType = @(x)any(validatestring(x,{'Equation','Expression'}));
validDependendVarName = @(x)ischar(x);
validIndependentVarName = @(x)ischar(x);
addParameter(p,"Order",defaultOrder,validOrder);
addParameter(p,"Type",defaultType,validType);
addParameter(p,"DependendVarName",defaultDependendVarName,validDependendVarName);
addParameter(p,"IndependentVarName",defaultIndependentVarName,validIndependentVarName);

parse(p,Parameters,varargin{:});

% Assign parsing results to named variables
Parameters = p.Results.Parameters;
Order = p.Results.Order;
Type = p.Results.Type;
DependendVarName = p.Results.DependendVarName;
IndependentVarName = p.Results.IndependentVarName;

if isequal(Order,'Descending')
    Parameters = Parameters(end:-1:1);
end
if isequal(Type,'Equation')
    EquationString = [DependendVarName ' == '];
elseif isequal(Type,'Expression')
    EquationString = '';
end
for i=1:length(Parameters)
    if Parameters(i)==0
        continue
    end
    if (sign(Parameters(i))+1)
        Sign = '+';
    else
        Sign = '-';
    end
    if i==1
        EquationString = [EquationString num2str(Parameters(i))];
    elseif i==2
        EquationString = [EquationString Sign num2str(abs(Parameters(i))) '*' IndependentVarName];
    else
        EquationString = [EquationString Sign num2str(abs(Parameters(i))) '*' IndependentVarName '^' num2str(i-1)];
    end
end

end