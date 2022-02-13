function LookupTable = numerically_solve_bivariate_equation_into_lookup_table(EquationString,DependentName,IndependentName,varargin)
% function LookupTable = numerically_solve_bivariate_equation_into_lookup_table(EquationString,DependentName,IndependentName,varargin)
%
% Solve a bivariate equation numerically within 'RangeIndependent' for 'NumSamples'
% points. DependendName and IndependentName refer to variables purpose
% AFTER solving, so referring to the output, e.g.:
% EquationString='y=x*sin(x)' and you want to solve for x, then
% DependentName='x' and IndependentName='y'.
% Equations can have arbitrarily many solutions: this function chooses the
% one closest to initial guess on the first loop and the does the following:
% After each individual solution, 'InitialGuess' is replaced with
% the previous solution and given the equation be continuous and
% 'NumSamples' sufficiently large, the resulting curve will stay within one
% branch of solutions as long as there exist solutions in the branch and
% jump to other branches otherwise.
%
% Solutions are not guaranteed and may lead to (partially) NaN outputs.
%
%
% Required inputs
% EquationString ... Character string containing the equation to be solved.
%                   Form should contain a '=' or'==' and exactly two named
%                   variables e.g.: '1/(1-exp(x))=log(y)'
% DependentName ... character string of the dependend variable e.g.:
%                   'y','x','C', etc.
% IndependentName ... character string of the independend variable e.g.:
%                   'y','x','C', etc.
%
% Name-Value pairs
% "RangeIndependent" ... Numeric range of the independent variable as 2x1
%                       or 1x2 real valued array with
%                       RangeIndependent(1)<RangeIndependent(2)
% "RangeDependent" ... Numeric range of the independent variable as 2x1
%                       or 1x2 real valued array with
%                       RangeIndependent(1)<RangeIndependent(2). being
%                       specific about this range can help with staying on
%                       the desired solution branch.
% "NumSamples" ... Number of evenly spaced points in RangeIndependend that
%                   the equation is solved for. More points take longer to
%                   compute but provide greater accuracy and can help with
%                   staying on one consistent solution branch
% "InitialGuess" ... Initial guess for the solution of the point
%                   RangeIndependent(1)
% "Verbose" ... Additionally draw some figures for the results

p = inputParser;
p.FunctionName = "numerically_solve_bivariate_equation_into_lookup_table";
p.CaseSensitive = false;
p.PartialMatching = true;

% Required inputs
validEquationString = @(x)ischar(x);
validDependentName = @(x)ischar(x);
validIndependentName = @(x)ischar(x);
addRequired(p,"EquationString",validEquationString);
addRequired(p,"DependentName",validDependentName);
addRequired(p,"IndependentName",validIndependentName);

% NameValue inputs
defaultRangeIndependent = [0 1];
defaultRangeDependent = [-inf inf];
defaultNumSamples = 100;
defaultInitialGuess = 0;
defaultVerbose = false;
validRangeIndependent = @(x)isnumeric(x)&&isequal(size(x),[1 2]);
validRangeDependent = @(x)isnumeric(x)&&isequal(size(x),[1 2]);
validNumSamples = @(x)isnumeric(x)&&isscalar(x);
validInitialGuess = @(x)isnumeric(x)&&isscalar(x);
validVerbose = @(x)islogical(x);
addParameter(p,"RangeIndependent",defaultRangeIndependent,validRangeIndependent);
addParameter(p,"RangeDependent",defaultRangeDependent,validRangeDependent);
addParameter(p,"NumSamples",defaultNumSamples,validNumSamples);
addParameter(p,"InitialGuess",defaultInitialGuess,validInitialGuess);
addParameter(p,"Verbose",defaultVerbose,validVerbose);

parse(p,EquationString,DependentName,IndependentName,varargin{:});

% Assign parsing results to named variables
EquationString = p.Results.EquationString;
DependentName = p.Results.DependentName;
IndependentName = p.Results.IndependentName;
RangeIndependent = p.Results.RangeIndependent;
RangeDependent = p.Results.RangeDependent;
NumSamples = p.Results.NumSamples;
InitialGuess = p.Results.InitialGuess;
Verbose = p.Results.Verbose;


% Define numerical values of independent var range
NumericalIndependent = linspace(RangeIndependent(1),RangeIndependent(2),NumSamples)';

% Check whether independent var is on one side only e.g.
% Independent==g(Dependent) and if so, switch the variables, insert the
% range of values
FoundAllSolutions = true;
if 0
    % TODO
else
    % Loop over numerical values, insert them into the equation and solve
    % numerically
    NumericalDependent = zeros(size(NumericalIndependent));
    for i=1:length(NumericalIndependent)
        EquationToSolve = strrep(EquationString,IndependentName,num2str(NumericalIndependent(i)));
        SymEquation = str2sym(EquationToSolve);
        SymToSolve = sym(DependentName);
        NumericalSolutions = vpasolve(SymEquation,SymToSolve);
        NumericalSolutions = double(NumericalSolutions);
        if isempty(NumericalSolutions)
            NumericalDependent(i) = nan;
            FoundAllSolutions = false;
            continue
        end
        NumericalSolutions((real(NumericalSolutions)<RangeDependent(1))|...
            (real(NumericalSolutions)>RangeDependent(2))) = [];
        if isempty(NumericalSolutions)
            NumericalDependent(i) = nan;
            FoundAllSolutions = false;
            continue
        end
        [~,NearInitIndex] = min(abs(NumericalSolutions-InitialGuess));
        NumericalDependent(i) = real(NumericalSolutions(NearInitIndex));
        InitialGuess = NumericalSolutions(NearInitIndex);
    end
end

LookupTable = [NumericalIndependent NumericalDependent];

if Verbose
    figure('Color','w','Name',p.FunctionName)
    plot(LookupTable(:,1),LookupTable(:,2))
    xlabel(IndependentName)
    ylabel(DependentName)
    legend(EquationString)
end

if ~FoundAllSolutions
    warning([newline ...
        'Could not find numerical solutions for some or all of the requested points.'...
        newline 'Try adjusting requested ranges and the initial guess.' ...
        newline 'Replaced missing solutions with NaN'])
end

end