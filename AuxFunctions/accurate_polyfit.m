function Parameters = accurate_polyfit(X,Y,n,varargin)
% function Parameters = accurate_polyfit(X,Y,n,varargin)
%
% Solves the polynomial fit directly by solving the ordinary least squares
% matrix equation. Parameters are ordered from highest polynomial degree to
% lowest. If isAffine is true or the polynomial is even or uneven the
% output will still be of length (n+1) and pad empty polynomials with zeros
% in order to work with the MATLAB built-in polyval.
%
% This function uses pinv() to solve the Vandermonde matrix, which makes it
% slower, than the MATLAB built-ins but overall more accurate.
%
% Required inputs
% X ... nx1 column vector
% Y ... nx1 column vector
% n ... Degree of fit
%
% Name-Value pairs
% "Symmetry" ... {'none','even','uneven'}
% "isAffine" ... Logical; if true there is only x-dependend terms: 
%                                               y = p1*x + p2*x^2 + ...

p = inputParser;
p.FunctionName = "simple_polyfit";
p.CaseSensitive = false;
p.PartialMatching = true;

% Required inputs
validX = @(x)ismatrix(x)&&size(x,2)==1;
validY = @(x)ismatrix(x)&&size(x,2)==1;
validn = @(x)isscalar(x)&&x>=0&&round(x)==x;
addRequired(p,"X",validX);
addRequired(p,"Y",validY);
addRequired(p,"n",validn);

% NameValue inputs
defaultSymmetry = 'none';
defaultisAffine = false;
validSymmetry = @(x)any(validatestring(x,{'none','even','uneven'}));
validisAffine = @(x)islogical(x);
addParameter(p,"Symmetry",defaultSymmetry,validSymmetry);
addParameter(p,"isAffine",defaultisAffine,validisAffine);

parse(p,X,Y,n,varargin{:});

% Assign parsing results to named variables
X = p.Results.X;
Y = p.Results.Y;
n = p.Results.n;
Symmetry = p.Results.Symmetry;
isAffine = p.Results.isAffine;

if length(X) ~= length(Y)
    error('X and Y need to be the same length');
end
% if length(X)<(n+1)
%     error('You need at least n+1 data points for an n-polyfit')
% end

% Set up Vandermonde Matrix
B = ones(length(X),n+1);
Parameters = zeros(n+1,1);
ParameterIndizes = [];
for i=n:-1:0
    if isAffine && i==0
        B(:,1) = [];
        continue
    end
    if isequal(Symmetry,'even') && mod(i,2)
        B(:,i+1) = [];
        continue
    elseif isequal(Symmetry,'uneven') && ~mod(i,2)
        B(:,i+1) = [];
        continue
    end
    B(:,i+1) = X.^i;
    ParameterIndizes = [ParameterIndizes i+1];
end

% Solve the Vandermonde Matrix
DenseParams = pinv(B'*B)*B'*Y;

Parameters(ParameterIndizes) = DenseParams(end:-1:1);

Parameters = Parameters(end:-1:1);

end