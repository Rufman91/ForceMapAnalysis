function [mu_star,sigma_star] = predictGP_mean(X,X_star,sigma,lambda,F,noise, varargin)
% Set default value for the Name-Value input argument
CalculateSigma = false;

% Process Name-Value input arguments
for i = 1:2:length(varargin)
    name = varargin{i};
    value = varargin{i+1};
    if strcmpi(name, 'CalculateSigma')
        CalculateSigma = value;
    end
end

FMean = mean(F);
FSTD = std(F);
FNorm = (F - FMean)/(FSTD+1e-10);
N = size(X,1);
K = kernel(X,X,sigma,lambda) + noise*eye(N);
K_star = kernel(X,X_star,sigma,lambda);
mu_star_norm = K_star'*(K\FNorm);

% Check if sigma_star should be calculated
if CalculateSigma
    sigma_star = kernel(X_star,X_star,sigma,lambda) - K_star'*(K\K_star);
else
    sigma_star = [];
end

mu_star = (mu_star_norm*FSTD + FMean);
end