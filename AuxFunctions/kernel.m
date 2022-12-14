function K = kernel(a,b,sigma,lambda)



% Check wheter expected matrix size is gonna make the device explode and
% partition accordingly
M = size(a,1);
N = size(b,1);
O = size(a,2);

if O~=size(b,2)
    error('Vector a and b need to have the same amount of columns!')
end

% % Vectorize that shit!
for i=1:O
    C(:,:,i) = a(:,i) - b(:,i)';
end
K = sigma*exp(-1/(2*lambda)*vecnorm(C,2,3).^2);

end
