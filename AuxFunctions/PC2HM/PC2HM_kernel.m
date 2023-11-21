function K = PC2HM_kernel(a,b,sigma,lambda)
    % Check whether expected matrix size is gonna make the device explode and
    % partition accordingly
    M = size(a,1);
    N = size(b,1);
    O = size(a,2);
    
    if O~=size(b,2)
        error('Vector a and b need to have the same amount of columns!')
    end
    
    % Vectorize that shit!
    C = zeros(M,N,O);
    for i=1:O
        C(:,:,i) = a(:,i) - b(:,i)';
    end
    K = sigma*exp(-1/(2*lambda)*sum(C.^2, 3));
end