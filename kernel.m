function K = kernel(a,b,sigma,lambda)
N = length(a);
M = length(b);
K = zeros(N,M);
for i=1:N
    for j=1:M
        K(i,j) = sigma*exp(-1/(2*lambda)*norm(a(i)-b(j))^2);
    end
end
end
