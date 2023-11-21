function K = kernel3D(a,b,z,sigma,lambda)
N = length(a);
M = length(b);
O = length(z);
K = zeros(N,M);
for i=1:N
    for j=1:M
        for k=1:O
            K(i,j,k) = sigma*exp(-1/(2*lambda)*norm(a(i)-b(j))^2);
        end
    end
end
end