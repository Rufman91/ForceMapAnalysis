p = 0.9;
A = X(:,:,:,1);

Dim = size(A);
A = reshape(A,[],1);
Z = zeros(size(A));
for i=1:length(A)
    if p < rand
        Z(i) = A(i);
    else
        Z(i) = 0;
    end
end
Z = reshape(Z,Dim);