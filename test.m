p = 0.5;

% Dim = size(A);
% A = reshape(A,[],1);
% Z = zeros(size(A));
% for i=1:length(A)
%     if p < rand
%         Z(i) = A(i);
%     else
%         Z(i) = 0;
%     end
% end
% Z = reshape(Z,Dim);

Dim = size(X);
Z = zeros(size(X));
for i=1:Dim(1)
    for j=1:Dim(2)
        for k=1:Dim(3)
            if p < rand
                Z(i,j,k) = X(i,j,k)*1/(1-p);
            end
        end
    end
end