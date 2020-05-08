times = 100;
YError = zeros(times,2);
for i=1:times
    YError(i,:) = predict(DropoutNet,X(:,:,:,1));
end
mean(YError(:,1))
std(YError(:,1))