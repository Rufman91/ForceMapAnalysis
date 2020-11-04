function show_images_with_cp(X,Y,Idxs)
% show_images_with_cp(X,Y,Idxs)
%
% X    ... 4-D Array with dimensions [ImgHeight ImgWidth ChannelNumber NImages]
% Y    ... NImages x 2 x- and y- coordinates of the contact point on the image
% Idxs ... vector containing 3 indizes
if nargin < 3
    RandPos = randi(length(X),3,1);
else
    RandPos = Idxs;
end
ImgSize = length(X(:,1,1,1));

for j = 1:3
    subplot(3,3,(1+3*(j-1)))
    imshow(X(:,:,1,RandPos(j)))
    drawpoint('Position',[Y(RandPos(j),1) (1 - Y(RandPos(j),2))]*ImgSize);
    title(sprintf('Image Nr. %i',RandPos(j)))
    subplot(3,3,(2+3*(j-1)))
    imshow(X(:,:,2,RandPos(j)))
    subplot(3,3,(3+3*(j-1)))
    imshow(X(:,:,3,RandPos(j)))
end

end