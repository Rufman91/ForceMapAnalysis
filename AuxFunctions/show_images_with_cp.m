function NextIdxs = show_images_with_cp(X,Y,M,N,Idxs)
% show_images_with_cp(X,Y,M,N,Idxs)
%
% X    ... 4-D Array with dimensions [ImgHeight ImgWidth ChannelNumber NImages]
% Y    ... NImages x 2 x NPredictions X and Y coordinates of the contact point on
%          the image 
% M    ... Number of rows
% N    ... Number of columns
% Idxs ... vector containing M*N indizes
if nargin < 5
    RandPos = randi(length(X),M*N,1);
else
    RandPos = Idxs;
end
ImgSize = length(X(:,1,1,1));
try
    NPredictions = size(Y,3);
catch
    NPredictions = 1;
end

Color = {'r','g','b','y'};

tiledlayout(M,N,'TileSpacing','none','Padding','none')

for j = 1:M*N
    nexttile
    imshow(X(:,:,1,RandPos(j)))
    for k = 1:NPredictions
        drawpoint('Position',[Y(RandPos(j),1,k) (1 - Y(RandPos(j),2,k))]*ImgSize,'Color',Color{k});
    end
    title(sprintf('Image Nr. %i',RandPos(j)))
end

NextIdxs = RandPos + N*M;

end