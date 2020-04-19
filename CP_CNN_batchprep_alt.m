function [X,Y] = CP_CNN_batchprep_alt(objcell,ImgSize)
% This function takes as input a 1xN objcell of N objects of the class
% 'ForceMap' and combines them into the matrizes needed for NN-training
if nargin < 2
    ImgSize = 128;
end
Nmaps = length(objcell);
% First check, if manual CPs are already given in the input class instance. If
% not, start manual CP selection.
for i=1:Nmaps
    if isempty(objcell{i}.Man_CP)
        objcell{i}.manual_CP();
    end
end

% Preallocate sizes of the outputvariables. The trainNetwork() function
% requires the batch X of predictor varables (in this case the force-images)
% to be a height-by-width-by-channelnumber-by-Numberofimages array and the
% prelabeled regression responses Y to be a
% Numberofimages-by-Numberofresponses array
Nimgs = 0;
for i=1:Nmaps
    Nimgs = Nimgs + sum(objcell{i}.selected_curves);
end
X = zeros(ImgSize,ImgSize,1,Nimgs);
Y = zeros(Nimgs,2);

% Transform the values for the CP such that they correspond to the correct
% positions in the force-images.
Norm_CP = cell(Nmaps,1);
for i=1:Nmaps
    jRange = find(objcell{i}.selected_curves);
    Norm_CP{i} = zeros(length(objcell{i}.basedapp),2);
    for j=jRange'
        Norm_CP{i}(j,1) = (objcell{i}.Man_CP(j,1)-min(objcell{i}.thapp{j}))/...
            range(objcell{i}.thapp{j});
        Norm_CP{i}(j,2) = (objcell{i}.Man_CP(j,2)-min(objcell{i}.basedapp{j}))/...
            range(objcell{i}.basedapp{j});
    end
end


fig = figure('Color','w');
k = 1;
for i=1:Nmaps
    jRange = find(objcell{i}.selected_curves);
    for j=jRange'
        % Save the plots as images and convert them into cropped [0
        % 1]-range grayscale images
        fig.Name = sprintf('%s curve nr.%i',objcell{i}.name,j);
        plot(objcell{i}.thapp{j},objcell{i}.basedapp{j},'color','black');
        axis([min(objcell{i}.thapp{j}) max(objcell{i}.thapp{j}) min(objcell{i}.basedapp{j}) max(objcell{i}.basedapp{j})])
        axis off
        graytest = imcomplement(rgb2gray(frame2im(getframe(gcf))));
        grayscaled = double(graytest)/255;
        % Crop off the rows and columns of the image only containing zeros
        crop = [1 size(grayscaled,2) 1 size(grayscaled,1)];
        while sum(grayscaled(:,crop(1)))== 0
            crop(1) = crop(1) + 1;
        end
        crop(1) = crop(1) - 1;
        while sum(grayscaled(:,crop(2)))== 0
            crop(2) = crop(2) - 1;
        end
        crop(2) = crop(2) + 1;
        while sum(grayscaled(crop(3),:))== 0
            crop(3) = crop(3) + 1;
        end
        crop(3) = crop(3) - 1;
        while sum(grayscaled(crop(4),:))== 0
            crop(4) = crop(4) - 1;
        end
        crop(4) = crop(4) + 1;
        grayscaled(:,[1:crop(1) crop(2):size(grayscaled,2)]) = [];
        grayscaled([1:crop(3) crop(4):size(grayscaled,1)],:) = [];
        grayfinal = imresize(grayscaled,[ImgSize ImgSize],'bilinear');
        % Fill the output varables X and Y
        X(:,:,1,k) = grayfinal;
        Y(k,1) = Norm_CP{i}(j,1);
        Y(k,2) = Norm_CP{i}(j,2);
        k = k + 1;
    end
end
close(fig);

end