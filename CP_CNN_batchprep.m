function [X,Y] = CP_CNN_batchprep(objcell,ImgSize)
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
    k = 1;
    Norm_CP{i} = zeros(sum(objcell{i}.selected_curves),2);
    for j=jRange'
        Norm_CP{i}(k,1) = (objcell{i}.Man_CP(j,1)-min(objcell{i}.thapp{j}))/...
            range(objcell{i}.thapp{j})*length(objcell{i}.thapp{j})/objcell{i}.header{2,2};
        Norm_CP{i}(k,2) = (objcell{i}.Man_CP(j,2)-min(objcell{i}.basedapp{j}))/...
            range(objcell{i}.thapp{j})*length(objcell{i}.basedapp{j})/objcell{i}.header{2,2};
        k = k + 1;
    end
end
        
% Fill the outputvariables X and Y
k = 1;
for i=1:Nmaps
    img = objcell{i}.force2img(ImgSize);
    for j=1:length(img)
        X(:,:,1,k) = img{j};
        Y(k,1) = Norm_CP{i}(j,1);
        Y(k,2) = Norm_CP{i}(j,2);
        k = k + 1;
    end
end

end