function [X,Y] = CP_CNN_batchprep(objcell,ImgSizeFinal,ImgSize,CutPercent)
% This function takes as input a 1xN objcell of N objects of the class
% 'ForceMap' and combines them into the matrizes needed for NN-training
if nargin < 2
    ImgSizeFinal = 128;
end
Nmaps = length(objcell);
% First check, if manual CPs are already given in the input class instance. If
% not, start manual CP selection.
% for i=1:Nmaps
%     if isempty(objcell{i}.Man_CP)
%         objcell{i}.estimate_cp_manually();
%     end
% end

% Preallocate sizes of the outputvariables. The trainNetwork() function
% requires the batch X of predictor varables (in this case the force-images)
% to be a height-by-width-by-channelnumber-by-Numberofimages array and the
% prelabeled regression responses Y to be a
% Numberofimages-by-Numberofresponses array
Nimgs = 0;
for i=1:Nmaps
    Nimgs = Nimgs + sum(objcell{i}.SelectedCurves);
end
X = objcell{1}.CP_batchprep_3_channel(objcell,ImgSizeFinal,ImgSize,CutPercent);
Y = zeros(Nimgs,2);



% Transform the values for the CP such that they correspond to the correct
% positions in the force-images.
k = 1;
Norm_CP = cell(Nmaps,1);
for i=1:Nmaps
    jRange = find(objcell{i}.SelectedCurves);
    Norm_CP{i} = zeros(length(objcell{i}.BasedApp),2);
    for j=jRange'
        Norm_CP{i}(j,1) = (objcell{i}.Man_CP(j,1)-min(objcell{i}.HHApp{j}))/...
            range(objcell{i}.HHApp{j});
        Norm_CP{i}(j,2) = (objcell{i}.Man_CP(j,2)-min(objcell{i}.BasedApp{j}))/...
            range(objcell{i}.BasedApp{j});
        Y(k,1) = Norm_CP{i}(j,1);
        Y(k,2) = Norm_CP{i}(j,2);
        k = k + 1;
    end
end


end