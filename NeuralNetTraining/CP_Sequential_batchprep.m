function [X,Y] = CP_Sequential_batchprep(objcell,SeqLen)
% This function takes as input a 1xN objcell of N objects of the class
% 'ForceMap' and combines them into the matrizes needed for NN-training
if nargin < 2
    SeqLen = 1024;
end
h = waitbar(0,'setting up');
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
X = zeros(SeqLen,1,1,Nimgs);
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

k = 1;
for i=1:Nmaps
    waitbar(i/Nmaps,h,sprintf('Batch Preparation for Force Map %i/%i',i,Nmaps));
    jRange = find(objcell{i}.selected_curves);
    for j=jRange'
        % Map the basedapp sequences to a length of SeqLen regardless of
        % initial length and normalize its values
        Norm_Force = (objcell{i}.basedapp{j}-min(objcell{i}.basedapp{j}))/range(objcell{i}.basedapp{j});
        Sequence = zeros(SeqLen,1);
        for m=1:SeqLen
            Idx_Prev = floor((m-1)/(SeqLen-1)*(length(objcell{i}.basedapp{j})-1)+1);
            Idx_Next = ceil((m-1)/(SeqLen-1)*(length(objcell{i}.basedapp{j})-1)+1);
            Dist_Prev = 1-(Idx_Prev-(m-1)/(SeqLen-1)*(length(objcell{i}.basedapp{j})-1));
            Dist_Next = -((m-1)/(SeqLen-1)*(length(objcell{i}.basedapp{j})-1)-(Idx_Next));
            Sequence(m) = (Dist_Next*Norm_Force(Idx_Prev)+Dist_Prev*Norm_Force(Idx_Next))/(1+(Idx_Next-Idx_Prev));
        end
        % Fill the output varables X and Y
        X(:,1,1,k) = Sequence;
        Y(k,1) = Norm_CP{i}(j,1);
        Y(k,2) = Norm_CP{i}(j,2);
        k = k + 1;
    end
end
close(h);
end