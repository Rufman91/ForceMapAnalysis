function OutMat = data_cell_to_mat(InCell)

if ~any(size(InCell)==1)
    error('InCell needs to be a nx1 or 1xn sized cell with nx1 or 1xn sized vectors as elements');
end
for i=1:length(InCell)
    if ~any(size(InCell{i})==1)
        error('InCell needs to be a nx1 or 1xn sized cell with nx1 or 1xn sized vectors as elements');
    end
end

L = length(InCell);
VecSizes = ones(1,L);
for i=1:L
    VecSizes(i) = length(InCell{i});
end

NRows = max(VecSizes);
OutMat = nan(NRows,L);

for i=1:L
    OutMat(:,i) = InCell{i};
end

end