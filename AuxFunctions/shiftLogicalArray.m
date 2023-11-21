function shiftedArray = shiftLogicalArray(logicalArray, ShiftPixX, ShiftPixY)
    % Initialize the resulting shifted array to all false
    shiftedArray = false(size(logicalArray));

    % Find the indices of all true entries in the logical array
    [rows, cols] = find(logicalArray);

    % Apply the shift
    shiftedRows = rows + ShiftPixX;
    shiftedCols = cols + ShiftPixY;

    % Filter out indices that are out of the bounds after the shift
    validInd = (shiftedRows > 0) & (shiftedRows <= size(logicalArray, 1)) & ...
               (shiftedCols > 0) & (shiftedCols <= size(logicalArray, 2));

    % Update the shifted array with the valid shifted indices
    shiftedArray(sub2ind(size(logicalArray), shiftedRows(validInd), shiftedCols(validInd))) = true;
end
