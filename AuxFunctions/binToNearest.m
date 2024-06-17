function binnedValues = binToNearest(values, binValues)
%function binnedValues = binToNearest(values, binValues)

    % Check the dimensions of binValues and ensure it's a row vector
    if size(binValues, 1) > 1
        binValues = binValues';
    end

    % Flatten the values array to a vector
    valuesVec = values(:);

    % Vectorized computation of the differences
    differences = abs(valuesVec - binValues);
    
    % Find the index of the minimum difference for each value
    [~, minIdx] = min(differences, [], 2);
    
    % Map the indices back to the bin values
    binnedValuesVec = binValues(minIdx);
    
    % Reshape the binned values back to the original shape of values
    binnedValues = reshape(binnedValuesVec, size(values));
end