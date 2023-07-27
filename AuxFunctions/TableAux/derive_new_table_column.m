function T = derive_new_table_column(T, newColumnName, columnCell, expression)

    % Check if new column name already exists
    if ismember(newColumnName, T.Properties.VariableNames)
        error('The new column name already exists in the table.');
    end

    % Convert cell array to table
    for i = 1:length(columnCell)
        ReplacedName = ['c' num2str(i)];
        NewName = ['columnCell' '{' num2str(i) '}'];
        expression = strrep(expression,ReplacedName,NewName);
    end
    
%     % Check if the columns in the expression are numeric
%     for i = 1:length(columnList)
%         columnName = ['c' num2str(i)];
%         if ~isnumeric(calcTable.(columnName))
%             error('One or more of the columns are not of numeric data type.');
%         end
%     end

    % Calculate the new column
    newColumn = eval(expression);

    % Add the new column to the table
    T.(newColumnName) = newColumn;
end
