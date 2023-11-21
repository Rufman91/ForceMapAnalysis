function clearAllExcept(varargin)
    % clearAllExcept Clears all variables in the base workspace except for the specified variables.
    %   clearAllExcept(varName1, varName2, ...) takes variable names as input and clears all other
    %   variables in the base workspace.

    allVars = evalin('base', 'who');
    varsToKeep = varargin;
    varsToClear = setdiff(allVars, varsToKeep);
    for i = 1:length(varsToClear)
        evalin('base', ['clear ' varsToClear{i}]);
    end
end
