function [ContainsSame,IsSame] = check_structs_for_same_fields_recursively(CheckedStruct,ComparedStruct,FirstCall)
% 

if nargin < 3
    FirstCall = true;
end

if FirstCall
    CheckedContainsCompared = check_structs_for_same_fields_recursively(CheckedStruct,ComparedStruct,false);
    ComparedContainsChecked = check_structs_for_same_fields_recursively(ComparedStruct,CheckedStruct,false);
    ContainsSame = CheckedContainsCompared;
    if ContainsSame && ComparedContainsChecked
        IsSame = true;
    else
        IsSame = false;
    end
    return
end

% first check if dimensions are the same and flatten structs to nx1 vector
if ~isequal(size(CheckedStruct),size(ComparedStruct))
    ContainsSame = false;
    IsSame = false;
    return
end
CheckedStruct = reshape(CheckedStruct,[],1);
ComparedStruct = reshape(ComparedStruct,[],1);

CheckedFieldNames = fieldnames(CheckedStruct);
ComparedFieldNames = fieldnames(ComparedStruct);
L = length(ComparedFieldNames);
StructL = length(ComparedStruct);

IsSameField = false(L,1);
for k=1:StructL
    for i=1:L
        for j=1:length(CheckedFieldNames)
            if isequal(ComparedFieldNames{i},CheckedFieldNames{j})
                if isstruct(ComparedStruct(k).(ComparedFieldNames{i})) &&...
                        isstruct(CheckedStruct(k).(CheckedFieldNames{j}))
                    IsSameField(i) = check_structs_for_same_fields_recursively(...
                        CheckedStruct(k).(CheckedFieldNames{j}),ComparedStruct(k).(ComparedFieldNames{i}));
                    break
                elseif isstruct(ComparedStruct(k).(ComparedFieldNames{i})) &&...
                        ~isstruct(CheckedStruct(k).(CheckedFieldNames{j}))
                    break
                else
                    IsSameField(i) = true;
                    break
                end
            end
        end
    end
end

ContainsSame = all(IsSameField);
IsSame = true;

end