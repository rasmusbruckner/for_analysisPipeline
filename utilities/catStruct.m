function mainStruct = catStruct(mainStruct, addStruct)
% Extract field names from structure
allNames = fieldnames(mainStruct);

% Cycle over field names and concatenate
for f = 1:length(allNames)
    fieldName = allNames{f};
    mainStruct.(fieldName) = cat(1, mainStruct.(fieldName), addStruct.(fieldName));
end
end