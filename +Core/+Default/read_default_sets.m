function [myStruct] = read_default_sets(setstxt,default)
    % Puts data from setstxt directly to structure

    if nargin < 2
        default = 1;
    end

    myStruct = struct();

    if default
        mFilePath = mfilename('fullpath');
        [threeLevelsUpDir, ~] = fileparts(fileparts(fileparts(mFilePath)));
        setstxt  =fullfile(threeLevelsUpDir,setstxt);
    end

    setsTable  = readtable(setstxt,'Format','%s%s%s');
    
    cellNames = {};
    for i = 1:size(setsTable,1)
        validFieldName = split(matlab.lang.makeValidName(setsTable.Var3{i}),'_');
        number = str2num(setsTable.Var1{i});
        if ~isnan(number)
            val = number;
        else
            val = strtrim(setsTable.Var1{i});
        end
        if length(validFieldName)>=3
            cellNames = [cellNames, validFieldName{1}];
            myStruct.(validFieldName{1}){str2num(validFieldName{2})} = val; % Example value assignment
        else
            myStruct.(validFieldName{1}) = val; % Example value assignment
        end
    end
    uniqueNames = unique(cellNames);
    for i=1:length(uniqueNames)
        myStruct.(uniqueNames{i}) = cell2mat(myStruct.(uniqueNames{i}));
    end

    



end

