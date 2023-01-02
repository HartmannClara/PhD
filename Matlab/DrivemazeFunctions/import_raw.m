function raw_f_res=import_raw(filename,cells)
    %UNTITLED10 Import raw
    %   Imports raw results file (made with imagej); csv file with Area
    %   Mean Min Max and 1st column frame
    %
    opts = delimitedTextImportOptions("NumVariables", (cells*4)+1);
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";
    %cells=cells;
    VariableNames = ["VarName1"];
    VariableTypes = ["double"];
    for i = 1:cells
        newVariableNames= ["Area"+ i , "Mean"+ i, "Min"+ i, "Max"+ i];
        B = convertStringsToChars(newVariableNames);
        VariableNames = [VariableNames, B];
    
        newVariableTypes = ["double", "double", "double", "double"];
        C = convertStringsToChars(newVariableTypes);
        VariableTypes = [VariableTypes, C];
    end    
    opts.VariableNames = VariableNames;
    opts.VariableTypes = VariableTypes;
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    Results = readtable(filename, opts);
    Results = table2array(Results);
    raw_f_res=Results(:,3:4:end);
end