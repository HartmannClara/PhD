function events=import_events(filename)
    %UNTITLED10 Import csv
    %   Imports events file (csv) from drivemaze as a table
    opts = delimitedTextImportOptions("NumVariables", 7);
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";
    opts.VariableNames = ["Date_Time", "amount_consumed", "latency_to_consumption", "Type", "frame", "hardware_time", "experiment"];
    opts.VariableTypes = ["string", "double", "double", "categorical", "double", "double", "string"];
    opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
    opts = setvaropts(opts, [1, 4], "EmptyFieldRule", "auto");
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    events = readtable(filename, opts);
    clear opts
end