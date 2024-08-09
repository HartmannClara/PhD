function events=import_DM_events(filename)
    %UNTITLED10 Import csv
    %   Imports events file (csv) from drivemaze as a table
    opts = delimitedTextImportOptions("NumVariables", 11);
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";
    opts.VariableNames = ["Date_Time","amount consumed","latency to consumption", "Type", "hardware_time", "experiment", "Pellet_number", "Latency", "Laser","Drop_Number", "Drop_latency"];
    opts.VariableTypes = ["string", "double", "double" "categorical", "double", "categorical", "double", "double", "double", "double", "double"];
    opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
    opts = setvaropts(opts, [1, 4], "EmptyFieldRule", "auto");
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    events = readtable(filename, opts);
    clear opts
end