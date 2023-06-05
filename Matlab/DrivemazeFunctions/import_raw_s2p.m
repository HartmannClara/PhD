function raw_f_res=import_raw_s2p(filename)
    %UNTITLED10 Import raw
    %   Imports F.csv file or Fneu.csv file from s2p
    %
    Results = readtable(filename);
    Results = table2array(Results);
    raw_f_res=Results';
end