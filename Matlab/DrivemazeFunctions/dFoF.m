function subneuro_trace = subneuro(raw_f_res,Fneu)
%UNTITLED11 dFoF
%   calculates deltaF over F0 with F0 calculated from a moving median (window size 120s)
%   takes in raw traces and gives out dFoF traces(df_f_trace)
medians = [];
    for ROI = 1:size(raw_f_res,2)
        cell_n = [];
        cell_n = raw_f_res(:,ROI); 
        medians_cell = movmedian(cell_n, 600);%window size in frames %120s(1200)
        medians = [medians, medians_cell];     
    end
    subneuro_trace = [];
    for cell_median = 1: size(medians,2) % Loop df/f
       v = [];
       d = [];
       v = medians(:,cell_median);%select column from medians
       d = raw_f_res(:,cell_median);%select column from raw_f_res
       for f = 1: size(d,1)
           a = [];
           a =  v(f) ;
           b = (d(f) - 0.9*a)/0.9*a;%0.7
           subneuro_trace(f,cell_median) = b;
       end
    end
end