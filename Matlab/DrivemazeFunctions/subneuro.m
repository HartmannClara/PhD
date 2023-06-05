function subneuro_trace = subneuro(raw_f_res,Fneu,factor)
%UNTITLED11 subneuro
%   subtracts Fneuro traces ( multiplied by factor ) from raw F traces, 
%   takes in raw traces and gives out neuropil subtracted traces


for i = 1:size(raw_f_res,2)
    subneuro_trace(:,i) = raw_f_res(:,i) - (factor* Fneu(:,i));
end    

