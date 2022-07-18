clear all;
traces = open('C:\Users\cha206\Anaconda3\envs\minian\raw_trace_0_1000.csv');
A = traces.data;
B = NaN([12 1000]);
id = [0:15];
id(:,4) = [];
id(:,2) = [];
id(:,12) = [];
id(:,12) = [];
z = 0;
for y=1:size(id,2)
    for x=1:size(A,1)
        if A(x,2) == id(y)
        B(y,x-z) = A(x,4);
        end
    end  
    z = z + 1000;   
end    

raw = B';
stackedplot(raw);   