clear all;
traces = open('C:\Users\cha206\Data\Pellet\g2\0804\YrA_afterfirstupdate.csv');
A = traces.data;
B = NaN([40 6040]);
id = [0:40];
%id(:,4) = [];
%id(:,2) = [];
%id(:,12) = [];
%id(:,12) = [];
z = 0;
for y=1:size(id,2)
    for x=1:size(A,1)
        if A(x,2) == id(y)
        B(y,x-z) = A(x,4);
        end
    end  
    z = z + 6040;   
end    

raw = B';
stackedplot(raw);   

%%
traces2 = open('C:\Users\cha206\Data\Pellet\g2\0604\YrA_aftersecupdate.csv');
C = traces2.data;
D = NaN([23 6042]);
id2 = [1:30];
id2(:,23) = [];
id2(:,20) = [];
id2(:,14) = [];
id2(:,5) = [];
id2(:,3) = [];
id2(:,2) = [];
id2(:,1) = [];

z = 0;
for y=1:size(id2,2)
    for x=1:size(C,1)
        if C(x,2) == id2(y)
        D(y,x-z) = C(x,4);
        end
    end  
    z = z + 6042;   
end    

raw2 = D';
stackedplot(raw2);

%%

stackedplot(F_raw');