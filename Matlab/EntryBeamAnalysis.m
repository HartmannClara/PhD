stack_path=['E:\Exp_2_DM_exploration\PFC-LH\EntryTimes\g4_0109_140350'];
cells = 6;
raw_f_res=import_raw(stack_path,cells);%df,ef,ds,es,dd,ed
plot(raw_f_res(:,1:2));




A= raw_f_res(:,1);
B = raw_f_res(:,2);
TF = islocalmin(A,'MinProminence',60);
TE = islocalmin(B,'MinProminence',60);
x=1:size(A,1);
plot(x,A,x,B,x(TF),A(TF),'r*',x(TE),B(TE),'b*')

dd =find(TF == 1);


plot(x,B,x(TE),B(TE),'r*')

