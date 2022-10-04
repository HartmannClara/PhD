clear all, close all
fa = figure
fb = figure
fc = figure
set(0,'defaultAxesFontSize',20)
%% only bgr
folder = ['13_27_41'; '13_48_56'; '14_08_17'; '14_31_46'; '14_50_06'];

for i = 1 : size(folder) 
   filename=['C:\Users\cha206\Data\Exp_2_DM_exploration\1_Raw\2022_02_08\' folder(i,:) '\My_V4_Miniscope\Results.csv'];
   output_dir = ['C:\Users\cha206\Data\Exp_2_DM_exploration\1_Raw\2022_02_08\' folder(i,:) '\My_V4_Miniscope'];
   raw = readtable(filename);
   raw_arr=table2array(raw);
   raw_arr(:,3:4:95) = raw_arr(:,3:4:95)-raw_arr(:,99);%background subtraction
   raw{:,:}= raw_arr(:,:);%put back into table
   raw_bgr = removevars(raw, [98 99 100 101] );%remove background columns
   save_format = fullfile(output_dir, 'Results_bgr.csv');
   writetable(raw_bgr, save_format);

end


%% import, background subtraction and dff of results files loop 
%always remember to index in the columms and rows. ex: folder(i,:) (with char variables, always index the whole column otherwise u only get the first letter)


folder = ['13_27_41'; '13_48_56'; '14_08_17'];

for i = 1 : size(folder) 
   filename=['C:\Users\cha206\Data\Exp_2_DM_exploration\1_Raw\2022_02_08\' folder(i,:) '\My_V4_Miniscope\Results.csv'];
   output_dir = ['C:\Users\cha206\Data\Exp_2_DM_exploration\1_Raw\2022_02_08\' folder(i,:) '\My_V4_Miniscope']; 
   raw = readtable(filename);
   raw_arr=table2array(raw);
   raw_arr(:,3:4:95) = raw_arr(:,3:4:95)-raw_arr(:,99); %background subtraction
   raw_arr = raw_arr(:,3:4:99);

   medians = [];
    for ROI = 1:size(raw_arr,2)
        cell_n = [];
        cell_n = raw_arr(:,ROI); 
        medians_cell = movmedian(cell_n, 600);%window size in frames %30s
        medians = [medians, medians_cell];     
    end

   df_f_trace = [];
    for cell_median = 1: size(medians,2) % Loop df/f
        v = [];
        d = [];
        v = medians(:,cell_median);%select column from medians
        d = raw_arr(:,cell_median);%select column from activity_bgr
        for f = 1: size(d,1)
            a = [];
            a =  v(f) ;
            b = (d(f) - a)/a;
            df_f_trace(f,cell_median) = b;
        end
    end
   
   raw{:,3:4:95} = df_f_trace(:,1:24);%put back into table
   raw_dff = removevars(raw, [98 99 100 101] );%remove background columns
   save_format = fullfile(output_dir, 'Results_dff.csv');
   writetable(raw_dff, save_format);

end



%% import and background subtraction of results files
clear all;
output_dir = 'C:\Users\cha206\Data\Exp_2_DM_exploration\1_Raw\2022_02_08\13_27_41\My_V4_Miniscope';
raw = readtable('C:\Users\cha206\Data\Exp_2_DM_exploration\1_Raw\2022_02_08\13_27_41\My_V4_Miniscope\Results.csv');
raw_arr=table2array(raw);
raw_arr(:,3:4:95) = raw_arr(:,3:4:95)-raw_arr(:,99);%background subtraction
raw{:,:}= raw_arr(:,:);%put back into table
raw_bgr = removevars(raw, [98 99 100 101] );%remove background columns
save_format = fullfile(output_dir, 'Results_bgr.csv');
writetable(raw_bgr, save_format);


%dff
medians = [];
for ROI = 1:size(raw_arr,2)
    cell_n = [];
    cell_n = raw_arr(:,ROI); 
    medians_cell = movmedian(cell_n, 300);%window size in frames %30s
    medians = [medians, medians_cell];     
end

df_f_trace = [];
for cell_median = 1: size(medians,2) % Loop df/f
    v = [];
    d = [];
    v = medians(:,cell_median);%select column from medians
    d = activity_bgr(:,cell_median);%select column from activity_bgr
    for f = 1: size(d,1)
        a = [];
        a =  v(f) ;
        b = (d(f) - a)/a;
        df_f_trace(f,cell_median) = b;
    end
end


figure(fc);
stackedplot(raw_bg)





