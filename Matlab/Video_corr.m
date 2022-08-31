clear all; close all;
%script to find missing frames and save them in a csv file 

%% import data
day = cellstr(['2022_07_07']);
datapath =['C:\Users\cha206\Data\DM\g2\' day{1, 1} ]; 
listing = dir(datapath);
isub = [listing(:).isdir];nameFolds = {listing(isub).name}';nameFolds(ismember(nameFolds,{'.','..'})) = [];%makes list with folder names

imaging_datasets=table(nameFolds);
f3=figure;
for s=1:3
%name paths and import
stack=cellstr(table2cell(imaging_datasets(s,1)));
list_path=['C:\Users\cha206\Data\DM\g2\' day{1, 1} '\' stack{1, 1} '\My_V4_Miniscope\timeStamps.csv'];
frame_list=import_framelist(list_path);

%% find dropped frames
frame=frame_list.(1);
time_stamp=frame_list.(2);
figure;subplot(2,1,1),plot(frame(2:end),diff(time_stamp),'ko-');%plot time stamps against frame list to visualize dropped frames
title('before');

fd=median(diff(time_stamp));
fake_time_stamps=time_stamp;
missed_frames=[];

miss=find(diff(fake_time_stamps)>1.5*fd,1);% if time difference is 1.5 times larger than median(fd)
while isnan(miss)==0 %if there are missed frames isnan(miss)= false(0)
    f=miss+1; %select the frame nr after the missed one
    missed_frames=[missed_frames f];
    new_stamp=fake_time_stamps(f-1)+fd;% makes new time stamp for missing frame
    fake_time_stamps=[fake_time_stamps(1:f-1); new_stamp; fake_time_stamps(f:end)];%adds new time stamp in fake time stamp list
    miss=find(diff(fake_time_stamps)>1.5*fd,1);%checks for remaining missing frames in fake_time_stamps list, this list gets corrected and miss=[] if all are corrected
end
subplot(2,1,2),plot(diff(fake_time_stamps),'ko-');
title('after');
missed_frames
filename =['C:\Users\cha206\Data\DM\g2\2022_07_07\' stack{1, 1} '\My_V4_Miniscope\missed_frames.csv'];
writematrix(missed_frames, filename);


end







function frame_list=import_framelist(list_path)
    %% Import csv
    opts = delimitedTextImportOptions("NumVariables", 3);
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";
    opts.VariableNames = ["FrameNumber", "TimeStampms", "BufferIndex"];
    opts.VariableTypes = ["double", "double", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    frame_list = readtable(list_path, opts);
    clear opts
end





