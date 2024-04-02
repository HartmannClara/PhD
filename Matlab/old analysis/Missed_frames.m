clear all; close all;
%% import data
imaging_datasets=table(['17_28_40'; '14_05_20';'14_09_21';'14_40_06';'15_17_30']);%;'14_48_41';'14_56_33']);%;'';'';''
f3=figure;
for s=1
%name paths and import
stack=cellstr(table2cell(imaging_datasets(s,1)));
%stack_path=['C:\Users\cha206\Data\DM\g2\2022_07_07\' stack{1, 1} '\My_V4_Miniscope\Results.csv'];
list_path=['E:\Exp_2_DM_exploration\PFC-LH\g5\2022_09_08\' stack{1, 1} '\My_V4_Miniscope\timeStamps.csv'];
events_path=['E:\Exp_2_DM_exploration\PFC-LH\3354903011_events.csv'];
frame_list=import_framelist(list_path);
events=import_events(events_path);
%raw_f=import_raw(stack_path);

%% refine events
%chop
if s==1
    chop_from=datetime('2022-09-08 17:28:00.000','InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    chop_to=datetime('2022-09-08 17:30:00.000','InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    elseif s==2
    chop_from=datetime('2023-10-27 14:05:00.000','InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    chop_to=datetime('2023-10-27 14:08:00.000','InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    elseif s==3
    chop_from=datetime('2023-10-27 14:09:00.000','InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    chop_to=datetime('2023-10-27 14:39:00.000','InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    elseif s==4
    chop_from=datetime('2023-10-27 14:40:00.000','InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    chop_to=datetime('2023-10-27 15:16:00.000','InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
end

events=sortrows(events,1); %spliced in some events from accidental save in different subject
x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
events(find(isnat(x)==1),:)=[];
x(find(isnat(x)==1))=[];
event_type=events.(4);
x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
events(find(x<chop_from),:)=[];
x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
events(find(x>chop_to),:)=[];
x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
events(find(isnat(x)==1),:)=[];
x(find(isnat(x)==1))=[];
event_type=events.(4);

%% find dropped frames
frame=frame_list.(1);
time_stamp=frame_list.(2);
figure;subplot(2,1,1),plot(frame(2:end),diff(time_stamp),'ko-');
title('before');

fd=median(diff(time_stamp));
fake_time_stamps=time_stamp;
missed_frames=[];

miss=find(diff(fake_time_stamps)>1.5*fd,1);% if time difference is 1.5 times larger than median(fd)
while isnan(miss)==0
    f=miss+1;
    missed_frames=[missed_frames f];
    new_stamp=fake_time_stamps(f-1)+fd;
    fake_time_stamps=[fake_time_stamps(1:f-1); new_stamp; fake_time_stamps(f:end)];
    miss=find(diff(fake_time_stamps)>1.5*fd,1);
end
subplot(2,1,2),plot(diff(fake_time_stamps),'ko-');
title('after');
missed_frames

%% correct event list for dropped frames
event_frames=events.(5);
new_event_frames=events.(5);
for i=1:length(missed_frames)
    to_shift=find(event_frames>missed_frames(i));
    new_event_frames(to_shift)=new_event_frames(to_shift)-1;
end
events.(5)=new_event_frames;


end

function raw_f=import_raw(filename)
    %% Import csv
    opts = delimitedTextImportOptions("NumVariables", 4);
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";
    opts.VariableNames = ["VarName1", "unit_id", "frame", "YrA"];
    opts.VariableTypes = ["double", "double", "double", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    Results = readtable(filename, opts);
    Results = table2array(Results);
    clear opts
    cells=unique(Results(:,2));
    for i=1:length(cells)
        raw_f(:,i)=Results(find(Results(:,2)==cells(i)),4);
    end
end
function events=import_events(filename)
    %% Import csv
    opts = delimitedTextImportOptions("NumVariables", 6);
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";
    opts.VariableNames = ["Date_Time", "amount_consumed", "latency_to_consumption", "Type", "frame", "hardware_time"];
    opts.VariableTypes = ["string", "double", "double", "categorical", "double", "double"];
    opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
    opts = setvaropts(opts, [1, 4], "EmptyFieldRule", "auto");
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    events = readtable(filename, opts);
    clear opts
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
