%% analzying missing frames and TTL pulse jitter
clear all ; 
close all;

%%
FrameInfo = [];
fn=fieldnames(ci_data.bsl);
%loop through the fields
for ii=1: numel(fn)
    fn1=fieldnames(ci_data.bsl.(fn{ii}));
    for j=1: numel(fn1)%day
        fn2=fieldnames(ci_data.bsl.(fn{ii}).(fn1{j}));
        for k=1: numel(fn2)
            for s=1:size(ci_data.bsl.(fn{ii}).(fn1{j}),2)%access the data sessions
                fn3= fieldnames(ci_data.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}));

                missed_frames = ci_data.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{4}); 
                NrMiss = size(missed_frames,2);% number of missed frames according to time stamps of miniscope

                %new_event_frames= ci_data.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{3}).(5);% in events column 5 (frame)
                day = fn1{j};
                day = erase(day,"d_");
                time = ci_data.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}).imaging_session;time = time{1,1}; time= time{1,1};
                list_path=['E:\Exp_2_DM_exploration\PFC-LH\' fn{ii} '\' day '\' time '\My_V4_Miniscope\timeStamps.csv'];
                frameList = import_framelist(list_path);
                img_stop_frames = frameList(end,1);img_stop_frames=img_stop_frames{1,1};

                raw_f_res=ci_data.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{1});
                Nr_frames = size(raw_f_res,1);%miniscope aquired

                Frames(1,1)=Nr_frames;%miniscope aquired
                Frames(1,2)=img_stop_frames+1; %from framelist +1 because python and matlab 0 difference
                Frames(1,3)= NrMiss;
                Frames(1,4)= (Nr_frames - (img_stop_frames+1))  + NrMiss; 
                FrameInfo = [FrameInfo; Frames];
            end
        end
    end
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