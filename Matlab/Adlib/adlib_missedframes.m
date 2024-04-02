% code for missed frame correction.
%something is wierd/not working with this. 09.01.24
clear all ;
close all;
load('ci_data_adlib.mat');
%% import data
%g10=load('3354916686_g10_femP150.mat');
%g12=load('3354909574_g12_femP150.mat');
%g15=load('1251667485196_g15_mP200.mat');

% manually transfer mahesh meta files to my format
%for f = 1:6
%    ci_data.adlib.g15.d_2023_09_01(f).session.raw_f_res =  g15.meta(1).imaging(f).cell_f;
%    ci_data.adlib.g15.d_2023_09_01(f).session.Fneu = g15.meta(1).imaging(f).np_f;
%    ci_data.adlib.g15.d_2023_09_01(f).session.events = g15.meta(1).imaging(f).events;
%    ci_data.adlib.g15.d_2023_09_01(f).session.imaging_session = g15.meta(1).imaging(f).stack_id;
%end    
%%
FrameInfo = [];
fn=fieldnames(ci_data.adlib);
%loop through the fields
for ii=1: numel(fn)
    fn1=fieldnames(ci_data.adlib.(fn{ii}));
    animal= fn{ii};
    if animal == 'g10'    
        RFID ='3354916686';
    elseif animal == 'g12'    
        RFID ='3354909574';
    elseif animal == 'g15'    
        RFID ='1251667485196';    
    end
    events_path=['D:\Drivemaze\PFC_LH\' RFID '_events.csv'];
    events=import_events(events_path);

    for j=1: numel(fn1)%day
        fn2=fieldnames(ci_data.adlib.(fn{ii}).(fn1{j}));
        for k=1: numel(fn2)
            for s=1:size(ci_data.adlib.(fn{ii}).(fn1{j}),2)%access the data sessions
                fn3= fieldnames(ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}));
                new_event_frames= ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{3}).(5);% in events column 5 (frame)these are already corrected

                day = fn1{j};
                day = erase(day,"d_");
                time = ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).imaging_session;
                if fn{ii} == 'g15'
                    time = time{1,1};
                end    

                list_path=['D:\Drivemaze\PFC_LH\' fn{ii} '\customEntValHere\' day '\' time '\My_V4_Miniscope\timeStamps.csv'];
                frameList = import_framelist(list_path);
                ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).frame_list = frameList;
                %check for missed frames
                frame=frameList.(1);
                time_stamp=frameList.(2);
                %figure;subplot(2,1,1),plot(frame(2:end),diff(time_stamp),'ko-');
                %title('before');
                %%
%                     figure;subplot(2,1,1),plot(frame(2:end),diff(time_stamp),'ko-');
%                     title('before');
%                     
%                     fd=median(diff(time_stamp));
%                     fake_time_stamps=time_stamp;
%                     missed_frames=[];
%                     
%                     miss=find(diff(fake_time_stamps)>1.4*fd,1);% if time difference is 1.5 times larger than median(fd)
%                     while isnan(miss)==0
%                         f=miss+1;
%                         missed_frames=[missed_frames f];
%                         new_stamp=fake_time_stamps(f-1)+fd;
%                         fake_time_stamps=[fake_time_stamps(1:f-1); new_stamp; fake_time_stamps(f:end)];
%                         miss=find(diff(fake_time_stamps)>1.4*fd,1);
%                     end
%                     subplot(2,1,2),plot(diff(fake_time_stamps),'ko-');
%                     title('after');
                
                fd=median(diff(time_stamp));
                fake_time_stamps=time_stamp;
                missed_frames=[];
                
                miss=find(diff(fake_time_stamps)>1.7*fd,1);% if time difference is 1.5 times larger than median(fd)
                while isnan(miss)==0
                    f=miss+1;
                    missed_frames=[missed_frames f];
                    new_stamp=fake_time_stamps(f-1)+fd;
                    fake_time_stamps=[fake_time_stamps(1:f-1); new_stamp; fake_time_stamps(f:end)];
                    miss=find(diff(fake_time_stamps)>1.7*fd,1);
                end
                ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).missed_frames = missed_frames;

                
                %% take python frame nr out of raw uncorrected events list
                DT_corr = ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).events.(1)(end);
                DT_raw = find(events.Date_Time == DT_corr); 
                event_img_stop = events.frame (DT_raw); % raw nr of timestamps pi aquired
                %correct saved events list in meta back to raw values
                DT_start = ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).events.(1)(1);
                DT_end = ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).events.(1)(end);

                raw_events_list = events(find(events.Date_Time == DT_start):find(events.Date_Time == DT_end),:);
                events_list = ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).events;
                events_list_frames = events_list.("frame");
                for entry =1:numel(events_list.(1))% iterate over altered events list and correct the timestamps back to raw values
                    events_list_frames(entry) = raw_events_list.frame(find(raw_events_list.Date_Time == events_list.Date_Time(entry))); 
                end
                
                ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).events.(5) = events_list_frames; % frame timestamps are now uncorrected
                %% correct raw_f_res and neuropil for dropped frame
                raw_f_res = ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{1});
                new_raw_f = ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{1});
                for i=1:length(missed_frames)
                    rowtoinsert = new_raw_f(missed_frames(i)-1,:);
                    new_raw_f=vertcat(new_raw_f((1:missed_frames(i)-1),:),rowtoinsert,new_raw_f(missed_frames(i):end,:));
                end
                ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).raw_f_res = new_raw_f; %save corrected data
                ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).raw_f_uncorrected = raw_f_res;

                Fneu = ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).Fneu;
                new_Fneu = ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).Fneu;
                for jj=1:length(missed_frames)
                    rowtoinsert2 =  new_Fneu(missed_frames(jj)-1,:);
                    new_Fneu=vertcat(new_Fneu((1:missed_frames(jj)-1),:),rowtoinsert2,new_Fneu(missed_frames(jj):end,:));
                end

                ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).Fneu = new_Fneu; %save corrected data
                ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).Fneu_uncorrected = Fneu;
                % now the missing frames were added to raw_f and not 
                % corrected in events list

                %% collect frame nr to compare
                NrMiss = size(missed_frames,2);% number of missed frames according to time stamps of miniscope
                Nr_frames = size(raw_f_res,1); % number of recorded frames uncorrected
                Nr_frames_corr = size(new_raw_f,1);% number of recorded frames corrected for missed frames
                framelist = frameList(end,1);framelist=framelist{1,1}; % last value from framelist
                
                
                Frames(1,1)= event_img_stop;%value from python events list
                Frames(1,2)= framelist; %from framelist 
                Frames(1,3)= (event_img_stop - framelist); % do frames aquired and python logged match? if positive, some frames dropped? 
                Frames(1,4)= NrMiss;% number of missed frames according to time stamps of miniscope
                Frames(1,5)= Nr_frames_corr-1;%nr frames after mf correction -1 because of matlab value
              
                FrameInfo = [FrameInfo; Frames];
            end
        end
    end
end



