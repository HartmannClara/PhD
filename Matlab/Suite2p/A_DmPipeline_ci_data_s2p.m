% Drivemaze pipeline 
% only analysis, creates missed frame corrected data in ci_data struct
% sorted by animal, day, session
%30.11.2022

%g5 27,08 aug
%g2 25 jul
%g4 01 02 sept
%g5, g4, g2 ,g12

clearvars -except ci_data; 
%% import data
animal = 4;
if animal == 2 
    RFID = '328340226232';
elseif animal == 4
    RFID = '335490167178';
elseif animal == 5    
    RFID ='3354903011';
elseif animal == 6    
    RFID ='8501114126';
elseif animal == 7    
    RFID ='8502199200';
elseif animal == 9    
    RFID ='8501181185';
elseif animal == 10    
    RFID ='3354916686';
elseif animal == 11    
    RFID ='8500188177';
elseif animal == 12    
    RFID ='3354909574';
elseif animal == 14    
    RFID ='8500142131';    
end    
events_path=['E:\Exp_2_DM_exploration\PFC-LH\' RFID '_events.csv'];
FrameDiff = [];
topLevelFolder = uigetdir();
% Get a list of all files and folders in this folder.
folders = getFolders(topLevelFolder);
ko=0;bad=[0]; %bad imaging dataset to be excluded later
imaging_datasets = folders(ko+1:(size(folders,1)),:);%imaging_datasets
if bad ~= 0 % removes bad datasets
    imaging_datasets(bad,:) = []; 
end    
%to add to linear imaging starts index, if imaging files start at 20, ko=19
for s=1:(size(imaging_datasets,1)-1)%:5%:3%change which imaging sets here !!!!!!!!!!!!!! remove faulty ones from list
%name paths and import
stack=cellstr(table2cell(imaging_datasets(s,1)));
stack_path=[topLevelFolder '\' stack{1, 1} '\My_V4_Miniscope\F.csv'];
stack2_path=[topLevelFolder '\' stack{1, 1} '\My_V4_Miniscope\Fneu.csv'];
list_path=[topLevelFolder '\' stack{1, 1} '\My_V4_Miniscope\timeStamps.csv'];
frame_list=import_framelist(list_path);


raw_f_res=import_raw_s2p(stack_path);%raw data from F.csv file: (filename)
Fneu = import_raw_s2p(stack2_path);
%raw_f= zscore(df_f_trace,0,1);% zscored using sample sd

events=import_events(events_path);
Og_events = events; 
%% chop
level = wildcardPattern + "\";pat = asManyOfPattern(level);
day = extractAfter(topLevelFolder,pat);day2 = strrep(day,'_','-');% get the day
chop_from=datetime([day2 ' 09:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');%chop the day
chop_to=datetime([day2 ' 22:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
events(find(x<chop_from),:)=[];%chops events to match the day
x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
events(find(x>chop_to),:)=[];

% make sure all blocks who have no imaging stop only block_end get one
type= table2array(events(:,4));
onlyBE = find(type == 'block_end');
for t=1:size(onlyBE,1)
    if type(onlyBE(t)+1) ~= "imaging_stop"
       newIS = events(onlyBE(t),:);newIS(1,"Type") = {"imaging_stop"};
       events= [events(1:onlyBE(t),:);newIS; events(onlyBE(t)+1:end,:)]; 
    end
    type= table2array(events(:,4));
    onlyBE = find(type == 'block_end');
end 

x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS'); %make the updated x axis
events(find(isnat(x)==1),:)=[];
x(find(isnat(x)==1))=[];
event_type=events.(4);

%% cleanup
% find imaging stops
u1=find(event_type=='imaging_stop');
frames =events.(5);
%abc = u1 - 1; 
imaging_ends=u1(find(frames(u1) > 5 ));%frames were aquired
%find imaging starts
u2=find(event_type=='imaging_start');
imaging_starts=[];%accumulator
for i=1:size(imaging_ends,1)
    u3=u2(find(u2 < imaging_ends(i)));%finds the previous imaging_starts
    imaging_starts=[imaging_starts; u3(end)];%adds the closest imaging start to list
end
%remove bad imaging datasets from imaging starts and stops 
if bad ~= 0
    imaging_starts(bad,:) = [];
    imaging_ends(bad,:) = [];
else
    imaging_starts = imaging_starts;
end

%% refine events
%chop sessions %maybe do this automatic?
chop_from = x(imaging_starts(s+ko));
chop_to=x(imaging_ends(s+ko));

x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
events(find(x<chop_from),:)=[];
x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
events(find(x>chop_to),:)=[];
x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
events(find(isnat(x)==1),:)=[];
x(find(isnat(x)==1))=[];
event_type=events.(4);

%% check if frames aquired matches event list
%loop through sessions and compare Fneu to events list frames last entry 
if size(raw_f_res,1) ~= size(frame_list,1)% miniscope frames aquired always one more than in timeStamps list because frame 0!!
   nomatch = ["no match in session", s, day, animal ]; 
   disp(nomatch)
end  
%% find dropped frames
frame=frame_list.(1);
time_stamp=frame_list.(2);

fd=median(diff(time_stamp));% find median difference between frames (10 Hz)
fake_time_stamps=time_stamp;
missed_frames=[];

miss=find(diff(fake_time_stamps)>1.5*fd,1);%checks if there are missing frames by comparing diff between timestamps to median
while isnan(miss)==0
    f=miss+1;%makes a new frame nr
    missed_frames=[missed_frames f];
    new_stamp=fake_time_stamps(f-1)+fd;%makes a new time stamp
    fake_time_stamps=[fake_time_stamps(1:f-1); new_stamp; fake_time_stamps(f:end)];
    miss=find(diff(fake_time_stamps)>1.5*fd,1);
end
%missed_frames;

%% correct event list for dropped frames
event_frames=events.(5);
new_event_frames=events.(5);
for i=1:length(missed_frames)
    to_shift=find(event_frames>missed_frames(i));
    new_event_frames(to_shift)=new_event_frames(to_shift)-1;
end
events.(5)=new_event_frames;  

%% check correct sorting
Fdiff = size(raw_f_res,1)- events.frame(end);
FrameDiff = [FrameDiff;Fdiff];

%% saving to struct
date = char(['d_', day]);
animalNr = char(['g', num2str(animal)]);
ci_data.bsl.(animalNr).(date)(s).session.raw_f_res = raw_f_res;
ci_data.bsl.(animalNr).(date)(s).session.Fneu = Fneu;
ci_data.bsl.(animalNr).(date)(s).session.events = events;
ci_data.bsl.(animalNr).(date)(s).session.missed_frames = missed_frames;
ci_data.bsl.(animalNr).(date)(s).session.imaging_session = imaging_datasets(s,1);
end

FrameDiff;
save("ci_data_s2p.mat", "ci_data");