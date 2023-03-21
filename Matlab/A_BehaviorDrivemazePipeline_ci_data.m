% Drivemaze pipeline without frames! just behavior
% only analysis, creates missed frame corrected data in ci_data struct
% sorted by animal, day, sessiohn
%08.03.2023

%g5 27,08 aug
%g2 25 jul
%g4 01 02 sept


clearvars -except ci_data;
%% import data
animal = 5;%;cells =23;%CHANGE nr of cells in FOV!!!!!!!!!!
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
elseif animal == 10    
    RFID ='3354916686';
elseif animal == 12    
    RFID ='3354909574';    
end    
%events_path=['D:\Exp_2_DM_exploration\PFC-LH\' RFID '_events.csv'];
events_path=['C:\Users\cha206\Data\DM\animals2023\' RFID '_events.csv'];
%topLevelFolder = uigetdir();
% Get a list of all files and folders in this folder.
%folders = getFolders(topLevelFolder);
%ko=0;bad=0; %bad imaging dataset to be excluded later
%imaging_datasets = folders(ko+1:(size(folders,1)),:);%imaging_datasets=table(['11_25_35'; '11_34_12']);
%to add to linear imaging starts index, if imaging files start at 20, ko=19

%for s=1:size(imaging_datasets,1)%:5%:3%change which imaging sets here
%name paths and import
%stack=cellstr(table2cell(imaging_datasets(s,1)));
%stack_path=[topLevelFolder '\' stack{1, 1} '\My_V4_Miniscope\Results.csv'];
%list_path=[topLevelFolder '\' stack{1, 1} '\My_V4_Miniscope\timeStamps.csv'];
%frame_list=import_framelist(list_path);
events=import_events(events_path);
%raw_f_res=import_raw(stack_path,cells); %raw data from results file: (filename,cells)
%raw_f= zscore(df_f_trace,0,1);% zscored using sample sd

%% chop
level = wildcardPattern + "\";pat = asManyOfPattern(level);
%day = extractAfter(topLevelFolder,pat);day2 = strrep(day,'_','-');% get the day
day2 = '2023-03-06';day = strrep(day2,'-','_');maze = char("Swap");
chop_from=datetime([day2 ' 12:22:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');%chop the day
chop_to=datetime([day2 ' 14:50:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
events(find(x<chop_from),:)=[];%chops events to match the day
x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
events(find(x>chop_to),:)=[];
x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS'); %make the updated x axis
events(find(isnat(x)==1),:)=[];
x(find(isnat(x)==1))=[];
event_type=events.(4);


%% cleanup
% find imaging stops
u1=find(event_type=='imaging_stop');
%frames =events.(5);
%abc = u1 - 1; 
imaging_ends=u1;%(find(frames(u1) > 5 ));%frames were aquired
%find imaging starts
u2=find(event_type=='imaging_start');

imaging_starts=[];%accumulator
for i=1:size(imaging_ends,1)
    u3=u2(find(u2 < imaging_ends(i)));%finds the previous imaging_starts
    imaging_starts=[imaging_starts; u3(end)];%adds the closest imaging start to list
end
% make fake imaging datasets list
for i = 1:numel(imaging_starts)
    id = table2array(events(imaging_starts(i),1));id= extractBetween(id,' ','.');id=strrep(id,':','_');
    imaging_datasets(i,1) = id;
end    
x2 =x;
events2=events;

for s=1:size(imaging_datasets,1)    
    
%% refine events
x=x2;
events=events2;
%chop sessions %maybe do this automatic?
%chop_from = x(imaging_starts(s+ko));
%chop_to=x(imaging_ends(s+ko));
chop_from = x(imaging_starts(s));
chop_to=x(imaging_ends(s));

x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
events(find(x<chop_from),:)=[];
x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
events(find(x>chop_to),:)=[];
x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
events(find(isnat(x)==1),:)=[];
x(find(isnat(x)==1))=[];
event_type=events.(4);

%saving to struct
date = char(['d_', day]);
animalNr = char(['g', num2str(animal)]);
%ci_data.bsl.(animalNr).(date)(s).session.raw_f_res = raw_f_res;
ci_data.(maze).(animalNr).(date)(s).session.events = events;
%ci_data.bsl.(animalNr).(date)(s).session.missed_frames = missed_frames;
%ci_data.bsl.(animalNr).(date)(s).session.imaging_session = imaging_datasets(s,1);
end
ci_data_beh = ci_data;
save("ci_data_beh.mat", "ci_data_beh");

%%
ci_data_beh.bsl.g5 = ci_data.bsl.g5;
ci_data_beh.bsl.g2 = ci_data.bsl.g2;
ci_data_beh.bsl.g4 = ci_data.bsl.g4;


