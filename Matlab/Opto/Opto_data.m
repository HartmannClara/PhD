% collection of OPTO data files
%% loop for all animals
clear all; close all;
animals = ['A','B', 'C', 'D', 'E', 'F'];
%intervals = [29, 43, 22, 23, 20, 17];%batch1
intervals = [41, 25, 31, 18, 14, 18];%batch2
%animals = ['A', 'C', 'D', 'E', 'F'];

allMt = table;
perc = [];mz_perc = [];
Opto_data_all =[];
for i=1:size(animals,2)
    animal = animals(i);
    events_path=['C:\Users\cha206\Desktop\OptoData\animal_' animal '_events'];
    events=import_opto(events_path);
    %correct events list for errors?




    % cut events into hab and the different experiments
    %batch1
    starthab= datetime(['2024-02-29 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    startinh_4=datetime(['2024-03-04 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    start_RF=datetime(['2024-03-12 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    endRF = datetime(['2024-03-13 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    startinh_1=datetime(['2024-03-14 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    %batch2
    % first adjust in events file for time change at 2:00 am on 31.03.24 time goes one
    % hour forward
%     starthab= datetime(['2024-03-21 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
%     start_RF=datetime(['2024-04-03 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
%     endRF = datetime(['2024-04-04 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
%     startinh_1=datetime(['2024-04-05 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
% 
%     if animal == "E" || animal == "D"
%         startinh_4=datetime(['2024-03-27 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
%     else
%         startinh_4=datetime(['2024-03-25 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
%     end

    % chop event data for experiments
    events_hab = events;
    x=datetime(events_hab.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS'); events_hab(find(x<starthab),:)=[];
    x=datetime(events_hab.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');events_hab(find(x>startinh_4),:)=[];
    events_inh4 = events;
    x=datetime(events_inh4.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS'); events_inh4(find(x<startinh_4),:)=[];
    x=datetime(events_inh4.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');events_inh4(find(x>start_RF),:)=[];
    events_RF = events;
    x=datetime(events_RF.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS'); events_RF(find(x<start_RF),:)=[];
    x=datetime(events_RF.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');events_RF(find(x>endRF),:)=[];
    events_inh1 = events;
    x=datetime(events_inh1.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS'); events_inh1(find(x<startinh_1),:)=[];
    

    experiment = table2array(events(:,"experiment"));
    % remove all inh_1 or hab experiments
    I1 = find(experiment == "inh_1");
    %I1 = find(experiment == "inh_4");
    HAB = find(experiment == "hab");
    excl = vertcat(I1, HAB);
    events(excl, :) = [];


    %start= datetime(['2024-03-04 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');%start day
    start= datetime(['2024-03-27 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');%start day
    days = 2;%how many days to plot
    hab = 0; %change if not habitutation

    
    
    %chop_from=datetime(['2024-03-01 09:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');%chop the day
    %chop_to=datetime(['2024-03-03 09:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    chop_from=start - hours(1);%chop the day
    chop_to=stop -hours(1);

    x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    events(find(x<chop_from),:)=[];%chops events to match the day
    x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    events(find(x>chop_to),:)=[];
    events_raw = events;
    
    Type = table2array(events(:,"Type"));
    Pellets = table2array(events(:,"Pellet_number"));
    Time = datetime(table2array(events(:,"Date_Time")),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    latency= table2array(events(:,'Latency'));
    
   

end 