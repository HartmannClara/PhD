% collection of OPTO data files into Optofiles struct
%errors already removed and timechange accounted for
%b1 and b2 are red in separetly and saved
% some manual steps to adjust for the different start days
%% loop for all animals
close all;
%clear all;
animals = ['A','B', 'C', 'D', 'E', 'F'];
%Labels = ['C','A', 'A', 'A', 'C', 'A'];
%intervals = [29, 43, 22, 23, 20, 17];%batch1
%intervals = [41, 25, 31, 18, 14, 18];%batch2
intervals = [15, 23, 15, 21, 12, 18];%batch3
%Labels = ['A','A', 'A', 'C', 'C', 'A'];%b2
Labels = ['A','C', 'A', 'C', 'C', 'C'];%b3
%b=0;
%b=6;
b=12;
allMt = table;
perc = [];mz_perc = [];
Opto_data_all =[];
for i=1:size(animals,2)
    animal = animals(i);
    Label = Labels(i);
    events_path=['C:\Users\cha206\Desktop\OptoData\animal_' animal '3_events'];
    events=import_opto(events_path);
    %correct events list for errors?
    Type = table2array(events(:,"Type"));
    FEDerror = find(Type == "FED error");
    events(FEDerror,:) = [];Type = table2array(events(:,"Type"));
    retrievals = find(Type == "Pellet_Retrieved");
    retrievals = events(retrievals,:);
    ret_lat = retrievals(:,"Latency");errors = find(table2array(ret_lat) < 1.5);
    events(errors,:) = [];


    % cut events into hab and the different experiments
    %batch1
% %     starthab= datetime(['2024-02-29 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
% %     startinh_4=datetime(['2024-03-04 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
% %     start_RF=datetime(['2024-03-12 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
% %     endRF = datetime(['2024-03-13 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
% %     startinh_1=datetime(['2024-03-14 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');

    %batch2
    % first adjust in events file for time change at 2:00 am on 31.03.24 time goes one
    % hour forward
% %     timechange= datetime(['2024-03-31 02:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
% %     x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
% %     for j =1:size(x,1)
% %         if x(j) > timechange
% %             x(j) = x(j)-hours(1);
% %         end    
% %     end
% %     events.(1) =   x;
% %     starthab= datetime(['2024-03-21 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
% %     start_RF=datetime(['2024-04-03 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
% %     endRF = datetime(['2024-04-04 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
% %     startinh_1=datetime(['2024-04-06 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
% % 
% %     if animal == "E" 
% %         startinh_4=datetime(['2024-03-27 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
% %     elseif animal == "D"
% %         startinh_4=datetime(['2024-03-27 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
% %     else
% %         startinh_4=datetime(['2024-03-25 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
% %     end

%batch3
    starthab= datetime(['2024-06-11 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    startinh_4=datetime(['2024-06-14 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    start_RF=datetime(['2024-06-22 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    endRF = datetime(['2024-06-22 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    startinh_1=datetime(['2024-06-22 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');




    
    % chop event data for experiments
    events_hab = events;
    x=datetime(events_hab.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS'); events_hab(find(x<starthab),:)=[];
    x=datetime(events_hab.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');events_hab(find(x>startinh_4),:)=[];
    events_inh4 = events;
    x=datetime(events_inh4.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS'); events_inh4(find(x<startinh_4),:)=[];
    x=datetime(events_inh4.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');events_inh4(find(x>start_RF),:)=[];
    events_inh4_2 = events;
    x=datetime(events_inh4_2.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS'); events_inh4_2(find(x<endRF),:)=[];
    x=datetime(events_inh4_2.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');events_inh4_2(find(x>startinh_1),:)=[];
    events_inh4 = vertcat(events_inh4,events_inh4_2);
    events_RF = events;
    x=datetime(events_RF.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS'); events_RF(find(x<start_RF),:)=[];
    x=datetime(events_RF.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');events_RF(find(x>endRF),:)=[];
    events_inh1 = events;
    x=datetime(events_inh1.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS'); events_inh1(find(x<startinh_1),:)=[];
    
    %save data in struct
    Optodata(i+b).events = events;
    Optodata(i+b).hab = events_hab;
    Optodata(i+b).inh4 = events_inh4;
    Optodata(i+b).RF = events_RF;
    Optodata(i+b).inh1 = events_inh1;
    Optodata(i+b).Label = Label;
    Optodata(i+b).setup = animal;
end 


