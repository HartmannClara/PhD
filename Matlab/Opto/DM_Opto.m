%% analysis code for drivemaze animals
% first part collects data an makes a struct DM_opto_data 
%second part plots the data separated into animals 
%02.08.2024
%first collect data, separate into hab and inh4 and then into days
clear all; close all;
%animals = ['F'];% 90 IPIpellet = 18s 30.04.24 IDI 75th perc 5s
%animals = ['E'];% 90 IPIpellet = 16s 12.06.24 IDI 75th perc 3.6s (4s) Eb1m
%animals = ['B'];%90 IPIpellet = 15s 08.07.24 IDI 75th perc 6s Bb1f
animals = {'Fb1m','Eb1m', 'Bb1f'};
Group= [A C A];%arch or control
Pintervals= [18 16 15];
Dintervals= [5 4 6];
for i =1:size(animals,2)
    animal= animals{1,i};
    events_path=['C:\Users\cha206\Desktop\Opto_DM\' animal '_events'];%name of the animal you want to check
    events=import_DM_opto(events_path);%imports events file from drivemaze
    events_raw = events;
    experiment = table2array(events(:,"experiment"));
    Time = datetime(table2array(events(:,"Date_Time")),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS'); 
    % find starts of hab and of inh_4
    hab= find(experiment == "hab");hab=hab(1,1);
    inh_4= find(experiment == "inh_4");inh4=inh_4(1,1);

    starthab= Time(hab);
    starthab.Hour=9;starthab.Minute=0;starthab.Second=0;% set the start time of the first day to 9 
    stophab= Time(inh4-1);
    habdays=caldays(between(starthab,stophab,"days"));
    startinh4=Time(inh4);
    startinh4.Hour=9;startinh4.Minute=0;startinh4.Second=0;
    stopinh4= Time(inh_4(end,1)-1);
    inh4days=caldays(between(startinh4,stopinh4,"days"));
   %split habituation days and save into struct, remove days that have no
   %events (weekend in between)
   da=0;
        for d=1:habdays+1
            da=da+1;
            events= events_raw;
            start = starthab + day(d-1);%adjust start date
            stop = start + day(1); 

            chop_from=start - hours(1);%chop the time
            chop_to=stop -hours(1);
            x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
            events(find(x<chop_from),:)=[];%chops events to match the time
            x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
            events(find(x>chop_to),:)=[];
            empt = isempty(events);
            if empt ~= 1;
                DM_Opto_data(i).hab(da).day = events;
                DM_Opto_data(i).ID = animal;
            else
                da=da-1;% reset placement in struct bc day is skipped
            end    
        end 
        %split inh_4 days and save into struct, remove days that have no
   %events (weekend in between)
   da=0;
        for d=1:inh4days+1% +1 to finish the final day
            da=da+1;
            events= events_raw;
            start = startinh4 + day(d-1);%adjust start date
            stop = start + day(1); 

            chop_from=start - hours(1);%chop the time
            chop_to=stop -hours(1);
            x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
            events(find(x<chop_from),:)=[];%chops events to match the time
            x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
            events(find(x>chop_to),:)=[];
            empt = isempty(events);
            if empt ~= 1;
                DM_Opto_data(i).inh4(da).day = events;
            else
                da=da-1;% reset placement in struct bc day is skipped
            end    
        end        
end

%% analysis of DM_Opto_data struct
close all;
clearvars -except DM_Opto_data;
%
f1= figure;
animals = {'Fb1m','Eb1m', 'Bb1f'};
Group= ['A', 'C', 'A'];%arch or control
Pintervals= [18 16 15];%pellet intervals
Dintervals= [5 4 6];%drop intervals

for i =1:size(animals,2)
    animal= animals{1,i};
    Meals = [];
    Drinks = [];
    for ii =1:size(DM_Opto_data(i).inh4,2)
        clear Omls Cmls mls Odr Cdr dri
        events= DM_Opto_data(i).inh4(ii);events=events.day; 
        Type = table2array(events(:,"Type"));
        Time = datetime(table2array(events(:,"Date_Time")),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');

        %plot Opto vs Ctrl meals and drinks for each animal seperately
        Optomeals = find(Type == "Opto meal ended");
        Ctrlmeals = find(Type == "Control meal ended");
        Shortmeals = find(Type == "Short meal ended");
        Otimes = Time(Optomeals);
        Ctimes = Time(Ctrlmeals);
        Stimes = Time(Shortmeals);
        Osizes = table2array(events(Optomeals,"Pellet_number"));
        Csizes = table2array(events(Ctrlmeals,"Pellet_number"));
        Ssizes = table2array(events(Shortmeals,"Pellet_number"));
        %save it in an array, concatenating across days 
        Omls(:,1) = Osizes;  Omls(1:size(Osizes),2) = 1;
        Cmls(:,1) = Csizes;  Cmls(1:size(Csizes),2) = 0;
        mls= vertcat(Omls,Cmls);
        Meals = vertcat(Meals, mls);%both meal types together over all days 

        %drops
        Optodrinks = find(Type == "Opto drink ended");
        Ctrldrinks = find(Type == "Control drink ended");
        Shortdrinks = find(Type == "Short drink ended");
        Otimes = Time(Optodrinks);
        Ctimes = Time(Ctrldrinks);
        Stimes = Time(Shortdrinks);
        Osizes = table2array(events(Optodrinks,"Drop_Number"));
        Csizes = table2array(events(Ctrldrinks,"Drop_Number"));
        Ssizes = table2array(events(Shortdrinks,"Drop_Number"));
        %save it in an array, concatenating across days 
        Odr(:,1) = Osizes;  Odr(1:size(Osizes),2) = 1;
        Cdr(:,1) = Csizes;  Cdr(1:size(Csizes),2) = 0;
        dri= vertcat(Odr,Cdr);
        Drinks = vertcat(Drinks, dri);%both meal types together over all days

    end 
    %plot comparison
    if Group(i) == "A"
        col= "#77AC30";
    else 
        col= "#4DBEEE";
    end

    figure(f1);
    subplot(2,size(animals,2),i);boxchart(Meals(:,2),Meals(:,1),BoxFaceColor=col);hold on;
            swarmchart(Meals(:,2),Meals(:,1),".");
            ylabel('meal size');
            title(animal);
            xticks([0 1]);xticklabels({'control', 'laser on'});
            ylim([0 30]);
            ms= Meals(:,1);
            ctrl_mean = mean(ms(find(Meals(:,2) == 0)));
            las_mean = mean(ms(find(Meals(:,2) == 1)));
            plot([0 1],[ctrl_mean las_mean], 'k*');%xy 
     subplot(2,size(animals,2),i+size(animals,2));boxchart(Drinks(:,2),Drinks(:,1),BoxFaceColor=col);hold on;
            swarmchart(Drinks(:,2),Drinks(:,1),".");
            ylabel('drink size');
            xticks([0 1]);xticklabels({'control', 'laser on'});
            ylim([0 30]);
            ds= Drinks(:,1);
            ctrld_mean = mean(ds(find(Drinks(:,2) == 0)));
            lasd_mean = mean(ds(find(Drinks(:,2) == 1)));
            plot([0 1],[ctrld_mean lasd_mean], 'k*');%xy 
end


    