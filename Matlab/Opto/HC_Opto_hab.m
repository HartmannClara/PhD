%% code to check homecage setup habituation data to get pellet intervals for inh experiments
clear all; close all;
animals = ['A' 'B' 'C' 'D' 'E' 'F'];
%pick the intervals best for the experiment once code is run
intvl = 0;%set to 0 if only looking for intervals not checking
intervals = [15, 23, 15, 21, 12, 18];
%% 
f1 = figure;
f2 = figure;
f3 = figure;
allMt = table;
perc = [];mz_perc = [];
Opto_data_all =[];
k=0;
for i=1:size(animals,2)
    animal = animals(i);
    events_path=['C:\Users\cha206\Desktop\OptoData\animal_' animal '3_events'];% change batch nr here
    events=import_opto(events_path);
    events_raw= events;
    experiment = table2array(events(:,"experiment"));
    % remove all inh_1 or inh_2 experiments
    I1 = find(experiment == "inh_1");
    I4 = find(experiment == "inh_4");
    excl = vertcat(I1, I4);
    events(excl, :) = [];

    %start= datetime(['2024-06-10 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');%start day of habituation, HC experiments start at 10! except B2 bc of time change
    days = 4;%how many days of hab to check
    for d=1:days
            start= datetime(['2024-06-10 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
            events= events_raw;
            start = start + day(d-1);%adjust start date
            stop = start + day(1); 

            chop_from=start - hours(1);%chop the time
            chop_to=stop -hours(1);
            x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
            events(find(x<chop_from),:)=[];%chops events to match the time
            x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
            events(find(x>chop_to),:)=[];
            
            Type = table2array(events(:,"Type"));
            Time = datetime(table2array(events(:,"Date_Time")),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
            
            %meals
            Pellets = table2array(events(:,"Pellet_number"));
            Latency= table2array(events(:,'Latency'));
            retrievals = find(Type == "Pellet_Retrieved");
            retrievals = events(retrievals,:);
            ret_lat = retrievals(:,"Latency");
            ret_las = retrievals(:,"Laser");
            ret_dat = retrievals(:,"Date_Time");

            if intvl == 1
                    sz = 1; 
                    meal_size = table;
                    for j = 1:size(ret_lat,1)
                        lat = table2array(ret_lat(j,:));
                        if lat <= intervals(i) % smaller or equal to selected IPI in seconds
                            sz = sz + 1; % size of the meal
                        elseif j == 1
                            meal = [ret_dat(j,:), array2table(sz) , ret_las(j,:) ];
                            meal_size = [meal_size; meal];
                            sz = 1;
                        else
                            meal = [ret_dat(j-1,:), array2table(sz) , ret_las(j-1,:) ];
                            meal_size = [meal_size; meal];
                            sz = 1;
                        end
                    end

    % find meal_size errors
    sizes= table2array(meal_size(:,2));
    err= find(sizes > 50);
    sizes(err,:) = [];
    meal_size(err,:) = [];
    edges2 = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
    figure(f3);
        subplot(6,days,d+k); histogram(sizes, edges2);
        xline(3,LineWidth=2);
        xlabel('IPI (5s bins)');
        ylabel(animal);
        title(string(start));
        ylim([0 90]);
        


            else  
            %%  Inter pellet intervals
                IPI =  table2array(ret_lat);
                er= find(IPI < 1.5);% look for errors, like a broken fed
                IPI(er) = [];
                IPI = round(IPI);
                big = find(IPI > 90);%exclude anything bigger than 90 seconds,
                IPI(big) = [];
                edges = [0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90];%bins to plot IPI histogram
                emptP = nanmean(Pellets,1);
                if emptP > 0.3
                    figure(f1);
                    subplot(6,days,d+k);histogram(IPI,edges);hold on;
                    sf = prctile(IPI,75);b=150;%plot([sf sf],[0 b]);
                    n0 = prctile(IPI,90);b=150;;xline(n0, "-", num2str(n0), LineWidth=2);%plot 90th percentile; 90% of all values are included in that nr
                    n5 = prctile(IPI,95);b=150;%plot([n5 n5],[0 b]);
                    xlabel('IPI (5s bins)');
                    ylabel(animal);
                    title(string(start));
                    xlim([0 90]);
                    ylim([0 100]);
    
                end
            end
    end

    k=k+days;
end
%%

