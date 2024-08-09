%% code to check IPI and IDI of new animal in Drivemaze HBITUATION
%02.08.2024
clear all; close all;
%animals = ['F'];% 90 IPIpellet = 18s 30.04.24 IDI 75th perc 5s
%animals = ['E'];% 90 IPIpellet = 16s 12.06.24 IDI 75th perc 3.6s (4s)
%animals = ['B'];%90 IPIpellet = 15s 08.07.24 IDI 75th perc 6s
%Cb1f; pellet = 23s drink = 4 7.8.2024
    events_path=['C:\Users\cha206\Desktop\Opto_DM\Cb1f_events'];%enter name of the animal you want to check
    events=import_DM_opto(events_path);%imports events file from drivemaze
    events_raw = events;
    experiment = table2array(events(:,"experiment"));
    % remove all inh_1 or inh_4 experiments
    I1 = find(experiment == "inh_1");
    I4 = find(experiment == "inh_4");
    excl = vertcat(I1, I4);
    start= datetime(['2024-08-05 09:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');%first proper habituation day (take from notes)
    days = 3;%how many days are the first habituation days, how many days to look at
    f1=figure;
    f2=figure;
    f3=figure;
    f4=figure;
        for d=1:days
            start= datetime(['2024-08-05 09:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
            events= events_raw;
            start = start + day(d-1);%adjust start date
            stop = start + day(1); 

            mid = start + day(1);
            daysplit = mid;
            day1lo= start + hours(12);%lights on
            day2lo= day1lo + day(1);

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
            latency= table2array(events(:,'Latency'));
            retrievals = find(Type == "Pellet_retrieved");
            retrievals = events(retrievals,:);
            ret_lat = retrievals(:,"Latency");
            ret_las = retrievals(:,"Laser");
            ret_dat = retrievals(:,"Date_Time");
        
            %drinks
            D_retrievals = find(Type == "drink");
            D_retrievals = events(D_retrievals,:);
            D_ret_lat = D_retrievals(:,"Drop_latency");
            D_ret_las = D_retrievals(:,"Laser");
            D_ret_dat = D_retrievals(:,"Date_Time");
        
            Drops = table2array(events(:,"Drop_Number"));
            D_latency= table2array(events(:,'Drop_latency'));
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
                subplot(1,days,d);histogram(IPI,edges);hold on;
                    sf = prctile(IPI,75);b=150;%plot([sf sf],[0 b]);
                    n0 = prctile(IPI,90);b=150;;xline(n0, "-", num2str(n0), LineWidth=2);%plot([n0 n0],[0 b]);%plot 90th percentile; 90% of all values are included in that nr
                    n5 = prctile(IPI,95);b=150;%plot([n5 n5],[0 b]);
                    xlabel('IPI (5s bins)');
                    title('Pellets')
                    xlim([0 90]);
                    %ylim([0 100]);  
            end
        
        %% IDIs
        emptD=nanmean(Drops);
        if emptD > 0.3
            IDI =  table2array(D_ret_lat);
            IDI = round(IDI,1);
            big = find(IDI > 90);
            IDI(big) = [];
            edges = [0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8 3 3.2 3.3 3.4 3.6 3.8 4];
            % find meal_size errors       
            figure(f2); 
                subplot(1,days,d);histogram(IDI,edges);hold on;
                n0D = prctile(IDI,90);xline(n0D, "-", num2str(n0D), LineWidth=2);
                title("Drinks");
                xlabel('IDI (0.2s bins)');
                %xlim([0 n0D+1]);
                %ylim([0 50]);
        end 

            %% plot pellets and drinks over time to check if animals eat regularly
            ret_dat = datetime(table2array(ret_dat),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
            figure(f3);
            subplot(days,1,d);histogram(ret_dat, "BinWidth", minutes(10));hold on;
            lim = 25;
            fill([start start+hours(12)+minutes(45) start+hours(12)+minutes(45) start],[0 0 40 40], 'k', FaceAlpha=0.1);
            ylim([0 lim]);
            ylabel('pellets retrieved');
            xlim([start stop]);
        if emptD > 0.3
            D_ret_dat = datetime(table2array(D_ret_dat),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
            figure(f4);
            subplot(days,1,d);histogram(D_ret_dat, "BinWidth", minutes(10));hold on;
            lim = 25;
            fill([start start+hours(12)+minutes(45) start+hours(12)+minutes(45) start],[0 0 40 40], 'k', FaceAlpha=0.1);
            ylim([0 lim]);
            ylabel('drops drunk');
            xlim([start stop]);
        end
        end
    