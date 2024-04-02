% analysis of OPTo files

%% loop for all animals
clear all; close all;
animals = ['A','B', 'C', 'D', 'E', 'F'];
%intervals = [29, 43, 22, 23, 20, 17];%batch1
intervals = [41, 25, 31, 18, 14, 18];%batch2
%animals = ['A', 'C', 'D', 'E', 'F'];
f0 = figure;
f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;
f5 = figure;
f6 = figure;
allMt = table;
perc = [];mz_perc = [];
Opto_data_all =[];
for i=1:size(animals,2)
    animal = animals(i);
    events_path=['C:\Users\cha206\Desktop\OptoData\animal_' animal '2_events'];
    events=import_opto(events_path);
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

    mid = start + day(1);
    stop = start + day (days);
    daysplit = mid;
    day1lo= start + hours(12);%lights on
    day2lo= day1lo + day(1);
    
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
    %create meals independent of code
    
    retrievals = find(Type == "Pellet_Retrieved");
    retrievals = events(retrievals,:);
    ret_lat = retrievals(:,"Latency");errors = find(table2array(ret_lat) < 1.5);
    events(errors,:) = [];
    Type(errors,:) = [];
    Pellets(errors,:) = [];
    Time(errors,:) = [];
    latency(errors,:) = [];
    retrievals(errors,:) = [];
    ret_lat = retrievals(:,"Latency");
    ret_las = retrievals(:,"Laser");
    ret_dat = retrievals(:,"Date_Time");
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
    
    %Inter pellet intervals
    IPI =  table2array(ret_lat);
    %er= find(IPI < 2);% look for errors???
    %IPI(er) = [];
    IPI = round(IPI);

    big = find(IPI > 90);
    IPI(big) = [];
    edges = [0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90];

%     figure(f1); 
%     subplot(2,6,i); histogram(IPI,edges);hold on;
%         sf = prctile(IPI,75);b=150;%plot([sf sf],[0 b]);
%         n0 = prctile(IPI,90);b=150;%plot([n0 n0],[0 b]);
%         n5 = prctile(IPI,95);b=150;%plot([n5 n5],[0 b]);
%         %collect percentiles
%         pt = [i,sf, n0, n5];
%         perc = vertcat(perc, pt);
%         title(animal);
%         xlabel('IPI (5s bins)');
%         xlim([0 90]);
%         ylim([0 100]);

    %plot latency to 4th, 5th or 6th pellet
    if animal == "D"| animal == "E"
         col = "#4DBEEE";
    else
        col = "#77AC30";
    end 

    ret_latency = table2array(retrievals(:,"Latency"));
    ret_pel = table2array(retrievals(:,"Pellet_number"));
    ret_laser = table2array(retrievals(:,"Laser"));
    lat_4 = ret_latency((find(ret_pel == 4)));short = find(lat_4 < 2);lat_4(short) = [];
    las_4 = ret_laser((find(ret_pel == 4)));las_4(short) = [];
    lat_5 = ret_latency((find(ret_pel == 5)));short = find(lat_5 < 2);lat_5(short) = [];
    las_5 = ret_laser((find(ret_pel == 5)));las_5(short) = [];
    lat_6 = ret_latency((find(ret_pel == 6)));short = find(lat_6 < 2);lat_6(short) = [];
    las_6 = ret_laser((find(ret_pel == 6)));las_6(short) = [];
    
    figure(f5);
    subplot(3,6,i);boxchart(las_4,lat_4,BoxFaceColor=col);hold on;
        swarmchart(las_4,lat_4,".");
        title(animal);
        ylabel('latency to 4th pellet');
        xticks([0 1]);xticklabels({'control', 'laser on'});
        ylim([0 30]);
        ctrl_mean = mean(lat_4(find(las_4 == 0)));
        las_mean = mean(lat_4(find(las_4 == 1)));
        plot([0 1],[ctrl_mean las_mean], 'k*');%xy

    subplot(3,6,i+6);boxchart(las_5,lat_5,BoxFaceColor=col);hold on;
        swarmchart(las_5,lat_5,".");
        title(animal);
        ylabel('latency to 5th pellet');
        xticks([0 1]);xticklabels({'control', 'laser on'});
        ylim([0 30]);
        ctrl_mean = mean(lat_5(find(las_5 == 0)));
        las_mean = mean(lat_5(find(las_5 == 1)));
        plot([0 1],[ctrl_mean las_mean], 'k*');%xy

    subplot(3,6,i+12);boxchart(las_6,lat_6,BoxFaceColor=col);hold on;
        swarmchart(las_6,lat_6,".");
        title(animal);
        ylabel('latency to 6th pellet');
        xticks([0 1]);xticklabels({'control', 'laser on'});
        ylim([0 30]);
        ctrl_mean = mean(lat_6(find(las_6 == 0)));
        las_mean = mean(lat_6(find(las_6 == 1)));
        plot([0 1],[ctrl_mean las_mean], 'k*');%xy
    
    
    % find meal_size errors
    %sizes= table2array(meal_size(:,2));
    %err= find(sizes > 30);
    %sizes(err,:) = [];
    %meal_size(err,:) = [];

    % classify satiety level according ot latency to last meal

    %cut meals and log single pellet trials if trial type consistent:

    Opto_meal_ends = find(Type == "Opto meal ended");
    Ctrl_meal_ends = find(Type == "Control meal ended");
    Short_meal_ends = find(Type == "Short meal ended");
    Shortmeals = Short_meal_ends;
    Optomeals = Opto_meal_ends;
    Ctrlmeals = Ctrl_meal_ends;
    
    for e = 1:size(Opto_meal_ends,1)
        en = Opto_meal_ends(e);
        ms = Pellets(en);
        Opto_meal_ends(e,3) = ms;
        Opto_meal_ends(e,2) = latency(en-ms);  
    end    
    for e = 1:size(Ctrl_meal_ends,1)
        en = Ctrl_meal_ends(e);
        ms = Pellets(en);
        Ctrl_meal_ends(e,3) = ms;
        Ctrl_meal_ends(e,2) = latency(en-ms);  
    end   

%plot meal size against time since last meal
%     figure(f0);
%     subplot(2,6,i); scatter(Ctrl_meal_ends(:,2),Ctrl_meal_ends(:,3)); 
%         title(animal); xlabel({'latency since last meal'});
%         ylabel({'ctrl meal size'});
%         xlim([0 7200]);ylim([2 10]);
%     subplot(2,6,i+6); scatter(Opto_meal_ends(:,2),Opto_meal_ends(:,3)); 
%         xlabel({'latency since last meal'});
%         ylabel({'opto meal size'});
%         xlim([0 7200]);ylim([2 10]);
%satiety index: nr of pellets consumed in the last 20 min from start of
%meal
    
    for m = 1:size(Opto_meal_ends,1)
        om = Opto_meal_ends(m,1);
        omtime = Time(om-(Opto_meal_ends(m,3)+1));
        last30 = omtime - minutes(20);%%%%%%%%%%%%%%%%%%%%%%% time
        last30time = events;
        x=datetime(last30time.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
        last30time(find(x>omtime),:)=[];%chops events to match the day
        x=datetime(last30time.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
        last30time(find(x<last30),:)=[];
        pell30 = size(find(table2array(last30time(:,2)) == "Pellet_Retrieved"),1);
        Opto_meal_ends(m,4) = pell30; %save satiety status
        if pell30 < 6
            Opto_meal_ends(m,5) = 4;%hungry
        elseif 5 < pell30 & pell30 < 8
            Opto_meal_ends(m,5) = 2;
        else 
            Opto_meal_ends(m,5) = 2;%sated
        end        
    end

    for m = 1:size(Ctrl_meal_ends,1)
        om = Ctrl_meal_ends(m,1);
        omtime = Time(om-Ctrl_meal_ends(m,3));;
        last30 = omtime - minutes(20);
        last30time = events;
        x=datetime(last30time.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
        last30time(find(x>omtime),:)=[];%chops events to match the day
        x=datetime(last30time.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
        last30time(find(x<last30),:)=[];
        pell30 = size(find(table2array(last30time(:,2)) == "Pellet_Retrieved"),1);
        Ctrl_meal_ends(m,4) = pell30; %save satiety status
        if pell30 < 6
            Ctrl_meal_ends(m,5) = 3;%hungry
        elseif 5 < pell30 & pell30 < 8
            Ctrl_meal_ends(m,5) = 1;
        else 
            Ctrl_meal_ends(m,5) = 1; %sated
        end    
            
    end
    % plot histogram of satiety measures
    figure(f0);
        subplot(2,6,i);histogram(Ctrl_meal_ends(:,4),[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26]);hold on;
        title(animal);
        xlabel('Nr pellets 20 min prior (Ctrl)')
        ylim([0 40]);
        subplot(2,6,i+6);histogram(Opto_meal_ends(:,4),[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26]);hold on;
        xlabel('Nr pellets 20 min prior (ArchT)');
        ylim([0 40]);

%     figure(f0);
%     subplot(2,6,i); boxchart(Ctrl_meal_ends(:,5),Ctrl_meal_ends(:,3));hold on;
%     swarmchart(Ctrl_meal_ends(:,5),Ctrl_meal_ends(:,3),".");
%     boxchart(Opto_meal_ends(:,5),Opto_meal_ends(:,3),BoxFaceColor=col);
%     swarmchart(Opto_meal_ends(:,5),Opto_meal_ends(:,3),".");
%     title(animal);
%     ylim([0 25]);
%     xlim([0 5]);
%     ylabel("meal size")
%     xticks([1 2 3 4]);
%     xticklabels({'sated ctrl','sated lon', 'hungry ctrl', 'hungry lon'});
%     edges = [0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 60];
%     subplot(2,6,i+6);histogram(Opto_meal_ends(:,4),edges);

%     figure(f1);
%     subplot(2,6,i+ 6); histogram(sizes);hold on;
%         mz_sf = prctile(sizes,75);b=100;%plot([mz_sf mz_sf],[0 b]);
%         mz_n0 = prctile(sizes,90);b=100;%plot([mz_n0 mz_n0],[0 b]);
%         mz_n5 = prctile(sizes,95);b=100;%plot([mz_n5 mz_n5],[0 b]);
%         %collect percentiles
%         mz_pt = [i,mz_sf, mz_n0, mz_n5];
%         mz_perc = vertcat(mz_perc, mz_pt);
%     
%         xlabel('meal size');
%         xlim([0 20]);
%         ylim([0 60]);

    % meal time
    %ret_dat(er,:) = [];
    %ret_dat_2 = ret_dat;
    ret_dat = datetime(table2array(ret_dat),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');

    %plot pellets or meals over time
    figure(f2);
            subplot(6,1,i);histogram(ret_dat, "BinWidth", minutes(10));hold on;
            lim = 25;
            fill([start start+hours(12) start+hours(12) start],[0 0 40 40], 'k', FaceAlpha=0.1);
            n2= start+hours(24);n3= start+hours(24*2);n4= start+hours(24*3);n5= start+hours(24*4); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            n6= start+hours(24*5);n7= start+hours(24*6);n8= start+hours(24*7);n9= start+hours(24*8);
            fill([n2 n2+hours(12) n2+hours(12) n2],[0 0 lim lim], 'k', FaceAlpha=0.1);
            fill([n3 n3+hours(12) n3+hours(12) n3],[0 0 lim lim], 'k', FaceAlpha=0.1);
            fill([n4 n4+hours(12) n4+hours(12) n4],[0 0 lim lim], 'k', FaceAlpha=0.1); 
            fill([n5 n5+hours(12) n5+hours(12) n5],[0 0 lim lim], 'k', FaceAlpha=0.1);
            fill([n6 n6+hours(12) n6+hours(12) n6],[0 0 lim lim], 'k', FaceAlpha=0.1);
            fill([n7 n7+hours(12) n7+hours(12) n7],[0 0 lim lim], 'k', FaceAlpha=0.1);
            fill([n8 n8+hours(12) n8+hours(12) n8],[0 0 lim lim], 'k', FaceAlpha=0.1);
            fill([n9 n9+hours(12) n9+hours(12) n9],[0 0 lim lim], 'k', FaceAlpha=0.1);
            ylim([0 lim]);
            ylabel('pellets retrieved');
            title(animal);
            xlim([start stop]);

    %meals
    mtimes = datetime(table2array(meal_size(:,"Date_Time")));
    meal_size.Date_Time = datetime(meal_size.Date_Time);
    msize = table2array(meal_size(:,"sz"));
    mlas = table2array(meal_size(:,"Laser"));

    Otimes=Time(Optomeals);
    Ctimes =Time(Ctrlmeals);
    Osizes =table2array(events(Optomeals,"Pellet_number"));
    Csizes = table2array(events(Ctrlmeals,"Pellet_number"));
    

    figure(f3);
    subplot(6,1,i); 
        %scatter(mtimes,msize, ".");hold on;
        stem(Otimes, Osizes,'-g',LineWidth=2 );hold on;
        stem(Ctimes,Csizes,'-b',LineWidth=2);
        lim = 40;
        fill([start start+hours(12) start+hours(12) start],[0 0 40 40], 'k', FaceAlpha=0.1);
            n2= start+hours(24);n3= start+hours(24*2);n4= start+hours(24*3);n5= start+hours(24*4); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            n6= start+hours(24*5);n7= start+hours(24*6);n8= start+hours(24*7);n9= start+hours(24*8);
            fill([n2 n2+hours(12) n2+hours(12) n2],[0 0 lim lim], 'k', FaceAlpha=0.1);
            fill([n3 n3+hours(12) n3+hours(12) n3],[0 0 lim lim], 'k', FaceAlpha=0.1);
            fill([n4 n4+hours(12) n4+hours(12) n4],[0 0 lim lim], 'k', FaceAlpha=0.1); 
            fill([n5 n5+hours(12) n5+hours(12) n5],[0 0 lim lim], 'k', FaceAlpha=0.1);
            fill([n6 n6+hours(12) n6+hours(12) n6],[0 0 lim lim], 'k', FaceAlpha=0.1);
            fill([n7 n7+hours(12) n7+hours(12) n7],[0 0 lim lim], 'k', FaceAlpha=0.1);
            fill([n8 n8+hours(12) n8+hours(12) n8],[0 0 lim lim], 'k', FaceAlpha=0.1);
            fill([n9 n9+hours(12) n9+hours(12) n9],[0 0 lim lim], 'k', FaceAlpha=0.1);
            ylim([0 lim]);
            %ylabel('pellets retrieved');
            ylabel('meal size');
            title(animal);
            xlim([start stop]);



% Opto plots using code meal ends
    O = Opto_meal_ends(:,3);O(:,2) = 1;
    C = Ctrl_meal_ends(:,3);C(:,2) = 0;
    mlas = vertcat(O(:,2), C(:,2));
    bigm = vertcat(O(:,1), C(:,1));
%     if animal == 'A'
%             small_meals = find(bigm < 4);%%%%%%%%%%%%%%change meal size here
%     else 
%             small_meals = find(bigm < 5);
%     end 
%     %bigm(small_meals)=[];
%     %mlas(small_meals)= [];
    if hab == 1
       mlas(find(mlas == 1)) = 0; 
    end 
    if animal == "D"| animal == "E"
         col = "#4DBEEE";
    else
        col = "#77AC30";
    end 

    figure(f4);
    subplot(1,6,i);boxchart(mlas,bigm,BoxFaceColor=col);hold on;
        swarmchart(mlas,bigm,".");
        title(animal);
        ylabel('meal size');
        xticks([0 1]);xticklabels({'control', 'laser on'});
        ylim([0 30]);
        ctrl_mean = mean(bigm(find(mlas == 0)));
        las_mean = mean(bigm(find(mlas == 1)));
        plot([0 1],[ctrl_mean las_mean], 'k*');%xy 
%         edges= [3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
%         subplot(2,6,i+6);histogram(bigm(find(mlas == 0)),edges); hold on;
%         histogram(bigm(find(mlas == 1)),edges,FaceColor="#77AC30", FaceAlpha=0.3);


    %collect data for analysis
    opt = [bigm,mlas];
    for o = 1:size(opt,1)
        if i == 1 || i == 5
            opt(o,4) = 0;
        else
            opt(o,4) = 1; 
        end
        opt(o,3) = i;
    end

    Opto_data_all = vertcat(Opto_data_all,opt);

end 

%% power anlyses

OD_A = Opto_data_all(find(Opto_data_all(:,3) == 1),[1:2 4]);
OD_B = Opto_data_all(find(Opto_data_all(:,3) == 2),[1:2 4]);
OD_C = Opto_data_all(find(Opto_data_all(:,3) == 3),[1:2 4]);
OD_D = Opto_data_all(find(Opto_data_all(:,3) == 4),[1:2 4]);
OD_E = Opto_data_all(find(Opto_data_all(:,3) == 5),[1:2 4]);
OD_F = Opto_data_all(find(Opto_data_all(:,3) == 6),[1:2 4]);
Opto_data_all(:,3) = [];
%plot means
Means(1,1) = mean(OD_A(find(OD_A(:,2) ==1)));Means(7,1) = mean(OD_A(find(OD_A(:,2) ==0)));
Means(2,1) = mean(OD_B(find(OD_B(:,2) ==1)));Means(8,1) = mean(OD_B(find(OD_B(:,2) ==0)));
Means(3,1) = mean(OD_C(find(OD_C(:,2) ==1)));Means(9,1) = mean(OD_C(find(OD_C(:,2) ==0)));
Means(4,1) = mean(OD_D(find(OD_D(:,2) ==1)));Means(10,1) = mean(OD_D(find(OD_D(:,2) ==0)));
Means(5,1) = mean(OD_E(find(OD_E(:,2) ==1)));Means(11,1) = mean(OD_E(find(OD_E(:,2) ==0)));
Means(6,1) = mean(OD_F(find(OD_F(:,2) ==1)));Means(12,1) = mean(OD_F(find(OD_F(:,2) ==0)));
%std
Means(1,2) = std(OD_A(find(OD_A(:,2) ==1)));Means(7,2) = std(OD_A(find(OD_A(:,2) ==0)));
Means(2,2) = std(OD_B(find(OD_B(:,2) ==1)));Means(8,2) = std(OD_B(find(OD_B(:,2) ==0)));
Means(3,2) = std(OD_C(find(OD_C(:,2) ==1)));Means(9,2) = std(OD_C(find(OD_C(:,2) ==0)));
Means(4,2) = std(OD_D(find(OD_D(:,2) ==1)));Means(10,2) = std(OD_D(find(OD_D(:,2) ==0)));
Means(5,2) = std(OD_E(find(OD_E(:,2) ==1)));Means(11,2) = std(OD_E(find(OD_E(:,2) ==0)));
Means(6,2) = std(OD_F(find(OD_F(:,2) ==1)));Means(12,2) = std(OD_F(find(OD_F(:,2) ==0)));
Means(1:6,3) = 1;Means(7:12,3) = 0;

figure;scatter(Means(:,3),Means(:,1));hold on;
xlim([-1 2]);xticks([0 1]); xticklabels({"ctrl", "laser"});
plot([Means(7,3) Means(1,3)],[Means(7,1) Means(1,1)],"b");
plot([Means(8,3) Means(2,3)],[Means(8,1) Means(2,1)],"g");
plot([Means(9,3) Means(3,3)],[Means(9,1) Means(3,1)],"g");
plot([Means(10,3) Means(4,3)],[Means(10,1) Means(4,1)],"g");
plot([Means(11,3) Means(5,3)],[Means(11,1) Means(5,1)],"b");
plot([Means(12,3) Means(6,3)],[Means(12,1) Means(6,1)],"g");







Ctrl = vertcat(OD_A, OD_E);
C_las = Ctrl(find(Ctrl(:,2)==1));
C_nl = Ctrl(find(Ctrl(:,2)==0));

%only D and F
Arch = vertcat(OD_D, OD_F);
A_las = Arch(find(Arch(:,2)==1));
A_nl = Arch(find(Arch(:,2)==0));
% B C
Arch2 = vertcat(OD_B, OD_C);
A2_las = Arch2(find(Arch2(:,2)==1));
A2_nl = Arch2(find(Arch2(:,2)==0));
figure; boxchart(Arch2(:,2),Arch2(:,1));hold on;
swarmchart(Arch2(:,2),Arch2(:,1),".");

%all Arch animals
Arch_all = vertcat(OD_D, OD_F,OD_B, OD_C);
Aa_las = Arch_all(find(Arch_all(:,2)==1));
Aa_nl = Arch_all(find(Arch_all(:,2)==0));

%power analysis all ArchT animals
pwr = sampsizepwr('t',[mean(Aa_nl), std(Aa_nl)], mean(Aa_las), [], 282)
n = sampsizepwr('t',[mean(Aa_nl), std(Aa_nl)], mean(Aa_las), 0.8, [])
%pwr only D and F
pwr = sampsizepwr('t',[mean(A_nl), std(A_nl)], mean(A_las), [], 157)
n = sampsizepwr('t',[mean(A_nl), std(A_nl)], mean(A_las), 0.9, [])
%pwr only B and C
pwr = sampsizepwr('t',[mean(A2_nl), std(A2_nl)], mean(A2_las), [], 124)
n = sampsizepwr('t',[mean(A2_nl), std(A2_nl)], mean(A2_las), 0.8, [])


%[h,p] = ttest2(log10(Aa_nl),log10(Aa_las)) 

Aa_las = OD_C(find(OD_C(:,2)==1));
Aa_nl = OD_C(find(OD_C(:,2)==0));
pwr = sampsizepwr('t',[mean(Aa_nl), std(Aa_nl)], mean(Aa_las), [], 56)
n = sampsizepwr('t',[mean(Aa_nl), std(Aa_nl)], mean(Aa_las), 0.8, [])




%% loop for just data collection no plots
