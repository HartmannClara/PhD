% analysis of OPTo files from Optofiles struct
% plotting

%% loop for all animals
clear all; close all;
intervals = [29, 43, 22, 23, 20, 17, 43, 22, 31, 18, 14, 18, 15, 23, 15, 21, 12, 18];%batch1%batch2%batch3
inh= [3, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4];
f0 = figure;
f1 = figure;
f2 = figure;
allMt = table;
perc = [];mz_perc = [];
Opto_data_all =[];
load("Optodata.mat");
fn=fieldnames(Optodata);

%loop through hab data of batch one and correct experiment name to "hab"
%and meals to short meals, in b1 hab was performed with inh 1 code but
%fibers where unplugged for testing
for j = 1:6
    events = Optodata(j).hab;
    exp = table2array(events(:,"experiment")); i1 = find(exp == "inh_1");
    Optodata(j).hab.experiment(i1) = "hab";
    Optodata(j).hab.Laser(i1) = 0;
    cmeals = table2array(events(:,"Type")); Om = find(cmeals == "Opto meal ended");Cm = find(cmeals == "Control meal ended");
    Optodata(j).hab.Type(Om) = "Short meal ended";Optodata(j).hab.Type(Cm) = "Short meal ended";
end

%for i= [1:6 9:12]%inh_1
for i = 1:size(Optodata,2)  %inh_4
    i
    animal = Optodata(i).Label;
    %events = Optodata(i).inh1; %%%%%%%%%%%%%%%%%%%%%%%change here!!
    events = Optodata(i).inh4;
    habt = Optodata(i).hab;

    % if refeeding should also be analyzed with the inh_4 data:
    %RF = Optodata(i).RF;
    %events2 = vertcat(events,RF);

    %for animals that have both inh4 and inh3 expeiments, select which data
    %to use:
    if i == 7
        exp = table2array(events(:,"experiment"));
        %inh4 = find(exp == "inh_4");
        inh4 = find(exp == "inh_3");%exclude inh_3 data
        events(inh4,:) = [];
    elseif i == 8  
        exp = table2array(events(:,"experiment"));
        %inh4 = find(exp == "inh_4");
        inh4 = find(exp == "inh_3");%change%here!!!!!!!!!!!!!!!!!!!!!!
        events(inh4,:) = [];
    end

    if i == 1 || i == 2 || i == 3 || i == 4 || i == 5 || i == 6 
        start= datetime(['2024-03-04 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');%start day
    elseif i == 10 || i == 11   
        start= datetime(['2024-03-27 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');%start day
    elseif i == 7 || i == 8 || i == 9 || i == 12
        start= datetime(['2024-03-25 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');%start day
    else 
        start= datetime(['2024-06-14 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');%start day
    end
    days = 8;%how many days to plot 
    mid = start + day(1);
    stop = start + day (days);

    % events analysis
    Type = table2array(events(:,"Type"));
    Pellets = table2array(events(:,"Pellet_number"));
    Time = datetime(table2array(events(:,"Date_Time")),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    latency= table2array(events(:,'Latency'));
    %create meals independent of code 
    retrievals = find(Type == "Pellet_Retrieved");
    retrievals = events(retrievals,:);
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

%habituation data
    habType = table2array(habt(:,"Type"));
    habPellets = table2array(habt(:,"Pellet_number"));
    habTime = datetime(table2array(habt(:,"Date_Time")),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    hablatency= table2array(habt(:,'Latency'));
    %create meals independent of code 
    habretrievals = find(habType == "Pellet_Retrieved");
    habretrievals = habt(habretrievals,:);
    habret_lat = habretrievals(:,"Latency");
    habret_las = habretrievals(:,"Laser");
    habret_dat = habretrievals(:,"Date_Time");
    sz = 1; 
    habmeal_size = table;
    for j = 1:size(habret_lat,1)
        hablat = table2array(habret_lat(j,:));
        if hablat <= intervals(i) % smaller or equal to selected IPI in seconds
            sz = sz + 1; % size of the meal
        elseif j == 1
            habmeal = [habret_dat(j,:), array2table(sz) , habret_las(j,:) ];
            habmeal_size = [habmeal_size; habmeal];
            sz = 1;
        else 
            habmeal = [habret_dat(j-1,:), array2table(sz) , habret_las(j-1,:) ];
            habmeal_size = [habmeal_size; habmeal];
            sz = 1;
        end      
    end  

%     %Inter pellet intervals
%     IPI =  table2array(habret_lat);
%     %er= find(IPI < 2);% look for errors???
%     %IPI(er) = [];
%     IPI = round(IPI);
%     big = find(IPI > 90);
%     IPI(big) = [];
%     edges = [0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90];
% 
%     figure(f1); 
%     subplot(2,10,i); histogram(IPI,edges);hold on;
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
% 

    %plot latency to 4th, 5th or 6th pellet
    if Optodata(i).Label =='C'
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
    %collect data for plotting
    LatencyInh4(i).pel4 = lat_4;
    LatencyInh4(i).pel5 = lat_5;
    LatencyInh4(i).pel6 = lat_6;
    LatencyInh4(i).las4 = las_4;
    LatencyInh4(i).las5 = las_5;
    LatencyInh4(i).las6 = las_6;
    

% classify satiety level according ot latency to last meal
% cut meals and log single pellet trials if trial type consistent:

    Opto_meal_ends = find(Type == "Opto meal ended");
    Ctrl_meal_ends = find(Type == "Control meal ended");
    Short_meal_ends = find(Type == "Short meal ended");

    Hab_meal_ends = find(habType == "Short meal ended");

    Shortmeals = Short_meal_ends;
    Optomeals = Opto_meal_ends;
    Ctrlmeals = Ctrl_meal_ends;
    Hab_meals = Hab_meal_ends;

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
    for e = 1:size(Hab_meal_ends,1)
        en = Hab_meal_ends(e);
        ms = habPellets(en);
        Hab_meal_ends(e,3) = ms;
        Hab_meal_ends(e,2) = hablatency(en-ms);  
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
%     for m = 1:size(Opto_meal_ends,1)
%         om = Opto_meal_ends(m,1);
%         omtime = Time(om-(Opto_meal_ends(m,3)+1));
%         last30 = omtime - minutes(20);%%%%%%%%%%%%%%%%%%%%%%% time
%         last30time = events;
%         x=datetime(last30time.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
%         last30time(find(x>omtime),:)=[];%chops events to match the day
%         x=datetime(last30time.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
%         last30time(find(x<last30),:)=[];
%         pell30 = size(find(table2array(last30time(:,2)) == "Pellet_Retrieved"),1);
%         Opto_meal_ends(m,4) = pell30; %save satiety status
%         if pell30 < 6
%             Opto_meal_ends(m,5) = 4;%hungry
%         elseif 5 < pell30 & pell30 < 8
%             Opto_meal_ends(m,5) = 2;
%         else 
%             Opto_meal_ends(m,5) = 2;%sated
%         end        
%     end
% 
%     for m = 1:size(Ctrl_meal_ends,1)
%         om = Ctrl_meal_ends(m,1);
%         omtime = Time(om-Ctrl_meal_ends(m,3));;
%         last30 = omtime - minutes(20);
%         last30time = events;
%         x=datetime(last30time.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
%         last30time(find(x>omtime),:)=[];%chops events to match the day
%         x=datetime(last30time.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
%         last30time(find(x<last30),:)=[];
%         pell30 = size(find(table2array(last30time(:,2)) == "Pellet_Retrieved"),1);
%         Ctrl_meal_ends(m,4) = pell30; %save satiety status
%         if pell30 < 6
%             Ctrl_meal_ends(m,5) = 3;%hungry
%         elseif 5 < pell30 & pell30 < 8
%             Ctrl_meal_ends(m,5) = 1;
%         else 
%             Ctrl_meal_ends(m,5) = 1; %sated
%         end    
%             
%     end
    % plot histogram of satiety measures
%     figure(f0);
%         subplot(2,12,i);histogram(Ctrl_meal_ends(:,4),[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26]);hold on;
%         %sf = prctile(Ctrl_meal_ends(:,4),75);b=30;plot([sf sf],[0 b]);
%         %n0 = prctile(Ctrl_meal_ends(:,4),90);b=30;plot([n0 n0],[0 b]);
%         %n5 = prctile(Ctrl_meal_ends(:,4),95);b=30;plot([n5 n5],[0 b]);
%         %title(animal);
%         %xlabel('Nr pellets 20 min prior (Ctrl)')
%         ylim([0 40]);
%         subplot(2,12,i+12);histogram(Opto_meal_ends(:,4),[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26]);hold on;
%         %sf = prctile(Opto_meal_ends(:,4),75);b=30;plot([sf sf],[0 b]);
%         %n0 = prctile(Opto_meal_ends(:,4),90);b=30;plot([n0 n0],[0 b]);
%         %n5 = prctile(Opto_meal_ends(:,4),95);b=30;plot([n5 n5],[0 b]);
%         %xlabel('Nr pellets 20 min prior (ArchT)');
%         ylim([0 40]);

%     figure(f2);
%         subplot(2,12,i); boxchart(Ctrl_meal_ends(:,5),Ctrl_meal_ends(:,3));hold on;
%         swarmchart(Ctrl_meal_ends(:,5),Ctrl_meal_ends(:,3),".");
%         boxchart(Opto_meal_ends(:,5),Opto_meal_ends(:,3),BoxFaceColor=col);
%         swarmchart(Opto_meal_ends(:,5),Opto_meal_ends(:,3),".");
%         title(animal);
%         ylim([0 25]);
%         xlim([0 5]);
%         ylabel("meal size")
%         xticks([1 2 3 4]);
%         xticklabels({'sated ctrl','sated lon', 'hungry ctrl', 'hungry lon'});
%         edges = [0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 60];
%         subplot(2,6,i+6);histogram(Opto_meal_ends(:,4),edges);

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
%     figure(f2);
%             subplot(12,1,i);histogram(ret_dat, "BinWidth", minutes(10));hold on;
%             lim = 25;
%             fill([start start+hours(12) start+hours(12) start],[0 0 40 40], 'k', FaceAlpha=0.1);
%             n2= start+hours(24);n3= start+hours(24*2);n4= start+hours(24*3);n5= start+hours(24*4); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             n6= start+hours(24*5);n7= start+hours(24*6);n8= start+hours(24*7);n9= start+hours(24*8);
%             fill([n2 n2+hours(12) n2+hours(12) n2],[0 0 lim lim], 'k', FaceAlpha=0.1);
%             fill([n3 n3+hours(12) n3+hours(12) n3],[0 0 lim lim], 'k', FaceAlpha=0.1);
%             fill([n4 n4+hours(12) n4+hours(12) n4],[0 0 lim lim], 'k', FaceAlpha=0.1); 
%             fill([n5 n5+hours(12) n5+hours(12) n5],[0 0 lim lim], 'k', FaceAlpha=0.1);
%             fill([n6 n6+hours(12) n6+hours(12) n6],[0 0 lim lim], 'k', FaceAlpha=0.1);
%             fill([n7 n7+hours(12) n7+hours(12) n7],[0 0 lim lim], 'k', FaceAlpha=0.1);
%             fill([n8 n8+hours(12) n8+hours(12) n8],[0 0 lim lim], 'k', FaceAlpha=0.1);
%             fill([n9 n9+hours(12) n9+hours(12) n9],[0 0 lim lim], 'k', FaceAlpha=0.1);
%             ylim([0 lim]);
%             ylabel('pellets retrieved');
%             title(animal);
%             xlim([start stop]);

    %meals
    mtimes = datetime(table2array(meal_size(:,"Date_Time")));
    meal_size.Date_Time = datetime(meal_size.Date_Time);
    msize = table2array(meal_size(:,"sz"));
    mlas = table2array(meal_size(:,"Laser"));

    Otimes=Time(Optomeals);
    Ctimes =Time(Ctrlmeals);
    Osizes =table2array(events(Optomeals,"Pellet_number"));
    Csizes = table2array(events(Ctrlmeals,"Pellet_number"));
    

%     figure(f3);
%     subplot(12,1,i); 
%         %scatter(mtimes,msize, ".");hold on;
%         stem(Otimes, Osizes,'-g',LineWidth=2 );hold on;
%         stem(Ctimes,Csizes,'-b',LineWidth=2);
%         lim = 40;
%         fill([start start+hours(12) start+hours(12) start],[0 0 40 40], 'k', FaceAlpha=0.1);
%             n2= start+hours(24);n3= start+hours(24*2);n4= start+hours(24*3);n5= start+hours(24*4); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             n6= start+hours(24*5);n7= start+hours(24*6);n8= start+hours(24*7);n9= start+hours(24*8);
%             fill([n2 n2+hours(12) n2+hours(12) n2],[0 0 lim lim], 'k', FaceAlpha=0.1);
%             fill([n3 n3+hours(12) n3+hours(12) n3],[0 0 lim lim], 'k', FaceAlpha=0.1);
%             fill([n4 n4+hours(12) n4+hours(12) n4],[0 0 lim lim], 'k', FaceAlpha=0.1); 
%             fill([n5 n5+hours(12) n5+hours(12) n5],[0 0 lim lim], 'k', FaceAlpha=0.1);
%             fill([n6 n6+hours(12) n6+hours(12) n6],[0 0 lim lim], 'k', FaceAlpha=0.1);
%             fill([n7 n7+hours(12) n7+hours(12) n7],[0 0 lim lim], 'k', FaceAlpha=0.1);
%             fill([n8 n8+hours(12) n8+hours(12) n8],[0 0 lim lim], 'k', FaceAlpha=0.1);
%             fill([n9 n9+hours(12) n9+hours(12) n9],[0 0 lim lim], 'k', FaceAlpha=0.1);
%             ylim([0 lim]);
%             %ylabel('pellets retrieved');
%             ylabel('meal size');
%             title(animal);
%             xlim([start stop]);



% Opto plots using code meal ends
    H = Hab_meal_ends(:,3);H(:,2) = -1;
    O = Opto_meal_ends(:,3);O(:,2) = 1;
    C = Ctrl_meal_ends(:,3);C(:,2) = 0;
    mlas = vertcat(O(:,2), C(:,2), H(:,2));
    bigm = vertcat(O(:,1), C(:,1), H(:,1));


    %collect data for analysis
    opt = [bigm,mlas];
    for o = 1:size(opt,1)
        if Optodata(i).Label == "C" 
            opt(o,4) = 0;
        else
            opt(o,4) = 1; 
        end
        opt(o,3) = i;
    end

    Opto_data_all = vertcat(Opto_data_all,opt);% meal size, meal type (-1 hab, 0 ctrl, 1 laser), animal, ctrl or arch(1) 
    Opto_data_n(i,1) = size(find(mlas == 0),1);%control
    Opto_data_n(i,2) = size(find(mlas == 1),1);%laser

end 

%% manual corrections data
%remove error meals (really high numbers)
Opto_data_all(Opto_data_all(:,1)>50,:) = [];

% for inh_4
%if inh_4 of animal B and C are used combine with data from b1:
Opto_data_n(2,:)= Opto_data_n(2,:) + Opto_data_n(7,:);
Opto_data_n(3,:)= Opto_data_n(3,:) + Opto_data_n(8,:);
Opto_data_n(7:8,:) = [];
mean(Opto_data_n(:,1)); std(Opto_data_n(:,1));
mean(Opto_data_n(:,2)); std(Opto_data_n(:,2));

%plot opto boxplots outside of loop to combine data of animals measured in both
%batches
Opto_data_cmb = Opto_data_all;

%if inh_3 for B and C in batch 2 is used
% Opto_data_cmb(find(Opto_data_cmb(:,3)== 2),:) =[];
% Opto_data_cmb(find(Opto_data_cmb(:,3)== 3),:) =[];
% Opto_data_all(find(Opto_data_all(:,3)== 2),:) =[];
% Opto_data_all(find(Opto_data_all(:,3)== 3),:) =[];

%correct meal sizes by replacing 7 and 8  with new labels 2 and 3
animB2 = find(Opto_data_cmb(:,3)== 7);Opto_data_cmb(animB2,3) = 2;
animC2 = find(Opto_data_cmb(:,3)== 8);Opto_data_cmb(animC2,3) = 3;

%habituation data
% correct hab data to remove small meals, inh_4 or inh_1 already only has
% meals eaqual or bigger than 3
smallhab = find(Opto_data_cmb(:,1) < 3);%rem all meals smaller than 3
Opto_data_cmb(smallhab,:) = [];
a1end = (find(Opto_data_cmb(:,3) == 1));

%smallhab2 = find(Opto_data_cmb(a1end(end)+1:end,1) < 4);%rem all meals smaller than 4 from all but animal 1
%Opto_data_cmb(smallhab2+a1end(end),:) = [];



%%
%remove hab data
Opto_data_cmb(Opto_data_cmb(:,2)==-1,:) = [];

%remove error meals (really high numbers)
Opto_data_cmb(Opto_data_cmb(:,1)>50,:) = [];

f4 = figure;
figure(f4);
    %anim  = [1 2 3 4 5 6 7 8 9 10]
%         pl = [1,5,6,7,2,8,9,3,4,10];
%         name = {'Ab1f','Eb1m','Db2m','Eb2f','Bb1f','Cb1f','Db1f','Fb1m','Cb2m','Fb2f'};
% b3 [11,12,13,14,15,16] [Ab3m, Bb3m, Cb3f, Db3f, Eb3m, Fb3m] [A C A C C C]
        %Optodata_corr2 = Optodata; Optodata_corr(7:8) = [];
        %IDs = {'Ab1f','Bb1f','Cb1f','Db1f','Eb1m','Fb1m','Cb2m','Db2m','Eb2f','Fb2f','Ab3m','Bb3m','Cb3f','Db3f','Eb3m','Fb3m'};

        %Optodata_corr2.ID = {'Ab1f','Bb1f','Cb1f','Db1f','Eb1m','Fb1m','Cb2m','Db2m','Eb2f','Fb2f','Ab3m','Bb3m','Cb3f','Db3f','Eb3m','Fb3m'};
        IDs = {'Ab1f','Bb1f','Cb1f','Db1f','Eb1m','Fb1m','Cb2m','Db2m','Eb2f','Fb2f','Ab3m','Bb3m','Cb3f','Db3f','Eb3m','Fb3m'};%org order
        name = {'Ab1f','Eb1m','Db2m','Eb2f','Bb3m','Db3f','Eb3m','Fb3m','Bb1f','Cb1f','Db1f','Fb1m','Cb2m','Fb2f', 'Ab3m','Cb3f'};% order i want to plot


        %ctrl = {'Ab1f','Eb1m','Db2m','Eb2f','Bb3m','Db3f','Eb3m','Fb3m'} % [1 5 8 9 12 14 15 16]
        %arch = {'Bb1f','Cb1f','Db1f','Fb1m','Cb2m','Fb2f', 'Ab3m','Cb3f'} % [2 3 4 6 7 10 11 13]

        i=1;
        for a = [1:6,9:18]
            mlas = Opto_data_cmb((find(Opto_data_cmb(:,3)==a)),2);%laser
            bigm = Opto_data_cmb((find(Opto_data_cmb(:,3)==a)),1);%meals
            Group = Opto_data_cmb((find(Opto_data_cmb(:,3) == a)),4);

            if Group(1,1) == 0
                col = "#4DBEEE";
            else
                col = "#77AC30";
            end 

            pl = find(contains(name,IDs(i)));
            tit = name(pl);

            subplot(1,16,pl);boxchart(mlas,bigm,BoxFaceColor=col);hold on;
            swarmchart(mlas,bigm,".");
            ylabel('meal size');

            title(tit);
            
            xticks([-1 0 1]);xticklabels({'hab','control', 'laser on'});
            ylim([0 30]);
            ctrl_mean = mean(bigm(find(mlas == 0)));
            las_mean = mean(bigm(find(mlas == 1)));
            hab_mean = mean(bigm(find(mlas == -1)));
            plot([-1 0 1],[hab_mean ctrl_mean las_mean], 'k*');%xy 
            i=i+1;
        end   
%% make one plot for control and one for arch
Optoanimals = Opto_data_cmb(find(Opto_data_cmb(:,4)== 1),:);
Ctrlanimals = Opto_data_cmb(find(Opto_data_cmb(:,4)== 0),:);

figure(f4);
subplot(1,2,1); boxchart(Ctrlanimals(:,2),Ctrlanimals(:,1),BoxFaceColor="#4DBEEE");hold on;
ylim([3 25]);xticks([0 1]);xticklabels({'control', 'laser on'});
subplot(1,2,2); boxchart(Optoanimals(:,2),Optoanimals(:,1),BoxFaceColor="#77AC30");hold on;
ylim([3 25]);xticks([0 1]);xticklabels({'control', 'laser on'});


%% just means with std connected to matched no laser trials
f4 = figure;
figure(f4);
    %anim  = [1 2 3 4 5 6 7 8 9 10]
        %pl = [1,5,6,7,2,8,9,3,4,10];
        %name = {'Ab1f','Eb1m','Db2m','Eb2f','Bb1f','Cb1f','Db1f','Fb1m','Cb2m','Fb2f'};
        i=1;
        for a = [1:6,9:18]
            mlas = Opto_data_cmb((find(Opto_data_cmb(:,3)==a)),2);%laser
            bigm = Opto_data_cmb((find(Opto_data_cmb(:,3)==a)),1);%meals
            Group = Opto_data_cmb((find(Opto_data_cmb(:,3) == a)),4);

            ctrl_mean = mean(bigm(find(mlas == 0)));
            las_mean = mean(bigm(find(mlas == 1)));

            optomice(i,1)=ctrl_mean;
            optomice(i,2)=las_mean;
            optomice(i,3)=Group(1,1);
            i=i+1;
            if Group(1,1) == 0
                plt = 1;%plot 1 control
                subplot(1,2,1);plot([0 1],[ctrl_mean las_mean], 'o',"Color","#4DBEEE",'LineWidth',2);hold on;
                plot([0 1],[ctrl_mean las_mean], '-',"Color","#4DBEEE",'LineWidth',2);
                ylabel('meal size');
                title('Control');
                xticks([0 1]);xticklabels({'control', 'laser on'});
                ylim([3 10]);
                xlim([-0.5 1.5]);
            else
                plt = 2;%plot 2 opto
                subplot(1,2,2);plot([0 1],[ctrl_mean las_mean], 'o', "Color","#77AC30", 'LineWidth',2); hold on;
                plot([0 1],[ctrl_mean las_mean], '-', "Color","#77AC30",'LineWidth',2);
                ylabel('meal size');
                title('ArchT');
                xticks([0 1]);xticklabels({'control', 'laser on'});
                ylim([3 10]);
                xlim([-0.5 1.5]);
            end    

        end   

%% TRY QUICK LME

DD = Opto_data_cmb(Opto_data_cmb(:,2)~=(-1),:);%remove all hab data
%correct meal sizes by replacing 7 and 8  with new labels 2 and 3
animB2 = find(Opto_data_cmb(:,3)== 7);Opto_data_cmb(animB2,3) = 2;
animC2 = find(Opto_data_cmb(:,3)== 8);Opto_data_cmb(animC2,3) = 3;

TT = table(DD(:,1),categorical(DD(:,2)),categorical(DD(:,3)),categorical(DD(:,4)),'VariableNames',{'MealSize','Opto','Mouse','Opsin'});

lme = fitlme(TT,'MealSize ~ Opsin*Opto + (1|Opsin:Mouse)');
disp(lme);

% vizualize?
%% 
DD = Opto_data_cmb(Opto_data_cmb(:,2)~=(-1),:);%remove all hab data

%correct meal sizes by replacing 7 and 8  with new labels 2 and 3
animB2 = find(Opto_data_cmb(:,3)== 7);Opto_data_cmb(animB2,3) = 2;
animC2 = find(Opto_data_cmb(:,3)== 8);Opto_data_cmb(animC2,3) = 3;

TT = table(DD(:,1),categorical(DD(:,2)),categorical(DD(:,3)),categorical(DD(:,4)),'VariableNames',{'MealSize','Opto','Mouse','Opsin'});

lme = fitlme(TT,'MealSize ~ Opsin*Opto + (1|Opsin:Mouse)');
disp(lme);
%% anova
% Assuming you have a table named 'data' with columns: 'MealSize', 'Opto', 'Opsin', and 'Mouse'
DD = Opto_data_cmb(Opto_data_cmb(:,2)~=(-1),:);%remove all hab data
data = table(DD(:,1),categorical(DD(:,2)),categorical(DD(:,3)),categorical(DD(:,4)),'VariableNames',{'MealSize','Opto','Mouse','Opsin'});

% Compute mean meal size for each animal in control and laser conditions
controlMeans = grpstats(data(data.Opto == categorical(0), :), 'Mouse', {'mean'}, 'DataVars', 'MealSize');
controlMeans.Properties.VariableNames = {'Mouse','N' ,'ControlMean'};

laserMeans = grpstats(data(data.Opto == categorical(1), :), 'Mouse', {'mean'}, 'DataVars', 'MealSize');
laserMeans.Properties.VariableNames = {'Mouse','N',  'LaserMean'};

%
p = anovan(data.MealSize, {data.Opto,data.Opsin})
p2 = anovan(data.MealSize, {data.Opto,data.Opsin},"model","interaction",'varnames',{'data.Opto','data.Opsin'})

%%
% Perform n-way ANOVA to compare the effect of 'Opto' and 'Opsin' on 'MealSize'
[p_anova, tbl_anova, stats_anova] = anovan(data.MealSize, {data.Opto, data.Opsin}, ...
    'varnames', {'Opto', 'Opsin'}, 'model', 'full');

% Display ANOVA table
disp(tbl_anova);

% post hoc multiple comparison tests (Tukey's HSD)
%comparison_results = multcompare(stats_anova, 'Dimension', [1 2], 'CType', "hsd");
%disp(comparison_results);

% post hoc multiple comparison tests (Bonferroni correction)
[results,means,~,gnames] = multcompare(stats_anova, 'Dimension', [1 2], 'CType', 'bonferroni');
tbl = array2table(results,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl.("Group A")=gnames(tbl.("Group A"));
tbl.("Group B")=gnames(tbl.("Group B"))


%% ttests
GFPmice= optomice(optomice(:,3) == 0,1:2);% column 1 is off 2 is laser on
Archmice= optomice(optomice(:,3) == 1,1:2);%ArchT


% Perform paired t-test on Opsin 0 group
[~, p_val_opsin0] = ttest(GFPmice(:,2), GFPmice(:,1));

% Perform paired t-test on Opsin 1 group
[~, p_val_opsin1] = ttest(Archmice(:,2), Archmice(:,1));

% Display p-values for Opsin 0 and Opsin 1 groups
disp(['P-value for GFP: ', num2str(p_val_opsin0)]);
disp(['P-value for Arch: ', num2str(p_val_opsin1)]);

%% normality?
% Sample paired data
x = Archmice(:,1);
y = Archmice(:,2);

% Step 1: Calculate Differences
differences = x - y;
% Step 2: Test for Normality of Differences
[h, p] = adtest(differences); %anderson darling test

% Output results of normality test
if h == 0
    disp('The differences are normally distributed (fail to reject null hypothesis)');
else
    disp('The differences are not normally distributed (reject null hypothesis)');
end
%% nested anova?
% Perform nested ANOVA
[p_anova, tbl_anova, stats_anova] = anovan(data.MealSize, {data.Opsin, data.Opto}, ...%Opsin is between subjects, Opto is within, Mouse is random
    'varnames', {'Opto', 'Opsin'}, 'model',2,'nested', [0 0; 1 0]);

% Display ANOVA results
disp(tbl_anova);


%%
DD2 = Opto_data_cmb(Opto_data_cmb(:,2)~=(-1),:);%remove all hab data
for i = 1:size(DD2,1)
    if DD2(i,3) == 5 | DD2(i,3) == 6 | DD2(i,3) == 9 | DD2(i,3) == 10
        DD2(i,5) = 1;%male
    else
        DD2(i,5) = 0;%female
    end    
end
for i = 1:size(DD2,1)
    if DD2(i,3) == 1 | DD2(i,3) == 2 | DD2(i,3) == 3 | DD2(i,3) == 4 | DD2(i,3) == 5 | DD2(i,3) == 6
        DD2(i,6) = 1;%batch1
    else
        DD2(i,6) = 2;%batch2
    end   
end


TT2 = table(DD2(:,1),categorical(DD2(:,2)),categorical(DD2(:,3)),categorical(DD2(:,4)),categorical(DD2(:,5)),categorical(DD2(:,6)),'VariableNames',{'MealSize','Opto','Mouse','Opsin','Sex','Batch'});

lme2 = fitlme(TT2,'MealSize ~ Opsin*Opto + (1|Mouse) + (1|Sex) + (1|Batch)');
disp(lme2);


% compare
results = compare(lme,lme2)

%% hab ipi plots
f1=figure;
f2=figure;
perc=[];
%order: A B C D E F C2 D2 E2 F2; C O O O C O   O C C O
intervals = [29, 43, 22, 23, 20, 17, 43, 22, 31, 18, 14, 18];%batch1%batch2
inh= [3, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4];

for i = [1:6 9:12]
    if i == 1 || i == 2 || i == 3 || i == 4 || i == 5 || i == 6 
        chop_from= datetime(['2024-03-02 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');%start day
    else
        chop_from= datetime(['2024-03-23 10:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');%start day
    end
   
    habt = Optodata(i).hab;
    chop_to = chop_from + day(1);

    x=datetime(habt.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    habt(find(x<chop_from),:)=[];%chops events to match the day
    x=datetime(habt.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    habt(find(x>chop_to),:)=[];

    habType = table2array(habt(:,"Type"));
    habPellets = table2array(habt(:,"Pellet_number"));
    habTime = datetime(table2array(habt(:,"Date_Time")),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    hablatency= table2array(habt(:,'Latency'));
    %create meals independent of code 
    habretrievals = find(habType == "Pellet_Retrieved");
    habretrievals = habt(habretrievals,:);
    habret_lat = habretrievals(:,"Latency");
    habret_las = habretrievals(:,"Laser");
    habret_dat = habretrievals(:,"Date_Time");
    sz = 1; 
    habmeal_size = table;
    for j = 1:size(habret_lat,1)
        hablat = table2array(habret_lat(j,:));
        if hablat <= intervals(i) % smaller or equal to selected IPI in seconds
            sz = sz + 1; % size of the meal
        elseif j == 1
            habmeal = [habret_dat(j,:), array2table(sz) , habret_las(j,:) ];
            habmeal_size = [habmeal_size; habmeal];
            sz = 1;
        else 
            habmeal = [habret_dat(j-1,:), array2table(sz) , habret_las(j-1,:) ];
            habmeal_size = [habmeal_size; habmeal];
            sz = 1;
        end      
    end  

    %Inter pellet intervals
    IPI =  table2array(habret_lat);
    er= find(IPI < 2);% look for errors???
    IPI(er) = [];
    IPI = round(IPI);
    big = find(IPI > 90);
    IPI(big) = [];
    edges = [0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90];

    figure(f1); 
    subplot(2,10,i); histogram(IPI,edges);hold on;
        sf = prctile(IPI,75);b=150;%plot([sf sf],[0 b]);
        n0 = prctile(IPI,90);b=150;plot([n0 n0],[0 b]);
        n5 = prctile(IPI,95);b=150;%plot([n5 n5],[0 b]);
        %collect percentiles
        pt = [i,sf, n0, n5];
        perc = vertcat(perc, pt);
        title(animal);
        xlabel('IPI (5s bins)');
        xlim([0 90]);
        ylim([0 160]);

    figure(f2); 
    edges2 = [1:30];
    subplot(2,10,i); histogram((habmeal_size.sz),edges2);hold on;
        title(animal);
        xlabel('IPI (5s bins)');
        xlim([0 30]);
        ylim([0 50]);

end    


%% latency to next pellet
LatencyInh4raw = LatencyInh4;
fn = fieldnames(LatencyInh4);%combine values of animals measured in both batches(2,3 is 7,8)
for i = 1:size(fn,1)
    LatencyInh4(2).(fn{i})= vertcat(LatencyInh4(2).(fn{i}),LatencyInh4(7).(fn{i}));
    LatencyInh4(3).(fn{i})= vertcat(LatencyInh4(3).(fn{i}),LatencyInh4(8).(fn{i}));
end    
LatencyInh4(7:8) = [];

%% plot
f5 = figure;
pl = [1,5,6,7,2,8,9,3,4,10];
name = {'Ab1f','Eb1m','Db2m','Eb2f','Bb1f','Cb1f','Db1f','Fb1m','Cb2m','Fb2f'};
for i = [1:10]
    lat_4 = LatencyInh4(i).pel4;las_4 = LatencyInh4(i).las4;
    lat_5 = LatencyInh4(i).pel5;las_5 = LatencyInh4(i).las5;
    lat_6 = LatencyInh4(i).pel6;las_6 = LatencyInh4(i).las6;

    if pl(i) == 1 ||pl(i) == 2||pl(i) == 3||pl(i) == 4
        col = "#4DBEEE";

    else
        col = "#77AC30";
    end 

    figure(f5);
    subplot(3,10,pl(i));boxchart(las_4,lat_4,BoxFaceColor=col);hold on;
        %swarmchart(las_4,lat_4,".");
        title(name(pl(i)));
        ylabel('latency to 4th pellet');
        xticks([0 1]);xticklabels({'control', 'laser on'});
        ylim([0 30]);
        ctrl_mean = mean(lat_4(find(las_4 == 0)));
        las_mean = mean(lat_4(find(las_4 == 1)));
        plot([0 1],[ctrl_mean las_mean], 'k*');%xy

    subplot(3,10,pl(i)+10);boxchart(las_5,lat_5,BoxFaceColor=col);hold on;
        %swarmchart(las_5,lat_5,".");
        ylabel('latency to 5th pellet');
        xticks([0 1]);xticklabels({'control', 'laser on'});
        ylim([0 30]);
        ctrl_mean = mean(lat_5(find(las_5 == 0)));
        las_mean = mean(lat_5(find(las_5 == 1)));
        plot([0 1],[ctrl_mean las_mean], 'k*');%xy

    subplot(3,10,pl(i)+20);boxchart(las_6,lat_6,BoxFaceColor=col);hold on;
        %swarmchart(las_6,lat_6,".");
        ylabel('latency to 6th pellet');
        xticks([0 1]);xticklabels({'control', 'laser on'});
        ylim([0 30]);
        ctrl_mean = mean(lat_6(find(las_6 == 0)));
        las_mean = mean(lat_6(find(las_6 == 1)));
        plot([0 1],[ctrl_mean las_mean], 'k*');%xy
    
end
%% plot only averg for each animal
f5 = figure;
pl = [1,5,6,7,2,8,9,3,4,10];
name = {'Ab1f','Eb1m','Db2m','Eb2f','Bb1f','Cb1f','Db1f','Fb1m','Cb2m','Fb2f'};
figure(f5);
for i = [1:10]
    lat_4 = LatencyInh4(i).pel4;las_4 = LatencyInh4(i).las4;
    lat_5 = LatencyInh4(i).pel5;las_5 = LatencyInh4(i).las5;
    lat_6 = LatencyInh4(i).pel6;las_6 = LatencyInh4(i).las6;
%means
    ctrl_mean4 = mean(lat_4(find(las_4 == 0)));
    las_mean4 = mean(lat_4(find(las_4 == 1)));
    ctrl_mean5 = mean(lat_4(find(las_5 == 0)));
    las_mean5 = mean(lat_4(find(las_5 == 1)));
    ctrl_mean6 = mean(lat_4(find(las_6 == 0)));
    las_mean6 = mean(lat_4(find(las_6 == 1)));


    if pl(i) == 1 ||pl(i) == 2||pl(i) == 3||pl(i) == 4
        col = "#4DBEEE";
        subplot(2,3,1);plot([0 1],[ctrl_mean4 las_mean4], 'o',"Color","#4DBEEE",'LineWidth',2);hold on;
                plot([0 1],[ctrl_mean4 las_mean4], '-',"Color","#4DBEEE",'LineWidth',2);
                ylabel('latency to 4th pellet');
                title('Control');
                xticks([0 1]);xticklabels({'laser off', 'laser on'});
                ylim([0 15]);
                xlim([-0.5 1.5]);
                xline(0);xline(1);
        subplot(2,3,2);plot([0 1],[ctrl_mean5 las_mean5], 'o',"Color","#4DBEEE",'LineWidth',2);hold on;
                plot([0 1],[ctrl_mean5 las_mean5], '-',"Color","#4DBEEE",'LineWidth',2);
                ylabel('latency to 5th pellet');
                %title('Control');
                xticks([0 1]);xticklabels({'laser off', 'laser on'});
                ylim([0 15]);
                xlim([-0.5 1.5]);  xline(0);xline(1);      
        subplot(2,3,3);plot([0 1],[ctrl_mean6 las_mean6], 'o',"Color","#4DBEEE",'LineWidth',2);hold on;
                plot([0 1],[ctrl_mean6 las_mean6], '-',"Color","#4DBEEE",'LineWidth',2);
                ylabel('latency to 6th pellet');
                %title('Control');
                xticks([0 1]);xticklabels({'laser off', 'laser on'});
                ylim([0 15]);
                xlim([-0.5 1.5]);xline(0);xline(1);
    else
        col = "#77AC30";
        subplot(2,3,4);plot([0 1],[ctrl_mean4 las_mean4], 'o', "Color","#77AC30", 'LineWidth',2); hold on;
                plot([0 1],[ctrl_mean4 las_mean4], '-', "Color","#77AC30",'LineWidth',2);
                ylabel('latency to 4th pellet');
                title('ArchT');
                xticks([0 1]);xticklabels({'laser off', 'laser on'});
                ylim([0 15]);
                xlim([-0.5 1.5]);xline(0);xline(1);
        subplot(2,3,5);plot([0 1],[ctrl_mean5 las_mean5], 'o', "Color","#77AC30", 'LineWidth',2); hold on;
                plot([0 1],[ctrl_mean5 las_mean5], '-', "Color","#77AC30",'LineWidth',2);
                ylabel('latency to 5th pellet');
                %title('ArchT');
                xline(0);xline(1);
                xticks([0 1]);xticklabels({'laser off', 'laser on'});
                ylim([0 15]);
                xlim([-0.5 1.5]);xline(0);xline(1);
        subplot(2,3,6);plot([0 1],[ctrl_mean6 las_mean6], 'o', "Color","#77AC30", 'LineWidth',2); hold on;
                plot([0 1],[ctrl_mean6 las_mean6], '-', "Color","#77AC30",'LineWidth',2);
                ylabel('latency to 6th pellet');
                %title('ArchT');
                xticks([0 1]);xticklabels({'laser off', 'laser on'});
                ylim([0 15]);
                xlim([-0.5 1.5]); xline(0);xline(1);       
    end     
end
%% averg latency to next pellet plot

fn = fieldnames(LatencyInh4);%combine values of animals measured in both batches(2,3 is 7,8)
for i = 1:size(fn,1)
    LatencyInh4MeanC.(fn{i}) = vertcat(LatencyInh4(1).(fn{i}),LatencyInh4(5).(fn{i}),LatencyInh4(8).(fn{i}),LatencyInh4(4).(fn{i}),LatencyInh4(9).(fn{i}));
    LatencyInh4MeanA.(fn{i}) = vertcat(LatencyInh4(2).(fn{i}),LatencyInh4(3).(fn{i}),LatencyInh4(4).(fn{i}),LatencyInh4(6).(fn{i}),LatencyInh4(7).(fn{i}),LatencyInh4(10).(fn{i}));
end    

lat_4 = LatencyInh4MeanC.pel4;las_4 = LatencyInh4MeanC.las4;
lat_5 = LatencyInh4MeanC.pel5;las_5 = LatencyInh4MeanC.las5;
lat_6 = LatencyInh4MeanC.pel6;las_6 = LatencyInh4MeanC.las6;

f6 = figure;
figure(f6);
col = "#4DBEEE";
    subplot(3,2,1);boxchart(las_4,lat_4,BoxFaceColor=col);hold on;
        %swarmchart(las_4,lat_4,".");
        ylabel('latency to 4th pellet');
        xticks([0 1]);xticklabels({'control', 'laser on'});
        ylim([0 30]);
        ctrl_mean = mean(lat_4(find(las_4 == 0)));
        las_mean = mean(lat_4(find(las_4 == 1)));
        plot([0 1],[ctrl_mean las_mean], 'k*');%xy

    subplot(3,2,3);boxchart(las_5,lat_5,BoxFaceColor=col);hold on;
        %swarmchart(las_5,lat_5,".");
        ylabel('latency to 5th pellet');
        xticks([0 1]);xticklabels({'control', 'laser on'});
        ylim([0 30]);
        ctrl_mean = mean(lat_5(find(las_5 == 0)));
        las_mean = mean(lat_5(find(las_5 == 1)));
        plot([0 1],[ctrl_mean las_mean], 'k*');%xy

    subplot(3,2,5);boxchart(las_6,lat_6,BoxFaceColor=col);hold on;
        %swarmchart(las_6,lat_6,".");
        ylabel('latency to 6th pellet');
        xticks([0 1]);xticklabels({'control', 'laser on'});
        ylim([0 30]);
        ctrl_mean = mean(lat_6(find(las_6 == 0)));
        las_mean = mean(lat_6(find(las_6 == 1)));
        plot([0 1],[ctrl_mean las_mean], 'k*');%xy


%archT

lat_4 = LatencyInh4MeanA.pel4;las_4 = LatencyInh4MeanA.las4;
lat_5 = LatencyInh4MeanA.pel5;las_5 = LatencyInh4MeanA.las5;
lat_6 = LatencyInh4MeanA.pel6;las_6 = LatencyInh4MeanA.las6;


figure(f6);
col = "#77AC30";
    subplot(3,2,2);boxchart(las_4,lat_4,BoxFaceColor=col);hold on;
        %swarmchart(las_4,lat_4,".");
        ylabel('latency to 4th pellet');
        xticks([0 1]);xticklabels({'control', 'laser on'});
        ylim([0 30]);
        ctrl_mean = mean(lat_4(find(las_4 == 0)));
        las_mean = mean(lat_4(find(las_4 == 1)));
        plot([0 1],[ctrl_mean las_mean], 'k*');%xy

    subplot(3,2,4);boxchart(las_5,lat_5,BoxFaceColor=col);hold on;
        %swarmchart(las_5,lat_5,".");
        ylabel('latency to 5th pellet');
        xticks([0 1]);xticklabels({'control', 'laser on'});
        ylim([0 30]);
        ctrl_mean = mean(lat_5(find(las_5 == 0)));
        las_mean = mean(lat_5(find(las_5 == 1)));
        plot([0 1],[ctrl_mean las_mean], 'k*');%xy

    subplot(3,2,6);boxchart(las_6,lat_6,BoxFaceColor=col);hold on;
        %swarmchart(las_6,lat_6,".");
        ylabel('latency to 6th pellet');
        xticks([0 1]);xticklabels({'control', 'laser on'});
        ylim([0 30]);
        ctrl_mean = mean(lat_6(find(las_6 == 0)));
        las_mean = mean(lat_6(find(las_6 == 1)));
        plot([0 1],[ctrl_mean las_mean], 'k*');%xy
